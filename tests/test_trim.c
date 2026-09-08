/*
 * AudioFader - unit tests (trim, sample I/O, fade, pad, WAV validation)
 *
 * Build: make check   (no external dependencies, C99 only)
 */
#include <stdio.h>
#include <string.h>

#include "../common.h"
#include "../audio_processing.h"
#include "../options.h"
#include "../wav_io.h"

static int failures = 0;
static int passes = 0;

#define CHECK(cond, ...) do { \
    if (cond) { passes++; } \
    else { \
        failures++; \
        fprintf(stderr, "FAIL %s:%d: ", __FILE__, __LINE__); \
        fprintf(stderr, __VA_ARGS__); \
        fprintf(stderr, "\n"); \
    } \
} while (0)

/* Write helpers (little-endian) for building fixtures */
static void put16le(unsigned char *p, int16_t v) {
    p[0] = (unsigned char)(v & 0xFF);
    p[1] = (unsigned char)((v >> 8) & 0xFF);
}

static void put24le(unsigned char *p, int32_t v) {
    p[0] = (unsigned char)(v & 0xFF);
    p[1] = (unsigned char)((v >> 8) & 0xFF);
    p[2] = (unsigned char)((v >> 16) & 0xFF);
}

static void put32le(unsigned char *p, int32_t v) {
    memcpy(p, &v, 4);
}

static void test_sample_io_roundtrip(void) {
    unsigned char buf[16];
    struct int24 b24;
    memset(buf, 0, sizeof(buf));

    /* 16-bit incl. negative */
    write_sample(buf, 0, 2, sizeof(buf), -10000, &b24);
    write_sample(buf, 1, 2, sizeof(buf), 10000, &b24);
    CHECK(read_sample(buf, 0, 2, sizeof(buf), &b24) == -10000, "16-bit negative roundtrip");
    CHECK(read_sample(buf, 1, 2, sizeof(buf), &b24) == 10000, "16-bit positive roundtrip");

    /* 24-bit sign extension: raw bytes, no bitfield */
    {
        unsigned char raw[3] = {0x00, 0x00, 0x80}; /* 0x800000 -> -8388608 */
        CHECK(read_sample(raw, 0, 3, 3, NULL) == -8388608, "24-bit min sign extension got %d",
              (int)read_sample(raw, 0, 3, 3, NULL));
        unsigned char raw2[3] = {0xFF, 0xFF, 0x7F}; /* 0x7FFFFF -> max */
        CHECK(read_sample(raw2, 0, 3, 3, NULL) == 8388607, "24-bit max sign extension");
        unsigned char raw3[3] = {0xFF, 0xFF, 0xFF}; /* -1 */
        CHECK(read_sample(raw3, 0, 3, 3, NULL) == -1, "24-bit -1 sign extension");
        /* write path must be portable too */
        unsigned char out[3];
        write_sample(out, 0, 3, 3, -8388608, NULL);
        CHECK(out[0] == 0x00 && out[1] == 0x00 && out[2] == 0x80, "24-bit min write bytes");
        write_sample(out, 0, 3, 3, -1, NULL);
        CHECK(out[0] == 0xFF && out[1] == 0xFF && out[2] == 0xFF, "24-bit -1 write bytes");
    }

    /* 32-bit incl. INT32_MIN (magnitude helper must not overflow) */
    {
        unsigned char b32[8];
        put32le(b32, INT32_MIN);
        put32le(b32 + 4, 123456);
        CHECK(read_sample(b32, 0, 4, 8, NULL) == INT32_MIN, "32-bit min roundtrip");
        CHECK(read_sample(b32, 1, 4, 8, NULL) == 123456, "32-bit pos roundtrip");
    }

    /* 8-bit is unsigned */
    {
        unsigned char b8[2] = {128, 255};
        CHECK(read_sample(b8, 0, 1, 2, NULL) == 128, "8-bit silence raw value");
        CHECK(read_sample(b8, 1, 1, 2, NULL) == 255, "8-bit max raw value");
    }

    /* NULL buf24 must be tolerated (portable API) */
    CHECK(read_sample(buf, 0, 2, sizeof(buf), NULL) == -10000, "read_sample NULL buf24");
}

static void test_negative_peak_stops_trim(void) {
    /* 16-bit mono: silence, then a LOUD NEGATIVE peak, then silence */
    int16_t vals[] = {0, 0, -12000, 0, 0, 0};
    unsigned char buf[sizeof(vals)];
    for (size_t i = 0; i < sizeof(vals) / sizeof(vals[0]); i++) {
        put16le(buf + i * 2, vals[i]);
    }
    int32_t threshold = calculate_threshold(2, 1.0); /* ~327 */
    int first = find_first_above_threshold(buf, sizeof(buf), 6, 2, threshold, 8000, 1, NULL);
    int last = find_last_above_threshold(buf, sizeof(buf), 6, 2, threshold, 8000, 1, NULL);
    CHECK(first == 2, "negative peak: first==2 got %d", first);
    CHECK(last == 3, "negative peak: last exclusive==3 got %d", last);

    /* Same for 24-bit and 32-bit */
    {
        unsigned char b24[6 * 3];
        memset(b24, 0, sizeof(b24));
        put24le(b24 + 2 * 3, -500000);
        int32_t th24 = calculate_threshold(3, 1.0);
        CHECK(find_first_above_threshold(b24, sizeof(b24), 6, 3, th24, 8000, 1, NULL) == 2,
              "24-bit negative peak detected");
        CHECK(find_last_above_threshold(b24, sizeof(b24), 6, 3, th24, 8000, 1, NULL) == 3,
              "24-bit negative peak end exclusive");
    }
    {
        unsigned char b32[6 * 4];
        memset(b32, 0, sizeof(b32));
        put32le(b32 + 2 * 4, -50000000);
        int32_t th32 = calculate_threshold(4, 0.01);
        CHECK(find_first_above_threshold(b32, sizeof(b32), 6, 4, th32, 8000, 1, NULL) == 2,
              "32-bit negative peak detected");
    }
}

static void test_8bit_centering(void) {
    /* Silence = 128. Old code compared raw 0..255 > threshold so it never trimmed. */
    unsigned char buf[] = {128, 128, 0, 128, 128}; /* index 2 is loud (-128) */
    int32_t threshold = calculate_threshold(1, 1.0); /* ~1 */
    int first = find_first_above_threshold(buf, sizeof(buf), 5, 1, threshold, 8000, 1, NULL);
    int last = find_last_above_threshold(buf, sizeof(buf), 5, 1, threshold, 8000, 1, NULL);
    CHECK(first == 2, "8-bit: first loud==2 got %d", first);
    CHECK(last == 3, "8-bit: last exclusive==3 got %d", last);

    /* All silence must not trim */
    unsigned char silence[] = {128, 128, 128, 128};
    CHECK(find_first_above_threshold(silence, sizeof(silence), 4, 1, threshold, 8000, 1, NULL) == 0,
          "8-bit all-silence first==0");
    CHECK(find_last_above_threshold(silence, sizeof(silence), 4, 1, threshold, 8000, 1, NULL) == 4,
          "8-bit all-silence last==total");

    /* 255 is +127 -> loud; 127/128/129 around silence -> quiet at 1% threshold */
    unsigned char edge[] = {127, 128, 129, 255};
    int f = find_first_above_threshold(edge, sizeof(edge), 4, 1, threshold, 8000, 1, NULL);
    CHECK(f == 3, "8-bit near-silence ignored, 255 loud: first==3 got %d", f);
}

static void test_stereo_frame_alignment(void) {
    /* Stereo frames (L,R): frame0 quiet, frame1 loud on R only, frame2 quiet */
    int16_t vals[] = {0, 0, 0, -9000, 0, 0};
    unsigned char buf[sizeof(vals)];
    for (size_t i = 0; i < 6; i++) {
        put16le(buf + i * 2, vals[i]);
    }
    int32_t th = calculate_threshold(2, 1.0);
    int first = find_first_above_threshold(buf, sizeof(buf), 6, 2, th, 8000, 2, NULL);
    int last = find_last_above_threshold(buf, sizeof(buf), 6, 2, th, 8000, 2, NULL);
    CHECK(first == 2, "stereo: first aligned to frame start 2, got %d", first);
    CHECK(last == 4, "stereo: last exclusive end of frame (2+2)=4, got %d", last);

    /* Loud sample in last frame must yield total (clamped), not drop the frame */
    int16_t tail[] = {0, 0, 0, 0, 0, 7000};
    unsigned char buf2[sizeof(tail)];
    for (size_t i = 0; i < 6; i++) {
        put16le(buf2 + i * 2, tail[i]);
    }
    int last2 = find_last_above_threshold(buf2, sizeof(buf2), 6, 2, th, 8000, 2, NULL);
    CHECK(last2 == 6, "stereo: loud tail frame end==total 6, got %d", last2);
}

static void test_long_file_sampling_path(void) {
    /* Exactness regression: single-sample peaks in a >10s file must be
     * found. The old coarse-stride scan missed these. */
    int total = 90000;
    int bps = 2;
    unsigned char *buf = (unsigned char *)calloc((size_t)total, (size_t)bps);
    if (!buf) {
        CHECK(0, "%s", "calloc failed for long-path test");
        return;
    }
    put16le(buf + (size_t)100 * 2, -15000);   /* single loud sample */
    put16le(buf + (size_t)80000 * 2, -15000); /* single loud sample */
    int32_t th = calculate_threshold(bps, 1.0);
    int first = find_first_above_threshold(buf, total * bps, total, bps, th, 8000, 1, NULL);
    int last = find_last_above_threshold(buf, total * bps, total, bps, th, 8000, 1, NULL);
    CHECK(first == 100, "long path: single-sample first==100, got %d", first);
    CHECK(last == 80001, "long path: single-sample last exclusive 80001, got %d", last);
    free(buf);
}

static void test_8bit_fade_preserves_silence(void) {
    /* Fade-in over 1ms at 8kHz = 8 samples: first frame factor ~0 must
     * stay 128 (unsigned silence), and loud audio must be attenuated. */
    wav_header_t hdr;
    memset(&hdr, 0, sizeof(hdr));
    hdr.sample_rate = 8000;
    hdr.num_channels = 1;
    hdr.bits_per_sample = 8;
    unsigned char buf8[8];
    memset(buf8, 128, sizeof(buf8));
    apply_fade(&hdr, buf8, sizeof(buf8), 1, FADE_DIRECTION_IN, FADE_CURVE_LINEAR, NULL);
    CHECK(buf8[0] == 128, "8-bit fade-in silence stays 128, got %d", buf8[0]);
    CHECK(buf8[7] == 128, "8-bit fade-in end silence stays 128, got %d", buf8[7]);

    /* Loud max (255 = +127) faded in: first sample ~0 gain -> ~128. */
    unsigned char loud[8];
    memset(loud, 255, sizeof(loud));
    apply_fade(&hdr, loud, sizeof(loud), 1, FADE_DIRECTION_IN, FADE_CURVE_LINEAR, NULL);
    CHECK(loud[0] == 128, "8-bit fade-in attenuates start to 128, got %d", loud[0]);
    CHECK(loud[7] == 255, "8-bit fade-in keeps end at 255, got %d", loud[7]);
}

static void test_trim_boundaries_integration(void) {
    /* End-to-end via find_trim_boundaries: leading/trailing silence removed, loud kept. */
    wav_header_t hdr;
    memset(&hdr, 0, sizeof(hdr));
    hdr.sample_rate = 8000;
    hdr.num_channels = 1;
    hdr.bits_per_sample = 16;
    int16_t vals[] = {0, 0, 5000, -6000, 0, 0};
    unsigned char buf[sizeof(vals)];
    for (size_t i = 0; i < 6; i++) {
        put16le(buf + i * 2, vals[i]);
    }
    int first = -1, last = -1;
    find_trim_boundaries(&hdr, buf, sizeof(buf), 1.0, 1.0, &first, &last, NULL);
    CHECK(first == 2, "boundaries: first==2 got %d", first);
    CHECK(last == 4, "boundaries: last exclusive==4 got %d", last);
}

static void test_sample_io_oob(void) {
    unsigned char buf[4] = {0, 0, 0, 0};
    /* Out-of-bounds reads return 0; writes return 0 and leave buffer alone. */
    CHECK(read_sample(NULL, 0, 2, 4, NULL) == 0, "read NULL buffer");
    CHECK(read_sample(buf, -1, 2, 4, NULL) == 0, "read negative index");
    CHECK(read_sample(buf, 2, 2, 4, NULL) == 0, "read past end");
    CHECK(read_sample(buf, 0, 5, 4, NULL) == 0, "read invalid bps");
    CHECK(write_sample(NULL, 0, 2, 4, 1, NULL) == 0, "write NULL buffer");
    CHECK(write_sample(buf, 2, 2, 4, 0x1234, NULL) == 0, "write past end");
    CHECK(write_sample(buf, 0, 5, 4, 1, NULL) == 0, "write invalid bps");
    CHECK(buf[0] == 0 && buf[1] == 0, "failed write leaves buffer intact");
}

static void test_threshold_and_duration(void) {
    CHECK(calculate_threshold(5, 1.0) == 0, "unknown bps threshold==0");
    CHECK(calculate_threshold(1, 1.0) == 1, "8-bit 1%% threshold==1, got %d",
          (int)calculate_threshold(1, 1.0));
    /* Overflow: 10M ms at 192kHz = 1.92B samples fits, but INT_MAX+1 must fail. */
    CHECK(duration_to_samples(10000000, 192000) == 1920000000,
          "duration max fits");
    CHECK(duration_to_samples(INT_MAX, 192000) < 0, "duration overflow returns -1");
    CHECK(create_fade_lut(0, 0, 0) == NULL, "fade lut 0 samples NULL");
    CHECK(create_fade_lut(1, 0, 0) == NULL, "fade lut 1 sample NULL");
    {
        double *lut = create_fade_lut(8, FADE_CURVE_LINEAR, FADE_DIRECTION_IN);
        CHECK(lut != NULL, "fade lut allocated");
        if (lut != NULL) {
            CHECK(lut[0] == 0.0, "fade-in starts at 0");
            CHECK(lut[7] == 1.0, "fade-in ends at 1");
            free(lut);
        }
    }
    {
        double *lut = create_fade_lut(8, FADE_CURVE_LINEAR, FADE_DIRECTION_OUT);
        CHECK(lut != NULL, "fade-out lut allocated");
        if (lut != NULL) {
            CHECK(lut[0] == 1.0, "fade-out starts at 1");
            CHECK(lut[7] == 0.0, "fade-out ends at 0");
            free(lut);
        }
    }
}

static void test_header_rejects_nonpcm(void) {
    wav_header_t hdr;
    memset(&hdr, 0, sizeof(hdr));
    memcpy(hdr.chunk_id, "RIFF", 4);
    memcpy(hdr.format, "WAVE", 4);
    memcpy(hdr.subchunk1_id, "fmt ", 4);
    hdr.subchunk1_size = 16;
    hdr.audio_format = 1;
    hdr.num_channels = 1;
    hdr.sample_rate = 8000;
    hdr.bits_per_sample = 16;
    hdr.block_align = 2;
    hdr.byte_rate = 16000;
    CHECK(validate_wav_header(&hdr) == 0, "valid PCM header accepted");

    hdr.audio_format = 3; /* IEEE float */
    CHECK(validate_wav_header(&hdr) != 0, "float format rejected");
    hdr.audio_format = 1;

    hdr.subchunk1_size = 18; /* extended fmt */
    CHECK(validate_wav_header(&hdr) != 0, "extended fmt rejected");
    hdr.subchunk1_size = 16;

    hdr.block_align = 99;
    CHECK(validate_wav_header(&hdr) != 0, "bad block_align rejected");
    hdr.block_align = 2;

    hdr.byte_rate = 1;
    CHECK(validate_wav_header(&hdr) != 0, "bad byte_rate rejected");
}

static void test_pad_edge_cases(void) {
    wav_header_t hdr;
    memset(&hdr, 0, sizeof(hdr));
    hdr.sample_rate = 192000;
    hdr.num_channels = 2;
    hdr.bits_per_sample = 32;
    {
        unsigned char in[8] = {1, 2, 3, 4, 5, 6, 7, 8};
        int out_size = 0;
        /* Zero-length input is rejected, not memcpy'd. */
        CHECK(create_padded_buffer(in, 0, &hdr, 10, 10, &out_size) == NULL,
              "pad zero-length input NULL");
        /* 10M ms at 192kHz stereo 32-bit = ~15 GB -> overflow guard NULL. */
        CHECK(create_padded_buffer(in, sizeof(in), &hdr, 10000000, 0,
                                   &out_size) == NULL,
              "pad huge start NULL (overflow guard)");
    }
}

static void test_same_file_rejected(void) {
    /* Identical argv strings are always the same file. */
    audio_fader_context_t ctx;
    memset(&ctx, 0, sizeof(ctx));
    char *argv1[] = {"prog", "same.wav", "same.wav"};
    CHECK(parse_options(3, argv1, &ctx) != 0, "identical in/out rejected");
}

static void test_write_output_exists_guard(void) {
    const char *path = "tests/tmp_exists_guard.wav";
    FILE *f = fopen(path, "wb");
    if (f != NULL) {
        fwrite("x", 1, 1, f);
        fclose(f);
    }
    audio_fader_context_t ctx;
    memset(&ctx, 0, sizeof(ctx));
    ctx.options.output_filename = path;
    ctx.options.force_overwrite = 0;
    wav_header_t hdr;
    memset(&hdr, 0, sizeof(hdr));
    unsigned char data[4] = {0, 0, 0, 0};
    CHECK(write_output(&ctx, &hdr, data, sizeof(data)) != 0,
          "existing output without --force rejected");
    remove(path);
}

static void write_test_wav(const char *path, const unsigned char *audio, int audio_size) {
    wav_header_t hdr;
    memset(&hdr, 0, sizeof(hdr));
    memcpy(hdr.chunk_id, "RIFF", 4);
    memcpy(hdr.format, "WAVE", 4);
    memcpy(hdr.subchunk1_id, "fmt ", 4);
    hdr.subchunk1_size = 16;
    hdr.audio_format = 1;
    hdr.num_channels = 1;
    hdr.sample_rate = 8000;
    hdr.bits_per_sample = 16;
    hdr.block_align = 2;
    hdr.byte_rate = 16000;
    memcpy(hdr.subchunk2_id, "data", 4);
    hdr.subchunk2_size = audio_size;
    hdr.chunk_size = audio_size + 36;
    FILE *f = fopen(path, "wb");
    if (!f) {
        return;
    }
    fwrite(&hdr, sizeof(hdr), 1, f);
    fwrite(audio, 1, (size_t)audio_size, f);
    fclose(f);
}

static void test_wav_load_data_offset(void) {
    /* Regression test: loader must return audio bytes, not header bytes.
     * Previously fread() ran from file offset 0 for standard WAVs, so the
     * "RIFF" header appeared as loud samples and broke trimming. */
    const char *path = "tests/tmp_io_test.wav";
    unsigned char audio[16];
    for (int i = 0; i < 8; i++) {
        put16le(audio + (size_t)i * 2, (int16_t)(i * 1000 - 4000));
    }
    write_test_wav(path, audio, sizeof(audio));

    audio_fader_context_t ctx;
    memset(&ctx, 0, sizeof(ctx));
    ctx.options.input_filename = path;
    wav_header_t hdr;
    memset(&hdr, 0, sizeof(hdr));
    unsigned char *data = NULL;
    int size = 0;
    int rc = load_wav_file(&ctx, &hdr, &data, &size);
    CHECK(rc == 0, "load_wav_file succeeds, got %d", rc);
    if (rc == 0) {
        CHECK(size == (int)sizeof(audio), "loaded size %d == audio size %d", size, (int)sizeof(audio));
        CHECK(memcmp(data, audio, sizeof(audio)) == 0, "loaded bytes match audio payload (no header)");
        /* First sample must be -4000, not 'RI' (0x4952 = 18770) from "RIFF" */
        struct int24 b24;
        memset(&b24, 0, sizeof(b24));
        int32_t s0 = read_sample(data, 0, 2, size, &b24);
        CHECK(s0 == -4000, "first sample -4000, got %d", (int)s0);
        free(data);
    }
    remove(path);
}

int main(void) {
    const char *names[] = {
        "sample_io_roundtrip", "negative_peak", "8bit_centering",
        "stereo_alignment", "long_exact", "8bit_fade", "boundaries",
        "sample_io_oob", "threshold_duration", "header_nonpcm",
        "pad_edge", "same_file", "exists_guard", "wav_offset"
    };
    void (*tests[])(void) = {
        test_sample_io_roundtrip, test_negative_peak_stops_trim,
        test_8bit_centering, test_stereo_frame_alignment,
        test_long_file_sampling_path, test_8bit_fade_preserves_silence,
        test_trim_boundaries_integration, test_sample_io_oob,
        test_threshold_and_duration, test_header_rejects_nonpcm,
        test_pad_edge_cases, test_same_file_rejected,
        test_write_output_exists_guard, test_wav_load_data_offset
    };
    size_t n = sizeof(tests) / sizeof(tests[0]);
    g_verbosity = VERBOSITY_QUIET; /* silence progress/log output during tests */
    for (size_t i = 0; i < n; i++) {
        int before = failures;
        tests[i]();
        printf("%s: %s\n", names[i], failures == before ? "ok" : "FAILED");
    }

    printf("trim tests: %d passed, %d failed\n", passes, failures);
    return failures == 0 ? 0 : 1;
}
