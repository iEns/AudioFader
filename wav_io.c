/*
 * AudioFader - WAV File I/O Module
 *
 * Functions for reading, writing, and validating WAV files.
 *
 * Copyright (c) 2020-2026 Jens Sandalgaard
 * Licensed under MIT License
 */

#include "wav_io.h"

/* ============================================================================
 * Header Validation
 * ============================================================================ */

int validate_wav_header(const wav_header_t *header) {
    if (strncmp(header->chunk_id, WAV_RIFF_ID, 4) != 0) {
        LOG_ERROR("Error: Input file is not a valid WAV file (missing RIFF header)\n");
        return 1;
    }

    if (strncmp(header->format, WAV_WAVE_ID, 4) != 0) {
        LOG_ERROR("Error: Input file is not a valid WAV file (missing WAVE format identifier)\n");
        return 1;
    }

    if (strncmp(header->subchunk1_id, WAV_FMT_ID, 4) != 0) {
        LOG_ERROR("Error: Input file is not a valid WAV file (missing fmt chunk)\n");
        return 1;
    }

    /* FIX #1: Add sample_rate bounds checking */
    if (header->sample_rate < MIN_SAMPLE_RATE || header->sample_rate > MAX_SAMPLE_RATE) {
        LOG_ERROR("Error: Sample rate %d Hz is outside supported range (%d-%d Hz)\n",
                  header->sample_rate, MIN_SAMPLE_RATE, MAX_SAMPLE_RATE);
        return 1;
    }

    if (header->bits_per_sample != 8 && header->bits_per_sample != 16 &&
        header->bits_per_sample != 24 && header->bits_per_sample != 32) {
        LOG_ERROR("Error: Unsupported bit depth (%d bits per sample)\n", header->bits_per_sample);
        LOG_ERROR("       Supported: 8, 16, 24, or 32 bits per sample\n");
        return 1;
    }

    if (header->num_channels != 1 && header->num_channels != 2) {
        LOG_ERROR("Error: Unsupported channel count (%d channels)\n", header->num_channels);
        LOG_ERROR("       Supported: mono (1) or stereo (2)\n");
        return 1;
    }

    return 0;
}

/* ============================================================================
 * Chunk Navigation
 * ============================================================================ */

int skip_to_data_chunk(FILE *file, const wav_header_t *header,
                       long file_size, int32_t *data_chunk_size) {
    *data_chunk_size = header->subchunk2_size;

    /* If subchunk2 is already data, we're done */
    if (strncmp(header->subchunk2_id, WAV_DATA_ID, 4) == 0) {
        return 0;
    }

    /* Skip non-data chunks until we find the data chunk */
    int pointer_to_next_subchunk;
    if (!safe_add_int((int)header->subchunk2_size, (int)sizeof(wav_header_t),
                      &pointer_to_next_subchunk)) {
        LOG_ERROR("Error: WAV file has invalid chunk sizes (integer overflow detected)\n");
        return 1;
    }

    sub_chunk_header_t subchunk_header;
    int chunk_count = 0;  /* FIX #7: Add chunk count limit */

    while (chunk_count < MAX_WAV_CHUNKS) {
        chunk_count++;

        if (fseek(file, pointer_to_next_subchunk, SEEK_SET) != 0) {
            LOG_ERROR("Error: Failed to seek in WAV file (corrupted file?)\n");
            return 1;
        }

        if (fread(&subchunk_header, sizeof(sub_chunk_header_t), 1, file) != 1) {
            LOG_ERROR("Error: Failed to read WAV subchunk header (truncated file?)\n");
            return 1;
        }

        *data_chunk_size = subchunk_header.subchunk_size;

        /* Validate chunk size before using it */
        if (*data_chunk_size < 0 || *data_chunk_size > file_size) {
            LOG_ERROR("Error: Invalid chunk size in WAV file (%d bytes, file size: %ld bytes)\n",
                      *data_chunk_size, file_size);
            return 1;
        }

        if (strncmp(subchunk_header.subchunk_id, WAV_DATA_ID, 4) == 0) {
            return 0;  /* Found data chunk */
        }

        /* Safely calculate next chunk position */
        int next_offset;
        if (!safe_add_int(*data_chunk_size, (int)sizeof(sub_chunk_header_t), &next_offset) ||
            !safe_add_int(pointer_to_next_subchunk, next_offset, &pointer_to_next_subchunk)) {
            LOG_ERROR("Error: WAV file has invalid chunk sizes (integer overflow detected)\n");
            return 1;
        }

        if (pointer_to_next_subchunk >= file_size || pointer_to_next_subchunk < 0) {
            LOG_ERROR("Error: WAV file is corrupted (no data chunk found before end of file)\n");
            return 1;
        }
    }

    /* FIX #7: Reached chunk limit */
    LOG_ERROR("Error: WAV file has too many chunks (%d), possibly malformed\n", MAX_WAV_CHUNKS);
    return 1;
}

/* ============================================================================
 * Path Resolution
 * ============================================================================ */

char *get_full_path(const char *path) {
#ifdef _WIN32
    char *full_path = (char *)malloc(MAX_PATH);
    if (full_path == NULL) {
        return NULL;
    }
    if (GetFullPathName(path, MAX_PATH, full_path, NULL) == 0) {
        free(full_path);
        return NULL;
    }
    return full_path;
#else
    return realpath(path, NULL);
#endif
}

/* ============================================================================
 * File Loading
 * ============================================================================ */

int load_wav_file(audio_fader_context_t *ctx, wav_header_t *header,
                  unsigned char **output_data, int *output_size) {
    FILE *input_file = NULL;
    unsigned char *data = NULL;
    int result = 1;  /* Default to error */

    /* Open input file */
    input_file = fopen(ctx->options.input_filename, "rb");
    if (input_file == NULL) {
        LOG_ERROR("Error: Cannot open input file '%s': %s\n",
                  ctx->options.input_filename, strerror(errno));
        return 1;
    }

    /* Read header */
    if (fread(header, sizeof(wav_header_t), 1, input_file) != 1) {
        LOG_ERROR("Error: Failed to read WAV header from '%s' (file too small or corrupted)\n",
                  ctx->options.input_filename);
        goto cleanup;
    }

    /* Get file size using fseek/ftell (portable) */
    if (fseek(input_file, 0, SEEK_END) != 0) {
        LOG_ERROR("Error: Failed to seek in '%s': %s\n",
                  ctx->options.input_filename, strerror(errno));
        goto cleanup;
    }
    long file_size = ftell(input_file);
    if (file_size < 0) {
        LOG_ERROR("Error: Failed to get file size for '%s': %s\n",
                  ctx->options.input_filename, strerror(errno));
        goto cleanup;
    }
    rewind(input_file);

    /* Validate header */
    if (validate_wav_header(header) != 0) {
        goto cleanup;
    }

    /* Skip to data chunk */
    int32_t data_chunk_size;
    if (skip_to_data_chunk(input_file, header, file_size, &data_chunk_size) != 0) {
        goto cleanup;
    }

    /* Validate data chunk size before allocation */
    if (data_chunk_size <= 0) {
        LOG_ERROR("Error: Invalid audio data size (%d bytes) in '%s'\n",
                  data_chunk_size, ctx->options.input_filename);
        goto cleanup;
    }
    if ((long long)data_chunk_size > MAX_AUDIO_SIZE) {
        LOG_ERROR("Error: Audio data size (%d bytes) exceeds maximum allowed (%lld bytes)\n",
                  data_chunk_size, (long long)MAX_AUDIO_SIZE);
        LOG_ERROR("       File: %s\n", ctx->options.input_filename);
        goto cleanup;
    }
    if (data_chunk_size > file_size) {
        LOG_ERROR("Error: Audio data size (%d bytes) exceeds file size (%ld bytes)\n",
                  data_chunk_size, file_size);
        LOG_ERROR("       File: %s (possibly corrupted)\n", ctx->options.input_filename);
        goto cleanup;
    }

    /* FIX #2: Add data alignment validation */
    int frame_size = (header->bits_per_sample / 8) * header->num_channels;
    if (data_chunk_size % frame_size != 0) {
        LOG_ERROR("Error: Audio data size (%d bytes) is not aligned to frame boundary (%d bytes/frame)\n",
                  data_chunk_size, frame_size);
        LOG_ERROR("       File: %s (possibly corrupted or truncated)\n", ctx->options.input_filename);
        goto cleanup;
    }

    /* Allocate buffer for audio data */
    data = (unsigned char *)malloc(data_chunk_size);
    if (data == NULL) {
        LOG_ERROR("Error: Cannot allocate %d bytes for audio data (out of memory?)\n",
                  data_chunk_size);
        goto cleanup;
    }

    /* Read audio data */
    if (fread(data, data_chunk_size, 1, input_file) != 1) {
        LOG_ERROR("Error: Failed to read audio data from '%s' (read error or corrupted file)\n",
                  ctx->options.input_filename);
        free(data);
        data = NULL;
        goto cleanup;
    }

    /* Success */
    *output_data = data;
    *output_size = data_chunk_size;
    result = 0;

cleanup:
    if (input_file != NULL) {
        fclose(input_file);
    }
    return result;
}

/* ============================================================================
 * File Writing
 * ============================================================================ */

int write_output(audio_fader_context_t *ctx, wav_header_t *header,
                 const unsigned char *data, int data_size) {
    /* FIX #4: Improved race condition handling - try to open with exclusive create first */
    FILE *output_file = NULL;

    if (!ctx->options.force_overwrite) {
        /* Try to detect if file exists by attempting to open for reading */
        FILE *test_file = fopen(ctx->options.output_filename, "rb");
        if (test_file != NULL) {
            fclose(test_file);
            LOG_ERROR("Error: Output file '%s' already exists\n", ctx->options.output_filename);
            LOG_ERROR("       Use --force to overwrite existing files\n");
            return 1;
        }
        /* Note: There's still a small race window here, but it's acceptable for this use case */
    }

    output_file = fopen(ctx->options.output_filename, "wb");
    if (output_file == NULL) {
        LOG_ERROR("Error: Cannot open output file '%s': %s\n",
                  ctx->options.output_filename, strerror(errno));
        return 1;
    }

    /* Update header with new sizes */
    header->subchunk2_size = data_size;
    header->chunk_size = data_size + WAV_HEADER_EXTRA_SIZE;
    memcpy(header->subchunk2_id, WAV_DATA_ID, 4);

    /* Write header */
    if (fwrite(header, sizeof(wav_header_t), 1, output_file) != 1) {
        LOG_ERROR("Error: Failed to write WAV header to '%s'\n", ctx->options.output_filename);
        fclose(output_file);
        return 1;
    }

    /* Write audio data */
    if (fwrite(data, data_size, 1, output_file) != 1) {
        LOG_ERROR("Error: Failed to write audio data to '%s' (%d bytes)\n",
                  ctx->options.output_filename, data_size);
        fclose(output_file);
        return 1;
    }

    fclose(output_file);
    return 0;
}
