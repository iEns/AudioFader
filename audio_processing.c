/*
 * AudioFader - Audio Processing Module
 *
 * Functions for trimming, fading, and padding audio data.
 *
 * Copyright (c) 2020-2026 Jens Sandalgaard
 * Licensed under MIT License
 */

#include "audio_processing.h"

/* ============================================================================
 * Sample I/O Functions
 * ============================================================================ */

int32_t read_sample(const unsigned char *buffer, int index,
                    int bytes_per_sample, int buffer_size,
                    struct int24 *buf24) {
    (void)buf24; /* kept for API compatibility; 24-bit is now handled portably */
    /* FIX #5: Add buffer bounds checking */
    if (buffer == NULL || index < 0) {
        return 0;
    }

    /* Calculate byte offset and check bounds */
    int byte_offset = index * bytes_per_sample;
    if (byte_offset < 0 || byte_offset + bytes_per_sample > buffer_size) {
        return 0;  /* Out of bounds */
    }

    const unsigned char *ptr = buffer + byte_offset;
    int32_t sample = 0;

    switch (bytes_per_sample) {
        case 1:
            sample = *ptr;
            break;
        case 2:
            sample = (int16_t)(ptr[1] << 8 | ptr[0]);
            break;
        case 3: {
            /* Portable 24-bit little-endian sign extension (no GCC bitfield). */
            int32_t v = (int32_t)((uint32_t)ptr[0] |
                                  ((uint32_t)ptr[1] << 8) |
                                  ((uint32_t)ptr[2] << 16));
            if (v & 0x00800000) {
                v |= (int32_t)(~0x00FFFFFF);
            }
            sample = v;
            break;
        }
        case 4:
            memcpy(&sample, ptr, 4);
            break;
        default:
            return 0;
    }

    return sample;
}

int write_sample(unsigned char *buffer, int index, int bytes_per_sample,
                 int buffer_size, int32_t value, struct int24 *buf24) {
    (void)buf24; /* kept for API compatibility; 24-bit is now handled portably */
    /* FIX #5: Add buffer bounds checking */
    if (buffer == NULL || index < 0) {
        return 0;
    }

    /* Calculate byte offset and check bounds */
    int byte_offset = index * bytes_per_sample;
    if (byte_offset < 0 || byte_offset + bytes_per_sample > buffer_size) {
        return 0;  /* Out of bounds */
    }

    unsigned char *ptr = buffer + byte_offset;

    switch (bytes_per_sample) {
        case 1:
            *ptr = (unsigned char)value;
            break;
        case 2:
            ptr[0] = (unsigned char)(value & 0xFF);
            ptr[1] = (unsigned char)((value >> 8) & 0xFF);
            break;
        case 3:
            ptr[0] = (unsigned char)(value & 0xFF);
            ptr[1] = (unsigned char)((value >> 8) & 0xFF);
            ptr[2] = (unsigned char)((value >> 16) & 0xFF);
            break;
        case 4:
            memcpy(ptr, &value, 4);
            break;
        default:
            return 0;
    }

    return 1;  /* Success */
}

int32_t calculate_threshold(int bytes_per_sample, double threshold_percent) {
    int32_t max_value = 0;

    switch (bytes_per_sample) {
        case 1: max_value = INT8_MAX_SAMPLE; break;
        case 2: max_value = INT16_MAX_SAMPLE; break;
        case 3: max_value = INT24_MAX_SAMPLE; break;
        case 4: max_value = INT32_MAX_SAMPLE; break;
    }

    return (int32_t)((double)max_value * threshold_percent / 100.0);
}

int duration_to_samples(int duration_ms, int sample_rate) {
    /* Use 64-bit arithmetic to detect overflow */
    int64_t samples = ((int64_t)duration_ms * sample_rate) / 1000;

    if (samples > INT_MAX) {
        LOG_ERROR("Error: Duration too long for sample rate (%d ms at %d Hz)\n",
                  duration_ms, sample_rate);
        return -1;
    }

    return (int)samples;
}

/* ============================================================================
 * Trim Boundary Detection
 * ============================================================================ */

/**
 * Magnitude of a sample for silence detection.
 *
 * Signed formats (16/24/32-bit): absolute value (INT32_MIN-safe).
 * Unsigned 8-bit PCM: distance from silence center (128).
 */
static int32_t sample_magnitude(int32_t sample, int bytes_per_sample) {
    if (bytes_per_sample == 1) {
        int32_t centered = sample - 128;
        return centered < 0 ? -centered : centered;
    }
    if (sample == INT32_MIN) {
        return INT32_MAX;
    }
    return sample < 0 ? -sample : sample;
}

/**
 * Round a sample index down to its channel-frame start.
 */
static int align_frame_start(int index, int num_channels) {
    return index - (index % num_channels);
}

/**
 * Round a sample index up to the exclusive end of its channel frame.
 * Result is clamped to total_samples by the caller.
 */
static int align_frame_end_exclusive(int index, int num_channels, int total_samples) {
    int end = align_frame_start(index, num_channels) + num_channels;
    return end > total_samples ? total_samples : end;
}

int find_first_above_threshold(const unsigned char *buffer, int buffer_size,
                               int total_samples, int bytes_per_sample,
                               int32_t threshold, int sample_rate,
                               int num_channels, struct int24 *buf24) {
    /* For files < 10 seconds, linear search is fast enough */
    int samples_in_threshold = sample_rate * TRIM_OPTIMIZATION_THRESHOLD_SECONDS * num_channels;
    if (total_samples < samples_in_threshold) {
        /* Linear search */
        for (int i = 0; i < total_samples; i++) {
            int32_t sample = read_sample(buffer, i, bytes_per_sample, buffer_size, buf24);
            if (sample_magnitude(sample, bytes_per_sample) > threshold) {
                /* Align to channel boundary */
                return align_frame_start(i, num_channels);
            }
        }
        return 0;
    }

    /* For longer files, use sampling strategy (optimization for large files) */
    /* Check every Nth sample (TRIM_SAMPLING_DIVISOR samples per second) */
    int sample_stride = sample_rate * num_channels / TRIM_SAMPLING_DIVISOR;
    int approx_start = 0;

    for (int i = 0; i < total_samples; i += sample_stride) {
        int32_t sample = read_sample(buffer, i, bytes_per_sample, buffer_size, buf24);
        if (sample_magnitude(sample, bytes_per_sample) > threshold) {
            approx_start = i;
            break;
        }
    }

    /* Narrow search: check samples in window [approx_start - stride, approx_start + stride] */
    int search_start = (approx_start > sample_stride) ? (approx_start - sample_stride) : 0;
    int search_end = (approx_start + sample_stride < total_samples) ?
                     (approx_start + sample_stride) : total_samples;

    for (int i = search_start; i < search_end; i++) {
        int32_t sample = read_sample(buffer, i, bytes_per_sample, buffer_size, buf24);
        if (sample_magnitude(sample, bytes_per_sample) > threshold) {
            /* Align to channel boundary */
            return align_frame_start(i, num_channels);
        }
    }

    return 0;
}

int find_last_above_threshold(const unsigned char *buffer, int buffer_size,
                              int total_samples, int bytes_per_sample,
                              int32_t threshold, int sample_rate,
                              int num_channels, struct int24 *buf24) {
    /* For files < 10 seconds, linear search is fast enough */
    int samples_in_threshold = sample_rate * TRIM_OPTIMIZATION_THRESHOLD_SECONDS * num_channels;
    if (total_samples < samples_in_threshold) {
        /* Linear search from end */
        for (int i = total_samples - 1; i >= 0; i--) {
            int32_t sample = read_sample(buffer, i, bytes_per_sample, buffer_size, buf24);
            if (sample_magnitude(sample, bytes_per_sample) > threshold) {
                /* Exclusive end of the frame containing the last loud sample */
                return align_frame_end_exclusive(i, num_channels, total_samples);
            }
        }
        return total_samples;
    }

    /* For longer files, use sampling strategy */
    int sample_stride = sample_rate * num_channels / TRIM_SAMPLING_DIVISOR;
    int approx_end = total_samples;

    for (int i = total_samples - 1; i >= 0; i -= sample_stride) {
        int32_t sample = read_sample(buffer, i, bytes_per_sample, buffer_size, buf24);
        if (sample_magnitude(sample, bytes_per_sample) > threshold) {
            approx_end = i;
            break;
        }
    }

    /* Narrow search */
    int search_start = (approx_end > sample_stride) ? (approx_end - sample_stride) : 0;
    int search_end = (approx_end + sample_stride < total_samples) ?
                     (approx_end + sample_stride) : total_samples;

    for (int i = search_end - 1; i >= search_start; i--) {
        int32_t sample = read_sample(buffer, i, bytes_per_sample, buffer_size, buf24);
        if (sample_magnitude(sample, bytes_per_sample) > threshold) {
            /* Exclusive end of the frame containing the last loud sample */
            return align_frame_end_exclusive(i, num_channels, total_samples);
        }
    }

    return total_samples;
}

int find_trim_boundaries(const wav_header_t *header,
                         const unsigned char *buffer, int buffer_size,
                         double trim_start_percent, double trim_end_percent,
                         int *first_sample_index, int *last_sample_index,
                         struct int24 *buf24) {
    int bytes_per_sample = header->bits_per_sample / 8;
    int total_samples = buffer_size / bytes_per_sample;

    *first_sample_index = 0;
    *last_sample_index = total_samples;

    /* Find first sample above threshold */
    if (trim_start_percent > 0) {
        int32_t threshold = calculate_threshold(bytes_per_sample, trim_start_percent);
        *first_sample_index = find_first_above_threshold(buffer, buffer_size, total_samples,
                                                         bytes_per_sample, threshold,
                                                         header->sample_rate,
                                                         header->num_channels, buf24);
    }

    /* Find last sample above threshold */
    if (trim_end_percent > 0) {
        int32_t threshold = calculate_threshold(bytes_per_sample, trim_end_percent);
        *last_sample_index = find_last_above_threshold(buffer, buffer_size, total_samples,
                                                       bytes_per_sample, threshold,
                                                       header->sample_rate,
                                                       header->num_channels, buf24);
    }

    return 0;
}

/* ============================================================================
 * Fade Functions
 * ============================================================================ */

double *create_fade_lut(int fade_samples, int curve_type, int direction) {
    /* FIX #3: Handle edge case where fade_samples <= 1 */
    if (fade_samples <= 1) {
        return NULL;  /* No fade needed for 0 or 1 samples */
    }

    double *lut = (double *)malloc(fade_samples * sizeof(double));
    if (lut == NULL) {
        return NULL;
    }

    for (int i = 0; i < fade_samples; i++) {
        double fade_index = (double)i / (double)(fade_samples - 1);

        /* Invert for fade-out */
        if (direction == FADE_DIRECTION_OUT) {
            fade_index = 1.0 - fade_index;
        }

        /* Calculate curve value */
        switch (curve_type) {
            case FADE_CURVE_LINEAR:
                lut[i] = fade_index;
                break;
            case FADE_CURVE_BEZIER:
                lut[i] = fade_index * fade_index * (3.0 - 2.0 * fade_index);
                break;
            case FADE_CURVE_LOGARITHMIC:
                /* Simplified from fade_index / (1.0 + (1.0 - fade_index)) */
                lut[i] = fade_index / (2.0 - fade_index);
                break;
            default:
                lut[i] = fade_index;
                break;
        }
    }

    return lut;
}

void apply_fade(const wav_header_t *header, unsigned char *buffer,
                int buffer_size, int fade_duration, int direction,
                int curve_type, struct int24 *buf24) {
    int fade_samples = duration_to_samples(fade_duration, header->sample_rate);
    if (fade_samples < 0) {
        /* Overflow detected in duration_to_samples */
        LOG_WARN("Fade duration calculation failed, skipping fade\n");
        return;
    }

    int bytes_per_sample = header->bits_per_sample / 8;
    int total_samples = buffer_size / bytes_per_sample;

    /* Validate fade duration doesn't exceed audio length */
    int fade_samples_with_channels = fade_samples * header->num_channels;
    if (fade_samples_with_channels > total_samples) {
        LOG_WARN("Fade duration (%d ms) exceeds audio length, clamping to fit\n",
                 fade_duration);
        fade_samples_with_channels = total_samples;
        fade_samples = fade_samples_with_channels / header->num_channels;
    }

    /* Pre-calculate fade curve lookup table - major performance optimization */
    double *fade_lut = create_fade_lut(fade_samples, curve_type, direction);
    if (fade_lut == NULL) {
        /* FIX #3: Handle edge case gracefully - no fade needed or allocation failed */
        if (fade_samples > 1) {
            LOG_ERROR("Error: Cannot allocate fade lookup table, skipping fade\n");
        }
        return;
    }

    int fade_start, fade_end;
    if (direction == FADE_DIRECTION_IN) {
        fade_start = 0;
        fade_end = fade_samples_with_channels;
    } else {
        fade_start = total_samples - fade_samples_with_channels;
        if (fade_start < 0) {
            fade_start = 0;  /* Clamp to buffer start */
        }
        fade_end = total_samples;
    }

    /* Apply fade using pre-calculated lookup table - no per-sample math! */
    for (int i = fade_start; i < fade_end; i++) {
        /* Look up pre-calculated fade factor from table */
        int lut_index = (i - fade_start) / header->num_channels;
        if (lut_index >= fade_samples) {
            lut_index = fade_samples - 1;  /* Safety clamp */
        }
        double fade_factor = fade_lut[lut_index];

        /* Read, apply fade, and write back */
        int32_t original_sample = read_sample(buffer, i, bytes_per_sample, buffer_size, buf24);

        /* Handle 8-bit unsigned offset (silence = 128) */
        if (bytes_per_sample == 1) {
            original_sample -= 128;
            int32_t faded_sample = (int32_t)(original_sample * fade_factor);
            write_sample(buffer, i, bytes_per_sample, buffer_size, faded_sample + 128, buf24);
        } else {
            int32_t faded_sample = (int32_t)(original_sample * fade_factor);
            write_sample(buffer, i, bytes_per_sample, buffer_size, faded_sample, buf24);
        }
    }

    free(fade_lut);
}

/* ============================================================================
 * Padding Functions
 * ============================================================================ */

unsigned char *create_padded_buffer(const unsigned char *input, int input_size,
                                    const wav_header_t *header, int pad_start_ms,
                                    int pad_end_ms, int *output_size) {
    int pad_start_samples = duration_to_samples(pad_start_ms, header->sample_rate);
    if (pad_start_samples < 0) {
        LOG_ERROR("Error: Padding start duration calculation failed (overflow at %d ms)\n",
                  pad_start_ms);
        return NULL;
    }
    int pad_start_bytes = pad_start_samples * header->num_channels *
                          (header->bits_per_sample / 8);

    int pad_end_samples = duration_to_samples(pad_end_ms, header->sample_rate);
    if (pad_end_samples < 0) {
        LOG_ERROR("Error: Padding end duration calculation failed (overflow at %d ms)\n",
                  pad_end_ms);
        return NULL;
    }
    int pad_end_bytes = pad_end_samples * header->num_channels *
                        (header->bits_per_sample / 8);

    /* Check for integer overflow in output size calculation */
    int temp_size;
    if (!safe_add_int(input_size, pad_start_bytes, &temp_size) ||
        !safe_add_int(temp_size, pad_end_bytes, output_size)) {
        LOG_ERROR("Error: Total padding (%d + %d ms) would cause buffer size overflow\n",
                  pad_start_ms, pad_end_ms);
        return NULL;
    }

    /* Validate input_size > 0 before memcpy */
    if (input_size <= 0) {
        LOG_ERROR("Error: Cannot pad zero-length audio data\n");
        return NULL;
    }

    /* Allocate buffer initialized to zero (silence) */
    unsigned char *padded_buffer = (unsigned char *)calloc(1, *output_size);
    if (padded_buffer == NULL) {
        LOG_ERROR("Error: Cannot allocate %d bytes for padded buffer (out of memory?)\n",
                  *output_size);
        return NULL;
    }

    /* Copy input data into middle of padded buffer */
    memcpy(padded_buffer + pad_start_bytes, input, input_size);

    return padded_buffer;
}

/* ============================================================================
 * Main Processing Pipeline
 * ============================================================================ */

int validate_options_with_audio(const options_t *options, const wav_header_t *header) {
    int sample_rate = header->sample_rate;
    int total_samples = header->subchunk2_size / (header->bits_per_sample / 8) / header->num_channels;
    int duration_ms = (total_samples * 1000) / sample_rate;

    LOG_VERBOSE("Audio duration: %d ms (%d samples at %d Hz)\n",
                duration_ms, total_samples, sample_rate);

    /* Check fade durations don't exceed audio length */
    if (options->fade_in + options->fade_out > duration_ms) {
        LOG_ERROR("Error: Combined fade duration (%d ms) exceeds audio length (%d ms)\n",
                  options->fade_in + options->fade_out, duration_ms);
        return 1;
    }

    /* Warn about extreme trimming */
    if (options->trim_start >= 50 || options->trim_end >= 50) {
        LOG_WARN("Trim threshold %.2f%% is very high and may remove most/all audio\n",
                 (options->trim_start > options->trim_end) ? options->trim_start : options->trim_end);
    }

    /* Check padding doesn't exceed reasonable limits */
    int total_pad_ms = options->pad_start + options->pad_end;
    if (total_pad_ms > MAX_DURATION_MS / 2) {
        LOG_WARN("Total padding (%d ms) is very large\n", total_pad_ms);
    }

    int final_duration = duration_ms + total_pad_ms;
    if (final_duration > MAX_DURATION_MS) {
        LOG_ERROR("Error: Output duration (%d ms) would exceed maximum (%d ms)\n",
                  final_duration, MAX_DURATION_MS);
        return 1;
    }

    return 0;
}

int process_audio(audio_fader_context_t *ctx, const wav_header_t *header,
                  const unsigned char *input_data, int input_size,
                  unsigned char **output_data, int *output_size) {
    unsigned char *current_buffer = NULL;
    unsigned char *working_buffer = NULL;
    int current_size = input_size;
    int bytes_per_sample = header->bits_per_sample / 8;
    int total_samples = input_size / bytes_per_sample;

    /* Phase 1: Trim (only allocate if actually trimming) */
    int first_sample = 0, last_sample = 0;
    int needs_trim = (ctx->options.trim_start > 0 || ctx->options.trim_end > 0);

    if (needs_trim) {
        if (find_trim_boundaries(header, input_data, input_size,
                                ctx->options.trim_start, ctx->options.trim_end,
                                &first_sample, &last_sample,
                                &ctx->sample24_buffer) != 0) {
            return 1;
        }

        /* Default to whole buffer if no trimming */
        if (last_sample == 0) {
            last_sample = total_samples;
        }

        /* Validate trim didn't remove everything */
        int trimmed_samples = last_sample - first_sample;
        if (trimmed_samples <= 0) {
            LOG_ERROR("Error: Trimming removed all audio (threshold too high?)\n");
            LOG_ERROR("       Trim thresholds: start=%.2f%%, end=%.2f%%\n",
                      ctx->options.trim_start, ctx->options.trim_end);
            LOG_ERROR("       Detected range: sample %d to %d (of %d total)\n",
                      first_sample, last_sample, total_samples);
            return 1;
        }

        /* FIX #8: Check minimum audio length (at least one channel frame) */
        if (trimmed_samples < header->num_channels) {
            LOG_ERROR("Error: Trimmed audio (%d samples) is less than one channel frame\n",
                      trimmed_samples);
            return 1;
        }

        /* Warn if very little audio remains */
        if (trimmed_samples < 100) {
            LOG_WARN("Only %d samples remain after trimming (from %d)\n",
                     trimmed_samples, total_samples);
        }

        /* Only allocate and copy if we're actually trimming */
        if (first_sample != 0 || last_sample != total_samples) {
            int trimmed_size = trimmed_samples * bytes_per_sample;
            working_buffer = (unsigned char *)malloc(trimmed_size);
            if (working_buffer == NULL) {
                LOG_ERROR("Error: Cannot allocate %d bytes for trimmed buffer (out of memory?)\n",
                          trimmed_size);
                return 1;
            }
            memcpy(working_buffer, input_data + (first_sample * bytes_per_sample), trimmed_size);
            current_buffer = working_buffer;
            current_size = trimmed_size;
        } else {
            /* No actual trimming needed, use input buffer directly */
            current_buffer = (unsigned char *)input_data;
            current_size = input_size;
        }
    } else {
        /* No trimming at all, use input buffer directly */
        current_buffer = (unsigned char *)input_data;
        current_size = input_size;
    }

    /* Phase 2: Fade (done in-place, no new allocation) */
    if (ctx->options.fade_in > 0) {
        apply_fade(header, current_buffer, current_size, ctx->options.fade_in,
                   FADE_DIRECTION_IN, ctx->options.fade_curve, &ctx->sample24_buffer);
    }

    if (ctx->options.fade_out > 0) {
        apply_fade(header, current_buffer, current_size, ctx->options.fade_out,
                   FADE_DIRECTION_OUT, ctx->options.fade_curve, &ctx->sample24_buffer);
    }

    /* Phase 3: Pad (only allocate if actually padding) */
    int needs_pad = (ctx->options.pad_start > 0 || ctx->options.pad_end > 0);

    if (needs_pad) {
        int padded_size;
        unsigned char *padded_buffer = create_padded_buffer(current_buffer, current_size, header,
                                                            ctx->options.pad_start, ctx->options.pad_end,
                                                            &padded_size);
        if (padded_buffer == NULL) {
            if (working_buffer != NULL) {
                free(working_buffer);
            }
            return 1;
        }

        /* Free the working buffer if we allocated one */
        if (working_buffer != NULL) {
            free(working_buffer);
        }

        *output_data = padded_buffer;
        *output_size = padded_size;
    } else {
        /* No padding needed */
        if (working_buffer != NULL) {
            /* We allocated a trimmed buffer, use it */
            *output_data = working_buffer;
            *output_size = current_size;
        } else {
            /* No trim or pad - need to copy input to output */
            unsigned char *output_copy = (unsigned char *)malloc(current_size);
            if (output_copy == NULL) {
                LOG_ERROR("Error: Cannot allocate %d bytes for output buffer (out of memory?)\n",
                          current_size);
                return 1;
            }
            memcpy(output_copy, current_buffer, current_size);
            *output_data = output_copy;
            *output_size = current_size;
        }
    }

    return 0;
}
