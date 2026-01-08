/*
 * AudioFader - Audio Processing Module
 *
 * Functions for trimming, fading, and padding audio data.
 *
 * Copyright (c) 2020-2026 Jens Sandalgaard
 * Licensed under MIT License
 */

#ifndef AUDIOFADER_AUDIO_PROCESSING_H
#define AUDIOFADER_AUDIO_PROCESSING_H

#include "common.h"

/**
 * Read a sample value from buffer at the specified index.
 *
 * @param buffer Audio data buffer
 * @param index Sample index (not byte index)
 * @param bytes_per_sample Number of bytes per sample (1, 2, 3, or 4)
 * @param buffer_size Total buffer size in bytes (for bounds checking)
 * @param buf24 Pointer to int24 buffer for 24-bit operations
 * @return Sample value as signed 32-bit integer, or 0 if out of bounds
 */
int32_t read_sample(const unsigned char *buffer, int index,
                    int bytes_per_sample, int buffer_size,
                    struct int24 *buf24);

/**
 * Write a sample value to buffer at the specified index.
 *
 * @param buffer Audio data buffer
 * @param index Sample index (not byte index)
 * @param bytes_per_sample Number of bytes per sample (1, 2, 3, or 4)
 * @param buffer_size Total buffer size in bytes (for bounds checking)
 * @param value Sample value to write
 * @param buf24 Pointer to int24 buffer for 24-bit operations
 * @return 1 on success, 0 if out of bounds
 */
int write_sample(unsigned char *buffer, int index, int bytes_per_sample,
                 int buffer_size, int32_t value, struct int24 *buf24);

/**
 * Calculate threshold value for trim operations.
 *
 * @param bytes_per_sample Number of bytes per sample (1, 2, 3, or 4)
 * @param threshold_percent Threshold as percentage of maximum value
 * @return Threshold value as signed 32-bit integer
 */
int32_t calculate_threshold(int bytes_per_sample, double threshold_percent);

/**
 * Convert duration in milliseconds to number of samples.
 *
 * Uses 64-bit arithmetic to detect potential integer overflow before
 * it occurs, preventing undefined behavior with long durations and high sample rates.
 *
 * @param duration_ms Duration in milliseconds
 * @param sample_rate Sample rate in Hz
 * @return Number of samples, or -1 if overflow would occur
 */
int duration_to_samples(int duration_ms, int sample_rate);

/**
 * Find first sample above threshold using optimized sampling strategy.
 *
 * For short files, uses linear search. For long files (>10 seconds), uses a
 * sampling strategy to find approximate location, then narrows the search.
 *
 * @param buffer Audio data buffer
 * @param buffer_size Buffer size in bytes
 * @param total_samples Total number of samples
 * @param bytes_per_sample Bytes per sample
 * @param threshold Threshold value
 * @param sample_rate Sample rate
 * @param num_channels Number of channels
 * @param buf24 Pointer to int24 buffer
 * @return Index of first sample above threshold
 */
int find_first_above_threshold(const unsigned char *buffer, int buffer_size,
                               int total_samples, int bytes_per_sample,
                               int32_t threshold, int sample_rate,
                               int num_channels, struct int24 *buf24);

/**
 * Find last sample above threshold using optimized sampling strategy.
 *
 * @param buffer Audio data buffer
 * @param buffer_size Buffer size in bytes
 * @param total_samples Total number of samples
 * @param bytes_per_sample Bytes per sample
 * @param threshold Threshold value
 * @param sample_rate Sample rate
 * @param num_channels Number of channels
 * @param buf24 Pointer to int24 buffer
 * @return Index of last sample above threshold
 */
int find_last_above_threshold(const unsigned char *buffer, int buffer_size,
                              int total_samples, int bytes_per_sample,
                              int32_t threshold, int sample_rate,
                              int num_channels, struct int24 *buf24);

/**
 * Find trim boundaries in audio buffer.
 *
 * @param header Pointer to WAV header
 * @param buffer Audio data buffer
 * @param buffer_size Buffer size in bytes
 * @param trim_start_percent Start trim threshold percentage
 * @param trim_end_percent End trim threshold percentage
 * @param first_sample_index Pointer to store first sample index
 * @param last_sample_index Pointer to store last sample index
 * @param buf24 Pointer to int24 buffer
 * @return 0 on success
 */
int find_trim_boundaries(const wav_header_t *header,
                         const unsigned char *buffer, int buffer_size,
                         double trim_start_percent, double trim_end_percent,
                         int *first_sample_index, int *last_sample_index,
                         struct int24 *buf24);

/**
 * Create fade curve lookup table.
 *
 * Pre-calculating fade values eliminates expensive floating-point math
 * from the inner loop, providing ~10-20x speedup for fade operations.
 *
 * @param fade_samples Number of samples in fade
 * @param curve_type Fade curve type
 * @param direction FADE_DIRECTION_IN or FADE_DIRECTION_OUT
 * @return Allocated lookup table, or NULL on error/edge case. Caller must free.
 */
double *create_fade_lut(int fade_samples, int curve_type, int direction);

/**
 * Apply fade effect to audio buffer.
 *
 * @param header Pointer to WAV header
 * @param buffer Audio data buffer
 * @param buffer_size Buffer size in bytes
 * @param fade_duration Duration in milliseconds
 * @param direction FADE_DIRECTION_IN or FADE_DIRECTION_OUT
 * @param curve_type Fade curve type (0=linear, 1=Bezier, 2=logarithmic)
 * @param buf24 Pointer to int24 buffer
 */
void apply_fade(const wav_header_t *header, unsigned char *buffer,
                int buffer_size, int fade_duration, int direction,
                int curve_type, struct int24 *buf24);

/**
 * Create a padded audio buffer.
 *
 * @param input Input audio buffer
 * @param input_size Input buffer size in bytes
 * @param header Pointer to WAV header
 * @param pad_start_ms Padding at start in milliseconds
 * @param pad_end_ms Padding at end in milliseconds
 * @param output_size Pointer to store output buffer size
 * @return Allocated padded buffer, or NULL on error. Caller must free.
 */
unsigned char *create_padded_buffer(const unsigned char *input, int input_size,
                                    const wav_header_t *header, int pad_start_ms,
                                    int pad_end_ms, int *output_size);

/**
 * Process audio data (trim, fade, pad).
 *
 * @param ctx Pointer to application context
 * @param header Pointer to WAV header
 * @param input_data Input audio data buffer
 * @param input_size Input data size in bytes
 * @param output_data Pointer to store allocated output buffer
 * @param output_size Pointer to store output size
 * @return 0 on success, 1 on error
 */
int process_audio(audio_fader_context_t *ctx, const wav_header_t *header,
                  const unsigned char *input_data, int input_size,
                  unsigned char **output_data, int *output_size);

/**
 * Validate options against audio file properties.
 *
 * Checks that fade durations don't exceed audio length, warns about
 * extreme settings, and validates that operations make sense together.
 *
 * @param options Pointer to options
 * @param header Pointer to WAV header
 * @return 0 on success, 1 on error
 */
int validate_options_with_audio(const options_t *options, const wav_header_t *header);

#endif /* AUDIOFADER_AUDIO_PROCESSING_H */
