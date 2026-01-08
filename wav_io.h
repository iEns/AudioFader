/*
 * AudioFader - WAV File I/O Module
 *
 * Functions for reading, writing, and validating WAV files.
 *
 * Copyright (c) 2020-2026 Jens Sandalgaard
 * Licensed under MIT License
 */

#ifndef AUDIOFADER_WAV_IO_H
#define AUDIOFADER_WAV_IO_H

#include "common.h"

/**
 * Validate WAV header fields.
 *
 * Checks for:
 * - Valid RIFF/WAVE/fmt identifiers
 * - Supported bit depths (8, 16, 24, 32)
 * - Supported channel counts (1 or 2)
 * - Valid sample rate range (8000-192000 Hz)
 *
 * @param header Pointer to WAV header structure
 * @return 0 on success, 1 on error
 */
int validate_wav_header(const wav_header_t *header);

/**
 * Skip to the data chunk in a WAV file (handles non-data chunks).
 *
 * WAV files may contain non-standard chunks before the data chunk.
 * This function scans through them to find the actual audio data.
 *
 * @param file File pointer
 * @param header Pointer to WAV header
 * @param file_size Total file size
 * @param data_chunk_size Pointer to store data chunk size
 * @return 0 on success, 1 on error
 */
int skip_to_data_chunk(FILE *file, const wav_header_t *header,
                       long file_size, int32_t *data_chunk_size);

/**
 * Get the full path of a file (platform-specific).
 *
 * On Windows, uses GetFullPathName().
 * On Unix/Linux, uses realpath().
 *
 * @param path Relative or absolute path
 * @return Allocated string with full path, or NULL on error. Caller must free.
 */
char *get_full_path(const char *path);

/**
 * Load and validate WAV file.
 *
 * Reads the WAV header, validates it, and loads all audio data into memory.
 *
 * @param ctx Pointer to application context
 * @param header Pointer to store WAV header
 * @param output_data Pointer to store allocated audio data buffer
 * @param output_size Pointer to store audio data size
 * @return 0 on success, 1 on error
 */
int load_wav_file(audio_fader_context_t *ctx, wav_header_t *header,
                  unsigned char **output_data, int *output_size);

/**
 * Write processed audio to output file.
 *
 * Updates the WAV header with new sizes and writes the complete file.
 *
 * @param ctx Pointer to application context
 * @param header Pointer to WAV header (will be modified)
 * @param data Audio data buffer
 * @param data_size Data size in bytes
 * @return 0 on success, 1 on error
 */
int write_output(audio_fader_context_t *ctx, wav_header_t *header,
                 const unsigned char *data, int data_size);

#endif /* AUDIOFADER_WAV_IO_H */
