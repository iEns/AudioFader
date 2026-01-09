/*
 * AudioFader - Options Parsing Module
 *
 * Functions for parsing and validating command-line options.
 *
 * Copyright (c) 2020-2026 Jens Sandalgaard
 * Licensed under MIT License
 */

#ifndef AUDIOFADER_OPTIONS_H
#define AUDIOFADER_OPTIONS_H

#include "common.h"

/**
 * Validate and parse an integer argument.
 *
 * @param arg String argument to parse
 * @param out_value Pointer to store parsed value
 * @param min Minimum valid value
 * @param max Maximum valid value
 * @param name Argument name for error messages
 * @return 0 on success, 1 on error
 */
int validate_integer_arg(const char *arg, int *out_value, int min,
                         int max, const char *name);

/**
 * Validate and parse a double argument.
 *
 * @param arg String argument to parse
 * @param out_value Pointer to store parsed value
 * @param min Minimum valid value
 * @param max Maximum valid value
 * @param name Argument name for error messages
 * @return 0 on success, 1 on error
 */
int validate_double_arg(const char *arg, double *out_value, double min,
                        double max, const char *name);

/**
 * Parse and validate command-line options.
 *
 * @param argc Argument count
 * @param argv Argument vector
 * @param ctx Pointer to application context
 * @return 0 on success, 1 on error, 2 if help shown (not an error)
 */
int parse_options(int argc, char *argv[], audio_fader_context_t *ctx);

/**
 * Print usage help message.
 *
 * @param program_name Name of the executable
 */
void print_help(const char *program_name);

/**
 * Print the application title and version.
 */
void print_title(void);

/**
 * Print the options that are currently in use.
 *
 * @param ctx Pointer to application context
 */
void print_options_used(const audio_fader_context_t *ctx);

/**
 * Print WAV file header statistics.
 *
 * @param header Pointer to WAV header structure
 */
void print_header_stats(const wav_header_t *header);

#endif /* AUDIOFADER_OPTIONS_H */
