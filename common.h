/*
 * AudioFader - Common Definitions
 *
 * Shared types, constants, and macros used across all modules.
 *
 * Copyright (c) 2020-2026 Jens Sandalgaard
 * Licensed under MIT License
 */

#ifndef AUDIOFADER_COMMON_H
#define AUDIOFADER_COMMON_H

#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _WIN32
    #include <windows.h>
    #ifdef _MSC_VER
        #define strcasecmp _stricmp
    #endif
#else
    #include <stdlib.h>
    #include <strings.h>
    #include <sys/stat.h>
#endif

/* ============================================================================
 * Version Information
 * ============================================================================ */

#define AUDIOFADER_VERSION "0.96beta"
#define AUDIOFADER_COPYRIGHT "Copyright 2021-2026 iEns Labs"

/* ============================================================================
 * Constants
 * ============================================================================ */

#define MAX_DURATION_MS 1000000
#define MAX_TRIM_PERCENT 100.0
#define MIN_ARGC 3
#define MAX_AUDIO_SIZE (1024 * 1024 * 1024)  /* 1GB limit for safety */

/* Sample max values for each bit depth */
#define INT8_MAX_SAMPLE 0x7F
#define INT16_MAX_SAMPLE 0x7FFF
#define INT24_MAX_SAMPLE 0x7FFFFF
#define INT32_MAX_SAMPLE 0x7FFFFFFF

/* WAV file format identifiers */
#define WAV_RIFF_ID "RIFF"
#define WAV_WAVE_ID "WAVE"
#define WAV_FMT_ID "fmt "
#define WAV_DATA_ID "data"
#define WAV_HEADER_EXTRA_SIZE 36

/* Fade direction constants */
#define FADE_DIRECTION_IN 0
#define FADE_DIRECTION_OUT 1

/* Fade curve types */
#define FADE_CURVE_LINEAR 0
#define FADE_CURVE_BEZIER 1
#define FADE_CURVE_LOGARITHMIC 2
#define NUM_FADE_CURVES 3

/* Optimization constants (extracted from magic numbers) */
#define TRIM_OPTIMIZATION_THRESHOLD_SECONDS 10
#define TRIM_SAMPLING_DIVISOR 10  /* Check every 1/10th second */

/* Maximum number of WAV chunks to scan (prevents malformed file DoS) */
#define MAX_WAV_CHUNKS 100

/* Minimum sample rate and maximum supported */
#define MIN_SAMPLE_RATE 8000
#define MAX_SAMPLE_RATE 192000

/* ============================================================================
 * Type Definitions
 * ============================================================================ */

/* WAV file header structure */
typedef struct {
    char chunk_id[4];
    int32_t chunk_size;
    char format[4];
    char subchunk1_id[4];
    int32_t subchunk1_size;
    int16_t audio_format;
    int16_t num_channels;
    int32_t sample_rate;
    int32_t byte_rate;
    int16_t block_align;
    int16_t bits_per_sample;
    char subchunk2_id[4];
    int32_t subchunk2_size;
} wav_header_t;

/* Subchunk header structure */
typedef struct {
    char subchunk_id[4];
    int32_t subchunk_size;
} sub_chunk_header_t;

/* Verbosity levels */
typedef enum {
    VERBOSITY_QUIET = 0,   /* Errors only */
    VERBOSITY_NORMAL = 1,  /* Default behavior */
    VERBOSITY_VERBOSE = 2  /* Debug information */
} verbosity_level_t;

/* Command-line options */
typedef struct {
    const char *input_filename;   /* Original from argv (not owned) */
    const char *output_filename;  /* Original from argv (not owned) */
    char *input_path_resolved;    /* Heap-allocated full path (owned, must free) */
    char *output_path_resolved;   /* Heap-allocated full path (owned, must free) */
    int fade_curve;
    int fade_in;
    int fade_out;
    double trim_start;
    double trim_end;
    int pad_start;
    int pad_end;
    int force_overwrite;          /* Allow overwriting existing files */
    int dry_run;                  /* Preview mode - don't write output */
    verbosity_level_t verbosity;  /* Output verbosity level */
} options_t;

/* 24-bit integer type (GCC-specific bitfield extension) */
struct int24 {
    signed int24 : 24;
} __attribute__((packed));

/* Progress tracker for large file processing */
typedef struct {
    clock_t start_time;
    int total_steps;
    int current_step;
} progress_tracker_t;

/* Application context */
typedef struct {
    options_t options;
    struct int24 sample24_buffer;
} audio_fader_context_t;

/* ============================================================================
 * Global Variables (declared extern, defined in ui.c)
 * ============================================================================ */

extern verbosity_level_t g_verbosity;

/* Fade curve names */
extern const char *FADE_CURVE_NAMES[];

/* ============================================================================
 * Logging Macros
 * ============================================================================ */

#define LOG_ERROR(...) fprintf(stderr, __VA_ARGS__)

#define LOG_WARN(...) do { \
    if (g_verbosity >= VERBOSITY_NORMAL) \
        fprintf(stderr, "Warning: " __VA_ARGS__); \
} while(0)

#define LOG_INFO(...) do { \
    if (g_verbosity >= VERBOSITY_NORMAL) \
        printf(__VA_ARGS__); \
} while(0)

#define LOG_VERBOSE(...) do { \
    if (g_verbosity >= VERBOSITY_VERBOSE) \
        printf(__VA_ARGS__); \
} while(0)

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

/**
 * Safely add two integers, detecting overflow.
 *
 * @param a First integer
 * @param b Second integer
 * @param result Pointer to store result
 * @return 1 on success, 0 if overflow would occur
 */
static inline int safe_add_int(int a, int b, int *result) {
    if (b > 0 && a > INT_MAX - b) {
        return 0;  /* Overflow */
    }
    if (b < 0 && a < INT_MIN - b) {
        return 0;  /* Underflow */
    }
    *result = a + b;
    return 1;
}

/**
 * Get the path separator character for this platform.
 *
 * @return Path separator ('\\' for Windows, '/' for Unix)
 */
static inline char get_path_separator(void) {
#ifdef _WIN32
    return '\\';
#else
    return '/';
#endif
}

/**
 * Get the program name from argv[0].
 *
 * @param argv0 The argv[0] value
 * @return Pointer to program name within argv0 string (not a copy)
 */
static inline const char *get_program_name(const char *argv0) {
    char separator = get_path_separator();
    const char *name = strrchr(argv0, separator);
    return (name != NULL) ? name + 1 : argv0;
}

#endif /* AUDIOFADER_COMMON_H */
