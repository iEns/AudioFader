/*
 * AudioFader - Options Parsing Module
 *
 * Functions for parsing and validating command-line options.
 *
 * Copyright (c) 2020-2026 Jens Sandalgaard
 * Licensed under MIT License
 */

#include "options.h"
#include "wav_io.h"

/* ============================================================================
 * Argument Validation
 * ============================================================================ */

int validate_integer_arg(const char *arg, int *out_value, int min,
                         int max, const char *name) {
    char *endptr;
    errno = 0;
    long value = strtol(arg, &endptr, 10);

    if (errno != 0 || *endptr != '\0' || endptr == arg) {
        LOG_ERROR("Error: Invalid %s argument '%s' (must be an integer)\n", name, arg);
        return 1;
    }

    if (value < min || value > max) {
        LOG_ERROR("Error: %s value %ld out of range [%d, %d]\n", name, value, min, max);
        return 1;
    }

    *out_value = (int)value;
    return 0;
}

int validate_double_arg(const char *arg, double *out_value, double min,
                        double max, const char *name) {
    char *endptr;
    errno = 0;
    double value = strtod(arg, &endptr);

    if (errno != 0 || *endptr != '\0' || endptr == arg) {
        LOG_ERROR("Error: Invalid %s argument '%s' (must be a number)\n", name, arg);
        return 1;
    }

    if (value < min || value > max) {
        LOG_ERROR("Error: %s value %f out of range [%f, %f]\n", name, value, min, max);
        return 1;
    }

    *out_value = value;
    return 0;
}

/* ============================================================================
 * Output Functions
 * ============================================================================ */

void print_title(void) {
    printf("\nAudioFader v%s %s\n\n", AUDIOFADER_VERSION, AUDIOFADER_COPYRIGHT);
}

void print_help(const char *program_name) {
    printf("Trim, fade and pad WAV audio files.\n\n");
    printf("USAGE:\n");
    printf("  %s <input.wav> <output.wav> [OPTIONS]\n\n", program_name);

    printf("AUDIO PROCESSING OPTIONS:\n");
    printf("  Trimming (removes silence from ends):\n");
    printf("    --trimstart N     Trim silence from start (0-100%% threshold, default: 0)\n");
    printf("    --trimend N       Trim silence from end (0-100%% threshold, default: 0)\n");
    printf("                      Common values: 0.01 to 1.0\n\n");

    printf("  Fading:\n");
    printf("    --fadein MS       Fade-in duration in milliseconds\n");
    printf("    --fadeout MS      Fade-out duration in milliseconds\n");
    printf("    --fadecurve N     Fade curve: 0=Linear (default), 1=Bezier, 2=Logarithmic\n\n");

    printf("  Padding (adds silence):\n");
    printf("    --padstart MS     Add silence at start (milliseconds)\n");
    printf("    --padend MS       Add silence at end (milliseconds)\n\n");

    printf("  Other:\n");
    printf("    --dry-run         Preview mode - show what would be done without writing output\n");
    printf("    --force           Overwrite existing output file without prompting\n");
    printf("    --quiet           Quiet mode (errors only)\n");
    printf("    --verbose         Verbose output (debug information)\n");
    printf("    --help            Show this help message\n\n");

    printf("EXAMPLES:\n");
    printf("  Basic fade:\n");
    printf("    %s in.wav out.wav --fadein 500 --fadeout 500\n\n", program_name);

    printf("  Trim silence and fade with Bezier curve:\n");
    printf("    %s in.wav out.wav --trimstart 0.01 --fadein 1000 --fadeout 1000 --fadecurve 1\n\n",
           program_name);

    printf("  Add 1 second padding:\n");
    printf("    %s in.wav out.wav --padstart 1000 --padend 1000\n\n", program_name);

    printf("  Preview changes without writing:\n");
    printf("    %s in.wav out.wav --trimstart 0.01 --dry-run\n\n", program_name);

    printf("  Complete workflow:\n");
    printf("    %s in.wav out.wav --trimstart 0.01 --trimend 0.01 \\\n", program_name);
    printf("                      --fadein 2000 --fadeout 2000 --fadecurve 1 \\\n");
    printf("                      --padstart 500 --padend 500\n\n");

    printf("PROCESSING ORDER:\n");
    printf("  Operations are applied in this sequence: Trim -> Fade -> Pad\n\n");

    printf("SUPPORTED FORMATS:\n");
    printf("  - WAV PCM (uncompressed)\n");
    printf("  - Bit depths: 8, 16, 24, 32 bits\n");
    printf("  - Channels: mono or stereo\n\n");
}

void print_options_used(const audio_fader_context_t *ctx) {
    printf("Options in use:\n");
    printf("\tInput file: \t\t%s\n", ctx->options.input_filename);
    printf("\tOutput file: \t\t%s\n", ctx->options.output_filename);

    if (ctx->options.dry_run) {
        printf("\tMode: \t\t\tDRY-RUN (no output file will be written)\n");
    }

    if (ctx->options.fade_in != 0) {
        printf("\tFade-in duration: \t%d ms\n", ctx->options.fade_in);
    }
    if (ctx->options.fade_out != 0) {
        printf("\tFade-out duration: \t%d ms\n", ctx->options.fade_out);
    }
    if (ctx->options.fade_in != 0 || ctx->options.fade_out != 0) {
        printf("\tFade-curve: \t\t%s\n", FADE_CURVE_NAMES[ctx->options.fade_curve]);
    }
    if (ctx->options.trim_start != 0) {
        printf("\tTrim start: \t\t%f%%\n", ctx->options.trim_start);
    }
    if (ctx->options.trim_end != 0) {
        printf("\tTrim end: \t\t%f%%\n", ctx->options.trim_end);
    }
    if (ctx->options.pad_start != 0) {
        printf("\tPad start: \t\t%d ms\n", ctx->options.pad_start);
    }
    if (ctx->options.pad_end != 0) {
        printf("\tPad end: \t\t%d ms\n", ctx->options.pad_end);
    }
    printf("\n");
}

void print_header_stats(const wav_header_t *header) {
    printf("Input File Information:\n");
    printf("\tChannels: \t\t%d\n", header->num_channels);
    printf("\tBits Per Sample: \t%d\n", header->bits_per_sample);
    printf("\tSample Rate: \t\t%d\n", header->sample_rate);
    printf("\tAudioFormat: \t\t%d\n", header->audio_format);
    printf("\tBlockAlign: \t\t%d\n", header->block_align);
    printf("\n");
}

/* ============================================================================
 * Options Parsing
 * ============================================================================ */

/**
 * Best-effort same-file check.
 *
 * 1. If both paths stat OK, compare (st_dev, st_ino): catches
 *    "in.wav" vs "./in.wav" vs symlinks even when realpath() failed
 *    for the not-yet-existing output.
 * 2. Otherwise compare canonical strings when available. POSIX uses
 *    case-sensitive strcmp (the old strcasecmp wrongly equated
 *    "In.Wav" with "in.wav" on Linux); Windows keeps _stricmp.
 */
static int paths_refer_to_same_file(const char *input_raw,
                                    const char *output_raw,
                                    const char *input_resolved,
                                    const char *output_resolved) {
#ifndef _WIN32
    struct stat sta, stb;
    if (input_raw != NULL && output_raw != NULL &&
        stat(input_raw, &sta) == 0 && stat(output_raw, &stb) == 0) {
        return (sta.st_dev == stb.st_dev && sta.st_ino == stb.st_ino);
    }
#endif
    if (input_resolved != NULL && output_resolved != NULL) {
#ifdef _WIN32
        return strcasecmp(input_resolved, output_resolved) == 0;
#else
        return strcmp(input_resolved, output_resolved) == 0;
#endif
    }
    if (input_raw != NULL && output_raw != NULL) {
#ifdef _WIN32
        return strcasecmp(input_raw, output_raw) == 0;
#else
        return strcmp(input_raw, output_raw) == 0;
#endif
    }
    return 0;
}

static void free_resolved_paths(audio_fader_context_t *ctx) {
    free(ctx->options.input_path_resolved);
    free(ctx->options.output_path_resolved);
    ctx->options.input_path_resolved = NULL;
    ctx->options.output_path_resolved = NULL;
}

int parse_options(int argc, char *argv[], audio_fader_context_t *ctx) {
    /* Check for --help first (before minimum arg check) */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--help") == 0) {
            const char *program_name = get_program_name(argv[0]);
            print_help(program_name);
            return 2;  /* Help shown - not an error */
        }
    }

    /* Check minimum argument count */
    if (argc < MIN_ARGC) {
        const char *program_name = get_program_name(argv[0]);
        print_help(program_name);
        return 1;
    }

    /* Initialize options with defaults */
    memset(&ctx->options, 0, sizeof(options_t));
    ctx->options.verbosity = VERBOSITY_NORMAL;  /* Default verbosity */

    /* Store original argv pointers (not owned, don't free) */
    ctx->options.input_filename = argv[1];
    ctx->options.output_filename = argv[2];

    /* Expand filenames to full paths (these are heap-allocated, must free) */
    ctx->options.input_path_resolved = get_full_path(ctx->options.input_filename);
    ctx->options.output_path_resolved = get_full_path(ctx->options.output_filename);

    /* Check that input and output are different (stat-based, so
     * "in.wav" vs "./in.wav" is caught even when realpath() returns
     * NULL for the not-yet-existing output file). */
    if (paths_refer_to_same_file(ctx->options.input_filename,
                                 ctx->options.output_filename,
                                 ctx->options.input_path_resolved,
                                 ctx->options.output_path_resolved)) {
        LOG_ERROR("Error: Input and output file names cannot be the same\n");
        LOG_ERROR("       Input:  %s\n", ctx->options.input_filename);
        LOG_ERROR("       Output: %s\n", ctx->options.output_filename);
        free_resolved_paths(ctx);
        return 1;
    }

    /* Parse command-line options */
    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "--fadein") == 0) {
            if (i + 1 >= argc || argv[i + 1][0] == '-') {
                LOG_ERROR("Error: Option --fadein requires a numeric argument (milliseconds)\n");
                goto error_cleanup;
            }
            if (validate_integer_arg(argv[i + 1], &ctx->options.fade_in,
                                     0, MAX_DURATION_MS, "fade-in") != 0) {
                goto error_cleanup;
            }
            i++;
        } else if (strcmp(argv[i], "--fadeout") == 0) {
            if (i + 1 >= argc || argv[i + 1][0] == '-') {
                LOG_ERROR("Error: Option --fadeout requires a numeric argument (milliseconds)\n");
                goto error_cleanup;
            }
            if (validate_integer_arg(argv[i + 1], &ctx->options.fade_out,
                                     0, MAX_DURATION_MS, "fade-out") != 0) {
                goto error_cleanup;
            }
            i++;
        } else if (strcmp(argv[i], "--fadecurve") == 0) {
            if (i + 1 >= argc || argv[i + 1][0] == '-') {
                LOG_ERROR("Error: Option --fadecurve requires an argument (0=Linear, 1=Bezier, 2=Logarithmic)\n");
                goto error_cleanup;
            }
            if (validate_integer_arg(argv[i + 1], &ctx->options.fade_curve,
                                     0, NUM_FADE_CURVES - 1, "fade-curve") != 0) {
                goto error_cleanup;
            }
            i++;
        } else if (strcmp(argv[i], "--trimstart") == 0) {
            if (i + 1 >= argc || argv[i + 1][0] == '-') {
                LOG_ERROR("Error: Option --trimstart requires a numeric argument (0-100 percent)\n");
                goto error_cleanup;
            }
            if (validate_double_arg(argv[i + 1], &ctx->options.trim_start,
                                    0, MAX_TRIM_PERCENT, "trim-start") != 0) {
                goto error_cleanup;
            }
            i++;
        } else if (strcmp(argv[i], "--trimend") == 0) {
            if (i + 1 >= argc || argv[i + 1][0] == '-') {
                LOG_ERROR("Error: Option --trimend requires a numeric argument (0-100 percent)\n");
                goto error_cleanup;
            }
            if (validate_double_arg(argv[i + 1], &ctx->options.trim_end,
                                    0, MAX_TRIM_PERCENT, "trim-end") != 0) {
                goto error_cleanup;
            }
            i++;
        } else if (strcmp(argv[i], "--padstart") == 0) {
            if (i + 1 >= argc || argv[i + 1][0] == '-') {
                LOG_ERROR("Error: Option --padstart requires a numeric argument (milliseconds)\n");
                goto error_cleanup;
            }
            if (validate_integer_arg(argv[i + 1], &ctx->options.pad_start,
                                     0, MAX_DURATION_MS, "pad-start") != 0) {
                goto error_cleanup;
            }
            i++;
        } else if (strcmp(argv[i], "--padend") == 0) {
            if (i + 1 >= argc || argv[i + 1][0] == '-') {
                LOG_ERROR("Error: Option --padend requires a numeric argument (milliseconds)\n");
                goto error_cleanup;
            }
            if (validate_integer_arg(argv[i + 1], &ctx->options.pad_end,
                                     0, MAX_DURATION_MS, "pad-end") != 0) {
                goto error_cleanup;
            }
            i++;
        } else if (strcmp(argv[i], "--force") == 0) {
            ctx->options.force_overwrite = 1;
        } else if (strcmp(argv[i], "--dry-run") == 0) {
            /* FIX #6: Implement --dry-run flag */
            ctx->options.dry_run = 1;
        } else if (strcmp(argv[i], "--quiet") == 0) {
            if (ctx->options.verbosity == VERBOSITY_VERBOSE) {
                LOG_ERROR("Error: --quiet and --verbose cannot be combined\n");
                goto error_cleanup;
            }
            ctx->options.verbosity = VERBOSITY_QUIET;
        } else if (strcmp(argv[i], "--verbose") == 0) {
            if (ctx->options.verbosity == VERBOSITY_QUIET) {
                LOG_ERROR("Error: --quiet and --verbose cannot be combined\n");
                goto error_cleanup;
            }
            ctx->options.verbosity = VERBOSITY_VERBOSE;
        } else {
            /* --help is handled at start of parse_options */
            LOG_ERROR("Error: Unknown option: %s\n", argv[i]);
            LOG_ERROR("       Use --help to see available options\n");
            goto error_cleanup;
        }
    }

    /* Set global verbosity level for logging macros */
    g_verbosity = ctx->options.verbosity;

    return 0;

error_cleanup:
    free_resolved_paths(ctx);
    return 1;
}
