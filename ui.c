/*
 * AudioFader - User Interface Module
 *
 * Functions for progress reporting, logging, and cleanup.
 *
 * Copyright (c) 2020-2026 Jens Sandalgaard
 * Licensed under MIT License
 */

#include "ui.h"

/* ============================================================================
 * Global Variables (defined here, declared extern in common.h)
 * ============================================================================ */

verbosity_level_t g_verbosity = VERBOSITY_NORMAL;

const char *FADE_CURVE_NAMES[] = {"Linear", "Bezier", "Logarithmic"};

/* ============================================================================
 * Progress Tracking Functions
 * ============================================================================ */

void progress_start(progress_tracker_t *tracker, int total) {
    tracker->start_time = clock();
    tracker->total_steps = total;
    tracker->current_step = 0;
}

void progress_update(progress_tracker_t *tracker, const char *operation) {
    tracker->current_step++;
    double elapsed = (double)(clock() - tracker->start_time) / CLOCKS_PER_SEC;
    int percent = (tracker->current_step * 100) / tracker->total_steps;

    if (g_verbosity >= VERBOSITY_NORMAL) {
        fprintf(stderr, "\r[%3d%%] %s (%.1fs elapsed)", percent, operation, elapsed);
        fflush(stderr);
    }
}

/* FIX #9: Removed unused tracker parameter */
void progress_done(void) {
    if (g_verbosity >= VERBOSITY_NORMAL) {
        fprintf(stderr, "\n");
    }
}

/* ============================================================================
 * Cleanup Functions
 * ============================================================================ */

int cleanup_and_exit(unsigned char *input_data, unsigned char *output_data,
                     audio_fader_context_t *ctx, clock_t start_time,
                     int exit_code) {
    if (input_data != NULL) {
        free(input_data);
    }
    if (output_data != NULL) {
        free(output_data);
    }

    /* Free resolved paths (heap-allocated in parse_options) */
    if (ctx->options.input_path_resolved != NULL) {
        free(ctx->options.input_path_resolved);
        ctx->options.input_path_resolved = NULL;
    }
    if (ctx->options.output_path_resolved != NULL) {
        free(ctx->options.output_path_resolved);
        ctx->options.output_path_resolved = NULL;
    }

    /* Do NOT free input_filename/output_filename - they point to argv[] */

    if (start_time != 0 && exit_code == 0) {
        double elapsed = (double)(clock() - start_time) / CLOCKS_PER_SEC;
        printf("\nFinished processing file in %f seconds.\n\n", elapsed);
    }

    return exit_code;
}
