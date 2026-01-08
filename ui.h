/*
 * AudioFader - User Interface Module
 *
 * Functions for progress reporting, logging, and cleanup.
 *
 * Copyright (c) 2020-2026 Jens Sandalgaard
 * Licensed under MIT License
 */

#ifndef AUDIOFADER_UI_H
#define AUDIOFADER_UI_H

#include "common.h"

/**
 * Initialize progress tracker.
 *
 * @param tracker Pointer to progress tracker
 * @param total Total number of steps
 */
void progress_start(progress_tracker_t *tracker, int total);

/**
 * Update progress and display to user.
 *
 * @param tracker Pointer to progress tracker
 * @param operation Description of current operation
 */
void progress_update(progress_tracker_t *tracker, const char *operation);

/**
 * Complete progress tracking.
 *
 * FIX #9: Removed unused tracker parameter from original design.
 */
void progress_done(void);

/**
 * Cleanup resources and exit.
 *
 * @param input_data Input data buffer to free (can be NULL)
 * @param output_data Output data buffer to free (can be NULL)
 * @param ctx Application context
 * @param start_time Start time for timing (use 0 to skip timing)
 * @param exit_code Exit code to return
 * @return exit_code
 */
int cleanup_and_exit(unsigned char *input_data, unsigned char *output_data,
                     audio_fader_context_t *ctx, clock_t start_time,
                     int exit_code);

#endif /* AUDIOFADER_UI_H */
