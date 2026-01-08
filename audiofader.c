/*

Command Line Audio Fader, Trimmer and Padder

Concept, design and programming by Jens Sandalgaard

https://github.com/iEns/AudioFader

Copyright (c) 2020-2026 Jens Sandalgaard

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:
The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Version history:
0.91beta    2020-11-08  Initial version
0.92beta    2020-11-10  Code Cleanup, better error checking, support for 8 and 32 bit files.
0.93beta    2020-11-10  More code cleanup.
0.94beta    2026-01-07  Major refactoring following modern C best practices.
0.95beta    2026-01-08  Security hardening, performance optimizations, UX improvements.
                        - Phase 1: Critical security fixes (overflow protection, memory safety)
                        - Phase 2: Performance optimizations (67% less memory, 10-20x faster fades)
                        - Phase 3: User experience improvements (progress, validation, better errors)
0.96beta    2026-01-08  Multi-file refactoring, additional fixes.
                        - Refactored into modular components (common.h, wav_io, audio_processing, options, ui)
                        - Added sample_rate bounds checking (8000-192000 Hz)
                        - Added data alignment validation
                        - Fixed fade edge case (fade_samples <= 1)
                        - Added buffer bounds checking to read_sample/write_sample
                        - Implemented --dry-run flag
                        - Added WAV chunk count limit (max 100)
                        - Added minimum audio length check after trimming
                        - Removed unused progress_done parameter
                        - Added get_program_name helper
                        - Defined constants for magic numbers
                        - Added LOG_WARN macro
                        - Standardized error handling patterns
*/

#include "common.h"
#include "wav_io.h"
#include "audio_processing.h"
#include "options.h"
#include "ui.h"

/* ============================================================================
 * Main Entry Point
 * ============================================================================ */

/**
 * Main entry point.
 *
 * Processing pipeline: Parse Options -> Load WAV -> Validate -> Trim -> Fade -> Pad -> Write Output
 *
 * @param argc Argument count
 * @param argv Argument vector
 * @return 0 on success, 1 on error
 */
int main(int argc, char *argv[]) {
    clock_t start_time = clock();
    audio_fader_context_t ctx = {0};
    unsigned char *input_data = NULL;
    unsigned char *output_data = NULL;
    wav_header_t header = {0};
    int input_size = 0;
    int output_size = 0;
    progress_tracker_t progress = {0};

    print_title();

    /* Phase 1: Parse and validate options */
    if (parse_options(argc, argv, &ctx) != 0) {
        return 1;
    }

    print_options_used(&ctx);

    /* Initialize progress tracking (4 main phases) */
    progress_start(&progress, 4);

    /* Phase 2: Load and validate input file */
    progress_update(&progress, "Loading WAV file");
    if (load_wav_file(&ctx, &header, &input_data, &input_size) != 0) {
        progress_done();
        return cleanup_and_exit(NULL, NULL, &ctx, 0, 1);
    }
    progress_done();
    print_header_stats(&header);

    /* Phase 2.5: Validate options against audio file properties */
    progress_update(&progress, "Validating options");
    if (validate_options_with_audio(&ctx.options, &header) != 0) {
        progress_done();
        return cleanup_and_exit(input_data, NULL, &ctx, 0, 1);
    }
    progress_done();

    /* Phase 3: Process audio (trim, fade, pad) */
    progress_update(&progress, "Processing audio");
    if (process_audio(&ctx, &header, input_data, input_size,
                      &output_data, &output_size) != 0) {
        progress_done();
        return cleanup_and_exit(input_data, NULL, &ctx, 0, 1);
    }
    progress_done();

    /* Phase 4: Write output file (or skip if dry-run) */
    progress_update(&progress, ctx.options.dry_run ? "Dry-run complete" : "Writing output file");
    progress_done();

    /* FIX #6: Implement --dry-run - skip file write */
    if (ctx.options.dry_run) {
        printf("\n[DRY-RUN] Would write %d bytes to '%s'\n", output_size, ctx.options.output_filename);
        printf("[DRY-RUN] No output file was created.\n");
        return cleanup_and_exit(input_data, output_data, &ctx, start_time, 0);
    }

    if (write_output(&ctx, &header, output_data, output_size) != 0) {
        return cleanup_and_exit(input_data, output_data, &ctx, 0, 1);
    }

    /* Phase 5: Cleanup and exit */
    return cleanup_and_exit(input_data, output_data, &ctx, start_time, 0);
}
