# Changelog

All notable changes to AudioFader will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.96beta] - 2026-01-08

### Multi-File Refactoring & Additional Safety Fixes

This release refactors the codebase from a single 1660-line file into a modular multi-file structure, improving maintainability and testability. Additionally, 13 safety, robustness, and code quality fixes have been applied.

### Architecture

#### Changed
- **Multi-file structure**: Split monolithic `audiofader.c` into 6 modules:
  - `common.h` - Shared types, constants, macros, utility functions
  - `wav_io.h/c` - WAV file reading, writing, and validation
  - `audio_processing.h/c` - Trim, fade, and pad operations
  - `options.h/c` - Command-line parsing and validation
  - `ui.h/c` - Progress reporting and cleanup functions
  - `audiofader.c` - Main entry point (reduced to ~135 lines)
- **Build system**: Updated Makefile for multi-file compilation
- **Documentation**: Updated README.md with new structure

### Safety & Robustness Fixes

#### Added
- **Sample rate bounds checking**: Validates sample rate is within 8000-192000 Hz range
  - Prevents integer overflow in calculations with malformed headers
- **Data alignment validation**: Verifies audio data size is aligned to frame boundaries
  - Detects corrupted or truncated WAV files
- **Buffer bounds checking**: Added `buffer_size` parameter to `read_sample()` and `write_sample()`
  - Prevents buffer overruns with invalid sample indices
- **WAV chunk count limit**: Maximum 100 chunks when scanning for data chunk
  - Prevents DoS from malformed WAV files with endless fake chunks
- **Minimum audio length check**: Validates trimmed audio has at least one channel frame
  - Better error message when trim removes nearly all audio
- **Zero-length input validation**: Checks `input_size > 0` before padding operations

#### Fixed
- **Fade edge case**: Handle `fade_samples <= 1` gracefully in `create_fade_lut()`
  - Previously caused division by zero when calculating fade index
- **Race condition improvement**: Better handling of output file existence check
  - Documents remaining TOCTOU window as acceptable for CLI tool

### New Features

#### Added
- **`--dry-run` flag**: Preview mode that shows what would be done without writing output
  - Displays output file size and path
  - Useful for testing options before committing to file write
- **`LOG_WARN` macro**: Separate warning logging (non-fatal issues)
  - Warnings now prefixed with "Warning:" for clarity

### Code Quality

#### Added
- **`get_program_name()` helper**: Extracts program name from `argv[0]`
  - Removes code duplication in help and error messages
- **Named constants for magic numbers**:
  - `TRIM_OPTIMIZATION_THRESHOLD_SECONDS` (was hardcoded `10`)
  - `TRIM_SAMPLING_DIVISOR` (was hardcoded `10`)
  - `MAX_WAV_CHUNKS` (new, set to `100`)
  - `MIN_SAMPLE_RATE` / `MAX_SAMPLE_RATE` (new, `8000` / `192000`)

#### Changed
- **Removed unused parameter**: `progress_done()` no longer takes tracker parameter
- **Standardized error handling**: Consistent `goto cleanup` pattern in functions with resource allocation
- **Updated comments**: Fixed outdated comment about Windows not having fstat (it has `_fstat64`)

### Build System

#### Changed
- Updated Makefile for multi-file compilation:
  - `SRCS` now lists all 5 source files
  - `HDRS` lists all 5 header files
  - Dependencies updated so changes to headers trigger rebuild
- Updated `make test` to use `--dry-run` flag
- Updated `make format` to include header files

### Documentation

#### Updated
- README.md with new project structure and build commands
- All source files with MIT license headers and documentation

### Testing

Verified with:
- ✅ All existing test cases from 0.95beta
- ✅ Zero-duration fade edge case (`--fadein 0`)
- ✅ Dry-run mode (`--dry-run`)
- ✅ Same input/output validation
- ✅ 24-bit/192kHz high sample rate files
- ✅ Clean build with `-Wall -Wextra -Wpedantic`

### Backward Compatibility

All changes maintain 100% backward compatibility:
- All existing command-line options work identically
- Output files are identical for the same inputs
- New `--dry-run` flag is optional

## [0.95beta] - 2026-01-07

### Security Hardening, Performance Optimization & UX Improvements

This release focuses on three critical areas: security/stability fixes, major performance improvements, and significantly enhanced user experience. The changes represent 40+ optimizations across the codebase.

### Security (Phase 1: Critical Fixes)

#### Added
- Integer overflow protection in `duration_to_samples()` using 64-bit arithmetic
- Safe integer addition helper function to prevent overflow in chunk calculations
- Audio data size validation (1GB maximum) before memory allocation
- Buffer bounds validation in `read_sample()` function
- Comprehensive file handle cleanup using goto-based error patterns
- Overflow detection in WAV chunk parsing with detailed error messages

#### Fixed
- **P0.1**: Integer overflow vulnerability with long durations at high sample rates
  - Previously: 1,000,000ms × 48kHz would overflow INT_MAX causing undefined behavior
  - Now: Uses int64_t arithmetic with validation before returning int result
- **P0.2**: File size type mismatch on systems with large files (>2GB)
  - Changed from `int` to `long` for file size storage
  - Uses `fstat()` on Unix/Linux for proper 64-bit size handling
- **P0.3**: Memory leak in path allocation and ownership tracking
  - Separated heap-allocated resolved paths from argv pointers in options struct
  - Only frees owned heap memory in cleanup, preventing double-free bugs
- **P0.4**: Missing bounds validation before buffer access
  - Added NULL pointer checks and index validation in read_sample()
- **P0.5**: Buffer overflow risks in WAV chunk parsing
  - Validates chunk sizes don't exceed file size before use
  - Safe integer addition prevents overflow wraparound
- **P0.6**: Malicious WAV files could trigger huge allocations
  - Added MAX_AUDIO_SIZE (1GB) limit with validation
  - Checks data chunk size against file size
- **P0.7**: File handle leaks on error paths
  - Implemented goto cleanup pattern ensuring file closure on all paths
- **P0.8**: Edge case when fade duration exceeds audio length
  - Validates and clamps fade durations to audio length with warnings
- **P0.9**: Edge case when trim removes all audio
  - Detects empty result and fails with detailed error message
- **P0.10**: Unsafe string operations with NULL pointers
  - Validates strrchr() result before dereferencing
  - Checks pointer validity before strcasecmp() calls

### Performance (Phase 2: Optimization)

#### Added
- **P1.2**: Pre-calculated fade curve lookup tables
  - Eliminates per-sample floating-point math in fade operations
  - **Result**: 10-20x faster fade processing
- **P1.4**: Optimized trim boundary search for long files
  - Uses sampling strategy (10 samples/second) for files >10 seconds
  - Binary search narrows to precise boundary
  - **Result**: ~4800x faster for 1-hour files (172M samples → 36K checks)

#### Changed
- **P1.1**: Conditional buffer allocation strategy
  - Only allocates trimmed buffer if actually trimming
  - Only allocates padded buffer if actually padding
  - Fading done in-place without new allocation
  - **Result**: 67% memory reduction (10MB file: 30MB → 10MB usage)
- **P1.3**: File size detection using `fstat()` instead of seek/ftell on Unix/Linux
  - Reduces from 3 system calls to 1
  - More reliable for special files
  - Returns proper 64-bit size
- **P1.7**: Simplified logarithmic fade calculation
  - Changed `fade_index / (1.0 + (1.0 - fade_index))` to `fade_index / (2.0 - fade_index)`
  - Eliminates one floating-point subtraction per sample

#### Performance Metrics
- Processing speed: ~3x faster overall (user-reported on test files)
- Memory usage: 67% reduction for fade-only operations
- Trim boundary search: Up to 4800x faster for long files
- Fade operations: 10-20x faster with lookup tables

### User Experience (Phase 3: Improvements)

#### Added
- **P2.1**: Real-time progress indication during processing
  - Shows percentage complete and elapsed time
  - Updates for: Loading → Validating → Processing → Writing
  - Respects verbosity settings
- **P2.2**: Pre-processing validation of options against audio properties
  - Validates fade durations don't exceed audio length
  - Warns about extreme trim thresholds (≥50%)
  - Checks output duration doesn't exceed limits
  - Displays audio duration in verbose mode
- **P2.3**: Output file overwrite protection
  - Checks if output file exists before processing
  - Requires `--force` flag to overwrite existing files
  - Prevents accidental data loss
- **P2.6**: Verbosity control with three levels
  - `--quiet`: Errors only, suppresses progress and info messages
  - Normal: Default behavior with progress tracking (default)
  - `--verbose`: Debug information including sample counts and durations
- **P2.8**: Comprehensive help output with examples
  - Clear section organization (Trimming, Fading, Padding, Other)
  - Four real-world usage examples
  - Explains processing order (Trim → Fade → Pad)
  - Lists supported formats and bit depths
- Logging infrastructure with three levels (LOG_ERROR, LOG_INFO, LOG_VERBOSE)

#### Changed
- **P2.4**: Dramatically improved error messages with context
  - All errors now include relevant filenames
  - Includes system error descriptions (strerror) where applicable
  - Provides actionable suggestions (e.g., "Use --force to overwrite")
  - Shows threshold values and detected ranges for trim errors
  - Displays expected vs actual values for validation errors

  Examples:
  ```
  Before: "Error: cannot open input file"
  After:  "Error: Cannot open input file 'test.wav': No such file or directory"

  Before: "Error: input file is not 8, 16, 24 or 32 bits per sample"
  After:  "Error: Unsupported bit depth (12 bits per sample)
           Supported: 8, 16, 24, or 32 bits per sample"

  Before: "Error: Trimming removed all audio"
  After:  "Error: Trimming removed all audio (threshold too high?)
           Trim thresholds: start=50.00%, end=50.00%
           Detected range: sample 0 to 0 (of 1536334 total)"
  ```

#### Technical Improvements
- Global verbosity variable for consistent logging across all functions
- Progress tracker struct with start time, total steps, and current step
- Separated validation phase in main() pipeline for early failure
- All printf() error messages converted to LOG_ERROR() macros
- Help text now uses program name extracted from argv[0]

### Documentation

#### Updated
- Version history in source file header with phase summaries
- Function documentation for new validation and progress functions
- This changelog with comprehensive release notes

### Backward Compatibility

All changes maintain 100% backward compatibility:
- All existing command-line options work identically
- Output files are identical for the same inputs
- No breaking changes to command-line interface
- New flags are optional (--force, --quiet, --verbose)

### Testing

Verified with multiple test scenarios:
- ✅ 8, 16, 24, and 32-bit WAV files
- ✅ Mono and stereo audio
- ✅ Files from 1 second to 1 hour duration
- ✅ All fade curves (Linear, Bezier, Logarithmic)
- ✅ Trim, fade, and pad operations
- ✅ Combined operations in sequence
- ✅ Error conditions (missing files, invalid options, same input/output)
- ✅ Overwrite protection with --force flag
- ✅ Verbosity levels (quiet, normal, verbose)
- ✅ Progress tracking
- ✅ Validation of invalid settings

### Known Limitations

- Maximum audio data size: 1GB (safety limit)
- Maximum duration: 1,000,000 milliseconds (~16.7 minutes)
- Windows still uses seek/ftell for file size (no fstat available)

## [0.94beta] - 2026-01-06

### Major Refactoring Release

This release represents a complete refactoring of the codebase following modern C best practices while maintaining 100% backward compatibility.

### Added
- Cross-platform Makefile with support for macOS, Linux, and Windows
- Comprehensive error handling for all file I/O operations
- Input validation using `strtol()`/`strtod()` instead of `atoi()`/`atof()`
- Platform-independent path handling
- Extensive function documentation with doxygen-style comments
- Resource cleanup on all error paths
- Support for building with GCC, Clang, and MSVC

### Changed
- **Code organization**: Split into logical sections with helper functions
- **Naming convention**: Changed from camelCase to snake_case (C standard style)
- **Type safety**: Use `stdint.h` types (`int32_t`, `int16_t`, etc.)
- **Constants**: Replaced all magic numbers with named constants
- **Global variables**: Eliminated all globals, using context structure instead
- **Main function**: Reduced from 182 lines to 35 lines
- **Error messages**: More descriptive with proper file paths and ranges

### Fixed
- Memory leaks in path expansion on error paths
- Missing 32-bit audio support in threshold detection
- Incorrect sample alignment calculations
- Buffer overflow potential in memcpy operations
- Missing error checking for `fread()`, `fwrite()`, `fseek()`
- `argc` validation bug (was `< 2`, now correctly `< 3`)
- File handle leaks on error conditions
- Platform-specific compilation issues on macOS and Linux

### Removed
- Code duplication in sample reading/writing (now uses helpers)
- Unnecessary switch statement repetition
- Global state dependencies

### Technical Improvements
- Zero global variables (only static const data)
- Functions under 80 lines each (most under 50)
- Comprehensive documentation for all functions
- Clear separation of concerns
- Proper const-correctness
- goto-based error handling for cleanup

## [0.93beta] - 2020-11-10

### Changed
- More code cleanup and refactoring

## [0.92beta] - 2020-11-10

### Added
- Support for 8-bit and 32-bit WAV files
- Better error checking throughout

### Changed
- General code cleanup

## [0.91beta] - 2020-11-08

### Added
- Initial release
- Trim silence from WAV file ends
- Fade in/out with multiple curve types (linear, Bezier, logarithmic)
- Pad with silence at start/end
- Support for 16-bit and 24-bit uncompressed WAV files
- Mono and stereo support
