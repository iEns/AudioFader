# AudioFader

#### Simple command line tool for fading and trimming of WAV files
AudioFader will trim silence from the ends of a WAV file, fade in/out and pad the ends with silence.
Fades can be done linearly, logarithmically or with a Bezier curve.

Works with 8, 16, 24 and 32 bit uncompressed stereo and mono WAV files.

### Platform Support

AudioFader is cross-platform and works on:
- **Windows** (tested with GCC/MinGW and MSVC)
- **macOS** (tested on Apple Silicon and Intel)
- **Linux** (any modern distribution with GCC)

Note: Code assumes little-endian architecture (not tested on big-endian systems).

### Project Structure

```
AudioFader/
├── audiofader.c         # Main entry point
├── common.h             # Shared types, constants, macros
├── wav_io.h/c           # WAV file reading/writing, validation
├── audio_processing.h/c # Trim, fade, pad functions
├── options.h/c          # CLI parsing, validation
└── ui.h/c               # Progress reporting, cleanup
```

### Building

#### Using Make (Recommended)

A cross-platform Makefile is provided for easy building:

```bash
make              # Build release version
make debug        # Build with debug symbols
make clean        # Clean build artifacts
make install      # Install to system (may require sudo)
make help         # Show all available targets
```

**Platform-specific builds:**
```bash
make CC=gcc       # Use GCC (default)
make CC=clang     # Use Clang
```

#### Manual Compilation

If you prefer to compile manually:

**On macOS/Linux:**
```bash
gcc -std=c99 -Wall -O2 -o audiofader audiofader.c wav_io.c audio_processing.c options.c ui.c -lm
```

**On Windows:**
```bash
gcc -std=c99 -Wall -O2 -o audiofader.exe audiofader.c wav_io.c audio_processing.c options.c ui.c
```

### Usage

```
AudioFader v0.96beta Copyright 2021-2026 iEns Labs

Trim, fade and pad WAV audio files.

USAGE:
  audiofader <input.wav> <output.wav> [OPTIONS]

AUDIO PROCESSING OPTIONS:
  Trimming (removes silence from ends):
    --trimstart N     Trim silence from start (0-100% threshold, default: 0)
    --trimend N       Trim silence from end (0-100% threshold, default: 0)
                      Common values: 0.01 to 1.0

  Fading:
    --fadein MS       Fade-in duration in milliseconds
    --fadeout MS      Fade-out duration in milliseconds
    --fadecurve N     Fade curve: 0=Linear (default), 1=Bezier, 2=Logarithmic

  Padding (adds silence):
    --padstart MS     Add silence at start (milliseconds)
    --padend MS       Add silence at end (milliseconds)

  Other:
    --dry-run         Preview mode - show what would be done without writing output
    --force           Overwrite existing output file without prompting
    --quiet           Quiet mode (errors only)
    --verbose         Verbose output (debug information)
    --help            Show this help message

Operations are always done in this order: Trim -> Fade -> Pad
WAV files must be uncompressed 8, 16, 24 or 32 bit, stereo or mono
Sample rates supported: 8000-192000 Hz
```

### Examples

**Basic fade:**
```bash
./audiofader input.wav output.wav --fadein 500 --fadeout 500
```

**Trim silence and fade with Bezier curve:**
```bash
./audiofader input.wav output.wav --trimstart 0.01 --fadein 1000 --fadeout 1000 --fadecurve 1
```

**Add 1 second padding:**
```bash
./audiofader input.wav output.wav --padstart 1000 --padend 1000
```

**Preview changes without writing (dry-run):**
```bash
./audiofader input.wav output.wav --trimstart 0.01 --dry-run
```

**Complete workflow:**
```bash
./audiofader input.wav output.wav --trimstart 0.01 --trimend 0.01 \
    --fadein 2000 --fadeout 2000 --fadecurve 1 \
    --padstart 500 --padend 500
```

### License

MIT License - see source files for details.
