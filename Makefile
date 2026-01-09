# =============================================================================
# AudioFader Makefile
# Cross-platform build system for macOS, Linux, and Windows
# =============================================================================

# Project information
PROJECT = audiofader
VERSION = 1.0.0

# Source files (multi-file structure)
SRCS = audiofader.c wav_io.c audio_processing.c options.c ui.c
HDRS = common.h wav_io.h audio_processing.h options.h ui.h

# Platform detection
UNAME_S := $(shell uname -s 2>/dev/null || echo Windows)
UNAME_M := $(shell uname -m 2>/dev/null || echo unknown)

# =============================================================================
# Compiler Configuration
# =============================================================================

# Default compiler
CC ?= gcc

# Common compiler flags
CFLAGS_COMMON = -std=c99 -Wall -Wextra -Wpedantic
CFLAGS_RELEASE = -O2 -DNDEBUG
CFLAGS_DEBUG = -O0 -g -DDEBUG

# Platform-specific settings
ifeq ($(UNAME_S),Darwin)
    # macOS
    PLATFORM = macOS
    TARGET = $(PROJECT)
    CFLAGS_PLATFORM =
    LDFLAGS = -lm
    INSTALL_DIR = /usr/local/bin
    RM = rm -f
    MKDIR = mkdir -p
endif

ifeq ($(UNAME_S),Linux)
    # Linux
    PLATFORM = Linux
    TARGET = $(PROJECT)
    CFLAGS_PLATFORM = -D_POSIX_C_SOURCE=200809L
    LDFLAGS = -lm
    INSTALL_DIR = /usr/local/bin
    RM = rm -f
    MKDIR = mkdir -p
endif

ifeq ($(UNAME_S),Windows)
    # Windows (MinGW)
    PLATFORM = Windows
    TARGET = $(PROJECT).exe
    CFLAGS_PLATFORM = -D_WIN32
    LDFLAGS =
    INSTALL_DIR = C:/Program Files/AudioFader
    RM = del /Q
    MKDIR = mkdir
endif

# Build type (default: release)
BUILD_TYPE ?= release

ifeq ($(BUILD_TYPE),debug)
    CFLAGS = $(CFLAGS_COMMON) $(CFLAGS_DEBUG) $(CFLAGS_PLATFORM)
    BUILD_DIR = build/debug
else
    CFLAGS = $(CFLAGS_COMMON) $(CFLAGS_RELEASE) $(CFLAGS_PLATFORM)
    BUILD_DIR = build/release
endif

# =============================================================================
# Targets
# =============================================================================

.PHONY: all clean debug release info install uninstall help

# Default target
all: info $(TARGET)

# Build the executable
$(TARGET): $(SRCS) $(HDRS)
	@echo "Building $(PROJECT) v$(VERSION) for $(PLATFORM) ($(BUILD_TYPE))..."
	$(CC) $(CFLAGS) -o $(TARGET) $(SRCS) $(LDFLAGS)
	@echo "Build complete: $(TARGET)"

# Debug build
debug:
	@$(MAKE) BUILD_TYPE=debug

# Release build (explicit)
release:
	@$(MAKE) BUILD_TYPE=release

# Display build information
info:
	@echo "==================================================================="
	@echo "  AudioFader Build System"
	@echo "==================================================================="
	@echo "  Version:      $(VERSION)"
	@echo "  Platform:     $(PLATFORM) ($(UNAME_M))"
	@echo "  Compiler:     $(CC)"
	@echo "  Build Type:   $(BUILD_TYPE)"
	@echo "  Target:       $(TARGET)"
	@echo "  Source Files: $(SRCS)"
	@echo "  CFLAGS:       $(CFLAGS)"
	@echo "  LDFLAGS:      $(LDFLAGS)"
	@echo "==================================================================="
	@echo ""

# Clean build artifacts
clean:
	@echo "Cleaning build artifacts..."
	-$(RM) $(PROJECT) $(PROJECT).exe 2>/dev/null || true
	-$(RM) *.o *.obj 2>/dev/null || true
	-$(RM) -r build 2>/dev/null || true
	@echo "Clean complete."

# Install to system
install: $(TARGET)
	@echo "Installing $(TARGET) to $(INSTALL_DIR)..."
ifeq ($(UNAME_S),Windows)
	@if not exist "$(INSTALL_DIR)" $(MKDIR) "$(INSTALL_DIR)"
	copy $(TARGET) "$(INSTALL_DIR)"
else
	$(MKDIR) $(INSTALL_DIR)
	install -m 755 $(TARGET) $(INSTALL_DIR)/$(TARGET)
endif
	@echo "Installation complete."
	@echo "Run 'audiofader' from anywhere to use the tool."

# Uninstall from system
uninstall:
	@echo "Uninstalling $(PROJECT)..."
ifeq ($(UNAME_S),Windows)
	-$(RM) "$(INSTALL_DIR)\$(TARGET)" 2>/dev/null || true
else
	-$(RM) $(INSTALL_DIR)/$(TARGET)
endif
	@echo "Uninstall complete."

# Static analysis (if available)
analyze:
	@echo "Running static analysis..."
	@which cppcheck >/dev/null 2>&1 && cppcheck --enable=all --suppress=missingIncludeSystem $(SRCS) || echo "cppcheck not found, skipping..."

# Format code (if clang-format available)
format:
	@echo "Formatting code..."
	@which clang-format >/dev/null 2>&1 && clang-format -i $(SRCS) $(HDRS) || echo "clang-format not found, skipping..."

# Run with example (requires test file)
test: $(TARGET)
	@echo "Running basic functionality test..."
	@if [ -f "test16sine.wav" ]; then \
		./$(TARGET) test16sine.wav output_test.wav --fadein 500 --fadeout 500 --dry-run; \
	else \
		echo "No test16sine.wav file found. Create one to test functionality."; \
	fi

# Help message
help:
	@echo "AudioFader Makefile - Available targets:"
	@echo ""
	@echo "  make              - Build release version (default)"
	@echo "  make all          - Same as 'make'"
	@echo "  make release      - Build optimized release version"
	@echo "  make debug        - Build debug version with symbols"
	@echo "  make clean        - Remove build artifacts"
	@echo "  make install      - Install to system (may require sudo)"
	@echo "  make uninstall    - Remove from system (may require sudo)"
	@echo "  make info         - Display build configuration"
	@echo "  make analyze      - Run static analysis (requires cppcheck)"
	@echo "  make format       - Format source code (requires clang-format)"
	@echo "  make test         - Run basic functionality test"
	@echo "  make help         - Display this help message"
	@echo ""
	@echo "Build options:"
	@echo "  CC=<compiler>     - Set compiler (default: gcc)"
	@echo "  BUILD_TYPE=debug  - Build with debug symbols"
	@echo ""
	@echo "Examples:"
	@echo "  make                          # Build release version"
	@echo "  make debug                    # Build debug version"
	@echo "  make CC=clang                 # Use clang compiler"
	@echo "  make clean && make            # Clean rebuild"
	@echo "  sudo make install             # Install to /usr/local/bin"
	@echo ""
