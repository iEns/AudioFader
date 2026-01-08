# Contributing to AudioFader

Thank you for your interest in contributing to AudioFader! This document provides guidelines and instructions for contributing.

## Code of Conduct

This project adheres to a Code of Conduct that all contributors are expected to follow. Please read [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md) before contributing.

## How to Contribute

### Reporting Bugs

Before creating a bug report, please check existing issues to avoid duplicates.

When filing a bug report, include:
- **Clear title and description**
- **Steps to reproduce** the issue
- **Expected behavior** vs. **actual behavior**
- **Platform information** (OS, compiler, version)
- **Sample WAV file** if relevant (or describe audio characteristics)
- **Command-line arguments** used
- **Error messages** (full output)

### Suggesting Features

Feature requests are welcome! Please:
- **Check existing issues** to avoid duplicates
- **Describe the use case** and why it would be valuable
- **Provide examples** of how the feature would be used
- **Consider scope** - keep features focused on audio processing

### Pull Requests

We actively welcome pull requests!

#### Before Starting

1. **Check existing issues** - someone may already be working on it
2. **Open an issue** to discuss major changes before implementing
3. **Fork the repository** and create a branch for your work

#### Development Setup

```bash
# Clone your fork
git clone https://github.com/iEns/AudioFader.git
cd AudioFader

# Build and test
make
./audiofader --help

# Run with test file
./audiofader input.wav output.wav --fadein 1000 --fadeout 1000
```

#### Code Standards

**Style Guide:**
- Follow C99 standard
- Use `snake_case` for all identifiers
- Use `_t` suffix for type names
- Maximum line length: 100 characters
- Use 4 spaces for indentation (no tabs)
- Format code with: `make format` (requires clang-format)

**Best Practices:**
- **No global variables** - pass state through function parameters
- **Use named constants** - no magic numbers
- **Check all I/O operations** - validate return values
- **Document functions** - include purpose, parameters, and return values
- **Keep functions small** - under 80 lines, ideally under 50
- **Use stdint.h types** - `int32_t`, `int16_t`, etc.
- **Handle errors** - clean up resources on all error paths

**Memory Management:**
- Always free allocated memory
- Use `cleanup_and_exit()` pattern for resource cleanup
- Check `malloc()`/`calloc()` return values
- Avoid memory leaks - test with valgrind or similar tools

#### Testing

Before submitting:

```bash
# Clean build
make clean && make

# Build with different compilers
make CC=gcc
make CC=clang

# Test debug build
make debug

# Run static analysis (if available)
make analyze

# Test with various WAV files
./audiofader test_8bit.wav out.wav --fadein 500
./audiofader test_16bit.wav out.wav --trimstart 0.01
./audiofader test_24bit.wav out.wav --fadeout 1000
./audiofader test_32bit.wav out.wav --padstart 1000
```

**Test coverage should include:**
- 8, 16, 24, and 32-bit audio
- Mono and stereo files
- All fade curve types
- Edge cases (very short files, very long fade times)
- Invalid inputs (error handling)

#### Commit Messages

Write clear, descriptive commit messages:

```
Short (50 chars or less) summary

More detailed explanation if needed. Wrap at 72 characters.
Explain the problem this commit solves and why you chose
this solution.

- Bullet points are fine
- Use present tense ("Add feature" not "Added feature")
- Reference issues: "Fixes #123" or "Relates to #456"
```

Examples:
```
Fix memory leak in path expansion error handling

Add support for 48kHz sample rate validation

Refactor sample I/O into helper functions
```

#### Pull Request Process

1. **Update documentation** - README, CHANGELOG, code comments
2. **Ensure all tests pass** - build on multiple platforms
3. **Update CHANGELOG.md** - add your changes under "Unreleased"
4. **Create pull request** with:
   - Clear title describing the change
   - Reference to related issue(s)
   - Description of what changed and why
   - Any breaking changes or migration notes
   - Test results from different platforms

5. **Respond to review feedback** - be patient and constructive

#### What We Look For

✅ **Good Pull Requests:**
- Solve a specific problem
- Follow the existing code style
- Include tests/verification
- Are well-documented
- Have clean commit history

❌ **Avoid:**
- Mixing unrelated changes
- Breaking existing functionality
- Adding unnecessary dependencies
- Changing code style throughout (do this separately)
- Incomplete implementations

## Development Workflow

### Setting Up Development Environment

**macOS:**
```bash
# Install Xcode Command Line Tools
xcode-select --install

# Or install GCC via Homebrew
brew install gcc
```

**Linux:**
```bash
# Debian/Ubuntu
sudo apt-get install build-essential

# Fedora/RHEL
sudo dnf install gcc make
```

**Windows:**
```bash
# Install MinGW-w64 or use MSVC
```

### Building Different Configurations

```bash
# Release build (optimized)
make release

# Debug build (with symbols)
make debug

# With specific compiler
make CC=clang

# Clean rebuild
make clean && make

# Install to system
sudo make install

# Show build info
make info
```

## Questions?

- **Open an issue** for questions about contributing
- **Check existing issues** for similar questions
- **Be respectful** and patient - this is a volunteer project

## License

By contributing, you agree that your contributions will be licensed under the MIT License.

Thank you for contributing to AudioFader! 🎵
