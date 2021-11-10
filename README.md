# AudioFader

#### Simple command line tool for fading and trimming of WAV files
Audiofader will trim silence from the ends of a WAV file, fade in/out and pad the ends with silence.
Fades can be done linearly, logarithmically or with a Bezier curve.

Works with 8, 16, 24 and 32 bit uncompressed stereo and mono WAV files.

Tested on Windows 11 (x64), compiled with GCC (MinGW32).
I have not tested on other architectures (Could potentially have problems with endianness).

### Usage:

	AudioFader v0.92beta Copyright 2021 iEns Labs

	Trim, fade and pad WAV-files.

	Usage: audiofader.exe <input file> <output file> <options>

	Options:
			--trim n        Trim start and end
			--trimstart n   Trim start
			--trimend n     Trim end
							Argument for any trim option is a percentage threshold
							n=0 means no trimming
							n=0.0001 to 0.01 are commonly used values
			--fadein n      Fade-in duration in milliseconds
			--fadeout n     Fade-out duration in milliseconds
			--fadecurve n   Choose the applied curve to the fade
							0 = Linear (default)
							1 = Bezier
							2 = Logarithmic
			--padstart n    Pad the start with silence, in milliseconds
							Padding is done after trimming and fading
			--padend n      Pad the end with silence, in milliseconds
							Padding is done after trimming and fading

	Operations are always done in this order: Trim->Fade->Pad
	WAV-files must be uncompressed 8, 16, 24 and 32 bit, stereo or mono

	Example: audiofader.exe input.wav output.wav --trim 0.01 --fadein 1000 --fadeout 1000 --padstart 1000 --padend 1000

### Compile:
	
	gcc audiofader.c -o audiofader.exe
