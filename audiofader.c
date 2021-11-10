/*

Audio Fader, Trimmer and Padder

Concept, design and programming by Jens Sandalgaard

https://github.com/iens/AudioFader

Copyright (c) 2021 Jens Sandalgaard

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
0.91beta	2020-11-08	Initial version
0.92beta	2020-11-10	Code Cleanup, better error checking, support for 8 and 32 bit files.
0.93beta	2020-11-10	More code cleanup.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <windows.h>

// define wav file header
typedef struct {
    char chunkID[4];
    int chunkSize;
    char format[4];
    char subchunk1ID[4];
    int subchunk1Size;
    short int audioFormat;
    short int numChannels;
    int sampleRate;
    int byteRate;
    short int blockAlign;
    short int bitsPerSample;
    char subchunk2ID[4];
    int subchunk2Size;
} WAV_HEADER;

// define subChunkHeader
typedef struct {
    char subchunkID[4];
    int subchunkSize;
} SUB_CHUNK_HEADER;

// define fadeCurveNames
char *fadeCurveNames[] = {"Linear", "Bezier", "Logarithmic"};

// define struct to hold the options
typedef struct {
    char *inputFileName;
    char *outputFileName;
    int fadeCurve;
    int fadeIn;
    int fadeOut;
    double trimStart;
    double trimEnd;
    int padStart;
    int padEnd;
} OPTIONS;

// initialize the options
OPTIONS options = {NULL, NULL, 0, 0, 0, 0, 0.0, 0.0, 0};

// struct for basic 24 bit compatibility
struct int24 {
    signed int24 : 24;
} __attribute__((packed));
struct int24 inputSample24;

void printTitle() { printf("\nAudioFader v0.93beta Copyright 2021 iEns Labs\n\n"); }

void printHelp(unsigned char *programName) {
    printf("Trim, fade and pad WAV-files.\n\n");
    printf("Usage: %s <input file> <output file> <options>\n\n", programName);
    printf("Options:\n");
    printf("\t--trimstart n\tTrim start\n");
    printf("\t--trimend n\tTrim end\n");
    printf("\t\t\tArgument for any trim option is a percentage threshold\n");
    printf("\t\t\tn=0 means no trimming\n");
    printf("\t\t\tn=0.0001 to 0.01 are commonly used values\n");
    printf("\t--fadein n\tFade-in duration in milliseconds\n");
    printf("\t--fadeout n\tFade-out duration in milliseconds\n");
    printf("\t--fadecurve n\tChoose the applied curve to the fade\n");
    printf("\t\t\t0 = Linear (default)\n");
    printf("\t\t\t1 = Bezier\n");
    printf("\t\t\t2 = Logarithmic\n");
    printf("\t--padstart n\tPad the start with silence, in milliseconds\n");
    printf("\t\t\tPadding is done after trimming and fading\n");
    printf("\t--padend n\tPad the end with silence, in milliseconds\n");
    printf("\t\t\tPadding is done after trimming and fading\n");
    printf("\n");
    printf("Operations are always done in this order: Trim->Fade->Pad\n");
    printf("WAV-files must be uncompressed 8, 16 , 24 or 32 bit, stereo or mono\n\n");
    printf("Example: %s input.wav output.wav --trimstart 0.01 --fadein 1000 --fadeout 1000 --padstart 1000 --padend 1000\n\n", programName);
}

int getAndCheckOptions(int argc, char *argv[]) {
    // check if the number of arguments is correct
    if (argc < 2) {
        // get the name of the program
        char *programName = strrchr(argv[0], '\\');
        if (programName == NULL) {
            programName = argv[0];
        } else {
            programName++;  // skip the last backslash
        }

        printHelp(programName);
        return 1;
    }

    options.inputFileName = argv[1];
    options.outputFileName = argv[2];

    // expand input file name to full path
    char *inputFileNameFullPath = (char *)malloc(MAX_PATH);
    GetFullPathName(options.inputFileName, MAX_PATH, inputFileNameFullPath, NULL);
    options.inputFileName = inputFileNameFullPath;

    // expand output file name to full path
    char *outputFileNameFullPath = (char *)malloc(MAX_PATH);
    GetFullPathName(options.outputFileName, MAX_PATH, outputFileNameFullPath, NULL);
    options.outputFileName = outputFileNameFullPath;

    // check that input file and output file are not the same (case insensitive)
    if (strcasecmp(options.inputFileName, options.outputFileName) == 0) {
        printf("Input and output file names cannot be the same\n\n");
        return 1;
    }

    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "--fadein") == 0) {
            // check that an argument is supplied and is not the next option
            if (i + 1 >= argc || argv[i + 1][0] == '-') {
                printf("Error: Option --fadein requires an argument\n");
                return 1;
            }
            options.fadeIn = atoi(argv[i + 1]);
            i++;
        } else if (strcmp(argv[i], "--fadeout") == 0) {
            if (i + 1 >= argc || argv[i + 1][0] == '-') {
                printf("Error: Option --fadeout requires an argument\n");
                return 1;
            }
            options.fadeOut = atoi(argv[i + 1]);
            i++;
        } else if (strcmp(argv[i], "--fadecurve") == 0) {
            if (i + 1 >= argc || argv[i + 1][0] == '-') {
                printf("Error: Option --fadecurve requires an argument\n");
                return 1;
            }
            options.fadeCurve = atoi(argv[i + 1]);
            i++;
        } else if (strcmp(argv[i], "--trimstart") == 0) {
            if (i + 1 >= argc || argv[i + 1][0] == '-') {
                printf("Error: Option --trimstart requires an argument\n");
                return 1;
            }
            options.trimStart = atof(argv[i + 1]);
            i++;
        } else if (strcmp(argv[i], "--trimend") == 0) {
            if (i + 1 >= argc || argv[i + 1][0] == '-') {
                printf("Error: Option --trimend requires an argument\n");
                return 1;
            }
            options.trimEnd = atof(argv[i + 1]);
            i++;
        } else if (strcmp(argv[i], "--padstart") == 0) {
            if (i + 1 >= argc || argv[i + 1][0] == '-') {
                printf("Error: Option --padstart requires an argument\n");
                return 1;
            }
            options.padStart = atoi(argv[i + 1]);
            i++;
        } else if (strcmp(argv[i], "--padend") == 0) {
            if (i + 1 >= argc || argv[i + 1][0] == '-') {
                printf("Error: Option --padend requires an argument\n");
                return 1;
            }
            options.padEnd = atoi(argv[i + 1]);
            i++;
        } else {
            printf("Unknown option: %s\n\n", argv[i]);
            return 1;
        }
    }

    // check for valid fade curve
    if (options.fadeCurve < 0 || options.fadeCurve > 2) {
        printf("Invalid fade curve: %d\n\n", options.fadeCurve);
        return 1;
    }

    // check for valid trim start
    if (options.trimStart < 0 || options.trimStart > 100) {
        printf("Invalid trim start: %f\n\n", options.trimStart);
        return 1;
    }

    // check for valid trim end
    if (options.trimEnd < 0 || options.trimEnd > 100) {
        printf("Invalid trim end: %f\n\n", options.trimEnd);
        return 1;
    }

    // check for valid pad start
    if (options.padStart < 0 || options.padStart > 1000000) {
        printf("Invalid pad start: %d\n\n", options.padStart);
        return 1;
    }

    // check for valid pad end
    if (options.padEnd < 0 || options.padEnd > 1000000) {
        printf("Invalid pad end: %d\n\n", options.padEnd);
        return 1;
    }

    // check for valid fade in
    if (options.fadeIn < 0 || options.fadeIn > 1000000) {
        printf("Invalid fade in: %d\n\n", options.fadeIn);
        return 1;
    }

    // check for valid fade out
    if (options.fadeOut < 0 || options.fadeOut > 1000000) {
        printf("Invalid fade out: %d\n\n", options.fadeOut);
        return 1;
    }
    return 0;
}

void printOptionsUsed() {
    // print all options
    printf("Options in use:\n");
    printf("\tInput file: \t\t%s\n", options.inputFileName);
    printf("\tOutput file: \t\t%s\n", options.outputFileName);
    if (options.fadeIn != 0) {
        printf("\tFade-in duration: \t%d ms\n", options.fadeIn);
    }
    if (options.fadeOut != 0) {
        printf("\tFade-out duration: \t%d ms\n", options.fadeOut);
    }
    if (options.fadeIn != 0 || options.fadeOut != 0) {
        printf("\tFade-curve: \t\t%s\n", fadeCurveNames[options.fadeCurve]);
    }
    if (options.trimStart != 0) {
        printf("\tTrim start: \t\t%f%%\n", options.trimStart);
    }
    if (options.trimEnd != 0) {
        printf("\tTrim end: \t\t%f%%\n", options.trimEnd);
    }
    if (options.padStart != 0) {
        printf("\tPad start: \t\t%d ms\n", options.padStart);
    }
    if (options.padEnd != 0) {
        printf("\tPad end: \t\t%d ms\n", options.padEnd);
    }
    printf("\n");
}

void printHeaderStats(WAV_HEADER header) {
    // print numchannels
    printf("Input File Information:\n");
    printf("\tChannels: \t\t%d\n", header.numChannels);
    printf("\tBits Per Sample: \t%d\n", header.bitsPerSample);
    printf("\tSample Rate: \t\t%d\n", header.sampleRate);
    printf("\tAudioFormat: \t\t%d\n", header.audioFormat);
    printf("\tBlockAlign: \t\t%d\n", header.blockAlign);
    printf("\n");
}
// convert duration to samples, aligned to sample rate
int durationToSamples(int duration, int sampleRate) { return (duration * sampleRate) / 1000; }

// buffer length in bytes, fade duration in milliseconds, direction: 0 = in, 1 = out
int fadeBuffer(WAV_HEADER header, unsigned char *buffer, int bufferSize, int fadeDuration, int direction, int fadeCurve) {
    int fadeSamples = durationToSamples(fadeDuration, header.sampleRate);

    int bytesPerFrame = header.numChannels * header.bitsPerSample / 8;
    int bytesPerSample = header.bitsPerSample / 8;

    int fadeInSample = 0;
    int inputSample;

    // fade buffer
    int fadeStart;
    int fadeEnd;

    if (direction == 0) {
        fadeStart = 0;
        fadeEnd = fadeSamples * header.numChannels;
    } else {
        fadeStart = bufferSize / bytesPerSample - fadeSamples * header.numChannels;
        fadeEnd = bufferSize / bytesPerSample;
    }

    for (int i = fadeStart; i < fadeEnd; i++) {
        // calculate fadein factor
        double fadeFactor;
        double fadeIndex;
        if (direction == 0) {
            fadeIndex = (double)(i - fadeStart) / (double)(fadeEnd - fadeStart);
        } else {
            fadeIndex = 1.0 - (double)(i - fadeStart) / (double)(fadeEnd - fadeStart);
        }
        switch (fadeCurve) {
            case 0:
                // Linear Curve
                if (direction == 0) {
                    fadeFactor = (double)i / (double)(fadeEnd);
                } else {
                    fadeFactor = 1.0 - ((double)i - fadeStart) / (double)(fadeEnd - fadeStart);
                }
                break;
            case 1:
                // Bezier curve
                fadeFactor = fadeIndex * fadeIndex * (3 - 2 * fadeIndex);
                break;
            case 2:
                // Logarithmic curve
                fadeFactor = (fadeIndex / (1 + (1 - fadeIndex)));
                break;
        }

        switch (bytesPerSample) {
            case 1:
                inputSample = (buffer[i * bytesPerSample]) - 127;
                fadeInSample = (inputSample * fadeFactor);
                buffer[i * bytesPerSample] = (unsigned char)(fadeInSample + 127);
                break;
            case 2:
                inputSample = (signed short int)(buffer[i * bytesPerSample + 1] << 8 | buffer[i * bytesPerSample]);
                fadeInSample = (inputSample * fadeFactor);
                buffer[i * bytesPerSample + 1] = (unsigned char)(fadeInSample >> 8);
                buffer[i * bytesPerSample] = (unsigned char)(fadeInSample & 0xff);
                break;
            case 3:
                inputSample24.int24 = ((buffer[i * bytesPerSample + 2] << 16) | (buffer[i * bytesPerSample + 1] << 8) | buffer[i * bytesPerSample]);
                inputSample = inputSample24.int24;
                fadeInSample = (inputSample * fadeFactor);

                inputSample24.int24 = (signed int)(fadeInSample);
                buffer[i * bytesPerSample + 2] = (unsigned char)(inputSample24.int24 >> 16);
                buffer[i * bytesPerSample + 1] = (unsigned char)(inputSample24.int24 >> 8);
                buffer[i * bytesPerSample] = (unsigned char)(inputSample24.int24 & 0xff);
                break;
            case 4:
                memcpy(&inputSample, &buffer[i * bytesPerSample], 4);
                inputSample = (inputSample * fadeFactor);
                memcpy(&buffer[i * bytesPerSample], &inputSample, 4);
                break;
        }
    }
}

int findFirstSampleOverThreshold(WAV_HEADER header, unsigned char *buffer, int bufferSize, double thresholdPercent) {
    int bytesPerFrame = header.numChannels * header.bitsPerSample / 8;
    int bytesPerSamplePerChannel = header.bitsPerSample / 8;

    int inputSample;
    int firstSample = 0;
    int threshold = 0;

    // set threshold to thresholdpercent of max value
    switch (bytesPerSamplePerChannel) {
        case 1:
            threshold = (int)((double)0x7f * (double)thresholdPercent / 100.0);
            break;
        case 2:
            threshold = (int)((double)0x7fff * (double)thresholdPercent / 100.0);
            break;
        case 3:
            threshold = (int)((double)0x7fffff * (double)thresholdPercent / 100.0);
            break;
    }

    for (int i = 0; i < bufferSize / bytesPerSamplePerChannel; i++) {
        switch (bytesPerSamplePerChannel) {
            case 1:
                inputSample = buffer[i * bytesPerSamplePerChannel];
                break;
            case 2:
                inputSample = (signed short int)(buffer[i * bytesPerSamplePerChannel + 1] << 8 | buffer[i * bytesPerSamplePerChannel]);
                break;
            case 3:
                inputSample24.int24 = ((buffer[i * bytesPerSamplePerChannel + 2] << 16) | (buffer[i * bytesPerSamplePerChannel + 1] << 8) | buffer[i * bytesPerSamplePerChannel]);
                inputSample = inputSample24.int24;
                break;
        }
        if (inputSample > threshold) {
            firstSample = i;
            // align to number of channels
            firstSample = firstSample - (firstSample % (header.numChannels * (header.bitsPerSample / 8)));
            break;
        }
    }

    return firstSample;
}

int findLastSampleOverThreshold(WAV_HEADER header, unsigned char *buffer, int bufferSize, double thresholdPercent) {
    int bytesPerFrame = header.numChannels * header.bitsPerSample / 8;
    int bytesPerSamplePerChannel = header.bitsPerSample / 8;

    int inputSample;
    int lastSample = 0;
    int threshold = 0;

    // set threshold to thresholdpercent of max value
    switch (bytesPerSamplePerChannel) {
        case 1:
            threshold = (int)((double)0x7f * (double)thresholdPercent / 100.0);
            break;
        case 2:
            threshold = (int)((double)0x7fff * (double)thresholdPercent / 100.0);
            break;
        case 3:
            threshold = (int)((double)0x7fffff * (double)thresholdPercent / 100.0);
            break;
    }

    for (int i = bufferSize / bytesPerSamplePerChannel - 1; i >= 0; i--) {
        switch (bytesPerSamplePerChannel) {
            case 1:
                inputSample = buffer[i * bytesPerSamplePerChannel];
                break;
            case 2:
                inputSample = (signed short int)(buffer[i * bytesPerSamplePerChannel + 1] << 8 | buffer[i * bytesPerSamplePerChannel]);
                break;
            case 3:
                inputSample24.int24 = ((buffer[i * bytesPerSamplePerChannel + 2] << 16) | (buffer[i * bytesPerSamplePerChannel + 1] << 8) | buffer[i * bytesPerSamplePerChannel]);
                inputSample = inputSample24.int24;
                break;
        }

        if (inputSample > threshold) {
            lastSample = i;
            // align to number of channels
            lastSample = lastSample - (lastSample % (header.numChannels * (header.bitsPerSample / 8)));
            break;
        }
    }

    return lastSample;
}

int main(int argc, char *argv[]) {
    // start a timer to measure the time it takes to process the file
    clock_t startTime = clock();

    printTitle();

    int retVal = getAndCheckOptions(argc, argv);
    if (retVal != 0) {
        return retVal;
    }

    // open input file
    FILE *inputFile = fopen(options.inputFileName, "rb");
    if (inputFile == NULL) {
        printf("Error: cannot open input file %s\n", argv[1]);
        return 1;
    }

    printOptionsUsed();

    // open output file
    FILE *outputFile = fopen(options.outputFileName, "wb");
    if (outputFile == NULL) {
        printf("Error: cannot open output file %s\n", argv[2]);
        return 1;
    }

    // read input file header
    WAV_HEADER inputHeader;
    SUB_CHUNK_HEADER inputSubChunkHeader;

    fread(&inputHeader, sizeof(WAV_HEADER), 1, inputFile);

    int nextChunkSize = inputHeader.subchunk2Size;

    // store current seek position, then get the filesize of open file called inputFile and then seek back to current position
    int currentPosition = ftell(inputFile);
    fseek(inputFile, 0, SEEK_END);
    int inputFileSize = ftell(inputFile);
    fseek(inputFile, currentPosition, SEEK_SET);

    // check if the input file is a valid wav file
    if (strncmp(inputHeader.chunkID, "RIFF", 4) != 0) {
        printf("Error: input file is not a valid wav file\n");
        return 1;
    }
    if (strncmp(inputHeader.format, "WAVE", 4) != 0) {
        printf("Error: input file is not a valid wav file\n");
        return 1;
    }
    if (strncmp(inputHeader.subchunk1ID, "fmt ", 4) != 0) {
        printf("Error: input file is not a valid wav file\n");
        return 1;
    }

    // check that the file is 8, 16, 24 or 32 bits
    if (inputHeader.bitsPerSample != 8 && inputHeader.bitsPerSample != 16 && inputHeader.bitsPerSample != 24 && inputHeader.bitsPerSample != 32) {
        printf("Error: input file is not 8, 16, 24 or 32 bits\n");
        return 1;
    }

    // check that the file is mono or stereo
    if (inputHeader.numChannels != 1 && inputHeader.numChannels != 2) {
        printf("Error: input file is not mono or stereo\n");
        return 1;
    }

    // skip all non-data chunks
    // TODO: Do we need to preserve the other chunks?
    if (strncmp(inputHeader.subchunk2ID, "data", 4) != 0) {
        // inputHeader.subchunk2ID); 		printf("Skipping...\n");
        int pointerToNextSubchunk = inputHeader.subchunk2Size + sizeof(WAV_HEADER);
        while (1) {
            fseek(inputFile, pointerToNextSubchunk, SEEK_SET);
            fread(&inputSubChunkHeader, sizeof(SUB_CHUNK_HEADER), 1, inputFile);
            nextChunkSize = inputSubChunkHeader.subchunkSize;
            if (strncmp(inputSubChunkHeader.subchunkID, "data", 4) == 0) {
                break;
            }
            pointerToNextSubchunk += nextChunkSize + sizeof(SUB_CHUNK_HEADER);
            if (pointerToNextSubchunk >= inputFileSize || pointerToNextSubchunk < 0) {
                printf("Error: input file is corrupted\n");
                return 1;
            }
        }
    }
    int inputDataChunkSize = nextChunkSize;

    printHeaderStats(inputHeader);

    // convert input file length to number of samples
    int numSamples = inputDataChunkSize / inputHeader.numChannels / (inputHeader.bitsPerSample / 8);

    // create buffer for input file data
    unsigned char *inputData = (unsigned char *)malloc(inputDataChunkSize);
    if (inputData == NULL) {
        printf("Error: cannot allocate memory for input file data\n");
        return 1;
    }

    // read input file data into buffer
    fread(inputData, inputDataChunkSize, 1, inputFile);

    // if options.trimEnd is set, find the last sample to be trimmed
    int lastSample = numSamples;
    if (options.trimEnd != 0) {
        lastSample = findLastSampleOverThreshold(inputHeader, inputData, inputDataChunkSize, options.trimEnd);
    }

    // if options.trimStart is set, find the first sample to be trimmed
    int firstSample = 0;
    if (options.trimStart != 0) {
        firstSample = findFirstSampleOverThreshold(inputHeader, inputData, inputDataChunkSize, options.trimStart);
    }

    int trimmedBufferLength = (lastSample - firstSample) * (inputHeader.bitsPerSample / 8);

    // create buffer for output file data
    unsigned char *trimmedBuffer = (unsigned char *)malloc(trimmedBufferLength);
    if (trimmedBuffer == NULL) {
        printf("Error: cannot allocate memory for output file data\n");
        return 1;
    }

    // copy input file data to output file data
    memcpy(trimmedBuffer, inputData + (firstSample * (inputHeader.bitsPerSample / 8)), trimmedBufferLength);
    free(inputData);

    if (options.fadeIn > 0) {
        fadeBuffer(inputHeader, trimmedBuffer, trimmedBufferLength, options.fadeIn, 0, options.fadeCurve);
    }

    if (options.fadeOut > 0) {
        fadeBuffer(inputHeader, trimmedBuffer, trimmedBufferLength, options.fadeOut, 1, options.fadeCurve);
    }

    int padStartSamples = durationToSamples(options.padStart, inputHeader.sampleRate);
    int padStartBytes = padStartSamples * inputHeader.numChannels * inputHeader.bitsPerSample / 8;
    int padEndSamples = durationToSamples(options.padEnd, inputHeader.sampleRate);
    int padEndBytes = padEndSamples * inputHeader.numChannels * inputHeader.bitsPerSample / 8;

    int paddedBufferLength = trimmedBufferLength + padStartBytes + padEndBytes;

    // create buffer for output file data, initialized to 0
    unsigned char *paddedBuffer = (unsigned char *)calloc(1, paddedBufferLength);
    if (paddedBuffer == NULL) {
        printf("Error: cannot allocate memory for output file data\n");
        return 1;
    }

    // copy trimmed buffer to padded buffer
    memcpy(paddedBuffer + padStartBytes, trimmedBuffer, trimmedBufferLength);
    free(trimmedBuffer);

    // modify inputheader with new length
    inputHeader.subchunk2Size = paddedBufferLength;

    // modify inputheader with new chunksize
    inputHeader.chunkSize = inputHeader.subchunk2Size + 36;

    // force subchunk2ID to data
    strncpy(inputHeader.subchunk2ID, "data", 4);

    // write output file header
    fwrite(&inputHeader, sizeof(WAV_HEADER), 1, outputFile);

    // write output file data
    fwrite(paddedBuffer, paddedBufferLength, 1, outputFile);

    // close input and output files
    fclose(inputFile);
    fclose(outputFile);

    // free memory
    free(paddedBuffer);

    // print that the program has finished processing the file, and include the time it took
    printf("\nFinished processing file in %f seconds.\n\n", (double)(clock() - startTime) / CLOCKS_PER_SEC);

    return 0;
}
