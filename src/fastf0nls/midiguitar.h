/// Pitch tracker.
#pragma once
#include <stdint.h>

/// FFT static dispatch based on size.
#define FFT_INIT arm_rfft_fast_init_1024_f32

/// FFT size.
enum { N_FFT_GRID = 1024 };

/// Downsampling factor.
enum { LOG_SAMPLE_DIVISOR = 0 };

/// Sample chunk size.
enum { AUDIO_CAP = 1024 };

/// Number of samples after downsampling.
enum { SAMPLES = AUDIO_CAP >> LOG_SAMPLE_DIVISOR };

/// Midi TX buffer size.
enum { MIDI_CAP = 16 };

/// Max number of harmonics.
enum { MAX_MODEL_ORDER = 3 };

/// Midpoint.
enum { MP = (N_FFT_GRID >> 1) + 1 };

/// FFT bounds.
enum { MIN_FFT_INDEX = 1, MAX_FFT_INDEX = N_FFT_GRID >> 1 };

/// Logarithm of ADC midpoint.
enum { LOG_OFFSET = 15 };

/// ADC midpoint.
enum { OFFSET = 1 << LOG_OFFSET };

/// Pitch tracking state.
struct midiguitar {
  /// Most recent note.
  uint8_t note;

  /// Bend for most recent note.
  uint16_t bend;

  /// Samples from previous calls, in case of overlap.
  float input[SAMPLES];
};

/// Converts audio to midi.
uint8_t midiguitar(struct midiguitar *midiguitar,
                   const uint16_t input[AUDIO_CAP], uint8_t output[MIDI_CAP]);
