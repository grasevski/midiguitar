/// Pitch tracker.
#pragma once
#include <stdint.h>

/// FFT static dispatch based on size.
#define FFT_INIT arm_rfft_fast_init_4096_f32

/// FFT size.
enum { N_FFT_GRID = 4096 };

/// Sample chunk size.
enum { SAMPLES = 512 };

/// Midi TX buffer size.
enum { MIDI_CAP = 16 };

/// Max number of harmonics.
enum { MAX_MODEL_ORDER = 4 };

/// Midpoint.
enum { MP = (N_FFT_GRID >> 1) + 1 };

/// FFT bounds.
enum { MIN_FFT_INDEX = 1, MAX_FFT_INDEX = N_FFT_GRID >> 1 };

/// Pitch tracking state.
struct midiguitar {
  /// Most recent note.
  uint8_t note;

  /// Bend for most recent note.
  uint16_t bend;
};

/// Converts audio to midi.
uint8_t midiguitar(struct midiguitar *midiguitar, const float input[SAMPLES],
                   uint8_t output[MIDI_CAP]);
