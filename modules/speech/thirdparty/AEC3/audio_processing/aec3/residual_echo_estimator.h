/*
 *  Copyright (c) 2017 The WebRTC project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef MODULES_AUDIO_PROCESSING_AEC3_RESIDUAL_ECHO_ESTIMATOR_H_
#define MODULES_AUDIO_PROCESSING_AEC3_RESIDUAL_ECHO_ESTIMATOR_H_

#include <array>
#include <memory>

#include "absl/types/optional.h"
#include "api/echo_canceller3_config.h"
#include "audio_processing/aec3/aec3_common.h"
#include "audio_processing/aec3/aec_state.h"
#include "audio_processing/aec3/render_buffer.h"
#include "audio_processing/aec3/reverb_model.h"
#include "audio_processing/aec3/spectrum_buffer.h"
#include "rtc_base/checks.h"

namespace webrtc {

class ResidualEchoEstimator {
 public:
  ResidualEchoEstimator(const EchoCanceller3Config& config,
                        size_t num_render_channels);
  ~ResidualEchoEstimator();

  ResidualEchoEstimator(const ResidualEchoEstimator&) = delete;
  ResidualEchoEstimator& operator=(const ResidualEchoEstimator&) = delete;

  void Estimate(
      const AecState& aec_state,
      const RenderBuffer& render_buffer,
      rtc::ArrayView<const std::array<float, kFftLengthBy2Plus1>> S2_linear,
      rtc::ArrayView<const std::array<float, kFftLengthBy2Plus1>> Y2,
      rtc::ArrayView<std::array<float, kFftLengthBy2Plus1>> R2);

 private:
  enum class ReverbType { kLinear, kNonLinear };

  // Resets the state.
  void Reset();

  // Updates estimate for the power of the stationary noise component in the
  // render signal.
  void UpdateRenderNoisePower(const RenderBuffer& render_buffer);

  // Adds the estimated unmodelled echo power to the residual echo power
  // estimate.
  void AddReverb(ReverbType reverb_type,
                 const AecState& aec_state,
                 const RenderBuffer& render_buffer,
                 rtc::ArrayView<std::array<float, kFftLengthBy2Plus1>> R2);

  const EchoCanceller3Config config_;
  const size_t num_render_channels_;
  std::array<float, kFftLengthBy2Plus1> X2_noise_floor_;
  std::array<int, kFftLengthBy2Plus1> X2_noise_floor_counter_;
  ReverbModel echo_reverb_;
};

}  // namespace webrtc

#endif  // MODULES_AUDIO_PROCESSING_AEC3_RESIDUAL_ECHO_ESTIMATOR_H_
