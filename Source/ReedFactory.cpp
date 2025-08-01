#include "ReedModel.h"
#include "LumpedReed.h"
#include "ClampedReed.h"
#include <memory>

std::unique_ptr<ReedModel> createReedModel(ReedModelType type, double targetFreq, double sampleRate) {
    switch (type) {
        case ReedModelType::LUMPED_PARAMETER:
            return std::make_unique<LumpedReed>(targetFreq, sampleRate);
            
        case ReedModelType::CLAMPED_BEAM:
            return std::make_unique<ClampedReed>(targetFreq, sampleRate);
            
        default:
            return std::make_unique<LumpedReed>(targetFreq, sampleRate);
    }
}