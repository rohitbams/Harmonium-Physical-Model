#pragma once
#include <memory>

enum class ReedModelType {
    LUMPED_PARAMETER,
    CLAMPED_BEAM
};

class ReedModel {
public:
    struct PhysicsState {
        double position = 0.0;      // Reed tip displacement (m)
        double velocity = 0.0;      // Reed tip velocity (m/s)
        double frequency = 0.0;     // Fundamental frequency (Hz)
        double amplitude = 0.0;     // Current amplitude (m)
    };

    virtual ~ReedModel() = default;
    
    // Core interface that all reed models must implement
    virtual void updateMotion(double pressureDifference, double dt) = 0;
    virtual void reset() = 0;
    virtual PhysicsState getState() const = 0;
    virtual void setFrequency(double frequency) = 0;
    virtual double getTipDisplacement() const = 0;
    virtual double getTipVelocity() const = 0;
    
    // Optional interface for advanced models
    virtual bool supportsModalAnalysis() const { return false; }
    virtual int getNumModes() const { return 1; }
    virtual double getModalAmplitude(int mode) const { return 0.0; }
};

// Factory function declaration
std::unique_ptr<ReedModel> createReedModel(ReedModelType type, double targetFreq, double sampleRate);
