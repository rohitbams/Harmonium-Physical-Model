// Voice.h
#pragma once
#include <JuceHeader.h>
#include <vector>
#include "HarmoniumState.h"
#include "PhysicsEngine.h"

/*
 * The Voice class.
 * This class handles polyphonic voice stealing, read instances
 */
class Voice {
public:
    enum ReedMode {
        SINGLE_REED = 1,
        DOUBLE_REED = 2,
        TRIPLE_REED = 3
    };

private:
    struct ReedInstance {
        HarmoniumState physicsState;
        HarmoniumPhysicsConfig config;
        double detuning_cents;
        
        ReedInstance(double detune_cents = 0.0) : detuning_cents(detune_cents) {
            physicsState.reset();
        }
        
        float processSample(const HarmoniumPhysicsConfig& config, double bellowsFlow, double bellowsPressure, double baseFrequency) {
            // TODO: calculate detuned pitches
            double detunedFreq = PhysicsEngine::calculateDetunedFrequency(baseFrequency, detuning_cents);
            double omega0 = PhysicsEngine::calculateOmega0(detunedFreq);
//            double frequency = currentFrequency_;
//            double omega0 = PhysicsEngine::calculateOmega0(baseFrequency);
            
            physicsState = PhysicsEngine::main(physicsState,config,bellowsFlow,bellowsPressure,omega0);
            
            
            return PhysicsEngine::generateAudio(physicsState, config);
        }
    };

    HarmoniumPhysicsConfig physicsConfig_;
    std::vector<ReedInstance> reedInstances_;
    ReedMode reedMode_;
    bool isActive_ = false;
    double currentFrequency_ = 0;

public:
    explicit Voice(double sampleRate) : reedMode_(SINGLE_REED) {
        physicsConfig_.sampleRate = sampleRate;
        physicsConfig_.dt = 1.0 / sampleRate;
        reedInstances_.emplace_back(0.0);
        DBG("Voice initialized with Physics Engine");
    }
    
    void setReedMode(ReedMode mode) {
        reedMode_ = mode;
        
        if (isActive_) {
            createReedsForMode(currentFrequency_);
        }
        
        DBG("Reed mode: " + juce::String(mode == SINGLE_REED ? "Single" :
                                       mode == DOUBLE_REED ? "Double" : "Triple"));
    }
    
    void startNote(double frequency) {
        if (isActive_) return;
        
        currentFrequency_ = frequency;
        isActive_ = true;
        
        createReedsForMode(frequency);
        DBG("Voice startNote: " + juce::String(frequency) + " Hz, stored as: " + juce::String(currentFrequency_));
    }
    
    void stopNote() {
        isActive_ = false;
        
        // reset all reed
        for (auto& reed : reedInstances_) {
            reed.physicsState.reset();
        }
        
        DBG("Note stopped");
    }
    
    float processSample(double bellowsFlowRate, double bellowsPressure) {
        if (!isActive_ || reedInstances_.empty()) return 0.0f;
        
        float mixedOutput = 0.0f;
        double airflowPerReed = bellowsFlowRate / reedInstances_.size();
        
        // process each reed
        for (auto& reedInstance : reedInstances_) {
            float reedOutput = reedInstance.processSample(
                physicsConfig_,
                airflowPerReed,
                bellowsPressure,
                currentFrequency_
            );
            mixedOutput += reedOutput;
        }
        
        // Natural mixing for multiple reeds
//        if (reedInstances_.size() > 1) {
//            mixedOutput /= std::sqrt(reedInstances_.size());
//        }
//        DBG("Processing with frequency: " + juce::String(currentFrequency_));

        return mixedOutput;
    }
    
    bool isAvailable() const { return !isActive_; }
    bool getIsActive() const { return isActive_; }
    bool matchesFrequency(double freq) const {
        return isActive_ && std::abs(currentFrequency_ - freq) < 1.0;
    }
    double getCurrentFrequency() const { return currentFrequency_; }
    ReedMode getReedMode() const { return reedMode_; }
    
    double getChamberPressure() const {
        return isActive_ && !reedInstances_.empty() ? reedInstances_[0].physicsState.p1 : 0.0;
    }
    
    double getJetPressure() const {
        return isActive_ && !reedInstances_.empty() ? reedInstances_[0].physicsState.p2 : 0.0;
    }
    
    double getReedPosition() const {
        return isActive_ && !reedInstances_.empty() ? reedInstances_[0].physicsState.reedPosition : 0.0;
    }
    
    double getReedVelocity() const {
        return isActive_ && !reedInstances_.empty() ? reedInstances_[0].physicsState.reedVelocity : 0.0;
    }
    
    double getJetVelocity() const {
        return isActive_ && !reedInstances_.empty() ? reedInstances_[0].physicsState.vj : 0.0;
    }
    
    double getTotalVolumeFlow() const {
        if (!isActive_ || reedInstances_.empty()) return 0.0;
        double totalFlow = 0.0;
        for (const auto& reed : reedInstances_) {
            totalFlow += reed.physicsState.u;
        }
        return totalFlow;
    }
    
    double getDampingForce() const {
        return isActive_ && !reedInstances_.empty() ? reedInstances_[0].physicsState.dampingForce : 0.0;
    }
    
    void setAmplitudeScale(int amp) {
        physicsConfig_.amplitudeScale = amp;
    }

//    void setLows(int value) {
//        physicsConfig_.low = value;
//    }
    
    void setMu(double mu) {
        physicsConfig_.mu = std::clamp(mu, 1e-6, 1e-1);
    }
    
    bool isOscillating() const {
        for (const auto& reed : reedInstances_) {
            if (reed.physicsState.isOscillating()) {
                return true;
            }
        }
        return false;
    }
    
    int getOscillatingReedCount() const {
        int count = 0;
        for (const auto& reed : reedInstances_) {
            if (reed.physicsState.isOscillating()) {
                count++;
            }
        }
        return count;
    }
    
//    double getAverageAmplitude() const {
//        if (!isActive_ || reedInstances_.empty()) return 0.0;
//
//        double totalAmplitude = 0.0;
//        for (const auto& reed : reedInstances_) {
//            totalAmplitude += reed.physicsState.getAmplitude();
//        }
//        return totalAmplitude / reedInstances_.size();
//    }

private:
    void createReedsForMode(double frequency) {
        reedInstances_.clear();
        
        switch (reedMode_) {
            case SINGLE_REED:
                reedInstances_.emplace_back(0.0); // No detuning
                break;
                
            case DOUBLE_REED:
                reedInstances_.emplace_back(0.0); // Fundamental
                reedInstances_.emplace_back(-12.1); // One octave lower
                break;
                    
            case TRIPLE_REED:
                reedInstances_.emplace_back(-12.5); // Slightly flat octave
                reedInstances_.emplace_back(0.0); // Fundamental
                reedInstances_.emplace_back(+12.2); // Slightly sharp octave
                break;
        }
        
        DBG("Created " + juce::String(reedInstances_.size()) + " reed instances");
    }
};
