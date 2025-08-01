// Harmonium.h - Complete implementation with simple bellows
#pragma once
#include <JuceHeader.h>
#include <memory>
#include <array>
#include "Voice.h"
#include "Harmonium.h"

class Harmonium {
public:
    struct PhysicsState {
        double p0 = 0.0;                // Bellows tank pressure
        double p1 = 0.0;                // Reed chamber pressure
        double p2 = 0.0;                // Jet pressure
        double u0 = 0.0;                // Bellows flow rate
        double u = 0.0;                 // Total air volume inflow
        double vj = 0.0;                // Jet velocity
        double reedPosition = 0.0;      // Reed displacement
        double reedVelocity = 0.0;      // Reed velocity
        int oscillatingReeds = 0;       // Number of oscillating reeds
        double averageAmplitude = 0.0;  // Average reed amplitude
        

    };

private:
    /*
     *
     */
    class Bellows {
    private:
        double sampleRate_;
        double dt_;
        double airMass_ = 0.003;        // kg - initial air in bellows
        double maxAirMass_ = 0.005;     // kg - maximum capacity
        double p0_ = 0.0;               // Bellows pressure (Pa)
        double u0_ = 0.0;               // Flow rate (m³/s)
        double modWheelValue_ = 0.0;    // 0.0 to 1.0
        bool keyPressed_ = false;
        double baseConsumptionRate_ = 0.003;  // kg/s
        double maxPumpingRate_ = 0.01;
        
        // dynamic pumping state
        double basePumpingRate = 0.005;
        double discretePumpAmount = 0.002;
        double pumpingSensitivity = 0.05;
        double currentPumpingRate = 0.0;
        double dynamicPressure = 0.0;
        double basePressure = 0.0;
        
    public:
        Bellows(double sampleRate) : sampleRate_(sampleRate), dt_(1.0/sampleRate) {}
        
        // TODO: add realistic pumping mechanism
        // TODO: noteStart should depend on the p1 (pressure in reed chamber, not bellow chamber)
        void updateBellows(double chamberPressure) {
            // Bellow pumping mechanism mod wheel adds air
            double pumpingIntensity = 1.0 - modWheelValue_;  // Invert the relationship
            if (pumpingIntensity > 0.1) {  // Only pump when mod wheel is down
                double pumpingRate = pumpingIntensity * maxPumpingRate_ * dt_;

                airMass_ += pumpingRate;
                airMass_ = std::min(airMass_, maxAirMass_);
            }
            
            DBG("modWHeelValue: " + juce::String(modWheelValue_));
            
            
            // Air consumption when key pressed
            if (keyPressed_ && airMass_ > 0.0) {
                double consumptionRate = baseConsumptionRate_ * dt_;
                airMass_ -= consumptionRate;
                airMass_ = std::max(0.0, airMass_);
            }
            
            // Calculate pressure from air mass (simple ideal gas)
            double airDensity = airMass_ / 0.01;  // Assume 0.01 m³ bellows volume
            p0_ = airDensity * 287.0 * 293.0;     // P = ρRT
            
            // Calculate flow rate (simple Bernoulli)
            if (keyPressed_ && airMass_ > 0.0 && p0_ > chamberPressure) {
                double pressureDiff = p0_ - chamberPressure;
                double velocity = std::sqrt(2.0 * pressureDiff / 1.225);  // v = √(2ΔP/ρ)
                u0_ = 0.001 * velocity;  // Q = A·v (assume 0.001 m² outlet area)
            } else {
                u0_ = 0.0;
            }
        }
        
        void setModWheelValue(double value) {
            modWheelValue_ = std::clamp(value, 0.0, 1.0);
        }
        
        void setKeyPressed(bool pressed) {
            keyPressed_ = pressed;
        }
        
        double getBellowsPressure() const { return p0_; }
        double getFlowRate() const { return u0_; }
        double getAirMass() const { return airMass_; }
        double getRelativePressure() const { return std::min(1.0, p0_ / 1000.0); }
        
        void setMaxAirMass(double mass) {
            maxAirMass_ = std::clamp(mass, 0.001, 0.02);
            if (airMass_ > maxAirMass_) airMass_ = maxAirMass_ * 0.8;
        }
        
        void setBaseConsumptionRate(double rate) {
            baseConsumptionRate_ = std::clamp(rate, 0.001, 0.02);
        }
    };
    
    // Voice management
    static const int NUM_VOICES = 4;
    std::array<std::unique_ptr<Voice>, NUM_VOICES> voices_;
    std::unique_ptr<Bellows> bellows_;
    
    // State
    PhysicsState currentState_;
    bool isPolyphonic_ = true;

public:
    explicit Harmonium(double sampleRate) {
        DBG("Harmonium initializing with Physics Engine architecture...");
        
        // Create voices
        for (int i = 0; i < NUM_VOICES; ++i) {
            voices_[i] = std::make_unique<Voice>(sampleRate);
        }
        
        // Create simple bellows
        bellows_ = std::make_unique<Bellows>(sampleRate);
        
        DBG("Harmonium initialized successfully");
    }
    
    // Main processing
    float processSample() {
        updateCurrentState();
        
        // Calculate average chamber pressure from active voices
        double avgChamberPressure = 0.0;
        int activeCount = 0;
        for (auto& voice : voices_) {
            if (voice->getIsActive()) {
                avgChamberPressure += voice->getChamberPressure();
                activeCount++;
            }
        }
        if (activeCount > 0) {
            avgChamberPressure /= activeCount;
        }
        
        // Update bellows
        bellows_->updateBellows(avgChamberPressure);
        
        // Process all voices
        float mixedOutput = 0.0f;
        double bellowsFlow = bellows_->getFlowRate();
        double bellowsPressure = bellows_->getBellowsPressure();
        
        for (auto& voice : voices_) {
            float voiceOutput = voice->processSample(bellowsFlow, bellowsPressure);
            mixedOutput += voiceOutput;
        }
        
        // Scale output for multiple voices
        if (activeCount > 0) {
            mixedOutput /= std::sqrt(activeCount);
        }
        
        return mixedOutput;
    }
    
    // Musical interface
    // TODO: why is the new note being played when an old note is held???
    void startNote(double frequency) {
        // if (bellows_->getAirMass() < 0.001) {
        if (currentState_.p1 > 0.001) {
            DBG("No air in reservoir - note blocked");
            return;
        }
        
        Voice* selectedVoice = selectVoice();
        if (selectedVoice) {
            selectedVoice->startNote(frequency);
            bellows_->setKeyPressed(true);
            DBG("Note started: " + juce::String(frequency) + " Hz");
        }
    }
    
    void stopNote(double frequency) {
        for (auto& voice : voices_) {
            if (voice->matchesFrequency(frequency)) {
                voice->stopNote();
                break;
            }
        }
        updateBellowsKeyState();
        DBG("Note stopped: " + juce::String(frequency) + " Hz");
    }
    
    void stopAllNotes() {
        for (auto& voice : voices_) {
            voice->stopNote();
        }
        bellows_->setKeyPressed(false);
        DBG("All notes stopped");
    }
    
    // Control interface
    void setModWheelValue(double value) {
        if (bellows_) {
            bellows_->setModWheelValue(value);
        }
    }
    
    double getModWheelValue() const {
        return bellows_ ? bellows_->getRelativePressure() : 0.0;
    }
    
    void setAirCapacity(float value) {
        if (bellows_) {
            double capacity = 0.001 + value * 0.01;
            bellows_->setMaxAirMass(capacity);
        }
    }
    
    void setAirConsumption(float value) {
        if (bellows_) {
            double consumption = 0.001 + value * 0.01;
            bellows_->setBaseConsumptionRate(consumption);
        }
    }
    
    void setPolyphonyMode(bool polyphonic) {
        isPolyphonic_ = polyphonic;
        DBG("Polyphony: " + juce::String(polyphonic ? "ON" : "OFF"));
    }
    
    // Reed configuration
    void setReedMode(Voice::ReedMode mode) {
        for (auto& voice : voices_) {
            voice->setReedMode(mode);
        }
        DBG("Reed mode: " + juce::String(mode == Voice::SINGLE_REED ? "Single" :
                                       mode == Voice::DOUBLE_REED ? "Double" : "Triple"));
    }
    
    Voice::ReedMode getReedMode() const {
        return voices_[0]->getReedMode();
    }
    
    void setAmplitudeScale(int amp) {
        for (auto& voice : voices_) {
            voice->setAmplitudeScale(amp);
        }
        DBG("Amplitude scale: " + juce::String(amp));
    }
    
    // Physics diagnostics
    const PhysicsState& getPhysicsState() const {
        return currentState_;
    }
    
    bool isNoteActive() const {
        for (const auto& voice : voices_) {
            if (voice->getIsActive()) {
                return true;
            }
        }
        return false;
    }
    
    int getTotalOscillatingReeds() const {
        int total = 0;
        for (const auto& voice : voices_) {
            total += voice->getOscillatingReedCount();
        }
        return total;
    }
    
    double getAverageReedAmplitude() const {
        double totalAmplitude = 0.0;
        int activeVoices = 0;
        
        for (const auto& voice : voices_) {
            if (voice->getIsActive()) {
                totalAmplitude += voice->getAverageAmplitude();
                activeVoices++;
            }
        }
        
        return activeVoices > 0 ? totalAmplitude / activeVoices : 0.0;
    }

private:
    std::unique_ptr<Harmonium> harmonium;
    
    Voice* selectVoice() {
        if (isPolyphonic_) {
            // Find available voice
            for (auto& voice : voices_) {
                if (voice->isAvailable()) {
                    return voice.get();
                }
            }
            // Voice stealing
            return voices_[0].get();
        } else {
            // Monophonic mode
            for (int i = 1; i < NUM_VOICES; ++i) {
                voices_[i]->stopNote();
            }
            return voices_[0].get();
        }
    }
    
    void updateBellowsKeyState() {
        bool anyActive = false;
        for (auto& voice : voices_) {
            if (voice->getIsActive()) {
                anyActive = true;
                break;
            }
        }
        
        if (!anyActive && bellows_) {
            bellows_->setKeyPressed(false);
        }
    }
    
    void updateCurrentState() {
        // Find first active voice for diagnostics
        Voice* activeVoice = nullptr;
        for (auto& voice : voices_) {
            if (voice->getIsActive()) {
                activeVoice = voice.get();
                break;
            }
        }
        
        // Bellows state
        if (bellows_) {
            currentState_.p0 = bellows_->getBellowsPressure();
            currentState_.u0 = bellows_->getFlowRate();
        }
        
        // Voice state
        if (activeVoice) {
            currentState_.p1 = activeVoice->getChamberPressure();
            currentState_.p2 = activeVoice->getJetPressure();
            currentState_.vj = activeVoice->getJetVelocity();
            currentState_.reedPosition = activeVoice->getReedPosition();
            currentState_.reedVelocity = activeVoice->getReedVelocity();
            currentState_.u = activeVoice->getTotalVolumeFlow();
        } else {
            // No active voice - zero out
            currentState_.p1 = 0.0;
            currentState_.p2 = 0.0;
            currentState_.vj = 0.0;
            currentState_.reedPosition = 0.0;
            currentState_.reedVelocity = 0.0;
            currentState_.u = 0.0;
        }
        
        // Global diagnostics
        currentState_.oscillatingReeds = getTotalOscillatingReeds();
        currentState_.averageAmplitude = getAverageReedAmplitude();
    }
};

