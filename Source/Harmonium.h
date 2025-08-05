// Harmonium.h
#pragma once
#include <JuceHeader.h>
#include <memory>
#include <array>
#include "Voice.h"
#include "Harmonium.h"
#include "ConvolutionProcessor.h"

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
     * The Bellows class.
     * This class handles the bellowing mechanism.
     * It calculates the air inflow and the pressure variation in the bellows chamber.
     * It also handels modulation wheen interaction on the MIDI keyboard.
     * The Bellows mode can be turned off to simulate a constant inifinite air supply.
     */
    class Bellows {
    private:
        double sampleRate_;
        double dt_;
        double airMass_ = 0.003;        // kg - initial air in bellows
        double maxAirMass_ = 0.02;     // kg - maximum capacity
        double reedChamberMass_ = 0;        // kg - initial air in Reed chamber
        double maxReedChamberMass_ = 0;     // kg - maximum capacity in reed chamber
        double narrowJetMass_ = 0;        // kg - initial air in NJ
        double maxNarrowJetMass_ = 0;     // kg - maximum NJ capacity
        double dampingForce = 0;     // damping force
        double p0_ = 0.0;               // Bellows pressure (Pa)
        double u0_ = 0.0;               // Flow rate (m³/s)
        double modWheelValue_ = 0.0;    // 0.0 to 1.0
        bool keyPressed_ = false;
        double baseConsumptionRate_ = 0.001;  // kg/s
        
        // dynamic pumping state
        double pumpAmount_ = 0.05;  // Air added per pump action
        double previousModWheelValue_ = 0.0;
        double minMovementThreshold_ = 0.01;  // Minimum movement to trigger pump
        double continuousPumpingRate_ = 0.0; // Current pumping rate (kgs/s)
        
        double springConstant_ = 5000000.0; //N/m - spring stiffness
        double restVolume_ = 0.004; // m³ bellows volume at rest (spring uncompressed)
        double effectiveBellowsArea_ = 0.05; // m² effective area for bellows calculation
        
    public:
        Bellows(double sampleRate) : sampleRate_(sampleRate), dt_(1.0/sampleRate) {}
        
        // simpleset version
        void updateBellows(double chamberPressure, int oscillatingReeds) {
            
            // Air mass management
            double modWheelMovement = modWheelValue_ - previousModWheelValue_;
            double movementSpeed = std::abs(modWheelMovement) / dt_;
            continuousPumpingRate_ = std::min(movementSpeed * 10.0, 20.0);
            
            if (continuousPumpingRate_ > 0.001) {
                airMass_ += continuousPumpingRate_ * dt_;
            }
            
            if (keyPressed_ && oscillatingReeds > 0) {
                double consumptionRate = baseConsumptionRate_ * oscillatingReeds * dt_;
                airMass_ -= consumptionRate;
            }
            
            airMass_ = std::clamp(airMass_, 0.0, maxAirMass_);
            
//            // SIMPLEST PRESSURE: direct proportion
//            p0_ = (airMass_ / maxAirMass_) * 2000.0;  // 0 to 4000 Pa
//            
//            // Flow
//            if (keyPressed_ && p0_ > chamberPressure) {
//                double pressureDiff = p0_ - chamberPressure;
//                u0_ = std::sqrt(pressureDiff) * 0.001;
//                u0_ = std::min(u0_, 0.1);
//            } else {
//                u0_ = 0.0;
//            }
            double airRatio = airMass_ / maxAirMass_;
            p0_ = airRatio * 4000.0;  // Your current calculation

            
            double pressureDiff = p0_ - chamberPressure;  // p0 - p1
            
            // Always calculate flow when there's pressure difference
            if (std::abs(pressureDiff) > 2.0) {  // Only if meaningful difference
                
                double flowDirection = (pressureDiff > 0) ? 1.0 : -1.0;
                double velocity = std::sqrt(std::abs(pressureDiff) / 1.225);
                
                if (keyPressed_) {
                    // Keys pressed: full connection (valve fully open)
                    u0_ = flowDirection * 0.001 * velocity;
                } else {
                    // Keys not pressed: restricted connection (small leak)
                    u0_ = flowDirection * 0.0003 * velocity;  // 30% of full flow
                }
                
                u0_ = std::clamp(u0_, -0.05, 0.05);  // Safety limits
                
            } else {
                // Pressures are equal: no flow (equilibrium)
                u0_ = 0.0;
            }

            
            previousModWheelValue_ = modWheelValue_;
        }


//        void updateBellows(double chamberPressure, int oscillatingReeds) {
//            // Mod wheel pumping - adds air
//                double modWheelMovement = modWheelValue_ - previousModWheelValue_;
//                double movementSpeed = std::abs(modWheelMovement) / dt_;
//                continuousPumpingRate_ = std::min(movementSpeed * 1, 12.0);
//                
//                if (continuousPumpingRate_ > 0.0001) {
//                    airMass_ += continuousPumpingRate_ * dt_;
//                }
//                // Air consumption - removes air
//                if (keyPressed_ && oscillatingReeds > 0) {
//                    double consumptionRate = baseConsumptionRate_ * oscillatingReeds * dt_;
//                    airMass_ -= consumptionRate;
//                }
//                // Keep air mass in bounds
//                airMass_ = std::clamp(airMass_, 0.0, maxAirMass_);
//    
//                // === PRESSURE CALCULATION ===
//                // Simple direct relationship: more air = more pressure
//                double airRatio = airMass_ / maxAirMass_;  // 0.0 to 1.0
//                double basePressure = airRatio * 3500.0;  // 0 to 3500 Pa
//                
//                // OPTIONAL: Add spring pressure (your enhancement)
//                double springPressure = 0.0;
//                if (airRatio > 0.5) {  // Only when bellows is significantly filled
//                    double compressionAmount = (airRatio - 0.5) * 2.0;  // 0.0 to 1.0
//                    springPressure = compressionAmount * compressionAmount * 500.0;  // 0 to 500 Pa
//                }
//                p0_ = basePressure + springPressure;
//                p0_ = std::clamp(p0_, 0.0, 4000.0);  // Your requested limit
//            
//                // === FLOW CALCULATION ===
//                if (keyPressed_ && airMass_ > 0.001 && p0_ > chamberPressure) {
//                    double pressureDiff = p0_ - chamberPressure;
//                    u0_ = std::sqrt(pressureDiff) * 0.001;  // Simple flow relationship
//                    u0_ = std::min(u0_, 0.1);
//                } else {
//                    u0_ = 0.0;
//                }
//        previousModWheelValue_ = modWheelValue_;
//    }
        
//        void updateBellows(double chamberPressure, int oscillatingReeds) {
//            // Pumping and consumption (unchanged)
//            double modWheelMovement = modWheelValue_ - previousModWheelValue_;
//            double movementSpeed = std::abs(modWheelMovement) / dt_;
//            continuousPumpingRate_ = std::min(movementSpeed * 2.5, 80.0);
//            
//            if (continuousPumpingRate_ > 0.0001) {
//                airMass_ += continuousPumpingRate_ * dt_;
//            }
//            
//            if (keyPressed_ && oscillatingReeds > 0) {
//                double consumptionRate = baseConsumptionRate_ * oscillatingReeds * dt_;
//                airMass_ -= consumptionRate;
//            }
//            
//            airMass_ = std::clamp(airMass_, 0.001, maxAirMass_);
//            
//            // === PHYSICS-BASED PRESSURE CALCULATION ===
//            // 1. Ideal gas pressure (dominant component)
//            const double AIR_DENSITY = 1.225;      // kg/m³
////            const double GAS_CONSTANT = 287.0;     // J/(kg·K)
////            const double TEMPERATURE = 293.0;      // K (20°C)
////            
////            double baseGasPressure = AIR_DENSITY * GAS_CONSTANT * TEMPERATURE;  // ≈103,011 Pa
////            double airMassRatio = airMass_ / maxAirMass_;  // 0.0 to 1.0
////            double gasPressure = baseGasPressure * airMassRatio;
//            
//            // 2. Spring pressure (your enhancement)
//            double springPressure = 0.0;
//            if (restVolume_ > 0.0 && effectiveBellowsArea_ > 1e-6) {
//                double currentVolume = airMass_ / AIR_DENSITY;
//                double volumeExpansion = std::max(0.0, currentVolume - restVolume_);
//                
//                if (volumeExpansion > 1e-6) {
//                    double springDisplacement = volumeExpansion / effectiveBellowsArea_;
//                    double springForce = springConstant_ * springDisplacement;
//                    springPressure = springForce / effectiveBellowsArea_;
//                    springPressure = std::clamp(springPressure, 0.0, 10000.0);  // Safety limit
//                }
//            }
//            
//            // 3. Combine both effects
//            p0_ = springPressure;
//            p0_ = std::clamp(p0_, 0.0, 200.0);  // Max ~2 atm
//            
//            // Flow calculation (unchanged)
//            if (keyPressed_ && airMass_ > 0.001 && p0_ > chamberPressure) {
//                double pressureDiff = std::max(0.0, p0_ - chamberPressure);
//                double velocity = std::sqrt(2.0 * pressureDiff / AIR_DENSITY);
//                u0_ = 0.001 * velocity;
//                u0_ = std::min(u0_, 0.1);
//            } else {
//                u0_ = 0.0;
//            }
//            
//            previousModWheelValue_ = modWheelValue_;
//        }

//        void updateBellows(double chamberPressure, int oscillatingReeds) {
//            // modwheel pumping calculation
//            double modWheelMovement = modWheelValue_ - previousModWheelValue_;
//            double movementSpeed = std::abs(modWheelMovement) / dt_; // capture movement per sample
//            continuousPumpingRate_ = std::min(movementSpeed * 2.0, 12.0);// cap to max rate
//            
//            if (continuousPumpingRate_ > 0.1) {
//                airMass_ += continuousPumpingRate_ * dt_; // increase airMass_ with pumping
//            }
//            
//             // air consumption
//            if (keyPressed_ && oscillatingReeds > 0) {
//                double consumptionRate = baseConsumptionRate_ * oscillatingReeds * dt_; // increase air consumption with more numbr of active reeds
//                airMass_ -= consumptionRate;
//            }
//            
//            airMass_ = std::clamp(airMass_, 0.001, maxAirMass_); // physical limit to never allow zero/negative value in airMass_
//            
//            // pressure calculation
//            double currentVolume = airMass_ / 1.225; // current air volume (ρ = 1.225 kg/m³)
//            double volumeExpansion = currentVolume - restVolume_;
//            volumeExpansion = std::max(0.0, volumeExpansion); // no negative compression
//            
//
//            // gas pressure (ideal gas law)
//            double airDensity = airMass_ / currentVolume;
//            double gasPressure = airDensity * 287.0 * 293.0;
//            
//            // Spring Pressure (resists expansion)
//            double springDisplacement = volumeExpansion / effectiveBellowsArea_; // Linear displacement
//            double springForce = springConstant_ * springDisplacement;
//            double springPressure = springForce / effectiveBellowsArea_;
//            
//            p0_ = gasPressure + springPressure;
//            p0_ = std::clamp(p0_, 0.0, gasPressure + springPressure);
//            
//            
//            if (keyPressed_ && airMass_ > 0.0 && p0_ > chamberPressure) {
//                double pressureDifference = p0_ - chamberPressure;
//                double velocity = std::sqrt(2.0 * pressureDifference / 1.225);
//                u0_ = 0.001 * velocity;
//            } else {
//                u0_ = 0.0;
//            }
//        }
        
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
            maxAirMass_ = std::clamp(mass, 0.1, 0.4);
            if (airMass_ > maxAirMass_) {
                airMass_ = maxAirMass_;
            }
        }

        void setMaxReedChamberAirMass(double mass) {
            maxReedChamberMass_ = std::clamp(mass, 0.1, 0.2);
            if (reedChamberMass_ > maxReedChamberMass_) reedChamberMass_ = maxReedChamberMass_ * 0.8;
        }
        
        void setMaxNarrowJetAirMass(double mass) {
            narrowJetMass_ = std::clamp(mass, 0.1, 0.2);
            if (narrowJetMass_ > maxNarrowJetMass_) narrowJetMass_ = maxNarrowJetMass_ * 0.8;
        }
        
        void setBaseConsumptionRate(double rate) {
            baseConsumptionRate_ = std::clamp(rate, 0.001, 0.02);
        }
    };
    
    // voice management
    static const int NUM_VOICES = 6;
    std::array<std::unique_ptr<Voice>, NUM_VOICES> voices_;
    std::unique_ptr<Bellows> bellows_;
    std::unique_ptr<ConvolutionProcessor> convolutionProcessor_;
    
    // state
    PhysicsState currentState_;
    bool isPolyphonic_ = true;

public:
    explicit Harmonium(double sampleRate) {
        DBG("Harmonium initializing with Physics Engine architecture...");
        
        // create voices
        for (int i = 0; i < NUM_VOICES; ++i) {
            voices_[i] = std::make_unique<Voice>(sampleRate);
        }
        
        bellows_ = std::make_unique<Bellows>(sampleRate);
        convolutionProcessor_ = std::make_unique<ConvolutionProcessor>();
        
        DBG("Harmonium initialized successfully");
    }
        
    // main processing
    float processSample() {
        updateCurrentState();
        
        double chamberPressure = 0.0;
        int oscillatingReeds = getTotalOscillatingReeds();
        bellows_->updateBellows(chamberPressure, oscillatingReeds);
        
        // process all voices
        float mixedOutput = 0.0f;
        double bellowsFlow = bellows_->getFlowRate();
        double bellowsPressure = bellows_->getBellowsPressure();
        
        for (auto& voice : voices_) {
            float voiceOutput = voice->processSample(bellowsFlow, bellowsPressure);
            mixedOutput += voiceOutput;
        }

        mixedOutput = convolutionProcessor_->processSample(mixedOutput);
        return mixedOutput;
    }
    
    void prepareToPlay(double sampleRate, int blockSize) {
        if (convolutionProcessor_) {
            convolutionProcessor_->prepare(sampleRate, blockSize);
        }
    }
    
    // ir loader
    bool loadImpulseResponse(const juce::File& irFile) {
        if (convolutionProcessor_) {
            bool success = convolutionProcessor_->loadIR(irFile);
            if(success) {
                DBG("Loaded IR: " + irFile.getFileName());
            } else {
                DBG("Failed to load IR: " + irFile.getFileName());
            }
            return success;
        }
        return false;
    }
    
    void setConvolutionEnabled(bool enabled) {
        if (convolutionProcessor_) {
            convolutionProcessor_->setEnabled(enabled);
            DBG("Convolution " + juce::String(enabled ? "enabled" : "disables"));
        }
    }
    
    bool isConvolutionEnabled() const {
        return convolutionProcessor_ ? convolutionProcessor_->getEnabled() : false;
    }
    
    // TODO: why is the new note being played when an old note is held???
    void startNote(double frequency) {
        
        Voice* selectedVoice = selectVoice();
        if (selectedVoice) {
            selectedVoice->startNote(frequency);
            updateBellowsKeyState();
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

    void setReedChamberCapacity(float value) {
        if (bellows_) {
            double reedChamberCapacity = 0.001 + value * 0.01;
            bellows_->setMaxReedChamberAirMass(reedChamberCapacity);
        }
    }
    
    void setNarrowJetCapacity(float value) {
        if (bellows_) {
            double narrowJetCapacity = 0.001 + value * 0.01;
            bellows_->setMaxNarrowJetAirMass(narrowJetCapacity);
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
    
    // reedMode config
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


private:
    std::unique_ptr<Harmonium> harmonium;
    
    Voice* selectVoice() {
        if (isPolyphonic_) {
            // find available voice
            for (auto& voice : voices_) {
                if (voice->isAvailable()) {
                    return voice.get();
                }
            }
            // voice stealing
            return voices_[0].get();
        } else {
            // monophonic mode
            for (int i = 1; i < NUM_VOICES; ++i) {
                voices_[i]->stopNote();
            }
            return voices_[0].get();
        }
    }
    
    int getActiveVoiceCount() const {
        int count = 0;
        for (const auto& voice : voices_) {
            if (voice->getIsActive()) {
                count++;
            }
        }
        return count;
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
        // find first active voice for diagnostics
        Voice* activeVoice = nullptr;
        for (auto& voice : voices_) {
            if (voice->getIsActive()) {
                activeVoice = voice.get();
                break;
            }
        }
        
        // bellows state
        if (bellows_) {
            currentState_.p0 = bellows_->getBellowsPressure();
            currentState_.u0 = bellows_->getFlowRate();
        }
        
        // voice state
        if (activeVoice) {
            currentState_.p1 = activeVoice->getChamberPressure();
            currentState_.p2 = activeVoice->getJetPressure();
            currentState_.vj = activeVoice->getJetVelocity();
            currentState_.reedPosition = activeVoice->getReedPosition();
            currentState_.reedVelocity = activeVoice->getReedVelocity();
            currentState_.u = activeVoice->getTotalVolumeFlow();
        } else {
            // no active voice - zero out
            currentState_.p1 = 0.0;
            currentState_.p2 = 0.0;
            currentState_.vj = 0.0;
            currentState_.reedPosition = 0.0;
            currentState_.reedVelocity = 0.0;
            currentState_.u = 0.0;
        }
        
        // global diagnostics
        currentState_.oscillatingReeds = getTotalOscillatingReeds();
    }
};



/*
 * my updateBellows() copy for reference
 *
 */

// // TODO: add realistic pumping mechanism
// // TODO: noteStart should depend on the p1 (pressure in reed chamber, not bellow chamber)
// void updateBellows(double chamberPressure, int oscillatingReeds) {
//     
//     double modWheelMovement = modWheelValue_ - previousModWheelValue_;
//     // continuous pumping
//     double movementSpeed = std::abs(modWheelMovement) / dt_; // movement per sample
////            continuousPumpingRate_ = std::min(movementSpeed * 0.1, 1.1); // cap to max rate
//     continuousPumpingRate_ = std::min(movementSpeed * 2.5, 80.0);
//     
//     // apply continuous pumping
////            if (continuousPumpingRate_ > 0.0001) {
//////                         if (continuousPumpingRate_ > baseConsumptionRate_) {
////                airMass_ += continuousPumpingRate_ * dt_; // air mass increasing per sample
////                airMass_ = std::max(0.1, airMass_);
////            }
//     
//     if (continuousPumpingRate_ > 0.0001) {
//         airMass_ += continuousPumpingRate_ * dt_; // air mass increasing per sample
//     }
//     
////            previousModWheelValue_ = modWheelValue_;
//     
//      // 2. SAFE AIR CONSUMPTION
//     if (keyPressed_ && oscillatingReeds > 0) {
//         double consumptionRate = baseConsumptionRate_ * oscillatingReeds * dt_;
//         airMass_ -= consumptionRate * 2;
//     }
//     
//     // 3. ENFORCE PHYSICAL LIMITS (CRITICAL!)
//     airMass_ = std::clamp(airMass_, 0.001, maxAirMass_); // Never allow zero/negative
//     
//     // consumption
//     if (keyPressed_ && airMass_ > 0.0 && oscillatingReeds > 0) {
//         double consumptionRate = baseConsumptionRate_ * oscillatingReeds * dt_;
//         airMass_ -= consumptionRate;
//         airMass_ = std::max(0.1, airMass_);
//     }
//     
//     // pressure calculation
//     double currentVolume = airMass_ / 1.225; // current air volume (ρ = 1.225 kg/m³)
//     double volumeExpansion = currentVolume - restVolume_;
//     volumeExpansion = std::max(0.0, volumeExpansion); // no negative compression
//     
////            double airDensity = airMass_ / 0.01;
////            p0_ = airDensity * 287.0 * 293.0;
//     
//     // gas pressure (ideal gas law)
//     double airDensity = airMass_ / currentVolume;
//     double gasPressure = airDensity * 287.0 * 293.0;
//     
//     // Spring Pressure (resists expansion)
//     double springDisplacement = volumeExpansion / effectiveBellowsArea_; // Linear displacement
//     double springForce = springConstant_ * springDisplacement;
//     double springPressure = springForce / effectiveBellowsArea_;
//     
//     p0_ = gasPressure + springPressure /** 10*/;
//     p0_ = std::clamp(p0_, 0.0, gasPressure + springPressure /** 10*/);
////                DBG("p0_ value: " + juce::String(p0_));
//     
//     
//     if (keyPressed_ && airMass_ > 0.0 && p0_ > chamberPressure) {
//         double pressureDifference = p0_ - chamberPressure;
//         double velocity = std::sqrt(2.0 * pressureDifference / 1.225);
//         u0_ = 0.001 * velocity;
////                DBG("Bellows flow: u0=" + juce::String(u0_, 6) + ", p0=" + juce::String(p0_, 1) + ", keyPressed=" + juce::String(keyPressed_));
//     } else {
//         u0_ = 0.0;
////                DBG("No bellows flow - keyPressed=" + juce::String(keyPressed_) + ", airMass=" + juce::String(airMass_, 4));
//
//     }
// }
//


//        new
//        void updateBellows(double chamberPressure, int oscillatingReeds) {
//            // 1. SAFE MODWHEEL PUMPING
//            double modWheelMovement = modWheelValue_ - previousModWheelValue_;
//            double movementSpeed = std::abs(modWheelMovement) / dt_;
//            continuousPumpingRate_ = std::min(movementSpeed * 2.5, 80.0); // Cap pumping rate
//
//            if (continuousPumpingRate_ > 0.0001) {
//                airMass_ += continuousPumpingRate_ * dt_;
//            }
//
//            // 2. SAFE AIR CONSUMPTION
//            if (keyPressed_ && oscillatingReeds > 0) {
//                double consumptionRate = baseConsumptionRate_ * oscillatingReeds * dt_;
//                airMass_ -= consumptionRate * 2;
//            }
//
//            // 3. ENFORCE PHYSICAL LIMITS (CRITICAL!)
//            airMass_ = std::clamp(airMass_, 0.001, maxAirMass_); // Never allow zero/negative
//
//            // 5. SIMPLIFIED PRESSURE MODEL (Avoid complex spring physics for now)
//            // Use direct pressure relationship instead of complex spring model
//            double basePressure = (airMass_ / maxAirMass_) * 3000.0; // Scale to reasonable pressure
//
//            // 6. PRESSURE ENHANCEMENT FROM PUMPING
//            double pumpingBoost = continuousPumpingRate_ * 500.0; // Immediate pressure from pumping
//
//            p0_ = basePressure + pumpingBoost;
//            p0_ = std::clamp(p0_, 0.0, 8000.0); // Physical limits
//
//            // 7. SAFE FLOW CALCULATION
//            if (keyPressed_ && airMass_ > 0.001 && p0_ > chamberPressure) {
//                double pressureDifference = std::max(0.0, p0_ - chamberPressure);
//                double velocity = std::sqrt(2.0 * pressureDifference / 1.225);
//                u0_ = std::min(0.001 * velocity, 0.1); // Cap flow rate
//            } else {
//                u0_ = 0.0;
//            }
//
//            // 8. UPDATE HISTORY
//            previousModWheelValue_ = modWheelValue_;
//
//            // 9. DEBUG OUTPUT (temporary)
//            if (std::isnan(p0_) || std::isinf(p0_)) {
//                DBG("ERROR: p0 is NaN/Inf! airMass=" + juce::String(airMass_) +
//                    ", basePressure=" + juce::String(basePressure));
//                p0_ = 100.0; // Emergency fallback
//            }
//        }
