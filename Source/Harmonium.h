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
        double p0 = 0.0; // bellows tank pressure
        double p1 = 0.0; // reed chamber pressure
        double p2 = 0.0; // jet pressure
        double u0 = 0.0; // bellows flow rate
        double u = 0.0; // total air volume inflow pumping
        double vj = 0.0; // narrow jet velocity
        double reedPosition = 0.0; // reed displacement
        double reedVelocity = 0.0; // reed velocity
        int oscillatingReeds = 0; // number of oscillating reeds
//        double averageAmplitude = 0.0; // average reed amplitude

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
        double airMass_ = 0.003; // initial air in bellows
        double maxAirMass_ = 0.02; //  maximum capacity
        double reedChamberMass_ = 0; // initial air in Reed chamber
        double maxReedChamberMass_ = 0; //  maximum capacity in reed chamber
        double narrowJetMass_ = 0; // initial air in narrow jet
        double maxNarrowJetMass_ = 0; // maximum NJ capacity
        double dampingForce = 0; // damping force
        double p0_ = 0.0; // bellows pressure (Pa)
        double u0_ = 0.0; // flow rate
        double modWheelValue_ = 0.0;
        bool keyPressed_ = false;
        double baseConsumptionRate_ = 0.01;
        
        // dynamic pumping state
        double pumpAmount_ = 0.05; // air mass added per pumping action
        double previousModWheelValue_ = 0.0;
        double minMovementThreshold_ = 0.01;  // minimum movement to trigger pump
        double continuousPumpingRate_ = 0.0; // current pumping rate (kgs/s)
        
        // bellows chamber spring pressure variables
//        double springConstant_ = 500.0; // spring stiffness
//        double restVolume_ = 0.004; // bellows volume at rest (spring uncompressed)
//        double effectiveBellowsArea_ = 0.05; // effective area for bellows calculation
        
    public:
        Bellows(double sampleRate) : sampleRate_(sampleRate), dt_(1.0/sampleRate) {}
        

        /*
         * this method handles bellowing (pumping mechanism) and the air pressure p0 buildup
         * in the bellows chameber V0, and air flow from the bellows chamber into the
         * reed chamber V2.
         */
        void updateBellows(double chamberPressure, int oscillatingReeds) {
            
            // pumping mechanism and air mass management
            double modWheelMovement = modWheelValue_ - previousModWheelValue_;
            double movementSpeed = std::abs(modWheelMovement) / dt_;
            continuousPumpingRate_ = std::min(movementSpeed * 10.0, 20.0);
            
            // add air mass while pumping
            if (continuousPumpingRate_ > 0.001) {
                airMass_ += continuousPumpingRate_ * dt_;
            }
            
            // reduce air mass while pressing note keys
            if (keyPressed_ && oscillatingReeds > 0) {
                double consumptionRate = baseConsumptionRate_ * oscillatingReeds * dt_;
                airMass_ -= consumptionRate;
            }
            
            airMass_ = std::clamp(airMass_, 0.0, maxAirMass_); // limit air mass to avoid negative value division
            
            double airRatio = airMass_ / maxAirMass_;
            p0_ = airRatio * 300.0;

            double pressureDiff = p0_ - chamberPressure;  // p0 - p1
            
            // calculate flow when there's significant pressure difference
            if (std::abs(pressureDiff) > 2.0) {
                
                double flowDirection = (pressureDiff > 0) ? 1.0 : -1.0;
                double velocity = std::sqrt(std::abs(pressureDiff) / 1.225);
                
                if (keyPressed_) {
                    // keys pressed, full connection (valve fully open)
                    u0_ = flowDirection * 0.001 * velocity;
                } else {
                    // keys not pressed, restricted connection (small leak)
                    u0_ = flowDirection * 0.0003 * velocity;  // 30% of full flow
                }
                
                u0_ = std::clamp(u0_, -0.05, 0.05);
            } else {
                // pressures are equal, no flow
                u0_ = 0.0;
            }

            
            previousModWheelValue_ = modWheelValue_;
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
//        DBG("Harmonium initialising");
        
        // create voices
        for (int i = 0; i < NUM_VOICES; ++i) {
            voices_[i] = std::make_unique<Voice>(sampleRate);
        }
        
        bellows_ = std::make_unique<Bellows>(sampleRate);
        convolutionProcessor_ = std::make_unique<ConvolutionProcessor>();
        
//        DBG("Harmonium initialised successfully");
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
//            if(success) {
//                DBG("Loaded IR: " + irFile.getFileName());
//            } else {
//                DBG("Failed to load IR: " + irFile.getFileName());
//            }
            return success;
        }
        return false;
    }
    
    void setConvolutionEnabled(bool enabled) {
        if (convolutionProcessor_) {
            convolutionProcessor_->setEnabled(enabled);
//            DBG("Convolution " + juce::String(enabled ? "enabled" : "disables"));
        }
    }
    
    bool isConvolutionEnabled() const {
        return convolutionProcessor_ ? convolutionProcessor_->getEnabled() : false;
    }
    
    // TODO: fix new note being played when an old note is held
    void startNote(double frequency) {
        
        Voice* selectedVoice = selectVoice();
        if (selectedVoice) {
            selectedVoice->startNote(frequency);
            updateBellowsKeyState();
            bellows_->setKeyPressed(true);
//            DBG("Note started: " + juce::String(frequency) + " Hz");
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
//        DBG("Note stopped: " + juce::String(frequency) + " Hz");
    }
    
    void stopAllNotes() {
        for (auto& voice : voices_) {
            voice->stopNote();
        }
        bellows_->setKeyPressed(false);
//        DBG("All notes stopped");
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
    
    // reed mode
    void setReedMode(Voice::ReedMode mode) {
        for (auto& voice : voices_) {
            voice->setReedMode(mode);
        }
//        DBG("Reed mode: " + juce::String(mode == Voice::SINGLE_REED ? "Single" :
//                                       mode == Voice::DOUBLE_REED ? "Double" : "Triple"));
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

    void setLows(int value) {
        for (auto& voice : voices_) {
            voice->setLows(value);
        }
        DBG("Lows: " + juce::String(value));
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
