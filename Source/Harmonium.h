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
//    struct PhysicsState {
//        double p0 = 0.0; // bellows tank pressure
//        double p1 = 0.0; // reed chamber pressure
//        double p2 = 0.0; // jet pressure
//        double u0 = 0.0; // bellows flow rate
//        double u = 0.0; // total air volume inflow pumping
//        double vj = 0.0; // narrow jet velocity
//        double reedPosition = 0.0; // reed displacement
//        double reedVelocity = 0.0; // reed velocity
//        int oscillatingReeds = 0; // number of oscillating reeds
////        double averageAmplitude = 0.0; // average reed amplitude
//
//    };

private:

    /*
     * The Bellows class.
     * This private class handles the novel bellowing mechanism.
     * It calculates the air inflow and the pressure variation in the bellows chamber.
     * It also handels modulation wheen interaction on the MIDI keyboard.
     * The Bellows mode can be turned off to simulate a constant inifinite air supply.
     */
    class Bellows {
    private:
        double sampleRate;
        double dt;
        double airMass = 0.003; // initial air in bellows
        double maxAirMass = 0.02; //  maximum capacity
        double reedChamberMass = 0; // initial air in Reed chamber
        double maxReedChamberMass = 0; //  maximum capacity in reed chamber
        double narrowJetMass = 0; // initial air in narrow jet
        double maxNarrowJetMass = 0; // maximum NJ capacity
        double dampingForce = 0; // damping force
        double p0 = 0.0; // bellows pressure (Pa)
        double u0 = 0.0; // flow rate
        double modWheelValue = 0.0;
        bool keyPressed = false;
        double baseConsumptionRate = 0.01;
        
        // dynamic pumping state
        double pumpAmount = 0.05; // air mass added per pumping action
        double previousModWheelValue = 0.0;
        double minMovementThreshold = 0.01;  // minimum movement to trigger pump
        double continuousPumpingRate = 0.0; // current pumping rate (kgs/s)
        double pumpingRate = 0.8;
        
        // bellows chamber spring pressure variables
//        double springConstant_ = 500.0; // spring stiffness
//        double restVolume_ = 0.004; // bellows volume at rest (spring uncompressed)
//        double effectiveBellowsArea_ = 0.05; // effective area for bellows calculation
        
    public:
        Bellows(double sampleRate) : sampleRate(sampleRate), dt(1.0/sampleRate) {}
        

        /*
         * this method handles bellowing (pumping mechanism) and the air pressure p0 buildup
         * in the bellows chameber V0, and air flow from the bellows chamber into the
         * reed chamber V2.
         */
        void updateBellows(double chamberPressure, int oscillatingReeds) {
            
            // pumping mechanism and air mass management
            double modWheelMovement = modWheelValue - previousModWheelValue;
            double movementSpeed = std::abs(modWheelMovement) / dt;
            continuousPumpingRate = std::min(movementSpeed * 10.0, 20.0);
            
            // add air mass while pumping
            if (continuousPumpingRate > 0.001) {
                airMass += continuousPumpingRate * dt;
            }
            
            // reduce air mass while pressing note keys
            if (keyPressed && oscillatingReeds > 0) {
                double consumptionRate = baseConsumptionRate * oscillatingReeds * dt;
                airMass -= consumptionRate;
            }
            
            airMass = std::clamp(airMass, 0.0, maxAirMass); // limit air mass to avoid negative value division
            
            double airRatio = airMass / maxAirMass;
            p0 = airRatio * 300.0;

            double pressureDiff = p0 - chamberPressure;  // p0 - p1
            
            // calculate flow when there's significant pressure difference
            if (std::abs(pressureDiff) > 2.0) {
                
                double flowDirection = (pressureDiff > 0) ? 1.0 : -1.0;
//                double velocity = std::sqrt(std::abs(pressureDiff) / 1.225); // 1.255 is air density (rho0)
                double velocity = std::sqrt(2.0 * std::abs(pressureDiff) / 1.225); // 1.255 is air density (rho0)
                
                if (keyPressed) {
                    // keys pressed, full connection (valve fully open)
                    u0 = flowDirection * 0.001 * velocity;
                } else {
                    // keys not pressed, restricted connection (small leak)
                    u0 = flowDirection * 0.0003 * velocity;  // 30% of full flow
                }
                
                u0 = std::clamp(u0, -0.05, 0.05);
            } else {
                // pressures are equal, no flow
                u0 = 0.0;
            }

            
            previousModWheelValue = modWheelValue;
        }
        
        void setModWheelValue(double value) {
            modWheelValue = std::clamp(value, 0.0, 1.0);
        }
        
        void setKeyPressed(bool pressed) {
            keyPressed = pressed;
        }
        
        double getBellowsPressure() const { return p0; }
        double getFlowRate() const { return u0; }
        double getAirMass() const { return airMass; }
        double getRelativePressure() const { return std::min(1.0, p0 / 1000.0); }
        
        void setMaxAirMass(double mass) {
            maxAirMass = std::clamp(mass, 0.1, 0.4);
            if (airMass > maxAirMass) {
                airMass = maxAirMass;
            }
        }

        void setMaxReedChamberAirMass(double mass) {
            maxReedChamberMass = std::clamp(mass, 0.1, 0.2);
            if (reedChamberMass > maxReedChamberMass) reedChamberMass = maxReedChamberMass * pumpingRate;
        }
        
        void setMaxNarrowJetAirMass(double mass) {
            narrowJetMass = std::clamp(mass, 0.1, 0.2);
            if (narrowJetMass > maxNarrowJetMass) narrowJetMass = maxNarrowJetMass * pumpingRate;
        }
        
        void setBaseConsumptionRate(double rate) {
            baseConsumptionRate = std::clamp(rate, 0.001, 0.02);
        }
    };
    
    // voice management
    static const int NUM_VOICES = 6;
    std::array<std::unique_ptr<Voice>, NUM_VOICES> voices;
    std::unique_ptr<Bellows> bellows;
    std::unique_ptr<ConvolutionProcessor> convolutionProcessor;
    
    // state
    HarmoniumState currentState;
    bool isPolyphonic = true;

public:
    explicit Harmonium(double sampleRate) {
//        DBG("Harmonium initialising");
        
        // create voices
        for (int i = 0; i < NUM_VOICES; ++i) {
            voices[i] = std::make_unique<Voice>(sampleRate);
        }
        
        bellows = std::make_unique<Bellows>(sampleRate);
        convolutionProcessor = std::make_unique<ConvolutionProcessor>();
        
//        DBG("Harmonium initialised successfully");
    }
        
    // main processing per sample
    float processSample() {
        updateCurrentState();
        
        double chamberPressure = 0.0;
        int oscillatingReeds = getTotalOscillatingReeds();
        bellows->updateBellows(chamberPressure, oscillatingReeds);
        
        // process all voices
        float mixedOutput = 0.0f;
        double bellowsFlow = bellows->getFlowRate();
        double bellowsPressure = bellows->getBellowsPressure();
        
        for (auto& voice : voices) {
            float voiceOutput = voice->processSample(bellowsFlow, bellowsPressure);
            mixedOutput += voiceOutput;
        }

        mixedOutput = convolutionProcessor->processSample(mixedOutput);
        return mixedOutput;
    }
    
    void prepareToPlay(double sampleRate, int blockSize) {
        if (convolutionProcessor) {
            convolutionProcessor->prepare(sampleRate, blockSize);
        }
    }
    
    // ir loader
    bool loadImpulseResponse(const juce::File& irFile) {
        if (convolutionProcessor) {
            bool success = convolutionProcessor->loadIR(irFile);
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
        if (convolutionProcessor) {
            convolutionProcessor->setEnabled(enabled);
//            DBG("Convolution " + juce::String(enabled ? "enabled" : "disables"));
        }
    }
    
    bool isConvolutionEnabled() const {
        return convolutionProcessor ? convolutionProcessor->getEnabled() : false;
    }
    
    
    void startNote(double frequency) {
        
        Voice* selectedVoice = selectVoice();
        if (selectedVoice) {
            selectedVoice->startNote(frequency);
            updateBellowsKeyState();
            bellows->setKeyPressed(true);
//            DBG("Note started: " + juce::String(frequency) + " Hz");
        }
    }
    
    void stopNote(double frequency) {
        for (auto& voice : voices) {
            if (voice->matchesFrequency(frequency)) {
                voice->stopNote();
                break;
            }
        }
        updateBellowsKeyState();
//        DBG("Note stopped: " + juce::String(frequency) + " Hz");
    }
    
    void stopAllNotes() {
        for (auto& voice : voices) {
            voice->stopNote();
        }
        bellows->setKeyPressed(false);
//        DBG("All notes stopped");
    }
    
    void setModWheelValue(double value) {
        if (bellows) {
            bellows->setModWheelValue(value);
        }
    }
    
    double getModWheelValue() const {
        return bellows ? bellows->getRelativePressure() : 0.0;
    }
    
    void setAirCapacity(float value) {
        if (bellows) {
            double capacity = 0.001 + value * 0.01;
            bellows->setMaxAirMass(capacity);
        }
    }

    void setReedChamberCapacity(float value) {
        if (bellows) {
            double reedChamberCapacity = 0.001 + value * 0.01;
            bellows->setMaxReedChamberAirMass(reedChamberCapacity);
        }
    }
    
    void setNarrowJetCapacity(float value) {
        if (bellows) {
            double narrowJetCapacity = 0.001 + value * 0.01;
            bellows->setMaxNarrowJetAirMass(narrowJetCapacity);
        }
    }
    
    void setAirConsumption(float value) {
        if (bellows) {
            double consumption = 0.001 + value * 0.01;
            bellows->setBaseConsumptionRate(consumption);
        }
    }
    
    void setPolyphonyMode(bool polyphonic) {
        isPolyphonic = polyphonic;
        DBG("Polyphony: " + juce::String(polyphonic ? "ON" : "OFF"));
    }
    
    // reed mode
    void setReedMode(Voice::ReedMode mode) {
        for (auto& voice : voices) {
            voice->setReedMode(mode);
        }
//        DBG("Reed mode: " + juce::String(mode == Voice::SINGLE_REED ? "Single" :
//                                       mode == Voice::DOUBLE_REED ? "Double" : "Triple"));
    }
    
    Voice::ReedMode getReedMode() const {
        return voices[0]->getReedMode();
    }
    
    void setAmplitudeScale(int amp) {
        for (auto& voice : voices) {
            voice->setAmplitudeScale(amp);
        }
        DBG("Amplitude scale: " + juce::String(amp));
    }

//    void setLows(int value) {
//        for (auto& voice : voices) {
//            voice->setLows(value);
//        }
//        DBG("Lows: " + juce::String(value));
//    }


    // physics diagnostics
    const HarmoniumState& getPhysicsState() const {
        return currentState;
    }
    
    bool isNoteActive() const {
        for (const auto& voice : voices) {
            if (voice->getIsActive()) {
                return true;
            }
        }
        return false;
    }
    
    int getTotalOscillatingReeds() const {
        int total = 0;
        for (const auto& voice : voices) {
            total += voice->getOscillatingReedCount();
        }
        return total;
    }


private:
    std::unique_ptr<Harmonium> harmonium;
    
    Voice* selectVoice() {
        if (isPolyphonic) {
            // find available voice
            for (auto& voice : voices) {
                if (voice->isAvailable()) {
                    return voice.get();
                }
            }
            // voice stealing
            return voices[0].get();
        } else {
            // monophonic mode
            for (int i = 1; i < NUM_VOICES; ++i) {
                voices[i]->stopNote();
            }
            return voices[0].get();
        }
    }
    
    int getActiveVoiceCount() const {
        int count = 0;
        for (const auto& voice : voices) {
            if (voice->getIsActive()) {
                count++;
            }
        }
        return count;
    }

    
    void updateBellowsKeyState() {
        bool anyActive = false;
        for (auto& voice : voices) {
            if (voice->getIsActive()) {
                anyActive = true;
                break;
            }
        }
        
        if (!anyActive && bellows) {
            bellows->setKeyPressed(false);
        }
    }
        
    /*
     * this method updates the physcics state per processed sample
     */
    void updateCurrentState() {
        // find first active voice for diagnostics
        Voice* activeVoice = nullptr;
        for (auto& voice : voices) {
            if (voice->getIsActive()) {
                activeVoice = voice.get();
                break;
            }
        }
        
        // bellows state
        if (bellows) {
            currentState.p0 = bellows->getBellowsPressure();
            currentState.u0 = bellows->getFlowRate();
        }
        
        // voice state
        if (activeVoice) {
            currentState.p1 = activeVoice->getChamberPressure();
            currentState.p2 = activeVoice->getJetPressure();
            currentState.vj = activeVoice->getJetVelocity();
            currentState.reedPosition = activeVoice->getReedPosition();
            currentState.reedVelocity = activeVoice->getReedVelocity();
            currentState.u = activeVoice->getTotalVolumeFlow();
        } else {
            // no active voice - zero out
            currentState.p1 = 0.0;
            currentState.p2 = 0.0;
            currentState.vj = 0.0;
            currentState.reedPosition = 0.0;
            currentState.reedVelocity = 0.0;
            currentState.u = 0.0;
        }
        
        // global diagnostics
        currentState.oscillatingReeds = getTotalOscillatingReeds();
    }
};
