#pragma once
#include <JuceHeader.h>
#include <memory>
#include <array>
#include "Voice.h"
#include "Harmonium.h"
#include "Convolution.h"


/*
 * Harmonium class
 */
class Harmonium {
private:
    
    // voice management
    static const int NUM_VOICES = 6;
    std::array<std::unique_ptr<Voice>, NUM_VOICES> voices;
    std::unique_ptr<Bellows> bellows;
    std::unique_ptr<Convolution> convolution;
    
    // state
    HarmoniumState currentState;
//    bool isPolyphonic = true;

public:
    explicit Harmonium(double sampleRate) {
        
        // create voices
        for (int i = 0; i < NUM_VOICES; ++i) {
            voices[i] = std::make_unique<Voice>(sampleRate);
        }
        
        bellows = std::make_unique<Bellows>(sampleRate);
        convolution = std::make_unique<Convolution>();
        
    }
        
    // main per sample processing
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

        mixedOutput = convolution->processSample(mixedOutput);
        return mixedOutput;
    }
    
    void prepareToPlay(double sampleRate, int blockSize) {
        if (convolution) {
            convolution->prepare(sampleRate, blockSize);
        }
    }
    
    // ir loader
    bool loadImpulseResponse(const juce::File& irFile) {
        if (convolution) {
            bool success = convolution->loadIR(irFile);
            return success;
        }
        return false;
    }
    
    void setConvolutionEnabled(bool enabled) {
        if (convolution) {
            convolution->setEnabled(enabled);
        }
    }
    
    bool isConvolutionEnabled() const {
        return convolution ? convolution->getEnabled() : false;
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
    
    void setReedMode(Voice::ReedMode mode) {
        for (auto& voice : voices) {
            voice->setReedMode(mode);
        }
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
    // voice stealing
    Voice* selectVoice() {
        for (auto& voice : voices) {
            if (voice->isAvailable()) {
                return voice.get();
            }
        }
        return voices[0].get();
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
        
    void updateCurrentState() {
        Voice* activeVoice = nullptr;
        for (auto& voice : voices) {
            if (voice->getIsActive()) {
                activeVoice = voice.get();
                break;
            }
        }
        
        if (bellows) {
            currentState.p0 = bellows->getBellowsPressure();
            currentState.u0 = bellows->getFlowRate();
        }
        
        if (activeVoice) {
            currentState.p1 = activeVoice->getChamberPressure();
            currentState.p2 = activeVoice->getJetPressure();
            currentState.vj = activeVoice->getJetVelocity();
            currentState.reedPosition = activeVoice->getReedPosition();
            currentState.reedVelocity = activeVoice->getReedVelocity();
            currentState.u = activeVoice->getTotalVolumeFlow();
        }
        
        else {
            currentState.p1 = 0.0;
            currentState.p2 = 0.0;
            currentState.vj = 0.0;
            currentState.reedPosition = 0.0;
            currentState.reedVelocity = 0.0;
            currentState.u = 0.0;
        }
        
        currentState.oscillatingReeds = getTotalOscillatingReeds();
    }
};
