// PluginProcessor.h

#pragma once
#include <JuceHeader.h>
#include "Harmonium.h"

class HarmoniumPhysicsEngineAudioProcessor : public juce::AudioProcessor {
public:
    HarmoniumPhysicsEngineAudioProcessor();
    ~HarmoniumPhysicsEngineAudioProcessor() override;

    void prepareToPlay(double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;
    void processBlock(juce::AudioBuffer<float>&, juce::MidiBuffer&) override;

    juce::AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override;

    const juce::String getName() const override;
    bool acceptsMidi() const override;
    bool producesMidi() const override;
    bool isMidiEffect() const override;
    double getTailLengthSeconds() const override;

    int getNumPrograms() override;
    int getCurrentProgram() override;
    void setCurrentProgram(int index) override;
    const juce::String getProgramName(int index) override;
    void changeProgramName(int index, const juce::String& newName) override;

    void getStateInformation(juce::MemoryBlock& destData) override;
    void setStateInformation(const void* data, int sizeInBytes) override;

    
    void startNote(double frequency);
    void stopNote(double frequency);
    void stopAllNotes();
    bool isNoteCurrentlyOn() const;
    
    void setLows(int value);
    
    void setAmplitudeScaling(int amp);
    void setAirCapacity(float value);
    void setAirConsumption(float value);
    void setMasterGain(float value);
    void setReedChamberCapacity(float value);
    void setNarrowJetCapacity(float value);
    void setPolyphonyMode(bool polyMode);
    void setReedMode(int mode);
    int getReedMode() const;

    double getBellowsPressure() const;
    double getChamberPressure() const;
    double getJetPressure() const;
    double getReedPosition() const;
    double getReedVelocity() const;
    double getJetVelocity() const;
    int getTotalOscillatingReeds() const;
    double getAverageReedAmplitude() const;

    void handleMidiCC(int ccNumber, int ccValue);
    void setModWheelValue(double value);

    bool loadImpulseResponse(const juce::File& irFile) {
        return harmonium ? harmonium->loadImpulseResponse(irFile) : false;
    }

    void setConvolutionEnabled(bool enabled) {
        if (harmonium) harmonium->setConvolutionEnabled(enabled);
    }

    bool isConvolutionEnabled() const {
        return harmonium ? harmonium->isConvolutionEnabled() : false;
    }
    
#ifndef JucePlugin_PreferredChannelConfigurations
    bool isBusesLayoutSupported(const BusesLayout& layouts) const override;
#endif

private:
    std::unique_ptr<Harmonium> harmonium;
    juce::MidiMessageCollector midiMessageCollector;
    std::atomic<float> masterGain{0.7f};
    std::atomic<float> airCapacity{0.5f};
    std::atomic<float> reedChamberCapacity{0.5f};
    std::atomic<float> narrowJetCapacity{0.5f};
    std::atomic<float> qFactor{0.5f};
    std::atomic<float> airConsumption{0.5f};
    std::atomic<bool> polyphonic{true};
    std::atomic<int> reedMode{1};
    
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(HarmoniumPhysicsEngineAudioProcessor)
};
