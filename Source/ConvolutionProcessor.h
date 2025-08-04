#pragma once
#include <JuceHeader.h>

class ConvolutionProcessor {
private:
    juce::dsp::Convolution convolution;
    juce::dsp::ProcessSpec spec;
    bool isEnabled = false;
    
    // Buffer for block processing
    juce::AudioBuffer<float> tempBuffer;
    
public:
    void prepare(double sampleRate, int blockSize) {
        spec.sampleRate = sampleRate;
        spec.maximumBlockSize = blockSize;
        spec.numChannels = 1;
        
        convolution.prepare(spec);
        
        // Prepare temp buffer for single sample processing
        tempBuffer.setSize(1, 1);
    }
    
    bool loadIR(const juce::File& irFile) {
        if (!irFile.existsAsFile()) return false;
        
        convolution.loadImpulseResponse(irFile,
                                      juce::dsp::Convolution::Stereo::no,
                                      juce::dsp::Convolution::Trim::yes,
                                      0,
                                      juce::dsp::Convolution::Normalise::yes);
        
        isEnabled = true;
        return true;
    }
    
    float processSample(float inputSample) {
        if (!isEnabled) return inputSample;
        
        // Use temp buffer for single sample
        tempBuffer.setSample(0, 0, inputSample);
        
        auto block = juce::dsp::AudioBlock<float>(tempBuffer);
        auto context = juce::dsp::ProcessContextReplacing<float>(block);
        convolution.process(context);
        
        return tempBuffer.getSample(0, 0);
    }
    
    void setEnabled(bool enabled) { isEnabled = enabled; }
    bool getEnabled() const { return isEnabled; }
};
