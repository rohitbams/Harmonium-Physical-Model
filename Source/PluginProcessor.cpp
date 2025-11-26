#include "PluginProcessor.h"
#include "PluginEditor.h"

/*
 * PluginProcessor class handles audio and MIDI processing.
 */
HarmoniumPhysicsEngineAudioProcessor::HarmoniumPhysicsEngineAudioProcessor()
#ifndef JucePlugin_PreferredChannelConfigurations
     : AudioProcessor (BusesProperties()
                     #if ! JucePlugin_IsMidiEffect
                      #if ! JucePlugin_IsSynth
                       .withInput  ("Input",  juce::AudioChannelSet::stereo(), true)
                      #endif
                       .withOutput ("Output", juce::AudioChannelSet::stereo(), true)
                     #endif
                       )
#endif
{
    DBG("HarmoniumPhysicsEngine processor initialized");
}

HarmoniumPhysicsEngineAudioProcessor::~HarmoniumPhysicsEngineAudioProcessor() {
}

const juce::String HarmoniumPhysicsEngineAudioProcessor::getName() const {
    return JucePlugin_Name;
}

bool HarmoniumPhysicsEngineAudioProcessor::acceptsMidi() const {
#if JucePlugin_WantsMidiInput
    return true;
#else
    return false;
#endif
}

bool HarmoniumPhysicsEngineAudioProcessor::producesMidi() const {
#if JucePlugin_ProducesMidiOutput
    return true;
#else
    return false;
#endif
}

bool HarmoniumPhysicsEngineAudioProcessor::isMidiEffect() const {
#if JucePlugin_IsMidiEffect
    return true;
#else
    return false;
#endif
}

double HarmoniumPhysicsEngineAudioProcessor::getTailLengthSeconds() const {
    return 0.0;
}

int HarmoniumPhysicsEngineAudioProcessor::getNumPrograms() {
    return 1;
}

int HarmoniumPhysicsEngineAudioProcessor::getCurrentProgram() {
    return 0;
}

void HarmoniumPhysicsEngineAudioProcessor::setCurrentProgram(int index) {
}

const juce::String HarmoniumPhysicsEngineAudioProcessor::getProgramName(int index) {
    return {};
}

void HarmoniumPhysicsEngineAudioProcessor::changeProgramName(int index, const juce::String& newName) {
}

void HarmoniumPhysicsEngineAudioProcessor::getStateInformation(juce::MemoryBlock& destData) {
}

void HarmoniumPhysicsEngineAudioProcessor::setStateInformation(const void* data, int sizeInBytes) {
}

#ifndef JucePlugin_PreferredChannelConfigurations
bool HarmoniumPhysicsEngineAudioProcessor::isBusesLayoutSupported(const BusesLayout& layouts) const {
#if JucePlugin_IsMidiEffect
    juce::ignoreUnused(layouts);
    return true;
#else
    if (layouts.getMainOutputChannelSet() != juce::AudioChannelSet::mono()
     && layouts.getMainOutputChannelSet() != juce::AudioChannelSet::stereo())
        return false;

#if ! JucePlugin_IsSynth
    if (layouts.getMainOutputChannelSet() != layouts.getMainInputChannelSet())
        return false;
#endif

    return true;
#endif
}
#endif

void HarmoniumPhysicsEngineAudioProcessor::prepareToPlay(double sampleRate, int samplesPerBlock) {
    DBG("Preparing to play: SR=" + juce::String(sampleRate) + ", Block=" + juce::String(samplesPerBlock));
    
    harmonium = std::make_unique<Harmonium>(sampleRate);
    
    if (harmonium) {
          harmonium->prepareToPlay(sampleRate, samplesPerBlock);
      }
      
    midiMessageCollector.reset(sampleRate);
    DBG("Harmonium Physics Engine ready for real-time processing");
}

void HarmoniumPhysicsEngineAudioProcessor::releaseResources() {
    DBG("Releasing resources");
}

void HarmoniumPhysicsEngineAudioProcessor::processBlock(juce::AudioBuffer<float>& buffer, juce::MidiBuffer& midiMessages) {
    juce::ScopedNoDenormals noDenormals;
    
    for (const auto metadata : midiMessages) {
        auto message = metadata.getMessage();
        
        if (message.isNoteOn()) {
            double frequency = juce::MidiMessage::getMidiNoteInHertz(message.getNoteNumber());
            startNote(frequency);
        }
        else if (message.isNoteOff()) {
            double frequency = juce::MidiMessage::getMidiNoteInHertz(message.getNoteNumber());
            stopNote(frequency);
        }
        else if (message.isController() && message.getControllerNumber() == 1) {
            
            double modWheelValue = message.getControllerValue() / 127.0;
            setModWheelValue(modWheelValue);
        }
    }
    
    buffer.clear();
    
    if (!harmonium) return;
    
    auto numChannels = buffer.getNumChannels();
    auto numSamples = buffer.getNumSamples();
    
    for (int sample = 0; sample < numSamples; ++sample) {
        float harmoniumOutput = harmonium->processSample();
        harmoniumOutput *= masterGain.load();

        for (int channel = 0; channel < numChannels; ++channel) {
            buffer.setSample(channel, sample, harmoniumOutput);
        }
    }
}

void HarmoniumPhysicsEngineAudioProcessor::startNote(double frequency) {
    if (harmonium) {
        harmonium->startNote(frequency);
        DBG("Note started: " + juce::String(frequency) + " Hz");
    }
}

void HarmoniumPhysicsEngineAudioProcessor::stopNote(double frequency) {
    if (harmonium) {
        harmonium->stopNote(frequency);
        DBG("Note stopped: " + juce::String(frequency) + " Hz");
    }
}

void HarmoniumPhysicsEngineAudioProcessor::stopAllNotes() {
    if (harmonium) {
        harmonium->stopAllNotes();
        DBG("All notes stopped");
    }
}

bool HarmoniumPhysicsEngineAudioProcessor::isNoteCurrentlyOn() const {
    return harmonium && harmonium->isNoteActive();
}

void HarmoniumPhysicsEngineAudioProcessor::setAmplitudeScaling(int amp) {
    if (harmonium) {
        harmonium->setAmplitudeScale(amp);
        DBG("Amplitude scale: " + juce::String(amp));
    }
}


// EQ
//void HarmoniumPhysicsEngineAudioProcessor::setLows(int value) {
//    if (harmonium) {
//        harmonium->setLows(value);
//        DBG("lows: " + juce::String(value));
//    }
//}


void HarmoniumPhysicsEngineAudioProcessor::setAirCapacity(float value) {
    airCapacity = value;
    if (harmonium) {
        harmonium->setAirCapacity(value);
        DBG("Air capacity: " + juce::String(value, 2));
    }
}

void HarmoniumPhysicsEngineAudioProcessor::setReedChamberCapacity(float value) {
    reedChamberCapacity = value;
    if (harmonium) {
        harmonium->setReedChamberCapacity(value);
        DBG("airCapacitySlider: " + juce::String(value, 2));
    }
}

void HarmoniumPhysicsEngineAudioProcessor::setNarrowJetCapacity(float value) {
    narrowJetCapacity = value;
    if (harmonium) {
        harmonium->setNarrowJetCapacity(value);
        DBG("NJCapacitySlider: " + juce::String(value, 2));
    }
}


void HarmoniumPhysicsEngineAudioProcessor::setAirConsumption(float value) {
    airConsumption = value;
    if (harmonium) {
        harmonium->setAirConsumption(value);
        DBG("Air consumption: " + juce::String(value, 2));
    }
}

void HarmoniumPhysicsEngineAudioProcessor::setMasterGain(float value) {
    masterGain = value;
    DBG("Master gain: " + juce::String(value, 2));
}

void HarmoniumPhysicsEngineAudioProcessor::setReedMode(int mode) {
    reedMode = juce::jlimit(1, 3, mode);
    if (harmonium) {
        harmonium->setReedMode(static_cast<Voice::ReedMode>(mode));
        DBG("Reed mode: " + juce::String(mode == 1 ? "Single" :
                                       mode == 2 ? "Double" : "Triple"));
    }
}

int HarmoniumPhysicsEngineAudioProcessor::getReedMode() const {
    return reedMode.load();
}

double HarmoniumPhysicsEngineAudioProcessor::getBellowsPressure() const {
    return harmonium ? harmonium->getPhysicsState().p0 : 0.0;
}

double HarmoniumPhysicsEngineAudioProcessor::getChamberPressure() const {
    return harmonium ? harmonium->getPhysicsState().p1 : 0.0;
}

double HarmoniumPhysicsEngineAudioProcessor::getJetPressure() const {
    return harmonium ? harmonium->getPhysicsState().p2 : 0.0;
}

double HarmoniumPhysicsEngineAudioProcessor::getReedPosition() const {
    return harmonium ? harmonium->getPhysicsState().reedPosition : 0.0;
}

double HarmoniumPhysicsEngineAudioProcessor::getReedVelocity() const {
    return harmonium ? harmonium->getPhysicsState().reedVelocity : 0.0;
}

double HarmoniumPhysicsEngineAudioProcessor::getJetVelocity() const {
    return harmonium ? harmonium->getPhysicsState().vj : 0.0;
}

int HarmoniumPhysicsEngineAudioProcessor::getTotalOscillatingReeds() const {
    return harmonium ? harmonium->getTotalOscillatingReeds() : 0;
}

//double HarmoniumPhysicsEngineAudioProcessor::getAverageReedAmplitude() const {
//    return harmonium ? harmonium->getAverageReedAmplitude() : 0.0;
//}

void HarmoniumPhysicsEngineAudioProcessor::handleMidiCC(int ccNumber, int ccValue) {
    if (ccNumber == 1 && ccValue >= 0 && ccValue <= 127) {
        double modWheelValue = ccValue / 127.0;
        setModWheelValue(modWheelValue);
    }
}

void HarmoniumPhysicsEngineAudioProcessor::setModWheelValue(double value) {
    if (harmonium && std::isfinite(value) && value >= 0.0 && value <= 1.0) {
        harmonium->setModWheelValue(value);
    }
}

juce::AudioProcessorEditor* HarmoniumPhysicsEngineAudioProcessor::createEditor() {
    return new HarmoniumPhysicsEngineAudioProcessorEditor(*this);
}

bool HarmoniumPhysicsEngineAudioProcessor::hasEditor() const {
    return true;
}

juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter() {
    return new HarmoniumPhysicsEngineAudioProcessor();
}
