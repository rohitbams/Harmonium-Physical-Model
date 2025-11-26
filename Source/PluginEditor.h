#pragma once
#include <JuceHeader.h>
#include "PluginProcessor.h"

/*
 * PluginEditor class handles GUI implementation.
 */
class HarmoniumPhysicsEngineAudioProcessorEditor : public juce::AudioProcessorEditor,
                                                   public juce::Slider::Listener,
                                                   public juce::MidiKeyboardStateListener,
                                                   public juce::Timer,
                                                   public juce::ComboBox::Listener
{
public:
    HarmoniumPhysicsEngineAudioProcessorEditor(HarmoniumPhysicsEngineAudioProcessor&);
    ~HarmoniumPhysicsEngineAudioProcessorEditor() override;

    void paint(juce::Graphics&) override;
    void resized() override;
    void timerCallback() override;

private:
    
    std::unique_ptr<juce::FileChooser> fileChooser;
    HarmoniumPhysicsEngineAudioProcessor& audioProcessor;
    
    juce::Slider amplitudeSlider;
    juce::Label amplitudeLabel;
    
    juce::Slider lowSlider;
    juce::Label lowLabel;
    juce::Slider lowMidSlider;
    juce::Label lowMidLabel;
    juce::Slider midSlider;
    juce::Label midLabel;
    juce::Slider highMidSlider;
    juce::Label highMidLabel;
    juce::Slider highSlider;
    juce::Label highLabel;
    
    juce::Slider airCapacitySlider;
    juce::Label airCapacityLabel;
    
    juce::Slider airConsumptionSlider;
    juce::Label airConsumptionLabel;

    juce::Slider reedChamberCapacitySlider;
    juce::Label reedChamberCapacityLabel;
    
    juce::Slider narrowJetCapacitySlider;
    juce::Label narrowJetCapacityLabel;

//    juce::Slider qSlider;
//    juce::Label qLabel;
    
    juce::Slider masterGainSlider;
    juce::Label masterGainLabel;
    
    juce::ComboBox reedModeSelector;
    juce::Label reedModeLabel;
    
    juce::TextButton loadIRButton;
    juce::ToggleButton convolutionToggle;
    juce::Label irStatusLabel;
    juce::Label convolutionLabel;
    
    juce::Label physicsTitle;
    juce::Label p0Display;
    juce::Label p1Display;
    juce::Label p2Display;
    juce::Label reedPositionDisplay;
    juce::Label reedVelocityDisplay;
    juce::Label oscillatingReedsDisplay;
    juce::Label modWheelDisplay;

    void handleNoteOn(juce::MidiKeyboardState*, int midiChannel, int midiNoteNumber, float velocity) override;
    void handleNoteOff(juce::MidiKeyboardState*, int midiChannel, int midiNoteNumber, float velocity) override;
    void sliderValueChanged(juce::Slider* slider) override;
    void comboBoxChanged(juce::ComboBox* comboBoxThatHasChanged) override;
    void loadIRFile();

    
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(HarmoniumPhysicsEngineAudioProcessorEditor)
};

// PluginEditor.cpp

#include "PluginEditor.h"

HarmoniumPhysicsEngineAudioProcessorEditor::HarmoniumPhysicsEngineAudioProcessorEditor(HarmoniumPhysicsEngineAudioProcessor& p) : AudioProcessorEditor(&p), audioProcessor(p) {
    
    // Amplitude control
    amplitudeSlider.setRange(100, 3000);
    amplitudeSlider.setValue(1000);
    amplitudeSlider.setTextBoxStyle(juce::Slider::TextBoxRight, false, 60, 20);
    amplitudeSlider.addListener(this);
    addAndMakeVisible(amplitudeSlider);
    
    amplitudeLabel.setText("Amplitude Scale", juce::dontSendNotification);
    amplitudeLabel.attachToComponent(&amplitudeSlider, true);
    addAndMakeVisible(amplitudeLabel);

    // EQ control
    //low
    lowSlider.setRange(100, 1000);
    lowSlider.setValue(500);
    lowSlider.setTextBoxStyle(juce::Slider::TextBoxRight, false, 60, 20);
    lowSlider.addListener(this);
    addAndMakeVisible(lowSlider);
    
    lowLabel.setText("Lows", juce::dontSendNotification);
    lowLabel.attachToComponent(&lowSlider, true);
    addAndMakeVisible(lowLabel);
    
    // Air capacity control
    airCapacitySlider.setSliderStyle(juce::Slider::LinearHorizontal);
    airCapacitySlider.setRange(0.0, 1.0, 0.01);
    airCapacitySlider.setValue(0.5);
    airCapacitySlider.addListener(this);
    addAndMakeVisible(airCapacitySlider);

    airCapacityLabel.setText("Bellows Air Capacity", juce::dontSendNotification);
    airCapacityLabel.attachToComponent(&airCapacitySlider, true);
    addAndMakeVisible(airCapacityLabel);

    // Air consumption control
    airConsumptionSlider.setSliderStyle(juce::Slider::LinearHorizontal);
    airConsumptionSlider.setRange(0.0, 1.0, 0.01);
    airConsumptionSlider.setValue(0.5);
    airConsumptionSlider.addListener(this);
    addAndMakeVisible(airConsumptionSlider);

    airConsumptionLabel.setText("Bellows Air Consumption", juce::dontSendNotification);
    airConsumptionLabel.attachToComponent(&airConsumptionSlider, true);
    addAndMakeVisible(airConsumptionLabel);

    // reedChamber capacity control
    reedChamberCapacitySlider.setSliderStyle(juce::Slider::LinearHorizontal);
    reedChamberCapacitySlider.setRange(0.0, 1.0, 0.01);
    reedChamberCapacitySlider.setValue(0.5);
    reedChamberCapacitySlider.addListener(this);
    addAndMakeVisible(reedChamberCapacitySlider);

    reedChamberCapacityLabel.setText("Reed Chamber Capacity", juce::dontSendNotification);
    reedChamberCapacityLabel.attachToComponent(&reedChamberCapacitySlider, true);
    addAndMakeVisible(reedChamberCapacityLabel);
    
    // narrowJet capacity control
    narrowJetCapacitySlider.setSliderStyle(juce::Slider::LinearHorizontal);
    narrowJetCapacitySlider.setRange(0.0, 1.0, 0.01);
    narrowJetCapacitySlider.setValue(0.5);
    narrowJetCapacitySlider.addListener(this);
    addAndMakeVisible(narrowJetCapacitySlider);

    narrowJetCapacityLabel.setText("Narrow Jet Capacity", juce::dontSendNotification);
    narrowJetCapacityLabel.attachToComponent(&narrowJetCapacitySlider, true);
    addAndMakeVisible(narrowJetCapacityLabel);

    // q control
//    qSlider.setSliderStyle(juce::Slider::LinearHorizontal);
//    qSlider.setRange(0.0, 1.0, 0.01);
//    qSlider.setValue(0.5);
//    qSlider.addListener(this);
//    addAndMakeVisible(qSlider);
//
//    qLabel.setText("Q factor", juce::dontSendNotification);
//    qLabel.attachToComponent(&qSlider, true);
//    addAndMakeVisible(qLabel);
    
    // Master gain control
    masterGainSlider.setSliderStyle(juce::Slider::LinearHorizontal);
    masterGainSlider.setRange(0.0, 1.0, 0.01);
    masterGainSlider.setValue(0.2);
    masterGainSlider.addListener(this);
    addAndMakeVisible(masterGainSlider);

    masterGainLabel.setText("Master Gain", juce::dontSendNotification);
    masterGainLabel.attachToComponent(&masterGainSlider, true);
    addAndMakeVisible(masterGainLabel);

    // reed mode selector
    reedModeSelector.addItem("Single Reed", 1);
    reedModeSelector.addItem("Double Reed", 2);
    reedModeSelector.addItem("Triple Reed", 3);
    reedModeSelector.setSelectedId(1, juce::dontSendNotification);
    reedModeSelector.addListener(this);
    addAndMakeVisible(reedModeSelector);

    reedModeLabel.setText("Reed Mode:", juce::dontSendNotification);
    reedModeLabel.attachToComponent(&reedModeSelector, true);
    addAndMakeVisible(reedModeLabel);
    
    // Physics displays
//    physicsTitle.setText("State variables", juce::dontSendNotification);
    physicsTitle.setFont(juce::FontOptions(16.0f, juce::Font::bold));
    addAndMakeVisible(physicsTitle);
    
    p0Display.setText("p0 (Bellows): -- Pa", juce::dontSendNotification);
    addAndMakeVisible(p0Display);
    
    p1Display.setText("p1 (Reed chamber): -- Pa", juce::dontSendNotification);
    addAndMakeVisible(p1Display);
    
    p2Display.setText("p2 (Narrow jet): -- Pa", juce::dontSendNotification);
    addAndMakeVisible(p2Display);
    
    reedPositionDisplay.setText("Reed Position: -- m", juce::dontSendNotification);
    addAndMakeVisible(reedPositionDisplay);
    
    reedVelocityDisplay.setText("Reed Velocity: -- m/s", juce::dontSendNotification);
    addAndMakeVisible(reedVelocityDisplay);
    
    oscillatingReedsDisplay.setText("Oscillating Reeds: 0", juce::dontSendNotification);
    addAndMakeVisible(oscillatingReedsDisplay);
    
    loadIRButton.setButtonText("Load IR");
    loadIRButton.onClick = [this] { loadIRFile(); };
    addAndMakeVisible(loadIRButton);

    convolutionToggle.setButtonText("Convolution");
    convolutionToggle.onStateChange = [this] {
        audioProcessor.setConvolutionEnabled(convolutionToggle.getToggleState());
    };
    addAndMakeVisible(convolutionToggle);

    irStatusLabel.setText("No IR loaded", juce::dontSendNotification);
    addAndMakeVisible(irStatusLabel);

    convolutionLabel.setText("Resonant Body:", juce::dontSendNotification);
    addAndMakeVisible(convolutionLabel);

    
    startTimer(50);
    
    setSize(900, 500);
    
}

HarmoniumPhysicsEngineAudioProcessorEditor::~HarmoniumPhysicsEngineAudioProcessorEditor() {}

void HarmoniumPhysicsEngineAudioProcessorEditor::paint(juce::Graphics& g) {
    g.fillAll(getLookAndFeel().findColour(juce::ResizableWindow::backgroundColourId));
    
    g.setColour(juce::Colours::white);
    g.setFont(20.0f);
    g.drawFittedText("Harmonium",
                     getLocalBounds().removeFromTop(30),
                     juce::Justification::centred, 1);
}

void HarmoniumPhysicsEngineAudioProcessorEditor::resized() {
    auto area = getLocalBounds();
    area.removeFromTop(35);
    
    const int sliderHeight = 30;
    const int labelWidth = 140;
    const int margin = 5;
    
    // controls section (left side)
    auto controlsArea = area.removeFromLeft(400);
    controlsArea.removeFromTop(5);
    
    // row 1: Amplitude
    auto ampArea = controlsArea.removeFromTop(sliderHeight);
    ampArea.removeFromLeft(labelWidth);
    amplitudeSlider.setBounds(ampArea);
    controlsArea.removeFromTop(margin);
    
    // row 2: Air capacity
    auto airCapArea = controlsArea.removeFromTop(sliderHeight);
    airCapArea.removeFromLeft(labelWidth);
    airCapacitySlider.setBounds(airCapArea);
    controlsArea.removeFromTop(margin);

    // row 3: Air consumption
    auto airConsArea = controlsArea.removeFromTop(sliderHeight);
    airConsArea.removeFromLeft(labelWidth);
    airConsumptionSlider.setBounds(airConsArea);
    controlsArea.removeFromTop(margin);
    
    // row 4: Master gain
    auto gainArea = controlsArea.removeFromTop(sliderHeight);
    gainArea.removeFromLeft(labelWidth);
    masterGainSlider.setBounds(gainArea);
    controlsArea.removeFromTop(margin);
    
    // row 5: Reed mode
    auto reedModeArea = controlsArea.removeFromTop(sliderHeight);
    reedModeArea.removeFromRight(labelWidth);
    reedModeArea.removeFromLeft(labelWidth);
    reedModeSelector.setBounds(reedModeArea);
    controlsArea.removeFromTop(margin);
    
    // row 7: reedChamber capacity
    auto reedChamberCapArea = controlsArea.removeFromTop(sliderHeight);
    reedChamberCapArea.removeFromLeft(labelWidth);
    reedChamberCapacitySlider.setBounds(reedChamberCapArea);
    controlsArea.removeFromTop(margin);

    // row 8: narrow jet capacity
    auto narrowJetCapArea = controlsArea.removeFromTop(sliderHeight);
    narrowJetCapArea.removeFromLeft(labelWidth);
    narrowJetCapacitySlider.setBounds(narrowJetCapArea);
    controlsArea.removeFromTop(margin);

//    // row 9: q factor
//    auto qArea = controlsArea.removeFromTop(sliderHeight);
//    qArea.removeFromLeft(labelWidth);
//    qSlider.setBounds(qArea);
//    controlsArea.removeFromTop(margin);
    
    // row 9: EQ
    auto lowArea = controlsArea.removeFromTop(sliderHeight);
    lowArea.removeFromLeft(labelWidth);
    lowSlider.setBounds(lowArea);
    controlsArea.removeFromTop(margin);

    
    // display section (right side)
    auto physicsArea = area.removeFromLeft(400);
    physicsArea.removeFromTop(5);
    
    physicsTitle.setBounds(physicsArea.removeFromTop(25));
    physicsArea.removeFromTop(5);
    
    p0Display.setBounds(physicsArea.removeFromTop(20));
    p1Display.setBounds(physicsArea.removeFromTop(20));
    p2Display.setBounds(physicsArea.removeFromTop(20));
    reedPositionDisplay.setBounds(physicsArea.removeFromTop(20));
    reedVelocityDisplay.setBounds(physicsArea.removeFromTop(20));
    oscillatingReedsDisplay.setBounds(physicsArea.removeFromTop(20));
    physicsArea.removeFromTop(10);
    modWheelDisplay.setBounds(physicsArea.removeFromTop(40));
    
    auto irArea = area.removeFromLeft(250);
    irArea.removeFromTop(5);

    convolutionLabel.setBounds(irArea.removeFromTop(20));
    loadIRButton.setBounds(irArea.removeFromTop(30));
    irArea.removeFromTop(5);
    convolutionToggle.setBounds(irArea.removeFromTop(25));
    irStatusLabel.setBounds(irArea.removeFromTop(40));

}

void HarmoniumPhysicsEngineAudioProcessorEditor::timerCallback() {
    double p0 = audioProcessor.getBellowsPressure();
    double p1 = audioProcessor.getChamberPressure();
    double p2 = audioProcessor.getJetPressure();
    double reedPosition = audioProcessor.getReedPosition();
    double reedVelocity = audioProcessor.getReedVelocity();
    int oscillatingReeds = audioProcessor.getTotalOscillatingReeds();
//    double low = audioProcessor.getLows();
//    double lowMids = audioProcessor.getLowMids();
//    double mids = audioProcessor.getMids();
//    double highMids = audioProcessor.getHighMids();
//    double highs = audioProcessor.getHighs();

    p0Display.setText("AirMass (Bellows): " + juce::String(p0, 1) + " Pa", juce::dontSendNotification);
    p1Display.setText("p1 (Reed chamber): " + juce::String(p1, 1) + " Pa", juce::dontSendNotification);
    p2Display.setText("p2 (Narrow jet): " + juce::String(p2, 1) + " Pa", juce::dontSendNotification);
    reedPositionDisplay.setText("Reed Position: " + juce::String(reedPosition, 6) + " m", juce::dontSendNotification);
    reedVelocityDisplay.setText("Reed Velocity: " + juce::String(reedVelocity, 3) + " m/s", juce::dontSendNotification);
    oscillatingReedsDisplay.setText("Oscillating Reeds: " + juce::String(oscillatingReeds), juce::dontSendNotification);
}

void HarmoniumPhysicsEngineAudioProcessorEditor::handleNoteOn(juce::MidiKeyboardState*, int midiChannel, int midiNoteNumber, float velocity) {
    double frequency = juce::MidiMessage::getMidiNoteInHertz(midiNoteNumber);
    audioProcessor.startNote(frequency);
    DBG("GUI Note ON: " + juce::String(midiNoteNumber) + " (" + juce::String(frequency, 1) + " Hz)");
}

void HarmoniumPhysicsEngineAudioProcessorEditor::handleNoteOff(juce::MidiKeyboardState*, int midiChannel, int midiNoteNumber, float velocity) {
    double frequency = juce::MidiMessage::getMidiNoteInHertz(midiNoteNumber);
    audioProcessor.stopNote(frequency);
    DBG("GUI Note OFF: " + juce::String(midiNoteNumber));
}

void HarmoniumPhysicsEngineAudioProcessorEditor::sliderValueChanged(juce::Slider* slider) {
    if (slider == &amplitudeSlider) {
        audioProcessor.setAmplitudeScaling(static_cast<int>(amplitudeSlider.getValue()));
    }
//    if (slider == &lowSlider) {
//        audioProcessor.setLows(static_cast<int>(lowSlider.getValue()));
//    }
//    else if (slider == &qSlider) {
//        audioProcessor.setQFactor(static_cast<float>(qSlider.getValue()));
//    }
    else if (slider == &airCapacitySlider) {
        audioProcessor.setAirCapacity(static_cast<float>(airCapacitySlider.getValue()));
    }
    else if (slider == &airConsumptionSlider) {
        audioProcessor.setAirConsumption(static_cast<float>(airConsumptionSlider.getValue()));
    }
    else if (slider == &masterGainSlider) {
        audioProcessor.setMasterGain(static_cast<float>(masterGainSlider.getValue()));
    }
    else if (slider == &narrowJetCapacitySlider) {
        audioProcessor.setNarrowJetCapacity(static_cast<float>(narrowJetCapacitySlider.getValue()));
    }
}

void HarmoniumPhysicsEngineAudioProcessorEditor::loadIRFile() {
    auto chooserFlags = juce::FileBrowserComponent::openMode | juce::FileBrowserComponent::canSelectFiles;
    
    fileChooser = std::make_unique<juce::FileChooser>("Select IR file", juce::File(), "*.wav");
    
    fileChooser->launchAsync(chooserFlags, [this](const juce::FileChooser& fc) {
        auto file = fc.getResult();
        if (file != juce::File{}) {
            if (audioProcessor.loadImpulseResponse(file)) {
                irStatusLabel.setText("Loaded: " + file.getFileNameWithoutExtension(),
                                     juce::dontSendNotification);
                convolutionToggle.setToggleState(true, juce::sendNotification);
            } else {
                irStatusLabel.setText("Failed to load IR", juce::dontSendNotification);
            }
        }
    });
}

void HarmoniumPhysicsEngineAudioProcessorEditor::comboBoxChanged(juce::ComboBox* comboBoxThatHasChanged) {
    if (comboBoxThatHasChanged == &reedModeSelector) {
        int selectedMode = reedModeSelector.getSelectedId();
        audioProcessor.setReedMode(selectedMode);
    }
}
