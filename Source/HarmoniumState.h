// HarmoniumState.h
#pragma once
#include <cmath>

/*
 * HarmoniumState class.
 *
 */
 struct HarmoniumState {
     // core physics state
     double p0 = 0.0; // bellows chamber pressure (Pa)
     double p1 = 0.0; // reed chamber pressure (Pa)
     double p2 = 0.0; // narrow jet pressure (Pa)
     double vj = 0.0; // narrow jet velocity (m/s)
     double reedPosition = 0.0; // reed displacement ζ (m)
     double reedVelocity = 0.0; // reed velocity dζ/dt (m/s)
     double u = 0.0; // total volume flow (m³/s)
     double u_previous = 0.0; // previous flow for derivatives
     double drivingForce = 0.0; // μ(p₂ - p_atm) (N)
     double springForce = 0.0; // ω₀²ζ (N)
     double dampingForce = 0.0; // (ω₀/Q)(dζ/dt) (N)
     double acceleration = 0.0; // Total acceleration (m/s²)
     double aperture = 0.0; // current reed aperture (m²)
     double netFlow = 0.0; // u₀ - u (m³/s)
    
    void reset() {
        *this = HarmoniumState{};
//        reedVelocity = 1e-4;
    }
    
    bool isOscillating() const {
        return std::abs(reedPosition) > 1e-5 || std::abs(reedVelocity) > 1e-2;
    }
  
};

struct HarmoniumPhysicsConfig {
    double sampleRate = 44100.0;
    double dt = 1.0 / sampleRate;
    
    double rho0 = 1.225; // air density (kg/m³)
    double c0 = 343.0; // sound speed (m/s)
    double p_atm = 0.0; // atmospheric pressure (Pa, gauge)
    double V1 = 0.5; // reed hamber volume (m³)
    double L2 = 0.005; // jet length (m)
    double S2 = 0.0001; // jet cross-section (m²)
    double Q = 50.0; // quality factor
    double mu = 5e-2; // pressure-to-force coupling
    double Sr = 4e-4; // eeed cross-section area (m²)
    double alpha = 0.8; // flow coefficient
    double reedWidth = 0.02; // reed width (m)
    int amplitudeScale = 50;
    double minOscillationLevel = 1e-6;
    double targetAmplitude = 5e-5;
    double restGap = 3e-5;
    double saturationFactor = 0.1;
    
    double low = 2.0;       // 15 - 200
    double lowMid = 1.0;    // 200 - 400
    double mid = 1.0;       // 400 - 900
    double highMid = 1.0;   // 900 - 4000
    double high = 1.0;      // 4000 - 20000
    
    // new small tuning params
    double nonlinK = 1.0;        // cubic nonlinearity strength (0.. few)
    double asymmetryK = 0.5;     // asymmetry gain producing even harmonics
    double hfBoost = 5.0;        // multiplier for HF differentiator (brightness)
    double noiseLevel = 1e-5;    // breath/noise level (added to output)
    double outputGain = 1.0;     // final gain applied to output
};
