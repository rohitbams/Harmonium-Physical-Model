// HarmoniumPhysicsState.h
#pragma once
#include <cmath>

// Pure data structure - all Puranik-Scavone variables in one place
struct HarmoniumPhysicsState {
    // Core physics state (Puranik-Scavone equations)
    double p0 = 0.0; // Bellows chamber pressure
    double p1 = 0.0;                    // reed chamber pressure (Pa)
    double p2 = 0.0;                    // jet pressure (Pa)
    double vj = 0.0;                    // jet velocity (m/s)
    double reedPosition = 0.0;          // reed displacement ζ (m)
    double reedVelocity = 0.0;          // reed velocity dζ/dt (m/s)
    double u = 0.0;                     // total volume flow (m³/s)
    double u_previous = 0.0;            // previous flow for derivatives
    
    // Computed values (for diagnostics)
    double drivingForce = 0.0;          // μ(p₂ - p_atm) (N)
    double springForce = 0.0;           // ω₀²ζ (N)
    double dampingForce = 0.0;          // (ω₀/Q)(dζ/dt) (N)
    double acceleration = 0.0;          // Total acceleration (m/s²)
    double aperture = 0.0;              // current reed aperture (m²)
    double netFlow = 0.0;               // u₀ - u (m³/s)
    
    // Reset to initial state
    void reset() {
        *this = HarmoniumPhysicsState{};
        reedVelocity = 1e-4;  // Small initial velocity to seed oscillation
    }
    
    // Utility methods
    bool isOscillating() const {
        return std::abs(reedPosition) > 1e-5 || std::abs(reedVelocity) > 1e-2;
    }
    
    //double getAmplitude() const {
    //    return std::abs(reedPosition);
    //}
};

// Physics configuration - single source of truth for all constants
struct HarmoniumPhysicsConfig {
    // Simulation
    double sampleRate = 44100.0;
    double dt = 1.0 / sampleRate;
    
    // Physical constants (Puranik-Scavone paper values)
    double rho0 = 1.225;                // Air density (kg/m³)
    double c0 = 343.0;                  // Sound speed (m/s)
    double p_atm = 0.0;                 // Atmospheric pressure (Pa, gauge)
    
    // Component parameters
    double V1 = 0.1;                   // Chamber volume (m³)
    double L2 = 0.005;                  // Jet length (m)
    double S2 = 0.0001;                 // Jet cross-section (m²)
    double Q = 50.0;                    // Quality factor
    double mu = 5e-2;                   // Pressure-to-force coupling
    double Sr = 4e-4;                   // Reed cross-section area (m²)
    double alpha = 0.8;                 // Flow coefficient
    double reedWidth = 0.02;            // Reed width (m)
    
    // Audio
    int amplitudeScaleFactor = 1000;    // Audio scaling
    
    // Advanced parameters
    double minOscillationLevel = 1e-6;
    double targetAmplitude = 5e-5;
    double rest_gap = 3e-5;
};
