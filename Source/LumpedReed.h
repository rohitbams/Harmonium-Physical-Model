#pragma once
#include "ReedModel.h"
#include <cmath>
#include <algorithm>

class LumpedReed : public ReedModel {
private:
    struct ReedParams {
        double omega0;              // Angular frequency (rad/s)
        double Q = 25.0;           // Quality factor
        double mu = 5e-2;          // Coupling coefficient
        double Sr = 4e-4;          // Reed area (m²)
    } params;
    
    double sampleRate;
    double dt;
    
    // State variables (your existing implementation)
    double reedPosition = 0.0;     // ζ (m)
    double reedVelocity = 0.0;     // dζ/dt (m/s)
    double reedAcceleration = 0.0; // d²ζ/dt² (m/s²)
    
    PhysicsState currentState;

public:
    LumpedReed(double targetFreq, double sampleRate)
        : sampleRate(sampleRate), dt(1.0/sampleRate) {
        // Initialize with your existing parameter values
        params.omega0 = 2.0 * M_PI * targetFreq;
        reset();
    }
    
    void updateMotion(double pressureDifference, double dt_override = 0.0) override {
        double time_step = (dt_override > 0.0) ? dt_override : dt;
        
        // Your existing reed equation: d²ζ/dt² + Q⁻¹ω₀ dζ/dt + ω₀²ζ = μ(p₂ - pₐₜₘ)
        double dampingTerm = (1.0 / params.Q) * params.omega0 * reedVelocity;
        double springTerm = params.omega0 * params.omega0 * reedPosition;
        double drivingForce = params.mu * pressureDifference;
        
        reedAcceleration = drivingForce - dampingTerm - springTerm;
        
        // Limit acceleration for stability
        reedAcceleration = std::clamp(reedAcceleration, -1000.0, 1000.0);
        
        // Euler integration (your existing method)
        reedVelocity += reedAcceleration * time_step;
        reedPosition += reedVelocity * time_step;
        
        // Apply physical limits
        reedPosition = std::clamp(reedPosition, -0.01, 0.01);
        reedVelocity = std::clamp(reedVelocity, -10.0, 10.0);
        
        // Update state
        updatePhysicsState();
    }
    
    void reset() override {
        reedPosition = 0.0;
        reedVelocity = 0.0;
        reedAcceleration = 0.0;
        currentState = PhysicsState{};
    }
    
    PhysicsState getState() const override {
        return currentState;
    }
    
    void setFrequency(double frequency) override {
        params.omega0 = 2.0 * M_PI * frequency;
    }
    
    double getTipDisplacement() const override {
        return reedPosition;
    }
    
    double getTipVelocity() const override {
        return reedVelocity;
    }
    
    // For backward compatibility with your existing code
    double getPosition() const { return reedPosition; }
    double getVelocity() const { return reedVelocity; }
    void setQ(double Q) { params.Q = Q; }
    void setMu(double mu) { params.mu = mu; }
    double getOmega0() const { return params.omega0; }

private:
    void updatePhysicsState() {
        currentState.position = reedPosition;
        currentState.velocity = reedVelocity;
        currentState.frequency = params.omega0 / (2.0 * M_PI);
        currentState.amplitude = std::abs(reedPosition);
    }
};
