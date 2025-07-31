// HarmoniumPhysicsEngine.h
#pragma once
#include <JuceHeader.h>
#include <cmath>
#include <algorithm>
#include "HarmoniumPhysicsState.h"

// Pure computational engine - all physics in one place
class HarmoniumPhysicsEngine {
public:
    // Main physics update - implements all 5 Puranik-Scavone equations
    static HarmoniumPhysicsState updatePhysics(
        const HarmoniumPhysicsState& currentState,
        const HarmoniumPhysicsConfig& config,
        double bellowsFlowRate,      // u₀ - flow from bellows
        double bellowsPressure,      // p₀ - bellows pressure
        double omega0                // ω₀ - reed natural frequency
    ) {
        HarmoniumPhysicsState newState = currentState;
        
        // Step 1: Calculate pressure-dependent damping
        calculatePressureDependentDamping(newState, config, omega0);
        
        // Step 2: Calculate reed forces (Equation 4 setup)
        calculateReedForces(newState, config, omega0);
        
        // Step 3: Handle oscillation startup
        handleOscillationStartup(newState, config);
        
        // Step 4: Apply energy input for sustained oscillation
        applyEnergyInput(newState, config);
        
        // Step 5: Integrate reed motion (Equation 4)
        integrateReedMotion(newState, config);
        
        // Step 6: Handle low amplitude restart
        handleLowAmplitudeRestart(newState, config);
        
        // Step 7: Calculate volume flow (Equation 5)
        calculateVolumeFlow(newState, config);
        
        // Step 8: Update chamber pressure (Equation 1)
        updateChamberPressure(newState, config, bellowsFlowRate);
        
        // Step 9: Update jet dynamics (Equations 2 & 3)
        updateJetDynamics(newState, config);
        
        // Step 10: Apply physical limits
        applyPhysicalLimits(newState);
        
        // Update flow history
        newState.u_previous = currentState.u;
        
        return newState;
    }
    
    // Audio generation
    static float generateAudio(const HarmoniumPhysicsState& state, int amplitudeScale = 1000) {
        double totalFlowFundamental = state.u * amplitudeScale;
        double reedPosFundamental = state.reedPosition * amplitudeScale;
        double richHarmonic = totalFlowFundamental + reedPosFundamental;
        richHarmonic = std::tanh(richHarmonic * 2.5);
        
        double mixWave = (totalFlowFundamental) + (reedPosFundamental);
        mixWave = std::tanh(mixWave * 3) * 0.5;
        return static_cast<float>(mixWave);
    }
    
    // Utility functions
    static double calculateOmega0(double frequency) {
        return 2.0 * M_PI * frequency;
    }
    
    static double calculateDetunedFrequency(double baseFreq, double cents) {
        return baseFreq * std::pow(2.0, cents / 1200.0);
    }
    
    static double getAperture(const HarmoniumPhysicsState& state, const HarmoniumPhysicsConfig& config) {
        double aperture = config.reedWidth * std::abs(state.reedPosition);
        double minAperture = config.reedWidth * config.rest_gap;
        return std::max(aperture, minAperture);
    }

private:
    // Pressure-dependent damping (key enhancement for realism)
    static void calculatePressureDependentDamping(HarmoniumPhysicsState& state, const HarmoniumPhysicsConfig& config, double omega0) {
        double effectiveQ = config.Q;
        double airPressureRatio = std::min(1.0, std::abs(state.p2) / 500.0);

        if (airPressureRatio > 0.2) {
            effectiveQ = config.Q * 2.0;  // High sustain with air pressure
        } else {
            effectiveQ = config.Q * (0.2 + 0.8 * airPressureRatio);  // Natural decay
        }
        
        state.dampingForce = (omega0 / effectiveQ) * state.reedVelocity;
    }
    
    // Reed forces calculation (Equation 4 components)
    static void calculateReedForces(HarmoniumPhysicsState& state, const HarmoniumPhysicsConfig& config, double omega0) {
        state.drivingForce = config.mu * state.p2;
        state.springForce = omega0 * omega0 * state.reedPosition;
        state.acceleration = state.drivingForce - state.springForce - state.dampingForce;
    }
    
    // Oscillation startup enhancement
    static void handleOscillationStartup(HarmoniumPhysicsState& state, const HarmoniumPhysicsConfig& config) {
        if (std::abs(state.p2) > 100.0) {
            double minAmplitude = 1e-5;
            if (std::abs(state.reedPosition) < minAmplitude && std::abs(state.reedVelocity) < 0.01) {
                double phaseDirection = (state.reedPosition >= 0) ? 1.0 : -1.0;
                state.reedPosition = minAmplitude * phaseDirection;
                state.reedVelocity = 0.1 * phaseDirection;
            }
        }
    }
    
    // Energy input mechanism (key breakthrough for sustained oscillation)
    static void applyEnergyInput(HarmoniumPhysicsState& state, const HarmoniumPhysicsConfig& config) {
        double currentAmplitude = std::abs(state.reedPosition);
        double amplitudeRatio = std::min(1.0, currentAmplitude / config.targetAmplitude);
        double energyReduction = 1.0 - amplitudeRatio;
        double energyInputStrength = 10.0 * energyReduction;
        double airflowEnergyInput = config.mu * state.p2 * config.dt * energyInputStrength;

        if (currentAmplitude < config.targetAmplitude && std::abs(state.reedVelocity) > 0.001) {
            double velocityDirection = (state.reedVelocity > 0) ? 1.0 : -1.0;
            state.reedVelocity += airflowEnergyInput * velocityDirection * 0.1;
        }
    }
    
    // Numerical integration (Euler method)
    static void integrateReedMotion(HarmoniumPhysicsState& state, const HarmoniumPhysicsConfig& config) {
        state.reedVelocity += state.acceleration * config.dt;
        state.reedPosition += state.reedVelocity * config.dt;
    }
    
    // Low amplitude restart mechanism
    static void handleLowAmplitudeRestart(HarmoniumPhysicsState& state, const HarmoniumPhysicsConfig& config) {
        if (std::abs(state.reedPosition) < config.minOscillationLevel &&
            std::abs(state.reedVelocity) < config.minOscillationLevel) {
            
            if (std::abs(state.p2) > 50.0) {
                double pressureRatio = std::clamp(std::abs(state.p2) / 1000.0, 0.01, 0.05);
                
                state.reedPosition += (state.p2 > 0 ? 1.0 : -1.0) * config.minOscillationLevel * 10.0 * pressureRatio;
                state.reedVelocity += config.mu * state.p2 * config.dt * pressureRatio;
            }
        }
    }
    
    // Volume flow calculation (Equation 5)
    static void calculateVolumeFlow(HarmoniumPhysicsState& state, const HarmoniumPhysicsConfig& config) {
        // Equation 5: u = Sr·(dζ/dt) + αSu·vj
        double reedPumpingFlow = config.Sr * state.reedVelocity;
        state.aperture = getAperture(state, config);
        double jetFlow = config.alpha * state.aperture * state.vj;
        state.u = reedPumpingFlow + jetFlow;
        state.u = std::clamp(state.u, -0.1, 0.1);
    }
    
    // Chamber pressure update (Equation 1 - mass conservation)
    static void updateChamberPressure(HarmoniumPhysicsState& state, const HarmoniumPhysicsConfig& config, double u0) {
        // V₁/c₀² × d(p₁-p_atm)/dt = ρ₀(u₀ - u)
        state.netFlow = u0 - state.u;
        double denominator = config.V1 / (config.c0 * config.c0);
        double dp1_dt = (config.rho0 * state.netFlow) / denominator;
        
        state.p1 += dp1_dt * config.dt;
        
        // Natural pressure decay when no net flow
        if (std::abs(state.netFlow) < 1e-6) {
            state.p1 *= 0.995;
            state.p1 = std::max(0.0, state.p1);
        }
    }
    
    // Jet dynamics update (Equations 2 & 3)
    static void updateJetDynamics(HarmoniumPhysicsState& state, const HarmoniumPhysicsConfig& config) {
        // Equation 2: p₁ = p₂ + ρ₀L₂/S₂ × du/dt
        double du_dt = (state.u - state.u_previous) / config.dt;
        double inertialDrop = (config.rho0 * config.L2 / config.S2) * du_dt;
        state.p2 = state.p1 - inertialDrop;
        
        // Equation 3: p₂ = p_atm + ½ρ₀vⱼ²
        double pressureDiff = state.p2 - config.p_atm;
        if (pressureDiff > 0.0) {
            state.vj = std::sqrt(2.0 * pressureDiff / config.rho0);
        } else {
            state.vj = 0.0;
        }
    }
    
    // Physical limits (essential for numerical stability)
    static void applyPhysicalLimits(HarmoniumPhysicsState& state) {
        state.reedPosition = std::clamp(state.reedPosition, -0.01, 0.01);
        state.reedVelocity = std::clamp(state.reedVelocity, -20.0, 20.0);
        state.p1 = std::clamp(state.p1, 0.0, 2000.0);
        state.p2 = std::clamp(state.p2, -100.0, 1500.0);
        state.vj = std::clamp(state.vj, 0.0, 100.0);
        state.u = std::clamp(state.u, -0.1, 0.1);
    }
};
