// PhysicsEngine.h
#pragma once
#include <JuceHeader.h>
#include <cmath>
#include <algorithm>
#include "HarmoniumState.h"
#include "PhysicsEngine.h"
#include "Voice.h"
#include "Harmonium.h"
/*
 * The PhysicsEngine class.
 * This class contains all physics values
 */
class PhysicsEngine {
public:
    
    
    /*
     * This method calls and implements all main physics equations
     *
     */
    static HarmoniumState updatePhysics(
        const HarmoniumState& currentState,
        const HarmoniumPhysicsConfig& config,
        double bellowsFlowRate,
        double bellowsPressure,
        double omega0
    ) {
        HarmoniumState newState = currentState;
        
        // calculate pressure dependent damping
        calculatePressureDependentDamping(newState, config, omega0);
        
        // calculate forces (equation 5)
        calculateReedForces(newState, config, omega0);
        
        // apply oscillation startup
        handleOscillationStartup(newState, config);
        
        // apply energy input for sustained oscillation
        applyEnergyInput(newState, config);
        
        // integrate reed motion (Equation 4)
        integrateReedMotion(newState, config);
        
        // handle low amplitude restart
        handleLowAmplitudeRestart(newState, config);
        
        // calculate volume flow (Equation 5)
        calculateVolumeFlow(newState, config);
        
        // update reed chamber pressure (Equation 1)
        updateChamberPressure(newState, config, bellowsFlowRate);
        
        // update jet dynamics (Equations 2 & 3)
        updateJetDynamics(newState, config);
        
        // apply physical limits
        applyPhysicalLimits(newState);
        
        // update flow history
        newState.u_previous = currentState.u;
        newState.p0 = bellowsPressure;
        
        return newState;
    }
    
    /*
     * This method adds harmonic richness to the main outputs and prepares audio samples to be
     */
    static float generateAudio(const HarmoniumState& state, const HarmoniumPhysicsConfig& config) {
        double totalFlowFundamental = state.reedPosition * config.amplitudeScale;
        double fundamental = state.reedPosition;
//        double reedPosFundamental = state.reedPosition * amplitudeScale;
        
        /*
         Sub-bass = 15-50 Hz
         Bass = 50-200 Hz
         Mid-bass = 200-400Hz
         Lower mids 400-900Hz
         Upper mids = 900-4.000Hz
         Treble – 4.000-20.000Hz
         */
            
        double low = config.low;           // 15 - 200
        double lowMid = config.lowMid;     // 200 - 400
        double mid = config.mid;           // 400 - 900
        double highMid = config.highMid;   // 900 - 4000
        double high = config.high;         // 4000 - 20000
        
        // -- rich harmonic profile based on hinge og recording -- //
        double harmonic2 = fundamental * 2 * (config.amplitudeScale * 0.96 /   low);     // 698.4 / 2
        double harmonic3 = fundamental * 5 * (config.amplitudeScale * 0.63 /   highMid); // 1746.0 / 5
        double harmonic4 = fundamental * 4 * (config.amplitudeScale * 0.60 /   highMid); // 1396.8 / 4
        double harmonic5 = fundamental * 3 * (config.amplitudeScale * 0.57 /   highMid); // 1047.6 / 3
        double harmonic6 = fundamental * 9 * (config.amplitudeScale * 0.34 /   highMid); // 3142.8 / 9
        double harmonic7 = fundamental * 6 * (config.amplitudeScale * 0.32 /   highMid); // 2095.2 / 6
        double harmonic8 = fundamental * 11 * (config.amplitudeScale * 0.32 /  highMid); // 3841.4 / 11
        double harmonic9 = fundamental * 13 * (config.amplitudeScale * 0.21 /  high);    // 4539.8 / 13
        double harmonic10 = fundamental * 17 * (config.amplitudeScale * 0.20 / high);    // 5936.6 / 17
        double harmonic11 = fundamental * 15 * (config.amplitudeScale * 0.15 / high);    // 5238.2 / 15
        double harmonic12 = fundamental * 19 * (config.amplitudeScale * 0.18 / high);    // 6635.0 / 19
        double harmonic13 = fundamental * 7 * (config.amplitudeScale * 0.14 /  high);    // 2444.4 / 7
        double harmonic14 = fundamental * 12 * (config.amplitudeScale * 0.10 / high);    // 4190.6 / 12
        double harmonic15 = fundamental * 21 * (config.amplitudeScale * 0.08 / high);    // 7333.3 / 21
        
        
        double richHarmonic = totalFlowFundamental
                            + harmonic2
                            + harmonic3
                            + harmonic4
                            + harmonic5
                            + harmonic6
                            + harmonic7
                            + harmonic8
                            + harmonic9
                            + harmonic10
                            + harmonic11
                            + harmonic12
                            + harmonic13
                            + harmonic14
                            + harmonic15;

        return static_cast<float>(richHarmonic);
        //        richHarmonic = std::tanh(richHarmonic * 2.5);
        
        // -- mix wave with saturation -- //
//        double mixWave = (totalFlowFundamental) + (reedPosFundamental) + (richHarmonic);
//        mixWave = std::tanh(mixWave * 3) * 0.5;
//        mixWave = std::tanh((std::tanh(mixWave * 3) * 0.5) * 1.5) * 0.7;
        
//        const HarmoniumState& saturationFactor;
        
//        totalFlowFundamental = std::tanh(totalFlowFundamental * 0.3);
//        return static_cast<float>(totalFlowFundamental);
//        return static_cast<float>(reedPosFundamental);
        
        
        
//        double finalOutput = std::tanh(finalOutput * 1.5) * 0.7;
//        return static_cast<float>(finalOutput);
        

        
    }
    
    
    static double calculateOmega0(double frequency) {
        return 2.0 * M_PI * frequency;
    }

     // TODO: calculate detuned frequencies for realistic variation
    static double calculateDetunedFrequency(double baseFreq, double cents) {
        return baseFreq * std::pow(2.0, cents / 1200.0);
    }
    
    static double getAperture(const HarmoniumState& state, const HarmoniumPhysicsConfig& config) {
        double aperture = config.reedWidth * std::abs(state.reedPosition);
        double minAperture = config.reedWidth * config.restGap;
        return std::max(aperture, minAperture);
    }

private:
    // pressure-dependent damping
    static void calculatePressureDependentDamping(HarmoniumState& state, const HarmoniumPhysicsConfig& config, double omega0) {
        double effectiveQ = config.Q;
        // TODO: WHY ARE WE USING airPressureRatio
        double airPressureRatio = std::min(1.0, std::abs(state.p2) / 500.0);

        if (airPressureRatio > 0.2) {
            effectiveQ = config.Q * 2.0;  // high sustain with air pressure
        } else {
            effectiveQ = config.Q * (0.2 + 0.8 * airPressureRatio);  // natural decay
        }
        
        state.dampingForce = (omega0 / effectiveQ) * state.reedVelocity;
    }
    
    // reed forces calculation (equation 4)
    static void calculateReedForces(HarmoniumState& state, const HarmoniumPhysicsConfig& config, double omega0) {
        state.drivingForce = config.mu * state.p2;
        state.springForce = omega0 * omega0 * state.reedPosition;
        state.acceleration = state.drivingForce - state.springForce - state.dampingForce;
    }
    
    // oscillation startup
    static void handleOscillationStartup(HarmoniumState& state, const HarmoniumPhysicsConfig& config) {
        if (std::abs(state.p2) > 50.0) {
            double minAmplitude = 1e-5;
            if (std::abs(state.reedPosition) < minAmplitude && std::abs(state.reedVelocity) < 0.001) {
//                double phaseDirection = (state.reedPosition >= 0) ? 1.0 : -1.0;
                state.reedPosition = minAmplitude/* * phaseDirection*/;
                state.reedVelocity = 0.1/* * phaseDirection*/;
            }
        }
    }
    
    // energy input mechanism for sustained oscillation
    static void applyEnergyInput(HarmoniumState& state, const HarmoniumPhysicsConfig& config) {
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
    
    // numerical integration Euler method
    static void integrateReedMotion(HarmoniumState& state, const HarmoniumPhysicsConfig& config) {
        state.reedVelocity += state.acceleration * config.dt;
        state.reedPosition += state.reedVelocity * config.dt;
    }
    
    // low amplitude restart mechanism
    static void handleLowAmplitudeRestart(HarmoniumState& state, const HarmoniumPhysicsConfig& config) {
        if (std::abs(state.reedPosition) < config.minOscillationLevel &&
            std::abs(state.reedVelocity) < config.minOscillationLevel) {
            
            if (std::abs(state.p2) > 100.0) {
//                double pressureRatio = std::clamp(std::abs(state.p2) / 1000.0, 0.01, 0.05);
                
                state.reedPosition += (state.p2 > 0 ? 1.0 : -1.0) * config.minOscillationLevel * 10.0 /** pressureRatio*/;
                state.reedVelocity += config.mu * state.p2 * config.dt /** pressureRatio*/;
            }
        }
    }
    
    // volume flow calculation (Equation 5 total volume flow)
    static void calculateVolumeFlow(HarmoniumState& state, const HarmoniumPhysicsConfig& config) {
        // equation 5: u = Sr·(dζ/dt) + αSu·vj
        double reedPumpingFlow = config.Sr * state.reedVelocity;
        state.aperture = getAperture(state, config);
        double jetFlow = config.alpha * state.aperture * state.vj;
        state.u = reedPumpingFlow + jetFlow;
        state.u = std::clamp(state.u, -0.1, 0.1);
    }
    
    // chamber pressure update (equation 1 mass conservation)
    static void updateChamberPressure(HarmoniumState& state, const HarmoniumPhysicsConfig& config, double u0) {
        // V₁/c₀² × d(p₁-p_atm)/dt = ρ₀(u₀ - u)
        state.netFlow = u0 - state.u;
        double denominator = config.V1 / (config.c0 * config.c0);
        double dp1_dt = (config.rho0 * state.netFlow) / denominator;

//        state.p1 += dp1_dt * config.dt;

        double p0 = state.p0;

        if (p0 < 3500) {
            state.p1 += dp1_dt * config.dt;
        } else {
            state.p1 += dp1_dt * config.dt * 5;
        }
        
        // natural pressure decay when no net flow
        if (std::abs(state.netFlow) < 1e-6) {
            state.p1 *= 0.995;
            state.p1 = std::max(0.0, state.p1);
        }
    }
    
    // Jet dynamics update (Equations 2 and 3)
    static void updateJetDynamics(HarmoniumState& state, const HarmoniumPhysicsConfig& config) {
    /// Equation 2:   p₁  =  p₂  +  ρ₀L₂/S₂  ×  du/dt
        double du_dt = (state.u - state.u_previous) / config.dt;
        double inertialDrop = (config.rho0 * config.L2 / config.S2) * du_dt;
        state.p2 = state.p1 - inertialDrop;
        
    /// Equation 3: p₂ = p_atm + ½ρ₀vⱼ²
        double pressureDiff = state.p2 - config.p_atm;
        if (pressureDiff > 0.0) {
            state.vj = std::sqrt(2.0 * pressureDiff / config.rho0);
        } else {
            state.vj = 0.0;
        }
    }
    
    // apply physical limits for numerical stability
    static void applyPhysicalLimits(HarmoniumState& state) {
        state.reedPosition = std::clamp(state.reedPosition, -0.01, 0.01);
        state.reedVelocity = std::clamp(state.reedVelocity, -20.0, 20.0);
        state.p1 = std::clamp(state.p1, 0.0, 2000.0);
        state.p2 = std::clamp(state.p2, -100.0, 1500.0);
        state.vj = std::clamp(state.vj, 0.0, 100.0);
        state.u = std::clamp(state.u, -0.1, 0.1);
//        state.saturationFactor = std::clamp(state.saturationFactor, 0.1, 5.0);
    }
};
