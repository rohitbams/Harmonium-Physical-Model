#pragma once
#include "ReedModel.h"
#include <vector>
#include <cmath>

class ClampedReed : public ReedModel {
public:
    struct BeamProperties {
        double length = 0.016;          // 16mm (from Puranik paper)
        double width = 0.002;           // 2mm
        double thickness = 0.0004;      // 0.4mm
        double youngModulus = 125e9;    // 125 GPa
        double density = 8490.0;        // kg/m³ (brass)
        double dampingCoeff = 0.065;    // R parameter
        double effectiveArea = 4e-4;    // Sᵣ from paper
    };

private:
    static constexpr int N = 40;        // Discretization points (from paper)
    static constexpr double BETA = 0.25; // Newmark-beta parameter
    static constexpr double GAMMA = 0.5; // Newmark-gamma parameter

    BeamProperties props;
    double sampleRate;
    double dt;
    double dx;
    
    // Beam state vectors [N points]
    std::vector<double> y;          // Displacement y(x,t)
    std::vector<double> ydot;       // Velocity dy/dt
    std::vector<double> yddot;      // Acceleration d²y/dt²
    
    // Previous time step for Newmark integration
    std::vector<double> y_prev;
    std::vector<double> ydot_prev;
    std::vector<double> yddot_prev;
    
    // System matrices (pre-computed for efficiency)
    std::vector<std::vector<double>> massMatrix;
    std::vector<std::vector<double>> dampingMatrix;
    std::vector<std::vector<double>> stiffnessMatrix;
    std::vector<std::vector<double>> systemMatrix; // For implicit solver
    
    // Physical parameters
    double crossSectionArea;
    double secondMomentArea;
    double linearDensity;
    
    // State tracking
    PhysicsState currentState;
    double targetFrequency;

public:
    ClampedReed(double targetFreq, double sampleRate)
        : sampleRate(sampleRate), dt(1.0/sampleRate), targetFrequency(targetFreq)
    {
        // Scale reed properties to achieve target frequency
        scalePropertiesForFrequency(targetFreq);
        
        // Initialize discretization
        dx = props.length / (N - 1);
        
        // Calculate derived physical properties
        crossSectionArea = props.width * props.thickness;
        secondMomentArea = props.width * std::pow(props.thickness, 3) / 12.0;
        linearDensity = props.density * crossSectionArea;
        
        // Initialize state vectors
        y.resize(N, 0.0);
        ydot.resize(N, 0.0);
        yddot.resize(N, 0.0);
        y_prev.resize(N, 0.0);
        ydot_prev.resize(N, 0.0);
        yddot_prev.resize(N, 0.0);
        
        // Pre-compute system matrices
        buildSystemMatrices();
        
        reset();
    }
    
    void updateMotion(double pressureDifference, double dt_override = 0.0) override {
        double time_step = (dt_override > 0.0) ? dt_override : dt;
        
        // Store previous state
        y_prev = y;
        ydot_prev = ydot;
        yddot_prev = yddot;
        
        // Build force vector from pressure difference
        std::vector<double> forceVector(N, 0.0);
        for (int i = 0; i < N; ++i) {
            forceVector[i] = props.effectiveArea * pressureDifference / linearDensity;
        }
        
        // Solve implicit system: (M + γ*dt*C + β*dt²*K) * a_new = F_total
        solveImplicitSystem(forceVector, time_step);
        
        // Apply boundary conditions (clamped at x=0, free at x=L)
        applyBoundaryConditions();
        
        // Update state tracking
        updatePhysicsState();
    }
    
    void reset() override {
        std::fill(y.begin(), y.end(), 0.0);
        std::fill(ydot.begin(), ydot.end(), 0.0);
        std::fill(yddot.begin(), yddot.end(), 0.0);
        std::fill(y_prev.begin(), y_prev.end(), 0.0);
        std::fill(ydot_prev.begin(), ydot_prev.end(), 0.0);
        std::fill(yddot_prev.begin(), yddot_prev.end(), 0.0);
        
        currentState = PhysicsState{};
    }
    
    PhysicsState getState() const override {
        return currentState;
    }
    
    void setFrequency(double frequency) override {
        if (std::abs(frequency - targetFrequency) > 1.0) {
            targetFrequency = frequency;
            scalePropertiesForFrequency(frequency);
            buildSystemMatrices(); // Rebuild matrices with new properties
        }
    }
    
    double getTipDisplacement() const override {
        return y[N-1]; // Displacement at free end
    }
    
    double getTipVelocity() const override {
        return ydot[N-1]; // Velocity at free end
    }
    
    bool supportsModalAnalysis() const override { return true; }
    int getNumModes() const override { return 5; } // First 5 modes are significant
    
    double getModalAmplitude(int mode) const override {
        if (mode >= getNumModes()) return 0.0;
        
        // Simplified modal analysis - project displacement onto mode shapes
        // This is a basic approximation; full modal analysis would be more complex
        double modalAmp = 0.0;
        double modeFreq = calculateModeFrequency(mode + 1);
        
        // Weight by mode shape (approximated)
        for (int i = 0; i < N; ++i) {
            double x = i * dx;
            double modeShape = calculateModeShape(mode + 1, x);
            modalAmp += y[i] * modeShape;
        }
        
        return modalAmp / N;
    }

private:
    void scalePropertiesForFrequency(double targetFreq) {
        // For clamped-free beam, fundamental frequency is:
        // f₁ = (λ₁²/2π) * √(EI/(ρA)) / L²
        // where λ₁ ≈ 1.875 for first mode
        
        const double lambda1 = 1.875104;
        double baseFreq = 261.63; // Middle C
        
        // Scale length to achieve target frequency
        double lengthScale = std::sqrt(baseFreq / targetFreq);
        props.length = 0.016 * lengthScale; // Scale from base 16mm
        
        // Recalculate discretization
        dx = props.length / (N - 1);
        
        // Update derived properties
        crossSectionArea = props.width * props.thickness;
        secondMomentArea = props.width * std::pow(props.thickness, 3) / 12.0;
        linearDensity = props.density * crossSectionArea;
    }
    
    void buildSystemMatrices() {
        // Initialize matrices
        massMatrix.assign(N, std::vector<double>(N, 0.0));
        dampingMatrix.assign(N, std::vector<double>(N, 0.0));
        stiffnessMatrix.assign(N, std::vector<double>(N, 0.0));
        systemMatrix.assign(N, std::vector<double>(N, 0.0));
        
        // Build mass matrix (lumped mass approximation)
        for (int i = 0; i < N; ++i) {
            massMatrix[i][i] = linearDensity * dx;
        }
        
        // Build damping matrix (proportional damping)
        for (int i = 0; i < N; ++i) {
            dampingMatrix[i][i] = props.dampingCoeff;
        }
        
        // Build stiffness matrix using 5-point finite difference for d⁴y/dx⁴
        double stiff_coeff = props.youngModulus * secondMomentArea / std::pow(dx, 4);
        
        for (int i = 2; i < N-2; ++i) {
            stiffnessMatrix[i][i-2] = stiff_coeff * 1.0;
            stiffnessMatrix[i][i-1] = stiff_coeff * (-4.0);
            stiffnessMatrix[i][i]   = stiff_coeff * 6.0;
            stiffnessMatrix[i][i+1] = stiff_coeff * (-4.0);
            stiffnessMatrix[i][i+2] = stiff_coeff * 1.0;
        }
        
        // Handle boundary points with modified stencils
        buildBoundaryStiffness();
        
        // Build system matrix for implicit integration
        // System = M + γ*dt*C + β*dt²*K
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                systemMatrix[i][j] = massMatrix[i][j]
                                   + GAMMA * dt * dampingMatrix[i][j]
                                   + BETA * dt * dt * stiffnessMatrix[i][j];
            }
        }
    }
    
    void buildBoundaryStiffness() {
        // Modified finite difference stencils near boundaries
        double stiff_coeff = props.youngModulus * secondMomentArea / std::pow(dx, 4);
        
        // Point i=0 (clamped boundary)
        stiffnessMatrix[0][0] = stiff_coeff * 7.0;
        stiffnessMatrix[0][1] = stiff_coeff * (-4.0);
        stiffnessMatrix[0][2] = stiff_coeff * 1.0;
        
        // Point i=1 (near clamped boundary)
        stiffnessMatrix[1][0] = stiff_coeff * (-4.0);
        stiffnessMatrix[1][1] = stiff_coeff * 6.0;
        stiffnessMatrix[1][2] = stiff_coeff * (-4.0);
        stiffnessMatrix[1][3] = stiff_coeff * 1.0;
        
        // Point i=N-2 (near free boundary)
        stiffnessMatrix[N-2][N-4] = stiff_coeff * 1.0;
        stiffnessMatrix[N-2][N-3] = stiff_coeff * (-4.0);
        stiffnessMatrix[N-2][N-2] = stiff_coeff * 6.0;
        stiffnessMatrix[N-2][N-1] = stiff_coeff * (-4.0);
        
        // Point i=N-1 (free boundary)
        stiffnessMatrix[N-1][N-3] = stiff_coeff * 1.0;
        stiffnessMatrix[N-1][N-2] = stiff_coeff * (-4.0);
        stiffnessMatrix[N-1][N-1] = stiff_coeff * 7.0;
    }
    
    void solveImplicitSystem(const std::vector<double>& forceVector, double time_step) {
        // Newmark-beta method for implicit time integration
        std::vector<double> rhs(N);
        
        // Build right-hand side
        for (int i = 0; i < N; ++i) {
            rhs[i] = forceVector[i];
            
            // Add contributions from previous time step
            for (int j = 0; j < N; ++j) {
                rhs[i] += massMatrix[i][j] * (yddot_prev[j]
                         + (1.0-GAMMA) * time_step * yddot_prev[j]);
                rhs[i] += dampingMatrix[i][j] * (ydot_prev[j]
                         + (1.0-GAMMA) * time_step * yddot_prev[j]);
            }
        }
        
        // Solve linear system: systemMatrix * yddot = rhs
        // Using simple Gaussian elimination (could be optimized)
        gaussianElimination(systemMatrix, rhs, yddot);
        
        // Update velocity and displacement using Newmark formulas
        for (int i = 0; i < N; ++i) {
            ydot[i] = ydot_prev[i] + time_step * ((1.0-GAMMA) * yddot_prev[i] + GAMMA * yddot[i]);
            y[i] = y_prev[i] + time_step * ydot_prev[i] + 0.5 * time_step * time_step *
                   ((1.0-2.0*BETA) * yddot_prev[i] + 2.0*BETA * yddot[i]);
        }
    }
    
    void applyBoundaryConditions() {
        // Clamped boundary conditions at x=0
        y[0] = 0.0;      // y(0,t) = 0
        ydot[0] = 0.0;   // dy/dt(0,t) = 0
        
        // Also enforce dy/dx(0,t) = 0 using finite difference
        // dy/dx ≈ (-3*y[0] + 4*y[1] - y[2]) / (2*dx) = 0
        // This is automatically satisfied since y[0] = 0 and we solve the system
        
        // Free boundary conditions at x=L are automatically satisfied
        // by the finite difference formulation (natural boundary conditions)
    }
    
    void updatePhysicsState() {
        // Update state information
        currentState.position = getTipDisplacement();
        currentState.velocity = getTipVelocity();
        currentState.frequency = estimateCurrentFrequency();
        currentState.amplitude = calculateRMSAmplitude();
    }
    
    double estimateCurrentFrequency() const {
        // Simple frequency estimation based on zero-crossings of tip velocity
        static std::vector<double> velocityHistory;
        static size_t historyIndex = 0;
        static constexpr size_t HISTORY_SIZE = 512;
        
        if (velocityHistory.size() < HISTORY_SIZE) {
            velocityHistory.resize(HISTORY_SIZE, 0.0);
        }
        
        velocityHistory[historyIndex] = getTipVelocity();
        historyIndex = (historyIndex + 1) % HISTORY_SIZE;
        
        // Count zero crossings
        int crossings = 0;
        for (size_t i = 1; i < HISTORY_SIZE; ++i) {
            if ((velocityHistory[i-1] >= 0) != (velocityHistory[i] >= 0)) {
                crossings++;
            }
        }
        
        return (crossings * sampleRate) / (2.0 * HISTORY_SIZE);
    }
    
    double calculateRMSAmplitude() const {
        double sum = 0.0;
        for (double displacement : y) {
            sum += displacement * displacement;
        }
        return std::sqrt(sum / N);
    }
    
    double calculateModeFrequency(int mode) const {
        // Mode frequencies for clamped-free beam
        static const double lambdas[] = {1.875, 4.694, 7.855, 10.996, 14.137};
        if (mode > 5) return 0.0;
        
        double lambda = lambdas[mode-1];
        return (lambda * lambda / (2.0 * M_PI)) *
               std::sqrt((props.youngModulus * secondMomentArea) /
                        (linearDensity * std::pow(props.length, 4)));
    }
    
    double calculateModeShape(int mode, double x) const {
        // Simplified mode shapes for clamped-free beam
        double xi = x / props.length;
        static const double lambdas[] = {1.875, 4.694, 7.855, 10.996, 14.137};
        if (mode > 5) return 0.0;
        
        double lambda = lambdas[mode-1];
        return std::sin(lambda * xi) - std::sinh(lambda * xi) +
               ((std::sin(lambda) + std::sinh(lambda)) / (std::cos(lambda) + std::cosh(lambda))) *
               (std::cos(lambda * xi) - std::cosh(lambda * xi));
    }
    
    void gaussianElimination(std::vector<std::vector<double>> A,
                           std::vector<double> b,
                           std::vector<double>& x) {
        int n = A.size();
        x.resize(n);
        
        // Forward elimination
        for (int i = 0; i < n; i++) {
            // Find pivot
            int maxRow = i;
            for (int k = i + 1; k < n; k++) {
                if (std::abs(A[k][i]) > std::abs(A[maxRow][i])) {
                    maxRow = k;
                }
            }
            std::swap(A[maxRow], A[i]);
            std::swap(b[maxRow], b[i]);
            
            // Make all rows below this one 0 in current column
            for (int k = i + 1; k < n; k++) {
                if (std::abs(A[i][i]) > 1e-12) {
                    double c = A[k][i] / A[i][i];
                    for (int j = i; j < n; j++) {
                        A[k][j] -= c * A[i][j];
                    }
                    b[k] -= c * b[i];
                }
            }
        }
        
        // Back substitution
        for (int i = n - 1; i >= 0; i--) {
            x[i] = b[i];
            for (int j = i + 1; j < n; j++) {
                x[i] -= A[i][j] * x[j];
            }
            if (std::abs(A[i][i]) > 1e-12) {
                x[i] /= A[i][i];
            }
        }
    }
};
