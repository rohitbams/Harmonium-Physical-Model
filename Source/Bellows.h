#pragma once
/*
 * The Bellows class.
 * This private class handles the novel bellowing mechanism.
 * It calculates the air inflow and the pressure variation in the bellows chamber.
 * It also handels modulation wheen interaction on the MIDI keyboard.
 * The Bellows mode can be turned off to simulate a constant inifinite air supply.
 */
class Bellows {
private:
    HarmoniumState state;
    double sampleRate;
    double dt;
    double airMass = 0.003;             // initial air in bellows
    double maxAirMass = 0.02;           // maximum capacity
    double reedChamberMass = 0;         // initial air in Reed chamber
    double maxReedChamberMass = 0;      // maximum capacity in reed chamber
    double narrowJetMass = 0;           // initial air in narrow jet
    double maxNarrowJetMass = 0;        // maximum narrow jet capacity
    double dampingForce = 0;            // damping force
    double p0 = state.p0;
    double u0 = state.u0;
    double modWheelValue = 0.0;
    bool keyPressed = false;
    double baseConsumptionRate = 0.01;
    double pumpAmount = 0.05;           // air mass added per pumping action
    double previousModWheelValue = 0.0;
    double minMovementThreshold = 0.01; // minimum movement to trigger pump
    double continuousPumpingRate = 0.0; // current pumping rate (kgs/s)
    double pumpingRate = 0.8;
    
public:
    Bellows(double sampleRate) : sampleRate(sampleRate), dt(1.0/sampleRate) {}
    

    /*
     * this method handles bellowing (pumping mechanism) and the air pressure p0 buildup
     * in the bellows chameber V0, and air flow from the bellows chamber into the
     * reed chamber V2.
     */
    void updateBellows(double chamberPressure, int oscillatingReeds) {
        
        // pumping mechanism and air mass management
        double modWheelMovement = modWheelValue - previousModWheelValue;
        double movementSpeed = std::abs(modWheelMovement) / dt;
        continuousPumpingRate = std::min(movementSpeed * 10.0, 20.0);
        
        // add air mass while pumping
        if (continuousPumpingRate > 0.001) {
            airMass += continuousPumpingRate * dt;
        }
        
        // reduce air mass while pressing note keys
        if (keyPressed && oscillatingReeds > 0) {
            double consumptionRate = baseConsumptionRate * oscillatingReeds * dt;
            airMass -= consumptionRate;
        }
        
        airMass = std::clamp(airMass, 0.0, maxAirMass); // limit air mass to avoid negative value division
        
        double airRatio = airMass / maxAirMass;
        p0 = airRatio * 300.0;

        double pressureDiff = p0 - chamberPressure;  // p0 - p1
        
        // calculate flow when there's significant pressure difference
        if (std::abs(pressureDiff) > 2.0) {
            
            double flowDirection = (pressureDiff > 0) ? 1.0 : -1.0;
            double velocity = std::sqrt(2.0 * std::abs(pressureDiff) / 1.225); // 1.255 is air density (rho0)
            
            if (keyPressed) {
                // keys pressed, full connection (valve fully open)
                u0 = flowDirection * 0.001 * velocity;
            } else {
                // keys not pressed, restricted connection
                u0 = flowDirection * 0.0003 * velocity;  // 30% of full flow
            }
            
            u0 = std::clamp(u0, -0.05, 0.05);
        } else {
            // pressures are equal, no flow
            u0 = 0.0;
        }

        
        previousModWheelValue = modWheelValue;
    }
    
    void setModWheelValue(double value) {
        modWheelValue = std::clamp(value, 0.0, 1.0);
    }
    
    void setKeyPressed(bool pressed) {
        keyPressed = pressed;
    }
    
    double getBellowsPressure() const { return p0; }
    double getFlowRate() const { return u0; }
    double getAirMass() const { return airMass; }
    double getRelativePressure() const { return std::min(1.0, p0 / 1000.0); }
    
    void setMaxAirMass(double mass) {
        maxAirMass = std::clamp(mass, 0.1, 0.4);
        if (airMass > maxAirMass) {
            airMass = maxAirMass;
        }
    }

    void setMaxReedChamberAirMass(double mass) {
        maxReedChamberMass = std::clamp(mass, 0.1, 0.2);
        if (reedChamberMass > maxReedChamberMass) reedChamberMass = maxReedChamberMass * pumpingRate;
    }
    
    void setMaxNarrowJetAirMass(double mass) {
        narrowJetMass = std::clamp(mass, 0.1, 0.2);
        if (narrowJetMass > maxNarrowJetMass) narrowJetMass = maxNarrowJetMass * pumpingRate;
    }
    
    void setBaseConsumptionRate(double rate) {
        baseConsumptionRate = std::clamp(rate, 0.001, 0.02);
    }
};
