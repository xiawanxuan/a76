#pragma once

#include "Types.h"
#include "Mesh2D.h"
#include <cmath>
#include <functional>

namespace RoadbedSim {

class ThermoHydro {
public:
    explicit ThermoHydro(const Mesh2D& mesh);
    ~ThermoHydro() = default;

    Scalar unfrozenWaterContent(Scalar temperature, Scalar totalWaterContent,
                                Scalar porosity) const;
    Scalar iceContent(Scalar temperature, Scalar totalWaterContent,
                      Scalar porosity) const;
    Scalar unfrozenWaterDerivative(Scalar temperature, Scalar totalWaterContent,
                                   Scalar porosity) const;

    Scalar equivalentThermalConductivity(Scalar temperature,
                                         const ZoneProperties& zone) const;
    Scalar equivalentHeatCapacity(Scalar temperature,
                                  const ZoneProperties& zone) const;

    Scalar hydraulicConductivity(Scalar temperature,
                                 const ZoneProperties& zone) const;
    Scalar waterPotential(Scalar temperature, Scalar waterContent,
                          const ZoneProperties& zone) const;

    Scalar frostHeaveStrain(Scalar temperature, Scalar initialWaterContent,
                            Scalar porosity) const;
    Scalar frostHeaveStrainDerivative(Scalar temperature,
                                      Scalar initialWaterContent,
                                      Scalar porosity) const;

    Scalar effectiveStressFactor(Scalar temperature, Scalar waterContent,
                                 Scalar porosity) const;

    void computeElementState(Index elemId, Scalar temperature,
                             Scalar waterContent,
                             Scalar& kThermal, Scalar& cEquiv,
                             Scalar& kHydraulic, Scalar& heaveStrain) const;

    VectorX computeInitialTemperatureField(Scalar surfaceTemp,
                                           Scalar geothermalGradient) const;
    VectorX computeInitialWaterContentField() const;

    Scalar annualAirTemperature(Scalar time, Scalar meanTemp,
                                Scalar amplitude, Scalar periodDays,
                                Scalar phaseShiftDays) const;

    void setFreezingPointDepression(Scalar d) { freezingPointDepression_ = d; }
    Scalar getFreezingPointDepression() const { return freezingPointDepression_; }

    void setVanGenuchtenParams(Scalar alpha, Scalar n, Scalar thetaR, Scalar thetaS) {
        vgAlpha_ = alpha;
        vgN_ = n;
        vgThetaR_ = thetaR;
        vgThetaS_ = thetaS;
    }

private:
    const Mesh2D& mesh_;
    Scalar freezingPointDepression_ = 0.5;
    Scalar vgAlpha_ = 0.01;
    Scalar vgN_ = 2.0;
    Scalar vgThetaR_ = 0.05;
    Scalar vgThetaS_ = 0.45;
};

} // namespace RoadbedSim
