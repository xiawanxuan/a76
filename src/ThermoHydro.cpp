#include "ThermoHydro.h"

namespace RoadbedSim {

ThermoHydro::ThermoHydro(const Mesh2D& mesh) : mesh_(mesh) {}

Scalar ThermoHydro::unfrozenWaterContent(Scalar temperature, Scalar totalWaterContent,
                                          Scalar porosity) const {
    if (temperature >= FREEZING_POINT - freezingPointDepression_) {
        return totalWaterContent;
    }
    Scalar deltaT = FREEZING_POINT - freezingPointDepression_ - temperature;
    if (deltaT <= 0) return totalWaterContent;
    Scalar maxIce = totalWaterContent * 0.85;
    Scalar unfrozen = totalWaterContent * std::exp(-0.5 * deltaT);
    return std::max(totalWaterContent - maxIce, unfrozen);
}

Scalar ThermoHydro::iceContent(Scalar temperature, Scalar totalWaterContent,
                                Scalar porosity) const {
    Scalar unfrozen = unfrozenWaterContent(temperature, totalWaterContent, porosity);
    return (totalWaterContent - unfrozen) * (WATER_DENSITY / ICE_DENSITY);
}

Scalar ThermoHydro::unfrozenWaterDerivative(Scalar temperature, Scalar totalWaterContent,
                                             Scalar porosity) const {
    if (temperature >= FREEZING_POINT - freezingPointDepression_) {
        return 0.0;
    }
    Scalar deltaT = FREEZING_POINT - freezingPointDepression_ - temperature;
    if (deltaT <= 0) return 0.0;
    Scalar maxIce = totalWaterContent * 0.85;
    Scalar unfrozen = totalWaterContent * std::exp(-0.5 * deltaT);
    Scalar lowerBound = totalWaterContent - maxIce;
    if (unfrozen <= lowerBound) return 0.0;
    return -0.5 * totalWaterContent * std::exp(-0.5 * deltaT);
}

Scalar ThermoHydro::equivalentThermalConductivity(Scalar temperature,
                                                   const ZoneProperties& zone) const {
    Scalar Tref = FREEZING_POINT - freezingPointDepression_;
    if (temperature >= Tref) {
        return zone.thermalCondUnfrozen;
    }
    Scalar deltaT = Tref - temperature;
    Scalar smoothFactor = std::tanh(deltaT / 2.0);
    Scalar kMix = 0.5 * (zone.thermalCondFrozen + zone.thermalCondUnfrozen) +
                  0.5 * (zone.thermalCondFrozen - zone.thermalCondUnfrozen) * smoothFactor;
    return kMix;
}

Scalar ThermoHydro::equivalentHeatCapacity(Scalar temperature,
                                            const ZoneProperties& zone) const {
    Scalar Tref = FREEZING_POINT - freezingPointDepression_;
    Scalar unfrozenW = unfrozenWaterContent(temperature, zone.waterContent, zone.porosity);
    Scalar w = zone.waterContent;
    Scalar rho = zone.density;
    Scalar cFrozen = zone.heatCapacityFrozen;
    Scalar cUnfrozen = zone.heatCapacityUnfrozen;
    Scalar cDry = 800.0;

    Scalar cMix = cDry + unfrozenW * 4180.0 + (w - unfrozenW) * 2100.0 * (WATER_DENSITY / ICE_DENSITY);

    Scalar smoothWidth = 1.0;
    Scalar factor = std::exp(-std::pow(temperature - Tref, 2) / (2.0 * smoothWidth * smoothWidth));
    Scalar cLatent = LATENT_HEAT * WATER_DENSITY * (-unfrozenWaterDerivative(temperature, w, zone.porosity));
    cMix += cLatent * factor / (std::sqrt(2.0 * M_PI) * smoothWidth);

    if (temperature >= Tref) {
        cMix = cDry + w * 4180.0;
    }
    return std::max(cMix, cDry);
}

Scalar ThermoHydro::hydraulicConductivity(Scalar temperature,
                                           const ZoneProperties& zone) const {
    Scalar Tref = FREEZING_POINT - freezingPointDepression_;
    Scalar kBase = zone.permeability;
    if (temperature >= Tref) {
        return kBase;
    }
    Scalar deltaT = Tref - temperature;
    Scalar reduction = std::exp(-3.0 * deltaT);
    return kBase * (0.0001 + 0.9999 * reduction);
}

Scalar ThermoHydro::waterPotential(Scalar temperature, Scalar waterContent,
                                    const ZoneProperties& zone) const {
    Scalar Tref = FREEZING_POINT - freezingPointDepression_;
    Scalar psiMatric = 0.0;
    if (waterContent > vgThetaR_ && vgThetaS_ > vgThetaR_) {
        Scalar Se = (waterContent - vgThetaR_) / (vgThetaS_ - vgThetaR_);
        Se = std::max(0.01, std::min(Se, 0.99));
        Scalar m = 1.0 - 1.0 / vgN_;
        psiMatric = -1.0 / vgAlpha_ * std::pow(std::pow(Se, -1.0 / m) - 1.0, 1.0 / vgN_);
    }

    Scalar psiThermal = 0.0;
    if (temperature < Tref) {
        Scalar deltaT = Tref - temperature;
        psiThermal = -1.0e5 * deltaT;
    }
    return psiMatric + psiThermal;
}

Scalar ThermoHydro::frostHeaveStrain(Scalar temperature, Scalar initialWaterContent,
                                      Scalar porosity) const {
    Scalar Tref = FREEZING_POINT - freezingPointDepression_;
    if (temperature >= Tref) return 0.0;
    Scalar deltaT = Tref - temperature;
    Scalar ice = iceContent(temperature, initialWaterContent, porosity);
    Scalar volumeExpansion = ice * (ICE_DENSITY / WATER_DENSITY - 1.0);
    Scalar strain = 0.09 * volumeExpansion * (1.0 - std::exp(-deltaT / 3.0));
    return strain;
}

Scalar ThermoHydro::frostHeaveStrainDerivative(Scalar temperature,
                                                Scalar initialWaterContent,
                                                Scalar porosity) const {
    Scalar Tref = FREEZING_POINT - freezingPointDepression_;
    if (temperature >= Tref) return 0.0;
    Scalar deltaT = Tref - temperature;
    Scalar eps = 1e-5;
    Scalar e1 = frostHeaveStrain(temperature + eps, initialWaterContent, porosity);
    Scalar e2 = frostHeaveStrain(temperature - eps, initialWaterContent, porosity);
    return (e2 - e1) / (2.0 * eps);
}

Scalar ThermoHydro::effectiveStressFactor(Scalar temperature, Scalar waterContent,
                                           Scalar porosity) const {
    Scalar Tref = FREEZING_POINT - freezingPointDepression_;
    Scalar ice = iceContent(temperature, waterContent, porosity);
    Scalar iceRatio = ice / (porosity + 1e-10);
    iceRatio = std::max(0.0, std::min(iceRatio, 1.0));
    return 1.0 - 0.7 * iceRatio;
}

void ThermoHydro::computeElementState(Index elemId, Scalar temperature,
                                       Scalar waterContent,
                                       Scalar& kThermal, Scalar& cEquiv,
                                       Scalar& kHydraulic, Scalar& heaveStrain) const {
    const auto& zone = mesh_.getElementZone(elemId);
    kThermal = equivalentThermalConductivity(temperature, zone);
    cEquiv = equivalentHeatCapacity(temperature, zone);
    kHydraulic = hydraulicConductivity(temperature, zone);
    heaveStrain = frostHeaveStrain(temperature, waterContent, zone.porosity);
}

VectorX ThermoHydro::computeInitialTemperatureField(Scalar surfaceTemp,
                                                     Scalar geothermalGradient) const {
    const Index N = mesh_.getNumNodes();
    VectorX T(N);
    Scalar maxY = mesh_.getBoundingBoxMaxY();
    for (Index i = 0; i < N; ++i) {
        Scalar depth = maxY - mesh_.getNode(i).y;
        T(i) = surfaceTemp + geothermalGradient * depth;
    }
    return T;
}

VectorX ThermoHydro::computeInitialWaterContentField() const {
    const Index N = mesh_.getNumNodes();
    VectorX W(N);
    const auto& nodeElems = mesh_.getNodeElements();
    for (Index i = 0; i < N; ++i) {
        if (!nodeElems.empty() && i < static_cast<Index>(nodeElems.size()) && !nodeElems[i].empty()) {
            Scalar avgWC = 0.0;
            for (Index eid : nodeElems[i]) {
                avgWC += mesh_.getElementZone(eid).waterContent;
            }
            W(i) = avgWC / static_cast<Scalar>(nodeElems[i].size());
        } else {
            W(i) = 0.25;
        }
    }
    return W;
}

Scalar ThermoHydro::annualAirTemperature(Scalar time, Scalar meanTemp,
                                          Scalar amplitude, Scalar periodDays,
                                          Scalar phaseShiftDays) const {
    Scalar periodSec = periodDays * 86400.0;
    Scalar phase = phaseShiftDays * 86400.0;
    return meanTemp + amplitude * std::sin(2.0 * M_PI * (time - phase) / periodSec);
}

} // namespace RoadbedSim
