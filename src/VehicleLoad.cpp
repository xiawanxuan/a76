#include "VehicleLoad.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <iostream>

namespace RoadbedSim {

VehicleLoad::VehicleLoad(const Mesh2D& mesh) : mesh_(mesh) {
    initializeSurfaceInfo();
}

void VehicleLoad::initializeSurfaceInfo() {
    surfaceNodeIds_ = mesh_.findSurfaceNodes();
    surfaceNodeX_.clear();
    surfaceNodeY_.clear();
    for (Index nid : surfaceNodeIds_) {
        const auto& n = mesh_.getNode(nid);
        surfaceNodeX_.push_back(n.x);
        surfaceNodeY_.push_back(n.y);
    }
    if (!surfaceNodeX_.empty()) {
        roadbedMinX_ = *std::min_element(surfaceNodeX_.begin(), surfaceNodeX_.end());
        roadbedMaxX_ = *std::max_element(surfaceNodeX_.begin(), surfaceNodeX_.end());
        roadbedWidth_ = roadbedMaxX_ - roadbedMinX_;
        if (!surfaceNodeY_.empty()) {
            roadbedSurfaceY_ = *std::max_element(surfaceNodeY_.begin(),
                                                   surfaceNodeY_.end());
        }
    }
    roughnessProfile_ = [](Scalar x) -> Scalar {
        return 0.002 * std::sin(2.0 * M_PI * x / 3.0) +
               0.001 * std::sin(2.0 * M_PI * x / 0.5);
    };
}

VehicleModel VehicleLoad::createStandardVehicle(VehicleClass vc, Scalar speedKmh,
                                                 Scalar laneCenterX,
                                                 Scalar entryTime) {
    Scalar speedMps = speedKmh / 3.6;
    VehicleModel vm{};
    vm.vClass = vc;
    vm.speed_mps = speedMps;
    vm.laneCenterX_m = laneCenterX;
    vm.entryTime_s = entryTime;
    vm.dynamicAmplificationFactor = 1.15;
    vm.impactFactor = 1.30;

    switch (vc) {
        case VehicleClass::LightCar:
            vm = createLightCar(speedKmh, laneCenterX, entryTime);
            break;
        case VehicleClass::MediumTruck: {
            vm.name = "MediumTruck_15t";
            vm.totalWeight_N = 15.0 * 1000.0 * GRAVITY;
            AxleLoad a1{AxleType::Single, 5.0 * 1000 * GRAVITY, 0.0, 700000.0, 0.10, 1.8};
            AxleLoad a2{AxleType::Single, 10.0 * 1000 * GRAVITY, 3.8, 700000.0, 0.10, 1.8};
            vm.axles = {a1, a2};
            vm.wheelBase_m = 3.8;
            vm.exitTime_s = entryTime + 10.0;
            vm.impactFactor = 1.40;
            break;
        }
        case VehicleClass::HeavyTruckBZZ100:
            vm = createBZZ100(speedKmh, laneCenterX, entryTime);
            break;
        case VehicleClass::OversizedLoad: {
            vm.name = "Oversized_55t";
            vm.totalWeight_N = 55.0 * 1000 * GRAVITY;
            AxleLoad a1{AxleType::Single, 8.0 * 1000 * GRAVITY, 0.0, 800000.0, 0.15, 2.0};
            AxleLoad a2{AxleType::Tandem, 23.5 * 1000 * GRAVITY, 3.6, 800000.0, 0.15, 2.0};
            AxleLoad a3{AxleType::Tridem, 23.5 * 1000 * GRAVITY, 9.5, 800000.0, 0.15, 2.0};
            vm.axles = {a1, a2, a3};
            vm.wheelBase_m = 9.5;
            vm.exitTime_s = entryTime + 15.0;
            vm.impactFactor = 1.50;
            vm.dynamicAmplificationFactor = 1.30;
            break;
        }
        default:
            vm = createBZZ100(speedKmh, laneCenterX, entryTime);
            break;
    }
    return vm;
}

VehicleModel VehicleLoad::createBZZ100(Scalar speedKmh, Scalar laneCenterX,
                                        Scalar entryTime) {
    VehicleModel vm{};
    vm.vClass = VehicleClass::HeavyTruckBZZ100;
    vm.name = "BZZ-100_Standard";
    Scalar speedMps = speedKmh / 3.6;
    vm.speed_mps = speedMps;
    vm.laneCenterX_m = laneCenterX;
    vm.entryTime_s = entryTime;
    vm.totalWeight_N = 100.0 * 1000.0 * GRAVITY;
    vm.wheelBase_m = 4.0 + 1.4;

    Scalar singleW = 60.0 * 1000.0 * GRAVITY;
    Scalar tandemW = 2 * 20.0 * 1000.0 * GRAVITY;

    AxleLoad steering{AxleType::Single, singleW, 0.0,
                       700000.0, 0.11, 1.82};
    AxleLoad drive{AxleType::Tandem, tandemW, 4.0,
                    700000.0, 0.11, 1.82};
    drive.spacing_m = 1.4;
    vm.axles = {steering, drive};

    Scalar roadLen = 100.0;
    vm.exitTime_s = entryTime + roadLen / std::max(0.1, speedMps) + 2.0;
    vm.dynamicAmplificationFactor = 1.20;
    vm.impactFactor = computeImpactFactor(speedKmh, 2.5);
    return vm;
}

VehicleModel VehicleLoad::createLightCar(Scalar speedKmh, Scalar laneCenterX,
                                          Scalar entryTime) {
    VehicleModel vm{};
    vm.vClass = VehicleClass::LightCar;
    vm.name = "LightCar_1.5t";
    Scalar speedMps = speedKmh / 3.6;
    vm.speed_mps = speedMps;
    vm.laneCenterX_m = laneCenterX;
    vm.entryTime_s = entryTime;
    vm.totalWeight_N = 1.5 * 1000.0 * GRAVITY;
    vm.wheelBase_m = 2.6;

    AxleLoad front{AxleType::Single, 0.55 * vm.totalWeight_N, 0.0,
                    250000.0, 0.08, 1.5};
    AxleLoad rear{AxleType::Single, 0.45 * vm.totalWeight_N, 2.6,
                   250000.0, 0.08, 1.5};
    vm.axles = {front, rear};

    Scalar roadLen = 100.0;
    vm.exitTime_s = entryTime + roadLen / std::max(0.1, speedMps) + 1.0;
    vm.dynamicAmplificationFactor = 1.10;
    vm.impactFactor = computeImpactFactor(speedKmh, 1.5);
    return vm;
}

Index VehicleLoad::addVehicle(const VehicleModel& vm) {
    vehicles_.push_back(vm);
    return static_cast<Index>(vehicles_.size() - 1);
}

Index VehicleLoad::addVehicleStream(VehicleClass vc, Scalar speedKmh,
                                     Scalar laneCenterX, Scalar startTime,
                                     Scalar endTime, Scalar headway_s,
                                     Index maxCount) {
    Index count = 0;
    for (Scalar t = startTime; t < endTime && count < maxCount; t += headway_s) {
        addVehicle(createStandardVehicle(vc, speedKmh, laneCenterX, t));
        ++count;
    }
    return count;
}

Index VehicleLoad::addVehicleStreamByAADT(VehicleClass vc, Scalar speedKmh,
                                           Scalar laneCenterX,
                                           Scalar aadtPerLane,
                                           Scalar startTime,
                                           Scalar truckPercent,
                                           Index maxCount) {
    Scalar flowRatePerSec = aadtPerLane / 86400.0;
    Scalar headwayMean = 1.0 / std::max(0.001, flowRatePerSec);

    if (vc != VehicleClass::LightCar) {
        flowRatePerSec *= truckPercent / 100.0;
        headwayMean = 1.0 / std::max(0.0001, flowRatePerSec);
    }

    Index count = 0;
    Scalar t = startTime;
    while (count < maxCount) {
        Scalar headway = headwayMean * (0.6 + 0.8 * std::fmod(sin(t * 1.7 + count * 0.3) + 1.0, 1.0));
        t += headway;
        addVehicle(createStandardVehicle(vc, speedKmh, laneCenterX, t));
        ++count;
        if (t > startTime + 86400.0 * 3) break;
    }
    return count;
}

void VehicleLoad::removeVehicle(Index vehicleId) {
    if (vehicleId < static_cast<Index>(vehicles_.size())) {
        vehicles_.erase(vehicles_.begin() + vehicleId);
    }
}

void VehicleLoad::clearAllVehicles() {
    vehicles_.clear();
    stepHistory_.clear();
}

Scalar VehicleLoad::computeVehiclePositionX(const VehicleModel& vm,
                                             Scalar time_s) const {
    if (time_s < vm.entryTime_s) return roadbedMinX_ - 50.0;
    if (time_s > vm.exitTime_s) return roadbedMaxX_ + 50.0;
    Scalar dt = time_s - vm.entryTime_s;
    Scalar accelPhase = std::min(1.0, dt / 3.0);
    Scalar effectiveSpeed = vm.speed_mps * accelPhase;
    return roadbedMinX_ - 5.0 + effectiveSpeed * dt;
}

bool VehicleLoad::isVehicleOnRoadbed(const VehicleModel& vm,
                                      Scalar time_s) const {
    if (time_s < vm.entryTime_s || time_s > vm.exitTime_s) return false;
    Scalar posX = computeVehiclePositionX(vm, time_s);
    Scalar margin = vm.wheelBase_m + 3.0;
    return (posX >= roadbedMinX_ - margin && posX <= roadbedMaxX_ + margin);
}

Scalar VehicleLoad::computeDynamicAmplification(Scalar temperature,
                                                 Scalar iceContentRatio,
                                                 Scalar baseDAF) const {
    Scalar deltaT = FREEZING_POINT - temperature;
    Scalar tempFactor = 1.0;
    if (deltaT > 0) {
        tempFactor = 1.0 + 0.015 * std::min(deltaT, 15.0);
    }
    Scalar iceFactor = 1.0 + 0.25 * std::max(0.0, std::min(iceContentRatio, 0.5));
    Scalar thawFactor = 1.0;
    if (temperature > FREEZING_POINT - 0.5 && temperature < FREEZING_POINT + 1.0) {
        thawFactor = 1.15 + 0.3 * std::exp(-std::pow(temperature - FREEZING_POINT, 2) / 0.5);
    }
    return baseDAF * tempFactor * iceFactor * thawFactor;
}

Scalar VehicleLoad::computeImpactFactor(Scalar speedKmh,
                                         Scalar roughnessIRI_mPerKm) const {
    Scalar v2 = speedKmh * speedKmh;
    Scalar iriFactor = 1.0 + roughnessIRI_mPerKm * 0.04;
    Scalar speedFactor = 1.0 + v2 / (80.0 * 80.0) * 0.5;
    return 1.0 + (iriFactor * speedFactor - 1.0) * 0.6;
}

Scalar VehicleLoad::getRoughnessAt(Scalar x) const {
    if (roughnessProfile_) return roughnessProfile_(x);
    return 0.0;
}

void VehicleLoad::generateWheelsForVehicle(const VehicleModel& vm, Scalar posX,
                                             std::vector<WheelPrint>& outWheels) const {
    Index vehicleId = 0;
    for (size_t i = 0; i < vehicles_.size(); ++i) {
        if (&vehicles_[i] == &vm) {
            vehicleId = static_cast<Index>(i);
            break;
        }
    }

    Scalar avgIce = 0.0;
    Scalar avgT = FREEZING_POINT - 1.0;
    for (size_t i = 0; i < surfaceNodeIds_.size() && i < 10; ++i) {
        Index nid = surfaceNodeIds_[i];
        const auto& zone = mesh_.getElementZone(
            mesh_.getNodeElements().empty() ? 0 :
            (mesh_.getNodeElements()[nid].empty() ? 0 : mesh_.getNodeElements()[nid][0]));
        avgIce += zone.waterContent;
        avgT -= 0.1;
    }
    Scalar iceRatio = surfaceNodeIds_.empty() ? 0.2 :
                      std::min(1.0, avgIce / std::max(size_t(1), surfaceNodeIds_.size()));
    Scalar DAF = computeDynamicAmplification(avgT, iceRatio,
                                              vm.dynamicAmplificationFactor);

    for (size_t ai = 0; ai < vm.axles.size(); ++ai) {
        const auto& axle = vm.axles[ai];
        Scalar axleX = posX - axle.spacing_m;
        int numWheelsPerSide = (axle.type == AxleType::Single) ? 1 :
                               (axle.type == AxleType::Tandem) ? 2 : 3;
        Scalar tandemSpacing = (axle.type == AxleType::Tandem) ? 1.2 :
                                (axle.type == AxleType::Tridem) ? 1.2 : 0.0;

        for (int side = 0; side < 2; ++side) {
            Scalar offsetY = (side == 0) ? -axle.wheelTrack_m * 0.5 :
                              axle.wheelTrack_m * 0.5;
            Scalar wheelY = vm.laneCenterX_m + offsetY;

            for (int w = 0; w < numWheelsPerSide; ++w) {
                Scalar wx = axleX - (numWheelsPerSide - 1) * tandemSpacing * 0.5
                           + w * tandemSpacing;
                Scalar loadPerWheel = axle.axleWeight_N / (2.0 * numWheelsPerSide);
                Scalar dynLoad = loadPerWheel * DAF * vm.impactFactor;

                WheelPrint wp{};
                wp.centerX_m = wx;
                wp.centerY_m = wheelY;
                wp.radius_m = axle.contactRadius_m;
                wp.load_N = dynLoad;
                wp.contactPressure_Pa = axle.tirePressure_Pa * DAF;
                wp.isActive = true;
                wp.vehicleId = vehicleId;
                wp.axleId = static_cast<Index>(ai);
                wp.wheelId = static_cast<Index>(side * numWheelsPerSide + w);
                outWheels.push_back(wp);
            }
        }
    }
}

std::vector<WheelPrint> VehicleLoad::computeActiveWheelsAtTime(
    Scalar time_s) const {
    std::vector<WheelPrint> wheels;
    wheels.reserve(vehicles_.size() * 8);
    for (const auto& vm : vehicles_) {
        if (!isVehicleOnRoadbed(vm, time_s)) continue;
        Scalar posX = computeVehiclePositionX(vm, time_s);
        generateWheelsForVehicle(vm, posX, wheels);
    }
    return wheels;
}

std::vector<WheelPrint> VehicleLoad::computeActiveWheelsInStep(
    Scalar timeBegin, Scalar timeEnd) const {
    std::vector<WheelPrint> allWheels;
    const int subSteps = 5;
    Scalar dt = (timeEnd - timeBegin) / subSteps;
    for (int s = 0; s <= subSteps; ++s) {
        Scalar t = timeBegin + s * dt;
        auto wheels = computeActiveWheelsAtTime(t);
        Scalar weight = (s == 0 || s == subSteps) ? 0.5 : 1.0;
        for (auto& wp : wheels) {
            wp.load_N *= weight / (subSteps + 1.0);
            wp.contactPressure_Pa *= weight / (subSteps + 1.0);
            allWheels.push_back(wp);
        }
    }
    return allWheels;
}

Scalar VehicleLoad::circularContactPressure(const WheelPrint& wp, Scalar dx,
                                              Scalar dy) const {
    Scalar r = wp.radius_m;
    Scalar dist2 = dx * dx + dy * dy;
    if (dist2 > r * r) return 0.0;
    return wp.contactPressure_Pa;
}

Scalar VehicleLoad::ellipticalContactPressure(const WheelPrint& wp, Scalar dx,
                                               Scalar dy) const {
    Scalar a = wp.radius_m * 1.2;
    Scalar b = wp.radius_m * 0.8;
    Scalar val = (dx * dx) / (a * a) + (dy * dy) / (b * b);
    if (val > 1.0) return 0.0;
    Scalar reduction = 1.0 - 0.3 * val;
    return wp.contactPressure_Pa * reduction;
}

Scalar VehicleLoad::interpolateContactPressure(const WheelPrint& wp,
                                                Scalar nodeX,
                                                Scalar nodeY) const {
    Scalar dx = nodeX - wp.centerX_m;
    Scalar dy = nodeY - wp.centerY_m;
    if (contactMode_ == 1) {
        return ellipticalContactPressure(wp, dx, dy);
    }
    return circularContactPressure(wp, dx, dy);
}

VectorX VehicleLoad::assembleNodalLoadVectorFromWheels(
    const std::vector<WheelPrint>& wheels) const {
    Index N = mesh_.getNumNodes();
    VectorX F = VectorX::Zero(4 * N);

    for (const auto& wp : wheels) {
        if (!wp.isActive) continue;
        Scalar r = wp.radius_m * 2.5;
        for (size_t si = 0; si < surfaceNodeIds_.size(); ++si) {
            Index nid = surfaceNodeIds_[si];
            Scalar nx = surfaceNodeX_[si];
            Scalar ny = roadbedSurfaceY_;
            Scalar dx = nx - wp.centerX_m;
            Scalar dy = ny - wp.centerY_m;
            if (std::abs(dx) > r || std::abs(dy) > r) continue;

            Scalar pressure = interpolateContactPressure(wp, nx, ny);
            if (pressure <= 0) continue;

            Scalar spacingX = (roadbedWidth_ > 0 && surfaceNodeIds_.size() > 1)
                ? roadbedWidth_ / (surfaceNodeIds_.size() - 1) : 1.0;
            Scalar spacingY = (wp.radius_m * 0.5 + spacingX * 0.5);
            Scalar area = spacingX * spacingY;
            Scalar force = -pressure * area;

            Index yDof = 2 * N + 2 * nid + 1;
            if (yDof >= 0 && yDof < 4 * N) {
                F(yDof) += force;
            }

            Index xDof = 2 * N + 2 * nid;
            if (xDof >= 0 && xDof < 4 * N && std::abs(dx) < wp.radius_m) {
                Scalar friction = 0.3;
                Scalar fricX = -force * friction * std::tanh(dx / std::max(0.05, wp.radius_m));
                F(xDof) += fricX;
            }
        }
    }
    return F;
}

VectorX VehicleLoad::assembleNodalLoadVector(Scalar time_s) const {
    auto wheels = computeActiveWheelsAtTime(time_s);
    return assembleNodalLoadVectorFromWheels(wheels);
}

VectorX VehicleLoad::computeDynamicStressContribution(
    Scalar time_s, const VectorX& dispY) const {
    Index N = mesh_.getNumNodes();
    VectorX dynStress = VectorX::Zero(N);
    auto wheels = computeActiveWheelsAtTime(time_s);
    if (wheels.empty()) return dynStress;

    for (size_t si = 0; si < surfaceNodeIds_.size(); ++si) {
        Index nid = surfaceNodeIds_[si];
        Scalar nx = surfaceNodeX_[si];
        Scalar totalP = 0.0;
        for (const auto& wp : wheels) {
            Scalar dx = nx - wp.centerX_m;
            Scalar dy = roadbedSurfaceY_ - wp.centerY_m;
            totalP += interpolateContactPressure(wp, dx, dy);
        }
        if (nid < N) dynStress(nid) = totalP;
    }
    return dynStress;
}

VectorX VehicleLoad::computeVehiclePositionField(Scalar time_s) const {
    Index N = mesh_.getNumNodes();
    VectorX field = VectorX::Zero(N);
    auto wheels = computeActiveWheelsAtTime(time_s);
    for (const auto& wp : wheels) {
        Scalar r = wp.radius_m * 2.0;
        for (Index i = 0; i < N; ++i) {
            const auto& n = mesh_.getNode(i);
            Scalar dx = n.x - wp.centerX_m;
            Scalar dy = n.y - wp.centerY_m;
            if (dx * dx + dy * dy <= r * r) {
                field(i) = std::max(field(i), static_cast<Scalar>(wp.vehicleId + 1));
            }
        }
    }
    return field;
}

VectorX VehicleLoad::computeWheelPressureField(Scalar time_s) const {
    Index N = mesh_.getNumNodes();
    VectorX field = VectorX::Zero(N);
    auto wheels = computeActiveWheelsAtTime(time_s);
    for (const auto& wp : wheels) {
        Scalar r = wp.radius_m * 3.0;
        for (Index i = 0; i < N; ++i) {
            const auto& n = mesh_.getNode(i);
            Scalar dx = n.x - wp.centerX_m;
            Scalar dy = n.y - wp.centerY_m;
            Scalar p = interpolateContactPressure(wp, n.x, n.y);
            if (dx * dx + dy * dy <= r * r && p > 0) {
                field(i) = std::max(field(i), p);
            }
        }
    }
    return field;
}

Scalar VehicleLoad::estimateTotalTrafficForce(Scalar time_s) const {
    auto wheels = computeActiveWheelsAtTime(time_s);
    Scalar total = 0.0;
    for (const auto& wp : wheels) total += wp.load_N;
    return total;
}

void VehicleLoad::recordStep(Index stepId, Scalar time_s) {
    MovingLoadRecord rec;
    rec.stepId = stepId;
    rec.time_s = time_s;
    rec.activeWheels = computeActiveWheelsAtTime(time_s);
    rec.equivalentStaticForce_N = 0.0;
    rec.totalDynamicForce_N = 0.0;
    for (const auto& wp : rec.activeWheels) {
        rec.totalDynamicForce_N += wp.load_N;
    }
    rec.dynamicFactor = rec.equivalentStaticForce_N > 1.0
        ? rec.totalDynamicForce_N / rec.equivalentStaticForce_N : 1.0;
    stepHistory_.push_back(rec);
}

} // namespace RoadbedSim
