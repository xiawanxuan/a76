#pragma once

#include "Types.h"
#include "Mesh2D.h"
#include <vector>
#include <array>
#include <string>
#include <memory>
#include <functional>
#include <cmath>

namespace RoadbedSim {

enum class AxleType : int {
    Single = 1,
    Tandem = 2,
    Tridem = 3
};

enum class VehicleClass : int {
    LightCar = 0,
    MediumTruck = 1,
    HeavyTruckBZZ100 = 2,
    OversizedLoad = 3,
    Custom = 99
};

struct AxleLoad {
    AxleType type;
    Scalar axleWeight_N;
    Scalar spacing_m;
    Scalar tirePressure_Pa;
    Scalar contactRadius_m;
    Scalar wheelTrack_m;
};

struct VehicleModel {
    VehicleClass vClass;
    std::string name;
    Scalar totalWeight_N;
    Scalar speed_mps;
    std::vector<AxleLoad> axles;
    Scalar wheelBase_m;
    Scalar laneCenterX_m;
    Scalar entryTime_s;
    Scalar exitTime_s;
    Scalar dynamicAmplificationFactor;
    Scalar impactFactor;
};

struct WheelPrint {
    Scalar centerX_m;
    Scalar centerY_m;
    Scalar radius_m;
    Scalar contactPressure_Pa;
    Scalar load_N;
    bool isActive;
    Index vehicleId;
    Index axleId;
    Index wheelId;
};

struct MovingLoadRecord {
    Index stepId;
    Scalar time_s;
    std::vector<WheelPrint> activeWheels;
    Scalar totalDynamicForce_N;
    Scalar equivalentStaticForce_N;
    Scalar dynamicFactor;
};

class VehicleLoad {
public:
    explicit VehicleLoad(const Mesh2D& mesh);
    ~VehicleLoad() = default;

    static VehicleModel createStandardVehicle(VehicleClass vc, Scalar speedKmh,
                                               Scalar laneCenterX, Scalar entryTime);
    static VehicleModel createBZZ100(Scalar speedKmh, Scalar laneCenterX,
                                     Scalar entryTime);
    static VehicleModel createLightCar(Scalar speedKmh, Scalar laneCenterX,
                                       Scalar entryTime);

    Index addVehicle(const VehicleModel& vm);
    Index addVehicleStream(VehicleClass vc, Scalar speedKmh, Scalar laneCenterX,
                           Scalar startTime, Scalar endTime, Scalar headway_s,
                           Index maxCount = 1000);
    Index addVehicleStreamByAADT(VehicleClass vc, Scalar speedKmh, Scalar laneCenterX,
                                 Scalar aadtPerLane, Scalar startTime,
                                 Scalar truckPercent = 10.0, Index maxCount = 500);

    void removeVehicle(Index vehicleId);
    void clearAllVehicles();

    Index getNumVehicles() const { return static_cast<Index>(vehicles_.size()); }
    Index getVehicleCount() const { return static_cast<Index>(vehicles_.size()); }
    bool hasVehicles() const { return !vehicles_.empty(); }
    const VehicleModel& getVehicle(Index id) const { return vehicles_.at(id); }
    const std::vector<VehicleModel>& getAllVehicles() const { return vehicles_; }

    std::vector<WheelPrint> computeActiveWheelsAtTime(Scalar time_s) const;
    std::vector<WheelPrint> computeActiveWheelsInStep(Scalar timeBegin, Scalar timeEnd) const;

    VectorX assembleNodalLoadVector(Scalar time_s) const;
    VectorX assembleNodalLoadVectorFromWheels(const std::vector<WheelPrint>& wheels) const;

    Scalar computeVehiclePositionX(const VehicleModel& vm, Scalar time_s) const;
    bool isVehicleOnRoadbed(const VehicleModel& vm, Scalar time_s) const;

    Scalar computeDynamicAmplification(Scalar temperature, Scalar iceContentRatio,
                                        Scalar baseDAF = 1.15) const;
    Scalar computeImpactFactor(Scalar speedKmh, Scalar roughnessIRI_mPerKm = 2.0) const;

    void computeRoadbedSurfaceNodes();
    const std::vector<Index>& getSurfaceNodeIds() const { return surfaceNodeIds_; }
    const std::vector<Scalar>& getSurfaceNodeXs() const { return surfaceNodeX_; }

    void setSurfaceWidth(Scalar width) { roadbedWidth_ = width; }
    Scalar getSurfaceWidth() const { return roadbedWidth_; }
    void setSurfaceY(Scalar y) { roadbedSurfaceY_ = y; }
    Scalar getSurfaceY() const { return roadbedSurfaceY_; }

    void setContactDistributionMode(int mode) { contactMode_ = mode; }
    int getContactDistributionMode() const { return contactMode_; }

    VectorX computeDynamicStressContribution(Scalar time_s,
                                              const VectorX& dispY) const;

    const std::vector<MovingLoadRecord>& getStepHistory() const { return stepHistory_; }
    void recordStep(Index stepId, Scalar time_s);

    void setRoughnessProfile(std::function<Scalar(Scalar x)> profile) {
        roughnessProfile_ = std::move(profile);
    }
    Scalar getRoughnessAt(Scalar x) const;

    VectorX computeVehiclePositionField(Scalar time_s) const;
    VectorX computeWheelPressureField(Scalar time_s) const;

    Scalar estimateTotalTrafficForce(Scalar time_s) const;

private:
    void initializeSurfaceInfo();
    Scalar interpolateContactPressure(const WheelPrint& wp, Scalar nodeX,
                                       Scalar nodeY) const;
    Scalar ellipticalContactPressure(const WheelPrint& wp, Scalar dx,
                                      Scalar dy) const;
    Scalar circularContactPressure(const WheelPrint& wp, Scalar dx,
                                    Scalar dy) const;

    void generateWheelsForVehicle(const VehicleModel& vm, Scalar posX,
                                  std::vector<WheelPrint>& outWheels) const;

    const Mesh2D& mesh_;
    std::vector<VehicleModel> vehicles_;
    std::vector<Index> surfaceNodeIds_;
    std::vector<Scalar> surfaceNodeX_;
    std::vector<Scalar> surfaceNodeY_;

    Scalar roadbedWidth_ = 20.0;
    Scalar roadbedSurfaceY_ = 0.0;
    Scalar roadbedMinX_ = 0.0;
    Scalar roadbedMaxX_ = 20.0;

    int contactMode_ = 0;
    Scalar defaultTirePressure_ = 700000.0;
    Scalar defaultContactRadius_ = 0.12;

    std::vector<MovingLoadRecord> stepHistory_;
    std::function<Scalar(Scalar x)> roughnessProfile_;
};

} // namespace RoadbedSim
