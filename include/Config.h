#pragma once

#include "Types.h"
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <string>
#include <map>
#include <vector>
#include <tuple>
#include <iostream>
#include <fstream>
#include <sstream>

namespace RoadbedSim {

struct SimulationConfig {
    std::string meshFile;
    std::string outputDir = "./vtk_output";
    std::string configFile;
    std::string logFile = "simulation.log";

    Scalar initialSurfaceTemp = 278.15;
    Scalar geothermalGradient = 0.025;
    Scalar meanAnnualTemp = 270.15;
    Scalar tempAmplitude = 20.0;
    Scalar periodDays = 365.0;
    Scalar phaseShiftDays = 15.0;

    Scalar totalTime = 3.154e7;
    Scalar timeStep = 3600.0;
    Index outputInterval = 24;

    Index maxIterations = 50;
    Scalar tolerance = 1e-6;
    bool useDirectSolver = true;

    Scalar bottomBoundaryTemp = 280.15;
    Scalar initialWaterContentDefault = 0.25;

    bool enableDamageAnalysis = true;
    Scalar damageThreshold = 0.5;

    Index blockSize = 1000;
    bool enableBlockLoading = true;

    bool enableTrafficCoupling = false;
    Scalar vehicleSpeedKmh = 60.0;
    Scalar truckCount = 50;
    Scalar lightCarCount = 200;
    std::vector<Scalar> laneCenterXs = {0.0, 3.5};
    Scalar AADT_perLane = 20000;
    Scalar truckPercent = 15.0;
    Scalar loadFrequencyHz = 10.0;
    bool useBZZ100Standard = true;
    Scalar startTime_s = 0.0;
    Scalar trafficDuration_s = 86400.0;
    Scalar headway_s = 6.0;
    Scalar surfaceElevation = 0.0;

    std::map<Index, ZoneProperties> customZones;
    std::vector<std::tuple<Scalar, Scalar, Scalar, Scalar, Index, ZoneProperties>> zoneRanges;
};

class Config {
public:
    Config();
    ~Config() = default;

    bool parseCommandLine(int argc, char* argv[]);
    bool loadFromFile(const std::string& filename);
    bool saveToFile(const std::string& filename) const;

    const SimulationConfig& get() const { return cfg_; }
    SimulationConfig& mutableConfig() { return cfg_; }

    void printSummary(std::ostream& os = std::cout) const;
    void validate() const;

    std::string getHelpMessage() const;

private:
    void setupDefaultZones();
    void setupDescription();

    SimulationConfig cfg_;
    boost::program_options::options_description desc_;
    boost::program_options::variables_map vm_;
};

} // namespace RoadbedSim
