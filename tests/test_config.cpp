#define BOOST_TEST_MODULE ConfigTest
#include <boost/test/included/unit_test.hpp>
#include "Config.h"
#include <fstream>
#include <filesystem>

using namespace RoadbedSim;

BOOST_AUTO_TEST_CASE(test_default_config) {
    Config config;
    const auto& cfg = config.get();
    BOOST_CHECK_CLOSE(cfg.totalTime, 3.154e7, 1e-3);
    BOOST_CHECK_CLOSE(cfg.timeStep, 3600.0, 1e-10);
    BOOST_CHECK_EQUAL(cfg.outputInterval, 24);
    BOOST_CHECK_EQUAL(cfg.maxIterations, 50);
    BOOST_CHECK_CLOSE(cfg.tolerance, 1e-6, 1e-15);
}

BOOST_AUTO_TEST_CASE(test_modify_config) {
    Config config;
    auto& cfg = config.mutableConfig();
    cfg.timeStep = 1800.0;
    cfg.maxIterations = 30;
    cfg.tolerance = 1e-7;
    BOOST_CHECK_CLOSE(config.get().timeStep, 1800.0, 1e-10);
    BOOST_CHECK_EQUAL(config.get().maxIterations, 30);
}

BOOST_AUTO_TEST_CASE(test_custom_zones) {
    Config config;
    auto& cfg = config.mutableConfig();

    ZoneProperties z{};
    z.waterContent = 0.5;
    z.youngModulus = 1e8;
    z.porosity = 0.45;
    cfg.customZones[5] = z;

    BOOST_CHECK_EQUAL(cfg.customZones.count(5), 1);
    BOOST_CHECK_CLOSE(cfg.customZones[5].waterContent, 0.5, 1e-10);
    BOOST_CHECK_CLOSE(cfg.customZones[5].youngModulus, 1e8, 1e-3);
}

BOOST_AUTO_TEST_CASE(test_save_load_config) {
    const std::string cfgFile = "test_config.json";

    Config config1;
    auto& cfg1 = config1.mutableConfig();
    cfg1.meshFile = "test.msh";
    cfg1.meanAnnualTemp = 271.0;
    cfg1.tempAmplitude = 18.0;
    cfg1.timeStep = 7200.0;

    ZoneProperties z{};
    z.waterContent = 0.33;
    z.porosity = 0.42;
    cfg1.customZones[7] = z;

    BOOST_CHECK(config1.saveToFile(cfgFile));

    Config config2;
    BOOST_CHECK(config2.loadFromFile(cfgFile));
    const auto& cfg2 = config2.get();

    BOOST_CHECK_EQUAL(cfg2.meshFile, "test.msh");
    BOOST_CHECK_CLOSE(cfg2.meanAnnualTemp, 271.0, 1e-3);
    BOOST_CHECK_CLOSE(cfg2.timeStep, 7200.0, 1e-10);

    std::remove(cfgFile.c_str());
}

BOOST_AUTO_TEST_CASE(test_validate_config) {
    Config config;
    auto& cfg = config.mutableConfig();
    cfg.meshFile = "dummy.msh";
    BOOST_CHECK_NO_THROW(config.validate());

    cfg.timeStep = -1;
    BOOST_CHECK_THROW(config.validate(), std::exception);
}

BOOST_AUTO_TEST_CASE(test_help_message) {
    Config config;
    std::string help = config.getHelpMessage();
    BOOST_CHECK(help.size() > 0);
    BOOST_CHECK(help.find("--help") != std::string::npos);
    BOOST_CHECK(help.find("--mesh") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(test_summary_output) {
    Config config;
    auto& cfg = config.mutableConfig();
    cfg.meshFile = "test.msh";

    std::ostringstream oss;
    config.printSummary(oss);
    BOOST_CHECK(!oss.str().empty());
}

BOOST_AUTO_TEST_CASE(test_zone_ranges) {
    Config config;
    auto& cfg = config.mutableConfig();

    ZoneProperties p{};
    p.waterContent = 0.4;
    cfg.zoneRanges.emplace_back(0.0, 10.0, 0.0, 2.0, 1, p);
    BOOST_CHECK_EQUAL(cfg.zoneRanges.size(), 1);
}
