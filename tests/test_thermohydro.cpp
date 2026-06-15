#define BOOST_TEST_MODULE ThermoHydroTest
#include <boost/test/included/unit_test.hpp>
#include "ThermoHydro.h"
#include "Mesh2D.h"

using namespace RoadbedSim;

struct ThermoHydroFixture {
    Mesh2D mesh;
    std::unique_ptr<ThermoHydro> th;

    ThermoHydroFixture() {
        mesh.addNode(0, 0); mesh.addNode(1, 0); mesh.addNode(0, 1); mesh.addNode(1, 1);
        mesh.addElement({0, 1, 2}, 0);
        mesh.addElement({1, 3, 2}, 0);
        mesh.computeElementAreas();

        ZoneProperties z{};
        z.waterContent = 0.25; z.porosity = 0.40;
        z.youngModulus = 30e6; z.poissonRatio = 0.35;
        z.thermalCondFrozen = 2.2; z.thermalCondUnfrozen = 1.8;
        z.heatCapacityFrozen = 1800; z.heatCapacityUnfrozen = 2200;
        z.density = 1800; z.permeability = 1e-8; z.frostHeaveCoeff = 0.09;
        mesh.addZone(0, z);

        th = std::make_unique<ThermoHydro>(mesh);
    }
};

BOOST_FIXTURE_TEST_CASE(test_unfrozen_water, ThermoHydroFixture) {
    Scalar unfrozenWarm = th->unfrozenWaterContent(280.0, 0.25, 0.40);
    BOOST_CHECK_CLOSE(unfrozenWarm, 0.25, 1e-6);

    Scalar unfrozenCold = th->unfrozenWaterContent(260.0, 0.25, 0.40);
    BOOST_CHECK_LT(unfrozenCold, 0.25);
    BOOST_CHECK_GT(unfrozenCold, 0.0);
}

BOOST_FIXTURE_TEST_CASE(test_ice_content, ThermoHydroFixture) {
    Scalar iceWarm = th->iceContent(280.0, 0.25, 0.40);
    BOOST_CHECK_CLOSE(iceWarm, 0.0, 1e-6);

    Scalar iceCold = th->iceContent(260.0, 0.25, 0.40);
    BOOST_CHECK_GT(iceCold, 0.0);
}

BOOST_FIXTURE_TEST_CASE(test_thermal_conductivity, ThermoHydroFixture) {
    const auto& zone = mesh.getZone(0);
    Scalar kFrozen = th->equivalentThermalConductivity(260.0, zone);
    Scalar kUnfrozen = th->equivalentThermalConductivity(280.0, zone);

    BOOST_CHECK_CLOSE(kUnfrozen, zone.thermalCondUnfrozen, 1e-3);
    BOOST_CHECK_CLOSE(kFrozen, zone.thermalCondFrozen, 5.0);
}

BOOST_FIXTURE_TEST_CASE(test_heat_capacity, ThermoHydroFixture) {
    const auto& zone = mesh.getZone(0);
    Scalar cWarm = th->equivalentHeatCapacity(280.0, zone);
    Scalar cCold = th->equivalentHeatCapacity(260.0, zone);
    BOOST_CHECK_GT(cWarm, 0);
    BOOST_CHECK_GT(cCold, 0);
}

BOOST_FIXTURE_TEST_CASE(test_hydraulic_conductivity, ThermoHydroFixture) {
    const auto& zone = mesh.getZone(0);
    Scalar kWarm = th->hydraulicConductivity(280.0, zone);
    Scalar kCold = th->hydraulicConductivity(260.0, zone);
    BOOST_CHECK_CLOSE(kWarm, zone.permeability, 1e-3);
    BOOST_CHECK_LT(kCold, kWarm);
}

BOOST_FIXTURE_TEST_CASE(test_frost_heave_strain, ThermoHydroFixture) {
    Scalar epsWarm = th->frostHeaveStrain(280.0, 0.25, 0.40);
    BOOST_CHECK_CLOSE(epsWarm, 0.0, 1e-10);

    Scalar epsCold = th->frostHeaveStrain(260.0, 0.25, 0.40);
    BOOST_CHECK_GT(epsCold, 0.0);
    BOOST_CHECK_LT(epsCold, 0.1);
}

BOOST_FIXTURE_TEST_CASE(test_annual_air_temp, ThermoHydroFixture) {
    Scalar mean = 270.0;
    Scalar amp = 20.0;
    Scalar T0 = th->annualAirTemperature(0, mean, amp, 365, 0);
    Scalar Tmid = th->annualAirTemperature(86400.0 * 365 * 0.25, mean, amp, 365, 0);
    BOOST_CHECK_CLOSE(T0, mean, 1e-3);
    BOOST_CHECK_CLOSE(Tmid, mean + amp, 1.0);
}

BOOST_FIXTURE_TEST_CASE(test_initial_fields, ThermoHydroFixture) {
    VectorX T = th->computeInitialTemperatureField(278.0, 0.025);
    VectorX W = th->computeInitialWaterContentField();

    BOOST_CHECK_EQUAL(T.size(), mesh.getNumNodes());
    BOOST_CHECK_EQUAL(W.size(), mesh.getNumNodes());
    BOOST_CHECK_GT(T(0), 200);
    BOOST_CHECK_LT(W(0), 1.0);
    BOOST_CHECK_GT(W(0), 0.0);
}

BOOST_FIXTURE_TEST_CASE(test_element_state, ThermoHydroFixture) {
    Scalar kT, cE, kH, hS;
    th->computeElementState(0, 270.0, 0.25, kT, cE, kH, hS);
    BOOST_CHECK_GT(kT, 0);
    BOOST_CHECK_GT(cE, 0);
    BOOST_CHECK_GT(kH, 0);
}
