#include "Types.h"
#include "Mesh2D.h"
#include "VehicleLoad.h"
#include "Assembler.h"
#include "ThermoHydro.h"
#include <boost/test/included/unit_test.hpp>
#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <memory>

using namespace RoadbedSim;

struct VehicleLoadFixture {
    Mesh2D mesh;
    VehicleLoad* vl;
    std::unique_ptr<VehicleLoad> vlPtr;

    VehicleLoadFixture() {
        std::vector<Node2D> nodes;
        for (int y = 0; y < 6; ++y) {
            for (int x = 0; x < 20; ++x) {
                nodes.push_back({static_cast<Scalar>(x * 0.5),
                                 static_cast<Scalar>(y * 0.5 - 2.5)});
            }
        }
        mesh.setNodes(nodes);

        std::vector<TriElement> elems;
        int Nx = 20, Ny = 6;
        for (int y = 0; y < Ny - 1; ++y) {
            for (int x = 0; x < Nx - 1; ++x) {
                int i0 = y * Nx + x;
                int i1 = y * Nx + x + 1;
                int i2 = (y + 1) * Nx + x;
                int i3 = (y + 1) * Nx + x + 1;
                elems.push_back({{i0, i1, i3}, 0, 0.0, {}});
                elems.push_back({{i0, i3, i2}, 0, 0.0, {}});
            }
        }
        for (auto& e : elems) {
            Scalar x0 = mesh.getNode(e.nodeIds[0]).x;
            Scalar y0 = mesh.getNode(e.nodeIds[0]).y;
            Scalar x1 = mesh.getNode(e.nodeIds[1]).x;
            Scalar y1 = mesh.getNode(e.nodeIds[1]).y;
            Scalar x2 = mesh.getNode(e.nodeIds[2]).x;
            Scalar y2 = mesh.getNode(e.nodeIds[2]).y;
            e.area = 0.5 * std::abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0));
            Scalar cx = (x0 + x1 + x2) / 3.0;
            (void)cx;
        }
        mesh.setElements(elems);
        mesh.recomputeAll();

        vlPtr = std::make_unique<VehicleLoad>(mesh);
        vl = vlPtr.get();
    }

    ~VehicleLoadFixture() = default;
};

BOOST_FIXTURE_TEST_SUITE(VehicleLoadTests, VehicleLoadFixture)

BOOST_AUTO_TEST_CASE(test_create_bzz100_standard) {
    auto vm = VehicleLoad::createBZZ100(60.0, 5.0, 0.0);
    BOOST_CHECK_EQUAL(vm.vClass, VehicleClass::HeavyTruckBZZ100);
    BOOST_CHECK_CLOSE(vm.speed_mps, 60.0 / 3.6, 1e-3);
    BOOST_CHECK_GT(vm.totalWeight_N, 98000.0 * 0.9);
    BOOST_CHECK_LT(vm.totalWeight_N, 102000.0 * 1.1);
    BOOST_CHECK_EQUAL(vm.axles.size(), 2u);
    BOOST_CHECK_EQUAL(vm.axles[0].type, AxleType::Single);
    BOOST_CHECK_EQUAL(vm.axles[1].type, AxleType::Tandem);
    BOOST_CHECK_CLOSE(vm.laneCenterX_m, 5.0, 1e-6);
    BOOST_CHECK_CLOSE(vm.entryTime_s, 0.0, 1e-6);
    BOOST_CHECK_GT(vm.axles[0].axleWeight_N, 50000.0);
    BOOST_CHECK_CLOSE(vm.axles[0].contactRadius_m, 0.11, 1e-3);
}

BOOST_AUTO_TEST_CASE(test_create_light_car) {
    auto vm = VehicleLoad::createLightCar(100.0, 7.5, 10.0);
    BOOST_CHECK_EQUAL(vm.vClass, VehicleClass::LightCar);
    BOOST_CHECK_EQUAL(vm.axles.size(), 2u);
    BOOST_CHECK_CLOSE(vm.speed_mps, 100.0 / 3.6, 1e-3);
    BOOST_CHECK_CLOSE(vm.laneCenterX_m, 7.5, 1e-6);
    BOOST_CHECK_CLOSE(vm.entryTime_s, 10.0, 1e-6);
    BOOST_CHECK_CLOSE(vm.axles[0].contactRadius_m, 0.08, 1e-3);
}

BOOST_AUTO_TEST_CASE(test_vehicle_add_and_count) {
    BOOST_CHECK(!vl->hasVehicles());
    BOOST_CHECK_EQUAL(vl->getVehicleCount(), 0u);

    auto vm1 = VehicleLoad::createBZZ100(60.0, 5.0, 0.0);
    auto vm2 = VehicleLoad::createLightCar(90.0, 7.5, 5.0);

    Index id1 = vl->addVehicle(vm1);
    Index id2 = vl->addVehicle(vm2);

    BOOST_CHECK_EQUAL(id1, 0);
    BOOST_CHECK_EQUAL(id2, 1);
    BOOST_CHECK_EQUAL(vl->getVehicleCount(), 2u);
    BOOST_CHECK(vl->hasVehicles());

    const auto& vm = vl->getVehicle(0);
    BOOST_CHECK_EQUAL(vm.vClass, VehicleClass::HeavyTruckBZZ100);
}

BOOST_AUTO_TEST_CASE(test_vehicle_stream) {
    Index count = vl->addVehicleStream(
        VehicleClass::HeavyTruckBZZ100,
        60.0, 5.0,
        0.0, 30.0, 5.0, 20);
    BOOST_CHECK_EQUAL(count, 7u);
    BOOST_CHECK_EQUAL(vl->getVehicleCount(), 7u);

    auto vehicles = vl->getAllVehicles();
    for (size_t i = 1; i < vehicles.size(); ++i) {
        Scalar dt = vehicles[i].entryTime_s - vehicles[i - 1].entryTime_s;
        BOOST_CHECK_CLOSE(dt, 5.0, 1e-3);
    }
}

BOOST_AUTO_TEST_CASE(test_vehicle_position_x) {
    auto vm = VehicleLoad::createBZZ100(36.0, 5.0, 0.0);
    Scalar speedMps = 36.0 / 3.6;

    Scalar p0 = vl->computeVehiclePositionX(vm, -1.0);
    Scalar p1 = vl->computeVehiclePositionX(vm, 0.0);
    Scalar p2 = vl->computeVehiclePositionX(vm, 1.0);
    Scalar p3 = vl->computeVehiclePositionX(vm, 1000.0);

    BOOST_CHECK_LT(p0, 0.0);
    BOOST_CHECK_LT(p1, p2);
    Scalar expected = p1 + speedMps * 1.0 * 1.0;
    BOOST_CHECK_CLOSE(p2, expected, speedMps * 0.5);
    BOOST_CHECK_GT(p3, p2);
}

BOOST_AUTO_TEST_CASE(test_is_vehicle_on_roadbed) {
    auto vm = VehicleLoad::createBZZ100(36.0, 5.0, 0.0);
    BOOST_CHECK(!vl->isVehicleOnRoadbed(vm, -100.0));
    BOOST_CHECK(vl->isVehicleOnRoadbed(vm, 2.0));
}

BOOST_AUTO_TEST_CASE(test_compute_dynamic_amplification) {
    Scalar baseDAF = 1.15;

    Scalar warmFrozen = vl->computeDynamicAmplification(272.0, 0.05, baseDAF);
    Scalar coldFrozen = vl->computeDynamicAmplification(268.0, 0.40, baseDAF);
    Scalar thawWeak = vl->computeDynamicAmplification(273.0, 0.20, baseDAF);
    Scalar warmUnfrozen = vl->computeDynamicAmplification(278.0, 0.00, baseDAF);

    BOOST_CHECK_GE(coldFrozen, warmFrozen);
    BOOST_CHECK_GE(thawWeak, baseDAF * 1.1);
    BOOST_CHECK_GE(warmFrozen, baseDAF);
    BOOST_CHECK_GE(warmUnfrozen, baseDAF * 0.9);
    BOOST_CHECK_LE(coldFrozen, baseDAF * 2.0);
}

BOOST_AUTO_TEST_CASE(test_impact_factor) {
    Scalar if1 = VehicleLoad::computeImpactFactor(40.0, 1.0);
    Scalar if2 = VehicleLoad::computeImpactFactor(80.0, 1.0);
    Scalar if3 = VehicleLoad::computeImpactFactor(80.0, 5.0);

    BOOST_CHECK_GE(if2, if1);
    BOOST_CHECK_GE(if3, if2);
    BOOST_CHECK_GE(if1, 1.0);
    BOOST_CHECK_LE(if1, 2.0);
}

BOOST_AUTO_TEST_CASE(test_active_wheels_at_time) {
    auto vm = VehicleLoad::createBZZ100(36.0, 5.0, 0.0);
    vl->addVehicle(vm);

    auto wheelsBefore = vl->computeActiveWheelsAtTime(-100.0);
    auto wheelsOn = vl->computeActiveWheelsAtTime(5.0);

    BOOST_CHECK_EQUAL(wheelsBefore.size(), 0u);
    BOOST_CHECK_GT(wheelsOn.size(), 0u);

    for (const auto& wp : wheelsOn) {
        BOOST_CHECK_GT(wp.radius_m, 0.0);
        BOOST_CHECK_GT(wp.load_N, 0.0);
        BOOST_CHECK(wp.isActive);
        BOOST_CHECK_GT(wp.contactPressure_Pa, 0.0);
    }
}

BOOST_AUTO_TEST_CASE(test_circular_contact_pressure) {
    WheelPrint wp{};
    wp.centerX_m = 5.0;
    wp.centerY_m = 0.1;
    wp.radius_m = 0.11;
    wp.contactPressure_Pa = 700000.0;
    wp.isActive = true;

    Scalar pCenter = vl->interpolateContactPressure(wp, 5.0, 0.1);
    Scalar pInside = vl->interpolateContactPressure(wp, 5.05, 0.1);
    Scalar pEdge = vl->interpolateContactPressure(wp, 5.0 + 0.109, 0.1);
    Scalar pOutside = vl->interpolateContactPressure(wp, 5.5, 0.1);

    BOOST_CHECK_CLOSE(pCenter, 700000.0, 1e-3);
    BOOST_CHECK_CLOSE(pInside, 700000.0, 1e-3);
    BOOST_CHECK_GT(pEdge, 0.0);
    BOOST_CHECK_EQUAL(pOutside, 0.0);
}

BOOST_AUTO_TEST_CASE(test_assemble_nodal_load_vector) {
    auto vm = VehicleLoad::createBZZ100(60.0, 5.0, 0.0);
    vl->addVehicle(vm);

    Index N = mesh.getNumNodes();
    VectorX F = vl->assembleNodalLoadVector(5.0);
    BOOST_CHECK_EQUAL(F.size(), 4 * N);

    Scalar sumFy = 0.0;
    for (Index i = 0; i < N; ++i) {
        sumFy += F(2 * N + 2 * i + 1);
    }

    BOOST_CHECK_LE(sumFy, 0.0);
    if (std::abs(sumFy) > 1e-10) {
        BOOST_CHECK_LT(sumFy, 0.0);
    }
}

BOOST_AUTO_TEST_CASE(test_position_and_pressure_fields) {
    auto vm = VehicleLoad::createBZZ100(60.0, 5.0, 0.0);
    vl->addVehicle(vm);

    Index N = mesh.getNumNodes();
    VectorX pos = vl->computeVehiclePositionField(5.0);
    VectorX press = vl->computeWheelPressureField(5.0);

    BOOST_CHECK_EQUAL(pos.size(), N);
    BOOST_CHECK_EQUAL(press.size(), N);
    BOOST_CHECK_GE(pos.minCoeff(), 0.0);
    BOOST_CHECK_GE(press.minCoeff(), 0.0);
}

BOOST_AUTO_TEST_CASE(test_record_and_history) {
    auto vm = VehicleLoad::createBZZ100(60.0, 5.0, 0.0);
    vl->addVehicle(vm);

    BOOST_CHECK_EQUAL(vl->getStepHistory().size(), 0u);
    vl->recordStep(0, 2.0);
    vl->recordStep(1, 5.0);
    BOOST_CHECK_EQUAL(vl->getStepHistory().size(), 2u);

    const auto& rec = vl->getStepHistory()[1];
    BOOST_CHECK_EQUAL(rec.stepId, 1);
    BOOST_CHECK_CLOSE(rec.time_s, 5.0, 1e-6);
}

BOOST_AUTO_TEST_CASE(test_clear_all) {
    auto vm = VehicleLoad::createBZZ100(60.0, 5.0, 0.0);
    vl->addVehicle(vm);
    vl->recordStep(0, 2.0);

    BOOST_CHECK_GT(vl->getVehicleCount(), 0u);
    BOOST_CHECK_GT(vl->getStepHistory().size(), 0u);

    vl->clearAllVehicles();
    BOOST_CHECK_EQUAL(vl->getVehicleCount(), 0u);
    BOOST_CHECK_EQUAL(vl->getStepHistory().size(), 0u);
}

BOOST_AUTO_TEST_CASE(test_aadt_stream) {
    Index count = vl->addVehicleStreamByAADT(
        VehicleClass::HeavyTruckBZZ100,
        60.0, 5.0,
        20000, 0.0, 20.0, 10);
    BOOST_CHECK_EQUAL(count, 10u);
    BOOST_CHECK_EQUAL(vl->getVehicleCount(), 10u);

    auto vehicles = vl->getAllVehicles();
    for (const auto& v : vehicles) {
        BOOST_CHECK(v.entryTime_s >= 0.0);
    }
}

BOOST_AUTO_TEST_CASE(test_remove_vehicle) {
    vl->addVehicle(VehicleLoad::createLightCar(80.0, 3.5, 0.0));
    vl->addVehicle(VehicleLoad::createBZZ100(60.0, 5.0, 0.0));
    BOOST_CHECK_EQUAL(vl->getVehicleCount(), 2u);

    vl->removeVehicle(0);
    BOOST_CHECK_EQUAL(vl->getVehicleCount(), 1u);
    const auto& remaining = vl->getVehicle(0);
    BOOST_CHECK_EQUAL(remaining.vClass, VehicleClass::HeavyTruckBZZ100);
}

BOOST_AUTO_TEST_CASE(test_assembler_dynamic_modulus) {
    ThermoHydro th(mesh);
    Assembler asmbler(mesh, th);

    Scalar E0 = 30e6;
    Scalar Efrozen = asmbler.computeElementDynamicModulus(0, 268.0, 0.30, E0, 10.0);
    Scalar Ephase = asmbler.computeElementDynamicModulus(0, 273.05, 0.20, E0, 10.0);
    Scalar EThaw = asmbler.computeElementDynamicModulus(0, 275.0, 0.00, E0, 10.0);
    Scalar EHiFreq = asmbler.computeElementDynamicModulus(0, 270.0, 0.20, E0, 50.0);
    Scalar ELowFreq = asmbler.computeElementDynamicModulus(0, 270.0, 0.20, E0, 1.0);

    BOOST_CHECK_GE(Efrozen, EThaw);
    BOOST_CHECK_GE(EHiFreq, ELowFreq);
    BOOST_CHECK_GT(E0, 0);
    BOOST_CHECK_GT(Ephase, 0);
}

BOOST_AUTO_TEST_CASE(test_add_vehicle_nodal_load) {
    ThermoHydro th(mesh);
    Assembler asmbler(mesh, th);

    Index N = mesh.getNumNodes();
    Index dofs = 4 * N;
    VectorX F = VectorX::Zero(dofs);
    VectorX VL = VectorX::LinSpaced(dofs, 1.0, static_cast<Scalar>(dofs));

    VectorX Fbefore = F;
    asmbler.addVehicleNodalLoad(F, VL);
    VectorX diff = F - Fbefore - VL;
    BOOST_CHECK_LT(diff.norm(), 1e-10);

    VectorX bad = VectorX::Zero(N);
    VectorX F2 = F;
    asmbler.addVehicleNodalLoad(F2, bad);
    BOOST_CHECK_CLOSE((F2 - F).norm(), 0.0, 1e-10);
}

BOOST_AUTO_TEST_SUITE_END()
