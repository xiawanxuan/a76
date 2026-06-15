#define BOOST_TEST_MODULE AssemblerTest
#include <boost/test/included/unit_test.hpp>
#include "Assembler.h"
#include "ThermoHydro.h"
#include "Mesh2D.h"

using namespace RoadbedSim;

struct AssemblerFixture {
    Mesh2D mesh;
    std::unique_ptr<ThermoHydro> th;
    std::unique_ptr<Assembler> asmb;

    AssemblerFixture() {
        for (int j = 0; j < 3; ++j) {
            for (int i = 0; i < 4; ++i) {
                mesh.addNode(static_cast<Scalar>(i), static_cast<Scalar>(j));
            }
        }
        for (int j = 0; j < 2; ++j) {
            for (int i = 0; i < 3; ++i) {
                Index n0 = j * 4 + i;
                Index n1 = j * 4 + i + 1;
                Index n2 = (j + 1) * 4 + i;
                Index n3 = (j + 1) * 4 + i + 1;
                mesh.addElement({n0, n1, n3}, 0);
                mesh.addElement({n0, n3, n2}, 0);
            }
        }
        mesh.computeElementAreas();
        mesh.computeNodeNeighbors();

        ZoneProperties z{};
        z.waterContent = 0.25; z.porosity = 0.40;
        z.youngModulus = 30e6; z.poissonRatio = 0.35;
        z.thermalCondFrozen = 2.2; z.thermalCondUnfrozen = 1.8;
        z.heatCapacityFrozen = 1800; z.heatCapacityUnfrozen = 2200;
        z.density = 1800; z.permeability = 1e-8; z.frostHeaveCoeff = 0.09;
        mesh.addZone(0, z);

        th = std::make_unique<ThermoHydro>(mesh);
        asmb = std::make_unique<Assembler>(mesh, *th);
    }
};

BOOST_FIXTURE_TEST_CASE(test_thermal_matrix_size, AssemblerFixture) {
    Index N = mesh.getNumNodes();
    VectorX T0 = VectorX::Constant(N, 275.0);
    VectorX W0 = VectorX::Constant(N, 0.25);
    VectorX U0 = VectorX::Zero(2 * N);

    SparseMatrixX K; VectorX F;
    asmb->assembleThermalSystem(T0, W0, U0, 3600.0, K, F);

    BOOST_CHECK_EQUAL(K.rows(), 4 * N);
    BOOST_CHECK_EQUAL(K.cols(), 4 * N);
    BOOST_CHECK_EQUAL(F.size(), 4 * N);
}

BOOST_FIXTURE_TEST_CASE(test_mechanical_matrix_size, AssemblerFixture) {
    Index N = mesh.getNumNodes();
    VectorX T0 = VectorX::Constant(N, 275.0);
    VectorX W0 = VectorX::Constant(N, 0.25);
    VectorX U0 = VectorX::Zero(2 * N);

    SparseMatrixX K; VectorX F;
    asmb->assembleMechanicalSystem(T0, W0, U0, 3600.0, K, F);

    BOOST_CHECK_EQUAL(K.rows(), 4 * N);
}

BOOST_FIXTURE_TEST_CASE(test_boundary_condition_apply, AssemblerFixture) {
    Index N = mesh.getNumNodes();
    VectorX T0 = VectorX::Constant(N, 275.0);
    VectorX W0 = VectorX::Constant(N, 0.25);
    VectorX U0 = VectorX::Zero(2 * N);

    SparseMatrixX K; VectorX F;
    asmb->assembleThermalSystem(T0, W0, U0, 3600.0, K, F);

    BoundaryCondition bc;
    bc.type = BoundaryCondition::Type::Temperature;
    auto surf = mesh.findSurfaceNodes();
    bc.nodeIds = surf;
    bc.valueFunc = [](Scalar) { return 280.0; };

    std::vector<BoundaryCondition> bcs{bc};
    asmb->applyBoundaryConditions(K, F, bcs, 0.0);

    for (Index nid : surf) {
        BOOST_CHECK_CLOSE(F(nid), 280.0, 1e-6);
    }
}

BOOST_FIXTURE_TEST_CASE(test_damage_computation, AssemblerFixture) {
    Index N = mesh.getNumNodes();
    Index nSteps = 10;
    VectorX hist = VectorX::Zero(nSteps * N);
    for (Index s = 0; s < nSteps; ++s) {
        for (Index i = 0; i < N; ++i) {
            Scalar cyclePos = std::sin(2 * M_PI * s / nSteps);
            hist(s * N + i) = 273.15 - 2.0 + cyclePos * 5.0;
        }
    }
    auto damage = asmb->computeElementFreezeThawDamage(hist, 0.2);
    BOOST_CHECK_EQUAL(damage.size(), mesh.getNumElements());
}

BOOST_FIXTURE_TEST_CASE(test_settlement_computation, AssemblerFixture) {
    Index N = mesh.getNumNodes();
    VectorX disp = VectorX::Zero(2 * N);
    for (Index i = 0; i < N; ++i) {
        disp(2 * i + 1) = -0.01;
    }
    VectorX settle = asmb->computeSettlementDisplacement(disp);
    BOOST_CHECK_EQUAL(settle.size(), N);
    BOOST_CHECK_CLOSE(settle(0), 0.01, 1e-10);
}

BOOST_FIXTURE_TEST_CASE(test_hydraulic_matrix, AssemblerFixture) {
    Index N = mesh.getNumNodes();
    VectorX T0 = VectorX::Constant(N, 275.0);
    VectorX W0 = VectorX::Constant(N, 0.25);
    VectorX U0 = VectorX::Zero(2 * N);

    SparseMatrixX K; VectorX F;
    asmb->assembleHydraulicSystem(T0, W0, U0, 3600.0, K, F);
    BOOST_CHECK_EQUAL(K.rows(), 4 * N);
    BOOST_CHECK_GT(K.nonZeros(), 0);
}
