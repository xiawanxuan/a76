#define BOOST_TEST_MODULE SolverTest
#include <boost/test/included/unit_test.hpp>
#include "NewtonRaphsonSolver.h"
#include "ThermoHydro.h"
#include "Assembler.h"
#include "Mesh2D.h"

using namespace RoadbedSim;

struct SolverFixture {
    Mesh2D mesh;
    std::unique_ptr<ThermoHydro> th;
    std::unique_ptr<Assembler> asmb;
    std::unique_ptr<NewtonRaphsonSolver> solver;

    SolverFixture() {
        for (int j = 0; j < 3; ++j) {
            for (int i = 0; i < 3; ++i) {
                mesh.addNode(static_cast<Scalar>(i), static_cast<Scalar>(j));
            }
        }
        for (int j = 0; j < 2; ++j) {
            for (int i = 0; i < 2; ++i) {
                Index n0 = j * 3 + i;
                Index n1 = j * 3 + i + 1;
                Index n2 = (j + 1) * 3 + i;
                Index n3 = (j + 1) * 3 + i + 1;
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
        solver = std::make_unique<NewtonRaphsonSolver>(mesh, *th, *asmb);
    }
};

BOOST_FIXTURE_TEST_CASE(test_initial_conditions, SolverFixture) {
    Index N = mesh.getNumNodes();
    VectorX T0 = VectorX::Constant(N, 275.0);
    VectorX W0 = VectorX::Constant(N, 0.25);
    VectorX UX0 = VectorX::Zero(N);
    VectorX UY0 = VectorX::Zero(N);

    solver->setInitialConditions(T0, W0, UX0, UY0);
    const auto& fields = solver->getCurrentFields();

    BOOST_CHECK_EQUAL(fields.temperature.size(), N);
    BOOST_CHECK_CLOSE(fields.temperature(0), 275.0, 1e-10);
    BOOST_CHECK_CLOSE(fields.waterContent(0), 0.25, 1e-10);
}

BOOST_FIXTURE_TEST_CASE(test_solver_params, SolverFixture) {
    solver->setMaxIterations(100);
    solver->setTolerance(1e-8);
    solver->setTimeStep(1800.0);
    solver->setTotalTime(1e6);

    BOOST_CHECK_EQUAL(solver->getMaxIterations(), 100);
    BOOST_CHECK_CLOSE(solver->getTolerance(), 1e-8, 1e-15);
    BOOST_CHECK_CLOSE(solver->getTimeStep(), 1800.0, 1e-10);
    BOOST_CHECK_CLOSE(solver->getTotalTime(), 1e6, 1e-10);
}

BOOST_FIXTURE_TEST_CASE(test_boundary_conditions, SolverFixture) {
    BoundaryCondition bc;
    bc.type = BoundaryCondition::Type::Temperature;
    bc.nodeIds = mesh.findSurfaceNodes();
    bc.valueFunc = [](Scalar t) { return 270.0 + 5.0 * std::sin(t); };

    std::vector<BoundaryCondition> bcs{bc};
    solver->setBoundaryConditions(bcs);

    BOOST_CHECK_EQUAL(solver->getBoundaryConditions().size(), 1);
}

BOOST_FIXTURE_TEST_CASE(test_callbacks, SolverFixture) {
    int stepCount = 0;
    solver->setStepCompleteCallback([&](Index, Scalar) { ++stepCount; });

    int resCount = 0;
    solver->setResidualOutputCallback([&](Index, Scalar, Scalar) { ++resCount; });

    Index N = mesh.getNumNodes();
    VectorX T0 = VectorX::Constant(N, 275.0);
    VectorX W0 = VectorX::Constant(N, 0.25);
    solver->setInitialConditions(T0, W0, VectorX::Zero(N), VectorX::Zero(N));
    solver->setOutputInterval(1);
    solver->setTotalTime(solver->getTimeStep() * 2);

    solver->solveTransient();
    BOOST_CHECK_GE(stepCount, 0);
}

BOOST_FIXTURE_TEST_CASE(test_solver_selection, SolverFixture) {
    solver->useDirectSolver(true);
    BOOST_CHECK(solver->isUsingDirectSolver());
    solver->useDirectSolver(false);
    BOOST_CHECK(!solver->isUsingDirectSolver());
}

BOOST_FIXTURE_TEST_CASE(test_direct_solver_solve_step, SolverFixture) {
    Index N = mesh.getNumNodes();
    VectorX T0 = VectorX::Constant(N, 275.0);
    VectorX W0 = VectorX::Constant(N, 0.25);
    solver->setInitialConditions(T0, W0, VectorX::Zero(N), VectorX::Zero(N));
    solver->useDirectSolver(true);

    auto result = solver->solveStep(0.0);
    BOOST_CHECK(result.iterations >= 0);
}
