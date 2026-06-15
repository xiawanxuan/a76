#define BOOST_TEST_MODULE VTKWriterTest
#include <boost/test/included/unit_test.hpp>
#include "VTKWriter.h"
#include "Mesh2D.h"
#include <fstream>
#include <filesystem>
#include <limits>
#include <cctype>

using namespace RoadbedSim;

struct VTKWriterFixture {
    Mesh2D mesh;
    std::unique_ptr<VTKWriter> writer;
    std::string outDir = "./test_vtk_out";

    VTKWriterFixture() {
        for (int j = 0; j < 2; ++j) {
            for (int i = 0; i < 2; ++i) {
                mesh.addNode(static_cast<Scalar>(i), static_cast<Scalar>(j));
            }
        }
        mesh.addElement({0, 1, 3}, 0);
        mesh.addElement({0, 3, 2}, 0);
        mesh.computeElementAreas();
        mesh.computeNodeNeighbors();

        ZoneProperties z{};
        z.waterContent = 0.25; z.youngModulus = 30e6; z.porosity = 0.40;
        z.thermalCondFrozen = 2.2; z.thermalCondUnfrozen = 1.8;
        z.heatCapacityFrozen = 1800; z.heatCapacityUnfrozen = 2200;
        z.density = 1800; z.permeability = 1e-8; z.frostHeaveCoeff = 0.09;
        z.poissonRatio = 0.35;
        mesh.addZone(0, z);

        writer = std::make_unique<VTKWriter>(mesh);
        writer->setOutputDirectory(outDir);
        std::filesystem::create_directories(outDir);
    }

    ~VTKWriterFixture() {
        try {
            std::filesystem::remove_all(outDir);
        } catch (...) {}
    }
};

BOOST_FIXTURE_TEST_CASE(test_write_mesh, VTKWriterFixture) {
    writer->writeMeshOnly("test_mesh.vtk");
    std::ifstream f(outDir + "/test_mesh.vtk");
    BOOST_CHECK(f.good());
}

BOOST_FIXTURE_TEST_CASE(test_write_step, VTKWriterFixture) {
    FieldVariables fv;
    Index N = mesh.getNumNodes();
    fv.temperature = VectorX::Constant(N, 273.15);
    fv.waterContent = VectorX::Constant(N, 0.25);
    fv.displaceX = VectorX::Zero(N);
    fv.displaceY = VectorX::Constant(N, -0.001);
    fv.iceContent = VectorX::Constant(N, 0.05);
    fv.stressXX = VectorX::Zero(N);
    fv.stressYY = VectorX::Constant(N, -1e5);
    fv.stressXY = VectorX::Zero(N);

    std::vector<Scalar> damage(mesh.getNumElements(), 0.1);
    writer->writeStep(0, 3600.0, fv, damage);

    std::ifstream f(outDir + "/roadbed_000000.vtk");
    BOOST_CHECK(f.good());
}

BOOST_FIXTURE_TEST_CASE(test_write_convergence_log, VTKWriterFixture) {
    ConvergenceLog log;
    log.timeSteps = {0, 1, 2};
    log.iterationCounts = {3, 2, 5};
    log.residuals = {1e-1, 1e-3, 1e-5, 1e-7, 1e-2, 1e-8, 1e-2, 1e-4, 1e-6, 1e-9};

    writer->writeConvergenceLog(log, "conv_log.txt");
    std::ifstream f(outDir + "/conv_log.txt");
    BOOST_CHECK(f.good());
}

BOOST_FIXTURE_TEST_CASE(test_write_settlement, VTKWriterFixture) {
    std::vector<Scalar> times = {0, 86400, 86400 * 2};
    std::vector<Scalar> sett = {0, 0.005, 0.012};

    writer->writeSettlementHistory(times, sett, "sett.txt");
    std::ifstream f(outDir + "/sett.txt");
    BOOST_CHECK(f.good());
}

BOOST_FIXTURE_TEST_CASE(test_filename_format, VTKWriterFixture) {
    auto name = VTKWriter::formatStepFilename("prefix", 42, 5);
    BOOST_CHECK_EQUAL(name, "prefix_00042.vtk");
}

BOOST_FIXTURE_TEST_CASE(test_precision, VTKWriterFixture) {
    writer->setPrecision(12);
    BOOST_CHECK_EQUAL(writer->getPrecision(), 12);
    writer->setOutputDirectory("my_dir");
    BOOST_CHECK_EQUAL(writer->getOutputDirectory(), "my_dir");
}

BOOST_FIXTURE_TEST_CASE(test_frozen_soil_displacement_serialization, VTKWriterFixture) {
    FieldVariables fv;
    Index N = mesh.getNumNodes();
    fv.temperature = VectorX::Constant(N, 260.0);
    fv.waterContent = VectorX::Constant(N, 0.30);
    fv.displaceX = VectorX::Constant(N, 1e-10);
    fv.displaceY = VectorX::Constant(N, -5e-7);
    fv.iceContent = VectorX::Constant(N, 0.15);
    fv.stressXX = VectorX::Constant(N, -2.5e5);
    fv.stressYY = VectorX::Constant(N, -5.0e5);
    fv.stressXY = VectorX::Constant(N, 1.0e4);

    fv.displaceX(1) = std::numeric_limits<Scalar>::quiet_NaN();
    fv.displaceY(2) = std::numeric_limits<Scalar>::infinity();
    fv.temperature(0) = -std::numeric_limits<Scalar>::infinity();

    writer->writeStep(999, 86400.0 * 30, fv, {});

    std::ifstream f(outDir + "/roadbed_000999.vtk");
    BOOST_CHECK(f.good());

    std::string line;
    int displacementLines = 0;
    bool inDisplacement = false;
    while (std::getline(f, line)) {
        if (line.find("VECTORS Displacement_m") != std::string::npos) {
            inDisplacement = true;
            continue;
        }
        if (inDisplacement && !line.empty() && isdigit(line[0])) {
            displacementLines++;
            BOOST_CHECK(line.find("nan") == std::string::npos);
            BOOST_CHECK(line.find("inf") == std::string::npos);
            BOOST_CHECK(line.find("-nan") == std::string::npos);
            BOOST_CHECK(line.find("-inf") == std::string::npos);
            if (displacementLines >= N) break;
        } else if (inDisplacement && line.find("SCALARS") != std::string::npos) {
            break;
        }
    }
    BOOST_CHECK_EQUAL(displacementLines, N);
}

BOOST_FIXTURE_TEST_CASE(test_all_nodes_written_regardless_of_field_size, VTKWriterFixture) {
    FieldVariables fv;
    Index N = mesh.getNumNodes();
    fv.temperature.resize(N - 1);
    fv.waterContent.resize(0);
    fv.displaceX.resize(1);
    fv.displaceY.resize(N + 10);
    fv.iceContent.resize(N);
    fv.stressXX.resize(N);
    fv.stressYY.resize(N);
    fv.stressXY.resize(N);

    fv.temperature.setConstant(265.0);
    fv.displaceX.setConstant(1e-8);
    fv.displaceY.setConstant(-2e-6);

    writer->writeStep(100, 1000.0, fv, {});

    std::ifstream f(outDir + "/roadbed_000100.vtk");
    BOOST_CHECK(f.good());

    std::string line;
    int tempLines = 0, dispLines = 0;
    bool inTemp = false, inDisp = false;
    while (std::getline(f, line)) {
        if (line.find("SCALARS Temperature_K") != std::string::npos) {
            inTemp = true; inDisp = false;
            std::getline(f, line);
            continue;
        }
        if (line.find("VECTORS Displacement_m") != std::string::npos) {
            inDisp = true; inTemp = false;
            continue;
        }
        if (inTemp && !line.empty() && (isdigit(line[0]) || line[0] == '-')) {
            tempLines++;
            if (tempLines >= N) inTemp = false;
        }
        if (inDisp && !line.empty() && (isdigit(line[0]) || line[0] == '-')) {
            dispLines++;
            if (dispLines >= N) inDisp = false;
        }
    }
    BOOST_CHECK_EQUAL(tempLines, N);
    BOOST_CHECK_EQUAL(dispLines, N);
}
