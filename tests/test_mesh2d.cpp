#define BOOST_TEST_MODULE Mesh2DTest
#include <boost/test/included/unit_test.hpp>
#include "Mesh2D.h"
#include <fstream>
#include <vector>

using namespace RoadbedSim;

BOOST_AUTO_TEST_CASE(test_add_nodes_elements) {
    Mesh2D mesh;
    mesh.addNode(0.0, 0.0);
    mesh.addNode(1.0, 0.0);
    mesh.addNode(0.0, 1.0);
    mesh.addNode(1.0, 1.0);

    mesh.addElement({0, 1, 2}, 0);
    mesh.addElement({1, 3, 2}, 0);

    BOOST_CHECK_EQUAL(mesh.getNumNodes(), 4);
    BOOST_CHECK_EQUAL(mesh.getNumElements(), 2);
}

BOOST_AUTO_TEST_CASE(test_compute_area) {
    Mesh2D mesh;
    mesh.addNode(0.0, 0.0);
    mesh.addNode(1.0, 0.0);
    mesh.addNode(0.0, 1.0);
    mesh.addElement({0, 1, 2}, 0);
    mesh.computeElementAreas();

    Scalar area = mesh.getElement(0).area;
    BOOST_CHECK_CLOSE(area, 0.5, 1e-10);
    BOOST_CHECK_CLOSE(mesh.computeElementArea(0), 0.5, 1e-10);
}

BOOST_AUTO_TEST_CASE(test_bounding_box) {
    Mesh2D mesh;
    mesh.addNode(0.0, 0.0);
    mesh.addNode(2.0, 0.0);
    mesh.addNode(0.0, 1.5);
    mesh.addNode(2.0, 1.5);

    BOOST_CHECK_CLOSE(mesh.getBoundingBoxMinX(), 0.0, 1e-10);
    BOOST_CHECK_CLOSE(mesh.getBoundingBoxMaxX(), 2.0, 1e-10);
    BOOST_CHECK_CLOSE(mesh.getBoundingBoxMinY(), 0.0, 1e-10);
    BOOST_CHECK_CLOSE(mesh.getBoundingBoxMaxY(), 1.5, 1e-10);
}

BOOST_AUTO_TEST_CASE(test_boundary_detection) {
    Mesh2D mesh;
    for (int i = 0; i <= 4; ++i) {
        for (int j = 0; j <= 3; ++j) {
            mesh.addNode(static_cast<Scalar>(i), static_cast<Scalar>(j));
        }
    }

    auto top = mesh.findSurfaceNodes();
    auto bottom = mesh.findBottomNodes();
    BOOST_CHECK_EQUAL(top.size(), 5);
    BOOST_CHECK_EQUAL(bottom.size(), 5);
}

BOOST_AUTO_TEST_CASE(test_mesh_save_load) {
    const std::string testFile = "test_mesh_temp.msh";
    Mesh2D mesh1;
    mesh1.addNode(0, 0); mesh1.addNode(1, 0); mesh1.addNode(0, 1); mesh1.addNode(1, 1);
    mesh1.addElement({0, 1, 2}, 0);
    mesh1.addElement({1, 3, 2}, 1);

    ZoneProperties z0{};
    z0.waterContent = 0.25; z0.youngModulus = 30e6;
    ZoneProperties z1{};
    z1.waterContent = 0.30; z1.youngModulus = 50e6;
    mesh1.addZone(0, z0); mesh1.addZone(1, z1);

    BOOST_CHECK(mesh1.saveToFile(testFile));

    Mesh2D mesh2;
    BOOST_CHECK(mesh2.loadFromFile(testFile));
    BOOST_CHECK_EQUAL(mesh2.getNumNodes(), 4);
    BOOST_CHECK_EQUAL(mesh2.getNumElements(), 2);
    BOOST_CHECK(mesh2.getNumZones() >= 2);

    std::remove(testFile.c_str());
}

BOOST_AUTO_TEST_CASE(test_node_neighbors) {
    Mesh2D mesh;
    mesh.addNode(0, 0); mesh.addNode(1, 0); mesh.addNode(0, 1); mesh.addNode(1, 1);
    mesh.addElement({0, 1, 2}, 0);
    mesh.addElement({1, 3, 2}, 0);
    mesh.computeNodeNeighbors();

    const auto& nn = mesh.getNodeNeighbors();
    BOOST_CHECK(!nn.empty());
    BOOST_CHECK(nn[0].size() >= 2);
}

BOOST_AUTO_TEST_CASE(test_zone_by_range) {
    Mesh2D mesh;
    mesh.addNode(0, 0); mesh.addNode(1, 0); mesh.addNode(0, 1); mesh.addNode(1, 1);
    mesh.addElement({0, 1, 2}, 0);
    mesh.addElement({1, 3, 2}, 0);
    mesh.computeElementAreas();

    ZoneProperties upper{};
    upper.waterContent = 0.5;
    mesh.setZonePropertyByRange(-1, 2, 0.4, 2, 1, upper);

    BOOST_CHECK_EQUAL(mesh.getElement(1).zoneId, 1);
}
