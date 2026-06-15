#pragma once

#include "Types.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <cmath>
#include <algorithm>

namespace RoadbedSim {

class Mesh2D {
public:
    Mesh2D() = default;
    ~Mesh2D() = default;

    bool loadFromFile(const std::string& filename);
    bool saveToFile(const std::string& filename) const;

    void addNode(Scalar x, Scalar y);
    void addElement(const std::array<Index, 3>& nodeIds, Index zoneId = 0);
    void addZone(Index zoneId, const ZoneProperties& props);

    void computeElementAreas();
    void computeNodeNeighbors();

    Index getNumNodes() const { return static_cast<Index>(nodes_.size()); }
    Index getNumElements() const { return static_cast<Index>(elements_.size()); }
    Index getNumZones() const { return static_cast<Index>(zones_.size()); }

    const Node2D& getNode(Index id) const { return nodes_[id]; }
    const TriElement& getElement(Index id) const { return elements_[id]; }
    const ZoneProperties& getZone(Index zoneId) const;
    const ZoneProperties& getElementZone(Index elemId) const;

    const std::vector<Node2D>& getNodes() const { return nodes_; }
    const std::vector<TriElement>& getElements() const { return elements_; }
    const std::map<Index, ZoneProperties>& getZones() const { return zones_; }

    const std::vector<std::vector<Index>>& getNodeElements() const { return nodeElements_; }
    const std::vector<std::vector<Index>>& getNodeNeighbors() const { return nodeNeighbors_; }

    std::array<Scalar, 3> computeShapeFunctionDerivatives(Index elemId) const;
    MatrixX computeJacobian(Index elemId) const;
    Scalar computeElementArea(Index elemId) const;

    std::vector<Index> findSurfaceNodes(Scalar tolerance = 1e-8) const;
    std::vector<Index> findBottomNodes(Scalar tolerance = 1e-8) const;
    std::vector<Index> findLeftBoundaryNodes(Scalar tolerance = 1e-8) const;
    std::vector<Index> findRightBoundaryNodes(Scalar tolerance = 1e-8) const;

    Scalar getBoundingBoxMinX() const { return bboxMinX_; }
    Scalar getBoundingBoxMaxX() const { return bboxMaxX_; }
    Scalar getBoundingBoxMinY() const { return bboxMinY_; }
    Scalar getBoundingBoxMaxY() const { return bboxMaxY_; }

    void setZonePropertyByRange(Scalar xMin, Scalar xMax, Scalar yMin, Scalar yMax,
                                Index zoneId, const ZoneProperties& props);

    void setNodes(const std::vector<Node2D>& nodes) { nodes_ = nodes; updateBoundingBox(); }
    void setElements(const std::vector<TriElement>& elements) { elements_ = elements; }
    void recomputeAll() {
        computeElementAreas();
        computeNodeNeighbors();
        updateBoundingBox();
    }

private:
    void updateBoundingBox();
    void parseMeshNode(std::istream& is);
    void parseMeshElement(std::istream& is);
    void parseMeshZone(std::istream& is);

    std::vector<Node2D> nodes_;
    std::vector<TriElement> elements_;
    std::map<Index, ZoneProperties> zones_;
    std::vector<std::vector<Index>> nodeElements_;
    std::vector<std::vector<Index>> nodeNeighbors_;

    Scalar bboxMinX_ = 0, bboxMaxX_ = 0;
    Scalar bboxMinY_ = 0, bboxMaxY_ = 0;
};

} // namespace RoadbedSim
