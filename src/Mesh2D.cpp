#include "Mesh2D.h"

namespace RoadbedSim {

static ZoneProperties MakeDefaultZone() {
    ZoneProperties z{};
    z.waterContent = 0.25;
    z.frostHeaveCoeff = 0.09;
    z.permeability = 1e-8;
    z.thermalCondFrozen = 2.2;
    z.thermalCondUnfrozen = 1.8;
    z.heatCapacityFrozen = 1800.0;
    z.heatCapacityUnfrozen = 2200.0;
    z.density = 1800.0;
    z.youngModulus = 30e6;
    z.poissonRatio = 0.35;
    z.porosity = 0.40;
    return z;
}

const ZoneProperties& Mesh2D::getZone(Index zoneId) const {
    static ZoneProperties def = MakeDefaultZone();
    auto it = zones_.find(zoneId);
    if (it != zones_.end()) return it->second;
    if (!zones_.empty()) return zones_.begin()->second;
    return def;
}

const ZoneProperties& Mesh2D::getElementZone(Index elemId) const {
    return getZone(elements_[elemId].zoneId);
}

void Mesh2D::addNode(Scalar x, Scalar y) {
    Node2D node{x, y, static_cast<Index>(nodes_.size())};
    nodes_.push_back(node);
    updateBoundingBox();
}

void Mesh2D::addElement(const std::array<Index, 3>& nodeIds, Index zoneId) {
    TriElement elem;
    elem.nodeIds = nodeIds;
    elem.id = static_cast<Index>(elements_.size());
    elem.zoneId = zoneId;
    elem.area = 0.0;
    elements_.push_back(elem);
}

void Mesh2D::addZone(Index zoneId, const ZoneProperties& props) {
    zones_[zoneId] = props;
}

void Mesh2D::updateBoundingBox() {
    if (nodes_.empty()) return;
    bboxMinX_ = bboxMaxX_ = nodes_[0].x;
    bboxMinY_ = bboxMaxY_ = nodes_[0].y;
    for (const auto& n : nodes_) {
        bboxMinX_ = std::min(bboxMinX_, n.x);
        bboxMaxX_ = std::max(bboxMaxX_, n.x);
        bboxMinY_ = std::min(bboxMinY_, n.y);
        bboxMaxY_ = std::max(bboxMaxY_, n.y);
    }
}

Scalar Mesh2D::computeElementArea(Index elemId) const {
    const auto& e = elements_[elemId];
    const auto& n0 = nodes_[e.nodeIds[0]];
    const auto& n1 = nodes_[e.nodeIds[1]];
    const auto& n2 = nodes_[e.nodeIds[2]];
    return 0.5 * std::abs((n1.x - n0.x) * (n2.y - n0.y) -
                          (n2.x - n0.x) * (n1.y - n0.y));
}

void Mesh2D::computeElementAreas() {
    for (auto& e : elements_) {
        e.area = computeElementArea(e.id);
    }
}

void Mesh2D::computeNodeNeighbors() {
    nodeElements_.assign(nodes_.size(), {});
    nodeNeighbors_.assign(nodes_.size(), {});

    for (const auto& e : elements_) {
        for (int i = 0; i < 3; ++i) {
            nodeElements_[e.nodeIds[i]].push_back(e.id);
            for (int j = 0; j < 3; ++j) {
                if (i != j) {
                    Index nb = e.nodeIds[j];
                    auto& vec = nodeNeighbors_[e.nodeIds[i]];
                    if (std::find(vec.begin(), vec.end(), nb) == vec.end()) {
                        vec.push_back(nb);
                    }
                }
            }
        }
    }
}

std::array<Scalar, 3> Mesh2D::computeShapeFunctionDerivatives(Index elemId) const {
    const auto& e = elements_[elemId];
    const auto& n0 = nodes_[e.nodeIds[0]];
    const auto& n1 = nodes_[e.nodeIds[1]];
    const auto& n2 = nodes_[e.nodeIds[2]];
    Scalar twoA = 2.0 * e.area;
    if (twoA < 1e-20) twoA = 1e-20;
    return { (n1.y - n2.y) / twoA, (n2.y - n0.y) / twoA, (n0.y - n1.y) / twoA };
}

MatrixX Mesh2D::computeJacobian(Index elemId) const {
    MatrixX J(2, 2);
    const auto& e = elements_[elemId];
    const auto& n0 = nodes_[e.nodeIds[0]];
    const auto& n1 = nodes_[e.nodeIds[1]];
    const auto& n2 = nodes_[e.nodeIds[2]];
    J(0, 0) = n1.x - n0.x; J(0, 1) = n2.x - n0.x;
    J(1, 0) = n1.y - n0.y; J(1, 1) = n2.y - n0.y;
    return J;
}

std::vector<Index> Mesh2D::findSurfaceNodes(Scalar tolerance) const {
    std::vector<Index> result;
    for (const auto& n : nodes_) {
        if (std::abs(n.y - bboxMaxY_) < tolerance) {
            result.push_back(n.id);
        }
    }
    return result;
}

std::vector<Index> Mesh2D::findBottomNodes(Scalar tolerance) const {
    std::vector<Index> result;
    for (const auto& n : nodes_) {
        if (std::abs(n.y - bboxMinY_) < tolerance) {
            result.push_back(n.id);
        }
    }
    return result;
}

std::vector<Index> Mesh2D::findLeftBoundaryNodes(Scalar tolerance) const {
    std::vector<Index> result;
    for (const auto& n : nodes_) {
        if (std::abs(n.x - bboxMinX_) < tolerance) {
            result.push_back(n.id);
        }
    }
    return result;
}

std::vector<Index> Mesh2D::findRightBoundaryNodes(Scalar tolerance) const {
    std::vector<Index> result;
    for (const auto& n : nodes_) {
        if (std::abs(n.x - bboxMaxX_) < tolerance) {
            result.push_back(n.id);
        }
    }
    return result;
}

void Mesh2D::setZonePropertyByRange(Scalar xMin, Scalar xMax, Scalar yMin, Scalar yMax,
                                     Index zoneId, const ZoneProperties& props) {
    zones_[zoneId] = props;
    for (auto& e : elements_) {
        const auto& n0 = nodes_[e.nodeIds[0]];
        const auto& n1 = nodes_[e.nodeIds[1]];
        const auto& n2 = nodes_[e.nodeIds[2]];
        Scalar cx = (n0.x + n1.x + n2.x) / 3.0;
        Scalar cy = (n0.y + n1.y + n2.y) / 3.0;
        if (cx >= xMin && cx <= xMax && cy >= yMin && cy <= yMax) {
            e.zoneId = zoneId;
        }
    }
}

void Mesh2D::parseMeshNode(std::istream& is) {
    Index id; Scalar x, y;
    while (is >> id >> x >> y) {
        addNode(x, y);
    }
}

void Mesh2D::parseMeshElement(std::istream& is) {
    Index id, z; std::array<Index, 3> n;
    while (is >> id >> n[0] >> n[1] >> n[2] >> z) {
        addElement(n, z);
    }
}

void Mesh2D::parseMeshZone(std::istream& is) {
    Index id;
    ZoneProperties p{};
    while (is >> id >> p.waterContent >> p.frostHeaveCoeff >> p.permeability
              >> p.thermalCondFrozen >> p.thermalCondUnfrozen
              >> p.heatCapacityFrozen >> p.heatCapacityUnfrozen
              >> p.density >> p.youngModulus >> p.poissonRatio >> p.porosity) {
        addZone(id, p);
    }
}

bool Mesh2D::loadFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;

    nodes_.clear();
    elements_.clear();
    zones_.clear();

    std::string section;
    while (file >> section) {
        if (section == "NODES") {
            Index count; file >> count;
            nodes_.reserve(count);
            parseMeshNode(file);
        } else if (section == "ELEMENTS") {
            Index count; file >> count;
            elements_.reserve(count);
            parseMeshElement(file);
        } else if (section == "ZONES") {
            parseMeshZone(file);
        }
    }
    file.close();

    computeElementAreas();
    computeNodeNeighbors();
    return true;
}

bool Mesh2D::saveToFile(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) return false;

    file << "NODES " << nodes_.size() << "\n";
    for (const auto& n : nodes_) {
        file << n.id << " " << n.x << " " << n.y << "\n";
    }
    file << "ELEMENTS " << elements_.size() << "\n";
    for (const auto& e : elements_) {
        file << e.id << " " << e.nodeIds[0] << " " << e.nodeIds[1] << " "
             << e.nodeIds[2] << " " << e.zoneId << "\n";
    }
    if (!zones_.empty()) {
        file << "ZONES\n";
        for (const auto& [id, p] : zones_) {
            file << id << " "
                 << p.waterContent << " " << p.frostHeaveCoeff << " "
                 << p.permeability << " "
                 << p.thermalCondFrozen << " " << p.thermalCondUnfrozen << " "
                 << p.heatCapacityFrozen << " " << p.heatCapacityUnfrozen << " "
                 << p.density << " " << p.youngModulus << " "
                 << p.poissonRatio << " " << p.porosity << "\n";
        }
    }
    file.close();
    return true;
}

} // namespace RoadbedSim
