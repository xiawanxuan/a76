#include "BlockLoader.h"

#ifdef USE_BLOCK_LOADING

namespace RoadbedSim {

BlockLoader::BlockLoader(const std::string& meshFile, Index blockSize)
    : meshFile_(meshFile), blockSize_(blockSize) {}

bool BlockLoader::requiresBlockLoading(Index numNodes, Index numElements,
                                        Index memoryLimitMB) {
    const Index bytesPerNode = 3 * sizeof(Scalar) + sizeof(Index);
    const Index bytesPerElem = 5 * sizeof(Index) + sizeof(Scalar);
    const Index estMatrixBytes = (numNodes * 4) * (numNodes * 4) * sizeof(Scalar) / 100;

    Index totalMB = (numNodes * bytesPerNode + numElements * bytesPerElem + estMatrixBytes)
                    / (1024 * 1024);
    return totalMB > memoryLimitMB;
}

bool BlockLoader::openMeshFile() {
    file_.open(meshFile_, std::ios::in);
    if (!file_.is_open()) return false;
    parseHeader();
    computeBlockLayout();
    return true;
}

void BlockLoader::closeMeshFile() {
    if (file_.is_open()) file_.close();
}

void BlockLoader::parseHeader() {
    if (!file_.is_open()) return;
    file_.clear();
    file_.seekg(0);

    std::string section;
    while (file_ >> section) {
        if (section == "NODES") {
            file_ >> totalNodes_;
            nodesStartPos_ = file_.tellg();
            for (Index i = 0; i < totalNodes_; ++i) {
                Index id; Scalar x, y;
                file_ >> id >> x >> y;
            }
        } else if (section == "ELEMENTS") {
            file_ >> totalElements_;
            elementsStartPos_ = file_.tellg();
            break;
        }
    }
}

void BlockLoader::computeBlockLayout() {
    blocks_.clear();
    if (blockSize_ == 0) blockSize_ = 1000;

    Index nodeBlockSize = blockSize_;
    Index elemBlockSize = blockSize_ * 2;

    for (Index ns = 0; ns < totalNodes_; ns += nodeBlockSize) {
        BlockInfo bi;
        bi.nodeStart = ns;
        bi.nodeCount = std::min(nodeBlockSize, totalNodes_ - ns);
        bi.elemStart = 0;
        bi.elemCount = 0;
        blocks_.push_back(bi);
    }

    if (!blocks_.empty()) {
        blocks_[0].elemStart = 0;
        blocks_[0].elemCount = totalElements_;
    }
}

bool BlockLoader::loadBlockNodes(Index blockId, std::vector<Node2D>& nodes) {
    if (!file_.is_open() || blockId >= static_cast<Index>(blocks_.size()))
        return false;

    const auto& bi = blocks_[blockId];
    file_.clear();
    file_.seekg(nodesStartPos_);

    for (Index i = 0; i < bi.nodeStart; ++i) {
        Index id; Scalar x, y;
        if (!(file_ >> id >> x >> y)) return false;
    }

    nodes.clear();
    nodes.reserve(bi.nodeCount);
    for (Index i = 0; i < bi.nodeCount; ++i) {
        Index id; Scalar x, y;
        if (!(file_ >> id >> x >> y)) return false;
        nodes.push_back(Node2D{x, y, bi.nodeStart + i});
    }
    return true;
}

bool BlockLoader::loadBlockElements(Index blockId, std::vector<TriElement>& elems) {
    if (!file_.is_open() || blockId >= static_cast<Index>(blocks_.size()))
        return false;

    const auto& bi = blocks_[blockId];
    file_.clear();
    file_.seekg(elementsStartPos_);

    for (Index i = 0; i < bi.elemStart; ++i) {
        Index id, n0, n1, n2, z;
        if (!(file_ >> id >> n0 >> n1 >> n2 >> z)) return false;
    }

    elems.clear();
    elems.reserve(bi.elemCount);
    for (Index i = 0; i < bi.elemCount; ++i) {
        Index id, z;
        std::array<Index, 3> n;
        if (!(file_ >> id >> n[0] >> n[1] >> n[2] >> z)) return false;
        TriElement e;
        e.id = bi.elemStart + i;
        e.nodeIds = n;
        e.zoneId = z;
        e.area = 0.0;
        elems.push_back(e);
    }
    return true;
}

bool BlockLoader::loadFullMesh(Mesh2D& mesh) {
    if (!openMeshFile()) return false;
    closeMeshFile();
    return mesh.loadFromFile(meshFile_);
}

bool BlockLoader::writeBlockedMesh(const Mesh2D& mesh,
                                    const std::string& outputDir) {
    std::filesystem::create_directories(outputDir);
    return mesh.saveToFile(outputDir + "/roadbed_mesh.msh");
}

} // namespace RoadbedSim

#endif
