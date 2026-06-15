#pragma once

#include "Types.h"
#include "Mesh2D.h"
#include <fstream>
#include <vector>
#include <string>
#include <memory>
#include <functional>
#include <filesystem>

#ifdef USE_BLOCK_LOADING

namespace RoadbedSim {

class BlockLoader {
public:
    explicit BlockLoader(const std::string& meshFile, Index blockSize = 1000);
    ~BlockLoader() = default;

    struct BlockInfo {
        Index nodeStart;
        Index nodeCount;
        Index elemStart;
        Index elemCount;
    };

    bool openMeshFile();
    void closeMeshFile();

    Index getTotalNodes() const { return totalNodes_; }
    Index getTotalElements() const { return totalElements_; }
    Index getNumBlocks() const { return static_cast<Index>(blocks_.size()); }
    Index getBlockSize() const { return blockSize_; }

    const std::vector<BlockInfo>& getBlockInfos() const { return blocks_; }

    bool loadBlockNodes(Index blockId, std::vector<Node2D>& nodes);
    bool loadBlockElements(Index blockId, std::vector<TriElement>& elems);

    bool loadFullMesh(Mesh2D& mesh);

    bool writeBlockedMesh(const Mesh2D& mesh,
                          const std::string& outputDir);

    void setProgressCallback(std::function<void(Index, Index)> cb) {
        progressCallback_ = std::move(cb);
    }

    static bool requiresBlockLoading(Index numNodes, Index numElements,
                                     Index memoryLimitMB = 2048);

private:
    void parseHeader();
    void computeBlockLayout();

    std::string meshFile_;
    mutable std::ifstream file_;
    Index blockSize_;
    Index totalNodes_ = 0;
    Index totalElements_ = 0;
    std::vector<BlockInfo> blocks_;

    std::streampos nodesStartPos_;
    std::streampos elementsStartPos_;

    std::function<void(Index, Index)> progressCallback_;
};

} // namespace RoadbedSim

#endif
