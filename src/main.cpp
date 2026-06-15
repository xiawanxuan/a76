#include "Types.h"
#include "Mesh2D.h"
#include "ThermoHydro.h"
#include "Assembler.h"
#include "NewtonRaphsonSolver.h"
#include "VTKWriter.h"
#include "Config.h"
#include "BlockLoader.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <chrono>
#include <iomanip>
#include <exception>
#include <ctime>
#include <cmath>
#include <filesystem>

using namespace RoadbedSim;

std::unique_ptr<std::ofstream> g_logFile;

void logMessage(const std::string& msg) {
    auto now = std::chrono::system_clock::now();
    auto timeT = std::chrono::system_clock::to_time_t(now);
    std::ostringstream oss;
    oss << std::put_time(std::localtime(&timeT), "[%Y-%m-%d %H:%M:%S] ") << msg;
    std::cout << oss.str() << std::endl;
    if (g_logFile && g_logFile->is_open()) {
        *g_logFile << oss.str() << std::endl;
        g_logFile->flush();
    }
}

std::vector<BoundaryCondition> buildBoundaryConditions(
    const Mesh2D& mesh, const ThermoHydro& th, const SimulationConfig& cfg) {

    std::vector<BoundaryCondition> bcs;

    auto surfaceNodes = mesh.findSurfaceNodes();
    auto bottomNodes = mesh.findBottomNodes();
    auto leftNodes = mesh.findLeftBoundaryNodes();
    auto rightNodes = mesh.findRightBoundaryNodes();

    logMessage("表面节点数: " + std::to_string(surfaceNodes.size()) +
               ", 底部节点数: " + std::to_string(bottomNodes.size()) +
               ", 左侧: " + std::to_string(leftNodes.size()) +
               ", 右侧: " + std::to_string(rightNodes.size()));

    {
        BoundaryCondition bc;
        bc.type = BoundaryCondition::Type::Temperature;
        bc.nodeIds = surfaceNodes;
        Scalar meanT = cfg.meanAnnualTemp;
        Scalar amp = cfg.tempAmplitude;
        Scalar period = cfg.periodDays;
        Scalar phase = cfg.phaseShiftDays;
        bc.valueFunc = [=](Scalar t) -> Scalar {
            return th.annualAirTemperature(t, meanT, amp, period, phase);
        };
        bcs.push_back(bc);
    }

    {
        BoundaryCondition bc;
        bc.type = BoundaryCondition::Type::Temperature;
        bc.nodeIds = bottomNodes;
        Scalar botT = cfg.bottomBoundaryTemp;
        bc.valueFunc = [=](Scalar) -> Scalar { return botT; };
        bcs.push_back(bc);
    }

    {
        BoundaryCondition bc;
        bc.type = BoundaryCondition::Type::DisplacementX;
        bc.nodeIds = leftNodes;
        bc.valueFunc = [](Scalar) -> Scalar { return 0.0; };
        bcs.push_back(bc);
    }
    {
        BoundaryCondition bc;
        bc.type = BoundaryCondition::Type::DisplacementX;
        bc.nodeIds = rightNodes;
        bc.valueFunc = [](Scalar) -> Scalar { return 0.0; };
        bcs.push_back(bc);
    }
    {
        BoundaryCondition bc;
        bc.type = BoundaryCondition::Type::DisplacementY;
        bc.nodeIds = bottomNodes;
        bc.valueFunc = [](Scalar) -> Scalar { return 0.0; };
        bcs.push_back(bc);
    }

    return bcs;
}

void applyCustomZonesAndRanges(Mesh2D& mesh, const SimulationConfig& cfg) {
    for (const auto& [zid, props] : cfg.customZones) {
        mesh.addZone(zid, props);
    }
    for (const auto& [xMin, xMax, yMin, yMax, zid, props] : cfg.zoneRanges) {
        mesh.setZonePropertyByRange(xMin, xMax, yMin, yMax, zid, props);
    }
}

int main(int argc, char* argv[]) {
    try {
        Config config;
        if (!config.parseCommandLine(argc, argv)) {
            return 0;
        }

        const auto& cfg = config.get();

        g_logFile = std::make_unique<std::ofstream>(cfg.logFile);
        if (!g_logFile->is_open()) {
            std::cerr << "警告: 无法打开日志文件 " << cfg.logFile << std::endl;
        }

        config.printSummary();
        logMessage("启动寒区道路冻融耦合仿真");

        Mesh2D mesh;
        bool loaded = false;

#ifdef USE_BLOCK_LOADING
        if (cfg.enableBlockLoading) {
            BlockLoader loader(cfg.meshFile, cfg.blockSize);
            if (loader.openMeshFile()) {
                logMessage("检测分块加载: 节点=" + std::to_string(loader.getTotalNodes()) +
                           ", 单元=" + std::to_string(loader.getTotalElements()));
                if (BlockLoader::requiresBlockLoading(loader.getTotalNodes(),
                                                      loader.getTotalElements())) {
                    logMessage("内存受限 - 启用分块加载模式");
                }
                loaded = loader.loadFullMesh(mesh);
            }
        }
#endif

        if (!loaded) {
            loaded = mesh.loadFromFile(cfg.meshFile);
        }

        if (!loaded) {
            logMessage("错误: 无法加载网格文件 " + cfg.meshFile);
            return 1;
        }

        logMessage("网格加载成功: 节点=" + std::to_string(mesh.getNumNodes()) +
                   ", 单元=" + std::to_string(mesh.getNumElements()) +
                   ", 分区=" + std::to_string(mesh.getNumZones()));

        applyCustomZonesAndRanges(mesh, cfg);

        const Scalar bboxW = mesh.getBoundingBoxMaxX() - mesh.getBoundingBoxMinX();
        const Scalar bboxH = mesh.getBoundingBoxMaxY() - mesh.getBoundingBoxMinY();
        logMessage("几何范围: 宽度=" + std::to_string(bboxW) + "m, 高度="
                   + std::to_string(bboxH) + "m");

        ThermoHydro thermoHydro(mesh);
        Assembler assembler(mesh, thermoHydro);
        assembler.setBlockSize(cfg.blockSize);
        NewtonRaphsonSolver solver(mesh, thermoHydro, assembler);

        solver.setMaxIterations(cfg.maxIterations);
        solver.setTolerance(cfg.tolerance);
        solver.setTimeStep(cfg.timeStep);
        solver.setTotalTime(cfg.totalTime);
        solver.setOutputInterval(cfg.outputInterval);
        solver.useDirectSolver(cfg.useDirectSolver);

        auto bcs = buildBoundaryConditions(mesh, thermoHydro, cfg);
        solver.setBoundaryConditions(bcs);
        logMessage("边界条件设置完成: 共 " + std::to_string(bcs.size()) + " 组");

        VectorX T0 = thermoHydro.computeInitialTemperatureField(
            cfg.initialSurfaceTemp, cfg.geothermalGradient);
        VectorX W0 = thermoHydro.computeInitialWaterContentField();
        VectorX UX0 = VectorX::Zero(mesh.getNumNodes());
        VectorX UY0 = VectorX::Zero(mesh.getNumNodes());

        solver.setInitialConditions(T0, W0, UX0, UY0);
        logMessage("初始条件设置完成");

        VTKWriter writer(mesh);
        writer.setOutputDirectory(cfg.outputDir);
        writer.writeMeshOnly("roadbed_mesh.vtk");
        logMessage("初始网格VTK已导出");

        Index outputStep = 0;
        Index totalSteps = static_cast<Index>(std::ceil(cfg.totalTime / cfg.timeStep));
        auto globalStart = std::chrono::high_resolution_clock::now();

        solver.setResidualOutputCallback([](Index it, Scalar res, Scalar t) {
            if (it % 5 == 0 || it < 3) {
                std::ostringstream oss;
                oss << "  NR迭代 it=" << it
                    << " | 残差=" << std::scientific << res
                    << " | t=" << std::fixed << t << "s";
                logMessage(oss.str());
            }
        });

        solver.setStepCompleteCallback([&](Index step, Scalar t) {
            if ((step + 1) % cfg.outputInterval == 0 || step == 0) {
                const auto& fields = solver.getCurrentFields();
                const auto& damageHist = solver.getDamageHistory();

                std::vector<Scalar> damage;
                if (!damageHist.empty()) {
                    damage = damageHist.back();
                }

                writer.writeStep(outputStep, t, fields, damage);

                Scalar settle = solver.getSettlementHistory().empty()
                    ? 0.0 : solver.getSettlementHistory().back();

                std::ostringstream oss;
                oss << "【步骤 " << step + 1 << "/" << totalSteps << "】"
                    << " 时间=" << std::fixed << std::setprecision(1) << t / 3600.0 << "h"
                    << " 融沉=" << std::scientific << settle << "m"
                    << " 输出ID=" << outputStep;
                logMessage(oss.str());
                ++outputStep;
            }

            if ((step + 1) % (totalSteps / 100 + 1) == 0) {
                auto now = std::chrono::high_resolution_clock::now();
                auto elapsed = std::chrono::duration<Scalar>(now - globalStart).count();
                Scalar pct = static_cast<Scalar>(step + 1) / totalSteps * 100.0;
                std::ostringstream oss;
                oss << "进度: " << std::fixed << std::setprecision(1) << pct << "%"
                    << " | 已用: " << std::setprecision(0) << elapsed << "s";
                logMessage(oss.str());
            }
        });

        logMessage("开始瞬态耦合求解...");
        solver.solveTransient();
        logMessage("瞬态求解完成");

        writer.writeConvergenceLog(solver.getConvergenceLog(), "convergence_log.txt");
        writer.writeSettlementHistory(solver.getTimeHistory(),
                                       solver.getSettlementHistory(),
                                       "settlement_history.txt");
        logMessage("收敛日志与沉降历史已导出");

        if (!solver.getSettlementHistory().empty()) {
            Scalar maxSettle = 0.0;
            Scalar avgSettle = 0.0;
            for (Scalar s : solver.getSettlementHistory()) {
                maxSettle = std::max(maxSettle, s);
                avgSettle += s;
            }
            avgSettle /= solver.getSettlementHistory().size();
            logMessage("最大竖向融沉位移: " + std::to_string(maxSettle) + " m");
            logMessage("平均竖向融沉位移: " + std::to_string(avgSettle) + " m");
        }

        if (!solver.getDamageHistory().empty()) {
            const auto& lastDamage = solver.getDamageHistory().back();
            Index highRisk = 0;
            for (Scalar d : lastDamage) {
                if (d >= cfg.damageThreshold) ++highRisk;
            }
            Scalar pct = static_cast<Scalar>(highRisk) / lastDamage.size() * 100.0;
            logMessage("冻融高风险单元数: " + std::to_string(highRisk)
                       + " (" + std::to_string(pct) + "%)");
        }

        auto globalEnd = std::chrono::high_resolution_clock::now();
        Scalar totalElapsed = std::chrono::duration<Scalar>(globalEnd - globalStart).count();
        logMessage("总计算时间: " + std::to_string(totalElapsed) + " 秒");
        logMessage("仿真全部完成, 结果保存在: " + cfg.outputDir);

        if (g_logFile) g_logFile->close();
        return 0;

    } catch (const std::exception& e) {
        logMessage("致命错误: " + std::string(e.what()));
        if (g_logFile) g_logFile->close();
        return 1;
    } catch (...) {
        logMessage("致命错误: 未知异常");
        if (g_logFile) g_logFile->close();
        return 1;
    }
}
