#include "Config.h"

namespace RoadbedSim {

Config::Config() {
    setupDescription();
    setupDefaultZones();
}

void Config::setupDefaultZones() {
    ZoneProperties base{};
    base.waterContent = 0.25;
    base.frostHeaveCoeff = 0.09;
    base.permeability = 1e-8;
    base.thermalCondFrozen = 2.2;
    base.thermalCondUnfrozen = 1.8;
    base.heatCapacityFrozen = 1800.0;
    base.heatCapacityUnfrozen = 2200.0;
    base.density = 1800.0;
    base.youngModulus = 30e6;
    base.poissonRatio = 0.35;
    base.porosity = 0.40;
    cfg_.customZones[0] = base;

    ZoneProperties subgrade = base;
    subgrade.waterContent = 0.30;
    subgrade.thermalCondFrozen = 2.0;
    subgrade.thermalCondUnfrozen = 1.6;
    subgrade.youngModulus = 50e6;
    subgrade.porosity = 0.38;
    cfg_.customZones[1] = subgrade;

    ZoneProperties embankment = base;
    embankment.waterContent = 0.20;
    embankment.thermalCondFrozen = 2.4;
    embankment.thermalCondUnfrozen = 2.0;
    embankment.youngModulus = 80e6;
    embankment.porosity = 0.35;
    cfg_.customZones[2] = embankment;
}

void Config::setupDescription() {
    using namespace boost::program_options;
    desc_.add_options()
        ("help,h", "显示帮助信息")
        ("mesh,m", value<std::string>(&cfg_.meshFile)->required(),
            "路基三角形网格文件路径")
        ("config,c", value<std::string>(&cfg_.configFile),
            "JSON配置文件路径")
        ("output,o", value<std::string>(&cfg_.outputDir)->default_value("./vtk_output"),
            "VTK输出目录")
        ("log,l", value<std::string>(&cfg_.logFile)->default_value("simulation.log"),
            "日志文件路径")
        ("mean-temp", value<Scalar>(&cfg_.meanAnnualTemp)->default_value(270.15),
            "年平均气温 (K)")
        ("amp-temp", value<Scalar>(&cfg_.tempAmplitude)->default_value(20.0),
            "气温年变幅 (K)")
        ("period", value<Scalar>(&cfg_.periodDays)->default_value(365.0),
            "冻融周期 (天)")
        ("phase", value<Scalar>(&cfg_.phaseShiftDays)->default_value(15.0),
            "相位偏移 (天)")
        ("surface-temp", value<Scalar>(&cfg_.initialSurfaceTemp)->default_value(278.15),
            "初始表面温度 (K)")
        ("geo-gradient", value<Scalar>(&cfg_.geothermalGradient)->default_value(0.025),
            "地温梯度 (K/m)")
        ("total-time", value<Scalar>(&cfg_.totalTime)->default_value(3.154e7),
            "总模拟时间 (秒)")
        ("time-step", value<Scalar>(&cfg_.timeStep)->default_value(3600.0),
            "时间步长 (秒)")
        ("output-interval", value<Index>(&cfg_.outputInterval)->default_value(24),
            "输出间隔 (步数)")
        ("max-iter", value<Index>(&cfg_.maxIterations)->default_value(50),
            "Newton迭代最大次数")
        ("tolerance", value<Scalar>(&cfg_.tolerance)->default_value(1e-6),
            "收敛容差")
        ("iterative-solver", value<bool>(&cfg_.useDirectSolver)->default_value(true),
            "使用直接求解器(true)或迭代求解器(false)")
        ("bottom-temp", value<Scalar>(&cfg_.bottomBoundaryTemp)->default_value(280.15),
            "底部温度边界 (K)")
        ("damage-thresh", value<Scalar>(&cfg_.damageThreshold)->default_value(0.5),
            "冻融损伤阈值")
        ("block-size", value<Index>(&cfg_.blockSize)->default_value(1000),
            "分块加载块大小")
        ("enable-block", value<bool>(&cfg_.enableBlockLoading)->default_value(true),
            "启用分块内存加载")
        ("summary,s", "打印配置摘要后退出")
    ;
}

std::string Config::getHelpMessage() const {
    std::ostringstream oss;
    oss << desc_;
    return oss.str();
}

bool Config::parseCommandLine(int argc, char* argv[]) {
    try {
        using namespace boost::program_options;
        store(parse_command_line(argc, argv, desc_), vm_);

        if (vm_.count("help")) {
            std::cout << getHelpMessage() << std::endl;
            return false;
        }

        if (!vm_["config"].as<std::string>().empty()) {
            if (!loadFromFile(vm_["config"].as<std::string>())) {
                std::cerr << "警告: 无法加载配置文件 "
                          << vm_["config"].as<std::string>() << std::endl;
            }
        }

        notify(vm_);

        if (vm_.count("summary")) {
            printSummary();
            return false;
        }

        validate();
        return true;
    } catch (const std::exception& e) {
        std::cerr << "命令行解析错误: " << e.what() << std::endl;
        std::cerr << getHelpMessage() << std::endl;
        return false;
    }
}

bool Config::loadFromFile(const std::string& filename) {
    try {
        boost::property_tree::ptree pt;
        boost::property_tree::read_json(filename, pt);

        cfg_.meshFile = pt.get("mesh_file", cfg_.meshFile);
        cfg_.outputDir = pt.get("output_dir", cfg_.outputDir);
        cfg_.logFile = pt.get("log_file", cfg_.logFile);

        cfg_.meanAnnualTemp = pt.get("mean_annual_temp_K", cfg_.meanAnnualTemp);
        cfg_.tempAmplitude = pt.get("temp_amplitude_K", cfg_.tempAmplitude);
        cfg_.periodDays = pt.get("period_days", cfg_.periodDays);
        cfg_.phaseShiftDays = pt.get("phase_shift_days", cfg_.phaseShiftDays);
        cfg_.initialSurfaceTemp = pt.get("initial_surface_temp_K", cfg_.initialSurfaceTemp);
        cfg_.geothermalGradient = pt.get("geothermal_gradient", cfg_.geothermalGradient);

        cfg_.totalTime = pt.get("total_time_s", cfg_.totalTime);
        cfg_.timeStep = pt.get("time_step_s", cfg_.timeStep);
        cfg_.outputInterval = pt.get("output_interval_steps", cfg_.outputInterval);
        cfg_.maxIterations = pt.get("max_iterations", cfg_.maxIterations);
        cfg_.tolerance = pt.get("tolerance", cfg_.tolerance);
        cfg_.useDirectSolver = pt.get("use_direct_solver", cfg_.useDirectSolver);
        cfg_.bottomBoundaryTemp = pt.get("bottom_temp_K", cfg_.bottomBoundaryTemp);
        cfg_.damageThreshold = pt.get("damage_threshold", cfg_.damageThreshold);
        cfg_.blockSize = pt.get("block_size", cfg_.blockSize);
        cfg_.enableBlockLoading = pt.get("enable_block_loading", cfg_.enableBlockLoading);

        if (auto zones = pt.get_child_optional("zones")) {
            for (const auto& [key, zoneNode] : *zones) {
                Index zid = std::stoi(key);
                ZoneProperties p{};
                p.waterContent = zoneNode.get("water_content", 0.25);
                p.frostHeaveCoeff = zoneNode.get("frost_heave_coeff", 0.09);
                p.permeability = zoneNode.get("permeability", 1e-8);
                p.thermalCondFrozen = zoneNode.get("thermal_cond_frozen", 2.2);
                p.thermalCondUnfrozen = zoneNode.get("thermal_cond_unfrozen", 1.8);
                p.heatCapacityFrozen = zoneNode.get("heat_cap_frozen", 1800.0);
                p.heatCapacityUnfrozen = zoneNode.get("heat_cap_unfrozen", 2200.0);
                p.density = zoneNode.get("density", 1800.0);
                p.youngModulus = zoneNode.get("young_modulus", 30e6);
                p.poissonRatio = zoneNode.get("poisson_ratio", 0.35);
                p.porosity = zoneNode.get("porosity", 0.40);
                cfg_.customZones[zid] = p;
            }
        }

        if (auto ranges = pt.get_child_optional("zone_ranges")) {
            for (const auto& [key, rNode] : *ranges) {
                Scalar xMin = rNode.get("x_min", 0.0);
                Scalar xMax = rNode.get("x_max", 0.0);
                Scalar yMin = rNode.get("y_min", 0.0);
                Scalar yMax = rNode.get("y_max", 0.0);
                Index zid = rNode.get("zone_id", 0);
                ZoneProperties p = cfg_.customZones.count(zid)
                    ? cfg_.customZones[zid] : ZoneProperties{};
                cfg_.zoneRanges.emplace_back(xMin, xMax, yMin, yMax, zid, p);
            }
        }

        return true;
    } catch (const std::exception& e) {
        std::cerr << "配置文件读取错误: " << e.what() << std::endl;
        return false;
    }
}

bool Config::saveToFile(const std::string& filename) const {
    try {
        boost::property_tree::ptree pt;
        pt.put("mesh_file", cfg_.meshFile);
        pt.put("output_dir", cfg_.outputDir);
        pt.put("log_file", cfg_.logFile);
        pt.put("mean_annual_temp_K", cfg_.meanAnnualTemp);
        pt.put("temp_amplitude_K", cfg_.tempAmplitude);
        pt.put("period_days", cfg_.periodDays);
        pt.put("phase_shift_days", cfg_.phaseShiftDays);
        pt.put("initial_surface_temp_K", cfg_.initialSurfaceTemp);
        pt.put("geothermal_gradient", cfg_.geothermalGradient);
        pt.put("total_time_s", cfg_.totalTime);
        pt.put("time_step_s", cfg_.timeStep);
        pt.put("output_interval_steps", cfg_.outputInterval);
        pt.put("max_iterations", cfg_.maxIterations);
        pt.put("tolerance", cfg_.tolerance);
        pt.put("use_direct_solver", cfg_.useDirectSolver);
        pt.put("bottom_temp_K", cfg_.bottomBoundaryTemp);
        pt.put("damage_threshold", cfg_.damageThreshold);
        pt.put("block_size", cfg_.blockSize);
        pt.put("enable_block_loading", cfg_.enableBlockLoading);

        boost::property_tree::ptree zonesNode;
        for (const auto& [zid, p] : cfg_.customZones) {
            boost::property_tree::ptree zNode;
            zNode.put("water_content", p.waterContent);
            zNode.put("frost_heave_coeff", p.frostHeaveCoeff);
            zNode.put("permeability", p.permeability);
            zNode.put("thermal_cond_frozen", p.thermalCondFrozen);
            zNode.put("thermal_cond_unfrozen", p.thermalCondUnfrozen);
            zNode.put("heat_cap_frozen", p.heatCapacityFrozen);
            zNode.put("heat_cap_unfrozen", p.heatCapacityUnfrozen);
            zNode.put("density", p.density);
            zNode.put("young_modulus", p.youngModulus);
            zNode.put("poisson_ratio", p.poissonRatio);
            zNode.put("porosity", p.porosity);
            zonesNode.put_child(std::to_string(zid), zNode);
        }
        pt.put_child("zones", zonesNode);

        boost::property_tree::write_json(filename, pt);
        return true;
    } catch (const std::exception& e) {
        std::cerr << "配置文件写入错误: " << e.what() << std::endl;
        return false;
    }
}

void Config::validate() const {
    if (cfg_.meshFile.empty()) {
        throw std::invalid_argument("必须指定网格文件路径");
    }
    if (cfg_.totalTime <= 0 || cfg_.timeStep <= 0) {
        throw std::invalid_argument("时间参数必须为正数");
    }
    if (cfg_.maxIterations <= 0) {
        throw std::invalid_argument("最大迭代次数必须为正整数");
    }
    if (cfg_.tolerance <= 0) {
        throw std::invalid_argument("收敛容差必须为正数");
    }
}

void Config::printSummary(std::ostream& os) const {
    os << "========== 寒区道路冻融仿真配置摘要 ==========\n";
    os << "网格文件:       " << cfg_.meshFile << "\n";
    os << "输出目录:       " << cfg_.outputDir << "\n";
    os << "日志文件:       " << cfg_.logFile << "\n";
    os << "----------------------------------------------\n";
    os << "年平均气温:     " << cfg_.meanAnnualTemp << " K\n";
    os << "气温年变幅:     " << cfg_.tempAmplitude << " K\n";
    os << "冻融周期:       " << cfg_.periodDays << " 天\n";
    os << "初始表面温度:   " << cfg_.initialSurfaceTemp << " K\n";
    os << "地温梯度:       " << cfg_.geothermalGradient << " K/m\n";
    os << "----------------------------------------------\n";
    os << "总模拟时间:     " << cfg_.totalTime << " s ("
       << cfg_.totalTime / 86400.0 << " 天)\n";
    os << "时间步长:       " << cfg_.timeStep << " s ("
       << cfg_.timeStep / 3600.0 << " 小时)\n";
    os << "总步数:         " << static_cast<Index>(cfg_.totalTime / cfg_.timeStep) << "\n";
    os << "输出间隔:       每 " << cfg_.outputInterval << " 步\n";
    os << "----------------------------------------------\n";
    os << "最大迭代次数:   " << cfg_.maxIterations << "\n";
    os << "收敛容差:       " << cfg_.tolerance << "\n";
    os << "求解器:         " << (cfg_.useDirectSolver ? "直接(SparseLU)" : "迭代(BiCGSTAB)") << "\n";
    os << "底部温度:       " << cfg_.bottomBoundaryTemp << " K\n";
    os << "损伤阈值:       " << cfg_.damageThreshold << "\n";
    os << "分块加载:       " << (cfg_.enableBlockLoading ? "启用" : "禁用") << "\n";
    os << "分区数量:       " << cfg_.customZones.size() << "\n";
    os << "==============================================\n";
}

} // namespace RoadbedSim
