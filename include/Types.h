#pragma once

#if !defined(_USE_MATH_DEFINES) && (defined(_WIN32) || defined(_WIN64))
#define _USE_MATH_DEFINES
#endif

#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <array>
#include <string>
#include <memory>
#include <map>
#include <cstdint>
#include <functional>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace RoadbedSim {

using Scalar = double;
using Index = Eigen::Index;
using VectorX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
using MatrixX = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
using SparseMatrixX = Eigen::SparseMatrix<Scalar>;
using Triplet = Eigen::Triplet<Scalar>;

struct Node2D {
    Scalar x;
    Scalar y;
    Index id;
};

struct TriElement {
    std::array<Index, 3> nodeIds;
    Index id;
    Index zoneId;
    Scalar area;
};

struct ZoneProperties {
    Scalar waterContent;
    Scalar frostHeaveCoeff;
    Scalar permeability;
    Scalar thermalCondFrozen;
    Scalar thermalCondUnfrozen;
    Scalar heatCapacityFrozen;
    Scalar heatCapacityUnfrozen;
    Scalar density;
    Scalar youngModulus;
    Scalar poissonRatio;
    Scalar porosity;
};

struct FieldVariables {
    VectorX temperature;
    VectorX waterContent;
    VectorX displaceX;
    VectorX displaceY;
    VectorX iceContent;
    VectorX stressXX;
    VectorX stressYY;
    VectorX stressXY;
};

struct BoundaryCondition {
    enum class Type { Temperature, WaterFlux, DisplacementX, DisplacementY, Stress };
    Type type;
    std::vector<Index> nodeIds;
    std::function<Scalar(Scalar time)> valueFunc;
};

struct SolverResult {
    bool converged;
    Index iterations;
    Scalar finalResidual;
    Scalar solveTime;
};

struct ConvergenceLog {
    std::vector<Scalar> residuals;
    std::vector<Index> iterationCounts;
    std::vector<Scalar> timeSteps;
};

inline constexpr Scalar FREEZING_POINT = 273.15;
inline constexpr Scalar GRAVITY = 9.81;
inline constexpr Scalar ICE_DENSITY = 917.0;
inline constexpr Scalar WATER_DENSITY = 1000.0;
inline constexpr Scalar LATENT_HEAT = 334000.0;

} // namespace RoadbedSim
