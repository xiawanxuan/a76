#include "VTKWriter.h"

namespace RoadbedSim {

namespace {
inline Scalar sanitizeScalar(Scalar val, Scalar fallback = 0.0) {
    if (!std::isfinite(val)) {
        return fallback;
    }
    constexpr Scalar MIN_ABS_DISP = 1e-12;
    if (std::abs(val) < MIN_ABS_DISP && val != 0.0) {
        return 0.0;
    }
    return val;
}
}

VTKWriter::VTKWriter(const Mesh2D& mesh) : mesh_(mesh) {
    ensureOutputDirectory();
}

void VTKWriter::ensureOutputDirectory() {
    std::filesystem::path p(outputDir_);
    if (!std::filesystem::exists(p)) {
        std::filesystem::create_directories(p);
    }
}

std::string VTKWriter::formatStepFilename(const std::string& prefix,
                                           Index stepIndex,
                                           Index totalDigits) {
    std::ostringstream oss;
    oss << prefix << "_" << std::setw(totalDigits) << std::setfill('0')
        << stepIndex << ".vtk";
    return oss.str();
}

void VTKWriter::setFieldNames(const std::vector<std::string>& scalars,
                               const std::vector<std::string>& vectors) {
    scalarNames_ = scalars;
    vectorNames_ = vectors;
}

void VTKWriter::writeHeader(std::ofstream& ofs, Scalar time) {
    ofs << std::setprecision(precision_) << std::scientific;
    ofs << "# vtk DataFile Version 3.0\n";
    ofs << "Roadbed Freeze-Thaw Simulation, time = " << time << " s\n";
    ofs << "ASCII\n";
    ofs << "DATASET UNSTRUCTURED_GRID\n";
}

void VTKWriter::writePoints(std::ofstream& ofs) {
    const auto& nodes = mesh_.getNodes();
    ofs << "POINTS " << nodes.size() << " double\n";
    for (const auto& n : nodes) {
        Scalar sx = sanitizeScalar(n.x, 0.0);
        Scalar sy = sanitizeScalar(n.y, 0.0);
        ofs << std::setprecision(precision_) << std::scientific
            << sx << " " << sy << " " << 0.0 << "\n";
    }
}

void VTKWriter::writeCells(std::ofstream& ofs) {
    const auto& elems = mesh_.getElements();
    Index numElems = static_cast<Index>(elems.size());
    Index totalSize = numElems * 4;
    ofs << "CELLS " << numElems << " " << totalSize << "\n";
    for (const auto& e : elems) {
        ofs << "3 " << e.nodeIds[0] << " " << e.nodeIds[1] << " " << e.nodeIds[2] << "\n";
    }
    ofs << "CELL_TYPES " << numElems << "\n";
    for (Index i = 0; i < numElems; ++i) {
        ofs << "5\n";
    }
}

void VTKWriter::writeCellData(std::ofstream& ofs, const std::vector<Scalar>& damage) {
    Index numElems = mesh_.getNumElements();
    ofs << "CELL_DATA " << numElems << "\n";
    if (!damage.empty() && static_cast<Index>(damage.size()) == numElems) {
        ofs << "SCALARS FreezeThawDamage double 1\n";
        ofs << "LOOKUP_TABLE default\n";
        for (Scalar d : damage) {
            Scalar sd = sanitizeScalar(d, 0.0);
            ofs << std::setprecision(precision_) << std::scientific << sd << "\n";
        }
        ofs << "SCALARS ZoneID int 1\n";
        ofs << "LOOKUP_TABLE default\n";
        for (const auto& e : mesh_.getElements()) {
            ofs << e.zoneId << "\n";
        }
    } else {
        ofs << "SCALARS ZoneID int 1\n";
        ofs << "LOOKUP_TABLE default\n";
        for (const auto& e : mesh_.getElements()) {
            ofs << e.zoneId << "\n";
        }
    }
}

void VTKWriter::writeScalarField(std::ofstream& ofs, const std::string& name,
                                  const VectorX& field) {
    const Index numNodes = mesh_.getNumNodes();
    ofs << "SCALARS " << name << " double 1\n";
    ofs << "LOOKUP_TABLE default\n";
    for (Index i = 0; i < numNodes; ++i) {
        Scalar val = (i < field.size()) ? field(i) : 0.0;
        val = sanitizeScalar(val, 0.0);
        ofs << std::setprecision(precision_) << std::scientific << val << "\n";
    }
}

void VTKWriter::writeVectorField(std::ofstream& ofs, const std::string& name,
                                  const VectorX& fieldX, const VectorX& fieldY) {
    const Index numNodes = mesh_.getNumNodes();
    ofs << "VECTORS " << name << " double\n";
    for (Index i = 0; i < numNodes; ++i) {
        Scalar vx = (i < fieldX.size()) ? fieldX(i) : 0.0;
        Scalar vy = (i < fieldY.size()) ? fieldY(i) : 0.0;
        vx = sanitizeScalar(vx, 0.0);
        vy = sanitizeScalar(vy, 0.0);
        ofs << std::setprecision(precision_) << std::scientific
            << vx << " " << vy << " " << 0.0 << "\n";
    }
}

void VTKWriter::writePointData(std::ofstream& ofs, const FieldVariables& fields,
                                const VectorX& vehiclePos,
                                const VectorX& wheelPress) {
    const Index numNodes = mesh_.getNumNodes();
    ofs << "POINT_DATA " << numNodes << "\n";

    writeScalarField(ofs, "Temperature_K", fields.temperature);
    writeScalarField(ofs, "WaterContent", fields.waterContent);
    writeScalarField(ofs, "IceContent", fields.iceContent);
    writeVectorField(ofs, "Displacement_m", fields.displaceX, fields.displaceY);

    {
        VectorX settlement(numNodes);
        for (Index i = 0; i < numNodes; ++i) {
            Scalar dy = (i < fields.displaceY.size()) ? fields.displaceY(i) : 0.0;
            settlement(i) = sanitizeScalar(-dy, 0.0);
        }
        writeScalarField(ofs, "Settlement_m", settlement);
    }

    writeScalarField(ofs, "Stress_XX_Pa", fields.stressXX);
    writeScalarField(ofs, "Stress_YY_Pa", fields.stressYY);
    writeScalarField(ofs, "Stress_XY_Pa", fields.stressXY);

    {
        VectorX tempC(numNodes);
        for (Index i = 0; i < numNodes; ++i) {
            Scalar tk = (i < fields.temperature.size()) ? fields.temperature(i) : FREEZING_POINT;
            tempC(i) = sanitizeScalar(tk - FREEZING_POINT, 0.0);
        }
        writeScalarField(ofs, "Temperature_C", tempC);
    }

    {
        VectorX freezeState(numNodes);
        for (Index i = 0; i < numNodes; ++i) {
            Scalar tk = (i < fields.temperature.size()) ? fields.temperature(i) : FREEZING_POINT;
            freezeState(i) = (tk < FREEZING_POINT - 0.5) ? 1.0 : 0.0;
        }
        writeScalarField(ofs, "FrozenState", freezeState);
    }

    writeScalarField(ofs, "VehiclePosition", vehiclePos);
    writeScalarField(ofs, "WheelPressure_Pa", wheelPress);
}

void VTKWriter::writeStep(Index stepIndex, Scalar time,
                           const FieldVariables& fields,
                           const std::vector<Scalar>& elemDamage,
                           const VectorX& vehiclePositionField,
                           const VectorX& wheelPressureField) {
    ensureOutputDirectory();
    std::string filename = outputDir_ + "/" + formatStepFilename("roadbed", stepIndex);
    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        throw std::runtime_error("Cannot open VTK file: " + filename);
    }

    writeHeader(ofs, time);
    writePoints(ofs);
    writeCells(ofs);
    writeCellData(ofs, elemDamage);
    writePointData(ofs, fields, vehiclePositionField, wheelPressureField);

    ofs.close();
}

void VTKWriter::writeConvergenceLog(const ConvergenceLog& log,
                                     const std::string& filename) {
    std::ofstream ofs(outputDir_ + "/" + filename);
    ofs << std::setprecision(precision_) << std::scientific;
    ofs << "# Step Iteration Residual\n";
    Index idx = 0;
    for (Index i = 0; i < static_cast<Index>(log.timeSteps.size()); ++i) {
        Index iters = (i < static_cast<Index>(log.iterationCounts.size()))
                      ? log.iterationCounts[i] : 0;
        for (Index j = 0; j < iters && idx < static_cast<Index>(log.residuals.size()); ++j) {
            ofs << i << " " << j << " " << log.residuals[idx++] << "\n";
        }
    }
    ofs.close();
}

void VTKWriter::writeSettlementHistory(const std::vector<Scalar>& times,
                                        const std::vector<Scalar>& settlements,
                                        const std::string& filename,
                                        const std::vector<Scalar>& trafficForceHistory) {
    std::ofstream ofs(outputDir_ + "/" + filename);
    ofs << std::setprecision(precision_) << std::scientific;
    if (trafficForceHistory.empty()) {
        ofs << "# Time_s Settlement_m\n";
        Index N = std::min(times.size(), settlements.size());
        for (Index i = 0; i < N; ++i) {
            ofs << times[i] << " " << settlements[i] << "\n";
        }
    } else {
        ofs << "# Time_s Settlement_m TotalTrafficForce_N\n";
        Index N = std::min({times.size(), settlements.size(), trafficForceHistory.size()});
        for (Index i = 0; i < N; ++i) {
            ofs << times[i] << " " << settlements[i] << " " << trafficForceHistory[i] << "\n";
        }
    }
    ofs.close();
}

void VTKWriter::writeMeshOnly(const std::string& filename) {
    ensureOutputDirectory();
    std::ofstream ofs(outputDir_ + "/" + filename);
    if (!ofs.is_open()) {
        throw std::runtime_error("Cannot open VTK file: " + filename);
    }
    writeHeader(ofs, 0.0);
    writePoints(ofs);
    writeCells(ofs);
    writeCellData(ofs, {});
    FieldVariables emptyFields;
    emptyFields.temperature = VectorX::Zero(mesh_.getNumNodes());
    emptyFields.waterContent = VectorX::Zero(mesh_.getNumNodes());
    emptyFields.iceContent = VectorX::Zero(mesh_.getNumNodes());
    emptyFields.displaceX = VectorX::Zero(mesh_.getNumNodes());
    emptyFields.displaceY = VectorX::Zero(mesh_.getNumNodes());
    emptyFields.stressXX = VectorX::Zero(mesh_.getNumNodes());
    emptyFields.stressYY = VectorX::Zero(mesh_.getNumNodes());
    emptyFields.stressXY = VectorX::Zero(mesh_.getNumNodes());
    writePointData(ofs, emptyFields, VectorX(), VectorX());
    ofs.close();
}

} // namespace RoadbedSim
