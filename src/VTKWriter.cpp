#include "VTKWriter.h"

namespace RoadbedSim {

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
        ofs << n.x << " " << n.y << " 0.0\n";
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
            ofs << d << "\n";
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
    ofs << "SCALARS " << name << " double 1\n";
    ofs << "LOOKUP_TABLE default\n";
    for (Index i = 0; i < field.size(); ++i) {
        ofs << field(i) << "\n";
    }
}

void VTKWriter::writeVectorField(std::ofstream& ofs, const std::string& name,
                                  const VectorX& fieldX, const VectorX& fieldY) {
    ofs << "VECTORS " << name << " double\n";
    Index N = std::min(fieldX.size(), fieldY.size());
    for (Index i = 0; i < N; ++i) {
        ofs << fieldX(i) << " " << fieldY(i) << " 0.0\n";
    }
}

void VTKWriter::writePointData(std::ofstream& ofs, const FieldVariables& fields) {
    Index N = mesh_.getNumNodes();
    ofs << "POINT_DATA " << N << "\n";

    if (fields.temperature.size() == N)
        writeScalarField(ofs, "Temperature_K", fields.temperature);

    if (fields.waterContent.size() == N)
        writeScalarField(ofs, "WaterContent", fields.waterContent);

    if (fields.iceContent.size() == N)
        writeScalarField(ofs, "IceContent", fields.iceContent);

    if (fields.displaceX.size() == N && fields.displaceY.size() == N)
        writeVectorField(ofs, "Displacement_m", fields.displaceX, fields.displaceY);

    if (fields.displaceY.size() == N) {
        VectorX settlement = -fields.displaceY;
        writeScalarField(ofs, "Settlement_m", settlement);
    }

    if (fields.stressXX.size() == N)
        writeScalarField(ofs, "Stress_XX_Pa", fields.stressXX);
    if (fields.stressYY.size() == N)
        writeScalarField(ofs, "Stress_YY_Pa", fields.stressYY);
    if (fields.stressXY.size() == N)
        writeScalarField(ofs, "Stress_XY_Pa", fields.stressXY);

    if (fields.temperature.size() == N) {
        VectorX tempC = fields.temperature.array() - FREEZING_POINT;
        writeScalarField(ofs, "Temperature_C", tempC);
    }
}

void VTKWriter::writeStep(Index stepIndex, Scalar time,
                           const FieldVariables& fields,
                           const std::vector<Scalar>& elemDamage) {
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
    writePointData(ofs, fields);

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
                                        const std::string& filename) {
    std::ofstream ofs(outputDir_ + "/" + filename);
    ofs << std::setprecision(precision_) << std::scientific;
    ofs << "# Time_s Settlement_m\n";
    Index N = std::min(times.size(), settlements.size());
    for (Index i = 0; i < N; ++i) {
        ofs << times[i] << " " << settlements[i] << "\n";
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
    ofs.close();
}

} // namespace RoadbedSim
