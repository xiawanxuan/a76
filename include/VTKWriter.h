#pragma once

#include "Types.h"
#include "Mesh2D.h"
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include <vector>

namespace RoadbedSim {

class VTKWriter {
public:
    explicit VTKWriter(const Mesh2D& mesh);
    ~VTKWriter() = default;

    void setOutputDirectory(const std::string& dir) { outputDir_ = dir; }
    const std::string& getOutputDirectory() const { return outputDir_; }

    void setPrecision(int p) { precision_ = p; }
    int getPrecision() const { return precision_; }

    void writeStep(Index stepIndex, Scalar time,
                   const FieldVariables& fields,
                   const std::vector<Scalar>& elemDamage = {},
                   const VectorX& vehiclePositionField = VectorX(),
                   const VectorX& wheelPressureField = VectorX());

    void writeConvergenceLog(const ConvergenceLog& log,
                             const std::string& filename);

    void writeSettlementHistory(const std::vector<Scalar>& times,
                                const std::vector<Scalar>& settlements,
                                const std::string& filename,
                                const std::vector<Scalar>& trafficForceHistory = {});

    void writeMeshOnly(const std::string& filename);

    static std::string formatStepFilename(const std::string& prefix,
                                          Index stepIndex,
                                          Index totalDigits = 6);

    void setFieldNames(const std::vector<std::string>& scalars,
                       const std::vector<std::string>& vectors);

private:
    void writeHeader(std::ofstream& ofs, Scalar time);
    void writePoints(std::ofstream& ofs);
    void writeCells(std::ofstream& ofs);
    void writeCellData(std::ofstream& ofs, const std::vector<Scalar>& damage);
    void writePointData(std::ofstream& ofs, const FieldVariables& fields,
                        const VectorX& vehiclePos, const VectorX& wheelPress);

    void writeScalarField(std::ofstream& ofs, const std::string& name,
                          const VectorX& field);
    void writeVectorField(std::ofstream& ofs, const std::string& name,
                          const VectorX& fieldX, const VectorX& fieldY);

    void ensureOutputDirectory();

    const Mesh2D& mesh_;
    std::string outputDir_ = "./vtk_output";
    int precision_ = 8;

    std::vector<std::string> scalarNames_;
    std::vector<std::string> vectorNames_;
};

} // namespace RoadbedSim
