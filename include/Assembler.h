#pragma once

#include "Types.h"
#include "Mesh2D.h"
#include "ThermoHydro.h"
#include <vector>
#include <functional>
#include <memory>

namespace RoadbedSim {

class Assembler {
public:
    Assembler(const Mesh2D& mesh, const ThermoHydro& th);
    ~Assembler() = default;

    enum class FieldType { Temperature, Water, Displacement };

    void assembleThermalSystem(const VectorX& tempPrev,
                               const VectorX& waterPrev,
                               const VectorX& dispPrev,
                               Scalar dt,
                               SparseMatrixX& K, VectorX& F);

    void assembleHydraulicSystem(const VectorX& tempPrev,
                                 const VectorX& waterPrev,
                                 const VectorX& dispPrev,
                                 Scalar dt,
                                 SparseMatrixX& K, VectorX& F);

    using Vector2 = Eigen::Matrix<Scalar, 2, 1>;
    using Matrix3 = Eigen::Matrix<Scalar, 3, 3>;
    using Matrix2 = Eigen::Matrix<Scalar, 2, 2>;
    using Matrix6 = Eigen::Matrix<Scalar, 6, 6>;
    using Vector6 = Eigen::Matrix<Scalar, 6, 1>;

    void assembleMechanicalSystem(const VectorX& tempPrev,
                                  const VectorX& waterPrev,
                                  const VectorX& dispPrev,
                                  Scalar dt,
                                  SparseMatrixX& K, VectorX& F);

    void assembleCoupledJacobian(const VectorX& temp,
                                 const VectorX& water,
                                 const VectorX& dispPrev,
                                 const VectorX& dispIncr,
                                 Scalar dt,
                                 SparseMatrixX& J);

    void assembleCoupledResidual(const VectorX& temp,
                                 const VectorX& tempPrev,
                                 const VectorX& water,
                                 const VectorX& waterPrev,
                                 const VectorX& dispPrev,
                                 const VectorX& dispIncr,
                                 Scalar dt,
                                 VectorX& R);

    void applyBoundaryConditions(SparseMatrixX& K, VectorX& F,
                                 const std::vector<BoundaryCondition>& bcs,
                                 Scalar time);

    void applyBCToJacobian(SparseMatrixX& J, VectorX& R,
                           const std::vector<BoundaryCondition>& bcs,
                           Scalar time);

    Index getThermalDofStart() const { return 0; }
    Index getHydraulicDofStart() const { return getNumNodes(); }
    Index getMechanicalDofStart() const { return 2 * getNumNodes(); }
    Index getTotalDofs() const { return 4 * getNumNodes(); }
    Index getNumNodes() const { return mesh_.getNumNodes(); }

    std::vector<Scalar> computeElementFreezeThawDamage(
        const VectorX& tempHistory,
        Scalar damageThreshold = 0.5) const;

    VectorX computeSettlementDisplacement(const VectorX& disp) const;

    void setBlockSize(Index size) { blockSize_ = size; }
    Index getBlockSize() const { return blockSize_; }

    void computeBMatrix(Index elemId, Eigen::Matrix<Scalar, 3, 6>& B) const;
    void computeDMatrix(Scalar E, Scalar nu, Eigen::Matrix<Scalar, 3, 3>& D) const;

#ifdef USE_BLOCK_LOADING
    void assembleThermalBlock(const VectorX& tempPrev,
                              const VectorX& waterPrev,
                              Scalar dt,
                              Index blockStart, Index blockEnd,
                              std::vector<Triplet>& triplets,
                              VectorX& F);
#endif

private:
    void computeElementThermalContribution(Index elemId,
                                           const std::array<Scalar, 3>& T,
                                           const std::array<Scalar, 3>& W,
                                           Scalar dt,
                                           std::vector<Triplet>& triplets,
                                           VectorX& F);
    void computeElementHydraulicContribution(Index elemId,
                                             const std::array<Scalar, 3>& T,
                                             const std::array<Scalar, 3>& W,
                                             Scalar dt,
                                             std::vector<Triplet>& triplets,
                                             VectorX& F);
    void computeElementMechanicalContribution(Index elemId,
                                              const std::array<Scalar, 3>& T,
                                              const std::array<Scalar, 3>& W,
                                              const std::array<Vector2, 3>& U,
                                              Scalar dt,
                                              std::vector<Triplet>& triplets,
                                              VectorX& F);

    Matrix3 computeThermalElementMatrix(Index elemId, Scalar kEff, Scalar cEff,
                                        Scalar rho, Scalar dt) const;
    Matrix3 computeHydraulicElementMatrix(Index elemId, Scalar kHyd, Scalar dt) const;
    Matrix6 computeMechanicalElementMatrix(Index elemId, Scalar E, Scalar nu,
                                           Scalar alphaFT) const;

    void computeGradientOperator(Index elemId, MatrixX& B) const;

    const Mesh2D& mesh_;
    const ThermoHydro& thermoHydro_;
    Index blockSize_ = 1000;
};

} // namespace RoadbedSim
