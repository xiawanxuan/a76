#pragma once

#include "Types.h"
#include "Assembler.h"
#include "Mesh2D.h"
#include "ThermoHydro.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>
#include <memory>
#include <functional>
#include <chrono>

namespace RoadbedSim {

class NewtonRaphsonSolver {
public:
    NewtonRaphsonSolver(const Mesh2D& mesh, const ThermoHydro& th, Assembler& asmb);
    ~NewtonRaphsonSolver() = default;

    void setMaxIterations(Index n) { maxIterations_ = n; }
    Index getMaxIterations() const { return maxIterations_; }

    void setTolerance(Scalar tol) { tolerance_ = tol; }
    Scalar getTolerance() const { return tolerance_; }

    void setTimeStep(Scalar dt) { dt_ = dt; }
    Scalar getTimeStep() const { return dt_; }

    void setTotalTime(Scalar t) { totalTime_ = t; }
    Scalar getTotalTime() const { return totalTime_; }

    void setOutputInterval(Index n) { outputInterval_ = n; }
    Index getOutputInterval() const { return outputInterval_; }

    void setBoundaryConditions(const std::vector<BoundaryCondition>& bcs) { bcs_ = bcs; }
    const std::vector<BoundaryCondition>& getBoundaryConditions() const { return bcs_; }

    void setInitialConditions(const VectorX& temp0, const VectorX& water0,
                              const VectorX& dispX0, const VectorX& dispY0);

    SolverResult solveStep(Scalar currentTime);
    void solveTransient();

    const FieldVariables& getCurrentFields() const { return current_; }
    const FieldVariables& getPreviousFields() const { return previous_; }
    const ConvergenceLog& getConvergenceLog() const { return convLog_; }

    const std::vector<Scalar>& getSettlementHistory() const { return settlementHistory_; }
    const std::vector<Scalar>& getTimeHistory() const { return timeHistory_; }
    const std::vector<std::vector<Scalar>>& getDamageHistory() const { return damageHistory_; }

    void setResidualOutputCallback(std::function<void(Index, Scalar, Scalar)> cb) {
        residualCallback_ = std::move(cb);
    }

    void setStepCompleteCallback(std::function<void(Index, Scalar)> cb) {
        stepCallback_ = std::move(cb);
    }

    void useDirectSolver(bool direct) { useDirectSolver_ = direct; }
    bool isUsingDirectSolver() const { return useDirectSolver_; }

    void setIterativeSolverTolerance(Scalar tol) { iterTol_ = tol; }

private:
    void packStateVector(const VectorX& T, const VectorX& W,
                         const VectorX& UX, const VectorX& UY,
                         VectorX& state) const;
    void unpackStateVector(const VectorX& state,
                           VectorX& T, VectorX& W,
                           VectorX& UX, VectorX& UY) const;

    SolverResult solveNonlinearStep(Scalar currentTime);

    bool solveLinearSystem(const SparseMatrixX& A, VectorX& x, const VectorX& b);

    void computeStressField();
    Scalar computeAverageSettlement() const;

    const Mesh2D& mesh_;
    const ThermoHydro& thermoHydro_;
    Assembler& assembler_;

    FieldVariables current_;
    FieldVariables previous_;
    FieldVariables increment_;

    Index maxIterations_ = 50;
    Scalar tolerance_ = 1e-6;
    Scalar dt_ = 3600.0;
    Scalar totalTime_ = 3.154e7;
    Index outputInterval_ = 24;
    Scalar iterTol_ = 1e-8;

    bool useDirectSolver_ = true;

    std::vector<BoundaryCondition> bcs_;
    ConvergenceLog convLog_;

    std::vector<Scalar> settlementHistory_;
    std::vector<Scalar> timeHistory_;
    std::vector<std::vector<Scalar>> damageHistory_;

    std::function<void(Index, Scalar, Scalar)> residualCallback_;
    std::function<void(Index, Scalar)> stepCallback_;

    Eigen::SparseLU<SparseMatrixX> directSolver_;
    Eigen::BiCGSTAB<SparseMatrixX> iterativeSolver_;
};

} // namespace RoadbedSim
