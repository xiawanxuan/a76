#include "NewtonRaphsonSolver.h"

namespace RoadbedSim {

NewtonRaphsonSolver::NewtonRaphsonSolver(const Mesh2D& mesh,
                                         const ThermoHydro& th,
                                         Assembler& asmb)
    : mesh_(mesh), thermoHydro_(th), assembler_(asmb) {
    iterativeSolver_.setTolerance(iterTol_);
}

void NewtonRaphsonSolver::setInitialConditions(const VectorX& temp0,
                                                const VectorX& water0,
                                                const VectorX& dispX0,
                                                const VectorX& dispY0) {
    Index N = mesh_.getNumNodes();
    current_.temperature = temp0;
    current_.waterContent = water0;
    current_.displaceX = dispX0;
    current_.displaceY = dispY0;
    current_.iceContent = VectorX::Zero(N);
    current_.stressXX = VectorX::Zero(N);
    current_.stressYY = VectorX::Zero(N);
    current_.stressXY = VectorX::Zero(N);
    previous_ = current_;
    increment_ = current_;

    vehiclePosField_ = VectorX::Zero(N);
    wheelPressField_ = VectorX::Zero(N);
    cachedNodalLoad_ = VectorX::Zero(4 * N);
}

void NewtonRaphsonSolver::packStateVector(const VectorX& T, const VectorX& W,
                                           const VectorX& UX, const VectorX& UY,
                                           VectorX& state) const {
    Index N = mesh_.getNumNodes();
    state.resize(4 * N);
    state << T, W, UX, UY;
}

void NewtonRaphsonSolver::unpackStateVector(const VectorX& state,
                                             VectorX& T, VectorX& W,
                                             VectorX& UX, VectorX& UY) const {
    Index N = mesh_.getNumNodes();
    T = state.segment(0, N);
    W = state.segment(N, N);
    UX = state.segment(2 * N, N);
    UY = state.segment(3 * N, N);
}

bool NewtonRaphsonSolver::solveLinearSystem(const SparseMatrixX& A,
                                             VectorX& x, const VectorX& b) {
    if (useDirectSolver_) {
        directSolver_.compute(A);
        if (directSolver_.info() != Eigen::Success) return false;
        x = directSolver_.solve(b);
        return directSolver_.info() == Eigen::Success;
    } else {
        if (!iterativeSolver_.isInitialized()) {
            iterativeSolver_.compute(A);
        } else {
            iterativeSolver_.analyzePattern(A);
            iterativeSolver_.factorize(A);
        }
        if (iterativeSolver_.info() != Eigen::Success) return false;
        x = iterativeSolver_.solve(b);
        return iterativeSolver_.error() < iterTol_ * 10.0;
    }
}

SolverResult NewtonRaphsonSolver::solveNonlinearStep(Scalar currentTime) {
    SolverResult result{false, 0, 0.0, 0.0};
    auto startTime = std::chrono::high_resolution_clock::now();

    Index N = mesh_.getNumNodes();
    Index totalDofs = 4 * N;

    if (trafficCouplingEnabled_ && vehicleLoad_) {
        cachedNodalLoad_ = vehicleLoad_->assembleNodalLoadVector(currentTime);
        vehiclePosField_ = vehicleLoad_->computeVehiclePositionField(currentTime);
        wheelPressField_ = vehicleLoad_->computeWheelPressureField(currentTime);
        Scalar totalForce = vehicleLoad_->estimateTotalTrafficForce(currentTime);
        trafficForceHistory_.push_back(totalForce);
        vehicleLoad_->recordStep(static_cast<Index>(trafficForceHistory_.size() - 1),
                                   currentTime);
        if (vehicleFieldCallback_) {
            vehicleFieldCallback_(static_cast<Index>(trafficForceHistory_.size() - 1),
                                  vehiclePosField_, wheelPressField_);
        }
    }

    VectorX state, statePrev;
    packStateVector(current_.temperature, current_.waterContent,
                    current_.displaceX, current_.displaceY, state);
    packStateVector(previous_.temperature, previous_.waterContent,
                    previous_.displaceX, previous_.displaceY, statePrev);

    VectorX delta = VectorX::Zero(totalDofs);
    Scalar residualNorm = 0.0;

    for (result.iterations = 0; result.iterations < maxIterations_; ++result.iterations) {
        VectorX Tk, Wk, UXk, UYk;
        unpackStateVector(state + delta, Tk, Wk, UXk, UYk);
        VectorX dispFlat(2 * N);
        for (Index i = 0; i < N; ++i) {
            dispFlat(2 * i) = UXk(i);
            dispFlat(2 * i + 1) = UYk(i);
        }

        SparseMatrixX J;
        VectorX R;
        assembler_.assembleCoupledJacobian(Tk, Wk, dispFlat, VectorX::Zero(2 * N),
                                           dt_, J);
        assembler_.assembleCoupledResidual(Tk, previous_.temperature,
                                           Wk, previous_.waterContent,
                                           dispFlat, VectorX::Zero(2 * N),
                                           dt_, R);

        if (trafficCouplingEnabled_ && cachedNodalLoad_.size() == totalDofs) {
            assembler_.addVehicleNodalLoad(R, cachedNodalLoad_);
        }

        assembler_.applyBCToJacobian(J, R, bcs_, currentTime);

        residualNorm = R.norm();
        convLog_.residuals.push_back(residualNorm);
        convLog_.iterationCounts.push_back(result.iterations);

        if (residualCallback_) {
            residualCallback_(result.iterations, residualNorm, currentTime);
        }

        if (residualNorm < tolerance_) {
            result.converged = true;
            break;
        }

        VectorX deltaStep;
        if (!solveLinearSystem(J, deltaStep, -R)) {
            break;
        }

        Scalar maxStep = deltaStep.cwiseAbs().maxCoeff();
        if (maxStep > 10.0) {
            deltaStep *= 10.0 / maxStep;
        }
        delta += deltaStep;
    }

    result.finalResidual = residualNorm;

    if (result.converged) {
        unpackStateVector(state + delta, current_.temperature, current_.waterContent,
                          current_.displaceX, current_.displaceY);

        for (Index i = 0; i < N; ++i) {
            const auto& zone = (i < static_cast<Index>(mesh_.getNodeElements().size())
                && !mesh_.getNodeElements()[i].empty())
                ? mesh_.getElementZone(mesh_.getNodeElements()[i][0])
                : mesh_.getZone(0);
            current_.iceContent(i) = thermoHydro_.iceContent(
                current_.temperature(i), zone.waterContent, zone.porosity);
        }

        computeStressField();
    }

    auto endTime = std::chrono::high_resolution_clock::now();
    result.solveTime = std::chrono::duration<Scalar>(endTime - startTime).count();
    return result;
}

void NewtonRaphsonSolver::computeStressField() {
    Index N = mesh_.getNumNodes();
    VectorX dispFlat(2 * N);
    for (Index i = 0; i < N; ++i) {
        dispFlat(2 * i) = current_.displaceX(i);
        dispFlat(2 * i + 1) = current_.displaceY(i);
    }

    current_.stressXX.setZero(N);
    current_.stressYY.setZero(N);
    current_.stressXY.setZero(N);
    VectorX nodeCounts = VectorX::Zero(N);

    for (Index eid = 0; eid < mesh_.getNumElements(); ++eid) {
        const auto& e = mesh_.getElement(eid);
        const auto& zone = mesh_.getElementZone(eid);

        Eigen::Matrix<Scalar, 3, 6> B;
        Eigen::Matrix<Scalar, 3, 3> D;
        using BMatType = Eigen::Matrix<Scalar, 3, 6>;
        B = BMatType::Zero();
        assembler_.computeBMatrix(eid, B);

        Scalar avgT = 0;
        for (int k = 0; k < 3; ++k) avgT += current_.temperature(e.nodeIds[k]);
        avgT /= 3.0;
        Scalar ice = thermoHydro_.iceContent(avgT, zone.waterContent, zone.porosity);
        Scalar Eeff = zone.youngModulus * (1.0 + 5.0 * ice / (zone.porosity + 1e-10));
        Eeff = std::min(Eeff, 10.0 * zone.youngModulus);

        assembler_.computeDMatrix(Eeff, zone.poissonRatio, D);

        Eigen::Matrix<Scalar, 6, 1> ue;
        for (int k = 0; k < 3; ++k) {
            ue(2 * k) = dispFlat(2 * e.nodeIds[k]);
            ue(2 * k + 1) = dispFlat(2 * e.nodeIds[k] + 1);
        }

        Scalar heave = thermoHydro_.frostHeaveStrain(avgT, zone.waterContent, zone.porosity);
        Eigen::Matrix<Scalar, 3, 1> eps0;
        eps0 << heave, heave, 0.0;

        Eigen::Matrix<Scalar, 3, 1> sigma = D * (B * ue - eps0);

        for (int k = 0; k < 3; ++k) {
            Index nid = e.nodeIds[k];
            current_.stressXX(nid) += sigma(0);
            current_.stressYY(nid) += sigma(1);
            current_.stressXY(nid) += sigma(2);
            nodeCounts(nid) += 1.0;
        }
    }

    for (Index i = 0; i < N; ++i) {
        if (nodeCounts(i) > 0) {
            current_.stressXX(i) /= nodeCounts(i);
            current_.stressYY(i) /= nodeCounts(i);
            current_.stressXY(i) /= nodeCounts(i);
        }
    }
}

Scalar NewtonRaphsonSolver::computeAverageSettlement() const {
    auto topNodes = mesh_.findSurfaceNodes();
    if (topNodes.empty()) return 0.0;
    Scalar sum = 0.0;
    for (Index nid : topNodes) {
        sum += -current_.displaceY(nid);
    }
    return sum / static_cast<Scalar>(topNodes.size());
}

SolverResult NewtonRaphsonSolver::solveStep(Scalar currentTime) {
    SolverResult result = solveNonlinearStep(currentTime);
    if (result.converged) {
        previous_ = current_;
    }
    return result;
}

void NewtonRaphsonSolver::solveTransient() {
    Index numSteps = static_cast<Index>(std::ceil(totalTime_ / dt_));
    settlementHistory_.reserve(numSteps / outputInterval_ + 1);
    timeHistory_.reserve(numSteps / outputInterval_ + 1);

    auto surfaceNodes = mesh_.findSurfaceNodes();
    Index numElems = mesh_.getNumElements();

    for (Index step = 0; step < numSteps; ++step) {
        Scalar currentTime = step * dt_;
        SolverResult result = solveStep(currentTime);

        if (stepCallback_) {
            stepCallback_(step, currentTime);
        }

        if ((step + 1) % outputInterval_ == 0 || step == 0) {
            settlementHistory_.push_back(computeAverageSettlement());
            timeHistory_.push_back(currentTime + dt_);

            std::vector<Scalar> elemDamage(numElems, 0.0);
            for (Index eid = 0; eid < numElems; ++eid) {
                const auto& e = mesh_.getElement(eid);
                Scalar avgT = 0;
                for (int k = 0; k < 3; ++k) {
                    avgT += current_.temperature(e.nodeIds[k]);
                }
                avgT /= 3.0;
                if (avgT < FREEZING_POINT - 0.25 && avgT > FREEZING_POINT - 1.0) {
                    elemDamage[eid] = 0.3 + 0.7 * std::abs(avgT - FREEZING_POINT + 0.5);
                }
            }
            damageHistory_.push_back(std::move(elemDamage));
        }

        convLog_.timeSteps.push_back(currentTime);
    }
}

} // namespace RoadbedSim
