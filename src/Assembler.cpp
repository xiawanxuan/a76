#include "Assembler.h"

namespace RoadbedSim {

Assembler::Assembler(const Mesh2D& mesh, const ThermoHydro& th)
    : mesh_(mesh), thermoHydro_(th) {}

void Assembler::computeGradientOperator(Index elemId, MatrixX& B) const {
    B.resize(2, 3);
    const auto& e = mesh_.getElement(elemId);
    const auto& n0 = mesh_.getNode(e.nodeIds[0]);
    const auto& n1 = mesh_.getNode(e.nodeIds[1]);
    const auto& n2 = mesh_.getNode(e.nodeIds[2]);
    Scalar twoA = 2.0 * e.area;
    if (twoA < 1e-20) twoA = 1e-20;

    B(0, 0) = (n1.y - n2.y) / twoA;
    B(0, 1) = (n2.y - n0.y) / twoA;
    B(0, 2) = (n0.y - n1.y) / twoA;
    B(1, 0) = (n2.x - n1.x) / twoA;
    B(1, 1) = (n0.x - n2.x) / twoA;
    B(1, 2) = (n1.x - n0.x) / twoA;
}

Assembler::Matrix3 Assembler::computeThermalElementMatrix(Index elemId, Scalar kEff,
                                                           Scalar cEff, Scalar rho,
                                                           Scalar dt) const {
    MatrixX B(2, 3);
    computeGradientOperator(elemId, B);
    Scalar area = mesh_.getElement(elemId).area;

    Matrix3 Kcond = area * kEff * (B.transpose() * B);

    Matrix3 Mcap;
    Mcap << 2.0, 1.0, 1.0,
            1.0, 2.0, 1.0,
            1.0, 1.0, 2.0;
    Mcap *= (rho * cEff * area) / 12.0;

    return Mcap / dt + Kcond;
}

Assembler::Matrix3 Assembler::computeHydraulicElementMatrix(Index elemId, Scalar kHyd,
                                                             Scalar dt) const {
    MatrixX B(2, 3);
    computeGradientOperator(elemId, B);
    Scalar area = mesh_.getElement(elemId).area;

    Matrix3 Kperm = area * kHyd * (B.transpose() * B);

    Matrix3 Mstorage;
    Mstorage << 2.0, 1.0, 1.0,
                1.0, 2.0, 1.0,
                1.0, 1.0, 2.0;
    const auto& zone = mesh_.getElementZone(elemId);
    Scalar Ss = zone.porosity * 1e-4;
    Mstorage *= (Ss * area) / 12.0;

    return Mstorage / dt + Kperm;
}

void Assembler::computeBMatrix(Index elemId, Eigen::Matrix<Scalar, 3, 6>& B) const {
    const auto& e = mesh_.getElement(elemId);
    const auto& n0 = mesh_.getNode(e.nodeIds[0]);
    const auto& n1 = mesh_.getNode(e.nodeIds[1]);
    const auto& n2 = mesh_.getNode(e.nodeIds[2]);
    Scalar twoA = 2.0 * e.area;
    if (twoA < 1e-20) twoA = 1e-20;

    Scalar b0 = n1.y - n2.y, b1 = n2.y - n0.y, b2 = n0.y - n1.y;
    Scalar c0 = n2.x - n1.x, c1 = n0.x - n2.x, c2 = n1.x - n0.x;

    B.setZero();
    B(0, 0) = b0 / twoA; B(0, 2) = b1 / twoA; B(0, 4) = b2 / twoA;
    B(1, 1) = c0 / twoA; B(1, 3) = c1 / twoA; B(1, 5) = c2 / twoA;
    B(2, 0) = c0 / twoA; B(2, 1) = b0 / twoA;
    B(2, 2) = c1 / twoA; B(2, 3) = b1 / twoA;
    B(2, 4) = c2 / twoA; B(2, 5) = b2 / twoA;
}

void Assembler::computeDMatrix(Scalar E, Scalar nu, Eigen::Matrix<Scalar, 3, 3>& D) const {
    Scalar factor = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
    D << (1.0 - nu), nu, 0.0,
         nu, (1.0 - nu), 0.0,
         0.0, 0.0, (1.0 - 2.0 * nu) / 2.0;
    D *= factor;
}

Assembler::Matrix6 Assembler::computeMechanicalElementMatrix(Index elemId, Scalar E,
                                                              Scalar nu,
                                                              Scalar alphaFT) const {
    Eigen::Matrix<Scalar, 3, 6> B;
    computeBMatrix(elemId, B);
    Eigen::Matrix<Scalar, 3, 3> D;
    computeDMatrix(E, nu, D);
    Scalar area = mesh_.getElement(elemId).area;
    return area * (B.transpose() * D * B);
}

void Assembler::computeElementThermalContribution(Index elemId,
                                                   const std::array<Scalar, 3>& T,
                                                   const std::array<Scalar, 3>& W,
                                                   Scalar dt,
                                                   std::vector<Triplet>& triplets,
                                                   VectorX& F) {
    const auto& e = mesh_.getElement(elemId);
    const auto& zone = mesh_.getElementZone(elemId);

    Scalar avgT = (T[0] + T[1] + T[2]) / 3.0;
    Scalar avgW = (W[0] + W[1] + W[2]) / 3.0;

    Scalar kEff, cEff, kHyd, heave;
    thermoHydro_.computeElementState(elemId, avgT, avgW, kEff, cEff, kHyd, heave);

    Matrix3 Ke = computeThermalElementMatrix(elemId, kEff, cEff, zone.density, dt);

    Matrix3 Mcap;
    Mcap << 2.0, 1.0, 1.0,
            1.0, 2.0, 1.0,
            1.0, 1.0, 2.0;
    Mcap *= (zone.density * cEff * e.area) / 12.0;

    Eigen::Matrix<Scalar, 3, 1> Tvec(T[0], T[1], T[2]);
    Eigen::Matrix<Scalar, 3, 1> Fe = (Mcap / dt) * Tvec;

    for (int i = 0; i < 3; ++i) {
        Index gi = e.nodeIds[i];
        F(gi) += Fe(i);
        for (int j = 0; j < 3; ++j) {
            Index gj = e.nodeIds[j];
            triplets.emplace_back(gi, gj, Ke(i, j));
        }
    }
}

void Assembler::computeElementHydraulicContribution(Index elemId,
                                                     const std::array<Scalar, 3>& T,
                                                     const std::array<Scalar, 3>& W,
                                                     Scalar dt,
                                                     std::vector<Triplet>& triplets,
                                                     VectorX& F) {
    const auto& e = mesh_.getElement(elemId);
    const auto& zone = mesh_.getElementZone(elemId);

    Scalar avgT = (T[0] + T[1] + T[2]) / 3.0;
    Scalar kHyd = thermoHydro_.hydraulicConductivity(avgT, zone);

    Scalar avgW = (W[0] + W[1] + W[2]) / 3.0;
    Scalar dWdT = thermoHydro_.unfrozenWaterDerivative(avgT, zone.waterContent, zone.porosity);

    Matrix3 Ke = computeHydraulicElementMatrix(elemId, kHyd, dt);

    const auto& n0 = mesh_.getNode(e.nodeIds[0]);
    const auto& n1 = mesh_.getNode(e.nodeIds[1]);
    const auto& n2 = mesh_.getNode(e.nodeIds[2]);
    Scalar avgY = (n0.y + n1.y + n2.y) / 3.0;
    Scalar gravityTerm = zone.porosity * kHyd * GRAVITY * e.area / 3.0;

    Matrix3 Mstorage;
    Mstorage << 2.0, 1.0, 1.0,
                1.0, 2.0, 1.0,
                1.0, 1.0, 2.0;
    Scalar Ss = zone.porosity * 1e-4;
    Mstorage *= (Ss * e.area) / 12.0;

    Eigen::Matrix<Scalar, 3, 1> Wvec(W[0], W[1], W[2]);
    Eigen::Matrix<Scalar, 3, 1> Fe = (Mstorage / dt) * Wvec;
    Eigen::Matrix<Scalar, 3, 1> gravVec(gravityTerm, gravityTerm, gravityTerm);
    Fe -= gravVec;

    Index N = mesh_.getNumNodes();
    for (int i = 0; i < 3; ++i) {
        Index gi = N + e.nodeIds[i];
        F(gi) += Fe(i);
        for (int j = 0; j < 3; ++j) {
            Index gj = N + e.nodeIds[j];
            triplets.emplace_back(gi, gj, Ke(i, j));
        }
    }
}

void Assembler::computeElementMechanicalContribution(Index elemId,
                                                      const std::array<Scalar, 3>& T,
                                                      const std::array<Scalar, 3>& W,
                                                      const std::array<Vector2, 3>& U,
                                                      Scalar dt,
                                                      std::vector<Triplet>& triplets,
                                                      VectorX& F) {
    const auto& e = mesh_.getElement(elemId);
    const auto& zone = mesh_.getElementZone(elemId);

    Scalar avgT = (T[0] + T[1] + T[2]) / 3.0;
    Scalar avgW = (W[0] + W[1] + W[2]) / 3.0;
    Scalar ice = thermoHydro_.iceContent(avgT, zone.waterContent, zone.porosity);
    Scalar E = zone.youngModulus * (1.0 + 5.0 * ice / (zone.porosity + 1e-10));
    Scalar Eeff = std::min(E, 10.0 * zone.youngModulus);
    Scalar heave = thermoHydro_.frostHeaveStrain(avgT, avgW, zone.porosity);

    Matrix6 Ke = computeMechanicalElementMatrix(elemId, Eeff, zone.poissonRatio, 0.0);

    Eigen::Matrix<Scalar, 3, 6> B;
    computeBMatrix(elemId, B);
    Eigen::Matrix<Scalar, 3, 3> D;
    computeDMatrix(Eeff, zone.poissonRatio, D);

    Eigen::Matrix<Scalar, 3, 1> eps0;
    eps0 << heave, heave, 0.0;
    Eigen::Matrix<Scalar, 6, 1> Fe = -e.area * (B.transpose() * (D * eps0));

    Scalar weight = zone.density * GRAVITY * e.area / 3.0;
    Index N = mesh_.getNumNodes();
    for (int i = 0; i < 3; ++i) {
        Index gx = 2 * N + 2 * e.nodeIds[i];
        Index gy = gx + 1;
        F(gx) += Fe(2 * i);
        F(gy) += Fe(2 * i + 1) - weight;
        for (int j = 0; j < 3; ++j) {
            Index gxj = 2 * N + 2 * e.nodeIds[j];
            Index gyj = gxj + 1;
            triplets.emplace_back(gx, gxj, Ke(2 * i, 2 * j));
            triplets.emplace_back(gx, gyj, Ke(2 * i, 2 * j + 1));
            triplets.emplace_back(gy, gxj, Ke(2 * i + 1, 2 * j));
            triplets.emplace_back(gy, gyj, Ke(2 * i + 1, 2 * j + 1));
        }
    }
}

void Assembler::assembleThermalSystem(const VectorX& tempPrev,
                                       const VectorX& waterPrev,
                                       const VectorX& dispPrev,
                                       Scalar dt,
                                       SparseMatrixX& K, VectorX& F) {
    Index N = mesh_.getNumNodes();
    Index totalDofs = 4 * N;
    std::vector<Triplet> triplets;
    triplets.reserve(mesh_.getNumElements() * 9 * 4);
    F = VectorX::Zero(totalDofs);

    for (Index eid = 0; eid < mesh_.getNumElements(); ++eid) {
        const auto& e = mesh_.getElement(eid);
        std::array<Scalar, 3> T{tempPrev(e.nodeIds[0]), tempPrev(e.nodeIds[1]), tempPrev(e.nodeIds[2])};
        std::array<Scalar, 3> W{waterPrev(e.nodeIds[0]), waterPrev(e.nodeIds[1]), waterPrev(e.nodeIds[2])};
        computeElementThermalContribution(eid, T, W, dt, triplets, F);
    }

    K.resize(totalDofs, totalDofs);
    K.setFromTriplets(triplets.begin(), triplets.end());
}

void Assembler::assembleHydraulicSystem(const VectorX& tempPrev,
                                         const VectorX& waterPrev,
                                         const VectorX& dispPrev,
                                         Scalar dt,
                                         SparseMatrixX& K, VectorX& F) {
    Index N = mesh_.getNumNodes();
    Index totalDofs = 4 * N;
    std::vector<Triplet> triplets;
    triplets.reserve(mesh_.getNumElements() * 9 * 4);
    F = VectorX::Zero(totalDofs);

    for (Index eid = 0; eid < mesh_.getNumElements(); ++eid) {
        const auto& e = mesh_.getElement(eid);
        std::array<Scalar, 3> T{tempPrev(e.nodeIds[0]), tempPrev(e.nodeIds[1]), tempPrev(e.nodeIds[2])};
        std::array<Scalar, 3> W{waterPrev(e.nodeIds[0]), waterPrev(e.nodeIds[1]), waterPrev(e.nodeIds[2])};
        computeElementHydraulicContribution(eid, T, W, dt, triplets, F);
    }

    K.resize(totalDofs, totalDofs);
    K.setFromTriplets(triplets.begin(), triplets.end());
}

void Assembler::assembleMechanicalSystem(const VectorX& tempPrev,
                                          const VectorX& waterPrev,
                                          const VectorX& dispPrev,
                                          Scalar dt,
                                          SparseMatrixX& K, VectorX& F) {
    Index N = mesh_.getNumNodes();
    Index totalDofs = 4 * N;
    std::vector<Triplet> triplets;
    triplets.reserve(mesh_.getNumElements() * 36);
    F = VectorX::Zero(totalDofs);

    for (Index eid = 0; eid < mesh_.getNumElements(); ++eid) {
        const auto& e = mesh_.getElement(eid);
        std::array<Scalar, 3> T{tempPrev(e.nodeIds[0]), tempPrev(e.nodeIds[1]), tempPrev(e.nodeIds[2])};
        std::array<Scalar, 3> W{waterPrev(e.nodeIds[0]), waterPrev(e.nodeIds[1]), waterPrev(e.nodeIds[2])};
        std::array<Vector2, 3> U{
            Vector2(dispPrev(2 * e.nodeIds[0]), dispPrev(2 * e.nodeIds[0] + 1)),
            Vector2(dispPrev(2 * e.nodeIds[1]), dispPrev(2 * e.nodeIds[1] + 1)),
            Vector2(dispPrev(2 * e.nodeIds[2]), dispPrev(2 * e.nodeIds[2] + 1))
        };
        computeElementMechanicalContribution(eid, T, W, U, dt, triplets, F);
    }

    K.resize(totalDofs, totalDofs);
    K.setFromTriplets(triplets.begin(), triplets.end());
}

void Assembler::assembleCoupledJacobian(const VectorX& temp,
                                         const VectorX& water,
                                         const VectorX& dispPrev,
                                         const VectorX& dispIncr,
                                         Scalar dt,
                                         SparseMatrixX& J) {
    assembleThermalSystem(temp, water, dispPrev, dt, J, VectorX());
    SparseMatrixX K2, K3;
    VectorX F2, F3;
    assembleHydraulicSystem(temp, water, dispPrev, dt, K2, F2);
    assembleMechanicalSystem(temp, water, dispPrev, dt, K3, F3);
    J += K2 + K3;
}

void Assembler::assembleCoupledResidual(const VectorX& temp,
                                         const VectorX& tempPrev,
                                         const VectorX& water,
                                         const VectorX& waterPrev,
                                         const VectorX& dispPrev,
                                         const VectorX& dispIncr,
                                         Scalar dt,
                                         VectorX& R) {
    SparseMatrixX Kt, Kh, Km;
    VectorX Ft, Fh, Fm;
    assembleThermalSystem(temp, water, dispPrev, dt, Kt, Ft);
    assembleHydraulicSystem(temp, water, dispPrev, dt, Kh, Fh);
    assembleMechanicalSystem(temp, water, dispPrev, dt, Km, Fm);

    Index N = mesh_.getNumNodes();
    VectorX state(4 * N);
    state << temp, water, dispPrev;
    R = (Kt + Kh + Km) * state - (Ft + Fh + Fm);
}

void Assembler::applyBoundaryConditions(SparseMatrixX& K, VectorX& F,
                                         const std::vector<BoundaryCondition>& bcs,
                                         Scalar time) {
    Index N = mesh_.getNumNodes();
    for (const auto& bc : bcs) {
        Scalar val = bc.valueFunc(time);
        for (Index nid : bc.nodeIds) {
            Index row;
            switch (bc.type) {
                case BoundaryCondition::Type::Temperature: row = nid; break;
                case BoundaryCondition::Type::WaterFlux: row = N + nid; break;
                case BoundaryCondition::Type::DisplacementX: row = 2 * N + 2 * nid; break;
                case BoundaryCondition::Type::DisplacementY: row = 2 * N + 2 * nid + 1; break;
                default: continue;
            }
            for (SparseMatrixX::InnerIterator it(K, row); it; ++it) {
                if (static_cast<Index>(it.index()) == row) {
                    F(row) = val * it.value();
                } else {
                    F(it.index()) -= it.value() * val;
                }
            }
        }
    }
    for (const auto& bc : bcs) {
        Scalar val = bc.valueFunc(time);
        for (Index nid : bc.nodeIds) {
            Index row;
            switch (bc.type) {
                case BoundaryCondition::Type::Temperature: row = nid; break;
                case BoundaryCondition::Type::WaterFlux: row = N + nid; break;
                case BoundaryCondition::Type::DisplacementX: row = 2 * N + 2 * nid; break;
                case BoundaryCondition::Type::DisplacementY: row = 2 * N + 2 * nid + 1; break;
                default: continue;
            }
            K.row(row) *= Scalar(0);
            K.col(row) *= Scalar(0);
            K.coeffRef(row, row) = 1.0;
            F(row) = val;
        }
    }
}

void Assembler::applyBCToJacobian(SparseMatrixX& J, VectorX& R,
                                   const std::vector<BoundaryCondition>& bcs,
                                   Scalar time) {
    applyBoundaryConditions(J, R, bcs, time);
}

std::vector<Scalar> Assembler::computeElementFreezeThawDamage(
    const VectorX& tempHistory, Scalar damageThreshold) const {
    Index numElems = mesh_.getNumElements();
    std::vector<Scalar> damage(numElems, 0.0);
    Index numSteps = tempHistory.size() / mesh_.getNumNodes();
    for (Index eid = 0; eid < numElems; ++eid) {
        const auto& e = mesh_.getElement(eid);
        int freezeThawCycles = 0;
        bool wasFrozen = false;
        for (Index step = 0; step < numSteps; ++step) {
            Scalar avgT = 0;
            for (int i = 0; i < 3; ++i) {
                avgT += tempHistory(step * mesh_.getNumNodes() + e.nodeIds[i]);
            }
            avgT /= 3.0;
            bool isFrozen = avgT < FREEZING_POINT - 0.5;
            if (wasFrozen && !isFrozen) ++freezeThawCycles;
            wasFrozen = isFrozen;
        }
        damage[eid] = std::min(1.0, freezeThawCycles * damageThreshold);
    }
    return damage;
}

VectorX Assembler::computeSettlementDisplacement(const VectorX& disp) const {
    Index N = mesh_.getNumNodes();
    VectorX settlement(N);
    auto topNodes = mesh_.findSurfaceNodes();
    for (Index i = 0; i < N; ++i) {
        settlement(i) = -disp(2 * i + 1);
    }
    return settlement;
}

void Assembler::addVehicleNodalLoad(VectorX& F, const VectorX& vehicleLoad) const {
    Index totalDofs = getTotalDofs();
    if (vehicleLoad.size() != totalDofs) return;
    F += vehicleLoad;
}

Scalar Assembler::computeElementDynamicModulus(Index elemId, Scalar temperature,
                                                Scalar iceContent,
                                                Scalar staticE,
                                                Scalar loadFrequencyHz) const {
    const auto& zone = mesh_.getElementZone(elemId);
    Scalar porosity = zone.porosity;
    Scalar iceRatio = porosity > 1e-6 ? iceContent / porosity : 0.0;
    iceRatio = std::max(0.0, std::min(1.0, iceRatio));

    Scalar freqFactor = 1.0 + 0.05 * std::log10(std::max(0.1, loadFrequencyHz));
    Scalar iceFactor = 1.0 + 8.0 * iceRatio;
    Scalar tempFactor = 1.0;
    if (temperature < FREEZING_POINT) {
        Scalar deltaT = FREEZING_POINT - temperature;
        tempFactor = 1.0 + 0.06 * std::min(deltaT, 20.0);
    } else {
        Scalar deltaT = temperature - FREEZING_POINT;
        tempFactor = std::max(0.5, 1.0 - 0.08 * std::min(deltaT, 10.0));
    }
    if (temperature > FREEZING_POINT - 0.5 && temperature < FREEZING_POINT + 1.0) {
        Scalar t = (temperature - FREEZING_POINT + 0.5) / 1.5;
        Scalar thawWeak = 1.0 - 0.6 * std::exp(-t * 3.0);
        tempFactor *= thawWeak;
    }

    Scalar baseE = staticE > 0 ? staticE : zone.youngModulus;
    Scalar dynE = baseE * freqFactor * iceFactor * tempFactor;
    return std::min(dynE, 50.0 * baseE);
}

void Assembler::assembleMechanicalWithDynamic(const VectorX& tempPrev,
                                                const VectorX& waterPrev,
                                                const VectorX& iceContent,
                                                const VectorX& dispPrev,
                                                Scalar dt,
                                                Scalar loadFrequency,
                                                SparseMatrixX& K, VectorX& F) {
    Index N = mesh_.getNumNodes();
    Index totalDofs = 4 * N;
    std::vector<Triplet> triplets;
    triplets.reserve(mesh_.getNumElements() * 36);
    F = VectorX::Zero(totalDofs);

    for (Index eid = 0; eid < mesh_.getNumElements(); ++eid) {
        const auto& e = mesh_.getElement(eid);
        std::array<Scalar, 3> T{tempPrev(e.nodeIds[0]), tempPrev(e.nodeIds[1]), tempPrev(e.nodeIds[2])};
        std::array<Scalar, 3> W{waterPrev(e.nodeIds[0]), waterPrev(e.nodeIds[1]), waterPrev(e.nodeIds[2])};
        std::array<Vector2, 3> U{
            Vector2(dispPrev(2 * e.nodeIds[0]), dispPrev(2 * e.nodeIds[0] + 1)),
            Vector2(dispPrev(2 * e.nodeIds[1]), dispPrev(2 * e.nodeIds[1] + 1)),
            Vector2(dispPrev(2 * e.nodeIds[2]), dispPrev(2 * e.nodeIds[2] + 1))
        };

        Scalar avgT = (T[0] + T[1] + T[2]) / 3.0;
        Scalar avgIce = (iceContent(e.nodeIds[0]) +
                         iceContent(e.nodeIds[1]) +
                         iceContent(e.nodeIds[2])) / 3.0;
        const auto& zone = mesh_.getElementZone(eid);
        Scalar EStatic = zone.youngModulus * (1.0 + 5.0 * avgIce / (zone.porosity + 1e-10));
        EStatic = std::min(EStatic, 10.0 * zone.youngModulus);
        Scalar EDyn = computeElementDynamicModulus(eid, avgT, avgIce, EStatic, loadFrequency);

        Scalar heave = thermoHydro_.frostHeaveStrain(avgT, (W[0]+W[1]+W[2])/3.0, zone.porosity);
        Scalar weight = zone.density * GRAVITY * e.area / 3.0;

        Eigen::Matrix<Scalar, 3, 6> B;
        Eigen::Matrix<Scalar, 3, 3> Dmat;
        computeBMatrix(eid, B);
        computeDMatrix(EDyn, zone.poissonRatio, Dmat);

        Matrix6 Ke = e.area * (B.transpose() * Dmat * B);
        Eigen::Matrix<Scalar, 3, 1> eps0;
        eps0 << heave, heave, 0.0;
        Eigen::Matrix<Scalar, 6, 1> Fe = -e.area * (B.transpose() * (Dmat * eps0));

        for (int i = 0; i < 3; ++i) {
            Index gx = 2 * N + 2 * e.nodeIds[i];
            Index gy = gx + 1;
            F(gx) += Fe(2 * i);
            F(gy) += Fe(2 * i + 1) - weight;
            for (int j = 0; j < 3; ++j) {
                Index gxj = 2 * N + 2 * e.nodeIds[j];
                Index gyj = gxj + 1;
                triplets.emplace_back(gx, gxj, Ke(2 * i, 2 * j));
                triplets.emplace_back(gx, gyj, Ke(2 * i, 2 * j + 1));
                triplets.emplace_back(gy, gxj, Ke(2 * i + 1, 2 * j));
                triplets.emplace_back(gy, gyj, Ke(2 * i + 1, 2 * j + 1));
            }
        }
    }

    K.resize(totalDofs, totalDofs);
    K.setFromTriplets(triplets.begin(), triplets.end());
}

void Assembler::assembleCoupledResidualWithVehicle(const VectorX& temp,
                                                    const VectorX& tempPrev,
                                                    const VectorX& water,
                                                    const VectorX& waterPrev,
                                                    const VectorX& dispPrev,
                                                    const VectorX& dispIncr,
                                                    const VectorX& vehicleLoad,
                                                    Scalar dt,
                                                    VectorX& R) {
    assembleCoupledResidual(temp, tempPrev, water, waterPrev,
                            dispPrev, dispIncr, dt, R);
    addVehicleNodalLoad(R, vehicleLoad);
}

#ifdef USE_BLOCK_LOADING
void Assembler::assembleThermalBlock(const VectorX& tempPrev,
                                      const VectorX& waterPrev,
                                      Scalar dt,
                                      Index blockStart, Index blockEnd,
                                      std::vector<Triplet>& triplets,
                                      VectorX& F) {
    for (Index eid = blockStart; eid < blockEnd && eid < mesh_.getNumElements(); ++eid) {
        const auto& e = mesh_.getElement(eid);
        std::array<Scalar, 3> T{tempPrev(e.nodeIds[0]), tempPrev(e.nodeIds[1]), tempPrev(e.nodeIds[2])};
        std::array<Scalar, 3> W{waterPrev(e.nodeIds[0]), waterPrev(e.nodeIds[1]), waterPrev(e.nodeIds[2])};
        computeElementThermalContribution(eid, T, W, dt, triplets, F);
    }
}
#endif

} // namespace RoadbedSim
