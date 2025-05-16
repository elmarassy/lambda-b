
#include <iostream>
#include <cstdlib>
#include <vector>

#include "FCCAnalyses/myUtils.h"
#include "FCCAnalyses/VertexFitterSimple.h"
#include "FCCAnalyses/ReconstructedParticle.h"
#include "FCCAnalyses/MCParticle.h"
#include "FCCAnalyses/Algorithms.h"

namespace FCCAnalyses {
    namespace myUtils {






testing testLb2LMuMu(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
                               ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop) {
    testing result;

    FCCAnalysesComposite2 oldHadrons;
    FCCAnalysesComposite2 oldMuons;
    FCCAnalysesComposite2 oldTotals;
    double_t oldPrimaryLbImpactParameters; //minimal distance between reconstructed lambdaB path and primary vertex
    double_t oldMuonHadronImpactParameters; //minimal distance between dihadron path and dimuon vertex
    double_t oldDisplacementProduct;
    int oldHadron1Types;
    int oldHadron2Types;

    FCCAnalysesComposite2 newHadrons;
    FCCAnalysesComposite2 newMuons;
    FCCAnalysesComposite2 newTotals;
    double_t newPrimaryLbImpactParameters = 1e9; //minimal distance between reconstructed lambdaB path and primary vertex
    double_t newMuonHadronImpactParameters; //minimal distance between dihadron path and dimuon vertex
    double_t newDisplacementProduct;
    int newHadron1Types;
    int newHadron2Types;
    //currently storing type of hadrons to use for debugging, will not be needed later

    std::cout << "Running...\n";
    VertexingUtils::FCCAnalysesVertex primary;
    int count = 0;
    for (auto &v: vertex) {
        //find the primary vertex
        if (v.vertex.primary) {
            primary = v;
            count++;
        }
    }
    int counterP = 0;
    for (auto &p: vertex) {
        if (p.vertex.primary || p.ntracks != 2) {
            counterP += 1;
            continue;
        }
        int chargeMuons = 0;
        int numMuons = 0;
        for (auto &potentialMuon: p.reco_ind) {
            if (recop.at(potentialMuon).type == 13) {
                numMuons += 1;
                chargeMuons += recop.at(potentialMuon).charge;
            }
        }
        if (chargeMuons == 0 && numMuons == 2) {
            TLorentzVector dimuonMomentum = build_tlv(recop, p.reco_ind);

            FCCAnalysesComposite2 muon1;
            muon1.particle = ReconstructedParticle::get_tlv(recop[p.reco_ind.at(0)]);
            muon1.charge = recop[p.reco_ind.at(0)].charge;

            FCCAnalysesComposite2 muon2;
            muon2.particle = ReconstructedParticle::get_tlv(recop[p.reco_ind.at(1)]);
            muon2.charge = recop[p.reco_ind.at(1)].charge;

            FCCAnalysesComposite2 dimuon;
            dimuon.vertex = counterP;
            dimuon.particle = dimuonMomentum;
            dimuon.charge = 0;

            int counterQ = 0;
            for (auto &q: vertex) {
                if (q.vertex.primary || q.ntracks != 2) {
                    counterQ += 1;
                    continue;
                }
                int chargeHadrons = 0;
                // int numHadrons = 0;
                int numP = 0;
                int numPi = 0;
                for (auto &potentialHadron: q.reco_ind) {
                    if (recop.at(potentialHadron).type == 2212) {
                        numP += 1;
                        chargeHadrons += recop.at(potentialHadron).charge;
                    }
                    if (recop.at(potentialHadron).type == 211) {
                        numPi += 1;
                        chargeHadrons += recop.at(potentialHadron).charge;
                    }
                    // if (recop.at(potentialHadron).type != 13 && recop.at(potentialHadron).type != 11) {
                    //     numHadrons += 1;
                    //     chargeHadrons += recop.at(potentialHadron).charge;
                    // }
                }

                // if (numHadrons == 2 && chargeHadrons == 0) {
                if (numP == 1 && numPi == 1 && chargeHadrons == 0) {
                    FCCAnalysesComposite2 hadron1;
                    hadron1.particle = ReconstructedParticle::get_tlv(recop[q.reco_ind.at(0)]);
                    hadron1.charge = recop[q.reco_ind.at(0)].charge;
                    int hadron1Type = recop[q.reco_ind.at(0)].type;

                    FCCAnalysesComposite2 hadron2;
                    hadron2.particle = ReconstructedParticle::get_tlv(recop[q.reco_ind.at(1)]);
                    hadron2.charge = recop[q.reco_ind.at(1)].charge;
                    int hadron2Type = recop[q.reco_ind.at(1)].type;

                    // if (hadron1.charge < 0) { //always have hadron1.charge > 0
                    //     FCCAnalysesComposite2 tempHadron = hadron1;
                    //     hadron1 = hadron2;
                    //     hadron2 = tempHadron;
                    // }

                    TLorentzVector dihadronMomentum = build_tlv(recop, q.reco_ind);
                    TVector3 dihadronDisplacement = TVector3(q.vertex.position[0] - p.vertex.position[0],
                                                             q.vertex.position[1] - p.vertex.position[1],
                                                             q.vertex.position[2] - p.vertex.position[2]);

                    TVector3 dihadron3Momentum = TVector3(dihadronMomentum.Px(), dihadronMomentum.Py(),
                                                          dihadronMomentum.Pz());

                    double_t flightDistanceL = dihadronDisplacement.Mag();
                    double_t muonHadronImpactParam = (dihadronDisplacement - (
                                           dihadronDisplacement.Dot(dihadron3Momentum) / dihadron3Momentum.Dot(
                                               dihadron3Momentum)) * dihadron3Momentum).Mag();
                    TLorentzVector LbMomentum = dimuonMomentum + dihadronMomentum;
                    TVector3 Lb3Momentum = TVector3(LbMomentum.Px(), LbMomentum.Py(), LbMomentum.Pz());
                    TVector3 displacementTotal = TVector3(p.vertex.position[0] - primary.vertex.position[0],
                                                          p.vertex.position[1] - primary.vertex.position[1],
                                                          p.vertex.position[2] - primary.vertex.position[2]);

                    double_t flightDistanceLb = displacementTotal.Mag();
                    double_t primaryLbImpactParam = (displacementTotal - (displacementTotal.Dot(Lb3Momentum) / Lb3Momentum.Dot(Lb3Momentum)) * Lb3Momentum).Mag();

                    double_t displacementProduct = displacementTotal.Dot(dihadronDisplacement);

                    std::cout << "Found new potential pair.\n";
                    std::cout << "Coordinates are: " << primary.vertex.position << " primary, " << p.vertex.position << " muon, " << q.vertex.position << " hadron.\n";
                    std::cout << "Displacements are: (" << displacementTotal[0] << ", " << displacementTotal[1] << ", " << displacementTotal[2] << ") total, " <<
                        "(" << dihadronDisplacement[0] << ", " << dihadronDisplacement[1] << ", " << dihadronDisplacement[2] << "), hadron.\n";
                    std::cout << "Displacement product is: " << displacementProduct << std::endl;
                    std::cout << "3-momenta are: (" << Lb3Momentum[0] << ", " << Lb3Momentum[1] << ", " << Lb3Momentum[2] << ") Lb, " <<
                        "(" << dihadron3Momentum[0] << ", " << dihadron3Momentum[1] << ", " << dihadron3Momentum[2] << "), hadron.\n";
                    std::cout << "Computed impact parameters: " << primaryLbImpactParam << " total, " << muonHadronImpactParam << "muonHadron.\n";
                    std::cout << "Computed Lb mass: " << LbMomentum.Mag() << "\n";

                    FCCAnalysesComposite2 dihadron;
                    dihadron.vertex = counterQ;
                    dihadron.particle = dihadronMomentum;
                    dihadron.charge = 0;

                    FCCAnalysesComposite2 total;
                    total.particle = LbMomentum;
                    total.charge = 0;
                    total.vertex = 0;

                    if (primaryLbImpactParam < newPrimaryLbImpactParameters) {
                        newMuons = dimuon;
                        newHadrons = dihadron;
                        newTotals = total;
                        newDisplacementProduct = displacementProduct;
                        newPrimaryLbImpactParameters = primaryLbImpactParam;
                        newMuonHadronImpactParameters = muonHadronImpactParam;
                        newHadron1Types = hadron1Type;
                        newHadron2Types = hadron2Type;
                    }

                }
                counterQ += 1;
            }
        }
        counterP += 1;
    }

    FCCAnalysesComposite2 compL;
    FCCAnalysesComposite2 compMu;
    FCCAnalysesComposite2 savedL;
    FCCAnalysesComposite2 savedMu;

    //   ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>
    // SelPrimaryTracks(ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind,
    //                  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,
    //                  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc,
    //                  TVector3 MC_EventPrimaryVertex)
    int counterL = 0;
    int pcount = 0;
    double max = 0;
    bool valid = false;
    for (auto &p: vertex) {
        if (p.vertex.primary) {
            counterL += 1;
            continue;
        }
        if (p.ntracks != 2) {
            counterL += 1;
            continue;
        }
        int charge_p = 0;
        int nobj_p = 0;
        for (auto &r: p.reco_ind) {
            if (recop.at(r).type == 2212) {
                nobj_p += 1;
                charge_p += recop.at(r).charge;
            }
        }
        int charge_pi = 0;
        int nobj_pi = 0;
        for (auto &r: p.reco_ind) {
            if (recop.at(r).type == 211) {
                nobj_pi += 1;
                charge_pi += recop.at(r).charge;
            }
        }
        if (nobj_pi == 1 && nobj_p == 1 && charge_pi + charge_p == 0) {
            compL.vertex = counterL;
            compL.particle = build_tlv(recop, p.reco_ind);
            compL.charge = charge_pi + charge_p;
            pcount += 1;
            int counterMu = 0;
            for (auto &q: vertex) {
                if (q.vertex.primary) {
                    counterMu += 1;
                    continue;
                }
                if (q.ntracks != 2) {
                    counterMu += 1;
                    continue;
                }
                int charge_mumu = 0;
                int nobj_mumu = 0;
                for (auto &s: q.reco_ind) {
                    if (recop.at(s).type == 13) {
                        nobj_mumu += 1;
                        charge_mumu += recop.at(s).charge;
                    }
                }
                if (nobj_mumu == 2 && charge_mumu == 0) {
                    compMu.vertex = counterMu;
                    compMu.particle = build_tlv(recop, q.reco_ind);
                    compMu.charge = charge_pi + charge_p;

                    double alignment = compMu.particle[0] * compL.particle[0] + compMu.particle[1] * compL.
                                       particle[1] + compMu.particle[2] * compL.particle[2];
                    if (alignment > max) {
                        max = alignment;
                        savedL = compL;
                        savedMu = compMu;
                        valid = true;
                    }
                }
                counterMu += 1;
            }
        }
        counterL += 1;
    }

    if (valid) {
                TLorentzVector dimuonMomentum = savedMu.particle;
                TLorentzVector dihadronMomentum = savedL.particle;
                TVector3 dihadronDisplacement = TVector3(vertex.at(savedL.vertex).vertex.position[0] - vertex.at(savedMu.vertex).vertex.position[0],
                                                         vertex.at(savedL.vertex).vertex.position[1] - vertex.at(savedMu.vertex).vertex.position[1],
                                                         vertex.at(savedL.vertex).vertex.position[2] - vertex.at(savedMu.vertex).vertex.position[2]);

                TVector3 dihadron3Momentum = TVector3(dihadronMomentum.Px(), dihadronMomentum.Py(),
                                                      dihadronMomentum.Pz());

                double_t flightDistanceL = dihadronDisplacement.Mag();
                double_t muonHadronImpactParam = (dihadronDisplacement - (
                                       dihadronDisplacement.Dot(dihadron3Momentum) / dihadron3Momentum.Dot(
                                           dihadron3Momentum)) * dihadron3Momentum).Mag();

                TLorentzVector LbMomentum = dimuonMomentum + dihadronMomentum;
                TVector3 Lb3Momentum = TVector3(LbMomentum.Px(), LbMomentum.Py(), LbMomentum.Pz());
                TVector3 displacementTotal = TVector3(primary.vertex.position[0] - vertex.at(savedMu.vertex).vertex.position[0],
                                                 primary.vertex.position[1] - vertex.at(savedMu.vertex).vertex.position[1],
                                                 primary.vertex.position[2] - vertex.at(savedMu.vertex).vertex.position[2]);

                double_t flightDistanceLb = displacementTotal.Mag();
                double_t primaryLbImpactParam = (displacementTotal - (displacementTotal.Dot(Lb3Momentum) / Lb3Momentum.Dot(Lb3Momentum)) * Lb3Momentum).Mag();

                double_t displacementProduct = displacementTotal.Dot(dihadronDisplacement);

                FCCAnalysesComposite2 total;
                total.particle = LbMomentum;
                total.charge = 0;
                total.vertex = 0;
                if (newPrimaryLbImpactParameters < 1e9) {
                    result.newHadrons.push_back(newHadrons);
                    result.newMuons.push_back(newMuons);
                    result.newPrimaryLbImpactParameters.push_back(newPrimaryLbImpactParameters);
                    result.newTotals.push_back(newTotals);
                    result.newMuonHadronImpactParameters.push_back(newMuonHadronImpactParameters);
                    result.newDisplacementProduct.push_back(newDisplacementProduct);
                    result.newHadron1Types.push_back(newHadron1Types);
                    result.newHadron2Types.push_back(newHadron2Types);

                    result.oldHadrons.push_back(savedL);
                    result.oldMuons.push_back(savedMu);
                    result.oldPrimaryLbImpactParameters.push_back(primaryLbImpactParam);
                    result.oldTotals.push_back(total);
                    result.oldMuonHadronImpactParameters.push_back(muonHadronImpactParam);
                    result.oldDisplacementProduct.push_back(displacementProduct);
                    result.oldHadron1Types.push_back(recop.at(vertex.at(savedL.vertex).reco_ind.at(0)).type);
                    result.oldHadron2Types.push_back(recop.at(vertex.at(savedL.vertex).reco_ind.at(1)).type);
                }
            }

            return result;
}



        testing2 test2Lb2LMuMu(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
                                       ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop) {

            testing2 result;

            double initial = 1e9;

            ROOT::VecOps::RVec<FCCAnalysesComposite2> newDimuons; //dimuon momentum
            ROOT::VecOps::RVec<FCCAnalysesComposite2> newMuons1; //muon1 momentum
            ROOT::VecOps::RVec<FCCAnalysesComposite2> newMuons2; //muon2 momentum

            ROOT::VecOps::RVec<FCCAnalysesComposite2> newDihadrons; //dihadron momentum
            ROOT::VecOps::RVec<FCCAnalysesComposite2> newHadrons1; //hadron1 momentum
            ROOT::VecOps::RVec<FCCAnalysesComposite2> newHadrons2; //hadron2 momentum

            ROOT::VecOps::RVec<FCCAnalysesComposite2> newTotals;
            //total momentum (dimuon + dihadron = lambdaB for valid selection)

            ROOT::VecOps::RVec<double_t> newFlightDistancesLb; //distance from dimuon vertex to primary vertex
            ROOT::VecOps::RVec<double_t> newFlightDistancesL; //distance from dihadron vertex to dimuon vertex

            ROOT::VecOps::RVec<double_t> newPrimaryLbImpactParams; //minimal distance between reconstructed lambdaB path and primary vertex
            ROOT::VecOps::RVec<double_t> newMuonHadronImpactParams; //minimal distance between dihadron path and dimuon vertex
            ROOT::VecOps::RVec<double_t> newDisplacementProducts;

            ROOT::VecOps::RVec<int> newHadron1Types;
            ROOT::VecOps::RVec<int> newHadron2Types;
            //currently storing type of hadrons to use for debugging, will not be needed later


            const double_t pionMass = 0.139570392;
            const double_t protonMass = 0.938272089;
            const double_t lambdaMass = 1.115683;

            std::cout << "Running...\n";
            VertexingUtils::FCCAnalysesVertex primary;
            int count = 0;
            for (auto &v: vertex) {
                //find the primary vertex
                if (v.vertex.primary) {
                    primary = v;
                    count++;
                }
            }
            int counterP = 0;
            for (auto &p: vertex) {
                if (p.vertex.primary || p.ntracks != 2) {
                    counterP += 1;
                    continue;
                }
                int chargeMuons = 0;
                int numMuons = 0;
                for (auto &potentialMuon: p.reco_ind) {
                    if (recop.at(potentialMuon).type == 13) {
                        numMuons += 1;
                        chargeMuons += recop.at(potentialMuon).charge;
                    }
                }
                if (chargeMuons == 0 && numMuons == 2) {
                    TLorentzVector dimuonMomentum = build_tlv(recop, p.reco_ind);

                    FCCAnalysesComposite2 muon1;
                    muon1.particle = ReconstructedParticle::get_tlv(recop[p.reco_ind.at(0)]);
                    muon1.charge = recop[p.reco_ind.at(0)].charge;

                    FCCAnalysesComposite2 muon2;
                    muon2.particle = ReconstructedParticle::get_tlv(recop[p.reco_ind.at(1)]);
                    muon2.charge = recop[p.reco_ind.at(1)].charge;

                    FCCAnalysesComposite2 dimuon;
                    dimuon.vertex = counterP;
                    dimuon.particle = dimuonMomentum;
                    dimuon.charge = 0;

                    int counterQ = 0;
                    for (auto &q: vertex) {
                        if (q.vertex.primary || q.ntracks != 2) {
                            counterQ += 1;
                            continue;
                        }
                        int chargeHadron1 = recop.at(q.reco_ind.at(0)).charge;
                        int chargeHadron2 = recop.at(q.reco_ind.at(1)).charge;
                        int numHadrons = 0;
                        for (auto &potentialHadron: q.reco_ind) {
                            if (recop.at(potentialHadron).type != 13 && recop.at(potentialHadron).type != 11) {
                                numHadrons += 1;
                            }
                        }
                        if (chargeHadron1 != 0 && chargeHadron1 + chargeHadron2 == 0 && numHadrons == 2) {
                            FCCAnalysesComposite2 hadron1;
                            // hadron1.particle = ReconstructedParticle::get_tlv(recop[q.reco_ind.at(0)]);
                            hadron1.charge = recop[q.reco_ind.at(0)].charge;
                            int hadron1Type;

                            FCCAnalysesComposite2 hadron2;
                            // hadron2.particle = ReconstructedParticle::get_tlv(recop[q.reco_ind.at(1)]);
                            hadron2.charge = recop[q.reco_ind.at(1)].charge;
                            int hadron2Type;

                            TLorentzVector tryPion1 = ReconstructedParticle::get_tlv(recop[q.reco_ind.at(0)]);
                            TLorentzVector tryPion2 = ReconstructedParticle::get_tlv(recop[q.reco_ind.at(1)]);

                            tryPion1.SetE(std::sqrt(pionMass * pionMass + tryPion1.Vect().Dot(tryPion1.Vect())));
                            tryPion2.SetE(std::sqrt(pionMass * pionMass + tryPion2.Vect().Dot(tryPion2.Vect())));
                            double_t norm = (tryPion1 + tryPion2).Mag();
                            if (norm >= 0.49 && norm <= 0.505) {
                                counterQ += 1;
                                continue;
                            }

                            TLorentzVector tryProton1 = ReconstructedParticle::get_tlv(recop[q.reco_ind.at(0)]);
                            TLorentzVector tryProton2 = ReconstructedParticle::get_tlv(recop[q.reco_ind.at(1)]);

                            tryProton1.SetE(std::sqrt(protonMass * protonMass + tryProton1.Vect().Dot(tryProton1.Vect())));
                            tryProton2.SetE(std::sqrt(protonMass * protonMass + tryProton2.Vect().Dot(tryProton2.Vect())));
                            double_t normProtonPion = (tryProton1 + tryPion2).Mag();
                            double_t normPionProton = (tryPion1 + tryProton2).Mag();
                            if (std::abs(normProtonPion - lambdaMass) < std::abs(normPionProton - lambdaMass)) {
                                hadron1.particle = tryProton1;
                                hadron2.particle = tryPion2;
                                hadron1Type = 2212;
                                hadron2Type = 211;
                            } else {
                                hadron1.particle = tryPion1;
                                hadron2.particle = tryProton2;
                                hadron1Type = 211;
                                hadron2Type = 2212;
                            }

                            if (std::abs((hadron1.particle + hadron2.particle).Mag() - lambdaMass) > 0.1) {
                                counterQ += 1;
                                continue;
                            }
                            TLorentzVector dihadronMomentum = hadron1.particle + hadron2.particle;
                            TVector3 dihadronDisplacement = TVector3(q.vertex.position[0] - p.vertex.position[0],
                                                                     q.vertex.position[1] - p.vertex.position[1],
                                                                     q.vertex.position[2] - p.vertex.position[2]);

                            TVector3 dihadron3Momentum = dihadronMomentum.Vect();

                            double_t flightDistanceL = dihadronDisplacement.Mag();
                            double_t muonHadronImpactParam = (dihadronDisplacement - (
                                                   dihadronDisplacement.Dot(dihadron3Momentum) / dihadron3Momentum.Dot(
                                                       dihadron3Momentum)) * dihadron3Momentum).Mag();

                            TLorentzVector LbMomentum = dimuonMomentum + dihadronMomentum;
                            TVector3 Lb3Momentum = LbMomentum.Vect();
                            TVector3 displacementTotal = TVector3(p.vertex.position[0] - primary.vertex.position[0],
                                                                  p.vertex.position[1] - primary.vertex.position[1],
                                                                  p.vertex.position[2] - primary.vertex.position[2]);

                            double_t flightDistanceLb = displacementTotal.Mag();
                            double_t primaryLbImpactParam = (displacementTotal - (displacementTotal.Dot(Lb3Momentum) / Lb3Momentum.Dot(Lb3Momentum)) * Lb3Momentum).Mag();
                            double_t displacementProduct = displacementTotal.Dot(dihadronDisplacement);

                            FCCAnalysesComposite2 dihadron;
                            // dihadron.vertex = counterQ;
                            dihadron.particle = dihadronMomentum;
                            dihadron.charge = 0;

                            FCCAnalysesComposite2 total;
                            total.particle = LbMomentum;
                            total.charge = 0;
                            // total.vertex = 0;

                            if (primaryLbImpactParam < initial) {
                                newDimuons.clear();
                                newMuons1.clear();
                                newMuons2.clear();
                                newDihadrons.clear();
                                newHadrons1.clear();
                                newHadrons2.clear();
                                newTotals.clear();
                                newFlightDistancesLb.clear();
                                newFlightDistancesL.clear();
                                newPrimaryLbImpactParams.clear();
                                newMuonHadronImpactParams.clear();
                                newDisplacementProducts.clear();
                                newHadron1Types.clear();
                                newHadron2Types.clear();

                                newDimuons.push_back(dimuon);
                                newMuons1.push_back(muon1);
                                newMuons2.push_back(muon2);

                                newDihadrons.push_back(dihadron);
                                newHadrons1.push_back(hadron1);
                                newHadrons2.push_back(hadron2);
                                newTotals.push_back(total);

                                newFlightDistancesL.push_back(flightDistanceL);
                                newFlightDistancesLb.push_back(flightDistanceLb);

                                newPrimaryLbImpactParams.push_back(primaryLbImpactParam);
                                newMuonHadronImpactParams.push_back(muonHadronImpactParam);
                                newDisplacementProducts.push_back(displacementProduct);

                                newHadron1Types.push_back(hadron1Type);
                                newHadron2Types.push_back(hadron2Type);
                                initial = primaryLbImpactParam;
                            }
                        }
                        counterQ += 1;
                    }
                }
                counterP += 1;
            }
            if (initial == 1e9) return result;

            FCCAnalysesComposite2 compL;
            FCCAnalysesComposite2 compMu;
            FCCAnalysesComposite2 savedL;
            FCCAnalysesComposite2 savedMu;

            //   ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>
            // SelPrimaryTracks(ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind,
            //                  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,
            //                  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc,
            //                  TVector3 MC_EventPrimaryVertex)
            int counterL = 0;
            int pcount = 0;
            double max = 0;
            bool valid = false;
            for (auto &p: vertex) {
                if (p.vertex.primary) {
                    counterL += 1;
                    continue;
                }
                if (p.ntracks != 2) {
                    counterL += 1;
                    continue;
                }
                int charge_p = 0;
                int nobj_p = 0;
                for (auto &r: p.reco_ind) {
                    if (recop.at(r).type == 2212) {
                        nobj_p += 1;
                        charge_p += recop.at(r).charge;
                    }
                }
                int charge_pi = 0;
                int nobj_pi = 0;
                for (auto &r: p.reco_ind) {
                    if (recop.at(r).type == 211) {
                        nobj_pi += 1;
                        charge_pi += recop.at(r).charge;
                    }
                }
                if (nobj_pi == 1 && nobj_p == 1 && charge_pi + charge_p == 0) {
                    compL.vertex = counterL;
                    compL.particle = build_tlv(recop, p.reco_ind);
                    compL.charge = charge_pi + charge_p;
                    pcount += 1;
                    int counterMu = 0;
                    for (auto &q: vertex) {
                        if (q.vertex.primary) {
                            counterMu += 1;
                            continue;
                        }
                        if (q.ntracks != 2) {
                            counterMu += 1;
                            continue;
                        }
                        int charge_mumu = 0;
                        int nobj_mumu = 0;
                        for (auto &s: q.reco_ind) {
                            if (recop.at(s).type == 13) {
                                nobj_mumu += 1;
                                charge_mumu += recop.at(s).charge;
                            }
                        }
                        if (nobj_mumu == 2 && charge_mumu == 0) {
                            compMu.vertex = counterMu;
                            compMu.particle = build_tlv(recop, q.reco_ind);
                            compMu.charge = charge_pi + charge_p;

                            double alignment = compMu.particle[0] * compL.particle[0] + compMu.particle[1] * compL.
                                               particle[1] + compMu.particle[2] * compL.particle[2];
                            if (alignment > max) {
                                max = alignment;
                                savedL = compL;
                                savedMu = compMu;
                                valid = true;
                            }
                        }
                        counterMu += 1;
                    }
                }
                counterL += 1;
            }

            if (valid) {
                TLorentzVector dimuonMomentum = savedMu.particle;
                TLorentzVector dihadronMomentum = savedL.particle;
                TVector3 dihadronDisplacement = TVector3(vertex.at(savedL.vertex).vertex.position[0] - vertex.at(savedMu.vertex).vertex.position[0],
                                                         vertex.at(savedL.vertex).vertex.position[1] - vertex.at(savedMu.vertex).vertex.position[1],
                                                         vertex.at(savedL.vertex).vertex.position[2] - vertex.at(savedMu.vertex).vertex.position[2]);

                TVector3 dihadron3Momentum = TVector3(dihadronMomentum.Px(), dihadronMomentum.Py(),
                                                      dihadronMomentum.Pz());

                double_t flightDistanceL = dihadronDisplacement.Mag();
                double_t muonHadronImpactParam = (dihadronDisplacement - (
                                       dihadronDisplacement.Dot(dihadron3Momentum) / dihadron3Momentum.Dot(
                                           dihadron3Momentum)) * dihadron3Momentum).Mag();

                TLorentzVector LbMomentum = dimuonMomentum + dihadronMomentum;
                TVector3 Lb3Momentum = TVector3(LbMomentum.Px(), LbMomentum.Py(), LbMomentum.Pz());
                TVector3 displacementTotal = TVector3(primary.vertex.position[0] - vertex.at(savedMu.vertex).vertex.position[0],
                                                 primary.vertex.position[1] - vertex.at(savedMu.vertex).vertex.position[1],
                                                 primary.vertex.position[2] - vertex.at(savedMu.vertex).vertex.position[2]);

                double_t flightDistanceLb = displacementTotal.Mag();
                double_t primaryLbImpactParam = (displacementTotal - (displacementTotal.Dot(Lb3Momentum) / Lb3Momentum.Dot(Lb3Momentum)) * Lb3Momentum).Mag();

                double_t displacementProduct = displacementTotal.Dot(dihadronDisplacement);

                FCCAnalysesComposite2 total;
                total.particle = LbMomentum;
                total.charge = 0;
                total.vertex = 0;

                result.newDimuons = newDimuons;
                result.newMuons1 = newMuons1;
                result.newMuons2 = newMuons2;

                result.newDihadrons = newDihadrons;
                result.newHadrons1 = newHadrons1;
                result.newHadrons2 = newHadrons2;

                result.newTotals = newTotals;

                result.newFlightDistancesLb = newFlightDistancesLb;
                result.newFlightDistancesL = newFlightDistancesL;

                result.newPrimaryLbImpactParameters = newPrimaryLbImpactParams;
                result.newMuonHadronImpactParameters = newMuonHadronImpactParams;
                result.newDisplacementProduct = newDisplacementProducts;

                result.newHadron1Types = newHadron1Types;
                result.newHadron2Types = newHadron2Types;

                result.oldDihadrons.push_back(savedL);
                FCCAnalysesComposite2 hadron1;
                hadron1.particle = ReconstructedParticle::get_tlv(recop[vertex.at(savedL.vertex).reco_ind.at(0)]);
                hadron1.charge = recop[vertex.at(savedL.vertex).reco_ind.at(0)].charge;
                int hadron1Type = recop[vertex.at(savedL.vertex).reco_ind.at(0)].type;
                FCCAnalysesComposite2 hadron2;
                hadron2.particle = ReconstructedParticle::get_tlv(recop[vertex.at(savedL.vertex).reco_ind.at(1)]);
                hadron2.charge = recop[vertex.at(savedL.vertex).reco_ind.at(1)].charge;
                int hadron2Type = recop[vertex.at(savedL.vertex).reco_ind.at(1)].type;
                result.oldHadrons1.push_back(hadron1);
                result.oldHadrons2.push_back(hadron2);
                result.oldHadron1Types.push_back(hadron1Type);
                result.oldHadron2Types.push_back(hadron2Type);

                result.oldDimuons.push_back(savedMu);
                FCCAnalysesComposite2 muon1;
                muon1.particle = ReconstructedParticle::get_tlv(recop[vertex.at(savedMu.vertex).reco_ind.at(0)]);
                muon1.charge = recop[vertex.at(savedMu.vertex).reco_ind.at(0)].charge;
                FCCAnalysesComposite2 muon2;
                muon2.particle = ReconstructedParticle::get_tlv(recop[vertex.at(savedMu.vertex).reco_ind.at(1)]);
                muon2.charge = recop[vertex.at(savedMu.vertex).reco_ind.at(1)].charge;
                result.oldMuons1.push_back(muon1);
                result.oldMuons2.push_back(muon2);

                result.oldPrimaryLbImpactParameters.push_back(primaryLbImpactParam);
                result.oldTotals.push_back(total);
                result.oldMuonHadronImpactParameters.push_back(muonHadronImpactParam);
                result.oldDisplacementProduct.push_back(displacementProduct);
                result.oldFlightDistancesL.push_back(flightDistanceL);
                result.oldFlightDistancesLb.push_back(flightDistanceLb);
            }

            return result;


        }



        buildLb2LMuMu Lb2LMuMu(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
                               ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop) {

            buildLb2LMuMu result;

            double initial = 1e9;

            ROOT::VecOps::RVec<FCCAnalysesComposite2> dimuons; //dimuon momentum
            ROOT::VecOps::RVec<FCCAnalysesComposite2> muons1; //muon1 momentum
            ROOT::VecOps::RVec<FCCAnalysesComposite2> muons2; //muon2 momentum

            ROOT::VecOps::RVec<FCCAnalysesComposite2> dihadrons; //dihadron momentum
            ROOT::VecOps::RVec<FCCAnalysesComposite2> hadrons1; //hadron1 momentum
            ROOT::VecOps::RVec<FCCAnalysesComposite2> hadrons2; //hadron2 momentum

            ROOT::VecOps::RVec<FCCAnalysesComposite2> totals;
            //total momentum (dimuon + dihadron = lambdaB for valid selection)

            ROOT::VecOps::RVec<double_t> flightDistancesLb; //distance from dimuon vertex to primary vertex
            ROOT::VecOps::RVec<double_t> flightDistancesL; //distance from dihadron vertex to dimuon vertex

            ROOT::VecOps::RVec<double_t> primaryLbImpactParams; //minimal distance between reconstructed lambdaB path and primary vertex
            ROOT::VecOps::RVec<double_t> muonHadronImpactParams; //minimal distance between dihadron path and dimuon vertex
            ROOT::VecOps::RVec<double_t> displacementProducts;

            ROOT::VecOps::RVec<int> hadron1Types;
            ROOT::VecOps::RVec<int> hadron2Types;
            //currently storing type of hadrons to use for debugging, will not be needed later


            const double_t pionMass = 0.139570392;
            const double_t protonMass = 0.938272089;
            const double_t lambdaMass = 1.115683;

            std::cout << "Running...\n";
            VertexingUtils::FCCAnalysesVertex primary;
            int count = 0;
            for (auto &v: vertex) {
                //find the primary vertex
                if (v.vertex.primary) {
                    primary = v;
                    count++;
                }
            }
            int counterP = 0;
            for (auto &p: vertex) {
                if (p.vertex.primary || p.ntracks != 2) {
                    counterP += 1;
                    continue;
                }
                int chargeMuons = 0;
                int numMuons = 0;
                for (auto &potentialMuon: p.reco_ind) {
                    if (recop.at(potentialMuon).type == 13) {
                        numMuons += 1;
                        chargeMuons += recop.at(potentialMuon).charge;
                    }
                }
                if (chargeMuons == 0 && numMuons == 2) {
                    TLorentzVector dimuonMomentum = build_tlv(recop, p.reco_ind);

                    FCCAnalysesComposite2 muon1;
                    muon1.particle = ReconstructedParticle::get_tlv(recop[p.reco_ind.at(0)]);
                    muon1.charge = recop[p.reco_ind.at(0)].charge;

                    FCCAnalysesComposite2 muon2;
                    muon2.particle = ReconstructedParticle::get_tlv(recop[p.reco_ind.at(1)]);
                    muon2.charge = recop[p.reco_ind.at(1)].charge;

                    FCCAnalysesComposite2 dimuon;
                    dimuon.vertex = counterP;
                    dimuon.particle = dimuonMomentum;
                    dimuon.charge = 0;

                    int counterQ = 0;
                    for (auto &q: vertex) {
                        if (q.vertex.primary || q.ntracks != 2) {
                            counterQ += 1;
                            continue;
                        }
                        int chargeHadron1 = recop.at(q.reco_ind.at(0)).charge;
                        int chargeHadron2 = recop.at(q.reco_ind.at(1)).charge;
                        int numHadrons = 0;
                        for (auto &potentialHadron: q.reco_ind) {
                            if (recop.at(potentialHadron).type != 13 && recop.at(potentialHadron).type != 11) {
                                numHadrons += 1;
                            }
                        }
                        if (chargeHadron1 != 0 && chargeHadron1 + chargeHadron2 == 0 && numHadrons == 2) {
                            FCCAnalysesComposite2 hadron1;
                            hadron1.particle = ReconstructedParticle::get_tlv(recop[q.reco_ind.at(0)]);
                            hadron1.charge = recop[q.reco_ind.at(0)].charge;
                            int hadron1Type;

                            FCCAnalysesComposite2 hadron2;
                            hadron2.particle = ReconstructedParticle::get_tlv(recop[q.reco_ind.at(1)]);
                            hadron2.charge = recop[q.reco_ind.at(1)].charge;
                            int hadron2Type;

                            TLorentzVector tryPion1 = ReconstructedParticle::get_tlv(recop[q.reco_ind.at(0)]);
                            TLorentzVector tryPion2 = ReconstructedParticle::get_tlv(recop[q.reco_ind.at(1)]);

                            tryPion1.SetE(std::sqrt(pionMass * pionMass + tryPion1.Vect().Dot(tryPion1.Vect())));
                            tryPion2.SetE(std::sqrt(pionMass * pionMass + tryPion2.Vect().Dot(tryPion2.Vect())));
                            double_t norm = (tryPion1 + tryPion2).Mag();
                            if (norm >= 0.49 && norm <= 0.505) {
                                counterQ += 1;
                                continue;
                            }

                            TLorentzVector tryProton1 = ReconstructedParticle::get_tlv(recop[q.reco_ind.at(0)]);
                            TLorentzVector tryProton2 = ReconstructedParticle::get_tlv(recop[q.reco_ind.at(1)]);

                            tryPion1.SetE(std::sqrt(protonMass * protonMass + tryPion1.Vect().Dot(tryPion1.Vect())));
                            tryPion2.SetE(std::sqrt(protonMass * protonMass + tryPion2.Vect().Dot(tryPion2.Vect())));
                            double_t normProtonPion = (tryProton1 + tryPion2).Mag();
                            double_t normPionProton = (tryPion1 + tryProton2).Mag();
                            if (std::abs(normProtonPion - lambdaMass) < std::abs(normPionProton - lambdaMass)) {
                                hadron1.particle = tryProton1;
                                hadron2.particle = tryPion2;
                                hadron1Type = 2212;
                                hadron2Type = 211;
                            } else {
                                hadron1.particle = tryPion1;
                                hadron2.particle = tryProton2;
                                hadron1Type = 211;
                                hadron2Type = 2212;
                            }
                            TLorentzVector dihadronMomentum = hadron1.particle + hadron2.particle;
                            TVector3 dihadronDisplacement = TVector3(q.vertex.position[0] - p.vertex.position[0],
                                                                     q.vertex.position[1] - p.vertex.position[1],
                                                                     q.vertex.position[2] - p.vertex.position[2]);

                            TVector3 dihadron3Momentum = dihadronMomentum.Vect();

                            double_t flightDistanceL = dihadronDisplacement.Mag();
                            double_t muonHadronImpactParam = (dihadronDisplacement - (
                                                   dihadronDisplacement.Dot(dihadron3Momentum) / dihadron3Momentum.Dot(
                                                       dihadron3Momentum)) * dihadron3Momentum).Mag();

                            TLorentzVector LbMomentum = dimuonMomentum + dihadronMomentum;
                            TVector3 Lb3Momentum = LbMomentum.Vect();
                            TVector3 displacementTotal = TVector3(p.vertex.position[0] - primary.vertex.position[0],
                                                                  p.vertex.position[1] - primary.vertex.position[1],
                                                                  p.vertex.position[2] - primary.vertex.position[2]);

                            double_t flightDistanceLb = displacementTotal.Mag();
                            double_t primaryLbImpactParam = (displacementTotal - (displacementTotal.Dot(Lb3Momentum) / Lb3Momentum.Dot(Lb3Momentum)) * Lb3Momentum).Mag();
                            double_t displacementProduct = displacementTotal.Dot(dihadronDisplacement);

                            FCCAnalysesComposite2 dihadron;
                            dihadron.vertex = counterQ;
                            dihadron.particle = dihadronMomentum;
                            dihadron.charge = 0;

                            FCCAnalysesComposite2 total;
                            total.particle = LbMomentum;
                            total.charge = 0;
                            total.vertex = 0;

                            if (primaryLbImpactParam < initial) {
                                dimuons.clear();
                                muons1.clear();
                                muons2.clear();
                                dihadrons.clear();
                                hadrons1.clear();
                                hadrons2.clear();
                                totals.clear();
                                flightDistancesLb.clear();
                                flightDistancesL.clear();
                                primaryLbImpactParams.clear();
                                muonHadronImpactParams.clear();
                                displacementProducts.clear();
                                hadron1Types.clear();
                                hadron2Types.clear();

                                dimuons.push_back(dimuon);
                                muons1.push_back(muon1);
                                muons2.push_back(muon2);

                                dihadrons.push_back(dihadron);
                                hadrons1.push_back(hadron1);
                                hadrons2.push_back(hadron2);
                                totals.push_back(total);

                                flightDistancesL.push_back(flightDistanceL);
                                flightDistancesLb.push_back(flightDistanceLb);

                                primaryLbImpactParams.push_back(primaryLbImpactParam);
                                muonHadronImpactParams.push_back(muonHadronImpactParam);
                                displacementProducts.push_back(displacementProduct);

                                hadron1Types.push_back(hadron1Type);
                                hadron2Types.push_back(hadron2Type);
                                initial = primaryLbImpactParam;
                            }
                        }
                        counterQ += 1;
                    }
                }
                counterP += 1;
            }
            result.dimuons = dimuons;
            result.muons1 = muons1;
            result.muons2 = muons2;

            result.dihadrons = dihadrons;
            result.hadrons1 = hadrons1;
            result.hadrons2 = hadrons2;

            result.totals = totals;

            result.flightDistancesLb = flightDistancesLb;
            result.flightDistancesL = flightDistancesL;

            result.primaryLbImpactParameters = primaryLbImpactParams;
            result.muonHadronImpactParameters = muonHadronImpactParams;
            result.displacementProduct = displacementProducts;

            result.hadron1Types = hadron1Types;
            result.hadron2Types = hadron2Types;

            return result;
        }


        std::pair<ROOT::VecOps::RVec<FCCAnalysesComposite2>, ROOT::VecOps::RVec<FCCAnalysesComposite2> > build_Lb2LMuMu(
            ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
            ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop, bool clean) {
            //Bd2KstNuNu
            ROOT::VecOps::RVec<FCCAnalysesComposite2> resultL;
            ROOT::VecOps::RVec<FCCAnalysesComposite2> resultMu;

            int events = 0;

            FCCAnalysesComposite2 compL;
            FCCAnalysesComposite2 compMu;
            FCCAnalysesComposite2 savedL;
            FCCAnalysesComposite2 savedMu;

            //   ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>
            // SelPrimaryTracks(ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind,
            //                  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,
            //                  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc,
            //                  TVector3 MC_EventPrimaryVertex)
            int counterL = 0;
            int pcount = 0;
            double max = 0;
            bool valid = false;
            for (auto &p: vertex) {
                if (p.vertex.primary) {
                    counterL += 1;
                    continue;
                }
                if (p.ntracks != 2) {
                    counterL += 1;
                    continue;
                }
                int charge_p = 0;
                int nobj_p = 0;
                for (auto &r: p.reco_ind) {
                    if (recop.at(r).type == 2212) {
                        nobj_p += 1;
                        charge_p += recop.at(r).charge;
                    }
                }
                int charge_pi = 0;
                int nobj_pi = 0;
                for (auto &r: p.reco_ind) {
                    if (recop.at(r).type == 211) {
                        nobj_pi += 1;
                        charge_pi += recop.at(r).charge;
                    }
                }
                if (nobj_pi == 1 && nobj_p == 1 && charge_pi + charge_p == 0) {
                    compL.vertex = counterL;
                    compL.particle = build_tlv(recop, p.reco_ind);
                    compL.charge = charge_pi + charge_p;
                    pcount += 1;
                    int counterMu = 0;
                    for (auto &q: vertex) {
                        if (q.vertex.primary) {
                            counterMu += 1;
                            continue;
                        }
                        if (q.ntracks != 2) {
                            counterMu += 1;
                            continue;
                        }
                        int charge_mumu = 0;
                        int nobj_mumu = 0;
                        for (auto &s: q.reco_ind) {
                            if (recop.at(s).type == 13) {
                                nobj_mumu += 1;
                                charge_mumu += recop.at(s).charge;
                            }
                        }
                        if (nobj_mumu == 2 && charge_mumu == 0) {
                            compMu.vertex = counterMu;
                            compMu.particle = build_tlv(recop, q.reco_ind);
                            compMu.charge = charge_pi + charge_p;

                            double alignment = compMu.particle[0] * compL.particle[0] + compMu.particle[1] * compL.
                                               particle[1] + compMu.particle[2] * compL.particle[2];
                            if (alignment > max) {
                                max = alignment;
                                savedL = compL;
                                savedMu = compMu;
                                valid = true;
                            }
                        }
                        counterMu += 1;
                    }
                }
                counterL += 1;
            }
            if (valid) {
                if (clean) {
                    if (pcount == 1) {
                        resultL.push_back(savedL);
                        resultMu.push_back(savedMu);
                    }
                } else {
                    resultL.push_back(savedL);
                    resultMu.push_back(savedMu);
                }
            }

            return std::make_pair(resultL, resultMu);
        }


        float get_d0(TVector3 x, TVector3 p) {
            float D = x[1] * p[0] / sqrt(p[0] * p[0] + p[1] * p[1]) - x[0] * p[1] / sqrt(p[0] * p[0] + p[1] * p[1]);
            return D;
        }

        float get_z0(TVector3 x, TVector3 p) {
            float Z = x[2] - p[2] / sqrt(p[0] * p[0] + p[1] * p[1]) * sqrt(
                          x[0] * x[0] + x[1] * x[1] - pow(get_d0(x, p), 2));
            return Z;
        }


        int get_Npos(ROOT::VecOps::RVec<float> in) {
            int result = 0;
            for (auto &p: in)
                if (p > 0.)result += 1;
            return result;
        }

        int get_Nneg(ROOT::VecOps::RVec<float> in) {
            int result = 0;
            for (auto &p: in)
                if (p < 0.)result += 1;
            return result;
        }

        float get_dPV2DV_min(ROOT::VecOps::RVec<float> in) {
            float result = 9999999.;
            for (auto &p: in) {
                if (p < result && p > 0)result = p;
            }
            return result;
        }

        float get_dPV2DV_max(ROOT::VecOps::RVec<float> in) {
            float result = -99999;
            for (auto &p: in) {
                if (p > result)result = p;
            }
            return result;
        }

        float get_dPV2DV_ave(ROOT::VecOps::RVec<float> in) {
            float result = 0.;
            for (auto &p: in) {
                result += p;
            }
            if (in.size() - 1 > 0)result = result / (in.size() - 1);
            return result;
        }

        int get_PV_ntracks(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex) {
            int result = 0;
            for (auto &p: vertex) {
                if (p.vertex.primary == 1) {
                    result = p.ntracks;
                    break;
                }
            }
            return result;
        }

        int hasPV(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex) {
            int result = 0;
            for (auto &p: vertex) {
                if (p.vertex.primary == 1) {
                    result = 1;
                    break;
                }
            }
            return result;
        }

        ROOT::VecOps::RVec<float> get_Vertex_mass(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
                                                  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco) {
            ROOT::VecOps::RVec<float> result;
            for (auto &p: vertex) {
                ROOT::VecOps::RVec<int> reco_ind = p.reco_ind;
                float mass = build_invmass(reco, reco_ind);
                result.push_back(mass);
            }
            return result;
        }

        ROOT::VecOps::RVec<float> get_Vertex_x(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex) {
            ROOT::VecOps::RVec<float> result;
            for (auto &p: vertex)
                result.push_back(p.vertex.position.x);
            return result;
        }

        ROOT::VecOps::RVec<float> get_Vertex_y(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex) {
            ROOT::VecOps::RVec<float> result;
            for (auto &p: vertex)
                result.push_back(p.vertex.position.y);
            return result;
        }

        ROOT::VecOps::RVec<float> get_Vertex_z(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex) {
            ROOT::VecOps::RVec<float> result;
            for (auto &p: vertex)
                result.push_back(p.vertex.position.z);
            return result;
        }

        ROOT::VecOps::RVec<float> get_Vertex_xErr(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex) {
            ROOT::VecOps::RVec<float> result;
            for (auto &p: vertex)
                result.push_back(sqrt(p.vertex.covMatrix[0]));
            return result;
        }

        ROOT::VecOps::RVec<float> get_Vertex_yErr(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex) {
            ROOT::VecOps::RVec<float> result;
            for (auto &p: vertex)
                result.push_back(sqrt(p.vertex.covMatrix[2]));
            return result;
        }

        ROOT::VecOps::RVec<float> get_Vertex_zErr(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex) {
            ROOT::VecOps::RVec<float> result;
            for (auto &p: vertex)
                result.push_back(sqrt(p.vertex.covMatrix[5]));
            return result;
        }

        ROOT::VecOps::RVec<float> get_Vertex_chi2(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex) {
            ROOT::VecOps::RVec<float> result;
            for (auto &p: vertex)
                result.push_back(p.vertex.chi2);
            return result;
        }

        ROOT::VecOps::RVec<int> get_Vertex_isPV(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex) {
            ROOT::VecOps::RVec<int> result;
            for (auto &p: vertex)
                result.push_back(p.vertex.primary);
            return result;
        }

        ROOT::VecOps::RVec<int> get_Vertex_ntracks(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex) {
            ROOT::VecOps::RVec<int> result;
            for (auto &p: vertex)
                result.push_back(p.ntracks);
            return result;
        }

        ROOT::VecOps::RVec<float> get_Vertex_d2PV(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
                                                  int comp) {
            ROOT::VecOps::RVec<float> result;
            VertexingUtils::FCCAnalysesVertex PV;
            for (auto &p: vertex)
                if (p.vertex.primary > 0) PV = p;

            for (auto &p: vertex) {
                if (p.vertex.primary > 0) result.push_back(0);
                else result.push_back(get_distanceVertex(PV.vertex, p.vertex, comp));
            }
            return result;
        }


        ROOT::VecOps::RVec<float> get_Vertex_d2PVError(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
                                                       int comp) {
            ROOT::VecOps::RVec<float> result;
            VertexingUtils::FCCAnalysesVertex PV;
            for (auto &p: vertex)
                if (p.vertex.primary > 0) PV = p;

            for (auto &p: vertex) {
                if (p.vertex.primary > 0) result.push_back(0);
                else result.push_back(get_distanceErrorVertex(PV.vertex, p.vertex, comp));
            }
            return result;
        }


        ROOT::VecOps::RVec<float> get_Vertex_d2MC(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
                                                  ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> mcver,
                                                  ROOT::VecOps::RVec<int> mcind,
                                                  int comp) {
            ROOT::VecOps::RVec<float> result;

            for (size_t i = 0; i < vertex.size(); ++i) {
                edm4hep::Vector3f recov = vertex.at(i).vertex.position;
                TVector3 mcv = mcver.at(mcind.at(i)).vertex;
                result.push_back(get_distance(recov, mcv, comp));
            }

            return result;
        }

        ROOT::VecOps::RVec<int> get_Vertex_indMC(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex) {
            ROOT::VecOps::RVec<int> result;
            for (size_t i = 0; i < vertex.size(); ++i)
                result.push_back(vertex.at(i).mc_ind);
            return result;
        }

        ROOT::VecOps::RVec<int> get_Vertex_indMC(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
                                                 ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> mcver) {
            ROOT::VecOps::RVec<int> result;
            for (auto &p: vertex) {
                float distance = 99999999999.;
                int index = -1;
                for (size_t i = 0; i < mcver.size(); ++i) {
                    float distance2 = sqrt(pow(p.vertex.position.x - mcver.at(i).vertex[0], 2) +
                                           pow(p.vertex.position.y - mcver.at(i).vertex[1], 2) +
                                           pow(p.vertex.position.z - mcver.at(i).vertex[2], 2));

                    if (distance2 < distance) {
                        distance = distance2;
                        index = i;
                    }
                }
                if (index > -1)result.push_back(index);
                else {
                    std::cout << "problem index get_Vertex_indMC " << index << "  distance " << distance << std::endl;
                    std::cout << "reco vtx chi2 " << p.vertex.chi2 << " ntrk " << p.ntracks << "x,y,z=" << p.vertex.
                            position.x << ", " << p.vertex.position.y << ", " << p.vertex.position.z << std::endl;
                    for (size_t i = 0; i < mcver.size(); ++i)
                        std::cout << "  mc vtx i=" << i << " ntrk " << mcver.at(i).mc_ind.size() << "x,y,z=" << mcver.
                                at(i).vertex[0] << ", " << mcver.at(i).vertex[1] << ", " << mcver.at(i).vertex[2] <<
                                std::endl;
                }
            }
            return result;
        }


        ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex>
        get_VertexObject(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> mcver,
                         ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,
                         ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
                         ROOT::VecOps::RVec<int> recin,
                         ROOT::VecOps::RVec<int> mcin) {
            ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> result;

            ROOT::VecOps::RVec<ROOT::VecOps::RVec<int> > rp2mc = ReconstructedParticle2MC::getRP2MC_indexVec(
                recin, mcin, reco);

            //for (auto &p: mcver){
            for (size_t v = 0; v < mcver.size(); ++v) {
                ROOT::VecOps::RVec<int> mc_indRVec = mcver.at(v).mc_ind;
                std::vector<int> mc_ind;

                for (size_t i = 0; i < mc_indRVec.size(); ++i)mc_ind.push_back(mc_indRVec.at(i));

                ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recoparticles;
                for (size_t i = 0; i < rp2mc.size(); ++i) {
                    for (size_t j = 0; j < rp2mc.at(i).size(); ++j) {
                        std::vector<int>::iterator it = std::find(mc_ind.begin(), mc_ind.end(), rp2mc.at(i).at(j));
                        if (it != mc_ind.end() && fabs(reco.at(i).charge) > 0)recoparticles.push_back(reco.at(i));
                    }
                }

                ROOT::VecOps::RVec<edm4hep::TrackState> tmptracks = ReconstructedParticle2Track::getRP2TRK(
                    recoparticles, tracks);

                if (recoparticles.size() < 2)continue;

                VertexingUtils::FCCAnalysesVertex TheVertex;
                if (v == 0) TheVertex = VertexFitterSimple::VertexFitter(
                                1, recoparticles, tracks, true, 4.5, 20e-3, 300);
                    //if (v==0) TheVertex = VertexFitterSimple::VertexFitter(1,recoparticles, tracks);
                else TheVertex = VertexFitterSimple::VertexFitter(0, recoparticles, tracks);
                if (std::isnan(TheVertex.vertex.chi2))continue;
                TheVertex.mc_ind = v;
                result.push_back(TheVertex);
            }
            return result;
        }

        int globalmm = 0;

        ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex>
        merge_VertexObjet(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> in) {
            ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> result;
            std::cout << "============================" << std::endl;
            for (size_t i = 0; i < in.size() - 1; ++i) {
                edm4hep::VertexData vi = in.at(i).vertex;
                std::array<float, 6> vi_covMatrix = vi.covMatrix;
                for (size_t j = i + 1; j < in.size(); ++j) {
                    edm4hep::VertexData vj = in.at(j).vertex;
                    std::array<float, 6> vj_covMatrix = vj.covMatrix;
                    float dist = get_distanceVertex(vi, vj, -1);
                    float err1 = sqrt(
                        vi_covMatrix[0] + vj_covMatrix[0] + vi_covMatrix[2] + vj_covMatrix[2] + vi_covMatrix[5] +
                        vj_covMatrix[5]);
                    float err2 = get_distanceErrorVertex(vi, vj, -1);
                    if (dist < err1 || dist < err2) {
                        globalmm += 1;
                        std::cout << globalmm << "  i,j=" << i << ", " << j << "  d= " << dist << "  err1  " << err1 <<
                                "  err2  " << err2 << " ntrk v1,v2=" << in.at(i).ntracks << ", " << in.at(j).ntracks <<
                                "  isPV=" << vi.primary << ", " << vj.primary << std::endl;
                    }
                }
            }
            //return result;
            return in;
        }


        std::vector<std::vector<int> > get_Vertex_ind(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex) {
            std::vector<std::vector<int> > result;
            for (auto &p: vertex) {
                std::vector<int> tmp;
                for (size_t i = 0; i < p.reco_ind.size(); ++i) tmp.push_back(p.reco_ind.at(i));
                result.push_back(tmp);
            }
            return result;
        }


        ROOT::VecOps::RVec<TVector3> get_MCVertex(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> vertex) {
            ROOT::VecOps::RVec<TVector3> result;
            for (auto &p: vertex)
                result.push_back(p.vertex);
            return result;
        }

        ROOT::VecOps::RVec<float> get_MCVertex_x(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> vertex) {
            ROOT::VecOps::RVec<float> result;
            for (auto &p: vertex)
                result.push_back(p.vertex[0]);
            return result;
        }

        ROOT::VecOps::RVec<float> get_MCVertex_y(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> vertex) {
            ROOT::VecOps::RVec<float> result;
            for (auto &p: vertex)
                result.push_back(p.vertex[1]);
            return result;
        }

        ROOT::VecOps::RVec<float> get_MCVertex_z(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> vertex) {
            ROOT::VecOps::RVec<float> result;
            for (auto &p: vertex)
                result.push_back(p.vertex[2]);
            return result;
        }

        ROOT::VecOps::RVec<int> get_NTracksMCVertex(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> vertex) {
            ROOT::VecOps::RVec<int> result;
            for (auto &p: vertex)
                result.push_back(p.mc_ind.size());
            return result;
        }

        std::vector<std::vector<int> >
        get_MCindMCVertex(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> vertex) {
            std::vector<std::vector<int> > result;
            for (auto &p: vertex) {
                std::vector<int> tmp;
                for (size_t i = 0; i < p.mc_ind.size(); ++i) tmp.push_back(p.mc_ind.at(i));
                result.push_back(tmp);
            }
            return result;
        }

        std::vector<std::vector<int> > get_MCpdgMCVertex(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> vertex,
                                                         ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) {
            std::vector<std::vector<int> > result;
            for (auto &p: vertex) {
                std::vector<int> tmp;
                for (size_t i = 0; i < p.mc_ind.size(); ++i) tmp.push_back(mc.at(p.mc_ind.at(i)).PDG);
                for (size_t i = 0; i < p.mc_indneutral.size(); ++i) tmp.push_back(mc.at(p.mc_indneutral.at(i)).PDG);
                result.push_back(tmp);
            }
            return result;
        }

        std::vector<std::vector<int> > get_MCpdgMotherMCVertex(
            ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> vertex,
            ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) {
            std::vector<std::vector<int> > result;
            for (auto &p: vertex) {
                std::vector<int> tmp;
                for (size_t i = 0; i < p.mother_ind.size(); ++i) tmp.push_back(mc.at(p.mother_ind.at(i)).PDG);
                result.push_back(tmp);
            }
            return result;
        }

        std::vector<std::vector<int> > get_MCpdgGMotherMCVertex(
            ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> vertex,
            ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) {
            std::vector<std::vector<int> > result;
            for (auto &p: vertex) {
                std::vector<int> tmp;
                for (size_t i = 0; i < p.gmother_ind.size(); ++i) tmp.push_back(mc.at(p.gmother_ind.at(i)).PDG);
                result.push_back(tmp);
            }
            return result;
        }


        ROOT::VecOps::RVec<int> get_MCMother1(ROOT::VecOps::RVec<edm4hep::MCParticleData> mc,
                                              ROOT::VecOps::RVec<int> ind) {
            ROOT::VecOps::RVec<int> result;
            for (size_t i = 0; i < mc.size(); ++i) {
                if (mc.at(i).parents_begin == mc.at(i).parents_end)result.push_back(-1);
                else result.push_back(ind.at(mc.at(i).parents_begin));
            }
            return result;
        }

        ROOT::VecOps::RVec<int> get_MCMother2(ROOT::VecOps::RVec<edm4hep::MCParticleData> mc,
                                              ROOT::VecOps::RVec<int> ind) {
            ROOT::VecOps::RVec<int> result;
            for (size_t i = 0; i < mc.size(); ++i) {
                if (mc.at(i).parents_begin == mc.at(i).parents_end)result.push_back(-1);
                else result.push_back(ind.at(mc.at(i).parents_end - 1));
            }
            return result;
        }

        ROOT::VecOps::RVec<int> get_MCDaughter1(ROOT::VecOps::RVec<edm4hep::MCParticleData> mc,
                                                ROOT::VecOps::RVec<int> ind) {
            ROOT::VecOps::RVec<int> result;
            for (size_t i = 0; i < mc.size(); ++i) {
                if (mc.at(i).daughters_begin == mc.at(i).daughters_end)result.push_back(-1);
                else result.push_back(ind.at(mc.at(i).daughters_begin));
            }
            return result;
        }

        ROOT::VecOps::RVec<int> get_MCDaughter2(ROOT::VecOps::RVec<edm4hep::MCParticleData> mc,
                                                ROOT::VecOps::RVec<int> ind) {
            ROOT::VecOps::RVec<int> result;
            for (size_t i = 0; i < mc.size(); ++i) {
                if (mc.at(i).daughters_begin == mc.at(i).daughters_end)result.push_back(-1);
                else result.push_back(ind.at(mc.at(i).daughters_end - 1));
            }
            return result;
        }


        ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> get_MCVertexObject(
            ROOT::VecOps::RVec<edm4hep::MCParticleData> mc,
            ROOT::VecOps::RVec<int> ind) {
            ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> result;
            ROOT::VecOps::RVec<TVector3> tmpvec;
            ROOT::VecOps::RVec<int> tmpvecint;
            for (size_t i = 0; i < mc.size(); ++i) {
                if (mc.at(i).charge == 0)continue;
                if (mc.at(i).generatorStatus == 1) {
                    TVector3 tmp;
                    tmp[0] = mc.at(i).vertex.x;
                    tmp[1] = mc.at(i).vertex.y;
                    tmp[2] = mc.at(i).vertex.z;
                    tmpvec.push_back(tmp);
                    tmpvecint.push_back(i);
                }
            }

            for (size_t i = 0; i < tmpvec.size(); ++i) {
                bool vertexfound = false;
                TVector3 vertexPos(tmpvec.at(i)[0], tmpvec.at(i)[1], tmpvec.at(i)[2]);

                if (result.size() == 0) {
                    VertexingUtils::FCCAnalysesVertexMC vertex;
                    ROOT::VecOps::RVec<int> ind;
                    ind.push_back(tmpvecint.at(i));
                    vertex.vertex = vertexPos;
                    vertex.mc_ind = ind;
                    result.push_back(vertex);
                } else {
                    for (size_t j = 0; j < result.size(); ++j) {
                        if (get_distance(result.at(j).vertex, vertexPos) < 0.00001) {
                            result.at(j).mc_ind.push_back(tmpvecint.at(i));
                            vertexfound = true;
                            break;
                        }
                    }
                    if (vertexfound == false) {
                        VertexingUtils::FCCAnalysesVertexMC vertex;
                        ROOT::VecOps::RVec<int> ind2;
                        ind2.push_back(tmpvecint.at(i));
                        vertex.vertex = vertexPos;
                        vertex.mc_ind = ind2;
                        result.push_back(vertex);
                    }
                }
            }


            //adding the neutrals
            for (auto &p: result) {
                ROOT::VecOps::RVec<int> mc_indneutral;
                for (size_t i = 0; i < mc.size(); ++i) {
                    if (mc.at(i).charge != 0)continue;
                    if (mc.at(i).generatorStatus == 1) {
                        TVector3 vertexPos(mc.at(i).vertex.x, mc.at(i).vertex.y, mc.at(i).vertex.z);
                        if (get_distance(p.vertex, vertexPos) < 0.000001) mc_indneutral.push_back(i);
                    }
                }
                p.mc_indneutral = mc_indneutral;
            }

            //adding the mother particles
            for (auto &p: result) {
                std::vector<int> mother_ind;

                //std::cout <<"n part="<<p.mc_ind.size()<<"  x=" << p.vertex[0] <<"  y=" << p.vertex[1] <<"  z=" << p.vertex[2] <<std::endl;
                ROOT::VecOps::RVec<int> mc_ind = p.mc_ind;
                for (size_t i = 0; i < mc_ind.size(); ++i) {
                    //std::cout << "i="<< i << "  mc_ind  "<< mc_ind.at(i) <<"  parent begin " << mc.at(mc_ind.at(i)).parents_begin <<"  parent end " << mc.at(mc_ind.at(i)).parents_end << " end-beg  "  << mc.at(mc_ind.at(i)).parents_end-mc.at(mc_ind.at(i)).parents_begin<< std::endl;

                    for (size_t j = mc.at(mc_ind.at(i)).parents_begin; j < mc.at(mc_ind.at(i)).parents_end; ++j) {
                        //std::cout << "     j="<<j << "  index  " << ind.at(j) << "  PDG ID " << mc.at(ind.at(j)).PDG << std::endl;
                        std::vector<int>::iterator it = std::find(mother_ind.begin(), mother_ind.end(), ind.at(j));
                        if (it == mother_ind.end())mother_ind.push_back(ind.at(j));
                    }
                }

                ROOT::VecOps::RVec<int> mother_indRVec;
                for (size_t i = 0; i < mother_ind.size(); ++i)mother_indRVec.push_back(mother_ind.at(i));
                p.mother_ind = mother_indRVec;
                //std::cout << "found n mothers="<<mother_ind.size()<<std::endl;
            }


            //adding the grand mother particles
            for (auto &p: result) {
                std::vector<int> gmother_ind;

                //std::cout <<"n part="<<p.mc_ind.size()<<"  x=" << p.vertex[0] <<"  y=" << p.vertex[1] <<"  z=" << p.vertex[2] <<std::endl;
                ROOT::VecOps::RVec<int> mc_ind = p.mother_ind;
                for (size_t i = 0; i < mc_ind.size(); ++i) {
                    //std::cout << "i="<< i << "  mc_ind  "<< mc_ind.at(i) <<"  parent begin " << mc.at(mc_ind.at(i)).parents_begin <<"  parent end " << mc.at(mc_ind.at(i)).parents_end << " end-beg  "  << mc.at(mc_ind.at(i)).parents_end-mc.at(mc_ind.at(i)).parents_begin<< std::endl;

                    for (size_t j = mc.at(mc_ind.at(i)).parents_begin; j < mc.at(mc_ind.at(i)).parents_end; ++j) {
                        //std::cout << "     j="<<j << "  index  " << ind.at(j) << "  PDG ID " << mc.at(ind.at(j)).PDG << std::endl;
                        std::vector<int>::iterator it = std::find(gmother_ind.begin(), gmother_ind.end(), ind.at(j));
                        if (it == gmother_ind.end())gmother_ind.push_back(ind.at(j));
                    }
                }

                ROOT::VecOps::RVec<int> gmother_indRVec;
                for (size_t i = 0; i < gmother_ind.size(); ++i)gmother_indRVec.push_back(gmother_ind.at(i));
                p.gmother_ind = gmother_indRVec;
                //std::cout << "found n mothers="<<mother_ind.size()<<std::endl;
            }


            /*std::cout <<"nvx MC  " << result.size()<<std::endl;
            for (size_t j = 0; j < result.size(); ++j)
            std::cout <<"n part="<<result.at(j).mc_ind.size()<<"  x=" << result.at(j).vertex[0] <<"  y=" << result.at(j).vertex[1] <<"  z=" << result.at(j).vertex[2] <<std::endl;*/

            return result;
        }

        float get_distance(edm4hep::Vector3f v1,
                           TVector3 v2,
                           int comp) {
            float result;
            if (comp == 0) result = v1.x - v2[0];
            else if (comp == 1) result = v1.y - v2[1];
            else if (comp == 2) result = v1.z - v2[2];

            else
                result = sqrt(pow(v1.x - v2[0], 2) +
                              pow(v1.y - v2[1], 2) +
                              pow(v1.z - v2[2], 2));
            return result;
        }

        float get_distance(TVector3 v1, TVector3 v2) {
            return sqrt(pow(v1[0] - v2[0], 2) +
                        pow(v1[1] - v2[1], 2) +
                        pow(v1[2] - v2[2], 2));
        }

        float get_distanceVertex(edm4hep::VertexData v1, edm4hep::VertexData v2, int comp) {
            float result;
            if (comp == 0) result = v1.position.x - v2.position.x;
            else if (comp == 1) result = v1.position.y - v2.position.y;
            else if (comp == 2) result = v1.position.z - v2.position.z;

            else
                result = sqrt(pow(v1.position.x - v2.position.x, 2) +
                              pow(v1.position.y - v2.position.y, 2) +
                              pow(v1.position.z - v2.position.z, 2));
            return result;
        }


        float get_distanceErrorVertex(edm4hep::VertexData v1, edm4hep::VertexData v2, int comp) {
            std::array<float, 6> v1_covMatrix = v1.covMatrix;
            std::array<float, 6> v2_covMatrix = v2.covMatrix;

            //when error on x, y, z only
            if (comp == 0) return sqrt(v1_covMatrix[0] + v2_covMatrix[0]);
            else if (comp == 1) return sqrt(v1_covMatrix[2] + v2_covMatrix[2]);
            else if (comp == 2) return sqrt(v1_covMatrix[5] + v2_covMatrix[5]);

            float result;

            edm4hep::Vector3f v1_position = v1.position;
            edm4hep::Vector3f v2_position = v2.position;

            float distance = get_distanceVertex(v1, v2, comp);

            float x =
                    (v1_position[0] - v2_position[0]) * (v1_covMatrix[0] + v2_covMatrix[0]) +
                    (v1_position[1] - v2_position[1]) * (v1_covMatrix[1] + v2_covMatrix[1]) +
                    (v1_position[2] - v2_position[2]) * (v1_covMatrix[3] + v2_covMatrix[3]);

            float y =
                    (v1_position[0] - v2_position[0]) * (v1_covMatrix[1] + v2_covMatrix[1]) +
                    (v1_position[1] - v2_position[1]) * (v1_covMatrix[2] + v2_covMatrix[2]) +
                    (v1_position[2] - v2_position[2]) * (v1_covMatrix[4] + v2_covMatrix[4]);

            float z =
                    (v1_position[0] - v2_position[0]) * (v1_covMatrix[3] + v2_covMatrix[3]) +
                    (v1_position[1] - v2_position[1]) * (v1_covMatrix[4] + v2_covMatrix[4]) +
                    (v1_position[2] - v2_position[2]) * (v1_covMatrix[5] + v2_covMatrix[5]);


            //\sigma_d = (\vec{x}_1-\vec{x}_2)^t\{C_1+C_2|}(\vec{x}_1-\vec{x}_2)/d^2
            //Where d is the distance between the two vertices
            //= \Sqrt((\vec{x}_1-\vec{x}_2)^t(\vec{x}_1-\vec{x}_2))

            x = x * (v1_position[0] - v2_position[0]);
            y = y * (v1_position[1] - v2_position[1]);
            z = z * (v1_position[2] - v2_position[2]);

            result = x + y + z;
            if (result > 0)
                result = sqrt((result) / pow(distance, 2));
            else result = 0;


            /*
            TVectorD x1(3);
            x1(0) = v1_position[0];
            x1(1) = v1_position[1];
            x1(2) = v1_position[2];

            TVectorD x2(3);
            x2(0) = v2_position[0];
            x2(1) = v2_position[1];
            x2(2) = v2_position[2];

            TVectorD xdiff = x1-x2;

            TMatrixDSym C1(3);
            C1[0][0] = v1_covMatrix[0];

            C1[1][0] = v1_covMatrix[1];
            C1[1][1] = v1_covMatrix[2];

            C1[0][1] = v1_covMatrix[1];

            C1[2][0] = v1_covMatrix[3];
            C1[2][1] = v1_covMatrix[4];
            C1[2][2] = v1_covMatrix[5];

            C1[0][2] = v1_covMatrix[3];
            C1[1][2] = v1_covMatrix[4];


            TMatrixDSym C2(3);
            C2[0][0] = v2_covMatrix[0];

            C2[1][0] = v2_covMatrix[1];
            C2[1][1] = v2_covMatrix[2];

            C2[0][1] = v2_covMatrix[1];

            C2[2][0] = v2_covMatrix[3];
            C2[2][1] = v2_covMatrix[4];
            C2[2][2] = v2_covMatrix[5];

            C2[0][2] = v2_covMatrix[3];
            C2[1][2] = v2_covMatrix[4];

            TMatrixDSym Csum = C1+C2;
            Double_t d = TMath::Sqrt(xdiff(0)*xdiff(0) +xdiff(1)*xdiff(1) +xdiff(2)*xdiff(2));
            TVectorD xDir = (1./d)*xdiff;
            //
            Double_t sig_d = TMath::Sqrt(Csum.Similarity(xDir));

            std::cout << "my result " << result << "  FB  " << sig_d <<std::endl;
            */

            return result;
        }


        sel_PV::sel_PV(bool arg_closest): m_closest(arg_closest) {
        };

        VertexingUtils::FCCAnalysesVertex
        sel_PV::operator()(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> pv) {
            VertexingUtils::FCCAnalysesVertex result;

            float min = 999999;
            for (auto &p: pv) {
                edm4hep::VertexData v = p.vertex;
                float dist = sqrt(pow(v.position.x, 2) + pow(v.position.y, 2) + pow(v.position.z, 2));
                if (dist < min) {
                    min = dist;
                    result = p;
                }
            }
            return result;
        }


        ROOT::VecOps::RVec<float> get_flightDistanceVertex(ROOT::VecOps::RVec<FCCAnalysesComposite> in,
                                                           edm4hep::VertexData pv) {
            ROOT::VecOps::RVec<float> result;

            for (auto &sv: in) {
                edm4hep::VertexData theSV = sv.vertex;
                result.push_back(get_distanceVertex(pv, theSV, -1));
            }

            return result;
        }


        ROOT::VecOps::RVec<float> get_flightDistanceVertex(ROOT::VecOps::RVec<FCCAnalysesComposite> in,
                                                           VertexingUtils::FCCAnalysesVertex pv) {
            ROOT::VecOps::RVec<float> result;
            edm4hep::VertexData thePV = pv.vertex;

            return get_flightDistanceVertex(in, thePV);
        }


        ROOT::VecOps::RVec<int> getMC_daughter(int daughterindex,
                                               ROOT::VecOps::RVec<edm4hep::MCParticleData> in,
                                               ROOT::VecOps::RVec<int> ind) {
            ROOT::VecOps::RVec<int> result;
            for (size_t i = 0; i < in.size(); ++i) {
                if (daughterindex + 1 > in.at(i).daughters_end - in.at(i).daughters_begin) {
                    result.push_back(-999);
                } else {
                    result.push_back(ind.at(in.at(i).daughters_begin + daughterindex));
                }
            }
            return result;
        }

        ROOT::VecOps::RVec<int> getMC_parent(int parentindex,
                                             ROOT::VecOps::RVec<edm4hep::MCParticleData> in,
                                             ROOT::VecOps::RVec<int> ind) {
            ROOT::VecOps::RVec<int> result;
            for (size_t i = 0; i < in.size(); ++i) {
                if (parentindex + 1 > in.at(i).parents_end - in.at(i).parents_begin) {
                    result.push_back(-999);
                } else {
                    result.push_back(ind.at(in.at(i).parents_begin + parentindex));
                }
            }
            return result;
        }

        int getMC_parent(int parentindex,
                         edm4hep::MCParticleData in,
                         ROOT::VecOps::RVec<int> ind) {
            int result;
            if (parentindex + 1 > in.parents_end - in.parents_begin)
                result = -999;
            else
                result = ind.at(in.parents_begin + parentindex);
            return result;
        }


        ROOT::VecOps::RVec<FCCAnalysesComposite> add_truthmatched(ROOT::VecOps::RVec<FCCAnalysesComposite> comp,
                                                                  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc,
                                                                  //ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> rp2mc){
                                                                  ROOT::VecOps::RVec<int> rp2mc,
                                                                  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>
                                                                  recop,
                                                                  ROOT::VecOps::RVec<int> ind) {
            for (size_t i = 0; i < comp.size(); ++i) {
                //std::cout << "compo " << i << "  charge " << comp.at(i).charge<< std::endl;
                ROOT::VecOps::RVec<int> index = comp.at(i).index;
                ROOT::VecOps::RVec<int> mother;
                ROOT::VecOps::RVec<int> motherPDG;

                for (size_t j = 0; j < index.size(); ++j) {
                    //ROOT::VecOps::RVec<int> mcassoc = rp2mc.at(index.at(j));
                    int mcassoc = rp2mc.at(index.at(j));
                    //if (mcassoc.size()==1){
                    //mother.push_back(mcassoc.at(0));
                    int mother1 = getMC_parent(0, mc.at(mcassoc), ind);
                    int mother2 = getMC_parent(1, mc.at(mcassoc), ind);

                    mother.push_back(mother1);
                    motherPDG.push_back(mc.at(mother1).PDG);
                    //std::cout << "mother 1 "<<mother1<<"  mother2  " << mother2<< std::endl;
                    //std::cout << " mc assoc j " << j << "  rp index  "<< index.at(j) << " mc index " << mcassoc <<  "  PDG ID  " << mc.at(mcassoc).PDG << "  charge "<< mc.at(mcassoc).charge<< " mc p="<<sqrt(pow(mc.at(mcassoc).momentum.z,2)+pow(mc.at(mcassoc).momentum.y,2)+pow(mc.at(mcassoc).momentum.x,2)) << "  rp p " << get_p(recop.at(index.at(j)))<<std::endl;

                    //std::cout << " mc assoc j " << j << "   " << mcassoc.at(0) <<  "  PDG ID  " << mc.at(mcassoc.at(0)).PDG << "  charge "<< mc.at(mcassoc.at(0)).charge<< "  px="<< mc.at(mcassoc.at(0)).momentum.x << "  py="<< mc.at(mcassoc.at(0)).momentum.y << "  pz="<< mc.at(mcassoc.at(0)).momentum.z << "  p="<<sqrt(pow(mc.at(mcassoc.at(0)).momentum.z,2)+pow(mc.at(mcassoc.at(0)).momentum.y,2)+pow(mc.at(mcassoc.at(0)).momentum.x,2))<<std::endl;
                    //}
                    //else std::cout <<"================================================================================================================MORE THAN 1 ASSOC"<<std::endl;
                }

                if (mother.size() > 0) {
                    bool truthmatched = true;
                    int tmpmother = mother.at(0);
                    int tmpmotherpdg = motherPDG.at(0);
                    for (size_t k = 1; k < mother.size(); ++k) {
                        if (tmpmother != mother.at(k) || tmpmotherpdg != motherPDG.at(k))truthmatched = false;
                    }

                    if (truthmatched == true) { comp.at(i).mc_index = mother.at(0); }
                    //std::cout <<"==============================mthchehhfewefkwefkwfwf"<<std::endl;}
                    else comp.at(i).mc_index = -999;
                }
            }

            /*  for (size_t i = 0; i < rp2mc.size(); ++i) {
             std::cout <<"RP " << i << std::endl;
             for (size_t j = 0; j < rp2mc.at(i).size(); ++j) {
               std::cout << "  MC  " << rp2mc.at(i).at(j) << std::endl;


             }
             }*/
            return comp;
        }


        ROOT::VecOps::RVec<int> get_trueVertex(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC> vertex,
                                               ROOT::VecOps::RVec<edm4hep::MCParticleData> mc,
                                               ROOT::VecOps::RVec<int> ind,
                                               int mother,
                                               int grandmother) {
            //std::cout<<"NEW EVENT"<<std::endl;
            ROOT::VecOps::RVec<int> result;
            for (size_t i = 0; i < vertex.size(); ++i) {
                ROOT::VecOps::RVec<int> mother_ind = vertex.at(i).mother_ind;
                for (size_t j = 0; j < mother_ind.size(); ++j) {
                    if (fabs(mc.at(mother_ind.at(j)).PDG) == mother) {
                        for (size_t k = mc.at(mother_ind.at(j)).parents_begin; k < mc.at(mother_ind.at(j)).parents_end;
                             ++k) {
                            if (fabs(mc.at(ind.at(k)).PDG) == grandmother) {
                                result.push_back(i);
                                //std::cout <<"i ="<< i << "  j=" << j << "  k=" << k << "  mother PDG " << mc.at(mother_ind.at(j)).PDG
                                //	      << "  gmother PDG " << mc.at(ind.at(k)).PDG<< std::endl;
                            }
                        }
                    }
                }
                //result.push_back(i);
            }
            return result;
        }

        ROOT::VecOps::RVec<FCCAnalysesComposite2> add_truthmatched2(ROOT::VecOps::RVec<FCCAnalysesComposite2> comp,
                                                                    ROOT::VecOps::RVec<edm4hep::MCParticleData> mc,
                                                                    ROOT::VecOps::RVec<
                                                                        VertexingUtils::FCCAnalysesVertex> vertex,
                                                                    ROOT::VecOps::RVec<int> rp2mc,
                                                                    ROOT::VecOps::RVec<
                                                                        edm4hep::ReconstructedParticleData> recop,
                                                                    ROOT::VecOps::RVec<int> ind) {
            for (size_t i = 0; i < comp.size(); ++i) {
                ROOT::VecOps::RVec<int> index = vertex.at(comp.at(i).vertex).reco_ind;
                ROOT::VecOps::RVec<int> mother;
                ROOT::VecOps::RVec<int> motherPDG;

                for (size_t j = 0; j < index.size(); ++j) {
                    int mcassoc = rp2mc.at(index.at(j));
                    int mother1 = getMC_parent(0, mc.at(mcassoc), ind);
                    int mother2 = getMC_parent(1, mc.at(mcassoc), ind);

                    mother.push_back(mother1);
                    motherPDG.push_back(mc.at(mother1).PDG);
                }

                if (mother.size() > 0) {
                    bool truthmatched = true;
                    int tmpmother = mother.at(0);
                    int tmpmotherpdg = motherPDG.at(0);
                    for (size_t k = 1; k < mother.size(); ++k) {
                        if (tmpmother != mother.at(k) || tmpmotherpdg != motherPDG.at(k))truthmatched = false;
                    }

                    if (truthmatched == true) comp.at(i).mc_index = mother.at(0);
                    else comp.at(i).mc_index = -999;
                }
            }
            return comp;
        }

        ROOT::VecOps::RVec<int> get_compmc(ROOT::VecOps::RVec<FCCAnalysesComposite> in) {
            ROOT::VecOps::RVec<int> result;
            for (size_t i = 0; i < in.size(); ++i)result.push_back(in.at(i).mc_index);
            return result;
        }


        /*ROOT::VecOps::RVec<FCCAnalysesComposite> add_truthmatched(ROOT::VecOps::RVec<FCCAnalysesComposite> comp,
                                           ROOT::VecOps::RVec<edm4hep::MCParticleData> mc,
                                           ROOT::VecOps::RVec<int> mcind,
                                           ROOT::VecOps::RVec<int> recoind){


          ROOT::VecOps::RVec<FCCAnalysesComposite> result;

          for (size_t i = 0; i < comp.size(); ++i) {
            ROOT::VecOps::RVec<int> index = i.index;
            for (size_t j = 0; j < index.size(); ++j) {



            }
          }

            if (parentindex+1>in.at(i).parents_end-in.at(i).parents_begin) {
              result.push_back(-999);
            }
            else {
              result.push_back(ind.at(in.at(i).parents_begin+parentindex));
            }
          }
          return result;
        }

        }*/

        ROOT::VecOps::RVec<FCCAnalysesComposite> build_Bu2D0Pi(
            ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop,
            ROOT::VecOps::RVec<FCCAnalysesComposite> d0,
            ROOT::VecOps::RVec<int> pions) {
            ROOT::VecOps::RVec<FCCAnalysesComposite> result;
            for (size_t i = 0; i < d0.size(); ++i) {
                ROOT::VecOps::RVec<int> d0index = d0.at(i).index;
                float kaoncharge = 0;
                //std::cout << "index sise "<< d0index.size() << " ind 0 " << d0index.at(0)<< " ind 1 " << d0index.at(1)  <<std::endl;
                //std::cout << " recop.at(index.at(0)) "<<recop.at(d0index.at(0)).type<< " recop.at(index.at(1)) "<<recop.at(d0index.at(1)).type<< std::endl;

                if (recop.at(d0index.at(0)).type == 321)kaoncharge = recop.at(d0index.at(0)).charge;
                else if (recop.at(d0index.at(1)).type == 321)kaoncharge = recop.at(d0index.at(1)).charge;
                else std::cout << "huston there iis a problem no kaon found build_Bu2D0Pi" << std::endl;
                for (size_t j = 0; j < pions.size(); ++j) {
                    if (get_p(recop.at(pions.at(j))) < 1.)continue;
                    if (kaoncharge != recop.at(pions.at(j)).charge)continue;

                    //Mass cut
                    TLorentzVector tlvpion = ReconstructedParticle::get_tlv(recop.at(pions.at(j)));
                    TLorentzVector tlvd0 = d0.at(i).particle;
                    TLorentzVector tlvB = tlvpion + tlvd0;

                    FCCAnalysesComposite B;
                    ROOT::VecOps::RVec<int> index;
                    index.push_back(d0index.at(0));
                    index.push_back(d0index.at(1));
                    index.push_back(pions.at(j));
                    B.particle = tlvB;
                    B.index = index;
                    result.push_back(B);

                    //if (fabs(tlvD0.M()-1.86483)>m_mass)continue;
                }
            }
            return result;
        }


        filter_PV::filter_PV(bool arg_pv): m_pv(arg_pv) {
        };

        ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>
        filter_PV::operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in,
                              ROOT::VecOps::RVec<int> index) {
            ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;

            for (auto &p: in) {
                bool found = false;
                for (size_t i = 0; i < index.size(); ++i) {
                    if (p.tracks_begin == index.at(i)) {
                        found = true;
                        break;
                    }
                }
                if (found == false && m_pv == false)result.push_back(p);
                else if (found == true && m_pv == true)result.push_back(p);
            }
            return result;
        }


        int getFCCAnalysesComposite_N(ROOT::VecOps::RVec<FCCAnalysesComposite> in) {
            return in.size();
        }

        ROOT::VecOps::RVec<float> getFCCAnalysesComposite_mass(ROOT::VecOps::RVec<FCCAnalysesComposite> in) {
            ROOT::VecOps::RVec<float> result;
            for (auto &p: in) {
                result.push_back(p.particle.M());
            }
            return result;
        }

        int getFCCAnalysesComposite_N(ROOT::VecOps::RVec<FCCAnalysesComposite2> in) {
            return in.size();
        }

        ROOT::VecOps::RVec<float> getFCCAnalysesComposite_mass(ROOT::VecOps::RVec<FCCAnalysesComposite2> in) {
            ROOT::VecOps::RVec<float> result;
            for (auto &p: in) {
                result.push_back(p.particle.M());
            }
            return result;
        }

        ROOT::VecOps::RVec<float> getFCCAnalysesComposite_mass(ROOT::VecOps::RVec<FCCAnalysesComposite2> in,
                                                               ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex>
                                                               vertex) {
            ROOT::VecOps::RVec<float> result;
            for (auto &p: in) {
                ROOT::VecOps::RVec<TVector3> updated_track_momentum_at_vertex = vertex.at(p.vertex).
                        updated_track_momentum_at_vertex;
                TLorentzVector tlv;
                for (size_t i = 0; i < updated_track_momentum_at_vertex.size(); ++i) {
                    TLorentzVector tlvtmp;
                    tlvtmp.SetXYZM(updated_track_momentum_at_vertex.at(i).Px(),
                                   updated_track_momentum_at_vertex.at(i).Py(),
                                   updated_track_momentum_at_vertex.at(i).Pz(),
                                   p.particle.M());
                    tlv += tlvtmp;
                }
                result.push_back(tlv.M());
            }
            return result;
        }

        ROOT::VecOps::RVec<float> getFCCAnalysesComposite_mass(ROOT::VecOps::RVec<FCCAnalysesComposite2> in,
                                                               ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex>
                                                               vertex,
                                                               ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>
                                                               recop,
                                                               int index) {
            ROOT::VecOps::RVec<float> result;
            for (auto &p: in) {
                int recoind = vertex.at(p.vertex).reco_ind.at(index);
                result.push_back(recop.at(recoind).mass);
            }
            return result;
        }

        ROOT::VecOps::RVec<int> getFCCAnalysesComposite_type(ROOT::VecOps::RVec<FCCAnalysesComposite2> in,
                                                             ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex>
                                                             vertex,
                                                             ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>
                                                             recop,
                                                             int index) {
            ROOT::VecOps::RVec<int> result;
            for (auto &p: in) {
                int recoind = vertex.at(p.vertex).reco_ind.at(index);
                result.push_back(recop.at(recoind).type);
            }
            return result;
        }


        ROOT::VecOps::RVec<int> getFCCAnalysesComposite_vertex(ROOT::VecOps::RVec<FCCAnalysesComposite2> in) {
            ROOT::VecOps::RVec<int> result;
            for (auto &p: in) {
                result.push_back(p.vertex);
            }
            return result;
        }

        ROOT::VecOps::RVec<int> getFCCAnalysesComposite_mcvertex(ROOT::VecOps::RVec<FCCAnalysesComposite2> in,
                                                                 ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex>
                                                                 vertex) {
            ROOT::VecOps::RVec<int> result;
            for (auto &p: in) {
                result.push_back(vertex.at(p.vertex).mc_ind);
            }
            return result;
        }


        ROOT::VecOps::RVec<float> getFCCAnalysesComposite_B(ROOT::VecOps::RVec<FCCAnalysesComposite2> in,
                                                            ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex>
                                                            vertex,
                                                            ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>
                                                            recop) {
            ROOT::VecOps::RVec<float> result;
            for (auto &p: in) {
                float energy = 0;
                ROOT::VecOps::RVec<int> reco_ind = vertex.at(p.vertex).reco_ind;
                std::vector<int> reco_indstd;
                for (size_t i = 0; i < reco_ind.size(); ++i) reco_indstd.push_back(reco_ind.at(i));

                for (size_t i = 0; i < recop.size(); ++i) {
                    std::vector<int>::iterator it = std::find(reco_indstd.begin(), reco_indstd.end(), i);
                    if (it != reco_indstd.end()) continue;
                    energy += recop.at(i).energy;
                }
                result.push_back(energy);
            }
            return result;
        }


        ROOT::VecOps::RVec<float> getFCCAnalysesComposite_p(ROOT::VecOps::RVec<FCCAnalysesComposite2> in,
                                                            int i) {
            ROOT::VecOps::RVec<float> result;
            for (auto &p: in) {
                if (i == 0)result.push_back(p.particle.Px());
                else if (i == 1)result.push_back(p.particle.Py());
                else if (i == 2)result.push_back(p.particle.Pz());
                else
                    result.push_back(sqrt(pow(p.particle.Px(), 2) +
                                          pow(p.particle.Py(), 2) +
                                          pow(p.particle.Pz(), 2)));
            }
            return result;
        }


        ROOT::VecOps::RVec<float> getFCCAnalysesComposite_p(ROOT::VecOps::RVec<FCCAnalysesComposite2> in,
                                                            ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex>
                                                            vertex,
                                                            ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>
                                                            recop,
                                                            int index,
                                                            int type) {
            ROOT::VecOps::RVec<float> result;
            for (auto &p: in) {
                int recoind = vertex.at(p.vertex).reco_ind.at(index);

                if (type == 0)result.push_back(recop.at(recoind).momentum.x);
                else if (type == 1)result.push_back(recop.at(recoind).momentum.y);
                else if (type == 2)result.push_back(recop.at(recoind).momentum.z);
                else
                    result.push_back(sqrt(pow(recop.at(recoind).momentum.x, 2) +
                                          pow(recop.at(recoind).momentum.y, 2) +
                                          pow(recop.at(recoind).momentum.z, 2)));
            }
            return result;
        }


        ROOT::VecOps::RVec<int> getFCCAnalysesComposite_q(ROOT::VecOps::RVec<FCCAnalysesComposite2> in,
                                                          ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
                                                          ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop,
                                                          int index) {
            ROOT::VecOps::RVec<int> result;
            for (auto &p: in) {
                int recoind = vertex.at(p.vertex).reco_ind.at(index);
                result.push_back(recop.at(recoind).charge);
            }
            return result;
        }


        ROOT::VecOps::RVec<edm4hep::TrackState> get_truetrack(ROOT::VecOps::RVec<int> in,
                                                              ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertexMC>
                                                              vertex,
                                                              ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) {
            ROOT::VecOps::RVec<edm4hep::TrackState> result;
            float charge = 0;
            float norm = 1e-3; // to convert from mm to meters
            for (auto &p: in) {
                ROOT::VecOps::RVec<int> mc_ind = vertex.at(p).mc_ind;
                TLorentzVector tlv;
                for (size_t i = 0; i < mc_ind.size(); ++i) {
                    TLorentzVector tlvtmp;
                    tlvtmp.SetXYZM(mc.at(mc_ind.at(i)).momentum.x,
                                   mc.at(mc_ind.at(i)).momentum.y,
                                   mc.at(mc_ind.at(i)).momentum.z,
                                   mc.at(mc_ind.at(i)).mass);
                    tlv += tlvtmp;
                    charge += mc.at(mc_ind.at(i)).charge;
                }

                edm4hep::TrackState track;
                TVector3 vertexFB(vertex.at(p).vertex.X() * norm,
                                  vertex.at(p).vertex.Y() * norm,
                                  vertex.at(p).vertex.Z() * norm);

                TVector3 momentum(tlv.Px(), tlv.Py(), tlv.Pz());

                TVectorD track_param = VertexingUtils::XPtoPar(vertexFB, momentum, charge);


                track.D0 = track_param[0] * 1e3; // from meters to mm
                track.phi = track_param[1];
                track.omega = track_param[2] / (0.5 * 1e3); // C from Franco = rho/2, and convert from m-1 to mm-1

                // need to change here, because the TracSTate of edm4heo currently use
                // the wrong sign !
                track.omega = -track.omega;

                track.Z0 = track_param[3] * 1e3; // from meters to mm
                track.tanLambda = track_param[4];

                track.referencePoint.x = vertexFB[0];
                track.referencePoint.y = vertexFB[1];
                track.referencePoint.z = vertexFB[2];
                result.push_back(track);
            }
            return result;
        }


        ROOT::VecOps::RVec<edm4hep::TrackState> get_pseudotrack(
            ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
            ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop) {
            ROOT::VecOps::RVec<edm4hep::TrackState> result;
            float norm = 1e-3; // to convert from mm to meters
            for (auto &p: vertex) {
                if (p.vertex.primary > 0)continue;
                edm4hep::TrackState track;
                TVector3 vertexFB(p.vertex.position.x * norm,
                                  p.vertex.position.y * norm,
                                  p.vertex.position.z * norm);
                TLorentzVector pseudop;
                float pseudopq = 0.;
                for (auto &q: p.reco_ind) {
                    TLorentzVector pseudotmp;
                    pseudotmp.SetXYZM(recop.at(q).momentum.x,
                                      recop.at(q).momentum.y,
                                      recop.at(q).momentum.z,
                                      recop.at(q).mass);
                    pseudop += pseudotmp;
                    pseudopq += recop.at(q).charge;
                }

                //std::cout<<"px, py, pz= " << pseudop.Px() << "  "<< pseudop.Py()<<"  "<<pseudop.Pz()<<std::endl;
                //std::cout<<"x, y, z=    " << vertexFB.X() << "  "<< vertexFB.Y()<<"  "<<vertexFB.Z()<<std::endl;
                //std::cout <<"charge="<<pseudopq<<std::endl;

                TVector3 momentum(pseudop.Px(), pseudop.Py(), pseudop.Pz());

                TVectorD track_param = VertexingUtils::XPtoPar(vertexFB, momentum, pseudopq);


                track.D0 = track_param[0] * 1e3; // from meters to mm
                track.phi = track_param[1];
                track.omega = track_param[2] / (0.5 * 1e3); // C from Franco = rho/2, and convert from m-1 to mm-1


                // need to change here, because the TracSTate of edm4heo currently use
                // the wrong sign !
                track.omega = -track.omega;

                track.Z0 = track_param[3] * 1e3; // from meters to mm
                track.tanLambda = track_param[4];

                track.referencePoint.x = vertexFB[0];
                track.referencePoint.y = vertexFB[1];
                track.referencePoint.z = vertexFB[2];


                //    std::cout <<"mag field D0, Z0    =  "<<track.D0<<"  "<<track.Z0<<std::endl;
                // std::cout <<"NO mag field D0, Z0 =  "<<get_d0(vertexFB, momentum)*1e3<<"  "<<get_z0(vertexFB, momentum)*1e3<<std::endl;

                track.D0 = get_d0(vertexFB, momentum) * 1e3; // from meters to mm
                track.Z0 = get_z0(vertexFB, momentum) * 1e3; // from meters to mm


                result.push_back(track);
            }
            return result;
        }


        ROOT::VecOps::RVec<edm4hep::TrackState> getFCCAnalysesComposite_track(
            ROOT::VecOps::RVec<FCCAnalysesComposite2> in,
            ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex) {
            ROOT::VecOps::RVec<edm4hep::TrackState> result;
            float norm = 1e-3; // to convert from mm to meters
            for (auto &p: in) {
                edm4hep::TrackState track;
                TVector3 vertexFB(vertex.at(p.vertex).vertex.position.x * norm,
                                  vertex.at(p.vertex).vertex.position.y * norm,
                                  vertex.at(p.vertex).vertex.position.z * norm);
                TVector3 momentum(p.particle.Px(), p.particle.Py(), p.particle.Pz());

                TVectorD track_param = VertexingUtils::XPtoPar(vertexFB, momentum, p.charge);


                track.D0 = track_param[0] * 1e3; // from meters to mm
                track.phi = track_param[1];
                track.omega = track_param[2] / (0.5 * 1e3); // C from Franco = rho/2, and convert from m-1 to mm-1

                // need to change here, because the TracSTate of edm4heo currently use
                // the wrong sign !
                track.omega = -track.omega;

                track.Z0 = track_param[3] * 1e3; // from meters to mm
                track.tanLambda = track_param[4];

                track.referencePoint.x = vertexFB[0];
                track.referencePoint.y = vertexFB[1];
                track.referencePoint.z = vertexFB[2];

                track.D0 = get_d0(vertexFB, momentum) * 1e3; // from meters to mm
                track.Z0 = get_z0(vertexFB, momentum) * 1e3; // from meters to mm


                result.push_back(track);
            }
            return result;
        }

        ROOT::VecOps::RVec<float> get_trackd0(ROOT::VecOps::RVec<edm4hep::TrackState> in) {
            ROOT::VecOps::RVec<float> result;
            for (auto &p: in) result.push_back(p.D0);
            return result;
        }

        ROOT::VecOps::RVec<float> get_trackz0(ROOT::VecOps::RVec<edm4hep::TrackState> in) {
            ROOT::VecOps::RVec<float> result;
            for (auto &p: in) result.push_back(p.Z0);
            return result;
        }


        ROOT::VecOps::RVec<float> getFCCAnalysesComposite_d0(ROOT::VecOps::RVec<FCCAnalysesComposite2> in,
                                                             ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex>
                                                             vertex,
                                                             int index) {
            ROOT::VecOps::RVec<float> result;
            for (auto &p: in) {
                float d0 = vertex.at(p.vertex).updated_track_parameters.at(index)[0];
                result.push_back(d0);
            }
            return result;
        }

        ROOT::VecOps::RVec<float> getFCCAnalysesComposite_z0(ROOT::VecOps::RVec<FCCAnalysesComposite2> in,
                                                             ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex>
                                                             vertex,
                                                             int index) {
            ROOT::VecOps::RVec<float> result;
            for (auto &p: in) {
                float z0 = vertex.at(p.vertex).updated_track_parameters.at(index)[3];
                result.push_back(z0);
            }
            return result;
        }

        ROOT::VecOps::RVec<int> getFCCAnalysesComposite_charge(ROOT::VecOps::RVec<FCCAnalysesComposite2> in) {
            ROOT::VecOps::RVec<int> result;
            for (auto &p: in) {
                result.push_back(p.charge);
            }
            return result;
        }

        ROOT::VecOps::RVec<int> getFCCAnalysesComposite_truthMatch(ROOT::VecOps::RVec<FCCAnalysesComposite2> in) {
            ROOT::VecOps::RVec<int> result;
            for (auto &p: in) {
                result.push_back(p.mc_index);
            }
            return result;
        }

        ROOT::VecOps::RVec<TLorentzVector>
        getFCCAnalysesComposite_particle(ROOT::VecOps::RVec<FCCAnalysesComposite> in) {
            ROOT::VecOps::RVec<TLorentzVector> result;
            for (auto &p: in) {
                result.push_back(p.particle);
            }
            return result;
        }


        ROOT::VecOps::RVec<ROOT::VecOps::RVec<int> > getFCCAnalysesComposite_index(
            ROOT::VecOps::RVec<FCCAnalysesComposite> in) {
            ROOT::VecOps::RVec<ROOT::VecOps::RVec<int> > result;
            for (auto &p: in) {
                result.push_back(p.index);
            }
            return result;
        }

        ROOT::VecOps::RVec<edm4hep::VertexData> getFCCAnalysesComposite_vertex(
            ROOT::VecOps::RVec<FCCAnalysesComposite> in) {
            ROOT::VecOps::RVec<edm4hep::VertexData> result;
            for (auto &p: in) {
                result.push_back(p.vertex);
            }
            return result;
        }

        bool isPV(edm4hep::ReconstructedParticleData recop, ROOT::VecOps::RVec<int> pvindex) {
            for (size_t i = 0; i < pvindex.size(); ++i) {
                if (recop.tracks_begin == pvindex.at(i))return true;
            }
            return false;
        }


        build_D0::build_D0(float arg_mass, float arg_p, bool arg_filterPV): m_mass(arg_mass), m_p(arg_p),
                                                                            m_filterPV(arg_filterPV) {
        };

        ROOT::VecOps::RVec<FCCAnalysesComposite>
        build_D0::operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop,
                             ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
                             ROOT::VecOps::RVec<int> pions,
                             ROOT::VecOps::RVec<int> kaons,
                             ROOT::VecOps::RVec<int> pvindex) {
            ROOT::VecOps::RVec<FCCAnalysesComposite> result;

            for (size_t i = 0; i < pions.size(); ++i) {
                //pion p cut
                if (get_p(recop.at(pions.at(i))) < m_p)continue;
                if (m_filterPV && isPV(recop.at(pions.at(i)), pvindex)) continue;
                for (size_t j = 0; j < kaons.size(); ++j) {
                    //kaon p cut
                    if (get_p(recop.at(kaons.at(j))) < m_p)continue;
                    if (m_filterPV && isPV(recop.at(kaons.at(j)), pvindex)) continue;
                    //std::cout << "i " << i << "  j " << j << "  pion a i " << pions.at(i) << "  kaons at j "<< kaons.at(j) << " nrecop  " << recop.size()<< std::endl;
                    //charge cut
                    int charge = recop.at(pions.at(i)).charge + recop.at(kaons.at(j)).charge;
                    if (charge != 0)continue;

                    //Mass cut
                    TLorentzVector tlvpion = ReconstructedParticle::get_tlv(recop.at(pions.at(i)));
                    TLorentzVector tlvkaon = ReconstructedParticle::get_tlv(recop.at(kaons.at(j)));
                    TLorentzVector tlvD0 = tlvpion + tlvkaon;
                    if (fabs(tlvD0.M() - 1.86483) > m_mass)continue;

                    //vertex cut
                    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> seltracks;
                    seltracks.push_back(recop.at(pions.at(i)));
                    seltracks.push_back(recop.at(kaons.at(j)));
                    VertexingUtils::FCCAnalysesVertex TheVertex =
                            VertexFitterSimple::VertexFitter(0, seltracks, tracks);
                    float chi2 = TheVertex.vertex.chi2;

                    if (chi2 < 0.01)continue;
                    if (chi2 > 10.)continue;

                    FCCAnalysesComposite D0;
                    ROOT::VecOps::RVec<int> index;
                    index.push_back(pions.at(i));
                    index.push_back(kaons.at(j));
                    D0.charge = charge;
                    //std::cout << "in d0 builder pion Index " << pions.at(i) << " type " << recop.at(pions.at(i)).type << "  kaon index " <<kaons.at(j)<<" type " << recop.at(kaons.at(j)).type  <<std::endl;
                    D0.vertex = TheVertex.vertex;
                    D0.particle = tlvD0;
                    D0.index = index;
                    result.push_back(D0);
                }
            }

            return result;
        }


        sel_PID::sel_PID(int arg_PDG): m_PDG(arg_PDG) {
        } ;

        ROOT::VecOps::RVec<int>
        sel_PID::operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop) {
            ROOT::VecOps::RVec<int> result;
            for (size_t i = 0; i < recop.size(); ++i) {
                if (recop.at(i).type == m_PDG)
                    result.push_back(i);
            }
            return result;
        }


        ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>
        PID(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop,
            ROOT::VecOps::RVec<int> recind,
            ROOT::VecOps::RVec<int> mcind,
            ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) {
            for (size_t i = 0; i < recind.size(); ++i) {
                //id a pion
                if (fabs(mc.at(mcind.at(i)).PDG) == 211) {
                    recop.at(recind.at(i)).type = 211;
                    recop.at(recind.at(i)).mass = 0.13957039;
                    recop.at(recind.at(i)).energy = sqrt(pow(recop.at(recind.at(i)).momentum.x, 2) +
                                                         pow(recop.at(recind.at(i)).momentum.y, 2) +
                                                         pow(recop.at(recind.at(i)).momentum.z, 2) +
                                                         pow(recop.at(recind.at(i)).mass, 2));
                }
                //id a kaon
                else if (fabs(mc.at(mcind.at(i)).PDG) == 321) {
                    recop.at(recind.at(i)).type = 321;
                    recop.at(recind.at(i)).mass = 0.493677;
                    recop.at(recind.at(i)).energy = sqrt(pow(recop.at(recind.at(i)).momentum.x, 2) +
                                                         pow(recop.at(recind.at(i)).momentum.y, 2) +
                                                         pow(recop.at(recind.at(i)).momentum.z, 2) +
                                                         pow(recop.at(recind.at(i)).mass, 2));
                }
                //id a proton
                else if (fabs(mc.at(mcind.at(i)).PDG) == 2212) {
                    recop.at(recind.at(i)).type = 2212;
                    recop.at(recind.at(i)).mass = 0.938272081;
                    recop.at(recind.at(i)).energy = sqrt(pow(recop.at(recind.at(i)).momentum.x, 2) +
                                                         pow(recop.at(recind.at(i)).momentum.y, 2) +
                                                         pow(recop.at(recind.at(i)).momentum.z, 2) +
                                                         pow(recop.at(recind.at(i)).mass, 2));
                }
                //id an electron
                else if (fabs(mc.at(mcind.at(i)).PDG) == 11) {
                    recop.at(recind.at(i)).type = 11;
                    recop.at(recind.at(i)).mass = 0.0005109989461;
                    recop.at(recind.at(i)).energy = sqrt(pow(recop.at(recind.at(i)).momentum.x, 2) +
                                                         pow(recop.at(recind.at(i)).momentum.y, 2) +
                                                         pow(recop.at(recind.at(i)).momentum.z, 2) +
                                                         pow(recop.at(recind.at(i)).mass, 2));
                }
                //id an muon
                else if (fabs(mc.at(mcind.at(i)).PDG) == 13) {
                    recop.at(recind.at(i)).type = 13;
                    recop.at(recind.at(i)).mass = 0.1056583745;
                    recop.at(recind.at(i)).energy = sqrt(pow(recop.at(recind.at(i)).momentum.x, 2) +
                                                         pow(recop.at(recind.at(i)).momentum.y, 2) +
                                                         pow(recop.at(recind.at(i)).momentum.z, 2) +
                                                         pow(recop.at(recind.at(i)).mass, 2));
                }
            }
            return recop;
        }


        ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> get_RP_atVertex(
            ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop,
            ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex) {
            for (auto &p: vertex) {
                ROOT::VecOps::RVec<int> reco_ind = p.reco_ind;
                ROOT::VecOps::RVec<TVector3> updated_track_momentum_at_vertex = p.updated_track_momentum_at_vertex;

                for (size_t i = 0; i < reco_ind.size(); ++i) {
                    recop.at(reco_ind.at(i)).momentum.x = updated_track_momentum_at_vertex.at(i).Px();
                    recop.at(reco_ind.at(i)).momentum.y = updated_track_momentum_at_vertex.at(i).Py();
                    recop.at(reco_ind.at(i)).momentum.z = updated_track_momentum_at_vertex.at(i).Pz();
                    //recop.at(reco_ind.at(i)).mass = 0.13957039;
                    recop.at(reco_ind.at(i)).energy = sqrt(pow(updated_track_momentum_at_vertex.at(i).Px(), 2) +
                                                           pow(updated_track_momentum_at_vertex.at(i).Py(), 2) +
                                                           pow(updated_track_momentum_at_vertex.at(i).Pz(), 2) +
                                                           pow(recop.at(reco_ind.at(i)).mass, 2));
                }
            }
            return recop;
        }


        float build_invmass(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop,
                            ROOT::VecOps::RVec<int> index) {
            float result = 0;
            TLorentzVector tlv;
            for (size_t i = 0; i < index.size(); i++) {
                TLorentzVector tmp_tlv = ReconstructedParticle::get_tlv(recop[index.at(i)]);
                tlv += tmp_tlv;
            }
            return tlv.M();
        }

        TLorentzVector build_tlv(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop,
                                 ROOT::VecOps::RVec<int> index) {
            float result = 0;
            TLorentzVector tlv;
            for (size_t i = 0; i < index.size(); i++) {
                TLorentzVector tmp_tlv = ReconstructedParticle::get_tlv(recop[index.at(i)]);
                tlv += tmp_tlv;
            }
            return tlv;
        }


        ROOT::VecOps::RVec<FCCAnalysesComposite2> build_tau23pi(
            ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
            ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop) {
            ROOT::VecOps::RVec<FCCAnalysesComposite2> result;

            int counter = 0;
            for (auto &p: vertex) {
                //not consider PV
                if (p.vertex.primary == 1) {
                    counter += 1;
                    continue;
                }
                //exactly 3 tracks
                if (p.ntracks != 3) {
                    counter += 1;
                    continue;
                }

                //3 tracks id as pions, |charge| is 1
                bool is3pi = true;
                int charge = 0;
                for (auto &r: p.reco_ind) {
                    if (recop.at(r).type != 211)is3pi = false;
                    charge += recop.at(r).charge;
                }
                if (is3pi == false) {
                    counter += 1;
                    continue;
                }
                if (fabs(charge) != 1) {
                    counter += 1;
                    continue;
                }

                FCCAnalysesComposite2 comp;
                comp.vertex = counter;
                comp.particle = build_tlv(recop, p.reco_ind);
                comp.charge = charge;

                result.push_back(comp);
                counter += 1;
            }
            return result;
        }


        ROOT::VecOps::RVec<FCCAnalysesComposite2> build_B2Kstee(
            ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
            ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop) {
            ROOT::VecOps::RVec<FCCAnalysesComposite2> result;

            int counter = 0;
            for (auto &p: vertex) {
                //not consider PV
                if (p.vertex.primary == 1) {
                    counter += 1;
                    continue;
                }
                //exactly 4 tracks
                if (p.ntracks != 4) {
                    counter += 1;
                    continue;
                }

                //2 tracks id as e+ e-
                int charge_ee = 0;
                int nobj_ee = 0;
                for (auto &r: p.reco_ind) {
                    if (recop.at(r).type == 11) {
                        nobj_ee += 1;
                        charge_ee += recop.at(r).charge;
                    }
                }
                //1 tracks id as kaon
                int charge_k = 0;
                int nobj_k = 0;
                for (auto &r: p.reco_ind) {
                    if (recop.at(r).type == 321) {
                        nobj_k += 1;
                        charge_k += recop.at(r).charge;
                    }
                }

                //1 tracks id as pion
                int charge_pi = 0;
                int nobj_pi = 0;
                for (auto &r: p.reco_ind) {
                    if (recop.at(r).type == 211) {
                        nobj_pi += 1;
                        charge_pi += recop.at(r).charge;
                    }
                }


                if (nobj_ee != 2) {
                    counter += 1;
                    continue;
                }
                if (nobj_pi != 1) {
                    counter += 1;
                    continue;
                }
                if (nobj_k != 1) {
                    counter += 1;
                    continue;
                }

                if (charge_ee != 0) {
                    counter += 1;
                    continue;
                }
                if (charge_pi + charge_k != 0) {
                    counter += 1;
                    continue;
                }
                FCCAnalysesComposite2 comp;
                comp.vertex = counter;
                comp.particle = build_tlv(recop, p.reco_ind);
                comp.charge = charge_ee + charge_pi + charge_k;

                result.push_back(comp);
                counter += 1;
            }
            return result;
        }

        ROOT::VecOps::RVec<FCCAnalysesComposite2> build_B2Kstmumu(
            ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
            ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop) {
            ROOT::VecOps::RVec<FCCAnalysesComposite2> result;

            int counter = 0;
            for (auto &p: vertex) {
                //not consider PV
                if (p.vertex.primary == 1) {
                    counter += 1;
                    continue;
                }
                //exactly 4 tracks
                if (p.ntracks != 4) {
                    counter += 1;
                    continue;
                }

                //2 tracks id as mu+ mu-
                int charge_mumu = 0;
                int nobj_mumu = 0;
                for (auto &r: p.reco_ind) {
                    if (recop.at(r).type == 13) {
                        nobj_mumu += 1;
                        charge_mumu += recop.at(r).charge;
                    }
                }

                //1 tracks id as kaon
                int charge_k = 0;
                int nobj_k = 0;
                for (auto &r: p.reco_ind) {
                    if (recop.at(r).type == 321) {
                        nobj_k += 1;
                        charge_k += recop.at(r).charge;
                    }
                }

                //1 tracks id as pion
                int charge_pi = 0;
                int nobj_pi = 0;
                for (auto &r: p.reco_ind) {
                    if (recop.at(r).type == 211) {
                        nobj_pi += 1;
                        charge_pi += recop.at(r).charge;
                    }
                }


                if (nobj_mumu != 2) {
                    counter += 1;
                    continue;
                }
                if (nobj_pi != 1) {
                    counter += 1;
                    continue;
                }
                if (nobj_k != 1) {
                    counter += 1;
                    continue;
                }

                if (charge_mumu != 0) {
                    counter += 1;
                    continue;
                }
                if (charge_pi + charge_k != 0) {
                    counter += 1;
                    continue;
                }

                FCCAnalysesComposite2 comp;
                comp.vertex = counter;
                comp.particle = build_tlv(recop, p.reco_ind);
                comp.charge = charge_mumu + charge_pi + charge_k;

                result.push_back(comp);
                counter += 1;
            }
            return result;
        }

        ROOT::VecOps::RVec<FCCAnalysesComposite2> build_Bd2KstNuNu(
            ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
            ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop) {
            ROOT::VecOps::RVec<FCCAnalysesComposite2> result;
            int counter = 0;
            for (auto &p: vertex) {
                //not consider PV
                if (p.vertex.primary == 1) {
                    counter += 1;
                    continue;
                }
                //exactly 2 tracks
                if (p.ntracks != 2) {
                    counter += 1;
                    continue;
                }

                //1 tracks id as kaon
                int charge_k = 0;
                int nobj_k = 0;
                for (auto &r: p.reco_ind) {
                    if (recop.at(r).type == 321) {
                        nobj_k += 1;
                        charge_k += recop.at(r).charge;
                    }
                }

                //1 tracks id as pion
                int charge_pi = 0;
                int nobj_pi = 0;
                for (auto &r: p.reco_ind) {
                    if (recop.at(r).type == 211) {
                        nobj_pi += 1;
                        charge_pi += recop.at(r).charge;
                    }
                }
                if (nobj_pi != 1) {
                    counter += 1;
                    continue;
                }
                if (nobj_k != 1) {
                    counter += 1;
                    continue;
                }
                if (charge_pi + charge_k != 0) {
                    counter += 1;
                    continue;
                }

                FCCAnalysesComposite2 comp;
                comp.vertex = counter;
                comp.particle = build_tlv(recop, p.reco_ind);
                comp.charge = charge_pi + charge_k;

                result.push_back(comp);
                counter += 1;
            }
            return result;
        }


        ROOT::VecOps::RVec<FCCAnalysesComposite2> build_Bs2PhiNuNu(
            ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
            ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop) {
            ROOT::VecOps::RVec<FCCAnalysesComposite2> result;
            //loop over the reconstructed vertex collection
            for (size_t i = 0; i < vertex.size(); i++) {
                //not consider PV, exactly 2 tracks
                if (vertex.at(i).vertex.primary == 1)continue;
                if (vertex.at(i).ntracks != 2) continue;

                //2 tracks id as kaons
                int charge_phi = 0;
                int nobj_phi = 0;
                for (auto &r: vertex.at(i).reco_ind) {
                    if (recop.at(r).type == 321) {
                        nobj_phi += 1;
                        charge_phi += recop.at(r).charge;
                    }
                }
                //select candidates with exactly 2 kaons and charge 0
                if (nobj_phi != 2) continue;
                if (charge_phi != 0) continue;

                //build a composite vertex
                FCCAnalysesComposite2 comp;
                comp.vertex = i;
                comp.particle = build_tlv(recop, vertex.at(i).reco_ind);
                comp.charge = charge_phi;

                //add the composite vertex to the collection
                result.push_back(comp);
            }
            return result;
        }


        ROOT::VecOps::RVec<FCCAnalysesComposite2> build_Bd2MuMu(
            ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
            ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop) {
            ROOT::VecOps::RVec<FCCAnalysesComposite2> result;
            //loop over the reconstructed vertex collection
            for (size_t i = 0; i < vertex.size(); i++) {
                //not consider PV, exactly 2 tracks
                if (vertex.at(i).vertex.primary == 1)continue;
                if (vertex.at(i).ntracks != 2) continue;

                //2 tracks id as muons
                int charge_Bd = 0;
                int nobj_Bd = 0;
                for (auto &r: vertex.at(i).reco_ind) {
                    if (recop.at(r).type == 13) {
                        nobj_Bd += 1;
                        charge_Bd += recop.at(r).charge;
                    }
                }
                //select candidates with exactly 2 muons and charge 0
                if (nobj_Bd != 2) continue;
                if (charge_Bd != 0) continue;

                //build a composite vertex
                FCCAnalysesComposite2 comp;
                comp.vertex = i;
                comp.particle = build_tlv(recop, vertex.at(i).reco_ind);
                comp.charge = charge_Bd;

                //add the composite vertex to the collection
                result.push_back(comp);
            }
            return result;
        }

        build_tau23pi::build_tau23pi(float arg_masslow, float arg_masshigh, float arg_p, float arg_angle,
                                     bool arg_rho): m_masslow(arg_masslow), m_masshigh(arg_masshigh), m_p(arg_p),
                                                    m_angle(arg_angle), m_rho(arg_rho) {
        };

        ROOT::VecOps::RVec<FCCAnalysesComposite2>
        build_tau23pi::operator()(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
                                  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop) {
            ROOT::VecOps::RVec<FCCAnalysesComposite2> result;
            //std::cout <<"n reco V " << vertex.size()<<std::endl;
            int counter = 0;
            for (auto &p: vertex) {
                //exactly 3 tracks
                if (p.ntracks != 3) {
                    counter += 1;
                    continue;
                }

                //3 tracks id as pions, |charge| is 1, p cut
                bool is3pi = true;
                bool pcut = true;
                bool hasrho = true;
                int charge = 0;
                float angle = -9999999.;
                for (auto &r: p.reco_ind) {
                    if (recop.at(r).type != 211)is3pi = false;
                    if (get_p(recop.at(r)) < m_p) pcut = false;
                    charge += recop.at(r).charge;
                    TVector3 p1(recop.at(r).momentum.x, recop.at(r).momentum.y, recop.at(r).momentum.z);

                    for (auto &r2: p.reco_ind) {
                        if (r == r2)continue;
                        TVector3 p2(recop.at(r2).momentum.x, recop.at(r2).momentum.y, recop.at(r2).momentum.z);
                        float delta = fabs(p1.Angle(p2));
                        if (delta > angle)angle = delta;
                        if (recop.at(r).charge != recop.at(r2).charge) {
                            ROOT::VecOps::RVec<int> tmpvec_rho;
                            tmpvec_rho.push_back(r);
                            tmpvec_rho.push_back(r2);
                            float mass_rho = build_invmass(recop, tmpvec_rho);
                            if (mass_rho > 0.9)hasrho = false;
                            if (mass_rho < 0.6)hasrho = false;
                        }
                    }
                }
                if (is3pi == false) {
                    counter += 1;
                    continue;
                }
                if (pcut == false) {
                    counter += 1;
                    continue;
                }
                if (fabs(charge) != 1) {
                    counter += 1;
                    continue;
                }
                if (m_rho && hasrho == false) {
                    counter += 1;
                    continue;
                }
                if (angle > m_angle)continue;

                float mass = build_invmass(recop, p.reco_ind);
                if (mass > m_masshigh || mass < m_masslow)continue;


                //std::cout << "SELECTED----------------"<<std::endl;
                //std::cout << "charge " << charge << std::endl;
                //std::cout << "mass   " << mass << std::endl;
                //std::cout << "chi2   " << p.vertex.chi2 << std::endl;
                //std::cout << "ntrk   " << p.ntracks << std::endl;
                //std::cout << "recoin   " << p.reco_ind.size() << std::endl;

                FCCAnalysesComposite2 comp;
                comp.vertex = counter;
                comp.particle = build_tlv(recop, p.reco_ind);
                comp.charge = charge;

                result.push_back(comp);
                counter += 1;
            }
            return result;
        }


        ROOT::VecOps::RVec<float> get_Vertex_thrusthemis_angle(
            ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
            ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop,
            ROOT::VecOps::RVec<float> thrust) {
            ROOT::VecOps::RVec<float> result;

            for (auto &p: vertex) {
                ROOT::VecOps::RVec<int> reco_ind = p.reco_ind;
                TLorentzVector tlv;
                for (auto &i: reco_ind) {
                    TLorentzVector tlvtmp = ReconstructedParticle::get_tlv(recop.at(i));
                    tlv += tlvtmp;
                }

                float angle = Algorithms::getAxisCosTheta(thrust, tlv.Px(), tlv.Py(), tlv.Pz());
                result.push_back(angle);
            }
            return result;
        }

        ROOT::VecOps::RVec<float> get_DVertex_thrusthemis_angle(
            ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
            ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop,
            ROOT::VecOps::RVec<float> thrust) {
            ROOT::VecOps::RVec<float> result;

            for (auto &p: vertex) {
                if (p.vertex.primary > 0)continue;
                ROOT::VecOps::RVec<int> reco_ind = p.reco_ind;
                TLorentzVector tlv;
                for (auto &i: reco_ind) {
                    TLorentzVector tlvtmp = ReconstructedParticle::get_tlv(recop.at(i));
                    tlv += tlvtmp;
                }

                float angle = Algorithms::getAxisCosTheta(thrust, tlv.Px(), tlv.Py(), tlv.Pz());
                result.push_back(angle);
            }
            return result;
        }


        ROOT::VecOps::RVec<int> get_Vertex_thrusthemis(ROOT::VecOps::RVec<float> angle,
                                                       int index) {
            ROOT::VecOps::RVec<int> result;
            for (auto &p: angle) {
                //positive angle == min energy
                if (index == 1) {
                    if (p > 0.) result.push_back(1);
                    else result.push_back(0);
                } else if (index == 0) {
                    if (p > 0.) result.push_back(0);
                    else result.push_back(1);
                } else std::cout << "unidentified index in get_Vertex_thrusthemis" << std::endl;
            }
            return result;
        }

        ROOT::VecOps::RVec<int> get_Vertex_thrusthemis_emin(ROOT::VecOps::RVec<float> angle,
                                                            float eneg,
                                                            float epos) {
            ROOT::VecOps::RVec<int> result;
            for (auto &p: angle) {
                if (eneg < epos && p < 0.) result.push_back(1);
                else if (eneg > epos && p > 0.) result.push_back(1);
                else result.push_back(0);
            }
            return result;
        }

        ROOT::VecOps::RVec<edm4hep::MCParticleData> build_truerho(ROOT::VecOps::RVec<int> vertexind,
                                                                  ROOT::VecOps::RVec<
                                                                      VertexingUtils::FCCAnalysesVertexMC> vertex,
                                                                  ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) {
            ROOT::VecOps::RVec<edm4hep::MCParticleData> result;
            int index = vertexind.at(0);
            VertexingUtils::FCCAnalysesVertexMC mcver = vertex.at(index);
            ROOT::VecOps::RVec<int> mc_ind = mcver.mc_ind;

            //std::cout <<"new MC cand"<<std::endl;
            //for (size_t i=0; i<mc_ind.size();++i)std::cout <<"MC i="<<i<< " px="<<mc.at(mc_ind.at(i)).momentum.x<< " py="<<mc.at(mc_ind.at(i)).momentum.y<< " pz="<<mc.at(mc_ind.at(i)).momentum.z<<std::endl;

            if (mc.at(mc_ind.at(0)).charge != mc.at(mc_ind.at(1)).charge) {
                edm4hep::MCParticleData cand;
                TLorentzVector rho1;
                rho1.SetXYZM(mc.at(mc_ind.at(0)).momentum.x, mc.at(mc_ind.at(0)).momentum.y,
                             mc.at(mc_ind.at(0)).momentum.z, mc.at(mc_ind.at(0)).mass);
                TLorentzVector rho2;
                rho2.SetXYZM(mc.at(mc_ind.at(1)).momentum.x, mc.at(mc_ind.at(1)).momentum.y,
                             mc.at(mc_ind.at(1)).momentum.z, mc.at(mc_ind.at(1)).mass);
                TLorentzVector rho = rho1 + rho2;
                cand.momentum.x = rho.Px();
                cand.momentum.y = rho.Py();
                cand.momentum.z = rho.Pz();
                cand.mass = rho.M();
                result.push_back(cand);
            }
            if (mc.at(mc_ind.at(0)).charge != mc.at(mc_ind.at(2)).charge) {
                edm4hep::MCParticleData cand;
                TLorentzVector rho1;
                rho1.SetXYZM(mc.at(mc_ind.at(0)).momentum.x, mc.at(mc_ind.at(0)).momentum.y,
                             mc.at(mc_ind.at(0)).momentum.z, mc.at(mc_ind.at(0)).mass);
                TLorentzVector rho2;
                rho2.SetXYZM(mc.at(mc_ind.at(2)).momentum.x, mc.at(mc_ind.at(2)).momentum.y,
                             mc.at(mc_ind.at(2)).momentum.z, mc.at(mc_ind.at(2)).mass);
                TLorentzVector rho = rho1 + rho2;
                cand.momentum.x = rho.Px();
                cand.momentum.y = rho.Py();
                cand.momentum.z = rho.Pz();
                cand.mass = rho.M();
                result.push_back(cand);
            }
            if (mc.at(mc_ind.at(1)).charge != mc.at(mc_ind.at(2)).charge) {
                edm4hep::MCParticleData cand;
                TLorentzVector rho1;
                rho1.SetXYZM(mc.at(mc_ind.at(1)).momentum.x, mc.at(mc_ind.at(1)).momentum.y,
                             mc.at(mc_ind.at(1)).momentum.z, mc.at(mc_ind.at(1)).mass);
                TLorentzVector rho2;
                rho2.SetXYZM(mc.at(mc_ind.at(2)).momentum.x, mc.at(mc_ind.at(2)).momentum.y,
                             mc.at(mc_ind.at(2)).momentum.z, mc.at(mc_ind.at(2)).mass);
                TLorentzVector rho = rho1 + rho2;
                cand.momentum.x = rho.Px();
                cand.momentum.y = rho.Py();
                cand.momentum.z = rho.Pz();
                cand.mass = rho.M();
                result.push_back(cand);
            }

            return result;
        }

        ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> > build_rho(
            ROOT::VecOps::RVec<FCCAnalysesComposite2> in,
            ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
            ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop) {
            ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> > result;

            for (auto &p: in) {
                ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> tmp;
                ROOT::VecOps::RVec<int> reco_ind = vertex.at(p.vertex).reco_ind;
                if (reco_ind.size() != 3)continue;
                //std::cout <<"new reco cand"<<std::endl;
                //for (size_t i=0; i<reco_ind.size();++i)std::cout <<"reco i="<<i<< " px="<<recop.at(reco_ind.at(i)).momentum.x<< " py="<<recop.at(reco_ind.at(i)).momentum.y<< " pz="<<recop.at(reco_ind.at(i)).momentum.z<<std::endl;

                if (recop[reco_ind.at(0)].charge != recop[reco_ind.at(1)].charge) {
                    edm4hep::ReconstructedParticleData cand;
                    TLorentzVector rho1;
                    rho1.SetXYZM(recop.at(reco_ind.at(0)).momentum.x, recop.at(reco_ind.at(0)).momentum.y,
                                 recop.at(reco_ind.at(0)).momentum.z, recop.at(reco_ind.at(0)).mass);
                    TLorentzVector rho2;
                    rho2.SetXYZM(recop.at(reco_ind.at(1)).momentum.x, recop.at(reco_ind.at(1)).momentum.y,
                                 recop.at(reco_ind.at(1)).momentum.z, recop.at(reco_ind.at(1)).mass);
                    TLorentzVector rho = rho1 + rho2;
                    cand.momentum.x = rho.Px();
                    cand.momentum.y = rho.Py();
                    cand.momentum.z = rho.Pz();
                    cand.mass = rho.M();
                    cand.energy = rho.E();
                    tmp.push_back(cand);
                }
                if (recop[reco_ind.at(0)].charge != recop[reco_ind.at(2)].charge) {
                    edm4hep::ReconstructedParticleData cand;
                    TLorentzVector rho1;
                    rho1.SetXYZM(recop.at(reco_ind.at(0)).momentum.x, recop.at(reco_ind.at(0)).momentum.y,
                                 recop.at(reco_ind.at(0)).momentum.z, recop.at(reco_ind.at(0)).mass);
                    TLorentzVector rho2;
                    rho2.SetXYZM(recop.at(reco_ind.at(2)).momentum.x, recop.at(reco_ind.at(2)).momentum.y,
                                 recop.at(reco_ind.at(2)).momentum.z, recop.at(reco_ind.at(2)).mass);
                    TLorentzVector rho = rho1 + rho2;
                    cand.momentum.x = rho.Px();
                    cand.momentum.y = rho.Py();
                    cand.momentum.z = rho.Pz();
                    cand.mass = rho.M();
                    cand.energy = rho.E();
                    tmp.push_back(cand);
                }
                if (recop[reco_ind.at(1)].charge != recop[reco_ind.at(2)].charge) {
                    edm4hep::ReconstructedParticleData cand;
                    TLorentzVector rho1;
                    rho1.SetXYZM(recop.at(reco_ind.at(1)).momentum.x, recop.at(reco_ind.at(1)).momentum.y,
                                 recop.at(reco_ind.at(1)).momentum.z, recop.at(reco_ind.at(1)).mass);
                    TLorentzVector rho2;
                    rho2.SetXYZM(recop.at(reco_ind.at(2)).momentum.x, recop.at(reco_ind.at(2)).momentum.y,
                                 recop.at(reco_ind.at(2)).momentum.z, recop.at(reco_ind.at(2)).mass);
                    TLorentzVector rho = rho1 + rho2;
                    cand.momentum.x = rho.Px();
                    cand.momentum.y = rho.Py();
                    cand.momentum.z = rho.Pz();
                    cand.mass = rho.M();
                    cand.energy = rho.E();
                    tmp.push_back(cand);
                }

                result.push_back(tmp);
            }
            return result;
        }


        ROOT::VecOps::RVec<float> get_mass(
            ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> > in,
            int index) {
            ROOT::VecOps::RVec<float> result;
            for (auto &p: in) {
                if (index >= p.size())continue;
                result.push_back(p.at(index).mass);
            }
            return result;
        }

        ROOT::VecOps::RVec<float> get_px(ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> > in,
                                         int index) {
            ROOT::VecOps::RVec<float> result;

            for (auto &p: in)
                result.push_back(p.at(index).momentum.x);
            return result;
        }

        ROOT::VecOps::RVec<float> get_py(ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> > in,
                                         int index) {
            ROOT::VecOps::RVec<float> result;

            for (auto &p: in)
                result.push_back(p.at(index).momentum.y);
            return result;
        }


        ROOT::VecOps::RVec<float> get_pz(ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> > in,
                                         int index) {
            ROOT::VecOps::RVec<float> result;

            for (auto &p: in)
                result.push_back(p.at(index).momentum.z);
            return result;
        }


        float get_p(edm4hep::ReconstructedParticleData in) {
            TLorentzVector tlv;
            tlv.SetXYZM(in.momentum.x, in.momentum.y, in.momentum.z, in.mass);
            return tlv.P();
        }


        ROOT::VecOps::RVec<float> getFCCAnalysesComposite_anglethrust(ROOT::VecOps::RVec<FCCAnalysesComposite2> in,
                                                                      ROOT::VecOps::RVec<float> thrust) {
            ROOT::VecOps::RVec<float> result;
            TVector3 thrustvec(thrust.at(1), thrust.at(3), thrust.at(5));
            for (auto &p: in) {
                TVector3 cand(p.particle.X(), p.particle.Y(), p.particle.Z());
                result.push_back(cand.Angle(thrustvec));
            }
            return result;
        }


        int has_anglethrust_emin(ROOT::VecOps::RVec<float> angle) {
            for (auto &p: angle)
                if (cos(p) > 0.)return 1;
            return -1;
        }
    } //end NS myUtils
} //end NS FCCAnalyses
