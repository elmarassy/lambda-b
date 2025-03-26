#Mandatory: List of processes
processList = {
    # 'p8_ee_Zbb_ecm91_EvtGen_Bd2MuMu':{},
    # 'p8_ee_Zbb_ecm91_EvtGen_Lb2LMuMu':{'fraction':0.05}
    'p8_ee_Zbb_ecm91_EvtGen_Lb2LMuMuphsp':{'fraction':1/23}
    # 'p8_ee_Zss_ecm91':{'fraction':1/4388}
}

#Mandatory: Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics
prodTag     = "FCCee/winter2023/IDEA/"#"/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"

#Optional: output directory, default is local running directory
outputDir   = "output"
# outputDirEos   = "/eos/experiment/fcc/ee/analyses/case-studies/flavour/Bd2MuMu/flatNtuples/spring2021/analysis_stage1/"
#Optional
nCPUS       = 64
runBatch    = False
#batchQueue = "longlunch"
#compGroup = "group_u_FCC.local_gen"

#Mandatory: RDFanalysis class where the use defines the operations on the TTree
class RDFanalysis():

    #__________________________________________________________
    #Mandatory: analysers funtion to define the analysers to process, please make sure you return the last dataframe, in this example it is df2
    def analysers(df):
        df2 = (df

        .Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
        .Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")
        .Alias("Particle0", "Particle#0.index")
        .Alias("Particle1", "Particle#1.index")
        .Alias("muon0", "Muon#0.index")





        #############################################
        ##MC record to study the Z->bb events types##
        #############################################
        .Define("MC_PDG", "FCCAnalyses::MCParticle::get_pdg(Particle)")
        .Define("MC_n",   "int(MC_PDG.size())")

        #.Define("MC_M1",  "FCCAnalyses::myUtils::get_MCMother1(Particle,Particle0)")
        #.Define("MC_M2",  "FCCAnalyses::myUtils::get_MCMother2(Particle,Particle0)")
        #.Define("MC_D1",  "FCCAnalyses::myUtils::get_MCDaughter1(Particle,Particle1)")
        #.Define("MC_D2",  "FCCAnalyses::myUtils::get_MCDaughter2(Particle,Particle1)")
        .Define("MC_M1",  "FCCAnalyses::myUtils::getMC_parent(0,Particle,Particle0)")
        .Define("MC_M2",  "FCCAnalyses::myUtils::getMC_parent(1,Particle,Particle0)")
        .Define("MC_D1",  "FCCAnalyses::myUtils::getMC_daughter(0,Particle,Particle1)")
        .Define("MC_D2",  "FCCAnalyses::myUtils::getMC_daughter(1,Particle,Particle1)")
        .Define("MC_D3",  "FCCAnalyses::myUtils::getMC_daughter(2,Particle,Particle1)")
        .Define("MC_D4",  "FCCAnalyses::myUtils::getMC_daughter(3,Particle,Particle1)")
        .Define("MC_D5",  "FCCAnalyses::myUtils::getMC_daughter(4,Particle,Particle1)")
        .Define("MC_D6",  "FCCAnalyses::myUtils::getMC_daughter(5,Particle,Particle1)")
        .Define("MC_orivtx_x",   "FCCAnalyses::MCParticle::get_vertex_x(Particle)")
        .Define("MC_orivtx_y",   "FCCAnalyses::MCParticle::get_vertex_y(Particle)")
        .Define("MC_orivtx_z",   "FCCAnalyses::MCParticle::get_vertex_z(Particle)")
        .Define("MC_endvtx_x",   "FCCAnalyses::MCParticle::get_endPoint_x(Particle)")
        .Define("MC_endvtx_y",   "FCCAnalyses::MCParticle::get_endPoint_y(Particle)")
        .Define("MC_endvtx_z",   "FCCAnalyses::MCParticle::get_endPoint_z(Particle)")
        .Define("MC_p",   "FCCAnalyses::MCParticle::get_p(Particle)")
        .Define("MC_pt",  "FCCAnalyses::MCParticle::get_pt(Particle)")
        .Define("MC_px",  "FCCAnalyses::MCParticle::get_pt(Particle)")
        .Define("MC_py",  "FCCAnalyses::MCParticle::get_pt(Particle)")
        .Define("MC_pz",  "FCCAnalyses::MCParticle::get_pt(Particle)")
        .Define("MC_e",   "FCCAnalyses::MCParticle::get_e(Particle)")
        .Define("MC_m",   "FCCAnalyses::MCParticle::get_mass(Particle)")
        .Define("MC_q",   "FCCAnalyses::MCParticle::get_charge(Particle)")
        .Define("MC_eta", "FCCAnalyses::MCParticle::get_eta(Particle)")
        .Define("MC_phi", "FCCAnalyses::MCParticle::get_phi(Particle)"))

        #############################################
        ##               Build MC Vertex           ##
        #############################################
        df2 = (df2.Define("MCVertexObject", "FCCAnalyses::myUtils::get_MCVertexObject(Particle, Particle0)")
               .Define("MC_Vertex_x",    "FCCAnalyses::myUtils::get_MCVertex_x(MCVertexObject)")
               .Define("MC_Vertex_y",    "FCCAnalyses::myUtils::get_MCVertex_y(MCVertexObject)")
               .Define("MC_Vertex_z",    "FCCAnalyses::myUtils::get_MCVertex_z(MCVertexObject)")
               .Define("MC_Vertex_ind",  "FCCAnalyses::myUtils::get_MCindMCVertex(MCVertexObject)")
               .Define("MC_Vertex_ntrk", "FCCAnalyses::myUtils::get_NTracksMCVertex(MCVertexObject)")
               .Define("MC_Vertex_n",    "int(MC_Vertex_x.size())")
               .Define("MC_Vertex_PDG",  "FCCAnalyses::myUtils::get_MCpdgMCVertex(MCVertexObject, Particle)")
               .Define("MC_Vertex_PDGmother",  "FCCAnalyses::myUtils::get_MCpdgMotherMCVertex(MCVertexObject, Particle)")
               .Define("MC_Vertex_PDGgmother", "FCCAnalyses::myUtils::get_MCpdgGMotherMCVertex(MCVertexObject, Particle)")

               #############################################
               ##              Build Reco Vertex          ##
               #############################################
               .Define("VertexObject", "FCCAnalyses::myUtils::get_VertexObject(MCVertexObject,ReconstructedParticles,EFlowTrack_1,MCRecoAssociations0,MCRecoAssociations1)")

               #############################################
               ##          Build PV var and filter        ##
               #############################################
               .Define("EVT_hasPV",    "FCCAnalyses::myUtils::hasPV(VertexObject)")
               .Define("EVT_NtracksPV", "float(FCCAnalyses::myUtils::get_PV_ntracks(VertexObject))")
               .Define("EVT_NVertex",   "float(VertexObject.size())")
               .Filter("EVT_hasPV==1"))


        #############################################
        ##          Build RECO P with PID          ##
        #############################################
        df2 = (df2.Define("RecoPartPID" ,"FCCAnalyses::myUtils::PID(ReconstructedParticles, MCRecoAssociations0,MCRecoAssociations1,Particle)")

               #############################################
               ##    Build RECO P with PID at vertex      ##
               #############################################
               .Define("RecoPartPIDAtVertex" ,"FCCAnalyses::myUtils::get_RP_atVertex(RecoPartPID, VertexObject)")

               #############################################
               ##         Build vertex variables          ##
               #############################################
               .Define("Vertex_x",        "FCCAnalyses::myUtils::get_Vertex_x(VertexObject)")
               .Define("Vertex_y",        "FCCAnalyses::myUtils::get_Vertex_y(VertexObject)")
               .Define("Vertex_z",        "FCCAnalyses::myUtils::get_Vertex_z(VertexObject)")
               .Define("Vertex_xErr",     "FCCAnalyses::myUtils::get_Vertex_xErr(VertexObject)")
               .Define("Vertex_yErr",     "FCCAnalyses::myUtils::get_Vertex_yErr(VertexObject)")
               .Define("Vertex_zErr",     "FCCAnalyses::myUtils::get_Vertex_zErr(VertexObject)")

               .Define("Vertex_chi2",     "FCCAnalyses::myUtils::get_Vertex_chi2(VertexObject)")
               .Define("Vertex_mcind",    "FCCAnalyses::myUtils::get_Vertex_indMC(VertexObject)")
               .Define("Vertex_ind",      "FCCAnalyses::myUtils::get_Vertex_ind(VertexObject)")
               .Define("Vertex_isPV",     "FCCAnalyses::myUtils::get_Vertex_isPV(VertexObject)")
               .Define("Vertex_ntrk",     "FCCAnalyses::myUtils::get_Vertex_ntracks(VertexObject)")
               .Define("Vertex_n",        "int(Vertex_x.size())")
               .Define("Vertex_mass",     "FCCAnalyses::myUtils::get_Vertex_mass(VertexObject,RecoPartPIDAtVertex)")

               .Define("Vertex_d2PV",     "FCCAnalyses::myUtils::get_Vertex_d2PV(VertexObject,-1)")
               .Define("Vertex_d2PVx",    "FCCAnalyses::myUtils::get_Vertex_d2PV(VertexObject,0)")
               .Define("Vertex_d2PVy",    "FCCAnalyses::myUtils::get_Vertex_d2PV(VertexObject,1)")
               .Define("Vertex_d2PVz",    "FCCAnalyses::myUtils::get_Vertex_d2PV(VertexObject,2)")

               .Define("Vertex_d2PVErr",  "FCCAnalyses::myUtils::get_Vertex_d2PVError(VertexObject,-1)")
               .Define("Vertex_d2PVxErr", "FCCAnalyses::myUtils::get_Vertex_d2PVError(VertexObject,0)")
               .Define("Vertex_d2PVyErr", "FCCAnalyses::myUtils::get_Vertex_d2PVError(VertexObject,1)")
               .Define("Vertex_d2PVzErr", "FCCAnalyses::myUtils::get_Vertex_d2PVError(VertexObject,2)")

               .Define("Vertex_d2PVSig",  "Vertex_d2PV/Vertex_d2PVErr")
               .Define("Vertex_d2PVxSig", "Vertex_d2PVx/Vertex_d2PVxErr")
               .Define("Vertex_d2PVySig", "Vertex_d2PVy/Vertex_d2PVyErr")
               .Define("Vertex_d2PVzSig", "Vertex_d2PVz/Vertex_d2PVzErr")

               .Define("Vertex_d2MC",     "FCCAnalyses::myUtils::get_Vertex_d2MC(VertexObject,MCVertexObject,Vertex_mcind,-1)")
               .Define("Vertex_d2MCx",    "FCCAnalyses::myUtils::get_Vertex_d2MC(VertexObject,MCVertexObject,Vertex_mcind,0)")
               .Define("Vertex_d2MCy",    "FCCAnalyses::myUtils::get_Vertex_d2MC(VertexObject,MCVertexObject,Vertex_mcind,1)")
               .Define("Vertex_d2MCz",    "FCCAnalyses::myUtils::get_Vertex_d2MC(VertexObject,MCVertexObject,Vertex_mcind,2)")

               .Define("EVT_dPV2DVmin",   "FCCAnalyses::myUtils::get_dPV2DV_min(Vertex_d2PV)")
               .Define("EVT_dPV2DVmax",   "FCCAnalyses::myUtils::get_dPV2DV_max(Vertex_d2PV)")
               .Define("EVT_dPV2DVave",   "FCCAnalyses::myUtils::get_dPV2DV_ave(Vertex_d2PV)")

               #############################################
               ##        Build Kstz -> KPi  candidates      ##
               #############################################
               .Define("TotalCands",         "FCCAnalyses::myUtils::build_Lb2LMuMu(VertexObject,RecoPartPIDAtVertex,true)") #true for clean (only take vertices with exactly one proton-pion pair)

               #.Define("LCands",         "FCCAnalyses::myUtils::Lb2LMuMuHelper(FCCAnalyses::myUtils::Lb2LMuMuHelper(TotalCands,0), 0)")
               #.Define("LCands",         "FCCAnalyses::myUtils::Lb2LMuMuHelper(FCCAnalyses::myUtils::Lb2LMuMuHelper(TotalCands,0), 1)")
               .Define("LCands",         "(TotalCands.first).first")
               .Define("anti",         "(TotalCands.first).second")


               .Define("MuCands", "TotalCands.second")
               .Define("LCandidates",    "FCCAnalyses::myUtils::add_truthmatched2(LCands, Particle, VertexObject, MCRecoAssociations1, ReconstructedParticles, Particle0)")
               #.Define("LCandidates2",    "FCCAnalyses::myUtils::add_truthmatched2(LCands2, Particle, VertexObject, MCRecoAssociations1, ReconstructedParticles, Particle0)")

               .Define("Muons",    "FCCAnalyses::myUtils::add_truthmatched2(MuCands, Particle, VertexObject, MCRecoAssociations1, ReconstructedParticles, Particle0)")

               #############################################
               ##       Filter Kstz -> KPi candidates      ##
               #############################################
               .Define("EVT_NKPi",              "float(FCCAnalyses::myUtils::getFCCAnalysesComposite_N(LCandidates))")
               .Filter("EVT_NKPi>0"))


        #############################################
        ##              Build the thrust           ##
        #############################################
        df2 = (df2.Define("RP_e",          "FCCAnalyses::ReconstructedParticle::get_e(RecoPartPIDAtVertex)")
               .Define("RP_px",         "FCCAnalyses::ReconstructedParticle::get_px(RecoPartPIDAtVertex)")
               .Define("RP_py",         "FCCAnalyses::ReconstructedParticle::get_py(RecoPartPIDAtVertex)")
               .Define("RP_pz",         "FCCAnalyses::ReconstructedParticle::get_pz(RecoPartPIDAtVertex)")
               .Define("RP_charge",     "FCCAnalyses::ReconstructedParticle::get_charge(RecoPartPIDAtVertex)")

               .Define("EVT_thrustNP",      'FCCAnalyses::Algorithms::minimize_thrust("Minuit2","Migrad")(RP_px, RP_py, RP_pz)')
               .Define("RP_thrustangleNP",  'FCCAnalyses::Algorithms::getAxisCosTheta(EVT_thrustNP, RP_px, RP_py, RP_pz)')
               .Define("EVT_thrust",        'FCCAnalyses::Algorithms::getThrustPointing(1.)(RP_thrustangleNP, RP_e, EVT_thrustNP)') # changed from 'Algorithms::getThrustPointing(RP_thrustangleNP, RP_e, EVT_thrustNP, 1.)' because of https://github.com/HEP-FCC/FCCAnalyses/commit/e9c4787f82505115be0c084da4453031c3cf8fdf
               .Define("RP_thrustangle",    'FCCAnalyses::Algorithms::getAxisCosTheta(EVT_thrust, RP_px, RP_py, RP_pz)')


               #############################################
               ##        Get thrust related values        ##
               #############################################
               ##hemis0 == negative angle == max energy hemisphere if pointing
               ##hemis1 == positive angle == min energy hemisphere if pointing
               .Define("EVT_thrusthemis0_n",    "FCCAnalyses::Algorithms::getAxisN(0)(RP_thrustangle, RP_charge)")
               .Define("EVT_thrusthemis1_n",    "FCCAnalyses::Algorithms::getAxisN(1)(RP_thrustangle, RP_charge)")
               .Define("EVT_thrusthemis0_e",    "FCCAnalyses::Algorithms::getAxisEnergy(0)(RP_thrustangle, RP_charge, RP_e)")
               .Define("EVT_thrusthemis1_e",    "FCCAnalyses::Algorithms::getAxisEnergy(1)(RP_thrustangle, RP_charge, RP_e)")

               .Define("EVT_ThrustEmax_E",         "EVT_thrusthemis0_e.at(0)")
               .Define("EVT_ThrustEmax_Echarged",  "EVT_thrusthemis0_e.at(1)")
               .Define("EVT_ThrustEmax_Eneutral",  "EVT_thrusthemis0_e.at(2)")
               .Define("EVT_ThrustEmax_N",         "float(EVT_thrusthemis0_n.at(0))")
               .Define("EVT_ThrustEmax_Ncharged",  "float(EVT_thrusthemis0_n.at(1))")
               .Define("EVT_ThrustEmax_Nneutral",  "float(EVT_thrusthemis0_n.at(2))")

               .Define("EVT_ThrustEmin_E",         "EVT_thrusthemis1_e.at(0)")
               .Define("EVT_ThrustEmin_Echarged",  "EVT_thrusthemis1_e.at(1)")
               .Define("EVT_ThrustEmin_Eneutral",  "EVT_thrusthemis1_e.at(2)")
               .Define("EVT_ThrustEmin_N",         "float(EVT_thrusthemis1_n.at(0))")
               .Define("EVT_ThrustEmin_Ncharged",  "float(EVT_thrusthemis1_n.at(1))")
               .Define("EVT_ThrustEmin_Nneutral",  "float(EVT_thrusthemis1_n.at(2))"))


        df2 = (df2.Define("Vertex_thrust_angle",   "FCCAnalyses::myUtils::get_Vertex_thrusthemis_angle(VertexObject, RecoPartPIDAtVertex, EVT_thrust)")
               .Define("DVertex_thrust_angle",  "FCCAnalyses::myUtils::get_DVertex_thrusthemis_angle(VertexObject, RecoPartPIDAtVertex, EVT_thrust)")
               ###0 == negative angle==max energy , 1 == positive angle == min energy
               .Define("Vertex_thrusthemis_emin",    "FCCAnalyses::myUtils::get_Vertex_thrusthemis(Vertex_thrust_angle, 1)")
               .Define("Vertex_thrusthemis_emax",    "FCCAnalyses::myUtils::get_Vertex_thrusthemis(Vertex_thrust_angle, 0)")

               .Define("EVT_ThrustEmin_NDV", "float(FCCAnalyses::myUtils::get_Npos(DVertex_thrust_angle))")
               .Define("EVT_ThrustEmax_NDV", "float(FCCAnalyses::myUtils::get_Nneg(DVertex_thrust_angle))")

               .Define("EVT_Thrust_Mag",  "EVT_thrust.at(0)")
               .Define("EVT_Thrust_X",    "EVT_thrust.at(1)")
               .Define("EVT_Thrust_XErr", "EVT_thrust.at(2)")
               .Define("EVT_Thrust_Y",    "EVT_thrust.at(3)")
               .Define("EVT_Thrust_YErr", "EVT_thrust.at(4)")
               .Define("EVT_Thrust_Z",    "EVT_thrust.at(5)")
               .Define("EVT_Thrust_ZErr", "EVT_thrust.at(6)")


               .Define("DV_tracks", "FCCAnalyses::myUtils::get_pseudotrack(VertexObject,RecoPartPIDAtVertex)")

               .Define("DV_d0",            "FCCAnalyses::myUtils::get_trackd0(DV_tracks)")
               .Define("DV_z0",            "FCCAnalyses::myUtils::get_trackz0(DV_tracks)")

               # Build MVA
               #.Define("MVAVec", ROOT.computeModel, ("EVT_ThrustEmin_E",        "EVT_ThrustEmax_E",
               #                                      "EVT_ThrustEmin_Echarged", "EVT_ThrustEmax_Echarged",
               #                                      "EVT_ThrustEmin_Eneutral", "EVT_ThrustEmax_Eneutral",
               #                                      "EVT_ThrustEmin_Ncharged", "EVT_ThrustEmax_Ncharged",
               #                                      "EVT_ThrustEmin_Nneutral", "EVT_ThrustEmax_Nneutral",
               #                                      "EVT_NtracksPV",           "EVT_NVertex",
               #                                      "EVT_NKPi",                "EVT_ThrustEmin_NDV",
               #                                      "EVT_ThrustEmax_NDV",      "EVT_dPV2DVmin",
               #                                      "EVT_dPV2DVmax",           "EVT_dPV2DVave"))
               #.Define("EVT_MVA1", "MVAVec.at(0)")
               #.Filter(MVAFilter)

               .Define("LCandidates_mass",    "FCCAnalyses::myUtils::getFCCAnalysesComposite_mass(LCandidates)")
               # .Define("LCandidates_mass2",    "FCCAnalyses::myUtils::getFCCAnalysesComposite_mass(LCandidates2)")

               .Define("LCandidates_q",       "FCCAnalyses::myUtils::getFCCAnalysesComposite_charge(LCandidates)")
               .Define("LCandidates_vertex",  "FCCAnalyses::myUtils::getFCCAnalysesComposite_vertex(LCandidates)")
               .Define("LCandidates_mcvertex","FCCAnalyses::myUtils::getFCCAnalysesComposite_mcvertex(LCandidates,VertexObject)")
               .Define("LCandidates_truth",   "FCCAnalyses::myUtils::getFCCAnalysesComposite_truthMatch(LCandidates)")
               .Define("LCandidates_px",      "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(LCandidates,0)")
               # .Define("LCandidates_px2",      "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(LCandidates2,0)")

               .Define("LCandidates_py",      "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(LCandidates,1)")
               # .Define("LCandidates_py2",      "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(LCandidates2,1)")

               .Define("LCandidates_pz",      "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(LCandidates,2)")
               # .Define("LCandidates_pz2",      "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(LCandidates2,2)")

               .Define("LCandidates_p",       "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(LCandidates,-1)")
               # .Define("LCandidates_p2",       "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(LCandidates2,-1)")
               .Define("LCandidates_B",       "FCCAnalyses::myUtils::getFCCAnalysesComposite_B(LCandidates, VertexObject, RecoPartPIDAtVertex)")

               .Define("LCandidates_track",   "FCCAnalyses::myUtils::getFCCAnalysesComposite_track(LCandidates, VertexObject)")
               .Define("LCandidates_d0",      "FCCAnalyses::myUtils::get_trackd0(LCandidates_track)")
               .Define("LCandidates_z0",      "FCCAnalyses::myUtils::get_trackz0(LCandidates_track)")

               .Define("LCandidates_anglethrust", "FCCAnalyses::myUtils::getFCCAnalysesComposite_anglethrust(LCandidates, EVT_thrust)")
               .Define("CUT_hasCandEmin",           "FCCAnalyses::myUtils::has_anglethrust_emin(LCandidates_anglethrust)")
               .Filter("CUT_hasCandEmin>0")

               .Define("LCandidates_h1px",   "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(LCandidates, VertexObject, RecoPartPIDAtVertex, 0, 0)")
               .Define("LCandidates_h1py",   "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(LCandidates, VertexObject, RecoPartPIDAtVertex, 0, 1)")
               .Define("LCandidates_h1pz",   "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(LCandidates, VertexObject, RecoPartPIDAtVertex, 0, 2)")
               .Define("LCandidates_h1p",    "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(LCandidates, VertexObject, RecoPartPIDAtVertex, 0, -1)")
               .Define("LCandidates_h1q",    "FCCAnalyses::myUtils::getFCCAnalysesComposite_q(LCandidates, VertexObject, RecoPartPIDAtVertex, 0)")
               .Define("LCandidates_h1m",    "FCCAnalyses::myUtils::getFCCAnalysesComposite_mass(LCandidates, VertexObject, RecoPartPIDAtVertex, 0)")
               .Define("LCandidates_h1type", "FCCAnalyses::myUtils::getFCCAnalysesComposite_type(LCandidates, VertexObject, RecoPartPIDAtVertex, 0)")
               .Define("LCandidates_h1d0",   "FCCAnalyses::myUtils::getFCCAnalysesComposite_d0(LCandidates, VertexObject, 0)")
               .Define("LCandidates_h1z0",   "FCCAnalyses::myUtils::getFCCAnalysesComposite_z0(LCandidates, VertexObject, 0)")

               .Define("LCandidates_h2px",   "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(LCandidates, VertexObject, RecoPartPIDAtVertex, 1, 0)")
               .Define("LCandidates_h2py",   "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(LCandidates, VertexObject, RecoPartPIDAtVertex, 1, 1)")
               .Define("LCandidates_h2pz",   "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(LCandidates, VertexObject, RecoPartPIDAtVertex, 1, 2)")
               .Define("LCandidates_h2p",    "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(LCandidates, VertexObject, RecoPartPIDAtVertex, 1, -1)")
               .Define("LCandidates_h2q",    "FCCAnalyses::myUtils::getFCCAnalysesComposite_q(LCandidates, VertexObject, RecoPartPIDAtVertex, 1)")
               .Define("LCandidates_h2m",    "FCCAnalyses::myUtils::getFCCAnalysesComposite_mass(LCandidates, VertexObject, RecoPartPIDAtVertex, 1)")
               .Define("LCandidates_h2type", "FCCAnalyses::myUtils::getFCCAnalysesComposite_type(LCandidates, VertexObject, RecoPartPIDAtVertex, 1)")
               .Define("LCandidates_h2d0",   "FCCAnalyses::myUtils::getFCCAnalysesComposite_d0(LCandidates, VertexObject, 1)")
               .Define("LCandidates_h2z0",   "FCCAnalyses::myUtils::getFCCAnalysesComposite_z0(LCandidates, VertexObject, 1)")

               .Define("TrueKPiBd_vertex",        "FCCAnalyses::myUtils::get_trueVertex(MCVertexObject,Particle,Particle0, 313, 511)")
               .Define("TrueKPiBd_track",         "FCCAnalyses::myUtils::get_truetrack(TrueKPiBd_vertex, MCVertexObject, Particle)")
               .Define("TrueKPiBd_d0",            "FCCAnalyses::myUtils::get_trackd0(TrueKPiBd_track)")
               .Define("TrueKPiBd_z0",            "FCCAnalyses::myUtils::get_trackz0(TrueKPiBd_track)"))

        #    .Define("muons",  "FCCAnalyses::ReconstructedParticle::get(muon0, ReconstructedParticles)")
        #    .Define("n_muons",  "FCCAnalyses::ReconstructedParticle::get_n(muons)") #count how many muons are in the event in total
        #    .Define("muon_px",      "FCCAnalyses::ReconstructedParticle::get_px(muons)")
        #    .Define("muon_py",      "FCCAnalyses::ReconstructedParticle::get_py(muons)")
        #    .Define("muon_pz",      "FCCAnalyses::ReconstructedParticle::get_pz(muons)")
        #    .Define("muon_p",       "FCCAnalyses::ReconstructedParticle::get_p(muons)")
        #    .Define("muon_mass",       "FCCAnalyses::ReconstructedParticle::get_mass(muons)")
        #    .Define("muon_charge",       "FCCAnalyses::ReconstructedParticle::get_charge(muons)")

        df2 = (df2.Define("Muons_px",      "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(Muons,0)")
               .Define("Muons_py",      "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(Muons,1)")
               .Define("Muons_pz",      "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(Muons,2)")
               .Define("Muons_p",       "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(Muons,-1)")
               .Define("Muons_B",       "FCCAnalyses::myUtils::getFCCAnalysesComposite_B(Muons, VertexObject, RecoPartPIDAtVertex)")
               .Define("Muons_mass",    "FCCAnalyses::myUtils::getFCCAnalysesComposite_mass(Muons)")
               .Define("Muons_q",       "FCCAnalyses::myUtils::getFCCAnalysesComposite_charge(Muons)")

               .Define("Muons_mu1px",   "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(Muons, VertexObject, RecoPartPIDAtVertex, 0, 0)")
               .Define("Muons_mu1py",   "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(Muons, VertexObject, RecoPartPIDAtVertex, 0, 1)")
               .Define("Muons_mu1pz",   "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(Muons, VertexObject, RecoPartPIDAtVertex, 0, 2)")
               .Define("Muons_mu1p",    "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(Muons, VertexObject, RecoPartPIDAtVertex, 0, -1)")
               .Define("Muons_mu1q",    "FCCAnalyses::myUtils::getFCCAnalysesComposite_q(Muons, VertexObject, RecoPartPIDAtVertex, 0)")
               .Define("Muons_mu1m",    "FCCAnalyses::myUtils::getFCCAnalysesComposite_mass(Muons, VertexObject, RecoPartPIDAtVertex, 0)")
               .Define("Muons_mu1type", "FCCAnalyses::myUtils::getFCCAnalysesComposite_type(Muons, VertexObject, RecoPartPIDAtVertex, 0)")
               .Define("Muons_mu1d0",   "FCCAnalyses::myUtils::getFCCAnalysesComposite_d0(Muons, VertexObject, 0)")
               .Define("Muons_mu1z0",   "FCCAnalyses::myUtils::getFCCAnalysesComposite_z0(Muons, VertexObject, 0)")

               .Define("Muons_mu2px",   "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(Muons, VertexObject, RecoPartPIDAtVertex, 1, 0)")
               .Define("Muons_mu2py",   "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(Muons, VertexObject, RecoPartPIDAtVertex, 1, 1)")
               .Define("Muons_mu2pz",   "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(Muons, VertexObject, RecoPartPIDAtVertex, 1, 2)")
               .Define("Muons_mu2p",    "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(Muons, VertexObject, RecoPartPIDAtVertex, 1, -1)")
               .Define("Muons_mu2q",    "FCCAnalyses::myUtils::getFCCAnalysesComposite_q(Muons, VertexObject, RecoPartPIDAtVertex, 1)")
               .Define("Muons_mu2m",    "FCCAnalyses::myUtils::getFCCAnalysesComposite_mass(Muons, VertexObject, RecoPartPIDAtVertex, 1)")
               .Define("Muons_mu2type", "FCCAnalyses::myUtils::getFCCAnalysesComposite_type(Muons, VertexObject, RecoPartPIDAtVertex, 1)")
               .Define("Muons_mu2d0",   "FCCAnalyses::myUtils::getFCCAnalysesComposite_d0(Muons, VertexObject, 1)")
               .Define("Muons_mu2z0",   "FCCAnalyses::myUtils::getFCCAnalysesComposite_z0(Muons, VertexObject, 1)")
               #.Define("Total_mass", "cmath::sqrt((cmath::sqrt(LCandidates_p**2 + LCandidates_mass**2) + cmath::sqrt(Muons_p**2 + Muons_mass**2))**2 - (Muons_px+LCandidates_px)**2 - (Muons_py+LCandidates_py)**2 - (Muons_pz + LCandidates_pz)**2)")

               #############################################
               ##          Aliases for # in python        ##
               #############################################
               # .Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
               # .Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")
               # .Alias("Particle0", "Particle#0.index")
               # .Alias("Particle1", "Particle#1.index")
               #
               # #############################################
               # ##               Build MC Vertex           ##
               # #############################################
               # .Define("MCVertexObject", "myUtils::get_MCVertexObject(Particle, Particle0)")
               # .Define("MC_PDG", "FCCAnalyses::MCParticle::get_pdg(Particle)")
               # .Define("MC_n",   "int(MC_PDG.size())")
               #
               # #############################################
               # ##              Build Reco Vertex          ##
               # #############################################
               # .Define("VertexObject", "myUtils::get_VertexObject(MCVertexObject,ReconstructedParticles,EFlowTrack_1,MCRecoAssociations0,MCRecoAssociations1)")
               #
               # #############################################
               # ##          Build PV var and filter        ##
               # #############################################
               # .Define("EVT_hasPV",     "myUtils::hasPV(VertexObject)")
               # .Define("EVT_NtracksPV", "myUtils::get_PV_ntracks(VertexObject)")
               # .Define("EVT_NVertex",   "VertexObject.size()")
               # .Filter("EVT_hasPV==1")
               #
               # #############################################
               # ##          Build new RecoP with PID       ##
               # #############################################
               # .Define("RecoPartPID" ,"myUtils::PID(ReconstructedParticles, MCRecoAssociations0,MCRecoAssociations1,Particle)")
               #
               # #############################################
               # ##  Build new RecoP with PID at vertex     ##
               # #############################################
               # .Define("RecoPartPIDAtVertex" ,"myUtils::get_RP_atVertex(RecoPartPID, VertexObject)")
               #
               # #############################################
               # ##        Build Lb -> LMuMu candidates      ##
               # #############################################
               # .Define("Lb2LMuMu",          "FCCAnalyses::myUtils::Lb2LMuMu(VertexObject, RecoPartPIDAtVertex)")
               #
               # .Define("Lb",                "Lb2LMuMu.totals")
               # .Define("dimuons",           "Lb2LMuMu.dimuons")
               # .Define("mu1",               "Lb2LMuMu.muons1")
               # .Define("mu2",               "Lb2LMuMu.muons2")
               # .Define("dihadrons",         "Lb2LMuMu.dihadrons")
               # .Define("h1",                "Lb2LMuMu.hadrons1")
               # .Define("h2",                "Lb2LMuMu.hadrons2")
               #
               # .Define("Lb_mass",           "FCCAnalyses::myUtils::getFCCAnalysesComposite_mass(Lb)")
               # .Define("Lb_P",              "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(Lb, -1)")
               # .Define("Lb_Px",             "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(Lb, 0)")
               # .Define("Lb_Py",             "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(Lb, 1)")
               # .Define("Lb_Pz",             "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(Lb, 2)")
               #
               #
               # .Define("dimuon_mass",       "FCCAnalyses::myUtils::getFCCAnalysesComposite_mass(dimuons)")
               # .Define("dimuon_P",          "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(dimuons, -1)")
               # .Define("dimuon_Px",         "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(dimuons, 0)")
               # .Define("dimuon_Py",         "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(dimuons, 1)")
               # .Define("dimuon_Pz",         "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(dimuons, 2)")
               #
               # .Define("mu1_mass",        "FCCAnalyses::myUtils::getFCCAnalysesComposite_mass(mu1)")
               # .Define("mu1_P",           "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(mu1, -1)")
               # .Define("mu1_Px",          "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(mu1, 0)")
               # .Define("mu1_Py",          "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(mu1, 1)")
               # .Define("mu1_Pz",          "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(mu1, 2)")
               # .Define("mu1_E",           "sqrt(pow(mu1_mass, 2) + pow(mu1_P, 2))")
               # .Define("mu1_q",           "FCCAnalyses::myUtils::getFCCAnalysesComposite_charge(mu1)")
               #
               # .Define("mu2_mass",        "FCCAnalyses::myUtils::getFCCAnalysesComposite_mass(mu2)")
               # .Define("mu2_P",           "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(mu2, -1)")
               # .Define("mu2_Px",          "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(mu2, 0)")
               # .Define("mu2_Py",          "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(mu2, 1)")
               # .Define("mu2_Pz",          "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(mu2, 2)")
               # .Define("mu2_E",           "sqrt(pow(mu2_mass, 2) + pow(mu2_P, 2))")
               # .Define("mu2_q",           "FCCAnalyses::myUtils::getFCCAnalysesComposite_charge(mu2)")
               #
               #
               # .Define("dihadron_mass",     "FCCAnalyses::myUtils::getFCCAnalysesComposite_mass(dihadrons)")
               # .Define("dihadron_P",        "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(dihadrons, -1)")
               # .Define("dihadron_Px",       "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(dihadrons, 0)")
               # .Define("dihadron_Py",       "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(dihadrons, 1)")
               # .Define("dihadron_Pz",       "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(dihadrons, 2)")
               #
               # .Define("h1_mass",         "FCCAnalyses::myUtils::getFCCAnalysesComposite_mass(h1)")
               # .Define("h1_P",            "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(h1, -1)")
               # .Define("h1_Px",           "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(h1, 0)")
               # .Define("h1_Py",           "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(h1, 1)")
               # .Define("h1_Pz",           "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(h1, 2)")
               # .Define("h1_E",            "sqrt(pow(h1_mass, 2) + pow(h1_P, 2))")
               # .Define("h1_q",            "FCCAnalyses::myUtils::getFCCAnalysesComposite_charge(h1)")
               #
               # .Define("h2_mass",         "FCCAnalyses::myUtils::getFCCAnalysesComposite_mass(h2)")
               # .Define("h2_P",            "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(h2, -1)")
               # .Define("h2_Px",           "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(h2, 0)")
               # .Define("h2_Py",           "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(h2, 1)")
               # .Define("h2_Pz",           "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(h2, 2)")
               # .Define("h2_E",            "sqrt(pow(h2_mass, 2) + pow(h2_P, 2))")
               # .Define("h2_q",            "FCCAnalyses::myUtils::getFCCAnalysesComposite_charge(h2)")
               #
               # .Define("h1_type",           "Lb2LMuMu.hadron1Types")
               # .Define("h2_type",           "Lb2LMuMu.hadron2Types")
               #
               # .Define("flightDistancesLb",           "Lb2LMuMu.flightDistancesLb")
               # .Define("flightDistancesL",           "Lb2LMuMu.flightDistancesL")
               #
               # .Define("errorsMu",           "Lb2LMuMu.errorsMu")
               # .Define("errorsH",           "Lb2LMuMu.errorsH")

               )
        return df2
    #__________________________________________________________
    #Mandatory: output function, please make sure you return the branchlist as a python list
    def output():

        branchList = ["MC_PDG","MC_M1","MC_M2","MC_n","MC_D1","MC_D2","MC_D3","MC_D4","MC_D5","MC_D6",
        "MC_p","MC_pt","MC_px","MC_py","MC_pz","MC_eta","MC_phi",
        "MC_orivtx_x","MC_orivtx_y","MC_orivtx_z",
        "MC_endvtx_x", "MC_endvtx_y", "MC_endvtx_z", "MC_e","MC_m",
        "EVT_ThrustEmin_E",          "EVT_ThrustEmax_E",
        "EVT_ThrustEmin_Echarged",   "EVT_ThrustEmax_Echarged",
        "EVT_ThrustEmin_Eneutral",   "EVT_ThrustEmax_Eneutral",
        "EVT_ThrustEmin_N",          "EVT_ThrustEmax_N",
        "EVT_ThrustEmin_Ncharged",   "EVT_ThrustEmax_Ncharged",
        "EVT_ThrustEmin_Nneutral",   "EVT_ThrustEmax_Nneutral",
        "EVT_ThrustEmin_NDV",        "EVT_ThrustEmax_NDV",
        "EVT_Thrust_Mag",
        "EVT_Thrust_X",  "EVT_Thrust_XErr",
        "EVT_Thrust_Y",  "EVT_Thrust_YErr",
        "EVT_Thrust_Z",  "EVT_Thrust_ZErr",

        "EVT_NtracksPV", "EVT_NVertex", "EVT_NKPi",

        "EVT_dPV2DVmin","EVT_dPV2DVmax","EVT_dPV2DVave",
        #"EVT_MVA1",

        "MC_Vertex_x", "MC_Vertex_y", "MC_Vertex_z",
        "MC_Vertex_ntrk", "MC_Vertex_n",

        "MC_Vertex_PDG","MC_Vertex_PDGmother","MC_Vertex_PDGgmother",

        "Vertex_x", "Vertex_y", "Vertex_z",
        "Vertex_xErr", "Vertex_yErr", "Vertex_zErr",
        "Vertex_isPV", "Vertex_ntrk", "Vertex_chi2", "Vertex_n",
        "Vertex_thrust_angle", "Vertex_thrusthemis_emin", "Vertex_thrusthemis_emax",

        "Vertex_d2PV", "Vertex_d2PVx", "Vertex_d2PVy", "Vertex_d2PVz",
        "Vertex_d2PVErr", "Vertex_d2PVxErr", "Vertex_d2PVyErr", "Vertex_d2PVzErr",
        "Vertex_mass",
        "DV_d0","DV_z0",

        "TrueKPiBd_vertex", "TrueKPiBd_d0", "TrueKPiBd_z0",

        "LCandidates_mass", "LCandidates_vertex", "LCandidates_mcvertex", "LCandidates_B",
        "LCandidates_truth",
        "LCandidates_px", "LCandidates_py", "LCandidates_pz", "LCandidates_p", "LCandidates_q",
        "LCandidates_d0",  "LCandidates_z0","LCandidates_anglethrust", "anti",

        # "LCandidates_mass2", "LCandidates_px2", "LCandidates_py2","LCandidates_pz2","LCandidates_p2",

        "LCandidates_h1px", "LCandidates_h1py", "LCandidates_h1pz",
        "LCandidates_h1p", "LCandidates_h1q", "LCandidates_h1m", "LCandidates_h1type",
        "LCandidates_h1d0", "LCandidates_h1z0",
        "LCandidates_h2px", "LCandidates_h2py", "LCandidates_h2pz",
        "LCandidates_h2p", "LCandidates_h2q", "LCandidates_h2m", "LCandidates_h2type",
        "LCandidates_h2d0", "LCandidates_h2z0",

        "Muons_px", "Muons_py", "Muons_pz", "Muons_p", "Muons_B", "Muons_mass", "Muons_q",
        "Muons_mu1px", "Muons_mu1py", "Muons_mu1pz", "Muons_mu1p", "Muons_mu1q", "Muons_mu1m", "Muons_mu1type", "Muons_mu1d0", "Muons_mu1z0",
        "Muons_mu2px", "Muons_mu2py", "Muons_mu2pz", "Muons_mu2p", "Muons_mu2q", "Muons_mu2m", "Muons_mu2type", "Muons_mu2d0", "Muons_mu2z0", #"Total_mass"
    ]
        # branchList = ["MC_PDG",
        #
        #               "Lb_mass", "Lb_P", "Lb_Px", "Lb_Py", "Lb_Pz",
        #
        #               "dimuon_mass", "dimuon_P", "dimuon_Px", "dimuon_Py", "dimuon_Pz",
        #               "mu1_mass", "mu1_P", "mu1_Px", "mu1_Py", "mu1_Pz", "mu1_E", "mu1_q",
        #               "mu2_mass", "mu2_P", "mu2_Px", "mu2_Py", "mu2_Pz", "mu2_E", "mu2_q",
        #
        #               "dihadron_mass", "dihadron_P", "dihadron_Px", "dihadron_Py", "dihadron_Pz",
        #               "h1_mass", "h1_P", "h1_Px", "h1_Py", "h1_Pz", "h1_E", "h1_q",
        #               "h2_mass", "h2_P", "h2_Px", "h2_Py", "h2_Pz", "h2_E", "h2_q",
        #               "h1_type", "h2_type",
        #
        #               "flightDistancesLb", "flightDistancesL",
        #               "errorsMu", "errorsH"
        #               ]
        return branchList