'''
ROOT::VecOps::RVec<FCCAnalysesComposite2> oldHadrons;
            ROOT::VecOps::RVec<FCCAnalysesComposite2> oldMuons;
            ROOT::VecOps::RVec<FCCAnalysesComposite2> oldTotals;
            ROOT::VecOps::RVec<double_t> oldPrimaryLbImpactParameters; //minimal distance between reconstructed lambdaB path and primary vertex
            ROOT::VecOps::RVec<double_t> oldMuonHadronImpactParameters; //minimal distance between dihadron path and dimuon vertex
            ROOT::VecOps::RVec<double_t> oldDisplacementProduct;
            ROOT::VecOps::RVec<int> oldHadron1Types;
            ROOT::VecOps::RVec<int> oldHadron2Types;

            ROOT::VecOps::RVec<FCCAnalysesComposite2> newHadrons;
            ROOT::VecOps::RVec<FCCAnalysesComposite2> newMuons;
            ROOT::VecOps::RVec<FCCAnalysesComposite2> newTotals;
            ROOT::VecOps::RVec<double_t> newPrimaryLbImpactParameters; //minimal distance between reconstructed lambdaB path and primary vertex
            ROOT::VecOps::RVec<double_t> newMuonHadronImpactParameters; //minimal distance between dihadron path and dimuon vertex
            ROOT::VecOps::RVec<double_t> newDisplacementProduct;
            ROOT::VecOps::RVec<int> newHadron1Types;
            ROOT::VecOps::RVec<int> newHadron2Types;
'''

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

               #############################################
               ##          Aliases for # in python        ##
               #############################################
               .Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
               .Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")
               .Alias("Particle0", "Particle#0.index")
               .Alias("Particle1", "Particle#1.index")

               #############################################
               ##               Build MC Vertex           ##
               #############################################
               .Define("MCVertexObject", "myUtils::get_MCVertexObject(Particle, Particle0)")
               .Define("MC_PDG", "FCCAnalyses::MCParticle::get_pdg(Particle)")
               .Define("MC_n",   "int(MC_PDG.size())")

               #############################################
               ##              Build Reco Vertex          ##
               #############################################
               .Define("VertexObject", "myUtils::get_VertexObject(MCVertexObject,ReconstructedParticles,EFlowTrack_1,MCRecoAssociations0,MCRecoAssociations1)")

               #############################################
               ##          Build PV var and filter        ##
               #############################################
               .Define("EVT_hasPV",     "myUtils::hasPV(VertexObject)")
               .Define("EVT_NtracksPV", "myUtils::get_PV_ntracks(VertexObject)")
               .Define("EVT_NVertex",   "VertexObject.size()")
               .Filter("EVT_hasPV==1")

               #############################################
               ##          Build new RecoP with PID       ##
               #############################################
               .Define("RecoPartPID" ,"myUtils::PID(ReconstructedParticles, MCRecoAssociations0,MCRecoAssociations1,Particle)")

               #############################################
               ##  Build new RecoP with PID at vertex     ##
               #############################################
               .Define("RecoPartPIDAtVertex" ,"myUtils::get_RP_atVertex(RecoPartPID, VertexObject)")

               #############################################
               ##        Build Lb -> LMuMu candidates      ##
               #############################################
               .Define("Lb2LMuMu",          "FCCAnalyses::myUtils::Lb2LMuMu(VertexObject, RecoPartPIDAtVertex)")

               .Define("newLb",                "Lb2LMuMu.newTotals")
               .Define("newDimuons",           "Lb2LMuMu.newMuons")
               .Define("newDihadrons",         "Lb2LMuMu.newHadrons")

               .Define("newLb_mass",           "FCCAnalyses::myUtils::getFCCAnalysesComposite_mass(newLb)")
               .Define("newLb_p",              "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(newLb, -1)")
               .Define("newLb_px",             "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(newLb, 0)")
               .Define("newLb_py",             "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(newLb, 1)")
               .Define("newLb_pz",             "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(newLb, 2)")


               .Define("newDimuon_mass",       "FCCAnalyses::myUtils::getFCCAnalysesComposite_mass(newDimuons)")
               .Define("newDimuon_p",          "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(newDimuons, -1)")
               .Define("newDimuon_px",         "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(newDimuons, 0)")
               .Define("newDimuon_py",         "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(newDimuons, 1)")
               .Define("newDimuon_pz",         "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(newDimuons, 2)")


               .Define("newDihadron_mass",     "FCCAnalyses::myUtils::getFCCAnalysesComposite_mass(newDihadrons)")
               .Define("newDihadron_p",        "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(newDihadrons, -1)")
               .Define("newDihadron_px",       "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(newDihadrons, 0)")
               .Define("newDihadron_py",       "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(newDihadrons, 1)")
               .Define("newDihadron_pz",       "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(newDihadrons, 2)")

               .Define("newH1_type",           "Lb2LMuMu.newHadron1Types")
               .Define("newH2_type",           "Lb2LMuMu.newHadron2Types")

               .Define("newFlightDistancesLb",           "Lb2LMuMu.newFlightDistanceLb")
               .Define("newFlightDistancesL",           "Lb2LMuMu.newFlightDistanceL")

               .Define("newPrimaryLbImpactParams",           "Lb2LMuMu.newPrimaryLbImpactParameters")
               .Define("newMuonHadronImpactParams",           "Lb2LMuMu.newMuonHadronImpactParameters")
               .Define("newDisplacementProducts",              "Lb2LMuMu.newDisplacementProduct")

               .Define("oldLb",                "Lb2LMuMu.oldTotals")
               .Define("oldDimuons",           "Lb2LMuMu.oldMuons")
               .Define("oldDihadrons",         "Lb2LMuMu.oldHadrons")

               .Define("oldLb_mass",           "FCCAnalyses::myUtils::getFCCAnalysesComposite_mass(oldLb)")
               .Define("oldLb_p",              "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(oldLb, -1)")
               .Define("oldLb_px",             "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(oldLb, 0)")
               .Define("oldLb_py",             "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(oldLb, 1)")
               .Define("oldLb_pz",             "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(oldLb, 2)")


               .Define("oldDimuon_mass",       "FCCAnalyses::myUtils::getFCCAnalysesComposite_mass(oldDimuons)")
               .Define("oldDimuon_p",          "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(oldDimuons, -1)")
               .Define("oldDimuon_px",         "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(oldDimuons, 0)")
               .Define("oldDimuon_py",         "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(oldDimuons, 1)")
               .Define("oldDimuon_pz",         "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(oldDimuons, 2)")


               .Define("oldDihadron_mass",     "FCCAnalyses::myUtils::getFCCAnalysesComposite_mass(oldDihadrons)")
               .Define("oldDihadron_p",        "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(oldDihadrons, -1)")
               .Define("oldDihadron_px",       "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(oldDihadrons, 0)")
               .Define("oldDihadron_py",       "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(oldDihadrons, 1)")
               .Define("oldDihadron_pz",       "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(oldDihadrons, 2)")

               .Define("oldH1_type",           "Lb2LMuMu.oldHadron1Types")
               .Define("oldH2_type",           "Lb2LMuMu.oldHadron2Types")

               .Define("oldFlightDistancesLb",           "Lb2LMuMu.oldFlightDistanceLb")
               .Define("oldFlightDistancesL",           "Lb2LMuMu.oldFlightDistanceL")

               .Define("oldPrimaryLbImpactParams",           "Lb2LMuMu.oldPrimaryLbImpactParameters")
               .Define("oldMuonHadronImpactParams",           "Lb2LMuMu.oldMuonHadronImpactParameters")
               .Define("oldDisplacementProducts",              "Lb2LMuMu.oldDisplacementProduct"))
        return df2
    #__________________________________________________________
    #Mandatory: output function, please make sure you return the branchlist as a python list
    def output():


        branchList = [
            "newLb_mass", "newLb_p", "newLb_px", "newLb_py", "newLb_pz",

            "newDimuon_mass", "newDimuon_p", "newDimuon_px", "newDimuon_py", "newDimuon_pz",

            "newDihadron_mass", "newDihadron_p", "newDihadron_px", "newDihadron_py", "newDihadron_pz",
            "newH1_type", "newH2_type",

            "newFlightDistancesLb", "newFlightDistancesL",
            "newPrimaryLbImpactParams", "newMuonHadronImpactParams", "newDisplacementProducts",

            "oldLb_mass", "oldLb_p", "oldLb_px", "oldLb_py", "oldLb_pz",

            "oldDimuon_mass", "oldDimuon_p", "oldDimuon_px", "oldDimuon_py", "oldDimuon_pz",

            "oldDihadron_mass", "oldDihadron_p", "oldDihadron_px", "oldDihadron_py", "oldDihadron_pz",
            "oldH1_type", "oldH2_type",

            "oldFlightDistancesLb", "oldFlightDistancesL",
            "oldPrimaryLbImpactParams", "oldMuonHadronImpactParams", "oldDisplacementProducts"
        ]
        return branchList