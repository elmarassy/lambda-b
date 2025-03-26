#Mandatory: List of processes
processList = {
    # 'p8_ee_Zbb_ecm91_EvtGen_Bd2MuMu':{},
    # 'p8_ee_Zbb_ecm91_EvtGen_Lb2LMuMu':{'fraction':0.05}
    # 'p8_ee_Zbb_ecm91_EvtGen_Lb2LMuMuphsp':{'fraction':1/22}
    'p8_ee_Zss_ecm91':{'fraction':1/4388}
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

               .Define("Lb",                "Lb2LMuMu.totals")
               .Define("dimuons",           "Lb2LMuMu.dimuons")
               .Define("mu1",               "Lb2LMuMu.muons1")
               .Define("mu2",               "Lb2LMuMu.muons2")
               .Define("dihadrons",         "Lb2LMuMu.dihadrons")
               .Define("h1",                "Lb2LMuMu.hadrons1")
               .Define("h2",                "Lb2LMuMu.hadrons2")

               .Define("Lb_mass",           "FCCAnalyses::myUtils::getFCCAnalysesComposite_mass(Lb)")
               .Define("Lb_P",              "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(Lb, -1)")
               .Define("Lb_Px",             "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(Lb, 0)")
               .Define("Lb_Py",             "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(Lb, 1)")
               .Define("Lb_Pz",             "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(Lb, 2)")


               .Define("dimuon_mass",       "FCCAnalyses::myUtils::getFCCAnalysesComposite_mass(dimuons)")
               .Define("dimuon_P",          "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(dimuons, -1)")
               .Define("dimuon_Px",         "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(dimuons, 0)")
               .Define("dimuon_Py",         "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(dimuons, 1)")
               .Define("dimuon_Pz",         "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(dimuons, 2)")

               .Define("mu1_mass",        "FCCAnalyses::myUtils::getFCCAnalysesComposite_mass(mu1)")
               .Define("mu1_P",           "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(mu1, -1)")
               .Define("mu1_Px",          "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(mu1, 0)")
               .Define("mu1_Py",          "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(mu1, 1)")
               .Define("mu1_Pz",          "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(mu1, 2)")
               .Define("mu1_E",           "sqrt(pow(mu1_mass, 2) + pow(mu1_P, 2))")
               .Define("mu1_q",           "FCCAnalyses::myUtils::getFCCAnalysesComposite_charge(mu1)")

               .Define("mu2_mass",        "FCCAnalyses::myUtils::getFCCAnalysesComposite_mass(mu2)")
               .Define("mu2_P",           "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(mu2, -1)")
               .Define("mu2_Px",          "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(mu2, 0)")
               .Define("mu2_Py",          "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(mu2, 1)")
               .Define("mu2_Pz",          "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(mu2, 2)")
               .Define("mu2_E",           "sqrt(pow(mu2_mass, 2) + pow(mu2_P, 2))")
               .Define("mu2_q",           "FCCAnalyses::myUtils::getFCCAnalysesComposite_charge(mu2)")


               .Define("dihadron_mass",     "FCCAnalyses::myUtils::getFCCAnalysesComposite_mass(dihadrons)")
               .Define("dihadron_P",        "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(dihadrons, -1)")
               .Define("dihadron_Px",       "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(dihadrons, 0)")
               .Define("dihadron_Py",       "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(dihadrons, 1)")
               .Define("dihadron_Pz",       "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(dihadrons, 2)")

               .Define("h1_mass",         "FCCAnalyses::myUtils::getFCCAnalysesComposite_mass(h1)")
               .Define("h1_P",            "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(h1, -1)")
               .Define("h1_Px",           "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(h1, 0)")
               .Define("h1_Py",           "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(h1, 1)")
               .Define("h1_Pz",           "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(h1, 2)")
               .Define("h1_E",            "sqrt(pow(h1_mass, 2) + pow(h1_P, 2))")
               .Define("h1_q",            "FCCAnalyses::myUtils::getFCCAnalysesComposite_charge(h1)")

               .Define("h2_mass",         "FCCAnalyses::myUtils::getFCCAnalysesComposite_mass(h2)")
               .Define("h2_P",            "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(h2, -1)")
               .Define("h2_Px",           "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(h2, 0)")
               .Define("h2_Py",           "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(h2, 1)")
               .Define("h2_Pz",           "FCCAnalyses::myUtils::getFCCAnalysesComposite_p(h2, 2)")
               .Define("h2_E",            "sqrt(pow(h2_mass, 2) + pow(h2_P, 2))")
               .Define("h2_q",            "FCCAnalyses::myUtils::getFCCAnalysesComposite_charge(h2)")

               .Define("h1_type",           "Lb2LMuMu.hadron1Types")
               .Define("h2_type",           "Lb2LMuMu.hadron2Types")

               .Define("flightDistancesLb",           "Lb2LMuMu.flightDistancesLb")
               .Define("flightDistancesL",           "Lb2LMuMu.flightDistancesL")

               .Define("errorsMu",           "Lb2LMuMu.errorsMu")
               .Define("errorsH",           "Lb2LMuMu.errorsH")

               )
        return df2
    #__________________________________________________________
    #Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        branchList = ["MC_PDG",

                      "Lb_mass", "Lb_P", "Lb_Px", "Lb_Py", "Lb_Pz",

                      "dimuon_mass", "dimuon_P", "dimuon_Px", "dimuon_Py", "dimuon_Pz",
                      "mu1_mass", "mu1_P", "mu1_Px", "mu1_Py", "mu1_Pz", "mu1_E", "mu1_q",
                      "mu2_mass", "mu2_P", "mu2_Px", "mu2_Py", "mu2_Pz", "mu2_E", "mu2_q",

                      "dihadron_mass", "dihadron_P", "dihadron_Px", "dihadron_Py", "dihadron_Pz",
                      "h1_mass", "h1_P", "h1_Px", "h1_Py", "h1_Pz", "h1_E", "h1_q",
                      "h2_mass", "h2_P", "h2_Px", "h2_Py", "h2_Pz", "h2_E", "h2_q",
                      "h1_type", "h2_type",

                      "flightDistancesLb", "flightDistancesL",
                      "errorsMu", "errorsH"
                      ]
        return branchList