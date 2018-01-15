#include "interface/SUSY_Upgrade_Skimmer.h"

void SUSY_Upgrade_Skimmer::addBranches(){

    // Event variables
    myskim->Branch("genWeight", &genWeight);
    myskim->Branch("nTot", &nTot);
    myskim->Branch("xs", &xs);

    // Electron variables
    myskim->Branch("el1_pt", &el1_pt);
    myskim->Branch("el1_eta", &el1_eta);
    myskim->Branch("el1_phi", &el1_phi);
    myskim->Branch("el1_q", &el1_q);
    myskim->Branch("el2_pt", &el2_pt);
    myskim->Branch("el2_eta", &el2_eta);
    myskim->Branch("el2_phi", &el2_phi);
    myskim->Branch("el2_q", &el2_q);

    //// Matched truth electron variables
    //myskim->Branch("el1_pt_truth_matched", &el1_pt_truth_matched);
    //myskim->Branch("el1_eta_truth_matched", &el1_eta_truth_matched);
    //myskim->Branch("el1_phi_truth_matched", &el1_phi_truth_matched);
    //myskim->Branch("el1_q_truth_matched", &el1_q_truth_matched);
    //myskim->Branch("el2_pt_truth_matched", &el2_pt_truth_matched);
    //myskim->Branch("el2_eta_truth_matched", &el2_eta_truth_matched);
    //myskim->Branch("el2_phi_truth_matched", &el2_phi_truth_matched);
    //myskim->Branch("el2_q_truth_matched", &el2_q_truth_matched);

    //// Unmatched truth electron variables
    //myskim->Branch("el1_pt_truth", &el1_pt_truth);
    //myskim->Branch("el1_eta_truth", &el1_eta_truth);
    //myskim->Branch("el1_phi_truth", &el1_phi_truth);
    //myskim->Branch("el1_q_truth", &el1_q_truth);
    //myskim->Branch("el2_pt_truth", &el2_pt_truth);
    //myskim->Branch("el2_eta_truth", &el2_eta_truth);
    //myskim->Branch("el2_phi_truth", &el2_phi_truth);
    //myskim->Branch("el2_q_truth", &el2_q_truth);

    // Muon variables
    myskim->Branch("mu1_pt", &mu1_pt);
    myskim->Branch("mu1_eta", &mu1_eta);
    myskim->Branch("mu1_phi", &mu1_phi);
    myskim->Branch("mu1_q", &mu1_q);
    myskim->Branch("mu2_pt", &mu2_pt);
    myskim->Branch("mu2_eta", &mu2_eta);
    myskim->Branch("mu2_phi", &mu2_phi);
    myskim->Branch("mu2_q", &mu2_q);

    //// Matched truth muon variables
    //myskim->Branch("mu1_pt_truth_matched", &mu1_pt_truth_matched);
    //myskim->Branch("mu1_eta_truth_matched", &mu1_eta_truth_matched);
    //myskim->Branch("mu1_phi_truth_matched", &mu1_phi_truth_matched);
    //myskim->Branch("mu1_q_truth_matched", &mu1_q_truth_matched);
    //myskim->Branch("mu2_pt_truth_matched", &mu2_pt_truth_matched);
    //myskim->Branch("mu2_eta_truth_matched", &mu2_eta_truth_matched);
    //myskim->Branch("mu2_phi_truth_matched", &mu2_phi_truth_matched);
    //myskim->Branch("mu2_q_truth_matched", &mu2_q_truth_matched);

    //// Matched truth muon variables
    //myskim->Branch("mu1_pt_truth", &mu1_pt_truth);
    //myskim->Branch("mu1_eta_truth", &mu1_eta_truth);
    //myskim->Branch("mu1_phi_truth", &mu1_phi_truth);
    //myskim->Branch("mu1_q_truth", &mu1_q_truth);
    //myskim->Branch("mu2_pt_truth", &mu2_pt_truth);
    //myskim->Branch("mu2_eta_truth", &mu2_eta_truth);
    //myskim->Branch("mu2_phi_truth", &mu2_phi_truth);
    //myskim->Branch("mu2_q_truth", &mu2_q_truth);

    // Lepton variables
    myskim->Branch("lep1_pt", &lep1_pt);
    myskim->Branch("lep1_eta", &lep1_eta);
    myskim->Branch("lep1_phi", &lep1_phi);
    myskim->Branch("lep1_mass", &lep1_mass);
    myskim->Branch("lep2_pt", &lep2_pt);
    myskim->Branch("lep2_eta", &lep2_eta);
    myskim->Branch("lep2_phi", &lep2_phi);
    myskim->Branch("lep2_mass", &lep2_mass);

    //// Matched truth lepton variables
    //myskim->Branch("lep1_pt_truth_matched", &lep1_pt_truth_matched);
    //myskim->Branch("lep1_eta_truth_matched", &lep1_eta_truth_matched);
    //myskim->Branch("lep1_phi_truth_matched", &lep1_phi_truth_matched);
    //myskim->Branch("lep1_mass_truth_matched", &lep1_mass_truth_matched);
    //myskim->Branch("lep2_pt_truth_matched", &lep2_pt_truth_matched);
    //myskim->Branch("lep2_eta_truth_matched", &lep2_eta_truth_matched);
    //myskim->Branch("lep2_phi_truth_matched", &lep2_phi_truth_matched);
    //myskim->Branch("lep2_mass_truth_matched", &lep2_mass_truth_matched);

    //// Unmatched truth lepton variables
    //myskim->Branch("lep1_pt_truth", &lep1_pt_truth);
    //myskim->Branch("lep1_eta_truth", &lep1_eta_truth);
    //myskim->Branch("lep1_phi_truth", &lep1_phi_truth);
    //myskim->Branch("lep1_mass_truth", &lep1_mass_truth);
    //myskim->Branch("lep2_pt_truth", &lep2_pt_truth);
    //myskim->Branch("lep2_eta_truth", &lep2_eta_truth);
    //myskim->Branch("lep2_phi_truth", &lep2_phi_truth);
    //myskim->Branch("lep2_mass_truth", &lep2_mass_truth);

    // Jet variables
    myskim->Branch("jet1_puppi_pt", &jet1_puppi_pt);
    myskim->Branch("jet1_puppi_eta", &jet1_puppi_eta);
    myskim->Branch("jet1_puppi_phi", &jet1_puppi_phi);
    myskim->Branch("jet1_puppi_q", &jet1_puppi_q);

    // Truth jet variables
    //myskim->Branch("jet1_pt_truth_matched", &jet1_pt_truth_matched);
    //myskim->Branch("jet1_eta_truth_matched", &jet1_eta_truth_matched);
    //myskim->Branch("jet1_phi_truth_matched", &jet1_phi_truth_matched);
    //myskim->Branch("jet1_q_truth_matched", &jet1_q_truth_matched);

    // MET variables
    myskim->Branch("met", &met);
    myskim->Branch("met_eta", &met_eta);
    myskim->Branch("met_phi", &met_phi);
    myskim->Branch("mht", &mht);
    myskim->Branch("mht_eta", &mht_eta);
    myskim->Branch("mht_phi", &mht_phi);
    myskim->Branch("mlt", &mlt);
    myskim->Branch("mlt_eta", &mlt_eta);
    myskim->Branch("mlt_phi", &mlt_phi);
    myskim->Branch("mhlt", &mhlt);
    myskim->Branch("mhlt_eta", &mhlt_eta);
    myskim->Branch("mhlt_phi", &mhlt_phi);

    // Other variables
    myskim->Branch("nLep", &nLep);
    myskim->Branch("nLep_truth", &nLep_truth);
    myskim->Branch("nEl", &nEl);
    myskim->Branch("nMu", &nMu);
    myskim->Branch("nSoftLep", &nSoftLep);
    myskim->Branch("nSoftEl", &nSoftEl);
    myskim->Branch("nSoftMu", &nSoftMu);
    myskim->Branch("nJet", &nJet);
    myskim->Branch("nBJet", &nBJet);
    myskim->Branch("nW", &nW);
    myskim->Branch("nZ", &nZ);
    myskim->Branch("ht", &ht);
    myskim->Branch("hasSFOS", &hasSFOS);
    myskim->Branch("hasSFOS_truth", &hasSFOS_truth);
    myskim->Branch("hasSoftSFOS", &hasSoftSFOS);
    myskim->Branch("hasSoftSFOS_truth", &hasSoftSFOS_truth);
    myskim->Branch("mllMin", &mllMin);
    myskim->Branch("mllMax", &mllMax);
    myskim->Branch("mt1", &mt1);
    myskim->Branch("mt2", &mt2);
    myskim->Branch("pt2l", &pt2l);
    myskim->Branch("drLep1Lep2", &drLep1Lep2);
    myskim->Branch("drLep1Jet1", &drLep1Jet1);
    myskim->Branch("drLep2Jet1", &drLep2Jet1);
    myskim->Branch("mu_pt5to10_origin_nghbr", &mu_pt5to10_origin_nghbr);
    myskim->Branch("mu_pt5to10_origin_cone", &mu_pt5to10_origin_cone);
    myskim->Branch("mu_pt10to30_origin_nghbr", &mu_pt10to30_origin_nghbr);
    myskim->Branch("mu_pt10to30_origin_cone", &mu_pt10to30_origin_cone);
}

void SUSY_Upgrade_Skimmer::clearVectors(){

    // Clear vectors
    el1_pt.clear();
    el1_eta.clear();
    el1_phi.clear();
    el1_q.clear();
    el2_pt.clear();
    el2_eta.clear();
    el2_phi.clear();
    el2_q.clear();
    //el1_pt_truth_matched.clear();
    //el1_eta_truth_matched.clear();
    //el1_phi_truth_matched.clear();
    //el1_q_truth_matched.clear();
    //el2_pt_truth_matched.clear();
    //el2_eta_truth_matched.clear();
    //el2_phi_truth_matched.clear();
    //el2_q_truth_matched.clear();
    //el1_pt_truth.clear();
    //el1_eta_truth.clear();
    //el1_phi_truth.clear();
    //el1_q_truth.clear();
    //el2_pt_truth.clear();
    //el2_eta_truth.clear();
    //el2_phi_truth.clear();
    //el2_q_truth.clear();
    mu1_pt.clear();
    mu1_eta.clear();
    mu1_phi.clear();
    mu1_q.clear();
    mu2_pt.clear();
    mu2_eta.clear();
    mu2_phi.clear();
    mu2_q.clear();
    //mu1_pt_truth_matched.clear();
    //mu1_eta_truth_matched.clear();
    //mu1_phi_truth_matched.clear();
    //mu1_q_truth_matched.clear();
    //mu2_pt_truth_matched.clear();
    //mu2_eta_truth_matched.clear();
    //mu2_phi_truth_matched.clear();
    //mu2_q_truth_matched.clear();
    //mu1_pt_truth.clear();
    //mu1_eta_truth.clear();
    //mu1_phi_truth.clear();
    //mu1_q_truth.clear();
    //mu2_pt_truth.clear();
    //mu2_eta_truth.clear();
    //mu2_phi_truth.clear();
    //mu2_q_truth.clear();
    lep1_pt.clear();
    lep1_eta.clear();
    lep1_phi.clear();
    lep1_mass.clear();
    lep2_pt.clear();
    lep2_eta.clear();
    lep2_phi.clear();
    lep2_mass.clear();
    //lep1_pt_truth_matched.clear();
    //lep1_eta_truth_matched.clear();
    //lep1_phi_truth_matched.clear();
    //lep1_mass_truth_matched.clear();
    //lep2_pt_truth_matched.clear();
    //lep2_eta_truth_matched.clear();
    //lep2_phi_truth_matched.clear();
    //lep2_mass_truth_matched.clear();
    //lep1_pt_truth.clear();
    //lep1_eta_truth.clear();
    //lep1_phi_truth.clear();
    //lep1_mass_truth.clear();
    //lep2_pt_truth.clear();
    //lep2_eta_truth.clear();
    //lep2_phi_truth.clear();
    //lep2_mass_truth.clear();
    jet1_puppi_pt.clear();
    jet1_puppi_eta.clear();
    jet1_puppi_phi.clear();
    jet1_puppi_q.clear();
    //jet1_pt_truth_matched.clear();
    //jet1_eta_truth_matched.clear();
    //jet1_phi_truth_matched.clear();
    //jet1_q_truth_matched.clear();
    mllMin.clear();
    mllMax.clear();
    mt1.clear();
    mt2.clear();
    pt2l.clear();
    drLep1Lep2.clear();
    drLep1Jet1.clear();
    drLep2Jet1.clear();
    mu_pt5to10_origin_nghbr.clear();
    mu_pt5to10_origin_cone.clear();
    mu_pt10to30_origin_nghbr.clear();
    mu_pt10to30_origin_cone.clear();
}

template <typename T> bool SUSY_Upgrade_Skimmer::isIsolated(const T particle){
    return true;
    if (particle->IsolationVarRhoCorr > iso_cut_rel){
        return false;
    }
    if (particle->SumPt > iso_cut_abs){
        return false;
    }
    return true;
}

bool SUSY_Upgrade_Skimmer::isMatched(const GenParticle* truthParticle, const double pt, const double eta, const double phi){
    if (fabs((truthParticle->PT - pt)/truthParticle->PT) > truth_match_diff_pt_rel){ return false; }
    if (DeltaR(truthParticle->Eta, eta, truthParticle->Phi, phi) > truth_match_diff_dr){ return false; }
    return true;
}

template <typename T> bool SUSY_Upgrade_Skimmer::isMatched(const GenParticle* truthParticle, const T particle){
    return isMatched(truthParticle, particle->PT, particle->Eta, particle->Phi);
}

double SUSY_Upgrade_Skimmer::DeltaR(double eta1, double eta2, double phi1, double phi2){
    double dEta = eta1 - eta2;
    double dPhi = DeltaPhi(phi1, phi2);
    return TMath::Sqrt(dEta*dEta + dPhi*dPhi);
}

double SUSY_Upgrade_Skimmer::DeltaPhi(double phi1, double phi2){
    double dPhi = phi1 - phi2;
    while (dPhi  >  TMath::Pi()) dPhi -= 2*TMath::Pi();
    while (dPhi <= -TMath::Pi()) dPhi += 2*TMath::Pi();
    return fabs(dPhi);
}

bool SUSY_Upgrade_Skimmer::isOverlap(const Jet* jet, d_ana::dBranchHandler<Electron>& elecs, d_ana::dBranchHandler<Muon>& muons){
    return false;
    for (size_t i=0; i<elecs.size(); ++i){
        if (DeltaR(jet->Eta, elecs.at(i)->Eta, jet->Phi, elecs.at(i)->Phi) < jet_or_dr){ return true; }
    }
    for (size_t i=0; i<muons.size(); ++i){
        if (DeltaR(jet->Eta, muons.at(i)->Eta, jet->Phi, muons.at(i)->Phi) < jet_or_dr){ return true; }
    }
    return false;
}

void SUSY_Upgrade_Skimmer::passRandomEfficiency(double eff, Float_t*& ppt){
    // Throw random number to decide if lepton can be reconstructed
    // (to mimic real world detector effects)
    if ((rand() % 1000) >= 1000*eff){
        // Oops! We lost that lepton!
        // Instead of deleting this object, we set pT to -1, it will then
        // fail any selection
        *ppt = -1.;
    }
}

void SUSY_Upgrade_Skimmer::effOnTopMuon(d_ana::dBranchHandler<Muon>& muontight){
    // Efficiencies are multiplied ON TOP of efficiencies in Delphes cards
    // In an optimal world, all values here would be 1
    // https://github.com/delphes/delphes/blob/3.4.2pre05/cards/CMS_PhaseII/muonTightId.tcl
    for (size_t i=0; i<muontight.size(); ++i){

        double pt = muontight.at(i)->PT;
        double eta = fabs(muontight.at(i)->Eta);
        Float_t* ppt = &muontight.at(i)->PT;

        if (eta >= 2.5){ passRandomEfficiency(.000, ppt); }
        if (pt < 2.){ passRandomEfficiency(.000, ppt); }
        else if (pt < 4.){
            if      (eta < 0.5){ passRandomEfficiency(.040/.9, ppt); }
            else if (eta < 1.0){ passRandomEfficiency(.050/.9, ppt); }
            else if (eta < 1.5){ passRandomEfficiency(.160/.9, ppt); }
            else if (eta < 2.0){ passRandomEfficiency(.240/.9, ppt); }
            else if (eta < 2.5){ passRandomEfficiency(.230/.9, ppt); }
        }else if (pt < 6.){
            if      (eta < 0.5){ passRandomEfficiency(.430/.9, ppt); }
            else if (eta < 1.0){ passRandomEfficiency(.470/.9, ppt); }
            else if (eta < 1.5){ passRandomEfficiency(.480/.9, ppt); }
            else if (eta < 2.0){ passRandomEfficiency(.440/.9, ppt); }
            else if (eta < 2.5){ passRandomEfficiency(.350/.9, ppt); }
        }else if (pt < 8.){
            if      (eta < 0.5){ passRandomEfficiency(.530/.9, ppt); }
            else if (eta < 1.0){ passRandomEfficiency(.560/.9, ppt); }
            else if (eta < 1.5){ passRandomEfficiency(.530/.9, ppt); }
            else if (eta < 2.0){ passRandomEfficiency(.510/.9, ppt); }
            else if (eta < 2.5){ passRandomEfficiency(.430/.9, ppt); }
        }else if (pt < 10.){
            if      (eta < 0.5){ passRandomEfficiency(.680/.9, ppt); }
            else if (eta < 1.0){ passRandomEfficiency(.690/.9, ppt); }
            else if (eta < 1.5){ passRandomEfficiency(.660/.9, ppt); }
            else if (eta < 2.0){ passRandomEfficiency(.710/.9, ppt); }
            else if (eta < 2.5){ passRandomEfficiency(.570/.9, ppt); }
        }else if (pt < 20.){
            if      (eta < 0.5){ passRandomEfficiency(.810/.9, ppt); }
            else if (eta < 1.0){ passRandomEfficiency(.790/.9, ppt); }
            else if (eta < 1.5){ passRandomEfficiency(.790/.9, ppt); }
            else if (eta < 2.0){ passRandomEfficiency(.770/.9, ppt); }
            else if (eta < 2.5){ passRandomEfficiency(.630/.9, ppt); }
        }else if (pt < 35.){
            if      (eta < 0.5){ passRandomEfficiency(.910/.9, ppt); }
            else if (eta < 1.0){ passRandomEfficiency(.930/.9, ppt); }
            else if (eta < 1.5){ passRandomEfficiency(.890/.9, ppt); }
            else if (eta < 2.0){ passRandomEfficiency(.910/.9, ppt); }
            else if (eta < 2.5){ passRandomEfficiency(.710/.9, ppt); }
        }else if (pt < 50.){
            if      (eta < 0.5){ passRandomEfficiency(.960/.9, ppt); }
            else if (eta < 1.0){ passRandomEfficiency(.940/.9, ppt); }
            else if (eta < 1.5){ passRandomEfficiency(.950/.9, ppt); }
            else if (eta < 2.0){ passRandomEfficiency(.920/.9, ppt); }
            else if (eta < 2.5){ passRandomEfficiency(.760/.9, ppt); }
        }else if (pt < 14000.){
            if      (eta < 0.5){ passRandomEfficiency(.910/.9, ppt); }
            else if (eta < 1.0){ passRandomEfficiency(.910/.9, ppt); }
            else if (eta < 1.5){ passRandomEfficiency(.910/.9, ppt); }
            else if (eta < 2.0){ passRandomEfficiency(.910/.9, ppt); }
            else if (eta < 2.5){ passRandomEfficiency(.810/.9, ppt); }
        }
    }
}

void SUSY_Upgrade_Skimmer::effOnTopElec(d_ana::dBranchHandler<Electron>& elecs){
    // Efficiencies are multiplied ON TOP of efficiencies in Delphes cards
    // In an optimal world, all values here would be 1
    // https://github.com/delphes/delphes/blob/3.4.2pre05/cards/CMS_PhaseII/CMS_PhaseII_Substructure_PIX4022_200PU.tcl#L1137-L1179
    for (size_t i=0; i<elecs.size(); ++i){

        double pt = elecs.at(i)->PT;
        double eta = fabs(elecs.at(i)->Eta);
        Float_t* ppt = &elecs.at(i)->PT;

        if (eta >= 2.5){ passRandomEfficiency(.000, ppt); }
        if (pt < 2.){ passRandomEfficiency(.000, ppt); }
        else if (pt < 4.){
            if      (eta < 0.5){ passRandomEfficiency(.100/.9*.55, ppt); }
            else if (eta < 1.0){ passRandomEfficiency(.130/.9*.55, ppt); }
            else if (eta < 1.5){ passRandomEfficiency(.440/.9*.55, ppt); }
            else if (eta < 2.0){ passRandomEfficiency(.880/.9*.55, ppt); }
            else if (eta < 2.5){ passRandomEfficiency(.740/.9*.55, ppt); }
        }else if (pt < 6.){
            if      (eta < 0.5){ passRandomEfficiency(.850/.9*.55, ppt); }
            else if (eta < 1.0){ passRandomEfficiency(.950/.9*.55, ppt); }
            else if (eta < 1.5){ passRandomEfficiency(.980/.9*.55, ppt); }
            else if (eta < 2.0){ passRandomEfficiency(.940/.9*.55, ppt); }
            else if (eta < 2.5){ passRandomEfficiency(.810/.9*.55, ppt); }
        }else if (pt < 8.){
            if      (eta < 0.5){ passRandomEfficiency(.950/.9*.55, ppt); }
            else if (eta < 1.0){ passRandomEfficiency(.940/.9*.55, ppt); }
            else if (eta < 1.5){ passRandomEfficiency(.950/.9*.55, ppt); }
            else if (eta < 2.0){ passRandomEfficiency(.980/.9*.55, ppt); }
            else if (eta < 2.5){ passRandomEfficiency(.910/.9*.55, ppt); }
        }else if (pt < 10.){
            if      (eta < 0.5){ passRandomEfficiency(.960/.9*.55, ppt); }
            else if (eta < 1.0){ passRandomEfficiency(.950/.9*.55, ppt); }
            else if (eta < 1.5){ passRandomEfficiency(.950/.9*.55, ppt); }
            else if (eta < 2.0){ passRandomEfficiency(.980/.9*.55, ppt); }
            else if (eta < 2.5){ passRandomEfficiency(.790/.9*.55, ppt); }
        }else if (pt < 20.){
            if      (eta < 0.5){ passRandomEfficiency(.960/.9*.55, ppt); }
            else if (eta < 1.0){ passRandomEfficiency(.950/.9*.55, ppt); }
            else if (eta < 1.5){ passRandomEfficiency(.980/.9*.55, ppt); }
            else if (eta < 2.0){ passRandomEfficiency(.980/.9*.55, ppt); }
            else if (eta < 2.5){ passRandomEfficiency(.810/.9*.55, ppt); }
        }else if (pt < 35.){
            if      (eta < 0.5){ passRandomEfficiency(.970/.9*.55, ppt); }
            else if (eta < 1.0){ passRandomEfficiency(.990/.9*.55, ppt); }
            else if (eta < 1.5){ passRandomEfficiency(.980/.9*.55, ppt); }
            else if (eta < 2.0){ passRandomEfficiency(.999/.9*.55, ppt); }
            else if (eta < 2.5){ passRandomEfficiency(.830/.9*.55, ppt); }
        }else if (pt < 50.){
            if      (eta < 0.5){ passRandomEfficiency(.960/.9*.55, ppt); }
            else if (eta < 1.0){ passRandomEfficiency(.970/.9*.55, ppt); }
            else if (eta < 1.5){ passRandomEfficiency(.980/.9*.55, ppt); }
            else if (eta < 2.0){ passRandomEfficiency(.980/.9*.55, ppt); }
            else if (eta < 2.5){ passRandomEfficiency(.870/.9*.55, ppt); }
        }else if (pt < 14000.){
            if      (eta < 0.5){ passRandomEfficiency(.970/.9*.55, ppt); }
            else if (eta < 1.0){ passRandomEfficiency(.970/.9*.55, ppt); }
            else if (eta < 1.5){ passRandomEfficiency(.980/.9*.55, ppt); }
            else if (eta < 2.0){ passRandomEfficiency(.980/.9*.55, ppt); }
            else if (eta < 2.5){ passRandomEfficiency(.870/.9*.55, ppt); }
        }
    }
}

int SUSY_Upgrade_Skimmer::getNghbr(int pid){
    return pid;
    //return fabs(pid);
    if (fabs(pid) <= 3 || pid == 21){
        // Light flavor
        return 0;
    }else if (fabs(pid) == 4){
        // c
        return 1;
    }else if (fabs(pid) == 5){
        // Heavy flavor
        return 2;
    }else if (fabs(pid) == 15){
        // Tau
        return 3;
    }else{
        // Others
        return 4;
    }
}

double SUSY_Upgrade_Skimmer::coneVeto(double pt, double eta, double phi, d_ana::dBranchHandler<GenParticle>& genpart){

    int nghbr = 99;
    double drMin = 99.;
    double drHfMin = 99.;
    double drTauMin = 99.;

    double iso = 0.;

    const double cone = .9;

    // Vectors with particles that have already been filled
    // These are used to not fill the same particle twice
    std::vector<int> filledPid;
    std::vector<double> filledPt;
    std::vector<double> filledEta;
    std::vector<double> filledPhi;

    for (size_t j=0; j<genpart.size(); ++j){

        // Remove muon itself
        if (isMatched(genpart.at(j), pt, eta, phi)){ continue; }

        // Check if this particle has been filled before
        // This needs to be checked since similar copies of truth
        // particles are stored
        bool skipParticle = false;
        for (size_t k=0; k<filledPid.size(); ++k){
            // PID needs to be the same and the particles need to match
            if ((genpart.at(j)->PID == filledPid.at(k)) && isMatched(genpart.at(j), filledPt.at(k), filledEta.at(k), filledPhi.at(k))){
                skipParticle = true;
                break;
            }
        }
        if (skipParticle){ continue; }

        double dr = DeltaR(eta, genpart.at(j)->Eta, phi, genpart.at(j)->Phi);

        // Check b's
        if ((fabs(genpart.at(j)->PID) == 4 || fabs(genpart.at(j)->PID == 5)) && dr < drHfMin){
            drHfMin = dr;
        }

        // Check taus
        if (fabs(genpart.at(j)->PID) == 15 && dr < drTauMin){
            drTauMin = dr;
        }

        if (dr < drMin){
            drMin = dr;
            nghbr = getNghbr(genpart.at(j)->PID);
        }
        if (dr < cone){
            if (pt > 5 && pt < 10){
                mu_pt5to10_origin_cone.push_back(getNghbr(genpart.at(j)->PID));
            }else if (pt < 30){
                mu_pt10to30_origin_cone.push_back(getNghbr(genpart.at(j)->PID));
            }
            filledPid.push_back(genpart.at(j)->PID);
            filledPt.push_back(genpart.at(j)->PT);
            filledEta.push_back(genpart.at(j)->Eta);
            filledPhi.push_back(genpart.at(j)->Phi);
            iso += genpart.at(j)->PT;
        }
    }

    if (pt > 5 && pt < 10){
        mu_pt5to10_origin_nghbr.push_back(nghbr);
    }else if (pt < 30){
        mu_pt10to30_origin_nghbr.push_back(nghbr);
    }

    // Weight to be returned
    // Here I finally sell my scientific soul to the gods of publish or perish
    //double wght = 1./(1. + wght_gen_iso*iso/pt);
    // Or maybe not?
    double wght = 1.;
    if (drHfMin > cone && drTauMin > cone){
        ;
    }else if (drHfMin < drTauMin){
        wght *= wght_hf_veto;
    }else{
        wght *= wght_tau_veto;
    }

    return wght;
}

void SUSY_Upgrade_Skimmer::analyze(size_t childid /* this info can be used for printouts */){

    d_ana::dBranchHandler<HepMCEvent> event(tree(),"Event");
    d_ana::dBranchHandler<Electron> elecs(tree(),"Electron");
    d_ana::dBranchHandler<Muon> muontight(tree(),"MuonTight");
    d_ana::dBranchHandler<Jet> jetpuppi(tree(), "JetPUPPI");
    d_ana::dBranchHandler<MissingET> puppimet(tree(), "PuppiMissingET");
    d_ana::dBranchHandler<GenParticle> genpart(tree(),"Particle");
    d_ana::dBranchHandler<Jet> genjet(tree(),"GenJet");
    //d_ana::dBranchHandler<Weight> rwgt(tree(), "Rwgt");
    //d_ana::dBranchHandler<ScalarHT> scalarht(tree(), "ScalarHT");
    //d_ana::dBranchHandler<Photon> photon(tree(),"Photon");

    //TH1* histo = addPlot(new TH1D("histo", "histo", 24, 0., 60.), "p_{T} (e_{2})", "Events");

    myskim=addTree();
    addBranches();

    // Event loop
    size_t nevents=tree()->entries();
    if(isTestMode()){
        nevents/=100;
    }
    for(size_t eventno=0;eventno<nevents;eventno++){

        // Report status and set event link
        reportStatus(eventno,nevents);
        tree()->setEntry(eventno);

        //if (elecs.size() >= 2){
        //    histo->Fill(elecs.at(1)->PT);
        //}
        //continue;

        clearVectors();

        // Fill event variables
        try{
            genWeight = event.at(0)->Weight;
        }catch (const std::out_of_range& oor){
            std::cerr << "Out of range error when accessing event vector: " << oor.what() << std::endl;
            return;
        }
        nTot = (getXsec()*3000.)/getNorm();
        xs = getXsec();

        //// If there is a tau in the event, skip the event 50 % of the time
        //// This mimics a tau veto, which cannot be implemented in Delphes
        //for (size_t i=0; i<genpart.size(); ++i){
        //    if (fabs(genpart.at(i)->PID) == 15){
        //        genWeight *= wght_tau_veto;
        //        // If we found already a tau, break here, otherwise we get
        //        // multiple shots at rejecting an event if the same tau is
        //        // stored multiple times in truth
        //        break;
        //    }
        //}

        bool evtFnd = false;
        if (event_by_event_comparison){
            // The hardcoded numbers are to be used with the following sample:
            // DYJetsToLL_M-5to50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_200PU/
            for (size_t i=0; i<genpart.size(); ++i){
                if ((fabs(genpart.at(i)->PT - 13.970113) < 1.e-6 && fabs(genpart.at(i)->Eta - 2.467708) < 1.e-6) ||
                    (fabs(genpart.at(i)->PT - 20.706135) < 1.e-6 && fabs(genpart.at(i)->Eta + 2.117317) < 1.e-6) ||
                    (fabs(genpart.at(i)->PT - 19.735157) < 1.e-6 && fabs(genpart.at(i)->Eta - 0.642992) < 1.e-6) ||
                    (fabs(genpart.at(i)->PT - 29.951921) < 1.e-6 && fabs(genpart.at(i)->Eta + 0.918342) < 1.e-6) ||
                    (fabs(genpart.at(i)->PT - 22.212555) < 1.e-6 && fabs(genpart.at(i)->Eta + 9.916969) < 1.e-6) ||
                    (fabs(genpart.at(i)->PT -  6.281808) < 1.e-6 && fabs(genpart.at(i)->Eta + 0.363280) < 1.e-6) ||
                    (fabs(genpart.at(i)->PT - 11.641006) < 1.e-6 && fabs(genpart.at(i)->Eta - 0.820936) < 1.e-6) ||
                    (fabs(genpart.at(i)->PT - 16.675346) < 1.e-6 && fabs(genpart.at(i)->Eta + 1.859797) < 1.e-6) ||
                    (fabs(genpart.at(i)->PT - 10.828916) < 1.e-6 && fabs(genpart.at(i)->Eta + 1.919281) < 1.e-6) ||
                    (fabs(genpart.at(i)->PT - 21.267307) < 1.e-6 && fabs(genpart.at(i)->Eta + 0.721485) < 1.e-6)){
                    evtFnd = true;
                    printf("Foo99 Match: Idx: %lu/%lu; ID: %d; Status: %d; pt: %f; eta: %f; phi:%f\n",
                           i+1, genpart.size(), genpart.at(i)->PID, genpart.at(i)->Status, genpart.at(i)->PT, genpart.at(i)->Eta, genpart.at(i)->Phi);
                }
            }
            if (evtFnd){
                printf("NEW EVENT\n");
                for (size_t i=0; i<genjet.size(); ++i){
                    printf("Foo00 Genjet: Idx: %lu/%lu; pt: %f; eta: %f\n",
                           i+1, genjet.size(), genjet.at(i)->PT, genjet.at(i)->Eta);
                }
                for (size_t i=0; i<genpart.size(); ++i){
                    //if (genpart.at(i)->Status != 1 && genpart.at(i)->Status != 23){ continue; }
                    //if ((fabs(genpart.at(i)->PID) != 11) && (fabs(genpart.at(i)->PID) != 13)){ continue; }
                    printf("Foo01 GenPart: Idx: %lu/%lu; ID: %d; Status: %d; pt: %f; eta: %f\n",
                           i+1, genpart.size(), genpart.at(i)->PID, genpart.at(i)->Status, genpart.at(i)->PT, genpart.at(i)->Eta);
                }
                for (size_t i=0; i<muontight.size(); ++i){
                    printf("Foo03 RecoMu: Idx: %lu/%lu; pt: %f; eta: %f; iso: %d; sumPt: %f\n",
                           i+1, muontight.size(), muontight.at(i)->PT, muontight.at(i)->Eta,
                           isIsolated(muontight.at(i)), muontight.at(i)->SumPt);
                }
                for (size_t i=0; i<elecs.size(); ++i){
                    printf("Foo04 RecoEl: Idx: %lu/%lu; pt: %f; eta: %f; iso: %d; sumPt: %f\n",
                           i+1, elecs.size(), elecs.at(i)->PT, elecs.at(i)->Eta,
                           isIsolated(elecs.at(i)), elecs.at(i)->SumPt);
                }
            }
        }

        // Lepton on-the-fly efficiencies
        effOnTopElec(elecs);
        effOnTopMuon(muontight);

        if (event_by_event_comparison){
            if (evtFnd){
                for (size_t i=0; i<muontight.size(); ++i){
                    printf("Foo05 RecoMu: Idx: %lu/%lu; pt: %f; eta: %f; iso: %d; sumPt: %f\n",
                           i+i, muontight.size(), muontight.at(i)->PT, muontight.at(i)->Eta,
                           isIsolated(muontight.at(i)), muontight.at(i)->SumPt);
                }
                for (size_t i=0; i<elecs.size(); ++i){
                    printf("Foo06 RecoEl: Idx: %lu/%lu; pt: %f; eta: %f; iso: %d; sumPt: %f\n",
                           i+1, elecs.size(), elecs.at(i)->PT, elecs.at(i)->Eta,
                           isIsolated(elecs.at(i)), elecs.at(i)->SumPt);
                }
                std::cout << "==================================" << std::endl;
            }
        }

        // Cutflow variables
        nLep = nEl = nMu = 0;
        nSoftLep = nSoftEl = nSoftMu = 0;
        nJet = nBJet = 0;
        nW = nZ = 0;
        for (size_t i=0; i<elecs.size(); ++i){
            if (elecs.at(i)->PT < el_pt_lo || !isIsolated(elecs.at(i))){ continue; }
            nLep++;
            nEl++;
            if (elecs.at(i)->PT < el_pt_hi){
                nSoftLep++;
                nSoftEl++;
            }
        }
        for (size_t i=0; i<muontight.size(); ++i){
            if (muontight.at(i)->PT < mu_pt_lo || !isIsolated(muontight.at(i))){ continue; }
            nLep++;
            nMu++;
            if (muontight.at(i)->PT < mu_pt_hi){
                nSoftLep++;
                nSoftMu++;
            }
        }
        for (size_t i=0; i<jetpuppi.size(); ++i){
            if (isOverlap(jetpuppi.at(i), elecs, muontight)){ continue; }
            if (jetpuppi.at(i)->PT > jet_pt_lo){
                nJet++;
                if (jetpuppi.at(i)->BTag){
                    nBJet++;
                }
            }
        }

        // Count particles with status 23 for real lepton efficiency histograms
        // Store also the max lepton pT's for these particles, since the genpart
        // vector is not sorted (there's an ambuigity here when several
        // particles have exactly the same pT, but oh well...)
        unsigned int nGenElStatus23 = 0;
        unsigned int nGenMuStatus23 = 0;
        std::vector<double> nGenLepPts;
        if (fill_rle){
            for (size_t i=0; i<genpart.size(); ++i){
                // Fill pT's of all status 1 particles
                if (genpart.at(i)->Status == 1){
                    nGenLepPts.push_back(genpart.at(i)->PT);
                }
                // Count status 23 particles with pT > 2 GeV
                if (genpart.at(i)->Status == 23 && genpart.at(i)->PT > 2){
                    if (fabs(genpart.at(i)->PID) == 11){
                        nGenElStatus23++;
                    }else if (fabs(genpart.at(i)->PID) == 13){
                        nGenMuStatus23++;
                    }
                }
            }
            // Sort vector
            std::sort(nGenLepPts.begin(), nGenLepPts.end());
            // Trim vector to number of status 23 particles
            nGenLepPts.resize(nGenElStatus23+nGenMuStatus23);
        }

        for (size_t i=0; i<genpart.size(); ++i){
            // Only count particles with status between 21 and 29 (to be revised?)
            if (genpart.at(i)->Status >= 21 && genpart.at(i)->Status <= 29){
                if (fabs(genpart.at(i)->PID) == 23){
                    nZ++;
                }else if (fabs(genpart.at(i)->PID) == 24){
                    nW++;
                }
            }

            if (fill_rle){
                // Fill real lepton efficiency histograms as many times as there are
                // Status 23 particles, but use Status 1 particles for this
                // Assume hardest leptons are the ones we are looking for
                if (fabs(genpart.at(i)->PID) == 11){
                    if (!nGenElStatus23){ continue; }
                    if (genpart.at(i)->Status != 1){ continue; }
                    if (std::find(nGenLepPts.begin(), nGenLepPts.end(), genpart.at(i)->PT) == nGenLepPts.end()){ continue; }

                    // Fill denominator histogram
                    rle_el_den->Fill(genpart.at(i)->PT, fabs(genpart.at(i)->Eta));
                    nGenElStatus23--;

                    // Check if we can match that particle
                    for (size_t j=0; j<elecs.size(); ++j){
                        if (!isIsolated(elecs.at(j))){ continue; }
                        // If we make it here, this is a proper reco electron
                        // Fill numerator histogram if it can be matched
                        if (isMatched(genpart.at(i), elecs.at(j))){
                            rle_el_num->Fill(genpart.at(i)->PT, fabs(genpart.at(i)->Eta));
                            break;
                        }
                    }

                }else if (fabs(genpart.at(i)->PID) == 13){
                    if (!nGenMuStatus23){ continue; }
                    if (genpart.at(i)->Status != 1){ continue; }
                    if (std::find(nGenLepPts.begin(), nGenLepPts.end(), genpart.at(i)->PT) == nGenLepPts.end()){ continue; }

                    // Fill denominator histogram
                    rle_mu_den->Fill(genpart.at(i)->PT, fabs(genpart.at(i)->Eta));
                    nGenMuStatus23--;

                    // Check if we can match that particle
                    for (size_t j=0; j<muontight.size(); ++j){
                        if (!isIsolated(muontight.at(j))){ continue; }
                        // If we make it here, this is a proper reco muon
                        // Fill numerator histogram if it can be matched
                        if (isMatched(genpart.at(i), muontight.at(j))){
                            rle_mu_num->Fill(genpart.at(i)->PT, fabs(genpart.at(i)->Eta));
                            break;
                        }
                    }
                }
            }
        }

        // Is a same flavour opposite sign lepton pair present?
        hasSFOS = hasSoftSFOS = false;
        for (size_t i=0; i<elecs.size(); ++i){
            // Only consider isolated particles with minimum pT
            if (elecs.at(i)->PT < el_pt_lo || !isIsolated(elecs.at(i))){ continue; }
            for (size_t j=i+1; j<elecs.size(); ++j){
                if (elecs.at(j)->PT < el_pt_lo || !isIsolated(elecs.at(j))){ continue; }
                if (elecs.at(i)->Charge*elecs.at(j)->Charge < 0){
                    hasSFOS = true;
                    // Check if both particles are soft
                    if (elecs.at(i)->PT < el_pt_hi && elecs.at(j)->PT < el_pt_hi){
                        hasSoftSFOS = true;

                        // Mll for soft SFOS
                        TLorentzVector l1, l2;
                        l1.SetPtEtaPhiM(elecs.at(i)->PT, elecs.at(i)->Eta, elecs.at(i)->Phi, mass_el);
                        l2.SetPtEtaPhiM(elecs.at(j)->PT, elecs.at(j)->Eta, elecs.at(j)->Phi, mass_el);
                        double mll = (l1+l2).M();
                        if (mllMin.size() == 0){
                            mllMin.push_back(mll);
                        }else if (mllMin.at(0) > mll){
                            mllMin.at(0) = mll;
                        }
                        if (mllMax.size() == 0){
                            mllMax.push_back(mll);
                        }else if (mllMax.at(0) < mll){
                            mllMax.at(0) = mll;
                        }
                    }
                }
            }
        }
        for (size_t i=0; i<muontight.size(); ++i){
            // Only consider isolated particles with minimum pT
            if (muontight.at(i)->PT < mu_pt_lo || !isIsolated(muontight.at(i))){ continue; }
            for (size_t j=i+1; j<muontight.size(); ++j){
                if (muontight.at(j)->PT < mu_pt_lo || !isIsolated(muontight.at(j))){ continue; }
                if (muontight.at(i)->Charge*muontight.at(j)->Charge < 0){
                    hasSFOS = true;
                    // Check if both particles are soft
                    if (muontight.at(i)->PT < mu_pt_hi && muontight.at(j)->PT < mu_pt_hi){
                        hasSoftSFOS = true;

                        // Mll for SFOS (only first found pair considered; this could potentially overwrite e+e- Mll)
                        TLorentzVector l1, l2;
                        l1.SetPtEtaPhiM(muontight.at(i)->PT, muontight.at(i)->Eta, muontight.at(i)->Phi, mass_mu);
                        l2.SetPtEtaPhiM(muontight.at(j)->PT, muontight.at(j)->Eta, muontight.at(j)->Phi, mass_mu);
                        double mll = (l1+l2).M();
                        if (mllMin.size() == 0){
                            mllMin.push_back(mll);
                        }else if (mllMin.at(0) > mll){
                            mllMin.at(0) = mll;
                        }
                        if (mllMax.size() == 0){
                            mllMax.push_back(mll);
                        }else if (mllMax.at(0) < mll){
                            mllMax.at(0) = mll;
                        }
                    }
                }
            }
        }

        // Skim
        if (nLep < 2){ continue; }
        if (nSoftLep < 2){ continue; }
        if (!hasSoftSFOS){ continue; }

        // Fill electrons
        for (size_t i=0; i<elecs.size(); ++i){
            if (elecs.at(i)->PT < el_pt_lo || !isIsolated(elecs.at(i))){ continue; }
            // Fill it if it hasn't been filled
            if (el1_pt.size() == 0){
                el1_pt.push_back(elecs.at(i)->PT);
                el1_eta.push_back(elecs.at(i)->Eta);
                el1_phi.push_back(elecs.at(i)->Phi);
                el1_q.push_back(elecs.at(i)->Charge);
            }else{
                // Otherwise fill second vector and break
                el2_pt.push_back(elecs.at(i)->PT);
                el2_eta.push_back(elecs.at(i)->Eta);
                el2_phi.push_back(elecs.at(i)->Phi);
                el2_q.push_back(elecs.at(i)->Charge);
                break;
            }
        }

        //// Fill truth electrons
        //for (size_t i=0; i<el1_pt.size(); ++i){
        //    for (size_t j=0; j<genpart.size(); ++j){
        //        if (genpart.at(j)->Status != 1){ continue; }
        //        if (fabs(genpart.at(j)->PID) != 11){ continue; }
        //        // Truth matching
        //        if (!isMatched(genpart.at(j), el1_pt.at(i), el1_eta.at(i), el1_phi.at(i))){ continue; }
        //        // If we make it here, the particle has been matched
        //        el1_pt_truth_matched.push_back(genpart.at(j)->PT);
        //        el1_eta_truth_matched.push_back(genpart.at(j)->Eta);
        //        el1_phi_truth_matched.push_back(genpart.at(j)->Phi);
        //        el1_q_truth_matched.push_back(genpart.at(j)->Charge);
        //        break;
        //    }
        //}
        //for (size_t i=0; i<el2_pt.size(); ++i){
        //    for (size_t j=0; j<genpart.size(); ++j){
        //        if (genpart.at(j)->Status != 1){ continue; }
        //        if (fabs(genpart.at(j)->PID) != 11){ continue; }
        //        // Truth matching
        //        if (!isMatched(genpart.at(j), el2_pt.at(i), el2_eta.at(i), el2_phi.at(i))){ continue; }
        //        // If we make it here, the particle has been matched
        //        el2_pt_truth_matched.push_back(genpart.at(j)->PT);
        //        el2_eta_truth_matched.push_back(genpart.at(j)->Eta);
        //        el2_phi_truth_matched.push_back(genpart.at(j)->Phi);
        //        el2_q_truth_matched.push_back(genpart.at(j)->Charge);
        //        break;
        //    }
        //}

        // Fill muons
        for (size_t i=0; i<muontight.size(); ++i){
            if (muontight.at(i)->PT < mu_pt_lo || !isIsolated(muontight.at(i))){ continue; }
            // Fill it if it hasn't been filled
            if (mu1_pt.size() == 0){
                mu1_pt.push_back(muontight.at(i)->PT);
                mu1_eta.push_back(muontight.at(i)->Eta);
                mu1_phi.push_back(muontight.at(i)->Phi);
                mu1_q.push_back(muontight.at(i)->Charge);
            }else{
                // Otherwise fill second vector and break
                mu2_pt.push_back(muontight.at(i)->PT);
                mu2_eta.push_back(muontight.at(i)->Eta);
                mu2_phi.push_back(muontight.at(i)->Phi);
                mu2_q.push_back(muontight.at(i)->Charge);
                break;
            }
        }

        //// Fill truth muons
        //for (size_t i=0; i<mu1_pt.size(); ++i){
        //    for (size_t j=0; j<genpart.size(); ++j){
        //        if (genpart.at(j)->Status != 1){ continue; }
        //        if (fabs(genpart.at(j)->PID) != 13){ continue; }
        //        // Truth matching
        //        if (!isMatched(genpart.at(j), mu1_pt.at(i), mu1_eta.at(i), mu1_phi.at(i))){ continue; }
        //        // If we make it here, the particle has been matched
        //        mu1_pt_truth_matched.push_back(genpart.at(j)->PT);
        //        mu1_eta_truth_matched.push_back(genpart.at(j)->Eta);
        //        mu1_phi_truth_matched.push_back(genpart.at(j)->Phi);
        //        mu1_q_truth_matched.push_back(genpart.at(j)->Charge);
        //        break;
        //    }
        //}
        //for (size_t i=0; i<mu2_pt.size(); ++i){
        //    for (size_t j=0; j<genpart.size(); ++j){
        //        if (genpart.at(j)->Status != 1){ continue; }
        //        if (fabs(genpart.at(j)->PID) != 13){ continue; }
        //        // Truth matching
        //        if (!isMatched(genpart.at(j), mu2_pt.at(i), mu2_eta.at(i), mu2_phi.at(i))){ continue; }
        //        // If we make it here, the particle has been matched
        //        mu2_pt_truth_matched.push_back(genpart.at(j)->PT);
        //        mu2_eta_truth_matched.push_back(genpart.at(j)->Eta);
        //        mu2_phi_truth_matched.push_back(genpart.at(j)->Phi);
        //        mu2_q_truth_matched.push_back(genpart.at(j)->Charge);
        //        break;
        //    }
        //}

        // Guess origin of leptons
        for (size_t i=0; i<mu1_pt.size(); ++i){
            genWeight *= coneVeto(mu1_pt.at(i), mu1_eta.at(i), mu1_phi.at(i), genpart);
        }
        for (size_t i=0; i<mu2_pt.size(); ++i){
            genWeight *= coneVeto(mu2_pt.at(i), mu2_eta.at(i), mu2_phi.at(i), genpart);
        }

        // Fill leptons
        // Put pT and eta into vector of vector for sorting
        std::vector<std::vector<double>> lepvec;
        if (el1_pt.size() != 0){
            lepvec.push_back({el1_pt.at(0), el1_eta.at(0), el1_phi.at(0), mass_el});
        }
        if (el2_pt.size() != 0){
            lepvec.push_back({el2_pt.at(0), el2_eta.at(0), el2_phi.at(0), mass_el});
        }
        if (mu1_pt.size() != 0){
            lepvec.push_back({mu1_pt.at(0), mu1_eta.at(0), mu1_phi.at(0), mass_mu});
        }
        if (mu2_pt.size() != 0){
            lepvec.push_back({mu2_pt.at(0), mu2_eta.at(0), mu2_phi.at(0), mass_mu});
        }
        // By definition, this sorts by the first element of the vector (in this case pT)
        if (lepvec.size() > 1){
            std::sort(begin(lepvec), end(lepvec));
            std::reverse(begin(lepvec), end(lepvec));
        }
        if (lepvec.size() >= 1){
            lep1_pt.push_back(lepvec[0][0]);
            lep1_eta.push_back(lepvec[0][1]);
            lep1_phi.push_back(lepvec[0][2]);
            lep1_mass.push_back(lepvec[0][3]);
        }
        if (lepvec.size() >= 2){
            lep2_pt.push_back(lepvec[1][0]);
            lep2_eta.push_back(lepvec[1][1]);
            lep2_phi.push_back(lepvec[1][2]);
            lep2_mass.push_back(lepvec[1][3]);
        }
        lepvec.clear();

        //// Fill matched truth leptons
        //// Put pT and eta into vector of vector for sorting
        //std::vector<std::vector<double>> lepvec_truth_matched;
        //if (el1_pt_truth_matched.size() != 0){
        //    lepvec_truth_matched.push_back({el1_pt_truth_matched.at(0), el1_eta_truth_matched.at(0), el1_phi_truth_matched.at(0), mass_el});
        //}
        //if (el2_pt_truth_matched.size() != 0){
        //    lepvec_truth_matched.push_back({el2_pt_truth_matched.at(0), el2_eta_truth_matched.at(0), el2_phi_truth_matched.at(0), mass_el});
        //}
        //if (mu1_pt_truth_matched.size() != 0){
        //    lepvec_truth_matched.push_back({mu1_pt_truth_matched.at(0), mu1_eta_truth_matched.at(0), mu1_phi_truth_matched.at(0), mass_mu});
        //}
        //if (mu2_pt_truth_matched.size() != 0){
        //    lepvec_truth_matched.push_back({mu2_pt_truth_matched.at(0), mu2_eta_truth_matched.at(0), mu2_phi_truth_matched.at(0), mass_mu});
        //}
        //// By definition, this sorts by the first element of the vector (in this case pT)
        //if (lepvec_truth_matched.size() > 1){
        //    std::sort(begin(lepvec_truth_matched), end(lepvec_truth_matched));
        //    std::reverse(begin(lepvec_truth_matched), end(lepvec_truth_matched));
        //}
        //if (lepvec_truth_matched.size() >= 1){
        //    lep1_pt_truth_matched.push_back(lepvec_truth_matched[0][0]);
        //    lep1_eta_truth_matched.push_back(lepvec_truth_matched[0][1]);
        //    lep1_phi_truth_matched.push_back(lepvec_truth_matched[0][2]);
        //    lep1_mass_truth_matched.push_back(lepvec_truth_matched[0][3]);
        //}
        //if (lepvec_truth_matched.size() >= 2){
        //    lep2_pt_truth_matched.push_back(lepvec_truth_matched[1][0]);
        //    lep2_eta_truth_matched.push_back(lepvec_truth_matched[1][1]);
        //    lep2_phi_truth_matched.push_back(lepvec_truth_matched[1][2]);
        //    lep2_mass_truth_matched.push_back(lepvec_truth_matched[1][3]);
        //}
        //lepvec_truth_matched.clear();

        // FIXME: The genpart vector is not sorted by pt, hence I'm filling here
        // random particles instead of the hardest ones; later I do sort them,
        // but at this point I might already have missed the hardest leptons
        //// Fill unmatched truth leptons
        //for (size_t j=0; j<genpart.size(); ++j){
        //    // Select only particles from hard process
        //    if (genpart.at(j)->Status != 1){ continue; }
        //    // Electrons
        //    if (fabs(genpart.at(j)->PID) == 11){
        //        if (el1_pt_truth.size() == 0){
        //            el1_pt_truth.push_back(genpart.at(j)->PT);
        //            el1_eta_truth.push_back(genpart.at(j)->Eta);
        //            el1_phi_truth.push_back(genpart.at(j)->Phi);
        //            el1_q_truth.push_back(genpart.at(j)->Charge);
        //        }else if (el2_pt_truth.size() == 0){
        //            el2_pt_truth.push_back(genpart.at(j)->PT);
        //            el2_eta_truth.push_back(genpart.at(j)->Eta);
        //            el2_phi_truth.push_back(genpart.at(j)->Phi);
        //            el2_q_truth.push_back(genpart.at(j)->Charge);
        //        }
        //    }
        //    //Muons
        //    else if (fabs(genpart.at(j)->PID) == 13){
        //        if (mu1_pt_truth.size() == 0){
        //            mu1_pt_truth.push_back(genpart.at(j)->PT);
        //            mu1_eta_truth.push_back(genpart.at(j)->Eta);
        //            mu1_phi_truth.push_back(genpart.at(j)->Phi);
        //            mu1_q_truth.push_back(genpart.at(j)->Charge);
        //        }else if (mu2_pt_truth.size() == 0){
        //            mu2_pt_truth.push_back(genpart.at(j)->PT);
        //            mu2_eta_truth.push_back(genpart.at(j)->Eta);
        //            mu2_phi_truth.push_back(genpart.at(j)->Phi);
        //            mu2_q_truth.push_back(genpart.at(j)->Charge);
        //        }
        //    }
        //}

        //// Fill unmatched truth leptons
        //// Put pT and eta into vector of vector for sorting
        //hasSFOS_truth = hasSoftSFOS_truth = false;
        //std::vector<std::vector<double>> lepvec_truth;
        //if (el1_pt_truth.size() != 0){
        //    lepvec_truth.push_back({el1_pt_truth.at(0), el1_eta_truth.at(0), el1_phi_truth.at(0), mass_el});
        //}
        //if (el2_pt_truth.size() != 0){
        //    lepvec_truth.push_back({el2_pt_truth.at(0), el2_eta_truth.at(0), el2_phi_truth.at(0), mass_el});
        //    // Do we have a SFOS in truth?
        //    if (el1_q_truth.at(0) * el2_q_truth.at(0) < 0){
        //        hasSFOS_truth = true;
        //        // Is the SFOS soft?
        //        if (el1_pt_truth.at(0) > 5 && el1_pt_truth.at(0) < 30 && el2_pt_truth.at(0) > 5 && el2_pt_truth.at(0) < 30){
        //            hasSoftSFOS_truth = true;
        //        }
        //    }
        //}
        //if (mu1_pt_truth.size() != 0){
        //    lepvec_truth.push_back({mu1_pt_truth.at(0), mu1_eta_truth.at(0), mu1_phi_truth.at(0), mass_mu});
        //}
        //if (mu2_pt_truth.size() != 0){
        //    lepvec_truth.push_back({mu2_pt_truth.at(0), mu2_eta_truth.at(0), mu2_phi_truth.at(0), mass_mu});
        //    // Do we have a SFOS in truth?
        //    if (mu1_q_truth.at(0) * mu2_q_truth.at(0) < 0){
        //        hasSFOS_truth = true;
        //        // Is the SFOS soft?
        //        if (mu1_pt_truth.at(0) > 5 && mu1_pt_truth.at(0) < 30 && mu2_pt_truth.at(0) > 5 && mu2_pt_truth.at(0) < 30){
        //            hasSoftSFOS_truth = true;
        //        }
        //    }
        //}
        //nLep_truth = lepvec_truth.size();
        //// By definition, this sorts by the first element of the vector (in this case pT)
        //if (lepvec_truth.size() > 1){
        //    std::sort(begin(lepvec_truth), end(lepvec_truth));
        //    std::reverse(begin(lepvec_truth), end(lepvec_truth));
        //}
        //if (lepvec_truth.size() >= 1){
        //    lep1_pt_truth.push_back(lepvec_truth[0][0]);
        //    lep1_eta_truth.push_back(lepvec_truth[0][1]);
        //    lep1_phi_truth.push_back(lepvec_truth[0][2]);
        //    lep1_mass_truth.push_back(lepvec_truth[0][3]);
        //}
        //if (lepvec_truth.size() >= 2){
        //    lep2_pt_truth.push_back(lepvec_truth[1][0]);
        //    lep2_eta_truth.push_back(lepvec_truth[1][1]);
        //    lep2_phi_truth.push_back(lepvec_truth[1][2]);
        //    lep2_mass_truth.push_back(lepvec_truth[1][3]);
        //}
        //lepvec_truth.clear();

        // Fill jets
        for (size_t i=0; i<jetpuppi.size(); ++i){
            if (jetpuppi.at(i)->PT < jet_pt_lo){ continue; }
            if (isOverlap(jetpuppi.at(i), elecs, muontight)){ continue; }
            jet1_puppi_pt.push_back(jetpuppi.at(0)->PT);
            jet1_puppi_eta.push_back(jetpuppi.at(0)->Eta);
            jet1_puppi_phi.push_back(jetpuppi.at(0)->Phi);
            jet1_puppi_q.push_back(jetpuppi.at(0)->Charge);
            break;
        }

        //// Fill truth jets
        //for (size_t i=0; i<jet1_puppi_pt.size(); ++i){
        //    for (size_t j=0; j<genjet.size(); ++j){
        //        // Truth matching
        //        if (fabs(genjet.at(j)->PT - jet1_puppi_pt.at(i)) > truth_match_diff_pt){ continue; }
        //        if (fabs(genjet.at(j)->Eta - jet1_puppi_eta.at(i) > truth_match_diff_eta)){ continue; }
        //        // If we make it here, the particle has been matched
        //        jet1_pt_truth_matched.push_back(genjet.at(j)->PT);
        //        jet1_eta_truth_matched.push_back(genjet.at(j)->Eta);
        //        jet1_phi_truth_matched.push_back(genjet.at(j)->Phi);
        //        jet1_q_truth_matched.push_back(genjet.at(j)->Charge);
        //        break;
        //    }
        //}

        // Fill MET
        try{
            met = puppimet.at(0)->MET;
            met_eta = puppimet.at(0)->Eta;
            met_phi = puppimet.at(0)->Phi;
        }catch (const std::out_of_range& oor){
            std::cerr << "Out of range error when accessing MET vector: " << oor.what() << std::endl;
            return;
        }

        // Fill poor man's MET
        TLorentzVector mht4, mlt4, mhlt4;
        for (size_t i=0; i<jetpuppi.size(); ++i){
            if (jetpuppi.at(i)->PT < jet_pt_lo){ continue; }
            if (isOverlap(jetpuppi.at(i), elecs, muontight)){ continue; }
            TLorentzVector j4;
            j4.SetPtEtaPhiM(jetpuppi.at(i)->PT, jetpuppi.at(i)->Eta, jetpuppi.at(i)->Phi, jetpuppi.at(i)->Mass);
            mht4 += j4;
            mhlt4 += j4;
        }
        for (size_t i=0; i<muontight.size(); ++i){
            if (muontight.at(i)->PT < mu_pt_lo || !isIsolated(muontight.at(i))){ continue; }
            TLorentzVector m4;
            m4.SetPtEtaPhiM(muontight.at(i)->PT, muontight.at(i)->Eta, muontight.at(i)->Phi, mass_mu);
            mlt4 += m4;
            mhlt4 += m4;
        }
        for (size_t i=0; i<elecs.size(); ++i){
            if (elecs.at(i)->PT < el_pt_lo || !isIsolated(elecs.at(i))){ continue; }
            TLorentzVector e4;
            e4.SetPtEtaPhiM(elecs.at(i)->PT, elecs.at(i)->Eta, elecs.at(i)->Phi, mass_el);
            mlt4 += e4;
            mhlt4 += e4;
        }
        mht = mht4.Pt();
        mht_eta = mht4.Eta();
        mht_phi = mht4.Phi();
        mlt = mlt4.Pt();
        mlt_eta = mlt4.Eta();
        mlt_phi = mlt4.Phi();
        mhlt = mhlt4.Pt();
        mhlt_eta = mhlt4.Eta();
        mhlt_phi = mhlt4.Phi();

        // Fill HT
        ht = 0.;
        for (size_t i=0; i<jetpuppi.size(); ++i){
            if (isOverlap(jetpuppi.at(i), elecs, muontight)){ continue; }
            if (jetpuppi.at(i)->PT > jet_pt_lo){
                ht += jetpuppi.at(i)->PT;
            }
        }

        // dr variables
        if (mu2_pt.size() != 0 && jet1_puppi_pt.size() != 0){
            drLep1Lep2.push_back(DeltaR(mu1_eta.at(0), mu2_eta.at(0), mu1_phi.at(0), mu2_phi.at(0)));
            drLep1Jet1.push_back(DeltaR(mu1_eta.at(0), jet1_puppi_eta.at(0), mu1_phi.at(0), jet1_puppi_phi.at(0)));
            drLep2Jet1.push_back(DeltaR(mu1_eta.at(0), jet1_puppi_eta.at(0), mu1_phi.at(0), jet1_puppi_phi.at(0)));
        }
        if (el2_pt.size() != 0 && jet1_puppi_pt.size() != 0){
            drLep1Lep2.push_back(DeltaR(el1_eta.at(0), el2_eta.at(0), el1_phi.at(0), el2_phi.at(0)));
            drLep1Jet1.push_back(DeltaR(el1_eta.at(0), jet1_puppi_eta.at(0), el1_phi.at(0), jet1_puppi_phi.at(0)));
            drLep2Jet1.push_back(DeltaR(el1_eta.at(0), jet1_puppi_eta.at(0), el1_phi.at(0), jet1_puppi_phi.at(0)));
        }

        // Transverse mass of leading two leptons
        if (lep1_pt.size() >= 1){
            mt1.push_back(std::sqrt(2*lep1_pt.at(0)*met*(1-std::cos(lep1_phi.at(0)-met_phi))));
        }
        if (lep2_pt.size() >= 1){
            mt2.push_back(std::sqrt(2*lep2_pt.at(0)*met*(1-std::cos(lep2_phi.at(0)-met_phi))));
            // Fill pt of two leptons
            TLorentzVector l1, l2;
            l1.SetPtEtaPhiM(lep1_pt[0], lep1_eta[0], lep1_phi[0], lep1_mass[0]);
            l2.SetPtEtaPhiM(lep2_pt[0], lep2_eta[0], lep2_phi[0], lep2_mass[0]);
            pt2l.push_back((l1 + l2).Pt());
        }

        myskim->Fill();
    }

    //mu1_pt_origin_nghbr->Write();
    //mu1_pt_origin_cone->Write();
    //mu2_pt_origin_nghbr->Write();
    //mu2_pt_origin_cone->Write();

    if (fill_rle){
        rle_el_num->Write();
        rle_el_den->Write();
        rle_mu_num->Write();
        rle_mu_den->Write();
        TH2D* rle_el = (TH2D*)rle_el_num->Clone();
        TH2D* rle_mu = (TH2D*)rle_mu_num->Clone();
        rle_el->SetNameTitle("rle_el", "rle_el");
        rle_mu->SetNameTitle("rle_mu", "rle_mu");
        rle_el->Divide(rle_el_den);
        rle_mu->Divide(rle_mu_den);
        rle_el->Write();
        rle_mu->Write();
    }

    /*
     * Must be called in the end, takes care of thread-safe writeout and
     * call-back to the parent process
     */
    processEndFunction();
}



void SUSY_Upgrade_Skimmer::postProcess(){
    /*
     * This function can be used to analyse the output histograms, e.g. extract a signal contribution etc.
     * The function can also be called directly on an output file with the histograms, if
     * RunOnOutputOnly = true
     * is set in the analyser's config file
     *
     * This function also represents an example of how the output of the analyser can be
     * read-back in an external program.
     * Just include the sampleCollection.h header and follow the procedure below
     *
     */

    /*
     * Here, the input file to the extraction of parameters from the histograms is the output file
     * of the parallelised analysis.
     * The sampleCollection class can also be used externally for accessing the output consistently
     */
    //d_ana::sampleCollection samplecoll;
    //samplecoll.readFromFile(getOutPath());

    //std::vector<TString> alllegends = samplecoll.listAllLegends();

    /*
     * Example how to process the output.
     * Usually, one would define the legendname of the histogram to be used here
     * by hand, e.g. "signal" or "background".
     * To make this example run in any case, I am using alllegends.at(0), which
     * could e.g. be the signal legend.
     *
     * So in practise, the following would more look like
     * samplecoll.getHistos("signal");
     */
    //if(alllegends.size()>0){
    //    d_ana::histoCollection histos=samplecoll.getHistos(alllegends.at(0));

    //    /*
    //     * here, the histogram created in the analyze() function is selected and evaluated
    //     * The histoCollection maintains ownership (you don't need to delete the histogram)
    //     */
    //    const TH1* myplot=histos.getHisto("histoname1");

    //    std::cout << "(example output): the integral is " << myplot->Integral() <<std::endl;

    //    /*
    //     * If the histogram is subject to changes, please clone it and take ownership
    //     */

    //    TH1* myplot2=histos.cloneHisto("histoname1");

    //    /*
    //     * do something with the histogram
    //     */

    //    delete myplot2;
    //}

    /*
     * do the extraction here.
     */



}



