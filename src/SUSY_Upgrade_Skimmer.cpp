#include "interface/SUSY_Upgrade_Skimmer.h"

void SUSY_Upgrade_Skimmer::addBranches(){

    // Event variables
    myskim->Branch("genWeight", &genWeight);
    myskim->Branch("nTot", &nTot);
    myskim->Branch("xs", &xs);
    myskim->Branch("metSF", &metSF);

    // Electron vectors
    myskim->Branch("el_pt", &el_pt);
    myskim->Branch("el_eta", &el_eta);
    myskim->Branch("el_phi", &el_phi);
    myskim->Branch("el_q", &el_q);
    myskim->Branch("el_sumPt", &el_sumPt);
    myskim->Branch("el_woIso_pt", &el_woIso_pt);
    myskim->Branch("el_woIso_eta", &el_woIso_eta);
    myskim->Branch("el_woIso_phi", &el_woIso_phi);
    myskim->Branch("el_woIso_q", &el_woIso_q);
    myskim->Branch("el_woIso_sumPt", &el_woIso_sumPt);

    // Electron truth vectors
    myskim->Branch("el_pt_truth", &el_pt_truth);
    myskim->Branch("el_eta_truth", &el_eta_truth);
    myskim->Branch("el_phi_truth", &el_phi_truth);

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

    // Muon vectors
    myskim->Branch("mu_pt", &mu_pt);
    myskim->Branch("mu_eta", &mu_eta);
    myskim->Branch("mu_phi", &mu_phi);
    myskim->Branch("mu_q", &mu_q);
    myskim->Branch("mu_sumPt", &mu_sumPt);
    myskim->Branch("mu_matched", &mu_matched);
    myskim->Branch("mu_st20to30", &mu_st20to30);
    myskim->Branch("mu_mother", &mu_mother);
    myskim->Branch("mu_woIso_pt", &mu_woIso_pt);
    myskim->Branch("mu_woIso_eta", &mu_woIso_eta);
    myskim->Branch("mu_woIso_phi", &mu_woIso_phi);
    myskim->Branch("mu_woIso_q", &mu_woIso_q);
    myskim->Branch("mu_woIso_sumPt", &mu_woIso_sumPt);

    // Muon truth vectors
    myskim->Branch("mu_pt_truth", &mu_pt_truth);
    myskim->Branch("mu_eta_truth", &mu_eta_truth);
    myskim->Branch("mu_phi_truth", &mu_phi_truth);

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

    // Lepton vectors
    myskim->Branch("lep_pt", &lep_pt);
    myskim->Branch("lep_eta", &lep_eta);
    myskim->Branch("lep_phi", &lep_phi);
    myskim->Branch("lep_mass", &lep_mass);

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
    myskim->Branch("jet_pt", &jet_pt);
    myskim->Branch("jet_eta", &jet_eta);
    myskim->Branch("jet_phi", &jet_phi);
    myskim->Branch("jet_q", &jet_q);

    // Truth jet variables
    //myskim->Branch("jet1_pt_truth_matched", &jet1_pt_truth_matched);
    //myskim->Branch("jet1_eta_truth_matched", &jet1_eta_truth_matched);
    //myskim->Branch("jet1_phi_truth_matched", &jet1_phi_truth_matched);
    //myskim->Branch("jet1_q_truth_matched", &jet1_q_truth_matched);

    // MET variables
    myskim->Branch("met", &met);
    myskim->Branch("met_eta", &met_eta);
    myskim->Branch("met_phi", &met_phi);
    //myskim->Branch("mht", &mht);
    //myskim->Branch("mht_eta", &mht_eta);
    //myskim->Branch("mht_phi", &mht_phi);
    myskim->Branch("mht25", &mht25);
    myskim->Branch("mht40", &mht40);
    myskim->Branch("mht60", &mht60);
    myskim->Branch("mht100", &mht100);
    myskim->Branch("mht150", &mht150);
    myskim->Branch("mlt", &mlt);
    myskim->Branch("mlt_eta", &mlt_eta);
    myskim->Branch("mlt_phi", &mlt_phi);
    myskim->Branch("mhlt25", &mhlt25);
    myskim->Branch("mhlt25_eta", &mhlt25_eta);
    myskim->Branch("mhlt25_phi", &mhlt25_phi);
    myskim->Branch("mhlt40", &mhlt40);
    myskim->Branch("mhlt40_eta", &mhlt40_eta);
    myskim->Branch("mhlt40_phi", &mhlt40_phi);
    myskim->Branch("PFmet", &PFmet);
    myskim->Branch("PFmet_eta", &PFmet_eta);
    myskim->Branch("PFmet_phi", &PFmet_phi);
    myskim->Branch("genmet", &genmet);
    myskim->Branch("genmet_eta", &genmet_eta);
    myskim->Branch("genmet_phi", &genmet_phi);
    myskim->Branch("genpumet", &genpumet);
    myskim->Branch("genpumet_eta", &genpumet_eta);
    myskim->Branch("genpumet_phi", &genpumet_phi);

    // Other variables
    myskim->Branch("nLep", &nLep);
    //myskim->Branch("nLep_truth", &nLep_truth);
    myskim->Branch("nEl", &nEl);
    myskim->Branch("nMu", &nMu);
    myskim->Branch("nSoftLep", &nSoftLep);
    myskim->Branch("nSoftEl", &nSoftEl);
    myskim->Branch("nSoftMu", &nSoftMu);
    myskim->Branch("nBJet", &nBJet);
    myskim->Branch("nW", &nW);
    myskim->Branch("nZ", &nZ);
    myskim->Branch("nJet25", &nJet25);
    myskim->Branch("nJet40", &nJet40);
    myskim->Branch("nJet60", &nJet60);
    myskim->Branch("nJet100", &nJet100);
    myskim->Branch("nJet150", &nJet150);
    myskim->Branch("ht25", &ht25);
    myskim->Branch("ht40", &ht40);
    myskim->Branch("ht60", &ht60);
    myskim->Branch("ht100", &ht100);
    myskim->Branch("ht150", &ht150);
    myskim->Branch("genht25", &genht25);
    myskim->Branch("genht40", &genht40);
    myskim->Branch("hasSFOS", &hasSFOS);
    myskim->Branch("hasSFOS_truth", &hasSFOS_truth);
    myskim->Branch("hasSoftSFOS", &hasSoftSFOS);
    myskim->Branch("hasSoftSFOS_truth", &hasSoftSFOS_truth);
    myskim->Branch("mllMin", &mllMin);
    myskim->Branch("mllMax", &mllMax);
    myskim->Branch("mt1", &mt1);
    myskim->Branch("mt2", &mt2);
    myskim->Branch("pt2l", &pt2l);
    myskim->Branch("ZtoLL", &ZtoLL);
    myskim->Branch("crazyMuon50", &crazyMuon50);
    myskim->Branch("crazyMuon200", &crazyMuon200);
    myskim->Branch("crazyMuon500", &crazyMuon500);
    //myskim->Branch("mu_pt5to10_origin_nghbr", &mu_pt5to10_origin_nghbr);
    //myskim->Branch("mu_pt5to10_origin_cone", &mu_pt5to10_origin_cone);
    //myskim->Branch("mu_pt10to20_origin_nghbr", &mu_pt10to20_origin_nghbr);
    //myskim->Branch("mu_pt10to20_origin_cone", &mu_pt10to20_origin_cone);
    //myskim->Branch("mu_pt20to30_origin_nghbr", &mu_pt20to30_origin_nghbr);
    //myskim->Branch("mu_pt20to30_origin_cone", &mu_pt20to30_origin_cone);
}

void SUSY_Upgrade_Skimmer::clearVectors(){

    // Clear vectors
    el_pt.clear();
    el_eta.clear();
    el_phi.clear();
    el_q.clear();
    el_sumPt.clear();
    el_woIso_pt.clear();
    el_woIso_eta.clear();
    el_woIso_phi.clear();
    el_woIso_q.clear();
    el_woIso_sumPt.clear();
    el_pt_truth.clear();
    el_eta_truth.clear();
    el_phi_truth.clear();
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
    mu_pt.clear();
    mu_eta.clear();
    mu_phi.clear();
    mu_q.clear();
    mu_sumPt.clear();
    mu_matched.clear();
    mu_st20to30.clear();
    mu_mother.clear();
    mu_woIso_pt.clear();
    mu_woIso_eta.clear();
    mu_woIso_phi.clear();
    mu_woIso_q.clear();
    mu_woIso_sumPt.clear();
    mu_pt_truth.clear();
    mu_eta_truth.clear();
    mu_phi_truth.clear();
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
    lep_pt.clear();
    lep_eta.clear();
    lep_phi.clear();
    lep_mass.clear();
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
    jet_pt.clear();
    jet_eta.clear();
    jet_phi.clear();
    jet_q.clear();
    //jet1_pt_truth_matched.clear();
    //jet1_eta_truth_matched.clear();
    //jet1_phi_truth_matched.clear();
    //jet1_q_truth_matched.clear();
    mllMin.clear();
    mllMax.clear();
    mt1.clear();
    mt2.clear();
    pt2l.clear();
    //mu_pt5to10_origin_nghbr.clear();
    //mu_pt5to10_origin_cone.clear();
    //mu_pt10to20_origin_nghbr.clear();
    //mu_pt10to20_origin_cone.clear();
    //mu_pt20to30_origin_nghbr.clear();
    //mu_pt20to30_origin_cone.clear();
}

template <typename T> bool SUSY_Upgrade_Skimmer::isIsolated(const T particle){
    //// For all samples except ttbar, don't apply isolation
    //if (getSampleFile()(0, 6) != "tt-4p-"){
    //    return true;
    //}
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
    for (size_t i=0; i<elecs.size(); ++i){
        //if (!isIsolated(elecs.at(i))){ continue; }
        if (DeltaR(jet->Eta, elecs.at(i)->Eta, jet->Phi, elecs.at(i)->Phi) < jet_or_dr){ return true; }
    }
    for (size_t i=0; i<muons.size(); ++i){
        //if (!isIsolated(muons.at(i))){ continue; }
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
            if      (eta < 0.50){ passRandomEfficiency(0., ppt); }
            else if (eta < 1.00){ passRandomEfficiency(0., ppt); }
            else if (eta < 1.45){ passRandomEfficiency(0., ppt); }
            else if (eta < 1.55){ passRandomEfficiency(0., ppt); }
            else if (eta < 2.00){ passRandomEfficiency(0., ppt); }
            else if (eta < 2.50){ passRandomEfficiency(0., ppt); }
            else if (eta < 3.00){ passRandomEfficiency(0., ppt); }
        }else if (pt < 6.){
            if      (eta < 0.50){ passRandomEfficiency(.018/.50, ppt); }
            else if (eta < 1.00){ passRandomEfficiency(.016/.50, ppt); }
            else if (eta < 1.45){ passRandomEfficiency(.005/.50, ppt); }
            else if (eta < 1.55){ passRandomEfficiency(.026/.35, ppt); }
            else if (eta < 2.00){ passRandomEfficiency(.061/.75, ppt); }
            else if (eta < 2.50){ passRandomEfficiency(.100/.65, ppt); }
            else if (eta < 3.00){ passRandomEfficiency(.049/.65, ppt); }
        }else if (pt < 8.){
            if      (eta < 0.50){ passRandomEfficiency(.252/.50, ppt); }
            else if (eta < 1.00){ passRandomEfficiency(.198/.50, ppt); }
            else if (eta < 1.45){ passRandomEfficiency(.029/.50, ppt); }
            else if (eta < 1.55){ passRandomEfficiency(.045/.35, ppt); }
            else if (eta < 2.00){ passRandomEfficiency(.191/.75, ppt); }
            else if (eta < 2.50){ passRandomEfficiency(.223/.65, ppt); }
            else if (eta < 3.00){ passRandomEfficiency(.152/.65, ppt); }
        }else if (pt < 10.){
            if      (eta < 0.50){ passRandomEfficiency(.480/.50, ppt); }
            else if (eta < 1.00){ passRandomEfficiency(.446/.50, ppt); }
            else if (eta < 1.45){ passRandomEfficiency(.108/.50, ppt); }
            else if (eta < 1.55){ passRandomEfficiency(.133/.35, ppt); }
            else if (eta < 2.00){ passRandomEfficiency(.337/.75, ppt); }
            else if (eta < 2.50){ passRandomEfficiency(.427/.65, ppt); }
            else if (eta < 3.00){ passRandomEfficiency(.436/.65, ppt); }
        }else if (pt < 20.){
            if      (eta < 0.50){ passRandomEfficiency(.681/.94, ppt); }
            else if (eta < 1.00){ passRandomEfficiency(.598/.94, ppt); }
            else if (eta < 1.45){ passRandomEfficiency(.289/.94, ppt); }
            else if (eta < 1.55){ passRandomEfficiency(.411/.40, ppt); }
            else if (eta < 2.00){ passRandomEfficiency(.475/.85, ppt); }
            else if (eta < 2.50){ passRandomEfficiency(.590/.75, ppt); }
            else if (eta < 3.00){ passRandomEfficiency(.679/.75, ppt); }
        }else if (pt < 35.){
            if      (eta < 0.50){ passRandomEfficiency(.792/.94, ppt); }
            else if (eta < 1.00){ passRandomEfficiency(.759/.94, ppt); }
            else if (eta < 1.45){ passRandomEfficiency(.570/.94, ppt); }
            else if (eta < 1.55){ passRandomEfficiency(.629/.40, ppt); }
            else if (eta < 2.00){ passRandomEfficiency(.605/.85, ppt); }
            else if (eta < 2.50){ passRandomEfficiency(.720/.75, ppt); }
            else if (eta < 3.00){ passRandomEfficiency(.778/.75, ppt); }
        }else if (pt < 50.){
            if      (eta < 0.50){ passRandomEfficiency(.862/.97, ppt); }
            else if (eta < 1.00){ passRandomEfficiency(.847/.97, ppt); }
            else if (eta < 1.45){ passRandomEfficiency(.743/.97, ppt); }
            else if (eta < 1.55){ passRandomEfficiency(.761/.45, ppt); }
            else if (eta < 2.00){ passRandomEfficiency(.713/.95, ppt); }
            else if (eta < 2.50){ passRandomEfficiency(.800/.90, ppt); }
            else if (eta < 3.00){ passRandomEfficiency(.830/.90, ppt); }
        }else if (pt < 14000.){
            if      (eta < 0.50){ passRandomEfficiency(.859/.99, ppt); }
            else if (eta < 1.00){ passRandomEfficiency(.872/.99, ppt); }
            else if (eta < 1.45){ passRandomEfficiency(.828/.99, ppt); }
            else if (eta < 1.55){ passRandomEfficiency(.752/.50, ppt); }
            else if (eta < 2.00){ passRandomEfficiency(.794/.98, ppt); }
            else if (eta < 2.50){ passRandomEfficiency(.840/.90, ppt); }
            else if (eta < 3.00){ passRandomEfficiency(.919/.90, ppt); }
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

//double SUSY_Upgrade_Skimmer::coneVeto(double pt, double eta, double phi, d_ana::dBranchHandler<GenParticle>& genpart){
//
//    int nghbr = 99;
//    double drMin = 99.;
//    double drHfMin = 99.;
//    double drTauMin = 99.;
//
//    double iso = 0.;
//
//    const double cone = .9;
//
//    // Vectors with particles that have already been filled
//    // These are used to not fill the same particle twice
//    std::vector<int> filledPid;
//    std::vector<double> filledPt;
//    std::vector<double> filledEta;
//    std::vector<double> filledPhi;
//
//    for (size_t j=0; j<genpart.size(); ++j){
//
//        // Remove muon itself
//        if (isMatched(genpart.at(j), pt, eta, phi)){ continue; }
//
//        // Check if this particle has been filled before
//        // This needs to be checked since similar copies of truth
//        // particles are stored
//        bool skipParticle = false;
//        for (size_t k=0; k<filledPid.size(); ++k){
//            // PID needs to be the same and the particles need to match
//            if ((genpart.at(j)->PID == filledPid.at(k)) && isMatched(genpart.at(j), filledPt.at(k), filledEta.at(k), filledPhi.at(k))){
//                skipParticle = true;
//                break;
//            }
//        }
//        if (skipParticle){ continue; }
//
//        double dr = DeltaR(eta, genpart.at(j)->Eta, phi, genpart.at(j)->Phi);
//
//        // Check b's
//        if ((fabs(genpart.at(j)->PID) == 4 || fabs(genpart.at(j)->PID == 5)) && dr < drHfMin){
//            drHfMin = dr;
//        }
//
//        // Check taus
//        if (fabs(genpart.at(j)->PID) == 15 && dr < drTauMin){
//            drTauMin = dr;
//        }
//
//        if (dr < drMin){
//            drMin = dr;
//            nghbr = getNghbr(genpart.at(j)->PID);
//        }
//        if (dr < cone && pt > 5){
//            if (pt < 10){
//                mu_pt5to10_origin_cone.push_back(getNghbr(genpart.at(j)->PID));
//            }else if (pt < 20){
//                mu_pt10to20_origin_cone.push_back(getNghbr(genpart.at(j)->PID));
//            }else if (pt < 30){
//                mu_pt20to30_origin_cone.push_back(getNghbr(genpart.at(j)->PID));
//            }
//            filledPid.push_back(genpart.at(j)->PID);
//            filledPt.push_back(genpart.at(j)->PT);
//            filledEta.push_back(genpart.at(j)->Eta);
//            filledPhi.push_back(genpart.at(j)->Phi);
//            iso += genpart.at(j)->PT;
//        }
//    }
//
//    if (pt > 5){
//        if (pt < 10){
//            mu_pt5to10_origin_nghbr.push_back(nghbr);
//        }else if (pt < 20){
//            mu_pt10to20_origin_nghbr.push_back(nghbr);
//        }else if (pt < 30){
//            mu_pt20to30_origin_nghbr.push_back(nghbr);
//        }
//    }
//
//    // Weight to be returned
//    // Here I finally sell my scientific soul to the gods of publish or perish
//    //double wght = 1./(1. + wght_gen_iso*iso/pt);
//    // Or maybe not?
//    double wght = 1.;
//    if (drHfMin > cone && drTauMin > cone){
//        ;
//    }else if (drHfMin < drTauMin){
//        wght *= wght_hf_veto;
//    }else{
//        wght *= wght_tau_veto;
//    }
//
//    return wght;
//}

// Print properties of particle
template <typename T> void SUSY_Upgrade_Skimmer::pppWpidWstatus(const char* text, const size_t idx, const size_t noParticles, const T particle, const char* addText) const {
    ppp(text, idx, noParticles, particle, particle->PID, particle->Status, addText);
}

template <typename T> void SUSY_Upgrade_Skimmer::ppp(const char* text, const size_t idx, const size_t noParticles, const T particle, const int pid, const int status, const char* addText) const {
    printf("%20s: Idx: %3lu/%3lu; ID: %8d; Status: %3d; pt: %8.3f; eta: %6.3f; phi: %6.3f; %s\n",
            text, idx, noParticles, pid, status, particle->PT, particle->Eta, particle->Phi, addText);
    fflush(stdout);
}

void SUSY_Upgrade_Skimmer::analyze(size_t childid /* this info can be used for printouts */){

    d_ana::dBranchHandler<HepMCEvent> event(tree(),"Event");
    d_ana::dBranchHandler<Electron> elecs(tree(),"Electron");
    d_ana::dBranchHandler<Muon> muontight(tree(),"MuonTight");
    d_ana::dBranchHandler<Jet> jetpuppi(tree(), "JetPUPPI");
    d_ana::dBranchHandler<MissingET> puppimet(tree(), "PuppiMissingET");
    d_ana::dBranchHandler<GenParticle> genpart(tree(),"Particle");
    d_ana::dBranchHandler<Jet> genjet(tree(),"GenJet");
    d_ana::dBranchHandler<MissingET> PFmeth(tree(), "MissingET");
    d_ana::dBranchHandler<MissingET> genmeth(tree(), "GenMissingET");
    d_ana::dBranchHandler<MissingET> genpumeth(tree(), "GenPileUpMissingET");
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

        if (dump_genpart){
            std::cout << "\nNEW EVENT!" << std::endl;
            for (size_t i=0; i<genpart.size(); ++i){
                std::cout << "N: " << i
                    << ", St: " << genpart.at(i)->Status
                    << ", PID: " << genpart.at(i)->PID
                    << ", E: " << genpart.at(i)->E
                    << ", Px: " << genpart.at(i)->Px
                    << ", Py: " << genpart.at(i)->Py
                    << ", Pz: " << genpart.at(i)->Pz
                    << ", M: " << genpart.at(i)->Mass
                    << ", M1: " << genpart.at(i)->M1
                    << ", M2: " << genpart.at(i)->M2
                    << ", D1: " << genpart.at(i)->D1
                    << ", D2: " << genpart.at(i)->D2 << std::endl;
            }
            std::cout << std::endl;
        }

        if (event_by_event_comparison){
            bool evtFound = false;
            // The hardcoded numbers are to be used with the following sample:
            // root://cmsxrootd.fnal.gov//store/mc/PhaseIITDRFall17MiniAOD/DYJetsToLL_M-10to50_TuneCUETP8M1_14TeV-madgraphMLM-pythia8/MINIAODSIM/PU200_93X_upgrade2023_realistic_v2-v2/150000/02D196DB-60C0-E711-8C46-24BE05CE2D41.root
            for (size_t i=0; i<genpart.size(); ++i){
                if ((fabs(genpart.at(i)->PT -  9.343750) < 1.e-6 && fabs(genpart.at(i)->Eta - 3.948608) < 1.e-6) ||
                        (fabs(genpart.at(i)->PT - 11.757812) < 1.e-6 && fabs(genpart.at(i)->Eta - 0.839775) < 1.e-6) ||
                        (fabs(genpart.at(i)->PT -  7.476562) < 1.e-6 && fabs(genpart.at(i)->Eta - 2.270211) < 1.e-6) ||
                        (fabs(genpart.at(i)->PT -  7.230469) < 1.e-6 && fabs(genpart.at(i)->Eta + 0.100905) < 1.e-6) ||
                        (fabs(genpart.at(i)->PT - 10.335938) < 1.e-6 && fabs(genpart.at(i)->Eta - 3.225501) < 1.e-6)){

                    // Found FullSim truth particle! Event matched!
                    printf("Event by event comparison. Compare event %lld with weight %f:\n", event.at(0)->Number, event.at(0)->Weight);
                    pppWpidWstatus("Matched with FullSim", i, genpart.size(), genpart.at(i));

                    evtFound = true;
                    break;
                }
            }

            if (evtFound){

                // Truth objects
                printf("%s\n", "Truth objects");
                // Truth electrons
                for (size_t i=0; i<genpart.size(); ++i){
                    if (fabs(genpart.at(i)->PID) != 11){ continue; }
                    pppWpidWstatus("Truth electrons", i, genpart.size(), genpart.at(i));
                }
                // Truth muons
                for (size_t i=0; i<genpart.size(); ++i){
                    if (fabs(genpart.at(i)->PID) != 13){ continue; }
                    pppWpidWstatus("Truth muons", i, genpart.size(), genpart.at(i));
                }
                // Truth particles (all)
                for (size_t i=0; i<genpart.size(); ++i){
                    pppWpidWstatus("Truth particles", i, genpart.size(), genpart.at(i));
                }
                // Truth jets
                for (size_t i=0; i<genjet.size(); ++i){
                    ppp("Truth jets", i, genjet.size(), genjet.at(i), -1, -1);
                }
                // Truth MET
                printf("%s: Idx: %3d/%3d; ID: %8s; Status: %3s; pt: %8.3f; eta: %6.3f; phi: %6.3f\n",
                        "Truth MET", 1, 1, "-", "-", genmeth.at(0)->MET, genmeth.at(0)->Eta, genmeth.at(0)->Phi);

                // Reco objects
                printf("%s\n", "Reco objects");
                // Reco electrons
                for (size_t i=0; i<elecs.size(); ++i){
                    // Additionally print SumPt
                    char SumPt[20];
                    snprintf(SumPt, sizeof SumPt, "%f", elecs.at(i)->SumPt);
                    char addText[60] = "sumPt: ";
                    strcat(addText, SumPt);
                    ppp("Reco electrons", i, elecs.size(), elecs.at(i), elecs.at(i)->Charge>0 ? -11: 11, 1, addText);
                }
                // Reco muons
                for (size_t i=0; i<muontight.size(); ++i){
                    // Additionally print SumPt
                    char SumPt[20];
                    snprintf(SumPt, sizeof SumPt, "%f", muontight.at(i)->SumPt);
                    char addText[60] = "sumPt: ";
                    strcat(addText, SumPt);
                    ppp("Reco muons", i, muontight.size(), muontight.at(i), muontight.at(i)->Charge>0 ? -13: 13, 1, addText);
                }
                // Reco jets
                for (size_t i=0; i<jetpuppi.size(); ++i){
                    ppp("Reco jets", i, jetpuppi.size(), jetpuppi.at(i), -1, -1);
                }
                // Reco MET
                printf("%20s: Idx: %3d/%3d; ID: %8s; Status: %3s; pt: %8.3f; eta: %6.3f; phi: %6.3f\n",
                        "Reco MET", 1, 1, "-", "-", puppimet.at(0)->MET, puppimet.at(0)->Eta, puppimet.at(0)->Phi);
            }
        }

        // Lepton on-the-fly efficiencies
        effOnTopElec(elecs);
        effOnTopMuon(muontight);

        // Cutflow variables
        nLep = nEl = nMu = 0;
        nSoftLep = nSoftEl = nSoftMu = 0;
        nBJet = 0;
        nW = nZ = 0;
        nJet25 = nJet40 = nJet60 = nJet100 = nJet150 = 0;
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
                if (jetpuppi.at(i)->BTag){ nBJet++; }
            }
            if (jetpuppi.at(i)->PT > 25.){ nJet25++; }
            if (jetpuppi.at(i)->PT > 40.){ nJet40++; }
            if (jetpuppi.at(i)->PT > 60.){ nJet60++; }
            if (jetpuppi.at(i)->PT > 100.){ nJet100++; }
            if (jetpuppi.at(i)->PT > 150.){ nJet150++; }
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

                    // Check if both particles are soft
                    if (elecs.at(i)->PT < el_pt_hi && elecs.at(j)->PT < el_pt_hi){
                        hasSoftSFOS = true;
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
        //if (nLep < 2){ continue; }
        //if (nSoftLep < 2){ continue; }
        //if (!hasSoftSFOS){ continue; }

        // Fill electrons
        for (size_t i=0; i<elecs.size(); ++i){
            if (elecs.at(i)->PT < el_pt_lo){ continue; }

            el_woIso_pt.push_back(elecs.at(i)->PT);
            el_woIso_eta.push_back(elecs.at(i)->Eta);
            el_woIso_phi.push_back(elecs.at(i)->Phi);
            el_woIso_q.push_back(elecs.at(i)->Charge);
            el_woIso_sumPt.push_back(elecs.at(i)->SumPt);

            if (!isIsolated(elecs.at(i))){ continue; }

            el_pt.push_back(elecs.at(i)->PT);
            el_eta.push_back(elecs.at(i)->Eta);
            el_phi.push_back(elecs.at(i)->Phi);
            el_q.push_back(elecs.at(i)->Charge);
            el_sumPt.push_back(elecs.at(i)->SumPt);
        }

        // Fill truth electrons
        for (size_t i=0; i<genpart.size(); ++i){
            if (genpart.at(i)->Status != 1){ continue; }
            if (genpart.at(i)->PID != 11){ continue; }
            if (genpart.at(i)->PT < el_pt_lo){ continue; }
            el_pt_truth.push_back(genpart.at(i)->PT);
            el_eta_truth.push_back(genpart.at(i)->Eta);
            el_phi_truth.push_back(genpart.at(i)->Phi);
        }

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

            if (logdebug){
                fprintf(stderr, "Inspect reco muon (pT: %f; eta: %f; phi: %f) [Muon #%lu of %lu].\n",
                        muontight.at(i)->PT, muontight.at(i)->Eta, muontight.at(i)->Phi, i+1, muontight.size());
            }

            if (muontight.at(i)->PT < mu_pt_lo){ continue; }

            // Fill muon vector ignoring isolation
            mu_woIso_pt.push_back(muontight.at(i)->PT);
            mu_woIso_eta.push_back(muontight.at(i)->Eta);
            mu_woIso_phi.push_back(muontight.at(i)->Phi);
            mu_woIso_q.push_back(muontight.at(i)->Charge);
            mu_woIso_sumPt.push_back(muontight.at(i)->SumPt);

            if (!isIsolated(muontight.at(i))){ continue; }

            // Default variables to eventually be changed and filled into vectors
            bool match = false;
            bool st20to30 = false;
            int mother = -999;

            // vector to store *index* of truth matched final state muon and all its ancestors
            std::vector<unsigned int> fsMu;

            // If readonly is true, no variables that are written to the output n-tuple are going to be changed;
            // this is helpful for better understanding an event, without changing the output;
            // for example, a b-quark is a final particle and what happened before that will not change the output,
            // but it can still be helpful to further inspect the event
            bool readonly = false;

            // Absolute values of PDGID's that are considered as final ancestors from muons;
            // in other words: once one of these ancestors is found, no more ancestors are checked
            const std::vector<int> mu_mother_final = {4, 5, 2212};

            // Check if you can match the muon
            for (size_t j=0; j<genpart.size(); ++j){
                if (genpart.at(j)->Status != 1){ continue; }
                if (fabs(genpart.at(j)->PID) != 13){ continue; }
                // Truth matching
                if (isMatched(genpart.at(j), muontight.at(i)->PT, muontight.at(i)->Eta, muontight.at(i)->Phi)){
                    if (!readonly){ match = true; }
                    fsMu.push_back(j);

                    if (logdebug){
                        fprintf(stderr, "Matched with truth muon (pT: %f; eta: %f; phi: %f)\n",
                                genpart.at(j)->PT, genpart.at(j)->Eta, genpart.at(j)->Phi);
                    }

                    // Loop through vector until it's empty
                    while (fsMu.size() > 0){

                        if (logdebug){
                            fprintf(stderr, "Content of muon vector: ");
                            for (size_t k=0; k<fsMu.size(); ++k){
                                fprintf(stderr, "%u; ", fsMu[k]);
                            }
                            fprintf(stderr, "now inspecting first element: %u (PDGID: %d; Status: %d).\n",
                                    fsMu[0], genpart.at(fsMu[0])->PID, genpart.at(fsMu[0])->Status);
                            //fprintf(stderr, "Ancestor of muon: %d.\n", fsMu[0]);
                        }

                        // Check if it is a muon with status between 20 and 30
                        if (fabs(genpart.at(fsMu[0])->PID) == 13 && genpart.at(fsMu[0])->Status >= 20 && genpart.at(fsMu[0])->Status <= 30){
                            if (!readonly){ st20to30 = true; }

                            if (logdebug){
                                fprintf(stderr, "Muon with status between 20 and 30 found.\n");
                            }
                        }

                        // Check if ancestor of muon is final
                        if (std::find(mu_mother_final.begin(), mu_mother_final.end(), fabs(genpart.at(fsMu[0])->PID)) != std::end(mu_mother_final)){
                            if (!readonly){ mother = genpart.at(fsMu[0])->PID; }

                            if (logdebug){
                                fprintf(stderr, "Found final mother of muon: %d.\n", genpart.at(fsMu[0])->PID);
                            }

                            if (logdebug){
                                if (!readonly){
                                    fprintf(stderr, "Enable read-only mode.\n");
                                    readonly = true;
                                }
                            }else{
                                break;
                            }
                        }

                        // Add mothers from particle to vector (and inspect them
                        // in the next iteration of the loop
                        int m1 = genpart.at(fsMu[0])->M1;
                        int m2 = genpart.at(fsMu[0])->M2;
                        if (m1 >= 0. && (std::find(fsMu.begin(), fsMu.end(), m1)) == std::end(fsMu)){
                            fsMu.push_back(m1);
                        }
                        if (m2 >= 0. && (std::find(fsMu.begin(), fsMu.end(), m2)) == std::end(fsMu)){
                            fsMu.push_back(m2);
                        }
                        fsMu.erase(fsMu.begin());
                    }

                    // break free from genpart loop, since reco particle has been matched
                    break;
                }
            }

            if (logdebug){
                fprintf(stderr, "Write to output: match: %d; st20to30: %d; mother: %d\n", match, st20to30, mother);
            }

            // Fill isolated muons
            mu_pt.push_back(muontight.at(i)->PT);
            mu_eta.push_back(muontight.at(i)->Eta);
            mu_phi.push_back(muontight.at(i)->Phi);
            mu_q.push_back(muontight.at(i)->Charge);
            mu_sumPt.push_back(muontight.at(i)->SumPt);
            mu_matched.push_back(match);
            mu_st20to30.push_back(st20to30);
            mu_mother.push_back(mother);
        }

        // Fill truth muons
        for (size_t i=0; i<genpart.size(); ++i){
            if (genpart.at(i)->Status != 1){ continue; }
            if (genpart.at(i)->PID != 13){ continue; }
            if (genpart.at(i)->PT < mu_pt_lo){ continue; }
            mu_pt_truth.push_back(genpart.at(i)->PT);
            mu_eta_truth.push_back(genpart.at(i)->Eta);
            mu_phi_truth.push_back(genpart.at(i)->Phi);
        }

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

        //// Guess origin of leptons
        //for (size_t i=0; i<mu1_pt.size(); ++i){
        //    genWeight *= coneVeto(mu1_pt.at(i), mu1_eta.at(i), mu1_phi.at(i), genpart);
        //}
        //for (size_t i=0; i<mu2_pt.size(); ++i){
        //    genWeight *= coneVeto(mu2_pt.at(i), mu2_eta.at(i), mu2_phi.at(i), genpart);
        //}

        // Fill leptons
        // Put pT and eta into vector of vector for simultaneous sorting
        std::vector<std::vector<double>> lepvec;
        for (size_t i=0; i<el_pt.size(); ++i){
            lepvec.push_back({el_pt.at(i), el_eta.at(i), el_phi.at(i), mass_el});
        }
        for (size_t i=0; i<mu_pt.size(); ++i){
            lepvec.push_back({mu_pt.at(i), mu_eta.at(i), mu_phi.at(i), mass_mu});
        }

        // By definition, this sorts by the first element of the vector (in this case pT)
        if (lepvec.size() > 1){
            std::sort(begin(lepvec), end(lepvec));
            std::reverse(begin(lepvec), end(lepvec));
        }

        // Fill specific vectors
        for (size_t i=0; i<lepvec.size(); ++i){
            lep_pt.push_back(lepvec[i][0]);
            lep_eta.push_back(lepvec[i][1]);
            lep_phi.push_back(lepvec[i][2]);
            lep_mass.push_back(lepvec[i][3]);
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
            jet_pt.push_back(jetpuppi.at(i)->PT);
            jet_eta.push_back(jetpuppi.at(i)->Eta);
            jet_phi.push_back(jetpuppi.at(i)->Phi);
            jet_q.push_back(jetpuppi.at(i)->Charge);
        }

        //// Fill truth jets
        //for (size_t i=0; i<jet1_pt.size(); ++i){
        //    for (size_t j=0; j<genjet.size(); ++j){
        //        // Truth matching
        //        if (fabs(genjet.at(j)->PT - jet1_pt.at(i)) > truth_match_diff_pt){ continue; }
        //        if (fabs(genjet.at(j)->Eta - jet1_eta.at(i) > truth_match_diff_eta)){ continue; }
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
        try{
            PFmet = PFmeth.at(0)->MET;
            PFmet_eta = PFmeth.at(0)->Eta;
            PFmet_phi = PFmeth.at(0)->Phi;
        }catch (const std::out_of_range& oor){
            std::cerr << "Out of range error when accessing PF MET vector: " << oor.what() << std::endl;
            return;
        }
        try{
            genmet = genmeth.at(0)->MET;
            genmet_eta = genmeth.at(0)->Eta;
            genmet_phi = genmeth.at(0)->Phi;
        }catch (const std::out_of_range& oor){
            std::cerr << "Out of range error when accessing Gen MET vector: " << oor.what() << std::endl;
            return;
        }
        try{
            genpumet = genpumeth.at(0)->MET;
            genpumet_eta = genpumeth.at(0)->Eta;
            genpumet_phi = genpumeth.at(0)->Phi;
        }catch (const std::out_of_range& oor){
            std::cerr << "Out of range error when accessing Gen PU MET vector: " << oor.what() << std::endl;
            return;
        }

        //// Fill poor man's MET
        TLorentzVector mlt4;
        TLorentzVector mht4v25, mht4v40, mht4v60, mht4v100, mht4v150;
        TLorentzVector mhlt4v25, mhlt4v40;
        for (size_t i=0; i<muontight.size(); ++i){
            if (muontight.at(i)->PT < mu_pt_lo || !isIsolated(muontight.at(i))){ continue; }
            TLorentzVector m4;
            m4.SetPtEtaPhiM(muontight.at(i)->PT, muontight.at(i)->Eta, muontight.at(i)->Phi, mass_mu);
            mlt4 += m4;
            mhlt4v25 += m4;
            mhlt4v40 += m4;
        }
        for (size_t i=0; i<elecs.size(); ++i){
            if (elecs.at(i)->PT < el_pt_lo || !isIsolated(elecs.at(i))){ continue; }
            TLorentzVector e4;
            e4.SetPtEtaPhiM(elecs.at(i)->PT, elecs.at(i)->Eta, elecs.at(i)->Phi, mass_el);
            mlt4 += e4;
            mhlt4v25 += e4;
            mhlt4v40 += e4;
        }
        for (size_t i=0; i<jetpuppi.size(); ++i){
            if (jetpuppi.at(i)->PT < jet_pt_lo){ continue; }
            if (isOverlap(jetpuppi.at(i), elecs, muontight)){ continue; }
            TLorentzVector j4;
            j4.SetPtEtaPhiM(jetpuppi.at(i)->PT, jetpuppi.at(i)->Eta, jetpuppi.at(i)->Phi, jetpuppi.at(i)->Mass);
            if (jetpuppi.at(i)->PT > 25.){ mht4v25 += j4; mhlt4v25 += j4; }
            if (jetpuppi.at(i)->PT > 40.){ mht4v40 += j4; mhlt4v40 += j4; }
            if (jetpuppi.at(i)->PT > 60.){ mht4v60 += j4; }
            if (jetpuppi.at(i)->PT > 100.){ mht4v100 += j4; }
            if (jetpuppi.at(i)->PT > 150.){ mht4v150 += j4; }
        }
        mlt = mlt4.Pt();
        mlt_eta = mlt4.Eta();
        mlt_phi = mlt4.Phi();
        mht25 = mht4v25.Pt();
        mht40 = mht4v40.Pt();
        mht60 = mht4v60.Pt();
        mht100 = mht4v100.Pt();
        mht150 = mht4v150.Pt();
        mhlt25 = mhlt4v25.Pt();
        mhlt25_eta = mhlt4v25.Eta();
        mhlt25_phi = mhlt4v25.Phi();
        mhlt40 = mhlt4v40.Pt();
        mhlt40_eta = mhlt4v40.Eta();
        mhlt40_phi = mhlt4v40.Phi();

        // Fill HT
        ht25 = ht40 = ht60 = ht100 = ht150 = 0.;
        for (size_t i=0; i<jetpuppi.size(); ++i){
            if (isOverlap(jetpuppi.at(i), elecs, muontight)){ continue; }
            if (jetpuppi.at(i)->PT > 25.){ ht25 += jetpuppi.at(i)->PT; }
            if (jetpuppi.at(i)->PT > 40.){ ht40 += jetpuppi.at(i)->PT; }
            if (jetpuppi.at(i)->PT > 60.){ ht60 += jetpuppi.at(i)->PT; }
            if (jetpuppi.at(i)->PT > 100.){ ht100 += jetpuppi.at(i)->PT; }
            if (jetpuppi.at(i)->PT > 150.){ ht150 += jetpuppi.at(i)->PT; }
        }

        // Fill GenHT
        genht25 = genht40 = 0.;
        for (size_t i=0; i<genjet.size(); ++i){
            if (genjet.at(i)->PT > 25.){
                genht25 += genjet.at(i)->PT;
            }
            if (genjet.at(i)->PT > 40.){
                genht40 += genjet.at(i)->PT;
            }
        }

        // In DYtoLL events, figure out what LL is
        unsigned int n11 = 0;
        unsigned int n13 = 0;
        unsigned int n15 = 0;
        for (size_t i=0; i<genpart.size(); ++i){
            // Only consider particles with a pT of at least 2 GeV
            if (genpart.at(i)->PID < 2){ continue; }
            // Count leptons (disregard status)
            if (fabs(genpart.at(i)->PID) == 11){
                n11++;
            }else if (fabs(genpart.at(i)->PID) == 13){
                n13++;
            }else if (fabs(genpart.at(i)->PID) == 15){
                n15++;
            }
        }
        // If there's at least one tau, it's a Z --> tautau event
        // Else, if there are more electrons than muons, it's a Z --> ee event
        // Else, if there are more muons than electrons, it's a Z --> mm event
        // Else, I don't know
        if (n15 >= 1){
            ZtoLL = 15;
        }else if (n11 > n13){
            ZtoLL = 11;
        }else if (n13 > n11){
            ZtoLL = 13;
        }else{
            ZtoLL = 999;
        }

        // Crazy muon filters (not used, for compatibility with FS)
        crazyMuon50 = crazyMuon200 = crazyMuon500 = false;

        // MET HT scale factors
        {
            TFile* fSF = new TFile("sf/met_ht_fs.root");
            TH1D* hSF = nullptr;
            if (getSampleFile()(0, 6) == "tt-4p-"){
                hSF = (TH1D*)fSF->Get("sf_met_ht200_coarse_varbin_tt_ratio");
            }else if (getSampleFile()(0, 10) == "DYJetsToLL"){
                hSF = (TH1D*)fSF->Get("sf_met_ht200_coarse_varbin_DY2_ratio");
            }else if (getSampleFile()(0, 10) == "WJetsToLNu"){
                hSF = (TH1D*)fSF->Get("sf_met_ht200_coarse_varbin_Wj2_ratio");
            }else if (getSampleFile()(0, 9) == "TChiWZOff"){
                hSF = (TH1D*)fSF->Get("sf_met_ht200_coarse_varbin_TChiWZ_400_375_Delphes_v09_ratio");
            }
            metSF = 1.;
            if (hSF){
                metSF = hSF->GetBinContent(hSF->FindBin(met));
            }
            fSF->Close();
        }

        // Transverse mass of leading two leptons
        if (lep_pt.size() >= 1){
            mt1.push_back(std::sqrt(2*lep_pt.at(0)*met*(1-std::cos(lep_phi.at(0)-met_phi))));
        }
        if (lep_pt.size() >= 2){
            mt2.push_back(std::sqrt(2*lep_pt.at(1)*met*(1-std::cos(lep_phi.at(1)-met_phi))));
            // Fill pt of two leptons
            TLorentzVector l1, l2;
            l1.SetPtEtaPhiM(lep_pt[0], lep_eta[0], lep_phi[0], lep_mass[0]);
            l2.SetPtEtaPhiM(lep_pt[1], lep_eta[1], lep_phi[1], lep_mass[1]);
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



