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

    // Truth electron variables
    myskim->Branch("el1_pt_truth", &el1_pt_truth);
    myskim->Branch("el1_eta_truth", &el1_eta_truth);
    myskim->Branch("el1_phi_truth", &el1_phi_truth);
    myskim->Branch("el1_q_truth", &el1_q_truth);
    myskim->Branch("el2_pt_truth", &el2_pt_truth);
    myskim->Branch("el2_eta_truth", &el2_eta_truth);
    myskim->Branch("el2_phi_truth", &el2_phi_truth);
    myskim->Branch("el2_q_truth", &el2_q_truth);

    // Muon variables
    myskim->Branch("mu1_tight_pt", &mu1_tight_pt);
    myskim->Branch("mu1_tight_eta", &mu1_tight_eta);
    myskim->Branch("mu1_tight_phi", &mu1_tight_phi);
    myskim->Branch("mu1_tight_q", &mu1_tight_q);
    myskim->Branch("mu2_tight_pt", &mu2_tight_pt);
    myskim->Branch("mu2_tight_eta", &mu2_tight_eta);
    myskim->Branch("mu2_tight_phi", &mu2_tight_phi);
    myskim->Branch("mu2_tight_q", &mu2_tight_q);

    // Truth muon variables
    myskim->Branch("mu1_pt_truth", &mu1_pt_truth);
    myskim->Branch("mu1_eta_truth", &mu1_eta_truth);
    myskim->Branch("mu1_phi_truth", &mu1_phi_truth);
    myskim->Branch("mu1_q_truth", &mu1_q_truth);
    myskim->Branch("mu2_pt_truth", &mu2_pt_truth);
    myskim->Branch("mu2_eta_truth", &mu2_eta_truth);
    myskim->Branch("mu2_phi_truth", &mu2_phi_truth);
    myskim->Branch("mu2_q_truth", &mu2_q_truth);

    // Lepton variables
    myskim->Branch("lep1_pt", &lep1_pt);
    myskim->Branch("lep1_eta", &lep1_eta);
    myskim->Branch("lep1_phi", &lep1_phi);
    myskim->Branch("lep1_mass", &lep1_mass);
    myskim->Branch("lep2_pt", &lep2_pt);
    myskim->Branch("lep2_eta", &lep2_eta);
    myskim->Branch("lep2_phi", &lep2_phi);
    myskim->Branch("lep2_mass", &lep2_mass);

    // Truth lepton variables
    myskim->Branch("lep1_pt_truth", &lep1_pt_truth);
    myskim->Branch("lep1_eta_truth", &lep1_eta_truth);
    myskim->Branch("lep1_phi_truth", &lep1_phi_truth);
    myskim->Branch("lep1_mass_truth", &lep1_mass_truth);
    myskim->Branch("lep2_pt_truth", &lep2_pt_truth);
    myskim->Branch("lep2_eta_truth", &lep2_eta_truth);
    myskim->Branch("lep2_phi_truth", &lep2_phi_truth);
    myskim->Branch("lep2_mass_truth", &lep2_mass_truth);

    // Jet variables
    myskim->Branch("jet1_puppi_pt", &jet1_puppi_pt);
    myskim->Branch("jet1_puppi_eta", &jet1_puppi_eta);
    myskim->Branch("jet1_puppi_phi", &jet1_puppi_phi);
    myskim->Branch("jet1_puppi_q", &jet1_puppi_q);

    // Truth jet variables
    myskim->Branch("jet1_pt_truth", &jet1_pt_truth);
    myskim->Branch("jet1_eta_truth", &jet1_eta_truth);
    myskim->Branch("jet1_phi_truth", &jet1_phi_truth);
    myskim->Branch("jet1_q_truth", &jet1_q_truth);

    // MET variables
    myskim->Branch("met", &met);
    myskim->Branch("met_eta", &met_eta);
    myskim->Branch("met_phi", &met_phi);

    // Other variables
    myskim->Branch("nLep", &nLep);
    myskim->Branch("nEl", &nEl);
    myskim->Branch("nMu", &nMu);
    myskim->Branch("nSoftLep", &nSoftLep);
    myskim->Branch("nSoftEl", &nSoftEl);
    myskim->Branch("nSoftMu", &nSoftMu);
    myskim->Branch("nJet", &nJet);
    myskim->Branch("nBJet", &nBJet);
    myskim->Branch("ht", &ht);
    myskim->Branch("hasSFOS", &hasSFOS);
    myskim->Branch("hasSoftSFOS", &hasSoftSFOS);
    myskim->Branch("mll", &mll);
    myskim->Branch("mt1", &mt1);
    myskim->Branch("mt2", &mt2);
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
    el1_pt_truth.clear();
    el1_eta_truth.clear();
    el1_phi_truth.clear();
    el1_q_truth.clear();
    el2_pt_truth.clear();
    el2_eta_truth.clear();
    el2_phi_truth.clear();
    el2_q_truth.clear();
    mu1_tight_pt.clear();
    mu1_tight_eta.clear();
    mu1_tight_phi.clear();
    mu1_tight_q.clear();
    mu2_tight_pt.clear();
    mu2_tight_eta.clear();
    mu2_tight_phi.clear();
    mu2_tight_q.clear();
    mu1_pt_truth.clear();
    mu1_eta_truth.clear();
    mu1_phi_truth.clear();
    mu1_q_truth.clear();
    mu2_pt_truth.clear();
    mu2_eta_truth.clear();
    mu2_phi_truth.clear();
    mu2_q_truth.clear();
    lep1_pt.clear();
    lep1_eta.clear();
    lep1_phi.clear();
    lep1_mass.clear();
    lep2_pt.clear();
    lep2_eta.clear();
    lep2_phi.clear();
    lep2_mass.clear();
    lep1_pt_truth.clear();
    lep1_eta_truth.clear();
    lep1_phi_truth.clear();
    lep1_mass_truth.clear();
    lep2_pt_truth.clear();
    lep2_eta_truth.clear();
    lep2_phi_truth.clear();
    lep2_mass_truth.clear();
    jet1_puppi_pt.clear();
    jet1_puppi_eta.clear();
    jet1_puppi_phi.clear();
    jet1_puppi_q.clear();
    jet1_pt_truth.clear();
    jet1_eta_truth.clear();
    jet1_phi_truth.clear();
    jet1_q_truth.clear();
    mll.clear();
    mt1.clear();
    mt2.clear();
}

template <typename T> bool SUSY_Upgrade_Skimmer::isIsolated(T particle){
    if (particle->IsolationVarRhoCorr/particle->PT < iso_cut){
        return true;
    }
    return false;
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

    myskim=addTree();
    addBranches();

    // Event loop
    size_t nevents=tree()->entries();
    if(isTestMode()){
        nevents/=1000;
    }
    for(size_t eventno=0;eventno<nevents;eventno++){

        // Report status and set event link
        reportStatus(eventno,nevents);
        tree()->setEntry(eventno);

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

        // Cutflow variables
        nLep = nEl = nMu = 0;
        nSoftLep = nSoftEl = nSoftMu = 0;
        nJet = nBJet = 0;
        for (size_t i=0; i<elecs.size(); ++i){
            if (elecs.at(i)->PT > el_pt_lo && isIsolated(elecs.at(i))){
                nLep++;
                nEl++;
                if (elecs.at(i)->PT < el_pt_hi){
                    nSoftLep++;
                    nSoftEl++;
                }
            }
        }
        for (size_t i=0; i<muontight.size(); ++i){
            if (muontight.at(i)->PT > mu_pt_lo && isIsolated(muontight.at(i))){
                nLep++;
                nMu++;
                if (muontight.at(i)->PT < mu_pt_hi){
                    nSoftLep++;
                    nSoftMu++;
                }
            }
        }
        for (size_t i=0; i<jetpuppi.size(); ++i){
            if (jetpuppi.at(i)->PT > jet_pt_lo){
                nJet++;
                if (jetpuppi.at(i)->BTag){
                    nBJet++;
                }
            }
        }
        //histo->Fill(nJet);
        //histo2->Fill(nJetN);
        //histo3->Fill(jetpuppi.size());
        //histo4->Fill(jetnotpuppi.size());
        //continue;

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

                        // Mll for soft SFOS (only first found pair considered)
                        TLorentzVector l1, l2;
                        l1.SetPtEtaPhiM(elecs.at(i)->PT, elecs.at(i)->Eta, elecs.at(i)->Phi, mass_el);
                        l2.SetPtEtaPhiM(elecs.at(j)->PT, elecs.at(j)->Eta, elecs.at(j)->Phi, mass_el);
                        mll.push_back((l1+l2).M());
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
                        mll.push_back((l1+l2).M());
                    }
                }
            }
        }

        // Skim
        if (nLep < 2){ continue; }
        if (nSoftLep < 2){ continue; }
        if (!hasSoftSFOS){ continue; }

        // Fill electrons
        if (nEl >= 1){
            el1_pt.push_back(elecs.at(0)->PT);
            el1_eta.push_back(elecs.at(0)->Eta);
            el1_phi.push_back(elecs.at(0)->Phi);
            el1_q.push_back(elecs.at(0)->Charge);
        }
        if (nEl >= 2){
            el2_pt.push_back(elecs.at(1)->PT);
            el2_eta.push_back(elecs.at(1)->Eta);
            el2_phi.push_back(elecs.at(1)->Phi);
            el2_q.push_back(elecs.at(1)->Charge);
        }

        // Fill truth electrons
        for (size_t i=0; i<genpart.size(); ++i){
            if (fabs(genpart.at(i)->PID) != 11){ continue; }
            // Fill it if it hasn't been filled
            if (el1_pt_truth.size() == 0){
                el1_pt_truth.push_back(genpart.at(i)->PT);
                el1_eta_truth.push_back(genpart.at(i)->Eta);
                el1_phi_truth.push_back(genpart.at(i)->Phi);
                el1_q_truth.push_back(genpart.at(i)->Charge);
            }else{
                el2_pt_truth.push_back(genpart.at(i)->PT);
                el2_eta_truth.push_back(genpart.at(i)->Eta);
                el2_phi_truth.push_back(genpart.at(i)->Phi);
                el2_q_truth.push_back(genpart.at(i)->Charge);
                // When second particle has been filled, break
                break;
            }
        }

        // Fill muons
        if (nMu >= 1){
            mu1_tight_pt.push_back(muontight.at(0)->PT);
            mu1_tight_eta.push_back(muontight.at(0)->Eta);
            mu1_tight_phi.push_back(muontight.at(0)->Phi);
            mu1_tight_q.push_back(muontight.at(0)->Charge);
        }
        if (nMu >= 2){
            mu2_tight_pt.push_back(muontight.at(1)->PT);
            mu2_tight_eta.push_back(muontight.at(1)->Eta);
            mu2_tight_phi.push_back(muontight.at(1)->Phi);
            mu2_tight_q.push_back(muontight.at(1)->Charge);
        }

        // Fill truth muons
        for (size_t i=0; i<genpart.size(); ++i){
            if (fabs(genpart.at(i)->PID) != 13){ continue; }
            // Fill it if it hasn't been filled
            if (mu1_pt_truth.size() == 0){
                mu1_pt_truth.push_back(genpart.at(i)->PT);
                mu1_eta_truth.push_back(genpart.at(i)->Eta);
                mu1_phi_truth.push_back(genpart.at(i)->Phi);
                mu1_q_truth.push_back(genpart.at(i)->Charge);
            }else{
                mu2_pt_truth.push_back(genpart.at(i)->PT);
                mu2_eta_truth.push_back(genpart.at(i)->Eta);
                mu2_phi_truth.push_back(genpart.at(i)->Phi);
                mu2_q_truth.push_back(genpart.at(i)->Charge);
                // When second particle has been filled, break
                break;
            }
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
        if (mu1_tight_pt.size() != 0){
            lepvec.push_back({mu1_tight_pt.at(0), mu1_tight_eta.at(0), mu1_tight_phi.at(0), mass_mu});
        }
        if (mu2_tight_pt.size() != 0){
            lepvec.push_back({mu2_tight_pt.at(0), mu2_tight_eta.at(0), mu2_tight_phi.at(0), mass_mu});
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

        // Fill truth leptons
        // Put pT and eta into vector of vector for sorting
        std::vector<std::vector<double>> lepvec_truth;
        if (el1_pt_truth.size() != 0){
            lepvec_truth.push_back({el1_pt_truth.at(0), el1_eta_truth.at(0), el1_phi_truth.at(0), mass_el});
        }
        if (el2_pt_truth.size() != 0){
            lepvec_truth.push_back({el2_pt_truth.at(0), el2_eta_truth.at(0), el2_phi_truth.at(0), mass_el});
        }
        if (mu1_pt_truth.size() != 0){
            lepvec_truth.push_back({mu1_pt_truth.at(0), mu1_eta_truth.at(0), mu1_phi_truth.at(0), mass_mu});
        }
        if (mu2_pt_truth.size() != 0){
            lepvec_truth.push_back({mu2_pt_truth.at(0), mu2_eta_truth.at(0), mu2_phi_truth.at(0), mass_mu});
        }
        // By definition, this sorts by the first element of the vector (in this case pT)
        if (lepvec_truth.size() > 1){
            std::sort(begin(lepvec_truth), end(lepvec_truth));
            std::reverse(begin(lepvec_truth), end(lepvec_truth));
        }
        if (lepvec_truth.size() >= 1){
            lep1_pt_truth.push_back(lepvec_truth[0][0]);
            lep1_eta_truth.push_back(lepvec_truth[0][1]);
            lep1_phi_truth.push_back(lepvec_truth[0][2]);
            lep1_mass_truth.push_back(lepvec_truth[0][3]);
        }
        if (lepvec_truth.size() >= 2){
            lep2_pt_truth.push_back(lepvec_truth[1][0]);
            lep2_eta_truth.push_back(lepvec_truth[1][1]);
            lep2_phi_truth.push_back(lepvec_truth[1][2]);
            lep2_mass_truth.push_back(lepvec_truth[1][3]);
        }
        lepvec_truth.clear();

        // Fill jets
        if (nJet >= 1){
            jet1_puppi_pt.push_back(jetpuppi.at(0)->PT);
            jet1_puppi_eta.push_back(jetpuppi.at(0)->Eta);
            jet1_puppi_phi.push_back(jetpuppi.at(0)->Phi);
            jet1_puppi_q.push_back(jetpuppi.at(0)->Charge);
        }

        // Fill truth jets
        for (size_t i=0; i<genjet.size(); ++i){
            if (genjet.at(i)->PT < jet_pt_lo){ continue; }
            jet1_pt_truth.push_back(genjet.at(i)->PT);
            jet1_eta_truth.push_back(genjet.at(i)->Eta);
            jet1_phi_truth.push_back(genjet.at(i)->Phi);
            jet1_q_truth.push_back(genjet.at(i)->Charge);
        }

        // Fill MET
        try{
            met = puppimet.at(0)->MET;
            met_eta = puppimet.at(0)->Eta;
            met_phi = puppimet.at(0)->Phi;
        }catch (const std::out_of_range& oor){
            std::cerr << "Out of range error when accessing MET vector: " << oor.what() << std::endl;
            return;
        }

        // Fill HT
        ht = 0.;
        for (size_t i=0; i<jetpuppi.size(); ++i){
            if (jetpuppi.at(i)->PT > jet_pt_lo){
                ht += jetpuppi.at(i)->PT;
            }
        }

        // Transverse mass of leading two leptons
        //if (lep1_pt.size() >= 1){
        //    mt1.push_back(std::sqrt(2*lep1_pt.at(0)*met*(1-std::cos(lep1_phi.at(0)-met_phi))));
        //}
        //if (lep1_pt.size() >= 2){
        //    mt2.push_back(std::sqrt(2*lep2_pt.at(1)*met*(1-std::cos(lep2_phi.at(1)-met_phi))));
        //}

        myskim->Fill();
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



