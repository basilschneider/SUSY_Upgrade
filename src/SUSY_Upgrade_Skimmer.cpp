#include "interface/SUSY_Upgrade_Skimmer.h"



void SUSY_Upgrade_Skimmer::analyze(size_t childid /* this info can be used for printouts */){

    /*
     * This skeleton analyser runs directly on the Delphes output.
     * It can be used to create histograms directly or a skim.
     * If a skim is created, a new input configuration will be written automatically
     * and stored in the output directory together with the ntuples.
     * The skim can contain delphes objects again or can be flat. This is up
     * to the user.
     * Examples for both are given here.
     *
     * The same skeleton can be used to read the skim. Please refer to the comments
     * marked with "==SKIM=="
     *
     * These parts are commented, since the code is supposed to work as an example without
     * modifications on Delphes output directly.
     */



    /*
     * Define the branches that are to be considered for the analysis
     * This branch handler (notice the "d")
     * is used to run directly in Delphes output.
     * For skimmed ntuples, see below
     */
    d_ana::dBranchHandler<Electron> elecs(tree(),"Electron");
    d_ana::dBranchHandler<Muon> muontight(tree(),"MuonTight");
    d_ana::dBranchHandler<Jet> jetpuppi(tree(), "JetPUPPI");
    d_ana::dBranchHandler<MissingET> puppimet(tree(), "PuppiMissingET");
    /*
     * Other branches might be the following
     * (for a full list, please inspect the Delphes sample root file with root)
     * For the Delphes class description, see $DELPHES_PATH/classes/DelphesClasses.h
     */
    //
    //d_ana::dBranchHandler<HepMCEvent>  event(tree(),"Event");
    //d_ana::dBranchHandler<Muon>        muontight(tree(),"MuonTight");
    //d_ana::dBranchHandler<Muon>        muonloose(tree(),"MuonLoose");
    //d_ana::dBranchHandler<Jet>         jetpuppi(tree(), "JetPUPPI");
    //d_ana::dBranchHandler<MissingET>   puppimet(tree(), "PuppiMissingET");
    ////d_ana::dBranchHandler<Weight>        rwgt(tree(), "Rwgt");
    //d_ana::dBranchHandler<ScalarHT>    scalarht(tree(), "ScalarHT");
    //d_ana::dBranchHandler<GenParticle> genpart(tree(),"Particle");
    //d_ana::dBranchHandler<Jet>         genjet(tree(),"GenJet");
    //d_ana::dBranchHandler<Jet>         jet(tree(),"Jet");
    //d_ana::dBranchHandler<Muon>        muontight(tree(),"MuonTight");
    //d_ana::dBranchHandler<Muon>        muonloose(tree(),"MuonLoose");
    //d_ana::dBranchHandler<Photon>      photon(tree(),"Photon");
    //d_ana::dBranchHandler<MissingET>   met(tree(),"MissingET");


    /* ==SKIM==
     *
     * If a skim of the Delphes outout was created in a way indicated
     * further below, use the tBranchHandler (please notice the "t")
     * to access vectors of objects...
     *
     */
    // d_ana::tBranchHandler<std::vector<Electron> > electrons(tree(),"Electrons");

    /*==SKIM==
     *
     * Or an object directly
     *
     */
    //d_ana::tBranchHandler<MissingET> met(tree(),"MET");



    /*
     * Always use this function to add a new histogram (can also be 2D)!
     * Histograms created this way are automatically added to the output file
     */
    //TH1* histo=addPlot(new TH1D("histoname1","histotitle1",100,0,100),"p_{T} [GeV]","N_{e}");


    /*
     * If (optionally) a skim or a flat ntuple is to be created, please use the following function to initialize
     * the tree.
     * The output files will be written automatically, and a config file will be created.
     */
    TTree* myskim=addTree();
    /*
     * Add a simple branch to the skim
     */
    //Double_t elecPt=0;
    //myskim->Branch("elecPt", &elecPt);
    /*
     * Or store a vector of objects (also possible to store only one object)
     */
    //std::vector<Electron> skimmedelecs;
    //myskim->Branch("Electrons",&skimmedelecs);

    // Cut variables
    double el_pt_lo = 5.;
    double el_pt_hi = 30.;
    double mu_pt_lo = 5.;
    double mu_pt_hi = 30.;
    double jet_pt_lo = 25.;

    // Event variables
    double genWeight, nTot;
    myskim->Branch("genWeight", &genWeight);
    myskim->Branch("nTot", &nTot);

    // Electron variables
    double el1_pt, el1_eta, el1_phi, el2_pt, el2_eta, el2_phi;
    int el1_q, el2_q;
    myskim->Branch("el1_pt", &el1_pt);
    myskim->Branch("el1_eta", &el1_eta);
    myskim->Branch("el1_phi", &el1_phi);
    myskim->Branch("el1_q", &el1_q);
    myskim->Branch("el2_pt", &el1_pt);
    myskim->Branch("el2_eta", &el1_eta);
    myskim->Branch("el2_phi", &el1_phi);
    myskim->Branch("el2_q", &el1_q);

    // Muon variables
    double mu1_tight_pt, mu1_tight_eta, mu1_tight_phi, mu2_tight_pt, mu2_tight_eta, mu2_tight_phi;
    int mu1_tight_q, mu2_tight_q;
    myskim->Branch("mu1_tight_pt", &mu1_tight_pt);
    myskim->Branch("mu1_tight_eta", &mu1_tight_eta);
    myskim->Branch("mu1_tight_phi", &mu1_tight_phi);
    myskim->Branch("mu1_tight_q", &mu1_tight_q);
    myskim->Branch("mu2_tight_pt", &mu2_tight_pt);
    myskim->Branch("mu2_tight_eta", &mu2_tight_eta);
    myskim->Branch("mu2_tight_phi", &mu2_tight_phi);
    myskim->Branch("mu2_tight_q", &mu2_tight_q);

    // Jet variables
    double jet1_puppi_pt, jet1_puppi_eta, jet1_puppi_phi;
    int jet1_puppi_q;
    myskim->Branch("jet1_puppi_pt", &jet1_puppi_pt);
    myskim->Branch("jet1_puppi_eta", &jet1_puppi_eta);
    myskim->Branch("jet1_puppi_phi", &jet1_puppi_phi);
    myskim->Branch("jet1_puppi_q", &jet1_puppi_q);

    // MET variables
    double met, met_eta, met_phi;
    myskim->Branch("met", &met);
    myskim->Branch("met_eta", &met_eta);
    myskim->Branch("met_phi", &met_phi);

    // HT variables
    double ht;
    myskim->Branch("ht", &ht);

    // Cutflow variables
    int nLep, nEl, nMu, nSoftLep, nSoftEl, nSoftMu, nJet, nBJet;
    bool hasSFOS;
    double mll;
    myskim->Branch("nLep", &nLep);
    myskim->Branch("nEl", &nEl);
    myskim->Branch("nMu", &nMu);
    myskim->Branch("nSoftLep", &nSoftLep);
    myskim->Branch("nSoftEl", &nSoftEl);
    myskim->Branch("nSoftMu", &nSoftMu);
    myskim->Branch("nJet", &nJet);
    myskim->Branch("nBJet", &nBJet);
    myskim->Branch("hasSFOS", &hasSFOS);
    myskim->Branch("mll", &mll);

    // Event loop
    size_t nevents=tree()->entries();
    if(isTestMode())
        nevents/=1000;
    for(size_t eventno=0;eventno<nevents;eventno++){
        /*
         * The following two lines report the status and set the event link
         * Do not remove!
         */
        reportStatus(eventno,nevents);
        tree()->setEntry(eventno);

        // Fill event variables
        genWeight = 1.;
        nTot = nevents;

        // Fill electrons
        el1_pt = el1_eta = el1_phi = el2_pt = el2_eta = el2_phi = -99.;
        el1_q = el2_q = -99;
        if (elecs.size() >= 1){
            el1_pt = elecs.at(0)->PT;
            el1_eta = elecs.at(0)->Eta;
            el1_phi = elecs.at(0)->Phi;
            el1_q = elecs.at(0)->Charge;
        }
        if (elecs.size() >= 2){
            el2_pt = elecs.at(1)->PT;
            el2_eta = elecs.at(1)->Eta;
            el2_phi = elecs.at(1)->Phi;
            el2_q = elecs.at(1)->Charge;
        }

        // Fill muons
        mu1_tight_pt = mu1_tight_eta = mu1_tight_phi = mu2_tight_pt = mu2_tight_eta = mu2_tight_phi = -99.;
        mu1_tight_q = mu2_tight_q = -99;
        if (muontight.size() >= 1){
            mu1_tight_pt = muontight.at(0)->PT;
            mu1_tight_eta = muontight.at(0)->Eta;
            mu1_tight_phi = muontight.at(0)->Phi;
            mu1_tight_q = muontight.at(0)->Charge;
        }
        if (muontight.size() >= 2){
            mu2_tight_pt = muontight.at(1)->PT;
            mu2_tight_eta = muontight.at(1)->Eta;
            mu2_tight_phi = muontight.at(1)->Phi;
            mu2_tight_q = muontight.at(1)->Charge;
        }

        // Fill jets
        jet1_puppi_pt = jet1_puppi_eta = jet1_puppi_phi = -99.;
        jet1_puppi_q = -99;
        if (jetpuppi.size() >= 1){
            jet1_puppi_pt = jetpuppi.at(0)->PT;
            jet1_puppi_eta = jetpuppi.at(0)->Eta;
            jet1_puppi_phi = jetpuppi.at(0)->Phi;
            jet1_puppi_q = jetpuppi.at(0)->Charge;
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

        // Cutflow variables
        nLep = nEl = nMu = 0;
        nSoftLep = nSoftEl = nSoftMu = 0;
        nJet = nBJet = 0;
        for (size_t i=0; i<elecs.size(); ++i){
            if (elecs.at(i)->PT > el_pt_lo){
                nLep++;
                nEl++;
                if (elecs.at(i)->PT < el_pt_hi){
                    nSoftEl++;
                    nSoftLep++;
                }
            }
        }
        for (size_t i=0; i<muontight.size(); ++i){
            if (muontight.at(i)->PT > mu_pt_lo){
                nLep++;
                nMu++;
                if (muontight.at(i)->PT < mu_pt_hi){
                    nSoftMu++;
                    nSoftLep++;
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

        // Is a same flavour opposite sign lepton pair present?
        hasSFOS = false;
        if (elecs.size() == 2 && elecs.at(0)->Charge*elecs.at(1)->Charge < 0){
            hasSFOS = true;
        }
        if (muontight.size() == 2 && muontight.at(0)->Charge*muontight.at(1)->Charge < 0){
            hasSFOS = true;
        }

        // Invariant mass of same flavour lepton pair
        mll = -99.;
        if (elecs.size() == 2){
            TLorentzVector el1, el2;
            el1.SetPtEtaPhiM(elecs.at(0)->PT,
                             elecs.at(0)->Eta,
                             elecs.at(0)->Phi,
                             0.000511);
            el2.SetPtEtaPhiM(elecs.at(1)->PT,
                             elecs.at(1)->Eta,
                             elecs.at(1)->Phi,
                             0.000511);
            mll = el1*el2;
        }
        if (muontight.size() == 2){
            TLorentzVector el1, el2;
            el1.SetPtEtaPhiM(muontight.at(0)->PT,
                             muontight.at(0)->Eta,
                             muontight.at(0)->Phi,
                             0.105658);
            el2.SetPtEtaPhiM(muontight.at(1)->PT,
                             muontight.at(1)->Eta,
                             muontight.at(1)->Phi,
                             0.105658);
            mll = el1*el2;
        }

        // Skim
        if (nSoftLep != 2){ continue; }
        if (!hasSFOS){ continue; }

        myskim->Fill();


        /*==SKIM==
         * Access the branches of the skim
         */
        //std::vector<Electron> * skimelecs=electrons.content();
        //for(size_t i=0;i<skimelecs->size();i++){
        //	histo->Fill(skimelecs->at(i).PT);
        //}
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



