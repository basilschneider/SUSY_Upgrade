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

    // Electron variables
    double el_pt, el_eta, el_phi;
    int el_q;
    myskim->Branch("el_pt", &el_pt);
    myskim->Branch("el_eta", &el_eta);
    myskim->Branch("el_phi", &el_phi);
    myskim->Branch("el_q", &el_q);

    // Muon variables
    double mu_tight_pt, mu_tight_eta, mu_tight_phi;
    int mu_tight_q;
    myskim->Branch("mu_tight_pt", &mu_tight_pt);
    myskim->Branch("mu_tight_eta", &mu_tight_eta);
    myskim->Branch("mu_tight_phi", &mu_tight_phi);
    myskim->Branch("mu_tight_q", &mu_tight_q);

    // Jet variables
    double jet_puppi_pt, jet_puppi_eta, jet_puppi_phi;
    int jet_puppi_q, jet_puppi_b;
    myskim->Branch("jet_puppi_pt", &jet_puppi_pt);
    myskim->Branch("jet_puppi_eta", &jet_puppi_eta);
    myskim->Branch("jet_puppi_phi", &jet_puppi_phi);
    myskim->Branch("jet_puppi_q", &jet_puppi_q);
    myskim->Branch("jet_puppi_b", &jet_puppi_b);

    // Cutflow variables
    int nLep;
    myskim->Branch("nLep", &nLep);
    bool hasSFOS;
    myskim->Branch("hasSFOS", &hasSFOS);
    double mll;
    myskim->Branch("mll", &mll);

    //std::vector<Event> skimmedevent;
    //myskim->Branch("Event",&skimmedevent);

    //std::vector<Muon> skimmedmuonstight;
    //myskim->Branch("MuonTight", &skimmedmuonstight);

    //std::vector<Muon> skimmedmuonsloose;
    //myskim->Branch("MuonLoose", &skimmedmuonsloose);

    //std::vector<Jet> skimmedjets;
    //myskim->Branch("JetPUPPI", &skimmedjets);

    //std::vector<MissingET> skimmedmet;
    //myskim->Branch("PuppiMissingET", &skimmedmet);

    ////std::vector<Weight> skimmedrwgt;
    ////myskim->Branch("Rwgt", &skimmedrwgt);

    //std::vector<ScalarHT> skimmedht;
    //myskim->Branch("ScalarHT", &skimmedht);

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

        // Fill electrons
        for (size_t i=0; i<elecs.size(); ++i){
            //flat info
            //elecPt=elecs.at(i)->PT;
            //if(elecs.at(i)->PT < 20) continue;
            //
            //or objects
            //skimmedelecs.push_back(*elecs.at(i));
            el_pt = elecs.at(i)->PT;
            el_eta = elecs.at(i)->Eta;
            el_phi = elecs.at(i)->Phi;
            el_q = elecs.at(i)->Charge;
        }

        // Fill muons
        for (size_t i=0; i<muontight.size(); ++i){
            mu_tight_pt = muontight.at(i)->PT;
            mu_tight_eta = muontight.at(i)->Eta;
            mu_tight_phi = muontight.at(i)->Phi;
            mu_tight_q = muontight.at(i)->Charge;
        }

        // Fill jets
        for (size_t i=0; i<jetpuppi.size(); ++i){
            jet_puppi_pt = jetpuppi.at(i)->PT;
            jet_puppi_eta = jetpuppi.at(i)->Eta;
            jet_puppi_phi = jetpuppi.at(i)->Phi;
            jet_puppi_q = jetpuppi.at(i)->Charge;
            jet_puppi_b = jetpuppi.at(i)->BTag;
        }

        // Number of leptons
        nLep = elecs.size() + muontight.size();

        // Is a same flavour opposite sign lepton pair present?
        hasSFOS = false;
        if (elecs.size() == 2 && elecs.at(0)->Charge*elecs.at(1)->Charge < 0){
            hasSFOS = true;
        }
        if (muontight.size() == 2 && muontight.at(0)->Charge*muontight.at(1)->Charge < 0){
            hasSFOS = true;
        }

        // Invariant mass of same flavour lepton pair
        mll = 0.;
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


        //skimmedevent.clear();
        //for(size_t i=0;i<event.size();i++){
        //    skimmedevent.push_back(*event.at(i));
        //}

        //skimmedmuonstight.clear();
        //for(size_t i=0;i<muontight.size();i++){
        //    skimmedmuonstight.push_back(*muontight.at(i));
        //}

        //skimmedmuonsloose.clear();
        //for(size_t i=0;i<muonloose.size();i++){
        //    skimmedmuonsloose.push_back(*muonloose.at(i));
        //}

        //skimmedjets.clear();
        //for(size_t i=0;i<jetpuppi.size();i++){
        //    skimmedjets.push_back(*jetpuppi.at(i));
        //}

        //skimmedmet.clear();
        //for(size_t i=0;i<puppimet.size();i++){
        //    skimmedmet.push_back(*puppimet.at(i));
        //}

        ////skimmedrwgt.clear();
        ////for(size_t i=0;i<rwgt.size();i++){
        ////    skimmedrwgt.push_back(*rwgt.at(i));
        ////}

        //skimmedht.clear();
        //for(size_t i=0;i<scalarht.size();i++){
        //    skimmedht.push_back(*scalarht.at(i));
        //}


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



