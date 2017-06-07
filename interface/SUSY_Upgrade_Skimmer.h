/*
 * SUSY_Upgrade_Skimmer.h
 *
 *  Created on: 24 Aug 2016
 *      Author: jkiesele
 */

#ifndef SUSY_Upgrade_Skimmer_H_
#define SUSY_Upgrade_Skimmer_H_

#include "interface/basicAnalyzer.h"
#include "interface/sampleCollection.h"
#include "classes/DelphesClasses.h"


class SUSY_Upgrade_Skimmer: public d_ana::basicAnalyzer{
    public:
        SUSY_Upgrade_Skimmer():d_ana::basicAnalyzer(){}
        ~SUSY_Upgrade_Skimmer(){}


    private:
        void analyze(size_t id);

        void postProcess();

        void addBranches();
        void clearVectors();
        template <typename T> bool isIsolated(T particle);

        // Tree
        TTree* myskim;

        // Cut variables
        static constexpr double el_pt_lo = 5.;
        static constexpr double el_pt_hi = 30.;
        static constexpr double mu_pt_lo = 5.;
        static constexpr double mu_pt_hi = 30.;
        static constexpr double jet_pt_lo = 25.;
        static constexpr double mass_el = .000511;
        static constexpr double mass_mu = .105658;
        static constexpr double iso_cut = .1;
        static constexpr double truth_match_diff_pt = 3.;
        static constexpr double truth_match_diff_eta = .1;

        // Event variables
        double genWeight, nTot, xs;

        // Electron variables
        std::vector<double> el1_pt, el1_eta, el1_phi, el2_pt, el2_eta, el2_phi;
        std::vector<int> el1_q, el2_q;

        // Truth electron variables
        std::vector<double> el1_pt_truth, el1_eta_truth, el1_phi_truth, el2_pt_truth, el2_eta_truth, el2_phi_truth;
        std::vector<int> el1_q_truth, el2_q_truth;

        // Truth electron vectors
        std::vector<double> el_pt_truth, el_eta_truth, el_phi_truth;
        std::vector<int> el_q_truth;

        // Muon variables
        std::vector<double> mu1_tight_pt, mu1_tight_eta, mu1_tight_phi, mu2_tight_pt, mu2_tight_eta, mu2_tight_phi;
        std::vector<int> mu1_tight_q, mu2_tight_q;

        // Truth muon variables
        std::vector<double> mu1_pt_truth, mu1_eta_truth, mu1_phi_truth, mu2_pt_truth, mu2_eta_truth, mu2_phi_truth;
        std::vector<int> mu1_q_truth, mu2_q_truth;

        // Truth muon vectors
        std::vector<double> mu_pt_truth, mu_eta_truth, mu_phi_truth;
        std::vector<int> mu_q_truth;

        // Lepton variables
        std::vector<double> lep1_pt, lep1_eta, lep1_phi, lep1_mass, lep2_pt, lep2_eta, lep2_phi, lep2_mass;

        // Truth lepton variables
        std::vector<double> lep1_pt_truth, lep1_eta_truth, lep1_phi_truth, lep1_mass_truth, lep2_pt_truth, lep2_eta_truth, lep2_phi_truth, lep2_mass_truth;

        // Jet variables
        std::vector<double> jet1_puppi_pt, jet1_puppi_eta, jet1_puppi_phi;
        std::vector<int> jet1_puppi_q;

        // Truth jet variables
        std::vector<double> jet1_pt_truth, jet1_eta_truth, jet1_phi_truth;
        std::vector<int> jet1_q_truth;

        // MET variables
        double met, met_eta, met_phi;

        // Other variables
        int nLep, nEl, nMu, nSoftLep, nSoftEl, nSoftMu, nJet, nBJet, nW, nZ;
        double ht;
        bool hasSFOS, hasSoftSFOS;
        std::vector<double> mll, mt1, mt2;

};





#endif /* SUSY_Upgrade_Skimmer_H_ */
