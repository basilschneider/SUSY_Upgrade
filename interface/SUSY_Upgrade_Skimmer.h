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
        double DeltaR(double eta1, double eta2, double phi1, double phi2);
        double DeltaPhi(double phi1, double phi2);
        bool isOverlap(const Jet* jet, d_ana::dBranchHandler<Electron>& elecs, d_ana::dBranchHandler<Muon>& muons);
        void effOnTopElec(d_ana::dBranchHandler<Electron>& elecs);
        void effOnTopMuon(d_ana::dBranchHandler<Muon>& muontight);
        void passRandomEfficiency(double eff, Float_t*& ppt);
        int getNghbr(int pid);
        double coneVeto(double pt, double eta, double phi, d_ana::dBranchHandler<GenParticle>& genpart);
        template <typename T> bool isMatched(const GenParticle* truthParticle, const T particle);
        bool isMatched(const GenParticle* truthParticle, const double pt, const double eta, const double phi);

        // Tree
        TTree* myskim;

        // Cutflow control
        static constexpr bool fill_rle = false;
        static constexpr bool event_by_event_comparison = false;

        // Cut variables
        static constexpr double el_pt_lo = 2.;
        static constexpr double el_pt_hi = 30.;
        static constexpr double mu_pt_lo = 2.;
        static constexpr double mu_pt_hi = 30.;
        static constexpr double jet_pt_lo = 25.;
        static constexpr double mass_el = .000511;
        static constexpr double mass_mu = .105658;
        static constexpr double iso_cut_rel = .0;
        static constexpr double iso_cut_abs = 0.;
        static constexpr double jet_or_dr = .4;
        static constexpr double truth_match_diff_pt_rel = .5;
        static constexpr double truth_match_diff_dr = .1;
        static constexpr double wght_tau_veto = 1.;
        static constexpr double wght_hf_veto = 1.;
        static constexpr double wght_gen_iso = .25;

        // Event variables
        double genWeight, nTot, xs;

        // Electron variables
        std::vector<double> el1_pt, el1_eta, el1_phi, el2_pt, el2_eta, el2_phi;
        std::vector<int> el1_q, el2_q;

        // Matched truth electron variables
        std::vector<double> el1_pt_truth_matched, el1_eta_truth_matched, el1_phi_truth_matched, el2_pt_truth_matched, el2_eta_truth_matched, el2_phi_truth_matched;
        std::vector<int> el1_q_truth_matched, el2_q_truth_matched;

        // Unmatched truth electron vectors
        std::vector<double> el1_pt_truth, el1_eta_truth, el1_phi_truth, el2_pt_truth, el2_eta_truth, el2_phi_truth;
        std::vector<int> el1_q_truth, el2_q_truth;

        // Muon variables
        std::vector<double> mu1_pt, mu1_eta, mu1_phi, mu2_pt, mu2_eta, mu2_phi;
        std::vector<int> mu1_q, mu2_q;

        // Matched truth muon variables
        std::vector<double> mu1_pt_truth_matched, mu1_eta_truth_matched, mu1_phi_truth_matched, mu2_pt_truth_matched, mu2_eta_truth_matched, mu2_phi_truth_matched;
        std::vector<int> mu1_q_truth_matched, mu2_q_truth_matched;

        // Unmatched truth muon vectors
        std::vector<double> mu1_pt_truth, mu1_eta_truth, mu1_phi_truth, mu2_pt_truth, mu2_eta_truth, mu2_phi_truth;
        std::vector<int> mu1_q_truth, mu2_q_truth;

        // Lepton variables
        std::vector<double> lep1_pt, lep1_eta, lep1_phi, lep1_mass, lep2_pt, lep2_eta, lep2_phi, lep2_mass;

        // Matched truth lepton variables
        std::vector<double> lep1_pt_truth_matched, lep1_eta_truth_matched, lep1_phi_truth_matched, lep1_mass_truth_matched, lep2_pt_truth_matched, lep2_eta_truth_matched, lep2_phi_truth_matched, lep2_mass_truth_matched;

        // Unmatched truth lepton variables
        std::vector<double> lep1_pt_truth, lep1_eta_truth, lep1_phi_truth, lep1_mass_truth, lep2_pt_truth, lep2_eta_truth, lep2_phi_truth, lep2_mass_truth;

        // Jet variables
        std::vector<double> jet1_puppi_pt, jet1_puppi_eta, jet1_puppi_phi;
        std::vector<int> jet1_puppi_q;

        //// Matched truth jet variables
        //std::vector<double> jet1_pt_truth_matched, jet1_eta_truth_matched, jet1_phi_truth_matched;
        //std::vector<int> jet1_q_truth_matched;

        // MET variables
        double met, met_eta, met_phi;

        // Other variables
        int nLep, nEl, nMu, nSoftLep, nSoftEl, nSoftMu, nJet, nBJet, nW, nZ;
        int nLep_truth;
        double ht;
        bool hasSFOS, hasSoftSFOS, hasSFOS_truth, hasSoftSFOS_truth;
        std::vector<double> mllMin, mllMax, mt1, mt2, pt2l;
        std::vector<double> drLep1Lep2, drLep1Jet1, drLep2Jet1;

        //// Guess origin of leptons
        //TH2D* mu1_pt_origin_nghbr = new TH2D("mu1_pt_origin_nghbr", "mu1_pt_origin_nghbr", 6, 0., 30., 5, -.5, 4.5);
        //TH2D* mu1_pt_origin_cone = new TH2D("mu1_pt_origin_cone", "mu1_pt_origin_cone", 6, 0., 30., 5, -.5, 4.5);
        //TH2D* mu2_pt_origin_nghbr = new TH2D("mu2_pt_origin_nghbr", "mu2_pt_origin_nghbr", 6, 0., 30., 5, -.5, 4.5);
        //TH2D* mu2_pt_origin_cone = new TH2D("mu2_pt_origin_cone", "mu2_pt_origin_cone", 6, 0., 30., 5, -.5, 4.5);
        std::vector<double> mu_pt5to10_origin_nghbr, mu_pt5to10_origin_cone;
        std::vector<double> mu_pt10to30_origin_nghbr, mu_pt10to30_origin_cone;

        // Real lepton efficiency histograms
        TH2D* rle_el_num = new TH2D("rle_el_num", "rle_el_num", 6, 0., 30., 8, 0., 4.);
        TH2D* rle_el_den = new TH2D("rle_el_den", "rle_el_den", 6, 0., 30., 8, 0., 4.);
        TH2D* rle_mu_num = new TH2D("rle_mu_num", "rle_mu_num", 6, 0., 30., 8, 0., 4.);
        TH2D* rle_mu_den = new TH2D("rle_mu_den", "rle_mu_den", 6, 0., 30., 8, 0., 4.);
};

#endif /* SUSY_Upgrade_Skimmer_H_ */
