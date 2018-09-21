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
        SUSY_Upgrade_Skimmer():d_ana::basicAnalyzer(){
            TString pathJetPtCorr = "./data/cumulativeGraphs_phaseII.root";
            std::cout << "Loading file from " << pathJetPtCorr << "." << std::endl;
            f_SMEAR = new TFile(pathJetPtCorr, "READ");
        }
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
        double coneVeto(double pt, double eta, double phi, d_ana::dBranchHandler<GenParticle>& genpart, std::string particle);
        template <typename T> bool isMatched(const GenParticle* truthParticle, const T particle);
        bool isMatched(const GenParticle* truthParticle, const double pt, const double eta, const double phi);
        template <typename T> void ppp(const char* text, const size_t idx, const size_t noParticles, const T particle, const int pid, const int status, const char* addText="") const;
        template <typename T> void pppWpidWstatus(const char* text, const size_t idx, const size_t noParticles, const T particle, const char* addText="") const;
        double getCorrSigma(double sigma, double genpt, double geneta);

        // Tree
        TTree* myskim;

        // Cutflow control
        static constexpr bool fill_rle = false;
        static constexpr bool event_by_event_comparison_primary = false;
        static constexpr bool event_by_event_comparison_secondary = false;
        static constexpr bool dump_genpart = false;
        static constexpr bool use_full_truth = false;
        static constexpr bool logdebug = false;

        // Cut variables
        static constexpr double el_pt_lo = 5.;
        static constexpr double el_pt_hi = 30.;
        static constexpr double mu_pt_lo = 5.;
        static constexpr double mu_pt_hi = 30.;
        static constexpr double el_eta_hi = 1.6;
        static constexpr double mu_eta_hi = 2.4;
        static constexpr double jet_pt_lo = 25.;
        static constexpr double mass_el = .000511;
        static constexpr double mass_mu = .105658;
        static constexpr double iso_cut_rel = .5;
        static constexpr double iso_cut_abs = 5.;
        static constexpr double jet_or_dr = .4;
        static constexpr double truth_match_diff_pt_rel = .5;
        static constexpr double truth_match_diff_dr = .1;
        static constexpr double wght_tau_veto = 1.;
        static constexpr double wght_hf_veto = 1.;
        static constexpr double wght_gen_iso = 1.;

        // Event variables
        double genWeight, nTot, xs;

        // Signal variables
        double mN1, mN2;
        TH2D* susy_masses = new TH2D("susy_masses", "susy_masses", 401, -1.25, 1001.25, 401, -1.25, 1001.25);

        // Electron vectors
        std::vector<double> el_pt, el_eta, el_phi, el_sumPt;
        std::vector<int> el_q;
        std::vector<bool> el_matched, el_st20to30;
        std::vector<int> el_mother;

        // Electron vectors before isolation
        std::vector<double> el_recoId_pt, el_recoId_eta, el_recoId_phi, el_recoId_sumPt;
        std::vector<int> el_recoId_q;

        // Electron truth vectors
        std::vector<double> el_pt_truth, el_eta_truth, el_phi_truth;

        // Matched truth electron variables
        //std::vector<double> el1_pt_truth_matched, el1_eta_truth_matched, el1_phi_truth_matched, el2_pt_truth_matched, el2_eta_truth_matched, el2_phi_truth_matched;
        //std::vector<int> el1_q_truth_matched, el2_q_truth_matched;

        // Unmatched truth electron vectors
        //std::vector<double> el1_pt_truth, el1_eta_truth, el1_phi_truth, el2_pt_truth, el2_eta_truth, el2_phi_truth;
        //std::vector<int> el1_q_truth, el2_q_truth;

        // Muon vectors
        std::vector<double> mu_pt, mu_eta, mu_phi, mu_sumPt;
        std::vector<int> mu_q;
        std::vector<bool> mu_matched, mu_st20to30;
        std::vector<int> mu_mother;

        // Muon vectors before isolation
        std::vector<double> mu_recoId_pt, mu_recoId_eta, mu_recoId_phi, mu_recoId_sumPt;
        std::vector<int> mu_recoId_q;

        // Muon truth vectors
        std::vector<double> mu_pt_truth, mu_eta_truth, mu_phi_truth;

        // Matched truth muon variables
        //std::vector<double> mu1_pt_truth_matched, mu1_eta_truth_matched, mu1_phi_truth_matched, mu2_pt_truth_matched, mu2_eta_truth_matched, mu2_phi_truth_matched;
        //std::vector<int> mu1_q_truth_matched, mu2_q_truth_matched;

        // Unmatched truth muon vectors
        //std::vector<double> mu1_pt_truth, mu1_eta_truth, mu1_phi_truth, mu2_pt_truth, mu2_eta_truth, mu2_phi_truth;
        //std::vector<int> mu1_q_truth, mu2_q_truth;

        // Lepton vectors
        std::vector<double> lep_pt, lep_eta, lep_phi, lep_mass;

        // Matched truth lepton variables
        //std::vector<double> lep1_pt_truth_matched, lep1_eta_truth_matched, lep1_phi_truth_matched, lep1_mass_truth_matched, lep2_pt_truth_matched, lep2_eta_truth_matched, lep2_phi_truth_matched, lep2_mass_truth_matched;

        // Unmatched truth lepton variables
        //std::vector<double> lep1_pt_truth, lep1_eta_truth, lep1_phi_truth, lep1_mass_truth, lep2_pt_truth, lep2_eta_truth, lep2_phi_truth, lep2_mass_truth;

        // Jet vectors
        std::vector<double> jet_pt, jet_eta, jet_phi;
        std::vector<int> jet_q;

        //// Matched truth jet variables
        //std::vector<double> jet1_pt_truth_matched, jet1_eta_truth_matched, jet1_phi_truth_matched;
        //std::vector<int> jet1_q_truth_matched;

        // MET variables
        double met, met_eta, met_phi;
        //double mht, mht_eta, mht_phi;
        double mht25, mht40, mht50, mht60, mht100, mht150;
        double mlt, mlt_eta, mlt_phi;
        double mhlt25, mhlt25_eta, mhlt25_phi;
        double mhlt40, mhlt40_eta, mhlt40_phi;
        double mhlt50, mhlt50_eta, mhlt50_phi;
        double mhlt60, mhlt60_eta, mhlt60_phi;
        double mhlt100, mhlt100_eta, mhlt100_phi;
        double mhlt150, mhlt150_eta, mhlt150_phi;
        double PFmet, PFmet_eta, PFmet_phi;
        double genmet, genmet_eta, genmet_phi;
        double genpumet, genpumet_eta, genpumet_phi;

        // Other variables
        int nLep, nEl, nMu, nSoftLep, nSoftEl, nSoftMu, nBJet, nW, nZ;
        int nJet25, nJet40, nJet50, nJet60, nJet100, nJet150;
        //int nLep_truth;
        double ht25, ht40, ht50, ht60, ht100, ht150;
        double genht25, genht40;
        bool hasSFOS, hasSoftSFOS, hasSFOS_truth, hasSoftSFOS_truth;
        std::vector<double> mllMin, mllMax, mt1, mt2, pt2l;
        double mTauTau;
        unsigned int ZtoLL;
        bool crazyMuon50, crazyMuon200, crazyMuon500;

        //// Guess origin of leptons
        //TH2D* mu1_pt_origin_nghbr = new TH2D("mu1_pt_origin_nghbr", "mu1_pt_origin_nghbr", 6, 0., 30., 5, -.5, 4.5);
        //TH2D* mu1_pt_origin_cone = new TH2D("mu1_pt_origin_cone", "mu1_pt_origin_cone", 6, 0., 30., 5, -.5, 4.5);
        //TH2D* mu2_pt_origin_nghbr = new TH2D("mu2_pt_origin_nghbr", "mu2_pt_origin_nghbr", 6, 0., 30., 5, -.5, 4.5);
        //TH2D* mu2_pt_origin_cone = new TH2D("mu2_pt_origin_cone", "mu2_pt_origin_cone", 6, 0., 30., 5, -.5, 4.5);
        std::vector<double> mu_origin_nghbr, mu_origin_cone;
        std::vector<double> mu_pt5to10_origin_nghbr, mu_pt5to10_origin_cone;
        std::vector<double> mu_pt10to20_origin_nghbr, mu_pt10to20_origin_cone;
        std::vector<double> mu_pt20to30_origin_nghbr, mu_pt20to30_origin_cone;
        std::vector<double> el_origin_nghbr, el_origin_cone;
        std::vector<double> el_pt5to10_origin_nghbr, el_pt5to10_origin_cone;
        std::vector<double> el_pt10to20_origin_nghbr, el_pt10to20_origin_cone;
        std::vector<double> el_pt20to30_origin_nghbr, el_pt20to30_origin_cone;

        // Real lepton efficiency histograms
        TH2D* rle_el_num = new TH2D("rle_el_num", "rle_el_num", 6, 0., 30., 8, 0., 4.);
        TH2D* rle_el_den = new TH2D("rle_el_den", "rle_el_den", 6, 0., 30., 8, 0., 4.);
        TH2D* rle_mu_num = new TH2D("rle_mu_num", "rle_mu_num", 6, 0., 30., 8, 0., 4.);
        TH2D* rle_mu_den = new TH2D("rle_mu_den", "rle_mu_den", 6, 0., 30., 8, 0., 4.);

        TFile * f_SMEAR;
};

#endif /* SUSY_Upgrade_Skimmer_H_ */
