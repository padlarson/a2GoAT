#ifndef __AdlarsonPhysics_h__
#define __AdlarsonPhysics_h__

#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <string> 

#include <TCutG.h>
#include <TFile.h>
#include "GTreeManager.h"
#include "GHistBGSub2.h"
#include "TVector.h"
#include "GTrue.h"
#include <random>
#include <vector>
#include <map>
#include <cmath>
#include "TTree.h"

#include <APLCON.hpp>
template<typename T>
std::vector<T> operator+(const std::vector<T>& v1, const std::vector<T>& v2) {
    std::vector<T> v = v1;
    v.insert(v.end(),v2.begin(),v2.end());
    return v;
}

class	AdlarsonPhysics  : public GTreeManager{

private:
    TRandom3*       pRandoms;
    Bool_t          MC;
    Bool_t          MC_weight  = false;
    Bool_t          MC_etapr   = false;
    Double_t        MCw  = 1.0;
    Bool_t          MCJuly14   = true;
    Bool_t          MCOctDec14 = false;

    Bool_t          eight_clusters;
    // Histograms and scatterplots
    // True
    GH1*            true_BeamE;
    GH1*            true_BeamE_weight;
    TH1*            true_norm; // normalisation factor

    TH2*            true_th_p_v_th_etapr_CM;
    TH2*    true_th_v_E_p;
    TH2*    true_th_v_E_eta_g;
    TH2*    true_th_v_E_eta_6g;
    TH2*    true_th_v_E_pi0_g;
    // Kinfit tests
    TH2*    true_six_dth_vs_th_p;
    TH2*    true_six_dR_vs_R_p;
    TH2*    true_six_dR_vs_det_p;
    TH2*    true_six_fit_dz_v_z;
    TH2*    fit_six_fit_dz_v_z;
    TH2*    true_six_fit_dz_v_p_th;
    GH1*            true_z_after_fit;
    GH1*            true_z_after_final_fit;
    // Physics Results
    TH1*            true_etapr_diff_distr;
    TH2*            true_eta_pr_production;
    TH2*            true_DP;
    TH2*            true_DP_005;
    TH2*            true_DP_075;
    TH2*            true_DP_010;
    TH2*            true_DP_015;
    TH1*            true_phy_DP_020;
    TH1*            true_phy_DP_015;
    TH1*            true_phy_DP_010;
    TH1*            true_phy_DP_075;
    TH1*            true_phy_DP_005;
    TH1*            true_M_pi1pi2_e2p;
    TH1*            true_M_etapi_e2p;
    TH1*            true_imng;

    // Sergey Binning
    TH2*            true_DP_SergeyBin;
    TH1*            true_X_SergeyBin;
    TH1*            true_Y_SergeyBin;
    TH1*            true_m_pipi_SergeyBin;
    TH1*            true_m_pipi_sq_SergeyBin;
    TH1*            true_m_pieta_SergeyBin;
    TH1*            true_m_pieta_sq_SergeyBin;

    TH2*    true_six_phy_DP_020;
    TH2*    true_six_phy_DP_015;
    TH2*    true_six_phy_DP_010;
    TH2*    true_six_phy_DP_075;
    TH2*    true_six_phy_DP_005;
    TH2*    true_six_phy_DP_020_pr;
    TH2*    true_six_phy_DP_015_pr;
    TH2*    true_six_phy_DP_010_pr;
    TH2*    true_six_phy_DP_075_pr;
    TH2*    true_six_phy_DP_005_pr;

    // in 6g analaysis
    TH2*    true_six_phy_dX_v_DPbin;
    TH2*    true_six_phy_dY_v_DPbin;
    TH2*    true_six_phy_dMpipi_v_Mpipi;
    TH2*    true_phy_3pi_IMpipi_v_IMppi;
    TH2*    true_six_phy_dX_v_X;
    TH2*    true_six_phy_dY_v_Y;
    TH2*    true_six_phy_Xtr_v_Xfit;
    TH2*    true_six_phy_Ytr_v_Yfit;
    TH2*    true_six_phy_Xtr_v_Xfit_metapr;
    TH2*    true_six_phy_Ytr_v_Yfit_metapr;

    // for m_etapr forced in kinfit
    TH2*    true_six_phy_dX_v_DPbin_metapr;
    TH2*    true_six_phy_dY_v_DPbin_metapr;
    TH2*    true_six_phy_dMpipi_v_Mpipi_metapr;
    TH2*    true_six_phy_dX_v_X_metapr;
    TH2*    true_six_phy_dY_v_Y_metapr;

    // in 10g analysis
    TH2*    true_ten_phy_dMpipi_v_Mpipi;
    TH2*    true_ten_phy_dX_v_DPbin;
    TH2*    true_ten_phy_dY_v_DPbin;

    // Tagger related
    GH1*            tag_BeamE;
    Double_t        tagger_min, tagger_max;

    // Proton related
    Double_t        TOF_CB;
    Double_t        TOF_CB_proton;
    TLorentzVector  MMp_vec;
    TLorentzVector  g[3], h[3], rc[3], rc_sig[3];


    GHistBGSub*     p_MM;
    GHistBGSub2*    p_th_v_E;
    GHistBGSub2*    p_E_v_dE_all;   // hits
    GHistBGSub2*    p_E_v_dE_cc;    // for random times
    GHistBGSub2*    p_E_v_dE_pr;    // for proton sel
    GHistBGSub2*    p_E_v_TOF;
    GHistBGSub2*    p_E_v_TOF_All;
    GHistBGSub2*    p_E_v_TOF_TAPS_1cl;
    GHistBGSub2*    p_E_v_TOF_CB_All_proton;
    GHistBGSub2*    p_E_v_TOF_CB_2PrID;
    GHistBGSub2*    p_E_v_TOF_CB_best;
    GHistBGSub2*    p_E_v_TOF_after_kfit;
    GHistBGSub2*    p_Erec_v_TOF_after_kfit;
    GHistBGSub2*    p_E_v_TOF_after_kfit_2;

    GHistBGSub*     CB_EnergySum;
    GHistBGSub*     CB_EnergySum_3pi0;
    GHistBGSub*     CB_EnergySum_etapr;

    GH1*     two_rec_IM;
    GH1*     kfit_pdf_2g;
    GH1*     IM2g_fit;

    // Photon related
    TLorentzVector  IM2g_vec;
    TLorentzVector  IM6g_vec;
    TLorentzVector  IM8g_vec;
    TLorentzVector  IM10g_vec;

    // 10g analysis
    GH1*            ten_rec_IM;
    GHistBGSub2*    ten_rec_IM_v_MMp;
    GHistBGSub2*    ten_fit_EvTh_g;

    GH1*            fi_diff_TAPSCB;
    GHistBGSub2*    fi_TAPSvsCB;
    GHistBGSub2*    fi_th_diff_TAPSCB;

    GHistBGSub2*    IMgg_v_det_3pi0_CB;
    GHistBGSub2*    IMgg_v_det_3pi0_TAPS;

    GHistBGSub2*    IM6g_v_det_etaprfit_CB;
    GHistBGSub2*    IM6g_v_det_etaprfit_TAPS;
    GHistBGSub2*    IM6g_v_det_etaprrec_CB;
    GHistBGSub2*    IM6g_v_det_etaprrec_TAPS;

    TH2D            time_TOF;
    TH2D            time_clusters_TAPS;
    TH2D            time_clusters_CB;
    TH2D            time_clusters_CBavg_CBdiff;
    TH2D            time_clusters_CBavg_TAPSdiff;
    TH2D            time_clusters_CBavg_EPT;
    TH1D            time_nr_AllClusters;
    TH1D            time_nr_ClustersinTime;
    TH1D            time_nr_FinalClusterSel;
    TH1D            six_time_TaggedTime;

    // 6g analysis
    // test analysis to test that kinfit APLCON is working properly
    GH1*            test_six_rec_IM;
    GH1*            test_six_fit_chi2;
    GH1*            test_six_fit_pdf;
    GHistBGSub2*    test_six_fit_Pulls;

    // normal 6g analysis

    GH1*            six_rec_IM;
    GHistBGSub2*    six_rec_IM_v_MMp;
    GH1*            six_rec_IM_eta2pi;
    GHistBGSub2*    six_rec_EvTh_6g;
    GHistBGSub2*    six_rec_EvTh_7g;
    GHistBGSub2*    time_clusters_CB_3pi0;

    GH1*            six_fit_chi2;
    GH1*            six_fit_pdf;
    GH1*            six_fit_etaprfinal_pdf;

    GH1*            NI6g;
    GH1*            NIeta2pi0;
    GH1*            NItetapr;

    GHistBGSub2*    NI3pi0vPDF;
    GHistBGSub2*    NIeta2pi0vPDF;
    GHistBGSub2*    NItetaprvPDF;

    GHistBGSub2*    proton_fit_e_v_th;
    GHistBGSub2*    proton_fit_e_v_th_final;

    GH1*            six_fit_IM;
    GHistBGSub2*    six_fit_IM_ncl;
    GHistBGSub2*    six_fit_IM_vz;

    GH1*            six_fit_which_place_best_3pi_cand;
    GH1*            six_fit_which_place_best_etapr_cand;
    GHistBGSub2*    six_fit_PDF_eta2pi_v_3pi;
    GHistBGSub2*    six_fit_PDF_eta2pi_v_3pi_2;
    GHistBGSub2*    six_fit_PDF_eta2pi_v_3pi_4;

    GH1*            six_fit_IM_eta2pi0_b;
    GH1*            six_fit_IM_eta2pi0_c;
    GH1*            six_fit_IM_eta2pi0_d;
    GH1*            six_fit_IM_eta2pi0_e;
    GH1*            six_fit_IM_eta2pi0_f;
    GH1*            six_fit_IM_eta2pi0_g;

    GHistBGSub2*    six_fit_PDF_eta2pi_v_Meta2pi;
    GH1*            six_fit_IM_3pi;
    GH1*            six_fit_IM_eta2pi;
    GH1*            six_fit_best_eta;
    GH1*            six_fit_best_eta_rec;

    GHistBGSub2*    six_fit_IM_eta2pi_v_ncl;

    GHistBGSub2*    six_fit_EvTh_g;
    GHistBGSub2*    six_fit_EvTh_g_final;
    GHistBGSub2*    six_fit_best_etapr_eta_E_v_th;
    GHistBGSub2*    six_fit_best_etapr_pi_E_v_th;
    GHistBGSub2*    six_fit_best_3pi0_pi_E_v_th;
    GH1*            six_fit_best_2pi;
    GH1*            six_fit_best_2pi_rec;

    GHistBGSub2*    six_fit_best_etapr_eta_mgg_v_thth;
    GHistBGSub2*    six_fit_best_etapr_pi0_mgg_v_thth;
    GHistBGSub2*    six_fit_best_3pi0_mgg_v_thth;

    // to check the energy of the  eta 2pi0 system vs its inv mass
    GHistBGSub2*    six_fit_best_eta_IM_v_E;
    GHistBGSub2*    six_fit_best_2pi_IM_v_E;

    GHistBGSub2*    six_fit_mgg_v_eth;
    GHistBGSub2*    six_fit_mgg_v_eth_2;
    GHistBGSub2*    six_fit_mgg_v_eth_3;
    GHistBGSub2*    six_fit_mgg_v_CB;
    GHistBGSub2*    six_fit_mgg_v_CB_2;
    GHistBGSub2*    six_fit_mgg_v_CB_3;
    GHistBGSub2*    six_fit_mgg_v_CB_4;
    GHistBGSub2*    six_fit_mgg_v_TAPS;
    GHistBGSub2*    six_fit_mgg_v_TAPS_2;
    GHistBGSub2*    six_fit_mgg_v_TAPS_3;
    GHistBGSub2*    six_fit_mgg_v_TAPS_4;

    GHistBGSub2*    six_rec_m6g_sig_v_eth;
    GHistBGSub2*    six_fit_fitted_etapr_de_v_eth;
    GHistBGSub2*    six_fit_fitted_3pi0_de_v_eth;
    GHistBGSub2*    six_fit_fitted_de_v_eth;
    GHistBGSub2*    six_fit_fitted_dth_v_eth;
    GHistBGSub2*    six_fit_fitted_p_th_v_det;
    GHistBGSub2*    six_fit_fitted_p_fi_v_det;

    // to check the energy of the 3pi0 system vs its inv mass
    GHistBGSub2*    six_fit_best_3pi_IM_v_E;
    GHistBGSub2*    six_phy_3pi_IMpipi_v_IMppi;

    GHistBGSub2*    six_phy_etapr_v_BeamE;
    GHistBGSub2*    six_phy_etapr_eta2pi_v_BeamE;
    GHistBGSub2*    six_phy_etapr_v_EPT;
    GHistBGSub2*    six_phy_etapr_eta2pi_v_EPT;

    GHistBGSub2*    six_phy_etapr_prod_diff_distr;
    GH1*            six_phy_etapr_prod_diff_distr_metapr;

    GHistBGSub2*    six_phy_DP_020;
    GHistBGSub2*    six_phy_DP_015;
    GHistBGSub2*    six_phy_DP_010;
    GHistBGSub2*    six_phy_DP_075;
    GHistBGSub2*    six_phy_DP_005;
    GHistBGSub2*    six_phy_DP_020_P1;
    GHistBGSub2*    six_phy_DP_015_P1;
    GHistBGSub2*    six_phy_DP_010_P1;
    GHistBGSub2*    six_phy_DP_075_P1;
    GHistBGSub2*    six_phy_DP_005_P1;
    GHistBGSub2*    six_phy_DP_X;
    GHistBGSub2*    six_phy_DP_Y;
    GHistBGSub2*    six_phy_M_pi1pi2_v_etapr;
    GHistBGSub2*    six_phy_M_etapi_v_etapr;
    GHistBGSub2*    six_phy_M_pi1pi2_v_etapr2;
    GHistBGSub2*    six_phy_M_etapi_v_etapr2;

    GH1*            six_phy_DP_020_pr;
    GH1*            six_phy_DP_015_pr;
    GH1*            six_phy_DP_010_pr;
    GH1*            six_phy_DP_075_pr;
    GH1*            six_phy_DP_005_pr;
    GH1*            six_phy_DP_X_pr;
    GH1*            six_phy_DP_Y_pr;
    GH1*            six_phy_M_pi1pi2_v_etapr_fit;
    GH1*            six_phy_M_etapi_v_etapr_fit;
    GH1*            six_phy_M_pi1pi2_v_etapr_fit2;
    GH1*            six_phy_M_etapi_v_etapr_fit2;

    TH1*            six_fit_IM_eta2pi_prompt;
    TH1*            six_fit_IM_eta2pi_random;

    TString   tree_file_name;
    TFile*    f_tree2;
    TTree*    tree2;

    double    fWeight;
    int       fNclusters;
    double    fTaggerenergy;
    double    fProton_th_fit;
    double    fProton_E_fit;
    double    fZ_vx_fit;
    double    fX;
    double    fY;
    double    fPeta2pi;
    double    fP3pi;
    double    fPetapr;
    double    fCosth_epr_cm;
    int       fDP005;
    int       fDP075;
    int       fDP010;
    int       fDP015;
    double    fMeta2pi;
    double    fMpipi;
    double    fMetapi1;
    double    fMetapi2;
    double    fCosth_epr_cm_pr;
    double    fXpr;
    double    fYpr;
    int       fDP005pr;
    int       fDP075pr;
    int       fDP010pr;
    int       fDP015pr;
    double    fMpipipr;
    double    fMetapi1pr;
    double    fMetapi2pr;


    //    Analysis 8g
    GH1*      eight_rec_IM;
    GH1*      kfit_pdf_8g;
    GH1*      IM8g_fit;

    // Kinfit related variables 10g
    GH1*            kfit_chi2_10g;
    GH1*            kfit_pdf_10g;
    GHistBGSub2*    kfit_Pulls_10g;

    GH1*            IM10g_fit;  
    GHistBGSub2*    ten_fit_PDF_eta2pi_v_eta6g;
    GH1*            IM10g_fit_best_cand;

    GH1*            ten_fit_PDF_eta2pi;
    GHistBGSub2*    ten_fit_X_v_pdf_eta2pi;
    GHistBGSub2*    ten_fit_Y_v_pdf_eta2pi;

    GHistBGSub2*    ten_phy_etapr_prod_diff_distr;
    GHistBGSub2*    ten_phy_etapr_prod_diff_distr_metapr;

    // proton identified from TAPS_E vs VETO_dE

    TFile*          cutFile;        // File which contains EdE cut
    TCutG*          cutProtonTAPS;
    TFile*          cutFile2;       // File which contains EToF cut
    TCutG*          cutProtonETOF;
    TFile*          cutFile3;       // File which contains MinIon cut
    TCutG*          cutMinIon;

    TCutG*          cutPionTAPS;
    TCutG*          cutElectronTAPS;

    TFile*          Evth_g_sel;
    TCutG*          sixg_cand;
    TFile*          Evth_7g_sel;
    TCutG*          seveng_cand;
    TFile*          PDF_cut_file;
    TCutG*          PDF_cut;

    TFile*          p_ToF_kfit_file;
    TCutG*          p_ToF_kfit_cut;

    TFile*          g_unc;
    TFile*          p_unc;
    // histograms with unc

    TH2*            g_e;
    TH2*            g_th;
    TH2*            g_phi;
    TH1*            g_R;
//    TH2*            p_TAPS_e;
    TH1*            p_TAPS_th;
    TH1*            p_TAPS_fi;
    TH1*            p_TAPS_R;

    TFile*          g_unc_vz;
    TFile*          p_unc_vz;
    TH2*            g_e_vz;
    TH2*            g_th_vz;
    TH2*            g_phi_vz;
    TH1*            p_TAPS_th_vz;
    TH1*            p_TAPS_fi_vz;

    TFile*          weight_bkgd;
    TH2*            MCw_bkgd;
    TFile*          weight_bkgd2;
    TH2*            MCw_bkgd2;
    TFile*          weight_bkgd3;
    TH2*            MCw_bkgd3;

    TFile*          etapr_MC_unc;
    TH1*            CB_unc;

    // True LorentzVectors
    TLorentzVector  eta_true;
    TLorentzVector  pi01_true;
    TLorentzVector  pi02_true;
    TLorentzVector  etapr_true[3];

    Int_t           diffbin, diffbin_pr;
    Double_t        Xtrue, Ytrue, bw;
    Int_t           DPnrTrue020, DPnrTrue015, DPnrTrue010, DPnrTrue075, DPnrTrue005;
    Double_t        m_etapi01True, m_etapi02True, m_2pi0True;

    Double_t        Xfit, Yfit, Xfit_pr, Yfit_pr;

    // Reconstructed Lorentz Vectors
    std::vector<int> ClustersInTime;
    std::vector<int> IgnoreTAPSCluster;
    std::vector<int> FinalClusterSelection;

    std::vector<int> PbWO4;
    std::vector<int> ring3_or_ring4_CB;
    std::vector<int> Is_CB;
    Bool_t edge_CB;

    UInt_t nrprotons;
    UInt_t iprtrack;

    TLorentzVector proton_vec;
    Double_t MMp;

    TLorentzVector etapr_sixgam[6];
    TLorentzVector dir3pi0_sixgam[6];

    TLorentzVector etapr_tengam[10];

    static Int_t perm4g[9][4];
    static Int_t perm6g[15][6];
    static Int_t perm8g[28][8];
    static Int_t perm6outof10g[210][6];

    GTrue   etapr_2gTrue;
    GTrue   etapr_6gTrue;
    GTrue   etapr_10gTrue;
    GTrue   threepi_etapi;

    std::vector<int>  detnr;
    std::vector<bool> CB_region;
    std::vector<bool> Is_CB_6g;
    std::vector<bool> Is_CB_10g;

    std::vector<double> obs;
    std::vector<double> unc;

    std::vector<double> Legendre; //etapr photo pr Legendre values

    std::vector<TLorentzVector> photons_rec;
    std::vector<TLorentzVector> photons_rec_eight;
    std::vector<TLorentzVector> photons_rec_ten;
    std::vector<TLorentzVector> photons_fit;
    std::vector<TLorentzVector> photons_fit_final;
    std::vector<TLorentzVector> photons_fit_final_metapr;
    std::vector<TLorentzVector> photons_fit_8g;
    std::vector<TLorentzVector> photons_fit_10g;

    Double_t    sigma_eta;
    Double_t    sigma_pi0;


    typedef std::pair<Int_t, std::vector<Double_t>> corr_pair;
    std::map<UInt_t, std::vector<Double_t>> CB_Ecorr;

    std::vector<Double_t> TAPS_EPT_toff;
    std::vector<Double_t> TAPS_Ecorr;

    std::vector<Double_t> CBtime_corr;
    std::vector<Double_t> CBtime_sgm;
    std::vector<Double_t> TAPS_CB_toff;

    std::vector<Double_t> CBgain;         //gain corr factors norm pi0 from 3pi0 at pipeak
    std::vector<Double_t> CBsmear;
    std::vector<Double_t> TAPSgain;
    std::vector<Double_t> TAPSsmear;

    std::vector<Double_t> TAPSth_corr;
    std::vector<Double_t> R_TAPS;
    std::vector<Double_t> R_TAPS_corr;

    typedef std::pair<Int_t, std::vector<Double_t>> time_pair;
    std::map<UInt_t, std::vector<Double_t>> MC_time_jump;


    typedef std::pair<std::vector<Int_t>, Double_t> comb;
    std::vector<comb> threepi_comb;
    std::vector<comb> etatwopi_comb;
    std::vector<comb> teng_comb;

    TFile*          thcorr_CB;        // File which contains TProfile
    TProfile*       dthvth_CB;
    TFile*          thcorr_TAPS;        // File which contains TProfile
    TProfile*       dthvth_TAPS;
    TFile*          Ecorr_CB;         // File which contains TH2F
    TH2F*           EvdetCB;

    TFile*          Ecorr_TAPS;         // File which contains TH2F
    TH2F*           EvdetTAPS;

    TFile*          Ecorr_gamma;        // correction compated to MC 3pi0
    TH2F*           Eth_gamma;


protected:
//    virtual Bool_t  Init(const char* configFile);
    virtual Bool_t  Start();

    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();


    static constexpr bool includeIMconstraint = true;
    static constexpr bool includeVertexFit = true;
    static constexpr size_t nPhotons_two    = 2;
    static constexpr size_t nPhotons_six    = 6;
    static constexpr size_t nPhotons_eight  = 8;
    static constexpr size_t nPhotons_ten    = 10;


    // lightweight structure for linking to fitter
        struct FitParticle{
            void SetFromVector(const TLorentzVector& p_) {
                Ek = p_.E()-p_.M();
                Theta = p_.Theta();
                Phi = p_.Phi();
            }
            void SetFromValues(const Double_t& E, const Double_t& th, const Double_t& ph) {
                Ek = E;
                Theta = th; // in case of TAPS, theta is replaced by Radius R
                Phi = ph;
            }

            static TLorentzVector Make(const std::vector<double>& EkThetaPhi,
                                               const Double_t m);
            static TLorentzVector Make(const FitParticle& p,
                                       const Double_t m) {
            return Make(std::vector<double>{p.Ek, p.Theta, p.Phi}, m);
            }

            std::vector<double*> Link() {
                return {
                        std::addressof(Ek),
                        std::addressof(Theta),
                        std::addressof(Phi)};
            }
            std::vector<double*> LinkSigma() {
                return {
                        std::addressof(Ek_Sigma),
                        std::addressof(Theta_Sigma),
                        std::addressof(Phi_Sigma)};
            }
            std::vector<APLCON::Variable_Settings_t> LinkSettings()
            {
                return{Ek_Setting, Theta_Setting, Phi_Setting};
            }

            void Smear(std::vector<double> unc , int particle);
            void Smear_R(std::vector<double> unc , int particle);

            void APLCONSettings();

            double Ek;
            double Ek_Sigma;
            APLCON::Variable_Settings_t Ek_Setting;
            double Theta;
            double Theta_Sigma;
            APLCON::Variable_Settings_t Theta_Setting;
            double Phi;
            double Phi_Sigma;
            APLCON::Variable_Settings_t Phi_Setting;


        private:
            static std::default_random_engine generator;

        };

    APLCON kinfit2g;
    APLCON kinfit;
    APLCON kinfit_final;
    APLCON kinfit3pi;
    APLCON kinfiteta2pi;
    APLCON kinfit8g;
//    APLCON kinfit8g_etapi;
    APLCON kinfit10g;
    APLCON kinfit10g_eta2pi;

    FitParticle beam2g;
    FitParticle beam;
    FitParticle beam_3pi;
    FitParticle beam_eta2pi;
    FitParticle beam_final;
    FitParticle beam8g;
//    FitParticle beam8g_etapi;
    FitParticle beam10g;
    FitParticle beam10g_eta2pi;

    std::vector<FitParticle> Photons_two;
    std::vector<FitParticle> Photons_six;
    std::vector<FitParticle> Photons_six_3pi;
    std::vector<FitParticle> Photons_six_eta2pi;
    std::vector<FitParticle> Photons_six_final;
    std::vector<FitParticle> Photons_eight;
//    std::vector<FitParticle> Photons_eight_etapi;
    std::vector<FitParticle> Photons_ten;
    std::vector<FitParticle> Photons_ten_eta2pi;

    FitParticle proton2g;
    FitParticle proton;
    FitParticle proton_3pi;
    FitParticle proton_eta2pi;
    FitParticle proton_final;
    FitParticle proton8g;
//    FitParticle proton8g_etapi;
    FitParticle proton10g;
    FitParticle proton10g_eta2pi;

			
public:
    AdlarsonPhysics();
    virtual ~AdlarsonPhysics();

    Bool_t	Init(const char* configfile);
    void            RandomTime();
    void            Tagger_corr();      // corrects tagged energy by subtracting 1.0 MeV for exp data
    void            Time_corr();
    void            Energy_corr();      // corrects theta for CB and TAPS for all clusters (Tracks).
    void            theta_corr();      // corrects theta for CB and TAPS for all clusters (Tracks).

    TLorentzVector GetLVCorrForZ(std::vector<double> EkPThPhi, const double v_z, Int_t &idet, double mass);
    void           CB_TAPS_boundary(); // checks if there are double hits close to CB-TAPS boundary;

    // function where true analysis is done for eta' --> eta 2pi0 --> 6g
    void TrueAnalysis_etapr2g8g(TString s);
    void TrueAnalysis_etapr6g(TString s);
    void TrueAnalysis_etapr10g();
    Double_t TrueAnalysis_threepi_etapi(); // returns weight as function of IM(ng);
    Double_t Get_etapr_weight_MC(Double_t beame, TLorentzVector eta_pr[2]);

    Double_t Get_ESumMC(Double_t &ESum);
    Double_t ESum;
    Double_t ESum_MC;

    double GetWeight3pi1(Double_t M1sq, Double_t M2sq);
    double GetWeight3pi2(Double_t M1sq, Double_t M2sq);
    double GetWeight3pi3(Double_t M1sq, Double_t M2sq);
    void   Kinfit_test();

    double GetGain(Double_t E, Double_t detnr);

    std::vector<double> Get_unc(Int_t apparatus_nr, Int_t particle, std::vector<double>& obs);
    std::vector<double> Get_unc_R(Int_t apparatus_nr, Int_t particle, std::vector<double>& obs);

    void twogAnalysis( UInt_t ipr );

    // functions specifically related to 6g analysis
    void sixgAnalysis( UInt_t ipr );
    void GetBest6gCombination(Double_t& sigma_eta, Double_t& sigma_pi0, Double_t& chi2min_eta2pi, Double_t& chi2min_3pi, std::vector<int>& imin_eta2pi, std::vector<int>& imin_3pi, std::vector<comb>& etatwopi_comb, std::vector<comb>& threepi_comb );
    void test_correct_hypothesis(Double_t &prob_etapr, Double_t& chi2min_eta2pi, Double_t& chi2min_3pi, std::vector<Int_t>& set_min, std::vector<int>& imin_eta2pi, std::vector<int>& imin_3pi,  std::vector<comb>& etatwopi_comb, std::vector<comb>& threepi_comb);
    void FillTree();

    void eightgAnalysis(UInt_t ipr);

    // functions specifically related to 10g analysis
    void tengAnalysis(UInt_t ipr );
    void GetBest10gCombination(Double_t& sigma_eta, Double_t& sigma_pi0, Double_t& chi2min_eta3pi, std::vector<int>& imin_eta3pi2pi );

    const std::vector<TLorentzVector> ClusterEnergyCorr();


    Int_t diff_distr(const Double_t beam_e, TLorentzVector& fin );
    void DalitzPlot(const TLorentzVector g[3] , Double_t &X, Double_t &Y, Double_t bw, Int_t &DP_nr);
    void m2pi0_metapi0(  TLorentzVector g[3], Double_t &m_etapi01, Double_t &m_etapi02, Double_t &m_2pi0 );

};
#endif

