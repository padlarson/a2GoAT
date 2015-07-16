#ifndef __AdlarsonPhysics_h__
#define __AdlarsonPhysics_h__

#include <iostream>
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

#include <APLCON.hpp>
template<typename T>
std::vector<T> operator+(const std::vector<T>& v1, const std::vector<T>& v2) {
    std::vector<T> v = v1;
    v.insert(v.end(),v2.begin(),v2.end());
    return v;
}

class	AdlarsonPhysics  : public GTreeManager
{

private:
    // Histograms and scatterplots
    // True
    GH1*            true_BeamE;
    GHistBGSub2*    true_th_p_v_th_etapr_CM;
    GHistBGSub2*    true_th_v_E_p;
    GHistBGSub2*    true_th_v_E_eta_g;
    GHistBGSub2*    true_th_v_E_pi0_g;
    // Kinfit tests
    GHistBGSub2*    true_six_dth_vs_th_p;
    GHistBGSub2*    true_six_z_v_Ncl;
    GHistBGSub2*    true_six_fit_dz_v_z;
    // Physics Results
    TH2*            true_DP;
    TH1*            true_phy_DP;
    TH1*            true_M_pi1pi2_e2p;
    TH1*            true_M_etapi_e2p;

    // in 6g analaysis
    GHistBGSub2*    true_six_phy_dX_v_DPbin;
    GHistBGSub2*    true_six_phy_dY_v_DPbin;
    GHistBGSub2*    true_six_phy_dMpipi_v_Mpipi;

    // in 7g analaysis
    GHistBGSub2*    true_seven_phy_dX_v_DPbin;
    GHistBGSub2*    true_seven_phy_dY_v_DPbin;
    GHistBGSub2*    true_seven_phy_dMpipi_v_Mpipi;


    // Tagger related
    GH1*            tag_BeamE;

    // Proton related
    TLorentzVector  MMp_vec;

    GHistBGSub*     p_nr;
    GHistBGSub*     p_MM;
    GHistBGSub2*    p_th_v_E;
    GHistBGSub2*    p_E_v_dE_all;   // hits
    GHistBGSub2*    p_E_v_dE_cc;    // for random times
    GHistBGSub2*    p_E_v_dE_pr;    // for proton sel
    GHistBGSub2*    p_E_v_TOF;
    GHistBGSub2*    p_E_v_TOF_All;
    GHistBGSub2*    p_E_v_TOF_All_wVeto;
    GHistBGSub2*    p_E_v_TOF_norm_to_c;

    // Photon related
    TLorentzVector  IM6g_vec;
    TLorentzVector  IM10g_vec;

    // 10g analysis
    GH1*            ten_rec_IM;
    GHistBGSub2*    ten_rec_IM_v_MMp;

    GH1*            fi_diff_TAPSCB;
    GHistBGSub2*    fi_TAPSvsCB;
    GHistBGSub2*    fi_th_diff_TAPSCB;

    GHistBGSub2*    IMgg_v_det_b4corr_CB;
    GHistBGSub2*    IMgg_v_det_afcorr_CB;
    GHistBGSub2*    IMgg_v_det_b4corr_TAPS;
    GHistBGSub2*    IMgg_v_det_afcorr_TAPS;

    GHistBGSub2*    IMgg_v_det_2pi0_CB;
    GHistBGSub2*    IMgg_v_det_etapi0_CB;
    GHistBGSub2*    IMgg_v_det_3pi0_CB;
    GHistBGSub2*    IMgg_v_det_3pi0_CB_fit;
    GHistBGSub2*    IMgg_v_det_2pi0_TAPS;
    GHistBGSub2*    IMgg_v_det_etapi0_TAPS;
    GHistBGSub2*    IMgg_v_det_3pi0_TAPS;
    GHistBGSub2*    IMgg_v_det_3pi0_TAPS_fit;

    TH2F            time_clusters_TAPS;
    TH2F            time_clusters_CB;
    TH1D            time_nr_AllClusters;
    TH1D            time_nr_ClustersinTime;
    TH1D            time_nr_CltimeVeto;
    TH1D            six_time_TaggedTime;


   // Analysis 4g

    TLorentzVector  IM4g_vec;
    TLorentzVector  IM4g_fit;

    // 4g analysis
    GH1*            four_rec_IM;
    GH1*            four_fit_IM;

    GH1*            four_fit_chi2;
    GH1*            four_fit_pdf;

    GHistBGSub2*    four_fit_PDF_etapi_v_2pi;


    // 6g analysis
    // test analysis to test that kinfit APLCON is working properly
    GH1*            test_six_rec_IM;
    GH1*            test_six_fit_chi2;
    GH1*            test_six_fit_pdf;
    GHistBGSub2*    test_six_fit_Pulls;

    // normal 6g analysis

    GH1*            six_rec_IM;
    GHistBGSub2*    six_rec_IM_v_MMp;

    GH1*            six_fit_chi2;
    GH1*            six_fit_pdf;
    GH1*            six_fit_eta_pdf;
    GHistBGSub2*    six_fit_Pulls;
    GHistBGSub2*    six_fit_Pulls_g_E_vs_Eth;
    GHistBGSub2*    six_fit_Pulls_g_th_vs_Eth;
    GHistBGSub2*    six_fit_Pulls_g_fi_vs_Eth;

    GHistBGSub2*    six_fit_Pulls_p_th_vs_Eth;
    GHistBGSub2*    six_fit_Pulls_p_fi_vs_Eth;

    GH1*            six_fit_IM;
    GH1*            six_fit_IM_rec;     // rec IM(6g) for events which passed the fit
    GH1*            six_fit_cos_th_cm;

    GHistBGSub2*    six_fit_PDF_eta2pi_v_3pi;
    GHistBGSub2*    six_fit_PDF_eta2pi_v_3pi_2;

    GHistBGSub2*    six_fit_PDF_eta2pi_v_Meta2pi;
    GHistBGSub2*    six_fit_PDF_2_eta2pi_v_Meta2pi;
    GHistBGSub2*    six_fit_eta_PDF_v_Metapr;
    GH1*            six_fit_IM_3pi;
    GH1*            six_fit_IM_eta2pi;
    GH1*            six_fiteta_IM2g;
    GH1*            six_fit_best_eta;
    GHistBGSub2*    six_fit_best_eta_E_v_th;
    GH1*            six_fit_best_2pi;

    // to check the energy of the  eta 2pi0 system vs its inv mass
    GHistBGSub2*    six_fit_best_eta_IM_v_E;
    GHistBGSub2*    six_fit_best_2pi_IM_v_E;

    // to check the energy of the 3pi0 system vs its inv mass
    GHistBGSub2*    six_fit_best_3pi_IM_v_E;
    GHistBGSub2*    six_phy_3pi_IMpipi_v_IMppi;

    GHistBGSub2*    six_phy_etapr_v_BeamE;
    GHistBGSub2*    six_phy_etapr_eta_v_BeamE;

    GHistBGSub2*    six_phy_DP;
    GHistBGSub2*    six_phy_M_pi1pi2_v_etapr;
    GHistBGSub2*    six_phy_M_etapi_v_etapr;

    TH3F            Ekfit_v_Eg_v_detnrCB_4g;
    TH3F            Ekfit_v_Eg_v_detnrTAPS_4g;
    TH3F            Ekfit_v_Eg_v_detnrCB_6g;
    TH3F            Ekfit_v_Eg_v_detnrTAPS_6g;
    TH3F            Ekfit_v_Eg_v_detnrCB_3pi0;
    TH3F            Ekfit_v_Eg_v_detnrCB_eta2pi0;
    TH3F            Ekfit_v_Eg_v_detnrTAPS_3pi0;
    TH3F            Ekfit_v_Eg_v_detnrTAPS_eta2pi0;
    TH3F            Ekfit_v_Eg_v_detnrCB_ng;
    TH3F            Ekfit_v_Eg_v_detnrTAPS_ng;

    TH3F            Ekfit_v_Eg_v_detnrCB_4g_bg;
    TH3F            Ekfit_v_Eg_v_detnrTAPS_4g_bg;
    TH3F            Ekfit_v_Eg_v_detnrCB_6g_bg;
    TH3F            Ekfit_v_Eg_v_detnrTAPS_6g_bg;
    TH3F            Ekfit_v_Eg_v_detnrCB_3pi0_bg;
    TH3F            Ekfit_v_Eg_v_detnrCB_eta2pi0_bg;
    TH3F            Ekfit_v_Eg_v_detnrTAPS_3pi0_bg;
    TH3F            Ekfit_v_Eg_v_detnrTAPS_eta2pi0_bg;
    TH3F            Ekfit_v_Eg_v_detnrCB_ng_bg;
    TH3F            Ekfit_v_Eg_v_detnrTAPS_ng_bg;


    //seven gamma analysis

    GHistBGSub2*    seven_rec_IM_v_MMp;

    GH1*            seven_fit_chi2;
    GH1*            seven_fit_pdf;
    GH1*            seven_fit_eta_pdf;
    GHistBGSub2*    seven_fit_Pulls;

    GH1*            seven_fit_IM;
    GH1*            seven_fit_IM_rec;     // rec IM(6g) for events which passed the fit

    GHistBGSub2*    seven_fit_PDF_eta2pi_v_3pi;
    GH1*            seven_fit_IM_3pi;
    GH1*            seven_fit_IM_eta2pi;
    GH1*            seven_fiteta_IM2g;
    GH1*            seven_fit_best_eta;
    GHistBGSub2*    seven_fit_best_eta_E_v_th;
    GH1*            seven_fit_best_2pi;

    // to check the energy of the  eta 2pi0 system vs its inv mass
    GHistBGSub2*    seven_fit_best_eta_IM_v_E;
    GHistBGSub2*    seven_fit_best_2pi_IM_v_E;

    // to check the energy of the 3pi0 system vs its inv mass
    GHistBGSub2*    seven_fit_best_3pi_IM_v_E;

    GHistBGSub2*    seven_phy_etapr_v_BeamE;
    GHistBGSub2*    seven_phy_etapr_eta_v_BeamE;

    GHistBGSub2*    seven_phy_DP;
    GHistBGSub2*    seven_phy_M_pi1pi2_v_etapr;
    GHistBGSub2*    seven_phy_M_etapi_v_etapr;

    // Kinfit related variables 10g

    GH1*            kfit_chi2_10g;
    GH1*            kfit_pdf_10g;
    GHistBGSub2*    kfit_Pulls_10g;

    GH1*            IM10g_fit;

    // proton identified from TAPS_E vs VETO_dE

    TFile*          cutFile;        // File which contains EdE cut
    TCutG*          cutProtonTAPS;
    TFile*          cutFile2;       // File which contains EToF cut
    TCutG*          cutProtonETOF;
    TFile*          cutFile3;       // File which contains MinIon cut
    TCutG*          cutMinIon;

    TCutG*          cutPionTAPS;
    TCutG*          cutElectronTAPS;

    TFile*          g_unc;
    TFile*          p_unc;
    // histograms with unc
    TH2*            g_CB_e;
    TH2*            g_CB_th;
    TH2*            g_CB_fi;
    TH2*            g_TAPS_e;
    TH2*            g_TAPS_th;
    TH2*            g_TAPS_fi;
    TH2*            p_TAPS_e;
    TH2*            p_TAPS_th;
    TH2*            p_TAPS_fi;

    // True LorentzVectors
    TLorentzVector  eta_true;
    TLorentzVector  pi01_true;
    TLorentzVector  pi02_true;
    TLorentzVector  etapr_true[3];

    Double_t        Xtrue1, Xtrue2, Ytrue;
    Int_t           DPnrTrue1, DPnrTrue2;
    Double_t        m_etapi01True, m_etapi02True, m_2pi0True;

    // Reconstructed Lorentz Vectors
    std::vector<Int_t> ClustersInTime;

    UInt_t nrprotons;
    UInt_t iprtrack;

    TLorentzVector proton_vec;
    Double_t MMp;

    TLorentzVector etapr_sixgam[6];
    TLorentzVector dir3pi0_sixgam[6];

    TLorentzVector etapr_tengam[10];

    static Int_t perm4g[9][4];
    static Int_t perm6g[15][6];
    static Int_t perm6outof7g[7][6];
    static Int_t perm6outof8g[28][6];
    static Int_t perm6outof10g[210][6];

    GTrue   etapr_6gTrue;

    std::vector<int> detnr;
    std::vector<bool> CB_region;
    std::vector<bool> Is_CB_6g;

    std::vector<double> obs;
    std::vector<double> unc;

    std::vector<TLorentzVector> photons_rec;
    std::vector<TLorentzVector> photons_fit;
    std::vector<TLorentzVector> photons_fit_eta;

    Double_t    sigma_eta;
    Double_t    sigma_pi0;

    const Double_t    taggerTimeCut = 6.0;

    typedef std::pair<UInt_t, std::vector<Double_t>> corr_pair;
    std::map<UInt_t, std::vector<Double_t>> CB_Ecorr;

    std::vector<Double_t> TAPS_Ecorr;

    std::vector<Double_t> CBgainMC;         //gain corr factors norm pi0 from 3pi0 at pipeak
    std::vector<Double_t> CBgainEXP;

    typedef std::pair<UInt_t, std::vector<Double_t>> EPT_TAPS_pair;
    std::map<UInt_t, std::vector<Double_t>> TOF_corr;

    TFile*          thcorr_CB;        // File which contains TProfile
    TProfile*       dthvth_CB;
    TFile*          thcorr_TAPS;        // File which contains TProfile
    TProfile*       dthvth_TAPS;


protected:
 //   virtual Bool_t  Init(const char* configFile);
    virtual Bool_t  Start();

    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();


    // choose here what you want to do
    // please also provide GoAT trees with matching MC true information...
    static constexpr bool includeIMconstraint = true;
    static constexpr bool includeVertexFit = true;
    static constexpr size_t nPhotons_four = 4;
    static constexpr size_t nPhotons_six = 6;
    static constexpr size_t nPhotons_ten = 10;



    // lightweight structure for linking to fitter
        struct FitParticle{
            void SetFromVector(const TLorentzVector& p_) {
                Ek = p_.E()-p_.M();
                Theta = p_.Theta();
                Phi = p_.Phi();
            }

            static TLorentzVector Make(const std::vector<double>& EkThetaPhi,
                                               const Double_t m);
            static TLorentzVector Make(const FitParticle& p,
                                       const Double_t m) {
            return Make(std::vector<double>{p.Ek, p.Theta, p.Phi}, m);
            }

            std::vector<double*> Link() {
                return {std::addressof(Ek),
                        std::addressof(Theta),
                        std::addressof(Phi)};
            }
            std::vector<double*> LinkSigma() {
                return {std::addressof(Ek_Sigma),
                        std::addressof(Theta_Sigma),
                        std::addressof(Phi_Sigma)};
            }

            std::vector<APLCON::Variable_Settings_t> LinkSettings()
            {
                return{Ek_Setting, Theta_Setting, Phi_Setting};
            }

//            void Smear(int itr_nr, Bool_t measured );
            void Smear(int itr_nr, Bool_t measured);
            void Smear_tmp(std::vector<double> unc , int particle);

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

            bool isCB;

        private:
            static std::default_random_engine generator;

        };

    APLCON kinfit4g;
    APLCON kinfit;
    APLCON kinfit_eta;
    APLCON kinfit3pi;
    APLCON kinfiteta2pi;
    APLCON kinfit10g;

    FitParticle beam4g;
    FitParticle beam;
    FitParticle beam_3pi;
    FitParticle beam_eta2pi;
    FitParticle beam_eta;
    FitParticle beam10g;
    std::vector<FitParticle> Photons_four;
    std::vector<FitParticle> Photons_six;
    std::vector<FitParticle> Photons_six_3pi;
    std::vector<FitParticle> Photons_six_eta2pi;
    std::vector<FitParticle> Photons_six_eta;
    std::vector<FitParticle> Photons_ten;
    FitParticle proton4g;
    FitParticle proton;
    FitParticle proton_3pi;
    FitParticle proton_eta2pi;
    FitParticle proton_eta;
    FitParticle proton10g;
			
public:
    AdlarsonPhysics();
    virtual ~AdlarsonPhysics();

    virtual Bool_t	Init(const char* configfile);
    void            Energy_corr_MC();      // corrects theta for CB and TAPS for all clusters (Tracks).
    void            Energy_corr_EXP();      // corrects theta for CB and TAPS for all clusters (Tracks).
    void            theta_corr();      // corrects theta for CB and TAPS for all clusters (Tracks).

    // calculates IM for n photons
    Double_t    IM_Ng( UInt_t n );

    // function where true analysis is done for eta' --> eta 2pi0 --> 6g
    void TrueAnalysis_etapr6g();
    void Kinfit_test();

    std::vector<double> Get_unc(Int_t apparatus_nr, Int_t particle, std::vector<double>& obs);

    // functions specifically related to 4g analysis
    void fourgAnalysis( UInt_t ipr );
    void Best4g_comb(std::vector<TLorentzVector>& photons_rec, std::vector<int> &detnr , std::vector<bool> &CB_region);

    // functions specifically related to 6g analysis
    void sixgAnalysis( UInt_t ipr );
    void GetBest6gCombination(Double_t& sigma_eta, Double_t& sigma_pi0, Double_t& chi2min_eta2pi, Double_t& chi2min_3pi, std::vector<int>& imin_eta2pi, std::vector<int>& imin_3pi );
    void test_correct_hypothesis(Double_t& chi2min_eta2pi, Double_t& chi2min_3pi, std::vector<Int_t>& set_min, std::vector<int>& imin_eta2pi, std::vector<int>& imin_3pi);


    void sevengAnalysis( UInt_t ipr );

    // functions specifically related to 10g analysis
    void tengAnalysis(UInt_t ipr );
    void GetBest10gCombination(Double_t& sigma_eta, Double_t& sigma_pi0, Double_t& chi2min_eta3pi, std::vector<int>& imin_eta3pi2pi );

    const std::vector<TLorentzVector> ClusterEnergyCorr();


    void DalitzPlot(const TLorentzVector g[3] , Double_t &X1, Double_t &X2, Double_t &Y, Int_t &DP_nr1, Int_t &DP_nr2);
    void m2pi0_metapi0(  TLorentzVector g[3], Double_t &m_etapi01, Double_t &m_etapi02, Double_t &m_2pi0 );

};
#endif

