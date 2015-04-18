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
    GH1*            True_BeamEnergy;
    GHistBGSub2*    ThpvsThetaprCM;
    GHistBGSub2*    ThvE_p;
    GHistBGSub2*    ThvE_eta_g;
    GHistBGSub2*    ThvE_pi0_g;
    GHistBGSub2*    DP_true;
    GH1*            M_pi1pi2_true;
    GH1*            M_etapi_true;

    // Tagger related
    GH1*            Tagged_BeamEnergy;

    // Proton related
    TLorentzVector  MMp_vec;

    GHistBGSub*     Nrprotons;
    GHistBGSub*     MM_p;
    GHistBGSub2*    ThvEp_rec;


    // Photon related
    TLorentzVector  IM6g_vec;
    GH1*            IM_6g;
    GHistBGSub2*    IM6gvMMp;

    TLorentzVector  IM10g_vec;
    GH1*            IM_10g;
    GHistBGSub2*    IM10gvMMp;

    GH1*            fi_diff_TAPSCB;
    GHistBGSub2*    fi_TAPSvsCB;
    GHistBGSub2*    fi_th_diff_TAPSCB;

    GHistBGSub2*    EvdE_TAPS_all;
    GHistBGSub2*    EvdE_TAPS_proton;
    GHistBGSub2*    EvTOF;
    GHistBGSub2*    EvTOFAll;
    GHistBGSub2*    EvTOFAllVeto;

    TH1D    RootHist;
    TH2F    ClustersTAPSTime;
    TH2F    ClustersCBTime;
    TH1D    AllClusters;
    TH1D    ClustersinTime;

    // Kinfit related variables

    GH1*            kfit_chi2;
    GH1*            kfit_pdf;
    GHistBGSub2*    kfit_Pulls;

    GH1*            IM6g_fit;

    GHistBGSub2*    PDF_eta2pi_v_3pi;
    GH1*            IM6g_fit_3pi;
    GH1*            IM6g_fit_eta2pi;


    GH1*            best_eta;
    GH1*            best_2pi;
    GHistBGSub2*    best_eta_EvTh;

    GH1*            M_pi1pi2_fit;
    GH1*            M_etapi_fit;
    GHistBGSub2*    deltaMpipi_v_Mpipi_fit;
    GHistBGSub2*    etapr_v_BeamE;

    GHistBGSub2*    DP_fit;
    GHistBGSub2*    deltaX_v_DPbin;
    GHistBGSub2*    deltaY_v_DPbin;


    // proton identified from TAPS_E vs VETO_dE

    TFile*          cutFile;
    TCutG*          cutProtonTAPS;
    TCutG*          cutPionTAPS;
    TCutG*          cutElectronTAPS;

    // True LorentzVectors
    TLorentzVector  eta_true;
    TLorentzVector  pi01_true;
    TLorentzVector  pi02_true;
    TLorentzVector  etapr_true[3];

    Double_t    Xtrue, Ytrue;
    Int_t       DPnrTrue;
    Double_t    m_etapi01True, m_etapi02True, m_2pi0True;

    // Reconstructed Lorentz Vectors
    std::vector<Int_t> ClustersInTime;

    Int_t nrprotons;
    Int_t iprtrack;

    TLorentzVector proton_vec;
    Double_t MMp;


    TLorentzVector etapr_sixgam[6];
    TLorentzVector dir3pi0_sixgam[6];

    TLorentzVector etapr_tengam[10];

    static Int_t perm6g[15][6];
    static Int_t perm6outof7g[7][6];
    static Int_t perm6outof8g[28][6];
    static Int_t perm6outof10g[210][6];

    GTrue   etapr_6gTrue;

    std::vector<TLorentzVector> photons_rec;
    std::vector<TLorentzVector> photons_fit;

    Double_t    sigma_eta;
    Double_t    sigma_pi0;


protected:
 //   virtual Bool_t  Init(const char* configFile);
    virtual Bool_t  Start();

    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();


    // choose here what you want to do
    // please also provide GoAT trees with matching MC true information...
    static constexpr bool includeIMconstraint = false;
    static constexpr bool includeVertexFit = true;
    static constexpr size_t nPhotons = 6;

    const double IM = MASS_ETA;


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

            void Smear(int itr_nr, Bool_t measured );
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

    APLCON kinfit;
    FitParticle beam;
    std::vector<FitParticle> photons;
    FitParticle proton;
			
public:
    AdlarsonPhysics();
    virtual ~AdlarsonPhysics();

    virtual Bool_t	Init(const char* configfile);

    // calculates IM for n photons
    Double_t    IM_Ng( UInt_t n );

    // function where true analysis is done for eta' --> eta 2pi0 --> 6g
    void TrueAnalysis_etapr6g();
    void Kinfit_test();

    // functions specifically related to 6g analysis
    void sixgAnalysis( Int_t ipr );
    void GetBest6gCombination(Double_t& sigma_eta, Double_t& sigma_pi0, Double_t& chi2min_eta2pi, Double_t& chi2min_3pi, std::vector<int>& imin_eta2pi, std::vector<int>& imin_3pi );

    // functions specifically related to 10g analysis
    void tengAnalysis(Int_t ipr );
    void GetBest6gCombination10g(Double_t& sigma_eta, Double_t& chi2min_eta3pi, std::vector<int>& imin_eta3pi );



    void DalitzPlot( const TLorentzVector g[3] , Double_t &X, Double_t &Y, Int_t &DP_nr );
    void m2pi0_metapi0(  TLorentzVector g[3], Double_t &m_etapi01, Double_t &m_etapi02, Double_t &m_2pi0 );

};
#endif

