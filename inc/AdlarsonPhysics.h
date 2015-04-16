#ifndef __AdlarsonPhysics_h__
#define __AdlarsonPhysics_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string> 

#include <TCutG.h>
#include "GTreeManager.h"
#include "GHistBGSub2.h"
#include "TVector.h"
#include "GTrue.h"


class	AdlarsonPhysics  : public GTreeManager
{
private:
    // Histograms and scatterplots declared
    GH1*	IM_6g;
    GH1*	IM_10g;

    GH1*            True_BeamEnergy;
    GHistBGSub2*    ThpvsThetaprCM;
    GHistBGSub2*    ThvE_p;
    GHistBGSub2*    ThvE_eta_g;
    GHistBGSub2*    ThvE_pi0_g;
    GHistBGSub2*    DP_true;
    GH1*            M_pi1pi2_true;
    GH1*            M_etapi_true;

    GHistBGSub2*    EvdE_TAPS;

    // proton identified from TAPS_E vs VETO_dE
    TCutG*          OpenCutFile(Char_t* filename, Char_t* cutname);
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

    TLorentzVector etapr_sixgam[6];
    TLorentzVector dir3pi0_sixgam[6];

    TLorentzVector etapr_tengam[10];

    static Int_t perm6g[15][6];
    static Int_t perm6outof7g[7][6];
    static Int_t perm6outof8g[28][6];
    static Int_t perm6outof10g[210][6];

    GTrue   etapr_6gTrue;

protected:
    virtual Bool_t  Start();

    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();
			
public:
    AdlarsonPhysics();
    virtual ~AdlarsonPhysics();

    virtual Bool_t	Init(const char* configfile);

    // calculates IM for n photons
    Double_t    IM_Ng( UInt_t n );

    // function where true analysis is done for eta' --> eta 2pi0 --> 6g
    void TrueAnalysis_etapr6g();

    // functions specifically related to 6g analysis
    void sixgAnalysis();
    void GetBest6gCombination();

    // functions specifically related to 10g analysis
    void tengAnalysis();
    void GetBest10gCombination();


    void DalitzPlot( const TLorentzVector g[3] , Double_t &X, Double_t &Y, Int_t &DP_nr );
    void m2pi0_metapi0(  TLorentzVector g[3], Double_t &m_etapi01, Double_t &m_etapi02, Double_t &m_2pi0 );
};
#endif
