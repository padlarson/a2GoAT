#ifndef __ProtonChargeRadius_h__
#define __ProtonChargeRadius_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string> 

#include <TCutG.h>
#include "GTreeManager.h"
#include "GHistBGSub2.h"
#include "TVector.h"
#include "GTrue.h"

#define MASS_ELECTRON       0.510998928
#define MASS_MUON           105.6583715


class	ProtonChargeRadius  : public GTreeManager
{
private:
    // histograms and scatterplots
   GH1* True_BeamEnergy;
   // Phase space observables
   GHistBGSub2* ThvT_p;
   GHistBGSub2* ThvT_e;
   GHistBGSub2* ThvT_mu;

   GHistBGSub2* ThpvEg;
   GHistBGSub2* dsigma_v_mll;
   GHistBGSub2* t_v_mll;
   GHistBGSub2* thlabpr_v_mll;

   GHistBGSub2* thlabpr_v_t;

   GHistBGSub2* thlabpr_v_th_e_1;
   GHistBGSub2* thlabpr_v_th_mu_1;
   GHistBGSub2* thlabpr_v_Tp_1;
   GHistBGSub2* thlabpr_v_th_e_2;
   GHistBGSub2* thlabpr_v_th_mu_2;
   GHistBGSub2* thlabpr_v_Tp_2;
   GHistBGSub2* thlabpr_v_th_e_3;
   GHistBGSub2* thlabpr_v_th_mu_3;
   GHistBGSub2* thlabpr_v_Tp_3;
   GHistBGSub2* thlabpr_v_Tp_13;

   GH1*         Ntracks;

   GHistBGSub2* E_v_dE;
   GHistBGSub2* dt_vs_t_rec;
   GHistBGSub2* dth_p_v_th_p;
   GHistBGSub2* dE_p_v_E_p;
   GHistBGSub2* proton_dtheta_v_E;
   GHistBGSub2* proton_rE_v_E;
   GHistBGSub2* proton_rt_v_E;

   GH1* MC_weight;
   GH1* mll;

   GHistBGSub2*        proton_E_v_th;
   GHistBGSub2*        proton_t_v_th;
   GHistBGSub3*        proton_rE_v_2D;
   GHistBGSub3*        proton_rt_v_2D;
   GHistBGSub3*        proton_dtheta_v_2D;
   GHistBGSub3*        proton_dphi_v_2D;

   TFile*      cutFile;        // File which contains EdE cut
   TCutG*      cutProtonCB;

   // Where all true observables are stored
   GTrue TrueObs;

   Double_t BeamE;
   Double_t t_true;





protected:
    virtual Bool_t  Start();

    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();
			
public:
    ProtonChargeRadius();
    virtual ~ProtonChargeRadius();

    virtual Bool_t	Init(const char* configfile);

    void TrueAnalysis_ll(int reaction);
    Double_t GEp(Double_t Qsqr);
    Double_t GMp(Double_t Qsqr);



};
#endif

