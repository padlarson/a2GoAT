#ifndef __TrueRecParticleGun_h__
#define __TrueRecParticleGun_h__

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


class	TrueRecParticleGun  : public GTreeManager
{
private:

    // photon
    TH2*        photon_true_Evth;
    TH2*        photon_true_Evth_CB;
    TH2*        photon_true_Evth_TAPS;
    TH2*        photon_true_Evth_rec;

    TH3*        photon_secondary_Evth;
    TH3*        proton_secondary_Evth;

    TH3*        photon_rE_v_2D_CB;
    TH3*        photon_rE_v_2D_TAPS;
    TH3*        photon_dtheta_v_2D_CB;
    TH3*        photon_dtheta_v_2D_TAPS;
    TH3*        photon_dphi_v_2D_CB;
    TH3*        photon_dphi_v_2D_TAPS;

    TH3*        photon_rE_v_det_CB;
    TH3*        photon_dtheta_v_det_CB;
    TH3*        photon_dphi_v_det_CB;
    TH3*        photon_rE_v_det_TAPS;
    TH3*        photon_dtheta_v_det_TAPS;
    TH3*        photon_dphi_v_det_TAPS;

    TH3*        photon_EtEr_v_det_CB;  // Etrue/Erec CB
    TH3*        photon_EtEr_v_det_TAPS;// Etrue/Erec TAPS


    //proton
    TH2*        proton_true_Evth;
    TH3*        proton_rE_v_2D_TAPS;
    TH3*        proton_dtheta_v_2D_TAPS;
    TH3*        proton_dphi_v_2D_TAPS;

    TH3*        proton_rE_v_det_TAPS;
    TH3*        proton_dtheta_v_det_TAPS;
    TH3*        proton_dphi_v_det_TAPS;

    TH2*        proton_dth_v_detnr;
    TH2*        proton_dphi_v_detnr;
    TH2*        proton_re_v_theta;

    TH2*        photon_dtheta_v_theta_CB;
    TH2*        photon_dtheta_v_theta_TAPS;
    TH2*        proton_dtheta_v_theta_TAPS;


    TH2*        proton_EvDE_TAPS;
    TH1*        z_vertex;

    TFile*          thcorr_CB;        // File which contains TProfile
    TProfile*       dthvth_CB;
    TFile*          thcorr_TAPS;        // File which contains TProfile
    TProfile*       dthvth_TAPS;


protected:
    virtual Bool_t  Start();

    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();
			
public:
    TrueRecParticleGun();
    virtual ~TrueRecParticleGun();

    virtual Bool_t	Init(const char* configfile);

    void            theta_corr();



};
#endif

