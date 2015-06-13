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

    TH3*        photon_rE_v_2D_CB;
    TH3*        photon_rE_v_2D_TAPS;
    TH3*        photon_dtheta_v_2D_CB;
    TH3*        photon_dtheta_v_2D_TAPS;
    TH3*        photon_dphi_v_2D_CB;
    TH3*        photon_dphi_v_2D_TAPS;

    //proton
    TH3*        proton_rE_v_2D_TAPS;
    TH3*        proton_dtheta_v_2D_TAPS;
    TH3*        proton_dphi_v_2D_TAPS;

    TH2*        proton_EvDE_TAPS;


protected:
    virtual Bool_t  Start();

    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();
			
public:
    TrueRecParticleGun();
    virtual ~TrueRecParticleGun();

    virtual Bool_t	Init(const char* configfile);



};
#endif

