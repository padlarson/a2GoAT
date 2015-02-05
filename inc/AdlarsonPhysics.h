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

    GH1*            TrueE;
    GH1*            TrueTh;
    GH1*            TruePhi;

    GHistBGSub2*    TrueEvTh;

    GHistBGSub2*    dthvdetnr_TAPS;
    GHistBGSub2*    dfivdetnr_TAPS;
    GHistBGSub2*    revdetnr_TAPS;

    GHistBGSub2*    dthvdetnr_CB;
    GHistBGSub2*    dfivdetnr_CB;
    GHistBGSub2*    revdetnr_CB;

    GHistBGSub2*    Ncl_v_trueEtrueth;
    GHistBGSub2*    EvdE_TAPS;


    double_t        th_true;
    double_t        fi_true;
    double_t        E_true;
    double_t        th_rec;
    double_t        fi_rec;
    double_t        E_rec;
    double_t        deltath;
    double_t        deltafi;
    double_t        REnergy;



protected:
    virtual Bool_t  Start();

    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();
			
public:
    AdlarsonPhysics();
    virtual ~AdlarsonPhysics();

    virtual Bool_t	Init(const char* configfile);


};
#endif
