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

   GH1* MC_weight;
   GH1* mll;

   // Where all true observables are stored
   GTrue TrueObs;

   Double_t BeamE;





protected:
    virtual Bool_t  Start();

    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();
			
public:
    ProtonChargeRadius();
    virtual ~ProtonChargeRadius();

    virtual Bool_t	Init(const char* configfile);

    void TrueAnalysis_ll(int reaction);


};
#endif

