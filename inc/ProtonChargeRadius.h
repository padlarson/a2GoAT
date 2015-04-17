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
   TH1F* True_BeamEnergy;
   // Phase space observables
   TH2* ThvT_p;
   TH2* ThvT_e;
   TH2* ThvT_mu;

   TH2* ThpvEg;

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


};
#endif

