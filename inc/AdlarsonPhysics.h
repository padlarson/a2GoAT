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
    // histograms and scatterplots
   GH1* True_BeamEnergy;
   // Phase space observables
   GHistBGSub2* ThvT_p;
   GHistBGSub2* ThpvEg;

   // Where all true observables are stored
   GTrue TrueObs;

   Double_t BeamE;



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
