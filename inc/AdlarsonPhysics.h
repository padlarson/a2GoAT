#ifndef __AdlarsonPhysics_h__
#define __AdlarsonPhysics_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string> 

#include "GTreeManager.h"
#include "GH1.h"

class	AdlarsonPhysics  : public GTreeManager
{
private:
    GH1*	IM_6g;
    GH1*	IM_10g;

    TLorentzVector sixgam_v;
    TLorentzVector tengam_v;

protected:
    virtual Bool_t  Start();

    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();
			
public:
    AdlarsonPhysics();
    virtual ~AdlarsonPhysics();

    virtual Bool_t	Init(const char* configfile);

    // calculates IM for n photons
    double_t    IM_Ng( UInt_t n );

    // method for 6g analysis
    void sixgAnalysis();
    void tengAnalysis();

    void PhysicsResults( const TLorentzVector &g );

};
#endif
