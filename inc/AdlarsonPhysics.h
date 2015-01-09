#ifndef __AdlarsonPhysics_h__
#define __AdlarsonPhysics_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string> 

#include "GTreeManager.h"
#include "GH1.h"
#include "TVector.h"

class	AdlarsonPhysics  : public GTreeManager
{
private:
    GH1*	IM_6g;
    GH1*	IM_10g;

    TLorentzVector etapr_sixgam[6];
    TLorentzVector dir3pi0_sixgam[6];

    TLorentzVector etapr_tengam[10];

    static Int_t perm6g[15][6];
    static Int_t perm6outof7g[7][6];
    static Int_t perm6outof8g[28][6];
    static Int_t perm6outof10g[210][6];

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

    // functions specifically related to 6g analysis
    void sixgAnalysis();
    void GetBest6gCombination();

    // functions specifically related to 10g analysis
    void tengAnalysis();
    void GetBest10gCombination();


    void DalitzPlot( const TLorentzVector *g , Double_t &X, Double_t &Y, Int_t &DP_nr );
    void m2pi0_metapi0( const TLorentzVector *g, Double_t &m_etapi01, Double_t &m_etapi02, Double_t &m_2pi0 );

};
#endif
