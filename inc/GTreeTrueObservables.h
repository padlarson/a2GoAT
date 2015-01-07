#ifndef GTREETRUEOBSERVABLES_H
#define GTREETRUEOBSERVABLES_H

#include "GTree.h"

#define MASS_PI0    134.9766
#define MASS_ETA    547.853
#define MASS_ETAP   957.78
#define MASS_PROTON 938.272046

#define GTreeTrueObservables_MaxEntries  32

using namespace std;

class  GTreeTrueObservables : public GTree
{
private:
    UInt_t              nTrueParticles;
    // Energy, theta and phi for each particle
    Double_t            TrueEnergy[GTreeTrueObservables_MaxEntries];
    Double_t            TrueTheta[GTreeTrueObservables_MaxEntries];
    Double_t            TruePhi[GTreeTrueObservables_MaxEntries];
    // True vertex
    Double_t            TrueXvertex[GTreeTrueObservables_MaxEntries];
    Double_t            TrueYvertex[GTreeTrueObservables_MaxEntries];
    Double_t            TrueZvertex[GTreeTrueObservables_MaxEntries];
    // weight of each event, default weight = 1
    Double_t            weight;

    // X, Y, m2pi0, tagger, angular distribution i.e. analysis specific
protected:
    virtual void    SetBranchAdresses();
    virtual void    SetBranches();
//    virtual void    ProcessEvent();

public:
    GTreeTrueObservables(GTreeManager *Manager);
    virtual ~GTreeTrueObservables();

//    virtual void            Clear();
};



#endif // GTREETRUEOBSERVABLES_H

