#ifndef __GTreeTrueObservables_h__
#define __GTreeTrueObservables_h__

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
    Double_t            TrueXVertex;
    Double_t            TrueYVertex;
    Double_t            TrueZVertex;
    // weight of each event, default weight = 1
    Double_t            weight;


    // X, Y, m2pi0, tagger, angular distribution i.e. analysis specific
protected:
    virtual void    SetBranchAdresses();
    virtual void    SetBranches();

 //           Bool_t  ProcessEventWithoutFilling();
 //   virtual void    ProcessEvent();

public:
    GTreeTrueObservables(GTreeManager *Manager);
    virtual ~GTreeTrueObservables();

    virtual void            Clear() {nTrueParticles = 0;}
    UInt_t          GetNTrueParticles()    const   {return nTrueParticles;}
    Double_t        GetTrueEnergy(const Int_t index)    const	{return TrueEnergy[index];}
    Double_t        GetTrueTheta(const Int_t index)    const	{return TrueTheta[index];}
    Double_t        GetTruePhi(const Int_t index)    const	{return TruePhi[index];}
    // Vertex
    Double_t        GetTrueXVertex()     const   {return TrueXVertex;}
    Double_t        GetTrueYVertex()     const   {return TrueYVertex;}
    Double_t        GetTrueZVertex()     const   {return TrueZVertex;}


    void        SetWeight(const Double_t w )   { weight = w;}
    Double_t    GetWeight() const   { return weight;}
};



#endif // __GTreeTrueObservables_h__

