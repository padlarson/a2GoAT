#ifndef __GTreeTrueParticle_h__
#define __GTreeTrueParticle_h__

#include <TClonesArray.h>
#include <TLorentzVector.h>

#include "GTree.h"

#define MASS_PI0    134.9766
#define MASS_ETA    547.853
#define MASS_ETAP   957.78
#define MASS_PROTON 938.272046

#define GTreeTrueParticle_MaxEntries  64

using namespace std;

class  GTreeTrueParticle : public GTree
{
private:

    UInt_t              nTrueParticles;
    TClonesArray*       TrueParticles;
    // weight of each event, default weight = 1
    Double_t            weight;
    // True vertex
    Double_t            XVertex[GTreeTrueParticle_MaxEntries];
    Double_t            YVertex[GTreeTrueParticle_MaxEntries];
    Double_t            ZVertex[GTreeTrueParticle_MaxEntries];



    // X, Y, m2pi0, tagger, angular distribution i.e. analysis specific
protected:
    virtual void    SetBranchAdresses();
    virtual void    SetBranches();

public:
    GTreeTrueParticle(GTreeManager *Manager, const TString& _Name);
    virtual         ~GTreeTrueParticle();
    void            AddTrueParticle(const TLorentzVector& truevec, const Double_t _XVertex, const Double_t _YVertex, const Double_t _ZVertex);
    virtual void    Clear() {nTrueParticles = 0; TrueParticles->Clear();}

    UInt_t          GetNTrueParticles()    const   {return nTrueParticles;}
 /*   Double_t        GetTrueEnergy(const Int_t index)    const	{return TrueEnergy[index];}
    Double_t        GetTrueTheta(const Int_t index)    const	{return TrueTheta[index];}
    Double_t        GetTruePhi(const Int_t index)    const	{return TruePhi[index];}
    // Vertex
    Double_t        GetTrueXVertex()     const   {return TrueXVertex;}
    Double_t        GetTrueYVertex()     const   {return TrueYVertex;}
    Double_t        GetTrueZVertex()     const   {return TrueZVertex;}
*/

    void        SetWeight(const Double_t w )   { weight = w;}
    Double_t    GetWeight() const   { return weight;}
};



#endif // __GTreeTrueParticle_h__

