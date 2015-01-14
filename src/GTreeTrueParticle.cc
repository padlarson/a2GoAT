#include "GTreeTrueParticle.h"
#include "GTreeManager.h"


using namespace std;

GTreeTrueParticle::GTreeTrueParticle(GTreeManager *Manager, const TString& _Name)    :
    GTree(Manager,_Name),
    nTrueParticles(0),
    weight(1.0),
    TrueParticles(new TClonesArray("TLorentzVector", GTreeTrueParticle_MaxEntries))
{
    for(int i=0; i<GTreeTrueParticle_MaxEntries; i++)
    {
        // True vertex
        XVertex[i]      = 0;
        YVertex[i]      = 0;
        ZVertex[i]      = 0;
    }
}



GTreeTrueParticle:: ~GTreeTrueParticle()
{

}

void    GTreeTrueParticle::SetBranchAdresses()
{
    inputTree->SetBranchAddress("nTrueParticles",&nTrueParticles);
    inputTree->SetBranchAddress("XVertex", &XVertex);
    inputTree->SetBranchAddress("YVertex", &YVertex);
    inputTree->SetBranchAddress("ZVertex", &ZVertex);
    inputTree->SetBranchAddress("weight", &weight);
}

void    GTreeTrueParticle::SetBranches()
{
    outputTree->Branch("nTrueParticles",&nTrueParticles, "nTrueParticles/i");
    outputTree->Branch("XVertex", &XVertex, "XVertex/D");
    outputTree->Branch("YVertex", &YVertex, "YVertex/D");
    outputTree->Branch("ZVertex", &ZVertex, "ZVertex/D");
    outputTree->Branch("weight", &weight, "weight/D");
}

void    GTreeTrueParticle::AddTrueParticle(const TLorentzVector& truevec, const Double_t _XVertex, const Double_t _YVertex, const Double_t _ZVertex)
{
    XVertex[nTrueParticles]     = _XVertex;
    YVertex[nTrueParticles]     = _YVertex;
    ZVertex[nTrueParticles]     = _ZVertex;
    new((*TrueParticles)[nTrueParticles]) TLorentzVector(truevec);
    nTrueParticles++;
}


