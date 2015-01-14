#include "GTreeTrueObservables.h"
#include "GTreeManager.h"


using namespace std;

GTreeTrueObservables::GTreeTrueObservables(GTreeManager *Manager) :
    GTree(Manager,  TString("trueObservables")),
    nTrueParticles(0),
    weight(1.0),
    TrueXVertex(0.0),
    TrueYVertex(0.0),
    TrueZVertex(0.0)

{
    for(int i=0; i<GTreeTrueObservables_MaxEntries; i++)
    {
        TrueEnergy[i]   = 0;
        TrueTheta[i]    = 0;
        TruePhi[i]      = 0;
    }

}

GTreeTrueObservables:: ~GTreeTrueObservables()
{

}

void    GTreeTrueObservables::SetBranchAdresses()
{
    inputTree->SetBranchAddress("nTrueParticles",&nTrueParticles);
    inputTree->SetBranchAddress("TrueEnergy",&TrueEnergy);
    inputTree->SetBranchAddress("TrueTheta", &TrueTheta);
    inputTree->SetBranchAddress("TruePhi", &TruePhi);
    inputTree->SetBranchAddress("TrueXVertex", &TrueXVertex);
    inputTree->SetBranchAddress("TrueYVertex", &TrueYVertex);
    inputTree->SetBranchAddress("TrueZVertex", &TrueZVertex);
    inputTree->SetBranchAddress("weight", &weight);
}

void    GTreeTrueObservables::SetBranches()
{
    outputTree->Branch("nTrueParticles",&nTrueParticles, "nTrueParticles/i");
    outputTree->Branch("TrueEnergy", &TrueEnergy, "TrueEnergy[nTrueParticles]/D");
    outputTree->Branch("TrueTheta", &TrueTheta, "TrueTheta[nTrueParticles]/D");
    outputTree->Branch("TruePhi", &TruePhi, "TruePhi[nTrueParticles]/D");
    outputTree->Branch("TrueXvertex", &TrueXVertex, "TrueXVertex/D");
    outputTree->Branch("TrueYvertex", &TrueYVertex, "TrueYVertex/D");
    outputTree->Branch("TrueZvertex", &TrueZVertex, "TrueZVertex/D");
    outputTree->Branch("weight", &weight, "weight/D");
}

/*
Bool_t  GTreeTrueObservables::ProcessEventWithoutFilling()
{
    cout << "test" << endl;
    return kTRUE;
}

void GTreeTrueObservables::ProcessEvent()
{

}
*/

