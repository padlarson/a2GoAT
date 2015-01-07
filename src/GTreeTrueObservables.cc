#include "GTreeTrueObservables.h"
#include "GTreeManager.h"


using namespace std;

GTreeTrueObservables::GTreeTrueObservables(GTreeManager *Manager) :
    GTree(Manager,  TString("trueObservables")),
    nTrueParticles(0),
    weight(1.0)
{
    for(int i=0; i<GTreeTrueObservables_MaxEntries; i++)
    {
        TrueEnergy[i]   = 0;
        TrueTheta[i]    = 0;
        TruePhi[i]      = 0;
        TrueXvertex[i]  = 0;
        TrueYvertex[i]  = 0;
        TrueZvertex[i]  = 0;
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
    inputTree->SetBranchAddress("TrueXvertex", &TrueXvertex);
    inputTree->SetBranchAddress("TrueYvertex", &TrueYvertex);
    inputTree->SetBranchAddress("TrueZvertex", &TrueZvertex);
    inputTree->SetBranchAddress("weight", &weight);
}

void    GTreeTrueObservables::SetBranches()
{
    outputTree->Branch("nTrueParticles",&nTrueParticles, "nTrueParticles/i");
    outputTree->Branch("TrueEnergy", &TrueEnergy, "TrueEnergy[nTrueParticles]/D");
    outputTree->Branch("TrueTheta", &TrueTheta, "TrueTheta[nTrueParticles]/D");
    outputTree->Branch("TruePhi", &TruePhi, "TruePhi[nTrueParticles]/D");
    outputTree->Branch("TrueXvertex", &TrueXvertex, "TrueXvertex[nTrueParticles]/D");
    outputTree->Branch("TrueYvertex", &TrueYvertex, "TrueYvertex[nTrueParticles]/D");
    outputTree->Branch("TrueZvertex", &TrueZvertex, "TrueZvertex[nTrueParticles]/D");
    outputTree->Branch("weight", &weight, "weight/D");
}



//void GTreeTrueObservables:: ProcessEvent()
//{}


