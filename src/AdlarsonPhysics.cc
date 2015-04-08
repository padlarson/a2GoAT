#include "GTreeManager.h"
#include "AdlarsonPhysics.h"
#include "GTrue.h"
#include "GTreeA2Geant.h"



AdlarsonPhysics::AdlarsonPhysics()
{
    // Beam Energy
     True_BeamEnergy = new GH1("True_BeamEnergy", "True Beam Energy", 200, 0., 0.800);
    // Phase space observables
    ThvT_p = new GHistBGSub2("ThvT_p","#theta_{p} vs Energy p", 200, 0., 400, 250, 0., 100.);
    ThpvEg = new GHistBGSub2("ThpvEg","#theta_{p} vs Beam Energy #gamma", 100, 100., 600, 200, 0., 100.);

    GHistBGSub::InitCuts(-20, 20, -55, -35);
    GHistBGSub::AddRandCut(35, 55);
}

AdlarsonPhysics::~AdlarsonPhysics()
{


}

Bool_t	AdlarsonPhysics::Start()
{
    if(!IsGoATFile())
    {
        cout << "ERROR: Input File is not a GoAT file." << endl;
        return kFALSE;
    }
    SetAsPhysicsFile();

    TraverseValidEvents();

	return kTRUE;
}

void	AdlarsonPhysics::ProcessEvent()
{

    TrueObs.Start(*GetPluto(), *GetGeant());

    True_BeamEnergy->Fill( TrueObs.GetTrueBeamEnergy() );

    Double_t protonE = TrueObs.GetTrueProtonLV().E();

    ThvT_p->Fill( TrueObs.GetTrueProtonLV().E()*1000 - MASS_PROTON, TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg() );
    ThpvEg->Fill( TrueObs.GetTrueBeamEnergy()*1000, TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg()   );


}

void	AdlarsonPhysics::ProcessScalerRead()
{
//    hist_eta.ScalerReadCorrection(Double_t(scalers->GetScaler(0))/scalers->GetScaler(1));
//    hist_etap.ScalerReadCorrection(Double_t(scalers->GetScaler(0))/scalers->GetScaler(1));
}


Bool_t	AdlarsonPhysics::Init(const char* configfile)
{
    return kTRUE;
}






