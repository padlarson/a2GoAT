#include "GTreeManager.h"
#include "ProtonChargeRadius.h"
#include "GTrue.h"
#include "GTreeA2Geant.h"
#include "TH1.h"



ProtonChargeRadius::ProtonChargeRadius()
{
    // Beam Energy
     True_BeamEnergy = new TH1F("True_BeamEnergy", "True Beam Energy", 200, 0., 0.800);
    // Phase space observables
    ThvT_p = new TH2D("ThvT_p","#theta_{p} vs Energy p", 400, 0., 400, 250, 0., 100.);
    ThvT_e = new TH2F("ThvT_e","#theta_{e} vs Energy e^{+}/e^{-}", 400, 0., 400, 180, 0., 180.);
    ThvT_mu = new TH2F("ThvT_mu","#theta_{#mu} vs Energy #mu^{+}/#mu^{-}", 200, 0., 400, 250, 0., 100.);

    ThpvEg = new TH2F("ThpvEg","#theta_{p} vs Beam Energy #gamma", 100, 100., 600, 200, 0., 100.);





    GHistBGSub::InitCuts(-20, 20, -55, -35);
    GHistBGSub::AddRandCut(35, 55);
}

ProtonChargeRadius::~ProtonChargeRadius()
{


}

Bool_t	ProtonChargeRadius::Start()
{
    if(!IsGoATFile())
    {
        cout << "ERROR: Input File is not a GoAT file." << endl;
        return kFALSE;
    }
    True_BeamEnergy->Reset();
    ThvT_p->Reset();
    ThvT_e->Reset();
    ThvT_mu->Reset();
    ThpvEg->Reset();
\

    SetAsPhysicsFile();

    TraverseValidEvents();

    outputFile->cd();
    gDirectory->mkdir("MC")->cd();

    True_BeamEnergy->Write();
    ThvT_p->SetStats(0111111);
    ThvT_p->Write();
    ThvT_e->Write();
    ThvT_mu->Write();
    ThpvEg->Write();



	return kTRUE;
}

void	ProtonChargeRadius::ProcessEvent()
{

    TrueObs.Start(*GetPluto(), *GetGeant());

    Double_t proton_test = 200;
    Double_t proton_test_theta = 50;
    Double_t weight = 7.1;

    True_BeamEnergy->Fill( TrueObs.GetTrueBeamEnergy(), weight );

    Double_t protonE = TrueObs.GetTrueProtonLV().E();

//    ThvT_p->Fill( TrueObs.GetTrueProtonLV().E()*1000 - MASS_PROTON, TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg(), weight );
    ThvT_p->Fill( proton_test, proton_test_theta, weight );
    ThpvEg->Fill( TrueObs.GetTrueBeamEnergy()*1000, TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg(), weight);

    ThvT_e->Fill( TrueObs.GetTrueElectronLV().E()*1000 - MASS_ELECTRON, TrueObs.GetTrueElectronLV().Theta()*TMath::RadToDeg() );
    ThvT_e->Fill( TrueObs.GetTruePositronLV().E()*1000 - MASS_ELECTRON, TrueObs.GetTruePositronLV().Theta()*TMath::RadToDeg() );

    ThvT_mu->Fill( TrueObs.GetTrueMuonNegLV().E()*1000 - MASS_MUON, TrueObs.GetTrueMuonNegLV().Theta()*TMath::RadToDeg() );
    ThvT_mu->Fill( TrueObs.GetTrueMuonPosLV().E()*1000 - MASS_MUON, TrueObs.GetTrueMuonPosLV().Theta()*TMath::RadToDeg() );


}

void	ProtonChargeRadius::ProcessScalerRead()
{
//    hist_eta.ScalerReadCorrection(Double_t(scalers->GetScaler(0))/scalers->GetScaler(1));
//    hist_etap.ScalerReadCorrection(Double_t(scalers->GetScaler(0))/scalers->GetScaler(1));
}


Bool_t	ProtonChargeRadius::Init(const char* configfile)
{
    return kTRUE;
}






