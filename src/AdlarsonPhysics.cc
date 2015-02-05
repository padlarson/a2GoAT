#include "GTreeManager.h"
#include "AdlarsonPhysics.h"
#include "GTrue.h"
#include "GTreeA2Geant.h"



AdlarsonPhysics::AdlarsonPhysics()
{

    TrueE               = new GH1("TrueE", "E_{true}", 200, 0., 1000);
    TrueTh              = new GH1("TrueTh", "#theta_{true}", 360, 0., 180);
    TruePhi             = new GH1("TruePhi", "#phi_{true}", 400, -200., 200.);
    TrueEvTh            = new GHistBGSub2("TrueEvTh", "E_{true} vs #theta_{true}", 100, 0., 1000, 60, 0., 30.);

    dthvdetnr_TAPS      = new GHistBGSub2("dthvdetnr_TAPS", "#Delta#theta vs detector nr TAPS", 80, -20., 20, 500, 0, 500);
    dfivdetnr_TAPS      = new GHistBGSub2("dfivdetnr_TAPS", "#Delta#phi vs detector nr TAPS", 80, -20., 20, 500, 0, 500);
    revdetnr_TAPS       = new GHistBGSub2("revdetnr_TAPS", "#Delta E /E_{rec} vs detector nr TAPS", 100, -2., 2, 500, 0, 500);

    dthvdetnr_CB        = new GHistBGSub2("dthvdetnr_CB", "#Delta#theta vs detector nr CB", 20, -20., 20, 800, 0, 800);
    dfivdetnr_CB        = new GHistBGSub2("dfivdetnr_CB", "#Delta#phi vs detector nr CB", 20, -20., 20, 800, 0, 800);
    revdetnr_CB         = new GHistBGSub2("revdetnr_CB", "#Delta E/E_{rec} vs detector nr CB", 100, -2., 2., 800, 0., 800.);

    Ncl_v_trueEtrueth   = new GHistBGSub2("revdetnr_CB", "#Delta E/E_{rec} vs detector nr CB", 100, -2., 2, 800, 0, 800);


     // Rec. TAPS - proton analysis
     EvdE_TAPS      = new GHistBGSub2("EvdE_TAPS", "E vs dE", 100, 0., 1000, 100, 0., 20.);



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

    th_true = (GetGeant()->GetTrueVector(0)).Theta()*TMath::RadToDeg();
    fi_true = (GetGeant()->GetTrueVector(0)).Phi()*TMath::RadToDeg();
    E_true  = (GetGeant()->GetTrueVector(0)).E()*1000 - MASS_PROTON;
//    E_true  = (GetGeant()->GetTrueVector(0)).E()*1000;

    TrueTh->Fill(th_true);
    TruePhi->Fill(fi_true);
    TrueE->Fill(E_true);
    TrueEvTh->Fill(E_true,th_true);

    for (Int_t i = 0; i< GetTracks()->GetNTracks(); i++)
    {

        if( GetTracks()->GetApparatus(i) == GTreeTrack::APPARATUS_CB )
        {
             E_rec  = GetTracks()->GetClusterEnergy(i);
             th_rec = GetTracks()->GetTheta(i);
             fi_rec = GetTracks()->GetPhi(i);

             deltath = th_true - th_rec;
             deltafi = fi_true - fi_rec;
             if(E_rec > 0.0)
             REnergy = (E_true-E_rec)/E_rec;

             Double_t itest = GetTracks()->GetCentralCrystal(i);

             dthvdetnr_CB->Fill(deltath, GetTracks()->GetCentralCrystal(i));
             dfivdetnr_CB->Fill(deltafi, GetTracks()->GetCentralCrystal(i));
             revdetnr_CB->Fill(REnergy, GetTracks()->GetCentralCrystal(i));
        }

        if( GetTracks()->GetApparatus(i) == GTreeTrack::APPARATUS_TAPS )
        {
             E_rec  = GetTracks()->GetClusterEnergy(i);
             th_rec = GetTracks()->GetTheta(i);
             fi_rec = GetTracks()->GetPhi(i);

             deltath = th_true - th_rec;
             deltafi = fi_true - fi_rec;
             if(E_rec > 0.0)
             REnergy = (E_true-E_rec)/E_rec;

             dthvdetnr_TAPS->Fill(deltath, GetTracks()->GetCentralCrystal(i));
             dfivdetnr_TAPS->Fill(deltafi, GetTracks()->GetCentralCrystal(i));
             revdetnr_TAPS->Fill(REnergy, GetTracks()->GetCentralCrystal(i));

             EvdE_TAPS->Fill(E_rec, GetTracks()->GetVetoEnergy(i));
        }

    }

   // for int

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






