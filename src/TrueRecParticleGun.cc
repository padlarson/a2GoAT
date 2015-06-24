#include "GTreeManager.h"
#include "TrueRecParticleGun.h"
#include "GTrue.h"
#include "GTreeA2Geant.h"
#include "TH1.h"



TrueRecParticleGun::TrueRecParticleGun()
{
    //Relative Energy
    proton_true_Evth        = new TH2F("proton_true_Evth", "proton; E_{true} (MeV); #theta_{true} (^{o})", 30, 0., 600., 30, 0, 30);
    photon_true_Evth        = new TH2F("photon_true_Evth", "#gamma; E_{true} (MeV); #theta_{true} (^{o})", 75, 0., 1500., 36, 0, 180);
    photon_true_Evth_CB     = new TH2F("photon_true_Evth_CB", "#gamma; E_{true} (MeV); #theta_{true} (^{o})", 75, 0., 1500., 36, 0, 180);
    photon_true_Evth_TAPS   = new TH2F("photon_true_Evth_TAPS", "#gamma; E_{true} (MeV); #theta_{true} (^{o})",  75, 0., 1500., 30, 0, 30);
    photon_true_Evth_rec    = new TH2F("photon_true_Evth_rec", "#gamma nr tracks ; E_{true} (MeV); #theta_{true} (^{o})", 75, 0., 1500., 36, 0, 180);

    photon_secondary_Evth   = new TH3F("photon_secondary_Evth", "Split-off photon; E_{rec} (MeV); #theta_{rec} (^{o}); E_{true} MeV", 150, 0., 1500, 90, 0., 180., 120, 0., 1200.);
    proton_secondary_Evth   = new TH3F("proton_secondary_Evth", "Split-off proton; E_{rec} (MeV); #theta_{rec} (^{o}); E_{true} MeV", 150, 0., 1500, 40, 0., 40., 60, 0., 600.);

    photon_rE_v_2D_CB       = new TH3F("photon_rE_v_2D_CB", "#gamma CB; E_{rec} (MeV); #theta_{rec} (^{o}); relative (E_{rec} - E_{true}) / E_{rec} ", 75, 0., 1500., 36, 0, 180, 200, -5., 5.);
    photon_dtheta_v_2D_CB   = new TH3F("photon_dtheta_v_2D_CB", "#gamma CB; E_{rec} (MeV); #theta_{rec} (^{o});  #theta_{rec} - #theta_{true} (^{o})", 75, 0., 1500., 36, 0, 180, 20, -40., 40.);
    photon_dphi_v_2D_CB     = new TH3F("photon_dphi_v_2D_CB", "#gamma CB; E_{rec} (MeV); #theta_{rec} (^{o});  #phi_{rec} - #phi_{true} (^{o})", 75, 0., 1500., 36, 0, 180, 800, -40., 40.);

    photon_rE_v_2D_TAPS     = new TH3F("photon_rE_v_2D_TAPS", "#gamma TAPS; E_{rec} (MeV); #theta_{rec} (^{o}); relative (E_{rec} - E_{true}) / E_{rec} ", 75, 0., 1500., 30, 0, 30, 500, -5., 5.);
    proton_rE_v_2D_TAPS     = new TH3F("proton_rE_v_2D_TAPS", "proton TAPS; E_{rec} (MeV); #theta_{rec} (^{o}); relative (E_{rec} - E_{true}) / E_{rec} ", 25, 0., 500., 30, 0, 30, 500, -5., 5.);
    photon_dtheta_v_2D_TAPS = new TH3F("photon_dtheta_v_2D_TAPS", "#gamma TAPS; E_{rec} (MeV); #theta_{rec} (^{o}); #theta_{rec} - #theta_{true} (^{o})", 75, 0., 1500., 30, 0, 30, 160, -20., 20.);
    proton_dtheta_v_2D_TAPS = new TH3F("proton_dtheta_v_2D_TAPS", "proton TAPS; E_{rec} (MeV); #theta_{rec} (^{o}); #theta_{rec} - #theta_{true} (^{o})", 75, 0., 1500., 30, 0, 30, 160, -20., 20.);
    photon_dphi_v_2D_TAPS   = new TH3F("photon_dphi_v_2D_TAPS", "#gamma TAPS; E_{rec} (MeV); #theta_{rec} (^{o}); #phi_{rec} - #phi_{true} (^{o})", 75, 0., 1500., 30, 0, 30, 160, -40., 40.);
    proton_dphi_v_2D_TAPS   = new TH3F("proton_dphi_v_2D_TAPS", "proton TAPS; E_{rec} (MeV); #theta_{rec} (^{o}); #phi_{rec} - #phi_{true} (^{o})", 150, 0., 1500., 60, 0, 30, 160, -40., 40.);

    proton_dth_v_detnr      = new TH2F("proton_dth_v_detnr", "proton; #Delta #theta ; TAPS element nr ", 160, -40., 40., 440, 0, 440);
    proton_dphi_v_detnr     = new TH2F("proton_dphi_v_detnr", "proton; #Delta #phi ; TAPS element nr ", 160, -40., 40., 440, 0, 440);

    proton_re_v_theta     = new TH2F("proton_re_v_theta", "proton; rel err ; #theta", 60, 0., 30., 500, -5., 5.);


    proton_EvDE_TAPS        = new TH2F("proton_EvDE_TAPS", "proton; E (MeV); #Delta E (MeV)", 1000, 0., 1000., 400, 0., 20.);

    GHistBGSub::InitCuts(-20, 20, -55, -35);
    GHistBGSub::AddRandCut(35, 55);
}

TrueRecParticleGun::~TrueRecParticleGun()
{

}
Bool_t	TrueRecParticleGun::Start()
{
//    if(!IsGoATFile())
//    {
//        cout << "ERROR: Input File is not a GoAT file." << endl;
//        return kFALSE;
//    }

    SetAsPhysicsFile();

    TraverseValidEvents();

    outputFile->cd();

	return kTRUE;
}

void	TrueRecParticleGun::ProcessEvent()
{

//    TrueObs.Start(*GetGeant());
    TLorentzVector trueobs;
    trueobs = GetGeant()->GetTrueVector(0);

    Double_t E_true = trueobs.E()*1000;
    Double_t th_true = trueobs.Theta()*TMath::RadToDeg();
    Double_t fi_true = trueobs.Phi()*TMath::RadToDeg();

    // particle:
    // 0  photon
    // 1  proton
    int particle = 1;

    if( particle == 0 )
    {
        photon_true_Evth->Fill(E_true, th_true);
        UInt_t i = 999;
        UInt_t k = 999;
        for( UInt_t j = 0; j < GetTracks()->GetNTracks(); j++ )
        {
            photon_true_Evth_rec->Fill(E_true, th_true);

            if( GetTracks()->GetNTracks() == 1)
                i = j;
            else if( GetTracks()->GetNTracks() == 2 )
            {
                Double_t dth0 = (th_true - GetTracks()->GetTheta(0));
                Double_t dfi0 = (fi_true - GetTracks()->GetPhi(0));
                Double_t chi0 = TMath::Power(dth0,2) + TMath::Power(dfi0,2);

                Double_t dth1 = (th_true - GetTracks()->GetTheta(1));
                Double_t dfi1 = (fi_true - GetTracks()->GetPhi(1));
                Double_t chi1 = TMath::Power(dth1,2) + TMath::Power(dfi1,2);
                if( chi0 < chi1 )
                {
                    i = 0;
                    k = 1;
                }
                else
                {
                    i = 1;
                    k = 0;
                }

            }

        }
        if( i < 2 )
        {
            Double_t E_rec    =  GetTracks()->GetClusterEnergy(i);
            Double_t th_rec   =  GetTracks()->GetTheta(i);

            Double_t rel_E  =  ( GetTracks()->GetClusterEnergy(i) - E_true )/(GetTracks()->GetClusterEnergy(i));
            Double_t dth    =  ( GetTracks()->GetTheta(i) - th_true );
            Double_t dfi    =  ( GetTracks()->GetPhi(i) - fi_true );
            if( k < 2 ) // plot the split-off E_rec and angle vs E_true (theta assumed at CB edge)
            {
                photon_secondary_Evth->Fill( GetTracks()->GetClusterEnergy(k), GetTracks()->GetTheta(k), E_true );
            }
            if(GetTracks()->HasCB(i))
            {
                photon_true_Evth_CB->Fill(E_true, th_true);
                photon_rE_v_2D_CB->Fill(E_rec, th_rec, rel_E);
                photon_dtheta_v_2D_CB->Fill(E_rec, th_rec, dth);
                photon_dphi_v_2D_CB->Fill(E_rec, th_rec, dfi);
            }
            else if(GetTracks()->HasTAPS(i)) //TAPS
            {
                photon_true_Evth_TAPS->Fill(E_true, th_true);
                photon_rE_v_2D_TAPS->Fill(E_rec, th_rec, rel_E);
                photon_dtheta_v_2D_TAPS->Fill(E_rec, th_rec, dth);
                photon_dphi_v_2D_TAPS->Fill(E_rec, th_rec, dfi);
            }
        }
    }
    if( particle == 1 )
    {
        Double_t E_true = trueobs.E()*1000;
        Double_t th_true = trueobs.Theta()*TMath::RadToDeg();
        Double_t fi_true = trueobs.Phi()*TMath::RadToDeg();
        proton_true_Evth->Fill(E_true-MASS_PROTON, th_true);
        UInt_t i = 999;
        UInt_t k = 999;
        for( UInt_t j = 0; j < GetTracks()->GetNTracks(); j++ )
        {
            if( GetTracks()->GetNTracks() == 1)
                i = j;
            else if( GetTracks()->GetNTracks() == 2 )
            {
                Double_t dth0 = (th_true - GetTracks()->GetTheta(0));
                Double_t dfi0 = (fi_true - GetTracks()->GetPhi(0));
                Double_t chi0 = TMath::Power(dth0,2) + TMath::Power(dfi0,2);

                Double_t dth1 = (th_true - GetTracks()->GetTheta(1));
                Double_t dfi1 = (fi_true - GetTracks()->GetPhi(1));
                Double_t chi1 = TMath::Power(dth1,2) + TMath::Power(dfi1,2);
                if( chi0 < chi1 )
                {
                    i = 0;
                    k = 1;
                }
                else
                {
                    i = 1;
                    k = 0;
                }
            }
        }
        if( i < 2 )
        {
            Double_t E_rec  =  GetTracks()->GetClusterEnergy(i);
            Double_t th_rec =  GetTracks()->GetTheta(i);
            Double_t fi_rec =  GetTracks()->GetPhi(i);

            Double_t rel_E  =  ( GetTracks()->GetClusterEnergy(i) - (E_true-MASS_PROTON) )/(GetTracks()->GetClusterEnergy(i));
            Double_t dth    =  ( th_rec - th_true );

            Double_t dfi    =  ( fi_rec - fi_true );
            if( k < 2 ) // plot the split-off E_rec and angle vs E_true (theta assumed at CB edge)
            {
                proton_secondary_Evth->Fill( GetTracks()->GetClusterEnergy(k), GetTracks()->GetTheta(k), E_true-MASS_PROTON );
            }

            if(GetTracks()->HasTAPS(i))
            {
                proton_rE_v_2D_TAPS->Fill(E_rec, th_rec, rel_E);
                proton_dtheta_v_2D_TAPS->Fill(E_rec, th_rec, dth);
                proton_dphi_v_2D_TAPS->Fill(E_rec, th_rec, dfi);

                proton_re_v_theta->Fill(th_rec , rel_E);

                proton_dth_v_detnr->Fill(dth, GetTracks()->GetCentralCrystal(i));
                proton_dphi_v_detnr->Fill(dfi, GetTracks()->GetCentralCrystal(i));

                proton_EvDE_TAPS->Fill(GetTracks()->GetClusterEnergy(i), GetTracks()->GetVetoEnergy(i));

            }
        }
    }


//        E_v_dE->FillWeighted( GetTracks()->GetClusterEnergy(i), GetTracks()->GetVetoEnergy(i), TrueObs.GetWeight());


}

void	TrueRecParticleGun::ProcessScalerRead()
{
//    hist_eta.ScalerReadCorrection(Double_t(scalers->GetScaler(0))/scalers->GetScaler(1));
//    hist_etap.ScalerReadCorrection(Double_t(scalers->GetScaler(0))/scalers->GetScaler(1));
}


Bool_t	TrueRecParticleGun::Init(const char* configfile)
{
    return kTRUE;
}









