#include "GTreeManager.h"
#include "TrueRecParticleGun.h"
#include "GTrue.h"
#include "GTreeA2Geant.h"
#include "TH1.h"
#include "TProfile.h"



TrueRecParticleGun::TrueRecParticleGun()
{
    //Relative Energy
    proton_true_Evth        = new TH2F("proton_true_Evth", "proton; E_{true} (MeV); #theta_{true} (^{o})", 30 , 0., 600., 30, 0, 30);
    photon_true_Evth        = new TH2F("photon_true_Evth", "#gamma; E_{true} (MeV); #theta_{true} (^{o})", 150, 0., 1500., 180, 0, 180);
    photon_true_Evth_CB     = new TH2F("photon_true_Evth_CB", "#gamma; E_{true} (MeV); #theta_{true} (^{o})", 150, 0., 1500., 90, 0, 180);
    photon_true_Evth_TAPS   = new TH2F("photon_true_Evth_TAPS", "#gamma; E_{true} (MeV); #theta_{true} (^{o})",  150, 0., 1500., 30, 0, 30);
    photon_true_Evth_rec    = new TH2F("photon_true_Evth_rec", "#gamma nr tracks ; E_{true} (MeV); #theta_{true} (^{o})", 75, 0., 1500., 90, 0, 180);

    photon_secondary_Evth   = new TH3F("photon_secondary_Evth", "Split-off photon; E_{rec} (MeV); #theta_{rec} (^{o}); E_{true} MeV", 150, 0., 1500, 90, 0., 180., 120, 0., 1200.);
    proton_secondary_Evth   = new TH3F("proton_secondary_Evth", "Split-off proton; E_{rec} (MeV); #theta_{rec} (^{o}); E_{true} MeV", 150, 0., 1500, 40, 0., 40., 60, 0., 600.);

    photon_rE_v_2D_CB       = new TH3F("photon_rE_v_2D_CB", "#gamma CB; E_{rec} (MeV); #theta_{rec} (^{o}); relative (E_{rec} - E_{true}) / E_{rec} ", 150, 0., 1500., 90, 0, 180, 400, -2., 2.);
    photon_dtheta_v_2D_CB   = new TH3F("photon_dtheta_v_2D_CB", "#gamma CB; E_{rec} (MeV); #theta_{rec} (^{o});  #theta_{rec} - #theta_{true} (^{o})", 150, 0., 1500., 90, 0, 180, 160, -20., 20.);

    photon_theta_phi_CB   = new TH2F("photon_theta_phi_CB", "#gamma CB; #theta_{rec} (MeV); #phi_{rec} (^{o});", 720, 0., 180., 720, 0., 180.);

    photon_dphi_v_2D_CB     = new TH3F("photon_dphi_v_2D_CB", "#gamma CB; E_{rec} (MeV); #theta_{rec} (^{o});  #phi_{rec} - #phi_{true} (^{o})", 150, 0., 1500., 90, 0, 180, 160, -20., 20.);

    photon_rE_v_2D       = new TH3F("photon_rE_v_2D", "#gamma; E_{rec} (MeV); #theta_{rec} (^{o}); relative (E_{rec} - E_{true}) / E_{rec} ", 150, 0., 1500., 180, 0, 180, 400, -2., 2.);
    photon_dtheta_v_2D   = new TH3F("photon_dtheta_v_2D", "#gamma; E_{rec} (MeV); #theta_{rec} (^{o});  #theta_{rec} - #theta_{true} (^{o})", 150, 0., 1500., 180, 0, 180, 160, -20., 20.);
    photon_dphi_v_2D     = new TH3F("photon_dphi_v_2D", "#gamma; E_{rec} (MeV); #theta_{rec} (^{o});  #phi_{rec} - #phi_{true} (^{o})", 150, 0., 1500., 180, 0, 180, 160, -20., 20.);

    photon_rE_v_det_CB      = new TH3F("photon_rE_v_det_CB", "#gamma CB; E_{rec} (MeV); det.nr; relative (E_{rec} - E_{true}) / E_{rec} ", 75, 0., 1500., 720, 0, 720, 200, -5., 5.);
    photon_dtheta_v_det_CB  = new TH3F("photon_dtheta_v_det_CB", "#gamma CB; E_{rec} (MeV); det.nr;  #theta_{rec} - #theta_{true} (^{o})", 75, 0., 1500., 720, 0, 720, 20, -40., 40.);
    photon_dphi_v_det_CB    = new TH3F("photon_dphi_v_det_CB", "#gamma CB; E_{rec} (MeV); det.nr;  #phi_{rec} - #phi_{true} (^{o})", 75, 0., 1500., 720, 0, 720, 800, -40., 40.);

    photon_rE_v_2D_TAPS     = new TH3F("photon_rE_v_2D_TAPS", "#gamma TAPS; E_{rec} (MeV); #theta_{rec} (^{o}); relative (E_{rec} - E_{true}) / E_{rec} ", 150, 0., 1500., 30, 0, 30, 500, -5., 5.);
    proton_rE_v_2D_TAPS     = new TH3F("proton_rE_v_2D_TAPS", "proton TAPS; E_{rec} (MeV); #theta_{rec} (^{o}); relative (E_{rec} - E_{true}) / E_{rec} ", 30, 0., 600., 30, 0, 30, 500, -5., 5.);
    photon_dtheta_v_2D_TAPS = new TH3F("photon_dtheta_v_2D_TAPS", "#gamma TAPS; E_{rec} (MeV); #theta_{rec} (^{o}); #theta_{rec} - #theta_{true} (^{o})", 150, 0., 1500., 30, 0, 30, 160, -20., 20.);
    proton_dtheta_v_2D_TAPS = new TH3F("proton_dtheta_v_2D_TAPS", "proton TAPS; E_{rec} (MeV); #theta_{rec} (^{o}); #theta_{rec} - #theta_{true} (^{o})", 30, 0., 600., 30, 0, 30, 160, -20., 20.);
    photon_dphi_v_2D_TAPS   = new TH3F("photon_dphi_v_2D_TAPS", "#gamma TAPS; E_{rec} (MeV); #theta_{rec} (^{o}); #phi_{rec} - #phi_{true} (^{o})", 150, 0., 1500., 30, 0, 30, 160, -40., 40.);
    proton_dphi_v_2D_TAPS   = new TH3F("proton_dphi_v_2D_TAPS", "proton TAPS; E_{rec} (MeV); #theta_{rec} (^{o}); #phi_{rec} - #phi_{true} (^{o})", 30, 0., 600., 30, 0, 30, 160, -40., 40.);

    photon_rE_v_det_TAPS     = new TH3F("photon_rE_v_det_TAPS", "#gamma TAPS; E_{rec} (MeV); #theta_{rec} (^{o}); relative (E_{rec} - E_{true}) / E_{rec} ", 150, 0., 1500., 440, 0, 440, 500, -5., 5.);
    proton_rE_v_det_TAPS     = new TH3F("proton_rE_v_det_TAPS", "proton TAPS; E_{rec} (MeV); #theta_{rec} (^{o}); relative (E_{rec} - E_{true}) / E_{rec} ", 30, 0., 600., 440, 0, 440, 500, -5., 5.);
    photon_dtheta_v_det_TAPS = new TH3F("photon_dtheta_v_det_TAPS", "#gamma TAPS; E_{rec} (MeV); #theta_{rec} (^{o}); #theta_{rec} - #theta_{true} (^{o})", 150, 0., 1500., 440, 0, 440, 160, -20., 20.);
    proton_dtheta_v_det_TAPS = new TH3F("proton_dtheta_v_det_TAPS", "proton TAPS; E_{rec} (MeV); #theta_{rec} (^{o}); #theta_{rec} - #theta_{true} (^{o})", 30, 0., 600., 440, 0, 440, 160, -20., 20.);
    photon_dR_v_det_TAPS     = new TH2F("photon_dR_v_det_TAPS", "photon_dR_v_det_TAPS", 440, 0, 440, 400, -20., 20.);


    photon_dphi_v_det_TAPS   = new TH3F("photon_dphi_v_det_TAPS", "#gamma TAPS; E_{rec} (MeV); #theta_{rec} (^{o}); #phi_{rec} - #phi_{true} (^{o})", 150, 0., 1500., 440, 0, 440, 160, -40., 40.);
    proton_dphi_v_det_TAPS   = new TH3F("proton_dphi_v_det_TAPS", "proton TAPS; E_{rec} (MeV); #theta_{rec} (^{o}); #phi_{rec} - #phi_{true} (^{o})", 30, 0., 600., 440, 0, 440, 160, -40., 40.);

    photon_EtEr_v_det_CB     = new TH3F("photon_EtEr_v_det_CB", "#gamma CB; E_{rec} (MeV); det nr CB; E_{true} / E_{rec} ", 120, 0., 1200., 720, 0, 720, 200, 0.5, 1.5);
    photon_EtEr_v_detmod_CB  = new TH3F("photon_EtEr_v_detmod_CB", "#gamma CB; E_{rec} (MeV); det nr CB; E_{true} / E_{rec} ", 1200, 0., 1200., 45, 0, 45, 400, 0.8, 1.2);
    photon_EtEr_v_det_TAPS   = new TH3F("photon_EtEr_v_det_TAPS", "#gamma TAPS; E_{rec} (MeV); det nr TAPS; E_{true} / E_{rec} ", 120, 0., 1200., 440, 0, 440, 200, 0.5, 1.5);

    proton_dth_v_detnr      = new TH2F("proton_dth_v_detnr", "proton; #Delta #theta ; TAPS element nr ", 160, -40., 40., 440, 0, 440);
    proton_dphi_v_detnr     = new TH2F("proton_dphi_v_detnr", "proton; #Delta #phi ; TAPS element nr ", 160, -40., 40., 440, 0, 440);

    proton_re_v_theta     = new TH2F("proton_re_v_theta", "proton; rel err ; #theta", 60, 0., 30., 500, -5., 5.);


    proton_EvDE_TAPS        = new TH2F("proton_EvDE_TAPS", "proton; E (MeV); #Delta E (MeV)", 1000, 0., 1000., 400, 0., 20.);

    photon_dtheta_v_theta_CB   = new TH2F("photon_dtheta_v_theta_CB", "#gamma CB; #theta_{rec} (^{o});  #theta_{rec} - #theta_{true} (^{o})", 360, 0, 180, 200, -20., 20.);
    photon_dtheta_v_theta_TAPS   = new TH2F("photon_dtheta_v_theta_TAPS", "#gamma TAPS; #theta_{rec} (^{o});  #theta_{rec} - #theta_{true} (^{o})", 150, 0, 30, 200, -10., 10.);
    proton_dtheta_v_theta_TAPS   = new TH2F("proton_dtheta_v_theta_TAPS", "proton TAPS; #theta_{rec} (^{o});  #theta_{rec} - #theta_{true} (^{o})", 60, 0, 30, 100, -10., 10.);

    photon_rE_v_E_CB            = new TH2F("photon_rE_v_E_CB", "#gamma CB; E_{rec} (MeV); photon rel En", 1000, 0, 1000, 400, -2., 2.);

    z_vertex                    = new TH1F("z_vertex", "z_vertex", 200, -10., 10.);
    logRelE                     = new TH1F("logRelE", "logRelE", 500, 0, 5);

    thcorr_CB                   = new TFile("configfiles/corr/CB_th_corr.root");
    dthvth_CB                   = (TProfile*)thcorr_CB->Get("photon_dtheta_v_theta_CB_pfx");

    Ecorr_CB                    = new TFile("configfiles/corr/CB_e_corr.root");
    EvdetCB                     = (TH2F*)Ecorr_CB->Get("g_peak_E_CB");

    Ecorr_TAPS                  = new TFile("configfiles/corr/TAPS_e_corr.root");
    EvdetTAPS                   = (TH2F*)Ecorr_TAPS->Get("g_peak_E_TAPS");

    thcorr_TAPS                 = new TFile("configfiles/corr/TAPS_th_corr.root");
    dthvth_TAPS                 = (TProfile*)thcorr_TAPS->Get("photon_dtheta_v_theta_TAPS_pfx");




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
    z_vertex->Fill(GetGeant()->GetVertex().Z());

    Double_t E_true = trueobs.E()*1000;
    Double_t th_true = trueobs.Theta()*TMath::RadToDeg();
    Double_t fi_true = trueobs.Phi()*TMath::RadToDeg();


    Energy_corr();
    theta_corr();




    // particle:
    // 0  photon
    // 1  proton
    int particle = 0;

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

            photon_theta_phi_CB->Fill(GetTracks()->GetTheta(i), GetTracks()->GetPhi(i));

            Double_t E_rec    =  GetTracks()->GetClusterEnergy(i);
            Double_t th_rec   =  GetTracks()->GetTheta(i);

            Double_t rel_E  =  ( GetTracks()->GetClusterEnergy(i) - E_true )/(GetTracks()->GetClusterEnergy(i));
            Double_t dth    =  ( GetTracks()->GetTheta(i) - th_true );
            Double_t dfi    =  ( GetTracks()->GetPhi(i) - fi_true );
            if( k < 2 ) // plot the split-off E_rec and angle vs E_true (theta assumed at CB edge)
            {
                photon_secondary_Evth->Fill( GetTracks()->GetClusterEnergy(k), GetTracks()->GetTheta(k), E_true );
            }

            photon_rE_v_2D->Fill(E_rec, th_rec, rel_E);
            photon_dtheta_v_2D->Fill(E_rec, th_rec, dth);
            photon_dphi_v_2D->Fill(E_rec, th_rec, dfi);

            if(GetTracks()->HasCB(i))
            {
                photon_true_Evth_CB->Fill(E_true, th_true);
                photon_rE_v_2D_CB->Fill(E_rec, th_rec, rel_E);
                photon_dtheta_v_2D_CB->Fill(E_rec, th_rec, dth);
                photon_dphi_v_2D_CB->Fill(E_rec, th_rec, dfi);

                photon_rE_v_det_CB->Fill(E_rec, GetTracks()->GetCentralCrystal(i), rel_E);
                photon_dtheta_v_det_CB->Fill(E_rec, GetTracks()->GetCentralCrystal(i), dth);
                photon_dphi_v_det_CB->Fill(E_rec, GetTracks()->GetCentralCrystal(i), dfi);

                photon_EtEr_v_det_CB->Fill(E_rec,  GetTracks()->GetCentralCrystal(i), E_true/E_rec);

                int idet_mod = int(GetTracks()->GetCentralCrystal(i)/16);
                photon_EtEr_v_detmod_CB->Fill(E_rec, idet_mod, rel_E);


                photon_dtheta_v_theta_CB->Fill(th_rec,dth );
                double test = (TMath::Log(-rel_E));
                logRelE->Fill(test);

                photon_rE_v_E_CB->Fill(E_rec, rel_E);

            }
            else if(GetTracks()->HasTAPS(i)) //TAPS
            {
                photon_true_Evth_TAPS->Fill(E_true, th_true);
                photon_rE_v_2D_TAPS->Fill(E_rec, th_rec, rel_E);
                photon_dtheta_v_2D_TAPS->Fill(E_rec, th_rec, dth);
                photon_dphi_v_2D_TAPS->Fill(E_rec, th_rec, dfi);

                photon_dtheta_v_theta_TAPS->Fill(th_rec,dth );

                photon_rE_v_det_TAPS->Fill(E_rec, GetTracks()->GetCentralCrystal(i), rel_E);
                photon_dtheta_v_det_TAPS->Fill(E_rec, GetTracks()->GetCentralCrystal(i), dth);
                photon_dphi_v_det_TAPS->Fill(E_rec, GetTracks()->GetCentralCrystal(i), dfi);

                photon_EtEr_v_det_TAPS->Fill(E_rec,  GetTracks()->GetCentralCrystal(i), E_true/E_rec);

                Double_t ztrue   =  0.0;
                Double_t R_true  =  (145.7-ztrue)*(TMath::Tan(trueobs.Theta()));
                Double_t R_rec   =  R_TAPS[GetTracks()->GetCentralCrystal(i)];
                Double_t dR      =  R_rec - R_true;

                photon_dR_v_det_TAPS->Fill( GetTracks()->GetCentralCrystal(i)  ,dR );


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

                proton_rE_v_det_TAPS->Fill(E_rec, GetTracks()->GetCentralCrystal(i), rel_E);
                proton_dtheta_v_det_TAPS->Fill(E_rec, GetTracks()->GetCentralCrystal(i), dth);
                proton_dphi_v_det_TAPS->Fill(E_rec, GetTracks()->GetCentralCrystal(i), dfi);

                proton_re_v_theta->Fill(th_rec , rel_E);

                proton_dth_v_detnr->Fill(dth, GetTracks()->GetCentralCrystal(i));
                proton_dphi_v_detnr->Fill(dfi, GetTracks()->GetCentralCrystal(i));

                proton_EvDE_TAPS->Fill(GetTracks()->GetClusterEnergy(i), GetTracks()->GetVetoEnergy(i));
                proton_dtheta_v_theta_TAPS->Fill(th_rec, dth);

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
    std::string         line;

    std::ifstream fileCB("configfiles/data/CB_lingain.txt");
    std::getline(fileCB, line);
    std::string         buffer2;
    stringstream        ss;
    ss << line;
    while (std::getline(ss, buffer2, '\t'))
    {
        CBgain.push_back(std::stod(buffer2));
    }
    fileCB.close();

    std::ifstream fileTAPS("configfiles/data/TAPS_lingain.txt");
    std::getline(fileTAPS, line);
    std::string         buffer3;
    std::stringstream   ss2;
    ss2 << line;
    while (std::getline(ss2, buffer3, '\t'))
    {
        TAPSgain.push_back(std::stod(buffer3));
    }
    fileTAPS.close();

    R_TAPS_corr.resize(0);
    std::ifstream fileTAPSRcorr("configfiles/corr/TAPS_R_corr.txt");
    if(fileTAPSRcorr){
        std::getline(fileTAPSRcorr, line);
        std::string         buffer4;
        std::stringstream   ss3;
        ss3 << line;
        while (std::getline(ss3, buffer4, '\t'))
        {
            R_TAPS_corr.push_back(std::stod(buffer4));
        }

        fileTAPSRcorr.close();
    }
    else
        for(int i = 0; i < 438; i++)
            R_TAPS_corr.push_back(0.0);

    R_TAPS.resize(0);
    std::ifstream fileTAPScentrcoord("configfiles/BaF2-PbWO4.dat");
    if(fileTAPScentrcoord){
        int el;
        double x,y,z, R;
        for(int iel = 0; iel < 438; iel++){
             fileTAPScentrcoord >> el >> x >> y >> z;
             R = TMath::Sqrt(x*x + y*y);
             R_TAPS.push_back(R);
        }
    }
    else{
        std::cout << "Did not find the file BaF2-PbWO4.dat"<<std::endl;
        exit(0);
    }

    for(int iR = 0 ; iR < 438; iR++)
    {
        R_TAPS[iR] -=  R_TAPS_corr[iR];
    }


   std::cout << " CB gain correction applied" << std::endl;
   std::cout << " TAPS gain correction applied" << std::endl;
    return kTRUE;
}


void    TrueRecParticleGun::theta_corr()
{
    Double_t th_c , dth;
    for(int tr = 0; tr < GetTracks()->GetNTracks(); tr++)
    {
        if( GetTracks()->HasCB(tr) )
        {
            dth = (Double_t)dthvth_CB->GetBinContent(dthvth_CB->FindBin(GetTracks()->GetTheta(tr)));
            th_c = GetTracks()->GetTheta(tr) - dth;
            tracks->SetTheta(tr, th_c);
        }
        if( GetTracks()->HasTAPS(tr) && GetTracks()->GetTheta(tr) < 22.0 )
        {
            dth = (Double_t)dthvth_TAPS->GetBinContent(dthvth_TAPS->FindBin(GetTracks()->GetTheta(tr)));
            th_c = GetTracks()->GetTheta(tr) - dth;
            tracks->SetTheta(tr, th_c);
        }
    }
    return;
}

void TrueRecParticleGun::Energy_corr()
{
    Double_t Erec, Ec_temp, DeltaE, Ec;
    int det;
    for (int i = 0; i < GetTracks()->GetNTracks() ; i++)
    {
        if( GetTracks()->HasCB(i) )
        {
            Erec = GetTracks()->GetVector(i).E();
            Ec_temp = CBgain[GetTracks()->GetCentralCrystal(i)]*Erec;
            DeltaE = Ec_temp*(Double_t)EvdetCB->GetBinContent(EvdetCB->FindBin(Ec_temp, GetTracks()->GetTheta(i)));
//             ( GetTracks()->GetClusterEnergy(i) - E_true )/(GetTracks()->GetClusterEnergy(i));

            Ec = Ec_temp - DeltaE;
//            Ec = Ec_temp*(Double_t)EvdetCB->GetBinContent(EvdetCB->FindBin(Ec_temp, GetTracks()->GetCentralCrystal(i)));
//            Ec  = Ec_temp;
            tracks->SetClusterEnergy(i, Ec);
        }
        else if(GetTracks()->HasTAPS(i) )
        {
            Erec = GetTracks()->GetVector(i).E();
            Ec_temp = TAPSgain[GetTracks()->GetCentralCrystal(i)]*Erec;
            DeltaE = Ec_temp*(Double_t)EvdetTAPS->GetBinContent(EvdetTAPS->FindBin(Ec_temp, GetTracks()->GetTheta(i)));
            Ec = Ec_temp - DeltaE;

            tracks->SetClusterEnergy(i, Ec);
        }
    }
};
