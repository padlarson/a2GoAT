#include "GTreeManager.h"
#include "ProtonChargeRadius.h"
#include "GTrue.h"
#include "GTreeA2Geant.h"
#include "TH1.h"



ProtonChargeRadius::ProtonChargeRadius()
{
    // Beam Energy
    True_BeamEnergy    = new GHistBGSub("True_BeamEnergy", "True Beam Energy (MeV)", 400, 0., 1.6);
    // Phase space observables
    ThvT_p              = new GHistBGSub2("ThvT_p","proton true; Energy p (MeV); #theta_{p} (^{o})", 300, 0., 150, 100, 0., 100.);
    proton_p_v_th       = new GHistBGSub2("proton_p_v_th", "proton; p_{true} (MeV/c); #theta_{true} (^{o})", 700, 0., 1400., 100, 0, 100);

    proton_p_vs_t       = new GHistBGSub2("proton_p_vs_t", "proton; p_{true} (MeV/c); -t(^{o})", 1000, 0., 500., 100, 0, 0.1);

    BeamEvETh_p         = new GHistBGSub2("BeamEvETh_p","proton true; Beam Energy #gamma (MeV); bin = int(T_{p}/10) + 1 + 20*int(#theta_{p}/10)", 160, 0., 1600, 300, 0., 300.);

    ThvT_e              = new GHistBGSub2("ThvT_e","electron/positron true; Energy e^{+}/e^{-} (MeV); #theta_{e} (^{o})", 75, 0., 1500, 90, 0., 180.);
    ThvT_mu             = new GHistBGSub2("ThvT_mu","#mu true; Energy #mu^{+}/#mu^{-} (MeV); #theta_{#mu} (^{o})", 75, 0., 1500, 90, 0., 180.);

    e_acceptance          = new GH1("e_acceptance","0 - both, 1 - 1 lept, 2 - none (< 45 degrees); ", 4, -1, 3);
    mu_acceptance         = new GH1("mu_acceptance","0 - both, 1 - 1 lept, 2 - none (< 45 degrees); ",4, -1, 3);


    ThpvEg              = new GHistBGSub2("ThpvEg","proton angles vs Beam Energy; Beam Energy #gamma (MeV); #theta_{p}" , 100, 100., 600, 200, 0., 100.);

    dsigma_v_mll        = new GHistBGSub2("dsigma_v_mll", " differential cross-section; M_{ll}^{2}; d#sigma  (#mub /GeV^{4}) ", 1000, 0.0, 0.5, 1000, 0., 1000);
    t_v_mll             = new GHistBGSub2("t vs M_{ll}", ";M_{ll}^{2} (GeV^{2}); mom transfer (-t) (GeV^{2})", 125, 0.0, 0.25, 500, 0, 0.5);
    thlabpr_v_mll       = new GHistBGSub2("thlabpr_v_mll", "#theta_{lab, pr} (y) vs M_{ll}^{2} (x)", 400, 0.0, 0.2, 200, 0, 100);
    thlabpr_v_t         = new GHistBGSub2("thlabpr_v_t", ";-t (GeV^{2});#theta_{lab, proton} (^{o})", 25, 0.0, 0.05, 100, 0, 100);

    mll                 = new GH1("mll", "M_{l^{+}l^{-}} (GeV^{2})", 100, 0.0, 0.20);

    Ntracks             = new GH1("Ntracks", "Nr of reconstructed tracks", 10, 0, 10);

    E_v_dE              = new GHistBGSub2("E_v_dE", "proton; E (MeV); #DeltaE (MeV)", 400, 0., 400, 100, 0., 10.);

    dt_vs_t_rec         = new GHistBGSub2("dt_vs_t_rec","#Delta t vs -t; #-t_{rec} (GeV^{2}); #Delta#t", 50, 0.0, 0.2, 50, -1., 1.);
    dth_p_v_th_p        = new GHistBGSub2("dth_p_v_th_p","#Delta#theta vs #theta_{rec, p}; #theta_p; #Delta#theta_{p}", 90, 0., 90, 50, -5., 5.);
    dE_p_v_E_p          = new GHistBGSub2("dE_p_v_E_p","#DeltaE vs E_{rec, p}", 90, 0., 90, 50, -5., 5.);

    proton_E_v_th      =  new GHistBGSub2("proton_E_v_th", "proton; E_{rec} (MeV); #theta_{rec} (^{o})", 50, 0., 150., 50, 0, 100);
    proton_t_v_th      =  new GHistBGSub2("proton_t_v_th", "proton; -t_{rec} (GeV^{2}); #theta_{rec} (^{o})", 125, 0.0, 0.25, 50, 0, 100);

    proton_E_v_mllth   =  new GHistBGSub2("proton_E_v_mllth", "proton; E_{true} (MeV); #theta_{true,ll} (^{o})", 50, 0., 150., 180, 0, 180);


    proton_rE_v_2D     =  new GHistBGSub3("proton_rE_v_2D", "proton; E_{rec} (MeV); #theta_{rec} (^{o}); relative (E_{rec} - E_{true}) / E_{rec} ", 75, 0., 150., 40, 0, 80, 100, -5., 5.);
    proton_dtheta_v_2D =  new GHistBGSub3("proton_dtheta_v_2D", "proton; E_{rec} (MeV); #theta_{rec} (^{o}); #theta_{rec} - #theta_{true} (^{o})", 75, 0., 150., 40, 0, 80, 40, -40., 40.);
    proton_dphi_v_2D   =  new GHistBGSub3("proton_dphi_v_2D", "proton; E_{rec} (MeV); #theta_{rec} (^{o}); #phi_{rec} - #phi_{true} (^{o})", 75, 0., 150., 40, 0, 80, 40, -40., 40.);

    proton_dtheta_v_E  =  new GHistBGSub2("proton_dtheta_v_E", "proton; E_{rec} (MeV); #theta_{rec} - #theta_{true} (^{o})", 15, 0., 150., 160, -40., 40.);
    proton_rE_v_E      =  new GHistBGSub2("proton_rE_v_E", "proton; E_{rec} (MeV);relative (E_{rec} - E_{true}) / E_{rec} ", 40, 0., 200., 200, -1., 1.);
    proton_rt_v_E      =  new GHistBGSub2("proton_rt_v_E", "proton; E_{rec} (MeV);relative (t_{rec} - t_{true}) / t_{rec} ", 40, 0., 200., 200, -1., 1.);

    diff_msqll         =  new GH1("diff_msqll", " m^{2}_{l^{+}l^{-},rec}-m^{2}_{l^{+}l^{-},true}; GeV^{2}", 500, -.005, .005);
    ds                 =  new GH1("ds", " s_{rec}-s_{true}", 200, -0.05, .05);
    dtheta             =  new GH1("dtheta", " th_{rec}-th_{true}", 200, -1., 1.);
    dt_v_t             =  new GHistBGSub2("dt_v_t", "t_{rec}, t_{rec}-t_{true}", 500, 0., 0.05, 500, -0.005 , 0.005);
    dmll_v_t           =  new GHistBGSub2("dmll_v_t", "t_{rec}, m^{2}_{l^{+}l^{-},rec}-m^{2}_{l^{+}l^{-},true}", 50, 0., 0.05, 200, -.01, .01);
    dmll_v_theta       =  new GHistBGSub2("dmll_v_theta", "theta_{rec}, m^{2}_{l^{+}l^{-},rec}-m^{2}_{l^{+}l^{-},true}", 500, 0., 100, 200, -.01, .01);

    // -t = 0.01
    dmll_v_bin         =  new GHistBGSub2("dmll_v_bin", ";binnr; m^{2}_{l^{+}l^{-},rec}-m^{2}_{l^{+}l^{-},true}", 2200, 0., 2200, 100, -.005, .005);
    dt_v_bin           =  new GHistBGSub2("dt_v_bin", ";binnr; t_{rec}-t_{true}, t = -0.01", 2200, 0., 2200, 1000, -0.05 , 0.05);

    // -t = 0.02
    dmll_v_bin2         =  new GHistBGSub2("dmll_v_bin2", ";binnr; m^{2}_{l^{+}l^{-},rec}-m^{2}_{l^{+}l^{-},true}", 2200, 0., 2200, 100, -.005, .005);
    dt_v_bin2           =  new GHistBGSub2("dt_v_bin2", ";binnr; t_{rec}-t_{true}, t = -0.02", 2200, 0., 2200, 400, -0.05 , 0.05);

    // -t = 0.03
    dmll_v_bin3         =  new GHistBGSub2("dmll_v_bin3", ";binnr; m^{2}_{l^{+}l^{-},rec}-m^{2}_{l^{+}l^{-},true}", 2200, 0., 2200, 100, -.005, .005);
    dt_v_bin3           =  new GHistBGSub2("dt_v_bin3", ";binnr; t_{rec}-t_{true}, t = -0.03", 2200, 0., 2200, 400, -0.05 , 0.05);

    cutFile             = new TFile("configfiles/cuts/CB_DeltaE_E_prchrad.root");
    cutProtonCB         = (TCutG*)cutFile->Get("Proton");

    GHistBGSub::InitCuts(-20, 20, -55, -35);
    GHistBGSub::AddRandCut(35, 55);
}

ProtonChargeRadius::~ProtonChargeRadius()
{

}
Bool_t	ProtonChargeRadius::Start()
{
//    if(!IsGoATFile())
//    {
//        cout << "ERROR: Input File is not a GoAT file." << endl;
//        return kFALSE;
//    }
    pRandoms = new TRandom3;
    pRandoms->SetSeed(0); //'Randomize Timer'

    SetAsPhysicsFile();

    TraverseValidEvents();

    outputFile->cd();
    gDirectory->mkdir("MC")->cd();


	return kTRUE;
}

void	ProtonChargeRadius::ProcessEvent()
{

    TrueObs.Start(*GetPluto(), *GetGeant());

    int reaction = 1;
    // reaction:
    // 0  p e+e- production   : 1/E x Bethe-Heidler process.
    // 1  p mu+mu- production   : 1/E x Bethe-Heidler process.
    // 2  p gamma             : 1/E x Total cross section.
    // 3  p pi0               : 1/E x Total cross section.
    // 4  p pi0pi0            : 1/E x Total cross section.
    // 5  p pi+pi-            : 1/E x Total cross section.
    TrueAnalysis_ll(reaction); //lepton-lepton analysis

    Double_t weight = TrueObs.GetWeight();
    True_BeamEnergy->FillWeighted( TrueObs.GetTrueBeamEnergy(), TrueObs.GetWeight() );


    Double_t E_true = TrueObs.GetTrueProtonLV().E()*1000 - MASS_PROTON;
    Double_t p_true = TMath::Sqrt(TrueObs.GetTrueProtonLV().E()*1000*TrueObs.GetTrueProtonLV().E()*1000 - MASS_PROTON*MASS_PROTON);
    Double_t th_true = TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg();
    Double_t fi_true = TrueObs.GetTrueProtonLV().Phi()*TMath::RadToDeg();


    ThvT_p->FillWeighted( TrueObs.GetTrueProtonLV().E()*1000 - MASS_PROTON, TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight() );
    proton_p_v_th->FillWeighted( p_true , TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight() );

    proton_E_v_mllth->FillWeighted( TrueObs.GetTrueProtonLV().E()*1000 - MASS_PROTON, TrueObs.GetTrueElectronLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight() );
    proton_E_v_mllth->FillWeighted( TrueObs.GetTrueProtonLV().E()*1000 - MASS_PROTON, TrueObs.GetTruePositronLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight() );

    ThpvEg->FillWeighted( TrueObs.GetTrueBeamEnergy()*1000, TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight());

    ThvT_e->FillWeighted( TrueObs.GetTrueElectronLV().E()*1000 - MASS_ELECTRON, TrueObs.GetTrueElectronLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight() );
    ThvT_e->FillWeighted( TrueObs.GetTruePositronLV().E()*1000 - MASS_ELECTRON, TrueObs.GetTruePositronLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight() );

    ThvT_mu->FillWeighted( TrueObs.GetTrueMuonNegLV().E()*1000 - MASS_MUON, TrueObs.GetTrueMuonNegLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight() );
    ThvT_mu->FillWeighted( TrueObs.GetTrueMuonPosLV().E()*1000 - MASS_MUON, TrueObs.GetTrueMuonPosLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight() );

    Int_t bin_nr = Int_t(( TrueObs.GetTrueProtonLV().E()*1000 - MASS_PROTON ) / 10) +1 + 20*Int_t(TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg()/10.0);
    BeamEvETh_p->FillWeighted(TrueObs.GetTrueBeamEnergy()*1000, bin_nr, TrueObs.GetWeight() );

    if(p_true > 1200) return;

    double theta_lim = 45.0;
    double acceptance = -1;
    if((TrueObs.GetTrueElectronLV().Theta()*TMath::RadToDeg() < theta_lim) && (TrueObs.GetTruePositronLV().Theta()*TMath::RadToDeg() < theta_lim)){
        acceptance = 0;
    }
    if((TrueObs.GetTrueElectronLV().Theta()*TMath::RadToDeg() > theta_lim) && (TrueObs.GetTruePositronLV().Theta()*TMath::RadToDeg() < theta_lim))
        acceptance = 1;
    if((TrueObs.GetTrueElectronLV().Theta()*TMath::RadToDeg() < theta_lim) && (TrueObs.GetTruePositronLV().Theta()*TMath::RadToDeg() > theta_lim))
        acceptance = 1;
    if((TrueObs.GetTrueElectronLV().Theta()*TMath::RadToDeg() > theta_lim) && (TrueObs.GetTruePositronLV().Theta()*TMath::RadToDeg() > theta_lim))
        acceptance = 2;

    e_acceptance->FillWeighted(acceptance , TrueObs.GetWeight());


    acceptance = -1;
    if((TrueObs.GetTrueMuonNegLV().Theta()*TMath::RadToDeg() < theta_lim) && (TrueObs.GetTruePositronLV().Theta()*TMath::RadToDeg() < theta_lim)){
        acceptance = 0;
    }
    else if((TrueObs.GetTrueMuonNegLV().Theta()*TMath::RadToDeg() > theta_lim) && (TrueObs.GetTrueMuonPosLV().Theta()*TMath::RadToDeg() < theta_lim))
        acceptance = 1;
    else if((TrueObs.GetTrueMuonNegLV().Theta()*TMath::RadToDeg() < theta_lim) && (TrueObs.GetTrueMuonPosLV().Theta()*TMath::RadToDeg() > theta_lim))
        acceptance = 1;
    else
        acceptance = 2;

    mu_acceptance->FillWeighted(acceptance , TrueObs.GetWeight());


//    CB_analysis();

//    BoNus_analysis();


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

void ProtonChargeRadius::BoNus_analysis(){

    Double_t p_true = TMath::Sqrt(TrueObs.GetTrueProtonLV().E()*1000*TrueObs.GetTrueProtonLV().E()*1000 - MASS_PROTON*MASS_PROTON);
    Double_t th_true = TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg();
    Double_t th_true_rad = TrueObs.GetTrueProtonLV().Theta();

    // sigma_theta_p 50 mrad; sigma_phi = 25 mrad; sigma_R 0.53 mm; magnetic field 4 T
    // sigma_p = 10 %, 20 %

    Double_t s_true     = Get_s(TrueObs.GetTrueBeamEnergy());
    Double_t t_true     = GetT(p_true/1000.);
    Double_t m_ll_true  = Get_mll(s_true, t_true, th_true_rad);

    // do alternative analysis here. Check if proton is in the Goldilock region/ sweet spot, above 50 degrees
    // check if proton momentum is in either of the three intervals 90-110 MeV/c, 132-152, 164-184 MeV/c


    if(th_true > 50.0){
        if( (p_true > 95 && p_true < 105) || (p_true > 137 && p_true < 147) || (p_true > 169 && p_true < 179)){
            for(int iter = 0; iter < 100; iter++ ){
                Double_t beam_rec   = TrueObs.GetTrueBeamEnergy() + pRandoms->Uniform(0,1.0e-4);

                Double_t rand1      = pRandoms->Uniform(0, 0.020);
                Double_t rand2      = pRandoms->Uniform(0, 0.005);

                Double_t p_rec      = p_true + pRandoms->Gaus(0, rand1*p_true);
                Double_t th_rec_rad = th_true_rad + pRandoms->Gaus(0, rand2);

                Double_t s_rec      = Get_s(beam_rec);
                Double_t t_rec      = GetT(p_rec/1000.);
                Double_t m_ll_rec   = Get_mll(s_rec, t_rec, th_rec_rad);


                Double_t diff_mll = m_ll_rec - mll_sq;
                diff_msqll->FillWeighted(diff_mll, TrueObs.GetWeight());
                ds->FillWeighted(s_rec-s_true, TrueObs.GetWeight());
                dtheta->FillWeighted(th_rec_rad-th_true_rad, TrueObs.GetWeight());
                dt_v_t->FillWeighted(-t_rec, t_rec-t_true, TrueObs.GetWeight());
                dmll_v_t->FillWeighted(-t_rec,diff_mll,TrueObs.GetWeight());
                dmll_v_theta->FillWeighted(th_rec_rad*TMath::RadToDeg(), diff_mll, TrueObs.GetWeight());

                // 20 in momentum res
                // 50 in theta res

                int bin = int(rand1/0.001) + 40*int(rand2/0.0002);
                if((p_true > 95 && p_true < 105) ){
                    dmll_v_bin->FillWeighted(bin, diff_mll, TrueObs.GetWeight());
                    dt_v_bin->FillWeighted(bin, t_rec-t_true, TrueObs.GetWeight());
                }

                if(p_true > 137 && p_true < 147){
                    dmll_v_bin2->FillWeighted(bin, diff_mll, TrueObs.GetWeight());
                    dt_v_bin2->FillWeighted(bin, t_rec-t_true, TrueObs.GetWeight());
                }

                if(p_true > 169 && p_true < 179){
                    dmll_v_bin3->FillWeighted(bin, diff_mll, TrueObs.GetWeight());
                    dt_v_bin3->FillWeighted(bin, t_rec-t_true, TrueObs.GetWeight());
                }
          }
       }
    }
}

Double_t ProtonChargeRadius::GetT(Double_t pr_lab_mom){
    Double_t mp         =   MASS_PROTON/1000.;
    return (2*mp*mp -TMath::Sqrt(4*mp*mp*mp*mp + (pr_lab_mom*pr_lab_mom)*(4*mp*mp) ));
}

Double_t ProtonChargeRadius::Get_mll(Double_t s, Double_t t, Double_t theta){
    Double_t mp         =   MASS_PROTON/1000;
    Double_t tau        =   (-t)/(4*mp*mp);

    Double_t term_1     = 2*(s-mp*mp)*TMath::Sqrt(tau*(1+tau));
    Double_t term_2     = 2*(s+mp*mp)*tau;
    Double_t mll = TMath::Cos(theta)*term_1 - term_2;
    return mll;
}

Double_t ProtonChargeRadius::Get_s(Double_t beam){
    Double_t mp         =   MASS_PROTON/1000;
    return (mp*mp + 2*mp*beam);   // Mandelstam s;
}

void ProtonChargeRadius::TrueAnalysis_ll(int reaction) // implementation BH-process
{
    Double_t weight = 1.0;
    Double_t norm;
         if( reaction < 2 ){
            Double_t ml;
            if(reaction == 0)
                ml = MASS_ELECTRON/1000;
            else
                ml = MASS_MUON/1000;

            if( TMath::Abs(ml - MASS_ELECTRON/1000) < 1.0e-4) //electrons/positrons
                mll_sq = (TrueObs.GetTrueElectronLV() + TrueObs.GetTruePositronLV()).M2();
            else                                                    //muons
                mll_sq = (TrueObs.GetTrueMuonNegLV() + TrueObs.GetTrueMuonPosLV()).M2();

            Double_t mp         =   MASS_PROTON/1000;
            Double_t hbarc      =   0.197327;
            Double_t pr_lab_mom =   TrueObs.GetTrueProtonLV().P();
            Double_t alpha      =   1.0/137.036;

            Double_t beta       =   TMath::Sqrt( 1 - ( (4*ml*ml)/mll_sq ) );
            Double_t s          =   (MASS_PROTON*MASS_PROTON + 2*MASS_PROTON*TrueObs.GetTrueBeamEnergy()*1000.0)/1.0e6;   // Mandelstam s

            Double_t t          =   GetT(pr_lab_mom);
            Double_t tau        =   (-t)/(4*mp*mp);


            Double_t term1, term2, term3;

            term1               =   TMath::Power(hbarc,2)*TMath::Power(10,4)*(alpha*alpha*alpha)/( (s-mp*mp)*(s-mp*mp) );
            term2               =   (4*beta)/( t*t*TMath::Power((mll_sq-t),4) );
            term3               =   1./(1.+tau);

            Double_t CE1, CE2, CM1, CM2;
            Double_t CE, CM;

            CE1                 =   (t)*(s - mp*mp)*(s - mp*mp - mll_sq + t)*( TMath::Power(mll_sq,2) + 6*mll_sq*t + t*t + 4*ml*ml*mll_sq ) + (TMath::Power((mll_sq -t),2))*( t*t*mll_sq + mp*mp*( TMath::Power((mll_sq+t),2) + 4*ml*ml*mp*mp*mll_sq));
            CE2                 =  (-t)*(s - mp*mp)*(s - mp*mp - mll_sq + t)*( TMath::Power(mll_sq,2) + t*t + 4*ml*ml*(mll_sq+2*t-2*ml*ml))+ TMath::Power((mll_sq - t), 2)*(-1*mp*mp*(mll_sq*mll_sq + t*t ) + 2*ml*ml*(-1*t*t - 2*mp*mp*mll_sq + 4*ml*ml*mp*mp));
            CE                  =  CE1 + CE2*(1/beta)*TMath::Log((1 + beta)/(1 - beta)) ;

            CM1                 =  CE1 - 2*mp*mp*(1 + tau)*TMath::Power((mll_sq-t),2)*(mll_sq*mll_sq + t*t + 4*ml*ml*mll_sq);
            CM2                 =  CE2 + 2*mp*mp*(1 + tau)*TMath::Power((mll_sq-t),2)*(mll_sq*mll_sq + t*t + 4*ml*ml*(mll_sq - t -2*ml*ml));  // Timothy said it is + and not -, check
            CM                  =  CM1 + CM2*(1/beta)*TMath::Log((1 + beta)/(1 - beta)) ;

            Double_t diffxs     =   term1*term2*term3*(CE*GEp(-t)*GEp(-t) + CM*tau*GMp(-t)*GMp(-t));

            // end of implementation BH process. Below a normalisation factor is found to match the nr of entries to integrated nr
            // Different if l = mu or e, and different given the beam photon energy.

           //            // norm 500 MeV e and mu
//              if(reaction == 0)
//                  norm = 0.000693963*5.0/5.9297;
//              else
//                  norm = 14.115561277*2.0/2.007587;

           // norm all beam energies MeV e and mu
//            if(reaction == 0)
//              norm = 1.0;
//            else
//              norm = 1.0e6/24047.0*(1.0e6/1153485.0)*(1.0e6/1391015.0);

//            if(reaction == 0) // 1500
//                norm = 1./245697563.*5000000.;
//            else
//                norm = 5000000./101149.5;

//            if(reaction == 0) //800 MeV
//                norm = 1./1085492244.*5000000.;
//            else
//                norm = 5.0e6/286049.;

            if(reaction == 0) //400 MeV
                norm = 1./964.323e9*5000000.;
            else
                norm = 5.0e6/765445.;

            Double_t weight_beame;
            Double_t weight_diffxs;

            weight_diffxs = diffxs*norm;
            //

            if(weight_diffxs < 0 )
                weight_diffxs = 0;



            TrueObs.SetWeight(weight_diffxs);

            mll->FillWeighted(mll_sq, TrueObs.GetWeight() );
            dsigma_v_mll->FillWeighted(mll_sq,diffxs, weight);
            t_v_mll->FillWeighted(mll_sq, -t, TrueObs.GetWeight());
            thlabpr_v_mll->FillWeighted(mll_sq, TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg(),  TrueObs.GetWeight());

            Double_t p_true = TMath::Sqrt(TrueObs.GetTrueProtonLV().E()*1000*TrueObs.GetTrueProtonLV().E()*1000 - MASS_PROTON*MASS_PROTON);
            if( (p_true > 70.) && (p_true < 150.))
                thlabpr_v_t->FillWeighted(-t, TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg(),  TrueObs.GetWeight());

            dsigma_v_mll->FillWeighted(mll_sq, diffxs, TrueObs.GetWeight());

         }
         if(reaction == 3){
             norm = 1.0;
             weight = (1.0/TrueObs.GetTrueBeamEnergy())*norm;
             TrueObs.SetWeight(weight);
         }
         if(reaction == 4){
             norm = 1.0;
             weight = (1.0/TrueObs.GetTrueBeamEnergy())*norm;
             TrueObs.SetWeight(weight);
         }
}

Double_t ProtonChargeRadius::GEp(Double_t Qsqr)
    {
    return  1./TMath::Power(1. + Qsqr/0.71,2)*(1. - 0.4980*Qsqr+5.4592*TMath::Power(Qsqr,2) -
            34.7281 * TMath::Power(Qsqr,3) + 124.3173 * TMath::Power(Qsqr,4) - 262.9808 * TMath::Power(Qsqr,5) +
            329.1395 * TMath::Power(Qsqr,6) - 227.3306 * TMath::Power(Qsqr,7) + 66.6980 * TMath::Power(Qsqr,8));
    }

Double_t ProtonChargeRadius::GMp(Double_t Qsqr)
    {
    Double_t PMAGN      =   2.79285;
    return PMAGN/TMath::Power(1 + Qsqr/0.71,2) * (1 + 0.2472 * Qsqr - 4.9123 * TMath::Power(Qsqr,2) +
           29.7509 * TMath::Power(Qsqr,3) - 84.0430 * TMath::Power(Qsqr,4) + 129.3256 * TMath::Power(Qsqr,5) -
           111.1068 * TMath::Power(Qsqr,6) + 49.9753 * TMath::Power(Qsqr,7) - 9.1659 * TMath::Power(Qsqr,8));
    }



void ProtonChargeRadius::CB_analysis(){

    Double_t E_true = TrueObs.GetTrueProtonLV().E()*1000 - MASS_PROTON;
    Double_t th_true = TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg();
    Double_t fi_true = TrueObs.GetTrueProtonLV().Phi()*TMath::RadToDeg();

    Ntracks->FillWeighted(GetTracks()->GetNTracks(), TrueObs.GetWeight());

    for(UInt_t i = 0; i < GetTracks()->GetNTracks(); i++)
    {
        E_v_dE->FillWeighted( GetTracks()->GetClusterEnergy(i), GetTracks()->GetVetoEnergy(i), TrueObs.GetWeight());

        Double_t E    =  GetTracks()->GetClusterEnergy(i);
        if(cutProtonCB->IsInside(E, GetTracks()->GetVetoEnergy(i)))
        {
            Double_t E    =  GetTracks()->GetClusterEnergy(i);
            double p[5] = {-4.8794, 0.220078, -0.00402405, 3.27302e-05, -9.58383e-08};
            Double_t E_rec = E - E*(-4.8794 + 0.220078*E - 0.00402405*E*E + 3.27302e-05*E*E*E - 9.58383e-08*E*E*E*E) + 0.04*E;  // "calibration of proton energy"
            Double_t mp = MASS_PROTON/1000;
            Double_t pr_lab_mom = TMath::Sqrt( (E_rec + MASS_PROTON)*(E_rec + MASS_PROTON) - MASS_PROTON*MASS_PROTON )/1000;
            Double_t t_rec =   2*mp*mp -TMath::Sqrt(4*mp*mp*mp*mp + (pr_lab_mom*pr_lab_mom)*(4*mp*mp) );           // momentum transfer t

            Double_t rel_E  = ( E_rec- E_true )/(E_rec);
            Double_t rel_t    = ( t_rec - t_true )/t_rec;
            Double_t dth    =  ( GetTracks()->GetTheta(i) - th_true );
            Double_t dfi    =  ( GetTracks()->GetPhi(i) - fi_true );

            Double_t th_rec    =  GetTracks()->GetTheta(i);

            proton_E_v_th->FillWeighted(E_rec, th_rec, TrueObs.GetWeight());
            proton_dtheta_v_E->FillWeighted(E_rec, dth, TrueObs.GetWeight());
            proton_rE_v_E->FillWeighted(E_rec, rel_E, TrueObs.GetWeight());
            proton_rt_v_E->FillWeighted(E_rec, rel_t, TrueObs.GetWeight());
            proton_t_v_th->FillWeighted(-t_rec, th_rec, TrueObs.GetWeight());

            proton_rE_v_2D->FillWeighted(E_rec, th_rec, rel_E, TrueObs.GetWeight());
            proton_dtheta_v_2D->FillWeighted(E_rec, th_rec, dth, TrueObs.GetWeight());
            proton_dphi_v_2D->FillWeighted(E_rec, th_rec, dfi, TrueObs.GetWeight());
        }
    }
}




