#include "GTreeManager.h"
#include "ProtonChargeRadius.h"
#include "GTrue.h"
#include "GTreeA2Geant.h"
#include "TH1.h"



ProtonChargeRadius::ProtonChargeRadius()
{
    // Beam Energy
    True_BeamEnergy    = new GH1("True_BeamEnergy", "True Beam Energy (MeV)", 200, 0., 0.800);
    // Phase space observables
    ThvT_p              = new GHistBGSub2("ThvT_p","proton true; Energy p (MeV); #theta_{p} (^{o})", 300, 0., 150, 100, 0., 100.);
    ThvT_e              = new GHistBGSub2("ThvT_e","electron/positron true; Energy e^{+}/e^{-} (MeV); #theta_{e} (^{o})", 60, 0., 600, 18, 0., 180.);
    ThvT_mu             = new GHistBGSub2("ThvT_mu","#mu true; Energy #mu^{+}/#mu^{-} (MeV); #theta_{#mu} (^{o})", 200, 0., 400, 90, 0., 180.);

    ThpvEg              = new GHistBGSub2("ThpvEg","proton angles vs Beam Energy; Beam Energy #gamma (MeV); #theta_{p}" , 100, 100., 600, 200, 0., 100.);

//    dsigma_v_mll        = new GHistBGSub2("dsigma_v_mll", " differential cross-section; M_{ll}^{2}; d#sigma  (#mub /GeV^{4}) ", 1000, 0.0, 0.2, 1000, 0, 1000);
    dsigma_v_mll        = new GHistBGSub2("dsigma_v_mll", " differential cross-section; M_{ll}^{2}; d#sigma  (#mub /GeV^{4}) ", 1000, 0.0, 0.5, 1000, 0., 1000);
    t_v_mll             = new GHistBGSub2("t vs M_{ll}", ";M_{ll}^{2} (GeV^{2}); mom transfer (-t) (GeV^{2})", 125, 0.0, 0.25, 500, 0, 0.5);
    thlabpr_v_mll       = new GHistBGSub2("thlabpr_v_mll", "#theta_{lab, pr} (y) vs M_{ll}^{2} (x)", 100, 0.0, 0.2, 50, 0, 100);
    thlabpr_v_t         = new GHistBGSub2("thlabpr_v_t", ";-t (GeV^{2});#theta_{lab, proton} (^{o})", 500, 0.0, 0.5, 100, 0, 100);

    thlabpr_v_th_e_1    = new GHistBGSub2("thlabpr_v_th_e_1", "#theta p vs #theta l 0.005 -t 0.015; #theta_p; #theta_e", 50, 0., 100., 45, 0., 180. );
    thlabpr_v_th_e_2    = new GHistBGSub2("thlabpr_v_th_e_2", "#theta p vs #theta l 0.015 -t 0.025; #theta_p; #theta_e", 50, 0., 100., 45, 0., 180. );
    thlabpr_v_th_e_3    = new GHistBGSub2("thlabpr_v_th_e_3", "#theta p vs #theta l 0.025 -t 0.035; #theta_p; #theta_e", 50, 0., 100., 45, 0., 180. );

    thlabpr_v_th_mu_1    = new GHistBGSub2("thlabpr_v_th_mu_1", "#theta p vs #theta l 0.005 -t 0.015; #theta_p; #theta_mu", 50, 0., 100., 45, 0., 180. );
    thlabpr_v_th_mu_2    = new GHistBGSub2("thlabpr_v_th_mu_2", "#theta p vs #theta l 0.015 -t 0.025; #theta_p; #theta_mu", 50, 0., 100., 45, 0., 180. );
    thlabpr_v_th_mu_3    = new GHistBGSub2("thlabpr_v_th_mu_3", "#theta p vs #theta l 0.025 -t 0.035; #theta_p; #theta_mu", 50, 0., 100., 45, 0., 180. );


    thlabpr_v_Tp_1      = new GHistBGSub2("thlabpr_v_Tp_1", "0.005 -t 0.015; T_{p}; #theta_p", 100, 0., 100., 100, 0., 100.);
    thlabpr_v_Tp_2      = new GHistBGSub2("thlabpr_v_Tp_2", "0.015 -t 0.025; T_{p}; #theta_p", 100, 0., 100., 100, 0., 100.);
    thlabpr_v_Tp_3      = new GHistBGSub2("thlabpr_v_Tp_3", "0.025 -t 0.035; T_{p}; #theta_p", 100, 0., 100., 100, 0., 100.);
    thlabpr_v_Tp_13      = new GHistBGSub2("thlabpr_v_Tp_13", " 0.005 < (-t) < 0.035; T_{p}; #theta_{p} (^{o})", 40, 0., 20., 50, 0., 100.);

    MC_weight           = new GH1("MC_weight", "MC_weight", 10000, 0, 1000);
    mll                 = new GH1("mll", "M_{l^{+}l^{-}} (GeV^{2})", 100, 0.0, 0.20);

    Ntracks             = new GH1("Ntracks", "Nr of reconstructed tracks", 10, 0, 10);

    E_v_dE              = new GHistBGSub2("E_v_dE", "proton; E (MeV); #DeltaE (MeV)", 400, 0., 400, 100, 0., 10.);

    dt_vs_t_rec         = new GHistBGSub2("dt_vs_t_rec","#Delta t vs -t; #-t_{rec} (GeV^{2}); #Delta#t", 50, 0.0, 0.2, 50, -1., 1.);
    dth_p_v_th_p        = new GHistBGSub2("dth_p_v_th_p","#Delta#theta vs #theta_{rec, p}; #theta_p; #Delta#theta_{p}", 90, 0., 90, 50, -5., 5.);
    dE_p_v_E_p          = new GHistBGSub2("dE_p_v_E_p","#DeltaE vs E_{rec, p}", 90, 0., 90, 50, -5., 5.);

    proton_E_v_th      =  new GHistBGSub2("proton_E_v_th", "proton; E_{rec} (MeV); #theta_{rec} (^{o})", 50, 0., 150., 50, 0, 100);
    proton_t_v_th      =  new GHistBGSub2("proton_t_v_th", "proton; -t_{rec} (GeV^{2}); #theta_{rec} (^{o})", 125, 0.0, 0.25, 50, 0, 100);

    proton_rE_v_2D     =  new GHistBGSub3("proton_rE_v_2D", "proton; E_{rec} (MeV); #theta_{rec} (^{o}); relative (E_{rec} - E_{true}) / E_{rec} ", 75, 0., 150., 40, 0, 80, 100, -5., 5.);
    proton_dtheta_v_2D =  new GHistBGSub3("proton_dtheta_v_2D", "proton; E_{rec} (MeV); #theta_{rec} (^{o}); #theta_{rec} - #theta_{true} (^{o})", 75, 0., 150., 40, 0, 80, 40, -40., 40.);
    proton_dphi_v_2D   =  new GHistBGSub3("proton_dphi_v_2D", "proton; E_{rec} (MeV); #theta_{rec} (^{o}); #phi_{rec} - #phi_{true} (^{o})", 75, 0., 150., 40, 0, 80, 40, -40., 40.);

    proton_dtheta_v_E  =  new GHistBGSub2("proton_dtheta_v_E", "proton; E_{rec} (MeV); #theta_{rec} - #theta_{true} (^{o})", 15, 0., 150., 160, -40., 40.);
    proton_rE_v_E      =  new GHistBGSub2("proton_rE_v_E", "proton; E_{rec} (MeV);relative (E_{rec} - E_{true}) / E_{rec} ", 40, 0., 200., 200, -1., 1.);
    proton_rt_v_E      =  new GHistBGSub2("proton_rt_v_E", "proton; E_{rec} (MeV);relative (t_{rec} - t_{true}) / t_{rec} ", 40, 0., 200., 200, -1., 1.);


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

    SetAsPhysicsFile();//        if(cutProtonCB->IsInside(GetTracks()->GetClusterEnergy(i), GetTracks()->GetVetoEnergy(i)))
    //        {


    //            proton_rE_v_2D->Fill(E_rec, th_rec, rel_E);
    //            proton_dtheta_v_2D->Fill(E_rec, th_rec, dth);
    //            proton_dphi_v_2D->Fill(E_rec, th_rec, dfi);
    //        }

    TraverseValidEvents();

    outputFile->cd();
    gDirectory->mkdir("MC")->cd();


	return kTRUE;
}

void	ProtonChargeRadius::ProcessEvent()
{

    TrueObs.Start(*GetPluto(), *GetGeant());

//    Double_t weight = 1.0;
    True_BeamEnergy->Fill( TrueObs.GetTrueBeamEnergy(), TrueObs.GetWeight() );


    Double_t E_true = TrueObs.GetTrueProtonLV().E()*1000 - MASS_PROTON;
    Double_t th_true = TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg();
    Double_t fi_true = TrueObs.GetTrueProtonLV().Phi()*TMath::RadToDeg();

    int reaction = 1;
    // reaction:
    // 0  p e+e- production   : 1/E x Bethe-Heidler process.
    // 1  p mu+mu- production   : 1/E x Bethe-Heidler process.
    // 2  p gamma             : 1/E x Total cross section.
    // 3  p pi0               : 1/E x Total cross section.
    // 4  p pi0pi0            : 1/E x Total cross section.
    // 5  p pi+pi-            : 1/E x Total cross section.

    TrueAnalysis_ll(reaction); //lepton-lepton analysis

    // Double_t protonE = TrueObs.GetTrueProtonLV().E();


    ThvT_p->FillWeighted( TrueObs.GetTrueProtonLV().E()*1000 - MASS_PROTON, TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight() );

    ThpvEg->FillWeighted( TrueObs.GetTrueBeamEnergy()*1000, TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight());

    ThvT_e->FillWeighted( TrueObs.GetTrueElectronLV().E()*1000 - MASS_ELECTRON, TrueObs.GetTrueElectronLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight() );
    ThvT_e->FillWeighted( TrueObs.GetTruePositronLV().E()*1000 - MASS_ELECTRON, TrueObs.GetTruePositronLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight() );

    ThvT_mu->FillWeighted( TrueObs.GetTrueMuonNegLV().E()*1000 - MASS_MUON, TrueObs.GetTrueMuonNegLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight() );
    ThvT_mu->FillWeighted( TrueObs.GetTrueMuonPosLV().E()*1000 - MASS_MUON, TrueObs.GetTrueMuonPosLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight() );


    Ntracks->FillWeighted(GetTracks()->GetNTracks(), TrueObs.GetWeight());

    for(UInt_t i = 0; i < GetTracks()->GetNTracks(); i++)
    {
        E_v_dE->FillWeighted( GetTracks()->GetClusterEnergy(i), GetTracks()->GetVetoEnergy(i), TrueObs.GetWeight());

        Double_t E    =  GetTracks()->GetClusterEnergy(i);
        if(cutProtonCB->IsInside(E, GetTracks()->GetVetoEnergy(i)))
        {
            Double_t E    =  GetTracks()->GetClusterEnergy(i);
            double p[5] = {-4.8794, 0.220078, -0.00402405, 3.27302e-05, -9.58383e-08};
            Double_t E_rec = E - E*(-4.8794 + 0.220078*E - 0.00402405*E*E + 3.27302e-05*E*E*E - 9.58383e-08*E*E*E*E) + 0.04*E;
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

void	ProtonChargeRadius::ProcessScalerRead()
{
//    hist_eta.ScalerReadCorrection(Double_t(scalers->GetScaler(0))/scalers->GetScaler(1));
//    hist_etap.ScalerReadCorrection(Double_t(scalers->GetScaler(0))/scalers->GetScaler(1));
}


Bool_t	ProtonChargeRadius::Init(const char* configfile)
{
    return kTRUE;
}


//void ProtonChargeRadius::TrueAnalysis_ll(int reaction)
//{
//    Double_t weight = 1.0;
//    Double_t norm;

////    norm = 50/163.9*5000/5001; //310 MeV beam energy norm
////    norm = 50000.0/99012.0; //510 MeV beam ebnergy norm
//    if(reaction == 4)
//        norm = 1.0;

//    weight = (1.0/TrueObs.GetTrueBeamEnergy())*norm;

//    if(reaction < 2)
//    {
//        Double_t ml;
//        if(reaction == 0)
//            ml = MASS_ELECTRON/1000;
//        else
//            ml = MASS_MUON/1000;

//        Double_t mll_sq; // mass squared lepton pairs
//        if( TMath::Abs(ml - MASS_ELECTRON/1000) < 1.0e-4) //electrons/positrons
//            mll_sq = (TrueObs.GetTrueElectronLV() + TrueObs.GetTruePositronLV()).M2();
//        else                                                    //muons
//            mll_sq = (TrueObs.GetTrueMuonNegLV() + TrueObs.GetTrueMuonPosLV()).M2();

//        Double_t mp         =   MASS_PROTON/1000;
//        Double_t pr_lab_mom =   TrueObs.GetTrueProtonLV().P();
//        Double_t pr_lab_mom2 =   TMath::Sqrt((TrueObs.GetTrueProtonLV().E()*TrueObs.GetTrueProtonLV().E()- MASS_PROTON*MASS_PROTON/1.0e6));
//        Double_t alpha      =   1.0/137.0;
//        Double_t beta       =   TMath::Sqrt( 1 - ( (4*ml*ml)/mll_sq ) );
//        Double_t s          =   (MASS_PROTON*MASS_PROTON + 2*MASS_PROTON*TrueObs.GetTrueBeamEnergy()*1000.0)/1.0e6;   // Mandelstam s
//        Double_t t          =   2*mp*mp -TMath::Sqrt(4*mp*mp*mp*mp + (pr_lab_mom*pr_lab_mom)*(4*mp*mp) );           // momentum transfer t
//        t_true = t;
//        Double_t tau        =   (-t)/(4*mp*mp);

//        Double_t term1, term2, term3;

//        term1               =   (alpha*alpha*alpha)/( (s-mp*mp)*(s-mp*mp) );
//        term2               =   (4*beta)/( t*t*TMath::Power((mll_sq-t),4) );
//        term3               =   1/(1+tau);

//        Double_t CE1, CE2, CM1, CM2;
//        Double_t CE, GEP, CM, GMP;

//        CE1                 =   (t)*(s - mp*mp)*(s - mp*mp - mll_sq + t)*( TMath::Power(mll_sq,2) + 6*mll_sq*t + t*t + 4*ml*ml*mll_sq ) + (TMath::Power((mll_sq -t),2))*( t*t*mll_sq + mp*mp*( TMath::Power((mll_sq+t),2) + 4*ml*ml*mp*mp*mll_sq));
//        CE2                 =  (-t)*(s - mp*mp)*(s - mp*mp - mll_sq + t)*( TMath::Power(mll_sq,2) + t*t + 4*ml*ml*(mll_sq+2*t-2*ml*ml))+ TMath::Power((mll_sq - t), 2)*(-1*mp*mp*(mll_sq*mll_sq + t*t ) + 2*ml*ml*(-1*t*t - 2*mp*mp*mll_sq + 4*ml*ml*mp*mp));
//        CE                  =  CE1 + CE2*(1/beta)*TMath::Log((1 + beta)/(1 - beta)) ;

//        CM1                 =  CE1 - 2*mp*mp*(1 + tau)*TMath::Power((mll_sq-t),2)*(mll_sq*mll_sq + t*t + 4*ml*ml*mll_sq);
//        CM2                 =  CE2 - 2*mp*mp*(1 + tau)*TMath::Power((mll_sq-t),2)*(mll_sq*mll_sq + t*t + 4*ml*ml*(mll_sq - t -2*ml*ml));
//        CM                  =  CM1 + CM2*(1/beta)*TMath::Log((1 + beta)/(1 - beta)) ;

//        GEP                 = 1.0; // In static limit GE(0) = 1
//        GMP                 = 1.0; // In static limit GM(0) = mu_p = 2.79* mu_N

//        Double_t diffxs     =   term1*term2*term3*(CE*GEP*GEP + CM*tau*GMP*GMP);

//}


void ProtonChargeRadius::TrueAnalysis_ll(int reaction)
{

    Double_t weight = 1.0;
    Double_t norm;
         if( reaction < 2 )
         {

            Double_t ml;
            if(reaction == 0)
                ml = MASS_ELECTRON/1000;
            else
                ml = MASS_MUON/1000;

            Double_t mll_sq; // mass squared lepton pairs
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
            Double_t t          =   2*mp*mp -TMath::Sqrt(4*mp*mp*mp*mp + (pr_lab_mom*pr_lab_mom)*(4*mp*mp) );           // momentum transfer t
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
            CM2                 =  CE2 - 2*mp*mp*(1 + tau)*TMath::Power((mll_sq-t),2)*(mll_sq*mll_sq + t*t + 4*ml*ml*(mll_sq - t -2*ml*ml));
            CM                  =  CM1 + CM2*(1/beta)*TMath::Log((1 + beta)/(1 - beta)) ;

            Double_t diffxs     =   term1*term2*term3*(CE*GEp(-t)*GEp(-t) + CM*tau*GMp(-t)*GMp(-t));

            // norm 510 MeV e and mu
                   if(reaction == 0)
                       norm = 1.0e-3*200000/38219*2.0/2.01099;
                   else
                       norm = 200000/13837*200000/198659*2.00000/1.98659;  // norm = 50000./15.11*2/1.985*2.0/1.99002;

            weight = diffxs*norm;
            Double_t weight_beame;
            weight_beame = 200000./397369.;


            weight = (1.0/TrueObs.GetTrueBeamEnergy())*weight*weight_beame;
                    if(weight < 0.1 )
                        weight = 0;
                    TrueObs.SetWeight(weight);


            dsigma_v_mll->FillWeighted(mll_sq,diffxs, weight);
            MC_weight->FillWeighted(weight,weight);
            mll->FillWeighted(mll_sq, TrueObs.GetWeight() );

            t_v_mll->FillWeighted(mll_sq, -t, TrueObs.GetWeight());
            thlabpr_v_mll->FillWeighted(mll_sq, TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg(),  TrueObs.GetWeight());

            thlabpr_v_t->FillWeighted(-t, TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg(),  TrueObs.GetWeight());
            dsigma_v_mll->FillWeighted(mll_sq,diffxs, TrueObs.GetWeight());
            mll->FillWeighted(mll_sq, TrueObs.GetWeight() );

            t_v_mll->FillWeighted(mll_sq, -t, TrueObs.GetWeight());
            thlabpr_v_mll->FillWeighted(mll_sq, TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg(),  TrueObs.GetWeight());

            thlabpr_v_t->FillWeighted(-t, TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg(),  TrueObs.GetWeight());

            if( (-t > 0.005) && (-t < 0.015) )
            {
                 thlabpr_v_Tp_13->FillWeighted( TrueObs.GetTrueProtonLV().E()*1000 - MASS_PROTON, TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight());
                 thlabpr_v_th_e_1->FillWeighted( TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg(), TrueObs.GetTrueElectronLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight());
                 thlabpr_v_th_e_1->FillWeighted( TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg(), TrueObs.GetTruePositronLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight());
                 thlabpr_v_th_mu_1->FillWeighted( TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg(), TrueObs.GetTrueMuonNegLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight());
                 thlabpr_v_th_mu_1->FillWeighted( TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg(), TrueObs.GetTrueMuonPosLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight());
                 thlabpr_v_Tp_1->FillWeighted( TrueObs.GetTrueProtonLV().E()*1000 - MASS_PROTON, TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg(),TrueObs.GetWeight());
             }
             if( (-t > 0.015) && (-t < 0.025) )
             {
                 thlabpr_v_Tp_13->FillWeighted( TrueObs.GetTrueProtonLV().E()*1000 - MASS_PROTON, TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight());
                 thlabpr_v_th_e_2->FillWeighted( TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg(), TrueObs.GetTrueElectronLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight());
                 thlabpr_v_th_e_2->FillWeighted( TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg(), TrueObs.GetTruePositronLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight());
                 thlabpr_v_th_mu_2->FillWeighted( TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg(), TrueObs.GetTrueMuonNegLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight());
                 thlabpr_v_th_mu_2->FillWeighted( TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg(), TrueObs.GetTrueMuonPosLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight());
                 thlabpr_v_Tp_2->FillWeighted( TrueObs.GetTrueProtonLV().E()*1000 - MASS_PROTON, TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg(),TrueObs.GetWeight());
             }
             if( (-t > 0.025) && (-t < 0.035) )
             {
                 thlabpr_v_Tp_13->FillWeighted( TrueObs.GetTrueProtonLV().E()*1000 - MASS_PROTON, TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight());
                 thlabpr_v_th_e_3->FillWeighted( TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg(), TrueObs.GetTrueElectronLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight());
                 thlabpr_v_th_e_3->FillWeighted( TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg(), TrueObs.GetTruePositronLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight());
                 thlabpr_v_th_mu_3->FillWeighted( TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg(), TrueObs.GetTrueMuonNegLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight());
                 thlabpr_v_th_mu_3->FillWeighted( TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg(), TrueObs.GetTrueMuonPosLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight());
                 thlabpr_v_Tp_3->FillWeighted( TrueObs.GetTrueProtonLV().E()*1000 - MASS_PROTON, TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg(),TrueObs.GetWeight());
             }
         }

         if(reaction == 3)
         {
             norm = 1.0;
             weight = (1.0/TrueObs.GetTrueBeamEnergy())*norm;
             TrueObs.SetWeight(weight);

         }

         if(reaction == 4)
         {
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



