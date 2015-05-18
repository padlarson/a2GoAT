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
    ThvT_p              = new GHistBGSub2("ThvT_p","#theta_{p} (y) vs Energy p (x)", 100, 0., 100, 100, 0., 100.);
    ThvT_e              = new GHistBGSub2("ThvT_e","#theta_{e} (y) vs Energy e^{+}/e^{-} (x)", 60, 0., 600, 18, 0., 180.);
    ThvT_mu             = new GHistBGSub2("ThvT_mu","#theta_{#mu} (y) vs Energy #mu^{+}/#mu^{-} (x)", 200, 0., 400, 90, 0., 180.);

    ThpvEg              = new GHistBGSub2("ThpvEg","#theta_{p} (y) vs Beam Energy #gamma (x)" , 100, 100., 600, 200, 0., 100.);

    dsigma_v_mll        = new GHistBGSub2("dsigma_v_mll", "d#sigma (y) M_{ll}^{2} (x) (#mub /GeV^{4}) ", 1000, 0.0, 0.2, 1000, 0, 1000);
    t_v_mll             = new GHistBGSub2("t_v_mll", "mom transfer (-t) (y) vs M_{ll}^{2} (x)", 50, 0.0, 0.2, 50, 0, 0.2);
    thlabpr_v_mll       = new GHistBGSub2("thlabpr_v_mll", "#theta_{lab, pr} (y) vs M_{ll}^{2} (x)", 100, 0.0, 0.2, 50, 0, 100);


    MC_weight           = new GH1("MC_weight", "MC_weight", 1000, 0, 10);
    mll                 = new GH1("mll", "M_{l^{+}l^{-}}", 100, 0.0, 0.20);


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

    SetAsPhysicsFile();

    TraverseValidEvents();

    outputFile->cd();
    gDirectory->mkdir("MC")->cd();


	return kTRUE;
}

void	ProtonChargeRadius::ProcessEvent()
{

    TrueObs.Start(*GetPluto(), *GetGeant());

    Double_t weight = 1.0;
    True_BeamEnergy->Fill( TrueObs.GetTrueBeamEnergy(), TrueObs.GetWeight() );

//    int reaction = 0;
    // reaction:
    // 0  p e+e- production   : 1/E x Bethe-Heidler process.
    // 1  p mu+mu- production   : 1/E x Bethe-Heidler process.
    // 2  p gamma             : 1/E x Total cross section.
    // 3  p pi0               : 1/E x Total cross section.
    // 4  p pi0pi0            : 1/E x Total cross section.
    // 5  p pi+pi-            : 1/E x Total cross section.

//    TrueAnalysis_ll(reaction); //lepton-lepton analysis

    // Double_t protonE = TrueObs.GetTrueProtonLV().E();

    ThvT_p->FillWeighted( TrueObs.GetTrueProtonLV().E()*1000 - MASS_PROTON, TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight() );

    ThpvEg->FillWeighted( TrueObs.GetTrueBeamEnergy()*1000, TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight());

    ThvT_e->FillWeighted( TrueObs.GetTrueElectronLV().E()*1000 - MASS_ELECTRON, TrueObs.GetTrueElectronLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight() );
    ThvT_e->FillWeighted( TrueObs.GetTruePositronLV().E()*1000 - MASS_ELECTRON, TrueObs.GetTruePositronLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight() );

    ThvT_mu->FillWeighted( TrueObs.GetTrueMuonNegLV().E()*1000 - MASS_MUON, TrueObs.GetTrueMuonNegLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight() );
    ThvT_mu->FillWeighted( TrueObs.GetTrueMuonPosLV().E()*1000 - MASS_MUON, TrueObs.GetTrueMuonPosLV().Theta()*TMath::RadToDeg(), TrueObs.GetWeight() );

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


void ProtonChargeRadius::TrueAnalysis_ll(int reaction)
{
    Double_t weight = 1.0;
    Double_t norm;

    norm = 50/163.9*5000/5001; //310 MeV beam energy norm
//    norm = 50000.0/99012.0; //510 MeV beam ebnergy norm
    weight = (1.0/TrueObs.GetTrueBeamEnergy())*norm;

    if(reaction < 2)
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
        Double_t pr_lab_mom =   TrueObs.GetTrueProtonLV().P();
        Double_t alpha      =   1.0/137.0;
        Double_t beta       =   TMath::Sqrt( 1 - ( (4*ml*ml)/mll_sq ) );
        Double_t s          =   (MASS_PROTON*MASS_PROTON + 2*MASS_PROTON*TrueObs.GetTrueBeamEnergy()*1000.0)/1.0e6;   // Mandelstam s
        Double_t t          =   2*mp*mp -TMath::Sqrt(4*mp*mp*mp*mp + (pr_lab_mom*pr_lab_mom)*(4*mp*mp) );           // momentum transfer t
        Double_t tau        =   (-t)/(4*mp*mp);

        Double_t term1, term2, term3;

        term1               =   (alpha*alpha*alpha)/( (s-mp*mp)*(s-mp*mp) );
        term2               =   (4*beta)/( t*t*TMath::Power((mll_sq-t),4) );
        term3               =   1/(1+tau);

        Double_t CE1, CE2, CM1, CM2;
        Double_t CE, GEP, CM, GMP;

        CE1                 =   (t)*(s - mp*mp)*(s - mp*mp - mll_sq + t)*( TMath::Power(mll_sq,2) + 6*mll_sq*t + t*t + 4*ml*ml*mll_sq ) + (TMath::Power((mll_sq -t),2))*( t*t*mll_sq + mp*mp*( TMath::Power((mll_sq+t),2) + 4*ml*ml*mp*mp*mll_sq));
        CE2                 =  (-t)*(s - mp*mp)*(s - mp*mp - mll_sq + t)*( TMath::Power(mll_sq,2) + t*t + 4*ml*ml*(mll_sq+2*t-2*ml*ml))+ TMath::Power((mll_sq - t), 2)*(-1*mp*mp*(mll_sq*mll_sq + t*t ) + 2*ml*ml*(-1*t*t - 2*mp*mp*mll_sq + 4*ml*ml*mp*mp));
        CE                  =  CE1 + CE2*(1/beta)*TMath::Log((1 + beta)/(1 - beta)) ;

        CM1                 =  CE1 - 2*mp*mp*(1 + tau)*TMath::Power((mll_sq-t),2)*(mll_sq*mll_sq + t*t + 4*ml*ml*mll_sq);
        CM2                 =  CE2 - 2*mp*mp*(1 + tau)*TMath::Power((mll_sq-t),2)*(mll_sq*mll_sq + t*t + 4*ml*ml*(mll_sq - t -2*ml*ml));
        CM                  =  CM1 + CM2*(1/beta)*TMath::Log((1 + beta)/(1 - beta)) ;

        GEP                 = 1.0; // In static limit GE(0) = 1
        GMP                 = 1.0; // In static limit GM(0) = mu_p = 2.79* mu_N

        Double_t diffxs     =   term1*term2*term3*(CE*GEP*GEP + CM*tau*GMP*GMP);
// norm 310 MeV e and mu
        if(reaction == 0)
            norm = 50000./143500.;
        else
            norm = 50000./12.56;


// norm 510 MeV e and mu
//       if(reaction == 0)
//            norm = 50000./10209.*5000.0/5046.4;
//        else
//            norm = 50000./15.11;

        weight = weight*diffxs*norm;
        if(weight < 0.1 )
            weight = 0;
        TrueObs.SetWeight(weight);

        Double_t test = TrueObs.GetWeight();

        dsigma_v_mll->FillWeighted(mll_sq,diffxs, TrueObs.GetWeight());
        mll->FillWeighted(mll_sq, TrueObs.GetWeight() );

        t_v_mll->FillWeighted(mll_sq, -t, TrueObs.GetWeight());
        thlabpr_v_mll->FillWeighted(mll_sq, TrueObs.GetTrueProtonLV().Theta()*TMath::RadToDeg(),  TrueObs.GetWeight());


    }

    MC_weight->FillWeighted(weight,weight);


//        Double_t ml;
//        Double_t mll_sq; // mass squared lepton pairs


//        for(Int_t k = 0; k < 2; k++)
//        {
//            if(k == 0)
//                ml = MASS_ELECTRON/1000;
//            else
//               ml =  MASS_MUON/1000;
//            for(Int_t j = 0; j < 3 ;j++)
//            {
//                Double_t Eg         =   0.5000;
//                Double_t t          =   -0.03+j*(0.01);

//                Double_t mp         =   MASS_PROTON/1000;

//                Double_t alpha      =   1.0/137.0;

//                for(int i = 0; i < 100; i++)
//                {
//                    mll_sq = 0.03 + i*0.0005;

//                    Double_t beta       =   TMath::Sqrt( 1 - ( (4*ml*ml)/mll_sq ) );
//                    Double_t s          =   (MASS_PROTON*MASS_PROTON + 2*MASS_PROTON*Eg*1000.0)/1.0e6;   // Mandelstam s

//                    Double_t tau        =   (-t)/(4*mp*mp);

//                    Double_t term1, term2, term3;

//                    term1               =   ( TMath::Power(alpha,3) )/( TMath::Power((s-mp*mp),2) );
//                    term2               =   (4*beta)/( t*t*TMath::Power((mll_sq-t),4) );
//                    term3               =   1/(1+tau);

//                    Double_t CE1, CE2, CM1, CM2;
//                    Double_t CE, GEP, CM, GMP;

//                    CE1                 =   (t)*(s - mp*mp)*(s - mp*mp - mll_sq + t)*( TMath::Power(mll_sq,2) + 6*mll_sq*t + t*t + 4*ml*ml*mll_sq ) + (TMath::Power((mll_sq -t),2))*( t*t*mll_sq + mp*mp*( TMath::Power((mll_sq+t),2) + 4*ml*ml*mp*mp*mll_sq));
//                    CE2                 =  (-t)*(s - mp*mp)*(s - mp*mp - mll_sq + t)*( TMath::Power(mll_sq,2) + t*t + 4*ml*ml*(mll_sq + 2*t - 2*ml*ml))+ TMath::Power((mll_sq - t), 2)*(-1*mp*mp*( mll_sq*mll_sq + t*t ) + 2*ml*ml*(-1*t*t - 2*mp*mp*mll_sq + 4*ml*ml*mp*mp));
//                    CE                  =  CE1 + CE2*(1/beta)*TMath::Log((1 + beta)/(1 - beta)) ;

//                    CM1                 =  CE1 - 2*mp*mp*(1 + tau)*TMath::Power((mll_sq-t),2)*(mll_sq*mll_sq + t*t + 4*ml*ml*mll_sq);
//                    CM2                 =  CE2 + 2*mp*mp*(1 + tau)*TMath::Power((mll_sq-t),2)*(mll_sq*mll_sq + t*t + 4*ml*ml*(mll_sq - t -2*ml*ml));
//                    CM                  =  CM1 + CM2*(1/beta)*TMath::Log((1 + beta)/(1 - beta)) ;

//                    GEP                 = 1.0; // In static limit GE(0) = 1
//                    GMP                 = 1.0; // In static limit GM(0) = mu_p = 2.79* mu_N

//                    Double_t diffxs     =   term1*term2*term3*(CE*GEP*GEP + CM*tau*GMP*GMP)*100*3.654;



//                    dsigma_v_mll->Fill(mll_sq,diffxs);
//                }
//            }
//        }



}






