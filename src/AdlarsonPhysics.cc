#include "GParticleReconstruction.h"
#include "AdlarsonPhysics.h"
#include "GTrue.h"



AdlarsonPhysics::AdlarsonPhysics()
{
// TRUE OBSERVABLES
// Beam energy
    True_BeamEnergy     = new GH1("True_BeamEnergy", "True Beam Energy", 250, 1.400, 1.650);
    Tagged_BeamEnergy    = new GH1("Tagged_BeamEnergy", "Tagged Beam Energy", 1650, 0, 1650);
// Phase space final state particles
    ThpvsThetaprCM      = new GHistBGSub2("ThpvsThetaprCM", "#theta_{p} vs #theta_{#eta^{'}}", 360, 0., 180, 360, 0., 180.);
    ThvE_p              = new GHistBGSub2("ThvE_p", "True E_{p} vs #theta_{p}", 100, 0., 0.6, 100, 0., 25.);
    ThvE_eta_g          = new GHistBGSub2("ThvE_eta_g", "E_{#gamma, #eta} vs #theta_{#gamma, #eta}", 100, 0, 1.000, 36, 0, 180);
    ThvE_pi0_g          = new GHistBGSub2("ThvE_pi0_g", "E_{#gamma, #pi^{0}} vs #theta_{#gamma, #pi^{0}}", 60, 0, .600, 36, 0, 180);
// Physics result
    DP_true             = new GHistBGSub2("DP_true", "True Dalitz Plot distribution", 60, -1.5, 1.5, 60, -1.5, 1.5);
    M_pi1pi2_true       = new GH1("M_pi1pi2_true", "True M_{#pi#pi,true}^{2}", 60, 0.05, 0.20);
    M_etapi_true        = new GH1("M_etapi_true", "True M_{#eta#pi,true}^{2}", 100, 0.45, 0.70);

// RECONSTRUCTED OBSERVABLES
// Rec. TAPS - proton analysis
    EvdE_TAPS_all      = new GHistBGSub2("EvdE_TAPS_all", "All E vs dE TAPS", 1200, 0., 600., 200, 0.5, 10.);
    EvdE_TAPS_proton   = new GHistBGSub2("EvdE_TAPS_proton", "Best proton cand E vs dE TAPS", 1200, 0., 600, 200, 0.5, 10.);

    Nrprotons           = new GHistBGSub("Nrprotons", "nr of protons", 5,  0, 5);
    ThvEp_rec           = new GHistBGSub2("ThvEp_rec", "Rec E_{p} vs #theta_{p}", 120, 0., 600, 50, 0., 25.);
    MM_p                = new GHistBGSub("MM_p", "Missing Mass calculated for proton", 300, 800., 1100.);

// Rec. Photons
    IM_6g	= new GH1("IM_6g", 	"IM_6g/7g", 240,   200, 1400);
    IM_10g	= new GH1("IM_10g", "IM_10g", 	240,   200, 1400);

    IM6gvMMp = new GHistBGSub2("IM6gvMMp", "MM(p) vs IM(6#gamma/7#gamma)", 300,800., 1100., 240, 200., 1400.);
    IM10gvMMp = new GHistBGSub2("IM10gvMMp", "MM(p) vs IM(10#gamma)", 300,800., 1100., 240, 200., 1400.);

    cutFile             = new TFile("/home/adlarson/a2GoAT/configfiles/cuts/TAPS_DeltaE-E-PA.root");
    cutProtonTAPS       = (TCutG*)cutFile->Get("CutProton");

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

    etapr_6gTrue.Start(*GetPluto(), *GetGeant());   // (pluto tree, n part in pluto per event)
    TrueAnalysis_etapr6g();                         // obtains the true observables


        Int_t nrprotons = 0;
        Int_t iprtrack = -1;

        MMp_vec.SetPxPyPzE(0., 0., 0., 0.);
        IM6g_vec.SetPxPyPzE(0., 0., 0., 0.);
        IM10g_vec.SetPxPyPzE(0., 0., 0., 0.);

        for (int i = 0; i < GetTracks()->GetNTracks() ; i++)
        {            
            if( GetTracks()->GetApparatus(i) == GTreeTrack::APPARATUS_TAPS )
            {
               EvdE_TAPS_all->Fill(GetTracks()->GetClusterEnergy(i),GetTracks()->GetVetoEnergy(i));
               if(cutProtonTAPS->IsInside(GetTracks()->GetClusterEnergy(i),GetTracks()->GetVetoEnergy(i)))
               {
                    nrprotons++;
                    iprtrack = i;

                    if( nrprotons > 1 )
                    {
                      if( GetTracks()->GetVetoEnergy(iprtrack) < GetTracks()->GetVetoEnergy(i) )
                        iprtrack = i;
                    }
                    else
                        iprtrack = i;

               }
            }
        }

        Nrprotons->Fill(nrprotons);
        if( iprtrack == -1 ) return;

        EvdE_TAPS_proton->Fill(GetTracks()->GetClusterEnergy(iprtrack),GetTracks()->GetVetoEnergy(iprtrack));
        ThvEp_rec->Fill(GetTracks()->GetClusterEnergy(iprtrack),GetTracks()->GetTheta(iprtrack));


        proton_vec = GetTracks()->GetVector(iprtrack, pdgDB->GetParticle("proton")->Mass()*1000);
        // Now construct missing mass calc for proton with tagger energies.
        for ( Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
        {
            Tagged_BeamEnergy->Fill(GetTagger()->GetTaggedEnergy(i));
            double t = GetTagger()->GetTaggedEnergy(i);

            MMp_vec.SetPxPyPzE(-proton_vec.Px(), -proton_vec.Py(), GetTagger()->GetTaggedEnergy(i) - proton_vec.Pz(), GetTagger()->GetTaggedEnergy(i) + pdgDB->GetParticle("proton")->Mass()*1000 - proton_vec.E());

            MMp = ( MMp_vec ).M();
            MM_p->Fill(MMp);
        }


        if( (GetTracks()->GetNTracks() == 7) || (GetTracks()->GetNTracks() == 8) || (GetTracks()->GetNTracks() == 9) )
        {
            sixgAnalysis(iprtrack);
        }

        if(GetPhotons()->GetNParticles() == 10)
        {
            tengAnalysis(iprtrack);
        }



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



void AdlarsonPhysics::TrueAnalysis_etapr6g()
{
    True_BeamEnergy->Fill(etapr_6gTrue.GetTrueBeamEnergy());

    TLorentzVector k;
    k.SetPxPyPzE(0.,0.,0.,0.);

    for(Int_t i = 0; i < etapr_6gTrue.GetNgamma(); i++ )
        k += etapr_6gTrue.GetTrueGammaLV(i);

    double testtrue = k.M();

    double hej = 0;

    // calculate Physics

    etapr_true[0] = etapr_6gTrue.GetTrueEtaLV();
    etapr_true[1] = etapr_6gTrue.GetTrueNeutralPiLV(0);
    etapr_true[2] = etapr_6gTrue.GetTrueNeutralPiLV(1);
    DalitzPlot(etapr_true, Xtrue, Ytrue, DPnrTrue);
    m2pi0_metapi0(etapr_true, m_etapi01True, m_etapi02True, m_2pi0True);

    // calculate True Energy vs Theta of final state particles

    ThvE_p->Fill(etapr_6gTrue.GetTrueProtonLV().E() - MASS_PROTON/1000 ,etapr_6gTrue.GetTrueProtonLV().Theta()*TMath::RadToDeg());

    for(Int_t i = 0; i < etapr_6gTrue.GetNgamma(); i++ )
    {
        if ( i > 3 ) // first four gammas come from 2pi0
            ThvE_eta_g->Fill(etapr_6gTrue.GetTrueGammaLV(i).E(),etapr_6gTrue.GetTrueGammaLV(i).Theta()*TMath::RadToDeg());
        else
            ThvE_pi0_g->Fill(etapr_6gTrue.GetTrueGammaLV(i).E(),etapr_6gTrue.GetTrueGammaLV(i).Theta()*TMath::RadToDeg());
    }
    return;
}


void AdlarsonPhysics::sixgAnalysis(Int_t ipr)
{

    for ( Int_t i = 0; i <= GetTracks()->GetNTracks() ; i++ )
    {
        if(i != ipr)
            IM6g_vec += GetTracks()->GetVector(i);
    }

    IM_6g->Fill( IM6g_vec.M() );
    IM6gvMMp->Fill(MMp,IM6g_vec.M());

    for( Int_t i = 0; i < GetTagger()->GetNTagged(); i++ )
    {

        // Run kinfit:
        // For all (6g or 7g) + proton obtain errors for all observables
        // Run kinfit and select the tagger + 6g + proton with lowest chi2

        // Here run post-kinfit analysis:
        // Find the 15 combinations to form a 2pi0+eta and 3pi0. Determin which reaction it is and remove 3pi0
        // run kinfit with rec. observables now also enforcing eta mass.

        // Obtain physics observables

    }

    // select best combination eta 2pi0 and 3pi0:
    // GetBest6gCombination();

    // cut in the scatterplot eta 2pi0 3pi0 pdf:

    // run kinfit again now enforcing eta mass (and pi0 pi0 as well?)

    // final results: Dalitz plot, m_pi0pi0, eta prime production as fcn energy
    // DalitzPlot();
    // m2pi0_metapi0;

}

void AdlarsonPhysics::tengAnalysis(Int_t ipr)
{
    for ( Int_t i = 0; i < GetTracks()->GetNTracks() ; i++ )
    {
        if(i != iprtrack)
            IM10g_vec += GetTracks()->GetVector(i);
    }

    IM_10g->Fill(IM10g_vec.M());
    IM10gvMMp->Fill( MMp , IM10g_vec.M());

    for( Int_t i = 0; i < GetTagger()->GetNTagged(); i++ )
    {


    }

    // with the best candidate:
    // select best combination eta (6g out of 10g, 220 combinations)
    // cut in the scatterplot eta 2pi0 3pi0 pdf

    // run kinfit again now enforcing eta mass (and pi0 pi0 as well?)

    // final results: Dalitz plot, m_pi0pi0, eta prime production as fcn energy
    //
    // DalitzPlot();
    // m2pi0_metapi0;
    //

}


void AdlarsonPhysics::DalitzPlot(const TLorentzVector g[3] , Double_t &X, Double_t &Y , Int_t &DP_nr)
{

    // g is the array of TLorentzVectors containing the final event sample in order
    // 0 1 2 - eta pi01 pi02
    // Dalitz plot variables

    Double_t T_eta, T_pi1, T_pi2;
    Double_t Q;
    Double_t X_max = 1.5, Y_max = 1.5, bin_width = 0.2;

    TVector3 etaprime_rest      = -(g[0]+g[1]+g[2]).BoostVector();
    TLorentzVector eta_ep_rest  = g[0];
    TLorentzVector pi01_ep_rest = g[1];
    TLorentzVector pi02_ep_rest = g[2];

    eta_ep_rest.Boost(etaprime_rest);
    pi01_ep_rest.Boost(etaprime_rest);
    pi02_ep_rest.Boost(etaprime_rest);

    T_eta = eta_ep_rest.E() - eta_ep_rest.M();
    if( pi01_ep_rest.E() > pi02_ep_rest.E() ){
        T_pi1 = pi01_ep_rest.E() - pi01_ep_rest.M();
        T_pi2 = pi02_ep_rest.E() - pi02_ep_rest.M();
    }
    else{
        T_pi1 = pi02_ep_rest.E() - pi02_ep_rest.M();
        T_pi2 = pi01_ep_rest.E() - pi01_ep_rest.M();
    }

    Q = T_eta + T_pi1 + T_pi2;
    X = ( TMath::Sqrt(3.0) / Q ) * ( T_pi1 - T_pi2 );
    Y = ( MASS_ETA + 2 * MASS_PI0 ) / MASS_PI0 * ( T_eta / Q ) - 1.0;

    DP_nr = ( int( X ) + X_max ) / 0.2 + 40*( int( Y ) + Y_max );

    DP_true->Fill(X,Y);
}


Double_t    AdlarsonPhysics::IM_Ng( UInt_t n )
{
    TLorentzVector g_vec;

    g_vec.SetPxPyPzE(0,0,0,0);
    for ( Int_t i = 0 ; i < n; i++ )
        g_vec += GetPhotons()->Particle(i);

    return g_vec.M();
}

void AdlarsonPhysics::m2pi0_metapi0(  TLorentzVector g[3], Double_t &m_etapi01, Double_t &m_etapi02, Double_t &m_2pi0 )
{
    // g is the array of TLorentzVectors containing the final event sample in order
    // 0 1 2 - eta pi01 pi02
    // Dalitz plot variables

    Double_t meta, mpi01, mpi02, metapr;


    meta = g[0].M();
    mpi01 = g[1].M();
    mpi02 = g[2].M();
    metapr = ( g[0] + g[1] + g[2] ).M();

    m_etapi01   = ( g[0] + g[1] ).M2();
    m_etapi02   = ( g[0] + g[2] ).M2();
    m_2pi0      = ( g[1] + g[2] ).M2();

    M_pi1pi2_true->Fill(m_2pi0);
    M_etapi_true->Fill(m_etapi01);
    M_etapi_true->Fill(m_etapi02);
}




Int_t AdlarsonPhysics::perm6g[15][6]=
{
    {0,1, 2,3, 4,5},
    {0,1, 2,4, 3,5},
    {0,1, 2,5, 4,3},

    {0,2, 1,3, 4,5},
    {0,2, 1,4, 3,5},
    {0,2, 1,5, 4,3},

    {0,3, 2,1, 4,5},
    {0,3, 2,4, 1,5},
    {0,3, 2,5, 4,1},

    {0,4, 2,3, 1,5},
    {0,4, 2,1, 3,5},
    {0,4, 2,5, 1,3},

    {0,5, 2,3, 4,1},
    {0,5, 2,4, 3,1},
    {0,5, 2,1, 4,3}
};

Int_t AdlarsonPhysics::perm6outof7g[7][6]=
{
    { 0, 1, 2, 3, 4, 5 },
    { 0, 1, 2, 3, 4, 6 },
    { 0, 1, 2, 3, 5, 6 },
    { 0, 1, 2, 4, 5, 6 },
    { 0, 1, 3, 4, 5, 6 },
    { 0, 2, 3, 4, 5, 6 },
    { 1, 2, 3, 4, 5, 6 }
};

Int_t AdlarsonPhysics::perm6outof8g[28][6]=
{
    { 0, 1, 2, 3, 4, 5 },   { 0, 1, 2, 3, 4, 6 },   { 0, 1, 2, 3, 4, 7 },
    { 0, 1, 2, 3, 5, 6 },   { 0, 1, 2, 3, 5, 7 },
    { 0, 1, 2, 3, 6, 7 },

    { 0, 1, 2, 4, 5, 6 },   { 0, 1, 2, 4, 5, 7 },
    { 0, 1, 2, 4, 6, 7 },

    { 0, 1, 2, 5, 6, 7 },

    { 0, 1, 3, 4, 5, 6 },   { 0, 1, 3, 4, 5, 7 },
    { 0, 1, 3, 4, 6, 7 },

    { 0, 1, 3, 5, 6, 7 },
    { 0, 1, 4, 5, 6, 7 },

    { 0, 2, 3, 4, 5, 6 },   { 0, 2, 3, 4, 5, 7 },
    { 0, 2, 3, 4, 6, 7 },
    { 0, 2, 3, 5, 6, 7 },
    { 0, 2, 4, 5, 6, 7 },

    { 0, 3, 4, 5, 6, 7 },

    { 1, 2, 3, 4, 5, 6 },   { 1, 2, 3, 4, 5, 7 },
    { 1, 2, 3, 4, 6, 7 },
    { 1, 2, 3, 5, 6, 7 },
    { 1, 2, 4, 5, 6, 7 },
    { 1, 3, 4, 5, 6, 7 },

    { 2, 3, 4, 5, 6, 7 }
};

Int_t AdlarsonPhysics::perm6outof10g[210][6]=
{
    // all possible permutations of 6 from 10 possibilities
    { 0, 1, 2, 3, 4, 5 }, 	{ 0, 1, 2, 3, 4, 6 }, 	{ 0, 1, 2, 3, 4, 7 }, 	{ 0, 1, 2, 3, 4, 8 }, 	{ 0, 1, 2, 3, 4, 9 },
    { 0, 1, 2, 3, 5, 6 }, 	{ 0, 1, 2, 3, 5, 7 }, 	{ 0, 1, 2, 3, 5, 8 }, 	{ 0, 1, 2, 3, 5, 9 },
    { 0, 1, 2, 3, 6, 7 }, 	{ 0, 1, 2, 3, 6, 8 }, 	{ 0, 1, 2, 3, 6, 9 },
    { 0, 1, 2, 3, 7, 8 }, 	{ 0, 1, 2, 3, 7, 9 },
    { 0, 1, 2, 3, 8, 9 },
    { 0, 1, 2, 4, 5, 6 }, 	{ 0, 1, 2, 4, 5, 7 }, 	{ 0, 1, 2, 4, 5, 8 }, 	{ 0, 1, 2, 4, 5, 9 },
    { 0, 1, 2, 4, 6, 7 }, 	{ 0, 1, 2, 4, 6, 8 }, 	{ 0, 1, 2, 4, 6, 9 },
    { 0, 1, 2, 4, 7, 8 }, 	{ 0, 1, 2, 4, 7, 9 },
    { 0, 1, 2, 4, 8, 9 },
    { 0, 1, 2, 5, 6, 7 }, 	{ 0, 1, 2, 5, 6, 8 }, 	{ 0, 1, 2, 5, 6, 9 },
    { 0, 1, 2, 5, 7, 8 }, 	{ 0, 1, 2, 5, 7, 9 },
    { 0, 1, 2, 5, 8, 9 },
    { 0, 1, 2, 6, 7, 8 }, 	{ 0, 1, 2, 6, 7, 9 },
    { 0, 1, 2, 6, 8, 9 },
    { 0, 1, 2, 7, 8, 9 },
    { 0, 1, 3, 4, 5, 6 }, 	{ 0, 1, 3, 4, 5, 7 }, 	{ 0, 1, 3, 4, 5, 8 }, 	{ 0, 1, 3, 4, 5, 9 },
    { 0, 1, 3, 4, 6, 7 }, 	{ 0, 1, 3, 4, 6, 8 }, 	{ 0, 1, 3, 4, 6, 9 },
    { 0, 1, 3, 4, 7, 8 }, 	{ 0, 1, 3, 4, 7, 9 },
    { 0, 1, 3, 4, 8, 9 },
    { 0, 1, 3, 5, 6, 7 }, 	{ 0, 1, 3, 5, 6, 8 }, 	{ 0, 1, 3, 5, 6, 9 },
    { 0, 1, 3, 5, 7, 8 }, 	{ 0, 1, 3, 5, 7, 9 },
    { 0, 1, 3, 5, 8, 9 },
    { 0, 1, 3, 6, 7, 8 }, 	{ 0, 1, 3, 6, 7, 9 },
    { 0, 1, 3, 6, 8, 9 },
    { 0, 1, 3, 7, 8, 9 },
    { 0, 1, 4, 5, 6, 7 }, 	{ 0, 1, 4, 5, 6, 8 }, 	{ 0, 1, 4, 5, 6, 9 },
    { 0, 1, 4, 5, 7, 8 }, 	{ 0, 1, 4, 5, 7, 9 },
    { 0, 1, 4, 5, 8, 9 },
    { 0, 1, 4, 6, 7, 8 }, 	{ 0, 1, 4, 6, 7, 9 },
    { 0, 1, 4, 6, 8, 9 },
    { 0, 1, 4, 7, 8, 9 },
    { 0, 1, 5, 6, 7, 8 }, 	{ 0, 1, 5, 6, 7, 9 },
    { 0, 1, 5, 6, 8, 9 },
    { 0, 1, 5, 7, 8, 9 },
    { 0, 1, 6, 7, 8, 9 },

    { 0, 2, 3, 4, 5, 6 }, 	{ 0, 2, 3, 4, 5, 7 }, 	{ 0, 2, 3, 4, 5, 8 }, 	{ 0, 2, 3, 4, 5, 9 },
    { 0, 2, 3, 4, 6, 7 }, 	{ 0, 2, 3, 4, 6, 8 }, 	{ 0, 2, 3, 4, 6, 9 },
    { 0, 2, 3, 4, 7, 8 }, 	{ 0, 2, 3, 4, 7, 9 },
    { 0, 2, 3, 4, 8, 9 },
    { 0, 2, 3, 5, 6, 7 }, 	{ 0, 2, 3, 5, 6, 8 }, 	{ 0, 2, 3, 5, 6, 9 },
    { 0, 2, 3, 5, 7, 8 }, 	{ 0, 2, 3, 5, 7, 9 },
    { 0, 2, 3, 5, 8, 9 },
    { 0, 2, 3, 6, 7, 8 }, 	{ 0, 2, 3, 6, 7, 9 },
    { 0, 2, 3, 6, 8, 9 },
    { 0, 2, 3, 7, 8, 9 },
    { 0, 2, 4, 5, 6, 7 }, 	{ 0, 2, 4, 5, 6, 8 }, 	{ 0, 2, 4, 5, 6, 9 },
    { 0, 2, 4, 5, 7, 8 }, 	{ 0, 2, 4, 5, 7, 9 },
    { 0, 2, 4, 5, 8, 9 },
    { 0, 2, 4, 6, 7, 8 }, 	{ 0, 2, 4, 6, 7, 9 },
    { 0, 2, 4, 6, 8, 9 },
    { 0, 2, 4, 7, 8, 9 },
    { 0, 2, 5, 6, 7, 8 }, 	{ 0, 2, 5, 6, 7, 9 },
    { 0, 2, 5, 6, 8, 9 },
    { 0, 2, 5, 7, 8, 9 },
    { 0, 2, 6, 7, 8, 9 },

    { 0, 3, 4, 5, 6, 7 }, 	{ 0, 3, 4, 5, 6, 8 }, 	{ 0, 3, 4, 5, 6, 9 },
    { 0, 3, 4, 5, 7, 8 }, 	{ 0, 3, 4, 5, 7, 9 },
    { 0, 3, 4, 5, 8, 9 },
    { 0, 3, 4, 6, 7, 8 }, 	{ 0, 3, 4, 6, 7, 9 },
    { 0, 3, 4, 6, 8, 9 },
    { 0, 3, 4, 7, 8, 9 },
    { 0, 3, 5, 6, 7, 8 }, 	{ 0, 3, 5, 6, 7, 9 },
    { 0, 3, 5, 6, 8, 9 },
    { 0, 3, 5, 7, 8, 9 },
    { 0, 3, 6, 7, 8, 9 },
    { 0, 4, 5, 6, 7, 8 }, 	{ 0, 4, 5, 6, 7, 9 },
    { 0, 4, 5, 6, 8, 9 },
    { 0, 4, 5, 7, 8, 9 },
    { 0, 4, 6, 7, 8, 9 },
    { 0, 5, 6, 7, 8, 9 },

    { 1, 2, 3, 4, 5, 6 }, 	{ 1, 2, 3, 4, 5, 7 }, 	{ 1, 2, 3, 4, 5, 8 }, 	{ 1, 2, 3, 4, 5, 9 },
    { 1, 2, 3, 4, 6, 7 }, 	{ 1, 2, 3, 4, 6, 8 }, 	{ 1, 2, 3, 4, 6, 9 },
    { 1, 2, 3, 4, 7, 8 }, 	{ 1, 2, 3, 4, 7, 9 },
    { 1, 2, 3, 4, 8, 9 },
    { 1, 2, 3, 5, 6, 7 }, 	{ 1, 2, 3, 5, 6, 8 }, 	{ 1, 2, 3, 5, 6, 9 },
    { 1, 2, 3, 5, 7, 8 }, 	{ 1, 2, 3, 5, 7, 9 },
    { 1, 2, 3, 5, 8, 9 },
    { 1, 2, 3, 6, 7, 8 }, 	{ 1, 2, 3, 6, 7, 9 },
    { 1, 2, 3, 6, 8, 9 },
    { 1, 2, 3, 7, 8, 9 },
    { 1, 2, 4, 5, 6, 7 }, 	{ 1, 2, 4, 5, 6, 8 }, 	{ 1, 2, 4, 5, 6, 9 },
    { 1, 2, 4, 5, 7, 8 }, 	{ 1, 2, 4, 5, 7, 9 },
    { 1, 2, 4, 5, 8, 9 },
    { 1, 2, 4, 6, 7, 8 }, 	{ 1, 2, 4, 6, 7, 9 },
    { 1, 2, 4, 6, 8, 9 },
    { 1, 2, 4, 7, 8, 9 },
    { 1, 2, 5, 6, 7, 8 }, 	{ 1, 2, 5, 6, 7, 9 },
    { 1, 2, 5, 6, 8, 9 },
    { 1, 2, 5, 7, 8, 9 },
    { 1, 2, 6, 7, 8, 9 },
    { 1, 3, 4, 5, 6, 7 }, 	{ 1, 3, 4, 5, 6, 8 }, 	{ 1, 3, 4, 5, 6, 9 },
    { 1, 3, 4, 5, 7, 8 }, 	{ 1, 3, 4, 5, 7, 9 },
    { 1, 3, 4, 5, 8, 9 },
    { 1, 3, 4, 6, 7, 8 }, 	{ 1, 3, 4, 6, 7, 9 },
    { 1, 3, 4, 6, 8, 9 },
    { 1, 3, 4, 7, 8, 9 },
    { 1, 3, 5, 6, 7, 8 }, 	{ 1, 3, 5, 6, 7, 9 },
    { 1, 3, 5, 6, 8, 9 },
    { 1, 3, 5, 7, 8, 9 },
    { 1, 3, 6, 7, 8, 9 },
    { 1, 4, 5, 6, 7, 8 }, 	{ 1, 4, 5, 6, 7, 9 },
    { 1, 4, 5, 6, 8, 9 },
    { 1, 4, 5, 7, 8, 9 },
    { 1, 4, 6, 7, 8, 9 },
    { 1, 5, 6, 7, 8, 9 },

    { 2, 3, 4, 5, 6, 7 }, 	{ 2, 3, 4, 5, 6, 8 }, 	{ 2, 3, 4, 5, 6, 9 },
    { 2, 3, 4, 5, 7, 8 }, 	{ 2, 3, 4, 5, 7, 9 },
    { 2, 3, 4, 5, 8, 9 },
    { 2, 3, 4, 6, 7, 8 }, 	{ 2, 3, 4, 6, 7, 9 },
    { 2, 3, 4, 6, 8, 9 },
    { 2, 3, 4, 7, 8, 9 },
    { 2, 3, 5, 6, 7, 8 }, 	{ 2, 3, 5, 6, 7, 9 },
    { 2, 3, 5, 6, 8, 9 },
    { 2, 3, 5, 7, 8, 9 },
    { 2, 3, 6, 7, 8, 9 },
    { 2, 4, 5, 6, 7, 8 }, 	{ 2, 4, 5, 6, 7, 9 },
    { 2, 4, 5, 6, 8, 9 },
    { 2, 4, 5, 7, 8, 9 },
    { 2, 4, 6, 7, 8, 9 },
    { 2, 5, 6, 7, 8, 9 },

    { 3, 4, 5, 6, 7, 8 }, 	{ 3, 4, 5, 6, 7, 9 },
    { 3, 4, 5, 6, 8, 9 },
    { 3, 4, 5, 7, 8, 9 },
    { 3, 4, 6, 7, 8, 9 },
    { 3, 5, 6, 7, 8, 9 },

    { 4, 5, 6, 7, 8, 9 }

};


