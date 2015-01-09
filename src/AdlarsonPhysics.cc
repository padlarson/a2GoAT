#include "AdlarsonPhysics.h"
#include "GTreeA2Geant.h"



AdlarsonPhysics::AdlarsonPhysics()
{

    IM_6g	= new GH1("IM_6g", 	"IM_6g/7g", 240,   0, 1200);
    IM_10g	= new GH1("IM_10g", "IM_10g", 	240,   0, 1200);


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
 //   if(protons->GetNParticles() > 0)
        // Strategy is :
        //  Identify all TAPS hits and check from TOF if the particle is a
        // nucleon or not. As default all hits are considered to be photons.

        for (int i = 0; photons->GetNParticles() ; i++)
        {
            if( photons->GetApparatus(i) == GTreeRawParticle::APPARATUS_TAPS )
            {
              // loop through all events
            }
        }

        if( (photons->GetNParticles() == 6) || (photons->GetNParticles() == 7) )
        {
            double_t im6g = IM_Ng( photons->GetNParticles()  );
            IM_6g->Fill( im6g );

            sixgAnalysis();
        }

        if(photons->GetNParticles() == 10)
        {
            double_t im10g = IM_Ng( photons->GetNParticles()  );
            IM_10g->Fill( im10g );

            tengAnalysis();
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


void AdlarsonPhysics::sixgAnalysis()
{
    for( Int_t i = 0; i < tagger->GetNTagged(); i++ )
    {
        for ( Int_t j = 0; j < protons->GetNParticles(); j++ )
        {

            // run kinematical fit energy and momentum conservation:

            // cut on pdf:

        }
    }

    // select best combination eta 2pi0 and 3pi0:
    // GetBest6gCombination();

    // cut in the scatterplot eta 2pi0 3pi0 pdf:

    // run kinfit again now enforcing eta mass (and pi0 pi0 as well?)

    // final results: Dalitz plot, m_pi0pi0, eta prime production as fcn energy
    // DalitzPlot();
    // m2pi0_metapi0;

}

void AdlarsonPhysics::tengAnalysis()
{
    for( Int_t i = 0; i < tagger->GetNTagged(); i++ )
    {
        for ( Int_t j = 0; j < protons->GetNParticles(); j++ )
        {

            // run kinematical fit
            // cut on pdf


        }
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

void AdlarsonPhysics::DalitzPlot(const TLorentzVector *g , Double_t &X, Double_t &Y , Int_t &DP_nr)
{

    // g is the array of TLorentzVectors containing the final event sample in order
    // 0 1 2 - eta pi01 pi02
    // Dalitz plot variables

    Double_t T_eta, T_pi1, T_pi2;
    Double_t Q;
    Double_t X_max = 1.5, Y_max = 1.5, bin_width = 0.2;

    T_eta = g[0].E() - g[0].M();
    if( g[1].E() > g[2].E() ){
        T_pi1 = g[1].E() - g[1].M();
        T_pi2 = g[2].E() - g[2].M();
    }
    else{
        T_pi1 = g[2].E() - g[2].M();
        T_pi2 = g[1].E() - g[1].M();
    }

    Q = T_eta + T_pi1 + T_pi2;
    X = ( TMath::Sqrt(3.0) / Q ) * ( T_pi1 - T_pi2 );
    Y = ( MASS_ETA + 2 * MASS_PI0 ) / MASS_PI0 * ( T_eta / Q ) - 1.0;

    DP_nr = ( int( X ) + X_max ) / 0.2 + 40*( int( Y ) + Y_max );
}


Double_t    AdlarsonPhysics::IM_Ng( UInt_t n )
{
    TLorentzVector g_vec;

    g_vec.SetPxPyPzE(0,0,0,0);
    for ( Int_t i = 0 ; i < n; i++ )
        g_vec += photons->Particle(i);

    return g_vec.M();
}

void m2pi0_metapi0( const TLorentzVector *g, Double_t &m_etapi01, Double_t &m_etapi02, Double_t &m_2pi0 )
{
    // g is the array of TLorentzVectors containing the final event sample in order
    // 0 1 2 - eta pi01 pi02
    // Dalitz plot variables

    m_etapi01   = ( g[0] + g[1] ).M();
    m_etapi02   = ( g[0] + g[2] ).M();
    m_2pi0      = ( g[1] + g[2] ).M();

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

