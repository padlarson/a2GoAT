#include "AdlarsonPhysics.h"



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

        if( (photons->GetNParticles() == 6) || (photons->GetNParticles() == 7))
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


double_t    AdlarsonPhysics::IM_Ng(UInt_t n )
{
    TLorentzVector g_vec;

    g_vec.SetPxPyPzE(0,0,0,0);
    for ( Int_t i = 0 ; i < n; i++ )
        g_vec += photons->Particle(i);

    return g_vec.M();
}


void AdlarsonPhysics::sixgAnalysis()
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
    // select best combination eta 2pi0 and 3pi0
    // cut in the scatterplot eta 2pi0 3pi0 pdf

    // run kinfit again now enforcing eta mass (and pi0 pi0 as well?)

    // final results: Dalitz plot, m_pi0pi0, eta prime production as fcn energy

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

}

void AdlarsonPhysics::PhysicsResults( const TLorentzVector &g )
{
    // g is the array of TLorentzVectors containing the final and kinfitted event sample in order
    // 0 1 2 3 - p eta pi01 pi02
    // Dalitz plot variables
    Double_t T_eta, T_pi1, T_pi2;
    Double_t X, Y, Q;

    Int_t DPbin;
    Double_t Xgen, Ygen, Xdf, Ydf;


    // Cusp effect variables
    Double_t m_2pi0;
    Double_t m_2pi0gen, m_2pi0df;

    Double_t m_etapr;

    /*
    T_eta = g[1].E() - g[1].M();
    if( g[2].E() > g[3].E() ){
        T_pi1 = g[2].E() - g[2].M();
        T_pi2 = g[3].E() - g[3].M();
    }
    else{
        T_pi1 = g[3].E() - g[3].M();
        T_pi2 = g[2].E() - g[2].M();
    }
    */

    Q = T_eta + T_pi1 + T_pi2;
    X = ( TMath::Sqrt(3.0) / Q ) * ( T_pi1 - T_pi2 );
    Y = ( MASS_ETA + 2 * MASS_PI0 ) / MASS_PI0 * ( T_eta / Q ) - 1.0;







}
