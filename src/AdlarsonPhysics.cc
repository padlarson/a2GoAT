#include "GParticleReconstruction.h"
#include "AdlarsonPhysics.h"
#include "GTrue.h"

#include <APLCON.hpp>



AdlarsonPhysics::AdlarsonPhysics():
    kinfit("etaprime"),
    photons(nPhotons)
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
    // Correlation CB and TAPS
    fi_diff_TAPSCB     = new GH1("fi_diff_TAPSCB", "#Delta#phi TAPS - CB", 200, -200., 200.);
    fi_th_diff_TAPSCB  = new GHistBGSub2("fi_th_diff_TAPSCB", "#Delta#phi v #Delta#theta TAPS - CB", 200, -200., 200., 80, -20., 20.);
    fi_TAPSvsCB        = new GHistBGSub2("fi_TAPSvsCB", "#phi TAPS vs CB", 100, -200.,200.,100, -200.,200.);


// Rec. TAPS - proton analysis
    EvdE_TAPS_all      = new GHistBGSub2("EvdE_TAPS_all", "All E vs dE TAPS", 1200, 0., 600., 200, 0.5, 10.);
    EvdE_TAPS_proton   = new GHistBGSub2("EvdE_TAPS_proton", "Best proton cand E vs dE TAPS", 1200, 0., 600, 200, 0.5, 10.);
    EvTOF              = new GHistBGSub2("EvTOF", "Energy TAPS vs TOF", 200, -20., 20., 400, 0., 800.);

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

    const TLorentzVector target(0,0,0, MASS_PROTON );

    // do I want to initialise my TOF-TAPS timing corrections here?
    // where does Patrik Ott do it?

    kinfit.LinkVariable("Beam",    beam.Link(),       beam.LinkSigma());
    kinfit.LinkVariable("Proton",  proton.Link(),     proton.LinkSigma());

    vector<string> photon_names;
    for(size_t i=0;i<nPhotons;i++) {
        stringstream s_photon;
        s_photon << "Photon" << (i+1);
        photon_names.push_back(s_photon.str());
        kinfit.LinkVariable(s_photon.str(), photons[i].Link(), photons[i].LinkSigma());
    }
    vector<string> all_names = {"Beam", "Proton"};
    all_names.insert(all_names.end(),photon_names.begin(),photon_names.end());


    // Constraint: Incoming 4-vector = Outgoing 4-vector
    auto EnergyMomentumBalance = [] (const vector< vector<double> >& particles) -> vector<double>
    {
        const TLorentzVector target(0,0,0, MASS_PROTON );
        // assume first particle is beam photon
        TLorentzVector diff = target + FitParticle::Make(particles[0], 0.0 );
        // assume second particle outgoing proton
        diff -= FitParticle::Make(particles[1], MASS_PROTON );
        // subtract the rest, assumed to be photons
        for(size_t i=2;i<particles.size();i++) {
            diff -= FitParticle::Make(particles[i], 0.0 );
        }

        return {diff.X(), diff.Y(), diff.Z(), diff.T()};

    };
    kinfit.AddConstraint("EnergyMomentumBalance", all_names, EnergyMomentumBalance);

    // Constraint: Invariant mass of nPhotons equals constant IM,
    // make lambda catch also this with [&] specification
    auto RequireIM = [&] (const vector< vector<double> >& photons) -> double
    {
        TLorentzVector sum(0,0,0,0);
        for(const auto& p : photons) {
            sum += FitParticle::Make(p, 0.0);
        }
        return sum.M() - IM;
    };
    if(includeIMconstraint)
        kinfit.AddConstraint("RequireIM", photon_names, RequireIM);

    // Constraint: Vertex position in z direction: v_z (positive if upstream)
    // if the photon originated from (0,0,v_z) instead of origin,
    // the corrected angle theta' is given by
    // tan(theta') = (R sin(theta))/(R cos(theta) - v_z)
    // R is the CB radius, 10in aka 25.4cm

    auto VertexConstraint = [&] (vector< vector<double> >& photons) -> double
    {
        constexpr double R = 25.4;
        // last element in photons is vz (scalar has dimension 1)
        // see AddConstraint below
        const double v_z = photons.back()[0];
        photons.resize(photons.size()-1); // get rid of last element
        // correct each photon's theta angle,
        // then calculate invariant mass of all photons
        TLorentzVector sum(0,0,0,0);
        for(auto& p : photons) {
            const double theta = p[1]; // second element is theta
            const double theta_p = std::atan2( R*sin(theta), R*cos(theta) - v_z);
            p[1] = theta_p;
            sum += FitParticle::Make(p,0.0);
        }
        return sum.M() - 0.0;
    };

    if(includeVertexFit) {
        kinfit.AddUnmeasuredVariable("v_z"); // default value 0
        kinfit.AddConstraint("VertexConstraint", photon_names + std::vector<string>{"v_z"}, VertexConstraint);
    }

//        static_assert(!(includeIMconstraint && includeVertexFit), "Do not enable Vertex and IM Fit at the same time");

    APLCON::Fit_Settings_t settings = kinfit.GetSettings();
    settings.MaxIterations = 50;
    kinfit.SetSettings(settings);

    cout.precision(3);
    APLCON::PrintFormatting::Width = 11;

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
            for( int j = 0; j < GetTracks()->GetNTracks() ; j++ )
            {
                if( (GetTracks()->GetApparatus(i) == GTreeTrack::APPARATUS_TAPS)  && ( GetTracks()->GetTheta(i) > 18.0 ) )
                {
                    if( (GetTracks()->GetApparatus(j) == GTreeTrack::APPARATUS_CB) && ( GetTracks()->GetTheta(j) < 30.0 ))
                    {
                        fi_diff_TAPSCB->Fill(GetTracks()->GetPhi(i) - GetTracks()->GetPhi(j) );
                        fi_th_diff_TAPSCB->Fill(GetTracks()->GetPhi(i) - GetTracks()->GetPhi(j),  GetTracks()->GetTheta(i) - GetTracks()->GetTheta(j));
                        fi_TAPSvsCB->Fill( GetTracks()->GetPhi(i), GetTracks()->GetPhi(j) );
                    }
                }
            }
        }

        for (Int_t i = 0; i < GetTracks()->GetNTracks() ; i++)
        {            
            if( GetTracks()->GetApparatus(i) == GTreeTrack::APPARATUS_TAPS )
            {
               EvdE_TAPS_all->Fill(GetTracks()->GetClusterEnergy(i),GetTracks()->GetVetoEnergy(i));
               Double_t radnm = 1.45/TMath::Cos( GetTracks()->GetThetaRad(i) );
               EvTOF->Fill( (GetTracks()->GetTime(i))/(radnm), GetTracks()->GetClusterEnergy(i));

               if( cutProtonTAPS->IsInside(GetTracks()->GetClusterEnergy(i), GetTracks()->GetVetoEnergy(i)) )
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


        if( (GetTracks()->GetNTracks() == 7) ) //|| (GetTracks()->GetNTracks() == 8) || (GetTracks()->GetNTracks() == 9) )
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



    for( Int_t i = 0; i < GetTagger()->GetNTagged(); i++ )
    {
        TLorentzVector k = GetTagger()->GetVector(i);
        TLorentzVector l = GetTagger()->GetVectorProtonTarget(i);
        double test = k.E();
        double test2 = l.E();
        beam.SetFromVector(GetTagger()->GetVector(i));

        int n_photons = 0;
        for ( Int_t j = 0; j <= GetTracks()->GetNTracks() ; j++ )
        {
            if( j != ipr)
            {
                IM6g_vec += GetTracks()->GetVector(j);
                IM_6g->Fill( IM6g_vec.M() );
                photons[n_photons].SetFromVector(GetTracks()->GetVector(j));
//                n_photons++;
            }
            else
                proton.SetFromVector(GetTracks()->GetVector(ipr));
        }
//        std::cout << "nr of photons are " << n_photons << std::endl;
    }

/*
    proton.Smear();
    for(auto& photon : photons)
        photon.Smear();
    beam.Smear();
 */

//    IM6gvMMp->Fill(MMp,IM6g_vec.M());



        // Run kinfit:
        // For all (6g or 7g) + proton obtain errors for all observables
        // Run kinfit and select the tagger + 6g + proton with lowest chi2

        // Here run post-kinfit analysis:
        // Find the 15 combinations to form a 2pi0+eta and 3pi0. Determin which reaction it is and remove 3pi0
        // run kinfit with rec. observables now also enforcing eta mass.

        // Obtain physics observables


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

TLorentzVector AdlarsonPhysics::FitParticle::Make(const std::vector<double> &EkThetaPhi, const Double_t m) {
    const double E = EkThetaPhi[0] + m;
    const Double_t p = sqrt( E*E - m*m );
    TVector3 pv(1,0,0);
    pv.SetMagThetaPhi(p, EkThetaPhi[1], EkThetaPhi[2]);
    TLorentzVector l(pv, E);
    return l;
}

void AdlarsonPhysics::FitParticle::Smear() {
    // set the sigmas here,
    // then the fitter knows them as well (because they're linked)

    Ek_Sigma = 0.02*Ek*pow(Ek,-0.36);
    Theta_Sigma = 2.5*TMath::DegToRad();
    if(Theta>20*TMath::DegToRad() && Theta<160*TMath::DegToRad()) {
        Phi_Sigma = Theta_Sigma/sin(Theta);
    }
    else {
        Phi_Sigma = 1*TMath::DegToRad();
    }

}



