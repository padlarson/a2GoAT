#include "GParticleReconstruction.h"
#include "AdlarsonPhysics.h"
#include "GTrue.h"
#include <cmath>

#include <APLCON.hpp>

std::default_random_engine AdlarsonPhysics::FitParticle::generator;

AdlarsonPhysics::AdlarsonPhysics():
    kinfit("etaprime"),
    photons(nPhotons),
    photons_fit(nPhotons),
    photons_rec(nPhotons),
    ClustersTAPSTime("ClustersTAPSTime", "TAPS cluster time" ,200, -50., 50., 450, 0, 450),
    ClustersCBTime("ClustersCBTime", "CB cluster time", 200, -50., 50., 720, 0, 720),
    AllClusters("AllClusters", "nr of all detected clusters",15, 0, 15),
    ClustersinTime("ClustersinTime", "nr of detected clusters in time",15, 0, 15),
    BestTaggedTime("BestTaggedTime", "Tagger time for the events in final result ", 2000, -500, 500),
    etapr_v_BeamE_test1("etapr_v_BeamE_test1", "IM(6#gamma) vs Beam Energy, prompt", 20, 1400, 1600, 100, 600., 1100.),
    etapr_v_BeamE_test2("etapr_v_BeamE_test2", "IM(6#gamma) vs Beam Energy, random", 20, 1400, 1600, 100, 600., 1100.)
{
// TRUE OBSERVABLES
// Beam energy
    True_BeamEnergy     = new GH1("True_BeamEnergy", "True Beam Energy", 250, 1.400, 1.650);
    Tagged_BeamEnergy   = new GH1("Tagged_BeamEnergy", "Tagged Beam Energy", 250, 1400, 1650);
// Phase space final state particles
    ThpvsThetaprCM      = new GHistBGSub2("ThpvsThetaprCM", "#theta_{p} vs #theta_{#eta^{'}}", 360, 0., 180, 360, 0., 180.);
    ThvE_p              = new GHistBGSub2("ThvE_p", "True E_{p} vs #theta_{p}", 100, 0., 0.6, 100, 0., 25.);
    ThvE_eta_g          = new GHistBGSub2("ThvE_eta_g", "E_{#gamma, #eta} vs #theta_{#gamma, #eta}", 100, 0, 1.000, 36, 0, 180);
    ThvE_pi0_g          = new GHistBGSub2("ThvE_pi0_g", "E_{#gamma, #pi^{0}} vs #theta_{#gamma, #pi^{0}}", 60, 0, .600, 36, 0, 180);
// Physics result
    DP_true             = new GHistBGSub2("DP_true", "True Dalitz Plot distribution", 60, -1.5, 1.5, 60, -1.5, 1.5);
    M_pi1pi2_true       = new GH1("M_pi1pi2_true", "True M_{#pi#pi,true}^{2}", 60, 0.05, 0.20);
    M_etapi_true        = new GH1("M_etapi_true", "True M_{#eta#pi,true}^{2}", 160, 0.30, 0.70);

// RECONSTRUCTED OBSERVABLES
    // Correlation CB and TAPS
    fi_diff_TAPSCB     = new GH1("fi_diff_TAPSCB", "#Delta#phi TAPS - CB", 200, -200., 200.);
    fi_th_diff_TAPSCB  = new GHistBGSub2("fi_th_diff_TAPSCB", "#Delta#phi v #Delta#theta TAPS - CB", 200, -200., 200., 80, -20., 20.);
    fi_TAPSvsCB        = new GHistBGSub2("fi_TAPSvsCB", "#phi TAPS vs CB", 100, -200.,200.,100, -200.,200.);


// Rec. TAPS - proton analysis
    EvdE_TAPS_all      = new GHistBGSub2("EvdE_TAPS_all", "All E vs dE TAPS", 1200, 0., 600., 200, 0.5, 10.);
    EvdE_TAPS_proton   = new GHistBGSub2("EvdE_TAPS_proton", "Best proton cand E vs dE TAPS", 1200, 0., 600, 200, 0.5, 10.);
    EvTOF              = new GHistBGSub2("EvTOF", "Energy TAPS vs TOF/ns", 800, -20., 20., 800, 0., 800.);
    EvTOFAll           = new GHistBGSub2("EvTOFAll", "Energy TAPS vs TOF ", 800, -20., 20., 800, 0., 800.);
    EvTOFAllVeto       = new GHistBGSub2("EvTOFAllVeto", "Energy TAPS vs TOF w VETO > 1MeV", 800, -20., 20., 800, 0., 800.);

    Nrprotons           = new GHistBGSub("Nrprotons", "nr of protons", 5,  0, 5);
    ThvEp_rec           = new GHistBGSub2("ThvEp_rec", "Rec E_{p} vs #theta_{p}", 120, 0., 600, 50, 0., 25.);
    MM_p                = new GHistBGSub("MM_p", "Missing Mass calculated for proton", 300, 800., 1100.);

// Rec. Photons
    IM_6g	= new GH1("IM_6g", 	"IM_6g/7g", 240,   200, 1400);
    IM_10g	= new GH1("IM_10g", "IM_10g", 	240,   200, 1400);

    IM6gvMMp    = new GHistBGSub2("IM6gvMMp", "MM(p) vs IM(6#gamma/7#gamma)",240, 200., 1200., 240, 200., 1400.);
    IM10gvMMp   = new GHistBGSub2("IM10gvMMp", "MM(p) vs IM(10#gamma)", 300,800., 1100., 240, 200., 1400.);

// Kinfit related variables

    kfit_chi2           = new GH1("kfit_chi2", "#chi^{2} kinfit", 500, 0, 500.);
    kfit_pdf            = new GH1("kfit_pdf", "#pdf kinfit", 100, 0, 1.);
    kfit_Pulls          = new GHistBGSub2("Pulls", "Pulls", 50, -5., 5., 25, 0, 25);

    IM6g_fit_rec        = new GH1("IM6g_fit_rec", "IM(6#gamma) reconstructed values for events which pass fit", 500, 400., 1400.);
    IM6g_fit            = new GH1("IM6g_fit", "IM(6#gamma) after APLCON fit", 500, 400., 1400.);
    IM6g_fit_3pi        = new GH1("IM6g_fit_3pi", "IM(6#gamma) for 3#pi^{0} candidates", 500, 400., 1400.);
    IM6g_fit_eta2pi     = new GH1("IM6g_fit_eta2pi", "IM(6#gamma) for #eta2#pi^{0} candidates", 500, 400., 1400.);

    PDF_eta2pi_v_3pi    = new GHistBGSub2("PDF_eta2pi_v_3pi", "PDF_eta2pi_v_3pi", 100, 0., 1., 100, 0., 1.);
    best_eta            = new GH1("best_eta", "best #eta cand from comb", 500, 200, 700.);
    best_2pi            = new GH1("best_2pi", "best 2#pi^{0} cand from comb", 500, 0., 500.);
    best_eta_EvTh       = new GHistBGSub2("best_eta_EvTh", "E_{#eta, #gamma} vs #Theta_{#eta, #gamma}", 100, 0., 1000., 50, 0., 200.);


// Physics results

    DP_fit             = new GHistBGSub2("DP_fit", "Rec fitted Dalitz Plot distribution vs run nr", 800, 0, 800, 100, 600., 1100.);
    M_pi1pi2_fit       = new GH1("M_pi1pi2_fit", "Rec M_{#pi#pi,fit}^{2}", 60 , 0.05, 0.20 );
    M_etapi_fit        = new GH1("M_etapi_fit", "Rec M_{#eta#pi,fit}^{2}", 80 , 0.30, 0.70 );
    deltaMpipi_v_Mpipi_fit = new GHistBGSub2("deltaMpipi_v_Mpipi_fit", "fitted - true value M_{#pi#pi,fit}^{2}", 60, 0.05, 0.20, 120, -0.05, 0.05);
    etapr_v_BeamE      = new GHistBGSub2("etapr_v_BeamE", "IM(6#gamma) vs Beam Energy", 20, 1400, 1600, 100, 600., 1100.);


    deltaX_v_DPbin     = new GHistBGSub2("deltaX_v_DPbin", "X_{fit} - X_{true} vd DP bin nr", 600, 0, 600, 40, -1.0, 1.0);
    deltaY_v_DPbin     = new GHistBGSub2("deltaY_v_DPbin", "Y_{fit} - Y_{true} vd DP bin nr", 600, 0, 600, 40, -1.0, 1.0);

    cutFile             = new TFile("/home/adlarson/a2GoAT/configfiles/cuts/TAPSEdEPA.root");
    cutProtonTAPS       = (TCutG*)cutFile->Get("CutProton");

    GHistBGSub::InitCuts(-10, 10, -530, -30);
    GHistBGSub::AddRandCut(30, 530);


    kinfit.LinkVariable("Beam",    beam.Link(),       beam.LinkSigma(),  beam.LinkSettings() );
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
\
    // Constraint: Incoming 4-vector = Outgoing 4-vector
    auto EnergyMomentumBalance = [] (const vector< vector<double> >& particles) -> vector<double>
    {
        const TLorentzVector target(0,0,0, MASS_PROTON);
        // assume first particle is beam photon
        TLorentzVector diff = target + FitParticle::Make(particles[0], 0.0 );
        diff -= FitParticle::Make(particles[1], MASS_PROTON );
        // subtract the rest, assumed to be photons
        for(size_t i=2;i<particles.size();i++) {
            diff -= FitParticle::Make(particles[i], 0.0 );
        }
        return {diff.X(), diff.Y(), diff.Z(), diff.T()};

    };
    kinfit.AddConstraint("EnergyMomentumBalance", all_names, EnergyMomentumBalance);

    /*
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
*/

//        static_assert(!(includeIMconstraint && includeVertexFit), "Do not enable Vertex and IM Fit at the same time");

    APLCON::Fit_Settings_t settings = kinfit.GetSettings();
    settings.MaxIterations = 20;

 //   settings.DebugLevel = 5;
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
//        cout << "ERROR: Input File is not a GoAT file." << endl;
//        return kFALSE;
//        return kTRUE;
    }
    SetAsPhysicsFile();

    ClustersTAPSTime.Reset();
    ClustersCBTime.Reset();
    ClustersinTime.Reset();
    AllClusters.Reset();
    etapr_v_BeamE_test1.Reset();
    etapr_v_BeamE_test2.Reset();

    TraverseValidEvents();

    outputFile->cd();
//    gDirectory->mkdir("MC")->cd();


    ClustersTAPSTime.Write();
    ClustersCBTime.Write();
    ClustersinTime.Write();
    AllClusters.Write();
    etapr_v_BeamE_test1.Write();
    etapr_v_BeamE_test2.Write();

	return kTRUE;
}

void	AdlarsonPhysics::ProcessEvent()
{
//    if(useTrueObservables)
//    etapr_6gTrue.Start(*GetPluto(), *GetGeant());   // (pluto tree, n part in pluto per event)
//    TrueAnalysis_etapr6g();                         // obtains the true observables


//    Kinfit_test();                                  // runs kinematical fit with true observables- for testing purposes

    ClustersInTime.clear();
    ClustersInTime.resize(0);
    UInt_t InTime = 0;
    UInt_t AllTimes = 0;

    for( Int_t i = 0; i < GetTracks()->GetNTracks(); i++ )
    {
        if( GetTracks()->HasTAPS(i) )
        {
            AllTimes++;
            ClustersTAPSTime.Fill( GetTracks()->GetTime(i), GetTracks()->GetCentralCrystal(i) );
            if( ( GetTracks()->GetTime(i) > -8.0) &&  ( GetTracks()->GetTime(i) < 15.0 ) )
            {
                ClustersInTime.push_back( i );
                InTime++;
            }
        }
        if( GetTracks()->HasCB(i) )
        {
            AllTimes++;
            ClustersCBTime.Fill( GetTracks()->GetTime(i), GetTracks()->GetCentralCrystal(i) );
            if( TMath::Abs( GetTracks()->GetTime(i) ) < 15.0 )
            {
                ClustersInTime.push_back( i );
                InTime++;
            }
        }
    }
    AllClusters.Fill(AllTimes);
    ClustersinTime.Fill(InTime);

    Bool_t good_multiplicity = 0;

    if( InTime == 5 )
        good_multiplicity = 1;
    if( (InTime > 6) && (InTime < 9) )
        good_multiplicity = 1;
    if( InTime == 11 )
        good_multiplicity = 1;

    if(!good_multiplicity) return;


    MMp_vec.SetPxPyPzE(0., 0., 0., 0.);
    IM6g_vec.SetPxPyPzE(0., 0., 0., 0.);
    IM10g_vec.SetPxPyPzE(0., 0., 0., 0.);

    for (int i = 0; i < GetTracks()->GetNTracks() ; i++)
    {
        for( int j = 0; j < GetTracks()->GetNTracks() ; j++ )
        {
            if( GetTracks()->HasTAPS(i)  && ( GetTracks()->GetTheta(i) > 18.0 ) )
            {
                if( GetTracks()->HasCB(j) && ( GetTracks()->GetTheta(j) < 30.0 ))
                {
                    fi_diff_TAPSCB->Fill(GetTracks()->GetPhi(i) - GetTracks()->GetPhi(j) );
                    fi_th_diff_TAPSCB->Fill(GetTracks()->GetPhi(i) - GetTracks()->GetPhi(j),  GetTracks()->GetTheta(i) - GetTracks()->GetTheta(j));
                    fi_TAPSvsCB->Fill( GetTracks()->GetPhi(i), GetTracks()->GetPhi(j) );
                }
            }
        }
    }

    Double_t TOF = -100.0;
    for(int t = 0; t < GetTagger()->GetNTagged(); t++)
    {
        for (UInt_t i = 0; i < ClustersInTime.size() ; i++)
        {
            UInt_t j = ClustersInTime[i];
            if( GetTracks()->HasTAPS(j) )
            {
                TOF = GetTagger()->GetTaggedTime(t) - GetTracks()->GetTime(j);
                Double_t radnm = 1.45/TMath::Cos( GetTracks()->GetThetaRad(j) );

                EvTOFAll->Fill( TOF/radnm , GetTracks()->GetClusterEnergy(j));
                if(  GetTracks()->GetVetoEnergy(j) > ( 4.5 -(4.5/180)*GetTracks()->GetClusterEnergy(j)) )
                {
                    if(GetTracks()->GetVetoEnergy(j) > 1.0)
                    {
                        EvTOFAllVeto->Fill(TOF/radnm, GetTracks()->GetClusterEnergy(j));
                    }
                }
                EvTOF->Fill( TOF /(radnm), GetTracks()->GetClusterEnergy(j));
           }
        }
    }

    iprtrack = -1;
    for (UInt_t i = 0; i < ClustersInTime.size() ; i++)
    {
        UInt_t j = ClustersInTime[i];
        if( GetTracks()->HasTAPS(j) )
        {
            EvdE_TAPS_all->Fill(GetTracks()->GetClusterEnergy(j),GetTracks()->GetVetoEnergy(j));
            if( cutProtonTAPS->IsInside(GetTracks()->GetClusterEnergy(j), GetTracks()->GetVetoEnergy(j)))
            {
                EvdE_TAPS_proton->Fill(GetTracks()->GetClusterEnergy(j),GetTracks()->GetVetoEnergy(j));
                nrprotons++;
                iprtrack = j;
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
        Tagged_BeamEnergy->Fill(GetTagger()->GetTaggedEnergy(i),GetTagger()->GetTaggedTime(i));
        MMp_vec.SetPxPyPzE(-proton_vec.Px(), -proton_vec.Py(), GetTagger()->GetTaggedEnergy(i) - proton_vec.Pz(), GetTagger()->GetTaggedEnergy(i) + pdgDB->GetParticle("proton")->Mass()*1000 - proton_vec.E());

        MMp = ( MMp_vec ).M();
        MM_p->Fill(MMp, GetTagger()->GetTaggedTime(i));
    }

    if( (ClustersInTime.size() == 5) )
    {
       // fourgAnalysis(iprtrack);
    }


    if( (ClustersInTime.size() == 7)  || (ClustersInTime.size() == 8)  )// || (GetTracks()->GetNTracks() == 9) )
    {
        sixgAnalysis( iprtrack );
    }

    if( ClustersInTime.size() == 11)
    {
//        tengAnalysis(iprtrack);
    }
}

void	AdlarsonPhysics::ProcessScalerRead()
{
//    hist_eta.ScalerReadCorrection(Double_t(scalers->GetScaler(0))/scalers->GetScaler(1));
//    hist_etap.ScalerReadCorrection(Double_t(scalers->GetScaler(0))/scalers->GetScaler(1));
}

Bool_t	AdlarsonPhysics::Init(const char* configFile)
{
    /*

    ifstream  stream;
    stream.open("/home/adlarson/a2GoAT/configfiles/TAPSvsTAGGERoffsetsAug2012.dat");
    std::vector<double> EPT_TAPS_TOFoff;
    EPT_TAPS_TOFoff.resize(0);
    std::string   line;
    std::getline(stream, line);


    Double newVal1;
    Double newVal2;
    Double newVal3;
    sscanf(GConfigFile::ReadConfig("TOFcorrections", 0, configFile).c_str(), "%lf %lf %lf", &newVal1, &newVal2, &newVal3);

    
    printf("try to enable CBEnergy correction per Run\n");
    
    if (sscanf(line, "%s", tmp) == 1) 
    {
      FILE*	f = fopen(tmp,"r");	
      if(!f)
      {
        Char_t tmp2[140];
        sprintf(tmp2,"data/%s",tmp);
        f = fopen(tmp2,"r");	
      }
      if(!f)
      {
        printf("\nERROR: could not open %s as calibration per Run file! --> exiting\n", tmp);
        printf("   correct the filename after the 'Use-CaLib-CBEnergyPerRun:' keyword in the physics class config file\n");
        printf("   or comment it out. Then no CBEnergy Calibration per Run is done.\n\n");
        exit(1);
      }

      ifstream  stream;
      stream.open(tmp);
      std::vector<double> gainsCB;
      gainsCB.resize(0);

      while( (!feof(f)) && (!file_found) )
      {
          std::string   line;
          std::getline(stream, line);

          //std::string   runnr(line.substr(0,5));
          std::stringstream s_line(line);
          int runnumber;

          */
    
    std::cout << "No Init function specified for this class." << std::endl;
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
    DP_true->Fill( Xtrue, Ytrue );

    m2pi0_metapi0(etapr_true, m_etapi01True, m_etapi02True, m_2pi0True);
    M_pi1pi2_true->Fill(m_2pi0True);
    M_etapi_true->Fill(m_etapi01True);
    M_etapi_true->Fill(m_etapi02True);

    // calculate True Energy vs Theta of final state particles

    ThvE_p->Fill(etapr_6gTrue.GetTrueProtonLV().E() - MASS_PROTON/1000 ,etapr_6gTrue.GetTrueProtonLV().Theta()*TMath::RadToDeg());

    for(UInt_t i = 0; i < etapr_6gTrue.GetNgamma(); i++ )
    {
        if ( i > 3 ) // first four gammas come from 2pi0
            ThvE_eta_g->Fill(etapr_6gTrue.GetTrueGammaLV(i).E(),etapr_6gTrue.GetTrueGammaLV(i).Theta()*TMath::RadToDeg());
        else
            ThvE_pi0_g->Fill(etapr_6gTrue.GetTrueGammaLV(i).E(),etapr_6gTrue.GetTrueGammaLV(i).Theta()*TMath::RadToDeg());
    }
    return;
}

void AdlarsonPhysics::fourgAnalysis(Int_t ipr)
{
    Double_t chi2_min = 1.0e6;
    Double_t prob_min = 1.0e6;
    std::vector<Int_t> set;     // particle conf: tagger, proton, photon 1, 2, 3, 4, 5, 6;
    std::vector<Int_t> set_min; // best particle configuration

    for(int t = 0; t < GetTagger()->GetNTagged(); t++)
    {
        set.resize(0);
        set.push_back(t);
        // Do I need to reset the APLCON linker for every loop?
        beam.SetFromVector( GetTagger()->GetVector(t) );
        beam.Smear(3, true);
        beam.APLCONSettings();

        set.push_back(ipr);
        proton.SetFromVector( GetTracks()->GetVector(ipr) );
        proton.Smear(2, false);

        Int_t n_photons = 0;
        for ( UInt_t j = 0; j < ClustersInTime.size() ; j++ )
        {
            // if 7 neutral clusters run loop here...?
            // Here loop through if more than 7 clusters, but only one proton candidate(?);
            UInt_t k = ClustersInTime[j];
            if( k != ipr ) // the id proton cluster
            {
                set.push_back(k);
                photons_rec[n_photons] = GetTracks()->GetVector(k);
                photons[n_photons].SetFromVector( GetTracks()->GetVector(k));

                if( GetTracks()->HasCB(k) )
                    photons[n_photons].Smear(1, true);
                else // ( GetTracks()->HasTAPS(k) )
                    photons[n_photons].Smear(2, true);
                n_photons++;
            }
        }
        const APLCON::Result_t& result = kinfit.DoFit();
        if(result.Status == APLCON::Result_Status_t::Success)
        {
            if( result.ChiSquare <  chi2_min )
            {
                chi2_min = result.ChiSquare;
                prob_min = result.Probability;
                set_min = set;
            }
        }
    }

        if( (prob_min < 0.01) || ( TMath::Abs(prob_min - 1.0e6) < 1.0e-4) ) return;


}


void AdlarsonPhysics::sixgAnalysis(Int_t ipr)
{
    for(Int_t tag = 0; tag < GetTagger()->GetNTagged(); tag++)
    {
        std::vector<Int_t> set;     // particle conf: tagger, proton, photon 1, 2, 3, 4, 5, 6;
        std::vector<Int_t> set_min; // best particle configuration
        set.resize(0);
        set_min.resize(0);

        for(UInt_t i = 0; i < nPhotons; i++)
        {
            photons_rec[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
            photons_fit[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
        }

        Double_t chi2_min = 1.0e6;
        Double_t prob_min = 1.0e6;


        set.push_back(tag);

        beam.SetFromVector( GetTagger()->GetVector(tag) );
        beam.Smear(3, true);
        beam.APLCONSettings();

        set.push_back(ipr);
        proton.SetFromVector( GetTracks()->GetVector(ipr) );
        proton.Smear(2, false);

        UInt_t n_photons = 0;
        for ( UInt_t jgam = 0; jgam < ClustersInTime.size() ; jgam++ )
        {
            // Here loop through if more than 7 clusters, but only one proton candidate(?);
            Int_t kgam = ClustersInTime[jgam];
            if( kgam != ipr ) // the id proton cluster
            {
                set.push_back(kgam);
                photons_rec[n_photons] = GetTracks()->GetVector(kgam);
                photons[n_photons].SetFromVector( GetTracks()->GetVector(kgam) );

                if( GetTracks()->HasCB(kgam) )
                    photons[n_photons].Smear(1, true);
                else // ( GetTracks()->HasTAPS(kgam) )
                    photons[n_photons].Smear(2, true);
                n_photons++;
            }
        }

        const APLCON::Result_t& result = kinfit.DoFit();
        if(result.Status == APLCON::Result_Status_t::Success)
        {
            if( result.ChiSquare <  chi2_min )
            {
                chi2_min = result.ChiSquare;
                prob_min = result.Probability;
                set_min = set;
            }
        }

        // Here run kinfit with the best combination:

        if( set_min.size() == 0 ) continue;
        if( (prob_min < 0.01) || ( TMath::Abs(prob_min - 1.0e6) < 1.0e-4) ) continue;


        // For the best combination with pdf >= 0.01 obtain the result.

        IM6g_vec.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
        TLorentzVector etap_fit(0.0, 0.0, 0.0, 0.0);

        beam.SetFromVector( GetTagger()->GetVector(set_min[0]) );
        beam.Smear(3, true);
        beam.APLCONSettings();

        proton.SetFromVector( GetTracks()->GetVector(set_min[1]) );
        proton.Smear(2, false);

        Int_t n_photons_min = 0;
        for ( Int_t jgam_min = 0; jgam_min < 6 ; jgam_min++ )
        {
            photons_rec[n_photons_min] = GetTracks()->GetVector( set_min[jgam_min+2] );
            photons[n_photons_min].SetFromVector( GetTracks()->GetVector(set_min[jgam_min+2]));
            IM6g_vec += photons_rec[n_photons_min];

            if( GetTracks()->HasCB(set_min[jgam_min+2]) )
                photons[n_photons_min].Smear(1, true);
            else // ( GetTracks()->HasTAPS(jgam_min) )
                photons[n_photons_min].Smear(2, true);
            n_photons_min++;
        }
        IM6g_fit_rec->Fill(IM6g_vec.M(),GetTagger()->GetTaggedTime(tag));
        IM6gvMMp->Fill( MMp, IM6g_vec.M(),GetTagger()->GetTaggedTime(tag) );

        const APLCON::Result_t& result_min = kinfit.DoFit();
        if(result_min.Status == APLCON::Result_Status_t::Success)
        {
           kfit_chi2->Fill(result_min.ChiSquare, GetTagger()->GetTaggedTime(tag));
           kfit_pdf->Fill(result_min.Probability, GetTagger()->GetTaggedTime(tag));

           Int_t inr = 0;
           for(const auto& it_map : result.Variables) {
     //          const string& varname = it_map.first;
               const APLCON::Result_Variable_t& var = it_map.second;
                kfit_Pulls->Fill(var.Pull, inr);
                inr++;
           }
           for(UInt_t igam_fit = 0; igam_fit < photons.size(); igam_fit++)
           {
               photons_fit[igam_fit] = FitParticle::Make(photons[igam_fit], 0);
               etap_fit += photons_fit[igam_fit];
           }
           IM6g_fit->Fill( etap_fit.M(),GetTagger()->GetTaggedTime(tag) );

           sigma_eta = 17; sigma_pi0 = 9;

           Double_t    chi2min_eta2pi  = 1.0e6;
           Double_t    chi2min_3pi     = 1.0e6;
           std::vector<int> imin_eta2pi;
           std::vector<int> imin_3pi;
           GetBest6gCombination(sigma_eta, sigma_pi0, chi2min_eta2pi, chi2min_3pi, imin_eta2pi, imin_3pi );

           PDF_eta2pi_v_3pi->Fill(TMath::Prob(chi2min_eta2pi,2), TMath::Prob(chi2min_3pi,2), GetTagger()->GetTaggedTime(tag));

           TLorentzVector g[3];
           g[0] = photons_fit[imin_eta2pi[0]] + photons_fit[imin_eta2pi[1]];
           g[1] = photons_fit[imin_eta2pi[2]] + photons_fit[imin_eta2pi[3]];
           g[2] = photons_fit[imin_eta2pi[4]] + photons_fit[imin_eta2pi[5]];

           if( (TMath::Prob(chi2min_3pi,2) > 0.02) && (TMath::Prob(chi2min_eta2pi,2) < 0.02) )
           {
                IM6g_fit_3pi->Fill(etap_fit.M(), GetTagger()->GetTaggedTime(tag));
                best_2pi->Fill( g[1].M(), GetTagger()->GetTaggedTime(tag) );
                best_2pi->Fill( g[2].M(), GetTagger()->GetTaggedTime(tag) );
           }
           if( (TMath::Prob(chi2min_eta2pi,2) > 0.02) && (TMath::Prob(chi2min_3pi,2) < 0.02) )
           {
               IM6g_fit_eta2pi->Fill( etap_fit.M(), GetTagger()->GetTaggedTime(tag) );
               etapr_v_BeamE->Fill((GetTagger()->GetVector(set_min[0])).E(), etap_fit.M(), GetTagger()->GetTaggedTime(tag));
               if(TMath::Abs(GetTagger()->GetTaggedTime(tag)) < 5.0)
                   etapr_v_BeamE_test1.Fill((GetTagger()->GetVector(set_min[0])).E(), etap_fit.M());
               else
                   etapr_v_BeamE_test2.Fill((GetTagger()->GetVector(set_min[0])).E(), etap_fit.M());


               best_eta->Fill( (g[0]).M(), GetTagger()->GetTaggedTime(tag) );
               best_eta_EvTh->Fill( g[1].E(), g[1].Theta()*TMath::RadToDeg(), GetTagger()->GetTaggedTime(tag));
               best_eta_EvTh->Fill( g[2].E(), g[2].Theta()*TMath::RadToDeg(), GetTagger()->GetTaggedTime(tag));

               Double_t m_etapi01_fit = 0;
               Double_t m_etapi02_fit = 0;
               Double_t m_2pi0_fit = 0;
               m2pi0_metapi0( g ,  m_etapi01_fit, m_etapi02_fit, m_2pi0_fit);
               M_pi1pi2_fit->Fill(m_2pi0_fit / 1.0e6, GetTagger()->GetTaggedTime(tag));
               M_etapi_fit->Fill(m_etapi01_fit / 1.0e6, GetTagger()->GetTaggedTime(tag));
               M_etapi_fit->Fill(m_etapi02_fit / 1.0e6, GetTagger()->GetTaggedTime(tag));
               deltaMpipi_v_Mpipi_fit->Fill(m_2pi0_fit / 1.0e6, m_2pi0_fit/1.0e6-m_2pi0True);

               Double_t Xfit = -2.0;
               Double_t Yfit = -2.0;
               Int_t DP_binnr_fit = 0;
               DalitzPlot(g, Xfit, Yfit, DP_binnr_fit);

               DP_fit->Fill(DP_binnr_fit, etap_fit.M(),GetTagger()->GetTaggedTime(tag) );
               deltaX_v_DPbin->Fill(DP_binnr_fit, Xfit - Xtrue);
               deltaY_v_DPbin->Fill(DP_binnr_fit, Yfit - Ytrue);
           }


        }

    } // end loop over tagged time hits

    // run kinfit again now enforcing eta mass (and pi0 pi0 as well?)

    // final results: Dalitz plot, m_pi0pi0, eta prime production as fcn energy
    // DalitzPlot();
    // m2pi0_metapi0;

}

void AdlarsonPhysics::tengAnalysis(Int_t ipr)
{
    photons_rec.resize(0);
    for ( UInt_t i = 0; i < ClustersInTime.size() ; i++ )
    {
        UInt_t j = ClustersInTime[i];
        if(j != iprtrack)
        {
            IM10g_vec += GetTracks()->GetVector(j);
            photons_rec.push_back(GetTracks()->GetVector(j));
        }
    }

    IM_10g->Fill(IM10g_vec.M());
    IM10gvMMp->Fill( MMp , IM10g_vec.M());

 //   for( Int_t i = 0; i < GetTagger()->GetNTagged(); i++ )
 //   {
    // loop through all tagged + proton + 10g combinations --> select permutation with lowest chi2

 //   }
    // chi2 vs taggertime
   /* if( result.Probability < 0.01 ) return;

    Int_t inr = 0;
    for(const auto& it_map : result.Variables) {
//          const string& varname = it_map.first;
        const APLCON::Result_Variable_t& var = it_map.second;
 //        kfit_Pulls10g->Fill(var.Pull, inr);
         inr++;
    }
    for(size_t t=0; t<photons.size(); t++)
    {
        photons_fit[t] = FitParticle::Make(photons[t], 0);
        etap_fit += photons_fit[t];
    }
    IM6g_fit->Fill( etap_fit.M() );

    sigma_eta = 17;

    Double_t    chi2min_eta3pi     = 1.0e6;
    std::vector<int> imin_eta3pi;

    GetBest6gCombination10g(sigma_eta, chi2min_eta3pi, imin_eta3pi );


    // with the best candidate:
    // select best combination eta (6g out of 10g, 220 combinations)

    // run kinfit again now enforcing eta mass (and pi0 pi0 as well?)

    // final results: Dalitz plot, m_pi0pi0, eta prime production as fcn energy
    //
    // DalitzPlot();
    // m2pi0_metapi0;
    //
*/
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

    DP_nr = int ( (X  + X_max ) / bin_width ) + 40*int( (Y  + Y_max ) / bin_width );


}

Double_t    AdlarsonPhysics::IM_Ng( UInt_t n )
{
    TLorentzVector g_vec;

    g_vec.SetPxPyPzE(0,0,0,0);
    for ( UInt_t i = 0 ; i < n; i++ )
        g_vec += GetPhotons()->Particle(i);

    return g_vec.M();
}

void AdlarsonPhysics::m2pi0_metapi0(  TLorentzVector g[3], Double_t &m_etapi01, Double_t &m_etapi02, Double_t &m_2pi0 )
{
    // g is the array of TLorentzVectors containing the final event sample in order
    // 0 1 2 - eta pi01 pi02
    // Dalitz plot variables

    /*    Double_t meta, mpi01, mpi02, metapr; // use for double check


    meta = g[0].M();
    mpi01 = g[1].M();
    mpi02 = g[2].M();
    metapr = ( g[0] + g[1] + g[2] ).M();
    */

    m_etapi01   = ( g[0] + g[1] ).M2();
    m_etapi02   = ( g[0] + g[2] ).M2();
    m_2pi0      = ( g[1] + g[2] ).M2();


}

TLorentzVector AdlarsonPhysics::FitParticle::Make(const std::vector<double> &EkThetaPhi, const Double_t m) {
    const double E = EkThetaPhi[0] + m;
    const Double_t p = sqrt( E*E - m*m );
    TVector3 pv(1,0,0);
    pv.SetMagThetaPhi(p, EkThetaPhi[1], EkThetaPhi[2]);
    TLorentzVector l(pv, E);
    return l;
}

void AdlarsonPhysics::FitParticle::Smear(Int_t apparatus_nr, Bool_t measured) {
    // set the sigmas here,
    // then the fitter knows them as well (because they're linked)

    if( apparatus_nr == 1 ) //CB
    {
        Ek_Sigma = (0.02*(Ek/1.0e3)*pow((Ek/1.0e3),-0.36))*1.0e3*1.3;
        Theta_Sigma = 2.5*TMath::DegToRad()*1.3;
        Phi_Sigma = Theta_Sigma/sin(Theta)*1.3;
    }
    else if( apparatus_nr == 2 ) //TAPS
    {
        if(measured)
            Ek_Sigma = ((0.018 + 0.008*TMath::Sqrt(Ek/1.0e3))*(Ek/1.0e3))*1.0e3*1.3;
        else
            Ek_Sigma = 0;
        Theta_Sigma = 1.0*TMath::DegToRad()*1.3;
        Phi_Sigma = 1.0*TMath::DegToRad()*1.3;
    }
    else  // beam
    {
        Ek_Sigma = 2.0;
        Theta_Sigma = 1e-4;
        Phi_Sigma = 1e-4;
    }

/*
    using gauss_t = std::normal_distribution<double>;
    gauss_t gauss_Ek(0, Ek_Sigma);
    Ek += gauss_Ek(generator);
    gauss_t gauss_Theta(0, Theta_Sigma);
    Theta += gauss_Theta(generator);
    gauss_t gauss_Phi(0, Phi_Sigma);
    Phi += gauss_Phi(generator);
*/

}

void AdlarsonPhysics::FitParticle::APLCONSettings()
{
    Ek_Setting.StepSize = 0;
    Theta_Setting.StepSize = 0;
    Phi_Setting.StepSize = 0;
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

void AdlarsonPhysics::GetBest6gCombination(Double_t& sigma_eta, Double_t& sigma_pi0, Double_t& chi2min_eta2pi, Double_t& chi2min_3pi, std::vector<int>& imin_eta2pi, std::vector<int>& imin_3pi )
{
    imin_eta2pi.resize(6);
    imin_3pi.resize(6);
    Double_t imgg[3];
    Double_t chi2_eta2pi[3];
    Double_t chi2_3pi;



    for( Int_t i = 0; i < 15 ; i++)
    {
       for (Int_t j = 0; j < 3; j++)
       {
            imgg[j] = ( photons_fit[perm6g[i][0 + j*2]] + photons_fit[perm6g[i][1 + j*2]]).M();
            chi2_eta2pi[j] = 0;
       }


        chi2_eta2pi[0] = std::pow((imgg[0] - MASS_ETA )/(sigma_eta), 2.0) + std::pow((imgg[1] - MASS_PI0 )/(sigma_pi0), 2) + std::pow((imgg[2] - MASS_PI0 )/(sigma_pi0), 2.0);
        chi2_eta2pi[1] = std::pow((imgg[1] - MASS_ETA )/(sigma_eta), 2.0) + std::pow((imgg[2] - MASS_PI0 )/(sigma_pi0), 2) + std::pow((imgg[0] - MASS_PI0 )/(sigma_pi0), 2.0);
        chi2_eta2pi[2] = std::pow((imgg[2] - MASS_ETA )/(sigma_eta), 2.0) + std::pow((imgg[0] - MASS_PI0 )/(sigma_pi0), 2) + std::pow((imgg[1] - MASS_PI0 )/(sigma_pi0), 2.0);


        for( Int_t t = 0; t < 3; t++ )
        if( chi2_eta2pi[t] < chi2min_eta2pi )
        {
            chi2min_eta2pi = chi2_eta2pi[t];
            if( t == 0 )
            {
                for(Int_t k = 0; k < 6; k++)
                   imin_eta2pi[k] = perm6g[i][k];
            }
            else if( t == 1 )
            {
                imin_eta2pi[0] = perm6g[i][2];
                imin_eta2pi[1] = perm6g[i][3];
                imin_eta2pi[2] = perm6g[i][0];
                imin_eta2pi[3] = perm6g[i][1];
                imin_eta2pi[4] = perm6g[i][4];
                imin_eta2pi[5] = perm6g[i][5];
            }
            else // t == 2
            {
                imin_eta2pi[0] = perm6g[i][4];
                imin_eta2pi[1] = perm6g[i][5];
                imin_eta2pi[2] = perm6g[i][0];
                imin_eta2pi[3] = perm6g[i][1];
                imin_eta2pi[4] = perm6g[i][2];
                imin_eta2pi[5] = perm6g[i][3];
            }
        }

        chi2_3pi = 0;
        for( Int_t k = 0; k < 3; k++ )
            chi2_3pi += std::pow((imgg[k] - MASS_PI0 )/(sigma_pi0), 2);
        if( chi2_3pi < chi2min_3pi )
        {
            chi2min_3pi = chi2_3pi;
            for(Int_t t = 0; t < 6; t++)
                imin_3pi[t] = perm6g[i][t];
        }
    }
}

/*
void AdlarsonPhysics::GetBest6gCombination10g(Double_t& sigma_eta, Double_t& chi2min_eta3pi, std::vector<int>& imin_eta3pi );
{
    imin_eta2pi.resize(6);
    imin_3pi.resize(6);
    Double_t imgg[3];
    Double_t chi2_eta2pi[3];
    Double_t chi2_3pi;

    for( Int_t i = 0; i < 15 ; i++)
    {
       for (Int_t j = 0; j < 3; j++)
            imgg[j] = ( photons_fit[perm6g[i][0 + j*2]] + photons_fit[perm6g[i][1 + j*2]]).M();


        chi2_eta2pi[0] = std::pow((imgg[0] - MASS_ETA )/(sigma_eta), 2.0) + std::pow((imgg[1] - MASS_PI0 )/(sigma_pi0), 2) + std::pow((imgg[2] - MASS_PI0 )/(sigma_pi0), 2.0);
        chi2_eta2pi[1] = std::pow((imgg[1] - MASS_ETA )/(sigma_eta), 2.0) + std::pow((imgg[2] - MASS_PI0 )/(sigma_pi0), 2) + std::pow((imgg[0] - MASS_PI0 )/(sigma_pi0), 2.0);
        chi2_eta2pi[2] = std::pow((imgg[2] - MASS_ETA )/(sigma_eta), 2.0) + std::pow((imgg[0] - MASS_PI0 )/(sigma_pi0), 2) + std::pow((imgg[1] - MASS_PI0 )/(sigma_pi0), 2.0);


        for( Int_t t = 0; t < 3; t++ )
        if( chi2_eta2pi[t] < chi2min_eta2pi )
        {
            chi2min_eta2pi = chi2_eta2pi[t];
            if( t == 0 )
            {
                for(Int_t k = 0; k < 6; k++)
                   imin_eta2pi[k] = perm6g[i][k];
            }
            else if( t == 1 )
            {
                imin_eta2pi[0] = perm6g[i][2];
                imin_eta2pi[1] = perm6g[i][3];
                imin_eta2pi[2] = perm6g[i][0];
                imin_eta2pi[3] = perm6g[i][1];
                imin_eta2pi[4] = perm6g[i][4];
                imin_eta2pi[5] = perm6g[i][5];
            }
            else // t == 2
            {
                imin_eta2pi[0] = perm6g[i][4];
                imin_eta2pi[1] = perm6g[i][5];
                imin_eta2pi[2] = perm6g[i][0];
                imin_eta2pi[3] = perm6g[i][1];
                imin_eta2pi[4] = perm6g[i][2];
                imin_eta2pi[5] = perm6g[i][3];
            }
        }

        for( Int_t k = 0; k < 3; k++ )
            chi2_3pi += std::pow((imgg[k] - MASS_PI0 )/(sigma_pi0), 2);
        if( chi2_3pi < chi2min_3pi )
        {
            chi2min_3pi = chi2_3pi;
            for(Int_t t = 0; t < 6; t++)
                imin_3pi[t] = perm6g[i][t];
        }
    }
}

*/


void AdlarsonPhysics::Kinfit_test()
{
    TLorentzVector beam_true, proton_true;
    TLorentzVector photons_true[6];

    beam_true.SetPxPyPzE(0., 0.,  etapr_6gTrue.GetTrueBeamEnergy()*1000, etapr_6gTrue.GetTrueBeamEnergy()*1000);
//    beam.SetFromVector( beam_true );
//    beam.Smear(3, true);

    proton_true.SetPxPyPzE(etapr_6gTrue.GetTrueProtonLV().Px()*1000, etapr_6gTrue.GetTrueProtonLV().Py()*1000, etapr_6gTrue.GetTrueProtonLV().Pz()*1000,etapr_6gTrue.GetTrueProtonLV().E()*1000);
    proton.SetFromVector( proton_true );
    proton.Smear(2, true);

    for(UInt_t i = 0; i < etapr_6gTrue.GetNgamma(); i++ )
    {
        (photons_true[i]).SetPxPyPzE( etapr_6gTrue.GetTrueGammaLV(i).Px()*1000, etapr_6gTrue.GetTrueGammaLV(i).Py()*1000,etapr_6gTrue.GetTrueGammaLV(i).Pz()*1000, etapr_6gTrue.GetTrueGammaLV(i).E()*1000 );
        photons[i].SetFromVector( photons_true[i] );
        photons[i].Smear(1, true);
    }

    TLorentzVector etap_rec;
    for(size_t t=0;t<photons.size();t++)
        etap_rec += FitParticle::Make(photons[t], 0);
    IM_6g->Fill(etap_rec.M());

    const APLCON::Result_t& result = kinfit.DoFit();
    if(result.Status == APLCON::Result_Status_t::Success)
    {
//        cout << result << endl;

       kfit_chi2->Fill(result.ChiSquare);
       kfit_pdf->Fill(result.Probability);

       Int_t inr = 0;
       for(const auto& it_map : result.Variables)
       {
           //const string& varname = it_map.first;
           const APLCON::Result_Variable_t& var = it_map.second;
           kfit_Pulls->Fill(var.Pull, inr);
           inr++;

       }

    }
};



