#include "GParticleReconstruction.h"
#include "AdlarsonPhysics.h"
#include "GTrue.h"
#include <cmath>
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TTree.h"
#include <stdio.h>
#include <string.h>

#include <APLCON.hpp>

std::default_random_engine AdlarsonPhysics::FitParticle::generator;

AdlarsonPhysics::AdlarsonPhysics():
//    time_TOF("time_TOF", "TAPS time_TOF time" ,23000, 0, 23000, 1600, -40., 40.),
    time_clusters_TAPS("time_clusters_TAPS", "TAPS cluster time" ,200, -50., 50., 450, 0, 450),
    time_clusters_CB("time_clusters_CB", "CB cluster time", 250, -125., 125., 720, 0, 720),
    time_clusters_CBavg_CBdiff("time_clusters_CBavg_CBdiff","time_clusters_CBavg_CBdiff",400, -20., 20., 720, 0, 720),
    time_clusters_CBavg_TAPSdiff("time_clusters_CBavg_TAPSdiff", "CB avg time - TAPS cluster time" ,400, -50., 50., 450, 0, 450),
    time_clusters_TAPS_TAPSdiff("time_clusters_TAPS_TAPSdiff", "TAPS - TAPS cluster time (for > 1 TAPS clusters)" ,400, -50., 50., 450, 0, 450),
    time_clusters_CBavg_EPT("time_clusters_CBavg_EPT", "CB avg time - Tagger time" ,400, -50., 50., 50, 0, 50),
    time_nr_AllClusters("time_nr_AllClusters", "nr of all detected clusters",15, 0, 15),
    time_nr_ClustersinTime("time_nr_ClustersinTime", "nr of detected clusters in time window and above Cluster Energy cut",15, 0, 15),
    time_nr_FinalClusterSel("time_nr_FinalClusterSel", "nr of clusters fulfilling all cuts",15, 0, 15),
    six_time_TaggedTime("six_time_TaggedTime", "Tagger time for the events with 7 clusters", 500, -250, 250),
    detnr(6),
    CB_region(4),
    photons_rec(nPhotons_six),
    photons_rec_ten(nPhotons_ten),
    photons_fit(nPhotons_six),
    photons_fit_final(nPhotons_six),
    photons_fit_10g(nPhotons_ten),
    kinfit4g("4greactions"),
    kinfit("etaprime"),
    kinfit_final("etaprime6g_final"),
    kinfit3pi("3pi0hyp"),
    kinfiteta2pi("eta2pihyp"),
    kinfit10g("etaprime10g"),
    kinfit10g_eta2pi("etaprime10g_eta2pihyp"),
    kinfit_2g_pi0("2g_pi0_hyp"),
    kinfit_2g_eta("2g_eta_hyp"),
    Photons_four(nPhotons_four),
    Photons_six(nPhotons_six),
    Photons_six_3pi(nPhotons_six),
    Photons_six_eta2pi(nPhotons_six),
    Photons_six_final(nPhotons_six),
    Photons_ten(nPhotons_ten),
    Photons_ten_eta2pi(nPhotons_ten)
{

// TRUE OBSERVABLES
    true_norm                   = new TH1D("true_norm", "normalisation factor between weighted and generated events", 1000, 0,1000.);
// Beam energy
    true_BeamE                  = new GH1("true_BeamE", "True Beam Energy", 100, 1.400, 1.60);
    true_BeamE_weight           = new GH1("true_BeamE_weight", "True Beam Energy weighted", 200, 1.40, 1.60);
    tag_BeamE                   = new GH1("tag_BeamE", "Tagged Beam Energy", 75, 1400, 1600);

    true_imng                   = new TH1D("true_imng", "IM(Ng) true for 3#pi and #eta#pi", 400, 400., 1200.);
// Phase space final state particles
    true_th_v_E_p               = new TH2F("true_th_v_E_p", "True E_{p} vs #theta_{p}", 100, 0., 0.6, 100, 0., 25.);
    true_th_v_E_eta_g           = new TH2F("true_th_v_E_eta_g", "E_{#gamma, #eta} vs #theta_{#gamma, #eta}", 100, 0, 1000, 90, 0, 180);
    true_th_v_E_eta_6g          = new TH2F("true_th_v_E_eta_6g", "E_{#gamma, #eta --> 6#gamma} vs #theta_{#gamma, #eta}", 100, 0, 1000, 90, 0, 180);
    true_th_v_E_pi0_g           = new TH2F("true_th_v_E_pi0_g", "E_{#gamma, #pi^{0}} vs #theta_{#gamma, #pi^{0}}", 60, 0, 600, 90, 0, 180);
//  Kinfit tests
    true_six_dth_vs_th_p        = new TH2F("true_six_dth_vs_th_p", "proton; #theta_{p,rec}; #theta_{rec}-#theta_{true} (^{o})", 200, 0, 25, 80, -5., 5.);
    true_six_dR_vs_R_p          = new TH2F("true_six_dR_vs_R_p", "proton; #R_{p,rec}; #R_{rec}-#R_{true} (^{o})", 40, 0, 80, 500, -20., 20.);
    true_six_dR_vs_det_p        = new TH2F("true_six_dR_vs_det_p", "proton; TAPS det element; #R_{rec}-#R_{true} (^{o})", 440, 0, 440, 500, -20., 20.);
    true_six_fit_dz_v_z         = new TH2F("six_fit_dz_v_z", "#Delta z vs z_{true}; z_{true} (cm); z_{true} - z_{fit} (cm)", 100, -10., 10., 50, -10., 10.);
    fit_six_fit_dz_v_z          = new TH2F("fit_six_fit_dz_v_z", "#Delta z vs z_{fit}; z_{fit} (cm); z_{fit} - z_{true} (cm)", 100, -10., 10., 50, -10., 10.);

    true_six_fit_dz_v_p_th      = new TH2F("true_six_fit_dz_v_p_th", "#Delta z vs p_th_{true}; pr_{#theta,true}; z_{true} - z_{fit} (cm)", 25, 0., 25, 50, -10., 10.);
    true_z_after_fit            = new GH1("true_z_after_fit", "z_{true} after 4C fit; z_true", 100, -10., 10.);
    true_z_after_final_fit      = new GH1("true_z_after_final_fit", "z_{true} after 7C fit; z_true", 100, -10., 10.);

// Physics result
    true_etapr_diff_distr       = new TH1D("true_etapr_diff_distr", "True differential distribution for eta prime prod", 260, 0, 260);
    true_eta_pr_production      = new TH2D("true_eta_pr_production", "eta prime d#sigma as fcn of #theta_{cm} divided in 12 beam ranges", 500, 0. ,25., 800, 0., 0.08);
    true_eta_pr_gg_effcorr      = new TH2D("true_eta_pr_gg_effcorr", "eta prime --> #gamma#gamma eff corr as fcn of #theta_{cm} divided in 12 beam ranges", 500, 0. ,25., 200, 0., 2.);
    true_DP                     = new TH2D("true_DP", "True Dalitz Plot distribution", 600, -1.5, 1.5, 600, -1.5, 1.5);
    true_phy_DP_020             = new TH1D("true_phy_DP_020", "True Dalitz Plot distribution bw 0.20", 200, 0, 200);
    true_phy_DP_015             = new TH1D("true_phy_DP_015", "True Dalitz Plot distribution bw 0.15", 400, 0, 400);
    true_phy_DP_010             = new TH1D("true_phy_DP_010", "True Dalitz Plot distribution bw 0.10", 800, 0, 800);
    true_phy_DP_075             = new TH1D("true_phy_DP_075", "True Dalitz Plot distribution bw 0.075", 2000, 0, 2000);
    true_phy_DP_005             = new TH1D("true_phy_DP_005", "True Dalitz Plot distribution bw 0.05", 2000, 0, 2000);
    true_M_pi1pi2_e2p           = new TH1D("true_M_pi1pi2_e2p", "True M_{#pi#pi,true}^{2}", 200 , 0.0, 200.);
    true_M_etapi_e2p            = new TH1D("true_M_etapi_e2p", "True M_{#eta#pi,true}^{2}", 400, 400.0, 800.);

    true_six_phy_DP_020         = new TH2F("true_six_phy_DP_020", "Dalitz Plot true vs rec, m(#eta2#pi0) 0.20 bin", 200, 0, 200, 200, 0, 200);
    true_six_phy_DP_015         = new TH2F("true_six_phy_DP_015", "Dalitz Plot true vs rec, m(#eta2#pi0) 0.15 bin", 400, 0, 400, 400, 0, 400);
    true_six_phy_DP_010         = new TH2F("true_six_phy_DP_010", "Dalitz Plot true vs rec, m(#eta2#pi0) 0.10 bin", 800, 0, 800, 800, 0, 800);
    true_six_phy_DP_075         = new TH2F("true_six_phy_DP_075", "Dalitz Plot true vs rec, m(#eta2#pi0) 0.075 bin", 2000, 0, 2000, 2000, 0, 2000);
    true_six_phy_DP_005         = new TH2F("true_six_phy_DP_005", "Dalitz Plot true vs rec, m(#eta2#pi0) 0.05 bin", 2000, 0, 2000, 2000, 0, 2000);

    true_six_phy_DP_020_pr      = new TH2F("true_six_phy_DP_020_pr", "Dalitz Plot true vs rec, m(#eta pr) 0.20 bin", 200, 0, 200, 200, 0, 200);
    true_six_phy_DP_015_pr      = new TH2F("true_six_phy_DP_015_pr", "Dalitz Plot true vs rec, m(#eta pr) 0.15 bin", 400, 0, 400, 400, 0, 400);
    true_six_phy_DP_010_pr      = new TH2F("true_six_phy_DP_010_pr", "Dalitz Plot true vs rec, m(#eta pr) 0.10 bin", 800, 0, 800, 800, 0, 800);
    true_six_phy_DP_075_pr      = new TH2F("true_six_phy_DP_075_pr", "Dalitz Plot true vs rec, m(#eta pr) 0.075 bin", 2000, 0, 2000, 2000, 0, 2000);
    true_six_phy_DP_005_pr      = new TH2F("true_six_phy_DP_005_pr", "Dalitz Plot true vs rec, m(#eta pr) 0.05 bin", 2000, 0, 2000, 2000, 0, 2000);

    six_phy_DP_020              = new GHistBGSub2("six_phy_DP_020", "Dalitz Plot vs m(#eta2#pi0) with 0.20 bin", 200, 0, 200, 250, 800., 1050.);
    six_phy_DP_015              = new GHistBGSub2("six_phy_DP_015", "Dalitz Plot vs m(#eta2#pi0) with 0.15 bin", 400, 0, 400, 250, 800., 1050.);
    six_phy_DP_010              = new GHistBGSub2("six_phy_DP_010", "Dalitz Plot vs m(#eta2#pi0) with 0.10 bin", 800, 0, 800, 250, 800., 1050.);
    six_phy_DP_075              = new GHistBGSub2("six_phy_DP_075", "Dalitz Plot vs m(#eta2#pi0) with 0.075 bin", 2000, 0, 2000, 250, 800., 1050.);
    six_phy_DP_005              = new GHistBGSub2("six_phy_DP_005", "Dalitz Plot vs m(#eta2#pi0) with 0.05 bin", 2000, 0, 2000, 250, 800., 1050.);

    six_phy_DP_020_P1           = new GHistBGSub2("six_phy_DP_020_P1", "Dalitz Plot vs m(#eta2#pi0) 0.2 P>0.01", 200, 0, 200, 250, 800., 1050.);
    six_phy_DP_015_P1           = new GHistBGSub2("six_phy_DP_015_P1", "Dalitz Plot vs m(#eta2#pi0) 0.15 P>0.01", 400, 0, 400, 250, 800., 1050.);
    six_phy_DP_010_P1           = new GHistBGSub2("six_phy_DP_010_P1", "Dalitz Plot vs m(#eta2#pi0) 0.10 P>0.01", 800, 0, 800, 250, 800., 1050.);
    six_phy_DP_075_P1           = new GHistBGSub2("six_phy_DP_075_P1", "Dalitz Plot vs m(#eta2#pi0) 0.075 P>0.01", 2000, 0, 2000, 250, 800., 1050.);
    six_phy_DP_005_P1           = new GHistBGSub2("six_phy_DP_005_P1", "Dalitz Plot vs m(#eta2#pi0) 0.05 P>0.01", 2000, 0, 2000, 250, 800., 1050.);

    six_phy_DP_020_pr           = new GH1("six_phy_DP_020_pr", "Dalitz Plot 0.20 bin for m(#eta pr)", 200, 0, 200);
    six_phy_DP_015_pr           = new GH1("six_phy_DP_015_pr", "Dalitz Plot 0.15 bin for m(#eta pr)", 400, 0, 400);
    six_phy_DP_010_pr           = new GH1("six_phy_DP_010_pr", "Dalitz Plot 0.10 bin for m(#eta pr)", 800, 0, 800);
    six_phy_DP_075_pr           = new GH1("six_phy_DP_075_pr", "Dalitz Plot 0.075 bin for m(#eta pr)", 2000, 0, 2000);
    six_phy_DP_005_pr           = new GH1("six_phy_DP_005_pr", "Dalitz Plot 0.05 bin for m(#eta pr)", 2000, 0, 2000);

// In 6g analysis
    true_phy_3pi_IMpipi_v_IMppi         = new TH2F("true_phy_3pi_IMpipi_v_IMppi", "True M_{#pi0#pi0} vs M_{p#pi0} for 3#pi^{0} ", 100, 0., 1., 200, 1., 3.);

    true_six_phy_dMpipi_v_Mpipi         = new TH2F("true_six_phy_dMpipi_v_Mpipi", "fitted - true value M_{#pi#pi,fit}^{2}", 200, 0.0, 200, 200, -100, 100);
    true_six_phy_dMpipi_v_Mpipi_metapr  = new TH2F("true_six_phy_dMpipi_v_Mpipi_metapr", "metapr fitted - true value M_{#pi#pi,fit}^{2}", 200, 0.0, 200, 200, -100, 100);
    true_six_phy_dX_v_DPbin             = new TH2F("true_six_phy_dX_v_DPbin", "X_{fit} - X_{true} v DP bin nr", 800, 0, 800, 800, -2.0, 2.0);
    true_six_phy_dY_v_DPbin             = new TH2F("true_six_phy_dY_v_DPbin", "Y_{fit} - Y_{true} v DP bin nr", 800, 0, 800, 800, -2.0, 2.0);
    true_six_phy_dX_v_X                 = new TH2F("true_six_phy_dX_v_X", "X_{fit} - X_{true} v X", 200, -2., 2., 200, -1.0, 1.0);
    true_six_phy_dY_v_Y                 = new TH2F("true_six_phy_dY_v_Y", "Y_{fit} - Y_{true} v Y", 200, -2., 2., 200, -1.0, 1.0);
    true_six_phy_dX_v_X_metapr          = new TH2F("true_six_phy_dX_v_X_metapr", "X_{fit} - X_{true} v X m_{eta pr}", 200, -2., 2., 200, -2.0, 2.0);
    true_six_phy_dY_v_Y_metapr          = new TH2F("true_six_phy_dY_v_Y_metapr", "Y_{fit} - Y_{true} v Y m_{eta pr}", 200, -2., 2., 200, -2.0, 2.0);

    true_six_phy_dX_v_DPbin_metapr     = new TH2F("true_six_phy_dX_v_DPbin_metapr", "X_{fit} - X_{true} vd DP bin nr with forced eta prime", 800, 0, 800, 200, -2.0, 2.0);
    true_six_phy_dY_v_DPbin_metapr     = new TH2F("true_six_phy_dY_v_DPbin_metapr", "Y_{fit} - Y_{true} vd DP bin nr with forced eta prime", 800, 0, 800, 200, -2.0, 2.0);

    true_ten_phy_dMpipi_v_Mpipi = new TH2F("true_ten_phy_dMpipi_v_Mpipi", " ten #gamma: fitted - true value M_{#pi#pi,fit}^{2}", 200, 0.0, 200, 200, -100, 100);
    true_ten_phy_dX_v_DPbin     = new TH2F("true_ten_phy_dX_v_DPbin", "ten #gamma: X_{fit} - X_{true} vd DP bin nr", 800, 0, 800, 200, -2.0, 2.0);
    true_ten_phy_dY_v_DPbin     = new TH2F("true_ten_phy_dY_v_DPbin", "ten #gamma: Y_{fit} - Y_{true} vd DP bin nr", 800, 0, 800, 200, -2.0, 2.0);


// RECONSTRUCTED OBSERVABLES
    // Correlation CB and TAPS
//    fi_diff_TAPSCB              = new GH1("fi_diff_TAPSCB", "For clusters close to CB-TAPS border ; #Delta#phi TAPS - CB (^{o}}); Events / 5^{o} ", 80, -200., 200.);
//    fi_th_diff_TAPSCB           = new GHistBGSub2("fi_th_diff_TAPSCB", "#Delta#phi v #Delta#theta TAPS - CB", 200, -200., 200., 80, -20., 20.);
//    fi_TAPSvsCB                 = new GHistBGSub2("fi_TAPSvsCB", "#phi TAPS vs CB", 100, -200.,200.,100, -200.,200.);
// Rec. TAPS - proton analysis

    p_E_v_TOF_TAPS_1cl          = new GHistBGSub2("p_E_v_TOF_TAPS_1cl", "Energy TAPS vs TOF with 1 TAPS cl ", 800, -20., 20., 800, 0., 800.);
    p_E_v_TOF_CB_All_proton     = new GHistBGSub2("p_E_v_TOF_CB_All_proton", "Energy TAPS vs TOF (TAPS - CB) for all proton cand", 800, -20., 20., 800, 0., 800.);
    p_E_v_TOF_CB_2PrID          = new GHistBGSub2("p_E_v_TOF_CB_2PrID", "Energy TAPS vs TOF (TAPS - CB) when more than one proton cand ", 800, -20., 20., 800, 0., 800.);
    p_E_v_TOF_CB_best           = new GHistBGSub2("p_E_v_TOF_CB_best", "Energy TAPS vs TOF (TAPS - CB) best out of multiple proton cand ", 800, -20., 20., 800, 0., 800.);
    p_E_v_TOF_TAPS_2PrID        = new GHistBGSub2("p_E_v_TOF_TAPS_2PrID", "Energy TAPS vs TOF (TAPS - TAPS) if more than 1 cluster in TAPS ", 800, -20., 20., 800, 0., 800.);
    p_E_v_TOF_after_kfit        = new GHistBGSub2("p_E_v_TOF_after_kfit", "fitted proton energy v ToF after kinfit", 800, -20., 20., 800, 0., 800.);
    p_Erec_v_TOF_after_kfit     = new GHistBGSub2("p_Erec_v_TOF_after_kfit", "reconstructed proton energy v ToF after kinfit", 800, -20., 20., 800, 0., 800.);
    p_E_v_TOF_after_kfit_2      = new GHistBGSub2("p_Erec_v_TOF_after_kfit_2", "reconstructed proton energy v ToF after kinfit after cut", 800, -20., 20., 800, 0., 800.);

    p_th_v_E                    = new GHistBGSub2("p_th_v_E", "Rec E_{p} vs #theta_{p}", 120, 0., 600, 50, 0., 25.);
    p_MM                        = new GHistBGSub("MM_p", "Missing Mass calculated for proton", 300, 800., 1100.);

    CB_EnergySum                = new GHistBGSub("CB_EnergySum", "Crystal Ball Energy Sum", 200, 0., 2000.);
    CB_EnergySum_2                = new GHistBGSub("CB_EnergySum_2", "Crystal Ball Energy Sum sum energy content crystals CB", 200, 0., 2000.);
//    CB_EnergySum_2pi0           = new GHistBGSub("CB_EnergySum_2pi0", "Crystal Ball Energy Sum for 2pi0", 400, 0., 2000.);
//    CB_EnergySum_etapi0         = new GHistBGSub("CB_EnergySum_etapi0", "Crystal Ball Energy Sum for etapi0", 400, 0., 2000.);
    CB_EnergySum_3pi0           = new GHistBGSub("CB_EnergySum_3pi0", "Crystal Ball Energy Sum for 3pi0", 200, 0., 2000.);
    CB_EnergySum_etapr          = new GHistBGSub("CB_EnergySum_etapr", "Crystal Ball Energy Sum for eta prime", 200, 0., 2000.);
    CB_EnergySum_3pi0_v_CB      = new GHistBGSub2("CB_EnergySum_3pi0_v_CB", "Crystal Ball Energy Sum for 3pi0 vs CB", 200, 0., 2000., 90, 0, 90);
    CB_EnergySum_etapr_v_CB     = new GHistBGSub2("CB_EnergySum_etapr_v_CB", "Crystal Ball Energy Sum for eta prime vs_CB", 200, 0., 2000., 90, 0,90);

    IMgg_v_det_2pi0_CB          =   new GHistBGSub2("IMgg_v_det_2pi0_CB", "IM(gg) 2#pi^{0}, CB", 50, 0, 250, 720, 0, 720);
    IMgg_v_det_etapi0_pi0_CB    =   new GHistBGSub2("IMgg_v_det_etapi0_pi0_CB", "IM(gg) #eta#pi^{0}, CB", 50, 0, 250, 720, 0, 720);
    IMgg_v_det_etapi0_eta_CB    =   new GHistBGSub2("IMgg_v_det_etapi0_eta_CB", "IM(gg) #eta#pi^{0}, CB", 50, 400, 650, 720, 0, 720);
    IMgg_v_det_3pi0_CB          =   new GHistBGSub2("IMgg_v_det_3pi0_CB", "IM(gg) 3#pi^{0}, CB, CB", 50, 0, 250, 720, 0, 720);
    IMgg_v_det_2pi0_TAPS        =   new GHistBGSub2("IMgg_v_det_2pi0_TAPS", "IM(gg) 2#pi^{0}, TAPS", 50, 0, 250, 440, 0, 440);
    IMgg_v_det_etapi0_TAPS      =   new GHistBGSub2("IMgg_v_det_etapi0_TAPS", "IM(gg) #eta#pi^{0}, TAPS", 200, 0, 1000, 440, 0, 440);
    IMgg_v_det_3pi0_TAPS        =   new GHistBGSub2("IMgg_v_det_3pi0_TAPS", "IM(gg) 3#pi^{0}, TAPS", 50, 0, 250, 440, 0, 440);

    IM6g_v_det_etaprfit_CB          =   new GHistBGSub2("IM6g_v_det_etaprfit_CB", "IM(6g) fit fin.ev.s, CB, CB", 80, 800, 1200, 720, 0, 720);
    IM6g_v_det_etaprfit_TAPS        =   new GHistBGSub2("IM6g_v_det_etaprfit_TAPS", "IM(6g) fit fin.ev.s, TAPS", 80, 800, 1200, 440, 0, 440);
    IM6g_v_det_etaprrec_CB          =   new GHistBGSub2("IM6g_v_det_etaprrec_CB", "IM(6g) rec fin.ev.s, CB, CB", 80, 800, 1200, 720, 0, 720);
    IM6g_v_det_etaprrec_TAPS        =   new GHistBGSub2("IM6g_v_det_etaprrec_TAPS", "IM(6g) rec fin.ev.s, TAPS", 80, 800, 1200, 440, 0, 440);

// Rec. Photons

    ten_rec_IM                  = new GH1("ten_rec_IM", " rec IM(10#gamma)", 400,  400, 1400);
    ten_rec_IM_v_MMp            = new GHistBGSub2("ten_rec_IM_v_MMp", "MM(p) vs IM(10#gamma)", 300,800., 1100., 240, 200., 1400.);

//  Related variables 4g analysis

//    four_rec_IM                 = new GH1("four_rec_IM", "rec. IM(4#gamma)", 280, 0, 1400);
//    four_fit_chi2               = new GH1("four_fit_chi2", "#chi^{2} kinfit for 4g final state", 500, 0, 500.);
//    four_fit_pdf                = new GH1("four_fit_pdf", "#pdf kinfit 4g final state", 100, 0, 1.);

//    four_fit_PDF_etapi_v_2pi    = new GHistBGSub2("four_fit_PDF_etapi_v_2pi", "PDF_etapi_v_2pi 4#gamma", 100, 0., 1., 100, 0., 1.);
//    four_fit_IM                 = new GH1("four_fit_IM", "IM_4#gamma fitted", 260, 0, 1400);

//    four_fit_best_2pi_IM_v_E        = new GHistBGSub2("four_fit_best_2pi_IM_v_E", "E_{#pi^{0}} vs M_{#gamma#gamma} for 2#pi^{0}", 100, 0., 1000., 100, 0., 400.);
//    four_fit_best_etapi_IM_v_E      = new GHistBGSub2("four_fit_best_etapi_IM_v_E", "E_{#pi^{0}} vs M_{#gamma#gamma} for #eta#pi^{0}", 100, 0., 1000., 250, 0., 1000.);

//    // 2pi0 inside phase space
//    four_fit_mgg_v_eth             = new GHistBGSub2("four_fit_mgg_v_eth", "m_{#gamma#gamma} vs E, #theta", 12500, 0, 12500, 50, 50., 250.);

//    // etapi0 eta-->gg cand inside phase space
//    four_fit_m_eta_gg_v_eth        = new GHistBGSub2("four_fit_m_eta_gg_v_eth", "m_{#gamma#gamma} vs E, #theta #eta-->#gamma#gamma" , 12500, 0, 12500, 80, 300., 700.);

//    four_fit_best_2pi0_pi_E_v_th      = new GHistBGSub2("four_fit_best_2pi0_pi_E_v_th", "E_{#pi^{0}, #gamma} vs #Theta_{#pi^{0}, #gamma} 2#pi^{0}", 275, 0., 1200., 200, 0., 200.);
//    four_fit_best_etapi0_pi_E_v_th    = new GHistBGSub2("four_fit_best_etapi0_pi_E_v_th", "E_{#pi^{0}, #gamma} vs #Theta_{#pi^{0}, #gamma} #eta#pi^{0}", 275, 0., 1200., 200, 0., 200.);

    // Kinfit related variables 6g

//    test_six_rec_IM             = new GH1("test_six_rec_IM", " rec. IM(6#gamma)", 280, 0., 1400.);
//    test_six_fit_chi2           = new GH1("test_six_fit_chi2", "#chi^{2} kinfit for 6#gamma", 500, 0, 500.);
//    test_six_fit_pdf            = new GH1("test_six_fit_pdf", "pdf kinfit for 6#gamma", 100, 0, 1.);
//    test_six_fit_Pulls          = new GHistBGSub2("test_six_fit_Pulls", "Pulls 6#gamma", 50, -5., 5., 25, 0, 25);

    six_rec_IM                  = new GH1("six_rec_IM", " rec. IM(6#gamma)", 200, 400., 1400.);
    six_rec_IM_v_MMp            = new GHistBGSub2("six_rec_IM_v_MMp", "MM(p) vs IM(6#gamma)",240, 200., 1200., 240, 200., 1400.);
    six_rec_IM_eta2pi           = new GH1("six_rec_IM_eta2pi", " rec. IM(6#gamma) for eta2pi0 cand", 200, 400., 1400.);

    time_clusters_CB_3pi0       = new GHistBGSub2("time_clusters_CB_3pi0", "CB cluster time for 3pi0 events", 200, -50., 50., 720, 0, 720),
    six_rec_EvTh_6g             = new GHistBGSub2("six_rec_EvTh_6g", "E_{#gamma} vs #Theta_{#gamma} rec 6g cand", 100, 0., 1000., 100, 0., 200.);
    six_rec_EvTh_7g             = new GHistBGSub2("six_rec_EvTh_7g", "E_{#gamma} vs #Theta_{#gamma} rec 7g cand", 100, 0., 1000., 100, 0., 200.);

    six_fit_chi2                = new GH1("six_fit_chi2", "#chi^{2} kinfit for 6#gamma", 20, 0, 20.);
    six_fit_pdf                 = new GH1("six_fit_pdf", "pdf kinfit for 6#gamma", 100, 0, 1.);
    six_fit_etaprfinal_pdf      = new GH1("six_fit_etaprfinal_pdf", "pdf kinfit with eta and 2pi mass enforced", 100, 0, 1.);

    NI6g                        = new GH1("NI6g", "Nr NI6g eta2pi0", 100, 0, 100.);
    NIeta2pi0                   = new GH1("NIeta2pi0", "Nr iterations eta2pi0", 100, 0, 100.);
    NItetapr                    = new GH1("NItetapr", "Nr iterations final etapr", 100, 0, 100.);

    NI3pi0vPDF                  = new GHistBGSub2("NI3pi0vPDF", "Nr iterations 3pi0 vs prob", 100, 0, 100. ,100, 0., 1.);
    NIeta2pi0vPDF               = new GHistBGSub2("NIeta2pi0vPDF", "Nr iterations eta2pi0 vs prob", 100, 0, 100.,100, 0., 1.);
    NItetaprvPDF                = new GHistBGSub2("NItetaprvPDF", "Nr iterations final etapr vs prob", 100, 0, 100.,100, 0., 1.);

    proton_fit_e_v_th           = new GHistBGSub2("proton_fit_e_v_th", "proton E vs #theta", 200, 0., 1000., 25, 0, 25);
    proton_fit_e_v_th_final     = new GHistBGSub2("proton_fit_e_v_th_final", "proton E vs #theta final ev sample", 200, 0., 1000., 25, 0, 25);

    six_fit_mgg_v_eth               = new GHistBGSub2("six_fit_mgg_v_eth", "m_{#gamma#gamma} vs E, #theta, 20 MeV, 1 deg", 13500, 0, 13500, 50, 0., 250.);
    six_fit_mgg_v_eth_2             = new GHistBGSub2("six_fit_mgg_v_eth_2", "m_{#gamma#gamma} vs E, #theta 40 MeV 2 deg", 3600, 0, 3600, 50, 0., 250.);
    six_rec_m6g_sig_v_eth           = new GHistBGSub2("six_rec_m6g_sig_v_eth", "m_{6#gamma} for eta prime decay products; E (MeV); #theta", 1500, 0, 1500, 80, 600., 1200.);

    six_fit_mgg_v_CB               = new GHistBGSub2("six_fit_mgg_v_CB", "m_{#gamma#gamma} vs E, det, 20 MeV, 1 det", 36000, 0, 36000, 50, 0., 250.);
    six_fit_mgg_v_CB_2             = new GHistBGSub2("six_fit_mgg_v_CB_2", "m_{#gamma#gamma} vs E, det, 40 MeV, 1 det", 18000, 0, 18000, 50, 0., 250.);
    six_fit_mgg_v_CB_3             = new GHistBGSub2("six_fit_mgg_v_CB_3", "m_{#gamma#gamma} vs E, det, 20 MeV, 2 det", 18000, 0, 18000, 50, 0., 250.);
    six_fit_mgg_v_CB_4             = new GHistBGSub2("six_fit_mgg_v_CB_4", "m_{#gamma#gamma} vs E, det, 40 MeV, 2 det", 9000, 0, 9000, 50, 0., 250.);

    six_fit_mgg_v_TAPS               = new GHistBGSub2("six_fit_mgg_v_TAPS", "m_{#gamma#gamma} vs E, det, 20 MeV, 1 det", 22000, 0, 22000, 50, 0., 250.);
    six_fit_mgg_v_TAPS_2             = new GHistBGSub2("six_fit_mgg_v_TAPS_2", "m_{#gamma#gamma} vs E, det, 40 MeV, 1 det", 11000, 0, 11000, 50, 0., 250.);
    six_fit_mgg_v_TAPS_3             = new GHistBGSub2("six_fit_mgg_v_TAPS_3", "m_{#gamma#gamma} vs E, det, 20 MeV, 2 det", 11000, 0, 11000, 50, 0., 250.);
    six_fit_mgg_v_TAPS_4             = new GHistBGSub2("six_fit_mgg_v_TAPS_4", "m_{#gamma#gamma} vs E, det, 40 MeV, 2 det", 5500, 0, 5500, 50, 0., 250.);

    six_fit_fitted_p_th_v_det       = new GHistBGSub2("six_fit_fitted_p_th_v_det", "fitted proton theta vs theta rec", 440, 0, 440, 120, -3.0, 3.0);
    six_fit_fitted_p_fi_v_det       = new GHistBGSub2("six_fit_fitted_p_fi_v_det", "fitted proton fi vs phi rec", 440, 0, 440, 200, -10.0, 10.0);

    six_fit_true_p_th_v_det       = new TH2F("six_fit_true_p_th_v_det", "true proton theta vs theta rec", 440, 0, 440, 400, -20.0, 20.0);
    six_fit_true_p_fi_v_det       = new TH2F("six_fit_true_p_fi_v_det", "true proton fi vs phi rec", 440, 0, 440, 800, -200.0, 200.0);

    six_fit_IM                  = new GH1("six_fit_IM", "IM(6#gamma) after APLCON fit", 400, 400., 1200.);
    six_fit_IM_ncl              = new GHistBGSub2("six_fit_IM_ncl", "IM(6#gamma) vs ncl", 10, 0., 10., 400, 400., 1200.);
    six_fit_IM_vz               = new GHistBGSub2("six_fit_IM_vz", "IM(6#gamma) vs z_true", 400, 400., 1200., 100, -10.0, 10.0);
    six_fit_dthpr_vz            = new GHistBGSub2("six_fit_dthpr_vz", "dth proton true - fit vs z_true", 100, -5., 5., 100, -10.0, 10.0);

    six_fit_IM_3pi              = new GH1("six_fit_IM_3pi", "IM(6#gamma) for 3#pi^{0} candidates", 400, 400., 1200.);
    six_fit_IM_eta2pi           = new GH1("six_fit_IM_eta2pi", "IM(6#gamma) for #eta2#pi^{0} candidates", 120, 800., 1100.);
    six_fit_IM_eta2pi_v_ncl     = new GHistBGSub2("six_fit_IM_eta2pi_v_ncl", "IM(6#gamma) for #eta2#pi^{0} vs ncl", 10,0.,10., 120, 800., 1100.);

    six_fit_IM_eta2pi0_b        = new GH1("six_fit_IM_eta2pi0_b", "IM(6#gamma) for #eta2#pi^{0}, P_{#eta2#pi} > 0.01", 120, 800., 1100.);
    six_fit_IM_eta2pi0_c        = new GH1("six_fit_IM_eta2pi0_c", "IM(6#gamma) for #eta2#pi^{0}, P_{#eta2#pi} > 0.04 & P_{3#pi} < 0.0075", 120, 800., 1100.);
    six_fit_IM_eta2pi0_d        = new GH1("six_fit_IM_eta2pi0_d", "IM(6#gamma) for #eta2#pi^{0}, P_{#eta2#pi} > 0.04 & P_{#eta pr} > 0.04", 120, 800., 1100.);
    six_fit_IM_eta2pi0_e        = new GH1("six_fit_IM_eta2pi0_e", "IM(6#gamma) for #eta2#pi^{0}, P_{#eta2#pi} > 0.04 & P_{#eta pr} < 0.04", 120, 800., 1100.);
    six_fit_IM_eta2pi0_f        = new GH1("six_fit_IM_eta2pi0_f", "IM(6#gamma) for #eta2#pi^{0}, P_{#eta2#pi} > 0.02 & P_{#eta pr} > 0.02", 120, 800., 1100.);
    six_fit_IM_eta2pi0_g        = new GH1("six_fit_IM_eta2pi0_g", "IM(6#gamma) for #eta2#pi^{0}, P_{#eta2#pi} > 0.02 & P_{#eta pr} < 0.02", 120, 800., 1100.);

    six_fit_which_place_best_3pi_cand = new GH1("six_fit_which_place_best_3pi_cand", "from rough chi2 test which fits kinfit hypoth best 3pi0", 11, -1, 10);
    six_fit_which_place_best_etapr_cand = new GH1("six_fit_which_place_best_etapr_cand", "from rough chi2 test which fits kinfit hypoth best eta prime", 11, -1, 10);

    six_fit_PDF_eta2pi_v_3pi        = new GHistBGSub2("six_fit_PDF_eta2pi_v_3pi", "PDF_eta2pi_v_3pi 6#gamma using kinfit", 100, 0., 1., 100, 0., 1.);
    six_fit_PDF_eta2pi_v_3pi_2      = new GHistBGSub2("six_fit_PDF_eta2pi_v_3pi_2", "PDF_eta2pi_v_3pi, P(#eta pr) > 0.02 )", 100, 0., 1., 100, 0., 1.);
    six_fit_PDF_eta2pi_v_3pi_4      = new GHistBGSub2("six_fit_PDF_eta2pi_v_3pi_4", "PDF_eta2pi_v_3pi, P(#eta pr) > 0.04 )", 100, 0., 1., 100, 0., 1.);

    six_fit_PDF_eta2pi_v_Meta2pi    = new GHistBGSub2("six_fit_PDF_eta2pi_v_Meta2pi", "PDF_eta2pi_v_3pi vs mass", 100, 0., 1., 80, 800., 1200.);

    six_fit_best_eta                 = new GH1("six_fit_best_eta", "best #eta cand from comb", 80, 400, 800.);
    six_fit_best_2pi                 = new GH1("six_fit_best_2pi", "best 2#pi^{0} cand from comb", 50, 0., 250.);
    six_fit_best_eta_rec             = new GH1("six_fit_best_eta_rec", "best #eta cand from comb rec", 80, 400, 800.);
    six_fit_best_2pi_rec             = new GH1("six_fit_best_2pi_rec", "best 2#pi^{0} cand from comb rec", 50, 0., 250.);

    six_fit_EvTh_g                   = new GHistBGSub2("six_fit_EvTh_g", "E_{#gamma} vs #Theta_{#gamma}", 100, 0., 1000., 100, 0., 200.);
    six_fit_EvTh_g_final             = new GHistBGSub2("six_fit_EvTh_g_final ev sample", "E_{#gamma} vs #Theta_{#gamma} final", 100, 0., 1000., 100, 0., 200.);

    six_fit_best_etapr_eta_E_v_th    = new GHistBGSub2("six_fit_best_etapr_eta_E_v_th", "E_{#eta, #gamma} vs #Theta_{#eta, #gamma}", 100, 0., 1000., 50, 0., 200.);
    six_fit_best_etapr_pi_E_v_th     = new GHistBGSub2("six_fit_best_etapr_pi_E_v_th", "E_{#pi^{0}, #gamma} vs #Theta_{#pi^{0}, #gamma}", 100, 0., 1000., 50, 0., 200.);
    six_fit_best_3pi0_pi_E_v_th      = new GHistBGSub2("six_fit_best_3pi0_pi_E_v_th", "E_{#pi^{0}, #gamma} vs #Theta_{#pi^{0}, #gamma} 3#pi^{0}", 100, 0., 1000., 50, 0., 200.);

    six_fit_eta_PDF_v_Metapr         = new GHistBGSub2("six_fit_eta_PDF_v_Metapr", "PDF when eta mass enforced vs IM, 6#gamma cand", 100, 0., 1., 100, 600., 1100.);

    six_phy_etapr_v_BeamE               = new GHistBGSub2("six_phy_etapr_v_BeamE", "IM(6#gamma) vs Beam Energy", 200, 1400, 1600, 300, 750., 1050.);
    six_phy_etapr_eta2pi_v_BeamE        = new GHistBGSub2("six_phy_etapr_eta2pi_v_BeamE", "IM(6#gamma) with enforced eta and 2pi0 mass vs Beam Energy", 200, 1400, 1600, 300, 750., 1050.);
    six_phy_etapr_v_EPT                 = new GHistBGSub2("six_phy_etapr_v_EPT", "IM(6#gamma) v EPT channel", 48, 0, 48, 300, 750., 1050.);
    six_phy_etapr_eta2pi_v_EPT          = new GHistBGSub2("six_phy_etapr_eta2pi_v_EPT", "IM(6#gamma) with enforced eta and 2pi0 mass vs v EPT channel", 48, 0, 48, 300, 750., 1050.);
    six_phy_etapr_prod_diff_distr       = new GHistBGSub2("six_phy_etapr_prod_diff_distr", "IM(6#gamma) w eta and 2pi0 mass vs E#gamma and #theta_{CM}", 260, 0, 260, 200, 800., 1200.);
    six_phy_etapr_prod_diff_distr_metapr= new GH1("six_phy_etapr_prod_diff_distr_metapr", "IM(6#gamma) w eta2pi0, eta pr mass vs E#gamma and #theta_{CM}", 260, 0, 260);

    six_phy_DP_020              = new GHistBGSub2("six_phy_DP_020", "Dalitz Plot vs m(#eta2#pi0) with 0.20 bin", 200, 0, 200, 250, 800., 1050.);
    six_phy_DP_015              = new GHistBGSub2("six_phy_DP_015", "Dalitz Plot vs m(#eta2#pi0) with 0.15 bin", 400, 0, 400, 250, 800., 1050.);
    six_phy_DP_010              = new GHistBGSub2("six_phy_DP_010", "Dalitz Plot vs m(#eta2#pi0) with 0.10 bin", 800, 0, 800, 250, 800., 1050.);
    six_phy_DP_075              = new GHistBGSub2("six_phy_DP_075", "Dalitz Plot vs m(#eta2#pi0) with 0.075 bin", 2000, 0, 2000, 250, 800., 1050.);
    six_phy_DP_005              = new GHistBGSub2("six_phy_DP_005", "Dalitz Plot vs m(#eta2#pi0) with 0.05 bin", 2000, 0, 2000, 250, 800., 1050.);
    six_phy_DP_X                = new GHistBGSub2("six_phy_DP_X", "Dalitz Plot X vs m(#eta2#pi0)", 30, 0, 1.5, 250, 800., 1050.);
    six_phy_DP_Y                = new GHistBGSub2("six_phy_DP_Y", "Dalitz Plot Y vs m(#eta2#pi0)", 60, -1.5, 1.5, 250, 800., 1050.);

    six_phy_DP_020_pr           = new GH1("six_phy_DP_020_pr", "Dalitz Plot 0.20 bin for m(#eta pr)", 200, 0, 200);
    six_phy_DP_015_pr           = new GH1("six_phy_DP_015_pr", "Dalitz Plot 0.15 bin for m(#eta pr)", 400, 0, 400);
    six_phy_DP_010_pr           = new GH1("six_phy_DP_010_pr", "Dalitz Plot 0.10 bin for m(#eta pr)", 800, 0, 800);
    six_phy_DP_075_pr           = new GH1("six_phy_DP_075_pr", "Dalitz Plot 0.075 bin for m(#eta pr)", 2000, 0, 2000);
    six_phy_DP_005_pr           = new GH1("six_phy_DP_005_pr", "Dalitz Plot 0.05 bin for m(#eta pr)", 2000, 0, 2000);
    six_phy_DP_X_pr             = new GH1("six_phy_DP_X_pr", "Dalitz Plot X proj m(#eta pr)", 30, 0., 1.5);
    six_phy_DP_Y_pr             = new GH1("six_phy_DP_Y_pr", "Dalitz Plot Y proj m(#eta pr)", 60, -1.5, 1.5);


    six_phy_M_pi1pi2_v_etapr    = new GHistBGSub2("six_phy_M_pi1pi2_v_etapr", "Fitted M_{#pi#pi,fit}^{2} 6#gamma", 200 , 0.0, 200., 250, 800., 1050. );
    six_phy_M_etapi_v_etapr     = new GHistBGSub2("six_phy_M_etapi_v_etapr", "Fitted M_{#eta#pi,fit}^{2} 6#gamma", 200 , 400.0, 800., 250, 800., 1050. );

    six_phy_M_pi1pi2_v_etapr    = new GHistBGSub2("six_phy_M_pi1pi2_v_etapr", "Fitted M_{#pi#pi,fit}^{2} 6#gamma", 200 , 0.0, 200., 250, 800., 1050. );
    six_phy_M_etapi_v_etapr     = new GHistBGSub2("six_phy_M_etapi_v_etapr", "Fitted M_{#eta#pi,fit}^{2} 6#gamma", 200 , 400.0, 800., 250, 800., 1050. );

    six_phy_M_pi1pi2_v_etapr_fit= new GH1("six_phy_M_pi1pi2_v_etapr_fit", "Fitted M_{#pi#pi,fit}^{2} 6#gamma with metapr fit", 200 , 0.0, 200.);
    six_phy_M_etapi_v_etapr_fit = new GH1("six_phy_M_etapi_v_etapr_fit", "Fitted M_{#eta#pi,fit}^{2} 6#gamma with metapr fit", 400 , 400.0, 800.);

    six_phy_M_pi1pi2_v_etapr2    = new GHistBGSub2("six_phy_M_pi1pi2_v_etapr2", "Fitted M_{#pi#pi,fit} 6#gamma", 400 , 200, 600, 250, 800., 1050. );
    six_phy_M_etapi_v_etapr2     = new GHistBGSub2("six_phy_M_etapi_v_etapr2", "Fitted M_{#eta#pi,fit} 6#gamma", 400 , 600, 1000, 250, 800., 1050. );

    six_phy_M_pi1pi2_v_etapr_fit2= new GH1("six_phy_M_pi1pi2_v_etapr_fit2", "Fitted M_{#pi#pi,fit} 6#gamma with metapr fit", 400 , 200, 600);
    six_phy_M_etapi_v_etapr_fit2 = new GH1("six_phy_M_etapi_v_etapr_fit2", "Fitted M_{#eta#pi,fit} 6#gamma with metapr fit", 400 , 600, 1000);

    // to check the energy of the  eta pi0 system vs its inv mass
    six_fit_best_eta_IM_v_E     = new GHistBGSub2("six_fit_best_eta_IM_v_E", "E_{#eta} vs M_{#eta}", 80, 400., 1200., 100, 300., 800.);
    six_fit_best_2pi_IM_v_E     = new GHistBGSub2("six_fit_best_2pi_IM_v_E", "E_{#pi^{0}} vs M_{#gamma#gamma} for #eta2#pi^{0}", 100, 0., 1000., 100, 0., 400.);

    // to check the energy of the pi0 system vs its inv mass
    six_fit_best_3pi_IM_v_E     = new GHistBGSub2("six_fit_best_3pi_IM_v_E", "E_{#pi^{0}} vs M_{#gamma#gamma} for 3#pi^{0}", 100, 0., 1000., 100, 0., 400.);
    six_phy_3pi_IMpipi_v_IMppi  = new GHistBGSub2("six_phy_3pi_IMpipi_v_IMppi", "M_{#pi0#pi0} vs M_{p#pi0} for 3#pi^{0}", 100, 0., 1., 200, 1., 3.);

// Kinfit related variables 10g

    kfit_chi2_10g       = new GH1("kfit_chi2_10g", "#chi^{2} kinfit 10#gamma", 20, 0, 20.);
    kfit_pdf_10g        = new GH1("kfit_pdf_10g", "#pdf kinfit 10#gamma", 100, 0, 1.);
    kfit_Pulls_10g      = new GHistBGSub2("kfit_Pulls_10g", "Pulls 10#gamma", 50, -5., 5., 40, 0, 40);
    ten_fit_EvTh_g      = new GHistBGSub2("ten_fit_EvTh_g", "E_{#eta, #gamma} vs #Theta_{#eta, #gamma} 10 #gamma", 200, 0., 1000., 100, 0., 200.);
    ten_fit_EvTh_g_min  = new GHistBGSub2("ten_fit_EvTh_g_min", "E_{#eta, #gamma} vs #Theta_{#eta, #gamma} min E 10 #gamma", 200, 0., 1000., 100, 0., 200.);
    IM10g_fit           = new GH1("IM10g_fit", "IM(10#gamma) after APLCON fit", 200, 700., 1100.);
    IM10g_fit_best_cand = new GH1("IM10g_fit_best_cand", "IM(10#gamma) after APLCON fit with eta and 2pi0 mass enforced", 200, 700., 1100.);

    ten_fit_PDF_eta2pi_v_eta6g   = new GHistBGSub2("ten_fit_PDF_eta2pi_v_eta6g", "ten_fit_PDF_eta2pi_v_eta6g 6#gamma", 100, 0., 1., 100, 0., 1.);

    ten_fit_PDF_eta2pi          = new GH1("ten_fit_PDF_eta2pi", "ten_fit_PDF_eta2pi 10#gamma", 100, 0., 1.);
    ten_fit_X_v_pdf_eta2pi      = new GHistBGSub2("ten_fit_X_v_pdf_eta2pi", "ten_fit_true_v_pdf_eta2pi 10#gamma", 100, 0., 1., 100, -5, 5.);
    ten_fit_Y_v_pdf_eta2pi      = new GHistBGSub2("ten_fit_Y_v_pdf_eta2pi", "ten_fit_true_v_pdf_eta2pi 10#gamma", 100, 0., 1., 100, -5., 5.);

    ten_fit_nfit_v_pdf_5pi0_v_eta2pi0 = new GHistBGSub2("ten_fit_nfit_v_pdf_5pi0_v_eta2pi0", "ten #gamma", 10000, 0, 10000, 50, 0, 50);
    ten_fit_dX_v_pdf_5pi0_v_eta2pi0 = new GHistBGSub2("ten_fit_dX_v_pdf_5pi0_v_eta2pi0", "ten #gamma", 10000, 0, 10000, 50, 0, 50);
    ten_fit_dY_v_pdf_5pi0_v_eta2pi0 = new GHistBGSub2("ten_fit_dY_v_pdf_5pi0_v_eta2pi0", "ten #gamma", 10000, 0, 10000, 50, 0, 50);

    ten_phy_etapr_prod_diff_distr       = new GHistBGSub2("ten_phy_etapr_prod_diff_distr", "IM(10#gamma) w eta and 2pi0 mass vs E#gamma and #theta_{CM}", 260, 0, 260, 200, 800., 1200.);
    ten_phy_etapr_prod_diff_distr_metapr= new GHistBGSub2("ten_phy_etapr_prod_diff_distr_metapr", "IM(10#gamma) w eta2pi0, eta pr mass vs E#gamma and #theta_{CM}", 260, 0, 260, 200,800., 1200.);

    // Physics results eta'
    // APLCON kinfit uncertainties

    g_unc                   = new TFile("configfiles/APLCONunc/photon_uncertainties_0z_R_test.root");
    p_unc                   = new TFile("configfiles/APLCONunc/proton_uncertainties_0z_R.root");
    g_unc_vz                = new TFile("configfiles/APLCONunc/photon_uncertainties_vz.root");
    p_unc_vz                = new TFile("configfiles/APLCONunc/proton_uncertainties_vz.root");

//  uncertainties obtained from particle gun with target length 0 and using R
    g_e     = (TH2F*)g_unc->Get("E")->Clone();
    g_th    = (TH2F*)g_unc->Get("theta")->Clone();
    g_phi   = (TH2F*)g_unc->Get("phi")->Clone();
    g_R     = (TH1F*)g_unc->Get("R")->Clone();

    p_TAPS_th = (TH2F*)p_unc->Get("p_theta")->Clone();
    p_TAPS_fi = (TH2F*)p_unc->Get("p_phi")->Clone();
    p_TAPS_R  = (TH1F*)p_unc->Get("p_R")->Clone();

//  uncertainties obtained from particle gun with target length 10 cm
    g_e_vz     = (TH2F*)g_unc_vz->Get("E")->Clone();
    g_th_vz    = (TH2F*)g_unc_vz->Get("theta")->Clone();
    g_phi_vz   = (TH2F*)g_unc_vz->Get("phi")->Clone();

    p_TAPS_th_vz = (TH2F*)p_unc_vz->Get("p_theta")->Clone();
    p_TAPS_fi_vz = (TH2F*)p_unc_vz->Get("p_phi")->Clone();


    etapr_MC_unc               = new TFile("configfiles/corr/etapr_unc.root");
    CB_unc                     = (TH1F*)etapr_MC_unc->Get("plot_etapr_CB_sgm_ratio")->Clone();

    thcorr_CB                  = new TFile("configfiles/corr/CB_th_corr.root");
    dthvth_CB                  = (TProfile*)thcorr_CB->Get("photon_dtheta_v_theta_CB_pfx")->Clone();

    Ecorr_CB                   = new TFile("configfiles/corr/CB_e_corr.root");
    EvdetCB                    = (TH2F*)Ecorr_CB->Get("g_peak_E_CB");

    Ecorr_TAPS                 = new TFile("configfiles/corr/TAPS_e_corr.root");
    EvdetTAPS                  = (TH2F*)Ecorr_TAPS->Get("g_peak_E_TAPS");

    Ecorr_gamma                = new TFile("configfiles/data/MCEXP3pi0diff.root");
    Eth_gamma                  = (TH2F*)Ecorr_gamma->Get("MCtoEXP");

    Ecorr_gammaNLO             = new TFile("configfiles/data/MCEXP3pi0diffNLO.root");
    Eth_gammaNLO               = (TH2F*)Ecorr_gamma->Get("MCtoEXP");

    thcorr_TAPS                = new TFile("configfiles/corr/TAPS_th_corr.root");
    dthvth_TAPS                = (TProfile*)thcorr_TAPS->Get("photon_dtheta_v_theta_TAPS_pfx")->Clone();

//    weight_bkgd                = new TFile("configfiles/corr/weight_3pi0_etapi0_cocktail.root");
//    MCw_bkgd                   = (TH1D*)weight_bkgd->Get("six_fit_IM_3pi_Buffer")->Clone();
    weight_bkgd                = new TFile("configfiles/corr/threepi_bkgd.root");
    MCw_bkgd                   = (TH2F*)weight_bkgd->Get("EXP")->Clone();
    weight_bkgd2                = new TFile("configfiles/corr/threepi_bkgd2.root");
    MCw_bkgd2                   = (TH2F*)weight_bkgd2->Get("EXP")->Clone();
    weight_bkgd3                = new TFile("configfiles/corr/threepi_bkgd3.root");
    MCw_bkgd3                   = (TH2F*)weight_bkgd3->Get("EXP")->Clone();

    Evth_g_sel                 = new TFile("configfiles/cuts/E_v_th_g_cut.root");
    sixg_cand                  = (TCutG*)Evth_g_sel->Get("CutProton")->Clone();

    Evth_7g_sel                = new TFile("configfiles/cuts/E_vs_th_7g_cut.root");
    seveng_cand                = (TCutG*)Evth_7g_sel->Get("CutProton")->Clone();

    p_ToF_kfit_file            = new TFile("configfiles/cuts/p_ToF_v_pfit.root");
    p_ToF_kfit_cut             = (TCutG*)p_ToF_kfit_file->Get("Hadron")->Clone();

    PDF_cut_file              = new TFile("configfiles/cuts/PDF_cut.root");
    PDF_cut                   = (TCutG*)PDF_cut_file->Get("CUTG")->Clone();
//    GHistBGSub::InitCuts(-8., 8., -88., -8.);
//    GHistBGSub::AddRandCut(8., 88.);

    GHistBGSub::InitCuts(-4.5, 3.5, -84.5, -4.5);
    GHistBGSub::AddRandCut(3.5, 83.5);

//  For final states including 6g
    kinfit.LinkVariable("Beam",    beam.Link(),       beam.LinkSigma(),  beam.LinkSettings() );
    kinfit.LinkVariable("Proton",  proton.Link(),     proton.LinkSigma() );

    vector<string> photon_names;
    for(size_t i=0;i<nPhotons_six;i++) {
        stringstream s_photon;
        s_photon << "Photon" << (i+1);
        photon_names.push_back(s_photon.str());
        kinfit.LinkVariable(s_photon.str(), Photons_six[i].Link(), Photons_six[i].LinkSigma());
    }
    vector<string> all_names = {"Beam", "Proton"};
    all_names.insert(all_names.end(),photon_names.begin(),photon_names.end());


//  To test hypotheses that we have a 3pi0 candidate
       kinfit3pi.LinkVariable("Beam_3pi",    beam_3pi.Link(),   beam_3pi.LinkSigma(),  beam_3pi.LinkSettings() );
       kinfit3pi.LinkVariable("Proton_3pi",  proton_3pi.Link(), proton_3pi.LinkSigma());

       vector<string> photon_names_3pi;
       for(size_t i=0;i<nPhotons_six;i++) {
           stringstream s_photon_3pi;
           s_photon_3pi << "Photon_3pi" << (i+1);
           photon_names_3pi.push_back(s_photon_3pi.str());
           kinfit3pi.LinkVariable(s_photon_3pi.str(), Photons_six_3pi[i].Link(), Photons_six_3pi[i].LinkSigma());
       }
       vector<string> all_names_3pi = {"Beam_3pi", "Proton_3pi"};
       all_names_3pi.insert(all_names_3pi.end(),photon_names_3pi.begin(),photon_names_3pi.end());

//  To test hypotheses that we have eta2pi0
      kinfiteta2pi.LinkVariable("Beam_eta2pi",    beam_eta2pi.Link(),   beam_eta2pi.LinkSigma(),  beam_eta2pi.LinkSettings() );
      kinfiteta2pi.LinkVariable("Proton_eta2pi",  proton_eta2pi.Link(), proton_eta2pi.LinkSigma());

      vector<string> photon_names_eta2pi;
      for(size_t i=0;i<nPhotons_six;i++) {
          stringstream s_photon_eta2pi;
          s_photon_eta2pi << "Photon_eta2pi" << (i+1);
          photon_names_eta2pi.push_back(s_photon_eta2pi.str());
          kinfiteta2pi.LinkVariable(s_photon_eta2pi.str(), Photons_six_eta2pi[i].Link(), Photons_six_eta2pi[i].LinkSigma());
      }
      vector<string> all_names_eta2pi = {"Beam_eta2pi", "Proton_eta2pi"};
       all_names_eta2pi.insert(all_names_eta2pi.end(),photon_names_eta2pi.begin(),photon_names_eta2pi.end());


   //  To test hypotheses that we have etaprime hyp fulfilled
         kinfit_final.LinkVariable("Beam_final",    beam_final.Link(),   beam_final.LinkSigma(),  beam_final.LinkSettings() );
         kinfit_final.LinkVariable("Proton_final",  proton_final.Link(), proton_final.LinkSigma());

         vector<string> photon_names_final;
         for(size_t i=0;i<nPhotons_six;i++) {
             stringstream s_photon_final;
             s_photon_final << "Photons_six_final" << (i+1);
             photon_names_final.push_back(s_photon_final.str());
             kinfit_final.LinkVariable(s_photon_final.str(), Photons_six_final[i].Link(), Photons_six_final[i].LinkSigma());
         }
         vector<string> all_names_final = {"Beam_final", "Proton_final"};
          all_names_final.insert(all_names_final.end(),photon_names_final.begin(),photon_names_final.end());

       //  For final states including 4g
    kinfit4g.LinkVariable("Beam4g",    beam4g.Link(),       beam4g.LinkSigma(),  beam4g.LinkSettings() );
    kinfit4g.LinkVariable("Proton4g",  proton4g.Link(),     proton4g.LinkSigma());

    vector<string> photon_names4g;
    for(size_t i=0;i<nPhotons_four;i++) {
        stringstream s_photon;
        s_photon << "Photon4g" << (i+1);
        photon_names4g.push_back(s_photon.str());
        kinfit4g.LinkVariable(s_photon.str(), Photons_four[i].Link(), Photons_four[i].LinkSigma());
    }
    vector<string> all_names4g = {"Beam4g", "Proton4g"};
    all_names4g.insert(all_names4g.end(),photon_names4g.begin(),photon_names4g.end());

//  For final states including 10g

    kinfit10g.LinkVariable("Beam10g",    beam10g.Link(),       beam10g.LinkSigma(),  beam10g.LinkSettings() );
    kinfit10g.LinkVariable("Proton10g",  proton10g.Link(),     proton10g.LinkSigma());

    vector<string> photon_names10g;
    for(size_t i=0;i<nPhotons_ten;i++) {
        stringstream s_photon;
        s_photon << "Photon10g" << (i+1);
        photon_names10g.push_back(s_photon.str());
        kinfit10g.LinkVariable(s_photon.str(), Photons_ten[i].Link(), Photons_ten[i].LinkSigma());
    }
    vector<string> all_names10g = {"Beam10g", "Proton10g"};
    all_names10g.insert(all_names10g.end(),photon_names10g.begin(),photon_names10g.end());

    //  For final states including 10g with eta and 2pi

    kinfit10g_eta2pi.LinkVariable("Beam10g_eta2pi",    beam10g_eta2pi.Link(),       beam10g_eta2pi.LinkSigma(),  beam10g_eta2pi.LinkSettings() );
    kinfit10g_eta2pi.LinkVariable("Proton10g_eta2pi",  proton10g_eta2pi.Link(),     proton10g_eta2pi.LinkSigma());

        vector<string> photon_names10g_eta2pi;
        for(size_t i=0;i<nPhotons_ten;i++) {
            stringstream s_photon;
            s_photon << "Photon10g_eta2pi" << (i+1);
            photon_names10g_eta2pi.push_back(s_photon.str());
            kinfit10g_eta2pi.LinkVariable(s_photon.str(), Photons_ten_eta2pi[i].Link(), Photons_ten_eta2pi[i].LinkSigma());
        }
        vector<string> all_names10g_eta2pi = {"Beam10g_eta2pi", "Proton10g_eta2pi"};
        all_names10g_eta2pi.insert(all_names10g_eta2pi.end(),photon_names10g_eta2pi.begin(),photon_names10g_eta2pi.end());

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
//    kinfit.AddConstraint("EnergyMomentumBalance", all_names, EnergyMomentumBalance);

//    kinfiteta2pi.AddConstraint("EnergyMomentumBalance", all_names_eta2pi, EnergyMomentumBalance);
//    kinfit3pi.AddConstraint("EnergyMomentumBalance", all_names_3pi, EnergyMomentumBalance);

//    kinfit_final.AddConstraint("EnergyMomentumBalance", all_names_final, EnergyMomentumBalance);

    kinfit4g.AddConstraint("EnergyMomentumBalance", all_names4g, EnergyMomentumBalance);
    kinfit10g.AddConstraint("EnergyMomentumBalance", all_names10g, EnergyMomentumBalance);
    kinfit10g_eta2pi.AddConstraint("EnergyMomentumBalance", all_names10g_eta2pi, EnergyMomentumBalance);

    // Constraint: Invariant mass of nPhotons equals constant IM,
    // make lambda catch also this with [&] specification



    auto RequireIM6g = [&] (const vector< vector<double> >& photons) -> double
    {
        TLorentzVector sum(0,0,0,0);
        for(int i = 0; i < 6; i++) {
            sum += FitParticle::Make(photons[i], 0.0);
        }
        return sum.M() - MASS_ETA;
    };

    auto RequireIMpi1 = [&] (const vector< vector<double> >& photons) -> double
    {
        TLorentzVector sum(0,0,0,0);
        for(int i = 0; i < 2; i++) {
            sum += FitParticle::Make(photons[i], 0.0);
        }
        return sum.M() - MASS_PI0;
    };
    auto RequireIMpi2 = [&] (const vector< vector<double> >& photons) -> double
    {
        TLorentzVector sum(0,0,0,0);
        for(int i = 2; i < 4; i++) {
            sum += FitParticle::Make(photons[i], 0.0);
        }
        return sum.M() - MASS_PI0;
    };
    auto RequireIMpi3 = [&] (const vector< vector<double> >& photons) -> double
    {
        TLorentzVector sum(0,0,0,0);
        for(int i = 4; i < 6; i++) {
            sum += FitParticle::Make(photons[i], 0.0);
        }
        return sum.M() - MASS_PI0;
    };

    auto RequireIMpi4 = [&] (const vector< vector<double> >& photons) -> double
    {
        TLorentzVector sum(0,0,0,0);
        for(int i = 6; i < 8; i++) {
            sum += FitParticle::Make(photons[i], 0.0);
        }
        return sum.M() - MASS_PI0;
    };

    auto RequireIMpi5 = [&] (const vector< vector<double> >& photons) -> double
    {
        TLorentzVector sum(0,0,0,0);
        for(int i = 8; i < 10; i++) {
            sum += FitParticle::Make(photons[i], 0.0);
        }
        return sum.M() - MASS_PI0;
    };

    kinfit10g_eta2pi.AddConstraint("RequireIM6g",photon_names10g_eta2pi, RequireIM6g); //added req that eta' candidate has eta mass.
    kinfit10g_eta2pi.AddConstraint("RequireIMpi1",photon_names10g_eta2pi, RequireIMpi1); //added req that first pair has pi0 mass.
    kinfit10g_eta2pi.AddConstraint("RequireIMpi2",photon_names10g_eta2pi, RequireIMpi2); //added req that second pair has pi0 mass.
    kinfit10g_eta2pi.AddConstraint("RequireIMpi3",photon_names10g_eta2pi, RequireIMpi3); //added req that third pair has pi0 mass.
    kinfit10g_eta2pi.AddConstraint("RequireIMpi4",photon_names10g_eta2pi, RequireIMpi4); //added req that fourth pair has pi0 mass.
    kinfit10g_eta2pi.AddConstraint("RequireIMpi5",photon_names10g_eta2pi, RequireIMpi5); //added req that fifth pair has pi0 mass.

    auto FourMomVertexConstraint = [&] (vector< vector<double> >& args) -> vector<double>
    {
       TLorentzVector diff(0.0, 0.0, 0.0, 0.0);
       TLorentzVector beam(0.0, 0.0, args[0][0], args[0][0]);
       const TLorentzVector target(0,0,0, MASS_PROTON);
       diff = beam + target;

       // E, px, py, pz, make a vector of doubles
       // E  = Energy
       // Px = P*sin(theta)*cos(phi)
       // Py = P*sin(theta)*sin(phi)
       // Pz = P*cos(theta)

       double P, E;
       double mass;
       int  idet;
       TLorentzVector LV;
       std::vector<Double_t> obs;       // observables given as E, P, theta/R, phi, in that order
       obs.resize(0);
       const double v_z = args.back()[0];
       args.resize( args.size()-1 ); // get rid of last element
       for(int j = 1; j < args.size(); j++){
           obs.resize(0);
           if( j == 1 ){ // Proton
                E = args[j][0] + MASS_PROTON;
                P = TMath::Sqrt( E*E - MASS_PROTON*MASS_PROTON );
                obs = {E, P , args[j][1], args[j][2]};
                idet = Is_CB[j-1];
                mass = MASS_PROTON;
                LV = GetLVCorrForZ( obs, v_z, idet, mass );
           }
           else{
               idet = Is_CB[j-1];
               obs = {args[j][0], args[j][0] , args[j][1], args[j][2]};
               mass = 0.0;
               LV = GetLVCorrForZ( obs, v_z, idet ,mass );
           }
           diff -= LV;
       }

       return {diff.X(), diff.Y(), diff.Z(), diff.E()};

    };

  auto v_z_settings = APLCON::Variable_Settings_t::Default;
  v_z_settings.Limit.High = 10.;
  v_z_settings.Limit.Low = -10.;

  kinfit.AddMeasuredVariable("v_z", 0., 2.3, v_z_settings); // default value 0
  kinfiteta2pi.AddMeasuredVariable("v_z", 0., 2.3, v_z_settings); // default value 0
  kinfit3pi.AddMeasuredVariable("v_z", 0., 2.3, v_z_settings); // default value 0
  kinfit_final.AddMeasuredVariable("v_z", 0., 2.3, v_z_settings); // default value 0
//  kinfit10g.AddMeasuredVariable("v_z", 0.0, 2.3, v_z_settings); // default value 0
//   kinfit10g.AddUnmeasuredVariable("v_z", 0.0, v_z_settings); // default value 0
//  kinfit10g_eta2pi.AddMeasuredVariable("v_z", 0., 2.3, v_z_settings); // default value 0

  kinfit.AddConstraint("FourMomVertexConstraint", all_names + std::vector<string>{"v_z"}, FourMomVertexConstraint);
  kinfiteta2pi.AddConstraint("FourMomVertexConstraint", all_names_eta2pi + std::vector<string>{"v_z"}, FourMomVertexConstraint);
  kinfit3pi.AddConstraint("FourMomVertexConstraint", all_names_3pi + std::vector<string>{"v_z"}, FourMomVertexConstraint);
  kinfit_final.AddConstraint("FourMomVertexConstraint", all_names_final + std::vector<string>{"v_z"}, FourMomVertexConstraint);

//  kinfit10g.AddConstraint("FourMomVertexConstraint", all_names10g + std::vector<string>{"v_z"}, FourMomVertexConstraint);
//  kinfit10g_eta2pi.AddConstraint("FourMomVertexConstraint", all_names10g_eta2pi+std::vector<string>{"v_z"}, FourMomVertexConstraint);

  auto RequireIM6g_etapr_vx = [&] (const vector< vector<double> >& args) -> double
  {
      TLorentzVector sum(0,0,0,0);
      int  idet;
      TLorentzVector LV;
      std::vector<Double_t> obs;       // observables given as E, P, theta/R, phi, in that order
      obs.resize(0);
      const double v_z = args.back()[0];
      for(int j = 2; j < 8; j++){
          idet = Is_CB[j-1];
          obs = {args[j][0], args[j][0] , args[j][1], args[j][2]};
          sum += GetLVCorrForZ( obs, v_z, idet ,0.0 );
      }

      return sum.M() - MASS_ETAP;
  };

  auto RequireIM_vx = [&] (const vector< vector<double> >& args) -> double
  {
      TLorentzVector sum(0,0,0,0);
      int  idet;
      TLorentzVector LV;
      std::vector<Double_t> obs;       // observables given as E, P, theta/R, phi, in that order
      obs.resize(0);
      const double v_z = args.back()[0];
      for(int j = 2; j < 4; j++){
          idet = Is_CB[j-1];
          obs = {args[j][0], args[j][0] , args[j][1], args[j][2]};
          sum += GetLVCorrForZ( obs, v_z, idet ,0.0 );
      }
      return sum.M() - MASS_ETA;
  };

//  auto RequireIM6g_vx = [&] (const vector< vector<double> >& args) -> double
//  {
//      TLorentzVector sum(0,0,0,0);
//      int  idet;
//      TLorentzVector LV;
//      std::vector<Double_t> obs;       // observables given as E, P, theta/R, phi, in that order
//      obs.resize(0);
//      const double v_z = args.back()[0];
//      for(int j = 2; j < 8; j++){
//          idet = Is_CB[j-1];
//          obs = {args[j][0], args[j][0] , args[j][1], args[j][2]};
//          sum += GetLVCorrForZ( obs, v_z, idet ,0.0 );
//      }
//      return sum.M() - MASS_ETA;
//  };

  auto RequireIMpi1_vx = [&] (const vector< vector<double> >& args) -> double
  {
      TLorentzVector sum(0,0,0,0);
      int  idet;
      TLorentzVector LV;
      std::vector<Double_t> obs;       // observables given as E, P, theta/R, phi, in that order
      obs.resize(0);
      const double v_z = args.back()[0];
      for(int j = 2; j < 4; j++){
          idet = Is_CB[j-1];
          obs = {args[j][0], args[j][0] , args[j][1], args[j][2]};
          sum += GetLVCorrForZ( obs, v_z, idet ,0.0 );
      }
      return sum.M() - MASS_PI0;
  };
  auto RequireIMpi2_vx = [&] (const vector< vector<double> >& args) -> double
  {
      TLorentzVector sum(0,0,0,0);
      int  idet;
      TLorentzVector LV;
      std::vector<Double_t> obs;       // observables given as E, P, theta/R, phi, in that order
      obs.resize(0);
      const double v_z = args.back()[0];
      for(int j = 4; j < 6; j++){
          idet = Is_CB[j-1];
          obs = {args[j][0], args[j][0] , args[j][1], args[j][2]};
          sum += GetLVCorrForZ( obs, v_z, idet ,0.0 );
      }
      return sum.M() - MASS_PI0;
  };
  auto RequireIMpi3_vx = [&] (const vector< vector<double> >& args) -> double
  {
      TLorentzVector sum(0,0,0,0);
      int  idet;
      TLorentzVector LV;
      std::vector<Double_t> obs;       // observables given as E, P, theta/R, phi, in that order
      obs.resize(0);
      const double v_z = args.back()[0];
      for(int j = 6; j < 8; j++){
          idet = Is_CB[j-1];
          obs = {args[j][0], args[j][0] , args[j][1], args[j][2]};
          sum += GetLVCorrForZ( obs, v_z, idet ,0.0 );
      }
      return sum.M() - MASS_PI0;
  };

//  auto RequireIMpi4_vx = [&] (const vector< vector<double> >& args) -> double
//  {
//      TLorentzVector sum(0,0,0,0);
//      int  idet;
//      TLorentzVector LV;
//      std::vector<Double_t> obs;       // observables given as E, P, theta/R, phi, in that order
//      obs.resize(0);
//      const double v_z = args.back()[0];
//      for(int j = 8; j < 10; j++){
//          idet = Is_CB[j-1];
//          obs = {args[j][0], args[j][0] , args[j][1], args[j][2]};
//          sum += GetLVCorrForZ( obs, v_z, idet ,0.0 );
//      }
//      return sum.M() - MASS_PI0;
//  };

//  auto RequireIMpi5_vx = [&] (const vector< vector<double> >& args) -> double
//  {
//      TLorentzVector sum(0,0,0,0);
//      int  idet;
//      TLorentzVector LV;
//      std::vector<Double_t> obs;       // observables given as E, P, theta/R, phi, in that order
//      obs.resize(0);
//      const double v_z = args.back()[0];
//      for(int j = 10; j < 12; j++){
//          idet = Is_CB[j-1];
//          obs = {args[j][0], args[j][0] , args[j][1], args[j][2]};
//          sum += GetLVCorrForZ( obs, v_z, idet ,0.0 );
//      }
//      return sum.M() - MASS_PI0;
//  };

  kinfit3pi.AddConstraint("RequireIMpi1_vx", all_names_3pi + std::vector<string>{"v_z"}, RequireIMpi1_vx); //added req that first pair has pi0 mass.
  kinfit3pi.AddConstraint("RequireIMpi2_vx", all_names_3pi + std::vector<string>{"v_z"}, RequireIMpi2_vx); //added req that second pair has pi0 mass.
  kinfit3pi.AddConstraint("RequireIMpi3_vx", all_names_3pi + std::vector<string>{"v_z"}, RequireIMpi3_vx); //added req that third pair has pi0 mass.

  kinfiteta2pi.AddConstraint("RequireIM_vx", all_names_eta2pi + std::vector<string>{"v_z"}, RequireIM_vx); //added req that eta' candidate has eta mass.
  kinfiteta2pi.AddConstraint("RequireIMpi2_vx", all_names_eta2pi + std::vector<string>{"v_z"}, RequireIMpi2_vx); //added req that second pair has pi0 mass.
  kinfiteta2pi.AddConstraint("RequireIMpi3_vx", all_names_eta2pi + std::vector<string>{"v_z"}, RequireIMpi3_vx); //added req that third pair has pi0 mass.

  kinfit_final.AddConstraint("RequireIM_vx", all_names_final + std::vector<string>{"v_z"}, RequireIM_vx); //added req that eta' candidate has eta mass.
  kinfit_final.AddConstraint("RequireIMpi2_vx", all_names_final + std::vector<string>{"v_z"}, RequireIMpi2_vx); //added req that second pair has pi0 mass.
  kinfit_final.AddConstraint("RequireIMpi3_vx", all_names_final + std::vector<string>{"v_z"}, RequireIMpi3_vx); //added req that third pair has pi0 mass.
  kinfit_final.AddConstraint("RequireIM6g_etapr_vx", all_names_final + std::vector<string>{"v_z"}, RequireIM6g_etapr_vx); //added req that third pair has pi0 mass.

//  kinfit10g_eta2pi.AddConstraint("RequireIM6g_vx", all_names10g_eta2pi+std::vector<string>{"v_z"}, RequireIM6g_vx);
//  kinfit10g_eta2pi.AddConstraint("RequireIMpi4_vx", all_names10g_eta2pi+std::vector<string>{"v_z"}, RequireIMpi4_vx);
//  kinfit10g_eta2pi.AddConstraint("RequireIMpi5_vx", all_names10g_eta2pi+std::vector<string>{"v_z"}, RequireIMpi5_vx);

//    kinfit.AddConstraint("VertexConstraint", all_names + std::vector<string>{"v_z"}, VertexConstraint);
//    kinfit.AddUnmeasuredVariable("v_z"); // default value 0

    APLCON::Fit_Settings_t settings = kinfit.GetSettings();
    settings.MaxIterations = 30;


    // 25, 35, 35, 35
    // 30 40 40 40

    APLCON::Fit_Settings_t settings_eta2pi = kinfiteta2pi.GetSettings();
    settings_eta2pi.MaxIterations = 30;
    APLCON::Fit_Settings_t settings_3pi = kinfit3pi.GetSettings();
    settings_3pi.MaxIterations = 30;
    APLCON::Fit_Settings_t settings_metapr = kinfit_final.GetSettings();
    settings_metapr.MaxIterations = 30;

//    settings.DebugLevel = 5;
//    settings_eta2pi.DebugLevel = 5;
    kinfit.SetSettings(settings);
    kinfiteta2pi.SetSettings(settings_eta2pi);
    kinfit3pi.SetSettings(settings_3pi);
    kinfit_final.SetSettings(settings_metapr);

    APLCON::Fit_Settings_t settings10g = kinfit10g.GetSettings();
    settings10g.MaxIterations = 30;

    cout.precision(5);
    APLCON::PrintFormatting::Width = 11;

}
AdlarsonPhysics::~AdlarsonPhysics()
{}

Bool_t	AdlarsonPhysics::Start()
{
    pRandoms = new TRandom3;
    pRandoms->SetSeed(0); //'Randomize Timer'
    if(!IsGoATFile())
    {
        cout << "ERROR: Input File is not a GoAT file." << endl;
        return kFALSE;
//        return kTRUE;
    }

    const char* str = outputFile->GetName();
    char * my_str = strdup(str);
    char * pch;
    char * pch2;
    TString prefix = "configfiles/output/";
    TString suffix = "result_";
    pch = strtok (my_str,"/");
    while (pch != NULL)
    {
      pch2 = pch;
      pch = strtok (NULL, "/");
    }
    
    tree_file_name = prefix + suffix + pch2;

//    f_tree = new TFile(tree_file_name,"recreate");
//    std::cout << tree_file_name << std::endl;
//    if(f_tree)
//        std::cout << "Created output tree in linked output folder"<< std::endl;
//    tree = new TTree("tree","result");

//    branch_weight        = new double[MAX_LENGTH];
//    branch_nclusters     = new int[MAX_LENGTH];
//    branch_taggerenergy  = new double[MAX_LENGTH];
//    branch_fitted_pr_th  = new double[MAX_LENGTH];
//    branch_fitted_pr_e   = new double[MAX_LENGTH];
//    branch_fitted_zvx    = new double[MAX_LENGTH];
//    branch_pdf_eta2pi    = new double[MAX_LENGTH];
//    branch_pdf_3pi       = new double[MAX_LENGTH];
//    branch_pdf_etapr     = new double[MAX_LENGTH];
//    branch_costh_epr_cm  = new double[MAX_LENGTH];
//    branch_X             = new double[MAX_LENGTH];
//    branch_Y             = new double[MAX_LENGTH];
//    branch_DP_005        = new int[MAX_LENGTH];
//    branch_DP_010        = new int[MAX_LENGTH];
//    branch_DP_015        = new int[MAX_LENGTH];
//    branch_M_eta2pi      = new double[MAX_LENGTH];
//    branch_M_etapi1      = new double[MAX_LENGTH];
//    branch_M_etapi2      = new double[MAX_LENGTH];
//    branch_M_pipi        = new double[MAX_LENGTH];
//    branch_costh_epr_cm_pr = new double[MAX_LENGTH];
//    branch_X_pr          = new double[MAX_LENGTH];
//    branch_Y_pr          = new double[MAX_LENGTH];
//    branch_DP_005_pr     = new int[MAX_LENGTH];
//    branch_DP_010_pr     = new int[MAX_LENGTH];
//    branch_DP_015_pr     = new int[MAX_LENGTH];
//    branch_M_etapi_pr1   = new double[MAX_LENGTH];
//    branch_M_etapi_pr2   = new double[MAX_LENGTH];
//    branch_M_pipi_pr     = new double[MAX_LENGTH];

//    tree->Branch("branch_length",&branch_length);
//    tree->Branch("weight", branch_weight, "branch_weight[branch_length]/D" );
//    tree->Branch("nclusters", branch_nclusters, "branch_nclusters[branch_length]/I");
//    tree->Branch("beam_energy", branch_taggerenergy, "branch_taggerenergy[branch_length]/D" );
//    tree->Branch("proton_th_fit", branch_fitted_pr_th, "branch_fitted_pr_th[branch_length]/D" );
//    tree->Branch("branch_fitted_pr_e", branch_fitted_pr_e, "branch_fitted_pr_e[branch_length]/D" );
//    tree->Branch("z_vx_fit", branch_fitted_zvx, "branch_fitted_zvx[branch_length]/D" );
//    tree->Branch("pdf_eta2pi", branch_pdf_eta2pi, "branch_pdf_eta2pi[branch_length]/D" );
//    tree->Branch("pdf_3pi", branch_pdf_3pi, "branch_pdf_3pi[branch_length]/D" );
//    tree->Branch("pdf_etapr", branch_pdf_etapr, "branch_pdf_etapr[branch_length]/D" );
//    tree->Branch("costh_epr_cm", branch_costh_epr_cm, "branch_costh_epr_cm[branch_length]/D" );
//    tree->Branch("X", branch_X, "branch_X[branch_length]/D" );
//    tree->Branch("Y", branch_Y, "branch_Y[branch_length]/D" );
//    tree->Branch("DP_005bin", branch_DP_005, "branch_DP_005[branch_length]/I" );
//    tree->Branch("DP_010bin", branch_DP_010, "branch_DP_010[branch_length]/I" );
//    tree->Branch("DP_015bin", branch_DP_015, "branch_DP_015[branch_length]/I" );
//    tree->Branch("Meta2pi", branch_M_eta2pi, "branch_M_eta2pi[branch_length]/D" );
//    tree->Branch("Metapi1", branch_M_etapi1, "branch_M_etapi1[branch_length]/D" );
//    tree->Branch("Metapi2", branch_M_etapi2, "branch_M_etapi2[branch_length]/D" );
//    tree->Branch("Mpipi", branch_M_pipi, "branch_M_pipi[branch_length]/D" );
//    tree->Branch("costh_epr_cm_pr", branch_costh_epr_cm_pr, "branch_costh_epr_cm_pr[branch_length]/D" );
//    tree->Branch("X_prime", branch_X_pr, "branch_X_pr[branch_length]/D" );
//    tree->Branch("Y_prime", branch_Y_pr, "branch_Y_pr[branch_length]/D" );
//    tree->Branch("DP_005bin_prime", branch_DP_005_pr, "branch_DP_005_pr[branch_length]/I" );
//    tree->Branch("DP_010bin_prime", branch_DP_010_pr, "branch_DP_010_pr[branch_length]/I" );
//    tree->Branch("DP_015bin_prime", branch_DP_015_pr, "branch_DP_015_pr[branch_length]/I" );
//    tree->Branch("Metapi_pr1", branch_M_etapi_pr1, "branch_M_etapi_pr1[branch_length]/D" );
//    tree->Branch("Metapi_pr2", branch_M_etapi_pr2, "branch_M_etapi_pr2[branch_length]/D" );
//    tree->Branch("Mpipi_pr", branch_M_pipi_pr, "branch_M_pipi_pr[branch_length]/D" );

//    TFile*    f_tree2;
//    TTree*    tree2;

    f_tree2 = new TFile(tree_file_name,"recreate");
    std::cout << tree_file_name << std::endl;
    if(f_tree2)
        std::cout << "Created output tree in linked output folder"<< std::endl;
    tree2 = new TTree("tree2","result");

    tree2->Branch("fWeight",&fWeight);
    tree2->Branch("fNclusters",&fNclusters);
    tree2->Branch("fTaggerenergy",&fTaggerenergy);
    tree2->Branch("fProton_th_fit",&fProton_th_fit);
    tree2->Branch("fProton_E_fit",&fProton_E_fit);
    tree2->Branch("fZ_vx_fit",&fZ_vx_fit);
    tree2->Branch("fX",&fX);
    tree2->Branch("fY",&fY);
    tree2->Branch("fPeta2pi",&fPeta2pi);
    tree2->Branch("fP3pi",&fP3pi);
    tree2->Branch("fPetapr",&fPetapr);
    tree2->Branch("fCosth_epr_cm",&fCosth_epr_cm);
    tree2->Branch("fDP005",&fDP005);
    tree2->Branch("fDP075",&fDP075);
    tree2->Branch("fDP010",&fDP010);
    tree2->Branch("fDP015",&fDP015);
    tree2->Branch("fMeta2pi",&fMeta2pi);
    tree2->Branch("fMpipi",&fMpipi);
    tree2->Branch("fMetapi1",&fMetapi1);
    tree2->Branch("fMetapi2",&fMetapi2);
    tree2->Branch("fCosth_epr_cm_pr",&fCosth_epr_cm_pr);
    tree2->Branch("fXpr",&fXpr);
    tree2->Branch("fYpr",&fYpr);
    tree2->Branch("fDP005pr",&fDP005pr);
    tree2->Branch("fDP075pr",&fDP075pr);
    tree2->Branch("fDP010pr",&fDP010pr);
    tree2->Branch("fDP015pr",&fDP015pr);
    tree2->Branch("fMpipipr",&fMpipipr);
    tree2->Branch("fMetapi1pr",&fMetapi1pr);
    tree2->Branch("fMetapi2pr",&fMetapi2pr);

    if (!GetScalers()->IsOpenForInput())
        MC = true;
    else
        MC = false;

    SetAsPhysicsFile();

    time_TOF.Reset();
    time_clusters_TAPS.Reset();
    time_clusters_CB.Reset();
    time_clusters_CBavg_CBdiff.Reset();
    time_clusters_CBavg_TAPSdiff.Reset();
    time_clusters_TAPS_TAPSdiff.Reset();
    time_clusters_CBavg_EPT.Reset();
    time_nr_AllClusters.Reset();
    time_nr_ClustersinTime.Reset();
    time_nr_FinalClusterSel.Reset();
    six_time_TaggedTime.Reset();

    TraverseValidEvents();

    outputFile->cd();

    time_TOF.Write();
    time_clusters_TAPS.Write();
    time_clusters_CB.Write();
    time_clusters_CBavg_CBdiff.Write();
    time_clusters_CBavg_TAPSdiff.Write();
    time_clusters_TAPS_TAPSdiff.Write();
    time_clusters_CBavg_EPT.Write();
    time_nr_AllClusters.Write();
    time_nr_ClustersinTime.Write();
    time_nr_FinalClusterSel.Write();
    six_time_TaggedTime.Write();

    f_tree2->Write();
    f_tree2->Close();
    
//    close and write physics tree
    
	return kTRUE;
}

void	AdlarsonPhysics::ProcessEvent()
{
    if(MC){

       MC_weight = true;
       etapr_6gTrue.Start(*GetPluto(), *GetGeant());   // (pluto tree, n part in pluto per event)
       TrueAnalysis_etapr6g("Sergey");                    // obtains the true observables
       MCw = etapr_6gTrue.GetWeight();
       if(MCJuly14)
         MCw *= 595881./500000.;
       else if(MCOctDec14)
         MCw *= 1008880./1000000.;

//*****       for 3pi0 and etapi0 MC *****
//       MCw = 1.0;
//       threepi_etapi.Start(*GetPluto(), *GetGeant());   // (pluto tree, n part in pluto per event)
//       MCw = TrueAnalysis_threepi_etapi();
//       MCw *= 1.5;

//       MC_weight = true;
//        etapr_10gTrue.Start(*GetPluto(), *GetGeant());   // (pluto tree, n part in pluto per event)
//        TrueAnalysis_etapr10g();
//        MCw = etapr_10gTrue.GetWeight();
        RandomTime();
    }
    if(!MC)
        Tagger_corr();
    tagger_min =-84.5;
    tagger_max = 83.5;

    for ( Int_t tag = 0; tag < GetTagger()->GetNTagged(); tag++){
        if(MC_weight)
            tag_BeamE->FillWeighted(GetTagger()->GetTaggedEnergy(tag), MCw);
        else
            tag_BeamE->Fill(GetTagger()->GetTaggedEnergy(tag),GetTagger()->GetTaggedTime(tag));
    }
    Time_corr();
    theta_corr();
    Energy_corr();

    ClustersInTime.clear();
    ClustersInTime.resize(0);
    UInt_t InTime = 0;
    UInt_t AllTimes = 0;
    UInt_t TAPS_clusters = 0;
    UInt_t CB_clusters = 0;

    eight_clusters = false;


    Double_t CB_avg_time = -1000.;
    Double_t CB_avg_time_tmp = 0.;
    ESum    = 0;
    ESum_MC = 0;

    // Here change MC Energy Sum
    // First step: count clusters and add CB clusters inside time window. Obtain CB time average and Energy
    for( Int_t i = 0; i < GetTracks()->GetNTracks(); i++ ){
        AllTimes++; // counting nr of CB and TAPS clusters
        if(GetTracks()->HasCB(i)){
            time_clusters_CB.Fill( GetTracks()->GetTime(i), GetTracks()->GetCentralCrystal(i) );
            ESum += GetTracks()->GetClusterEnergy(i);
            if( TMath::Abs( GetTracks()->GetTime(i) ) < 30.0 ){                               
                    CB_avg_time_tmp += GetTracks()->GetClusterEnergy(i)*GetTracks()->GetTime(i);
            }
        }
    }

    time_nr_AllClusters.Fill(AllTimes);
    CB_avg_time = CB_avg_time_tmp/ESum;
    if(MC)
        CB_avg_time = 0.0;

    // Here put in the condition that CB Energy Sum has to be fulfilled. For MC one has to consider also which photons contribute since in July 2014 data last 16 channels were off
    if(MC_weight){
        Double_t ESum_threshold = pRandoms->Gaus(540.,52.);
        ESum_MC = Get_ESumMC(ESum);
        if( (ESum_MC < ESum_threshold)) return;

        CB_EnergySum->FillWeighted(ESum_MC, MCw);
        CB_EnergySum_2->Fill(ESum_MC,MCw);

    }
    else{
        CB_EnergySum->Fill(ESum);
    }

    Double_t dt;
    // Second step: count clusters and add CB clusters inside time window. Obtain CB time average and Energy
    for( Int_t i = 0; i < GetTracks()->GetNTracks(); i++ ){
        if(GetTracks()->HasCB(i)){
            dt = GetTracks()->GetTime(i) - CB_avg_time;
            time_clusters_CBavg_CBdiff.Fill(dt, GetTracks()->GetCentralCrystal(i) );
            if( TMath::Abs( dt ) < 10.0 ){
                if(!sixg_cand->IsInside(GetTracks()->GetClusterEnergy(i)+10., GetTracks()->GetTheta(i))){
                    CB_clusters++;
                    InTime++;
                    ClustersInTime.push_back(i);
                }
            }
        }
    }

    for( Int_t i = 0; i < GetTagger()->GetNTagged(); i++ ){
        dt = GetTagger()->GetTaggedTime(i) - CB_avg_time;
        GetTagger()->SetTaggedTime( i, dt );
        time_clusters_CBavg_EPT.Fill(dt, GetTagger()->GetTaggedChannel(i));
    }

    TOF_CB = -100.;

    Double_t radnm;
    std::vector<int> TAPS_cl;
    TAPS_cl.resize(0);

 // Second step: count and add TAPS clusters in time window
    for( Int_t j = 0; j < GetTracks()->GetNTracks(); j++ ){
        if(GetTracks()->HasTAPS(j)){
            time_clusters_TAPS.Fill( GetTracks()->GetTime(j), GetTracks()->GetCentralCrystal(j) );
            radnm = 1.45/TMath::Cos( GetTracks()->GetThetaRad(j) );
            TOF_CB = (( CB_avg_time -(PbWO4[GetTracks()->GetCentralCrystal(j)])*GetTracks()->GetTime(j) )/radnm) - TAPS_CB_toff[GetTracks()->GetCentralCrystal(j)];

            if(MC){
                TOF_CB = -1*(( CB_avg_time - GetTracks()->GetTime(j) )/radnm);
            }
            time_clusters_CBavg_TAPSdiff.Fill(TOF_CB, GetTracks()->GetCentralCrystal(j));
            if( ( TOF_CB > -8.0) &&  ( TOF_CB < 2.0 ) ){
                if(!sixg_cand->IsInside(GetTracks()->GetClusterEnergy(j)+10., GetTracks()->GetTheta(j))){
                    TAPS_clusters++;
                    InTime++;
                    ClustersInTime.push_back(j);
                    TAPS_cl.push_back(j);
                }
            }
        }
    }

    if(ClustersInTime.size() < 7) return;
    time_nr_ClustersinTime.Fill(ClustersInTime.size());

    //  Third step: Identify the proton cluster in TAPS

    nrprotons = 0;
    iprtrack = 10000;

    // third step- select proton cluster (only TAPS are allowed)
    if(TAPS_cl.size() == 1){
       Int_t q = TAPS_cl[0];
       radnm = 1.45/TMath::Cos( GetTracks()->GetThetaRad(q) );
       if(MC){
           TOF_CB = -1*(( CB_avg_time - GetTracks()->GetTime(q) )/radnm);
           TOF_CB_proton = TOF_CB;
       }
       else{
           TOF_CB = (( CB_avg_time -(PbWO4[GetTracks()->GetCentralCrystal(q)])*GetTracks()->GetTime(q) )/radnm) - TAPS_CB_toff[GetTracks()->GetCentralCrystal(q)];
           TOF_CB_proton = TOF_CB;
       }
       p_E_v_TOF_TAPS_1cl->Fill( TOF_CB , GetTracks()->GetClusterEnergy(q));
       nrprotons++;
       iprtrack = q;
    }
    else{ // more than one cluster
        UInt_t pr_cand_temp = 1000;
        Double_t time_pr_temp = 1000;
        for(UInt_t p = 0; p < TAPS_cl.size(); p++){
            UInt_t q = TAPS_cl[p];
            radnm = 1.45/TMath::Cos( GetTracks()->GetThetaRad(q) );
            if(MC){
                TOF_CB = -1*(( CB_avg_time - GetTracks()->GetTime(q) )/radnm);
            }
            else{
                TOF_CB = (( CB_avg_time -(PbWO4[GetTracks()->GetCentralCrystal(q)])*GetTracks()->GetTime(q) )/radnm) - TAPS_CB_toff[GetTracks()->GetCentralCrystal(q)];
            }
            p_E_v_TOF_CB_2PrID->Fill(TOF_CB, GetTracks()->GetClusterEnergy(q) );
            if(TOF_CB < time_pr_temp){
                pr_cand_temp = q;
                time_pr_temp = TOF_CB;
                TOF_CB_proton = TOF_CB;
            }
        }
        nrprotons++;
        iprtrack = pr_cand_temp;

        if(iprtrack != 1000)
            p_E_v_TOF_CB_best->Fill(TOF_CB_proton, GetTracks()->GetClusterEnergy(iprtrack) );
    }

    if(iprtrack != 1000){
        if(MC_weight)
            p_E_v_TOF_CB_All_proton->FillWeighted(TOF_CB_proton, GetTracks()->GetClusterEnergy(iprtrack),MCw);
        else
            p_E_v_TOF_CB_All_proton->Fill(TOF_CB_proton, GetTracks()->GetClusterEnergy(iprtrack) );
    }

    FinalClusterSelection.resize(0);
    FinalClusterSelection = ClustersInTime;

    time_nr_FinalClusterSel.Fill(FinalClusterSelection.size());

    TLorentzVector IM6_vec;
    IM6_vec.SetPxPyPzE(0., 0., 0., 0.);
    if(FinalClusterSelection.size() == 7){
        for(UInt_t p = 0; p < FinalClusterSelection.size() ; p++){
            UInt_t q = FinalClusterSelection[p];
            if(q != iprtrack){
                IM6_vec += GetTracks()->GetVector(q);
                if(MC_weight)
                    six_rec_EvTh_6g->FillWeighted(GetTracks()->GetClusterEnergy(q),GetTracks()->GetTheta(q), MCw);
                else
                    six_rec_EvTh_6g->Fill(GetTracks()->GetClusterEnergy(q),GetTracks()->GetTheta(q));
            }
        }
    }

    if(FinalClusterSelection.size() == 8){
        for(UInt_t p = 0; p < FinalClusterSelection.size() ; p++){
            UInt_t q = FinalClusterSelection[p];
            if(q != iprtrack){
                IM6_vec += GetTracks()->GetVector(q);
                if(MC_weight)
                    six_rec_EvTh_7g->FillWeighted(GetTracks()->GetClusterEnergy(q),GetTracks()->GetTheta(q), MCw);
                else
                    six_rec_EvTh_7g->Fill(GetTracks()->GetClusterEnergy(q),GetTracks()->GetTheta(q));
            }
        }
    }


   if((FinalClusterSelection.size() == 8) && (TAPS_cl.size() >= 2)){
       IM6_vec.SetPxPyPzE(0.,0.,0.,0.);
       eight_clusters = true;
       Int_t imin_E_cluster = 0;
       Double_t Emin = 1.0e6;
       for(UInt_t p = 0; p < TAPS_cl.size(); p++){
           UInt_t q = TAPS_cl[p];
           if(q != iprtrack){
               if(GetTracks()->GetClusterEnergy(q) < Emin){
                   Emin = GetTracks()->GetClusterEnergy(q);
                   imin_E_cluster = q;
               }
           }
       }
       // here re-set the FinalClusterSelection
       FinalClusterSelection.resize(0);
       for(int icl = 0; icl < ClustersInTime.size(); icl++){
           UInt_t q = FinalClusterSelection[icl];
           if(imin_E_cluster != q){
             FinalClusterSelection.push_back(q);
             if(q != iprtrack){
                IM6_vec += GetTracks()->GetVector(q);
             }
           }
       }

   }


    Double_t mass_hyp = IM6_vec.M();

    //    Kinfit_test();  // runs kinematical fit with true observables- for testing purposes

    MMp_vec.SetPxPyPzE(0., 0., 0., 0.);

    if( (iprtrack == 10000) || (nrprotons == 0) ) return;
    if( mass_hyp < 600. || mass_hyp > 1200.) return;

    // correction of proton from comparing true to reconstructed proton theta for MC eta prime'
    Double_t dth2 = TAPSth_corr[GetTracks()->GetCentralCrystal(iprtrack)];
    Double_t th_pr = GetTracks()->GetTheta(iprtrack) + dth2;
    tracks->SetTheta(iprtrack, th_pr);

    if(MC){
        Double_t ztrue   =  etapr_6gTrue.GetTrueVertex().Z();
        Double_t R_true  =  (145.7-ztrue)*(TMath::Tan(etapr_6gTrue.GetTrueProtonLV().Theta()));
        Double_t R_rec   =  R_TAPS[GetTracks()->GetCentralCrystal(iprtrack)];
        Double_t dR      =  R_rec - R_true;

        true_six_dth_vs_th_p->Fill(GetTracks()->GetTheta(iprtrack),GetTracks()->GetTheta(iprtrack) - etapr_6gTrue.GetTrueProtonLV().Theta()*TMath::RadToDeg(), MCw);
        true_six_dR_vs_R_p->Fill( R_rec, dR, MCw);
        true_six_dR_vs_det_p->Fill(GetTracks()->GetCentralCrystal(iprtrack), dR , MCw);
    }


    p_th_v_E->Fill(GetTracks()->GetClusterEnergy(iprtrack),GetTracks()->GetTheta(iprtrack));
    proton_vec = GetTracks()->GetVector(iprtrack, pdgDB->GetParticle("proton")->Mass()*1000);

    // Now construct missing mass calc for proton with tagger energies.
    for ( Int_t tag = 0; tag < GetTagger()->GetNTagged(); tag++){
        MMp_vec.SetPxPyPzE(-proton_vec.Px(), -proton_vec.Py(), GetTagger()->GetTaggedEnergy(tag) - proton_vec.Pz(), GetTagger()->GetTaggedEnergy(tag) + pdgDB->GetParticle("proton")->Mass()*1000 - proton_vec.E());
        if(MC_weight)
            p_MM->FillWeighted(MMp_vec.M(), MCw);
        else
            p_MM->Fill(MMp_vec.M(), GetTagger()->GetTaggedTime(tag));
    }

//    CB_TAPS_boundary();

    ClustersInTime.resize(0);
    ClustersInTime = FinalClusterSelection;
    FinalClusterSelection.resize(0);

//    if( (ClustersInTime.size() == 5) && (nrprotons > 0) )
//        fourgAnalysis(iprtrack);

    if( ClustersInTime.size() == 7 && (nrprotons > 0) )
        sixgAnalysis( iprtrack );

    if( ClustersInTime.size() == 11 && (nrprotons > 0) )
        tengAnalysis(iprtrack);
}

void	AdlarsonPhysics::ProcessScalerRead()
{
//    hist_eta.ScalerReadCorrection(Double_t(scalers->GetScaler(0))/scalers->GetScaler(1));
//    hist_etap.ScalerReadCorrection(Double_t(scalers->GetScaler(0))/scalers->GetScaler(1));
}


Bool_t	AdlarsonPhysics::Init(const char* configFile){

   std::string         line;

   std::ifstream file_etapr_leg("configfiles/data/etaprime_Legendrecoeff_effcorr_eta2g.txt");
//   std::ifstream file_etapr_leg("configfiles/data/etaprime_Legendrecoeff.txt");       // gives diff x-section
   if(file_etapr_leg){
       while(std::getline(file_etapr_leg, line)){
           std::string  buffer;
           std::stringstream ss;
           ss << line;
           while(std::getline(ss, buffer, '\t')){
                Legendre.push_back(std::stod(buffer));
           }
       }
   }
   else{
       std::cout << "etaprime_Legendrecoeff.txt not found !" << std::endl;
       exit(0);
   }
   CBtime_corr.resize(0);
   std::ifstream fileCBtime("configfiles/data/time_corr_CB.txt");
   if(fileCBtime){
        std::getline(fileCBtime, line);
        std::string         buffer1;
        std::stringstream   ss0;
        ss0 << line;
        while (std::getline(ss0, buffer1, '\t'))
        {
           CBtime_corr.push_back(std::stod(buffer1));
        }
        fileCBtime.close();
   }
    else
        for(int i = 0; i < 720; i++)
            CBtime_corr.push_back(0.0);


   CBtime_sgm.resize(0);
   std::ifstream fileCBsgm("configfiles/data/time_corr_CB_sgm.txt");
   if(fileCBsgm){
        std::getline(fileCBsgm, line);
        std::string         buffer1;
        std::stringstream   ss0;
        ss0 << line;
        while (std::getline(ss0, buffer1, '\t'))
        {
           CBtime_sgm.push_back(std::stod(buffer1));
        }
        fileCBsgm.close();
   }
    else
        for(int i = 0; i < 720; i++)
            CBtime_sgm.push_back(0.0);

    TAPS_CB_toff.resize(0);
    std::ifstream fileCBTAPS_off("configfiles/data/TAPS_CB_offsets.txt");
    if(fileCBTAPS_off){
    std::getline(fileCBTAPS_off, line); // first line gives explanation
    std::getline(fileCBTAPS_off, line);
    std::string         buffer1;
    std::stringstream   ss0;
    ss0 << line;
    while (std::getline(ss0, buffer1, '\t'))
    {
       TAPS_CB_toff.push_back(std::stod(buffer1));
    }
    fileCBTAPS_off.close();
   }
    else
        for(int i = 0; i < 438; i++)
            TAPS_CB_toff.push_back(0.0);

   CBgain.resize(0);
   std::ifstream fileCB("configfiles/data/CB_lingain.txt");
   std::getline(fileCB, line);
   std::string         buffer2;
   std::stringstream   ss;
   ss << line;
   while (std::getline(ss, buffer2, '\t'))
   {
       CBgain.push_back(std::stod(buffer2));
   }
   fileCB.close();

   TAPSgain.resize(0);
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

   TAPSth_corr.resize(0);

   std::ifstream fileTAPSthcorr("configfiles/corr/TAPS_theta_etapr_corr.txt");
   if(fileTAPSthcorr){
       std::getline(fileTAPSthcorr, line); // first line gives explanation
       std::getline(fileTAPSthcorr, line);
       std::string         buffer4;
       std::stringstream   ss3;
       ss3 << line;
       while (std::getline(ss3, buffer4, '\t'))
       {
           TAPSth_corr.push_back(std::stod(buffer4));
       }

       fileTAPSthcorr.close();
   }
   else
       for(int i = 0; i < 438; i++)
           TAPSth_corr.push_back(0.0);


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

   PbWO4.resize(0);
   for(int iPbWO = 0; iPbWO < 438; iPbWO++){
       if( iPbWO < 12 ) PbWO4.push_back(-1);
       else if( iPbWO > 72 && iPbWO < 85)
           PbWO4.push_back(-1);
       else if( iPbWO > 145 && iPbWO < 158)
           PbWO4.push_back(-1);
       else if( iPbWO > 218 && iPbWO < 231)
           PbWO4.push_back(-1);
       else if( iPbWO > 291 && iPbWO < 304)
           PbWO4.push_back(-1);
       else if( iPbWO > 364 && iPbWO < 377)
           PbWO4.push_back(-1);
       else
           PbWO4.push_back(1);
   }

   Int_t ring3_4_out[70] = {24, 25, 27, 28, 39, 41, 42, 43, 72, 73, 74, 77, 80, 81, 82, 636, 638, 639, 643, 644, 646, 647, 676, 677, 678, 680, 690, 693, 694, 695, 18, 21, 22, 23, 44, 45, 46, 66, 69, 70, 71, 75, 76, 78, 79, 83, 84, 85, 86, 88, 631, 633, 634, 635, 637, 640, 641, 642, 645, 648, 649, 651, 652, 673, 674, 675, 696, 697, 699, 700};

   ring3_or_ring4_CB.resize(0);
   bool el;
   for(Int_t iCB = 0; iCB < 720; iCB++){
       el = false;
       for(Int_t j = 0; j < 70; j++){
            if(iCB == ring3_4_out[j])
                el = true;
       }

       if(el)
           ring3_or_ring4_CB.push_back(1);
       else
           ring3_or_ring4_CB.push_back(0);
   }

   // structure of file: groups of 8 CB elements (0-89), cumulative distribution for -25, 0, 25, 50, 100. If elements -1 then no time jumps.
   // use of random value in RandomTime to determine if time jump should be implemented

   string sFile;
   if(MCJuly14)
        sFile = "configfiles/data/CB_time_jump_July14.txt";
   else // MCOctDec14
        sFile = "configfiles/data/CB_time_jump_OctDec14.txt";

   std::ifstream file_MC_t_jump(sFile);
   if(file_MC_t_jump){
       int group;
       std::vector<Double_t> t_jump;
       for(int igr = 0; igr < 90; igr++){
           t_jump.resize(5);
           file_MC_t_jump >> group >> t_jump[0] >> t_jump[1]  >> t_jump[2] >> t_jump[3] >> t_jump[4] ;
           MC_time_jump.insert(time_pair(group,t_jump));
       }
   }

      std::cout << "CB gain correction applied" << std::endl;
      std::cout << "TAPS gain correction applied" << std::endl;
//    std::cout << "No Init function specified for this class." << std::endl;

    return kTRUE;
}


void AdlarsonPhysics::fourgAnalysis(UInt_t ipr)
{
    for(int tag = 0; tag < GetTagger()->GetNTagged(); tag++){
        if((GetTagger()->GetTaggedTime(tag) < tagger_min) || (GetTagger()->GetTaggedTime(tag) > tagger_max)) continue;
        Double_t chi2_min = 1.0e6;
        Double_t prob_min = 1.0e6;
        std::vector<Int_t> set;     // particle conf: tagger, proton, photon 1, 2, 3, 4, 5, 6;
        std::vector<double> obs;
        std::vector<double> unc;

        detnr.resize(0);
        CB_region.resize(0);

//        IM4g_vec.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
//        IM4g_fit.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);

        photons_rec.resize(0);
        photons_fit.resize(0);

        for(UInt_t i = 0; i < nPhotons_four; i++){
            photons_rec[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
            photons_fit[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
        }

        // beam energy
        set.resize(0);
        set.push_back(tag);
        beam4g.SetFromVector( GetTagger()->GetVector(tag) );

        obs.resize(0);
        unc.resize(0);
        unc = Get_unc(3,0,obs);
        beam4g.Smear(unc, 2);
        beam4g.APLCONSettings();

        // proton
        set.push_back(ipr);
        proton4g.SetFromVector( GetTracks()->GetVector(ipr) );
        obs.resize(0);
        obs.push_back(GetTracks()->GetCentralCrystal(ipr));
        unc.resize(0);
        unc = Get_unc( 2, 2, obs);
        proton4g.Smear(unc, 1);

        Int_t n_photons = 0;
        for ( UInt_t jgam = 0; jgam < ClustersInTime.size() ; jgam++ ){
            UInt_t kgam = ClustersInTime[jgam];
            if( kgam != ipr ){ // the id proton cluster
                set.push_back(kgam);
                photons_rec[n_photons] = GetTracks()->GetVector(kgam);
                Photons_four[n_photons].SetFromVector( GetTracks()->GetVector(kgam));
                detnr.push_back(GetTracks()->GetCentralCrystal(kgam));

                if( GetTracks()->HasCB(kgam) ){
                    obs.resize(0);
                    obs.push_back(Photons_four[n_photons].Ek);
                    obs.push_back(Photons_four[n_photons].Theta);
                    unc.resize(0);
                    unc = Get_unc( 1, 1, obs);
                    Photons_four[n_photons].Smear(unc, 0);
                }
                else{ // ( GetTracks()->HasTAPS(k) )
                    obs.resize(0);
                    obs.push_back(Photons_four[n_photons].Ek);
                    obs.push_back(Photons_four[n_photons].Theta);
                    unc.resize(0);
                    unc = Get_unc( 2, 1, obs);
                    Photons_four[n_photons].Smear(unc, 0);
                }
                n_photons++;
            }
        }
        const APLCON::Result_t& result4g = kinfit4g.DoFit();
        if(result4g.Status == APLCON::Result_Status_t::Success  && (result4g.ChiSquare <  chi2_min)){
            chi2_min = result4g.ChiSquare;
            prob_min = result4g.Probability;
        }

        if( (prob_min < 0.01) || ( TMath::Abs(prob_min - 1.0e6) < 1.0e-4) ) continue;

//        four_fit_chi2->Fill( result4g.ChiSquare, GetTagger()->GetTaggedTime(tag));
//        four_fit_pdf->Fill( result4g.Probability, GetTagger()->GetTaggedTime(tag));

        for(UInt_t i = 0; i < nPhotons_four; i++){
            IM4g_vec += photons_rec[i];
            photons_fit[i] = FitParticle::Make(Photons_four[i], 0);
            IM4g_fit += photons_fit[i];
        }

//        four_rec_IM->Fill(IM4g_vec.M());
//        four_fit_IM->Fill(IM4g_fit.M());

        for ( Int_t jgam = 0; jgam < 4 ; jgam++ )
            if( GetTracks()->HasCB(set[jgam+2]) )
              CB_region.push_back(true);
            else
              CB_region.push_back(false);

        Double_t sigma_eta = 25.0;
        Double_t sigma_pi0 = 15.0;

        std::vector<Double_t> chi2_etapi0;
        std::vector<Double_t>  chi2_2pi0;
        chi2_etapi0.resize(0);
        chi2_2pi0.resize(0);


        Double_t chi2min_etapi0     = 1.0e6;
        Double_t chi2min_2pi0       = 1.0e6;

        Double_t chi2_etapi0_tmp    = 1.0e6;
        Double_t chi2_2pi0_tmp      = 1.0e6;

        int imin_etapi0 = -1;
        int imin_2pi0   = -1;
        int idet;

        for(int i = 0; i < 9; i++){
            if( i < 6 ){
                chi2_etapi0_tmp =  TMath::Power( ( ( photons_rec[perm4g[i][0]] + photons_rec[perm4g[i][1]] ).M() - MASS_ETA) / sigma_eta , 2 ) + TMath::Power( ( ( photons_rec[perm4g[i][2]] + photons_rec[perm4g[i][3]] ).M() - MASS_PI0 ) / sigma_pi0 , 2 );
                chi2_etapi0.push_back(chi2_etapi0_tmp);
                if( chi2_etapi0[i] < chi2min_etapi0 ){
                    chi2min_etapi0 = chi2_etapi0[i];
                    imin_etapi0 = i;
                }
            }
            else{
                chi2_2pi0_tmp = TMath::Power( ( (photons_rec[perm4g[i][0]] + photons_rec[perm4g[i][1]]).M() - MASS_PI0) / sigma_pi0 , 2 ) + TMath::Power( ( (photons_rec[perm4g[i][2]] + photons_rec[perm4g[i][3]]).M() - MASS_PI0) / sigma_pi0 , 2 );
                chi2_2pi0.push_back(chi2_2pi0_tmp);
                if( chi2_2pi0[i-6] < chi2min_2pi0 ){
                    chi2min_2pi0 = chi2_2pi0[i-6];
                    imin_2pi0 = i;
                }
            }
        }
//        four_fit_PDF_etapi_v_2pi->Fill(TMath::Prob(chi2min_etapi0,1), TMath::Prob(chi2min_2pi0,1), GetTagger()->GetTaggedTime(tag) );

        if( (chi2min_2pi0 < chi2min_etapi0) && ( TMath::Prob(chi2min_2pi0,1) > 0.01) ){ // 2pi0
            for(int j = 0; j < 4; j++){

//                four_fit_best_2pi0_pi_E_v_th->Fill(photons_fit[j].E(), photons_fit[j].Theta()*TMath::RadToDeg(), GetTagger()->GetTaggedTime(tag));

//                int imass = int(j/2);
//                idet = perm4g[imin_2pi0][j];
//                int nEN =  int((photons_rec[perm4g[ imin_2pi0 ][ j ]]).E()/20.0);
//                int nTH =  75*int((photons_rec[perm4g[ imin_2pi0 ][ j]]).Theta()*TMath::RadToDeg()/1.0);
//                int nBIN = nEN + nTH;
//                double mgg = (photons_rec[perm4g[ imin_2pi0 ][ imass*2 ]] + photons_rec[perm4g[imin_2pi0][imass*2+1]] ).M();

//                four_fit_mgg_v_eth->Fill( nBIN, mgg, GetTagger()->GetTaggedTime(tag) );
//                if(CB_region[perm4g[ imin_2pi0 ][ j ]]){
//                    IMgg_v_det_2pi0_CB->Fill( (photons_rec[perm4g[ imin_2pi0 ][ imass*2 ]] + photons_rec[perm4g[imin_2pi0][imass*2+1]] ).M(), detnr[idet], GetTagger()->GetTaggedTime(tag));
//                }
//                else{
//                    IMgg_v_det_2pi0_TAPS->Fill((photons_rec[perm4g[ imin_2pi0 ][ imass*2 ]] + photons_rec[perm4g[imin_2pi0][imass*2+1]] ).M(), detnr[idet], GetTagger()->GetTaggedTime(tag) );

//                }
            }
        }
        else if( (chi2min_etapi0 < chi2min_2pi0  ) && ( TMath::Prob(chi2min_etapi0,1) > 0.01)){   //etapi0

            // -------------------------------++++++++++++++++++++++++++++++++++++++++++++----------------------------------

            continue;

            // -------------------------------++++++++++++++++++++++++++++++++++++++++++++----------------------------------

            for(int j = 0; j < 4; j++){
//                four_fit_best_etapi0_pi_E_v_th->Fill(photons_fit[j].E(), photons_fit[j].Theta()*TMath::RadToDeg(), GetTagger()->GetTaggedTime(tag));

                idet = perm4g[imin_etapi0][j];
                if(CB_region[perm4g[ imin_etapi0 ][ j ]]){

                }
            }
        }
    }
}

void AdlarsonPhysics::sixgAnalysis(UInt_t ipr){

    // First reset the vectors
//    vec_weight.resize(0);
//    vec_nclusters.resize(0);
//    vec_X.resize(0);
//    vec_Y.resize(0);
//    vec_M_eta2pi.resize(0);
//    vec_taggerenergy.resize(0);
//    vec_fitted_pr_th.resize(0);
//    vec_fitted_pr_e.resize(0);
//    fitted_zvx.resize(0);
//    vec_X.resize(0);
//    vec_Y.resize(0);
//    vec_pdf_eta2pi.resize(0);
//    vec_pdf_3pi.resize(0);
//    vec_pdf_etapr.resize(0);
//    vec_costh_epr_cm.resize(0);
//    vec_DP_005.resize(0);
//    vec_DP_010.resize(0);
//    vec_DP_015.resize(0);
//    vec_M_eta2pi.resize(0);
//    vec_M_etapi1.resize(0);
//    vec_M_etapi2.resize(0);
//    vec_M_pipi.resize(0);
//    vec_costh_epr_cm_pr.resize(0);
//    vec_X_pr.resize(0);
//    vec_Y_pr.resize(0);
//    vec_DP_005_pr.resize(0);
//    vec_DP_010_pr.resize(0);
//    vec_DP_015_pr.resize(0);
//    vec_M_etapi_pr1.resize(0);
//    vec_M_etapi_pr2.resize(0);
//    vec_M_pipi_pr.resize(0);

    for(Int_t tag = 0; tag < GetTagger()->GetNTagged(); tag++){
        if((GetTagger()->GetTaggedTime(tag) < tagger_min) || (GetTagger()->GetTaggedTime(tag) > tagger_max)) continue;
        Double_t chi2_min = 1.0e6;
        Double_t prob_min = 1.0e6;

        IM6g_vec.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
        six_time_TaggedTime.Fill(GetTagger()->GetTaggedTime(tag));

        std::vector<Int_t> set;     // particle conf: tagger, proton, photon 1, 2, 3, 4, 5, 6;
        std::vector<Int_t> set_min; // best particle configuration
        set.resize(0);
        set_min.resize(0);
        Is_CB.resize(0);

        photons_rec.resize(0);
        photons_fit.resize(0);
        photons_fit_final.resize(0);
        for(UInt_t i = 0; i < nPhotons_six; i++){
            photons_rec[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
            photons_fit[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
            photons_fit_final[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
        }

        for(UInt_t i = 0; i < 3; i++){
            g[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
            h[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
            rc[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
            rc_sig[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
        }

        for ( UInt_t jgam = 0; jgam < ClustersInTime.size() ; jgam++ ){
            UInt_t kgam = ClustersInTime[jgam];
            if( kgam != ipr ){ // the id proton cluster
                IM6g_vec += GetTracks()->GetVector(kgam);
            }
        }

        set.push_back(tag);

        beam.SetFromVector( GetTagger()->GetVector(tag) );
        obs.resize(0);
        unc.resize(0);
        unc = Get_unc_R(3,0,obs);
        beam.Smear(unc, 2);
        beam.APLCONSettings();

        set.push_back(ipr);

        Double_t test2 = (GetTagger()->GetVector(tag) - IM6g_vec).E();

        proton.SetFromValues( test2, R_TAPS[GetTracks()->GetCentralCrystal(ipr)],GetTracks()->GetVector(ipr).Phi() );
        obs.resize(0);
        obs.push_back(GetTracks()->GetCentralCrystal(ipr));
        unc.resize(0);
        unc = Get_unc_R( 2, 2, obs);
        proton.Smear_R(unc, 1);
        Is_CB.push_back(0);

        UInt_t n_photons = 0;
        for ( UInt_t jgam = 0; jgam < ClustersInTime.size() ; jgam++ ){
            UInt_t kgam = ClustersInTime[jgam];
            if( kgam != ipr ){ // the id proton cluster
                set.push_back(kgam);
                photons_rec.push_back(GetTracks()->GetVector(kgam));
//                IM6g_vec += photons_rec[n_photons];

                if( GetTracks()->HasCB(kgam) ){
                    Photons_six[n_photons].SetFromValues( GetTracks()->GetVector(kgam).E(),GetTracks()->GetVector(kgam).Theta(),GetTracks()->GetVector(kgam).Phi() );
                    Is_CB.push_back(1);
                    obs.resize(0);
                    obs.push_back(Photons_six[n_photons].Ek);
                    obs.push_back(Photons_six[n_photons].Theta);
                    unc.resize(0);
                    unc = Get_unc_R( 1, 1, obs);
                    Photons_six[n_photons].Smear_R(unc, 0);
                }
                else{   // ( GetTracks()->HasTAPS(kgam) )
                    Is_CB.push_back(0);
                    Photons_six[n_photons].SetFromValues( GetTracks()->GetVector(kgam).E(), R_TAPS[GetTracks()->GetCentralCrystal(kgam)],GetTracks()->GetVector(kgam).Phi() );
                    obs.resize(0);
                    obs.push_back(Photons_six[n_photons].Ek);
                    obs.push_back(Photons_six[n_photons].Theta);
                    obs.push_back(GetTracks()->GetCentralCrystal(kgam));
                    unc.resize(0);
                    unc = Get_unc_R( 2, 1, obs);
                    Photons_six[n_photons].Smear_R(unc, 0);
                }
                n_photons++;
            }
        }

        if(MC_weight){
            six_rec_IM->FillWeighted(IM6g_vec.M(), MCw);
            six_rec_IM_v_MMp->FillWeighted(IM6g_vec.M(), MMp_vec.M(), MCw);
        }
        else{
            six_rec_IM->Fill(IM6g_vec.M(), GetTagger()->GetTaggedTime(tag));
            six_rec_IM_v_MMp->Fill(IM6g_vec.M(), MMp_vec.M(), GetTagger()->GetTaggedTime(tag));
        }

        const APLCON::Result_t& result_min = kinfit.DoFit();
        if(result_min.Status == APLCON::Result_Status_t::Success){
            if( result_min.ChiSquare <  chi2_min ){
                chi2_min = result_min.ChiSquare;
                prob_min = result_min.Probability;
                set_min = set;
                NI6g->Fill(result_min.NIterations);
            }
        }

        if( set_min.size() == 0 ) continue;
        if( (prob_min < 0.01) || ( TMath::Abs(prob_min - 1.0e6) < 1.0e-4) ) continue;
        detnr.resize(0);


//        TLorentzVector trueobs;
//      trueobs = GetGeant()->GetTrueVector(1);

//        Double_t p_th_true = trueobs.Theta()*TMath::RadToDeg();
//        Double_t p_phi_true = trueobs.Phi()*TMath::RadToDeg();

//        Double_t p_th_fit = GetTracks()->GetTheta(iprtrack);
//        Double_t p_th_rec = GetTracks()->GetTheta(iprtrack);
//        Double_t p_phi_rec = GetTracks()->GetPhi(iprtrack);

//         six_fit_true_p_th_v_det->Fill(GetTracks()->GetCentralCrystal(iprtrack), p_th_true - p_th_fit );
//        six_fit_true_p_fi_v_det->Fill(GetTracks()->GetCentralCrystal(iprtrack), p_phi_true- p_phi_rec );


        // For the best combination with pdf >= 0.01 obtain the result.

        std::vector<double> obs_kfit;
        obs_kfit.resize(0);
        IM6g_vec.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
        TLorentzVector etap_fit(0.0, 0.0, 0.0, 0.0);
        TLorentzVector etap_fit_tmp(0.0, 0.0, 0.0, 0.0);

        for ( Int_t n_photons_min = 0; n_photons_min < 6 ; n_photons_min++ ){
//            photons_fit.push_back( FitParticle::Make(Photons_six[n_photons_min], 0));
            detnr.push_back( GetTracks()->GetCentralCrystal( set_min[n_photons_min+2]) );
            IM6g_vec += photons_rec[n_photons_min];
            obs_kfit.resize(0);
            string s = Form("Photon%0i[0]", n_photons_min+1); // Energy
            string t = Form("Photon%0i[1]", n_photons_min+1); // th
            string u = Form("Photon%0i[2]", n_photons_min+1); // phi
            string z = Form("v_z"); // z

            double E = result_min.Variables.at(s).Value.After;
            double th = result_min.Variables.at(t).Value.After;
            double ph = result_min.Variables.at(u).Value.After;
            double vx_z = result_min.Variables.at(z).Value.After;
            if(Is_CB[n_photons_min+1] == 1){
                obs_kfit = {E, E, th, ph};
                etap_fit_tmp = GetLVCorrForZ(obs_kfit, vx_z, Is_CB[n_photons_min+1], 0.0);
                photons_fit.push_back(etap_fit_tmp);
                etap_fit += etap_fit_tmp;
            }
            else{
                obs_kfit = {E, E, th, ph}; // th is R in TAPS

                etap_fit_tmp = GetLVCorrForZ(obs_kfit, vx_z, Is_CB[n_photons_min+1], 0.0);
                photons_fit.push_back(etap_fit_tmp);
                etap_fit += etap_fit_tmp;
            }
        }
        if(result_min.Status == APLCON::Result_Status_t::Success ){

//            TLorentzVector proton_fit = FitParticle::Make(proton, MASS_PROTON);
            // proton LV:
            double E    = result_min.Variables.at("Proton[0]").Value.After + MASS_PROTON;
            double P    = TMath::Sqrt( E*E - MASS_PROTON*MASS_PROTON );
            double R    = result_min.Variables.at("Proton[1]").Value.After;
            double ph   = result_min.Variables.at("Proton[2]").Value.After;
            double vx_z = result_min.Variables.at("v_z").Value.After;
            int    idet = 0;
            double mass = MASS_PROTON;(GetTagger()->GetVector(tag) - IM6g_vec).E();


            obs_kfit = {E, P, R, ph};
            TLorentzVector proton_fit = GetLVCorrForZ(obs_kfit, vx_z, idet, mass);


            if(MC_weight){
//                Double_t p_th_fit = TMath::ATan(R/(145.7-vx_z))*TMath::RadToDeg();
                Double_t p_th_true = etapr_6gTrue.GetTrueProtonLV().Theta()*TMath::RadToDeg();
//                Double_t p_th_fit = proton_fit.Theta()*TMath::RadToDeg();
//                Double_t p_th_fit = TMath::ATan(proton_fit.Theta()/145.7)*TMath::RadToDeg();
//                Double_t p_E_true = etapr_6gTrue.GetTrueProtonLV().E()*1.0e3 - MASS_PROTON;

                Double_t ztrue = etapr_6gTrue.GetTrueVertex().Z();
                true_z_after_fit->Fill(ztrue, GetTagger()->GetTaggedTime(tag));

                Double_t diff = vx_z - ztrue;
                true_six_fit_dz_v_z->Fill(ztrue, diff);
                true_six_fit_dz_v_p_th->Fill(p_th_true, diff);
                fit_six_fit_dz_v_z->Fill(vx_z, diff);
            }

            if(etap_fit.M() < 830.)
                continue;
            if( GetTagger()->GetTaggedChannel(tag) > 39 )
                continue;

            if(MC_weight){
                six_fit_chi2->FillWeighted(result_min.ChiSquare, MCw );
                six_fit_pdf->FillWeighted( result_min.Probability, MCw );

                if(result_min.Probability >0.01 && result_min.Probability <1.0)
                    six_fit_IM->FillWeighted( etap_fit.M(), MCw );

                Double_t ztrue = etapr_6gTrue.GetTrueVertex().Z();
                six_fit_IM_vz->FillWeighted( etap_fit.M(), ztrue, MCw );
//                six_fit_dthpr_vz->FillWeighted(p_th_true - proton_fit.Theta()*TMath::RadToDeg(),ztrue,MCw);

                proton_fit_e_v_th->FillWeighted(proton_fit.E() - MASS_PROTON, proton_fit.Theta()*TMath::RadToDeg(), MCw );
                p_E_v_TOF_after_kfit->FillWeighted(TOF_CB_proton, proton_fit.E() - MASS_PROTON, MCw );
                p_Erec_v_TOF_after_kfit->FillWeighted(TOF_CB_proton, GetTracks()->GetClusterEnergy(ipr) , MCw );

                six_fit_fitted_p_th_v_det->FillWeighted(GetTracks()->GetCentralCrystal(ipr), proton_fit.Theta()*TMath::RadToDeg()-GetTracks()->GetTheta(ipr) , MCw );
                six_fit_fitted_p_fi_v_det->FillWeighted(GetTracks()->GetCentralCrystal(ipr), proton_fit.Phi()*TMath::RadToDeg()-GetTracks()->GetPhi(ipr) , MCw );

                for(UInt_t igam_fit = 0; igam_fit < Photons_six.size(); igam_fit++)
                     six_fit_EvTh_g->FillWeighted(photons_fit[igam_fit].E(), photons_fit[igam_fit].Theta()*TMath::RadToDeg(), MCw );

                // true observable only in MC
//                six_fit_true_p_th_v_det->FillWeighted(GetTracks()->GetCentralCrystal(iprtrack), p_th_true - proton_fit.Theta()*TMath::RadToDeg(), MCw );
            }
            else{
                six_fit_chi2->Fill( result_min.ChiSquare, GetTagger()->GetTaggedTime(tag));
                six_fit_pdf->Fill( result_min.Probability, GetTagger()->GetTaggedTime(tag) );

                if(result_min.Probability >0.01 && result_min.Probability <1.0){
                    six_fit_IM->Fill( etap_fit.M(),GetTagger()->GetTaggedTime(tag) );

                    int ncl = 7;
                    if(eight_clusters)
                        ncl = 8;

                    six_fit_IM_ncl->Fill( ncl, etap_fit.M(), GetTagger()->GetTaggedTime(tag)  );
                }

                proton_fit_e_v_th->Fill(proton_fit.E() - MASS_PROTON, proton_fit.Theta()*TMath::RadToDeg(), GetTagger()->GetTaggedTime(tag));
                p_E_v_TOF_after_kfit->Fill(TOF_CB_proton, proton_fit.E()-MASS_PROTON, GetTagger()->GetTaggedTime(tag));
                p_Erec_v_TOF_after_kfit->Fill(TOF_CB_proton, GetTracks()->GetClusterEnergy(ipr) ,GetTagger()->GetTaggedTime(tag));

                six_fit_fitted_p_th_v_det->Fill(GetTracks()->GetCentralCrystal(ipr), proton_fit.Theta()*TMath::RadToDeg()-GetTracks()->GetTheta(ipr) , GetTagger()->GetTaggedTime(tag) );
                six_fit_fitted_p_fi_v_det->Fill(GetTracks()->GetCentralCrystal(ipr), proton_fit.Phi()*TMath::RadToDeg()-GetTracks()->GetPhi(ipr) , GetTagger()->GetTaggedTime(tag) );

                for(UInt_t igam_fit = 0; igam_fit < Photons_six.size(); igam_fit++)
                     six_fit_EvTh_g->Fill(photons_fit[igam_fit].E(), photons_fit[igam_fit].Theta()*TMath::RadToDeg(), GetTagger()->GetTaggedTime(tag) );
            }


           if(p_ToF_kfit_cut->IsInside(TOF_CB_proton, proton_fit.E() - MASS_PROTON))
               return;

           p_E_v_TOF_after_kfit_2->Fill(TOF_CB_proton, proton_fit.E()-MASS_PROTON, GetTagger()->GetTaggedTime(tag));

           sigma_eta = 20; sigma_pi0 = 12;

           Double_t    chi2min_eta2pi  = 1.0e6;
           Double_t    chi2min_3pi     = 1.0e6;
           Double_t    probmin_etapr  = 0.;
           Double_t    probmin_eta2pi  = 1.0e6;
           Double_t    probmin_3pi     = 1.0e6;
           std::vector<int> imin_eta2pi;
           std::vector<int> imin_3pi;
           imin_eta2pi.resize(0);
           imin_3pi.resize(0);

           std::vector<comb> threepi_comb;
           std::vector<comb> etatwopi_comb;
           threepi_comb.resize(0);
           etatwopi_comb.resize(0);

           GetBest6gCombination(sigma_eta, sigma_pi0, chi2min_eta2pi, chi2min_3pi, imin_eta2pi, imin_3pi, etatwopi_comb, threepi_comb );
           test_correct_hypothesis(probmin_etapr, probmin_eta2pi, probmin_3pi, set_min, imin_eta2pi, imin_3pi, etatwopi_comb, threepi_comb);

           if(MC_weight){
               six_fit_PDF_eta2pi_v_3pi->FillWeighted(probmin_eta2pi, probmin_3pi, MCw);
               if(probmin_etapr > 0.02 && probmin_etapr < 1.0 )
                    six_fit_PDF_eta2pi_v_3pi_2->FillWeighted(probmin_eta2pi, probmin_3pi, MCw);
               if(probmin_etapr > 0.04 &&  probmin_etapr < 1.0 )
                    six_fit_PDF_eta2pi_v_3pi_4->FillWeighted(probmin_eta2pi, probmin_3pi, MCw);
           }
           else{
               six_fit_PDF_eta2pi_v_3pi->Fill(probmin_eta2pi, probmin_3pi, GetTagger()->GetTaggedTime(tag));

               if(probmin_etapr > 0.02 && probmin_etapr < 1.0 )
                    six_fit_PDF_eta2pi_v_3pi_2->Fill(probmin_eta2pi, probmin_3pi, GetTagger()->GetTaggedTime(tag));
               if(probmin_etapr > 0.04 &&  probmin_etapr < 1.0 )
                    six_fit_PDF_eta2pi_v_3pi_4->Fill(probmin_eta2pi, probmin_3pi, GetTagger()->GetTaggedTime(tag));
           }

           etatwopi_comb.resize(0);
           threepi_comb.resize(0);

           // photons after kinematical fit only under energy and momentum constraint

           g[0] = photons_fit[imin_eta2pi[0]] + photons_fit[imin_eta2pi[1]];
           g[1] = photons_fit[imin_eta2pi[2]] + photons_fit[imin_eta2pi[3]];
           g[2] = photons_fit[imin_eta2pi[4]] + photons_fit[imin_eta2pi[5]];

           h[0] = photons_fit[imin_3pi[0]] + photons_fit[imin_3pi[1]];
           h[1] = photons_fit[imin_3pi[2]] + photons_fit[imin_3pi[3]];
           h[2] = photons_fit[imin_3pi[4]] + photons_fit[imin_3pi[5]];

           // reconstructed photons
           rc[0] = photons_rec[imin_3pi[0]] + photons_rec[imin_3pi[1]];
           rc[1] = photons_rec[imin_3pi[2]] + photons_rec[imin_3pi[3]];
           rc[2] = photons_rec[imin_3pi[4]] + photons_rec[imin_3pi[5]];

           rc_sig[0] = photons_rec[imin_eta2pi[0]] + photons_rec[imin_eta2pi[1]];
           rc_sig[1] = photons_rec[imin_eta2pi[2]] + photons_rec[imin_eta2pi[3]];
           rc_sig[2] = photons_rec[imin_eta2pi[4]] + photons_rec[imin_eta2pi[5]];

           double mass_shift = 1.3;

           TLorentzVector etap_fit_final(0.0, 0.0, 0.0, 0.0);
           TLorentzVector fin[3];
           Xfit = -2.0;
           Yfit = -2.0;
           Double_t m_etapi01_fit = 0;
           Double_t m_etapi02_fit = 0;
           Double_t m_2pi0_fit = 0;
           Int_t DP_binnr_fit020 = -100; //0.20
           Int_t DP_binnr_fit015 = -100; //0.15
           Int_t DP_binnr_fit010 = -100; //0.10
           Int_t DP_binnr_fit075 = -100; //0.075
           Int_t DP_binnr_fit005 = -100; //0.05

           Xfit_pr = -2.0;
           Yfit_pr = -2.0;
           Double_t m_etapi01_fit_pr = 0;
           Double_t m_etapi02_fit_pr = 0;
           Double_t m_2pi0_fit_pr = 0;
           Int_t DP_binnr_fit020_pr = -100; //0.20
           Int_t DP_binnr_fit015_pr = -100; //0.15
           Int_t DP_binnr_fit010_pr = -100; //0.10
           Int_t DP_binnr_fit075_pr = -100; //0.075
           Int_t DP_binnr_fit005_pr = -100; //0.05

           // Sergeys comparison plots
           if(probmin_eta2pi > 0.01 && probmin_eta2pi < 1.0){
               for(UInt_t igam_fit = 0; igam_fit < Photons_six.size(); igam_fit++)
                   etap_fit_final += photons_fit_final[igam_fit];

               fin[0] = photons_fit_final[0] + photons_fit_final[1];
               fin[1] = photons_fit_final[2] + photons_fit_final[3];
               fin[2] = photons_fit_final[4] + photons_fit_final[5];

               // Calculation of final results

               TLorentzVector fin_all = fin[0] + fin[1] + fin[2];
               diffbin =  diff_distr( GetTagger()->GetTaggedEnergy(tag), fin_all );

               m2pi0_metapi0( fin ,  m_etapi01_fit, m_etapi02_fit, m_2pi0_fit);
               bw = 0.2;
               DalitzPlot(fin, Xfit, Yfit, bw, DP_binnr_fit020);
               bw = 0.15;
               DalitzPlot(fin, Xfit, Yfit, bw, DP_binnr_fit015);
               bw = 0.1;
               DalitzPlot(fin, Xfit, Yfit, bw, DP_binnr_fit010);
               bw = 0.075;
               DalitzPlot(fin, Xfit, Yfit, bw, DP_binnr_fit075);
               bw = 0.05;
               DalitzPlot(fin, Xfit, Yfit, bw, DP_binnr_fit005);

               int n_clusters = 7;
               if(eight_clusters)
                    n_clusters = 8;


               if(!MC)
                    if( (GetTagger()->GetTaggedTime(tag) > 3.5) || (GetTagger()->GetTaggedTime(tag) < -4.5))
                        fWeight = -8./160.;
                    else
                        fWeight = 1.0;
               else
                   fWeight = MCw;

                fNclusters      = n_clusters;
                fTaggerenergy   = GetTagger()->GetTaggedEnergy(tag);
                fProton_E_fit   = proton_fit.E()-MASS_PROTON;
                fProton_th_fit  = proton_fit.Theta()*TMath::RadToDeg();
                fZ_vx_fit       = vx_z;
                fPeta2pi        = probmin_eta2pi;
                fP3pi           = probmin_3pi;
                fPetapr         = probmin_etapr;
                fMeta2pi        = etap_fit_final.M();
                fX              = Xfit;
                fY              = Yfit;
                fCosth_epr_cm   = diffbin;
                fMpipi          = m_2pi0_fit/1.0e3;
                fMetapi1        = m_etapi01_fit/1.0e3;
                fMetapi2        = m_etapi02_fit/1.0e3;
                fDP005          = DP_binnr_fit005;
                fDP075          = DP_binnr_fit075;
                fDP010          = DP_binnr_fit010;
                fDP015          = DP_binnr_fit015;

                if(probmin_etapr > 0.01 && probmin_etapr < 1.0){

                    TLorentzVector fin_metapr[3];
                    fin_metapr[0] = photons_fit_final_metapr[0] + photons_fit_final_metapr[1];
                    fin_metapr[1] = photons_fit_final_metapr[2] + photons_fit_final_metapr[3];
                    fin_metapr[2] = photons_fit_final_metapr[4] + photons_fit_final_metapr[5];

                    TLorentzVector fin_all_metapr = fin_metapr[0] + fin_metapr[1] + fin_metapr[2];
                    diffbin_pr =  diff_distr( GetTagger()->GetTaggedEnergy(tag), fin_all_metapr );

                    m_etapi01_fit_pr = 0.;
                    m_etapi02_fit_pr = 0.;
                    m_2pi0_fit_pr    = 0.;
                    m2pi0_metapi0( fin_metapr ,  m_etapi01_fit_pr, m_etapi02_fit_pr, m_2pi0_fit_pr);

                    bw = 0.20;
                    DalitzPlot(fin_metapr, Xfit_pr, Yfit_pr, bw, DP_binnr_fit020_pr);
                    bw = 0.15;
                    DalitzPlot(fin_metapr, Xfit_pr, Yfit_pr, bw, DP_binnr_fit015_pr);
                    bw = 0.10;
                    DalitzPlot(fin_metapr, Xfit_pr, Yfit_pr, bw, DP_binnr_fit010_pr);
                    bw = 0.075;
                    DalitzPlot(fin_metapr, Xfit_pr, Yfit_pr, bw, DP_binnr_fit075_pr);
                    bw = 0.05;
                    DalitzPlot(fin_metapr, Xfit_pr, Yfit_pr, bw, DP_binnr_fit005_pr);

                }

                fCosth_epr_cm_pr    = diffbin_pr;
                fXpr                = Xfit_pr;
                fYpr                = Yfit_pr;
                fDP005pr            = DP_binnr_fit005_pr;
                fDP075pr            = DP_binnr_fit075_pr;
                fDP010pr            = DP_binnr_fit010_pr;
                fDP015pr            = DP_binnr_fit015_pr;
                fMpipipr            = m_2pi0_fit_pr/1.0e3;
                fMetapi1pr          = m_etapi01_fit_pr/1.0e3;
                fMetapi2pr          = m_etapi02_fit_pr/1.0e3;

                tree2->Fill();

               if(MC_weight){
                    six_fit_IM_eta2pi0_b->FillWeighted( etap_fit_final.M()+mass_shift, MCw );
                    six_rec_IM_eta2pi->FillWeighted(IM6g_vec.M(), MCw);
                    six_phy_DP_020_P1->FillWeighted(etap_fit_final.M()+mass_shift, DP_binnr_fit020_pr, MCw);
                    six_phy_DP_015_P1->FillWeighted(etap_fit_final.M()+mass_shift, DP_binnr_fit015_pr, MCw);
                    six_phy_DP_010_P1->FillWeighted(etap_fit_final.M()+mass_shift, DP_binnr_fit010_pr, MCw);
                    six_phy_DP_075_P1->FillWeighted(etap_fit_final.M()+mass_shift, DP_binnr_fit075_pr, MCw);
                    six_phy_DP_005_P1->FillWeighted(etap_fit_final.M()+mass_shift, DP_binnr_fit005_pr, MCw);
               }
               else{
                    six_fit_IM_eta2pi0_b->Fill( etap_fit_final.M(), GetTagger()->GetTaggedTime(tag) );
                    six_rec_IM_eta2pi->Fill(IM6g_vec.M(), GetTagger()->GetTaggedTime(tag) );
                    six_phy_DP_020_P1->Fill(etap_fit_final.M(), DP_binnr_fit020_pr, GetTagger()->GetTaggedTime(tag));
                    six_phy_DP_015_P1->Fill(etap_fit_final.M(), DP_binnr_fit015_pr, GetTagger()->GetTaggedTime(tag));
                    six_phy_DP_010_P1->Fill(etap_fit_final.M(), DP_binnr_fit010_pr, GetTagger()->GetTaggedTime(tag));
                    six_phy_DP_075_P1->Fill(etap_fit_final.M(), DP_binnr_fit075_pr, GetTagger()->GetTaggedTime(tag));
                    six_phy_DP_005_P1->Fill(etap_fit_final.M(), DP_binnr_fit005_pr, GetTagger()->GetTaggedTime(tag));
               }

               if(probmin_eta2pi > 0.04 && probmin_eta2pi < 1.0){
                   if( probmin_3pi < 0.0075 ){
                        if(MC_weight)
                            six_fit_IM_eta2pi0_c->FillWeighted( etap_fit_final.M()+mass_shift, MCw );
                        else
                            six_fit_IM_eta2pi0_c->Fill( etap_fit_final.M(), GetTagger()->GetTaggedTime(tag) );

                        if( (probmin_etapr > 0.04) && (probmin_etapr < 1.0) ){
                            if(MC_weight){
                                six_fit_IM_eta2pi0_d->FillWeighted( etap_fit_final.M()+mass_shift, MCw );
                                for(int ietapr = 0; ietapr < 6; ietapr++){
                                   int nEN =  int(photons_rec[imin_eta2pi[ietapr]].E()/40.);
                                   int nTH =  30*int(photons_rec[imin_eta2pi[ietapr]].Theta()*TMath::RadToDeg()/4.0);
                                   int nBIN = nEN + nTH;
                                   six_rec_m6g_sig_v_eth->FillWeighted(nBIN, IM6g_vec.M(), MCw );

                                   if( GetTracks()->HasCB( set_min[imin_eta2pi[ietapr]+2]) ){
                                            IM6g_v_det_etaprrec_CB->FillWeighted( IM6g_vec.M(), detnr[imin_eta2pi[ietapr]], MCw);
                                            IM6g_v_det_etaprfit_CB->FillWeighted( etap_fit_final.M()+mass_shift, detnr[imin_eta2pi[ietapr]], MCw);
                                   }
                                   else{
                                            IM6g_v_det_etaprrec_TAPS->FillWeighted( IM6g_vec.M(), detnr[imin_eta2pi[ietapr]], MCw);
                                            IM6g_v_det_etaprfit_TAPS->FillWeighted( etap_fit_final.M(), detnr[imin_eta2pi[ietapr]], MCw);
                                       }
                                 }
                            }
                            else{
                                six_fit_IM_eta2pi0_d->Fill( etap_fit_final.M(), GetTagger()->GetTaggedTime(tag) );
                                for(int ietapr = 0; ietapr < 6; ietapr++){
                                    int nEN =  int(photons_rec[imin_eta2pi[ietapr]].E()/40.);
                                    int nTH =  30*int(photons_rec[imin_eta2pi[ietapr]].Theta()*TMath::RadToDeg()/4.0);
                                    int nBIN = nEN + nTH;
                                    six_rec_m6g_sig_v_eth->Fill(nBIN, IM6g_vec.M(), GetTagger()->GetTaggedTime(tag) );
                                    if( GetTracks()->HasCB( set_min[imin_eta2pi[ietapr]+2]) ){
                                            IM6g_v_det_etaprrec_CB->Fill( IM6g_vec.M(), detnr[imin_eta2pi[ietapr]], GetTagger()->GetTaggedTime(tag));
                                            IM6g_v_det_etaprfit_CB->Fill( etap_fit_final.M(), detnr[imin_eta2pi[ietapr]], GetTagger()->GetTaggedTime(tag));
                                    }
                                    else{
                                            IM6g_v_det_etaprrec_TAPS->Fill( IM6g_vec.M(), detnr[imin_eta2pi[ietapr]], GetTagger()->GetTaggedTime(tag));
                                            IM6g_v_det_etaprfit_TAPS->Fill( etap_fit_final.M(), detnr[imin_eta2pi[ietapr]], GetTagger()->GetTaggedTime(tag));
                                        }
                                    }
                                }
                            }
                        else{
                            if(MC_weight)
                                six_fit_IM_eta2pi0_e->FillWeighted( etap_fit_final.M()+mass_shift, MCw );
                            else
                                six_fit_IM_eta2pi0_e->Fill( etap_fit_final.M(), GetTagger()->GetTaggedTime(tag) );
                        }
                   }
               }
           }
           if( probmin_eta2pi > 0.02 && probmin_eta2pi < 1.0 ){
              if( probmin_3pi < 0.01 ){
                 if( (probmin_etapr > 0.02) && (probmin_etapr < 1.0) ){
                    if(MC_weight)
                        six_fit_IM_eta2pi0_f->FillWeighted( etap_fit_final.M()+mass_shift, MCw );
                    else
                        six_fit_IM_eta2pi0_f->Fill( etap_fit_final.M(), GetTagger()->GetTaggedTime(tag) );
                }
                else{
                    if(MC_weight)
                        six_fit_IM_eta2pi0_g->FillWeighted( etap_fit_final.M()+mass_shift, MCw );
                    else
                        six_fit_IM_eta2pi0_g->Fill( etap_fit_final.M(), GetTagger()->GetTaggedTime(tag) );
                }
              }
           }
           if( (probmin_3pi > 0.01) && (probmin_eta2pi < 0.1) ){ // dir 3pi0
               if(etap_fit.M() > 650.0){
                   if(MC_weight)
                       CB_EnergySum_3pi0->Fill(ESum_MC, MCw);
                   else
                       CB_EnergySum_3pi0->Fill(GetTrigger()->GetEnergySum(),GetTagger()->GetTaggedTime(tag));
               }

               for(Int_t t = 0; t < GetTracks()->GetNTracks(); t++)
                   if(GetTracks()->HasCB(t))
                        time_clusters_CB_3pi0->Fill(  GetTracks()->GetTime(t), GetTracks()->GetCentralCrystal(t) ,GetTagger()->GetTaggedTime(tag));

               six_fit_best_3pi_IM_v_E->Fill(h[0].E(), h[0].M(), GetTagger()->GetTaggedTime(tag));
               six_fit_best_3pi_IM_v_E->Fill(h[1].E(), h[1].M(), GetTagger()->GetTaggedTime(tag));
               six_fit_best_3pi_IM_v_E->Fill(h[2].E(), h[2].M(), GetTagger()->GetTaggedTime(tag));

               if(etap_fit.M() > 650.0){
                   if(MC_weight){
                       six_phy_3pi_IMpipi_v_IMppi->FillWeighted((h[0]+h[1]).M2()/1.0e6, (h[2]+proton_fit).M2()/1.0e6, MCw);
                       six_phy_3pi_IMpipi_v_IMppi->FillWeighted((h[2]+h[0]).M2()/1.0e6, (h[1]+proton_fit).M2()/1.0e6, MCw);
                       six_phy_3pi_IMpipi_v_IMppi->FillWeighted((h[1]+h[2]).M2()/1.0e6, (h[0]+proton_fit).M2()/1.0e6, MCw);
                   }
                   else{
                       six_phy_3pi_IMpipi_v_IMppi->Fill((h[0]+h[1]).M2()/1.0e6, (h[2]+proton_fit).M2()/1.0e6, GetTagger()->GetTaggedTime(tag));
                       six_phy_3pi_IMpipi_v_IMppi->Fill((h[2]+h[0]).M2()/1.0e6, (h[1]+proton_fit).M2()/1.0e6, GetTagger()->GetTaggedTime(tag));
                       six_phy_3pi_IMpipi_v_IMppi->Fill((h[1]+h[2]).M2()/1.0e6, (h[0]+proton_fit).M2()/1.0e6, GetTagger()->GetTaggedTime(tag));
                   }
               }

               for(Int_t isix = 0; isix < 6; isix++){
                    six_fit_best_3pi0_pi_E_v_th->Fill( photons_fit[imin_3pi[isix]].E(), photons_fit[imin_3pi[isix]].Theta()*TMath::RadToDeg(), GetTagger()->GetTaggedTime(tag));

                    int imass = int(isix/2);
                    int nEN =  int(photons_rec[imin_3pi[isix]].E()/20.);
                    int nTH =  75*int(photons_rec[imin_3pi[isix]].Theta()*TMath::RadToDeg()/1.0);
                    int nBIN = nEN + nTH;

                    int nEN2 =  int(photons_rec[imin_3pi[isix]].E()/40.);
                    int nTH2 =  40*int(photons_rec[imin_3pi[isix]].Theta()*TMath::RadToDeg()/2.0);
                    int nBIN2 = nEN2 + nTH2;

                    if(etap_fit.M() > 650.0){
                        six_fit_mgg_v_eth->Fill(nBIN, rc[imass].M(),GetTagger()->GetTaggedTime(tag));
                        six_fit_mgg_v_eth_2->Fill(nBIN2, rc[imass].M(),GetTagger()->GetTaggedTime(tag));
                    }

                    if( GetTracks()->HasCB( set_min[imin_3pi[isix]+2]) ){
                        int idet1 = detnr[imin_3pi[isix]];
                        int idet2 = int(detnr[imin_3pi[isix]])/2.;

                        int inEn1 = int(photons_rec[imin_3pi[isix]].E()/20.);
                        int inEn2 = int(photons_rec[imin_3pi[isix]].E()/40.);

                        int nBIN  = idet1*50 + inEn1;
                        int nBIN2 = idet1*25 + inEn2;
                        int nBIN3 = idet2*50 + inEn1;
                        int nBIN4 = idet2*25 + inEn2;

                        six_fit_mgg_v_CB->Fill(nBIN, rc[imass].M(),GetTagger()->GetTaggedTime(tag));
                        six_fit_mgg_v_CB_2->Fill(nBIN2, rc[imass].M(),GetTagger()->GetTaggedTime(tag));
                        six_fit_mgg_v_CB_3->Fill(nBIN3, rc[imass].M(),GetTagger()->GetTaggedTime(tag));
                        six_fit_mgg_v_CB_4->Fill(nBIN4, rc[imass].M(),GetTagger()->GetTaggedTime(tag));

                        if(etap_fit.M() > 900.0)
                            IMgg_v_det_3pi0_CB->Fill( rc[imass].M(), detnr[imin_3pi[isix]], GetTagger()->GetTaggedTime(tag));
                    }
                    else{

                        int idet1 = detnr[imin_3pi[isix]];
                        int idet2 = int(detnr[imin_3pi[isix]])/2.;

                        int inEn1 = int(photons_rec[imin_3pi[isix]].E()/20.);
                        int inEn2 = int(photons_rec[imin_3pi[isix]].E()/40.);

                        int nBIN  = idet1*50 + inEn1;
                        int nBIN2 = idet1*25 + inEn2;
                        int nBIN3 = idet2*50 + inEn1;
                        int nBIN4 = idet2*25 + inEn2;

                        six_fit_mgg_v_TAPS->Fill(nBIN, rc[imass].M(),GetTagger()->GetTaggedTime(tag));
                        six_fit_mgg_v_TAPS_2->Fill(nBIN2, rc[imass].M(),GetTagger()->GetTaggedTime(tag));
                        six_fit_mgg_v_TAPS_3->Fill(nBIN3, rc[imass].M(),GetTagger()->GetTaggedTime(tag));
                        six_fit_mgg_v_TAPS_4->Fill(nBIN4, rc[imass].M(),GetTagger()->GetTaggedTime(tag));
                        if(etap_fit.M() > 900.0){
                            IMgg_v_det_3pi0_TAPS->Fill( rc[imass].M() , detnr[imin_3pi[isix]], GetTagger()->GetTaggedTime(tag));
                        }
                    }
               }

           if(MC_weight)
                six_fit_IM_3pi->FillWeighted(etap_fit.M(), MCw);
           else
                six_fit_IM_3pi->Fill(etap_fit.M(), GetTagger()->GetTaggedTime(tag));

             } // end of p 3pi0 final analysis

//             if(PDF_cut->IsInside(probmin_eta2pi, probmin_3pi)){
           if( (probmin_3pi < 0.08) && (probmin_eta2pi > 0.02) ){ //eta prime
//               if( (probmin_3pi < 0.1) && (probmin_eta2pi > 0.04) ){ //eta prime


             if(MC_weight)
                CB_EnergySum_etapr->FillWeighted(ESum_MC, MCw );
             else
                CB_EnergySum_etapr->Fill(GetTrigger()->GetEnergySum(), GetTagger()->GetTaggedTime(tag));

             if(MC_weight)
               proton_fit_e_v_th_final->FillWeighted(proton_fit.E() - MASS_PROTON, proton_fit.Theta()*TMath::RadToDeg(), MCw);
             else
               proton_fit_e_v_th_final->Fill(proton_fit.E() - MASS_PROTON, proton_fit.Theta()*TMath::RadToDeg(), GetTagger()->GetTaggedTime(tag));

               for(UInt_t t = 0; t < 6; t++){
                   if( MC_weight )
                       six_fit_EvTh_g_final->FillWeighted(photons_fit_final[t].E(), photons_fit_final[t].Theta()*TMath::RadToDeg(), MCw);
                   else
                       six_fit_EvTh_g_final->Fill(photons_fit_final[t].E(), photons_fit_final[t].Theta()*TMath::RadToDeg(),GetTagger()->GetTaggedTime(tag));
               }

               if(MC_weight){
                    six_fit_PDF_eta2pi_v_Meta2pi->FillWeighted(probmin_eta2pi, etap_fit_final.M(), MCw);
                    six_fit_best_eta_IM_v_E->FillWeighted(g[0].E(), g[0].M(), MCw );
                    six_fit_best_2pi_IM_v_E->FillWeighted(g[1].E(), g[1].M(), MCw);
                    six_fit_best_2pi_IM_v_E->FillWeighted(g[2].E(), g[2].M(), MCw);

                    six_phy_etapr_v_BeamE->FillWeighted((GetTagger()->GetVector(set_min[0])).E(), etap_fit.M(), MCw);
                    six_phy_etapr_eta2pi_v_BeamE->FillWeighted((GetTagger()->GetVector(set_min[0])).E(), etap_fit_final.M(), MCw);
                    six_phy_etapr_v_EPT->FillWeighted((GetTagger()->GetTaggedChannel(set_min[0])), etap_fit.M(), MCw);
                    six_phy_etapr_eta2pi_v_EPT->FillWeighted((GetTagger()->GetTaggedChannel(set_min[0])), etap_fit_final.M(), MCw);

                    six_fit_best_eta->FillWeighted( g[0].M(), MCw );
                    six_fit_best_2pi->FillWeighted( g[1].M(), MCw );
                    six_fit_best_2pi->FillWeighted( g[2].M(), MCw );

                    six_fit_best_eta_rec->FillWeighted( rc_sig[0].M(), MCw );
                    six_fit_best_2pi_rec->FillWeighted( rc_sig[1].M(), MCw );
                    six_fit_best_2pi_rec->FillWeighted( rc_sig[2].M(), MCw );

                    int ncl = 7;
                    if(eight_clusters)
                        ncl = 8 ;
                    six_fit_IM_eta2pi_v_ncl->Fill(ncl, etap_fit_final.M(), GetTagger()->GetTaggedTime(tag) );

                    six_fit_IM_eta2pi->FillWeighted( etap_fit_final.M()+mass_shift, MCw );
                    six_fit_etaprfinal_pdf->FillWeighted( probmin_eta2pi, MCw );
                    six_fit_eta_PDF_v_Metapr->FillWeighted( probmin_eta2pi, etap_fit_final.M(), MCw);

                    Double_t ztrue = etapr_6gTrue.GetTrueVertex().Z();
                    true_z_after_final_fit->FillWeighted(ztrue, MCw);
               }
               else{
                   six_fit_PDF_eta2pi_v_Meta2pi->Fill(probmin_eta2pi, etap_fit_final.M(),GetTagger()->GetTaggedTime(tag));
                   six_fit_best_eta_IM_v_E->Fill(g[0].E(), g[0].M(), GetTagger()->GetTaggedTime(tag));
                   six_fit_best_2pi_IM_v_E->Fill(g[1].E(), g[1].M(), GetTagger()->GetTaggedTime(tag));
                   six_fit_best_2pi_IM_v_E->Fill(g[2].E(), g[2].M(), GetTagger()->GetTaggedTime(tag));

                   six_phy_etapr_v_BeamE->Fill((GetTagger()->GetVector(set_min[0])).E(), etap_fit.M(), GetTagger()->GetTaggedTime(tag));
                   six_phy_etapr_eta2pi_v_BeamE->Fill((GetTagger()->GetVector(set_min[0])).E(), etap_fit_final.M(), GetTagger()->GetTaggedTime(tag));
                   six_phy_etapr_v_EPT->Fill((GetTagger()->GetTaggedChannel(set_min[0])), etap_fit.M(), GetTagger()->GetTaggedTime(tag));
                   six_phy_etapr_eta2pi_v_EPT->Fill((GetTagger()->GetTaggedChannel(set_min[0])), etap_fit_final.M(), GetTagger()->GetTaggedTime(tag));

                   six_fit_best_eta->Fill( (g[0]).M(), GetTagger()->GetTaggedTime(tag) );
                   six_fit_best_2pi->Fill( g[1].M(), GetTagger()->GetTaggedTime(tag) );
                   six_fit_best_2pi->Fill( g[2].M(), GetTagger()->GetTaggedTime(tag) );

                   six_fit_best_eta_rec->Fill( (rc_sig[0]).M(), GetTagger()->GetTaggedTime(tag)  );
                   six_fit_best_2pi_rec->Fill( rc_sig[1].M(), GetTagger()->GetTaggedTime(tag)  );
                   six_fit_best_2pi_rec->Fill( rc_sig[2].M(), GetTagger()->GetTaggedTime(tag)  );

                   int ncl = 7;
                   if(eight_clusters)
                       ncl = 8 ;


                   six_fit_IM_eta2pi_v_ncl->Fill(ncl, etap_fit_final.M(), GetTagger()->GetTaggedTime(tag) );

                   six_fit_IM_eta2pi->Fill( etap_fit_final.M(), GetTagger()->GetTaggedTime(tag) );

                   six_fit_etaprfinal_pdf->Fill( probmin_eta2pi, GetTagger()->GetTaggedTime(tag) );
                   six_fit_eta_PDF_v_Metapr->Fill( probmin_eta2pi, etap_fit_final.M(),GetTagger()->GetTaggedTime(tag));
               }

               for(int isix = 0; isix < 6; isix++){
                   if(isix < 2)
                       six_fit_best_etapr_eta_E_v_th->Fill( photons_fit[imin_eta2pi[isix]].E(), photons_fit[imin_eta2pi[isix]].Theta()*TMath::RadToDeg(), GetTagger()->GetTaggedTime(tag));
                   else
                       six_fit_best_etapr_pi_E_v_th->Fill( photons_fit[imin_eta2pi[isix]].E(), photons_fit[imin_eta2pi[isix]].Theta()*TMath::RadToDeg(), GetTagger()->GetTaggedTime(tag));
               }

               // Filling of final results

               if(MC_weight){
                   six_phy_DP_X->FillWeighted(Xfit, etap_fit_final.M()+mass_shift, MCw );
                   six_phy_DP_Y->FillWeighted(Yfit, etap_fit_final.M()+mass_shift, MCw );
                   six_phy_DP_020->FillWeighted(DP_binnr_fit020, etap_fit_final.M()+mass_shift, MCw );
                   six_phy_DP_015->FillWeighted(DP_binnr_fit015, etap_fit_final.M()+mass_shift, MCw );
                   six_phy_DP_010->FillWeighted(DP_binnr_fit010, etap_fit_final.M()+mass_shift, MCw );
                   six_phy_DP_075->FillWeighted(DP_binnr_fit075, etap_fit_final.M()+mass_shift, MCw );
                   six_phy_DP_005->FillWeighted(DP_binnr_fit005, etap_fit_final.M()+mass_shift, MCw );

                   six_phy_M_pi1pi2_v_etapr->FillWeighted(m_2pi0_fit / 1.0e3, etap_fit_final.M()+mass_shift, MCw);
                   six_phy_M_etapi_v_etapr->FillWeighted(m_etapi01_fit / 1.0e3, etap_fit_final.M()+mass_shift, MCw);
                   six_phy_M_etapi_v_etapr->FillWeighted(m_etapi02_fit / 1.0e3, etap_fit_final.M()+mass_shift, MCw);

                   six_phy_M_pi1pi2_v_etapr2->FillWeighted(TMath::Sqrt(m_2pi0_fit), etap_fit_final.M()+mass_shift, MCw);
                   six_phy_M_etapi_v_etapr2->FillWeighted(TMath::Sqrt(m_etapi01_fit), etap_fit_final.M()+mass_shift, MCw);
                   six_phy_M_etapi_v_etapr2->FillWeighted(TMath::Sqrt(m_etapi02_fit), etap_fit_final.M()+mass_shift, MCw);

                   if(etap_fit_final.M()> 930. && etap_fit_final.M() < 980.){
                       true_six_phy_dMpipi_v_Mpipi->Fill(m_2pi0_fit / 1.0e3, (m_2pi0_fit/1.0e3-m_2pi0True*1.0e3), MCw);
                       true_six_phy_DP_020->Fill(DPnrTrue020, DP_binnr_fit020, MCw);
                       true_six_phy_DP_015->Fill(DPnrTrue015, DP_binnr_fit015, MCw);
                       true_six_phy_DP_010->Fill(DPnrTrue010, DP_binnr_fit010, MCw);
                       true_six_phy_DP_075->Fill(DPnrTrue075, DP_binnr_fit075, MCw);
                       true_six_phy_DP_005->Fill(DPnrTrue005, DP_binnr_fit005, MCw);

                       true_six_phy_dX_v_DPbin->Fill(DP_binnr_fit015, Xfit - Xtrue, MCw);
                       true_six_phy_dY_v_DPbin->Fill(DP_binnr_fit015, Yfit - Ytrue, MCw);

                       true_six_phy_dX_v_X->Fill(Xfit, Xfit - Xtrue, MCw);
                       true_six_phy_dY_v_Y->Fill(Yfit, Yfit - Ytrue, MCw);
                   }
               }
               else{
                   six_phy_M_pi1pi2_v_etapr->Fill(m_2pi0_fit / 1.0e3, etap_fit_final.M(), GetTagger()->GetTaggedTime(tag));
                   six_phy_M_etapi_v_etapr->Fill(m_etapi01_fit / 1.0e3, etap_fit_final.M(), GetTagger()->GetTaggedTime(tag));
                   six_phy_M_etapi_v_etapr->Fill(m_etapi02_fit / 1.0e3, etap_fit_final.M(), GetTagger()->GetTaggedTime(tag));

                   six_phy_M_pi1pi2_v_etapr2->Fill(TMath::Sqrt(m_2pi0_fit), etap_fit_final.M(), GetTagger()->GetTaggedTime(tag));
                   six_phy_M_etapi_v_etapr2->Fill(TMath::Sqrt(m_etapi01_fit), etap_fit_final.M(), GetTagger()->GetTaggedTime(tag));
                   six_phy_M_etapi_v_etapr2->Fill(TMath::Sqrt(m_etapi02_fit), etap_fit_final.M(), GetTagger()->GetTaggedTime(tag));

                   six_phy_DP_X->Fill(Xfit, etap_fit_final.M(), GetTagger()->GetTaggedTime(tag) );
                   six_phy_DP_Y->Fill(Yfit, etap_fit_final.M(), GetTagger()->GetTaggedTime(tag) );
                   six_phy_DP_020->Fill(DP_binnr_fit020, etap_fit_final.M(), GetTagger()->GetTaggedTime(tag) );
                   six_phy_DP_015->Fill(DP_binnr_fit015, etap_fit_final.M(), GetTagger()->GetTaggedTime(tag) );
                   six_phy_DP_010->Fill(DP_binnr_fit010, etap_fit_final.M(), GetTagger()->GetTaggedTime(tag) );
                   six_phy_DP_075->Fill(DP_binnr_fit075, etap_fit_final.M(), GetTagger()->GetTaggedTime(tag) );
                   six_phy_DP_005->Fill(DP_binnr_fit005, etap_fit_final.M(), GetTagger()->GetTaggedTime(tag) );

               }

               if(MC_weight)
                    six_phy_etapr_prod_diff_distr->FillWeighted(diffbin, etap_fit_final.M()+mass_shift, MCw);
               else
                    six_phy_etapr_prod_diff_distr->Fill(diffbin, etap_fit_final.M(), GetTagger()->GetTaggedTime(tag));

               bool final_cond = false;
               if(probmin_eta2pi > 0.04 && probmin_eta2pi < 1.0)
                   if( probmin_3pi < 0.0075 )
                    if((probmin_etapr > 0.04) && (probmin_etapr < 1.0))
                        final_cond = true;

               if(final_cond){
                  if(MC_weight){
                       six_phy_etapr_prod_diff_distr_metapr->FillWeighted( diffbin_pr, MCw );
                       six_phy_M_pi1pi2_v_etapr_fit->FillWeighted(m_2pi0_fit_pr / 1.0e3, MCw );
                       six_phy_M_etapi_v_etapr_fit->FillWeighted(m_etapi01_fit_pr / 1.0e3, MCw );
                       six_phy_M_etapi_v_etapr_fit->FillWeighted(m_etapi02_fit_pr / 1.0e3, MCw );
                       six_phy_M_pi1pi2_v_etapr_fit2->FillWeighted(TMath::Sqrt(m_2pi0_fit_pr), MCw );
                       six_phy_M_etapi_v_etapr_fit2->FillWeighted(TMath::Sqrt(m_etapi01_fit_pr), MCw );
                       six_phy_M_etapi_v_etapr_fit2->FillWeighted(TMath::Sqrt(m_etapi02_fit_pr), MCw );
                       true_six_phy_dMpipi_v_Mpipi_metapr->Fill(m_2pi0_fit_pr / 1.0e3, (m_2pi0_fit_pr/1.0e3-m_2pi0True*1.0e3), MCw );

                       six_phy_DP_X_pr->FillWeighted(Xfit_pr, MCw );
                       six_phy_DP_Y_pr->FillWeighted(Yfit_pr, MCw );
                       six_phy_DP_020_pr->FillWeighted(DP_binnr_fit020_pr, MCw );
                       six_phy_DP_015_pr->FillWeighted(DP_binnr_fit015_pr, MCw );
                       six_phy_DP_010_pr->FillWeighted(DP_binnr_fit010_pr, MCw );
                       six_phy_DP_075_pr->FillWeighted(DP_binnr_fit075_pr, MCw );
                       six_phy_DP_005_pr->FillWeighted(DP_binnr_fit005_pr, MCw );

                       true_six_phy_DP_020_pr->Fill(DPnrTrue020, DP_binnr_fit020_pr, MCw);
                       true_six_phy_DP_015_pr->Fill(DPnrTrue015, DP_binnr_fit015_pr, MCw);
                       true_six_phy_DP_010_pr->Fill(DPnrTrue010, DP_binnr_fit010_pr, MCw);
                       true_six_phy_DP_075_pr->Fill(DPnrTrue075, DP_binnr_fit075_pr, MCw);
                       true_six_phy_DP_005_pr->Fill(DPnrTrue005, DP_binnr_fit005_pr, MCw);

                       true_six_phy_dX_v_DPbin_metapr->Fill(DP_binnr_fit015, Xfit_pr - Xtrue, MCw);
                       true_six_phy_dY_v_DPbin_metapr->Fill(DP_binnr_fit015, Yfit_pr - Ytrue, MCw);

                       true_six_phy_dX_v_X_metapr->Fill(Xfit_pr, Xfit_pr - Xtrue, MCw);
                       true_six_phy_dY_v_Y_metapr->Fill(Yfit_pr, Yfit_pr - Ytrue, MCw);

                   }
                   else{
                       six_phy_etapr_prod_diff_distr_metapr->Fill( diffbin_pr, GetTagger()->GetTaggedTime(tag) );
                       six_phy_M_pi1pi2_v_etapr_fit->Fill(m_2pi0_fit_pr / 1.0e3, GetTagger()->GetTaggedTime(tag) );
                       six_phy_M_etapi_v_etapr_fit->Fill(m_etapi01_fit_pr / 1.0e3, GetTagger()->GetTaggedTime(tag) );
                       six_phy_M_etapi_v_etapr_fit->Fill(m_etapi02_fit_pr / 1.0e3, GetTagger()->GetTaggedTime(tag) );
                       six_phy_M_pi1pi2_v_etapr_fit2->Fill(TMath::Sqrt(m_2pi0_fit_pr), GetTagger()->GetTaggedTime(tag) );
                       six_phy_M_etapi_v_etapr_fit2->Fill(TMath::Sqrt(m_etapi01_fit_pr), GetTagger()->GetTaggedTime(tag) );
                       six_phy_M_etapi_v_etapr_fit2->Fill(TMath::Sqrt(m_etapi02_fit_pr), GetTagger()->GetTaggedTime(tag) );

                       six_phy_DP_X_pr->Fill(Xfit_pr, GetTagger()->GetTaggedTime(tag) );
                       six_phy_DP_Y_pr->Fill(Yfit_pr, GetTagger()->GetTaggedTime(tag) );

                       six_phy_DP_020_pr->Fill(DP_binnr_fit020_pr, GetTagger()->GetTaggedTime(tag) );
                       six_phy_DP_015_pr->Fill(DP_binnr_fit015_pr, GetTagger()->GetTaggedTime(tag) );
                       six_phy_DP_010_pr->Fill(DP_binnr_fit010_pr, GetTagger()->GetTaggedTime(tag) );
                       six_phy_DP_075_pr->Fill(DP_binnr_fit075_pr, GetTagger()->GetTaggedTime(tag) );
                       six_phy_DP_005_pr->Fill(DP_binnr_fit005_pr, GetTagger()->GetTaggedTime(tag) );
                   }
               }
           }
        }
    }

//    if(vec_weight.size() > 0 )
//        FillTree();

}


void AdlarsonPhysics::GetBest6gCombination(Double_t& sigma_eta, Double_t& sigma_pi0, Double_t& chi2min_eta2pi, Double_t& chi2min_3pi, std::vector<int>& imin_eta2pi, std::vector<int>& imin_3pi, std::vector<comb>& etatwopi_comb, std::vector<comb>& threepi_comb)
{
    imin_eta2pi.resize(6);
    imin_3pi.resize(6);
    Double_t imgg[3];
    Double_t chi2_eta2pi[3];
    Double_t chi2_3pi;
    Double_t sgm_eta, sgm_pi0, sgm_3pi0;
    chi2min_eta2pi = 1.0e6;
    chi2min_3pi = 1.0e6;

    std::vector<Int_t> sig_perm, bkgd_perm;
    etatwopi_comb.resize(0);
    threepi_comb.resize(0);

    for( Int_t i = 0; i < 15 ; i++){
       for (Int_t j = 0; j < 3; j++){
//            imgg[j] = ( photons_fit[perm6g[i][0 + j*2]] + photons_fit[perm6g[i][1 + j*2]]).M();
            imgg[j] = ( photons_rec[perm6g[i][0 + j*2]] + photons_rec[perm6g[i][1 + j*2]]).M();

            sgm_eta = 40.0;
            sgm_pi0 = 14.0;
            chi2_eta2pi[j] = 0;
       }

       chi2_eta2pi[0] = TMath::Power((imgg[0] - MASS_ETA )/(sgm_eta), 2.0) + TMath::Power((imgg[1] - MASS_PI0 )/(sgm_pi0), 2) + TMath::Power((imgg[2] - MASS_PI0 )/(sgm_pi0), 2.0);
       chi2_eta2pi[1] = TMath::Power((imgg[1] - MASS_ETA )/(sgm_eta), 2.0) + TMath::Power((imgg[2] - MASS_PI0 )/(sgm_pi0), 2) + TMath::Power((imgg[0] - MASS_PI0 )/(sgm_pi0), 2.0);
       chi2_eta2pi[2] = TMath::Power((imgg[2] - MASS_ETA )/(sgm_eta), 2.0) + TMath::Power((imgg[0] - MASS_PI0 )/(sgm_pi0), 2) + TMath::Power((imgg[1] - MASS_PI0 )/(sgm_pi0), 2.0);

       for( Int_t t = 0; t < 3; t++ ){
           sig_perm.resize(0);
           if(t == 0){
               if(chi2_eta2pi[0] < chi2min_eta2pi){
                   imin_eta2pi.resize(0);
                   chi2min_eta2pi = chi2_eta2pi[0];
                   for(Int_t u = 0; u < 6; u++)
                        imin_eta2pi.push_back(perm6g[i][u]);
               }
               for(Int_t u = 0; u < 6; u++)
                    sig_perm.push_back(perm6g[i][u]);

               etatwopi_comb.push_back(comb(sig_perm, chi2_eta2pi[0]));

           }
           else if( t == 1){               
               if(chi2_eta2pi[1] < chi2min_eta2pi)
               {
                   imin_eta2pi.resize(0);
                   chi2min_eta2pi = chi2_eta2pi[1];
                   imin_eta2pi.push_back(perm6g[i][2]);
                   imin_eta2pi.push_back(perm6g[i][3]);
                   imin_eta2pi.push_back(perm6g[i][0]);
                   imin_eta2pi.push_back(perm6g[i][1]);
                   imin_eta2pi.push_back(perm6g[i][4]);
                   imin_eta2pi.push_back(perm6g[i][5]);

               }
               sig_perm.push_back(perm6g[i][2]);
               sig_perm.push_back(perm6g[i][3]);
               sig_perm.push_back(perm6g[i][0]);
               sig_perm.push_back(perm6g[i][1]);
               sig_perm.push_back(perm6g[i][4]);
               sig_perm.push_back(perm6g[i][5]);

               etatwopi_comb.push_back(comb(sig_perm, chi2_eta2pi[1]));
           }
           else{
               if(chi2_eta2pi[2] < chi2min_eta2pi)
               {
                   imin_eta2pi.resize(0);
                   chi2min_eta2pi = chi2_eta2pi[2];
                   imin_eta2pi.push_back(perm6g[i][4]);
                   imin_eta2pi.push_back(perm6g[i][5]);
                   imin_eta2pi.push_back(perm6g[i][0]);
                   imin_eta2pi.push_back(perm6g[i][1]);
                   imin_eta2pi.push_back(perm6g[i][2]);
                   imin_eta2pi.push_back(perm6g[i][3]);
               }
               sig_perm.push_back(perm6g[i][4]);
               sig_perm.push_back(perm6g[i][5]);
               sig_perm.push_back(perm6g[i][0]);
               sig_perm.push_back(perm6g[i][1]);
               sig_perm.push_back(perm6g[i][2]);
               sig_perm.push_back(perm6g[i][3]);

               etatwopi_comb.push_back(comb(sig_perm, chi2_eta2pi[2]));
           }
        }

        bkgd_perm.resize(0);
        sgm_3pi0 = 14.0;
        chi2_3pi = TMath::Power( (imgg[0] - MASS_PI0 )/(sgm_3pi0), 2) + TMath::Power( (imgg[1] - MASS_PI0 )/(sgm_3pi0), 2) + TMath::Power( (imgg[2] - MASS_PI0 )/(sgm_3pi0), 2);

        if(chi2_3pi < chi2min_3pi)
        {
            imin_3pi.resize(0);
            chi2min_3pi = chi2_3pi;
            for(Int_t u = 0; u < 6; u++)
                 imin_3pi.push_back(perm6g[i][u]);
        }
        for(Int_t u = 0; u < 6; u++)
            bkgd_perm.push_back(perm6g[i][u]);

        threepi_comb.push_back(comb(bkgd_perm,chi2_3pi));

    }

    auto vec_it = etatwopi_comb.cbegin();
    while( vec_it != etatwopi_comb.cend()){
        if(vec_it->second < chi2min_eta2pi){
            imin_eta2pi.resize(0);
            chi2min_eta2pi = vec_it->second;
            imin_eta2pi = vec_it->first;
        }
        ++vec_it;
    }

    std::sort(etatwopi_comb.begin(), etatwopi_comb.end(), [](const comb left, const comb right){ return left.second < right.second; });

    auto vec_it_3pi = threepi_comb.cbegin();
    while( vec_it_3pi != threepi_comb.cend()){
        if(vec_it_3pi->second < chi2min_3pi ){
            imin_3pi.resize(0);
            chi2min_3pi = vec_it_3pi->second;
            imin_3pi = vec_it_3pi->first;
        }
        ++vec_it_3pi;
    }

    std::sort(threepi_comb.begin(), threepi_comb.end(), [](const comb left, const comb right){ return left.second < right.second; });
    return;
}


void AdlarsonPhysics::test_correct_hypothesis(Double_t& prob_etapr, Double_t& prob_eta2pi, Double_t& prob_3pi, std::vector<Int_t>& set_min, std::vector<int>& imin_eta2pi, std::vector<int>& imin_3pi, std::vector<comb>& imap_eta2pi, std::vector<comb>& imap_3pi)
{
    Double_t prob_3pi_fit       = -1.0e6;
    Double_t prob_eta2pi_fit    = -1.0e6;
    Double_t chi2_3pi_fit       =  1.0e6;
    Double_t chi2_eta2pi_fit    =  1.0e6;

    int niter_eta2pi = 0;
    int niter_3pi = 0;

    Int_t n_photons_min;
    bool etapr_candidate = false;

    Int_t iteration_place_etatwopi = -1;
    Int_t iteration_place_threepi = -1;
    imin_eta2pi.resize(0);
    imin_3pi.resize(0);

    std::vector<int> imin_eta2pi_temp;
    std::vector<int> imin_3pi_temp;
    photons_fit_final.resize(0);
    photons_fit_final_metapr.resize(0);

    for( Int_t icand = 0; icand < 3; icand++){

        Is_CB.resize(0);

        beam_3pi.SetFromVector( GetTagger()->GetVector(set_min[0]));
        obs.resize(0);
        unc.resize(0);
        unc = Get_unc_R(3,0,obs);
        beam_3pi.Smear(unc, 2);
        beam_3pi.APLCONSettings();

        Double_t test = (GetTagger()->GetVector(set_min[0]) - IM6g_vec).E();
        // GetTracks()->GetVector(set_min[1]).E()

        proton_3pi.SetFromValues( test ,R_TAPS[GetTracks()->GetCentralCrystal(set_min[1])],GetTracks()->GetVector(set_min[1]).Phi() );
        obs.resize(0);
        obs.push_back(GetTracks()->GetCentralCrystal(set_min[1]));
        unc.resize(0);
        unc = Get_unc_R( 2, 2, obs);
        proton_3pi.Smear_R(unc, 1);
        Is_CB.push_back(0);

        n_photons_min = 0;
        imin_3pi_temp.resize(0);
        //        imin_3pi_temp = imin_3pi;
        imin_3pi_temp = imap_3pi.at(icand).first;
        Photons_six_3pi.resize(0);
        for ( Int_t jgam_min = 0; jgam_min < 6 ; jgam_min++ ){
            int tr3pi = set_min[ imin_3pi_temp[jgam_min] + 2];
            if( GetTracks()->HasCB(tr3pi)){
                Photons_six_3pi[n_photons_min].SetFromValues( GetTracks()->GetVector(tr3pi).E(),GetTracks()->GetVector(tr3pi).Theta(),GetTracks()->GetVector(tr3pi).Phi() );
                Is_CB.push_back(1);
                obs.resize(0);
                obs.push_back(Photons_six_3pi[n_photons_min].Ek);
                obs.push_back(Photons_six_3pi[n_photons_min].Theta);                
                unc.resize(0);
                unc = Get_unc_R( 1, 1, obs);
                Photons_six_3pi[n_photons_min].Smear_R(unc ,0);
            }
            else{ // ( GetTracks()->HasTAPS(jgam_min) )
                Is_CB.push_back(0);
                Photons_six_3pi[n_photons_min].SetFromValues( GetTracks()->GetVector(tr3pi).E(), R_TAPS[GetTracks()->GetCentralCrystal(tr3pi)],GetTracks()->GetVector(tr3pi).Phi() );
                obs.resize(0);
                obs.push_back(Photons_six_3pi[n_photons_min].Ek);
                obs.push_back(Photons_six_3pi[n_photons_min].Theta);
                obs.push_back(GetTracks()->GetCentralCrystal(tr3pi));
                unc.resize(0);
                unc = Get_unc_R( 2, 1, obs);
                Photons_six_3pi[n_photons_min].Smear_R(unc ,0);
            }
            n_photons_min++;
        }

            const APLCON::Result_t& result_3pi = kinfit3pi.DoFit();
            if(result_3pi.Status == APLCON::Result_Status_t::Success){
                if( (result_3pi.ChiSquare < chi2_3pi_fit) && (result_3pi.Probability < 1.0)  && (result_3pi.Probability > 1.0e-5)  ){
                    chi2_3pi_fit = result_3pi.ChiSquare;
                    prob_3pi_fit = result_3pi.Probability;
                    imin_3pi.resize(0);
                    imin_3pi = imin_3pi_temp;
                    iteration_place_threepi = icand;

                    niter_3pi = result_3pi.NIterations;
                }
            }

        Is_CB.resize(0);

        beam_eta2pi.SetFromVector( GetTagger()->GetVector( set_min[0]) );
        obs.resize(0);
        unc.resize(0);
        unc = Get_unc_R(3,0,obs);
        beam_eta2pi.Smear(unc, 2);
        beam_eta2pi.APLCONSettings();

        proton_eta2pi.SetFromValues( test ,R_TAPS[GetTracks()->GetCentralCrystal(set_min[1])],GetTracks()->GetVector(set_min[1]).Phi() );
        obs.resize(0);
        obs.push_back(GetTracks()->GetCentralCrystal(set_min[1]));
        unc.resize(0);
        unc = Get_unc_R( 2, 2, obs);
        proton_eta2pi.Smear_R(unc, 1);
        Is_CB.push_back(0);

        n_photons_min = 0;
        imin_eta2pi_temp.resize(0);
        Photons_six_eta2pi.resize(0);
        imin_eta2pi_temp = imap_eta2pi.at(icand).first;

//        imin_eta2pi_temp = imin_eta2pi;

        for ( Int_t jgam_min = 0; jgam_min < 6 ; jgam_min++ ){
            int treta2pi = set_min[ imin_eta2pi_temp[jgam_min] + 2 ];
            if( GetTracks()->HasCB( treta2pi ) )
            {
                Is_CB.push_back(1);
                Photons_six_eta2pi[n_photons_min].SetFromValues( GetTracks()->GetVector(treta2pi).E(),GetTracks()->GetVector(treta2pi).Theta(), GetTracks()->GetVector(treta2pi).Phi() );
                obs.resize(0);
                obs.push_back(Photons_six_eta2pi[n_photons_min].Ek);
                obs.push_back(Photons_six_eta2pi[n_photons_min].Theta);
                unc.resize(0);
                unc = Get_unc_R( 1, 1, obs);
                Photons_six_eta2pi[n_photons_min].Smear_R(unc, 0);
            }
            else // ( GetTracks()->HasTAPS(treta2pi) )
            {
                Is_CB.push_back(0);
                Photons_six_eta2pi[n_photons_min].SetFromValues( GetTracks()->GetVector(treta2pi).E(), R_TAPS[GetTracks()->GetCentralCrystal(treta2pi)], GetTracks()->GetVector(treta2pi).Phi() );
                obs.resize(0);
                obs.push_back(Photons_six_eta2pi[n_photons_min].Ek);
                obs.push_back(Photons_six_eta2pi[n_photons_min].Theta);
                obs.push_back(GetTracks()->GetCentralCrystal(treta2pi));
                unc.resize(0);
                unc = Get_unc_R( 2, 1, obs);
                Photons_six_eta2pi[n_photons_min].Smear_R(unc ,0);
            }
            n_photons_min++;
        }

        const APLCON::Result_t& result_eta2pi = kinfiteta2pi.DoFit();
        if(result_eta2pi.Status == APLCON::Result_Status_t::Success){
            if((result_eta2pi.ChiSquare < chi2_eta2pi_fit) && (result_eta2pi.Probability < 1.0) && (result_eta2pi.Probability > 1.0e-5)){
                chi2_eta2pi_fit = result_eta2pi.ChiSquare;
                prob_eta2pi_fit = result_eta2pi.Probability;

                niter_eta2pi = result_eta2pi.NIterations;

                imin_eta2pi.resize(0);
                imin_eta2pi = imin_eta2pi_temp;

                iteration_place_etatwopi = icand;
                photons_fit_final.resize(0);

                for ( Int_t n_photons_min = 0; n_photons_min < 6 ; n_photons_min++ ){
                    std::vector<double> obs_kfit;
                    obs_kfit.resize(0);
                    TLorentzVector eta2pi(0.0, 0.0, 0.0, 0.0);

                    string s = Form("Photon_eta2pi%0i[0]", n_photons_min+1); // Energy
                    string t = Form("Photon_eta2pi%0i[1]", n_photons_min+1); // th
                    string u = Form("Photon_eta2pi%0i[2]", n_photons_min+1); // phi
                    string z = Form("v_z"); // z

                    double E    = result_eta2pi.Variables.at(s).Value.After;
                    double th   = result_eta2pi.Variables.at(t).Value.After;
                    double ph   = result_eta2pi.Variables.at(u).Value.After;
                    double vx_z = result_eta2pi.Variables.at(z).Value.After;

                    if(Is_CB[n_photons_min+1] == 1){
                        obs_kfit = {E, E, th, ph};
                        eta2pi = GetLVCorrForZ(obs_kfit, vx_z, Is_CB[n_photons_min+1], 0.0);
                        photons_fit_final.push_back(eta2pi);
                    }
                    else{
                        obs_kfit = {E, E, th, ph}; // th is R in TAPS
                        eta2pi = GetLVCorrForZ(obs_kfit, vx_z, Is_CB[n_photons_min+1], 0.0);
                        photons_fit_final.push_back(eta2pi);
                    }
                }

//                photons_fit_final.resize(0);
//                for(UInt_t igam_fit = 0; igam_fit < Photons_six.size(); igam_fit++)
//                    photons_fit_final.push_back( FitParticle::Make( Photons_six_eta2pi[igam_fit], 0.0) );

            }
        }
    }

    if((prob_eta2pi_fit > 0.0) && (prob_eta2pi_fit < 1.0)){
        prob_eta2pi = prob_eta2pi_fit;
        NIeta2pi0->Fill(niter_eta2pi);
        NIeta2pi0vPDF->Fill(niter_eta2pi, prob_eta2pi);


        if(prob_eta2pi > 0.01)
            etapr_candidate = true;
    }
    else
        prob_eta2pi = 0.0;

    if((prob_3pi_fit > 0.0) && (prob_3pi_fit < 1.0)){
        prob_3pi = prob_3pi_fit;
        NI3pi0vPDF->Fill(niter_3pi, prob_3pi);
    }
    else
        prob_3pi = 0.0;

    // now choose the best found 3pi cand and eta2pi cand and get the final result for these

    six_fit_which_place_best_etapr_cand->Fill(iteration_place_etatwopi);
    six_fit_which_place_best_3pi_cand->Fill(iteration_place_threepi);

    if(etapr_candidate){

        Is_CB.resize(0);
        beam_final.SetFromVector( GetTagger()->GetVector(set_min[0]) );
        obs.resize(0);
        unc.resize(0);
        unc = Get_unc_R(3,0,obs);
        beam_final.Smear(unc, 2);
        beam_final.APLCONSettings();


        Double_t test = (GetTagger()->GetVector(set_min[0]) - IM6g_vec).E();
        proton_final.SetFromValues( test , R_TAPS[GetTracks()->GetCentralCrystal(set_min[1])],GetTracks()->GetVector(set_min[1]).Phi() );
        obs.resize(0);
        obs.push_back(GetTracks()->GetCentralCrystal(set_min[1]));
        unc.resize(0);
        unc = Get_unc_R( 2, 2, obs);
        proton_final.Smear_R(unc, 1);
        Is_CB.push_back(0);

        Photons_six_final.resize(0);
        Int_t n_photons_min = 0;
        for ( Int_t jgam_min = 0; jgam_min < 6 ; jgam_min++ ){
            int tretapr = set_min[ imin_eta2pi[jgam_min] + 2 ];
            if( GetTracks()->HasCB( tretapr ) )
            {
                Is_CB.push_back(1);
                Photons_six_final[n_photons_min].SetFromValues( GetTracks()->GetVector(tretapr).E(),GetTracks()->GetVector(tretapr).Theta(), GetTracks()->GetVector(tretapr).Phi() );
                obs.resize(0);
                obs.push_back(Photons_six_final[n_photons_min].Ek);
                obs.push_back(Photons_six_final[n_photons_min].Theta);
                unc.resize(0);
                unc = Get_unc_R( 1, 1, obs);
                Photons_six_final[n_photons_min].Smear_R(unc, 0);
            }
            else
            {
                Is_CB.push_back(0);
                Photons_six_final[n_photons_min].SetFromValues( GetTracks()->GetVector(tretapr).E(), R_TAPS[GetTracks()->GetCentralCrystal(tretapr)], GetTracks()->GetVector(tretapr).Phi() );
                obs.resize(0);
                obs.push_back(Photons_six_final[n_photons_min].Ek);
                obs.push_back(Photons_six_final[n_photons_min].Theta);
                obs.push_back(GetTracks()->GetCentralCrystal(tretapr));
                unc.resize(0);
                unc = Get_unc_R( 2, 1, obs);
                Photons_six_final[n_photons_min].Smear_R(unc ,0);
            }
            n_photons_min++;
        }
    }

    const APLCON::Result_t& result_final = kinfit_final.DoFit();
    if(result_final.Status == APLCON::Result_Status_t::Success){
        if( (result_final.Probability < 1.0 ) && ( result_final.Probability > 0.01 ) )
        {
            NItetapr->Fill(result_final.NIterations);
            prob_etapr = result_final.Probability;

            NItetaprvPDF->Fill(result_final.NIterations, prob_etapr);
            photons_fit_final_metapr.resize(0);
            for ( Int_t n_photons_min = 0; n_photons_min < 6 ; n_photons_min++ ){
                std::vector<double> obs_kfit;
                obs_kfit.resize(0);
                TLorentzVector etapr(0.0, 0.0, 0.0, 0.0);

                string s = Form("Photons_six_final%0i[0]", n_photons_min+1); // Energy
                string t = Form("Photons_six_final%0i[1]", n_photons_min+1); // th
                string u = Form("Photons_six_final%0i[2]", n_photons_min+1); // phi
                string z = Form("v_z"); // z

                double E    = result_final.Variables.at(s).Value.After;
                double th   = result_final.Variables.at(t).Value.After;
                double ph   = result_final.Variables.at(u).Value.After;
                double vx_z = result_final.Variables.at(z).Value.After;

                if(Is_CB[n_photons_min+1] == 1){
                    obs_kfit = {E, E, th, ph};
                    etapr = GetLVCorrForZ(obs_kfit, vx_z, Is_CB[n_photons_min+1], 0.0);
                    photons_fit_final_metapr.push_back(etapr);
                }
                else{
                    obs_kfit = {E, E, th, ph}; // th is R in TAPS
                    etapr = GetLVCorrForZ(obs_kfit, vx_z, Is_CB[n_photons_min+1], 0.0);
                    photons_fit_final_metapr.push_back(etapr);
                }
            }

//            for(UInt_t igam_fit = 0; igam_fit < Photons_six.size(); igam_fit++)
//                photons_fit_final_metapr.push_back( FitParticle::Make( Photons_six_final[igam_fit], 0.0) );

        }
    }

    return;
}

void AdlarsonPhysics::FillTree(){
//    if(vec_weight.size() >= MAX_LENGTH){
//        std::cout << "too many tagger hits!" <<  std::endl;
//        return;
//    }

//    branch_length = vec_weight.size();
//    for(int i = 0; i < branch_length ; i++){
//        branch_weight[i]        = vec_weight[i];
//        branch_nclusters[i]     = vec_nclusters[i];
//        branch_taggerenergy[i]  = vec_taggerenergy[i];
//        branch_fitted_pr_th[i]  = vec_fitted_pr_th[i];
//        branch_fitted_pr_e[i]   = vec_fitted_pr_e[i];
//        branch_fitted_zvx[i]    = fitted_zvx[i];
//        branch_pdf_eta2pi[i]    = vec_pdf_eta2pi[i];
//        branch_pdf_3pi[i]       = vec_pdf_3pi[i];
//        branch_pdf_etapr[i]     = vec_pdf_etapr[i];
//        branch_costh_epr_cm[i]  = vec_costh_epr_cm[i];
//        branch_X[i]             = vec_X[i];
//        branch_Y[i]             = vec_Y[i];
//        branch_DP_005[i]        = vec_DP_005[i];
//        branch_DP_010[i]        = vec_DP_010[i];
//        branch_DP_015[i]        = vec_DP_015[i];
//        branch_M_eta2pi[i]      = vec_M_eta2pi[i];
//        branch_M_etapi1[i]      = vec_M_etapi1[i];
//        branch_M_etapi2[i]      = vec_M_etapi2[i];
//        branch_M_pipi[i]        = vec_M_pipi[i];
//        branch_costh_epr_cm_pr[i] = vec_costh_epr_cm_pr[i];
//        branch_X_pr[i]          = vec_X_pr[i];
//        branch_Y_pr[i]          = vec_Y_pr[i];
//        branch_DP_005_pr[i]     = vec_DP_005_pr[i];
//        branch_DP_010_pr[i]     = vec_DP_010_pr[i];
//        branch_DP_015_pr[i]     = vec_DP_015_pr[i];
//        branch_M_etapi_pr1[i]   = vec_M_etapi_pr1[i];
//        branch_M_etapi_pr2[i]   = vec_M_etapi_pr2[i];
//        branch_M_pipi_pr[i]     = vec_M_pipi_pr[i];
//    }

//    tree->Fill();

}

void AdlarsonPhysics::tengAnalysis(UInt_t ipr)
{
    for(Int_t tag = 0; tag < GetTagger()->GetNTagged(); tag++)
    {
        if(TMath::Abs(GetTagger()->GetTaggedTime(tag))> 40.) continue;

        sigma_eta = 20; sigma_pi0 = 12;
        Is_CB.resize(0);
        photons_rec_ten.resize(0);
        photons_fit_10g.resize(0);
        IM10g_vec.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
        TLorentzVector etap10g_fit(0.0, 0.0, 0.0, 0.0);

        std::vector<Int_t> set, set_min;     // particle conf: tagger, proton, photon 1, 2, 3, 4, 5, 6, 7, 8, 9, 10;
        set.resize(0);
        set_min.resize(0);

        for(UInt_t i = 0; i < nPhotons_ten; i++)
        {
            photons_rec_ten[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
            photons_fit_10g[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
        }

        Double_t chi2_min = 1.0e6;
        Double_t prob_min = 1.0e6;

        set.push_back(tag);

        beam10g.SetFromVector( GetTagger()->GetVector(tag) );
        obs.resize(0);
        unc.resize(0);
        unc = Get_unc_R(3,0,obs);
        beam10g.Smear(unc, 2);
        beam10g.APLCONSettings();

        set.push_back(ipr);
        proton10g.SetFromVector( GetTracks()->GetVector(ipr) );
        obs.push_back(GetTracks()->GetCentralCrystal(ipr));
        unc.resize(0);
        unc = Get_unc( 2, 2, obs);
        proton10g.Smear(unc, 1);

        UInt_t n_photons = 0;
        for ( UInt_t jgam = 0; jgam < ClustersInTime.size() ; jgam++ ){
            UInt_t kgam = ClustersInTime[jgam];
            if( kgam != ipr ) // the id proton cluster
            {
                set.push_back(kgam);
                photons_rec_ten.push_back(GetTracks()->GetVector(kgam));
                Photons_ten[n_photons].SetFromVector( GetTracks()->GetVector(kgam) );
                IM10g_vec += GetTracks()->GetVector(kgam);

                if( GetTracks()->HasCB(kgam) )
                {
                    Is_CB.push_back(1);
                    obs.resize(0);
                    unc.resize(0);
                    obs.push_back(Photons_ten[n_photons].Ek);
                    obs.push_back(Photons_ten[n_photons].Theta);
                    unc = Get_unc( 1, 1, obs);
                    Photons_ten[n_photons].Smear(unc, 0);
                }
                else // ( GetTracks()->HasTAPS(kgam) )
                {
                    Is_CB.push_back(0);
                    obs.resize(0);
                    obs.push_back(Photons_ten[n_photons].Ek);
                    obs.push_back(Photons_ten[n_photons].Theta);
                    obs.push_back(GetTracks()->GetCentralCrystal(kgam));
                    unc.resize(0);
                    unc = Get_unc( 2, 1, obs);
                    Photons_ten[n_photons].Smear(unc, 0);
                }
                n_photons++;
            }
        }

        ten_rec_IM->Fill( IM10g_vec.M(),GetTagger()->GetTaggedTime(tag) );

        const APLCON::Result_t& result10g = kinfit10g.DoFit();
        if(result10g.Status == APLCON::Result_Status_t::Success)
        {
            if( result10g.ChiSquare <  chi2_min )
            {
                chi2_min = result10g.ChiSquare;
                prob_min = result10g.Probability;
                set_min = set;
            }
        }
        // Here run kinfit with the best combination:

        if( (prob_min < 0.01) || ( TMath::Abs(prob_min - 1.0e6) < 1.0e-4) ) continue;

        kfit_chi2_10g->Fill( result10g.ChiSquare, GetTagger()->GetTaggedTime(tag) );
        kfit_pdf_10g->Fill( result10g.Probability, GetTagger()->GetTaggedTime(tag) );

        Int_t inr = 0;
           for(const auto& it_map : result10g.Variables) {
                const APLCON::Result_Variable_t& var = it_map.second;
                kfit_Pulls_10g->Fill(var.Pull, inr, GetTagger()->GetTaggedTime(tag));
                inr++;
           }

        ten_rec_IM_v_MMp->Fill( MMp , IM10g_vec.M(),GetTagger()->GetTaggedTime(tag));

        std::vector<double> obs_kfit;
        obs_kfit.resize(0);
        photons_fit_10g.resize(0);
        IM6g_vec.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
        etap10g_fit.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
        TLorentzVector etap_fit_tmp(0.0, 0.0, 0.0, 0.0);

        double E_min = 1000.;
        int imin = -10;
        if(result10g.Status == APLCON::Result_Status_t::Success){
            for(UInt_t igam_fit = 0; igam_fit < Photons_ten.size(); igam_fit++){
                photons_fit[igam_fit] = FitParticle::Make(Photons_ten[igam_fit], 0);
                if(photons_fit[igam_fit].E() < E_min){
                    E_min = photons_fit[igam_fit].E();
                    imin = igam_fit;
                }

            }
        }
        ten_fit_EvTh_g_min->FillWeighted(photons_fit[imin].E(), photons_fit[imin].Theta()*TMath::RadToDeg(), MCw);

        if(result10g.Status == APLCON::Result_Status_t::Success){
            for(UInt_t igam_fit = 0; igam_fit < Photons_ten.size(); igam_fit++){
                photons_fit_10g[igam_fit] = FitParticle::Make(Photons_ten[igam_fit], 0);
                ten_fit_EvTh_g->Fill(photons_fit_10g[igam_fit].E(), photons_fit_10g[igam_fit].Theta()*TMath::RadToDeg(), GetTagger()->GetTaggedTime(tag) );
                etap10g_fit += photons_fit_10g[igam_fit];
            }
        }

        Int_t diffbin =  diff_distr( GetTagger()->GetTaggedEnergy(tag), etap10g_fit );

        if(MC){
            IM10g_fit->FillWeighted(etap10g_fit.M(), MCw );
            ten_phy_etapr_prod_diff_distr->FillWeighted(diffbin, etap10g_fit.M(), MCw);
        }
        else{
            IM10g_fit->Fill(etap10g_fit.M(), GetTagger()->GetTaggedTime(tag));
            ten_phy_etapr_prod_diff_distr->Fill(diffbin ,etap10g_fit.M(), GetTagger()->GetTaggedTime(tag));
        }

        return;

        Double_t chi2_eta2pi0 = 1.0e6;
        std::vector<int> imin_eta3pi2pi;
        imin_eta3pi2pi.resize(0);

        // find best 10g combination here
        GetBest10gCombination(sigma_eta, sigma_pi0, chi2_eta2pi0, imin_eta3pi2pi);


        // here fit eta and 2pi0 best combination(s?)

        beam10g_eta2pi.SetFromVector( GetTagger()->GetVector(tag) );
        obs.resize(0);
        unc.resize(0);
        unc = Get_unc(3,0,obs);
        beam10g_eta2pi.Smear(unc, 2);
        beam10g_eta2pi.APLCONSettings();

        set.push_back(ipr);
        proton10g_eta2pi.SetFromVector( GetTracks()->GetVector(ipr) );
        obs.resize(0);
//        obs.push_back(proton10g_eta2pi.Ek);
//        obs.push_back(proton10g_eta2pi.Theta);
        obs.push_back(GetTracks()->GetCentralCrystal(ipr));
        unc.resize(0);
        unc = Get_unc( 2, 2, obs);
        proton10g_eta2pi.Smear(unc, 1);

        for ( UInt_t jgam = 0; jgam < 10 ; jgam++ ){
            Photons_ten_eta2pi[jgam].SetFromVector( GetTracks()->GetVector( set_min[imin_eta3pi2pi[jgam] + 2]));
            if( GetTracks()->HasCB(set_min[imin_eta3pi2pi[jgam] + 2]) )
            {
                obs.resize(0);
                obs.push_back(Photons_ten_eta2pi[jgam].Ek);
                obs.push_back(Photons_ten_eta2pi[jgam].Theta);
                unc.resize(0);
                unc = Get_unc( 1, 1, obs);
                Photons_ten_eta2pi[jgam].Smear(unc, 0);
            }
            else // in TAPS
            {
                obs.resize(0);
                obs.push_back(Photons_ten_eta2pi[jgam].Ek);
                obs.push_back(Photons_ten_eta2pi[jgam].Theta);
                unc.resize(0);
                unc = Get_unc( 2, 1, obs);
                Photons_ten_eta2pi[jgam].Smear(unc, 0);
            }
        }
        const APLCON::Result_t& result10g_eta2pi = kinfit10g_eta2pi.DoFit();
        if(result10g_eta2pi.Status == APLCON::Result_Status_t::Success)
        {
            if(result10g_eta2pi.Probability > 0.01 && result10g_eta2pi.Probability < 1.0){

                ten_fit_PDF_eta2pi->Fill(result10g_eta2pi.Probability,GetTagger()->GetTaggedTime(tag));

                photons_fit.resize(0);
                for(UInt_t igam_fit = 0; igam_fit < Photons_ten.size(); igam_fit++)
                {
                    photons_fit[igam_fit] = FitParticle::Make(Photons_ten_eta2pi[igam_fit], 0);
                }

                TLorentzVector fin[3];
                for(int i6g = 0; i6g < 6; i6g++)
                    fin[0] += photons_fit[i6g];

                fin[1] = photons_fit[6] + photons_fit[7];
                fin[2] = photons_fit[8] + photons_fit[9];

                IM10g_fit_best_cand->Fill((fin[0]+fin[1]+fin[2]).M(),GetTagger()->GetTaggedTime(tag));

                Double_t Xfit = -2.0;
                Double_t Yfit = -2.0;
                Int_t DP_binnr_fit020 = -100;
                bw = 0.2;
                DalitzPlot(fin, Xfit, Yfit, bw, DP_binnr_fit020);

                true_ten_phy_dX_v_DPbin->Fill(DP_binnr_fit020, Xfit - Xtrue);
                true_ten_phy_dY_v_DPbin->Fill(DP_binnr_fit020, Yfit - Ytrue);

                ten_fit_X_v_pdf_eta2pi->Fill(result10g_eta2pi.Probability, Xfit - Xtrue);
                ten_fit_Y_v_pdf_eta2pi->Fill(result10g_eta2pi.Probability, Yfit - Ytrue);
            }
        }

    }
}

void AdlarsonPhysics::GetBest10gCombination( Double_t& sigma_eta, Double_t& sigma_pi0, Double_t& chi2min_eta3pi, std::vector<int>& imin_eta3pi2pi )
{

    // here also use sgm as fcn of theta_g2 and theta_g1. In addition use chi2 test for eta -->3pi0 -->6g
    TLorentzVector  eta_6g_cand;
    Double_t        im6g;
    Double_t        chi2;
    Double_t        chi2_min = 1.0e6;

    std::vector<Int_t> sixg;
    std::vector<Int_t> sixg_min;

    std::vector<Int_t> teng;
    std::vector<Int_t> teng_min;
    teng_comb.resize(0);

    sixg.resize(0);
    sixg_min.resize(0);
    for(Int_t iperm = 0; iperm < 210; iperm++){
        eta_6g_cand.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
        std::vector<int> pi0c;
        im6g = 0;
        Double_t im2g_1, im2g_2;

        sixg.resize(0);
        for(int jperm = 0; jperm < 6; jperm++){
            eta_6g_cand += ( photons_rec_ten[perm6outof10g[iperm][jperm]] );
            sixg.push_back( perm6outof10g[iperm][jperm] );
        }


        // here calculate the best 3pi0 combination
        double chi2_3pi0_min = 1.0e6;
        double chi2_3pi0     = 1.0e6;
        for(int isix = 0; isix < 15; isix++)
        {
             Double_t mgg_1=0., mgg_2=0., mgg_3 = 0.;

             mgg_1 = (photons_rec_ten[ sixg[ perm6g[isix][0] ] ] + photons_rec_ten[ sixg[ perm6g[isix][1] ] ]).M();
             mgg_2 = (photons_rec_ten[ sixg[ perm6g[isix][2] ] ] + photons_rec_ten[ sixg[ perm6g[isix][3] ] ]).M();
             mgg_3 = (photons_rec_ten[ sixg[ perm6g[isix][4] ] ] + photons_rec_ten[ sixg[ perm6g[isix][5] ] ]).M();

             chi2_3pi0 = TMath::Power( (mgg_1 - MASS_PI0 )/sigma_pi0 ,2) + TMath::Power( (mgg_2 - MASS_PI0 )/sigma_pi0 ,2) + TMath::Power( (mgg_3 - MASS_PI0 )/sigma_pi0 ,2);
             if(chi2_3pi0 < chi2_3pi0_min){
                 chi2_3pi0_min = chi2_3pi0;
                 sixg_min.resize(0);
                 for( auto i = 0; i < 6; i++ )
                     sixg_min.push_back( sixg[ perm6g[isix][i] ] );
             }
        }

        for(int ipi = 0; ipi < 10; ipi++){
            bool perm_available = true;
            for(Int_t kperm = 0; kperm < 6; kperm++){
                if(perm6outof10g[iperm][kperm] == ipi)
                    perm_available = false;
            }
            if(perm_available)
                pi0c.push_back(ipi);
        }

        im6g = eta_6g_cand.M();

        for(int i2pi0 = 0; i2pi0 < 3; i2pi0++){
            teng.resize(0);
            for(Int_t igam = 0; igam < 6; igam++){
//                teng.push_back(perm6outof10g[iperm][igam]);
                teng.push_back(sixg_min[igam]);
            }

            im6g = 0;
            im2g_1 = 0;
            im2g_2 = 0;
            if(i2pi0 == 0){
                for(int &r : pi0c)
                    teng.push_back(r);
            }
            else if(i2pi0 == 1){
                teng.push_back(pi0c[0]);
                teng.push_back(pi0c[2]);
                teng.push_back(pi0c[1]);
                teng.push_back(pi0c[3]);
            }
            else{ //(i2pi0 == 2)
                teng.push_back(pi0c[0]);
                teng.push_back(pi0c[3]);
                teng.push_back(pi0c[1]);
                teng.push_back(pi0c[2]);
            }
            im6g = eta_6g_cand.M();
            im2g_1 = ( photons_rec_ten[teng[6]] + photons_rec_ten[teng[7]] ).M();
            im2g_2 = ( photons_rec_ten[teng[8]] + photons_rec_ten[teng[9]] ).M();

            chi2 = chi2_3pi0_min + TMath::Power( (im6g - MASS_ETA )/sigma_eta ,2) + TMath::Power( (im2g_1 - MASS_PI0 )/sigma_pi0 ,2) + TMath::Power( (im2g_2 - MASS_PI0 )/sigma_pi0 ,2);
            teng_comb.push_back(comb(teng,chi2));
        }
    }

        auto vec_it = teng_comb.cbegin();
        while( vec_it != teng_comb.cend() ){
            if(vec_it->second < chi2_min){
                teng_min.resize(0);
                chi2_min = vec_it->second;
                teng_min = vec_it->first;
            }
            ++vec_it;
        }

        std::sort(teng_comb.begin(), teng_comb.end(), [](const comb left, const comb right){ return left.second < right.second; });

        chi2min_eta3pi = chi2_min;
        imin_eta3pi2pi = teng_min;

        return;

}

Int_t AdlarsonPhysics::diff_distr(const Double_t beam_e, TLorentzVector &fin ){

    Int_t diffbin = 0;
    TLorentzVector sqrt_s(0.,0., beam_e, beam_e + MASS_PROTON);
    TVector3 cm_vector  = -(sqrt_s).BoostVector();
    TLorentzVector etapr_cm  = fin;
    etapr_cm.Boost(cm_vector);
    Double_t angle = etapr_cm.Theta();
    Double_t x = TMath::Cos( angle );
    double width;

    Int_t ee = -100;
    if( beam_e < (Legendre[0]+3.25) )
         ee = 0;
    else if( beam_e >= (Legendre[66]-6.5) )
         ee = 11;
    else{

        if(beam_e < 1473.1)
            width = 3.25;
        else
            width = 6.5;

        bool energy_region = false;
        int k = 6;
        while(! energy_region && k < 61){
            if( ( beam_e >= (Legendre[k]-width) ) && ( beam_e < (Legendre[k+6]+width) ) ){
                 ee = k/6;
                 energy_region = true;
            }
            else
                k += 6;
        }
    }
    return diffbin = ee*20 + int((1. + x)/0.2);
}

void AdlarsonPhysics::DalitzPlot(const TLorentzVector g[3] , Double_t &X, Double_t &Y , Double_t bw, Int_t &DP_nr)
{
    // g is the array of TLorentzVectors containing the final event sample in order
    // 0 1 2 - eta pi01 pi02
    // Dalitz plot variables

    Double_t T_eta, T_pi1, T_pi2;
    Double_t Q;
    Double_t X1, X2;
    Double_t Y_max = 1.5;

    TVector3 etaprime_rest      = -(g[0]+g[1]+g[2]).BoostVector();
    TLorentzVector eta_ep_rest  = g[0];
    TLorentzVector pi01_ep_rest = g[1];
    TLorentzVector pi02_ep_rest = g[2];

    eta_ep_rest.Boost(etaprime_rest);
    pi01_ep_rest.Boost(etaprime_rest);
    pi02_ep_rest.Boost(etaprime_rest);

    T_eta = eta_ep_rest.E() - eta_ep_rest.M();
    T_pi1 = pi01_ep_rest.E() - pi01_ep_rest.M();
    T_pi2 = pi02_ep_rest.E() - pi02_ep_rest.M();

    Q = T_eta + T_pi1 + T_pi2;
    X1 = ( TMath::Sqrt(3.0) / Q ) * ( T_pi1 - T_pi2 );
    X2 = ( TMath::Sqrt(3.0) / Q ) * ( T_pi2 - T_pi1 );

    if(X1 > X2)
        X = X1;
    else
        X = X2;

    Y = ( MASS_ETA + 2 * MASS_PI0 ) / MASS_PI0 * ( T_eta / Q ) - 1.0;

    if(TMath::Abs(bw-0.2) < 1.0e-3)
        DP_nr = int ( X / bw ) + 10*int( (Y  + Y_max ) / bw );
    if(TMath::Abs(bw-0.15) < 1.0e-3)
        DP_nr = int ( X / bw ) + 15*int( (Y  + Y_max ) / bw );
    if(TMath::Abs(bw-0.10) < 1.0e-3)
        DP_nr = int ( X / bw ) + 20*int( (Y  + Y_max ) / bw );
    if(TMath::Abs(bw-0.075) < 1.0e-3)
        DP_nr = int ( X / bw ) + 30*int( (Y  + Y_max ) / bw );
    if(TMath::Abs(bw-0.05) < 1.0e-3)
        DP_nr = int ( X / bw ) + 30*int( (Y  + Y_max ) / bw );

    return;

}

void AdlarsonPhysics::m2pi0_metapi0(  TLorentzVector g[3], Double_t &m_etapi01, Double_t &m_etapi02, Double_t &m_2pi0 )
{
    // g is the array of TLorentzVectors containing the final event sample in order
    // 0 1 2 - eta pi01 pi02
    // Dalitz plot variables


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

std::vector<double> AdlarsonPhysics::Get_unc(Int_t apparatus_nr, Int_t particle, std::vector<double>& obs)
{
    double Ek_s, Theta_s, Phi_s;
    std::vector<double> unc;
    unc.resize(0);

    if( particle == 1 )         // gamma
    {
        double Ek = obs[0];
        double theta = obs[1]*TMath::RadToDeg();

        Ek_s    = g_e_vz->GetBinContent( g_e_vz->FindBin(Ek,theta) );
        Theta_s = g_th_vz->GetBinContent( g_th_vz->FindBin(Ek,theta) );
        Phi_s   = g_phi_vz->GetBinContent( g_phi_vz->FindBin(Ek,theta) );

        if(TMath::Abs(Ek_s) < 1.0e-5){
            if(apparatus_nr == 1 )
                Ek_s = 0.06; // CB
            else
                Ek_s = 0.09; // TAPS
        }
        if(TMath::Abs(Theta_s) < 1.0e-5){
            if(apparatus_nr == 1 )
                Theta_s = 3.9; // CB
            else
                Theta_s = 1.0; // TAPS
        }
        if(TMath::Abs(Phi_s) < 1.0e-5){
            if(apparatus_nr == 1 )
                Phi_s = 3.1; // CB
            else
                Phi_s = 2.5; // TAPS
        }
    }
    else if( particle == 2 ){    // proton
        Int_t inr = obs[0];
        Ek_s = 0.0;
        Theta_s = p_TAPS_th_vz->GetBinContent(p_TAPS_th_vz->FindBin( inr ));
        Phi_s = p_TAPS_fi_vz->GetBinContent(p_TAPS_fi_vz->FindBin( inr ));


        if(TMath::Abs(Theta_s) < 1.0e-5)
            Theta_s = 0.7;
        if(TMath::Abs(Phi_s) < 1.0e-5)
            Phi_s = 4.0;

    }
    else // beam energy photon
    {
        Ek_s = 2.0/TMath::Sqrt(3.0);
        Theta_s = 1e-3;
        Phi_s = 1e-3;
    }

    unc.push_back(Ek_s);
    unc.push_back(Theta_s);
    unc.push_back(Phi_s);

    return unc;
}

std::vector<double> AdlarsonPhysics::Get_unc_R(Int_t apparatus_nr, Int_t particle, std::vector<double>& obs)
{
    double Ek_s, Theta_s, Phi_s;
//    double e_g = 1.0;
    std::vector<double> unc;
    unc.resize(0);
    if( particle == 1 )         // gamma
    {
        double Ek = obs[0];
        double theta = obs[1]*TMath::RadToDeg();

        Ek_s    = g_e->GetBinContent( g_e->FindBin(Ek,theta) );

        Phi_s   = g_phi->GetBinContent( g_phi->FindBin(Ek,theta) );
        Phi_s   *= TMath::DegToRad();

        if(apparatus_nr == 2){ // 1 = CB, 2 = TAPS
            Int_t inr = obs[2]; // for photon, Ek, theta and Central Crystal element is in obs[0], obs[1] and obs[2]
            Theta_s = g_R->GetBinContent( g_R->FindBin(inr) );
        }
//        if(apparatus_nr == 2){ // 1 = CB, 2 = TAPS
//            Theta_s = g_th->GetBinContent( g_th->FindBin(Ek,theta) );
//            Theta_s *= TMath::DegToRad();
//        }
        else{ // in CB
            Theta_s = g_th->GetBinContent( g_th->FindBin(Ek,theta) );
            Theta_s *= TMath::DegToRad();
        }

        if(TMath::Abs(Ek_s) < 1.0e-5){
            if(apparatus_nr == 1 )
                Ek_s = 0.06; // CB
            else
                Ek_s = 0.09; // TAPS
        }
        if(TMath::Abs(Theta_s) < 1.0e-5){
            if(apparatus_nr == 1 )
                Theta_s = 3.9*TMath::DegToRad(); // CB
            else{
//                Theta_s = 2.0; // TAPS - R [cm] default
                Theta_s = 2.0; // TAPS - th [cm] default
            }
        }
        if(TMath::Abs(Phi_s) < 1.0e-5){
            if(apparatus_nr == 1 )
                Phi_s = 3.1*TMath::DegToRad(); // CB
            else
                Phi_s = 2.5*TMath::DegToRad(); // TAPS
        }
    }
    else if( particle == 2 ){    // proton
        Int_t inr = obs[0];
        Ek_s = 0.0;
                Theta_s = p_TAPS_R->GetBinContent(p_TAPS_R->FindBin( inr )); // sigmaR in TAPS
//                Theta_s = p_TAPS_th->GetBinContent(p_TAPS_th->FindBin( inr )); // sigma theta in TAPS
//                Theta_s *= TMath::DegToRad();
                Phi_s = p_TAPS_fi->GetBinContent(p_TAPS_fi->FindBin( inr ));
                Phi_s *= TMath::DegToRad();

        if(TMath::Abs(Theta_s) < 1.0e-5){
            Theta_s = 2.0;  // TAPS - R [cm] default
 //           Theta_s = 2.0*TMath::DegToRad();  // TAPS - theta [rad] default
        }
        if(TMath::Abs(Phi_s) < 1.0e-5)
            Phi_s = 4.0*TMath::DegToRad();

    }
    else // beam energy photon
    {
        Ek_s = 2.0/TMath::Sqrt(3.0);
        Theta_s = 1e-3;
        Phi_s = 1e-3;
    }

    unc.push_back(Ek_s);
    unc.push_back(Theta_s);
    unc.push_back(Phi_s);

    return unc;
}

void AdlarsonPhysics::FitParticle::Smear(std::vector<double> unc, int particle) {
    if( particle == 0 ) // gamma
    {
        Ek_Sigma = unc[0]*Ek;
        Theta_Sigma = unc[1]*TMath::DegToRad();
        Phi_Sigma = unc[2]*TMath::DegToRad();
    }
    else if(particle == 1) //proton
    {
        Ek_Sigma = 0;
        Theta_Sigma = unc[1]*TMath::DegToRad();
        Phi_Sigma = unc[2]*TMath::DegToRad();
    }
    else // beam
    {
        Ek_Sigma = unc[0];
        Theta_Sigma = unc[1];
        Phi_Sigma = unc[2];
    }
}

void AdlarsonPhysics::FitParticle::Smear_R(std::vector<double> unc, int particle) {
    if( particle == 0 ) // gamma
    {
        Ek_Sigma = unc[0]*Ek;
        Theta_Sigma = unc[1];
        Phi_Sigma = unc[2];
    }
    else if(particle == 1) //proton
    {
        Ek_Sigma = 0;
        Theta_Sigma = unc[1];
        Phi_Sigma = unc[2];
    }
    else // beam
    {
        Ek_Sigma = unc[0];
        Theta_Sigma = unc[1];
        Phi_Sigma = unc[2];
    }
}

void AdlarsonPhysics::FitParticle::APLCONSettings()
{
    Ek_Setting.StepSize = 0;
    Theta_Setting.StepSize = 0;
    Phi_Setting.StepSize = 0;
}

Int_t AdlarsonPhysics::perm4g[9][4]=
{
      {0, 1,  2, 3}, {2, 3,  0, 1}, {0, 2,  1, 3}, {1, 3,  0, 2}, {0, 3,  1, 2}, {1, 2,  0, 3},
      {0, 1,  2, 3}, {0, 2,  1, 3}, {0, 3,  1, 2}
};

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

void AdlarsonPhysics::Kinfit_test()
{
    TLorentzVector beam_true, proton_true;
    TLorentzVector photons_true[6];

    beam_true.SetPxPyPzE(0., 0., etapr_6gTrue.GetTrueBeamEnergy()*1000, etapr_6gTrue.GetTrueBeamEnergy()*1000);
//    beam.SetFromVector( beam_true );
//    beam.Smear(3, true);

    proton_true.SetPxPyPzE(etapr_6gTrue.GetTrueProtonLV().Px()*1000, etapr_6gTrue.GetTrueProtonLV().Py()*1000, etapr_6gTrue.GetTrueProtonLV().Pz()*1000,etapr_6gTrue.GetTrueProtonLV().E()*1000);
    proton.SetFromVector( proton_true );
//    proton.Smear(2, true);

    for(UInt_t i = 0; i < etapr_6gTrue.GetNgamma(); i++ )
    {
        (photons_true[i]).SetPxPyPzE( etapr_6gTrue.GetTrueGammaLV(i).Px()*1000, etapr_6gTrue.GetTrueGammaLV(i).Py()*1000,etapr_6gTrue.GetTrueGammaLV(i).Pz()*1000, etapr_6gTrue.GetTrueGammaLV(i).E()*1000 );
        Photons_six[i].SetFromVector( photons_true[i] );
//        Photons_six[i].Smear(1, true);
    }

    TLorentzVector etap_rec;
    for(size_t t=0;t<Photons_six.size();t++)
        etap_rec += FitParticle::Make(Photons_six[t], 0);
    test_six_rec_IM->Fill(etap_rec.M());

    const APLCON::Result_t& result = kinfit.DoFit();
    if(result.Status == APLCON::Result_Status_t::Success)
    {
//        cout << result << endl;

       test_six_fit_chi2->Fill(result.ChiSquare);
       test_six_fit_pdf->Fill(result.Probability);

       Int_t inr = 0;
       for(const auto& it_map : result.Variables)
       {
           //const string& varname = it_map.first;
           const APLCON::Result_Variable_t& var = it_map.second;
           test_six_fit_Pulls->Fill(var.Pull, inr);
           inr++;
       }
    }
};

void::AdlarsonPhysics::Tagger_corr(){
    Double_t new_e;
    for ( Int_t tag = 0; tag < GetTagger()->GetNTagged(); tag++){
        new_e = GetTagger()->GetTaggedEnergy(tag)-1.0;
        GetTagger()->SetTaggedEnergy( tag, new_e );
    }
}

void AdlarsonPhysics::Time_corr(){

    if(TMath::Abs(CBtime_corr[0]) < 1.0e-5 ) return;

    for (int i = 0; i < GetTracks()->GetNTracks() ; i++)
    {
        if( GetTracks()->HasCB(i) ){
            tracks->SetTime(i, GetTracks()->GetTime(i) - CBtime_corr[GetTracks()->GetCentralCrystal(i)]);
        }
    }
}

void AdlarsonPhysics::Energy_corr()
{
    Double_t Erec, Ec_temp, DeltaE, Ec, Ec2;
    Double_t smear;

//    int group;
    for (int i = 0; i < GetTracks()->GetNTracks() ; i++){
        if( GetTracks()->HasCB(i) ){
            Erec = GetTracks()->GetVector(i).E();
            Ec_temp = CBgain[GetTracks()->GetCentralCrystal(i)]*Erec;

            DeltaE = Ec_temp*(Double_t)EvdetCB->GetBinContent(EvdetCB->FindBin(Ec_temp, GetTracks()->GetTheta(i)));
            Ec = Ec_temp - DeltaE;

            tracks->SetClusterEnergy(i, Ec);

            if(!MC){
                double gain = GetGain(Ec, GetTracks()->GetCentralCrystal(i));
                Ec2 = Ec*gain;
                tracks->SetClusterEnergy(i, Ec2);
            }

            if(MC){
                smear = CB_unc->GetBinContent(CB_unc->FindBin(GetTracks()->GetCentralCrystal(i)));
                if(smear > 0.0)
                    Ec = pRandoms->Gaus(Ec, smear*Ec);
            }

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

void    AdlarsonPhysics::theta_corr()
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

void    AdlarsonPhysics::CB_TAPS_boundary(){
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
    return;
}

void    AdlarsonPhysics::RandomTime(){
    Double_t Trigger_t, CB_t, TAPS_t;
    Double_t sgm_CB;

    std::vector<Double_t> time_jump;
    Double_t thr;
    int CB_group;
    time_jump.resize(0);
    for(int it = 0; it < 90; it++)
    {
        if((MC_time_jump.find(it)->second[0]  +1.) < 1.0e-4)
            time_jump.push_back(0.0);
        else{
            thr = pRandoms->Uniform(0,1.);
            if(thr < MC_time_jump.find(it)->second[0] )
                time_jump.push_back(-25.0);
            else if(thr >= MC_time_jump.find(it)->second[0] && thr < MC_time_jump.find(it)->second[1])
                time_jump.push_back(0.);
            else if(thr >= MC_time_jump.find(it)->second[1] && thr < MC_time_jump.find(it)->second[2])
                time_jump.push_back(25.);
            else if(thr >= MC_time_jump.find(it)->second[2] && thr < MC_time_jump.find(it)->second[3])
                time_jump.push_back(50.);
            else
                time_jump.push_back(100.);
        }

    }

    Trigger_t   = pRandoms->Gaus(0, 2.2);
    CB_t        = Trigger_t;
    TAPS_t      = pRandoms->Gaus(0, 0.61);

    for( int j = 0; j < GetTracks()->GetNTracks() ; j++ )
    {
        if(GetTracks()->HasCB(j)){
            sgm_CB = TMath::Sqrt( CBtime_sgm[GetTracks()->GetCentralCrystal(j)]*CBtime_sgm[GetTracks()->GetCentralCrystal(j)] - 2.2*2.2 );
            CB_t += pRandoms->Gaus(0, sgm_CB);

            // time to add the offset
            CB_group = int((GetTracks()->GetCentralCrystal(j))/8);

            if( GetTracks()->GetCentralCrystal(j) == 719 )
                tracks->SetTime(j, 0. + CB_t + time_jump[CB_group]);
            else
                tracks->SetTime(j, GetTracks()->GetTime(j) + CB_t + time_jump[CB_group]);
        }
        if(GetTracks()->HasTAPS(j))
            tracks->SetTime(j,GetTracks()->GetTime(j) + TAPS_t);
    }
}

Double_t AdlarsonPhysics::Get_ESumMC(Double_t& ESum){
    Double_t ESum_shift = 17.11;
    for( int j = 0; j < GetTracks()->GetNTracks() ; j++ ){
        if(MCJuly14){
            if( GetTracks()->HasCB(j) ){
                if((GetTracks()->GetCentralCrystal(j) >= 352) && (GetTracks()->GetCentralCrystal(j) <= 415)) // July beam time only
                    ESum -= GetTracks()->GetClusterEnergy(j);
            }
        }
    }
    return ESum + ESum_shift;
}

TLorentzVector    AdlarsonPhysics::GetLVCorrForZ(std::vector<double> EkPThPhi, const double v_z, Int_t& idet, double mass)
{

    double X0       = 2.588;
    double X0_TAPS  = 2.026;
    double Ec       = 13.3;
    double Ec_TAPS  = 13.7;
    constexpr double R_CB = 25.4;
    constexpr double Z_TAPS = 145.7;
    double E, P, th, ph;
    double Px, Py, Pz;
    double R;
    TLorentzVector LVmod(0.0, 0.0, 0.0, 0.0);

    if(idet == 1){ // CB
        E = EkPThPhi[0];
        P = EkPThPhi[1];

        R  = R_CB + X0*std::log2(E/Ec)/std::pow(std::sin(EkPThPhi[2]), 3.0);
        th = TMath::ACos(( R*TMath::Cos(EkPThPhi[2]) - v_z)/ R );
        ph = EkPThPhi[3];

        Px = P*TMath::Sin(th)*TMath::Cos(ph);
        Py = P*TMath::Sin(th)*TMath::Sin(ph);
        Pz = P*TMath::Cos(th);

        LVmod.SetPxPyPzE(Px,Py,Pz,E);
    }
    else{   // TAPS
        E = EkPThPhi[0];
        P = EkPThPhi[1];

        if(mass > 900.)
            R  =    X0_TAPS*std::log2((E-MASS_PROTON)/Ec_TAPS);
        else
            R  =    X0_TAPS*std::log2(E/Ec_TAPS);
        th = TMath::ATan( EkPThPhi[2] / (Z_TAPS - v_z + R));

        ph = EkPThPhi[3];

        Px = P*TMath::Sin(th)*TMath::Cos(ph);
        Py = P*TMath::Sin(th)*TMath::Sin(ph);
        Pz = P*TMath::Cos(th);

        LVmod.SetPxPyPzE(Px,Py,Pz,E);
    }
    return LVmod;
}

void AdlarsonPhysics::TrueAnalysis_etapr6g(TString s){

    Double_t N, a, b, d;

    true_BeamE->Fill( etapr_6gTrue.GetTrueBeamEnergy() );
    Double_t beam_e = etapr_6gTrue.GetTrueBeamEnergy();

    etapr_true[0] = etapr_6gTrue.GetTrueEtaLV();
    etapr_true[1] = etapr_6gTrue.GetTrueNeutralPiLV(0);
    etapr_true[2] = etapr_6gTrue.GetTrueNeutralPiLV(1);

    Double_t weight = Get_etapr_weight_MC(beam_e, etapr_true);

    // calculate Physics

    // here calculate the boost
    TVector3 etaprime_rest  = -(etapr_true[0]+etapr_true[1]+etapr_true[2]).BoostVector();
    TLorentzVector etapr_rest  = etapr_6gTrue.GetTrueEtaPrimeLV();
    etapr_rest.Boost(etaprime_rest);

    // here calculate the eta prime theta as function of beam energy
    bw = 0.20;
    DalitzPlot(etapr_true, Xtrue, Ytrue, bw, DPnrTrue020);
    bw = 0.15;
    DalitzPlot(etapr_true, Xtrue, Ytrue, bw, DPnrTrue015);
    bw = 0.10;
    DalitzPlot(etapr_true, Xtrue, Ytrue, bw, DPnrTrue010);
    bw = 0.075;
    DalitzPlot(etapr_true, Xtrue, Ytrue, bw, DPnrTrue075);
    bw = 0.05;
    DalitzPlot(etapr_true, Xtrue, Ytrue, bw, DPnrTrue005);

    Double_t weight2 = 1.0;
    // Weigh the MC with the Dalitz plot parameters.
    if(s == "Sergey"){
        N = 1.0e5/94190;
        a = -0.071;
        b = -0.070;
        d = -0.061;
        Double_t weight2 = N*( 1+ a*Ytrue + b*Ytrue*Ytrue + d*Xtrue*Xtrue );
    }
    else // Phase Space
    {
        weight2 = 1.0;
    }

    weight *= weight2;
    etapr_6gTrue.SetWeight(weight);

    true_phy_DP_020->Fill(DPnrTrue020,etapr_6gTrue.GetWeight());
    true_phy_DP_015->Fill(DPnrTrue015,etapr_6gTrue.GetWeight());
    true_phy_DP_010->Fill(DPnrTrue010,etapr_6gTrue.GetWeight());
    true_phy_DP_075->Fill(DPnrTrue005,etapr_6gTrue.GetWeight());
    true_phy_DP_005->Fill(DPnrTrue005,etapr_6gTrue.GetWeight());


    true_BeamE_weight->FillWeighted(etapr_6gTrue.GetTrueBeamEnergy(),etapr_6gTrue.GetWeight());

    true_DP->Fill( Xtrue, Ytrue,etapr_6gTrue.GetWeight() );

    m2pi0_metapi0(etapr_true, m_etapi01True, m_etapi02True, m_2pi0True);
    true_M_pi1pi2_e2p->Fill(m_2pi0True*1.0e3);
    true_M_etapi_e2p->Fill(m_etapi01True*1.0e3);
    true_M_etapi_e2p->Fill(m_etapi02True*1.0e3);

    // calculate True Energy vs Theta of final state particles

    true_th_v_E_p->Fill(etapr_6gTrue.GetTrueProtonLV().E() - MASS_PROTON/1000 ,etapr_6gTrue.GetTrueProtonLV().Theta()*TMath::RadToDeg());

    for(UInt_t i = 0; i < etapr_6gTrue.GetNgamma(); i++ )
    {
        if ( i > 3 ) // first four gammas come from 2pi0
            true_th_v_E_eta_g->Fill(etapr_6gTrue.GetTrueGammaLV(i).E()*1000.,etapr_6gTrue.GetTrueGammaLV(i).Theta()*TMath::RadToDeg());
        else
            true_th_v_E_pi0_g->Fill(etapr_6gTrue.GetTrueGammaLV(i).E()*1000.,etapr_6gTrue.GetTrueGammaLV(i).Theta()*TMath::RadToDeg());

    }
    return;
}

void AdlarsonPhysics::TrueAnalysis_etapr10g(){
    true_BeamE->Fill( etapr_10gTrue.GetTrueBeamEnergy() );
    Double_t beam_e = etapr_10gTrue.GetTrueBeamEnergy();

    etapr_true[0] = etapr_10gTrue.GetTrueEtaLV();
    etapr_true[1] = etapr_10gTrue.GetTrueNeutralPiLV(0);
    etapr_true[2] = etapr_10gTrue.GetTrueNeutralPiLV(1);
    Double_t weight = Get_etapr_weight_MC( beam_e, etapr_true );

    etapr_10gTrue.SetWeight(weight);

    if(TMath::Abs(weight) < 1.0e-5)
        weight = 1.;

    // calculate Physics

    etapr_true[0] = etapr_10gTrue.GetTrueEtaLV();
    etapr_true[1] = etapr_10gTrue.GetTrueNeutralPiLV(0);
    etapr_true[2] = etapr_10gTrue.GetTrueNeutralPiLV(1);
    DalitzPlot(etapr_true, Xtrue, Ytrue, bw, DPnrTrue020);
    true_DP->Fill( Xtrue, Ytrue, weight);
    true_phy_DP_020->Fill(DPnrTrue020, weight);

    m2pi0_metapi0(etapr_true, m_etapi01True, m_etapi02True, m_2pi0True);
    true_M_pi1pi2_e2p->Fill(m_2pi0True*1.0e3, weight);
    true_M_etapi_e2p->Fill(m_etapi01True*1.0e3, weight);
    true_M_etapi_e2p->Fill(m_etapi02True*1.0e3, weight);

    // calculate True Energy vs Theta of final state particles

    true_th_v_E_p->Fill(etapr_10gTrue.GetTrueProtonLV().E() - MASS_PROTON/1000 ,etapr_10gTrue.GetTrueProtonLV().Theta()*TMath::RadToDeg());

    for(UInt_t i = 0; i < etapr_10gTrue.GetNgamma(); i++ )
    {
        if ( i > 3 ) // first four gammas come from 2pi0
            true_th_v_E_eta_6g->Fill(etapr_10gTrue.GetTrueGammaLV(i).E(),etapr_10gTrue.GetTrueGammaLV(i).Theta()*TMath::RadToDeg(), weight);
        else
            true_th_v_E_pi0_g->Fill(etapr_10gTrue.GetTrueGammaLV(i).E(),etapr_10gTrue.GetTrueGammaLV(i).Theta()*TMath::RadToDeg(), weight);
    }
    return;
}

Double_t AdlarsonPhysics::Get_etapr_weight_MC(Double_t beame, TLorentzVector etapr_true[2]){

    Double_t w = 1.0;
    Double_t N, c1, c2, c3, c4;
    Double_t leg0, leg1, leg2, leg3, leg4;
    Double_t weight, w1, w2;
    Double_t be = beame*1000;
    Int_t ee = -100;
    TLorentzVector sqrt_s(0.,0., be, be+MASS_PROTON);

    double width;

//    etapr_true[0] = etapr_6gTrue.GetTrueEtaLV();
//    etapr_true[1] = etapr_6gTrue.GetTrueNeutralPiLV(0);
//    etapr_true[2] = etapr_6gTrue.GetTrueNeutralPiLV(1);


//    etapr_true[0] = etapr_10gTrue.GetTrueEtaLV();
//    etapr_true[1] = etapr_10gTrue.GetTrueNeutralPiLV(0);
//    etapr_true[2] = etapr_10gTrue.GetTrueNeutralPiLV(1);

    // here calculate the boost
    TVector3 cm_vector  = -(sqrt_s).BoostVector();
    TLorentzVector etapr_cm  = etapr_6gTrue.GetTrueEtaPrimeLV();
//    TLorentzVector etapr_cm  = etapr_10gTrue.GetTrueEtaPrimeLV();
    etapr_cm.Boost(cm_vector);

    //in the CM frame!
    Double_t angle = etapr_cm.Theta();
    Double_t x = TMath::Cos( angle );

    leg0 = 1.;
    leg1 = x;
    leg2 = (0.5)*( 3.*TMath::Power(x,2) - 1. );
    leg3 = (0.5)*( 5.*TMath::Power(x,3) - 3.*x );
    leg4 = (0.125)*( 35.*TMath::Power(x,4) - 30.*TMath::Power(x,2.) + 3.);

    if(be <= Legendre[0]){
        width = 6.5;
        N   =   Legendre[0+1];
        c1  =   Legendre[0+2];
        c2  =   Legendre[0+3];
        c3  =   Legendre[0+4];
        c4  =   Legendre[0+5];

        w = N*leg0 + c1*leg1 + c2*leg2 + c3*leg3 + c4*leg4;
        w /= width;
        ee = 0;
    }
    else if(be >= (Legendre[66])){
        width = 13.;
        N   =   Legendre[66+1];
        c1  =   Legendre[66+2];
        c2  =   Legendre[66+3];
        c3  =   Legendre[66+4];
        c4  =   Legendre[66+5];

        w = N*leg0 + c1*leg1 + c2*leg2 + c3*leg3 + c4*leg4;
        w /= width;
        ee = 11;
    }
    else{
        bool be_range = true;
        int k = 0;
        while( be_range && ( k < 61 ) ){
            if( be >= Legendre[k] && be < Legendre[k+6] ){
                be_range = false;
                ee = k/6;
            }
            else
                k += 6;
        }

        N   =   Legendre[k+1];
        c1  =   Legendre[k+2];
        c2  =   Legendre[k+3];
        c3  =   Legendre[k+4];
        c4  =   Legendre[k+5];

        w1 =  N*leg0 + c1*leg1 + c2*leg2 + c3*leg3 + c4*leg4;

        N   =   Legendre[k+7];
        c1  =   Legendre[k+8];
        c2  =   Legendre[k+9];
        c3  =   Legendre[k+10];
        c4  =   Legendre[k+11];

        w2 =  N*leg0 + c1*leg1 + c2*leg2 + c3*leg3 + c4*leg4;

        //should add up to 1

        if(be < 1469.8){
            width = 6.5;
            double frac1 = (1- (be-Legendre[k+0])/width);
            double frac2 = (1- (Legendre[k+6]-be)/width);
            w = (frac1*w1 + frac2*w2)/width;
        }
        else if((be >= 1469.8) && (be < 1473.05 )){
            width = 6.5;
            double frac1 = (1- (be-Legendre[k+0])/width);
            double frac2 = (1- (Legendre[k+6]-be)/width);
//            w = (frac1*w1 + frac2*w2)/width;
            w = w1/width;
        }
        else if((be >= 1473.05) && (be < 1479.5 )){
            width = 13.;
//            double frac1 = (1- (be-Legendre[k+0])/width);
//            double frac2 = (1- (Legendre[k+6]-be)/width);
//            w = (frac1*w1 + frac2*w2)/width;
            w = w2/width;
        }

        else{
            width = 13.;
            double frac1 = (1- (be-Legendre[k+0])/width);
            double frac2 = (1- (Legendre[k+6]-be)/width);
            w = (frac1*w1 + frac2*w2)/width;
        }
    }

//        w /= width;

     w /= 104.440282;
     if(be < 1473.05)
         w *= 1.1;

    TLorentzVector eta_pr_all = etapr_true[0] + etapr_true[1] + etapr_true[2];
    Int_t bin =  diff_distr( be, eta_pr_all );
    bin = ee*20 + int((1. + x)/0.2);

    true_eta_pr_gg_effcorr->Fill(bin,w);
    true_etapr_diff_distr->Fill(bin,w);

    // here split up into different beams;

    true_eta_pr_production->Fill(bin,w);

    true_norm->Fill(1,w);
    MC_weight = true;

    return w;

}

Double_t AdlarsonPhysics::TrueAnalysis_threepi_etapi(){
    TLorentzVector true_im(0.0, 0.0, 0.0, 0.0);

    MCw = 1.0;
    for(UInt_t ig = 0; ig < threepi_etapi.GetNgamma(); ig++)
        true_im += threepi_etapi.GetTrueGammaLV(ig);

    if(threepi_etapi.GetNgamma() == 6){
     Double_t MCw_temp = 1.0; Double_t MCw_temp1 = 1.0; Double_t MCw_temp2 = 1.0; Double_t MCw_temp3 = 1.0;
     Double_t MCw_tempb = 1.0; Double_t MCw_temp4 = 1.0; Double_t MCw_temp5 = 1.0; Double_t MCw_temp6 = 1.0;
     Double_t MCw_tempc = 1.0; Double_t MCw_temp7 = 1.0; Double_t MCw_temp8 = 1.0; Double_t MCw_temp9 = 1.0;

     TLorentzVector three_pi0, three_pi1, three_pi2, proton_LV;
     three_pi0 = threepi_etapi.GetTrueNeutralPiLV(0);
     three_pi1 = threepi_etapi.GetTrueNeutralPiLV(1);
     three_pi2 = threepi_etapi.GetTrueNeutralPiLV(2);
     proton_LV = threepi_etapi.GetTrueProtonLV();

     MCw_temp1 = GetWeight3pi1( (three_pi0+three_pi1).M2(), (three_pi2+proton_LV).M2());
     if(MCw_temp1 > 0.){
        MCw_temp *= MCw_temp1;
     }
     else
        MCw_temp *=1.0;

     MCw_temp2 = GetWeight3pi1( (three_pi0+three_pi2).M2(), (three_pi1+proton_LV).M2());
     if(MCw_temp2 > 0.){
        MCw_temp *= MCw_temp2;
     }
     else
        MCw_temp *=1.0;

     MCw_temp3 = GetWeight3pi1( (three_pi1+three_pi2).M2(), (three_pi0+proton_LV).M2());
     if(MCw_temp3 > 0.){
        MCw_temp *= MCw_temp3;
     }
     else
        MCw_temp *=1.0;

     MCw = TMath::Power(MCw_temp, 1./3.);

     MCw_temp4 = GetWeight3pi2( (three_pi0+three_pi1).M2(), (three_pi2+proton_LV).M2());
     if(MCw_temp4 > 0.){
        MCw_tempb *= MCw_temp4;
     }
     else
        MCw_tempb *=1.0;

     MCw_temp5 = GetWeight3pi2( (three_pi0+three_pi2).M2(), (three_pi1+proton_LV).M2());
     if(MCw_temp5 > 0.){
        MCw_tempb *= MCw_temp5;
     }
     else
        MCw_tempb *=1.0;

     MCw_temp6 = GetWeight3pi2( (three_pi1+three_pi2).M2(), (three_pi0+proton_LV).M2());
     if(MCw_temp6 > 0.){
        MCw_tempb *= MCw_temp6;
     }
     else
        MCw_tempb *=1.0;

     MCw *= TMath::Power(MCw_tempb, 1./3.);

     MCw_temp7 = GetWeight3pi3( (three_pi0+three_pi1).M2(), (three_pi2+proton_LV).M2());
     if(MCw_temp7 > 0.){
        MCw_tempc *= MCw_temp7;
     }
     else
        MCw_tempc *=1.0;

     MCw_temp8 = GetWeight3pi3( (three_pi0+three_pi2).M2(), (three_pi1+proton_LV).M2());
     if(MCw_temp8 > 0.){
        MCw_tempc *= MCw_temp8;
     }
     else
        MCw_tempc *=1.0;

     MCw_temp9 = GetWeight3pi3( (three_pi1+three_pi2).M2(), (three_pi0+proton_LV).M2());
     if(MCw_temp9 > 0.){
        MCw_tempc *= MCw_temp9;
     }
     else
        MCw_tempc *=1.0;

     MCw *= TMath::Power(MCw_tempc, 1./3.);

     Double_t M = true_im.M()*1.0e3;
     Double_t Mw = M;
        if(M > 830.){
            true_phy_3pi_IMpipi_v_IMppi->Fill((three_pi0+three_pi1).M2(), (three_pi2+proton_LV).M2());
            true_phy_3pi_IMpipi_v_IMppi->Fill((three_pi0+three_pi2).M2(), (three_pi1+proton_LV).M2());
            true_phy_3pi_IMpipi_v_IMppi->Fill((three_pi1+three_pi2).M2(), (three_pi0+proton_LV).M2());
        }
    }
    return MCw;
}

double AdlarsonPhysics::GetWeight3pi1(Double_t M1sq, Double_t M2sq){

    double r_c, r_temp, r_sum;
    double ms1_c, ms2_c, ms1_temp, ms2_temp;
    double r_min = 1.0e5;
    double MCw_center, MCw_neighbor;
    int x_n, y_n;

    int xc = MCw_bkgd->GetXaxis()->FindBin(M1sq);
    int yc = MCw_bkgd->GetYaxis()->FindBin(M2sq);

    ms1_c = MCw_bkgd->GetXaxis()->GetBinCenter(xc);
    ms2_c = MCw_bkgd->GetYaxis()->GetBinCenter(yc);

    r_c = TMath::Sqrt(TMath::Abs((MCw_bkgd->GetXaxis()->GetBinCenter(xc)*MCw_bkgd->GetXaxis()->GetBinCenter(xc)-M1sq*M1sq)) + TMath::Abs((MCw_bkgd->GetYaxis()->GetBinCenter(yc)*MCw_bkgd->GetYaxis()->GetBinCenter(yc)-M2sq*M2sq)));
    // test all 8 possibilities, from theta 0 - 360 deg, counter clock-wise
    // 1) xbin + 1, ybin + 0
    r_temp = TMath::Sqrt(TMath::Abs((MCw_bkgd->GetXaxis()->GetBinCenter(xc+1)*MCw_bkgd->GetXaxis()->GetBinCenter(xc+1)-M1sq*M1sq)) + TMath::Abs((MCw_bkgd->GetYaxis()->GetBinCenter(yc)*MCw_bkgd->GetYaxis()->GetBinCenter(yc)-M2sq*M2sq)));
    ms1_temp =  MCw_bkgd->GetXaxis()->GetBinCenter(xc+1);
    ms2_temp = MCw_bkgd->GetYaxis()->GetBinCenter(yc);
    if(r_temp < r_min){
        r_min = r_temp;
        x_n = xc + 1;
        y_n = yc ;
    }
    // 2) xbin + 1, ybin + 1
    r_temp = TMath::Sqrt(TMath::Abs((MCw_bkgd->GetXaxis()->GetBinCenter(xc+1)*MCw_bkgd->GetXaxis()->GetBinCenter(xc+1)-M1sq*M1sq)) + TMath::Abs((MCw_bkgd->GetYaxis()->GetBinCenter(yc+1)*MCw_bkgd->GetYaxis()->GetBinCenter(yc+1)-M2sq*M2sq)));
    if(r_temp < r_min){
        r_min = r_temp;
        x_n = xc + 1;
        y_n = yc + 1 ;
    }
    // 3) xbin + 0, ybin + 1
    r_temp = TMath::Sqrt(TMath::Abs((MCw_bkgd->GetXaxis()->GetBinCenter(xc)*MCw_bkgd->GetXaxis()->GetBinCenter(xc)-M1sq*M1sq)) + TMath::Abs((MCw_bkgd->GetYaxis()->GetBinCenter(yc+1)*MCw_bkgd->GetYaxis()->GetBinCenter(yc+1)-M2sq*M2sq)));
    if(r_temp < r_min){
        r_min = r_temp;
        x_n = xc ;
        y_n = yc + 1 ;
    }
    // 4) xbin - 1, ybin + 1
    r_temp = TMath::Sqrt(TMath::Abs((MCw_bkgd->GetXaxis()->GetBinCenter(xc-1)*MCw_bkgd->GetXaxis()->GetBinCenter(xc-1)-M1sq*M1sq)) + TMath::Abs((MCw_bkgd->GetYaxis()->GetBinCenter(yc+1)*MCw_bkgd->GetYaxis()->GetBinCenter(yc+1)-M2sq*M2sq)));
    if(r_temp < r_min){
        r_min = r_temp;
        x_n = xc - 1;
        y_n = yc + 1;
    }
    // 5) xbin - 1, ybin + 0
    r_temp = TMath::Sqrt(TMath::Abs((MCw_bkgd->GetXaxis()->GetBinCenter(xc-1)*MCw_bkgd->GetXaxis()->GetBinCenter(xc-1)-M1sq*M1sq)) + TMath::Abs((MCw_bkgd->GetYaxis()->GetBinCenter(yc)*MCw_bkgd->GetYaxis()->GetBinCenter(yc)-M2sq*M2sq)));
    if(r_temp < r_min){
        r_min = r_temp;
        x_n = xc - 1;
        y_n = yc ;
    }
    // 6) xbin - 1, ybin - 1
    r_temp = TMath::Sqrt(TMath::Abs((MCw_bkgd->GetXaxis()->GetBinCenter(xc-1)*MCw_bkgd->GetXaxis()->GetBinCenter(xc-1)-M1sq*M1sq)) + TMath::Abs((MCw_bkgd->GetYaxis()->GetBinCenter(yc-1)*MCw_bkgd->GetYaxis()->GetBinCenter(yc-1)-M2sq*M2sq)));
    if(r_temp < r_min){
        r_min = r_temp;
        x_n = xc - 1;
        y_n = yc - 1;
    }
    // 7) xbin + 0, ybin - 1
    r_temp = TMath::Sqrt(TMath::Abs((MCw_bkgd->GetXaxis()->GetBinCenter(xc)*MCw_bkgd->GetXaxis()->GetBinCenter(xc)-M1sq*M1sq)) + TMath::Abs((MCw_bkgd->GetYaxis()->GetBinCenter(yc-1)*MCw_bkgd->GetYaxis()->GetBinCenter(yc-1)-M2sq*M2sq)));
    if(r_temp < r_min){
        r_min = r_temp;
        x_n = xc ;
        y_n = yc - 1;
    }
    // 8) xbin + 1, ybin - 1
    r_temp = TMath::Sqrt(TMath::Abs((MCw_bkgd->GetXaxis()->GetBinCenter(xc+1)*MCw_bkgd->GetXaxis()->GetBinCenter(xc+1)-M1sq*M1sq)) + TMath::Abs((MCw_bkgd->GetYaxis()->GetBinCenter(yc-1)*MCw_bkgd->GetYaxis()->GetBinCenter(yc-1)-M2sq*M2sq)));
    if(r_temp < r_min){
        r_min = r_temp;
        x_n = xc + 1;
        y_n = yc - 1;
    }

    MCw_center    = MCw_bkgd->GetBinContent(xc,yc);
    MCw_neighbor  = MCw_bkgd->GetBinContent(x_n,y_n);

    r_sum = r_min + r_c;

    return MCw_center*(r_c/r_sum) + MCw_neighbor*(r_min/r_sum);

}

double AdlarsonPhysics::GetWeight3pi2(Double_t M1sq, Double_t M2sq){

    double r_c, r_temp, r_sum;
    double ms1_c, ms2_c, ms1_temp, ms2_temp;
    double r_min = 1.0e5;
    double MCw_center, MCw_neighbor;
    int x_n, y_n;

    int xc = MCw_bkgd2->GetXaxis()->FindBin(M1sq);
    int yc = MCw_bkgd2->GetYaxis()->FindBin(M2sq);

    ms1_c = MCw_bkgd2->GetXaxis()->GetBinCenter(xc);
    ms2_c = MCw_bkgd2->GetYaxis()->GetBinCenter(yc);

    r_c = TMath::Sqrt(TMath::Abs((MCw_bkgd2->GetXaxis()->GetBinCenter(xc)*MCw_bkgd2->GetXaxis()->GetBinCenter(xc)-M1sq*M1sq)) + TMath::Abs((MCw_bkgd2->GetYaxis()->GetBinCenter(yc)*MCw_bkgd2->GetYaxis()->GetBinCenter(yc)-M2sq*M2sq)));
    // test all 8 possibilities, from theta 0 - 360 deg, counter clock-wise
    // 1) xbin + 1, ybin + 0
    r_temp = TMath::Sqrt(TMath::Abs((MCw_bkgd2->GetXaxis()->GetBinCenter(xc+1)*MCw_bkgd2->GetXaxis()->GetBinCenter(xc+1)-M1sq*M1sq)) + TMath::Abs((MCw_bkgd2->GetYaxis()->GetBinCenter(yc)*MCw_bkgd2->GetYaxis()->GetBinCenter(yc)-M2sq*M2sq)));
    ms1_temp =  MCw_bkgd2->GetXaxis()->GetBinCenter(xc+1);
    ms2_temp = MCw_bkgd2->GetYaxis()->GetBinCenter(yc);
    if(r_temp < r_min){
        r_min = r_temp;
        x_n = xc + 1;
        y_n = yc ;
    }
    // 2) xbin + 1, ybin + 1
    r_temp = TMath::Sqrt(TMath::Abs((MCw_bkgd2->GetXaxis()->GetBinCenter(xc+1)*MCw_bkgd2->GetXaxis()->GetBinCenter(xc+1)-M1sq*M1sq)) + TMath::Abs((MCw_bkgd2->GetYaxis()->GetBinCenter(yc+1)*MCw_bkgd2->GetYaxis()->GetBinCenter(yc+1)-M2sq*M2sq)));
    if(r_temp < r_min){
        r_min = r_temp;
        x_n = xc + 1;
        y_n = yc + 1 ;
    }
    // 3) xbin + 0, ybin + 1
    r_temp = TMath::Sqrt(TMath::Abs((MCw_bkgd2->GetXaxis()->GetBinCenter(xc)*MCw_bkgd2->GetXaxis()->GetBinCenter(xc)-M1sq*M1sq)) + TMath::Abs((MCw_bkgd2->GetYaxis()->GetBinCenter(yc+1)*MCw_bkgd2->GetYaxis()->GetBinCenter(yc+1)-M2sq*M2sq)));
    if(r_temp < r_min){
        r_min = r_temp;
        x_n = xc ;
        y_n = yc + 1 ;
    }
    // 4) xbin - 1, ybin + 1
    r_temp = TMath::Sqrt(TMath::Abs((MCw_bkgd2->GetXaxis()->GetBinCenter(xc-1)*MCw_bkgd2->GetXaxis()->GetBinCenter(xc-1)-M1sq*M1sq)) + TMath::Abs((MCw_bkgd2->GetYaxis()->GetBinCenter(yc+1)*MCw_bkgd2->GetYaxis()->GetBinCenter(yc+1)-M2sq*M2sq)));
    if(r_temp < r_min){
        r_min = r_temp;
        x_n = xc - 1;
        y_n = yc + 1;
    }
    // 5) xbin - 1, ybin + 0
    r_temp = TMath::Sqrt(TMath::Abs((MCw_bkgd2->GetXaxis()->GetBinCenter(xc-1)*MCw_bkgd2->GetXaxis()->GetBinCenter(xc-1)-M1sq*M1sq)) + TMath::Abs((MCw_bkgd2->GetYaxis()->GetBinCenter(yc)*MCw_bkgd2->GetYaxis()->GetBinCenter(yc)-M2sq*M2sq)));
    if(r_temp < r_min){
        r_min = r_temp;
        x_n = xc - 1;
        y_n = yc ;
    }
    // 6) xbin - 1, ybin - 1
    r_temp = TMath::Sqrt(TMath::Abs((MCw_bkgd2->GetXaxis()->GetBinCenter(xc-1)*MCw_bkgd2->GetXaxis()->GetBinCenter(xc-1)-M1sq*M1sq)) + TMath::Abs((MCw_bkgd2->GetYaxis()->GetBinCenter(yc-1)*MCw_bkgd2->GetYaxis()->GetBinCenter(yc-1)-M2sq*M2sq)));
    if(r_temp < r_min){
        r_min = r_temp;
        x_n = xc - 1;
        y_n = yc - 1;
    }
    // 7) xbin + 0, ybin - 1
    r_temp = TMath::Sqrt(TMath::Abs((MCw_bkgd2->GetXaxis()->GetBinCenter(xc)*MCw_bkgd2->GetXaxis()->GetBinCenter(xc)-M1sq*M1sq)) + TMath::Abs((MCw_bkgd2->GetYaxis()->GetBinCenter(yc-1)*MCw_bkgd2->GetYaxis()->GetBinCenter(yc-1)-M2sq*M2sq)));
    if(r_temp < r_min){
        r_min = r_temp;
        x_n = xc ;
        y_n = yc - 1;
    }
    // 8) xbin + 1, ybin - 1
    r_temp = TMath::Sqrt(TMath::Abs((MCw_bkgd2->GetXaxis()->GetBinCenter(xc+1)*MCw_bkgd2->GetXaxis()->GetBinCenter(xc+1)-M1sq*M1sq)) + TMath::Abs((MCw_bkgd2->GetYaxis()->GetBinCenter(yc-1)*MCw_bkgd2->GetYaxis()->GetBinCenter(yc-1)-M2sq*M2sq)));
    if(r_temp < r_min){
        r_min = r_temp;
        x_n = xc + 1;
        y_n = yc - 1;
    }

    MCw_center    = MCw_bkgd2->GetBinContent(xc,yc);
    MCw_neighbor  = MCw_bkgd2->GetBinContent(x_n,y_n);

    r_sum = r_min + r_c;

    return MCw_center*(r_c/r_sum) + MCw_neighbor*(r_min/r_sum);

}

double AdlarsonPhysics::GetWeight3pi3(Double_t M1sq, Double_t M2sq){

    double r_c, r_temp, r_sum;
    double ms1_c, ms2_c, ms1_temp, ms2_temp;
    double r_min = 1.0e5;
    double MCw_center, MCw_neighbor;
    int x_n, y_n;

    int xc = MCw_bkgd3->GetXaxis()->FindBin(M1sq);
    int yc = MCw_bkgd3->GetYaxis()->FindBin(M2sq);

    ms1_c = MCw_bkgd3->GetXaxis()->GetBinCenter(xc);
    ms2_c = MCw_bkgd3->GetYaxis()->GetBinCenter(yc);

    r_c = TMath::Sqrt(TMath::Abs((MCw_bkgd3->GetXaxis()->GetBinCenter(xc)*MCw_bkgd3->GetXaxis()->GetBinCenter(xc)-M1sq*M1sq)) + TMath::Abs((MCw_bkgd3->GetYaxis()->GetBinCenter(yc)*MCw_bkgd3->GetYaxis()->GetBinCenter(yc)-M2sq*M2sq)));
    // test all 8 possibilities, from theta 0 - 360 deg, counter clock-wise
    // 1) xbin + 1, ybin + 0
    r_temp = TMath::Sqrt(TMath::Abs((MCw_bkgd3->GetXaxis()->GetBinCenter(xc+1)*MCw_bkgd3->GetXaxis()->GetBinCenter(xc+1)-M1sq*M1sq)) + TMath::Abs((MCw_bkgd3->GetYaxis()->GetBinCenter(yc)*MCw_bkgd3->GetYaxis()->GetBinCenter(yc)-M2sq*M2sq)));
    ms1_temp =  MCw_bkgd3->GetXaxis()->GetBinCenter(xc+1);
    ms2_temp = MCw_bkgd3->GetYaxis()->GetBinCenter(yc);
    if(r_temp < r_min){
        r_min = r_temp;
        x_n = xc + 1;
        y_n = yc ;
    }
    // 2) xbin + 1, ybin + 1
    r_temp = TMath::Sqrt(TMath::Abs((MCw_bkgd3->GetXaxis()->GetBinCenter(xc+1)*MCw_bkgd3->GetXaxis()->GetBinCenter(xc+1)-M1sq*M1sq)) + TMath::Abs((MCw_bkgd3->GetYaxis()->GetBinCenter(yc+1)*MCw_bkgd3->GetYaxis()->GetBinCenter(yc+1)-M2sq*M2sq)));
    if(r_temp < r_min){
        r_min = r_temp;
        x_n = xc + 1;
        y_n = yc + 1 ;
    }
    // 3) xbin + 0, ybin + 1
    r_temp = TMath::Sqrt(TMath::Abs((MCw_bkgd3->GetXaxis()->GetBinCenter(xc)*MCw_bkgd3->GetXaxis()->GetBinCenter(xc)-M1sq*M1sq)) + TMath::Abs((MCw_bkgd3->GetYaxis()->GetBinCenter(yc+1)*MCw_bkgd3->GetYaxis()->GetBinCenter(yc+1)-M2sq*M2sq)));
    if(r_temp < r_min){
        r_min = r_temp;
        x_n = xc ;
        y_n = yc + 1 ;
    }
    // 4) xbin - 1, ybin + 1
    r_temp = TMath::Sqrt(TMath::Abs((MCw_bkgd3->GetXaxis()->GetBinCenter(xc-1)*MCw_bkgd3->GetXaxis()->GetBinCenter(xc-1)-M1sq*M1sq)) + TMath::Abs((MCw_bkgd3->GetYaxis()->GetBinCenter(yc+1)*MCw_bkgd3->GetYaxis()->GetBinCenter(yc+1)-M2sq*M2sq)));
    if(r_temp < r_min){
        r_min = r_temp;
        x_n = xc - 1;
        y_n = yc + 1;
    }
    // 5) xbin - 1, ybin + 0
    r_temp = TMath::Sqrt(TMath::Abs((MCw_bkgd3->GetXaxis()->GetBinCenter(xc-1)*MCw_bkgd3->GetXaxis()->GetBinCenter(xc-1)-M1sq*M1sq)) + TMath::Abs((MCw_bkgd3->GetYaxis()->GetBinCenter(yc)*MCw_bkgd3->GetYaxis()->GetBinCenter(yc)-M2sq*M2sq)));
    if(r_temp < r_min){
        r_min = r_temp;
        x_n = xc - 1;
        y_n = yc ;
    }
    // 6) xbin - 1, ybin - 1
    r_temp = TMath::Sqrt(TMath::Abs((MCw_bkgd3->GetXaxis()->GetBinCenter(xc-1)*MCw_bkgd3->GetXaxis()->GetBinCenter(xc-1)-M1sq*M1sq)) + TMath::Abs((MCw_bkgd3->GetYaxis()->GetBinCenter(yc-1)*MCw_bkgd3->GetYaxis()->GetBinCenter(yc-1)-M2sq*M2sq)));
    if(r_temp < r_min){
        r_min = r_temp;
        x_n = xc - 1;
        y_n = yc - 1;
    }
    // 7) xbin + 0, ybin - 1
    r_temp = TMath::Sqrt(TMath::Abs((MCw_bkgd3->GetXaxis()->GetBinCenter(xc)*MCw_bkgd3->GetXaxis()->GetBinCenter(xc)-M1sq*M1sq)) + TMath::Abs((MCw_bkgd3->GetYaxis()->GetBinCenter(yc-1)*MCw_bkgd3->GetYaxis()->GetBinCenter(yc-1)-M2sq*M2sq)));
    if(r_temp < r_min){
        r_min = r_temp;
        x_n = xc ;
        y_n = yc - 1;
    }
    // 8) xbin + 1, ybin - 1
    r_temp = TMath::Sqrt(TMath::Abs((MCw_bkgd3->GetXaxis()->GetBinCenter(xc+1)*MCw_bkgd3->GetXaxis()->GetBinCenter(xc+1)-M1sq*M1sq)) + TMath::Abs((MCw_bkgd3->GetYaxis()->GetBinCenter(yc-1)*MCw_bkgd3->GetYaxis()->GetBinCenter(yc-1)-M2sq*M2sq)));
    if(r_temp < r_min){
        r_min = r_temp;
        x_n = xc + 1;
        y_n = yc - 1;
    }

    MCw_center    = MCw_bkgd3->GetBinContent(xc,yc);
    MCw_neighbor  = MCw_bkgd3->GetBinContent(x_n,y_n);

    r_sum = r_min + r_c;

    return MCw_center*(r_c/r_sum) + MCw_neighbor*(r_min/r_sum);

}

double AdlarsonPhysics::GetGain(Double_t E, Double_t detnr){

    double  gain_c, gain_n, gain;
    double  E_c, E_n, de;
    double  frac_c, frac_n;

    de = 40.;

    int xc  = Eth_gamma->GetXaxis()->FindBin(E);
    int yc  = Eth_gamma->GetYaxis()->FindBin(detnr);

    E_c     = Eth_gamma->GetXaxis()->GetBinCenter(xc);
    gain_c = Eth_gamma->GetBinContent(xc,yc); // in the same bin
    if( E < E_c ){
        E_n     = Eth_gamma->GetXaxis()->GetBinCenter(xc-1);
        gain_n  = Eth_gamma->GetBinContent(xc-1,yc);
    }
    else{
        E_n     = Eth_gamma->GetXaxis()->GetBinCenter(xc+1);
        gain_n  = Eth_gamma->GetBinContent(xc+1,yc);
    }

    // 4 cases-
        // i)   c and n have content
        // ii)  c has content, not n
        // iii) n has content, not c
        // iv)  no content c nor n

    if((TMath::Abs(gain_c) > 1.0e-3) && (TMath::Abs(gain_n) > 1.0e-3)){
        frac_c = 1-TMath::Abs(E-E_c)/de;
        frac_n = 1-TMath::Abs(E-E_n)/de;
        gain = frac_c*gain_c + frac_n*gain_n;
    }
    else if((TMath::Abs(gain_c) > 1.0e-3) && (TMath::Abs(gain_n) < 1.0e-3)){
        gain = gain_c;
    }
    else if((TMath::Abs(gain_n) > 1.0e-3) && (TMath::Abs(gain_c) < 1.0e-3)){
        gain = gain_n;
    }
    else{
        int ibin = 1;
        while((TMath::Abs(Eth_gamma->GetBinContent(xc-ibin,yc))<1.0e-3) && (ibin < 10))
            ibin++;
        gain = Eth_gamma->GetBinContent(xc-ibin,yc);

        if(TMath::Abs(gain)<1.0e-3)
            gain =1.;

    }
    if(TMath::Abs(gain) < 1.0e-3)
            gain = 1.0;

    return gain;

}

double AdlarsonPhysics::GetGainNLO(Double_t E, Double_t detnr){

    double  gain_c, gain_n, gain;
    double  E_c, E_n, de;
    double  frac_c, frac_n;
    int     x_n, y_n;

    de = 40.;

    int xc  = Eth_gammaNLO->GetXaxis()->FindBin(E);
    int yc  = Eth_gammaNLO->GetYaxis()->FindBin(detnr);

    E_c     = Eth_gammaNLO->GetXaxis()->GetBinCenter(xc);
    gain_c = Eth_gammaNLO->GetBinContent(xc,yc); // in the same bin
    if( E < E_c ){
        E_n     = Eth_gammaNLO->GetXaxis()->GetBinCenter(xc-1);
        gain_n  = Eth_gammaNLO->GetBinContent(xc-1,yc);
    }
    else{
        E_n     = Eth_gammaNLO->GetXaxis()->GetBinCenter(xc+1);
        gain_n  = Eth_gammaNLO->GetBinContent(xc+1,yc);
    }

    // 4 cases-
        // i)   c and n have content
        // ii)  c has content, not n
        // iii) n has content, not c
        // iv)  no content c nor n

    if((TMath::Abs(gain_c) > 1.0e-3) && (TMath::Abs(gain_n) > 1.0e-3)){
        frac_c = 1-TMath::Abs(E-E_c)/de;
        frac_n = 1-TMath::Abs(E-E_n)/de;
        gain = frac_c*gain_c + frac_n*gain_n;
    }
    else if((TMath::Abs(gain_c) > 1.0e-3) && (TMath::Abs(gain_n) < 1.0e-3)){
        gain = gain_c;
    }
    else if((TMath::Abs(gain_n) > 1.0e-3) && (TMath::Abs(gain_c) < 1.0e-3)){
        gain = gain_n;
    }
    else{
        int ibin = 1;
        while((TMath::Abs(Eth_gammaNLO->GetBinContent(xc-ibin,yc))<1.0e-3) && (ibin < 10))
            ibin++;
        gain = Eth_gammaNLO->GetBinContent(xc-ibin,yc);

        if(TMath::Abs(gain)<1.0e-3)
            gain =1.;
    }

    if(TMath::Abs(gain) < 1.0e-3)
            gain = 1.0;

    return gain;

}


