#include "GParticleReconstruction.h"
#include "AdlarsonPhysics.h"
#include "GTrue.h"
#include <cmath>
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TRandom.h"

#include <APLCON.hpp>

std::default_random_engine AdlarsonPhysics::FitParticle::generator;

AdlarsonPhysics::AdlarsonPhysics():
    time_TOF("time_TOF", "TAPS time_TOF time" ,23000, 0, 23000, 1600, -40., 40.),
    time_clusters_TAPS("time_clusters_TAPS", "TAPS cluster time" ,200, -50., 50., 450, 0, 450),
    time_clusters_CB("time_clusters_CB", "CB cluster time", 200, -50., 50., 720, 0, 720),
    time_clusters_CBavg_CBdiff("time_clusters_CBavg_CBdiff","time_clusters_CBavg_CBdiff",400, -20., 20., 720, 0, 720),
    time_clusters_CBavg_TAPSdiff("time_clusters_CBavg_TAPSdiff", "CB avg time - TAPS cluster time" ,800, -50., 50., 450, 0, 450),
    time_clusters_TAPS_TAPSdiff("time_clusters_TAPS_TAPSdiff", "TAPS - TAPS cluster time (for > 1 TAPS clusters)" ,800, -50., 50., 450, 0, 450),
    time_clusters_CBavg_EPT("time_clusters_CBavg_EPT", "CB avg time - Tagger time" ,800, -50., 50., 50, 0, 50),
    time_nr_AllClusters("time_nr_AllClusters", "nr of all detected clusters",15, 0, 15),
    time_nr_ClustersinTime("time_nr_ClustersinTime", "nr of detected clusters in time window and above Cluster Energy cut",15, 0, 15),
    time_nr_FinalClusterSel("time_nr_FinalClusterSel", "nr of clusters fulfilling all cuts",15, 0, 15),
    six_time_TaggedTime("six_time_TaggedTime", "Tagger time for the events with 7 clusters", 2000, -500, 500),
    detnr(6),
    CB_region(4),
    photons_rec(nPhotons_six),
    photons_fit(nPhotons_six),
    photons_fit_final(nPhotons_six),
    kinfit4g("4greactions"),
    kinfit("etaprime"),
    kinfit_final("etaprime6g_final"),
    kinfit3pi("3pi0hyp"),
    kinfiteta2pi("eta2pihyp"),
    kinfit10g("etaprime10g"),
    kinfit10g_eta2pi("etaprime10g_eta2pihyp"),
    Photons_four(nPhotons_four),
    Photons_six(nPhotons_six),
    Photons_six_3pi(nPhotons_six),
    Photons_six_eta2pi(nPhotons_six),
    Photons_ten(nPhotons_ten),
    Photons_ten_eta2pi(nPhotons_ten)
{

// TRUE OBSERVABLES
// Beam energy
    true_BeamE                  = new GH1("true_BeamE", "True Beam Energy", 250, 1.400, 1.650);
    tag_BeamE                   = new GH1("tag_BeamE", "Tagged Beam Energy", 250, 1400, 1650);
// Phase space final state particles
    true_th_p_v_th_etapr_CM     = new GHistBGSub2("true_th_p_v_th_etapr_CM", "#theta_{p} vs #theta_{#eta^{'}}", 360, 0., 180, 360, 0., 180.);
    true_th_v_E_p               = new GHistBGSub2("true_th_v_E_p", "True E_{p} vs #theta_{p}", 100, 0., 0.6, 100, 0., 25.);
    true_th_v_E_eta_g           = new GHistBGSub2("true_th_v_E_eta_g", "E_{#gamma, #eta} vs #theta_{#gamma, #eta}", 1000, 0, 1000, 180, 0, 180);
    true_th_v_E_eta_6g           = new GHistBGSub2("true_th_v_E_eta_6g", "E_{#gamma, #eta --> 6#gamma} vs #theta_{#gamma, #eta}", 1000, 0, 1000, 180, 0, 180);
    true_th_v_E_pi0_g           = new GHistBGSub2("true_th_v_E_pi0_g", "E_{#gamma, #pi^{0}} vs #theta_{#gamma, #pi^{0}}", 600, 0, 600, 180, 0, 180);
//  Kinfit tests
    true_six_dth_vs_th_p        = new GHistBGSub2("true_six_dth_vs_th_p", "proton; #theta_{p,rec}; #theta_{rec}-#theta_{true} (^{o})", 200, 0, 25, 80, -5., 5.);
    true_six_z_v_Ncl            = new GHistBGSub2("true_six_z_v_Ncl", "Nr rec clusters vs z_{true}; z_{true} (cm); z_{true} - z_{fit} (cm)", 100, -10., 10., 15, 0, 15);
    true_six_fit_dz_v_z         = new GHistBGSub2("six_fit_dz_v_z", "#Delta z vs z_{true}; z_{true} (cm); z_{true} - z_{fit} (cm)", 100, -10., 10., 50, -10., 10.);

// Physics result
    true_DP                     = new TH2D("true_DP", "True Dalitz Plot distribution", 600, -1.5, 1.5, 600, -1.5, 1.5);
    true_phy_DP                 = new TH1D("true_phy_DP", "True Dalitz Plot distribution as bin nr 6#gamma", 800, 0, 800);
    true_M_pi1pi2_e2p           = new TH1D("true_M_pi1pi2_e2p", "True M_{#pi#pi,true}^{2}", 100 , 0.0, 200.);
    true_M_etapi_e2p            = new TH1D("true_M_etapi_e2p", "True M_{#eta#pi,true}^{2}", 350, 0.0, 700.);
// In 6g analysis
    true_six_phy_dMpipi_v_Mpipi = new GHistBGSub2("true_six_phy_dMpipi_v_Mpipi", "fitted - true value M_{#pi#pi,fit}^{2}", 200, 0.0, 200, 200, -100, 100);
    true_six_phy_dX_v_DPbin     = new GHistBGSub2("true_six_phy_dX_v_DPbin", "X_{fit} - X_{true} vd DP bin nr", 800, 0, 800, 200, -2.0, 2.0);
    true_six_phy_dY_v_DPbin     = new GHistBGSub2("true_six_phy_dY_v_DPbin", "Y_{fit} - Y_{true} vd DP bin nr", 800, 0, 800, 200, -2.0, 2.0);

    true_ten_phy_dMpipi_v_Mpipi = new GHistBGSub2("true_ten_phy_dMpipi_v_Mpipi", " ten #gamma: fitted - true value M_{#pi#pi,fit}^{2}", 200, 0.0, 200, 200, -100, 100);
    true_ten_phy_dX_v_DPbin     = new GHistBGSub2("true_ten_phy_dX_v_DPbin", "ten #gamma: X_{fit} - X_{true} vd DP bin nr", 800, 0, 800, 200, -2.0, 2.0);
    true_ten_phy_dY_v_DPbin     = new GHistBGSub2("true_ten_phy_dY_v_DPbin", "ten #gamma: Y_{fit} - Y_{true} vd DP bin nr", 800, 0, 800, 200, -2.0, 2.0);


// RECONSTRUCTED OBSERVABLES
    // Correlation CB and TAPS
    fi_diff_TAPSCB              = new GH1("fi_diff_TAPSCB", "For clusters close to CB-TAPS border ; #Delta#phi TAPS - CB (^{o}}); Events / 5^{o} ", 80, -200., 200.);
    fi_th_diff_TAPSCB           = new GHistBGSub2("fi_th_diff_TAPSCB", "#Delta#phi v #Delta#theta TAPS - CB", 200, -200., 200., 80, -20., 20.);
    fi_TAPSvsCB                 = new GHistBGSub2("fi_TAPSvsCB", "#phi TAPS vs CB", 100, -200.,200.,100, -200.,200.);
// Rec. TAPS - proton analysis
    p_E_v_dE_all                = new GHistBGSub2("p_E_v_dE_all", "All E vs dE TAPS", 1200, 0., 600., 200, 0.0, 10.);
    p_E_v_dE_cc                 = new GHistBGSub2("p_E_v_dE_cc", "All E vs dE TAPS chance coinc", 1200, 0., 600., 200, 0.0, 10.);
    p_E_v_dE_pr                 = new GHistBGSub2("p_E_v_dE_pr", "Best proton cand E vs dE TAPS", 1200, 0., 600, 200, 0.0, 10.);
    p_E_v_TOF                   = new GHistBGSub2("p_E_v_TOF", "Energy TAPS vs TOF/ns", 800, -20., 20., 800, 0., 800.);
    p_E_v_TOF_All               = new GHistBGSub2("p_E_v_TOF_All", "Energy TAPS vs TOF ", 800, -20., 20., 800, 0., 800.);
    p_E_v_TOF_TAPS_1cl          = new GHistBGSub2("p_E_v_TOF_TAPS_1cl", "Energy TAPS vs TOF with 1 TAPS cl ", 800, -20., 20., 800, 0., 800.);
    p_E_v_TOF_CB_All            = new GHistBGSub2("p_E_v_TOF_CB_All", "Energy TAPS vs TOF (TAPS - CB) ", 800, -20., 20., 800, 0., 800.);
    p_E_v_TOF_TAPS_All          = new GHistBGSub2("p_E_v_TOF_TAPS_All", "Energy TAPS vs TOF (TAPS - TAPS) when > 1 TAPS cluster ", 800, -20., 20., 800, 0., 800.);
    p_E_v_TOF_CB_2PrID          = new GHistBGSub2("p_E_v_TOF_CB_2PrID", "Energy TAPS vs TOF (TAPS - CB) when more than one proton cand ", 800, -20., 20., 800, 0., 800.);
    p_E_v_TOF_CB_best           = new GHistBGSub2("p_E_v_TOF_CB_best", "Energy TAPS vs TOF (TAPS - CB) best out of multiple proton cand ", 800, -20., 20., 800, 0., 800.);
    p_E_v_TOF_TAPS_2PrID        = new GHistBGSub2("p_E_v_TOF_TAPS_2PrID", "Energy TAPS vs TOF (TAPS - TAPS) if more than 1 cluster in TAPS ", 800, -20., 20., 800, 0., 800.);
    p_E_v_TOF_All_wVeto         = new GHistBGSub2("p_E_v_TOF_All_wVeto", "Energy TAPS vs TOF w VETO > 1MeV", 800, -20., 20., 800, 0., 800.);

    p_nr                        = new GHistBGSub("p_nr", "nr of protons", 5,  0, 5);
    p_th_v_E                    = new GHistBGSub2("p_th_v_E", "Rec E_{p} vs #theta_{p}", 120, 0., 600, 50, 0., 25.);
    p_MM                        = new GHistBGSub("MM_p", "Missing Mass calculated for proton", 300, 800., 1100.);

    CB_EnergySum                = new GHistBGSub("CB_EnergySum", "Crystal Ball Energy Sum", 400, 0., 2000.);
    CB_EnergySum_2pi0           = new GHistBGSub("CB_EnergySum_2pi0", "Crystal Ball Energy Sum for 2pi0", 400, 0., 2000.);
    CB_EnergySum_etapi0         = new GHistBGSub("CB_EnergySum_etapi0", "Crystal Ball Energy Sum for etapi0", 400, 0., 2000.);
    CB_EnergySum_3pi0           = new GHistBGSub("CB_EnergySum_3pi0", "Crystal Ball Energy Sum for 3pi0", 400, 0., 2000.);
    CB_EnergySum_etapr           = new GHistBGSub("CB_EnergySum_etapr", "Crystal Ball Energy Sum for eta prime", 400, 0., 2000.);

    IMgg_v_det_2pi0_CB          =   new GHistBGSub2("IMgg_v_det_2pi0_CB", "IM(gg) 2#pi^{0}, CB", 500, 0, 1000, 720, 0, 720);
    IMgg_v_det_etapi0_pi0_CB    =   new GHistBGSub2("IMgg_v_det_etapi0_pi0_CB", "IM(gg) #eta#pi^{0}, CB", 400, 0, 400, 720, 0, 720);
    IMgg_v_det_etapi0_eta_CB    =   new GHistBGSub2("IMgg_v_det_etapi0_eta_CB", "IM(gg) #eta#pi^{0}, CB", 400, 400, 800, 720, 0, 720);
    IMgg_v_det_3pi0_CB          =   new GHistBGSub2("IMgg_v_det_3pi0_CB", "IM(gg) 3#pi^{0}, CB, CB", 500, 0, 1000, 720, 0, 720);
    IMgg_v_det_3pi0_CB_fit      =   new GHistBGSub2("IMgg_v_det_3pi0_CB_fit", "IM(gg) 3#pi^{0} fit, CB, CB", 500, 0, 1000, 720, 0, 720);
    IMgg_v_det_2pi0_TAPS        =   new GHistBGSub2("IMgg_v_det_2pi0_TAPS", "IM(gg) 2#pi^{0}, TAPS", 500, 0, 1000, 440, 0, 440);
    IMgg_v_det_etapi0_TAPS      =   new GHistBGSub2("IMgg_v_det_etapi0_TAPS", "IM(gg) #eta#pi^{0}, TAPS", 500, 0, 1000, 440, 0, 440);
    IMgg_v_det_3pi0_TAPS        =   new GHistBGSub2("IMgg_v_det_3pi0_TAPS", "IM(gg) 3#pi^{0}, TAPS", 500, 0, 1000, 440, 0, 440);
    IMgg_v_det_3pi0_TAPS_fit    =   new GHistBGSub2("IMgg_v_det_3pi0_TAPS_fit", "IM(gg) 3#pi^{0} fit, TAPS", 500, 0, 1000, 440, 0, 440);

// Rec. Photons

    ten_rec_IM                  = new GH1("ten_rec_IM", " rec IM(10#gamma)", 240,  200, 1400);
    ten_rec_IM_v_MMp            = new GHistBGSub2("ten_rec_IM_v_MMp", "MM(p) vs IM(10#gamma)", 300,800., 1100., 240, 200., 1400.);

//  Related variables 4g analysis

    four_rec_IM                 = new GH1("four_rec_IM", "rec. IM(4#gamma)", 280, 0, 1400);
    four_fit_chi2               = new GH1("four_fit_chi2", "#chi^{2} kinfit for 4g final state", 500, 0, 500.);
    four_fit_pdf                = new GH1("four_fit_pdf", "#pdf kinfit 4g final state", 100, 0, 1.);

    four_fit_PDF_etapi_v_2pi    = new GHistBGSub2("four_fit_PDF_etapi_v_2pi", "PDF_etapi_v_2pi 4#gamma", 100, 0., 1., 100, 0., 1.);
    four_fit_IM                 = new GH1("four_fit_IM", "IM_4#gamma fitted", 260, 0, 1400);

//    four_fit_Pulls_g_E_vs_E_CB       = new GHistBGSub2("four_fit_Pulls_g_E_vs_E_CB", "4g Pulls #gamma E vs E CB", 100, 0, 1000, 50, -5., 5.);
//    four_fit_Pulls_g_th_vs_th_CB     = new GHistBGSub2("four_fit_Pulls_g_th_vs_th_CB", "4g Pulls #gamma #theta vs #theta CB", 90, 0, 180, 50, -5., 5.);
//    four_fit_Pulls_g_phi_vs_th_CB    = new GHistBGSub2("four_fit_Pulls_g_phi_vs_th_CB", "4g Pulls #gamma #phi vs #theta CB", 90, 0, 180., 50, -5., 5.);
//    four_fit_Pulls_g_E_vs_E_TAPS       = new GHistBGSub2("four_fit_Pulls_g_E_vs_E_TAPS", "4g Pulls #gamma E vs E TAPS", 100, 0, 1000, 50, -5., 5.);
//    four_fit_Pulls_g_th_vs_th_TAPS     = new GHistBGSub2("four_fit_Pulls_g_th_vs_th_TAPS", "4g Pulls #gamma #theta vs #theta TAPS", 25, 0, 50, 50, -5., 5.);
//    four_fit_Pulls_g_phi_vs_th_TAPS    = new GHistBGSub2("four_fit_Pulls_g_phi_vs_th_TAPS", "4g Pulls #gamma #phi vs #theta TAPS", 25, 0, 50, 50, -5., 5.);
    four_fit_Pulls_g_E_vs_eth_CB        = new GHistBGSub2("four_fit_Pulls_g_E_vs_eth_CB", "4g Pulls #gamma E vs E, #theta CB", 7000, 0, 7000, 50, -5., 5.);
    four_fit_Pulls_g_th_vs_eth_CB       = new GHistBGSub2("four_fit_Pulls_g_th_vs_eth_CB", "4g Pulls #gamma #theta vs E, #theta CB", 7000, 0, 7000, 50, -5., 5.);
    four_fit_Pulls_g_phi_vs_eth_CB      = new GHistBGSub2("four_fit_Pulls_g_phi_vs_eth_CB", "4gPulls #gamma #phi vs E, #theta CB", 7000, 0, 7000, 50, -5., 5.);

    four_fit_best_2pi_IM_v_E        = new GHistBGSub2("four_fit_best_2pi_IM_v_E", "E_{#pi^{0}} vs M_{#gamma#gamma} for 2#pi^{0}", 100, 0., 1000., 100, 0., 400.);
    four_fit_best_etapi_IM_v_E      = new GHistBGSub2("four_fit_best_etapi_IM_v_E", "E_{#pi^{0}} vs M_{#gamma#gamma} for #eta#pi^{0}", 100, 0., 1000., 250, 0., 1000.);

    // 2pi0 inside phase space
    four_fit_mgg_v_eth             = new GHistBGSub2("four_fit_mgg_v_eth", "m_{#gamma#gamma} vs E, #theta", 27000, 0, 27000, 200, 0., 1000.);
    four_fit_mgg_v_eth_BaF2        = new GHistBGSub2("four_fit_mgg_v_eth_BaF2", "m_{#gamma#gamma} vE, #theta BaF2", 4000, 0, 4000, 200., 0, 1000.);
    four_fit_mgg_v_eth_PbWO4       = new GHistBGSub2("four_fit_mgg_v_eth_PbWO4", "m_{#gamma#gamma} v E, #theta PbWO4", 1000, 0, 1000, 200., 0, 1000.);

    // etapi0 eta-->gg cand inside phase space
    four_fit_m_eta_gg_v_eth        = new GHistBGSub2("four_fit_m_eta_gg_v_eth", "m_{#gamma#gamma} vs E, #theta #eta-->#gamma#gamma" , 27000, 0, 27000, 200, 0., 1000.);
    four_fit_m_etagg_v_eth_BaF2    = new GHistBGSub2("four_fit_m_etagg_v_eth_BaF2", "m_{#gamma#gamma} vE, #theta BaF2 #eta-->#gamma#gamma", 4000, 0, 4000, 200., 0, 1000.);
    four_fit_m_etagg_v_eth_PbWO4   = new GHistBGSub2("four_fit_m_etagg_v_eth_PbWO4", "m_{#gamma#gamma} v E, #theta PbWO4 #eta-->#gamma#gamma", 1000, 0, 1000, 200., 0, 1000.);

//    four_fit_mgg_v_eth_2             = new GHistBGSub2("four_fit_mgg_v_eth_2", "m_{#gamma#gamma} vs E, #theta inside cut", 27000, 0, 27000, 200, 0., 1000.);
//    four_fit_mgg_v_eth_BaF2_2        = new GHistBGSub2("four_fit_mgg_v_eth_BaF2_2", "m_{#gamma#gamma} vE, #theta BaF2 inside cut", 4000, 0, 4000, 200., 0, 1000.);
//    four_fit_mgg_v_eth_PbWO4_2       = new GHistBGSub2("four_fit_mgg_v_eth_PbWO4_2", "m_{#gamma#gamma} v E, #theta PbWO4 inside cut", 1000, 0, 1000, 200., 0, 1000.);

    // etapi0 eta-->gg cand inside phase space
//    four_fit_m_eta_gg_v_eth_2        = new GHistBGSub2("four_fit_m_eta_gg_v_eth_2", "m_{#gamma#gamma} vs E, #theta #eta-->#gamma#gamma inside cut" , 27000, 0, 27000, 200, 0., 1000.);
//    four_fit_m_etagg_v_eth_BaF2_2    = new GHistBGSub2("four_fit_m_etagg_v_eth_BaF2_2", "m_{#gamma#gamma} vE, #theta BaF2 #eta-->#gamma#gamma inside cut", 4000, 0, 4000, 200., 0, 1000.);
//    four_fit_m_etagg_v_eth_PbWO4_2   = new GHistBGSub2("four_fit_m_etagg_v_eth_PbWO4_2", "m_{#gamma#gamma} v E, #theta PbWO4 #eta-->#gamma#gamma inside cut", 1000, 0, 1000, 200., 0, 1000.);

    four_fit_best_2pi0_pi_E_v_th      = new GHistBGSub2("four_fit_best_2pi0_pi_E_v_th", "E_{#pi^{0}, #gamma} vs #Theta_{#pi^{0}, #gamma} 2#pi^{0}", 600, 0., 1200., 200, 0., 200.);
    four_fit_best_etapi0_pi_E_v_th    = new GHistBGSub2("four_fit_best_etapi0_pi_E_v_th", "E_{#pi^{0}, #gamma} vs #Theta_{#pi^{0}, #gamma} #eta#pi^{0}", 600, 0., 1200., 200, 0., 200.);

    // Kinfit related variables 6g

//    test_six_rec_IM             = new GH1("test_six_rec_IM", " rec. IM(6#gamma)", 280, 0., 1400.);
//    test_six_fit_chi2           = new GH1("test_six_fit_chi2", "#chi^{2} kinfit for 6#gamma", 500, 0, 500.);
//    test_six_fit_pdf            = new GH1("test_six_fit_pdf", "pdf kinfit for 6#gamma", 100, 0, 1.);
//    test_six_fit_Pulls          = new GHistBGSub2("test_six_fit_Pulls", "Pulls 6#gamma", 50, -5., 5., 25, 0, 25);

    six_rec_IM                  = new GH1("six_rec_IM", " rec. IM(6#gamma)", 280, 0., 1400.);
    six_rec_IM_v_MMp            = new GHistBGSub2("six_rec_IM_v_MMp", "MM(p) vs IM(6#gamma)",240, 200., 1200., 240, 200., 1400.);

    time_clusters_CB_3pi0       = new GHistBGSub2("time_clusters_CB_3pi0", "CB cluster time for 3pi0 events", 200, -50., 50., 720, 0, 720),

    six_rec_EvTh_5g             = new GHistBGSub2("six_rec_EvTh_5g", "E_{#gamma} vs #Theta_{#gamma} rec 5g cand", 1000, 0., 1000., 200, 0., 200.);
    six_rec_EvTh_6g             = new GHistBGSub2("six_rec_EvTh_6g", "E_{#gamma} vs #Theta_{#gamma} rec 6g cand", 1000, 0., 1000., 200, 0., 200.);
    six_rec_EvTh_7g             = new GHistBGSub2("six_rec_EvTh_7g", "E_{#gamma} vs #Theta_{#gamma} rec 7g cand", 1000, 0., 1000., 200, 0., 200.);

    six_rec_EvTh_6g_from_5g_sample = new GHistBGSub2("six_rec_EvTh_6g_from_5g_sample", "E_{#gamma} vs #Theta_{#gamma} rec 5g cand", 1000, 0., 1000., 200, 0., 200.);

    six_fit_chi2                = new GH1("six_fit_chi2", "#chi^{2} kinfit for 6#gamma", 100, 0, 100.);
    six_fit_pdf                 = new GH1("six_fit_pdf", "pdf kinfit for 6#gamma", 100, 0, 1.);
    six_fit_etaprfinal_pdf      = new GH1("six_fit_etaprfinal_pdf", "pdf kinfit with eta and 2pi mass enforced", 100, 0, 1.);
    six_fit_Pulls               = new GHistBGSub2("six_fit_Pulls", "Pulls 6#gamma", 50, -5., 5., 25, 0, 25);

    proton_fit_e_v_th           = new GHistBGSub2("proton_fit_e_v_th", "proton E vs #theta", 200, 0., 1000., 50, 0, 25);
    proton_fit_e_v_th_final     = new GHistBGSub2("proton_fit_e_v_th_final", "proton E vs #theta final ev sample", 200, 0., 1000., 50, 0, 25);

    six_fit_Pulls_g_E_vs_E_CB       = new GHistBGSub2("six_fit_Pulls_g_E_vs_E_CB", "Pulls #gamma E vs E CB", 100, 0, 1000, 50, -5., 5.);
    six_fit_Pulls_g_th_vs_th_CB     = new GHistBGSub2("six_fit_Pulls_g_th_vs_th_CB", "Pulls #gamma #theta vs #theta CB", 90, 0, 180, 50, -5., 5.);
    six_fit_Pulls_g_phi_vs_th_CB    = new GHistBGSub2("six_fit_Pulls_g_phi_vs_th_CB", "Pulls #gamma #phi vs #theta CB", 90, 0, 180., 50, -5., 5.);
    six_fit_Pulls_g_E_vs_E_TAPS       = new GHistBGSub2("six_fit_Pulls_g_E_vs_E_TAPS", "Pulls #gamma E vs E TAPS", 100, 0, 1000, 50, -5., 5.);
    six_fit_Pulls_g_th_vs_th_TAPS     = new GHistBGSub2("six_fit_Pulls_g_th_vs_th_TAPS", "Pulls #gamma #theta vs #theta TAPS", 50, 0, 50, 50, -5., 5.);
    six_fit_Pulls_g_phi_vs_th_TAPS    = new GHistBGSub2("six_fit_Pulls_g_phi_vs_th_TAPS", "Pulls #gamma #phi vs #theta TAPS", 50, 0, 50, 50, -5., 5.);

    six_fit_Pulls_g_E_vs_eth        = new GHistBGSub2("six_fit_Pulls_g_E_vs_eth", "Pulls #gamma E vs E, #theta CB", 7000, 0, 7000, 50, -5., 5.);
    six_fit_Pulls_g_th_vs_eth       = new GHistBGSub2("six_fit_Pulls_g_th_vs_eth", "Pulls #gamma #theta vs E, #theta CB", 7000, 0, 7000, 50, -5., 5.);
    six_fit_Pulls_g_phi_vs_eth      = new GHistBGSub2("six_fit_Pulls_g_phi_vs_eth", "Pulls #gamma #phi vs E, #theta CB", 7000, 0, 7000, 50, -5., 5.);

    six_fit_Pulls_p_th_vs_det_TAPS     = new GHistBGSub2("six_fit_Pulls_p_th_vs_det_TAPS", "Pulls proton #theta vs #theta TAPS", 440, 0, 440, 50, -5., 5.);
    six_fit_Pulls_p_phi_vs_det_TAPS    = new GHistBGSub2("six_fit_Pulls_p_phi_vs_det_TAPS", "Pulls proton #phi vs #theta TAPS", 440, 0, 440, 50, -5., 5.);

    six_fit_mgg_v_eth           = new GHistBGSub2("six_fit_mgg_v_eth", "m_{#gamma#gamma} vs E, #theta", 27000, 0, 27000, 200, 0., 1000.);
    six_fit_mgg_v_eth_BaF2      = new GHistBGSub2("six_fit_mgg_v_eth_BaF2", "m_{#gamma#gamma} vs E, #theta for BaF2", 8000, 0, 8000, 200, 0., 1000.);
    six_fit_mgg_v_eth_PbWO4     = new GHistBGSub2("six_fit_mgg_v_eth_PbWO4", "m_{#gamma#gamma} vs E, #theta for PbWO4", 2000, 0, 2000, 200, 0., 1000.);
    six_fit_mgg_v_e_BaF2        = new GHistBGSub2("six_fit_mgg_v_e_BaF2", "m_{#gamma#gamma} vs E for BaF2", 100, 0, 1000, 200, 0., 1000.);
    six_fit_mgg_v_e_PbWO4       = new GHistBGSub2("six_fit_mgg_v_e_PbWO4", "m_{#gamma#gamma} vs E for PbWO4", 100, 0, 1000, 200, 0., 1000.);

//    six_fit_mgg_v_eth_2           = new GHistBGSub2("six_fit_mgg_v_eth_2", "m_{#gamma#gamma} vs E , #theta in cut", 42000, 0, 42000, 200, 0., 1000.);
//    six_fit_mgg_v_eth_BaF2_2      = new GHistBGSub2("six_fit_mgg_v_eth_BaF2_2", "m_{#gamma#gamma} vs E, #theta for BaF2 in cut", 8000, 0, 8000, 200, 0., 1000.);
//    six_fit_mgg_v_eth_PbWO4_2     = new GHistBGSub2("six_fit_mgg_v_eth_PbWO4_2", "m_{#gamma#gamma} vs E, #theta for PbWO4 in cut", 2000, 0, 2000, 200, 0., 1000.);
    fit_mgg_pi_v_CB_2            = new GHistBGSub2("fit_mgg_pi_v_CB_2", "m_{#gamma#gamma} vs r vs detnr in cut CB", 18500, 0, 18500, 400, 0., 2.0);
    fit_mgg_eta_v_CB_2            = new GHistBGSub2("fit_mgg_eta_v_CB_2", "m_{#gamma#gamma} vs r vs detnr in cut CB", 18500, 0, 18500, 400, 0., 2.0);
    fit_mgg_pi_v_TAPS_2          = new GHistBGSub2("fit_mgg_pi_v_TAPS_2", "m_{#gamma#gamma} vs r vs detnr in cut TAPS", 13500, 0, 13500, 400, 0., 2.0);
    fit_mgg_eta_v_TAPS_2          = new GHistBGSub2("fit_mgg_eta_v_TAPS_2", "m_{#gamma#gamma} vs r vs detnr in cut TAPS", 13500, 0, 13500, 400, 0., 2.0);

    six_fit_IM                  = new GH1("six_fit_IM", "IM(6#gamma) after APLCON fit", 500, 400., 1400.);
    six_fit_IM_3pi              = new GH1("six_fit_IM_3pi", "IM(6#gamma) for 3#pi^{0} candidates", 500, 400., 1400.);
    six_fit_IM_eta2pi           = new GH1("six_fit_IM_eta2pi", "IM(6#gamma) for #eta2#pi^{0} candidates", 500, 400., 1400.);

    six_fit_which_place_best_3pi_cand = new GH1("six_fit_which_place_best_3pi_cand", "from rough chi2 test which fits kinfit hypoth best 3pi0", 10, 0, 10);
    six_fit_which_place_best_etapr_cand = new GH1("six_fit_which_place_best_etapr_cand", "from rough chi2 test which fits kinfit hypoth best eta prime", 10, 0, 10);;

    six_fit_PDF_eta2pi_v_3pi    = new GHistBGSub2("six_fit_PDF_eta2pi_v_3pi_2", "PDF_eta2pi_v_3pi 6#gamma using kinfit", 100, 0., 1., 100, 0., 1.);
    six_fit_PDF_eta2pi_v_Meta2pi= new GHistBGSub2("six_fit_PDF_eta2pi_v_Meta2pi", "PDF_eta2pi_v_3pi vs mass", 100, 0., 1., 200, 800., 1200.);
    six_fit_best_eta            = new GH1("six_fit_best_eta", "best #eta cand from comb", 500, 200, 700.);

    six_fit_EvTh_g                   = new GHistBGSub2("six_fit_EvTh_g", "E_{#gamma} vs #Theta_{#gamma}", 1000, 0., 1000., 200, 0., 200.);
    six_fit_EvTh_g_final             = new GHistBGSub2("six_fit_EvTh_g_final ev sample", "E_{#gamma} vs #Theta_{#gamma} final", 10000, 0., 1000., 200, 0., 200.);
    six_fit_EvTh_g_final_removed     = new GHistBGSub2("six_fit_EvTh_g_final_removed ev sample", "E_{#gamma} vs #Theta_{#gamma} final", 500, 0., 1000., 100, 0., 200.);

    six_fit_best_etapr_eta_E_v_th    = new GHistBGSub2("six_fit_best_etapr_eta_E_v_th", "E_{#eta, #gamma} vs #Theta_{#eta, #gamma}", 100, 0., 1000., 50, 0., 200.);
    six_fit_best_etapr_pi_E_v_th     = new GHistBGSub2("six_fit_best_etapr_pi_E_v_th", "E_{#pi^{0}, #gamma} vs #Theta_{#pi^{0}, #gamma}", 100, 0., 1000., 50, 0., 200.);
    six_fit_best_3pi0_pi_E_v_th      = new GHistBGSub2("six_fit_best_3pi0_pi_E_v_th", "E_{#pi^{0}, #gamma} vs #Theta_{#pi^{0}, #gamma} 3#pi^{0}", 500, 0., 1000., 200, 0., 200.);
    six_fit_best_3pi0_mgg_v_thth     = new GHistBGSub2("six_fit_best_3pi0_mgg_v_thth", "m_{#gamma#gamma} vs #theta #gamma1, #theta #gamma2 for 3#pi^{0}", 1600, 0, 1600, 200, 0., 1000.);
    six_fit_best_2pi                 = new GH1("six_fit_best_2pi", "best 2#pi^{0} cand from comb", 500, 0., 500.);

    six_fit_final_coplanarity        = new GH1("final_coplanarity", "#phi difference proton and #eta prime", 400, -20., 20.);
    six_fit_final_theta_diff         = new GH1("six_fit_final_theta_diff", "#theta difference for vector beam target - proton and #eta prime}", 400, -20., 20.);

    six_fit_eta_PDF_v_Metapr        = new GHistBGSub2("six_fit_eta_PDF_v_Metapr", "PDF when eta mass enforced vs IM, 6#gamma cand", 100, 0., 1., 100, 600., 1100.);

    six_phy_etapr_v_BeamE       = new GHistBGSub2("six_phy_etapr_v_BeamE", "IM(6#gamma) vs Beam Energy", 20, 1400, 1600, 500, 600., 1100.);
    six_phy_etapr_eta_v_BeamE   = new GHistBGSub2("six_phy_etapr_eta_v_BeamE", "IM(6#gamma) with enforced eta mass vs Beam Energy", 20, 1400, 1600, 500, 600., 1100.);

    six_phy_DP                  = new GHistBGSub2("six_phy_DP", "Rec fitted Dalitz Plot distribution vs bin nr 6#gamma", 800, 0, 800, 500, 600., 1100.);
    six_phy_M_pi1pi2_v_etapr    = new GHistBGSub2("six_phy_M_pi1pi2_v_etapr", "Rec M_{#pi#pi,fit}^{2} 6#gamma", 400 , 0.0, 200, 200, 600., 1100. );
    six_phy_M_etapi_v_etapr     = new GHistBGSub2("six_phy_M_etapi_v_etapr", "Rec M_{#eta#pi,fit}^{2} 6#gamma", 700 , 0.0, 700, 200, 600., 1100. );

    // to check the energy of the  eta pi0 system vs its inv mass
    six_fit_best_eta_IM_v_E     = new GHistBGSub2("six_fit_best_eta_IM_v_E", "E_{#eta} vs M_{#eta}", 100, 0., 1000., 150, 200., 800.);
    six_fit_best_2pi_IM_v_E     = new GHistBGSub2("six_fit_best_2pi_IM_v_E", "E_{#pi^{0}} vs M_{#gamma#gamma} for #eta2#pi^{0}", 100, 0., 1000., 100, 0., 400.);

    // to check the energy of the pi0 system vs its inv mass
    six_fit_best_3pi_IM_v_E     = new GHistBGSub2("six_fit_best_3pi_IM_v_E", "E_{#pi^{0}} vs M_{#gamma#gamma} for 3#pi^{0}", 100, 0., 1000., 100, 0., 400.);
    six_phy_3pi_IMpipi_v_IMppi  = new GHistBGSub2("six_phy_3pi_IMpipi_v_IMppi", "M_{#pi0#pi0} vs M_{p#pi0} for 3#pi^{0}", 100, 0., 1., 200, 1., 3.);


// Kinfit related variables 10g

    kfit_chi2_10g       = new GH1("kfit_chi2_10g", "#chi^{2} kinfit 10#gamma", 100, 0, 100.);
    kfit_pdf_10g        = new GH1("kfit_pdf_10g", "#pdf kinfit 10#gamma", 100, 0, 1.);
    kfit_Pulls_10g      = new GHistBGSub2("kfit_Pulls_10g", "Pulls 10#gamma", 50, -5., 5., 40, 0, 40);
    ten_fit_EvTh_g      = new GHistBGSub2("ten_fit_EvTh_g", "E_{#eta, #gamma} vs #Theta_{#eta, #gamma} 10 #gamma", 200, 0., 1000., 100, 0., 200.);
    IM10g_fit           = new GH1("IM10g_fit", "IM(10#gamma) after APLCON fit", 500, 400., 1400.);
    IM10g_fit_best_cand = new GH1("IM10g_fit_best_cand", "IM(10#gamma) after APLCON fit with eta and 2pi0 mass enforced", 500, 400., 1400.);

    ten_fit_PDF_eta2pi_v_eta6g   = new GHistBGSub2("ten_fit_PDF_eta2pi_v_eta6g", "ten_fit_PDF_eta2pi_v_eta6g 6#gamma", 100, 0., 1., 100, 0., 1.);

    ten_fit_PDF_eta2pi          = new GH1("ten_fit_PDF_eta2pi", "ten_fit_PDF_eta2pi 10#gamma", 100, 0., 1.);
    ten_fit_X_v_pdf_eta2pi      = new GHistBGSub2("ten_fit_X_v_pdf_eta2pi", "ten_fit_true_v_pdf_eta2pi 10#gamma", 100, 0., 1., 100, -5, 5.);
    ten_fit_Y_v_pdf_eta2pi      = new GHistBGSub2("ten_fit_Y_v_pdf_eta2pi", "ten_fit_true_v_pdf_eta2pi 10#gamma", 100, 0., 1., 100, -5., 5.);

    ten_fit_nfit_v_pdf_5pi0_v_eta2pi0 = new GHistBGSub2("ten_fit_nfit_v_pdf_5pi0_v_eta2pi0", "ten #gamma", 10000, 0, 10000, 50, 0, 50);
    ten_fit_dX_v_pdf_5pi0_v_eta2pi0 = new GHistBGSub2("ten_fit_dX_v_pdf_5pi0_v_eta2pi0", "ten #gamma", 10000, 0, 10000, 50, 0, 50);
    ten_fit_dY_v_pdf_5pi0_v_eta2pi0 = new GHistBGSub2("ten_fit_dY_v_pdf_5pi0_v_eta2pi0", "ten #gamma", 10000, 0, 10000, 50, 0, 50);

    // Physics results eta'

    cutFile             = new TFile("configfiles/cuts/TAPSEdEPA.root");
    cutProtonTAPS       = (TCutG*)cutFile->Get("CutProton");

    cutFile2            = new TFile("configfiles/cuts/TAPS_TOF_PA4.root");
    cutProtonETOF       = (TCutG*)cutFile2->Get("Hadron");

    // APLCON kinfit uncertainties

    g_unc               = new TFile("configfiles/APLCONunc/photon_uncertainties_vz.root");
    p_unc               = new TFile("configfiles/APLCONunc/proton_uncertainties_vz.root");

//  uncertainties obtained from particle gun
    g_e     = (TH2F*)g_unc->Get("E")->Clone();
    g_th    = (TH2F*)g_unc->Get("theta")->Clone();
    g_phi   = (TH2F*)g_unc->Get("phi")->Clone();

    p_TAPS_th = (TH2F*)p_unc->Get("p_theta")->Clone();
    p_TAPS_fi = (TH2F*)p_unc->Get("p_phi")->Clone();

    // corrections to kinfit based on Pull results for eta'

    unc_corr                   = new TFile("configfiles/corr/unc_corr.root");
    g_e_c1                     = (TH2F*)unc_corr->Get("Eg_unc_corr")->Clone();
    g_th_c1                    = (TH2F*)unc_corr->Get("thg_unc_corr")->Clone();
    p_th_c1                    = (TH1F*)unc_corr->Get("th_pr_unc_corr")->Clone();
    p_fi_c1                    = (TH1F*)unc_corr->Get("phi_pr_unc_corr")->Clone();

    fin_sel                    = new TFile("configfiles/cuts/final_state_selection.root");
    etapr_eta                  = (TH2F*)fin_sel->Get("etapr_eta_gg")->Clone();
    etapr_pi                   = (TH2F*)fin_sel->Get("etapr_pi0_gg")->Clone();
    dir3pi_gg                  = (TH2F*)fin_sel->Get("3pi0_pi0_gg")->Clone();

    thcorr_CB                   = new TFile("configfiles/corr/CB_th_corr.root");
    dthvth_CB                   = (TProfile*)thcorr_CB->Get("photon_dtheta_v_theta_CB_pfx")->Clone();

    Ecorr_CB                    = new TFile("configfiles/corr/CB_e_corr.root");
    EvdetCB                     = (TH2F*)Ecorr_CB->Get("g_peak_E_CB");

    Ecorr_TAPS                  = new TFile("configfiles/corr/TAPS_e_corr.root");
    EvdetTAPS                   = (TH2F*)Ecorr_TAPS->Get("g_peak_E_TAPS");

    thcorr_TAPS                 = new TFile("configfiles/corr/TAPS_th_corr.root");
    dthvth_TAPS                 = (TProfile*)thcorr_TAPS->Get("photon_dtheta_v_theta_TAPS_pfx")->Clone();

    Evth_g_sel                  = new TFile("configfiles/cuts/E_v_th_g_cut.root");
    sixg_cand                   = (TCutG*)Evth_g_sel->Get("CutProton")->Clone();

    Evth_7g_sel                  = new TFile("configfiles/cuts/E_vs_th_7g_cut.root");
    seveng_cand                   = (TCutG*)Evth_7g_sel->Get("CutProton")->Clone();

    eta_cand                    = new TFile("configfiles/cuts/eta_cand.root");
    eta_g_cand                  = (TCutG*)eta_cand->Get("CutProton")->Clone();
    pi0_cand                    = new TFile("configfiles/cuts/pi0_cand.root");
    pi0_g_cand                  = (TCutG*)pi0_cand->Get("CutProton")->Clone();

    PDF_cut_file                = new TFile("configfiles/cuts/PDF_cut.root");
    PDF_cut                   = (TCutG*)PDF_cut_file->Get("CUTG")->Clone();

    // cuts for 3pi phase space, 2pi phase space and etapi0 phase space
    threepi_Evth                = new TFile("configfiles/cuts/threepi_Evth.root");
    threepi_g_cand              = (TCutG*)threepi_Evth->Get("CutProton")->Clone();
    twopi_Evth                  = new TFile("configfiles/cuts/twopi_Evth.root");
    twopi_g_cand                = (TCutG*)twopi_Evth->Get("CutProton")->Clone();
    etapi_Evth                  = new TFile("configfiles/cuts/etapi_Evth.root");
    etapi_g_cand                = (TCutG*)etapi_Evth->Get("CutProton")->Clone();


    GHistBGSub::InitCuts(-12., 12., -30., -20.);
    GHistBGSub::AddRandCut(20., 30.);

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
    kinfit.AddConstraint("EnergyMomentumBalance", all_names, EnergyMomentumBalance);

    kinfiteta2pi.AddConstraint("EnergyMomentumBalance", all_names_eta2pi, EnergyMomentumBalance);
    kinfit3pi.AddConstraint("EnergyMomentumBalance", all_names_3pi, EnergyMomentumBalance);

    kinfit4g.AddConstraint("EnergyMomentumBalance", all_names4g, EnergyMomentumBalance);
    kinfit10g.AddConstraint("EnergyMomentumBalance", all_names10g, EnergyMomentumBalance);
    kinfit10g_eta2pi.AddConstraint("EnergyMomentumBalance", all_names10g_eta2pi, EnergyMomentumBalance);

    // Constraint: Invariant mass of nPhotons equals constant IM,
    // make lambda catch also this with [&] specification

    // Constraint: Coplanarity between eta' and recoil proton
//    auto CoplanarityConstraint = [] (const vector<vector<double>>& particles) -> double
//    {
//        TLorentzVector proton = FitParticle::Make(particles[1], MASS_PROTON);
//        TLorentzVector etap(0., 0., 0., 0.);
//        for(size_t i=2;i<particles.size();i++) {
//            etap += FitParticle::Make(particles[i],0.0);
//        }

//        return abs(etap.Phi() - proton.Phi())*TMath::RadToDeg() - 180.;
//    };

//    kinfit.AddConstraint("CoplanarityConstraint", all_names, CoplanarityConstraint);

//    kinfiteta2pi.AddConstraint("CoplanarityConstraint", all_names_eta2pi, CoplanarityConstraint);
//    kinfit3pi.AddConstraint("CoplanarityConstraint", all_names_3pi, CoplanarityConstraint);

//    auto RequireEtaprime = [&] (const vector< vector<double> >& photons) -> double
//    {
//        TLorentzVector sum(0,0,0,0);
//        for(int i = 0; i < 6; i++) {
//            sum += FitParticle::Make(photons[i], 0.0);
//        }
//        return sum.M() - MASS_ETAP;
//    };
    auto RequireIM = [&] (const vector< vector<double> >& photons) -> double
    {
        TLorentzVector sum(0,0,0,0);
        for(int i = 0; i < 2; i++) {
            sum += FitParticle::Make(photons[i], 0.0);
        }
        return sum.M() - MASS_ETA;
    };

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

    kinfit3pi.AddConstraint("RequireIMpi1", photon_names_3pi, RequireIMpi1); //added req that first pair has pi0 mass.
    kinfit3pi.AddConstraint("RequireIMpi2", photon_names_3pi, RequireIMpi2); //added req that second pair has pi0 mass.
    kinfit3pi.AddConstraint("RequireIMpi3", photon_names_3pi, RequireIMpi3); //added req that third pair has pi0 mass.

    kinfiteta2pi.AddConstraint("RequireIM", photon_names_eta2pi, RequireIM); //added req that eta' candidate has eta mass.
    kinfiteta2pi.AddConstraint("RequireIMpi2", photon_names_eta2pi, RequireIMpi2); //added req that second pair has pi0 mass.
    kinfiteta2pi.AddConstraint("RequireIMpi3", photon_names_eta2pi, RequireIMpi3); //added req that third pair has pi0 mass.

    kinfit10g_eta2pi.AddConstraint("RequireIM6g",photon_names10g_eta2pi, RequireIM6g); //added req that eta' candidate has eta mass.
    kinfit10g_eta2pi.AddConstraint("RequireIMpi4",photon_names10g_eta2pi, RequireIMpi4); //added req that second pair has pi0 mass.
    kinfit10g_eta2pi.AddConstraint("RequireIMpi5",photon_names10g_eta2pi, RequireIMpi5); //added req that third pair has pi0 mass.

//    auto VertexConstraint = [&] (vector< vector<double> >& args) -> vector<double>
//    {
//       TLorentzVector diff(0.0, 0.0, 0.0, 0.0);
//       const TLorentzVector target(0,0,0, MASS_PROTON);
//       // assume first particle is beam photon
//       diff = target + FitParticle::Make(args[0], 0.0 ); // beam photon

//       double R;
//       double X0 = 2.59;
//       double Ec = 20;
//       constexpr double R_CB = 25.4;
//       constexpr double R_TAPS = 145.7;
//       double theta_p;

//       const double v_z = args.back()[0];
//       args.resize(args.size()-1); // get rid of last element
//       // correct each final state particle theta angle,

//       for( UInt_t i = 1; i < 8; i++)
//       {
//           double theta = args[i][1]; // 2nd element theta
//           if( !(Is_CB_6g[i-1]) )
//           {
//               R = R_TAPS/cos(theta) + TMath::Log10(args[i][0]/Ec);;
//            //   R= R_CB;
////               theta_p = theta;
//               theta_p = std::atan2( R*sin(theta), R*cos(theta) - v_z);
//           }
//           else{
//               R = R_CB + X0*TMath::Log10(args[i][0]/Ec);
//               theta_p = std::atan2( R*sin(theta), R*cos(theta) - v_z);
//           }

////           double theta_p = std::atan2( R*sin(theta), R*cos(theta) - v_z);
//           args[i][1] = theta_p;                // now set the re-defined theta angle as the new theta
//           if(i == 1)
//               diff -= FitParticle::Make(args[i], MASS_PROTON );
//           else
//               diff -= FitParticle::Make(args[i], 0.0 );
//       }
//       return {diff.X(), diff.Y(), diff.Z(), diff.T()};

//    };
//    kinfit.AddUnmeasuredVariable("v_z"); // default value 0
//    kinfit.AddConstraint("VertexConstraint", all_names + std::vector<string>{"v_z"}, VertexConstraint);
//    kinfiteta2pi.AddUnmeasuredVariable("v_z"); // default value 0
//    kinfiteta2pi.AddConstraint("VertexConstraint", all_names_eta2pi + std::vector<string>{"v_z"}, VertexConstraint);
//    kinfit3pi.AddUnmeasuredVariable("v_z"); // default value 0
//    kinfit3pi.AddConstraint("VertexConstraint", all_names_3pi + std::vector<string>{"v_z"}, VertexConstraint);

    APLCON::Fit_Settings_t settings = kinfit.GetSettings();
    settings.MaxIterations = 15;

    APLCON::Fit_Settings_t settings_eta2pi = kinfiteta2pi.GetSettings();
    settings_eta2pi.MaxIterations = 15;
    APLCON::Fit_Settings_t settings_3pi = kinfit3pi.GetSettings();
    settings_3pi.MaxIterations = 15;


//    settings.DebugLevel = 5;
    kinfit.SetSettings(settings);
    kinfiteta2pi.SetSettings(settings_eta2pi);
    kinfit3pi.SetSettings(settings_3pi);

    cout.precision(3);
    APLCON::PrintFormatting::Width = 11;
}

AdlarsonPhysics::~AdlarsonPhysics()
{

}

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
    
//    makeoutput tree

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

//    tree.write
    
//    close and write physics tree
    
	return kTRUE;
}

void	AdlarsonPhysics::ProcessEvent()
{
    if(MC){
//       etapr_6gTrue.Start(*GetPluto(), *GetGeant());   // (pluto tree, n part in pluto per event)
//       TrueAnalysis_etapr6g();                     // obtains the true observables

//        etapr_10gTrue.Start(*GetPluto(), *GetGeant());   // (pluto tree, n part in pluto per event)
//        TrueAnalysis_etapr10g();
        RandomTime();
    }
    for ( Int_t tag = 0; tag < GetTagger()->GetNTagged(); tag++)
        tag_BeamE->Fill(GetTagger()->GetTaggedEnergy(tag),GetTagger()->GetTaggedTime(tag));

    theta_corr();
    Energy_corr();


    ClustersInTime.clear();
    ClustersInTime.resize(0);
    UInt_t InTime = 0;
    UInt_t AllTimes = 0;
    UInt_t TAPS_clusters = 0;
    UInt_t CB_clusters = 0;

    Double_t CB_avg_time = -1000.;
    Double_t CB_avg_time_tmp = 0.;
    Double_t CB_en_sum = 0;



    // First step: count clusters and add CB clusters inside time window. Obtain CB time average and Energy
    for( Int_t i = 0; i < GetTracks()->GetNTracks(); i++ ){
        AllTimes++; // counting nr of CB and TAPS clusters
        if(GetTracks()->HasCB(i)){
            time_clusters_CB.Fill( GetTracks()->GetTime(i), GetTracks()->GetCentralCrystal(i) );
            CB_en_sum += GetTracks()->GetClusterEnergy(i);
            if( TMath::Abs( GetTracks()->GetTime(i) ) < 30.0 ){                               
                    CB_avg_time_tmp += GetTracks()->GetClusterEnergy(i)*GetTracks()->GetTime(i);
            }
        }
    }

    time_nr_AllClusters.Fill(AllTimes);
    CB_avg_time = CB_avg_time_tmp/CB_en_sum;
    CB_EnergySum->Fill(CB_en_sum);
    Double_t dt;
    // Second step: count clusters and add CB clusters inside time window. Obtain CB time average and Energy
    for( Int_t i = 0; i < GetTracks()->GetNTracks(); i++ ){
        if(GetTracks()->HasCB(i)){
            dt = GetTracks()->GetTime(i) - CB_avg_time;
            time_clusters_CBavg_CBdiff.Fill(dt, GetTracks()->GetCentralCrystal(i) );
            if( TMath::Abs( dt ) < 4.0 ){
                if(!sixg_cand->IsInside(GetTracks()->GetClusterEnergy(i), GetTracks()->GetTheta(i))){
//                    CB_avg_time_tmp += GetTracks()->GetClusterEnergy(i)*GetTracks()->GetTime(i);
                    CB_clusters++;
                    InTime++;
                    ClustersInTime.push_back(i);
                }
            }
        }
    }


//    Double_t TOF = -100.0;
    Double_t TOF_CB = -100.;

    Double_t radnm;
    std::vector<int> TAPS_cl;
    TAPS_cl.resize(0);

 // Second step: count and add TAPS clusters in time window
    for( Int_t j = 0; j < GetTracks()->GetNTracks(); j++ ){
        if(GetTracks()->HasTAPS(j)){
            time_clusters_TAPS.Fill( GetTracks()->GetTime(j), GetTracks()->GetCentralCrystal(j) );
            radnm = 1.45/TMath::Cos( GetTracks()->GetThetaRad(j) );
            TOF_CB = (( CB_avg_time -(PbWO4[GetTracks()->GetCentralCrystal(j)])*GetTracks()->GetTime(j) )/radnm);

            if(MC){
                TOF_CB = -1*(( CB_avg_time - GetTracks()->GetTime(j) )/radnm);
            }

            time_clusters_CBavg_TAPSdiff.Fill(TOF_CB, GetTracks()->GetCentralCrystal(j));
            if( ( TOF_CB > -15.0) &&  ( TOF_CB < 2.0 ) ){

                if(!sixg_cand->IsInside(GetTracks()->GetClusterEnergy(j), GetTracks()->GetTheta(j))){
                    TAPS_clusters++;
                    InTime++;
                    ClustersInTime.push_back(j);
                    TAPS_cl.push_back(j);
                }
            }
        }

    }

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
       }
       else
           TOF_CB = (( CB_avg_time -(PbWO4[GetTracks()->GetCentralCrystal(q)])*GetTracks()->GetTime(q) )/radnm);


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
                TOF_CB = (( CB_avg_time -(PbWO4[GetTracks()->GetCentralCrystal(q)])*GetTracks()->GetTime(q) )/radnm);
            }
            p_E_v_TOF_CB_2PrID->Fill(TOF_CB, GetTracks()->GetClusterEnergy(q) );
            if(TOF_CB < time_pr_temp){
                pr_cand_temp = q;
                time_pr_temp = TOF_CB;
            }
        }
        nrprotons++;
        iprtrack = pr_cand_temp;
        p_E_v_TOF_CB_best->Fill(time_pr_temp, GetTracks()->GetClusterEnergy(iprtrack) );
    }

    FinalClusterSelection.resize(0);
    IgnoreTAPSCluster.resize(0);
    // fourth step- go through good CB cluster. If CB cluster is in the edge, check if there is a corresponding TAPS hit in phi. If inside 40 deg then merge clusters
    for(UInt_t i = 0; i < ClustersInTime.size() ; i++){
        UInt_t j = ClustersInTime[i];
        if(GetTracks()->HasCB(j)){
            if( ring3_or_ring4_CB[ GetTracks()->GetCentralCrystal(j) ] && ClustersInTime.size() ==8 ){
                for(UInt_t p = 0; p < TAPS_cl.size(); p++){
                    UInt_t q = TAPS_cl[p];
                    if(q != iprtrack){
                        Double_t dfi = GetTracks()->GetPhi(j) - GetTracks()->GetPhi(q);
                        if(TMath::Abs(dfi) < 40. ){
                            Double_t E_CB   = GetTracks()->GetClusterEnergy(j);
                            Double_t E_TAPS = GetTracks()->GetClusterEnergy(q);
                            tracks->SetClusterEnergy(j, E_CB + E_TAPS );
                            IgnoreTAPSCluster.push_back(q);
                        }
                    }
                }
            }
            if(!sixg_cand->IsInside(GetTracks()->GetClusterEnergy(j), GetTracks()->GetTheta(j))){
                FinalClusterSelection.push_back(j);
            }
        }
    }
    for(UInt_t p = 0; p < TAPS_cl.size(); p++){
        UInt_t q = TAPS_cl[p];

        if( q == iprtrack )
            FinalClusterSelection.push_back(q);
        else{
            Bool_t ignore_el = false;
            for(uint s = 0; s < IgnoreTAPSCluster.size(); s++)
                if(IgnoreTAPSCluster[s] == q)
                    ignore_el = true;
            if(!ignore_el)
                    FinalClusterSelection.push_back(q);
            }
        }

    std::vector<Int_t> fin_temp;
    fin_temp.resize(0);
    if(FinalClusterSelection.size() == 8){
        for(UInt_t p = 0; p < FinalClusterSelection.size(); p++){
            UInt_t q = FinalClusterSelection[p];
            if(q == iprtrack)
                    fin_temp.push_back(q);
            else if(!seveng_cand->IsInside(GetTracks()->GetClusterEnergy(q), GetTracks()->GetTheta(q)))
                fin_temp.push_back(q);

        }

        FinalClusterSelection.resize(0);
        FinalClusterSelection = fin_temp;

    }

    time_nr_FinalClusterSel.Fill(FinalClusterSelection.size());

    for (Int_t ii = 0; ii < GetTracks()->GetNTracks(); ii++)
        if( GetTracks()->HasTAPS(ii) )
            p_E_v_dE_all->Fill(GetTracks()->GetClusterEnergy(ii),GetTracks()->GetVetoEnergy(ii));

    if(FinalClusterSelection.size() == 6){
        for(UInt_t p = 0; p < FinalClusterSelection.size() ; p++){
            UInt_t q = FinalClusterSelection[p];
                if(q != iprtrack)
                    six_rec_EvTh_5g->Fill(GetTracks()->GetClusterEnergy(q),GetTracks()->GetTheta(q));
        }
//        std::vector<Int_t> fin_temp;
//        fin_temp.resize(0);
//        int n_six = 0;
//        for(UInt_t p = 0; p < GetTracks()->GetNTracks(); p++){
//            if( GetTracks()->HasCB(p) ){
//                dt = GetTracks()->GetTime(p) - CB_avg_time;
//                if( TMath::Abs( dt ) < 4.0 ){
//                    n_six++;
//                    fin_temp.push_back(p);
//                }
//            else if(GetTracks()->HasTAPS(p) && (p!= iprtrack) ){
//                    radnm = 1.45/TMath::Cos( GetTracks()->GetThetaRad(p) );
//                    TOF_CB = (( CB_avg_time -(PbWO4[GetTracks()->GetCentralCrystal(p)])*GetTracks()->GetTime(p) )/radnm);

//                    if(MC){
//                        TOF_CB = -1*(( CB_avg_time - GetTracks()->GetTime(p) )/radnm);
//                    }
//                    if( ( TOF_CB > -5.0) &&  ( TOF_CB < 2.0 ) ){
//                        n_six++;
//                        fin_temp.push_back(p);
//                    }

//            }
//        }
//        if(n_six == 6)
//            for(UInt_t p = 0; p < GetTracks()->GetNTracks(); p++){
//                six_rec_EvTh_6g_from_5g_sample->Fill(GetTracks()->GetClusterEnergy(p),GetTracks()->GetTheta(p));
//            }
//        }
    }
    if(FinalClusterSelection.size() == 7){
        for(UInt_t p = 0; p < FinalClusterSelection.size() ; p++){
            UInt_t q = FinalClusterSelection[p];
                if(q != iprtrack)
                    six_rec_EvTh_6g->Fill(GetTracks()->GetClusterEnergy(q),GetTracks()->GetTheta(q));
        }
    }
    if(FinalClusterSelection.size() == 8){
        for(UInt_t p = 0; p < FinalClusterSelection.size() ; p++){
            UInt_t q = FinalClusterSelection[p];
                if(q != iprtrack){
                    six_rec_EvTh_7g->Fill(GetTracks()->GetClusterEnergy(q),GetTracks()->GetTheta(q));
                }
        }
    }

    //    Kinfit_test();  // runs kinematical fit with true observables- for testing purposes

    MMp_vec.SetPxPyPzE(0., 0., 0., 0.);



//    nrprotons = 0;
//    iprtrack = 10000;
//    Double_t TOF = -100.0;
//    Double_t TOF_CB = -100.; // obtain the CB time from the crystal with highest energy (least impact from TimeWalk) and compare this time to the TAPS time
//    for(int tag = 0; tag < GetTagger()->GetNTagged(); tag++){
//        time_clusters_CBavg_EPT.Fill(GetTagger()->GetTaggedTime(tag) -CB_avg_time,GetTagger()->GetTaggedChannel(tag));
//        if( TMath::Abs(GetTagger()->GetTaggedTime(tag)) > 12) continue;
//        for(UInt_t i = 0; i < ClustersInTime.size() ; i++){
//            UInt_t j = ClustersInTime[i];
//            if( GetTracks()->HasTAPS(j) ){

//                Double_t radnm = 1.45/TMath::Cos( GetTracks()->GetThetaRad(j) );
//                int nbin = GetTagger()->GetTaggedChannel(tag) + GetTracks()->GetCentralCrystal(j)*50;

//                if(MC){
//                    TOF = -1*(( GetTagger()->GetTaggedTime(tag) - GetTracks()->GetTime(j) )/radnm) - TAPS_EPT_toff[nbin];
//                    TOF_CB = -1*(( CB_avg_time - GetTracks()->GetTime(j) )/radnm);
//                    Double_t dTOF = pRandoms->Gaus(0, 0.35);
//                    TOF += dTOF;
//                    TOF_CB += dTOF;

//                }
//                else{
//                    TOF = (( GetTagger()->GetTaggedTime(tag) -(PbWO4[GetTracks()->GetCentralCrystal(j)])*GetTracks()->GetTime(j) )/radnm) - TAPS_EPT_toff[nbin];
//                    TOF_CB = (( CB_avg_time -(PbWO4[GetTracks()->GetCentralCrystal(j)])*GetTracks()->GetTime(j) )/radnm);
//                    time_clusters_CBavg_TAPSdiff.Fill(TOF_CB, GetTracks()->GetCentralCrystal(j));

//                }

//                if(TMath::Abs(TOF) < 40.1){
//                    int nbin = GetTagger()->GetTaggedChannel(tag) + GetTracks()->GetCentralCrystal(j)*50;
//                    time_TOF.Fill(nbin, TOF);

//                }
//                p_E_v_TOF_All->Fill( TOF , GetTracks()->GetClusterEnergy(j));
//                p_E_v_TOF_CB_All->Fill(TOF_CB, GetTracks()->GetClusterEnergy(j));
//                if( cutProtonETOF->IsInside( TOF, GetTracks()->GetClusterEnergy(j)) ){
//                    p_E_v_TOF->Fill( TOF, GetTracks()->GetClusterEnergy(j));
//                    nrprotons++ ;
//                    if(nrprotons > 1){
//                        Double_t TOF_champion = CB_avg_time -((PbWO4[GetTracks()->GetCentralCrystal(iprtrack)])*GetTracks()->GetTime(iprtrack) )/radnm;
//                        Double_t TOF_rocky = CB_avg_time  - ((PbWO4[GetTracks()->GetCentralCrystal(j)])*GetTracks()->GetTime(j) )/radnm;

//                        p_E_v_TOF_CB_2PrID->Fill(TOF_champion, GetTracks()->GetClusterEnergy(iprtrack) );
//                        p_E_v_TOF_CB_2PrID->Fill(TOF_rocky, GetTracks()->GetClusterEnergy(j) );
//                        if(TOF_rocky < TOF_champion){
//                            iprtrack = j;
//                            p_E_v_TOF_CB_best->Fill(CB_avg_time -((PbWO4[GetTracks()->GetCentralCrystal(iprtrack)])*GetTracks()->GetTime(iprtrack) )/radnm, GetTracks()->GetClusterEnergy(iprtrack));
//                        }
//                    }
//                    else
//                        iprtrack = j;
//                }

//                if(  GetTracks()->GetVetoEnergy(j) > ( 4.5 -(4.5/180)*GetTracks()->GetClusterEnergy(j)) )
//                    if(GetTracks()->GetVetoEnergy(j) > 1.0)
//                        p_E_v_TOF_All_wVeto->Fill(TOF, GetTracks()->GetClusterEnergy(j));

//           }
//        }
//    }


    for (Int_t ii = 0; ii < GetTracks()->GetNTracks(); ii++)
        if( GetTracks()->HasTAPS(ii) )
            p_E_v_dE_all->Fill(GetTracks()->GetClusterEnergy(ii),GetTracks()->GetVetoEnergy(ii));

    p_nr->Fill(nrprotons);
    if( (iprtrack == 10000) || (nrprotons == 0) ) return;

    true_six_dth_vs_th_p->Fill(GetTracks()->GetTheta(iprtrack),GetTracks()->GetTheta(iprtrack) - etapr_6gTrue.GetTrueProtonLV().Theta()*TMath::RadToDeg());

    p_E_v_dE_pr->Fill(GetTracks()->GetClusterEnergy(iprtrack),GetTracks()->GetVetoEnergy(iprtrack));
    p_th_v_E->Fill(GetTracks()->GetClusterEnergy(iprtrack),GetTracks()->GetTheta(iprtrack));
    proton_vec = GetTracks()->GetVector(iprtrack, pdgDB->GetParticle("proton")->Mass()*1000);

    // Now construct missing mass calc for proton with tagger energies.
    for ( Int_t tag = 0; tag < GetTagger()->GetNTagged(); tag++){
        MMp_vec.SetPxPyPzE(-proton_vec.Px(), -proton_vec.Py(), GetTagger()->GetTaggedEnergy(tag) - proton_vec.Pz(), GetTagger()->GetTaggedEnergy(tag) + pdgDB->GetParticle("proton")->Mass()*1000 - proton_vec.E());
        p_MM->Fill(MMp_vec.M(), GetTagger()->GetTaggedTime(tag));
    }

    CB_TAPS_boundary();

    ClustersInTime.resize(0);
    ClustersInTime = FinalClusterSelection;


//    if( (FinalClusterSelection.size() == 5) && (nrprotons > 0) )
//        fourgAnalysis(iprtrack);

    if( FinalClusterSelection.size() == 7 && (nrprotons > 0) )
        sixgAnalysis( iprtrack );

//    if( FinalClusterSelection.size() == 11 && (nrprotons > 0) )
//        tengAnalysis(iprtrack);
}

void	AdlarsonPhysics::ProcessScalerRead()
{
//    hist_eta.ScalerReadCorrection(Double_t(scalers->GetScaler(0))/scalers->GetScaler(1));
//    hist_etap.ScalerReadCorrection(Double_t(scalers->GetScaler(0))/scalers->GetScaler(1));
}


Bool_t	AdlarsonPhysics::Init(const char* configFile)
{
   std::string         line;

//   std::ifstream file_etapr_leg("configfiles/data/etaprime_Legendrecoeff.txt");
//   if(file_etapr_leg){
//       while(std::getline(file_etapr_leg, line)){
//           std::string  buffer;
//           std::stringstream ss;
//           ss << line;
//           while(std::getline(ss, buffer, '\t')){
//                Legendre.push_back(buffer);
//           }
//       }



//       while()


//   }


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

   std::ifstream file_EPT_TAPS_off("configfiles/data/TAPS_EPT_offsets.txt");
   if(file_EPT_TAPS_off){
   std::getline(file_EPT_TAPS_off, line);
   std::string         buffer4;
   std::stringstream   ss3;
   ss3 << line;
   while (std::getline(file_EPT_TAPS_off, line))
   {
       std::string         buffer;
       std::stringstream   ss;
       ss << line;

       while (std::getline(ss, buffer, '\t')) {
            TAPS_EPT_toff.push_back(std::stod(buffer));
       }
   }
   file_EPT_TAPS_off.close();
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
//     Int_t ring3_4_out[30] = {24, 25, 27, 28, 39, 41, 42, 43, 72, 73, 74, 77, 80, 81, 82, 636, 638, 639, 643, 644, 646, 647, 676, 677, 678, 680, 690, 693, 694, 695};

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




      std::cout << " CB gain correction applied" << std::endl;
      std::cout << " TAPS gain correction applied" << std::endl;
//    std::cout << "No Init function specified for this class." << std::endl;

    return kTRUE;
}


void AdlarsonPhysics::TrueAnalysis_etapr6g(){

    true_BeamE->Fill( etapr_6gTrue.GetTrueBeamEnergy() );

    // calculate Physics

    etapr_true[0] = etapr_6gTrue.GetTrueEtaLV();
    etapr_true[1] = etapr_6gTrue.GetTrueNeutralPiLV(0);
    etapr_true[2] = etapr_6gTrue.GetTrueNeutralPiLV(1);
    
    // here calculate the boost
    TVector3 etaprime_rest  = -(etapr_true[0]+etapr_true[1]+etapr_true[2]).BoostVector();
    TLorentzVector etapr_rest  = etapr_6gTrue.GetTrueEtaPrimeLV();
    etapr_rest.Boost(etaprime_rest);

    Double_t x = TMath::Cos(etapr_rest.Theta());
    Double_t leg1, leg2, leg3, leg4, leg5;
    Double_t weight, w1, w2;

    leg1 = 1;
    leg2 = x;
    leg3 = (1/2)*( 3*TMath::Power(x,2) - 1 );
    leg4 = (1/2)*( 5*TMath::Power(x,3) - 3*x );
    leg5 = (1/8)*( 35*TMath::Power(x,4) - 30*TMath::Power(x,2) + 3);


    
    // here calculate the eta prime theta as function of beam energy
    
    
    DalitzPlot(etapr_true, Xtrue1, Xtrue2, Ytrue, DPnrTrue1, DPnrTrue2);
    true_DP->Fill( Xtrue1, Ytrue );
    true_DP->Fill( Xtrue2, Ytrue );
    true_phy_DP->Fill(DPnrTrue1);
    true_phy_DP->Fill(DPnrTrue2);

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

void AdlarsonPhysics::TrueAnalysis_etapr10g()
{
    true_BeamE->Fill( etapr_10gTrue.GetTrueBeamEnergy() );

    // calculate Physics

    etapr_true[0] = etapr_10gTrue.GetTrueEtaLV();
    etapr_true[1] = etapr_10gTrue.GetTrueNeutralPiLV(0);
    etapr_true[2] = etapr_10gTrue.GetTrueNeutralPiLV(1);
    DalitzPlot(etapr_true, Xtrue1, Xtrue2, Ytrue, DPnrTrue1, DPnrTrue2);
    true_DP->Fill( Xtrue1, Ytrue );
    true_DP->Fill( Xtrue2, Ytrue );
    true_phy_DP->Fill(DPnrTrue1);
    true_phy_DP->Fill(DPnrTrue2);

    m2pi0_metapi0(etapr_true, m_etapi01True, m_etapi02True, m_2pi0True);
    true_M_pi1pi2_e2p->Fill(m_2pi0True*1.0e3);
    true_M_etapi_e2p->Fill(m_etapi01True*1.0e3);
    true_M_etapi_e2p->Fill(m_etapi02True*1.0e3);

    // calculate True Energy vs Theta of final state particles

    true_th_v_E_p->Fill(etapr_10gTrue.GetTrueProtonLV().E() - MASS_PROTON/1000 ,etapr_10gTrue.GetTrueProtonLV().Theta()*TMath::RadToDeg());

    for(UInt_t i = 0; i < etapr_10gTrue.GetNgamma(); i++ )
    {
        if ( i > 3 ) // first four gammas come from 2pi0
            true_th_v_E_eta_6g->Fill(etapr_10gTrue.GetTrueGammaLV(i).E(),etapr_10gTrue.GetTrueGammaLV(i).Theta()*TMath::RadToDeg());
        else
            true_th_v_E_pi0_g->Fill(etapr_10gTrue.GetTrueGammaLV(i).E(),etapr_10gTrue.GetTrueGammaLV(i).Theta()*TMath::RadToDeg());
    }
    return;
}

void AdlarsonPhysics::fourgAnalysis(UInt_t ipr)
{
    for(int tag = 0; tag < GetTagger()->GetNTagged(); tag++){
        Double_t chi2_min = 1.0e6;
        Double_t prob_min = 1.0e6;
        std::vector<Int_t> set;     // particle conf: tagger, proton, photon 1, 2, 3, 4, 5, 6;
        std::vector<double> obs;
        std::vector<double> unc;

        detnr.resize(0);
        CB_region.resize(0);

        IM4g_vec.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
        IM4g_fit.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);

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

        four_fit_chi2->Fill( result4g.ChiSquare, GetTagger()->GetTaggedTime(tag));
        four_fit_pdf->Fill( result4g.Probability, GetTagger()->GetTaggedTime(tag));

        for(UInt_t i = 0; i < nPhotons_four; i++){
            IM4g_vec += photons_rec[i];
            photons_fit[i] = FitParticle::Make(Photons_four[i], 0);
            IM4g_fit += photons_fit[i];
        }

        four_rec_IM->Fill(IM4g_vec.M());
        four_fit_IM->Fill(IM4g_fit.M());

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
        four_fit_PDF_etapi_v_2pi->Fill(TMath::Prob(chi2min_etapi0,1), TMath::Prob(chi2min_2pi0,1), GetTagger()->GetTaggedTime(tag) );

        if( (chi2min_2pi0 < chi2min_etapi0) && ( TMath::Prob(chi2min_2pi0,1) > 0.01) ){ // 2pi0
            for(int j = 0; j < 4; j++){

                four_fit_best_2pi0_pi_E_v_th->Fill(photons_fit[j].E(), photons_fit[j].Theta()*TMath::RadToDeg(), GetTagger()->GetTaggedTime(tag));

                int imass = int(j/2);
                idet = perm4g[imin_2pi0][j];
                if(CB_region[perm4g[ imin_2pi0 ][ j ]]){
                    double En = (photons_rec[perm4g[ imin_2pi0 ][ j ]]).E();
                    double Th = (photons_rec[perm4g[ imin_2pi0 ][ j]]).Theta()*TMath::RadToDeg();
                    double mgg = (photons_rec[perm4g[ imin_2pi0 ][ imass*2 ]] + photons_rec[perm4g[imin_2pi0][imass*2+1]] ).M();
                    int nEN =  int(En/20.0);
                    int nTH =  75*int(Th/2.0);
                    int nBIN = nEN + nTH;

                    IMgg_v_det_2pi0_CB->Fill( (photons_rec[perm4g[ imin_2pi0 ][ imass*2 ]] + photons_rec[perm4g[imin_2pi0][imass*2+1]] ).M(), detnr[idet], GetTagger()->GetTaggedTime(tag));
                    four_fit_mgg_v_eth->Fill( nBIN, mgg, GetTagger()->GetTaggedTime(tag) );

//                    if(twopi_g_cand->IsInside( photons_rec[perm4g[ imin_2pi0 ][ j ]].E(), photons_rec[perm4g[ imin_2pi0 ][ j ]].Theta()*TMath::RadToDeg()))
//                        four_fit_mgg_v_eth_2->Fill(nBIN, mgg, GetTagger()->GetTaggedTime(tag));

                }
                else{
                    IMgg_v_det_2pi0_TAPS->Fill((photons_rec[perm4g[ imin_2pi0 ][ imass*2 ]] + photons_rec[perm4g[imin_2pi0][imass*2+1]] ).M(), detnr[idet], GetTagger()->GetTaggedTime(tag) );
                    double En = (photons_rec[perm4g[ imin_2pi0 ][ j ]]).E();
                    double Th = (photons_rec[perm4g[ imin_2pi0 ][ j]]).Theta()*TMath::RadToDeg();
                    double mgg = (photons_rec[perm4g[ imin_2pi0 ][ imass*2 ]] + photons_rec[perm4g[imin_2pi0][imass*2+1]] ).M();
                    int nEN =  int(En/10.0);
                    int nTH =  150*int(Th/1.0);
                    int nBIN2 = nEN + nTH;
                    if( PbWO4[detnr[idet]] == -1 ){
                        four_fit_mgg_v_eth_PbWO4->Fill( nBIN2, mgg , GetTagger()->GetTaggedTime(tag));
//                        if(twopi_g_cand->IsInside(photons_rec[perm4g[ imin_2pi0 ][ j ]].E(), photons_rec[perm4g[ imin_2pi0 ][ j ]].Theta()*TMath::RadToDeg()))
//                            four_fit_mgg_v_eth_PbWO4_2->Fill(nBIN2, mgg, GetTagger()->GetTaggedTime(tag));
                    }
                    else{
                        four_fit_mgg_v_eth_BaF2->Fill( nBIN2, mgg , GetTagger()->GetTaggedTime(tag));
//                        if(twopi_g_cand->IsInside(photons_rec[perm4g[ imin_2pi0 ][ j ]].E(), photons_rec[perm4g[ imin_2pi0 ][ j ]].Theta()*TMath::RadToDeg()))
//                            four_fit_mgg_v_eth_BaF2_2->Fill(nBIN2, mgg, GetTagger()->GetTaggedTime(tag));

                    }

                }
            }
        }
        else if( (chi2min_etapi0 < chi2min_2pi0  ) && ( TMath::Prob(chi2min_etapi0,1) > 0.01)){   //etapi0

            Int_t nEn, nTh, nbin;
            for(UInt_t iPull = 0; iPull < 4; iPull++)
            {
                 string s = Form("Photon4g%0i[0]", iPull+1); // Energy
                 string t = Form("Photon4g%0i[1]", iPull+1); // theta
                 string u = Form("Photon4g%0i[2]", iPull+1); // phi

                 if( GetTracks()->HasCB( set[iPull+2]) ){
                     if( etapi_g_cand->IsInside( photons_rec[iPull].E(), photons_rec[iPull].Theta()*TMath::RadToDeg())){
                         nEn =  int(photons_rec[iPull].E()/10.);
                         nTh =  150*int(photons_rec[iPull].Theta()*TMath::RadToDeg()/1.0);
                         nbin = nEn + nTh;
                         four_fit_Pulls_g_E_vs_eth_CB->Fill(nbin, result4g.Variables.at(s).Pull, GetTagger()->GetTaggedTime(tag));
                         four_fit_Pulls_g_th_vs_eth_CB->Fill(nbin, result4g.Variables.at(t).Pull, GetTagger()->GetTaggedTime(tag));
                         four_fit_Pulls_g_phi_vs_eth_CB->Fill(nbin, result4g.Variables.at(u).Pull, GetTagger()->GetTaggedTime(tag));
                     }
                 }
            }
            for(int j = 0; j < 4; j++){
                four_fit_best_etapi0_pi_E_v_th->Fill(photons_fit[j].E(), photons_fit[j].Theta()*TMath::RadToDeg(), GetTagger()->GetTaggedTime(tag));

                int imass = int(j/2);
                idet = perm4g[imin_etapi0][j];
                if(CB_region[perm4g[ imin_etapi0 ][ j ]]){
                    if(imass == 0){
                        IMgg_v_det_etapi0_eta_CB->Fill( (photons_rec[perm4g[ imin_etapi0 ][ imass*2 ]] + photons_rec[perm4g[imin_etapi0][imass*2+1]] ).M(), detnr[idet], GetTagger()->GetTaggedTime(tag) );
                        double En = (photons_rec[perm4g[ imin_etapi0 ][ j ]]).E();
                        double Th = (photons_rec[perm4g[ imin_etapi0 ][ j]]).Theta()*TMath::RadToDeg();
                        double mgg = (photons_rec[perm4g[ imin_etapi0 ][ imass*2 ]] + photons_rec[perm4g[imin_etapi0][imass*2+1]] ).M();
                        int nEN =  int(En/10.0);
                        int nTH =  150*int(Th/1.0);
                        int nBIN = nEN + nTH;
                        four_fit_m_eta_gg_v_eth->Fill( nBIN, mgg, GetTagger()->GetTaggedTime(tag) );
                        if(etapi_g_cand->IsInside(photons_rec[perm4g[ imin_etapi0 ][ j ]].E(), photons_rec[perm4g[ imin_etapi0 ][ j ]].Theta()*TMath::RadToDeg())){
//                            four_fit_m_eta_gg_v_eth_2->Fill(nBIN, mgg, GetTagger()->GetTaggedTime(tag));
                            int nEN =  int(En/40.0);
                            int nBIN3 = 25*detnr[idet] + nEN;
                            fit_mgg_eta_v_CB_2->Fill(nBIN3, MASS_ETA/mgg, GetTagger()->GetTaggedTime(tag));
                        }                        
                    }
                    else{  //pi0 cand
//                        double En = (photons_rec[perm4g[ imin_etapi0 ][ j ]]).E();
//                        double mgg = (photons_rec[perm4g[ imin_etapi0 ][ imass*2 ]] + photons_rec[perm4g[imin_etapi0][imass*2+1]] ).M();

                        IMgg_v_det_etapi0_pi0_CB->Fill( (photons_rec[perm4g[ imin_etapi0 ][ imass*2 ]] + photons_rec[perm4g[imin_etapi0][imass*2+1]] ).M(), detnr[idet], GetTagger()->GetTaggedTime(tag) );
                    }
                }
                else{
                    IMgg_v_det_etapi0_TAPS->Fill( (photons_rec[perm4g[ imin_etapi0 ][ imass*2 ]] + photons_rec[perm4g[imin_etapi0][imass*2+1]] ).M(), detnr[idet], GetTagger()->GetTaggedTime(tag) );
                    if(imass == 0){
                        double En = (photons_rec[perm4g[ imin_etapi0 ][ j ]]).E();
                        double Th = (photons_rec[perm4g[ imin_etapi0 ][ j]]).Theta()*TMath::RadToDeg();
                        double mgg = (photons_rec[perm4g[ imin_etapi0 ][ imass*2 ]] + photons_rec[perm4g[imin_etapi0][imass*2+1]] ).M();
                        int nEN =  int(En/10.0);
                        int nTH =  150*int(Th/1.0);
                        int nBIN2 = nEN + nTH;
                        if( PbWO4[detnr[idet]] == -1 ){
                            four_fit_m_etagg_v_eth_PbWO4->Fill( nBIN2, mgg , GetTagger()->GetTaggedTime(tag));
                            if(etapi_g_cand->IsInside(photons_rec[perm4g[ imin_etapi0 ][ j ]].E(), photons_rec[perm4g[ imin_etapi0 ][ j ]].Theta()*TMath::RadToDeg())){
//                                four_fit_m_etagg_v_eth_PbWO4_2->Fill( nBIN2, mgg , GetTagger()->GetTaggedTime(tag));
                                int nEN =  int(En/40.0);
                                int nBIN3 = 30*detnr[idet] + nEN;
                                fit_mgg_eta_v_TAPS_2->Fill(nBIN3, MASS_ETA/mgg, GetTagger()->GetTaggedTime(tag));
                            }
                        }
                        else{
                            four_fit_m_etagg_v_eth_BaF2->Fill( nBIN2, mgg , GetTagger()->GetTaggedTime(tag));
                            if(etapi_g_cand->IsInside(photons_rec[perm4g[ imin_etapi0 ][ j ]].E(), photons_rec[perm4g[ imin_etapi0 ][ j ]].Theta()*TMath::RadToDeg())){
//                                four_fit_m_etagg_v_eth_BaF2_2->Fill( nBIN2, mgg , GetTagger()->GetTaggedTime(tag));
                                int nEN =  int(En/40.0);
                                int nBIN3 = 30*detnr[idet] + nEN;
                                fit_mgg_eta_v_TAPS_2->Fill(nBIN3, MASS_ETA/mgg, GetTagger()->GetTaggedTime(tag));
                            }
                        }

                    }
//                    else{}
                }
            }
        }
    }
}

void AdlarsonPhysics::sixgAnalysis(UInt_t ipr)
{
    for(Int_t tag = 0; tag < GetTagger()->GetNTagged(); tag++)
    {
        Double_t chi2_min = 1.0e6;
        Double_t prob_min = 1.0e6;

        IM6g_vec.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
        six_time_TaggedTime.Fill(GetTagger()->GetTaggedTime(tag));

        std::vector<Int_t> set;     // particle conf: tagger, proton, photon 1, 2, 3, 4, 5, 6;
        std::vector<Int_t> set_min; // best particle configuration
        set.resize(0);
        set_min.resize(0);
//        Is_CB_6g.resize(0);

        for(UInt_t i = 0; i < nPhotons_six; i++)
        {
            photons_rec[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
            photons_fit[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
            photons_fit_final[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
        }

        set.push_back(tag);

        beam.SetFromVector( GetTagger()->GetVector(tag) );
        obs.resize(0);
        unc.resize(0);
        unc = Get_unc(3,0,obs);
        beam.Smear(unc, 2);
        beam.APLCONSettings();

        set.push_back(ipr);
        proton.SetFromVector( GetTracks()->GetVector(ipr) );
        obs.resize(0);
        obs.push_back(GetTracks()->GetCentralCrystal(ipr));
        unc.resize(0);
        unc = Get_unc( 2, 2, obs);
        proton.Smear(unc, 1);
//        Is_CB_6g.push_back(false);

        UInt_t n_photons = 0;
        for ( UInt_t jgam = 0; jgam < ClustersInTime.size() ; jgam++ ){
            UInt_t kgam = ClustersInTime[jgam];
            if( kgam != ipr ){ // the id proton cluster
                set.push_back(kgam);
                photons_rec[n_photons] = GetTracks()->GetVector(kgam);
                IM6g_vec += photons_rec[n_photons];
                Photons_six[n_photons].SetFromVector( GetTracks()->GetVector(kgam) );

                if( GetTracks()->HasCB(kgam) ){
                    obs.resize(0);
                    obs.push_back(Photons_six[n_photons].Ek);
                    obs.push_back(Photons_six[n_photons].Theta);
                    unc.resize(0);
                    unc = Get_unc( 1, 1, obs);
                    Photons_six[n_photons].Smear(unc, 0);
//                    Is_CB_6g.push_back(true);
                }
                else{   // ( GetTracks()->HasTAPS(kgam) )
                    obs.resize(0);
                    obs.push_back(Photons_six[n_photons].Ek);
                    obs.push_back(Photons_six[n_photons].Theta);
                    unc.resize(0);
                    unc = Get_unc( 2, 1, obs);
                    Photons_six[n_photons].Smear(unc, 0);
//                    Is_CB_6g.push_back(false);
                }
                n_photons++;
            }
        }

        six_rec_IM->Fill(IM6g_vec.M());
        six_rec_IM_v_MMp->Fill(IM6g_vec.M(), MMp_vec.M());
        const APLCON::Result_t& result = kinfit.DoFit();
        if(result.Status == APLCON::Result_Status_t::Success){
            if( result.ChiSquare <  chi2_min ){
                chi2_min = result.ChiSquare;
                prob_min = result.Probability;
                set_min = set;
            }
        }
        // Here run kinfit with the best combination:

        if( set_min.size() == 0 ) continue;
        if( (prob_min < 0.01) || ( TMath::Abs(prob_min - 1.0e6) < 1.0e-4) ) continue;
        detnr.resize(0);

        // For the best combination with pdf >= 0.01 obtain the result.

        IM6g_vec.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
        TLorentzVector etap_fit(0.0, 0.0, 0.0, 0.0);

        beam.SetFromVector( GetTagger()->GetVector(set_min[0]) );
        obs.resize(0);
        unc.resize(0);
        unc = Get_unc(3,0,obs);
        beam.Smear(unc, 2);
        beam.APLCONSettings();

        proton.SetFromVector( GetTracks()->GetVector(set_min[1]) );
        obs.resize(0);
        obs.push_back(GetTracks()->GetCentralCrystal(set_min[1]));
        unc.resize(0);
        unc = Get_unc( 2, 2, obs);
        proton.Smear(unc, 1);

        Int_t n_photons_min = 0;
        for ( Int_t jgam_min = 0; jgam_min < 6 ; jgam_min++ )
        {
            photons_rec[n_photons_min] = GetTracks()->GetVector( set_min[jgam_min+2] );
            Photons_six[n_photons_min].SetFromVector( GetTracks()->GetVector(set_min[jgam_min+2]));
            detnr.push_back( GetTracks()->GetCentralCrystal( set_min[jgam_min+2]) );
            IM6g_vec += photons_rec[n_photons_min];

            if( GetTracks()->HasCB(set_min[jgam_min+2]) ){
                obs.resize(0);
                obs.push_back(Photons_six[n_photons_min].Ek);
                obs.push_back(Photons_six[n_photons_min].Theta);
                unc.resize(0);
                unc = Get_unc( 1, 1, obs);
                Photons_six[n_photons_min].Smear(unc, 0);
//                Is_CB_6g.push_back(true);
            }
            else{
                // ( GetTracks()->HasTAPS(jgam_min+2) )
                obs.resize(0);
                obs.push_back(Photons_six[n_photons_min].Ek);
                obs.push_back(Photons_six[n_photons_min].Theta);
                unc.resize(0);
                unc = Get_unc( 2, 1, obs);
                Photons_six[n_photons_min].Smear(unc, 0);
//                Is_CB_6g.push_back(false);
            }
            n_photons_min++;
        }

        six_rec_IM->Fill(IM6g_vec.M(),GetTagger()->GetTaggedTime(tag));
        six_rec_IM_v_MMp->Fill( MMp, IM6g_vec.M(),GetTagger()->GetTaggedTime(tag) );

        const APLCON::Result_t& result_min = kinfit.DoFit();
        if(result_min.Status == APLCON::Result_Status_t::Success)
        {
//            Double_t zrec = result_min.Variables.at("v_z").Value.After;
//            Double_t ztrue = etapr_6gTrue.GetTrueVertex().Z();
//            Double_t diff = zrec -ztrue;
//            true_six_fit_dz_v_z->Fill(ztrue, ztrue - zrec);

           six_fit_chi2->Fill( result_min.ChiSquare, GetTagger()->GetTaggedTime(tag) );
           six_fit_pdf->Fill( result_min.Probability, GetTagger()->GetTaggedTime(tag) );

           for(UInt_t igam_fit = 0; igam_fit < Photons_six.size(); igam_fit++){
               photons_fit[igam_fit] = FitParticle::Make(Photons_six[igam_fit], 0);
               etap_fit += photons_fit[igam_fit];
           }
           TLorentzVector proton_fit = FitParticle::Make(proton, MASS_PROTON);

           six_fit_IM->Fill( etap_fit.M(),GetTagger()->GetTaggedTime(tag) );
           proton_fit_e_v_th->Fill(proton_fit.E()-MASS_PROTON, proton_fit.Theta()*TMath::RadToDeg(), GetTagger()->GetTaggedTime(tag));

           for(UInt_t igam_fit = 0; igam_fit < Photons_six.size(); igam_fit++)
                six_fit_EvTh_g->Fill(photons_fit[igam_fit].E(), photons_fit[igam_fit].Theta()*TMath::RadToDeg(), GetTagger()->GetTaggedTime(tag) );

           Int_t inr = 0;
           for(const auto& it_map : result_min.Variables){
                const APLCON::Result_Variable_t& var = it_map.second;
                six_fit_Pulls->Fill(var.Pull, inr, GetTagger()->GetTaggedTime(tag));
                inr++;
           }


           Int_t nEn, nTh, nbin;
           for(UInt_t iPull = 0; iPull < nPhotons_six; iPull++){
                string s = Form("Photon%0i[0]", iPull+1); // Energy
                string t = Form("Photon%0i[1]", iPull+1); // theta
                string u = Form("Photon%0i[2]", iPull+1); // phi

                nEn =  int(photons_rec[iPull].E()/20.);
                nTh =  75*int(photons_rec[iPull].Theta()*TMath::RadToDeg()/2.0);
                nbin = nEn + nTh;

                six_fit_Pulls_g_E_vs_eth->Fill(nbin, result_min.Variables.at(s).Pull, GetTagger()->GetTaggedTime(tag));
                six_fit_Pulls_g_th_vs_eth->Fill(nbin, result_min.Variables.at(t).Pull, GetTagger()->GetTaggedTime(tag));
                six_fit_Pulls_g_phi_vs_eth->Fill(nbin, result_min.Variables.at(u).Pull, GetTagger()->GetTaggedTime(tag));

                if( GetTracks()->HasCB( set[iPull+2]) ){
                    six_fit_Pulls_g_E_vs_E_CB->Fill(photons_rec[iPull].E(), result_min.Variables.at(s).Pull, GetTagger()->GetTaggedTime(tag));
                    six_fit_Pulls_g_th_vs_th_CB->Fill(photons_rec[iPull].Theta()*TMath::RadToDeg(), result_min.Variables.at(t).Pull, GetTagger()->GetTaggedTime(tag));
                    six_fit_Pulls_g_phi_vs_th_CB->Fill(photons_rec[iPull].Theta()*TMath::RadToDeg(), result_min.Variables.at(u).Pull, GetTagger()->GetTaggedTime(tag));
                }
                else{
                    six_fit_Pulls_g_E_vs_E_TAPS->Fill(photons_rec[iPull].E(), result_min.Variables.at(s).Pull, GetTagger()->GetTaggedTime(tag));
                    six_fit_Pulls_g_th_vs_th_TAPS->Fill(photons_rec[iPull].Theta()*TMath::RadToDeg(), result_min.Variables.at(t).Pull, GetTagger()->GetTaggedTime(tag));
                    six_fit_Pulls_g_phi_vs_th_TAPS->Fill(photons_rec[iPull].Theta()*TMath::RadToDeg(), result_min.Variables.at(u).Pull, GetTagger()->GetTaggedTime(tag));
                }
           }

           six_fit_Pulls_p_th_vs_det_TAPS->Fill(GetTracks()->GetCentralCrystal(set_min[1]), result_min.Variables.at("Proton[1]").Pull, GetTagger()->GetTaggedTime(tag));
           six_fit_Pulls_p_phi_vs_det_TAPS->Fill(GetTracks()->GetCentralCrystal(set_min[1]), result_min.Variables.at("Proton[2]").Pull, GetTagger()->GetTaggedTime(tag));

           sigma_eta = 28; sigma_pi0 = 17;

           Double_t    chi2min_eta2pi  = 1.0e6;
           Double_t    chi2min_3pi     = 1.0e6;
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

           GetBest6gCombination(sigma_eta, sigma_pi0, chi2min_eta2pi, chi2min_3pi, imin_eta2pi, imin_3pi , etatwopi_comb, threepi_comb  );

           test_correct_hypothesis(probmin_eta2pi, probmin_3pi, set_min, imin_eta2pi, imin_3pi, etatwopi_comb, threepi_comb);
           six_fit_PDF_eta2pi_v_3pi->Fill(probmin_eta2pi, probmin_3pi, GetTagger()->GetTaggedTime(tag));


           // photons after kinematical fit only under energy and momentum constraint
           TLorentzVector g[3];
           g[0] = photons_fit[imin_eta2pi[0]] + photons_fit[imin_eta2pi[1]];
           g[1] = photons_fit[imin_eta2pi[2]] + photons_fit[imin_eta2pi[3]];
           g[2] = photons_fit[imin_eta2pi[4]] + photons_fit[imin_eta2pi[5]];

           TLorentzVector h[3];
           h[0] = photons_fit[imin_3pi[0]] + photons_fit[imin_3pi[1]];
           h[1] = photons_fit[imin_3pi[2]] + photons_fit[imin_3pi[3]];
           h[2] = photons_fit[imin_3pi[4]] + photons_fit[imin_3pi[5]];

           // reconstructed photons
           TLorentzVector rc[3];
           rc[0] = photons_rec[imin_3pi[0]] + photons_rec[imin_3pi[1]];
           rc[1] = photons_rec[imin_3pi[2]] + photons_rec[imin_3pi[3]];
           rc[2] = photons_rec[imin_3pi[4]] + photons_rec[imin_3pi[5]];

           if( (probmin_3pi > 0.01) && (probmin_eta2pi < 0.2) ){ // dir 3pi0
               if(etap_fit.M() > 650.0)
                    CB_EnergySum_3pi0->Fill(GetTrigger()->GetEnergySum(),GetTagger()->GetTaggedTime(tag));

               for(Int_t t = 0; t < GetTracks()->GetNTracks(); t++)
                   if(GetTracks()->HasCB(t))
                        time_clusters_CB_3pi0->Fill(  GetTracks()->GetTime(t), GetTracks()->GetCentralCrystal(t) ,GetTagger()->GetTaggedTime(tag));

               six_fit_best_3pi_IM_v_E->Fill(h[0].E(), h[0].M(), GetTagger()->GetTaggedTime(tag));
               six_fit_best_3pi_IM_v_E->Fill(h[1].E(), h[1].M(), GetTagger()->GetTaggedTime(tag));
               six_fit_best_3pi_IM_v_E->Fill(h[2].E(), h[2].M(), GetTagger()->GetTaggedTime(tag));

               six_phy_3pi_IMpipi_v_IMppi->Fill((h[0]+h[1]).M2()/1.0e6, (h[2]+proton_fit).M2()/1.0e6, GetTagger()->GetTaggedTime(tag));
               six_phy_3pi_IMpipi_v_IMppi->Fill((h[2]+h[0]).M2()/1.0e6, (h[1]+proton_fit).M2()/1.0e6, GetTagger()->GetTaggedTime(tag));
               six_phy_3pi_IMpipi_v_IMppi->Fill((h[1]+h[2]).M2()/1.0e6, (h[0]+proton_fit).M2()/1.0e6, GetTagger()->GetTaggedTime(tag));

               for(Int_t isix = 0; isix < 6; isix++){
                    six_fit_best_3pi0_pi_E_v_th->Fill(photons_fit[imin_3pi[isix]].E(), photons_fit[imin_3pi[isix]].Theta()*TMath::RadToDeg(), GetTagger()->GetTaggedTime(tag));

                    int imass = int(isix/2);
//                    double mgg = rc[imass].M();
                    int nEN =  int(photons_rec[imin_3pi[isix]].E()/10.);
                    int nTH =  150*int(photons_rec[imin_3pi[isix]].Theta()*TMath::RadToDeg()/1.0);
                    int nBIN = nEN + nTH;

                    if(GetTracks()->HasCB(set_min[imin_3pi[isix]+2]))
                    {
                        if(etap_fit.M() > 650.0){

                            IMgg_v_det_3pi0_CB->Fill( rc[imass].M(), detnr[imin_3pi[isix]], GetTagger()->GetTaggedTime(tag));
                            IMgg_v_det_3pi0_CB_fit->Fill( h[imass].M(), detnr[imin_3pi[isix]], GetTagger()->GetTaggedTime(tag));
                            six_fit_mgg_v_eth->Fill(nBIN, rc[imass].M(),GetTagger()->GetTaggedTime(tag));

                            if(threepi_g_cand->IsInside(photons_rec[imin_3pi[isix]].E(), photons_rec[imin_3pi[isix]].Theta()*TMath::RadToDeg())){
//                                six_fit_mgg_v_eth_2->Fill(nBIN, rc[imass].M(),GetTagger()->GetTaggedTime(tag));
                                int nEN =  int(photons_rec[imin_3pi[isix]].E()/40.0);
                                int nBIN3 = 25*detnr[imin_3pi[isix]] + nEN;
                                fit_mgg_pi_v_CB_2->Fill(nBIN3, MASS_PI0/rc[imass].M(), GetTagger()->GetTaggedTime(tag));
                            }
                        }
                    }
                    else{
                        if(etap_fit.M() > 650.0){
                            IMgg_v_det_3pi0_TAPS->Fill( rc[imass].M() , detnr[imin_3pi[isix]], GetTagger()->GetTaggedTime(tag));
                            IMgg_v_det_3pi0_TAPS_fit->Fill( h[imass].M() , detnr[imin_3pi[isix]], GetTagger()->GetTaggedTime(tag));


                            string s = Form("Photon%0i[0]", isix+1); // Energy
                            string t = Form("Photon%0i[1]", isix+1); // theta
                            string u = Form("Photon%0i[2]", isix+1); // phi
                            int nEN =  int(photons_rec[imin_3pi[isix]].E()/20.);
                            int nTH =  75*int(photons_rec[imin_3pi[isix]].Theta()*TMath::RadToDeg()/2.0);
                            int nBIN2 = nEN + nTH;

                            if( PbWO4[detnr[imin_3pi[isix]]] == -1 ){
                                six_fit_mgg_v_e_PbWO4->Fill(photons_rec[imin_3pi[isix]].E(), rc[imass].M());
                                six_fit_mgg_v_eth_PbWO4->Fill(nBIN2, rc[imass].M(), GetTagger()->GetTaggedTime(tag));

                                if(threepi_g_cand->IsInside(photons_rec[imin_3pi[isix]].E(), photons_rec[imin_3pi[isix]].Theta()*TMath::RadToDeg())){
//                                    six_fit_mgg_v_eth_PbWO4_2->Fill(nBIN2, rc[imass].M(),GetTagger()->GetTaggedTime(tag));
                                    int nEN =  int(photons_rec[imin_3pi[isix]].E()/40.0);
                                    int nBIN3 = 30*detnr[imin_3pi[isix]] + nEN;
                                    fit_mgg_pi_v_TAPS_2->Fill(nBIN3, MASS_PI0/rc[imass].M(), GetTagger()->GetTaggedTime(tag));
                                }
                            }
                            else{
                                six_fit_mgg_v_e_BaF2->Fill(photons_rec[imin_3pi[isix]].E(), rc[imass].M());
                                six_fit_mgg_v_eth_BaF2->Fill(nBIN2, rc[imass].M() , GetTagger()->GetTaggedTime(tag));

                                if(threepi_g_cand->IsInside(photons_rec[imin_3pi[isix]].E(), photons_rec[imin_3pi[isix]].Theta()*TMath::RadToDeg())){
//                                    six_fit_mgg_v_eth_BaF2_2->Fill(nBIN2, rc[imass].M(),GetTagger()->GetTaggedTime(tag));
                                    int nEN =  int(photons_rec[imin_3pi[isix]].E()/40.0);
                                    int nBIN3 = 30*detnr[imin_3pi[isix]] + nEN;
                                    fit_mgg_pi_v_TAPS_2->Fill(nBIN3, MASS_PI0/rc[imass].M(), GetTagger()->GetTaggedTime(tag));
                                }
                            }
                        }
                    }
               }
                    six_fit_IM_3pi->Fill(etap_fit.M(), GetTagger()->GetTaggedTime(tag));

                }
             if(PDF_cut->IsInside(probmin_eta2pi, probmin_3pi)){
//           if( (probmin_3pi < 0.02) && (probmin_eta2pi > 0.05) ){ //eta prime

               proton_fit_e_v_th_final->Fill(proton_fit.E()-MASS_PROTON, proton_fit.Theta()*TMath::RadToDeg(), GetTagger()->GetTaggedTime(tag));

               TLorentzVector etap_fit_final(0.0, 0.0, 0.0, 0.0);

               for(UInt_t igam_fit = 0; igam_fit < Photons_six.size(); igam_fit++)
                   etap_fit_final += photons_fit_final[igam_fit];

               TLorentzVector fin[3];
               fin[0] = photons_fit_final[0] + photons_fit_final[1];
               fin[1] = photons_fit_final[2] + photons_fit_final[3];
               fin[2] = photons_fit_final[4] + photons_fit_final[5];
               // calculate the theta and azimuthal angle. Perhaps this will reduce background from 2pi0 and etapi0

               Double_t dfi = (proton_fit.Phi() - (fin[0]+fin[1]+fin[2]).Phi())*TMath::RadToDeg();
               six_fit_final_coplanarity->Fill(TMath::Abs(dfi) - 180.0 , GetTagger()->GetTaggedTime(tag));

               TLorentzVector miss_proton(0., 0., 0., 0.);
               TLorentzVector target(0,0,0, MASS_PROTON);
               miss_proton = GetTagger()->GetVector(set_min[0]) + target - proton_fit;

               Double_t th_missproton   = miss_proton.Theta()*TMath::RadToDeg();
               Double_t th_etapr        = etap_fit_final.Theta()*TMath::RadToDeg();
               Double_t dth             = th_missproton - th_etapr;
               six_fit_final_theta_diff->Fill(dth , GetTagger()->GetTaggedTime(tag));

               bool eta_photon_split    = false;
               bool gamma_photon_split  = false;

               for(UInt_t t = 0; t < 2; t++)
                   if(photons_fit_final[t].Theta()*TMath::RadToDeg() > 18 && photons_fit_final[t].Theta()*TMath::RadToDeg() < 28)
                       eta_photon_split = true;

               for(UInt_t t = 0; t < 6; t++){
                   if(photons_fit_final[t].Theta()*TMath::RadToDeg() < 24 && photons_fit_final[t].E() < 100.)
                       gamma_photon_split = true;
               }

               if(gamma_photon_split && eta_photon_split){
                   for(UInt_t t = 0; t < 6; t++)
                        six_fit_EvTh_g_final_removed->Fill(photons_fit_final[t].E(), photons_fit_final[t].Theta()*TMath::RadToDeg() );
                   continue;
               }

               for(UInt_t t = 0; t < 6; t++)
                    six_fit_EvTh_g_final->Fill(photons_fit_final[t].E(), photons_fit_final[t].Theta()*TMath::RadToDeg() );


               six_fit_PDF_eta2pi_v_Meta2pi->Fill(probmin_eta2pi, etap_fit.M(),GetTagger()->GetTaggedTime(tag));
               six_fit_best_eta_IM_v_E->Fill(g[0].E(), g[0].M(), GetTagger()->GetTaggedTime(tag));
               six_fit_best_2pi_IM_v_E->Fill(g[1].E(), g[1].M(), GetTagger()->GetTaggedTime(tag));
               six_fit_best_2pi_IM_v_E->Fill(g[2].E(), g[2].M(), GetTagger()->GetTaggedTime(tag));

               // here fill Energy vs IM(gg) vs detector nr for eta and pi0 separately
               six_phy_etapr_v_BeamE->Fill((GetTagger()->GetVector(set_min[0])).E(), etap_fit.M(), GetTagger()->GetTaggedTime(tag));

               six_fit_best_eta->Fill( (g[0]).M(), GetTagger()->GetTaggedTime(tag) );
               six_fit_best_2pi->Fill( g[1].M(), GetTagger()->GetTaggedTime(tag) );
               six_fit_best_2pi->Fill( g[2].M(), GetTagger()->GetTaggedTime(tag) );

               six_fit_IM_eta2pi->Fill( etap_fit.M(), GetTagger()->GetTaggedTime(tag) );

               for(int isix = 0; isix < 6; isix++){
                   if(isix < 2)
                       six_fit_best_etapr_eta_E_v_th->Fill( photons_fit[imin_eta2pi[isix]].E(), photons_fit[imin_eta2pi[isix]].Theta()*TMath::RadToDeg(), GetTagger()->GetTaggedTime(tag));
                   else
                       six_fit_best_etapr_pi_E_v_th->Fill( photons_fit[imin_eta2pi[isix]].E(), photons_fit[imin_eta2pi[isix]].Theta()*TMath::RadToDeg(), GetTagger()->GetTaggedTime(tag));
               }

                   six_fit_etaprfinal_pdf->Fill( probmin_eta2pi, GetTagger()->GetTaggedTime(tag) );
                   six_fit_eta_PDF_v_Metapr->Fill(probmin_eta2pi, etap_fit_final.M(),GetTagger()->GetTaggedTime(tag));



                   six_phy_etapr_eta_v_BeamE->Fill((GetTagger()->GetVector(set_min[0])).E(), etap_fit_final.M(), GetTagger()->GetTaggedTime(tag));

                   Double_t m_etapi01_fit = 0;
                   Double_t m_etapi02_fit = 0;
                   Double_t m_2pi0_fit = 0;
                   m2pi0_metapi0( fin ,  m_etapi01_fit, m_etapi02_fit, m_2pi0_fit);
                   six_phy_M_pi1pi2_v_etapr->Fill(m_2pi0_fit / 1.0e3, etap_fit_final.M(), GetTagger()->GetTaggedTime(tag));
                   six_phy_M_etapi_v_etapr->Fill(m_etapi01_fit / 1.0e3, etap_fit_final.M(), GetTagger()->GetTaggedTime(tag));
                   six_phy_M_etapi_v_etapr->Fill(m_etapi02_fit / 1.0e3, etap_fit_final.M(), GetTagger()->GetTaggedTime(tag));
                   true_six_phy_dMpipi_v_Mpipi->Fill(m_2pi0_fit / 1.0e3, (m_2pi0_fit/1.0e3-m_2pi0True*1.0e3));

                   Double_t Xfit1 = -2.0;
                   Double_t Xfit2 = -2.0;
                   Double_t Yfit = -2.0;
                   Int_t DP_binnr_fit1 = -100;
                   Int_t DP_binnr_fit2 = -100;
                   DalitzPlot(fin, Xfit1, Xfit2, Yfit, DP_binnr_fit1, DP_binnr_fit2);

                   six_phy_DP->Fill(DP_binnr_fit1, etap_fit_final.M(), GetTagger()->GetTaggedTime(tag) );
                   six_phy_DP->Fill(DP_binnr_fit2, etap_fit_final.M(), GetTagger()->GetTaggedTime(tag) );
                   if( TMath::Abs(Xfit1 - Xtrue1) > TMath::Abs(Xfit1 - Xtrue2)){
                       Double_t Xtemp = Xtrue1;
                       Xtrue1 = Xtrue2;
                       Xtrue2 = Xtemp;
                   }

                   true_six_phy_dX_v_DPbin->Fill(DP_binnr_fit1, Xfit1 - Xtrue1);
                   true_six_phy_dY_v_DPbin->Fill(DP_binnr_fit1, Yfit - Ytrue);
                   true_six_phy_dX_v_DPbin->Fill(DP_binnr_fit2, Xfit2 - Xtrue2);
                   true_six_phy_dY_v_DPbin->Fill(DP_binnr_fit2, Yfit - Ytrue);
           }
        }
    }
}


void AdlarsonPhysics::GetBest6gCombination(Double_t& sigma_eta, Double_t& sigma_pi0, Double_t& chi2min_eta2pi, Double_t& chi2min_3pi, std::vector<int>& imin_eta2pi, std::vector<int>& imin_3pi, std::vector<comb>& etatwopi_comb, std::vector<comb>& threepi_comb )
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

    for( Int_t i = 0; i < 15 ; i++)
    {
       for (Int_t j = 0; j < 3; j++){
            imgg[j] = ( photons_fit[perm6g[i][0 + j*2]] + photons_fit[perm6g[i][1 + j*2]]).M();

            sgm_eta = 27.0;
            sgm_pi0 = 15.0;
            sgm_3pi0 = 15.0;

            chi2_eta2pi[j] = 0;
       }

       chi2_eta2pi[0] = std::pow((imgg[0] - MASS_ETA )/(sgm_eta), 2.0) + std::pow((imgg[1] - MASS_PI0 )/(sgm_pi0), 2) + std::pow((imgg[2] - MASS_PI0 )/(sgm_pi0), 2.0);
       chi2_eta2pi[1] = std::pow((imgg[1] - MASS_ETA )/(sgm_eta), 2.0) + std::pow((imgg[2] - MASS_PI0 )/(sgm_pi0), 2) + std::pow((imgg[0] - MASS_PI0 )/(sgm_pi0), 2.0);
       chi2_eta2pi[2] = std::pow((imgg[2] - MASS_ETA )/(sgm_eta), 2.0) + std::pow((imgg[0] - MASS_PI0 )/(sgm_pi0), 2) + std::pow((imgg[1] - MASS_PI0 )/(sgm_pi0), 2.0);

       for( Int_t t = 0; t < 3; t++ ){
           sig_perm.resize(0);

           if(t == 0){
               for(Int_t u = 0; u < 6; u++)
                    sig_perm.push_back(perm6g[i][u]);

               etatwopi_comb.push_back(comb(sig_perm, chi2_eta2pi[0]));
           }
           else if( t == 1){
               sig_perm.push_back(perm6g[i][2]);
               sig_perm.push_back(perm6g[i][3]);
               sig_perm.push_back(perm6g[i][0]);
               sig_perm.push_back(perm6g[i][1]);
               sig_perm.push_back(perm6g[i][4]);
               sig_perm.push_back(perm6g[i][5]);

               etatwopi_comb.push_back(comb(sig_perm, chi2_eta2pi[1]));
           }
           else{
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
        chi2_3pi = std::pow((imgg[0] - MASS_PI0 )/(sgm_3pi0), 2) + std::pow((imgg[1] - MASS_PI0 )/(sgm_3pi0), 2) + std::pow((imgg[2] - MASS_PI0 )/(sgm_3pi0), 2);
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
}


void AdlarsonPhysics::test_correct_hypothesis(Double_t& prob_eta2pi, Double_t& prob_3pi, std::vector<Int_t>& set_min, std::vector<int>& imin_eta2pi, std::vector<int>& imin_3pi, std::vector<comb>& imap_eta2pi, std::vector<comb>& imap_3pi)
{
    Double_t prob_3pi_fit       = -1.0e6;
    Double_t prob_eta2pi_fit    = -1.0e6;
    Double_t chi2_3pi_fit       = 1.0e6;
    Double_t chi2_eta2pi_fit    = 1.0e6;

    Int_t iteration_place_etatwopi = -100;
    Int_t iteration_place_threepi = -100;
    imin_eta2pi.resize(0);
    imin_3pi.resize(0);

    std::vector<int> imin_eta2pi_temp;
    std::vector<int> imin_3pi_temp;
    photons_fit_final.resize(0);

    for( Int_t icand = 0; icand < 2; icand++){

        beam_3pi.SetFromVector( GetTagger()->GetVector(set_min[0]) );
        obs.resize(0);
        unc.resize(0);
        unc = Get_unc(3,0,obs);
        beam_3pi.Smear(unc, 2);
        beam_3pi.APLCONSettings();

        proton_3pi.SetFromVector( GetTracks()->GetVector(set_min[1]) );
        obs.resize(0);
        obs.push_back(GetTracks()->GetCentralCrystal(set_min[1]));
        unc.resize(0);
        unc = Get_unc( 2, 2, obs);
        proton_3pi.Smear(unc, 1);

        beam_eta2pi.SetFromVector( GetTagger()->GetVector(set_min[0]) );
        obs.resize(0);
        unc.resize(0);
        unc = Get_unc(3,0,obs);
        beam_eta2pi.Smear(unc, 2);
        beam_eta2pi.APLCONSettings();

        proton_eta2pi.SetFromVector( GetTracks()->GetVector(set_min[1]) );
        obs.resize(0);
        obs.push_back(GetTracks()->GetCentralCrystal(set_min[1]));
        unc.resize(0);
        unc = Get_unc( 2, 2, obs);
        proton_eta2pi.Smear(unc, 1);

        Int_t n_photons_min = 0;
        imin_eta2pi_temp.resize(0);
        imin_3pi_temp.resize(0);
        Photons_six_3pi.resize(0);
        Photons_six_eta2pi.resize(0);
        imin_eta2pi_temp = imap_eta2pi.at(icand).first;
        imin_3pi_temp = imap_3pi.at(icand).first;
        for ( Int_t jgam_min = 0; jgam_min < 6 ; jgam_min++ ){
            Photons_six_3pi[n_photons_min].SetFromVector( GetTracks()->GetVector( set_min[ imin_3pi_temp[jgam_min] + 2] ) );
            Photons_six_eta2pi[n_photons_min].SetFromVector( GetTracks()->GetVector( set_min[ imin_eta2pi_temp[jgam_min] + 2] ) );
            if( GetTracks()->HasCB(set_min[ imin_3pi_temp[jgam_min] + 2] ) ){
                obs.resize(0);
                obs.push_back(Photons_six_3pi[n_photons_min].Ek);
                obs.push_back(Photons_six_3pi[n_photons_min].Theta);
                unc.resize(0);
                unc = Get_unc( 1, 1, obs);
                Photons_six_3pi[n_photons_min].Smear(unc ,0);
            }
            else{ // ( GetTracks()->HasTAPS(jgam_min) )
                obs.resize(0);
                obs.push_back(Photons_six_3pi[n_photons_min].Ek);
                obs.push_back(Photons_six_3pi[n_photons_min].Theta);
                unc.resize(0);
                unc = Get_unc( 2, 1, obs);
                Photons_six_3pi[n_photons_min].Smear(unc, 0);
            }

            if( GetTracks()->HasCB( set_min[ imin_eta2pi_temp[jgam_min] + 2 ] ) )
            {
                obs.resize(0);
                obs.push_back(Photons_six_eta2pi[n_photons_min].Ek);
                obs.push_back(Photons_six_eta2pi[n_photons_min].Theta);
                unc.resize(0);
                unc = Get_unc( 1, 1, obs);
                Photons_six_eta2pi[n_photons_min].Smear(unc ,0);
            }
            else // ( GetTracks()->HasTAPS(jgam_min) )
            {
                obs.resize(0);
                obs.push_back(Photons_six_eta2pi[n_photons_min].Ek);
                obs.push_back(Photons_six_eta2pi[n_photons_min].Theta);
                unc.resize(0);
                unc = Get_unc( 2, 1, obs);
                Photons_six_eta2pi[n_photons_min].Smear(unc ,0);
            }
            n_photons_min++;
        }

        const APLCON::Result_t& result_3pi = kinfit3pi.DoFit();
        if(result_3pi.Status == APLCON::Result_Status_t::Success){

            if( (result_3pi.ChiSquare < chi2_3pi_fit) && (result_3pi.Probability < 1.0) ){

                chi2_3pi_fit = result_3pi.ChiSquare;
                prob_3pi_fit = result_3pi.Probability;
                imin_3pi.resize(0);
                imin_3pi = imin_3pi_temp;

                iteration_place_threepi = icand;

            }
        }

        const APLCON::Result_t& result_eta2pi = kinfiteta2pi.DoFit();
        if(result_eta2pi.Status == APLCON::Result_Status_t::Success){
            if((result_eta2pi.ChiSquare < chi2_eta2pi_fit) && (result_eta2pi.Probability < 1.0)){
                chi2_eta2pi_fit = result_eta2pi.ChiSquare;
                prob_eta2pi_fit = result_eta2pi.Probability;
                imin_eta2pi = imin_eta2pi_temp;
                iteration_place_etatwopi = icand;
                photons_fit_final.resize(0);
                for(UInt_t igam_fit = 0; igam_fit < Photons_six.size(); igam_fit++)
                    photons_fit_final[igam_fit] = FitParticle::Make(Photons_six_eta2pi[igam_fit], 0.0);
            }
        }
    }

    if(prob_eta2pi_fit > 0.0 && prob_eta2pi_fit < 1.0)
        prob_eta2pi = prob_eta2pi_fit;
    else
        prob_eta2pi = 0.0;

    if(prob_3pi_fit > 0.0 && prob_3pi_fit < 1.0)
        prob_3pi = prob_3pi_fit;
    else
        prob_3pi = 0.0;


    six_fit_which_place_best_etapr_cand->Fill(iteration_place_etatwopi);
    six_fit_which_place_best_3pi_cand->Fill(iteration_place_threepi);

    return;
}


void AdlarsonPhysics::tengAnalysis(UInt_t ipr)
{
    for(Int_t tag = 0; tag < GetTagger()->GetNTagged(); tag++)
    {
        sigma_eta = 20; sigma_pi0 = 10;
        Is_CB_10g.resize(0);
        photons_rec.resize(0);
        IM10g_vec.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
        TLorentzVector etap10g_fit(0.0, 0.0, 0.0, 0.0);

        std::vector<Int_t> set, set_min;     // particle conf: tagger, proton, photon 1, 2, 3, 4, 5, 6, 7, 8, 9, 10;
        set.resize(0);
        set_min.resize(0);

        for(UInt_t i = 0; i < nPhotons_ten; i++)
        {
            photons_rec[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
            photons_fit[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
        }

        Double_t chi2_min = 1.0e6;
        Double_t prob_min = 1.0e6;

        set.push_back(tag);

        beam10g.SetFromVector( GetTagger()->GetVector(tag) );
        std::vector<double> obs;
        std::vector<double> unc;
        obs.resize(0);
        unc.resize(0);
        unc = Get_unc(3,0,obs);
        beam10g.Smear(unc, 2);
        beam10g.APLCONSettings();

        set.push_back(ipr);
        proton10g.SetFromVector( GetTracks()->GetVector(ipr) );
        obs.resize(0);
//        obs.push_back(proton10g.Ek);
//        obs.push_back(proton10g.Theta);
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
                photons_rec[n_photons] = GetTracks()->GetVector(kgam);
                Photons_ten[n_photons].SetFromVector( GetTracks()->GetVector(kgam) );
                IM10g_vec += GetTracks()->GetVector(kgam);

                if( GetTracks()->HasCB(kgam) )
                {
                    obs.resize(0);
                    obs.push_back(Photons_ten[n_photons].Ek);
                    obs.push_back(Photons_ten[n_photons].Theta);
                    unc.resize(0);
                    unc = Get_unc( 1, 1, obs);
                    Photons_ten[n_photons].Smear(unc, 0);
                    Is_CB_10g.push_back(true);
                }
                else // ( GetTracks()->HasTAPS(kgam) )
                {
                    obs.resize(0);
                    obs.push_back(Photons_ten[n_photons].Ek);
                    obs.push_back(Photons_ten[n_photons].Theta);
                    unc.resize(0);
                    unc = Get_unc( 2, 1, obs);
                    Photons_ten[n_photons].Smear(unc, 0);
                    Is_CB_10g.push_back(false);
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

        for(UInt_t igam_fit = 0; igam_fit < Photons_ten.size(); igam_fit++)
        {
            photons_fit[igam_fit] = FitParticle::Make(Photons_ten[igam_fit], 0);
            ten_fit_EvTh_g->Fill(photons_fit[igam_fit].E(), photons_fit[igam_fit].Theta()*TMath::RadToDeg(), GetTagger()->GetTaggedTime(tag) );
            
            etap10g_fit += photons_fit[igam_fit];
        }

        Double_t chi2_eta2pi0 = 1.0e6;
        std::vector<int> imin_eta3pi2pi;
        imin_eta3pi2pi.resize(0);
        IM10g_fit->Fill(etap10g_fit.M(),GetTagger()->GetTaggedTime(tag));

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
            if( Is_CB_10g[imin_eta3pi2pi[jgam]] )
            {
                obs.resize(0);
                obs.push_back(Photons_ten_eta2pi[jgam].Ek);
                obs.push_back(Photons_ten_eta2pi[jgam].Theta);
                unc.resize(0);
                unc = Get_unc( 1, 1, obs);
                Photons_ten_eta2pi[jgam].Smear(unc, 0);
            }
            else // Is_CB_10g[jgam] is false
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

                Double_t Xfit1 = -2.0;
                Double_t Xfit2 = -2.0;
                Double_t Yfit = -2.0;
                Int_t DP_binnr_fit1 = -100;
                Int_t DP_binnr_fit2 = -100;
                DalitzPlot(fin, Xfit1, Xfit2, Yfit, DP_binnr_fit1, DP_binnr_fit2);

                if( TMath::Abs(Xfit1 - Xtrue1) > TMath::Abs(Xfit1 - Xtrue2))
                {
                    Double_t Xtemp = Xtrue1;
                    Xtrue1 = Xtrue2;
                    Xtrue2 = Xtemp;
                }

                true_ten_phy_dX_v_DPbin->Fill(DP_binnr_fit1, Xfit1 - Xtrue1);
                true_ten_phy_dY_v_DPbin->Fill(DP_binnr_fit1, Yfit - Ytrue);
                true_ten_phy_dX_v_DPbin->Fill(DP_binnr_fit2, Xfit2 - Xtrue2);
                true_ten_phy_dY_v_DPbin->Fill(DP_binnr_fit2, Yfit - Ytrue);

                ten_fit_X_v_pdf_eta2pi->Fill(result10g_eta2pi.Probability, Xfit1 - Xtrue1);
                ten_fit_X_v_pdf_eta2pi->Fill(result10g_eta2pi.Probability, Xfit2 - Xtrue2);
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

    std::vector<Int_t> teng;
    std::vector<Int_t> teng_min;
    teng_comb.resize(0);

    for(Int_t iperm = 0; iperm < 210; iperm++)
    {
        eta_6g_cand.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
        std::vector<int> pi0c;
        im6g = 0;
        Double_t im2g_1, im2g_2;

        for(int jperm = 0; jperm < 6; jperm++)
            eta_6g_cand += ( photons_fit[perm6outof10g[iperm][jperm]] );

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
            for(Int_t igam = 0; igam < 6; igam++)
                teng.push_back(perm6outof10g[iperm][igam]);

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
            im2g_1 = ( photons_fit[teng[6]] + photons_fit[teng[7]] ).M();
            im2g_2 = ( photons_fit[teng[8]] + photons_fit[teng[9]] ).M();

            chi2 = TMath::Power( (im6g - MASS_ETA )/sigma_eta ,2) + TMath::Power( (im2g_1 - MASS_PI0 )/sigma_pi0 ,2) + TMath::Power( (im2g_2 - MASS_PI0 )/sigma_pi0 ,2);
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


void AdlarsonPhysics::DalitzPlot(const TLorentzVector g[3] , Double_t &X1, Double_t &X2, Double_t &Y , Int_t &DP_nr1, Int_t &DP_nr2)
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
    T_pi1 = pi01_ep_rest.E() - pi01_ep_rest.M();
    T_pi2 = pi02_ep_rest.E() - pi02_ep_rest.M();

    Q = T_eta + T_pi1 + T_pi2;
    X1 = ( TMath::Sqrt(3.0) / Q ) * ( T_pi1 - T_pi2 );
    X2 = ( TMath::Sqrt(3.0) / Q ) * ( T_pi2 - T_pi1 );
    Y = ( MASS_ETA + 2 * MASS_PI0 ) / MASS_PI0 * ( T_eta / Q ) - 1.0;

    DP_nr1 = int ( (X1  + X_max ) / bin_width ) + 40*int( (Y  + Y_max ) / bin_width );
    DP_nr2 = int ( (X2  + X_max ) / bin_width ) + 40*int( (Y  + Y_max ) / bin_width );

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
    double e_g, th_g, p_th, p_fi;

    if( particle == 1 )         // gamma
    {
        double Ek = obs[0];
        double theta = obs[1]*TMath::RadToDeg();

//        TH2*            g_e_c1;
//        TH2*            g_th_c1;
//        TH1*            p_th_c1;
//        TH1*            p_fi_c1;

        Ek_s    = g_e->GetBinContent( g_e->FindBin(Ek,theta) );
        Theta_s = g_th->GetBinContent( g_th->FindBin(Ek,theta) );
        Phi_s   = g_phi->GetBinContent( g_phi->FindBin(Ek,theta) );

//        e_g     = g_e_c1->GetBinContent(g_e_c1->FindBin(Ek,theta) );
//        th_g    = g_th_c1->GetBinContent(g_th_c1->FindBin(Ek,theta) );

//        if( e_g < 0.5 || e_g > 1.5 )
//            e_g = 1.0;
//        if( th_g < 0.5 || th_g > 1.5 )
//            th_g = 1.0;

//        Ek_s    = Ek_s*e_g*e_g;
//        Theta_s = Theta_s*th_g*th_g;

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
        Theta_s = p_TAPS_th->GetBinContent(p_TAPS_th->FindBin( inr ));
        Phi_s = p_TAPS_fi->GetBinContent(p_TAPS_fi->FindBin( inr ));


        p_th    = p_th_c1->GetBinContent(p_th_c1->FindBin(inr) );
        p_fi    = p_fi_c1->GetBinContent(p_fi_c1->FindBin(inr) );

//        if( p_th < 0.5 || p_th > 1.7 )
//            p_th = 1.0;
//        if( p_fi < 0.5 || p_fi > 1.7 )
//            p_fi = 1.0;

//        Theta_s = Theta_s*p_th*p_th;
//        Phi_s   = Phi_s*p_fi*p_fi;


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


void AdlarsonPhysics::Energy_corr()
{
    Double_t Erec, Ec_temp, DeltaE, Ec;
    for (int i = 0; i < GetTracks()->GetNTracks() ; i++)
    {
        if( GetTracks()->HasCB(i) )
        {
            Erec = GetTracks()->GetVector(i).E();
            Ec_temp = CBgain[GetTracks()->GetCentralCrystal(i)]*Erec;

            DeltaE = Ec_temp*(Double_t)EvdetCB->GetBinContent(EvdetCB->FindBin(Ec_temp, GetTracks()->GetTheta(i)));
            Ec = Ec_temp - DeltaE;

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
//    Double_t EPT_t;

    Trigger_t   = pRandoms->Gaus(0, 2.5);
//    EPT_t       = Trigger_t + pRandoms->Gaus(0, 0.2);
    CB_t        = Trigger_t + pRandoms->Gaus(0, 0.8);
    TAPS_t      = Trigger_t + pRandoms->Gaus(0, 0.1);

    for( int j = 0; j < GetTracks()->GetNTracks() ; j++ )
    {
        if(GetTracks()->HasCB(j))
            tracks->SetTime(j, GetTracks()->GetTime(j) + CB_t);
        if(GetTracks()->HasTAPS(j))
            tracks->SetTime(j,GetTracks()->GetTime(j) + TAPS_t);
    }

}
