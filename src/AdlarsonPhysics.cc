#include "GParticleReconstruction.h"
#include "AdlarsonPhysics.h"
#include "GTrue.h"
#include <cmath>
#include "TH1.h"
#include "TProfile.h"

#include <APLCON.hpp>

std::default_random_engine AdlarsonPhysics::FitParticle::generator;

AdlarsonPhysics::AdlarsonPhysics():
    time_TOF("time_TOF", "TAPS time_TOF time" ,23000, 0, 23000, 400, -10., 10.),
    time_clusters_TAPS("time_clusters_TAPS", "TAPS cluster time" ,200, -50., 50., 450, 0, 450),
    time_clusters_CB("time_clusters_CB", "CB cluster time", 200, -50., 50., 720, 0, 720),
    time_nr_AllClusters("time_nr_AllClusters", "nr of all detected clusters",15, 0, 15),
    time_nr_ClustersinTime("time_nr_ClustersinTime", "nr of detected clusters in time window",15, 0, 15),
    time_nr_CltimeVeto("CltimeVeto", "nr of detected clusters in time window and not e+/-",15, 0, 15),
    six_time_TaggedTime("six_time_TaggedTime", "Tagger time for the events with 7 clusters", 2000, -500, 500),
    detnr(6),
    CB_region(4),
    Is_CB_6g(7),
    photons_rec(nPhotons_six),
    photons_fit(nPhotons_six),
    photons_fit_eta(nPhotons_six),
    kinfit4g("4greactions"),
    kinfit("etaprime"),
    kinfit_eta("etaprime6g_eta_fix"),
    kinfit3pi("3pi0hyp"),
    kinfiteta2pi("eta2pihyp"),
    kinfit10g("etaprime10g"),
    Photons_four(nPhotons_four),
    Photons_six(nPhotons_six),
    Photons_six_3pi(nPhotons_six),
    Photons_six_eta2pi(nPhotons_six),
    Photons_six_eta(nPhotons_six),
    Photons_ten(nPhotons_ten)
{


// TRUE OBSERVABLES
// Beam energy
    true_BeamE                  = new GH1("true_BeamE", "True Beam Energy", 250, 1.400, 1.650);
    tag_BeamE                   = new GH1("tag_BeamE", "Tagged Beam Energy", 250, 1400, 1650);
// Phase space final state particles
    true_th_p_v_th_etapr_CM     = new GHistBGSub2("true_th_p_v_th_etapr_CM", "#theta_{p} vs #theta_{#eta^{'}}", 360, 0., 180, 360, 0., 180.);
    true_th_v_E_p               = new GHistBGSub2("true_th_v_E_p", "True E_{p} vs #theta_{p}", 100, 0., 0.6, 100, 0., 25.);
    true_th_v_E_eta_g           = new GHistBGSub2("true_th_v_E_eta_g", "E_{#gamma, #eta} vs #theta_{#gamma, #eta}", 100, 0, 1.000, 36, 0, 180);
    true_th_v_E_pi0_g           = new GHistBGSub2("true_th_v_E_pi0_g", "E_{#gamma, #pi^{0}} vs #theta_{#gamma, #pi^{0}}", 60, 0, .600, 36, 0, 180);
//  Kinfit tests
    true_six_dth_vs_th_p        = new GHistBGSub2("true_six_dth_vs_th_p", "proton; #theta_{p,rec}; #theta_{rec}-#theta_{true} (^{o})", 200, 0, 25, 80, -5., 5.);
    true_six_z_v_Ncl            = new GHistBGSub2("true_six_z_v_Ncl", "Nr rec clusters vs z_{true}; z_{true} (cm); z_{true} - z_{fit} (cm)", 100, -10., 10., 15, 0, 15);
    true_six_fit_dz_v_z         = new GHistBGSub2("six_fit_dz_v_z", "#Delta z vs z_{true}; z_{true} (cm); z_{true} - z_{fit} (cm)", 100, -10., 10., 50, -10., 10.);

// Physics result
    true_DP                     = new TH2D("true_DP", "True Dalitz Plot distribution", 600, -1.5, 1.5, 600, -1.5, 1.5);
    true_phy_DP                 = new TH1D("true_phy_DP", "True Dalitz Plot distribution as bin nr 6#gamma", 800, 0, 800);
    true_M_pi1pi2_e2p           = new TH1D("true_M_pi1pi2_e2p", "True M_{#pi#pi,true}^{2}", 100 , 0.0, 200);
    true_M_etapi_e2p            = new TH1D("true_M_etapi_e2p", "True M_{#eta#pi,true}^{2}", 350, 0.0, 700.);
// In 6g analysis
    true_six_phy_dMpipi_v_Mpipi = new GHistBGSub2("true_six_phy_dMpipi_v_Mpipi", "fitted - true value M_{#pi#pi,fit}^{2}", 200, 0.0, 200, 200, -100, 100);
    true_six_phy_dX_v_DPbin     = new GHistBGSub2("true_six_phy_dX_v_DPbin", "X_{fit} - X_{true} vd DP bin nr", 800, 0, 800, 200, -2.0, 2.0);
    true_six_phy_dY_v_DPbin     = new GHistBGSub2("true_six_phy_dY_v_DPbin", "Y_{fit} - Y_{true} vd DP bin nr", 800, 0, 800, 200, -2.0, 2.0);
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
    p_E_v_TOF_All_wVeto         = new GHistBGSub2("p_E_v_TOF_All_wVeto", "Energy TAPS vs TOF w VETO > 1MeV", 800, -20., 20., 800, 0., 800.);
    p_E_v_TOF_norm_to_c         = new GHistBGSub2("p_E_v_TOF_norm_to_c", "Energy TAPS vs TOF norm to c", 800, -20., 20., 800, 0., 800.);

    p_nr                        = new GHistBGSub("p_nr", "nr of protons", 5,  0, 5);
    p_th_v_E                    = new GHistBGSub2("p_th_v_E", "Rec E_{p} vs #theta_{p}", 120, 0., 600, 50, 0., 25.);
    p_MM                        = new GHistBGSub("MM_p", "Missing Mass calculated for proton", 300, 800., 1100.);

    IMgg_v_det_2pi0_CB          =   new GHistBGSub2("IMgg_v_det_2pi0_CB", "IM(gg) 2#pi^{0}, CB", 500, 0, 1000, 720, 0, 720);
    IMgg_v_det_etapi0_pi0_CB        =   new GHistBGSub2("IMgg_v_det_etapi0_pi0_CB", "IM(gg) #eta#pi^{0}, CB", 400, 0, 400, 720, 0, 720);
    IMgg_v_det_etapi0_eta_CB        =   new GHistBGSub2("IMgg_v_det_etapi0_eta_CB", "IM(gg) #eta#pi^{0}, CB", 400, 400, 800, 720, 0, 720);
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


    four_fit_Pulls_g_E_vs_E_CB       = new GHistBGSub2("four_fit_Pulls_g_E_vs_E_CB", "4g Pulls #gamma E vs E CB", 100, 0, 1000, 50, -5., 5.);
    four_fit_Pulls_g_th_vs_th_CB     = new GHistBGSub2("four_fit_Pulls_g_th_vs_th_CB", "4g Pulls #gamma #theta vs #theta CB", 90, 0, 180, 50, -5., 5.);
    four_fit_Pulls_g_phi_vs_th_CB    = new GHistBGSub2("four_fit_Pulls_g_phi_vs_th_CB", "4g Pulls #gamma #phi vs #theta CB", 90, 0, 180., 50, -5., 5.);
    four_fit_Pulls_g_E_vs_E_TAPS       = new GHistBGSub2("four_fit_Pulls_g_E_vs_E_TAPS", "4g Pulls #gamma E vs E TAPS", 100, 0, 1000, 50, -5., 5.);
    four_fit_Pulls_g_th_vs_th_TAPS     = new GHistBGSub2("four_fit_Pulls_g_th_vs_th_TAPS", "4g Pulls #gamma #theta vs #theta TAPS", 25, 0, 50, 50, -5., 5.);
    four_fit_Pulls_g_phi_vs_th_TAPS    = new GHistBGSub2("four_fit_Pulls_g_phi_vs_th_TAPS", "4g Pulls #gamma #phi vs #theta TAPS", 25, 0, 50, 50, -5., 5.);


    four_fit_best_2pi_IM_v_E        = new GHistBGSub2("four_fit_best_2pi_IM_v_E", "E_{#pi^{0}} vs M_{#gamma#gamma} for 2#pi^{0}", 100, 0., 1000., 100, 0., 400.);
    four_fit_best_etapi_IM_v_E      = new GHistBGSub2("four_fit_best_etapi_IM_v_E", "E_{#pi^{0}} vs M_{#gamma#gamma} for #eta#pi^{0}", 100, 0., 1000., 250, 0., 1000.);

    four_fit_re_v_eth              = new GHistBGSub2("four_fit_re_v_eth", "#DeltaE vs E, #theta", 5000, 0, 5000, 200, -1.0, 1.0);
    four_fit_mgg_v_eth             = new GHistBGSub2("four_fit_mgg_v_eth", "m_{#gamma#gamma} vs E, #theta", 5000, 0, 5000, 200, 0, 1000);
    four_fit_mgg_v_edet            = new GHistBGSub2("four_fit_mgg_v_edet", "m_{#gamma#gamma} vs E, #detnr", 5000, 0, 5000, 200, 0, 1000);
    four_fit_mgg_v_edet_TAPS       = new GHistBGSub2("four_fit_mgg_v_edet_TAPS", "m_{#gamma#gamma} vs E, #detnr, TAPS", 2500, 0, 2500, 200, 0, 1000);


    // Kinfit related variables 6g

    test_six_rec_IM             = new GH1("test_six_rec_IM", " rec. IM(6#gamma)", 280, 0., 1400.);
    test_six_fit_chi2           = new GH1("test_six_fit_chi2", "#chi^{2} kinfit for 6#gamma", 500, 0, 500.);
    test_six_fit_pdf            = new GH1("test_six_fit_pdf", "pdf kinfit for 6#gamma", 100, 0, 1.);
    test_six_fit_Pulls          = new GHistBGSub2("test_six_fit_Pulls", "Pulls 6#gamma", 50, -5., 5., 25, 0, 25);

    six_rec_IM                  = new GH1("six_rec_IM", " rec. IM(6#gamma)", 280, 0., 1400.);
    six_rec_IM_v_MMp            = new GHistBGSub2("six_rec_IM_v_MMp", "MM(p) vs IM(6#gamma)",240, 200., 1200., 240, 200., 1400.);

    six_fit_chi2                = new GH1("six_fit_chi2", "#chi^{2} kinfit for 6#gamma", 100, 0, 100.);
    six_fit_pdf                 = new GH1("six_fit_pdf", "pdf kinfit for 6#gamma", 100, 0, 1.);
    six_fit_eta_pdf             = new GH1("six_fit_eta_pdf", "pdf kinfit with eta mass enforced 6#gamma", 100, 0, 1.);
    six_fit_Pulls               = new GHistBGSub2("six_fit_Pulls", "Pulls 6#gamma", 50, -5., 5., 25, 0, 25);

    six_fit_Pulls_g_E_vs_E_CB       = new GHistBGSub2("six_fit_Pulls_g_E_vs_E_CB", "Pulls #gamma E vs E CB", 100, 0, 1000, 50, -5., 5.);
    six_fit_Pulls_g_th_vs_th_CB     = new GHistBGSub2("six_fit_Pulls_g_th_vs_th_CB", "Pulls #gamma #theta vs #theta CB", 90, 0, 180, 50, -5., 5.);
    six_fit_Pulls_g_phi_vs_th_CB    = new GHistBGSub2("six_fit_Pulls_g_phi_vs_th_CB", "Pulls #gamma #phi vs #theta CB", 90, 0, 180., 50, -5., 5.);
    six_fit_Pulls_g_E_vs_E_TAPS       = new GHistBGSub2("six_fit_Pulls_g_E_vs_E_TAPS", "Pulls #gamma E vs E TAPS", 100, 0, 1000, 50, -5., 5.);
    six_fit_Pulls_g_th_vs_th_TAPS     = new GHistBGSub2("six_fit_Pulls_g_th_vs_th_TAPS", "Pulls #gamma #theta vs #theta TAPS", 50, 0, 50, 50, -5., 5.);
    six_fit_Pulls_g_phi_vs_th_TAPS    = new GHistBGSub2("six_fit_Pulls_g_phi_vs_th_TAPS", "Pulls #gamma #phi vs #theta TAPS", 50, 0, 50, 50, -5., 5.);

    six_fit_Pulls_p_th_vs_th_TAPS     = new GHistBGSub2("six_fit_Pulls_p_th_vs_th_TAPS", "Pulls proton #theta vs #theta TAPS", 25, 0, 25, 50, -5., 5.);
    six_fit_Pulls_p_phi_vs_th_TAPS    = new GHistBGSub2("six_fit_Pulls_p_phi_vs_th_TAPS", "Pulls proton #phi vs #theta TAPS", 25, 0, 25, 50, -5., 5.);

    six_fit_re_v_eth                = new GHistBGSub2("six_fit_re_v_eth", "#DeltaE vs E, #theta", 5000, 0, 5000, 200, -1.0, 1.0);
    six_fit_mgg_v_eth               = new GHistBGSub2("six_fit_mgg_v_eth", "m_{#gamma#gamma} vs E, #theta", 5000, 0, 5000, 200, 0., 1000.);
    six_fit_mgg_v_edet               = new GHistBGSub2("six_fit_mgg_v_edet", "m_{#gamma#gamma} vs E, det mod", 5000, 0, 5000, 200, 0., 1000.);
    six_fit_mgg_v_edet_TAPS       = new GHistBGSub2("six_fit_mgg_v_edet_TAPS", "m_{#gamma#gamma} vs E, #detnr, TAPS", 2500, 0, 2500, 200, 0, 1000);

    n_fit_mgg_v_e              = new GHistBGSub2("n_fit_mgg_v_e", "m_{#gamma#gamma} vs E", 100, 0, 1000, 200, 0, 1000);
    n_fit_mgg_v_eth              = new GHistBGSub2("n_fit_mgg_v_eth", "m_{#gamma#gamma} vs E, #theta", 5000, 0, 5000, 200, 0, 1000);
    n_fit_mgg_v_edet              = new GHistBGSub2("n_fit_mgg_v_edet", "m_{#gamma#gamma} vs E, #detnr", 5000, 0, 5000, 200, 0, 1000);

    n_fit_mgg_v_e_TAPS       = new GHistBGSub2("n_fit_mgg_v_e_TAPS", "m_{#gamma#gamma} vs E, TAPS", 100, 0, 1000, 200, 0, 1000);
    n_fit_mgg_v_edet_TAPS       = new GHistBGSub2("n_fit_mgg_v_edet_TAPS", "m_{#gamma#gamma} vs E, #detnr TAPS", 2500, 0, 2500, 200, 0, 1000);


    six_fit_IM                  = new GH1("six_fit_IM", "IM(6#gamma) after APLCON fit", 500, 400., 1400.);
    six_fit_IM_rec              = new GH1("six_fit_IM_rec", "IM(6#gamma) reconstructed values for events which pass fit", 500, 400., 1400.);
    six_fit_IM_3pi              = new GH1("six_fit_IM_3pi", "IM(6#gamma) for 3#pi^{0} candidates", 500, 400., 1400.);
    six_fit_IM_eta2pi           = new GH1("six_fit_IM_eta2pi", "IM(6#gamma) for #eta2#pi^{0} candidates", 500, 400., 1400.);
    six_fiteta_IM2g             = new GH1("six_fiteta_IM2g", "IM(2g)_{#eta} where #eta mass enforced", 500, 500, 600.);

    six_fit_PDF_eta2pi_v_3pi    = new GHistBGSub2("six_fit_PDF_eta2pi_v_3pi", "PDF_eta2pi_v_3pi 6#gamma", 100, 0., 1., 100, 0., 1.);
    six_fit_PDF_eta2pi_v_3pi_2  = new GHistBGSub2("six_fit_PDF_eta2pi_v_3pi_2", "PDF_eta2pi_v_3pi 6#gamma using kinfit", 100, 0., 1., 100, 0., 1.);
    six_fit_best_eta            = new GH1("six_fit_best_eta", "best #eta cand from comb", 500, 200, 700.);
    six_fit_best_eta_E_v_th     = new GHistBGSub2("six_fit_best_eta_E_v_th", "E_{#eta, #gamma} vs #Theta_{#eta, #gamma}", 100, 0., 1000., 50, 0., 200.);
    six_fit_best_2pi            = new GH1("six_fit_best_2pi", "best 2#pi^{0} cand from comb", 500, 0., 500.);  

    six_fit_PDF_eta2pi_v_Meta2pi    = new GHistBGSub2("six_fit_PDF_eta2pi_v_Meta2pi", "PDF_eta2pi_sel vs IM(6#gamma), 6g cand", 100, 0., 1., 100, 600., 1100.);
    six_fit_PDF_2_eta2pi_v_Meta2pi  = new GHistBGSub2("six_fit_PDF_2_eta2pi_v_Meta2pi", "PDF_eta2pi kinfit vs IM(6#gamma), 6g cand", 100, 0., 1., 100, 600., 1100.);
    six_fit_eta_PDF_v_Metapr        = new GHistBGSub2("six_fit_eta_PDF_v_Metapr", "PDF when eta mass enforced vs IM, 6#gamma cand", 100, 0., 1., 100, 600., 1100.);

    six_phy_etapr_v_BeamE       = new GHistBGSub2("six_phy_etapr_v_BeamE", "IM(6#gamma) vs Beam Energy", 20, 1400, 1600, 500, 600., 1100.);
    six_phy_etapr_eta_v_BeamE   = new GHistBGSub2("six_phy_etapr_eta_v_BeamE", "IM(6#gamma) with enforced eta mass vs Beam Energy", 20, 1400, 1600, 500, 600., 1100.);

    six_phy_DP                  = new GHistBGSub2("six_phy_DP", "Rec fitted Dalitz Plot distribution vs bin nr 6#gamma", 800, 0, 800, 500, 600., 1100.);
    six_phy_M_pi1pi2_v_etapr    = new GHistBGSub2("six_phy_M_pi1pi2_v_etapr", "Rec M_{#pi#pi,fit}^{2} 6#gamma", 100 , 0.0, 200, 100, 600., 1100. );
    six_phy_M_etapi_v_etapr     = new GHistBGSub2("six_phy_M_etapi_v_etapr", "Rec M_{#eta#pi,fit}^{2} 6#gamma", 350 , 0.0, 700, 100, 600., 1100. );

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

    IM10g_fit           = new GH1("IM10g_fit", "IM(10#gamma) after APLCON fit", 500, 400., 1400.);

// Physics results eta'

    cutFile             = new TFile("configfiles/cuts/TAPSEdEPA.root");
    cutProtonTAPS       = (TCutG*)cutFile->Get("CutProton");

    cutFile2            = new TFile("configfiles/cuts/TAPS_TOF_PA3.root");
    cutProtonETOF       = (TCutG*)cutFile2->Get("Hadron");

    cutFile3            = new TFile("configfiles/cuts/minion_cut.root");
    cutMinIon           = (TCutG*)cutFile3->Get("CutElectron");


    // APLCON kinfit uncertainties
    g_unc               = new TFile("configfiles/APLCONunc/photon_uncertainties.root");
    p_unc               = new TFile("configfiles/APLCONunc/proton_uncertainties.root");

    g_CB_e = (TH2F*)g_unc->Get("CB_E")->Clone();
    g_CB_th = (TH2F*)g_unc->Get("CB_theta")->Clone();
    g_CB_fi = (TH2F*)g_unc->Get("CB_phi")->Clone();
    g_TAPS_e = (TH2F*)g_unc->Get("TAPS_E")->Clone();
    g_TAPS_th = (TH2F*)g_unc->Get("TAPS_theta")->Clone();
    g_TAPS_fi = (TH2F*)g_unc->Get("TAPS_phi")->Clone();

    p_TAPS_e = (TH2F*)p_unc->Get("TAPS_E")->Clone();
    p_TAPS_th = (TH2F*)p_unc->Get("TAPS_theta")->Clone();
    p_TAPS_fi = (TH2F*)p_unc->Get("TAPS_phi")->Clone();

    thcorr_CB                   = new TFile("configfiles/corr/CB_th_corr.root");
    dthvth_CB                   = (TProfile*)thcorr_CB->Get("photon_dtheta_v_theta_CB_pfx")->Clone();

    Ecorr_CB                    = new TFile("configfiles/corr/CB_e_corr.root");
    EvdetCB                     = (TH2F*)Ecorr_CB->Get("g_peak_E_CB");

    thcorr_TAPS                 = new TFile("configfiles/corr/TAPS_th_corr.root");
    dthvth_TAPS                 = (TProfile*)thcorr_TAPS->Get("photon_dtheta_v_theta_TAPS_pfx")->Clone();

    GHistBGSub::InitCuts(-10., 10., -60., -20.);
    GHistBGSub::AddRandCut(20., 60.);

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


 //  For final states including 6g and 2g is = eta mass
        kinfit_eta.LinkVariable("Beam_eta",    beam_eta.Link(),   beam_eta.LinkSigma(),  beam_eta.LinkSettings() );
        kinfit_eta.LinkVariable("Proton_eta",  proton_eta.Link(), proton_eta.LinkSigma());

        vector<string> photon_names_eta;
        for(size_t i=0;i<nPhotons_six;i++) {
            stringstream s_photon_eta;
            s_photon_eta << "Photon_eta" << (i+1);
            photon_names_eta.push_back(s_photon_eta.str());
            kinfit_eta.LinkVariable(s_photon_eta.str(), Photons_six_eta[i].Link(), Photons_six_eta[i].LinkSigma());
        }
        vector<string> all_names_eta = {"Beam_eta", "Proton_eta"};
        all_names_eta.insert(all_names_eta.end(),photon_names_eta.begin(),photon_names_eta.end());


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
    kinfit_eta.AddConstraint("EnergyMomentumBalance", all_names_eta, EnergyMomentumBalance);


    kinfiteta2pi.AddConstraint("EnergyMomentumBalance", all_names_eta2pi, EnergyMomentumBalance);
    kinfit3pi.AddConstraint("EnergyMomentumBalance", all_names_3pi, EnergyMomentumBalance);

    kinfit4g.AddConstraint("EnergyMomentumBalance", all_names4g, EnergyMomentumBalance);
    kinfit10g.AddConstraint("EnergyMomentumBalance", all_names10g, EnergyMomentumBalance);


    // Constraint: Invariant mass of nPhotons equals constant IM,
    // make lambda catch also this with [&] specification
    auto RequireIM = [&] (const vector< vector<double> >& photons) -> double
    {
        TLorentzVector sum(0,0,0,0);
        for(int i = 0; i < 2; i++) {
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

    kinfit3pi.AddConstraint("RequireIMpi1", photon_names_3pi, RequireIMpi1); //added req that first pair has pi0 mass.
    kinfit3pi.AddConstraint("RequireIMpi2", photon_names_3pi, RequireIMpi2); //added req that second pair has pi0 mass.
    kinfit3pi.AddConstraint("RequireIMpi3", photon_names_3pi, RequireIMpi3); //added req that third pair has pi0 mass.

    kinfiteta2pi.AddConstraint("RequireIM", photon_names_eta2pi, RequireIM); //added req that eta' candidate has eta mass.
    kinfiteta2pi.AddConstraint("RequireIMpi2", photon_names_eta2pi, RequireIMpi2); //added req that second pair has pi0 mass.
    kinfiteta2pi.AddConstraint("RequireIMpi3", photon_names_eta2pi, RequireIMpi3); //added req that third pair has pi0 mass.

    kinfit_eta.AddConstraint("RequireIM", photon_names_eta, RequireIM); //added req that eta' candidate has eta mass.
    kinfit_eta.AddConstraint("RequireIMpi2", photon_names_eta, RequireIMpi2); //added req that second pair has pi0 mass.
    kinfit_eta.AddConstraint("RequireIMpi3", photon_names_eta, RequireIMpi3); //added req that third pair has pi0 mass.


    auto VertexConstraint = [&] (vector< vector<double> >& args) -> vector<double>
    {
       TLorentzVector diff(0.0, 0.0, 0.0, 0.0);
       const TLorentzVector target(0,0,0, MASS_PROTON);
       // assume first particle is beam photon
       diff = target + FitParticle::Make(args[0], 0.0 ); // beam photon

       double R;
       constexpr double R_CB = 25.4;
       constexpr double R_TAPS = 145.7;

       const double v_z = args.back()[0];
       args.resize(args.size()-1); // get rid of last element
       // correct each final state particle theta angle,

       for( UInt_t i = 1; i < 8; i++)
       {
           double theta = args[i][1]; // 2nd element theta

           if( !(Is_CB_6g[i-1]) )
           {
               R = R_TAPS/cos(theta);
            //   R= R_CB;
           }
           else
               R = R_CB;

           double theta_p = std::atan2( R*sin(theta), R*cos(theta) - v_z);
           args[i][1] = theta_p;                // now set the re-defined theta angle as the new theta
           if(i == 1)
               diff -= FitParticle::Make(args[i], MASS_PROTON );
           else
               diff -= FitParticle::Make(args[i], 0.0 );
       }

       return {diff.X(), diff.Y(), diff.Z(), diff.T()};

    };
//    kinfit.AddUnmeasuredVariable("v_z"); // default value 0
//    kinfit.AddConstraint("VertexConstraint", all_names + std::vector<string>{"v_z"}, VertexConstraint);


    APLCON::Fit_Settings_t settings = kinfit.GetSettings();
    settings.MaxIterations = 20;

//    settings.DebugLevel = 5;
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

    time_TOF.Reset();
    time_clusters_TAPS.Reset();
    time_clusters_CB.Reset();
    time_nr_AllClusters.Reset();
    time_nr_ClustersinTime.Reset();
    time_nr_CltimeVeto.Reset();
    six_time_TaggedTime.Reset();

    TraverseValidEvents();

    outputFile->cd();

    time_TOF.Write();
    time_clusters_TAPS.Write();
    time_clusters_CB.Write();
    time_nr_AllClusters.Write();
    time_nr_ClustersinTime.Write();
    time_nr_CltimeVeto.Write();
    six_time_TaggedTime.Write();

	return kTRUE;
}

void	AdlarsonPhysics::ProcessEvent()
{
//    if(useTrueObservables)
//    {
//        std::cout << "useTrueObservables = " << useTrueObservables << std::endl;
//    }
//    etapr_6gTrue.Start(*GetPluto(), *GetGeant());   // (pluto tree, n part in pluto per event)
//    TrueAnalysis_etapr6g();                         // obtains the true observables

    Energy_corr();
    theta_corr();
//    Kinfit_test();                                  // runs kinematical fit with true observables- for testing purposes


    ClustersInTime.clear();
    ClustersInTime.resize(0);
    UInt_t InTime = 0;
    UInt_t AllTimes = 0;
    UInt_t InTime_minion = 0;

    for( Int_t i = 0; i < GetTracks()->GetNTracks(); i++ )
    {
        if( GetTracks()->HasTAPS(i) )
        {
            AllTimes++;
            time_clusters_TAPS.Fill( GetTracks()->GetTime(i), GetTracks()->GetCentralCrystal(i) );
            if( ( GetTracks()->GetTime(i) > -8.0) &&  ( GetTracks()->GetTime(i) < 15.0 ) )
            {

                InTime++;
                if(cutMinIon->IsInside(GetTracks()->GetClusterEnergy(i), GetTracks()->GetVetoEnergy(i)))
                {
                   p_E_v_dE_cc->Fill(GetTracks()->GetClusterEnergy(i),GetTracks()->GetVetoEnergy(i));
                }
                else
                {
                    if(GetTracks()->GetClusterEnergy(i) > 20.0)
                    {
                        InTime_minion++;
                        ClustersInTime.push_back( i );
                    }
                }
            }
        }
        if( GetTracks()->HasCB(i) )
        {
            AllTimes++;
            time_clusters_CB.Fill( GetTracks()->GetTime(i), GetTracks()->GetCentralCrystal(i) );
            if( TMath::Abs( GetTracks()->GetTime(i) ) < 30.0 )
            {
                if(GetTracks()->GetClusterEnergy(i) > 0.0)
                {
                    ClustersInTime.push_back( i );
                    InTime++;
                    InTime_minion++;
                }
            }
        }
    }
    time_nr_AllClusters.Fill(AllTimes);
    time_nr_ClustersinTime.Fill(InTime);
    time_nr_CltimeVeto.Fill(InTime_minion);

    MMp_vec.SetPxPyPzE(0., 0., 0., 0.);

    nrprotons = 0;
    iprtrack = 10000;
    Double_t TOF = -100.0;
    for(int tag = 0; tag < GetTagger()->GetNTagged(); tag++){
        for(UInt_t i = 0; i < ClustersInTime.size() ; i++){
            UInt_t j = ClustersInTime[i];
            if( GetTracks()->HasTAPS(j) ){
                Double_t radnm = 1.45/TMath::Cos( GetTracks()->GetThetaRad(j) );
                int nbin = GetTagger()->GetTaggedChannel(tag) + GetTracks()->GetCentralCrystal(j)*50;
                TOF = (( GetTagger()->GetTaggedTime(tag) - GetTracks()->GetTime(j))/radnm)-TAPS_EPT_toff[nbin];
                if(TMath::Abs(TOF) < 10.1){
                    int nbin = GetTagger()->GetTaggedChannel(tag) + GetTracks()->GetCentralCrystal(j)*50;
                    time_TOF.Fill(nbin, TOF);
                }
                p_E_v_TOF_All->Fill( TOF , GetTracks()->GetClusterEnergy(j));
                if( cutProtonETOF->IsInside( TOF, GetTracks()->GetClusterEnergy(j)) ){
                    p_E_v_TOF->Fill( TOF, GetTracks()->GetClusterEnergy(j));
                    // uncomment for EXP
                    nrprotons++;
                    if(nrprotons >= 1)
                    {
                        if( GetTracks()->GetVetoEnergy(j) > GetTracks()->GetVetoEnergy(iprtrack))
                            iprtrack = j;
                    }
                    else
                        iprtrack = j;
                }

                if(  GetTracks()->GetVetoEnergy(j) > ( 4.5 -(4.5/180)*GetTracks()->GetClusterEnergy(j)) )
                {
                    if(GetTracks()->GetVetoEnergy(j) > 1.0)
                    {
                        p_E_v_TOF_All_wVeto->Fill(TOF, GetTracks()->GetClusterEnergy(j));
                    }
                }

           }
        }
    }


    // test TOF in terms of ns per meter
    if(iprtrack != 10000)
    {
        Double_t TOF_corr = -1.0e6;
        Double_t offset = 0.625;
        Double_t c_ns_per_m = 0.2998;
        Double_t radnm = 1.45/TMath::Cos( GetTracks()->GetThetaRad(iprtrack) );
        for(int tag = 0; tag < GetTagger()->GetNTagged(); tag++)
        {

            TOF_corr = (GetTagger()->GetTaggedTime(tag) - GetTracks()->GetTime(iprtrack))/radnm + offset - c_ns_per_m;
            p_E_v_TOF_norm_to_c->Fill(-TOF_corr, GetTracks()->GetClusterEnergy(iprtrack),GetTagger()->GetTaggedTime(tag));
        }
    }


//    for (UInt_t i = 0; i < ClustersInTime.size() ; i++)
    for (UInt_t ii = 0; ii < GetTracks()->GetNTracks(); ii++)
    {
        if( GetTracks()->HasTAPS(ii) )
        {
            p_E_v_dE_all->Fill(GetTracks()->GetClusterEnergy(ii),GetTracks()->GetVetoEnergy(ii));
            if( cutProtonTAPS->IsInside(GetTracks()->GetClusterEnergy(ii), GetTracks()->GetVetoEnergy(ii)))
            {
//                nrprotons++;
//                if(nrprotons >= 1)
//                {
//                    if( GetTracks()->GetVetoEnergy(ii) > GetTracks()->GetVetoEnergy(iprtrack))
//                        iprtrack =ii;
//                }
//                else
//                    iprtrack = ii;

                p_E_v_dE_pr->Fill(GetTracks()->GetClusterEnergy(ii),GetTracks()->GetVetoEnergy(ii));

            }

        }
    }



    p_nr->Fill(nrprotons);
    if( iprtrack == 10000 ) return;
    if( nrprotons == 0 ) return;

    true_six_dth_vs_th_p->Fill(GetTracks()->GetTheta(iprtrack),GetTracks()->GetTheta(iprtrack) - etapr_6gTrue.GetTrueProtonLV().Theta()*TMath::RadToDeg());

    p_E_v_dE_pr->Fill(GetTracks()->GetClusterEnergy(iprtrack),GetTracks()->GetVetoEnergy(iprtrack));
    p_th_v_E->Fill(GetTracks()->GetClusterEnergy(iprtrack),GetTracks()->GetTheta(iprtrack));
    proton_vec = GetTracks()->GetVector(iprtrack, pdgDB->GetParticle("proton")->Mass()*1000);
    // Now construct missing mass calc for proton with tagger energies.
    for ( Int_t tag = 0; tag < GetTagger()->GetNTagged(); tag++)
    {
        if( TMath::Abs(GetTagger()->GetTaggedTime(tag)) > taggerTimeCut ) continue;
        tag_BeamE->Fill(GetTagger()->GetTaggedEnergy(tag),GetTagger()->GetTaggedTime(tag));
        MMp_vec.SetPxPyPzE(-proton_vec.Px(), -proton_vec.Py(), GetTagger()->GetTaggedEnergy(tag) - proton_vec.Pz(), GetTagger()->GetTaggedEnergy(tag) + pdgDB->GetParticle("proton")->Mass()*1000 - proton_vec.E());

        MMp = ( MMp_vec ).M();
        p_MM->Fill(MMp, GetTagger()->GetTaggedTime(tag));
    }



    Bool_t good_multiplicity = 0;
    if( InTime == 3 )
        good_multiplicity = 1;
    if( InTime == 5 )
        good_multiplicity = 1;
    if( (InTime > 6) && (InTime < 9) )
        good_multiplicity = 1;
    if( InTime == 11 )
        good_multiplicity = 1;

    if(!good_multiplicity) return;

    // do for all clusters energy check. If Energy of cluster less than x MeV reject?


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

    if(GetTracks()->GetTheta(iprtrack) > 18.0 )
    {
        true_six_z_v_Ncl->Fill(etapr_6gTrue.GetTrueVertex().Z(), InTime_minion);
    }

    if( (ClustersInTime.size() == 5) && (nrprotons > 0) )
    {
        fourgAnalysis(iprtrack);
    }

    if( ClustersInTime.size() == 7 && (nrprotons > 0) )
    {
        sixgAnalysis( iprtrack );
    }

    if( ClustersInTime.size() == 11 && (nrprotons > 0) )
    {
        tengAnalysis(iprtrack);
    }
}

void	AdlarsonPhysics::ProcessScalerRead()
{
//    hist_eta.ScalerReadCorrection(Double_t(scalers->GetScaler(0))/scalers->GetScaler(1));
//    hist_etap.ScalerReadCorrection(Double_t(scalers->GetScaler(0))/scalers->GetScaler(1));
}

Bool_t	AdlarsonPhysics::Init(const char* configFile)
{

   std::string         line;

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


   std::ifstream file1("configfiles/data/CB_gain_vs_erec_mod_2.txt");
   if(file1){
   std::getline(file1, line); // first line describes input


   while(std::getline(file1, line)){
       std::string         buffer;
       std::stringstream   ss;
       Int_t  module;
       std::vector<Double_t> vec;
       ss << line;
       std::getline(ss, buffer, '\t');
       module = atoi(buffer.c_str());
       while (std::getline(ss, buffer, '\t')) {
            vec.push_back(std::stod(buffer));
       }
       CB_Ecorr.insert(corr_pair(module, vec));
   }

   file1.close();
   }


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



//   typedef std::pair<UInt_t, std::vector<Double_t>> EPT_TAPS_pair;
//   std::map<UInt_t, std::vector<Double_t>> TOF_corr;


      std::cout << " CB gain correction applied" << std::endl;
      std::cout << " TAPS gain correction applied" << std::endl;
//    std::cout << "No Init function specified for this class." << std::endl;

    return kTRUE;

}


void AdlarsonPhysics::TrueAnalysis_etapr6g()
{
    true_BeamE->Fill( etapr_6gTrue.GetTrueBeamEnergy() );

    // calculate Physics

    etapr_true[0] = etapr_6gTrue.GetTrueEtaLV();
    etapr_true[1] = etapr_6gTrue.GetTrueNeutralPiLV(0);
    etapr_true[2] = etapr_6gTrue.GetTrueNeutralPiLV(1);
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
            true_M_etapi_e2p->Fill(etapr_6gTrue.GetTrueGammaLV(i).E(),etapr_6gTrue.GetTrueGammaLV(i).Theta()*TMath::RadToDeg());
        else
            true_M_pi1pi2_e2p->Fill(etapr_6gTrue.GetTrueGammaLV(i).E(),etapr_6gTrue.GetTrueGammaLV(i).Theta()*TMath::RadToDeg());
    }
    return;
}

void AdlarsonPhysics::fourgAnalysis(UInt_t ipr)
{
    for(int tag = 0; tag < GetTagger()->GetNTagged(); tag++)
    {
        Double_t chi2_min = 1.0e6;
        Double_t prob_min = 1.0e6;
        std::vector<Int_t> set;     // particle conf: tagger, proton, photon 1, 2, 3, 4, 5, 6;
        detnr.resize(0);
        CB_region.resize(0);

        IM4g_vec.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
        IM4g_fit.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);

        for(UInt_t i = 0; i < nPhotons_four; i++){
            photons_rec[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
            photons_fit[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
        }

        // beam energy
        set.resize(0);
        set.push_back(tag);
        beam4g.SetFromVector( GetTagger()->GetVector(tag) );
        std::vector<double> obs;
        std::vector<double> unc;
//        beam4g.Smear(3, true);
        obs.resize(0);
        unc.resize(0);
        unc = Get_unc(3,0,obs);
        beam4g.Smear_tmp(unc, 2);
        beam4g.APLCONSettings();

        // proton
        set.push_back(ipr);
        proton4g.SetFromVector( GetTracks()->GetVector(ipr) );
//        proton4g.Smear(2, false);
        obs.resize(0);
        obs.push_back(proton4g.Ek);
        obs.push_back(proton4g.Theta);
        unc.resize(0);
        unc = Get_unc( 2, 2, obs);
        proton4g.Smear_tmp(unc, 1);

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
                    Photons_four[n_photons].Smear_tmp(unc, 0);
//                    Photons_four[n_photons].Smear(1, true);
                }
                else{ // ( GetTracks()->HasTAPS(k) )
                    obs.resize(0);
                    obs.push_back(Photons_four[n_photons].Ek);
                    obs.push_back(Photons_four[n_photons].Theta);
                    unc.resize(0);
                    unc = Get_unc( 2, 1, obs);
                    Photons_four[n_photons].Smear_tmp(unc, 0);
//                    Photons_four[n_photons].Smear(2, true);
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

        Double_t sigma_eta = 30.0;
        Double_t sigma_pi0 = 20.0;

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


        // Here put in fcn which finds the best 4g combo...

            //            four_fit_Pulls_g_E_vs_E_CB
            //            four_fit_Pulls_g_th_vs_th_CB
            //            four_fit_Pulls_g_phi_vs_th_CB
            //            four_fit_Pulls_g_E_vs_E_TAPS
            //            four_fit_Pulls_g_th_vs_th_TAPS
            //            four_fit_Pulls_g_phi_vs_th_TAPS

            //    four_fit_best_2pi_IM_v_E
            //            four_fit_best_etapi_IM_v_E



        four_fit_PDF_etapi_v_2pi->Fill(TMath::Prob(chi2min_etapi0,1), TMath::Prob(chi2min_2pi0,1), GetTagger()->GetTaggedTime(tag) );

        if( (chi2min_2pi0 < chi2min_etapi0) && ( TMath::Prob(chi2min_2pi0,1) > 0.001) ){
            for(int j = 0; j < 4; j++){
                int imass = int(j/2);
                idet = perm4g[imin_2pi0][j];
                if(CB_region[j]){
                    double En = (photons_rec[perm4g[ imin_2pi0 ][ j ]]).E();
                    double Th = (photons_rec[perm4g[ imin_2pi0 ][ j]]).Theta()*TMath::RadToDeg();
                    double mgg = (photons_rec[perm4g[ imin_2pi0 ][ imass*2 ]] + photons_rec[perm4g[imin_2pi0][imass*2+1]] ).M();
                    int nEN =  int(En/10.0);
                    int nTH =  100*int(Th/5.0);
                    int nBIN = nEN + nTH;

                    int nDET = 100*(int( detnr[idet]/16 ));
                    int nBIN2 = nEN + nDET;

                    IMgg_v_det_2pi0_CB->Fill( (photons_rec[perm4g[ imin_2pi0 ][ imass*2 ]] + photons_rec[perm4g[imin_2pi0][imass*2+1]] ).M(), detnr[idet], GetTagger()->GetTaggedTime(tag));
                    four_fit_mgg_v_eth->Fill(nBIN, mgg, GetTagger()->GetTaggedTime(tag));
                    four_fit_mgg_v_edet->Fill(nBIN2, mgg, GetTagger()->GetTaggedTime(tag));

                    n_fit_mgg_v_e->Fill(En,mgg, GetTagger()->GetTaggedTime(tag));
                    n_fit_mgg_v_eth->Fill(nBIN, mgg, GetTagger()->GetTaggedTime(tag));
                    n_fit_mgg_v_edet->Fill(nBIN2,mgg, GetTagger()->GetTaggedTime(tag));
                }
                else{
                    IMgg_v_det_2pi0_TAPS->Fill((photons_rec[perm4g[ imin_2pi0 ][ imass*2 ]] + photons_rec[perm4g[imin_2pi0][imass*2+1]] ).M(), detnr[idet], GetTagger()->GetTaggedTime(tag) );
                    double En = (photons_rec[perm4g[ imin_2pi0 ][ j ]]).E();
                    double mgg = (photons_rec[perm4g[ imin_2pi0 ][ imass*2 ]] + photons_rec[perm4g[imin_2pi0][imass*2+1]] ).M();
                    int nEN =  int(En/10.0);
                    int nDET = 100*(int( detnr[idet]/16 ));
                    int nBIN2 = nEN + nDET;

                    four_fit_mgg_v_edet_TAPS->Fill(nBIN2, mgg, GetTagger()->GetTaggedTime(tag));
                    n_fit_mgg_v_edet_TAPS->Fill(nBIN2, mgg, GetTagger()->GetTaggedTime(tag));
                    n_fit_mgg_v_e_TAPS->Fill(En, mgg, GetTagger()->GetTaggedTime(tag));
                }
            }
        }
        else if( (chi2min_etapi0 < chi2min_2pi0  ) && ( TMath::Prob(chi2min_etapi0,1) > 0.001)){   //etapi0
            for(int j = 0; j < 4; j++){
                int imass = int(j/2);
//                idet = perm4g[imin_etapi0][j]-1;
                idet = perm4g[imin_etapi0][j];
                if(CB_region[j]){
                    if(imass == 0)
                        IMgg_v_det_etapi0_eta_CB->Fill( (photons_rec[perm4g[ imin_etapi0 ][ imass*2 ]] + photons_rec[perm4g[imin_etapi0][imass*2+1]] ).M(), detnr[idet], GetTagger()->GetTaggedTime(tag) );
                    else{
                        double En = (photons_rec[perm4g[ imin_etapi0 ][ j ]]).E();
                        double mgg = (photons_rec[perm4g[ imin_etapi0 ][ imass*2 ]] + photons_rec[perm4g[imin_etapi0][imass*2+1]] ).M();
                        int nEN =  int(En/10.0);
                        int nDET = 100*(int( detnr[idet]/16 ));
                        int nBIN2 = nEN + nDET;
                        four_fit_mgg_v_edet->Fill(nBIN2, mgg, GetTagger()->GetTaggedTime(tag));

                        IMgg_v_det_etapi0_pi0_CB->Fill( (photons_rec[perm4g[ imin_etapi0 ][ imass*2 ]] + photons_rec[perm4g[imin_etapi0][imass*2+1]] ).M(), detnr[idet], GetTagger()->GetTaggedTime(tag) );
                    }
                }
                else{
                    IMgg_v_det_etapi0_TAPS->Fill( (photons_rec[perm4g[ imin_etapi0 ][ imass*2 ]] + photons_rec[perm4g[imin_etapi0][imass*2+1]] ).M(), detnr[idet], GetTagger()->GetTaggedTime(tag) );
                    if(imass != 0){
                        double En = (photons_rec[perm4g[ imin_etapi0 ][ j ]]).E();
                        double mgg = (photons_rec[perm4g[ imin_etapi0 ][ imass*2 ]] + photons_rec[perm4g[imin_etapi0][imass*2+1]] ).M();
                        n_fit_mgg_v_e_TAPS->Fill(En, mgg, GetTagger()->GetTaggedTime(tag));
                    }
                }
            }
        }
    }
}

//void AdlarsonPhysics::GetBest6gCombination(Double_t& sigma_eta, Double_t& sigma_pi0, Double_t& chi2min_eta2pi, Double_t& chi2min_3pi, std::vector<int>& imin_eta2pi, std::vector<int>& imin_3pi )


//void AdlarsonPhysics::Best4g_comb(std::vector<TLorentzVector>& photons_rec , std::vector<int>& detnr, std::vector<bool>& CB_region)
//{

//}



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
        Is_CB_6g.resize(0);

        for(UInt_t i = 0; i < nPhotons_six; i++)
        {
            photons_rec[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
            photons_fit[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
            photons_fit_eta[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
        }

        set.push_back(tag);

        beam.SetFromVector( GetTagger()->GetVector(tag) );
//        beam.Smear(3, true);
        obs.resize(0);
        unc.resize(0);
        unc = Get_unc(3,0,obs);
        beam.Smear_tmp(unc, 2);
        beam.APLCONSettings();


        set.push_back(ipr);
        proton.SetFromVector( GetTracks()->GetVector(ipr) );
//        proton.Smear(2, false);
        obs.resize(0);
        obs.push_back(proton.Ek);
        obs.push_back(proton.Theta);
        unc.resize(0);
        unc = Get_unc( 2, 2, obs);
        proton.Smear_tmp(unc, 1);
        Is_CB_6g.push_back(false);

        UInt_t n_photons = 0;
        for ( UInt_t jgam = 0; jgam < ClustersInTime.size() ; jgam++ )
        {
            UInt_t kgam = ClustersInTime[jgam];
            if( kgam != ipr ) // the id proton cluster
            {
                set.push_back(kgam);
                photons_rec[n_photons] = GetTracks()->GetVector(kgam);
                IM6g_vec += photons_rec[n_photons];
                Photons_six[n_photons].SetFromVector( GetTracks()->GetVector(kgam) );

                if( GetTracks()->HasCB(kgam) )
                {
                    obs.resize(0);
                    obs.push_back(Photons_six[n_photons].Ek);
                    obs.push_back(Photons_six[n_photons].Theta);
                    unc.resize(0);
                    unc = Get_unc( 1, 1, obs);
                    Photons_six[n_photons].Smear_tmp(unc, 0);
//                    Photons_six[n_photons].Smear(1, true);
                    Is_CB_6g.push_back(true);
                }
                else
                {// ( GetTracks()->HasTAPS(kgam) )
                    obs.resize(0);
                    obs.push_back(Photons_six[n_photons].Ek);
                    obs.push_back(Photons_six[n_photons].Theta);
                    unc.resize(0);
                    unc = Get_unc( 2, 1, obs);
                    Photons_six[n_photons].Smear_tmp(unc, 0);
//                    Photons_six[n_photons].Smear(2, true);
                    Is_CB_6g.push_back(false);
                }
                n_photons++;
            }
        }

        six_rec_IM->Fill(IM6g_vec.M());
        const APLCON::Result_t& result = kinfit.DoFit();
        if(result.Status == APLCON::Result_Status_t::Success)
        {
            if( result.ChiSquare <  chi2_min )
            {
                chi2_min = result.ChiSquare;
                prob_min = result.Probability;
                set_min = set;
  //              true_six_fit_dz_v_z->Fill(etapr_6gTrue.GetTrueVertex().Z(), etapr_6gTrue.GetTrueVertex().Z() - result.Variables.at("v_z").Value.After);

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
//        beam.Smear(3, true);
        obs.resize(0);
        unc.resize(0);
        unc = Get_unc(3,0,obs);
        beam.Smear_tmp(unc, 2);
        beam.APLCONSettings();

//        proton.SetFromVector( GetTracks()->GetVector(set_min[1]) );
//        proton.Smear(2, false);

        proton.SetFromVector( GetTracks()->GetVector(set_min[1]) );
        obs.resize(0);
        obs.push_back(proton.Ek);
        obs.push_back(proton.Theta);
        unc.resize(0);
        unc = Get_unc( 2, 2, obs);
        proton.Smear_tmp(unc, 1);

        Int_t n_photons_min = 0;
        for ( Int_t jgam_min = 0; jgam_min < 6 ; jgam_min++ )
        {
            photons_rec[n_photons_min] = GetTracks()->GetVector( set_min[jgam_min+2] );
            Photons_six[n_photons_min].SetFromVector( GetTracks()->GetVector(set_min[jgam_min+2]));
            detnr.push_back( GetTracks()->GetCentralCrystal( set_min[jgam_min+2]) );
            IM6g_vec += photons_rec[n_photons_min];

            if( GetTracks()->HasCB(set_min[jgam_min+2]) )
            {

                obs.resize(0);
                obs.push_back(Photons_six[n_photons_min].Ek);
                obs.push_back(Photons_six[n_photons_min].Theta);
                unc.resize(0);
                unc = Get_unc( 1, 1, obs);
                Photons_six[n_photons_min].Smear_tmp(unc, 0);
//                    Photons_six[n_photons_min].Smear(1, true);
                Is_CB_6g.push_back(true);
            }
            else
            {// ( GetTracks()->HasTAPS(jgam_min+2) )
                obs.resize(0);
                obs.push_back(Photons_six[n_photons_min].Ek);
                obs.push_back(Photons_six[n_photons_min].Theta);
                unc.resize(0);
                unc = Get_unc( 2, 1, obs);
                Photons_six[n_photons_min].Smear_tmp(unc, 0);
//                    Photons_six[n_photons_min].Smear(2, true);
                Is_CB_6g.push_back(false);
            }

            n_photons_min++;
        }

        six_rec_IM->Fill(IM6g_vec.M(),GetTagger()->GetTaggedTime(tag));
        six_rec_IM_v_MMp->Fill( MMp, IM6g_vec.M(),GetTagger()->GetTaggedTime(tag) );

        const APLCON::Result_t& result_min = kinfit.DoFit();
        if(result_min.Status == APLCON::Result_Status_t::Success)
        {
           six_fit_chi2->Fill( result_min.ChiSquare, GetTagger()->GetTaggedTime(tag) );
           six_fit_pdf->Fill( result_min.Probability, GetTagger()->GetTaggedTime(tag) );

           Int_t inr = 0;
           for(const auto& it_map : result_min.Variables) {
                const string& varname = it_map.first;
                const APLCON::Result_Variable_t& var = it_map.second;
                six_fit_Pulls->Fill(var.Pull, inr, GetTagger()->GetTaggedTime(tag));
                inr++;
           }

           for(UInt_t iPull = 0; iPull < nPhotons_six; iPull++)
           {
                int nth = 100*int(photons_rec[iPull].Theta()*TMath::RadToDeg()/10);
                int nEn = int((photons_rec[iPull].E()/20));
                int ibin = nEn + nth;

                string s = Form("Photon%0i[0]", iPull+1); // Energy
                string t = Form("Photon%0i[1]", iPull+1); // theta
                string u = Form("Photon%0i[2]", iPull+1); // phi

                if( GetTracks()->HasCB( set[iPull+2]) ){

                    six_fit_Pulls_g_E_vs_E_CB->Fill(photons_rec[iPull].E(), result_min.Variables.at(s).Pull, GetTagger()->GetTaggedTime(tag));
                    six_fit_Pulls_g_th_vs_th_CB->Fill(photons_rec[iPull].Theta()*TMath::RadToDeg(), result_min.Variables.at(t).Pull, GetTagger()->GetTaggedTime(tag));
                    six_fit_Pulls_g_phi_vs_th_CB->Fill(photons_rec[iPull].Theta()*TMath::RadToDeg(), result_min.Variables.at(u).Pull, GetTagger()->GetTaggedTime(tag));
                }
                else
                {
                    six_fit_Pulls_g_E_vs_E_TAPS->Fill(photons_rec[iPull].E(), result_min.Variables.at(s).Pull, GetTagger()->GetTaggedTime(tag));
                    six_fit_Pulls_g_th_vs_th_TAPS->Fill(photons_rec[iPull].Theta()*TMath::RadToDeg(), result_min.Variables.at(t).Pull, GetTagger()->GetTaggedTime(tag));
                    six_fit_Pulls_g_phi_vs_th_TAPS->Fill(photons_rec[iPull].Theta()*TMath::RadToDeg(), result_min.Variables.at(u).Pull, GetTagger()->GetTaggedTime(tag));
                }

           }

           six_fit_Pulls_p_th_vs_th_TAPS->Fill(GetTracks()->GetVector(set_min[1]).Theta()*TMath::RadToDeg(), result_min.Variables.at("Proton[1]").Pull, GetTagger()->GetTaggedTime(tag));
           six_fit_Pulls_p_phi_vs_th_TAPS->Fill(GetTracks()->GetVector(set_min[1]).Theta()*TMath::RadToDeg(), result_min.Variables.at("Proton[2]").Pull, GetTagger()->GetTaggedTime(tag));

           TLorentzVector proton_fit = FitParticle::Make(proton, MASS_PROTON);

           for(UInt_t igam_fit = 0; igam_fit < Photons_six.size(); igam_fit++)
           {
               photons_fit[igam_fit] = FitParticle::Make(Photons_six[igam_fit], 0);
               etap_fit += photons_fit[igam_fit];
           }
           six_fit_IM->Fill( etap_fit.M(),GetTagger()->GetTaggedTime(tag) );

           sigma_eta = 25; sigma_pi0 = 15;

           Double_t    chi2min_eta2pi  = 1.0e6;
           Double_t    chi2min_3pi     = 1.0e6;
           Double_t    probmin_eta2pi_2  = 1.0e6;
           Double_t    probmin_3pi_2     = 1.0e6;
           std::vector<int> imin_eta2pi;
           std::vector<int> imin_3pi;
           imin_eta2pi.resize(0);
           imin_3pi.resize(0);
           GetBest6gCombination(sigma_eta, sigma_pi0, chi2min_eta2pi, chi2min_3pi, imin_eta2pi, imin_3pi );

           six_fit_PDF_eta2pi_v_3pi->Fill(TMath::Prob(chi2min_eta2pi,2), TMath::Prob(chi2min_3pi,2), GetTagger()->GetTaggedTime(tag));

           test_correct_hypothesis(probmin_eta2pi_2, probmin_3pi_2, set_min, imin_eta2pi, imin_3pi);
           six_fit_PDF_eta2pi_v_3pi_2->Fill(probmin_eta2pi_2, probmin_3pi_2,GetTagger()->GetTaggedTime(tag));

           TLorentzVector g[3];
           g[0] = photons_fit[imin_eta2pi[0]] + photons_fit[imin_eta2pi[1]];
           g[1] = photons_fit[imin_eta2pi[2]] + photons_fit[imin_eta2pi[3]];
           g[2] = photons_fit[imin_eta2pi[4]] + photons_fit[imin_eta2pi[5]];

           TLorentzVector h[3];
           h[0] = photons_fit[imin_3pi[0]] + photons_fit[imin_3pi[1]];
           h[1] = photons_fit[imin_3pi[2]] + photons_fit[imin_3pi[3]];
           h[2] = photons_fit[imin_3pi[4]] + photons_fit[imin_3pi[5]];

           TLorentzVector rc[3];
           rc[0] = photons_rec[imin_3pi[0]] + photons_rec[imin_3pi[1]];
           rc[1] = photons_rec[imin_3pi[2]] + photons_rec[imin_3pi[3]];
           rc[2] = photons_rec[imin_3pi[4]] + photons_rec[imin_3pi[5]];

//           if( (TMath::Prob(chi2min_3pi,2) > 0.0) && (TMath::Prob(chi2min_eta2pi,2) < 0.20) ) // dir 3pi0
           if( (probmin_3pi_2 > 0.0) && (probmin_eta2pi_2 < 0.20) ) // dir 3pi0
           {
               six_fit_best_3pi_IM_v_E->Fill(h[0].E(), h[0].M(), GetTagger()->GetTaggedTime(tag));
               six_fit_best_3pi_IM_v_E->Fill(h[1].E(), h[1].M(), GetTagger()->GetTaggedTime(tag));
               six_fit_best_3pi_IM_v_E->Fill(h[2].E(), h[2].M(), GetTagger()->GetTaggedTime(tag));

               six_phy_3pi_IMpipi_v_IMppi->Fill((h[0]+h[1]).M2()/1.0e6, (h[2]+proton_fit).M2()/1.0e6, GetTagger()->GetTaggedTime(tag));
               six_phy_3pi_IMpipi_v_IMppi->Fill((h[2]+h[0]).M2()/1.0e6, (h[1]+proton_fit).M2()/1.0e6, GetTagger()->GetTaggedTime(tag));
               six_phy_3pi_IMpipi_v_IMppi->Fill((h[1]+h[2]).M2()/1.0e6, (h[0]+proton_fit).M2()/1.0e6, GetTagger()->GetTaggedTime(tag));


               //            six_fit_re_v_eth
               //            six_fit_mgg_v_eth
               for(Int_t isix = 0; isix < 6; isix++)
               {
                    int imass = int(isix/2);
                    if(GetTracks()->HasCB(set_min[imin_3pi[isix]+2]))
                    {
                        double En = photons_rec[imin_3pi[isix]].E();
                        double Th = photons_rec[imin_3pi[isix]].Theta()*TMath::RadToDeg();
                        double mgg = rc[imass].M();
                        int nEN =  int(photons_rec[imin_3pi[isix]].E()/10.);
                        int nTH =  100*int(photons_rec[imin_3pi[isix]].Theta()*TMath::RadToDeg()/5.0);
                        int nBIN = nEN + nTH;

                        int nDET = 100*(int(detnr[imin_3pi[isix]]/16));
                        int nBIN2 = nEN + nDET;

                        IMgg_v_det_3pi0_CB->Fill( rc[imass].M() , detnr[imin_3pi[isix]], GetTagger()->GetTaggedTime(tag));
                        IMgg_v_det_3pi0_CB_fit->Fill( h[imass].M() , detnr[imin_3pi[isix]], GetTagger()->GetTaggedTime(tag));
                        six_fit_mgg_v_eth->Fill(nBIN, rc[imass].M(),GetTagger()->GetTaggedTime(tag));
                        six_fit_mgg_v_edet->Fill(nBIN2, rc[imass].M(),GetTagger()->GetTaggedTime(tag));

                        n_fit_mgg_v_e->Fill(En,mgg, GetTagger()->GetTaggedTime(tag));
                        n_fit_mgg_v_eth->Fill(nBIN,mgg, GetTagger()->GetTaggedTime(tag));
                        n_fit_mgg_v_edet->Fill(nBIN2,mgg, GetTagger()->GetTaggedTime(tag));

                    }
                    else
                    {
                        IMgg_v_det_3pi0_TAPS->Fill( rc[imass].M() , detnr[imin_3pi[isix]], GetTagger()->GetTaggedTime(tag));
                        IMgg_v_det_3pi0_TAPS_fit->Fill( h[imass].M() , detnr[imin_3pi[isix]], GetTagger()->GetTaggedTime(tag));

                        double En = photons_rec[imin_3pi[isix]].E();
                        double mgg = rc[imass].M();
                        int nEN =  int(En/10.0);
                        int nDET = 100*(int(detnr[imin_3pi[isix]]/16));
                        int nBIN2 = nEN + nDET;

                        six_fit_mgg_v_edet_TAPS->Fill(nBIN2,mgg, GetTagger()->GetTaggedTime(tag));
                        n_fit_mgg_v_edet_TAPS->Fill(nBIN2,mgg, GetTagger()->GetTaggedTime(tag));
                        n_fit_mgg_v_e_TAPS->Fill(En, mgg, GetTagger()->GetTaggedTime(tag));

                    }
               }
                six_fit_IM_3pi->Fill(etap_fit.M(), GetTagger()->GetTaggedTime(tag));

           }
//           if( (TMath::Prob(chi2min_eta2pi,2) > 0.01) && (TMath::Prob(chi2min_3pi,2) < 0.20) ) //eta prime
           if( (probmin_3pi_2 < 0.2) && (probmin_eta2pi_2 > 0.05) ) //eta prime
           {
               six_fit_PDF_eta2pi_v_Meta2pi->Fill(TMath::Prob(chi2min_eta2pi,2),etap_fit.M(),GetTagger()->GetTaggedTime(tag));
               six_fit_PDF_2_eta2pi_v_Meta2pi->Fill(probmin_eta2pi_2,etap_fit.M(),GetTagger()->GetTaggedTime(tag));
               six_fit_best_eta_IM_v_E->Fill(g[0].E(), g[0].M(), GetTagger()->GetTaggedTime(tag));
               six_fit_best_2pi_IM_v_E->Fill(g[1].E(), g[1].M(), GetTagger()->GetTaggedTime(tag));
               six_fit_best_2pi_IM_v_E->Fill(g[2].E(), g[2].M(), GetTagger()->GetTaggedTime(tag));

               // here fill Energy vs IM(gg) vs detector nr for eta and pi0 separately
               six_phy_etapr_v_BeamE->Fill((GetTagger()->GetVector(set_min[0])).E(), etap_fit.M(), GetTagger()->GetTaggedTime(tag));

               six_fit_best_eta->Fill( (g[0]).M(), GetTagger()->GetTaggedTime(tag) );
               six_fit_best_eta_E_v_th->Fill( g[1].E(), g[1].Theta()*TMath::RadToDeg(), GetTagger()->GetTaggedTime(tag));
               six_fit_best_eta_E_v_th->Fill( g[2].E(), g[2].Theta()*TMath::RadToDeg(), GetTagger()->GetTaggedTime(tag));
               six_fit_best_2pi->Fill( g[1].M(), GetTagger()->GetTaggedTime(tag) );
               six_fit_best_2pi->Fill( g[2].M(), GetTagger()->GetTaggedTime(tag) );

               six_fit_IM_eta2pi->Fill( etap_fit.M(), GetTagger()->GetTaggedTime(tag) );


 // now do a final kinematical fit with the two eta->gg candidates set to m_eta

               beam_eta.SetFromVector( GetTagger()->GetVector(set_min[0]) );
               obs.resize(0);
               unc.resize(0);
               unc = Get_unc(3,0,obs);
               beam_eta.Smear_tmp(unc, 2);
               beam_eta.APLCONSettings();

//               beam_eta.Smear(3, true);
               beam_eta.APLCONSettings();

               proton_eta.SetFromVector( GetTracks()->GetVector(set_min[1]) );
               obs.resize(0);
               obs.push_back(proton_eta.Ek);
               obs.push_back(proton_eta.Theta);
               unc.resize(0);
               unc = Get_unc( 2, 2, obs);
               proton_eta.Smear_tmp(unc, 1);
//               proton_eta.Smear(2, false);

               Int_t n_photons_min = 0;
               for ( Int_t jgam_min = 0; jgam_min < 6 ; jgam_min++ )
               {
                   Photons_six_eta[n_photons_min].SetFromVector( GetTracks()->GetVector( set_min[imin_eta2pi[jgam_min]+2] ) );
                   if( GetTracks()->HasCB(set_min[imin_eta2pi[jgam_min]+2]) )
                   {
//                       Photons_six_eta[n_photons_min].Smear(1, true);
                       obs.resize(0);
                       obs.push_back(Photons_six_eta[n_photons_min].Ek);
                       obs.push_back(Photons_six_eta[n_photons_min].Theta);
                       unc.resize(0);
                       unc = Get_unc( 1, 1, obs);
                       Photons_six_eta[n_photons_min].Smear_tmp(unc, 0);
                   }
                   else // ( GetTracks()->HasTAPS(jgam_min) )
                   {
//                       Photons_six_eta[n_photons_min].Smear(2, true);
                       obs.resize(0);
                       obs.push_back(Photons_six_eta[n_photons_min].Ek);
                       obs.push_back(Photons_six_eta[n_photons_min].Theta);
                       unc.resize(0);
                       unc = Get_unc( 2, 1, obs);
                       Photons_six_eta[n_photons_min].Smear_tmp(unc, 0);
                   }
                   n_photons_min++;
               }

               const APLCON::Result_t& result_eta = kinfit_eta.DoFit();
               if(result_eta.Status == APLCON::Result_Status_t::Success)
               {
                   six_fit_eta_pdf->Fill( result_eta.Probability, GetTagger()->GetTaggedTime(tag) );
//                   if(result_eta.Probability > 0.01)
//                   {
                       TLorentzVector etap_fit_eta(0.0, 0.0, 0.0, 0.0);
                       for(UInt_t igam_fit = 0; igam_fit < Photons_six.size(); igam_fit++)
                       {
                           photons_fit_eta[igam_fit] = FitParticle::Make(Photons_six_eta[igam_fit], 0);
                           etap_fit_eta += photons_fit_eta[igam_fit];
                       }

                       six_fit_eta_PDF_v_Metapr->Fill(result_eta.Probability,etap_fit_eta.M(),GetTagger()->GetTaggedTime(tag));

                       TLorentzVector fin[3];
                       fin[0] = photons_fit_eta[0] + photons_fit_eta[1];
                       fin[1] = photons_fit_eta[2] + photons_fit_eta[3];
                       fin[2] = photons_fit_eta[4] + photons_fit_eta[5];

                       six_fiteta_IM2g->Fill(fin[0].M(),GetTagger()->GetTaggedTime(tag) );

                       six_phy_etapr_eta_v_BeamE->Fill((GetTagger()->GetVector(set_min[0])).E(), etap_fit_eta.M(), GetTagger()->GetTaggedTime(tag));

                       Double_t m_etapi01_fit = 0;
                       Double_t m_etapi02_fit = 0;
                       Double_t m_2pi0_fit = 0;
                       m2pi0_metapi0( fin ,  m_etapi01_fit, m_etapi02_fit, m_2pi0_fit);
                       six_phy_M_pi1pi2_v_etapr->Fill(m_2pi0_fit / 1.0e3, etap_fit_eta.M(), GetTagger()->GetTaggedTime(tag));
                       six_phy_M_etapi_v_etapr->Fill(m_etapi01_fit / 1.0e3, etap_fit_eta.M(), GetTagger()->GetTaggedTime(tag));
                       six_phy_M_etapi_v_etapr->Fill(m_etapi02_fit / 1.0e3, etap_fit_eta.M(), GetTagger()->GetTaggedTime(tag));
                       true_six_phy_dMpipi_v_Mpipi->Fill(m_2pi0_fit / 1.0e3, (m_2pi0_fit/1.0e3-m_2pi0True*1.0e3));

                       Double_t Xfit1 = -2.0;
                       Double_t Xfit2 = -2.0;
                       Double_t Yfit = -2.0;
                       Int_t DP_binnr_fit1 = -100;
                       Int_t DP_binnr_fit2 = -100;
                       DalitzPlot(fin, Xfit1, Xfit2, Yfit, DP_binnr_fit1, DP_binnr_fit2);

                       six_phy_DP->Fill(DP_binnr_fit1, etap_fit.M(), GetTagger()->GetTaggedTime(tag) );
                       six_phy_DP->Fill(DP_binnr_fit2, etap_fit.M(), GetTagger()->GetTaggedTime(tag) );
                       if( TMath::Abs(Xfit1 - Xtrue1) > TMath::Abs(Xfit1 - Xtrue2))
                       {
                           Double_t Xtemp = Xtrue1;
                           Xtrue1 = Xtrue2;
                           Xtrue2 = Xtemp;
                       }

                       true_six_phy_dX_v_DPbin->Fill(DP_binnr_fit1, Xfit1 - Xtrue1);
                       true_six_phy_dY_v_DPbin->Fill(DP_binnr_fit1, Yfit - Ytrue);
                       true_six_phy_dX_v_DPbin->Fill(DP_binnr_fit2, Xfit2 - Xtrue2);
                       true_six_phy_dY_v_DPbin->Fill(DP_binnr_fit2, Yfit - Ytrue);
//                    }
               }
           }
        }
    }
}


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


void AdlarsonPhysics::test_correct_hypothesis(Double_t& prob_eta2pi, Double_t& prob_3pi, std::vector<Int_t>& set_min, std::vector<int>& imin_eta2pi, std::vector<int>& imin_3pi)
{
    Double_t prob_3pi_fit = 1.0e6;
    Double_t prob_eta2pi_fit = 1.0e6;

    beam_3pi.SetFromVector( GetTagger()->GetVector(set_min[0]) );
//    beam_3pi.Smear(3, true);
//    beam_3pi.APLCONSettings();
    obs.resize(0);
    unc.resize(0);
    unc = Get_unc(3,0,obs);
    beam_3pi.Smear_tmp(unc, 2);
    beam_3pi.APLCONSettings();

    proton_3pi.SetFromVector( GetTracks()->GetVector(set_min[1]) );
    obs.resize(0);
    obs.push_back(proton_3pi.Ek);
    obs.push_back(proton_3pi.Theta);
    unc.resize(0);
    unc = Get_unc( 2, 2, obs);
    proton_3pi.Smear_tmp(unc, 1);
//    proton_3pi.Smear(2, false);

    beam_eta2pi.SetFromVector( GetTagger()->GetVector(set_min[0]) );
//    beam_eta2pi.Smear(3, true);
//    beam_eta2pi.APLCONSettings();
    obs.resize(0);
    unc.resize(0);
    unc = Get_unc(3,0,obs);
    beam_eta2pi.Smear_tmp(unc, 2);
    beam_eta2pi.APLCONSettings();

    proton_eta2pi.SetFromVector( GetTracks()->GetVector(set_min[1]) );
//    proton_eta2pi.Smear(2, false);
    obs.resize(0);
    obs.push_back(proton_eta2pi.Ek);
    obs.push_back(proton_eta2pi.Theta);
    unc.resize(0);
    unc = Get_unc( 2, 2, obs);
    proton_eta2pi.Smear_tmp(unc, 1);

    Int_t n_photons_min = 0;
    for ( Int_t jgam_min = 0; jgam_min < 6 ; jgam_min++ )
    {
        Photons_six_3pi[n_photons_min].SetFromVector( GetTracks()->GetVector( set_min[ imin_3pi[jgam_min] + 2] ) );
        Photons_six_eta2pi[n_photons_min].SetFromVector( GetTracks()->GetVector( set_min[ imin_eta2pi[jgam_min] + 2] ) );
        if( GetTracks()->HasCB(set_min[ imin_3pi[jgam_min] + 2] ) )
        {
//            Photons_six_3pi[n_photons_min].Smear(1, true);
            obs.resize(0);
            obs.push_back(Photons_six_3pi[n_photons_min].Ek);
            obs.push_back(Photons_six_3pi[n_photons_min].Theta);
            unc.resize(0);
            unc = Get_unc( 1, 1, obs);
            Photons_six_3pi[n_photons_min].Smear_tmp(unc ,0);
        }
        else // ( GetTracks()->HasTAPS(jgam_min) )
        {
//            Photons_six_3pi[n_photons_min].Smear(2, true);
            obs.resize(0);
            obs.push_back(Photons_six_3pi[n_photons_min].Ek);
            obs.push_back(Photons_six_3pi[n_photons_min].Theta);
            unc.resize(0);
            unc = Get_unc( 2, 1, obs);
            Photons_six_3pi[n_photons_min].Smear_tmp(unc, 0);
        }

        if( GetTracks()->HasCB( set_min[ imin_eta2pi[jgam_min] + 2 ] ) )
        {
//            Photons_six_eta2pi[n_photons_min].Smear(1, true);
            obs.resize(0);
            obs.push_back(Photons_six_eta2pi[n_photons_min].Ek);
            obs.push_back(Photons_six_eta2pi[n_photons_min].Theta);
            unc.resize(0);
            unc = Get_unc( 1, 1, obs);
            Photons_six_eta2pi[n_photons_min].Smear_tmp(unc ,0);
        }
        else // ( GetTracks()->HasTAPS(jgam_min) )
        {
//            Photons_six_eta2pi[n_photons_min].Smear(2, true);
            obs.resize(0);
            obs.push_back(Photons_six_eta2pi[n_photons_min].Ek);
            obs.push_back(Photons_six_eta2pi[n_photons_min].Theta);
            unc.resize(0);
            unc = Get_unc( 2, 1, obs);
            Photons_six_eta2pi[n_photons_min].Smear_tmp(unc ,0);
        }
        n_photons_min++;
    }

    const APLCON::Result_t& result_3pi = kinfit3pi.DoFit();
    if(result_3pi.Status == APLCON::Result_Status_t::Success)
        prob_3pi_fit = result_3pi.Probability;
    else
        prob_3pi_fit = 0.0;


    const APLCON::Result_t& result_eta2pi = kinfiteta2pi.DoFit();
    if(result_eta2pi.Status == APLCON::Result_Status_t::Success)
        prob_eta2pi_fit = result_eta2pi.Probability;
    else
        prob_eta2pi_fit = 0.0;

    prob_eta2pi = prob_eta2pi_fit;
    prob_3pi = prob_3pi_fit;

    return;
}




void AdlarsonPhysics::tengAnalysis(UInt_t ipr)
{
    for(Int_t tag = 0; tag < GetTagger()->GetNTagged(); tag++)
    {
        sigma_eta = 25; sigma_pi0 = 15;
        photons_rec.resize(0);
        IM10g_vec.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
        TLorentzVector etap10g_fit(0.0, 0.0, 0.0, 0.0);

        std::vector<Int_t> set;     // particle conf: tagger, proton, photon 1, 2, 3, 4, 5, 6, 7, 8, 9, 10;
        set.resize(0);

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
        beam10g.Smear_tmp(unc, 2);
//        beam10g.Smear(3, true);
        beam10g.APLCONSettings();



        set.push_back(ipr);
        proton10g.SetFromVector( GetTracks()->GetVector(ipr) );
//        proton10g.Smear(2, false);
        obs.resize(0);
        obs.push_back(proton10g.Ek);
        obs.push_back(proton10g.Theta);
        unc.resize(0);
        unc = Get_unc( 2, 2, obs);
        proton10g.Smear_tmp(unc, 1);

        UInt_t n_photons = 0;
        for ( UInt_t jgam = 0; jgam < ClustersInTime.size() ; jgam++ )
        {
            UInt_t kgam = ClustersInTime[jgam];
            if( kgam != ipr ) // the id proton cluster
            {
                set.push_back(kgam);
                photons_rec[n_photons] = GetTracks()->GetVector(kgam);
                IM10g_vec += photons_rec[n_photons];
                Photons_ten[n_photons].SetFromVector( GetTracks()->GetVector(kgam) );

                if( GetTracks()->HasCB(kgam) )
                {
                    obs.resize(0);
                    obs.push_back(Photons_ten[n_photons].Ek);
                    obs.push_back(Photons_ten[n_photons].Theta);
                    unc.resize(0);
                    unc = Get_unc( 1, 1, obs);
                    Photons_ten[n_photons].Smear_tmp(unc, 0);
//                    Photons_ten[n_photons].Smear(1, true);
                }
                else // ( GetTracks()->HasTAPS(kgam) )
                {
                    obs.resize(0);
                    obs.push_back(Photons_ten[n_photons].Ek);
                    obs.push_back(Photons_ten[n_photons].Theta);
                    unc.resize(0);
                    unc = Get_unc( 2, 1, obs);
                    Photons_ten[n_photons].Smear_tmp(unc, 0);
//                    Photons_ten[n_photons].Smear(2, true);
                }
                n_photons++;
            }
        }

        const APLCON::Result_t& result10g = kinfit10g.DoFit();
        if(result10g.Status == APLCON::Result_Status_t::Success)
        {
            if( result10g.ChiSquare <  chi2_min )
            {
                chi2_min = result10g.ChiSquare;
                prob_min = result10g.Probability;
            }
        }
        // Here run kinfit with the best combination:

        if( (prob_min < 0.01) || ( TMath::Abs(prob_min - 1.0e6) < 1.0e-4) ) continue;

        kfit_chi2_10g->Fill( result10g.ChiSquare, GetTagger()->GetTaggedTime(tag) );
        kfit_pdf_10g->Fill( result10g.Probability, GetTagger()->GetTaggedTime(tag) );

        Int_t inr = 0;
           for(const auto& it_map : result10g.Variables) {
//                const string& varname = it_map.first;
                const APLCON::Result_Variable_t& var = it_map.second;
                kfit_Pulls_10g->Fill(var.Pull, inr, GetTagger()->GetTaggedTime(tag));
                inr++;
           }

        ten_rec_IM->Fill( IM10g_vec.M(),GetTagger()->GetTaggedTime(tag) );
        ten_rec_IM_v_MMp->Fill( MMp , IM10g_vec.M(),GetTagger()->GetTaggedTime(tag));

        for(UInt_t igam_fit = 0; igam_fit < Photons_ten.size(); igam_fit++)
        {
            photons_fit[igam_fit] = FitParticle::Make(Photons_ten[igam_fit], 0);
            etap10g_fit += photons_fit[igam_fit];
        }

        Double_t chi2_eta2pi0 = 1.0e6;
        std::vector<int> imin_eta3pi2pi;
        imin_eta3pi2pi.resize(0);
        IM10g_fit->Fill(etap10g_fit.M(),GetTagger()->GetTaggedTime(tag));
        GetBest10gCombination(sigma_eta, sigma_pi0, chi2_eta2pi0, imin_eta3pi2pi);
    }
}


void AdlarsonPhysics::GetBest10gCombination( Double_t& sigma_eta, Double_t& sigma_pi0, Double_t& chi2min_eta3pi, std::vector<int>& imin_eta3pi2pi )
{
    TLorentzVector  eta_6g_cand;
    Double_t        im6g, im6g_min;
    Double_t        chi2;
    Double_t        chi2_min = 1.0e6;
    UInt_t          imin = 211;
    for(Int_t iperm = 0; iperm < 210; iperm++)
    {
        eta_6g_cand.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
        im6g = 0;
        chi2 = 1.0e6;
        for(Int_t jperm = 0; jperm < 6; jperm++)
        {
            eta_6g_cand += ( photons_fit[perm6outof10g[iperm][jperm]] );
        }
        im6g = eta_6g_cand.M();
        chi2 = TMath::Power( (im6g - MASS_ETA )/sigma_eta ,2);
        if( chi2 < chi2_min )
        {
            chi2_min = chi2;
            imin     = iperm;
            im6g_min = im6g;

        }
    }

    double hej = im6g_min;
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

void AdlarsonPhysics::FitParticle::Smear(Int_t apparatus_nr, Bool_t measured) {
    // set the sigmas here,
    // then the fitter knows them as well (because they're linked)

    if( apparatus_nr == 1 ) //CB
    {
        if(measured)
        {
            Ek_Sigma = (0.02*(Ek/1.0e3)*pow((Ek/1.0e3),-0.36))*1.0e3*2.0;
        }
        else
            Ek_Sigma = 0;
//        Theta_Sigma = 2.5*TMath::DegToRad()*2.0;
        Theta_Sigma = 4.5*TMath::DegToRad();
        Phi_Sigma = 2.5*TMath::DegToRad()/sin(Theta);
//        Theta_Sigma = 2.0*TMath::DegToRad();
//        Phi_Sigma = TMath::DegToRad()/sin(Theta);
    }
    else if( apparatus_nr == 2 ) //TAPS
    {
        if(measured)
        {
            Ek_Sigma = ((0.018 + 0.008*TMath::Sqrt(Ek/1.0e3))*(Ek/1.0e3))*1.0e3*2.0;
        }
        else
            Ek_Sigma = 0;

        Theta_Sigma = 2.0*TMath::DegToRad();
        Phi_Sigma = 2.0*TMath::DegToRad();
 //       Theta_Sigma = 1.0*TMath::DegToRad();
 //       Phi_Sigma = 1.0*TMath::DegToRad();
    }
    else  // beam
    {
        Ek_Sigma = 2.0;
        Theta_Sigma = 1e-3;
        Phi_Sigma = 1e-3;
    }
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
        if(apparatus_nr == 1 )  // CB
        {
//            Ek_s = (0.02*(Ek/1.0e3)*pow((Ek/1.0e3),-0.36))*1.0e3*2.0/Ek;
            Ek_s = g_CB_e->GetBinContent( g_CB_e->FindBin(Ek,theta) );
            Theta_s = g_CB_th->GetBinContent( g_CB_th->FindBin(Ek,theta) );
            Phi_s = g_CB_fi->GetBinContent( g_CB_fi->FindBin(Ek,theta) );
        }
        else                    // TAPS
        {
//            Ek_s = ((0.018 + 0.008*TMath::Sqrt(Ek/1.0e3))*(Ek/1.0e3))*1.0e3*2.0/Ek;
            Ek_s = g_TAPS_e->GetBinContent( g_TAPS_e->FindBin( Ek,theta ) );
            Theta_s = g_TAPS_th->GetBinContent(g_TAPS_th->FindBin( Ek,theta ));
            Phi_s = g_TAPS_fi->GetBinContent(g_TAPS_fi->FindBin( Ek,theta ));
        }
    }
    else if( particle == 2 )    // proton
    {
        double Ek = obs[0];
        double theta = obs[1]*TMath::RadToDeg();
        if( apparatus_nr == 1 )  // CB
        {
            std::cout << "Proton not detected in CB, check your code! "<< std::endl;
        }
        if( apparatus_nr == 2 )// TAPS
        {
        //    Ek_s = p_TAPS_e->GetBinContent(p_TAPS_e->FindBin( Ek,theta ));
            Ek_s = 0.0;
            Theta_s = p_TAPS_th->GetBinContent(p_TAPS_th->FindBin( Ek,theta ));
            Phi_s = p_TAPS_fi->GetBinContent(p_TAPS_fi->FindBin( Ek,theta ));
        }
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


void AdlarsonPhysics::FitParticle::Smear_tmp(std::vector<double> unc, int particle) {
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

void AdlarsonPhysics::Kinfit_test()
{
    TLorentzVector beam_true, proton_true;
    TLorentzVector photons_true[6];

    beam_true.SetPxPyPzE(0., 0., etapr_6gTrue.GetTrueBeamEnergy()*1000, etapr_6gTrue.GetTrueBeamEnergy()*1000);
//    beam.SetFromVector( beam_true );
//    beam.Smear(3, true);

    proton_true.SetPxPyPzE(etapr_6gTrue.GetTrueProtonLV().Px()*1000, etapr_6gTrue.GetTrueProtonLV().Py()*1000, etapr_6gTrue.GetTrueProtonLV().Pz()*1000,etapr_6gTrue.GetTrueProtonLV().E()*1000);
    proton.SetFromVector( proton_true );
    proton.Smear(2, true);

    for(UInt_t i = 0; i < etapr_6gTrue.GetNgamma(); i++ )
    {
        (photons_true[i]).SetPxPyPzE( etapr_6gTrue.GetTrueGammaLV(i).Px()*1000, etapr_6gTrue.GetTrueGammaLV(i).Py()*1000,etapr_6gTrue.GetTrueGammaLV(i).Pz()*1000, etapr_6gTrue.GetTrueGammaLV(i).E()*1000 );
        Photons_six[i].SetFromVector( photons_true[i] );
        Photons_six[i].Smear(1, true);
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
    Double_t DeltaE;
    Double_t Erec, Ec1, Ec2, Efin;
    Double_t g1;
    std::vector<double> pol;
    int mod;
    for (int i = 0; i < GetTracks()->GetNTracks() ; i++)
    {
        if( GetTracks()->HasCB(i) )
        {

            //            DeltaE = Ec_temp*(Double_t)EvdetCB->GetBinContent(EvdetCB->FindBin(Ec_temp, GetTracks()->GetTheta(i)));
            //             ( GetTracks()->GetClusterEnergy(i) - E_true )/(GetTracks()->GetClusterEnergy(i));
            //            Ec = Ec_temp; //(Double_t)EvdetCB->GetBinContent(EvdetCB->FindBin(Ec_temp, GetTracks()->GetCentralCrystal(i)));

            Erec = GetTracks()->GetVector(i).E();
            Ec1 = CBgain[GetTracks()->GetCentralCrystal(i)]*Erec; // linear gain

//            mod = (GetTracks()->GetCentralCrystal(i)/16);             // energy dependence as fcn of module
//            pol = CB_Ecorr.at(mod);                                   // mod is the key to vector
////
//            if(pol[0] > 0){
//                if(Ec1 > pol[0])  // to have continuity as fcn of energy
//                    Ec2 = pol[0];
//                else
//                    Ec2 = Ec1;

//                g1 = 0;
//                for(int j = 1; j < 8; j++)
//                {
//                    g1 += pol[j]*TMath::Power(Ec2,j-1);
//                }
//                Efin = Ec1*(g1*g1);
//            }
//            else
//                Efin = Ec1;

            DeltaE = Ec1*(Double_t)EvdetCB->GetBinContent(EvdetCB->FindBin(Ec1, GetTracks()->GetTheta(i)));
            Efin = Ec1 - DeltaE;
            //             ( GetTracks()->GetClusterEnergy(i) - E_true )/(GetTracks()->GetClusterEnergy(i));
            //            Ec = Ec_temp; //(Double_t)EvdetCB->GetBinContent(EvdetCB->FindBin(Ec_temp, GetTracks()->GetCentralCrystal(i)));

            tracks->SetClusterEnergy(i, Efin);
        }
        else if(GetTracks()->HasTAPS(i) )
        {
            Erec = GetTracks()->GetVector(i).E();
            Ec1 = TAPSgain[GetTracks()->GetCentralCrystal(i)]*Erec;
            tracks->SetClusterEnergy(i, Ec1);
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

