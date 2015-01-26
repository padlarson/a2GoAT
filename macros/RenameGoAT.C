void RenameGoAT(TString sFile){
  
  TFile fOld(sFile,"READ");

  // Acqu Trees
  TTree *taggerOld       = (TTree*)fOld.Get("treeTagger");
  TTree *tracksOld       = (TTree*)fOld.Get("treeRawEvent");
  TTree *detectorHitsOld = (TTree*)fOld.Get("treeDetectorHits");
  TTree *linPolOld       = (TTree*)fOld.Get("treeLinPol");
  TTree *triggerOld      = (TTree*)fOld.Get("treeTrigger");
  TTree *scalersOld       = (TTree*)fOld.Get("treeScaler");

  // GoAT Trees
  TTree *eventParametersOld = (TTree*)fOld.Get("treeEventParameters");
  TTree *rootinosOld        = (TTree*)fOld.Get("rootino");
  TTree *photonsOld         = (TTree*)fOld.Get("gamma");
  TTree *electronsOld       = (TTree*)fOld.Get("e-");
  TTree *chargedPionsOld    = (TTree*)fOld.Get("pi+");
  TTree *protonsOld         = (TTree*)fOld.Get("proton");
  TTree *neutronsOld        = (TTree*)fOld.Get("neutron");
  TTree *neutralPionsOld    = (TTree*)fOld.Get("pi0");
  TTree *etasOld            = (TTree*)fOld.Get("eta");
  TTree *etaPrimesOld       = (TTree*)fOld.Get("eta'");

  TClonesArray *particles   = new TClonesArray("TLorentzVector", 64);

  sFile.ReplaceAll(".root","_renamed.root");
  TFile fNew(sFile,"RECREATE");

  if(taggerOld){
    printf("Renaming treeTagger to tagger\n");

    Int_t nTagged          = 0;
    Double_t *taggedEnergy = new Double_t[4096];
    Int_t *taggedChannel   = new Int_t[4096];
    Double_t *taggedTime   = new Double_t[4096];

    TTree *tagger = new TTree("tagger", "tagger");

    if(taggerOld->GetBranch("nTagged")){
      taggerOld->SetBranchAddress("nTagged", &nTagged);
      tagger->Branch("nTagged", &nTagged,"nTagged/I");
      if(taggerOld->GetBranch("tagged_ch")){
	taggerOld->SetBranchAddress("tagged_ch", taggedChannel);
	tagger->Branch("taggedChannel", taggedChannel, "taggedChannel[nTagged]/I");
      }
      if(taggerOld->GetBranch("tagged_t")){
	taggerOld->SetBranchAddress("tagged_t", taggedTime);
	tagger->Branch("taggedTime", taggedTime, "taggedTime[nTagged]/D");
      }
      if(taggerOld->GetBranch("photonbeam_E")){
	taggerOld->SetBranchAddress("photonbeam_E", taggedEnergy);
	tagger->Branch("taggedEnergy", taggedEnergy, "taggedEnergy[nTagged]/D");
      }
    }

    for(Int_t i=0; i<taggerOld->GetEntries(); i++){
      taggerOld->GetEvent(i);
      tagger->Fill();
    }
    tagger->Write();
  }

  Int_t nParticles      = 0;
  Double_t *time        = new Double_t[128];
  Int_t *clusterSize    = new Int_t[128];
  Int_t *apparatus      = new Int_t[128];
  Double_t *vetoEnergy  = new Double_t[128];
  Double_t *MWPC0Energy = new Double_t[128];
  Double_t *MWPC1Energy = new Double_t[128];

  if(tracksOld){
    printf("Renaming treeRawEvent to tracks\n");

    Int_t nTracks           = 0;
    Double_t *clusterEnergy = new Double_t[128];
    Double_t *theta         = new Double_t[128];
    Double_t *phi           = new Double_t[128];
    Int_t *centralCrystal   = new Int_t[128];
    Int_t *centralVeto      = new Int_t[128];

    TTree *tracks = new TTree("tracks", "tracks");

    if(tracksOld->GetBranch("nParticles")){
      tracksOld->SetBranchAddress("nParticles", &nTracks);
      tracks->Branch("nTracks", &nTracks, "nTracks/I");
      if(tracksOld->GetBranch("Ek")){
	tracksOld->SetBranchAddress("Ek", clusterEnergy);
	tracks->Branch("clusterEnergy", clusterEnergy, "clusterEnergy[nTracks]/D");
      }
      if(tracksOld->GetBranch("Theta")){
	tracksOld->SetBranchAddress("Theta", theta);
	tracks->Branch("theta", theta, "theta[nTracks]/D");
      }
      if(tracksOld->GetBranch("Phi")){
	tracksOld->SetBranchAddress("Phi", phi);
	tracks->Branch("phi", phi, "phi[nTracks]/D");
      }
      if(tracksOld->GetBranch("time")){
	tracksOld->SetBranchAddress("time", time);
	tracks->Branch("time", time, "time[nTracks]/D");
      }
      if(tracksOld->GetBranch("clusterSize")){
	tracksOld->SetBranchAddress("clusterSize", clusterSize);
	tracks->Branch("clusterSize", clusterSize, "clusterSize[nTracks]/I");
      }
      if(tracksOld->GetBranch("centralCrys")){
	tracksOld->SetBranchAddress("centralCrys", centralCrystal);
	tracks->Branch("centralCrystal", centralCrystal, "centralCrystal[nTracks]/I");
      }
      if(tracksOld->GetBranch("centralVeto")){
	tracksOld->SetBranchAddress("centralVeto", centralVeto);
	tracks->Branch("centralVeto", centralVeto, "centralVeto[nTracks]/I");
      }
      if(tracksOld->GetBranch("Apparatus")){
	tracksOld->SetBranchAddress("Apparatus", apparatus);
	tracks->Branch("apparatus", apparatus, "apparatus[nTracks]/I");
      }
      if(tracksOld->GetBranch("d_E")){
	tracksOld->SetBranchAddress("d_E", vetoEnergy);
	tracks->Branch("vetoEnergy", vetoEnergy, "vetoEnergy[nTracks]/D");
      }
      if(tracksOld->GetBranch("WC0_E")){
	tracksOld->SetBranchAddress("WC0_E", MWPC0Energy);
	tracks->Branch("MWPC0Energy", MWPC0Energy, "MWPC0Energy[nTracks]/D");
      }
      if(tracksOld->GetBranch("WC1_E")){
	tracksOld->SetBranchAddress("WC1_E", MWPC1Energy);
	tracks->Branch("MWPC1Energy", MWPC1Energy, "MWPC1Energy[nTracks]/D");
      }
    }

    for(Int_t i=0; i<tracksOld->GetEntries(); i++){
      tracksOld->GetEvent(i);
      tracks->Fill();
    }
    tracks->Write();
  }

  if(detectorHitsOld){
    printf("Renaming treeDetectorHits to detectorHits\n");

    Int_t nNaIHits          = 0;
    Int_t *NaIHits          = new Int_t[860];
    Int_t *NaICluster       = new Int_t[860];
    Int_t nPIDHits          = 0;
    Int_t *PIDHits          = new Int_t[860];
    Int_t nMWPCHits         = 0;
    Int_t *MWPCHits         = new Int_t[860];
    Int_t nBaF2Hits    = 0;
    Int_t *BaF2Hits    = new Int_t[860];
    Int_t *BaF2Cluster = new Int_t[860];
    Int_t nVetoHits         = 0;
    Int_t *VetoHits         = new Int_t[860];
	
    TTree *detectorHits = new TTree("detectorHits", "detectorHits");

    if(detectorHitsOld->GetBranch("nNaI_Hits")){
      detectorHitsOld->SetBranchAddress("nNaI_Hits", &nNaIHits);
      detectorHits->Branch("nNaIHits", &nNaIHits, "nNaIHits/I");
      if(detectorHitsOld->GetBranch("NaI_Hits")){
	detectorHitsOld->SetBranchAddress("NaI_Hits", NaIHits);
	detectorHits->Branch("NaIHits", NaIHits, "NaIHits[nNaIHits]/I");
      }
      if(detectorHitsOld->GetBranch("NaI_Cluster")){
	detectorHitsOld->SetBranchAddress("NaI_Cluster", NaICluster);
	detectorHits->Branch("NaICluster", NaICluster, "NaICluster[nNaIHits]/I");
      }
    }
    if(detectorHitsOld->GetBranch("nPID_Hits")){
      detectorHitsOld->SetBranchAddress("nPID_Hits", &nPIDHits);
      detectorHits->Branch("nPIDHits", &nPIDHits, "nPIDHits/I");
      if(detectorHitsOld->GetBranch("PID_Hits")){
	detectorHitsOld->SetBranchAddress("PID_Hits", PIDHits);
	detectorHits->Branch("PIDHits", PIDHits, "PIDHits[nPIDHits]/I");
      }
    }
    if(detectorHitsOld->GetBranch("nWC_Hits")){
      detectorHitsOld->SetBranchAddress("nWC_Hits", &nMWPCHits);
      detectorHits->Branch("nMWPCHits", &nMWPCHits, "nMWPCHits/I");
      if(detectorHitsOld->GetBranch("WC_Hits")){
	detectorHitsOld->SetBranchAddress("WC_Hits", MWPCHits);
	detectorHits->Branch("MWPCHits", MWPCHits, "MWPCHits[nMWPCHits]/I");
      }
    }
    if(detectorHitsOld->GetBranch("nBaF2_PbWO4_Hits")){
      detectorHitsOld->SetBranchAddress("nBaF2_PbWO4_Hits", &nBaF2Hits);
      detectorHits->Branch("nBaF2Hits", &nBaF2Hits, "nBaF2Hits/I");
      if(detectorHitsOld->GetBranch("BaF2_PbWO4_Hits")){
	detectorHitsOld->SetBranchAddress("BaF2_PbWO4_Hits", BaF2Hits);
	detectorHits->Branch("BaF2Hits", BaF2Hits, "BaF2Hits[nBaF2Hits]/I");
      }
      if(detectorHitsOld->GetBranch("BaF2_PbWO4_Cluster")){
	detectorHitsOld->SetBranchAddress("BaF2_PbWO4_Cluster", BaF2Cluster);
	detectorHits->Branch("BaF2Cluster", BaF2Cluster, "BaF2Cluster[nBaF2Hits]/I");
      }
    }
    if(detectorHitsOld->GetBranch("nVeto_Hits")){
      detectorHitsOld->SetBranchAddress("nVeto_Hits", &nVetoHits);
      detectorHits->Branch("nVetoHits", &nVetoHits, "nVetoHits/I");
      if(detectorHitsOld->GetBranch("Veto_Hits")){
	detectorHitsOld->SetBranchAddress("Veto_Hits", VetoHits);
	detectorHits->Branch("VetoHits", VetoHits, "VetoHits[nVetoHits]/I");
      }
    }

    for(Int_t i=0; i<detectorHitsOld->GetEntries(); i++){
      detectorHitsOld->GetEvent(i);
      detectorHits->Fill();
    }
    detectorHits->Write();
  }

  if(linPolOld){
    printf("Renaming treeLinPol to linPol\n");

    Int_t plane                 = 0;
    Double_t edge               = 0;
    Double_t edgeSetting        = 0;
    Double_t *polarizationTable = new Double_t[352];
    Double_t *enhancementTable  = new Double_t[352];

    TTree *linPol = new TTree("linPol", "linPol");

    if(linPolOld->GetBranch("plane")){
      linPolOld->SetBranchAddress("plane", &plane);
      linPol->Branch("plane", &plane, "plane/I");
    }
    if(linPolOld->GetBranch("edge")){
      linPolOld->SetBranchAddress("edge", &edge);
      linPol->Branch("edge", &edge, "edge/D");
    }
    if(linPolOld->GetBranch("edgeSetting")){
      linPolOld->SetBranchAddress("edgeSetting", &edgeSetting);
      linPol->Branch("edgeSetting", &edgeSetting, "edgeSetting/D");
    }
    if(linPolOld->GetBranch("polTable")){
      linPolOld->SetBranchAddress("polTable", polarizationTable);
      linPol->Branch("polarizationTable", polarizationTable, "polarizationTable[352]/D");
    }
    if(linPolOld->GetBranch("enhTable")){
      linPolOld->SetBranchAddress("enhTable", enhancementTable);
      linPol->Branch("enhancementTable", enhancementTable, "enhancementTable[352]/D");
    }

    for(Int_t i=0; i<linPolOld->GetEntries(); i++){
      linPolOld->GetEvent(i);
      linPol->Fill();
    }
    linPol->Write();
  }

  if(triggerOld){
    printf("Renaming treeTrigger to trigger\n");

    Double_t energySum      = 0;
    Int_t multiplicity      = 0;
    Bool_t helicity         = false;
    Int_t nErrors           = 0;
    Int_t nErrors           = 0;
    Int_t *errorModuleID    = new Int_t[256];
    Int_t *errorModuleIndex = new Int_t[256];
    Int_t *errorCode        = new Int_t[256];
    Int_t nTriggerPattern   = 0;
    Int_t *triggerPattern   = new Int_t[32];

    TTree *trigger = new TTree("trigger", "trigger");

    if(triggerOld->GetBranch("ESum")){
      triggerOld->SetBranchAddress("ESum", &energySum);
      trigger->Branch("energySum", &energySum, "energySum/D");
    }
    if(triggerOld->GetBranch("Mult")){
      triggerOld->SetBranchAddress("Mult", &multiplicity);
      trigger->Branch("multiplicity", &multiplicity, "multiplicity/I");
    }
    if(triggerOld->GetBranch("Helicity")){
      triggerOld->SetBranchAddress("Helicity", &helicity);
      trigger->Branch("helicity", &helicity, "helicity/O");
    }
    if(triggerOld->GetBranch("nError")){
      triggerOld->SetBranchAddress("nError", &nErrors);
      trigger->Branch("nErrors", &nErrors, "nErrors/I");
      if(triggerOld->GetBranch("ErrModID")){
	triggerOld->SetBranchAddress("ErrModID", errorModuleID);
	trigger->Branch("errorModuleID", errorModuleID, "errorModuleID[nErrors]/I");
      }
      if(triggerOld->GetBranch("ErrModIndex")){
	triggerOld->SetBranchAddress("ErrModIndex", errorModuleIndex);
	trigger->Branch("errorModuleIndex", errorModuleIndex, "errorModuleIndex[nErrors]/I");
      }
      if(triggerOld->GetBranch("ErrCode")){
	triggerOld->SetBranchAddress("ErrCode", errorCode);
	trigger->Branch("errorCode", errorCode, "errorCode[nErrors]/I");
      }
    }
    if(triggerOld->GetBranch("nTriggerPattern")){
      triggerOld->SetBranchAddress("nTriggerPattern", &nTriggerPattern);
      trigger->Branch("nTriggerPattern", &nTriggerPattern, "nTriggerPattern/I");
      if(triggerOld->GetBranch("TriggerPattern")){
	triggerOld->SetBranchAddress("TriggerPattern", triggerPattern);
	trigger->Branch("triggerPattern", triggerPattern, "triggerPattern[nTriggerPattern]/I");
      }
    }

    for(Int_t i=0; i<triggerOld->GetEntries(); i++){
      triggerOld->GetEvent(i);
      trigger->Fill();
    }
    trigger->Write();
  }

  Int_t eventNumber = 0;

  if(scalersOld){
    printf("Renaming treeScaler to scalers\n");

    TString sScaler = scalersOld->GetBranch("Scaler")->GetTitle();
    sScaler.Remove(0,sScaler.First("[")+1);
    sScaler.Remove(sScaler.First("]"));
    const Int_t nScaler = sScaler.Atoi();
    Char_t str[256];
    sprintf(str, "scalers[%d]/i", nScaler);

    Int_t eventID = 0;
    UInt_t *scalerArray = new UInt_t[nScaler];

    TTree *scalers = new TTree("scalers", "scalers");
    
    if(scalersOld->GetBranch("eventNumber")){
      scalersOld->SetBranchAddress("eventNumber", &eventNumber);
      scalers->Branch("eventNumber", &eventNumber, "eventNumber/I");
    }
    if(scalersOld->GetBranch("eventID")){
      scalersOld->SetBranchAddress("eventID", &eventID);
      scalers->Branch("eventID", &eventID, "eventID/I");
    }
    if(scalersOld->GetBranch("Scaler")){
      scalersOld->SetBranchAddress("Scaler", scalerArray);
      scalers->Branch("scalers", scalerArray, str);
    }

    for(Int_t i=0; i<scalersOld->GetEntries(); i++){
      scalersOld->GetEvent(i);
      scalers->Fill();
    }
    scalers->Write();
  }

  if(eventParametersOld){
    printf("Renaming treeEventParameters to eventParameters\n");

    Int_t nReconstructed = 0;

    TTree *eventParameters = new TTree("eventParameters", "eventParameters");

    if(eventParametersOld->GetBranch("EventNumber")){
      eventParametersOld->SetBranchAddress("EventNumber", &eventNumber);
      eventParameters->Branch("eventNumber", &eventNumber, "eventNumber/I");
    }
    if(eventParametersOld->GetBranch("nReconstructed")){
      eventParametersOld->SetBranchAddress("nReconstructed", &nReconstructed);
      eventParameters->Branch("nReconstructed", &nReconstructed, "nReconstructed/I");
    }

    for(Int_t i=0; i<eventParametersOld->GetEntries(); i++){
      eventParametersOld->GetEvent(i);
      eventParameters->Fill();
    }
    eventParameters->Write();
  }

  if(rootinosOld){
    printf("Renaming rootino to rootinos\n");

    TTree *rootinos = new TTree("rootinos", "rootinos");

    if(rootinosOld->GetBranch("nParticles")){
      rootinosOld->SetBranchAddress("nParticles", &nParticles);
      rootinos->Branch("nParticles", &nParticles, "nParticles/i");
      if(rootinosOld->GetBranch("particles.")){
	rootinosOld->SetBranchAddress("particles.", &particles);
	rootinos->Branch("particles", &particles, 32000, 0);
      }
      if(rootinosOld->GetBranch("time")){
	rootinosOld->SetBranchAddress("time", time);
	rootinos->Branch("time", time, "time[nParticles]/D");
      }
      if(rootinosOld->GetBranch("clusterSize")){
	rootinosOld->SetBranchAddress("clusterSize", clusterSize);
	rootinos->Branch("clusterSize", clusterSize, "clusterSize[nParticles]/b");
      }
      if(rootinosOld->GetBranch("Apparatus")){
	rootinosOld->SetBranchAddress("Apparatus", apparatus);
	rootinos->Branch("apparatus", apparatus, "apparatus[nParticles]/b");
      }
      if(rootinosOld->GetBranch("d_E")){
	rootinosOld->SetBranchAddress("d_E", vetoEnergy);
	rootinos->Branch("vetoEnergy", vetoEnergy, "vetoEnergy[nParticles]/D");
      }
      if(rootinosOld->GetBranch("WC0_E")){
	rootinosOld->SetBranchAddress("WC0_E", MWPC0Energy);
	rootinos->Branch("MWPC0Energy", MWPC0Energy, "MWPC0Energy[nParticles]/D");
      }
      if(rootinosOld->GetBranch("WC1_E")){
	rootinosOld->SetBranchAddress("WC1_E", MWPC1Energy);
	rootinos->Branch("MWPC1Energy", MWPC1Energy, "MWPC1Energy[nParticles]/D");
      }
    }

    for(Int_t i=0; i<rootinosOld->GetEntries(); i++){
      rootinosOld->GetEvent(i);
      rootinos->Fill();
    }
    rootinos->Write();
  }

  if(photonsOld){
    printf("Renaming gamma to photons\n");

    TTree *photons = new TTree("photons", "photons");

    if(photonsOld->GetBranch("nParticles")){
      photonsOld->SetBranchAddress("nParticles", &nParticles);
      photons->Branch("nParticles", &nParticles, "nParticles/I");
      if(photonsOld->GetBranch("particles.")){
	photonsOld->SetBranchAddress("particles.", &particles);
	photons->Branch("particles", &particles, 32000, 0);
      }
      if(photonsOld->GetBranch("time")){
	photonsOld->SetBranchAddress("time", time);
	photons->Branch("time", time, "time[nParticles]/D");
      }
      if(photonsOld->GetBranch("clusterSize")){
	photonsOld->SetBranchAddress("clusterSize", clusterSize);
	photons->Branch("clusterSize", clusterSize, "clusterSize[nParticles]/I");
      }
      if(photonsOld->GetBranch("Apparatus")){
	photonsOld->SetBranchAddress("Apparatus", apparatus);
	photons->Branch("apparatus", apparatus, "apparatus[nParticles]/I");
      }
      if(photonsOld->GetBranch("d_E")){
	photonsOld->SetBranchAddress("d_E", vetoEnergy);
	photons->Branch("vetoEnergy", vetoEnergy, "vetoEnergy[nParticles]/D");
      }
      if(photonsOld->GetBranch("WC0_E")){
	photonsOld->SetBranchAddress("WC0_E", MWPC0Energy);
	photons->Branch("MWPC0Energy", MWPC0Energy, "MWPC0Energy[nParticles]/D");
      }
      if(photonsOld->GetBranch("WC1_E")){
	photonsOld->SetBranchAddress("WC1_E", MWPC1Energy);
	photons->Branch("MWPC1Energy", MWPC1Energy, "MWPC1Energy[nParticles]/D");
      }
    }

    for(Int_t i=0; i<photonsOld->GetEntries(); i++){
      photonsOld->GetEvent(i);
      photons->Fill();
    }
    photons->Write();
  }

  if(electronsOld){
    printf("Renaming e- to electrons\n");

    TTree *electrons = new TTree("electrons", "electrons");

    if(electronsOld->GetBranch("nParticles")){
      electronsOld->SetBranchAddress("nParticles", &nParticles);
      electrons->Branch("nParticles", &nParticles, "nParticles/I");
      if(electronsOld->GetBranch("particles.")){
	electronsOld->SetBranchAddress("particles.", &particles);
	electrons->Branch("particles", &particles, 32000, 0);
      }
      if(electronsOld->GetBranch("time")){
	electronsOld->SetBranchAddress("time", time);
	electrons->Branch("time", time, "time[nParticles]/D");
      }
      if(electronsOld->GetBranch("clusterSize")){
	electronsOld->SetBranchAddress("clusterSize", clusterSize);
	electrons->Branch("clusterSize", clusterSize, "clusterSize[nParticles]/I");
      }
      if(electronsOld->GetBranch("Apparatus")){
	electronsOld->SetBranchAddress("Apparatus", apparatus);
	electrons->Branch("apparatus", apparatus, "apparatus[nParticles]/I");
      }
      if(electronsOld->GetBranch("d_E")){
	electronsOld->SetBranchAddress("d_E", vetoEnergy);
	electrons->Branch("vetoEnergy", vetoEnergy, "vetoEnergy[nParticles]/D");
      }
      if(electronsOld->GetBranch("WC0_E")){
	electronsOld->SetBranchAddress("WC0_E", MWPC0Energy);
	electrons->Branch("MWPC0Energy", MWPC0Energy, "MWPC0Energy[nParticles]/D");
      }
      if(electronsOld->GetBranch("WC1_E")){
	electronsOld->SetBranchAddress("WC1_E", MWPC1Energy);
	electrons->Branch("MWPC1Energy", MWPC1Energy, "MWPC1Energy[nParticles]/D");
      }
    }

    for(Int_t i=0; i<electronsOld->GetEntries(); i++){
      electronsOld->GetEvent(i);
      electrons->Fill();
    }
    electrons->Write();
  }

  if(chargedPionsOld){
    printf("Renaming pi+ to chargedPions\n");

    TTree *chargedPions = new TTree("chargedPions", "chargedPions");

    if(chargedPionsOld->GetBranch("nParticles")){
      chargedPionsOld->SetBranchAddress("nParticles", &nParticles);
      chargedPions->Branch("nParticles", &nParticles, "nParticles/I");
      if(chargedPionsOld->GetBranch("particles.")){
	chargedPionsOld->SetBranchAddress("particles.", &particles);
	chargedPions->Branch("particles", &particles, 32000, 0);
      }
      if(chargedPionsOld->GetBranch("time")){
	chargedPionsOld->SetBranchAddress("time", time);
	chargedPions->Branch("time", time, "time[nParticles]/D");
      }
      if(chargedPionsOld->GetBranch("clusterSize")){
	chargedPionsOld->SetBranchAddress("clusterSize", clusterSize);
	chargedPions->Branch("clusterSize", clusterSize, "clusterSize[nParticles]/I");
      }
      if(chargedPionsOld->GetBranch("Apparatus")){
	chargedPionsOld->SetBranchAddress("Apparatus", apparatus);
	chargedPions->Branch("apparatus", apparatus, "apparatus[nParticles]/I");
      }
      if(chargedPionsOld->GetBranch("d_E")){
	chargedPionsOld->SetBranchAddress("d_E", vetoEnergy);
	chargedPions->Branch("vetoEnergy", vetoEnergy, "vetoEnergy[nParticles]/D");
      }
      if(chargedPionsOld->GetBranch("WC0_E")){
	chargedPionsOld->SetBranchAddress("WC0_E", MWPC0Energy);
	chargedPions->Branch("MWPC0Energy", MWPC0Energy, "MWPC0Energy[nParticles]/D");
      }
      if(chargedPionsOld->GetBranch("WC1_E")){
	chargedPionsOld->SetBranchAddress("WC1_E", MWPC1Energy);
	chargedPions->Branch("MWPC1Energy", MWPC1Energy, "MWPC1Energy[nParticles]/D");
      }
    }

    for(Int_t i=0; i<chargedPionsOld->GetEntries(); i++){
      chargedPionsOld->GetEvent(i);
      chargedPions->Fill();
    }
    chargedPions->Write();
  }

  if(protonsOld){
    printf("Renaming proton to protons\n");

    TTree *protons = new TTree("protons", "protons");

    if(protonsOld->GetBranch("nParticles")){
      protonsOld->SetBranchAddress("nParticles", &nParticles);
      protons->Branch("nParticles", &nParticles, "nParticles/I");
      if(protonsOld->GetBranch("particles.")){
	protonsOld->SetBranchAddress("particles.", &particles);
	protons->Branch("particles", &particles, 32000, 0);
      }
      if(protonsOld->GetBranch("time")){
	protonsOld->SetBranchAddress("time", time);
	protons->Branch("time", time, "time[nParticles]/D");
      }
      if(protonsOld->GetBranch("clusterSize")){
	protonsOld->SetBranchAddress("clusterSize", clusterSize);
	protons->Branch("clusterSize", clusterSize, "clusterSize[nParticles]/I");
      }
      if(protonsOld->GetBranch("Apparatus")){
	protonsOld->SetBranchAddress("Apparatus", apparatus);
	protons->Branch("apparatus", apparatus, "apparatus[nParticles]/I");
      }
      if(protonsOld->GetBranch("d_E")){
	protonsOld->SetBranchAddress("d_E", vetoEnergy);
	protons->Branch("vetoEnergy", vetoEnergy, "vetoEnergy[nParticles]/D");
      }
      if(protonsOld->GetBranch("WC0_E")){
	protonsOld->SetBranchAddress("WC0_E", MWPC0Energy);
	protons->Branch("MWPC0Energy", MWPC0Energy, "MWPC0Energy[nParticles]/D");
      }
      if(protonsOld->GetBranch("WC1_E")){
	protonsOld->SetBranchAddress("WC1_E", MWPC1Energy);
	protons->Branch("MWPC1Energy", MWPC1Energy, "MWPC1Energy[nParticles]/D");
      }
    }

    for(Int_t i=0; i<protonsOld->GetEntries(); i++){
      protonsOld->GetEvent(i);
      protons->Fill();
    }
    protons->Write();
  }

  if(neutronsOld){
    printf("Renaming neutron to neutrons\n");

    TTree *neutrons = new TTree("neutrons", "neutrons");

    if(neutronsOld->GetBranch("nParticles")){
      neutronsOld->SetBranchAddress("nParticles", &nParticles);
      neutrons->Branch("nParticles", &nParticles, "nParticles/I");
      if(neutronsOld->GetBranch("particles.")){
	neutronsOld->SetBranchAddress("particles.", &particles);
	neutrons->Branch("particles", &particles, 32000, 0);
      }
      if(neutronsOld->GetBranch("time")){
	neutronsOld->SetBranchAddress("time", time);
	neutrons->Branch("time", time, "time[nParticles]/D");
      }
      if(neutronsOld->GetBranch("clusterSize")){
	neutronsOld->SetBranchAddress("clusterSize", clusterSize);
	neutrons->Branch("clusterSize", clusterSize, "clusterSize[nParticles]/I");
      }
      if(neutronsOld->GetBranch("Apparatus")){
	neutronsOld->SetBranchAddress("Apparatus", apparatus);
	neutrons->Branch("apparatus", apparatus, "apparatus[nParticles]/I");
      }
      if(neutronsOld->GetBranch("d_E")){
	neutronsOld->SetBranchAddress("d_E", vetoEnergy);
	neutrons->Branch("vetoEnergy", vetoEnergy, "vetoEnergy[nParticles]/D");
      }
      if(neutronsOld->GetBranch("WC0_E")){
	neutronsOld->SetBranchAddress("WC0_E", MWPC0Energy);
	neutrons->Branch("MWPC0Energy", MWPC0Energy, "MWPC0Energy[nParticles]/D");
      }
      if(neutronsOld->GetBranch("WC1_E")){
	neutronsOld->SetBranchAddress("WC1_E", MWPC1Energy);
	neutrons->Branch("MWPC1Energy", MWPC1Energy, "MWPC1Energy[nParticles]/D");
      }
    }

    for(Int_t i=0; i<neutronsOld->GetEntries(); i++){
      neutronsOld->GetEvent(i);
      neutrons->Fill();
    }
    neutrons->Write();
  }

  UChar_t *nSubParticles    = new UChar_t[64];
  UChar_t *nSubRootinos     = new UChar_t[64];
  UChar_t *nSubPhotons      = new UChar_t[64];
  UChar_t *nSubChargedPions = new UChar_t[64];
    
  TClonesArray *subRootinos     = new TClonesArray("TLorentzVector", 64);
  TClonesArray *subPhotons      = new TClonesArray("TLorentzVector", 64);
  TClonesArray *subChargedPions = new TClonesArray("TLorentzVector", 64);

  if(neutralPionsOld){
    printf("Renaming pi0 to neutralPions\n");

    TTree *neutralPions = new TTree("neutralPions", "neutralPions");

    if(neutralPionsOld->GetBranch("nParticles")){
      neutralPionsOld->SetBranchAddress("nParticles", &nParticles);
      neutralPions->Branch("nParticles", &nParticles, "nParticles/I");
      if(neutralPionsOld->GetBranch("particles.")){
	neutralPionsOld->SetBranchAddress("particles.", &particles);
	neutralPions->Branch("particles", &particles, 32000, 0);
      }
      if(neutralPionsOld->GetBranch("time")){
	neutralPionsOld->SetBranchAddress("time", time);
	neutralPions->Branch("time", time, "time[nParticles]/D");
      }
      if(neutralPionsOld->GetBranch("clusterSize")){
	neutralPionsOld->SetBranchAddress("clusterSize", clusterSize);
	neutralPions->Branch("clusterSize", clusterSize, "clusterSize[nParticles]/I");
      }
      if(neutralPionsOld->GetBranch("Apparatus")){
	neutralPionsOld->SetBranchAddress("Apparatus", apparatus);
	neutralPions->Branch("apparatus", apparatus, "apparatus[nParticles]/I");
      }
      if(neutralPionsOld->GetBranch("d_E")){
	neutralPionsOld->SetBranchAddress("d_E", vetoEnergy);
	neutralPions->Branch("vetoEnergy", vetoEnergy, "vetoEnergy[nParticles]/D");
      }
      if(neutralPionsOld->GetBranch("WC0_E")){
	neutralPionsOld->SetBranchAddress("WC0_E", MWPC0Energy);
	neutralPions->Branch("MWPC0Energy", MWPC0Energy, "MWPC0Energy[nParticles]/D");
      }
      if(neutralPionsOld->GetBranch("WC1_E")){
	neutralPionsOld->SetBranchAddress("WC1_E", MWPC1Energy);
	neutralPions->Branch("MWPC1Energy", MWPC1Energy, "MWPC1Energy[nParticles]/D");
      }
      if(neutralPionsOld->GetBranch("nSubParticles")){
	neutralPionsOld->SetBranchAddress("nSubParticles", nSubParticles);
	neutralPions->Branch("nSubParticles", nSubParticles, "nSubParticles[nParticles]/I");
      }
      if(neutralPionsOld->GetBranch("nSubRootinos")){
	neutralPionsOld->SetBranchAddress("nSubRootinos", nSubRootinos);
	neutralPions->Branch("nSubRootinos", nSubRootinos, "nSubRootinos[nParticles]/I");
      }
      if(neutralPionsOld->GetBranch("nSubPhotons")){
	neutralPionsOld->SetBranchAddress("nSubPhotons", nSubPhotons);
	neutralPions->Branch("nSubPhotons", nSubPhotons, "nSubPhotons[nParticles]/I");
      }
      if(neutralPionsOld->GetBranch("nSubChargedPi")){
	neutralPionsOld->SetBranchAddress("nSubChargedPi", nSubChargedPions);
	neutralPions->Branch("nSubChargedPions", nSubChargedPions, "nSubChargedPions[nParticles]/I");
      }
      if(neutralPionsOld->GetBranch("subRootinos")){
	neutralPionsOld->SetBranchAddress("subRootinos", &subRootinos);
	neutralPions->Branch("subRootinos", &subRootinos, 32, 0);
      }
      if(neutralPionsOld->GetBranch("subPhotons")){
	neutralPionsOld->SetBranchAddress("subPhotons", &subPhotons);
	neutralPions->Branch("subPhotons", &subPhotons, 32, 0);
      }
      if(neutralPionsOld->GetBranch("subChargedPi")){
	neutralPionsOld->SetBranchAddress("subChargedPi", &subChargedPions);
	neutralPions->Branch("subChargedPions", &subChargedPions, 32, 0);
      }
    }

    for(Int_t i=0; i<neutralPionsOld->GetEntries(); i++){
      neutralPionsOld->GetEvent(i);
      neutralPions->Fill();
    }
    neutralPions->Write();
  }

  if(etasOld){
    printf("Renaming eta to etas\n");

    TTree *etas = new TTree("etas", "etas");

    if(etasOld->GetBranch("nParticles")){
      etasOld->SetBranchAddress("nParticles", &nParticles);
      etas->Branch("nParticles", &nParticles, "nParticles/I");
      if(etasOld->GetBranch("particles.")){
	etasOld->SetBranchAddress("particles.", &particles);
	etas->Branch("particles", &particles, 32000, 0);
      }
      if(etasOld->GetBranch("time")){
	etasOld->SetBranchAddress("time", time);
	etas->Branch("time", time, "time[nParticles]/D");
      }
      if(etasOld->GetBranch("clusterSize")){
	etasOld->SetBranchAddress("clusterSize", clusterSize);
	etas->Branch("clusterSize", clusterSize, "clusterSize[nParticles]/I");
      }
      if(etasOld->GetBranch("Apparatus")){
	etasOld->SetBranchAddress("Apparatus", apparatus);
	etas->Branch("apparatus", apparatus, "apparatus[nParticles]/I");
      }
      if(etasOld->GetBranch("d_E")){
	etasOld->SetBranchAddress("d_E", vetoEnergy);
	etas->Branch("vetoEnergy", vetoEnergy, "vetoEnergy[nParticles]/D");
      }
      if(etasOld->GetBranch("WC0_E")){
	etasOld->SetBranchAddress("WC0_E", MWPC0Energy);
	etas->Branch("MWPC0Energy", MWPC0Energy, "MWPC0Energy[nParticles]/D");
      }
      if(etasOld->GetBranch("WC1_E")){
	etasOld->SetBranchAddress("WC1_E", MWPC1Energy);
	etas->Branch("MWPC1Energy", MWPC1Energy, "MWPC1Energy[nParticles]/D");
      }
      if(etasOld->GetBranch("nSubParticles")){
	etasOld->SetBranchAddress("nSubParticles", nSubParticles);
	etas->Branch("nSubParticles", nSubParticles, "nSubParticles[nParticles]/I");
      }
      if(etasOld->GetBranch("nSubRootinos")){
	etasOld->SetBranchAddress("nSubRootinos", nSubRootinos);
	etas->Branch("nSubRootinos", nSubRootinos, "nSubRootinos[nParticles]/I");
      }
      if(etasOld->GetBranch("nSubPhotons")){
	etasOld->SetBranchAddress("nSubPhotons", nSubPhotons);
	etas->Branch("nSubPhotons", nSubPhotons, "nSubPhotons[nParticles]/I");
      }
      if(etasOld->GetBranch("nSubChargedPi")){
	etasOld->SetBranchAddress("nSubChargedPi", nSubChargedPions);
	etas->Branch("nSubChargedPions", nSubChargedPions, "nSubChargedPions[nParticles]/I");
      }
      if(etasOld->GetBranch("subRootinos")){
	etasOld->SetBranchAddress("subRootinos", &subRootinos);
	etas->Branch("subRootinos", &subRootinos, 32, 0);
      }
      if(etasOld->GetBranch("subPhotons")){
	etasOld->SetBranchAddress("subPhotons", &subPhotons);
	etas->Branch("subPhotons", &subPhotons, 32, 0);
      }
      if(etasOld->GetBranch("subChargedPi")){
	etasOld->SetBranchAddress("subChargedPi", &subChargedPions);
	etas->Branch("subChargedPions", &subChargedPions, 32, 0);
      }
    }

    for(Int_t i=0; i<etasOld->GetEntries(); i++){
      etasOld->GetEvent(i);
      etas->Fill();
    }
    etas->Write();
  }

  if(etaPrimesOld){
    printf("Renaming eta' to etaPrimes\n");

    TTree *etaPrimes = new TTree("etaPrimes", "etaPrimes");

    if(etaPrimesOld->GetBranch("nParticles")){
      etaPrimesOld->SetBranchAddress("nParticles", &nParticles);
      etaPrimes->Branch("nParticles", &nParticles, "nParticles/I");
      if(etaPrimesOld->GetBranch("particles.")){
	etaPrimesOld->SetBranchAddress("particles.", &particles);
	etaPrimes->Branch("particles", &particles, 32000, 0);
      }
      if(etaPrimesOld->GetBranch("time")){
	etaPrimesOld->SetBranchAddress("time", time);
	etaPrimes->Branch("time", time, "time[nParticles]/D");
      }
      if(etaPrimesOld->GetBranch("clusterSize")){
	etaPrimesOld->SetBranchAddress("clusterSize", clusterSize);
	etaPrimes->Branch("clusterSize", clusterSize, "clusterSize[nParticles]/I");
      }
      if(etaPrimesOld->GetBranch("Apparatus")){
	etaPrimesOld->SetBranchAddress("Apparatus", apparatus);
	etaPrimes->Branch("apparatus", apparatus, "apparatus[nParticles]/I");
      }
      if(etaPrimesOld->GetBranch("d_E")){
	etaPrimesOld->SetBranchAddress("d_E", vetoEnergy);
	etaPrimes->Branch("vetoEnergy", vetoEnergy, "vetoEnergy[nParticles]/D");
      }
      if(etaPrimesOld->GetBranch("WC0_E")){
	etaPrimesOld->SetBranchAddress("WC0_E", MWPC0Energy);
	etaPrimes->Branch("MWPC0Energy", MWPC0Energy, "MWPC0Energy[nParticles]/D");
      }
      if(etaPrimesOld->GetBranch("WC1_E")){
	etaPrimesOld->SetBranchAddress("WC1_E", MWPC1Energy);
	etaPrimes->Branch("MWPC1Energy", MWPC1Energy, "MWPC1Energy[nParticles]/D");
      }
      if(etaPrimesOld->GetBranch("nSubParticles")){
	etaPrimesOld->SetBranchAddress("nSubParticles", nSubParticles);
	etaPrimes->Branch("nSubParticles", nSubParticles, "nSubParticles[nParticles]/I");
      }
      if(etaPrimesOld->GetBranch("nSubRootinos")){
	etaPrimesOld->SetBranchAddress("nSubRootinos", nSubRootinos);
	etaPrimes->Branch("nSubRootinos", nSubRootinos, "nSubRootinos[nParticles]/I");
      }
      if(etaPrimesOld->GetBranch("nSubPhotons")){
	etaPrimesOld->SetBranchAddress("nSubPhotons", nSubPhotons);
	etaPrimes->Branch("nSubPhotons", nSubPhotons, "nSubPhotons[nParticles]/I");
      }
      if(etaPrimesOld->GetBranch("nSubChargedPi")){
	etaPrimesOld->SetBranchAddress("nSubChargedPi", nSubChargedPions);
	etaPrimes->Branch("nSubChargedPions", nSubChargedPions, "nSubChargedPions[nParticles]/I");
      }
      if(etaPrimesOld->GetBranch("subRootinos")){
	etaPrimesOld->SetBranchAddress("subRootinos", &subRootinos);
	etaPrimes->Branch("subRootinos", &subRootinos, 32, 0);
      }
      if(etaPrimesOld->GetBranch("subPhotons")){
	etaPrimesOld->SetBranchAddress("subPhotons", &subPhotons);
	etaPrimes->Branch("subPhotons", &subPhotons, 32, 0);
      }
      if(etaPrimesOld->GetBranch("subChargedPi")){
	etaPrimesOld->SetBranchAddress("subChargedPi", &subChargedPions);
	etaPrimes->Branch("subChargedPions", &subChargedPions, 32, 0);
      }
    }

    for(Int_t i=0; i<etaPrimesOld->GetEntries(); i++){
      etaPrimesOld->GetEvent(i);
      etaPrimes->Fill();
    }
    etaPrimes->Write();
  }

  if(fOld.GetListOfKeys()->Contains("CheckCB")){
    fNew.mkdir("CheckCB");
    fNew.cd("CheckCB");

    TH2F *Check_CBdE_E = (TH2F*)fOld.Get("CheckCB/Check_CBdE_E");
    Check_CBdE_E->Write();
    TH2F *Check_CBPhiCorr = (TH2F*)fOld.Get("CheckCB/Check_CBPhiCorr");	
    Check_CBPhiCorr->Write();
    TH2F *Check_CBdE_E_1PID = (TH2F*)fOld.Get("CheckCB/Check_CBdE_E_1PID");	
    Check_CBdE_E_1PID->Write();
    TH2F *Check_CBPhiCorr_1PID = (TH2F*)fOld.Get("CheckCB/Check_CBPhiCorr_1PID");	
    Check_CBPhiCorr_1PID->Write();
    TH2F *Check_CBdE_E_pi0 = (TH2F*)fOld.Get("CheckCB/Check_CBdE_E_pi0");	
    Check_CBdE_E_pi0->Write();
    TH2F *Check_CBPhiCorr_pi0 = (TH2F*)fOld.Get("CheckCB/Check_CBPhiCorr_pi0");	
    Check_CBPhiCorr_pi0->Write();

    fNew.cd();
  }

  if(fOld.GetListOfKeys()->Contains("CheckTAPS")){
    fNew.mkdir("CheckTAPS");
    fNew.cd("CheckTAPS");

    TH2F *Check_TAPSdE_E = (TH2F*)fOld.Get("CheckTAPS/Check_TAPSdE_E");
    Check_TAPSdE_E->Write();
    TH2F *Check_TAPSPhiCorr = (TH2F*)fOld.Get("CheckTAPS/Check_TAPSPhiCorr");
    Check_TAPSPhiCorr->Write();
    TH2F *Check_TAPSdE_E_1Veto = (TH2F*)fOld.Get("CheckTAPS/Check_TAPSdE_E_1PID");
    Check_TAPSdE_E_1Veto->SetName("Check_TAPSdE_E_1Veto");
    Check_TAPSdE_E_1Veto->Write();
    TH2F *Check_TAPSPhiCorr_1Veto = (TH2F*)fOld.Get("CheckTAPS/Check_TAPSPhiCorr_1Veto");
    Check_TAPSPhiCorr_1Veto->Write();

    fNew.cd();
  }

  if(fOld.GetListOfKeys()->Contains("CheckHitPatterns")){
    fNew.mkdir("CheckHitPatterns");
    fNew.cd("CheckHitPatterns");	

    TH2F *Check_CBHits = (TH2F*)fOld.Get("CheckHitPatterns/Check_CBHits");
    Check_CBHits->Write();
    TH2F *Check_CBADCHits = (TH2F*)fOld.Get("CheckHitPatterns/Check_CBADCHits");
    Check_CBADCHits->Write();
    TH2F *Check_CBTDCHits = (TH2F*)fOld.Get("CheckHitPatterns/Check_CBTDCHits");
    Check_CBTDCHits->Write();

    TH2F *Check_PIDHits = (TH2F*)fOld.Get("CheckHitPatterns/Check_PIDHits");
    Check_PIDHits->Write();
    TH2F *Check_PIDADCHits = (TH2F*)fOld.Get("CheckHitPatterns/Check_PIDADCHits");
    Check_PIDADCHits->Write();
    TH2F *Check_PIDTDCHits = (TH2F*)fOld.Get("CheckHitPatterns/Check_PIDTDCHits");
    Check_PIDTDCHits->Write();

    TH2F *Check_TAPSHits = (TH2F*)fOld.Get("CheckHitPatterns/Check_TAPSHits");
    Check_TAPSHits->Write();
    TH2F *Check_TAPSADCHits = (TH2F*)fOld.Get("CheckHitPatterns/Check_TAPSADCHits");
    Check_TAPSADCHits->Write();
    TH2F *Check_TAPSTDCHits = (TH2F*)fOld.Get("CheckHitPatterns/Check_TAPSTDCHits");
    Check_TAPSTDCHits->Write();

    TH2F *Check_VetoHits = (TH2F*)fOld.Get("CheckHitPatterns/Check_VetoHits");
    Check_VetoHits->Write();
    TH2F *Check_VetoADCHits = (TH2F*)fOld.Get("CheckHitPatterns/Check_VetoADCHits");
    Check_VetoADCHits->Write();
    TH2F *Check_VetoTDCHits = (TH2F*)fOld.Get("CheckHitPatterns/Check_VetoTDCHits");
    Check_VetoTDCHits->Write();

    fNew.cd();
  }

  if(fOld.GetListOfKeys()->Contains("GoAT_File")){
    TNamed *GoAT_File = (TNamed*)fOld.Get("GoAT_File");
    GoAT_File->Write();
  }

  if(fOld.GetListOfKeys()->Contains("CountScalerValid")){
    TH1I *CountScalerValid = (TH1I*)fOld.Get("CountScalerValid");
    CountScalerValid->Write();
  }

  fNew.Close();
  fOld.Close();

}
