#include "GTreeDetectorHits.h"



GTreeDetectorHits::GTreeDetectorHits(GTreeManager *Manager)    :
    GTree(Manager, TString("detectorHits")),
    nNaIHits(0),
    nPIDHits(0),
    nMWPCHits(0),
    nBaF2Hits(0),
    nVetoHits(0)
{
}

GTreeDetectorHits::~GTreeDetectorHits()
{
}

void    GTreeDetectorHits::SetBranchAdresses()
{
    inputTree->SetBranchAddress("nNaIHits", &nNaIHits);
    inputTree->SetBranchAddress("NaIHits", NaIHits);
    if(inputTree->GetBranch("NaICluster")) inputTree->SetBranchAddress("NaICluster", NaICluster);
    if(inputTree->GetBranch("NaIEnergy")) inputTree->SetBranchAddress("NaIEnergy", NaIEnergy);
    if(inputTree->GetBranch("NaITime")) inputTree->SetBranchAddress("NaITime", NaITime);
    inputTree->SetBranchAddress("nPIDHits", &nPIDHits);
    inputTree->SetBranchAddress("PIDHits", PIDHits);
    if(inputTree->GetBranch("PIDEnergy")) inputTree->SetBranchAddress("PIDEnergy", PIDEnergy);
    if(inputTree->GetBranch("PIDTime")) inputTree->SetBranchAddress("PIDTime", PIDTime);
    inputTree->SetBranchAddress("nMWPCHits", &nMWPCHits);
    inputTree->SetBranchAddress("MWPCHits", MWPCHits);
    inputTree->SetBranchAddress("nBaF2Hits", &nBaF2Hits);
    inputTree->SetBranchAddress("BaF2Hits", BaF2Hits);
    if(inputTree->GetBranch("BaF2Cluster")) inputTree->SetBranchAddress("BaF2Cluster", BaF2Cluster);
    if(inputTree->GetBranch("BaF2Energy")) inputTree->SetBranchAddress("BaF2Energy", BaF2Energy);
    if(inputTree->GetBranch("BaF2Time")) inputTree->SetBranchAddress("BaF2Time", BaF2Time);
    inputTree->SetBranchAddress("nVetoHits", &nVetoHits);
    inputTree->SetBranchAddress("VetoHits", VetoHits);
    if(inputTree->GetBranch("VetoEnergy")) inputTree->SetBranchAddress("VetoEnergy", VetoEnergy);
    if(inputTree->GetBranch("VetoTime")) inputTree->SetBranchAddress("VetoTime", VetoTime);
}

void    GTreeDetectorHits::SetBranches()
{
    outputTree->Branch("nNaIHits", &nNaIHits, "nNaIHits/I");
    outputTree->Branch("NaIHits", NaIHits, "NaIHits[nNaIHits]/I");
    if(inputTree->GetBranch("NaICluster")) outputTree->Branch("NaICluster", NaICluster, "NaICluster[nNaIHits]/I");
    if(inputTree->GetBranch("NaIEnergy")) outputTree->Branch("NaIEnergy", NaIEnergy, "NaIEnergy[nNaIHits]/D");
    if(inputTree->GetBranch("NaITime")) outputTree->Branch("NaITime", NaITime, "NaITime[nNaIHits]/D");
    outputTree->Branch("nPIDHits", &nPIDHits, "nPIDHits/I");
    outputTree->Branch("PIDHits", PIDHits, "PIDHits[nPIDHits]/I");
    if(inputTree->GetBranch("PIDEnergy")) outputTree->Branch("PIDEnergy", PIDEnergy, "PIDEnergy[nPIDHits]/D");
    if(inputTree->GetBranch("PIDTime")) outputTree->Branch("PIDTime", PIDTime, "PIDTime[nPIDHits]/D");
    outputTree->Branch("nMWPCHits", &nMWPCHits, "nMWPCHits/I");
    outputTree->Branch("MWPCHits", MWPCHits, "MWPCHits[nMWPCHits]/I");
    outputTree->Branch("nBaF2Hits", &nBaF2Hits, "nBaF2Hits/I");
    outputTree->Branch("BaF2Hits", BaF2Hits, "BaF2Hits[nBaF2Hits]/I");
    if(inputTree->GetBranch("BaF2Cluster")) outputTree->Branch("BaF2Cluster", BaF2Cluster, "BaF2Cluster[nBaF2Hits]/I");
    if(inputTree->GetBranch("BaF2Energy")) outputTree->Branch("BaF2Energy", BaF2Energy, "BaF2Energy[nBaF2Hits]/D");
    if(inputTree->GetBranch("BaF2Time")) outputTree->Branch("BaF2Time", BaF2Time, "BaF2Time[nBaF2Hits]/D");
    outputTree->Branch("nVetoHits", &nVetoHits, "nVetoHits/I");
    outputTree->Branch("VetoHits", VetoHits, "VetoHits[nVetoHits]/I");
    if(inputTree->GetBranch("VetoEnergy")) outputTree->Branch("VetoEnergy", VetoEnergy, "VetoEnergy[nVetoHits]/D");
    if(inputTree->GetBranch("VetoTime")) outputTree->Branch("VetoTime", VetoTime, "VetoTime[nVetoHits]/D");
}


