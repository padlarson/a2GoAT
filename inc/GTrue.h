#ifndef __GTrue_h__
#define __GTrue_h__

#include "GTreeManager.h"
#include "GTreePluto.h"
#include "TLorentzVector.h"
#include "GH1.h"

#define GTrue_MaxEntries  16

class GTrue
{

private:
//    GH1                 vX;
    Double_t            TrueBeamEnergy;
    // weight of each event, default weight = 1
    Double_t            weight;
    // LorentzVectors of all produced true particles;
    TLorentzVector     neutron;
    TLorentzVector     proton;
    TLorentzVector     etaprime;
    TLorentzVector     omega;
    TLorentzVector     eta;
    TLorentzVector     chpi[GTrue_MaxEntries];
    TLorentzVector     pi0[GTrue_MaxEntries];
    TLorentzVector     gamma[GTrue_MaxEntries];
    TLorentzVector     electron;
    TLorentzVector     positron;
    // Incrementing nr of each particle type;
    UInt_t              nchpi;
    UInt_t              npi0;
    UInt_t              ngamma;
    // True vertex
    TVector3            vertex;

    void                SetWeight(const Double_t w )   { weight = w;}


public:
    GTrue();
    virtual ~GTrue();
    void Init(GTreePluto& pluto);
    void Start(GTreePluto& pluto, GTreeA2Geant& geant);

            Double_t    GetTrueBeamEnergy()                     { return TrueBeamEnergy;}
    const   Double_t    GetTrueBeamEnergy()     const           { return TrueBeamEnergy;}
    // proton LorentzVector
            TLorentzVector& GetTrueProtonLV()           { return proton;}
    const   TLorentzVector& GetTrueProtonLV()   const   { return proton;}
    // neutron LorentzVector
            TLorentzVector& GetTrueNeutronLV()          { return neutron;}
    const   TLorentzVector& GetTrueNeutronLV()  const   { return neutron;}
    // eta' LorentzVector
            TLorentzVector& GetTrueEtaPrimeLV()         { return etaprime;}
    const   TLorentzVector& GetTrueEtaPrimeLV() const   { return etaprime;}
    // omega LorentzVector
            TLorentzVector& GetTrueOmegaLV()            { return omega;}
    const   TLorentzVector& GetTrueOmegaLV()    const   { return omega;}
    // eta LorentzVector
            TLorentzVector& GetTrueEtaLV()              { return eta;}
    const   TLorentzVector& GetTrueEtaLV()      const   { return eta;}
    // pi-/+ LorentzVector
            TLorentzVector& GetTrueChargedPiLV(const Int_t index )       { return chpi[index];}
    const   TLorentzVector& GetTrueChargedPiLV(const Int_t index ) const { return chpi[index];}
    // pi0 LorentzVector
            TLorentzVector& GetTrueNeutralPiLV( Int_t index )       { return pi0[index];}
    const   TLorentzVector& GetTrueNeutralPiLV( Int_t index ) const { return pi0[index];}
    // gamma LorentzVector
            TLorentzVector& GetTrueGammaLV( Int_t index )       { return gamma[index];}
    const   TLorentzVector& GetTrueGammaLV( Int_t index ) const { return gamma[index];}
    // electron LorentzVector
            TLorentzVector& GetTrueElectronLV()       { return electron;}
    const   TLorentzVector& GetTrueElectronLV() const { return electron;}
    // positron LorentzVector
            TLorentzVector& GetTruePositronLV()       { return positron;}
    const   TLorentzVector& GetTruePositronLV() const { return positron;}
//
    Double_t    GetWeight() const   { return weight;}

    UInt_t      GetNchpi()  const   { return nchpi;}
    UInt_t      GetNpi0()   const   { return npi0;}
    UInt_t      GetNgamma() const   { return ngamma;}

    // Get Vertices
    const TVector3&    GetTrueVertex()   const   { return vertex;}

};

#endif //__GTrue_h__
