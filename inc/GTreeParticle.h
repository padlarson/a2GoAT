#ifndef __GTreeParticle_h__
#define __GTreeParticle_h__


#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "Rtypes.h"
#include "GTree.h"

#define GTreeParticle_MAX 128

class   GTreeMeson;

class  GTreeParticle    : public GTree
{
private:
    Int_t               nParticles;
    TClonesArray*       particles;	// reconstructed
    Int_t               apparatus[GTreeParticle_MAX];
    Double_t            time[GTreeParticle_MAX];
    Int_t               clusterSize[GTreeParticle_MAX];
    //Charged detector energies
    Double_t            vetoEnergy[GTreeParticle_MAX];
    Double_t            MWPC0Energy[GTreeParticle_MAX];
    Double_t            MWPC1Energy[GTreeParticle_MAX];
    Int_t               trackIndex[GTreeParticle_MAX];  // index of the corresponding tack in the track list, -1 => No track


protected:
    virtual void    SetBranchAdresses();
    virtual void    SetBranches();
            void    PrintParticle(const Int_t i) const;

public:
    GTreeParticle(GTreeManager *Manager, const TString& _Name);
    virtual ~GTreeParticle();

            void            AddParticle(const TLorentzVector& vec, const Int_t _Apparatus = 0, const Double_t _vetoEnergy = 0, const Double_t _MWPC0Energy = 0, const Double_t _MWPC1Energy = 0, const Double_t _Time = 0, const Int_t _ClusterSize = 0, const Int_t _trackIndex = -1);
    virtual void            Clear() {nParticles = 0; particles->Clear();}
            Int_t           GetApparatus(const Int_t index)     const	{return apparatus[index];}
            Int_t           GetClusterSize(const Int_t index)   const 	{return clusterSize[index];}
            Double_t        GetVetoEnergy(const Int_t index)    const	{return vetoEnergy[index];}
            Int_t           GetNParticles()                     const	{return nParticles;}
    const	Double_t*       GetTime()                           const	{return time;}
            Double_t        GetTime(const Int_t index)          const	{return time[index];}
    const	Double_t*       GetMWPC0Energy()                    const	{return MWPC0Energy;}
            Double_t        GetMWPC0Energy(const Int_t index)   const	{return MWPC0Energy[index];}
    const	Double_t*       GetMWPC1Energy()                    const	{return MWPC1Energy;}
            Double_t        GetMWPC1Energy(const Int_t index)   const	{return MWPC1Energy[index];}

            Int_t           GetTrackIndex(const Int_t index)    const   {return trackIndex[index];}
    const   Int_t*          GetTrackIndex()                     const   {return trackIndex; }

    inline          TLorentzVector& Particle(const Int_t particle);
    inline  const   TLorentzVector& Particle(const Int_t particle) const;
    virtual void            Print() const;
            void            RemoveParticles(const Int_t nIndices, const Int_t* indices);
            void            RemoveAllParticles();
    virtual Bool_t          Write();


    friend  class GTreeMeson;
};


TLorentzVector&         GTreeParticle::Particle(const Int_t particle)
{
    if(particle>=nParticles)
    {
        std::cout << "ERROR: Try to access particle " << particle << " in " << GetName() << ". Only " << nParticles << " in tree." << std::endl;
        static  TLorentzVector  helpVectorUnity(0,0,0,0);
        helpVectorUnity.SetPxPyPzE(0,0,0,0);
        return helpVectorUnity;
    }
    return *((TLorentzVector*)particles->At(particle));
}

const   TLorentzVector& GTreeParticle::Particle(const Int_t particle) const
{
    if(particle>=nParticles)
    {
        std::cout << "ERROR: Try to access particle " << particle << " in " << GetName() << ". Only " << nParticles << " in tree." << std::endl;
        static  TLorentzVector  helpVectorUnity(0,0,0,0);
        helpVectorUnity.SetPxPyPzE(0,0,0,0);
        return helpVectorUnity;
    }
    return *((TLorentzVector*)particles->At(particle));
}


#endif
