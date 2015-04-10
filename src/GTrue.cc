#include "GTrue.h"
#include "GTreeA2Geant.h"
#include "GTreePluto.h"
#include "PParticle.h"


GTrue::GTrue():
 //   vX("vertexX", "vertexX", 2050, -1025, 1025),
    weight(1.0),
    TrueBeamEnergy(-10000.0),
    vertex(-1000,-1000,-1000),
    ngamma(0),
    npi0(0),
    nchpi(0)
{
    weight = 1;
}

GTrue::~GTrue()
{
}


void GTrue::Start(GTreePluto& pluto, GTreeA2Geant& geant)
{
    nchpi = 0;
    npi0 = 0;
    ngamma = 0;
    TrueBeamEnergy = -10000.0;

    vertex  = geant.GetVertex();

    TrueBeamEnergy = pluto.GetTrueBeam().E();
//    TrueBeamEnergy = pluto.GetAllParticles().front()->P();

    /*for(std::list<PParticle*>::iterator i = pluto.GetAllParticles().begin(); i!=pluto.GetAllParticles().end(); i++)
    {
        if((*i)->Is("gamma")){
            //gamma[ngamma] = (const TLorentzVector*)(*i);
            ngamma++;
            continue;
        }
    }*/


    for(int i = 0; i < pluto.GetAllParticles().size(); i++)
    {
        if(pluto.GetMCTrue(i)->Is("n")){
            neutron = pluto.GetTrueP4(i);
            continue;
        }
        if(pluto.GetMCTrue(i)->Is("p")){
            proton = pluto.GetTrueP4(i);
            continue;
        }
        if(pluto.GetMCTrue(i)->Is("eta'")){
            etaprime = pluto.GetTrueP4(i);
            continue;
        }
        if(pluto.GetMCTrue(i)->Is("w")){
            omega = pluto.GetTrueP4(i);
            continue;
        }
        if(pluto.GetMCTrue(i)->Is("eta")){
            eta = pluto.GetTrueP4(i);
            continue;
        }
        if(pluto.GetMCTrue(i)->Is("pi-") || pluto.GetMCTrue(i)->Is("pi+")){
            chpi[nchpi] = pluto.GetTrueP4(i);
            nchpi++;
            continue;
        }
        if(pluto.GetMCTrue(i)->Is("pi0")){
            pi0[npi0] = pluto.GetTrueP4(i);
            npi0++;
            continue;
        }
        if(pluto.GetMCTrue(i)->Is("gamma")){
            gamma[ngamma] = pluto.GetTrueP4(i);
            ngamma++;
            continue;
        }
    }
}
