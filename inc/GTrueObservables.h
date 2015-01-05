#ifndef GTRUEOBSERVABLES_H
#define GTRUEOBSERVABLES_H

#endif // GTRUEOBSERVABLES_H

#include <TLorentzVector.h>

#define MASS_PI0    134.9766
#define MASS_ETA    547.853
#define MASS_ETAP   957.78
#define MASS_PROTON 938.272046

using namespace std;

class  GTrueObservables : public GTreeA2Geant
{
private:
    // here one needs all vertices, energies, theta, phi
    // X, Y, m2pi0, tagger, angular distribution
protected:
    virtual void    SetBranchAdresses();
    virtual void    SetBranches();

public:
    GTrueObservables(GTreeManager *Manager, const TString& _Name);
    virtual ~GTrueObservables();

    virtual void            Clear();
};



#endif

