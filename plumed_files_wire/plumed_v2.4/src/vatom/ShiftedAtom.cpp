/*
    created: Jan 3, 2018
    author: Chenghan Li, lch004218@gmail.com
    
    shift the input atom by offset, i.e.
    \vec{r}_\text{shifted} = \vec{r}_\text{input} + \vec{r}_\text{offset}
*/
#include "ActionWithVirtualAtom.h"
#include "ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"

#include <vector>
#include <cmath>

using namespace std;

namespace PLMD {
namespace vatom {

class ShiftedAtom : public ActionWithVirtualAtom {
    vector<double> offset;
public:
    explicit ShiftedAtom(const ActionOptions&ao);
    void calculate();
    static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(ShiftedAtom,"SHIFTEDATOM") 

void ShiftedAtom::registerKeywords(Keywords& keys) {
    ActionWithVirtualAtom::registerKeywords(keys);
    keys.add("atoms","ATOM","Atom to be rotated");
    keys.add("compulsory","OFFSET","The offset vector");
}

ShiftedAtom::ShiftedAtom(const ActionOptions&ao):
    Action(ao),
    ActionWithVirtualAtom(ao)
{
    vector<AtomNumber> atom;
    parseAtomList("ATOM",atom);
    if(atom.size()!=1) error("should only be one atom specified");

    parseVector("OFFSET", offset);
    checkRead();

    if(offset.size() != 3) error("OFFSET should be a 3-d vector");

    log.printf("  the offset vector is :\n");
    log<<"  "<<offset[0]<<" "<<offset[1]<<" "<<offset[2]<<"\n";
    log.printf("\n");
    requestAtoms(atom);
}

void ShiftedAtom::calculate() {
    Vector pos;
    vector<Tensor> deriv(getNumberOfAtoms());
    for(int i = 0; i < 3; ++i) 
       pos[i] = getPosition(0)[i] + offset[i];
    deriv[0] = Tensor::identity();

    setPosition(pos);
    setMass(getMass(0));
    //setCharge(getCharge(0));
    setAtomsDerivatives(deriv);
}

}
}
