import ase.io, os

# set of utility routines specific this this model/testing framework

# the current model
import model

ats = ase.io.read(os.path.join(os.path.dirname(__file__),"model-GAP-vacancy-path.neb.xyz"),":")
Es = []
for at in ats:
    at.set_calculator(model.calculator)
    e = at.get_potential_energy()
    Es.append(e)
properties = {'vacancy-path' : Es}