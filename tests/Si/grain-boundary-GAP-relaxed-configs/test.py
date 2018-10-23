# This script defines a test case which computes one or more physical
# properties with a given model
#
# INPUTS:
#   model.calculator -- an ase.calculator.Calculator instance
#     this script can assume the calculator is checkpointed.
#
# OUTPUTS:
#   properties -- dictionary of key/value pairs corresponding
#     to physical quantities computed by this test

# standard ASE structure generation routines

import ase.io, os

# set of utility routines specific this this model/testing framework

# the current model
import model

#ats = ase.io.read(os.path.join(os.path.dirname(__file__),"gb.0.25_0.0-relaxed.opt.xyz"),":")
#Es = []
#for at in ats:
    # e0 = at.get_potential_energy()
#    at.set_calculator(model.calculator)
#    e = at.get_potential_energy()
#    Es.append(e)
at = ase.io.read(os.path.join(os.path.dirname(__file__),"gb.0.25_0.0-relaxed.opt.xyz"))
at.set_calculator(model.calculator)
energy = at.get_potential_energy()

properties = {'grain-boundary-GAP-relaxed':energy}
