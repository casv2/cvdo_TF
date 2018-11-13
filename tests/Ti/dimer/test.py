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

from ase import Atoms
import model
import numpy as np
from ase.lattice.cubic import BodyCenteredCubic

single_atom = Atoms('Ti', cell=[[20,0,0],[0,20,0],[0,0,20]])
single_atom.set_calculator(model.calculator)
e0 = single_atom.get_potential_energy()

a0 = 3.32

bulk = BodyCenteredCubic(symbol='Ti', latticeconstant=a0)
rnn = np.linalg.norm(bulk.positions[1])

dimer = Atoms('Ti2',
               positions=[[0, 0, 0],
                          [0, 0, 1.5]], cell=[[20,0,0],[0,20,0],[0,0,20]])


dimer.set_calculator(model.calculator)

    
p = dimer.get_positions()
e = []
rr = []

for r in np.arange(2,7.0,0.1):
    p[1,2] = r
    dimer.set_positions(p)
    e.append(dimer.get_potential_energy() - (2*e0))
    rr.append(r)

properties = {'dimer distance':rr, 'dimer energy': e, 'rnn' : rnn}