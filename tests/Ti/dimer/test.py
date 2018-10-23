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

dimer = Atoms('Ti2',
               positions=[[0, 0, 0],
                          [0, 0, 1.2]], cell=[[20,0,0],[0,20,0],[0,0,20]])


dimer.set_calculator(model.calculator)
    
p = dimer.get_positions()
e = []
rr = []
for r in np.arange(1.2,6.0,0.1):
    p[1,2] = r
    dimer.set_positions(p)
    e.append(dimer.get_potential_energy())
    rr.append(r)

properties = {'dimer distance':rr, 'dimer energy': e}
