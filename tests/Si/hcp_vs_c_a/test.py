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
from ase.lattice.hexagonal import HexagonalClosedPacked
import ase.io
from math import sqrt
from utilities import relax_atoms_cell

# the current model
import model 
tol = 1.0e-3

data = []
for c_over_a_i in range(13):
    c_over_a = 0.6 + 0.1 * c_over_a_i
    a0 = (16.55*2*2/sqrt(3.0)/c_over_a)**(1.0/3.0)# initial guess at lattice constant, cell will be relaxed below

    bulk = HexagonalClosedPacked(symbol='Si', latticeconstant=(a0,a0*c_over_a))
    # bulk = relax_atoms_cell(bulk, tol=tol, traj_file="bulk.relax_traj.c_a_{}.extxyz".format(c_over_a), hydrostatic_strain=True, method='lbfgs')
    bulk = relax_atoms_cell(bulk, tol=tol, hydrostatic_strain=True, method='cg_n')
    data.append((c_over_a, bulk.get_potential_energy()/len(bulk), bulk.get_volume()/len(bulk)))
    print "got ", data[-1]
    ase.io.write("relaxed.c_a_{}.extxyz".format(c_over_a), bulk)

# dictionary of computed properties - this is output of this test, to
#   be compared with other models
properties = { 'c_over_a' : [ x[0] for x in data ], 
               'E' : [ x[1] for x in data ], 
               'V' : [ x[2] for x in data ] } 
