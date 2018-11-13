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
from math import sqrt

import lattice_tetragonal

# the current model
import model 

#c_over_a = 1.8
#a0 = (16.0*2*2/sqrt(3.0)/c_over_a)**(1.0/3.0)# initial guess at lattice constant, cell will be relaxed below

c_over_a = 1.588
a0 = 2.95
#a0 = 2.95

# set up the a
bulk = HexagonalClosedPacked(symbol='Ti', latticeconstant=(a0,a0*c_over_a)) 

(c11, c33, c12, c13, c44, c66, E_vs_V) = lattice_tetragonal.do_lattice(bulk, elastic=True)#, tol=1.0e-2)

# dictionary of computed properties - this is output of this test, to
#   be compared with other models
#properties = {'hcp_E_vs_V': E_vs_V }

properties = {'hcp_c11': c11, 'hcp_c33': c33, 'hcp_c12': c12, 'hcp_c13': c13, 'hcp_c44': c44, 'hcp_c66': c66, 'hcp_E_vs_V': E_vs_V }
