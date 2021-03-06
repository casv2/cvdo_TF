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

import numpy as np
import lattice_tetragonal

# the current model
import model 

#c_over_a = 1.8
#a0 = (16.0*2*2/sqrt(3.0)/c_over_a)**(1.0/3.0)# initial guess at lattice constant, cell will be relaxed below

#c_over_a = 1.588
#a0 = 2.95

def create_omega_custom(c_vs_a, a, z):
   import quippy

   #a = (2 * V / (3.0**(0.5) * c_vs_a))**(1.0/3.0)
   c = c_vs_a * a

   lattice = []
   lattice.append([3.0**(0.5) /2.0 * a,-a/2.0,0])
   lattice.append([3.0**(0.5) /2.0 * a, a/2.0,0])
   lattice.append([0,0,c])
   lattice = np.transpose(lattice)
   unitcell = quippy.Atoms(n=0, lattice=lattice)

   pos = []
   pos.append([3.0**(0.5) /6.0 * a,0,0.0])
   pos.append([3.0**(0.5) /2.0 * a,0,c/2.0])
   pos.append([3.0**(0.5) * 5.0/6.0 * a,0,c/2.0])

   for i in range(0,len(pos)):
      unitcell.add_atoms(pos[i],z)

   return unitcell

# set up the a
#bulk = HexagonalClosedPacked(symbol='Ti', latticeconstant=(a0,a0*c_over_a))

c_vs_a = 0.610
a = 4.630

bulk = create_omega_custom(c_vs_a, a, 22)

(c11, c33, c12, c13, c44, c66, E_vs_V) = lattice_tetragonal.do_lattice(bulk, elastic=True)#, tol=1.0e-2)

# dictionary of computed properties - this is output of this test, to
#   be compared with other models
#properties = {'hcp_E_vs_V': E_vs_V }

properties = {'omega_c11': c11, 'omega_c33': c33, 'omega_c12': c12, 'omega_c13': c13, 'omega_c44': c44, 'omega_c66': c66, 'omega_E_vs_V': E_vs_V }
