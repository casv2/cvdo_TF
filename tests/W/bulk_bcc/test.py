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
from ase.lattice.cubic import BodyCenteredCubic

import lattice_cubic

# the current model
import model 

a0 = 3.16 # initial guess at lattice constant, cell will be relaxed below

# set up the a
bulk = BodyCenteredCubic(symbol='W', latticeconstant=a0)

(c11, c12, c44, E_vs_V) = lattice_cubic.do_lattice(bulk, elastic=True)

#properties = {'bcc_E_vs_V': E_vs_V }

# dictionary of computed properties - this is output of this test, to
#   be compared with other models
properties = {'bcc_c11': c11, 'bcc_c12': c12, 'bcc_c44': c44, 'bcc_E_vs_V': E_vs_V }
