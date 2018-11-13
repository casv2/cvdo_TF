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

import os

# standard ASE structure generation routines
from ase.lattice.cubic import Diamond
from ase import Atom
from ase.optimize import FIRE, MDMin
#from quippy.neb import NEB
from ase.neb import NEB
import numpy as np
#from quippy.atoms import Atoms

import ase.io, sys

# set of utility routines specific this this model/testing framework
from utilities import relax_atoms, relax_atoms_cell

# the current model
import model 

a0 = 5.44 # initial guess at lattice constant, cell will be relaxed below
tol = 1e-3 # maximum force following relaxtion [eV/A]
N = 3 # number of unit cells in each direction
fmax = 0.05
n_images = 5
k = 1.0 # NEB spring constant, in eV/A^2

# set up the a
bulk = Diamond(symbol='Si', latticeconstant=a0)

# specify that we will use model.calculator to compute forces, energies and stresses
bulk.set_calculator(model.calculator)

# use one of the routines from utilities module to relax the initial
# unit cell and atomic positions
print 'starting bulk cell relaxation'
bulk = relax_atoms_cell(bulk, tol=tol, traj_file=None)
print 'bulk cell relaxation done, energy=', bulk.get_potential_energy()
print bulk.cell
e_bulk_per_atom = bulk.get_potential_energy()/len(bulk)

bulk *= (N, N, N)

Nat = bulk.get_number_of_atoms()
initial_struct = bulk.copy()
initial_struct.set_calculator(bulk.get_calculator())
# add an atom to introduce an interstitial
initial_struct.append(Atom('Si', (0.001, 0.002, 5.44/2.0+0.003)))

cell = initial_struct.get_cell()
print "shift ",-initial_struct.positions[len(initial_struct)-1,0]+sum(cell[:,0])/2.0, -initial_struct.positions[len(initial_struct)-1,1]+sum(cell[:,1])/2.0, -initial_struct.positions[len(initial_struct)-1,2]+sum(cell[:,2])/2.0
initial_struct.positions[:,0] += -initial_struct.positions[len(initial_struct)-1,0]+sum(cell[:,0])/2.0
initial_struct.positions[:,1] += -initial_struct.positions[len(initial_struct)-1,1]+sum(cell[:,1])/2.0
initial_struct.positions[:,2] += -initial_struct.positions[len(initial_struct)-1,2]+sum(cell[:,2])/2.0
initial_struct.wrap()

# start final from unrelaxed initial - easier to keep interstitial at high symmetry pos
final_struct = initial_struct.copy()

# relax atom positions, holding cell fixed
print 'starting initial interstitial cell relaxation'
initial_struct = relax_atoms(initial_struct, tol=tol, traj_file="model-"+model.name+"-test-interstitial-tetrahedral-initial.opt.xyz")
print 'initial interstitial cell relaxation finished, energy=', initial_struct.get_potential_energy()

final_struct.positions[len(initial_struct)-1,:] += [-a0/4.0, a0/4.0, -a0/4.0]
final_struct.set_calculator(bulk.get_calculator())

# relax atom positions, holding cell fixed
final_struct = relax_atoms(final_struct, tol=tol, traj_file="model-"+model.name+"-test-interstitial-tetrahedral-final.opt.xyz")

# do neb
# make chain
images = [ initial_struct ]
images += [ initial_struct.copy() for i in range(n_images) ]
images += [ final_struct ]

neb = NEB(images, k=k)
neb.interpolate()
for (i,img) in enumerate(images):
   img.set_calculator(model.calculator)
   print "NEB start energy ", img.get_potential_energy()
   ase.io.write(sys.stdout, img, 'extxyz')
   ase.io.write("model-"+model.name+"-test-interstitial.neb-start-%02d.xyz" % i, img, format='extxyz')

# optimizer = FIRE(neb)
optimizer = MDMin(neb)
optimizer.run(fmax=fmax, steps=300)

for (i,img) in enumerate(images):
   ase.io.write("model-"+model.name+"-test-interstitial.neb-end-%02d.xyz" % i, img, format='extxyz')

Es = []
for img in images:
   print "NEB end energy ", img.get_potential_energy()
   Es.append(img.get_potential_energy() - e_bulk_per_atom*len(img))
   ase.io.write(sys.stdout, img, 'extxyz')

# dictionary of computed properties - this is output of this test, to
#   be compared with other models
properties = {'tetrahedral_interstitial_path':
                Es}
