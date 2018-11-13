from ase.build import bulk
from ase.dft.kpoints import ibz_points, bandpath
from ase.phonons import Phonons
from ase.io import read, write
from ase.lattice.hexagonal import HexagonalClosedPacked
from ase.optimize import BFGS

import model
# Setup crystal and EMT calculator
c_over_a = 1.588
a0 = 2.95
#a0 = 2.95

# set up the a
atoms = HexagonalClosedPacked(symbol='Ti', latticeconstant=(a0,a0*c_over_a)) 

#atoms.rattle(0.001)

#atoms.set_calculator(model.calculator)

#dyn = BFGS(atoms)
#dyn.run(fmax=0.0001)

#atoms = read("/Users/Cas/gits/gap-testing-framework/tests/Ti/phonon_hcp/hcp_relaxed.xyz")

#atoms = read("/Users/Cas/gits/gap-testing-framework/tests/Ti/phonon_hcp/hcp_init_relaxed.xyz")

# Phonon calculator
N = 3
ph = Phonons(atoms, model.calculator, supercell=(N, N, N), delta=0.01)
ph.run()

# Read forces and assemble the dynamical matrix
ph.read(acoustic=True)

FCM = ph.get_force_constant()

# High-symmetry points in the Brillouin zone
points = ibz_points['hexagonal']
G = points['Gamma']
K = points['K']
M = points['M']
A = points['A']

point_names = ['$\Gamma$', 'K', 'M', '$\Gamma$', 'A']
path = [G, K, M, G, A]

path_kc, q, Q = bandpath(path, atoms.cell, 100)
omega_kn = 1000 * ph.band_structure(path_kc)

properties = {"omega_kn" : omega_kn.tolist(), "point_names" : point_names, "path" : path, "Q" : list(Q), "q" : list(q), "FCM" : FCM.tolist() }
