from ase.build import bulk
from ase.calculators.emt import EMT
from ase.dft.kpoints import ibz_points, bandpath
from ase.phonons import Phonons
from ase.optimize import BFGS
from ase.lattice.cubic import BodyCenteredCubic
import ase
import quippy
#from ase.io import read

import model
# Setup crystal and EMT calculator
atoms = BodyCenteredCubic("Ti", latticeconstant=3.32)

#atoms.rattle(0.001)

#atoms.set_calculator(model.calculator)

#dyn = BFGS(atoms)
#dyn.run(fmax=0.00001)

atoms = ase.Atoms(quippy.Atoms("/Users/Cas/gits/gap-testing-framework/tests/Ti/phonon_bcc/bcc_k4_relaxed.xyz"))

# Phonon calculator
N = 3
ph = Phonons(atoms, model.calculator, supercell=(N, N, N), delta=0.01)
ph.run()

# Read forces and assemble the dynamical matrix
ph.read(acoustic=True)

# High-symmetry points in the Brillouin zone
points = ibz_points['bcc']
G = points['Gamma']
H = points['H']
P = points['P']
N = points['N']

point_names = ['$\Gamma$', 'H', 'P', '$\Gamma$', 'N']
path = [G, H, P, G, N]

path_kc, q, Q = bandpath(path, atoms.cell, 100)
omega_kn = 1000 * ph.band_structure(path_kc)

properties = {"omega_kn" : omega_kn.tolist(), "point_names" : point_names, "path" : path, "Q" : list(Q), "q" : list(q) }

