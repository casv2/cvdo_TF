from ase.build import bulk
from ase.calculators.emt import EMT
from ase.dft.kpoints import ibz_points, bandpath
from ase.phonons import Phonons
from ase.lattice.cubic import Diamond

import model

atoms = Diamond("Si", latticeconstant=5.42)

N = 3
ph = Phonons(atoms, model.calculator, supercell=(N, N, N), delta=0.05)
ph.run()

ph.read(acoustic=True)

