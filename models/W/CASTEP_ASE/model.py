from ase.calculators.castep import Castep
import os
from distutils import spawn
import ase

mpirun = spawn.find_executable('mpirun')
mpirun_args = "-n 16"
castep = "/opt/castep/18.1/bin/castep.mpi"
print mpirun
print mpirun_args
print castep

os.environ['CASTEP_COMMAND'] = '{0} {1} {2}'.format(mpirun, mpirun_args, castep)

name = 'CASTEP'
Castep.name = name
Castep.todict = lambda self: {}

no_checkpoint = True

def start(test_name):
        global calculator
        calculator = ase.calculators.castep.Castep(directory=test_name,
                         cut_off_energy=600,
                         opt_strategy='speed',
                         max_scf_cycles=250,
                         smearing_width='0.1',
                         kpoints_mp_spacing='0.015')
                                                
