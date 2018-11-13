from ase.lattice.cubic import Diamond
from ase.md.langevin import Langevin
from ase import units
from ase.md.velocitydistribution import (MaxwellBoltzmannDistribution,
                                         Stationary, ZeroRotation)
# the current model
import model 
import sys
from cStringIO import StringIO
import numpy as np 
from ase.io.trajectory import Trajectory

#def MD(model, temp):
def MD():
	old_stdout = sys.stdout
	sys.stdout = mystdout = StringIO()
	

	orig_stdout = sys.stdout
	f = open('out.txt', 'w')
	sys.stdout = f

	#f = open('out_' + str(temp) + '.txt', 'w')
	#sys.stdout = f
	Tmin = 1000
	Tmax = 15000

	a0 = 5.43
	N = 1
	atoms = Diamond(symbol='Si', latticeconstant=a0)
	atoms *= (N,N,N)
	atoms.set_calculator(model.calculator)
	MaxwellBoltzmannDistribution(atoms, Tmin * units.kB) ###is this allowed?
	Stationary(atoms)  # zero linear momentum
	ZeroRotation(atoms)

	def printenergy(a=atoms):
		epot = a.get_potential_energy() / len(a)
		ekin = a.get_kinetic_energy() / len(a)
		print(ekin)

	traj = Trajectory('Si.traj', 'w', atoms)

 	for temp in np.linspace(Tmin,Tmax,5):
		dyn = Langevin(atoms, 5*units.fs, units.kB * temp, 0.002)
		dyn.attach(traj.write, interval=1)
		dyn.attach(printenergy, interval=1)
		dyn.run(100)

	sys.stdout = orig_stdout
	f.close()
	#sys.stdout = old_stdout
	
	#print "helo", mystdout.getvalue(), mystdout.getvalue().split("\n")[:-1]

	#ekin_list = mystdout.getvalue().split("\n")[:-1]
	#return ekin_list

properties = {'1000K': MD()}








#MD(model, 1000)