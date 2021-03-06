import model
import quippy
import numpy as np
#from ase.lattice.hexagonal import HexagonalClosedPacked

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

def force_calculator(calculator):
    lattice_displacements = [[0.001, 0, 0], [-0.001, 0, 0], [0, 0.001, 0], [0, -0.001, 0], [0,0,0.001], [0,0,-0.001]]
    data = []
    
    c_vs_a = 0.610
    a = 4.630
    
    num_atoms = len(create_omega_custom(c_vs_a, a, 22))
    
    for i in xrange(0,num_atoms):
        forces = []
        for disp in lattice_displacements:
            bulk = create_omega_custom(c_vs_a, a, 22)
            big_bulk = bulk * (2,2,2)
            big_bulk[i].position = big_bulk[i].position + disp
            big_bulk.set_calculator(calculator)
            force = big_bulk.get_forces()
            forces.append(force.tolist())
        data.append(forces)
    
    return data[0], data[1], data[2]


(config1, config2, config3) = force_calculator(model.calculator)

properties = {"config1":config1, "config2":config2, "config3":config3}
