import model
import quippy
import numpy as np
from ase.lattice.cubic import Diamond

def force_calculator(calculator):
    lattice_displacements = [[0.001, 0, 0], [-0.001, 0, 0], [0, 0.001, 0], [0, -0.001, 0], [0,0,0.001], [0,0,-0.001]]
    data = []
    
    a0 = (20.0*8)**(1.0/3.0)
    
    num_atoms = len(Diamond(symbol='Si', latticeconstant=a0))
    
    for i in xrange(0,num_atoms):
        forces = []
        for disp in lattice_displacements:
            bulk = Diamond(symbol='Si', latticeconstant=a0)
            big_bulk = bulk * (2,2,2)
            big_bulk[i].position = big_bulk[i].position + disp
            big_bulk.set_calculator(calculator)
            force = big_bulk.get_forces()
            forces.append(force.tolist())
        data.append(forces)
    
    return data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7]


(config1, config2, config3, config4, config5, config6, config7, config8) = force_calculator(model.calculator)

properties = {"config1":config1, "config2":config2, "config3":config3, "config4":config4, "config4":config4, "config5":config5, "config6":config6, "config7":config7}
