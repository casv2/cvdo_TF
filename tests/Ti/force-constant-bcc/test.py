import model
import quippy
import numpy as np
from ase.lattice.cubic import BodyCenteredCubic

def force_calculator(calculator):
    lattice_displacements = [[0.001, 0, 0], [-0.001, 0, 0], [0, 0.001, 0], [0, -0.001, 0], [0,0,0.001], [0,0,-0.001]]
    data = []
    
    a0 = 3.32
    
    num_atoms = len(BodyCenteredCubic(symbol='Ti', latticeconstant=a0))
    
    for i in xrange(0,num_atoms):
        forces = []
        for disp in lattice_displacements:
            bulk = BodyCenteredCubic(symbol='Ti', latticeconstant=a0)
            big_bulk = bulk * (2,2,2)
            big_bulk[i].position = big_bulk[i].position + disp
            big_bulk.set_calculator(calculator)
            force = big_bulk.get_forces()
            forces.append(force.tolist())
        data.append(forces)
    
    return data[0], data[1]


(config1, config2) = force_calculator(model.calculator)

properties = {"config1":config1, "config2":config2}
