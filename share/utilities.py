from ase.io import read, write
from ase.calculators.calculator import Calculator
#from ase.calculators.singlepoint import SinglePointCalculator
from ase.optimize import FIRE
from ase.optimize.precon import Exp, PreconLBFGS
from ase.constraints import UnitCellFilter, FixAtoms, voigt_6_to_full_3x3_stress, full_3x3_to_voigt_6_stress
from quippy.potential import Minim
import numpy as np
import os.path

import symmetrize

import os
import StringIO

from ase.io import read, write
from ase.calculators.singlepoint import SinglePointCalculator
from ase.optimize import FIRE
from ase.constraints import UnitCellFilter, voigt_6_to_full_3x3_stress, full_3x3_to_voigt_6_stress

try:
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
except:
    pass
from ase import Atoms
import numpy as np

#from quippy.io import AtomsWriter
#from quippy.cinoutput import CInOutput,OUTPUT
#from quippy.atoms import Atoms

def read_bulk_reference(model_name, test_name):
    log_file = 'model-{0}-test-{1}.txt'.format(model_name, test_name)
    log_lines = open(os.path.join('..', log_file), 'r').readlines()
    start_idx = log_lines.index('relaxed bulk\n')
    n_atoms = int(log_lines[start_idx+1])
    xyz_lines = log_lines[start_idx+1:start_idx+n_atoms+3]
    fh = StringIO.StringIO(''.join(xyz_lines))
    bulk = read(fh, format='extxyz')
    return bulk

def relax_atoms(atoms, tol=1e-3, method='lbfgs', max_steps=1000, traj_file=None, **kwargs):
    import model
    atoms.set_calculator(model.calculator)
    if hasattr(model, 'Optimizer'):
        method = 'model_optimizer'
        opt = model.Optimizer(atoms)
        opt.run(tol, max_steps)
    elif method.startswith('lbfgs') or method == 'fire' or method == 'cg_n':
        if method == 'lbfgs_ASE':
            from ase.optimize import LBFGS
            opt = LBFGS(atoms, **kwargs)
        elif method == 'cg_n':
            from quippy import Minim
            opt = Minim(atoms, relax_positions=True, relax_cell=False, method='cg_n')
        else:
            from ase.optimize.precon.precon import Exp
            from ase.optimize.precon.lbfgs import PreconLBFGS
            precon = None
            if method.endswith('precon'):
                precon = Exp(3.0, recalc_mu=True)
            if method.startswith('lbfgs'):
                opt = PreconLBFGS(atoms, precon=precon, **kwargs)
            else:
                opt = FIRE(atoms, **kwargs)
        if traj_file is not None and method != 'cg_n':
            traj = open(traj_file, 'w')
            def write_trajectory():
                write(traj, atoms, format='extxyz')
            opt.attach(write_trajectory)
        opt.run(tol, max_steps)
        try:
            traj.close()
        except:
            pass
    else:
        raise ValueError('unknown method %s!' % method)

    return atoms

import symmetrize
from ase.calculators.calculator import Calculator
class SymmetrizedCalculator(Calculator):
   implemented_properties = ['energy','forces','stress']
   def __init__(self, calc, atoms, *args, **kwargs):
      Calculator.__init__(self, *args, **kwargs)
      self.calc = calc
      (self.rotations, self.translations, self.symm_map) = symmetrize.prep(atoms, symprec=0.01)

   def calculate(self, atoms, properties, system_changes):
        Calculator.calculate(self, atoms, properties, system_changes)
        if system_changes:
            self.results = {}
        if 'energy' in properties and 'energy' not in self.results:
            self.results['energy'] = self.calc.get_potential_energy(atoms)
        if 'forces' in properties and 'forces' not in self.results:
            raw_forces = self.calc.get_forces(atoms)
            self.results['forces'] = symmetrize.forces(atoms.get_cell(), atoms.get_reciprocal_cell().T, raw_forces, 
                self.rotations, self.translations, self.symm_map)
        if 'stress' in properties and 'stress' not in self.results:
            raw_stress = self.calc.get_stress(atoms)
            if len(raw_stress) == 6: # Voigt
                raw_stress = voigt_6_to_full_3x3_stress(raw_stress)
            symmetrized_stress = symmetrize.stress(atoms.get_cell(), atoms.get_reciprocal_cell().T, raw_stress, self.rotations)
            self.results['stress'] = full_3x3_to_voigt_6_stress(symmetrized_stress)


def relax_atoms_cell(atoms, tol=1e-3, stol=None, method='lbfgs', max_steps=100, mask=None, traj_file=None,
                     hydrostatic_strain=False, constant_volume=False, precon_apply_positions=True,
                     precon_apply_cell=True, symmetrize = False, **kwargs):
    import model
    print "relax_atoms_cell using method",method
    if symmetrize:
        atoms.set_calculator(SymmetrizedCalculator(model.calculator, atoms))
    else:
        atoms.set_calculator(model.calculator)
    ## print "relax_atoms_cell initial e ", atoms.get_potential_energy()
    ## print "relax_atoms_cell initial f ", atoms.get_forces()
    ## print "relax_atoms_cell initial s ", atoms.get_stress()
    if hasattr(model, 'Optimizer'):
        method = 'model_optimizer'
    if method != 'cg_n':
        atoms = UnitCellFilter(atoms, mask=mask,
                               hydrostatic_strain=hydrostatic_strain,
                               constant_volume=constant_volume)
    if method.startswith('lbfgs') or method == 'fire' or method == 'cg_n':
        if method == 'cg_n':
            from quippy import Minim, fzeros
            atoms.info['Minim_Hydrostatic_Strain'] = hydrostatic_strain
            atoms.info['Minim_Constant_Volume'] = constant_volume
            if mask is not None:
                 atoms.info['Minim_Lattice_Fix'] = fzeros( (3,3) )
                 if not mask[0]:
                     atoms.info['Minim_Lattice_Fix'][1,1] = 1.0
                 if not mask[1]:
                     atoms.info['Minim_Lattice_Fix'][2,2] = 1.0
                 if not mask[2]:
                     atoms.info['Minim_Lattice_Fix'][3,3] = 1.0
                 if not mask[3]:
                     atoms.info['Minim_Lattice_Fix'][1,2] = 1.0
                     atoms.info['Minim_Lattice_Fix'][2,1] = 1.0
                 if not mask[4]:
                     atoms.info['Minim_Lattice_Fix'][2,3] = 1.0
                     atoms.info['Minim_Lattice_Fix'][3,2] = 1.0
                 if not mask[5]:
                     atoms.info['Minim_Lattice_Fix'][1,3] = 1.0
                     atoms.info['Minim_Lattice_Fix'][3,1] = 1.0
            opt = Minim(atoms, relax_positions=True, relax_cell=True, method='cg_n')
        else:
            from ase.optimize.precon.precon import Exp
            from ase.optimize.precon.lbfgs import PreconLBFGS
            precon = None
            if method.endswith('precon'):
                precon = Exp(3.0, apply_positions=precon_apply_positions,
                             apply_cell=precon_apply_cell, recalc_mu=True)
            if method.startswith('lbfgs'):
                opt = PreconLBFGS(atoms, precon=precon, **kwargs)
            else:
                opt = FIRE(atoms, **kwargs)
        if traj_file is not None:
            traj = open(traj_file, 'w')
            def write_trajectory():
                try:
                    #print atoms.atoms.info["map_shift"]
                    print atoms.atoms.__dict__
                    #print atoms.atoms.map_shift
                    #print atoms.atoms.n_neighb
                    #print atoms.atoms.orig_index
                    #print atoms.atoms.travel
                    write(traj, atoms.atoms, format='extxyz')
                except:
                    write(traj, atoms, format='extxyz')
            #opt.attach(write_trajectory)
        if method != 'cg_n' and isinstance(opt, PreconLBFGS):
            opt.run(tol, max_steps, smax=stol)
        else:
            opt.run(tol, max_steps)
        if traj_file is not None:
            traj.close()
    elif method == 'model_optimizer':
        opt = model.Optimizer(atoms.atoms)
        opt.run()
    else:
        raise ValueError('unknown method %s!' % method)

    if isinstance(atoms, UnitCellFilter):
        return atoms.atoms
    else:
        return atoms


# def relax_atoms(atoms, tol=1e-3, method='cg', max_steps=100, traj_file=None, **kwargs):
#     import model
#     if traj_file:
#        traj = CInOutput(traj_file, action=OUTPUT)
#     else:
#        traj=None
#     at=Atoms(atoms)
#     if hasattr(model.calculator, 'cutoff'):
#         at.set_cutoff(model.calculator.cutoff()+0.1)
#     model.calculator.minim(at, method=method, convergence_tol=tol, max_steps=max_steps,
#                            do_print=(traj is not None), print_cinoutput=traj, do_pos=True, do_lat=False)
#     atoms.set_positions(at.get_positions())
#     return atoms

# def relax_atoms_cell(atoms, tol=1e-3, method='cg', max_steps=100, mask=None, traj_file=None, **kwargs):
#     import model
#     at = Atoms(atoms)
#     if hasattr(model.calculator, 'cutoff'):
#         at.set_cutoff(model.calculator.cutoff()+0.1)
#     if traj_file is not None:
#         traj = CInOutput(traj_file, action=OUTPUT)
#     else:
#         traj=None
#     if mask is not None:
#         at.info['Minim_Lattice_Fix'] = quippy.fzeros( (3,3) )
#         if not mask[0]:
#             at.info['Minim_Lattice_Fix'][1,1] = 1.0
#         if not mask[1]:
#             at.info['Minim_Lattice_Fix'][2,2] = 1.0
#         if not mask[2]:
#             at.info['Minim_Lattice_Fix'][3,3] = 1.0
#         if not mask[3]:
#             at.info['Minim_Lattice_Fix'][1,2] = 1.0
#             at.info['Minim_Lattice_Fix'][2,1] = 1.0
#         if not mask[4]:
#             at.info['Minim_Lattice_Fix'][2,3] = 1.0
#             at.info['Minim_Lattice_Fix'][3,2] = 1.0
#         if not mask[5]:
#             at.info['Minim_Lattice_Fix'][1,3] = 1.0
#             at.info['Minim_Lattice_Fix'][3,1] = 1.0
#     model.calculator.minim(at, method=method, convergence_tol=tol, max_steps=max_steps,
#                            do_print=(traj is not None), do_lat=True, do_pos=True, print_cinoutput=traj)
#     atoms.set_calculator(model.calculator)
#     atoms.set_positions(at.get_positions())
#     atoms.set_cell(at.get_cell())
#     return atoms

def evaluate(atoms, do_energy=True, do_forces=True, do_stress=True):
    import model
    atoms.set_calculator(model.calculator)

    energy = None
    if do_energy:
        energy = atoms.get_potential_energy()

    forces = None
    if do_forces:
        forces = atoms.get_forces()

    stress = None
    if do_stress:
        stress = atoms.get_stress()

    spc = SinglePointCalculator(atoms,
                                energy=energy,
                                forces=forces,
                                stress=stress)
    atoms.set_calculator(spc)
    return atoms

def phonons(model,bulk,supercell,dx,mesh=None,points=None,n_points=50):
    
    import model
    
    unitcell = PhonopyAtoms(symbols=bulk.get_chemical_symbols(),
                            cell=bulk.get_cell(),
                            scaled_positions=bulk.get_scaled_positions())
    phonon = Phonopy(unitcell,supercell)
    phonon.generate_displacements(distance=dx)
    
    sets_of_forces = []
    
    for s in phonon.get_supercells_with_displacements():
        at = Atoms(cell=s.get_cell(),
                       symbols=s.get_chemical_symbols(),
                       scaled_positions=s.get_scaled_positions(),
                       pbc=3*[True])
        at.set_calculator(model.calculator)
        sets_of_forces.append(at.get_forces())

    phonon.set_forces(sets_of_forces=sets_of_forces)
    phonon.produce_force_constants()

    properties = {}
    
    if mesh is not None:
        phonon.set_mesh(mesh,is_gamma_center=True)
        qpoints, weights, frequencies, eigvecs = phonon.get_mesh()

        properties["frequencies"] = frequencies.tolist()
        properties["weights"] = weights.tolist()
    
    if points is not None:
        bands = []
        for i in range(len(points)-1):
            band = []
            for r in np.linspace(0,1,n_points):
                band.append(points[i]+(points[i+1]-points[i])*r)
            bands.append(band)
            
        phonon.set_band_structure(bands,is_eigenvectors=True,is_band_connection=False)
        band_q_points, band_distances, band_frequencies, band_eigvecs = phonon.get_band_structure()
        
        band_distance_max = np.max(band_distances)
        band_distances = [(_b / band_distance_max).tolist() for _b in band_distances]

        band_frequencies = [ _b.tolist() for _b in band_frequencies ]
        
        properties["band_q_points"] = band_q_points
        properties["band_distances"] = band_distances
        properties["band_frequencies"] = band_frequencies
        properties["band_eigvecs"] = band_eigvecs
    
    properties["phonopy"] = phonon

    return properties

def name_of_file(file):
    if '/' in file:
        name = os.path.basename(os.path.dirname(file))
    else:
        name = file
    return name

def model_test_root(u_model_name=None, u_test_name=None, base_model=False):
    if u_model_name is None:
        if base_model:
            u_model_name = base_model_name
        else:
            u_model_name = model_name
    if u_test_name is None:
        u_test_name = test_name
    if system_label != '':
        return '{0}-model-{1}-test-{2}'.format(system_label, u_model_name, u_test_name)
    else:
        return 'model-{0}-test-{1}'.format(model_name, test_name)

class SymmetrizedCalculator(Calculator):
   implemented_properties = ['energy','forces','stress']
   def __init__(self, calc, atoms, *args, **kwargs):
      Calculator.__init__(self, *args, **kwargs)
      self.calc = calc
      (self.rotations, self.translations, self.symm_map) = symmetrize.prep(atoms)

   def calculate(self, atoms, properties, system_changes):
        Calculator.calculate(self, atoms, properties, system_changes)
        #print self.calc.get_stress(atoms)
        if system_changes:
            self.results = {}
        if 'energy' in properties and 'energy' not in self.results:
            self.results['energy'] = self.calc.get_potential_energy(atoms)
        if 'forces' in properties and 'forces' not in self.results:
            raw_forces = self.calc.get_forces(atoms)
            print raw_forces
            self.results['forces'] = symmetrize.forces(atoms.get_cell(), atoms.get_reciprocal_cell().T, raw_forces, 
                self.rotations, self.translations, self.symm_map)
        if 'stress' in properties and 'stress' not in self.results:
            #raw_stress = voigt_6_to_full_3x3_stress(self.calc.get_stress(atoms))
            raw_stress = self.calc.get_stress(atoms)
            symmetrized_stress = symmetrize.stress(atoms.get_cell(), atoms.get_reciprocal_cell().T, raw_stress, self.rotations)
            self.results['stress'] = full_3x3_to_voigt_6_stress(symmetrized_stress)

def relax_config(atoms, relax_pos, relax_cell, tol=1e-3, method='fire', max_steps=200, traj_file=None, constant_volume=False,
    refine_symmetry_tol=None, keep_symmetry=False, strain_mask = None, config_label=None, from_base_model=False, save_config=False, **kwargs):

    # get from base model if requested
    import model
    if from_base_model:
        if config_label is None:
            raise ValueError('from_base_model is set but no config_label provided')
        try:
            base_run_file = os.path.join('..',base_run_root,base_run_root+'-'+config_label+'-relaxed.xyz')
            atoms_in = read(base_run_file, format='extxyz')
            # set positions from atoms_in rescaling to match current cell
            saved_cell = atoms.get_cell().copy()
            atoms.set_cell(atoms_in.get_cell())
            atoms.set_positions(atoms_in.get_positions())
            atoms.set_cell(saved_cell, scale_atoms=True)
            print "relax_config read config from ",base_run_file
        except:
            try:
                print "relax_config failed to read base run config from ",base_run_root+'-'+config_label+'-relaxed.xyz'
            except:
                print "relax_config failed to determined base_run_root"

    if refine_symmetry_tol is not None:
        symmetrize.refine(atoms, refine_symmetry_tol)
    if keep_symmetry:
        atoms.set_calculator(SymmetrizedCalculator(model.calculator, atoms))
    else:
        atoms.set_calculator(model.calculator)

    if method == 'lbfgs':
        if 'move_mask' in atoms.arrays:
            atoms.set_constraint(FixAtoms(np.where(atoms.arrays['move_mask'] == 0)[0]))
        if relax_cell:
            atoms = UnitCellFilter(atoms, mask=strain_mask, constant_volume=constant_volume)
        precon = Exp(3.0)
        opt = PreconLBFGS(atoms, precon=precon, **kwargs)
        if traj_file is not None:
            traj = open(traj_file, "w")
            def write_trajectory():
                write(traj, atoms, format='extxyz')
            opt.attach(write_trajectory)
    elif method == 'cg_n':
        if strain_mask is not None:
            raise(Exception("strain_mask not supported with method='cg_n'"))
        atoms.info['Minim_Constant_Volume'] = constant_volume
        opt = Minim(atoms, relax_positions=relax_pos, relax_cell=relax_cell, method='cg_n')
    else:
        raise ValueError('unknown method %s!' % method)

    opt.run(tol, max_steps)

    # in case we wrapped in a UnitCellFilter
    try:
        atoms = atoms.atoms
    except:
        pass

    symmetrize.check(atoms, refine_symmetry_tol)
    symmetrize.check(atoms, 1.0e-6)

    # in case we had a trajectory saved
    try:
        traj.close()
    except:
        pass

    if save_config:
        if config_label is None:
            raise ValueError('save_config is set but no config_label provided')
        write(run_root+'-'+config_label+'-relaxed.xyz', atoms, format='extxyz')

    if keep_symmetry:
        atoms.set_calculator(model.calculator)

    return atoms


def evaluate(atoms, do_energy=True, do_forces=True, do_stress=True):
    import model
    from ase.calculators.singlepoint import SinglePointCalculator
    atoms.set_calculator(model.calculator)

    energy = None
    if do_energy:
        energy = atoms.get_potential_energy()

    forces = None
    if do_forces:
        forces = atoms.get_forces()

    stress = None
    if do_stress:
        stress = atoms.get_stress()

    spc = SinglePointCalculator(atoms,
                                energy=energy,
                                forces=forces,
                                stress=stress)
    atoms.set_calculator(spc)
    return atoms
    
    import os
    import StringIO

    from ase.io import read, write
    from ase.calculators.singlepoint import SinglePointCalculator
    from ase.optimize import FIRE
    from ase.constraints import UnitCellFilter, voigt_6_to_full_3x3_stress, full_3x3_to_voigt_6_stress

    try:
        from phonopy import Phonopy
        from phonopy.structure.atoms import PhonopyAtoms
    except:
        pass
    from ase import Atoms
    import numpy as np

    #from quippy.io import AtomsWriter
    #from quippy.cinoutput import CInOutput,OUTPUT
    #from quippy.atoms import Atoms
