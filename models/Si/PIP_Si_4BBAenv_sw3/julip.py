"""
Reduced ASE to JuLIP.jl interface

Requires `pyjulia` package from https://github.com/JuliaPy/pyjulia
"""

import numpy as np
from ase.calculators.calculator import Calculator

from julia import Julia
julia = Julia()
julia.using("JuLIP")
julia.using("ASE")

# Workaround limitiation in PyCall that does not allow types to be called
#   https://github.com/JuliaPy/PyCall.jl/issues/319

ASEAtoms = julia.eval('ASEAtoms(a) = ASE.ASEAtoms(a)')
ASECalculator = julia.eval('ASECalculator(c) = ASE.ASECalculator(c)')
fixedcell = julia.eval('fixedcell(a) = JuLIP.Constraints.FixedCell(a)')
variablecell = julia.eval('variablecell(a) = JuLIP.Constraints.VariableCell(a)')

class JulipCalculator(Calculator):
    """
    ASE-compatible Calculator that calls JuLIP.jl for forces and energy
    """
    implemented_properties = ['forces', 'energy', 'stress']
    default_parameters = {}
    name = 'JulipCalculator'

    def __init__(self, julip_calculator):
        Calculator.__init__(self)
        self.julip_calculator = julia.eval(julip_calculator)

    def calculate(self, atoms, properties, system_changes):
        Calculator.calculate(self, atoms, properties, system_changes)
        julia_atoms = ASEAtoms(atoms)
        self.results = {}
        if 'energy' in properties:
            self.results['energy'] = julia.energy(self.julip_calculator, julia_atoms)
        if 'forces' in properties:
            self.results['forces'] = np.array(julia.forces(self.julip_calculator, julia_atoms))
        if 'stress' in properties:
            self.results['stress'] = np.array(julia.stress(self.julip_calculator, julia_atoms))
