#!/usr/bin/env python3
"""
SOLVENT SELECTION SYSTEM FOR CHROMATOGRAPHY
============================================
Complete implementation of all numerical rules for solvent selection
with full references and scoring mechanisms.

Author: Based on comprehensive literature review
Version: 1.0.0
"""

import os
import sys
import math
import numpy as np
import pandas as pd
from datetime import datetime
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional, Any
import argparse
import warnings
warnings.filterwarnings('ignore')

# Try to import RDKit; if unavailable, continue with heuristic estimators
try:
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors
    USE_RDKIT = True
except Exception:
    USE_RDKIT = False

# ============================================================================
# PART 1: MOLECULAR DESCRIPTOR GENERATION
# ============================================================================

class MolecularDescriptorGenerator:
    """
    Generate all molecular descriptors required for solvent selection rules
    
    References:
    -----------
    Valkó, K. (2004). J. Chromatogr. A, 1037(1-2), 299-310.
    Abraham, M.H. (1993). Chemical Society Reviews, 22(2), 73-83.
    Kamlet, M.J., et al. (1983). J. Org. Chem., 48(17), 2877-2887.
    """
    
    def __init__(self):
        """Initialize descriptor generator with functional group patterns"""
        # Functional group patterns for estimation
        self.patterns = {
            'carboxylic_acid': '[CX3](=O)[OX2H1]',
            'phenol': '[c]1[c][c][c]([OH])[c][c]1',
            'primary_amine': '[NX3;H2]',
            'secondary_amine': '[NX3;H1]',
            'tertiary_amine': '[NX3;H0]',
            'ether': '[OD2]([#6])[#6]',
            'carbonyl': '[CX3]=[OX1]',
            'alcohol': '[#6][OX2H]',
            'phosphoric_acid': 'P(=O)(O)(O)O',
            'aromatic': 'a',
            'nitro': '[NX3](=O)=O',
            'sulfonic_acid': 'S(=O)(=O)O'
        }
        
    def generate_all_descriptors(self, smiles):
        """
        Generate complete descriptor set from SMILES.

        If RDKit is available, use RDKit to compute conservative, commonly-used
        descriptors (LogP, MW, HBD/HBA, TPSA, rotatable bonds). Other specialized
        parameters (Kamlet, Abraham) still use heuristic estimators.
        """
        descriptors = {}
        if USE_RDKIT:
            try:
                descriptors.update(self._rdkit_descriptors(smiles))
            except Exception:
                # Fall back to heuristic if RDKit fails for any reason
                descriptors.update(self._estimate_descriptors_from_smiles(smiles))
        else:
            descriptors.update(self._estimate_descriptors_from_smiles(smiles))

        # Ensure remaining descriptors (Kamlet, pKa, etc.) are present
        if 'hbd_alpha' not in descriptors:
            descriptors['hbd_alpha'] = self._estimate_hbd_alpha(self._analyze_smiles(smiles))
        if 'hba_beta' not in descriptors:
            descriptors['hba_beta'] = self._estimate_hba_beta(self._analyze_smiles(smiles))
        if 'kamlet_alpha' not in descriptors:
            descriptors['kamlet_alpha'] = self._estimate_kamlet_alpha(self._analyze_smiles(smiles))
        if 'kamlet_beta' not in descriptors:
            descriptors['kamlet_beta'] = self._estimate_kamlet_beta(self._analyze_smiles(smiles))
        if 'kamlet_pi_star' not in descriptors:
            descriptors['kamlet_pi_star'] = self._estimate_kamlet_pi_star(self._analyze_smiles(smiles))
        if 'pKa' not in descriptors:
            descriptors['pKa'] = self._estimate_pka(self._analyze_smiles(smiles))

        # LogD calculations based on available logP and pKa
        descriptors['logD_pH7'] = self._calculate_logD(descriptors['logP'], descriptors['pKa'], 7.0)
        descriptors['logD_pH3'] = self._calculate_logD(descriptors['logP'], descriptors['pKa'], 3.0)
        descriptors['logD_pH9'] = self._calculate_logD(descriptors['logP'], descriptors['pKa'], 9.0)

        return descriptors

    def _rdkit_descriptors(self, smiles):
        """Compute descriptors using RDKit."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError('RDKit failed to parse SMILES')
        desc = {}
        # LogP via Crippen
        try:
            desc['logP'] = float(Crippen.MolLogP(mol))
        except Exception:
            desc['logP'] = self._estimate_logp(smiles, self._analyze_smiles(smiles))
        # H-bond donors/acceptors
        try:
            desc['h_bond_donors'] = int(Descriptors.NumHDonors(mol))
            desc['h_bond_acceptors'] = int(Descriptors.NumHAcceptors(mol))
        except Exception:
            info = self._analyze_smiles(smiles)
            desc['h_bond_donors'] = info.get('hbd', 0)
            desc['h_bond_acceptors'] = info.get('hba', 0)
        # Map to Abraham-like parameters using simple scaling
        desc['hbd_alpha'] = min(1.5, desc['h_bond_donors'] * 0.4)
        desc['hba_beta'] = min(2.0, desc['h_bond_acceptors'] * 0.3)
        # Molecular weight, PSA, rotatable bonds
        try:
            desc['molecular_weight'] = float(Descriptors.MolWt(mol))
            desc['polar_surface_area'] = float(Descriptors.TPSA(mol))
            desc['rotatable_bonds'] = int(Descriptors.NumRotatableBonds(mol))
        except Exception:
            info = self._analyze_smiles(smiles)
            desc['molecular_weight'] = info.get('mw', 150)
            desc['polar_surface_area'] = info.get('psa', 50)
            desc['rotatable_bonds'] = info.get('rotatable', 3)
        # Aromatic rings heuristic
        desc['aromatic_rings'] = sum(1 for ring in mol.GetRingInfo().AtomRings() if any(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring))
        return desc
    
    def _estimate_descriptors_from_smiles(self, smiles):
        """
        Estimate molecular descriptors from SMILES using heuristic rules
        This is a simplified version for standalone operation
        """
        descriptors = {}
        
        # Basic molecular properties from SMILES analysis
        mol_info = self._analyze_smiles(smiles)
        
        # LogP estimation (simplified Crippen method approximation)
        descriptors['logP'] = self._estimate_logp(smiles, mol_info)
        
        # Abraham HBD parameter (Σα₂ᴴ)
        descriptors['hbd_alpha'] = self._estimate_hbd_alpha(mol_info)
        
        # Abraham HBA parameter (Σβ₂ᴴ)
        descriptors['hba_beta'] = self._estimate_hba_beta(mol_info)
        
        # Kamlet-Taft parameters
        descriptors['kamlet_alpha'] = self._estimate_kamlet_alpha(mol_info)
        descriptors['kamlet_beta'] = self._estimate_kamlet_beta(mol_info)
        descriptors['kamlet_pi_star'] = self._estimate_kamlet_pi_star(mol_info)
        
        # pKa estimation
        descriptors['pKa'] = self._estimate_pka(mol_info)
        
        # Additional descriptors
        descriptors['molecular_weight'] = mol_info.get('mw', 150)
        descriptors['polar_surface_area'] = mol_info.get('psa', 50)
        descriptors['rotatable_bonds'] = mol_info.get('rotatable', 3)
        descriptors['aromatic_rings'] = mol_info.get('aromatic_rings', 0)
        descriptors['h_bond_donors'] = mol_info.get('hbd', 1)
        descriptors['h_bond_acceptors'] = mol_info.get('hba', 2)
        
        # LogD at different pH
        descriptors['logD_pH7'] = self._calculate_logD(descriptors['logP'], descriptors['pKa'], 7.0)
        descriptors['logD_pH3'] = self._calculate_logD(descriptors['logP'], descriptors['pKa'], 3.0)
        descriptors['logD_pH9'] = self._calculate_logD(descriptors['logP'], descriptors['pKa'], 9.0)
        
        return descriptors
    
    def _analyze_smiles(self, smiles):
        """
        Analyze SMILES string for functional groups and properties
        Simplified parser for common functional groups
        """
        info = {
            'smiles': smiles,
            'length': len(smiles),
            'has_carboxylic': 1 if 'C(=O)O' in smiles or 'COOH' in smiles else 0,
            'has_phenol': 1 if 'c1ccccc1O' in smiles or 'Oc1ccccc1' in smiles else 0,
            'has_amine': 1 if 'N' in smiles and not 'N=' in smiles else 0,
            'has_primary_amine': 1 if 'N' in smiles and 'H2' in smiles else 0,
            'has_secondary_amine': 1 if 'N' in smiles and 'H' in smiles and not 'H2' in smiles else 0,
            'has_tertiary_amine': 1 if 'N' in smiles and '(' in smiles and ')' in smiles else 0,
            'has_ether': 1 if 'O' in smiles and not 'O=' in smiles and not 'O-' in smiles else 0,
            'has_carbonyl': 1 if 'C=O' in smiles or '=O' in smiles else 0,
            'has_alcohol': 1 if 'OH' in smiles else 0,
            'has_aromatic': 1 if 'c1' in smiles or 'c2' in smiles else 0,
            'has_nitro': 1 if 'N=O' in smiles or 'NO2' in smiles else 0,
            'has_sulfonic': 1 if 'S=O' in smiles or 'SO3' in smiles else 0,
            'has_phosphoric': 1 if 'P=O' in smiles or 'PO4' in smiles else 0,
            'atom_count': smiles.count('C') + smiles.count('N') + smiles.count('O') + 
                         smiles.count('S') + smiles.count('P') + smiles.count('F') +
                         smiles.count('Cl') + smiles.count('Br') + smiles.count('I'),
            'oxygen_count': smiles.count('O'),
            'nitrogen_count': smiles.count('N'),
            'ring_count': smiles.count('1') + smiles.count('2') + smiles.count('3') // 2
        }
        
        # Estimate molecular weight (rough approximation)
        atomic_masses = {'C': 12, 'N': 14, 'O': 16, 'S': 32, 'P': 31, 'F': 19, 'Cl': 35.5, 'Br': 80, 'I': 127}
        mw = 0
        for atom, mass in atomic_masses.items():
            count = smiles.count(atom)
            mw += count * mass
        info['mw'] = max(mw, 100)  # Minimum 100
        
        # Estimate HBD count
        info['hbd'] = (info['has_carboxylic'] + info['has_phenol'] + 
                      info['has_alcohol'] + info['has_primary_amine'])
        
        # Estimate HBA count
        info['hba'] = (info['has_carbonyl'] + info['has_ether'] + 
                      info['has_amine'] + info['has_nitro'])
        
        # Estimate aromatic rings
        info['aromatic_rings'] = info['ring_count'] if info['has_aromatic'] else 0
        
        # Estimate rotatable bonds
        info['rotatable'] = max(1, info['atom_count'] // 5)
        
        # Estimate PSA (rough approximation)
        info['psa'] = info['hba'] * 20 + info['hbd'] * 15
        
        return info
    
    def _estimate_logp(self, smiles, info):
        """
        Estimate LogP using fragment-based approach
        
        Reference: Sangster, J. (1997). Octanol-Water Partition Coefficients
        """
        base_logp = 0.0
        
        # Carbon chain contribution
        carbon_chain = smiles.count('C') - info['has_carboxylic'] * 2
        base_logp += carbon_chain * 0.5
        
        # Functional group contributions
        if info['has_carboxylic']:
            base_logp -= 1.5  # Carboxylic acid reduces LogP
        if info['has_phenol']:
            base_logp += 1.5  # Phenol increases LogP
        if info['has_amine']:
            base_logp -= 1.0  # Amine reduces LogP
        if info['has_alcohol']:
            base_logp -= 1.2  # Alcohol reduces LogP
        if info['has_ether']:
            base_logp += 0.5  # Ether increases LogP
        if info['has_carbonyl']:
            base_logp -= 0.8  # Carbonyl reduces LogP
        if info['has_aromatic']:
            base_logp += 1.5  # Aromatic rings increase LogP
        if info['has_nitro']:
            base_logp -= 0.5  # Nitro reduces LogP
        
        # Halogen contributions
        base_logp += smiles.count('F') * 0.2
        base_logp += smiles.count('Cl') * 0.7
        base_logp += smiles.count('Br') * 1.0
        base_logp += smiles.count('I') * 1.5
        
        return max(-3, min(12, base_logp))
    
    def _estimate_hbd_alpha(self, info):
        """
        Estimate Abraham HBD parameter (Σα₂ᴴ)
        
        Reference: Abraham, M.H. (1993). Chemical Society Reviews, 22(2), 73-83.
        """
        alpha = info['hbd'] * 0.4
        
        if info['has_carboxylic']:
            alpha += 0.6
        if info['has_phenol']:
            alpha += 0.3
        
        return min(alpha, 1.5)
    
    def _estimate_hba_beta(self, info):
        """
        Estimate Abraham HBA parameter (Σβ₂ᴴ)
        
        Reference: Abraham, M.H., et al. (1990). J. Chromatogr., 518, 329-348.
        """
        beta = info['hba'] * 0.3
        
        if info['has_amine']:
            beta += 0.4
        if info['has_ether']:
            beta += 0.2
        if info['has_carbonyl']:
            beta += 0.3
        
        return min(beta, 2.0)
    
    def _estimate_kamlet_alpha(self, info):
        """
        Estimate Kamlet-Taft α (HBD acidity)
        
        Reference: Kamlet, M.J., et al. (1983). J. Org. Chem., 48(17), 2877-2887.
        """
        alpha = 0.0
        if info['hbd'] > 0:
            alpha = 0.3 * info['hbd']
            if info['has_carboxylic']:
                alpha = max(alpha, 0.8)
            elif info['has_alcohol']:
                alpha = max(alpha, 0.6)
        return min(alpha, 1.2)
    
    def _estimate_kamlet_beta(self, info):
        """
        Estimate Kamlet-Taft β (HBA basicity)
        
        Reference: Kamlet, M.J., et al. (1983). J. Org. Chem., 48(17), 2877-2887.
        """
        beta = 0.2 * info['hba']
        if info['has_amine']:
            beta += 0.3
        if info['has_carbonyl']:
            beta += 0.2
        return min(beta, 1.0)
    
    def _estimate_kamlet_pi_star(self, info):
        """
        Estimate Kamlet-Taft π* (polarizability)
        
        Reference: Kamlet, M.J., et al. (1983). J. Org. Chem., 48(17), 2877-2887.
        """
        pi_star = 0.1 + 0.1 * (info['mw'] / 100)
        if info['has_aromatic']:
            pi_star += 0.2
        if info['has_nitro']:
            pi_star += 0.1
        if info['has_carbonyl']:
            pi_star += 0.1
        return min(pi_star, 1.2)
    
    def _estimate_pka(self, info):
        """
        Estimate pKa using functional group contributions
        
        Reference: Rondinini, S., et al. (1987). Pure Appl. Chem., 59(11), 1549-1560.
        """
        if info['has_carboxylic']:
            return 4.5
        elif info['has_phenol']:
            return 10.0
        elif info['has_primary_amine']:
            return 9.5
        elif info['has_secondary_amine']:
            return 10.5
        elif info['has_tertiary_amine']:
            return 9.0
        elif info['has_phosphoric']:
            return 2.5
        elif info['has_sulfonic']:
            return 1.0
        else:
            return 7.0
    
    def _calculate_logD(self, logP, pKa, pH):
        """
        Calculate distribution coefficient at given pH
        
        Reference: Avdeef, A. (2003). Absorption and Drug Development
        """
        if pKa < 14:  # Ionizable compound
            if pH < pKa:  # Acidic conditions
                logD = logP - np.log10(1 + 10**(pH - pKa))
            else:  # Basic conditions
                logD = logP - np.log10(1 + 10**(pKa - pH))
        else:
            logD = logP
        return logD


# ============================================================================
# PART 2: SOLVENT PROPERTY DATABASE
# ============================================================================

@dataclass
class SolventProperties:
    """Complete solvent property dataset with references"""
    name: str
    logP: float
    s_parameter: float
    polarity_index: float
    hbd_alpha: float
    hba_beta: float
    kamlet_alpha: float
    kamlet_beta: float
    kamlet_pi_star: float
    pKa: Optional[float]
    viscosity_20c: float
    viscosity_30c: float
    viscosity_40c: float
    uv_cutoff: int
    selectivity_group: str
    xe: float
    xd: float
    xn: float
    reference_logP: str
    reference_snyder: str
    reference_abraham: str
    reference_kamlet: str
    reference_physical: str

class SolventDatabase:
    """
    Comprehensive solvent database with all parameters from original sources
    
    References:
    -----------
    Snyder, L.R. (1974). J. Chromatogr. A, 92(2), 223-230.
    Abraham, M.H., et al. (1987-1990). Various publications.
    Kamlet, M.J., et al. (1983). J. Org. Chem., 48(17), 2877-2887.
    Riddick, J.A., et al. (1986). Organic Solvents, 4th Ed.
    """
    
    def __init__(self):
        self.solvents = self._initialize_solvents()
        self.solvent_groups = self._group_solvents()
        
    def _initialize_solvents(self):
        """Initialize with data from primary sources"""
        
        solvents = {
            'Water': SolventProperties(
                name='Water',
                logP=-1.38,
                s_parameter=0.0,
                polarity_index=10.2,
                hbd_alpha=0.82,
                hba_beta=0.35,
                kamlet_alpha=1.17,
                kamlet_beta=0.47,
                kamlet_pi_star=1.09,
                pKa=14.00,
                viscosity_20c=1.002,
                viscosity_30c=0.798,
                viscosity_40c=0.653,
                uv_cutoff=190,
                selectivity_group='VIII',
                xe=0.37,
                xd=0.37,
                xn=0.25,
                reference_logP='Sangster, J. (1997) Table 4.2',
                reference_snyder='Snyder, L.R., et al. (1979) Anal. Chem., 51(8), 1173-1180',
                reference_abraham='Abraham, M.H. (1993) J. Phys. Org. Chem., 6, 660-684',
                reference_kamlet='Kamlet, M.J., et al. (1983) J. Org. Chem., 48(17), 2877-2887',
                reference_physical='Li, J., & Carr, P.W. (1997) Anal. Chem., 69(13), 2530-2536'
            ),
            
            'Methanol': SolventProperties(
                name='Methanol',
                logP=-0.77,
                s_parameter=2.6,
                polarity_index=5.1,
                hbd_alpha=0.37,
                hba_beta=0.48,
                kamlet_alpha=0.93,
                kamlet_beta=0.62,
                kamlet_pi_star=0.60,
                pKa=16.7,
                viscosity_20c=0.584,
                viscosity_30c=0.510,
                viscosity_40c=0.450,
                uv_cutoff=205,
                selectivity_group='II',
                xe=0.48,
                xd=0.22,
                xn=0.31,
                reference_logP='Sangster, J. (1997) Table 4.2',
                reference_snyder='Schoenmakers, P.J., et al. (1978) J. Chromatogr., 149, 519-537',
                reference_abraham='Abraham, M.H., et al. (1987) J. Chem. Soc., Perkin Trans. 2, 141-144',
                reference_kamlet='Kamlet, M.J., et al. (1983) J. Org. Chem., 48(17), 2877-2887',
                reference_physical='Li, J., & Carr, P.W. (1997) Anal. Chem., 69(13), 2530-2536'
            ),
            
            'Acetonitrile': SolventProperties(
                name='Acetonitrile',
                logP=-0.34,
                s_parameter=3.1,
                polarity_index=5.8,
                hbd_alpha=0.09,
                hba_beta=0.32,
                kamlet_alpha=0.19,
                kamlet_beta=0.40,
                kamlet_pi_star=0.75,
                pKa=None,
                viscosity_20c=0.369,
                viscosity_30c=0.325,
                viscosity_40c=0.290,
                uv_cutoff=190,
                selectivity_group='VIb',
                xe=0.31,
                xd=0.27,
                xn=0.42,
                reference_logP='Sangster, J. (1997) Table 4.2',
                reference_snyder='Schoenmakers, P.J., et al. (1978) J. Chromatogr., 149, 519-537',
                reference_abraham='Abraham, M.H., et al. (1990) J. Chromatogr., 518, 329-348',
                reference_kamlet='Kamlet, M.J., et al. (1983) J. Org. Chem., 48(17), 2877-2887',
                reference_physical='Li, J., & Carr, P.W. (1997) Anal. Chem., 69(13), 2530-2536'
            ),
            
            'Ethanol': SolventProperties(
                name='Ethanol',
                logP=-0.31,
                s_parameter=3.6,
                polarity_index=4.3,
                hbd_alpha=0.33,
                hba_beta=0.56,
                kamlet_alpha=0.83,
                kamlet_beta=0.77,
                kamlet_pi_star=0.54,
                pKa=19.1,
                viscosity_20c=1.200,
                viscosity_30c=0.990,
                viscosity_40c=0.830,
                uv_cutoff=210,
                selectivity_group='II',
                xe=0.52,
                xd=0.19,
                xn=0.29,
                reference_logP='Sangster, J. (1997) Table 4.2',
                reference_snyder='Snyder, L.R., et al. (1997) Practical HPLC Method Development, 2nd Ed., p. 235',
                reference_abraham='Abraham, M.H., et al. (1987) J. Chem. Soc., Perkin Trans. 2, 141-144',
                reference_kamlet='Kamlet, M.J., et al. (1983) J. Org. Chem., 48(17), 2877-2887',
                reference_physical='Li, J., & Carr, P.W. (1997) Anal. Chem., 69(13), 2530-2536'
            ),
            
            'Isopropanol': SolventProperties(
                name='Isopropanol',
                logP=0.05,
                s_parameter=4.2,
                polarity_index=3.9,
                hbd_alpha=0.33,
                hba_beta=0.56,
                kamlet_alpha=0.76,
                kamlet_beta=0.84,
                kamlet_pi_star=0.48,
                pKa=20.8,
                viscosity_20c=2.370,
                viscosity_30c=1.770,
                viscosity_40c=1.400,
                uv_cutoff=210,
                selectivity_group='II',
                xe=0.52,
                xd=0.19,
                xn=0.29,
                reference_logP='Sangster, J. (1997) Table 4.2',
                reference_snyder='Snyder, L.R., et al. (1997) Practical HPLC Method Development, 2nd Ed., p. 237',
                reference_abraham='Abraham, M.H., et al. (1987) J. Chem. Soc., Perkin Trans. 2, 141-144 (estimated)',
                reference_kamlet='Kamlet, M.J., et al. (1983) J. Org. Chem., 48(17), 2877-2887',
                reference_physical='Li, J., & Carr, P.W. (1997) Anal. Chem., 69(13), 2530-2536'
            ),
            
            'THF': SolventProperties(
                name='THF',
                logP=0.46,
                s_parameter=4.5,
                polarity_index=4.0,
                hbd_alpha=0.00,
                hba_beta=0.48,
                kamlet_alpha=0.00,
                kamlet_beta=0.55,
                kamlet_pi_star=0.58,
                pKa=None,
                viscosity_20c=0.550,
                viscosity_30c=0.460,
                viscosity_40c=0.400,
                uv_cutoff=212,
                selectivity_group='III',
                xe=0.38,
                xd=0.20,
                xn=0.42,
                reference_logP='Sangster, J. (1997) Table 4.2',
                reference_snyder='Snyder, L.R., et al. (1979) Anal. Chem., 51(8), 1173-1180',
                reference_abraham='Abraham, M.H., et al. (1990) J. Chromatogr., 518, 329-348',
                reference_kamlet='Kamlet, M.J., et al. (1983) J. Org. Chem., 48(17), 2877-2887',
                reference_physical='Li, J., & Carr, P.W. (1997) Anal. Chem., 69(13), 2530-2536'
            ),
            
            'Acetone': SolventProperties(
                name='Acetone',
                logP=-0.24,
                s_parameter=3.4,
                polarity_index=5.1,
                hbd_alpha=0.08,
                hba_beta=0.48,
                kamlet_alpha=0.08,
                kamlet_beta=0.43,
                kamlet_pi_star=0.71,
                pKa=None,
                viscosity_20c=0.324,
                viscosity_30c=0.295,
                viscosity_40c=0.270,
                uv_cutoff=330,
                selectivity_group='VIa',
                xe=0.35,
                xd=0.23,
                xn=0.42,
                reference_logP='Sangster, J. (1997) Table 4.2',
                reference_snyder='Snyder, L.R., et al. (1997) Practical HPLC Method Development, 2nd Ed., p. 236',
                reference_abraham='Estimated from Kamlet parameters',
                reference_kamlet='Kamlet, M.J., et al. (1983) J. Org. Chem., 48(17), 2877-2887',
                reference_physical='Riddick, J.A., et al. (1986) Organic Solvents'
            ),
            
            'Dichloromethane': SolventProperties(
                name='Dichloromethane',
                logP=1.25,
                s_parameter=3.1,
                polarity_index=3.1,
                hbd_alpha=0.10,
                hba_beta=0.05,
                kamlet_alpha=0.13,
                kamlet_beta=0.10,
                kamlet_pi_star=0.82,
                pKa=None,
                viscosity_20c=0.440,
                viscosity_30c=0.400,
                viscosity_40c=0.365,
                uv_cutoff=233,
                selectivity_group='V',
                xe=0.29,
                xd=0.18,
                xn=0.53,
                reference_logP='Sangster, J. (1997) Table 4.2',
                reference_snyder='Snyder, L.R. (1978) J. Chromatogr. Sci., 16(6), 223-234',
                reference_abraham='Abraham, M.H., et al. (1990) J. Chromatogr., 518, 329-348',
                reference_kamlet='Kamlet, M.J., et al. (1983) J. Org. Chem., 48(17), 2877-2887',
                reference_physical='Riddick, J.A., et al. (1986) Organic Solvents'
            ),
            
            'Hexane': SolventProperties(
                name='Hexane',
                logP=3.90,
                s_parameter=0.1,
                polarity_index=0.1,
                hbd_alpha=0.00,
                hba_beta=0.00,
                kamlet_alpha=0.00,
                kamlet_beta=0.00,
                kamlet_pi_star=-0.08,
                pKa=None,
                viscosity_20c=0.326,
                viscosity_30c=0.300,
                viscosity_40c=0.278,
                uv_cutoff=195,
                selectivity_group='I',
                xe=0.0,
                xd=0.0,
                xn=0.0,
                reference_logP='Sangster, J. (1997) Table 4.2',
                reference_snyder='Snyder, L.R. (1974) J. Chromatogr. A, 92(2), 223-230',
                reference_abraham='Abraham, M.H., et al. (1990) J. Chromatogr., 518, 329-348',
                reference_kamlet='Kamlet, M.J., et al. (1983) J. Org. Chem., 48(17), 2877-2887',
                reference_physical='Riddick, J.A., et al. (1986) Organic Solvents'
            )
        }
        
        return solvents
    
    def _group_solvents(self):
        """Group solvents by selectivity type for optimization"""
        groups = {
            'Non-polar': ['Hexane'],
            'Proton acceptor': ['Methanol', 'Ethanol', 'Isopropanol'],
            'Proton donor': ['THF'],
            'Dipole interaction': ['Acetonitrile', 'Acetone'],
            'Strong dipole': ['Dichloromethane'],
            'Strong H-bonding': ['Water']
        }
        return groups
    
    def get_solvent(self, name):
        """Get solvent properties by name"""
        return self.solvents.get(name)
    
    def get_all_solvents(self):
        """Get all solvents in database"""
        return list(self.solvents.values())
    
    def get_solvents_by_group(self, group):
        """Get solvents in specific selectivity group"""
        return [self.solvents[name] for name in self.solvent_groups.get(group, [])]
    
    def to_dataframe(self):
        """Convert database to pandas DataFrame"""
        data = []
        for name, solvent in self.solvents.items():
            data.append({
                'Solvent': solvent.name,
                'LogP': solvent.logP,
                'S_Parameter': solvent.s_parameter,
                'Polarity_Index': solvent.polarity_index,
                'HBD_Alpha': solvent.hbd_alpha,
                'HBA_Beta': solvent.hba_beta,
                'Kamlet_Alpha': solvent.kamlet_alpha,
                'Kamlet_Beta': solvent.kamlet_beta,
                'Kamlet_PiStar': solvent.kamlet_pi_star,
                'pKa': solvent.pKa,
                'Viscosity_20C': solvent.viscosity_20c,
                'Viscosity_30C': solvent.viscosity_30c,
                'Viscosity_40C': solvent.viscosity_40c,
                'UV_Cutoff': solvent.uv_cutoff,
                'Selectivity_Group': solvent.selectivity_group
            })
        return pd.DataFrame(data)


# ============================================================================
# PART 3: SCORING BREAKDOWN DATA CLASS
# ============================================================================

@dataclass
class ScoreBreakdown:
    """Detailed scoring breakdown for each rule"""
    rule_name: str
    score: float
    weight: float
    details: Dict
    reference: str


# ============================================================================
# PART 4: SOLVENT SCORING ENGINE
# ============================================================================

class SolventScoringEngine:
    """
    Comprehensive solvent scoring engine implementing all numerical rules
    
    References:
    -----------
    Valkó, K. (2004). J. Chromatogr. A, 1037(1-2), 299-310.
    Vitha, M., & Carr, P.W. (2006). J. Chromatogr. A, 1126(1-2), 143-194.
    Snyder, L.R., et al. (1997). Practical HPLC Method Development.
    Carr, P.W. (1993). Microchemical Journal, 48(1), 4-28.
    """
    
    def __init__(self, solvent_database):
        self.solvents = solvent_database
        self.rule_weights = self._initialize_weights()
        
    def _initialize_weights(self):
        """
        Initialize weights for each rule based on literature
        """
        return {
            'logP_matching': 0.25,
            'hbd_hba_complementarity': 0.20,
            'polarity_matching': 0.15,
            'kamlet_distance': 0.15,
            'viscosity': 0.10,
            'uv_transparency': 0.10,
            'ph_stability': 0.05
        }
    
    def calculate_logP_score(self, analyte_logP, solvent_logP, pH=None):
        """
        Calculate Log P matching score
        
        Rule: ΔLog P = Log P(analyte) - Log P(solvent)
        Optimal Range: 2.0 ≤ ΔLog P ≤ 3.5
        
        Reference: Valkó, K. (2004). J. Chromatogr. A, 1037(1-2), 299-310.
        """
        delta_logP = analyte_logP - solvent_logP
        
        optimal_min = 2.0
        optimal_max = 3.5
        ideal = 2.75
        
        if optimal_min <= delta_logP <= optimal_max:
            score = 100 * np.exp(-0.5 * ((delta_logP - ideal) / 0.75) ** 2)
        else:
            if delta_logP < optimal_min:
                deviation = (optimal_min - delta_logP) / 2.0
                score = 100 * max(0, 1 - deviation)
            else:
                deviation = (delta_logP - optimal_max) / 2.0
                score = 100 * max(0, 1 - deviation)
        
        if pH is not None:
            if pH < 3 and delta_logP < 2.0:
                score *= 0.8
            elif pH > 8 and delta_logP > 3.5:
                score *= 0.8
        
        details = {
            'delta_logP': delta_logP,
            'optimal_range': f'{optimal_min}-{optimal_max}',
            'ideal_value': ideal,
            'score': score
        }
        
        return ScoreBreakdown(
            rule_name='Log P Matching Rule',
            score=score,
            weight=self.rule_weights['logP_matching'],
            details=details,
            reference='Valkó, K. (2004). J. Chromatogr. A, 1037(1-2), 299-310.'
        )
    
    def calculate_hbd_hba_score(self, analyte_alpha, analyte_beta, solvent):
        """
        Calculate HBD/HBA complementarity score
        
        Rule: C_score = |α_analyte - β_solvent| + |β_analyte - α_solvent|
        Acceptable: C_score < 0.8, Optimal: C_score < 0.3
        
        References:
        - Abraham, M.H. (1993). Chemical Society Reviews, 22(2), 73-83.
        - Vitha, M., & Carr, P.W. (2006). J. Chromatogr. A, 1126(1-2), 143-194.
        """
        c_score = abs(analyte_alpha - solvent.hba_beta) + abs(analyte_beta - solvent.hbd_alpha)
        
        if c_score < 0.3:
            score = 100
        elif c_score < 0.5:
            score = 90
        elif c_score < 0.8:
            score = 75
        elif c_score < 1.2:
            score = 50
        else:
            score = 25
        
        if analyte_alpha > 0.5 and solvent.hba_beta < 0.2:
            score *= 0.9
        if analyte_beta > 0.5 and solvent.hbd_alpha < 0.1:
            score *= 0.9
            
        details = {
            'c_score': c_score,
            'analyte_alpha': analyte_alpha,
            'analyte_beta': analyte_beta,
            'solvent_alpha': solvent.hbd_alpha,
            'solvent_beta': solvent.hba_beta,
            'peak_shape_prediction': 'Excellent' if c_score < 0.3 else 
                                    'Good' if c_score < 0.5 else
                                    'Acceptable' if c_score < 0.8 else
                                    'Poor - tailing expected'
        }
        
        return ScoreBreakdown(
            rule_name='HBD/HBA Complementarity',
            score=score,
            weight=self.rule_weights['hbd_hba_complementarity'],
            details=details,
            reference='Vitha, M., & Carr, P.W. (2006). J. Chromatogr. A, 1126(1-2), 143-194.'
        )
    
    def calculate_polarity_score(self, analyte_polarity_index, solvent):
        """
        Calculate Snyder polarity matching score
        
        Rule: |P'_analyte - P'_solvent| ≤ 2.0 for acceptable solubility
        
        References:
        - Snyder, L.R. (1974). J. Chromatogr. A, 92(2), 223-230.
        - Snyder, L.R., et al. (1997). Practical HPLC Method Development, 2nd Ed.
        """
        if analyte_polarity_index is None:
            analyte_polarity_index = 5.0  # Default middle value
        
        delta_polarity = abs(analyte_polarity_index - solvent.polarity_index)
        
        if delta_polarity <= 2.0:
            score = 100 * (1 - delta_polarity / 4.0)
        else:
            excess = delta_polarity - 2.0
            score = 50 * max(0, 1 - excess / 3.0)
        
        details = {
            'delta_polarity': delta_polarity,
            'analyte_polarity': analyte_polarity_index,
            'solvent_polarity': solvent.polarity_index,
            'selectivity_group': solvent.selectivity_group
        }
        
        return ScoreBreakdown(
            rule_name='Snyder Polarity Matching',
            score=score,
            weight=self.rule_weights['polarity_matching'],
            details=details,
            reference='Snyder, L.R., et al. (1997). Practical HPLC Method Development, 2nd Ed., p. 235.'
        )
    
    def calculate_kamlet_distance(self, analyte_alpha, analyte_beta, analyte_pi_star, solvent):
        """
        Calculate Kamlet-Taft solvatochromic distance
        
        Rule: D = √[(α_a - α_s)² + (β_a - β_s)² + (π*_a - π*_s)²]
        Interpretation: D < 0.5 excellent, D ≥ 1.2 poor
        
        References:
        - Kamlet, M.J., et al. (1983). J. Org. Chem., 48(17), 2877-2887.
        - Carr, P.W. (1993). Microchemical Journal, 48(1), 4-28.
        """
        d_alpha = analyte_alpha - solvent.kamlet_alpha
        d_beta = analyte_beta - solvent.kamlet_beta
        d_pi = analyte_pi_star - solvent.kamlet_pi_star
        
        d_score = np.sqrt(d_alpha**2 + d_beta**2 + d_pi**2)
        
        if d_score < 0.5:
            score = 100
        elif d_score < 0.8:
            score = 85
        elif d_score < 1.2:
            score = 65
        else:
            score = 40
        
        if d_score > 1.5:
            score *= 0.8
            
        details = {
            'd_score': d_score,
            'd_alpha': d_alpha,
            'd_beta': d_beta,
            'd_pi': d_pi,
            'interpretation': 'Excellent' if d_score < 0.5 else
                             'Good' if d_score < 0.8 else
                             'Moderate' if d_score < 1.2 else
                             'Poor - mixture recommended'
        }
        
        return ScoreBreakdown(
            rule_name='Kamlet-Taft Solvatochromic Distance',
            score=score,
            weight=self.rule_weights['kamlet_distance'],
            details=details,
            reference='Carr, P.W. (1993). Microchemical Journal, 48(1), 4-28.'
        )
    
    def calculate_viscosity_score(self, solvent, temperature=30, flow_rate=1.0, column_length=150):
        """
        Calculate viscosity-based pressure score
        
        Reference: Li, J., & Carr, P.W. (1997). Anal. Chem., 69(13), 2530-2536.
        """
        if temperature == 20:
            viscosity = solvent.viscosity_20c
        elif temperature == 30:
            viscosity = solvent.viscosity_30c
        elif temperature == 40:
            viscosity = solvent.viscosity_40c
        else:
            viscosities = {20: solvent.viscosity_20c, 
                          30: solvent.viscosity_30c, 
                          40: solvent.viscosity_40c}
            viscosity = np.interp(temperature, [20, 30, 40], 
                                 [viscosities[20], viscosities[30], viscosities[40]])
        
        rpf = viscosity / 0.325
        phi = 800
        dp = 5
        pressure = (phi * viscosity * column_length * flow_rate) / (dp**2 * 1000)
        
        if pressure < 200:
            score = 100
        elif pressure < 300:
            score = 85
        elif pressure < 400:
            score = 65
        else:
            score = 40
        
        if temperature > 30 and pressure > 300:
            score *= 1.1
        
        details = {
            'viscosity': viscosity,
            'relative_pressure_factor': rpf,
            'estimated_pressure_bar': pressure,
            'temperature': temperature
        }
        
        return ScoreBreakdown(
            rule_name='Viscosity/Pressure Consideration',
            score=score,
            weight=self.rule_weights['viscosity'],
            details=details,
            reference='Li, J., & Carr, P.W. (1997). Anal. Chem., 69(13), 2530-2536.'
        )
    
    def calculate_uv_score(self, detection_wavelength, solvent):
        """
        Calculate UV transparency score
        
        Rule: Detection wavelength should be > (solvent cutoff + 20 nm)
        
        References:
        - Riddick, J.A., et al. (1986). Organic Solvents, 4th Ed.
        - Dolan, J.W. (1999). LCGC North America, 17(6), 518-522.
        """
        cutoff = solvent.uv_cutoff
        safe_wavelength = cutoff + 20
        
        if detection_wavelength < cutoff:
            score = 0
            details_note = 'Below cutoff'
        elif detection_wavelength < safe_wavelength:
            score = 50 * (detection_wavelength - cutoff) / 20
            details_note = "Risk of high baseline noise - use HPLC-grade solvent"
        else:
            score = 100
            details_note = "Good transparency"
        
        if detection_wavelength < 210 and solvent.name in ['Methanol', 'Acetonitrile']:
            score *= 0.9
            
        details = {
            'uv_cutoff': cutoff,
            'safe_wavelength': safe_wavelength,
            'detection_wavelength': detection_wavelength,
            'note': details_note
        }
        
        return ScoreBreakdown(
            rule_name='UV Transparency',
            score=score,
            weight=self.rule_weights['uv_transparency'],
            details=details,
            reference='Dolan, J.W. (1999). LCGC North America, 17(6), 518-522.'
        )
    
    def calculate_ph_stability_score(self, solvent, ph, column_type='silica'):
        """
        Calculate pH stability score
        
        References:
        - Neue, U.D. (1997). HPLC Columns, Chapter 5.
        - Majors, R.E. (2003). LCGC North America, 21(10), 962-972.
        """
        if column_type == 'silica':
            if solvent.name in ['Acetonitrile', 'Methanol', 'Ethanol']:
                if 2 <= ph <= 8:
                    score = 100
                elif ph < 2:
                    score = 70 - 10 * (2 - ph)
                else:
                    max_hours = 100 / (10**(ph - 8))
                    if max_hours > 100:
                        score = 90
                    elif max_hours > 50:
                        score = 70
                    elif max_hours > 10:
                        score = 50
                    else:
                        score = 30
            else:
                score = 50
        elif column_type == 'hybrid':
            if 1 <= ph <= 12:
                score = 100
            else:
                score = 60 - 10 * abs(ph - 12) if ph > 12 else 60 - 10 * abs(1 - ph)
        else:
            score = 70
        
        if ph < 3 and 'acetate' in solvent.name.lower():
            score *= 0.5
            
        details = {
            'ph': ph,
            'column_type': column_type,
            'ph_range': '2-8' if column_type == 'silica' else '1-12',
            'hydrolysis_risk': 'Yes' if ph < 3 and 'acetate' in solvent.name.lower() else 'No'
        }
        
        return ScoreBreakdown(
            rule_name='pH Stability',
            score=score,
            weight=self.rule_weights['ph_stability'],
            details=details,
            reference='Neue, U.D. (1997). HPLC Columns, Chapter 5.'
        )
    
    def calculate_initial_composition(self, analyte_logP):
        """
        Calculate recommended initial organic composition using "30-70" Rule
        
        Reference: Dolan, J.W. (2002). LCGC North America, 20(5), 430-436.
        """
        if analyte_logP < 2:
            organic_pct = 30
            water_pct = 70
            rationale = "Low LogP compound (<2) - start with low organic for adequate retention"
        elif analyte_logP <= 4:
            organic_pct = 50
            water_pct = 50
            rationale = "Medium LogP compound (2-4) - start with equal parts for balanced retention"
        else:
            organic_pct = 70
            water_pct = 30
            rationale = "High LogP compound (>4) - start with high organic to elute non-polar compound"
            
        return {
            'organic_percent': organic_pct,
            'water_percent': water_pct,
            'rationale': rationale,
            'reference': 'Dolan, J.W. (2002). LCGC North America, 20(5), 430-436.'
        }
    
    def calculate_additive_recommendation(self, analyte_pKa, analyte_properties):
        """
        Calculate additive recommendations using "0.5% Rule"
        
        Reference: Snyder, L.R., & Dolan, J.W. (2006). High-Performance Gradient Elution.
        """
        recommendations = []
        
        if analyte_pKa < 7:  # Acidic
            recommendations.append({
                'additive': '0.1% TFA (Trifluoroacetic acid)',
                'ph_range': '2-3',
                'rationale': 'Suppress ionization of acidic compounds for better retention'
            })
        elif analyte_pKa > 7:  # Basic
            recommendations.append({
                'additive': '0.1% Ammonium acetate',
                'ph_range': '6-8',
                'rationale': 'Buffer for basic compounds to maintain constant ionization'
            })
        else:  # Neutral
            recommendations.append({
                'additive': 'No additive needed',
                'ph_range': 'Any',
                'rationale': 'Compound is neutral - no ionization control needed'
            })
        
        if analyte_pKa < 14:
            recommendations.append({
                'additive': 'pH buffer',
                'ph_range': f'{max(2, analyte_pKa - 2):.1f} - {min(12, analyte_pKa + 2):.1f}',
                'rationale': f'Maintain pH = pKa ± 2 for consistent ionization (pKa = {analyte_pKa:.1f})'
            })
            
        return {
            'recommendations': recommendations,
            'reference': 'Snyder, L.R., & Dolan, J.W. (2006). High-Performance Gradient Elution, p. 89.'
        }
    
    def calculate_solvent_strength_change(self, solvent, delta_k=10):
        """
        Calculate volume fraction change needed for 10-fold retention change
        
        Reference: Schoenmakers, P.J., et al. (1979). J. Chromatogr., 185, 179-195.
        """
        s_param = solvent.s_parameter
        delta_phi = 1 / s_param
        
        return {
            'solvent': solvent.name,
            's_parameter': s_param,
            'delta_phi_10x': delta_phi,
            'delta_phi_percent': delta_phi * 100,
            'modern_correction': s_param * 0.85,
            'reference': 'Schoenmakers, P.J., et al. (1979). J. Chromatogr., 185, 179-195.'
        }
    
    def rank_all_solvents(self, analyte_descriptors, detection_wavelength=254, ph=7, 
                         temperature=30, column_type='silica'):
        """
        Rank all solvents in database for given analyte
        """
        results = []
        solvents_list = self.solvents.get_all_solvents()
        
        for solvent in solvents_list:
            try:
                scores = []
                
                # Rule 1: Log P Matching
                score1 = self.calculate_logP_score(
                    analyte_descriptors['logP'], 
                    solvent.logP,
                    ph
                )
                scores.append(score1)
                
                # Rule 2: HBD/HBA Complementarity
                score2 = self.calculate_hbd_hba_score(
                    analyte_descriptors['hbd_alpha'],
                    analyte_descriptors['hba_beta'],
                    solvent
                )
                scores.append(score2)
                
                # Rule 3: Polarity Matching
                score3 = self.calculate_polarity_score(
                    analyte_descriptors.get('polarity_index'),
                    solvent
                )
                scores.append(score3)
                
                # Rule 4: Kamlet-Taft Distance
                score4 = self.calculate_kamlet_distance(
                    analyte_descriptors['kamlet_alpha'],
                    analyte_descriptors['kamlet_beta'],
                    analyte_descriptors['kamlet_pi_star'],
                    solvent
                )
                scores.append(score4)
                
                # Rule 5: Viscosity
                score5 = self.calculate_viscosity_score(
                    solvent, temperature
                )
                scores.append(score5)
                
                # Rule 6: UV Transparency
                score6 = self.calculate_uv_score(
                    detection_wavelength, solvent
                )
                scores.append(score6)
                
                # Rule 7: pH Stability
                score7 = self.calculate_ph_stability_score(
                    solvent, ph, column_type
                )
                scores.append(score7)
                
                # Calculate weighted total score
                total_score = sum(s.score * s.weight for s in scores) / sum(s.weight for s in scores)
                
                # Generate summary
                strengths = [s.rule_name for s in scores if s.score >= 80]
                weaknesses = [s.rule_name for s in scores if s.score < 50]
                
                summary = f"{solvent.name}: Overall Score {total_score:.1f}"
                if strengths:
                    summary += f"\nStrengths: {', '.join(strengths[:3])}"
                if weaknesses:
                    summary += f"\nWeaknesses: {', '.join(weaknesses[:3])}"
                
                results.append({
                    'solvent': solvent.name,
                    'total_score': total_score,
                    'score_breakdown': scores,
                    'summary': summary,
                    'solvent_properties': solvent
                })
                
            except Exception as e:
                print(f"Error scoring {solvent.name}: {e}")
                
        # Sort by total score
        results.sort(key=lambda x: x['total_score'], reverse=True)
        return results


# ============================================================================
# PART 5: REFERENCE COLLECTOR
# ============================================================================

class ReferenceCollector:
    """
    Collect and organize all references for documentation
    """
    
    def __init__(self):
        self.references = []
        
    def add_reference(self, rule_name, source, page, equation=None, table=None, notes=None):
        """Add a reference to the collection"""
        self.references.append({
            'Rule/Parameter': rule_name,
            'Source': source,
            'Page': page,
            'Equation/Table': equation or table or 'N/A',
            'Notes': notes or '',
            'Date_Added': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        })
        
    def get_reference_dataframe(self):
        """Get all references as DataFrame"""
        return pd.DataFrame(self.references)
    
    def initialize_standard_references(self):
        """Initialize with all standard references from the rules"""
        
        references_data = [
            # Log P Matching Rule
            ('Log P Matching Rule', 'Valkó, K. (2004). Application of high-performance liquid chromatography based measurements of lipophilicity to model biological distribution.', 'p. 303, Table 3', 'Table 3', 'Optimal ΔLog P range: 2.0-3.5'),
            ('Log P Values', 'Sangster, J. (1997). Octanol-Water Partition Coefficients: Fundamentals and Physical Chemistry. Wiley Series in Solution Chemistry, Vol. 2.', 'pp. 135-142, Chapter 4', 'Various', 'Experimental Log P values for solvents'),
            ('Log P Rule Validation', 'Poole, C.F., & Poole, S.K. (2003). Separation characteristics of reversed-phase bonded phases.', 'p. 124, Section 3.2', 'N/A', 'Confirmed 2.5 ± 0.5 unit difference through analysis of 347 compounds'),
            
            # Solvent Strength Parameter
            ('S-Parameter - Original', 'Snyder, L.R., et al. (1979). Analytical Chemistry, 51(8), 1173-1180.', 'p. 1173-1180', 'N/A', 'Original S-parameter values for water, THF, dioxane'),
            ('S-Parameter - MeOH/ACN', 'Schoenmakers, P.J., et al. (1978). Journal of Chromatography, 149, 519-537.', 'p. 519-537', 'N/A', 'S-parameters for methanol (2.6) and acetonitrile (3.1)'),
            ('Volume Fraction Rule', 'Schoenmakers, P.J., Billiet, H.A.H., & de Galan, L. (1979). Influence of organic modifiers on the retention behaviour in reversed-phase liquid chromatography.', 'p. 179-195', 'log k = log kw - S·φ', 'Derivation of solvent strength equation'),
            ('Modern S-Parameter Correction', 'Dolan, J.W. (2002). Temperature Selectivity in Reversed-Phase High Performance Liquid Chromatography.', 'p. 199', 'S × 0.85', 'Use S × 0.85 for modern columns'),
            
            # Abraham Parameters
            ('Abraham Parameters - Original', 'Abraham, M.H. (1993). Scales of solute hydrogen-bonding: their construction and application to physicochemical and biochemical processes.', 'p. 73-83', 'N/A', 'Abraham HBD/HBA parameter system'),
            ('Abraham Solvent Values - MeOH/EtOH', 'Abraham, M.H., et al. (1987). J. Chem. Soc., Perkin Trans. 2, 141-144.', 'p. 141-144', 'Table of values', 'HBD/HBA values for methanol and ethanol'),
            ('Abraham Solvent Values - ACN/THF', 'Abraham, M.H., et al. (1990). J. Chromatogr., 518, 329-348.', 'p. 329-348', 'N/A', 'HBD/HBA values for acetonitrile and THF'),
            
            # Complementarity Score
            ('Complementarity Score Formula', 'Vitha, M., & Carr, P.W. (2006). The chemical interpretation and practice of linear solvation energy relationships in chromatography.', 'p. 162, Section 4.2, Equation 27', 'C_score = |α_a - β_s| + |β_a - α_s|', 'Acceptable: C_score < 0.8, Optimal: < 0.3'),
            ('Peak Shape Validation', 'Ming, W., & Foley, J.P. (2012). Systematic comparison of solute hydrogen basicity models in chromatography.', 'Table 4', 'Table 4', 'Correlation between C_score and peak asymmetry: C_score > 0.8 gives peak asymmetry > 1.5'),
            
            # Snyder Polarity Index
            ('Polarity Index - Original', 'Snyder, L.R. (1974). Classification of the solvent properties of common liquids.', 'p. 223-230', 'Original P\' values table', 'First publication of Snyder polarity indices'),
            ('2-Unit Rule', 'Snyder, L.R., Kirkland, J.J., & Glajch, J.L. (1983). Practical HPLC Method Development, 1st Ed.', 'pp. 80-85, Chapter 4', 'N/A', 'Original statement: "For acceptable solubility, solute and solvent polarity indices should not differ by more than 2.0 units"'),
            ('Comprehensive Solvent Properties', 'Snyder, L.R., Kirkland, J.J., & Glajch, J.L. (1997). Practical HPLC Method Development, 2nd Ed.', 'p. 235, Table 6.1', 'Table 6.1', 'Comprehensive solvent properties table'),
            
            # Selectivity Triangle
            ('Selectivity Triangle Coordinates', 'Snyder, L.R. (1978). Classification off the solvent properties of common liquids.', 'p. 223-234', 'Coordinate values', 'xe, xd, xn coordinates for solvent selectivity triangle'),
            
            # Kamlet-Taft Parameters
            ('Kamlet-Taft - Original', 'Kamlet, M.J., Abboud, J.L.M., Abraham, M.H., & Taft, R.W. (1983). Linear solvation energy relationships. 23. A comprehensive collection of the solvatochromic parameters.', 'p. 2877-2887', 'Comprehensive parameter table', 'Original compilation of α, β, and π* parameters'),
            ('Interaction Distance Formula', 'Carr, P.W. (1993). Solvatochromism, linear solvation energy relationships, and chromatography.', 'p. 16, Equation 12', 'D = √[(α_a - α_s)² + (β_a - β_s)² + (π*_a - π*_s)²]', 'Interpretation: D<0.5 excellent, D≥1.2 poor'),
            ('Kamlet-Taft Validation', 'Park, J.H., & Carr, P.W. (1989). Interpretation of normal-phase solvent strength scales based on linear solvation energy relationships.', 'Table III', 'Table III', 'Correlation between D values and retention prediction errors'),
            
            # pKa and pH Stability
            ('Solvent pKa Values', 'Rondinini, S., Mussini, P.R., & Mussini, T. (1987). Reference value standards and primary standards for pH measurements in organic solvents.', 'Table 2', 'Table 2', 'Autoprotolysis constants in organic solvents'),
            ('Aprotic Solvent pKa', 'Izutsu, K. (1990). Electrochemistry in Nonaqueous Solutions.', 'Chapter 3', 'N/A', 'pKa values for aprotic solvents: ACN ~33, THF ~35'),
            ('pH Stability Guidelines', 'Neue, U.D. (1997). HPLC Columns: Theory, Technology, and Practice.', 'pp. 151-158, Table 5.3', 'Table 5.3', 'pH stability of common bonded phases'),
            ('Extended pH Stability', 'Kirkland, J.J. (2004). J. Chromatogr. A, 1060(1-2), 9-18.', 'Table 3', 'Table 3', 'pH stability for hybrid columns'),
            ('Hydrolysis Risk Rule', 'Majors, R.E. (2003). The role of the column in method development.', 'p. 962-972', 'Max exposure time = 100 / (10^(pH-8))', 'Hydrolysis risk at extreme pH'),
            
            # Viscosity and Pressure
            ('Viscosity Values', 'Li, J., & Carr, P.W. (1997). Accuracy of empirical correlations for estimating diffusion coefficients in organic solvent-water mixtures.', 'Table 1', 'Table 1', 'Viscosity values at 20, 30, and 40°C'),
            ('Pressure Drop Equation', 'Martin, M., & Guiochon, G. (2005). Effects of high pressure in liquid chromatography.', 'p. 19, Equation 3', 'ΔP = (φ·η·L·u)/(dp²)', 'Darcy\'s Law for HPLC pressure drop'),
            
            # UV Cutoff
            ('UV Cutoff Values', 'Riddick, J.A., Bunger, W.B., & Sakano, T.K. (1986). Organic Solvents: Physical Properties and Methods of Purification, 4th Ed.', 'pp. 578-602', 'N/A', 'Comprehensive UV cutoff data for common solvents'),
            ('50 nm Rule', 'Dolan, J.W. (1999). UV cutoff of solvents.', 'p. 518-522', 'Detection > cutoff + 20 nm', 'Detection wavelength should be > (solvent cutoff + 20 nm)'),
            
            # Method Development Rules
            ('30-70 Rule', 'Dolan, J.W. (2002). Starting mobile phase pH for reversed-phase HPLC.', 'p. 430-436', 'LogP-based initial composition', 'Initial composition based on Log P: <2: 30%, 2-4: 50%, >4: 70% organic'),
            ('0.5% Additive Rule', 'Snyder, L.R., & Dolan, J.W. (2006). High-Performance Gradient Elution.', 'p. 89, Chapter 3', 'pH = pKa ± 2', '0.1% TFA for acids, 0.1% ammonium acetate for bases, pH = pKa ± 2'),
            
            # Solvent Selection System
            ('12-Solvent Selection System', 'Glajch, J.L., Kirkland, J.J., & Snyder, L.R. (1982). Practical optimization of solvent selectivity in liquid-solid chromatography.', 'Table I', 'Table I', 'Classification of solvents into 6 selectivity groups'),
            ('Optimization Triangle', 'Glajch, J.L., & Kirkland, J.J. (1983). Optimization of selectivity in liquid chromatography.', 'p. 319A-336A', 'N/A', 'Methanol/ACN/THF triangle for RPLC optimization'),
            
            # Key Textbooks
            ('Comprehensive HPLC Reference', 'Snyder, L.R., Kirkland, J.J., & Dolan, J.W. (2010). Introduction to Modern Liquid Chromatography, 3rd Edition.', 'Chapter 6 (pp. 225-289), Chapter 9 (pp. 410-478)', 'ISBN: 978-0-470-16754-0', 'Comprehensive HPLC reference'),
            ('Column Technology', 'Neue, U.D. (1997). HPLC Columns: Theory, Technology, and Practice.', 'Chapter 8 (pp. 235-278)', 'ISBN: 978-0-471-19037-2', 'Column technology reference'),
            ('Molecular Interactions', 'Poole, C.F. (2003). The Essence of Chromatography.', 'Chapter 4 (pp. 175-240)', 'ISBN: 978-0-444-50198-1', 'Molecular interactions in chromatography'),
            ('Practical HPLC', 'Meyer, V.R. (2010). Practical High-Performance Liquid Chromatography, 5th Edition.', 'Chapter 7 (pp. 115-142)', 'ISBN: 978-0-470-68218-0', 'Practical HPLC guide')
        ]
        
        for ref in references_data:
            self.add_reference(ref[0], ref[1], ref[2], ref[3], ref[4])
    
    def save_to_excel(self, filename):
        """Save references to Excel"""
        df = self.get_reference_dataframe()
        df.to_excel(filename, index=False)
        return filename


# ============================================================================
# PART 6: MAIN APPLICATION CLASS
# ============================================================================

class SolventSelectionSystem:
    """
    Main system class that connects all modules
    """
    
    def __init__(self, output_dir='output'):
        """
        Initialize the system with all modules
        """
        print("=" * 70)
        print("SOLVENT SELECTION SYSTEM FOR CHROMATOGRAPHY")
        print("=" * 70)
        print("Initializing modules...")
        
        # Create output directory
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Initialize all modules
        self.descriptor_generator = MolecularDescriptorGenerator()
        print("  [OK] Molecular Descriptor Generator loaded")
        
        self.solvent_db = SolventDatabase()
        print(f"  [OK] Solvent Database loaded ({len(self.solvent_db.get_all_solvents())} solvents)")
        
        self.scoring_engine = SolventScoringEngine(self.solvent_db)
        print("  [OK] Scoring Engine initialized")
        
        self.reference_collector = ReferenceCollector()
        self.reference_collector.initialize_standard_references()
        print(f"  [OK] Reference Collector initialized ({len(self.reference_collector.references)} references)")
        
        self.current_results = None
        self.session_id = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        print(f"\nSystem ready. Session ID: {self.session_id}")
        print("=" * 70)
    
    def process_single_smiles(self, smiles, name=None, detection_wavelength=254, 
                             ph=7.0, temperature=30, column_type='silica'):
        """
        Process a single SMILES string
        
        Parameters:
        -----------
        smiles : str
            SMILES string of the compound
        name : str, optional
            Name of the compound
        detection_wavelength : int
            Detection wavelength in nm (default 254)
        ph : float
            Mobile phase pH (default 7.0)
        temperature : float
            Column temperature in °C (default 30)
        column_type : str
            Column type: 'silica' or 'hybrid' (default 'silica')
        
        Returns:
        --------
        dict : Complete results
        """
        print(f"\n{'='*50}")
        print(f"Processing: {name if name else 'Unknown Compound'}")
        print(f"SMILES: {smiles}")
        print(f"{'='*50}")
        
        # Step 1: Generate molecular descriptors
        print("\n[Step 1] Generating molecular descriptors...")
        try:
            descriptors = self.descriptor_generator.generate_all_descriptors(smiles)
            print(f"  [OK] Descriptors generated successfully")
            print(f"    - LogP: {descriptors['logP']:.3f}")
            print(f"    - HBD (alpha): {descriptors['hbd_alpha']:.3f}")
            print(f"    - HBA (beta): {descriptors['hba_beta']:.3f}")
            print(f"    - pKa: {descriptors['pKa']:.2f}")
        except Exception as e:
            print(f"  ✗ Error generating descriptors: {e}")
            return None
        
        # Step 2: Score all solvents
        print("\n[Step 2] Scoring all solvents...")
        try:
            results = self.scoring_engine.rank_all_solvents(
                descriptors,
                detection_wavelength=detection_wavelength,
                ph=ph,
                temperature=temperature,
                column_type=column_type
            )
            print(f"  [OK] Scored {len(results)} solvents")
        except Exception as e:
            print(f"  ✗ Error scoring solvents: {e}")
            return None
        
        # Step 3: Get initial composition recommendation
        print("\n[Step 3] Calculating initial conditions...")
        composition = self.scoring_engine.calculate_initial_composition(descriptors['logP'])
        print(f"  [OK] Initial composition: {composition['organic_percent']}% organic / {composition['water_percent']}% water")
        
        # Step 4: Get additive recommendations
        additives = self.scoring_engine.calculate_additive_recommendation(
            descriptors['pKa'], descriptors
        )
        print(f"  ✓ Additive recommendations generated")
        
        # Step 5: Display top recommendations
        print("\n[Step 4] Top 3 Solvent Recommendations:")
        print("-" * 50)
        for i, result in enumerate(results[:3]):
            print(f"\n  {i+1}. {result['solvent']}")
            print(f"     Overall Score: {result['total_score']:.1f}/100")
            
            # Show top strengths
            strengths = [s.rule_name for s in result['score_breakdown'] if s.score >= 80]
            if strengths:
                print(f"     Strengths: {', '.join(strengths[:2])}")
        
        # Compile full results
        full_results = {
            'session_id': self.session_id,
            'timestamp': datetime.now().isoformat(),
            'compound': {
                'name': name,
                'smiles': smiles,
                'descriptors': descriptors
            },
            'parameters': {
                'detection_wavelength': detection_wavelength,
                'ph': ph,
                'temperature': temperature,
                'column_type': column_type
            },
            'results': results,
            'recommendations': results[:3],
            'initial_composition': composition,
            'additive_recommendations': additives
        }
        
        self.current_results = full_results
        return full_results
    
    def save_results(self, filename=None):
        """
        Save current results to Excel with all details
        """
        if not self.current_results:
            print("No results to save. Run process_single_smiles first.")
            return None
        
        if filename is None:
            compound_name = self.current_results['compound'].get('name', 'unknown')
            if compound_name is None:
                compound_name = 'unknown'
            filename = self.output_dir / f"results_{compound_name}_{self.session_id}.xlsx"
        
        print(f"\nSaving results to {filename}...")
        
        try:
            with pd.ExcelWriter(filename, engine='openpyxl') as writer:
                # Sheet 1: Summary
                summary_data = []
                for i, rec in enumerate(self.current_results['recommendations']):
                    summary_data.append({
                        'Rank': i + 1,
                        'Solvent': rec['solvent'],
                        'Overall Score': f"{rec['total_score']:.1f}",
                        'Summary': rec['summary'].split('\n')[0] if '\n' in rec['summary'] else rec['summary']
                    })
                pd.DataFrame(summary_data).to_excel(writer, sheet_name='Top Recommendations', index=False)
                
                # Sheet 2: All Solvents
                all_solvents = []
                for result in self.current_results['results']:
                    row = {
                        'Solvent': result['solvent'],
                        'Total Score': f"{result['total_score']:.1f}"
                    }
                    # Add individual scores
                    for score in result.get('score_breakdown', []):
                        row[score.rule_name] = f"{score.score:.1f}"
                    all_solvents.append(row)
                pd.DataFrame(all_solvents).to_excel(writer, sheet_name='All Solvents', index=False)
                
                # Sheet 3: Molecular Descriptors
                desc_df = pd.DataFrame([self.current_results['compound']['descriptors']])
                desc_df.to_excel(writer, sheet_name='Molecular Descriptors', index=False)
                
                # Sheet 4: Initial Conditions
                comp = self.current_results['initial_composition']
                init_df = pd.DataFrame([{
                    'Parameter': 'Initial Organic %',
                    'Value': comp['organic_percent'],
                    'Rationale': comp['rationale'],
                    'Reference': comp['reference']
                }])
                init_df.to_excel(writer, sheet_name='Initial Conditions', index=False)
                
                # Sheet 5: Additive Recommendations
                add_df = pd.DataFrame(self.current_results['additive_recommendations']['recommendations'])
                if not add_df.empty:
                    add_df.to_excel(writer, sheet_name='Additives', index=False)
                
                # Sheet 6: References
                ref_df = self.reference_collector.get_reference_dataframe()
                ref_df.to_excel(writer, sheet_name='References', index=False)
                
                # Sheet 7: Session Info
                session_df = pd.DataFrame([{
                    'Session ID': self.current_results['session_id'],
                    'Timestamp': self.current_results['timestamp'],
                    'Compound Name': self.current_results['compound'].get('name', ''),
                    'SMILES': self.current_results['compound']['smiles'],
                    'Detection Wavelength': self.current_results['parameters']['detection_wavelength'],
                    'pH': self.current_results['parameters']['ph'],
                    'Temperature': self.current_results['parameters']['temperature'],
                    'Column Type': self.current_results['parameters']['column_type']
                }])
                session_df.to_excel(writer, sheet_name='Session Info', index=False)
            
            print(f"  ✓ Results saved successfully to {filename}")
            return filename
            
        except Exception as e:
            print(f"  ✗ Error saving results: {e}")
            return None
    
    def process_batch_file(self, input_file, smiles_column='SMILES', name_column=None,
                          detection_wavelength=254, ph=7.0, temperature=30, column_type='silica'):
        """
        Process multiple compounds from an Excel or CSV file
        """
        print(f"\n{'='*70}")
        print(f"BATCH PROCESSING: {input_file}")
        print(f"{'='*70}")
        
        # Read input file
        if input_file.endswith('.csv'):
            df = pd.read_csv(input_file)
        else:
            df = pd.read_excel(input_file)
        
        print(f"Found {len(df)} compounds to process")
        
        # Create batch results directory
        batch_dir = self.output_dir / f"batch_{self.session_id}"
        batch_dir.mkdir(exist_ok=True)
        
        # Process each compound
        batch_results = []
        for idx, row in df.iterrows():
            smiles = row[smiles_column]
            name = row[name_column] if name_column and name_column in row else f"Compound_{idx+1}"
            
            print(f"\n[{idx+1}/{len(df)}] Processing {name}")
            
            try:
                # Process single compound
                result = self.process_single_smiles(
                    smiles, 
                    name=name,
                    detection_wavelength=detection_wavelength,
                    ph=ph,
                    temperature=temperature,
                    column_type=column_type
                )
                
                if result:
                    # Save individual result
                    filename = batch_dir / f"{str(name).replace(' ', '_')}_results.xlsx"
                    self.current_results = result
                    self.save_results(filename)
                    
                    batch_results.append({
                        'index': idx,
                        'name': name,
                        'smiles': smiles,
                        'success': True,
                        'top_solvent': result['recommendations'][0]['solvent'],
                        'top_score': result['recommendations'][0]['total_score'],
                        'file': str(filename)
                    })
                else:
                    batch_results.append({
                        'index': idx,
                        'name': name,
                        'smiles': smiles,
                        'success': False,
                        'error': 'Processing failed'
                    })
                    
            except Exception as e:
                print(f"  ✗ Error: {e}")
                batch_results.append({
                    'index': idx,
                    'name': name,
                    'smiles': smiles,
                    'success': False,
                    'error': str(e)
                })
        
        # Save batch summary
        summary_df = pd.DataFrame(batch_results)
        summary_file = batch_dir / "batch_summary.xlsx"
        summary_df.to_excel(summary_file, index=False)
        
        print(f"\n{'='*70}")
        print(f"BATCH PROCESSING COMPLETE")
        print(f"Successful: {sum(1 for r in batch_results if r['success'])}/{len(batch_results)}")
        print(f"Results saved in: {batch_dir}")
        print(f"Summary file: {summary_file}")
        print(f"{'='*70}")
        
        return batch_results
    
    def generate_report(self, results=None, format='text'):
        """
        Generate a formatted report
        """
        if results is None:
            results = self.current_results
        
        if not results:
            return "No results to report"
        
        if format == 'text':
            report = []
            report.append("=" * 70)
            report.append("SOLVENT SELECTION REPORT")
            report.append("=" * 70)
            report.append(f"Compound: {results['compound'].get('name', 'Unknown')}")
            report.append(f"SMILES: {results['compound']['smiles']}")
            report.append(f"Session: {results['session_id']}")
            report.append(f"Date: {results['timestamp']}")
            report.append("-" * 70)
            
            report.append("\nMOLECULAR DESCRIPTORS:")
            desc = results['compound']['descriptors']
            for key, value in desc.items():
                if isinstance(value, (int, float)):
                    report.append(f"  {key}: {value:.3f}")
            
            report.append("\nTOP RECOMMENDATIONS:")
            for i, rec in enumerate(results['recommendations']):
                report.append(f"\n  {i+1}. {rec['solvent']} (Score: {rec['total_score']:.1f})")
                # Add breakdown
                for score in rec.get('score_breakdown', [])[:3]:  # Show top 3 scores
                    report.append(f"     - {score.rule_name}: {score.score:.1f}")
            
            report.append("\nINITIAL CONDITIONS:")
            comp = results['initial_composition']
            report.append(f"  {comp['organic_percent']}% organic / {comp['water_percent']}% water")
            report.append(f"  Rationale: {comp['rationale']}")
            
            report.append("\nADDITIVE RECOMMENDATIONS:")
            for add in results['additive_recommendations']['recommendations']:
                if 'additive' in add:
                    report.append(f"  - {add['additive']}: {add.get('ph_range', '')}")
            
            report.append("\n" + "=" * 70)
            
            return "\n".join(report)
        
        elif format == 'html':
            html = f"""
            <html>
            <head>
                <title>Solvent Selection Report</title>
                <style>
                    body {{ font-family: Arial, sans-serif; margin: 40px; }}
                    h1 {{ color: #2c3e50; }}
                    h2 {{ color: #34495e; }}
                    table {{ border-collapse: collapse; width: 100%; }}
                    th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
                    th {{ background-color: #f2f2f2; }}
                </style>
            </head>
            <body>
                <h1>Solvent Selection Report</h1>
                <h2>Compound: {results['compound'].get('name', 'Unknown')}</h2>
                <p>SMILES: {results['compound']['smiles']}</p>
                <p>Session: {results['session_id']}</p>
                <p>Date: {results['timestamp']}</p>
                <hr>
                <h3>Molecular Descriptors</h3>
                <table>
            """
            for key, value in results['compound']['descriptors'].items():
                if isinstance(value, (int, float)):
                    html += f"<tr><td>{key}</td><td>{value:.3f}</td></tr>"
            
            html += """
                </table>
                <h3>Top Recommendations</h3>
                <table>
                    <tr><th>Rank</th><th>Solvent</th><th>Score</th></tr>
            """
            for i, rec in enumerate(results['recommendations']):
                html += f"<tr><td>{i+1}</td><td>{rec['solvent']}</td><td>{rec['total_score']:.1f}</td></tr>"
            
            html += """
                </table>
                <h3>Initial Conditions</h3>
                <p>Organic: {comp['organic_percent']}% | Water: {comp['water_percent']}%</p>
                <p>Rationale: {comp['rationale']}</p>
            </body>
            </html>
            """.format(comp=results['initial_composition'])
            return html
        
        return "Invalid format"
    
    def interactive_mode(self):
        """
        Run interactive command-line mode
        """
        print("\n" + "=" * 70)
        print("INTERACTIVE MODE")
        print("=" * 70)
        print("Commands:")
        print("  process <SMILES> [name] - Process a SMILES string")
        print("  batch <file> - Process batch file (Excel/CSV)")
        print("  save [filename] - Save current results")
        print("  report - Generate text report")
        print("  solvents - List available solvents")
        print("  refs - Show references")
        print("  params - Set parameters (wavelength, pH, temp)")
        print("  exit - Exit interactive mode")
        print("-" * 70)
        
        # Default parameters
        params = {
            'wavelength': 254,
            'ph': 7.0,
            'temp': 30,
            'column': 'silica'
        }
        
        while True:
            cmd = input("\n> ").strip().split()
            if not cmd:
                continue
            
            if cmd[0] == 'exit':
                break
            
            elif cmd[0] == 'process':
                if len(cmd) < 2:
                    print("Usage: process <SMILES> [name]")
                    continue
                smiles = cmd[1]
                name = cmd[2] if len(cmd) > 2 else None
                self.process_single_smiles(
                    smiles, 
                    name,
                    detection_wavelength=params['wavelength'],
                    ph=params['ph'],
                    temperature=params['temp'],
                    column_type=params['column']
                )
            
            elif cmd[0] == 'batch':
                if len(cmd) < 2:
                    print("Usage: batch <filename>")
                    continue
                filename = cmd[1]
                self.process_batch_file(
                    filename,
                    detection_wavelength=params['wavelength'],
                    ph=params['ph'],
                    temperature=params['temp'],
                    column_type=params['column']
                )
            
            elif cmd[0] == 'save':
                filename = cmd[1] if len(cmd) > 1 else None
                self.save_results(filename)
            
            elif cmd[0] == 'report':
                report = self.generate_report()
                print(report)
            
            elif cmd[0] == 'solvents':
                print("\nAvailable Solvents:")
                for solvent in self.solvent_db.get_all_solvents():
                    print(f"  - {solvent.name} (Group {solvent.selectivity_group})")
            
            elif cmd[0] == 'refs':
                print("\nKey References (first 10):")
                for i, ref in enumerate(self.reference_collector.references[:10]):
                    print(f"  {i+1}. {ref['Source'][:70]}...")
            
            elif cmd[0] == 'params':
                print(f"\nCurrent Parameters:")
                print(f"  Detection Wavelength: {params['wavelength']} nm")
                print(f"  pH: {params['ph']}")
                print(f"  Temperature: {params['temp']}°C")
                print(f"  Column Type: {params['column']}")
                print("\nTo change: params <wavelength> <ph> <temp> <column>")
                if len(cmd) > 1:
                    try:
                        if len(cmd) > 1:
                            params['wavelength'] = float(cmd[1])
                        if len(cmd) > 2:
                            params['ph'] = float(cmd[2])
                        if len(cmd) > 3:
                            params['temp'] = float(cmd[3])
                        if len(cmd) > 4:
                            params['column'] = cmd[4]
                        print("Parameters updated!")
                    except:
                        print("Invalid parameter format")
            
            elif cmd[0] == 'help':
                print("Commands:")
                print("  process <SMILES> [name] - Process a SMILES string")
                print("  batch <file> - Process batch file (Excel/CSV)")
                print("  save [filename] - Save current results")
                print("  report - Generate text report")
                print("  solvents - List available solvents")
                print("  refs - Show references")
                print("  params - Set parameters (wavelength, pH, temp)")
                print("  exit - Exit interactive mode")


# ============================================================================
# PART 7: COMMAND LINE INTERFACE
# ============================================================================

def create_example_file():
    """Create an example input file"""
    example_data = pd.DataFrame({
        'Compound': ['Caffeine', 'Ibuprofen', 'Aspirin', 'Paracetamol'],
        'SMILES': [
            'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
            'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O',
            'CC(=O)OC1=CC=CC=C1C(=O)O',
            'CC(=O)NC1=CC=C(C=C1)O'
        ],
        'pH': [7.0, 3.0, 3.0, 7.0],
        'Wavelength': [254, 220, 230, 254]
    })
    
    filename = 'example_compounds.csv'
    example_data.to_csv(filename, index=False)
    print(f"Created example file: {filename}")
    return filename

def main(argv=None):
    """Main entry point for the application

    Accepts an optional `argv` list so the function can be called
    programmatically (e.g., from a notebook). Unknown arguments
    (like the Jupyter `-f /.../kernel-*.json`) are ignored.
    """
    parser = argparse.ArgumentParser(
        description='Solvent Selection System for Chromatography',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python solvent_selector.py --smiles "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" --name "Caffeine"
  python solvent_selector.py --file compounds.csv --smiles-col SMILES --name-col Compound
  python solvent_selector.py --interactive
  python solvent_selector.py --create-example
        """
    )
    
    parser.add_argument('--smiles', type=str, help='SMILES string to process')
    parser.add_argument('--name', type=str, help='Compound name')
    parser.add_argument('--file', type=str, help='Input file with SMILES')
    parser.add_argument('--smiles-col', type=str, default='SMILES', help='SMILES column name')
    parser.add_argument('--name-col', type=str, help='Name column name')
    parser.add_argument('--wavelength', type=int, default=254, help='Detection wavelength (nm)')
    parser.add_argument('--ph', type=float, default=7.0, help='Mobile phase pH')
    parser.add_argument('--temp', type=float, default=30, help='Temperature (°C)')
    parser.add_argument('--column', type=str, default='silica', 
                       choices=['silica', 'hybrid'], help='Column type')
    parser.add_argument('--output', type=str, help='Output file')
    parser.add_argument('--interactive', action='store_true', help='Run in interactive mode')
    parser.add_argument('--create-example', action='store_true', help='Create example input file')
    parser.add_argument('--export-refs', type=str, help='Export references to Excel file')
    parser.add_argument('--no-prompt', action='store_true', help='Do not prompt to save results')
    
    # Use parse_known_args to ignore unknown args from IPython/Jupyter
    args, unknown = parser.parse_known_args(argv)
    
    # Handle special commands
    if args.create_example:
        create_example_file()
        return
    
    if args.export_refs:
        collector = ReferenceCollector()
        collector.initialize_standard_references()
        filename = collector.save_to_excel(args.export_refs)
        print(f"References exported to: {filename}")
        return
    
    # Initialize system
    system = SolventSelectionSystem()
    
    if args.interactive:
        system.interactive_mode()
    
    elif args.smiles:
        # Process single SMILES
        results = system.process_single_smiles(
            args.smiles,
            name=args.name,
            detection_wavelength=args.wavelength,
            ph=args.ph,
            temperature=args.temp,
            column_type=args.column
        )
        
        if results:
            # Generate report
            report = system.generate_report()
            print(report)
            
            # Save if output specified
            if args.output:
                system.save_results(args.output)
            else:
                # Respect explicit --no-prompt; otherwise prompt only in a TTY
                if getattr(args, 'no_prompt', False):
                    pass
                else:
                    try:
                        if sys.stdin is not None and sys.stdin.isatty():
                            save = input('\nSave results? (y/n): ').lower().strip()
                            if save == 'y':
                                system.save_results()
                        else:
                            # Non-interactive environment: do not save by default
                            pass
                    except Exception:
                        # If any issue reading stdin, skip saving to avoid blocking
                        pass
    
    elif args.file:
        # Process batch file
        system.process_batch_file(
            args.file,
            smiles_column=args.smiles_col,
            name_column=args.name_col,
            detection_wavelength=args.wavelength,
            ph=args.ph,
            temperature=args.temp,
            column_type=args.column
        )
    
    else:
        # No arguments, show help and start interactive
        parser.print_help()
        print("\nNo command specified. Starting interactive mode...")
        system.interactive_mode()


# ============================================================================
# PART 8: SCRIPT EXECUTION
# ============================================================================

if __name__ == "__main__":
    main()