"""
Buffer Selection Engine for HPLC Method Development
Based on 10 Comprehensive Rules with Direct References
"""

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen, MolSurf
from rdkit.Chem.rdMolDescriptors import CalcNumHBA, CalcNumHBD
import logging
from datetime import datetime
import openpyxl

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class MolecularDescriptorCalculator:
    """Calculate molecular descriptors needed for buffer selection rules"""
    
    def __init__(self, smiles):
        self.smiles = smiles
        self.mol = Chem.MolFromSmiles(smiles)
        if self.mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles}")
        
    def calculate_all_descriptors(self):
        """Calculate all molecular descriptors required for the 10 rules"""
        
        # Basic molecular properties
        descriptors = {
            'molecular_weight': Descriptors.MolWt(self.mol),
            'logP': Descriptors.MolLogP(self.mol),
            'tpsa': Descriptors.TPSA(self.mol),
            'num_hba': CalcNumHBA(self.mol),  # Hydrogen bond acceptors
            'num_hbd': CalcNumHBD(self.mol),  # Hydrogen bond donors
            'num_rotatable_bonds': Descriptors.NumRotatableBonds(self.mol),
            'num_rings': Descriptors.RingCount(self.mol),
            'num_aromatic_rings': Descriptors.NumAromaticRings(self.mol),
        }
        
        # Functional group detection (Rule 7 - Reactivity)
        descriptors.update(self._detect_functional_groups())
        
        # Metal binding sites (Rule 7.3 and Rule 10)
        descriptors.update(self._detect_metal_binding_sites())
        
        # UV chromophores (Rule 2 - UV detection)
        descriptors.update(self._estimate_uv_absorption())
        
        # Charge state prediction at different pH
        descriptors.update(self._predict_charge_states())
        
        # Stability indicators (Rule 9)
        descriptors.update(self._assess_chemical_stability())
        
        return descriptors
    
    def _detect_functional_groups(self):
        """Detect functional groups relevant to buffer compatibility rules"""
        mol = self.mol
        
        # SMARTS patterns for functional group detection
        patterns = {
            'has_aldehyde': '[CX3H1](=O)[#6]',  # Aldehyde group
            'has_ketone': '[#6][CX3](=O)[#6]',  # Ketone group
            'has_ester': '[#6][CX3](=O)O[#6]',  # Ester group
            'has_amide': '[#6][CX3](=O)N[#6]',  # Amide group
            'has_carboxylic_acid': '[CX3](=O)[OX2H1]',  # Carboxylic acid
            'has_primary_amine': '[NX3;H2]',  # Primary amine
            'has_secondary_amine': '[NX3;H1]',  # Secondary amine
            'has_tertiary_amine': '[NX3;H0]',  # Tertiary amine
            'has_phosphate': '[PX4](=O)([OX2H])([OX2H])[OX2H,OX1-]',  # Phosphate
            'has_vicinal_diol': '[CH1X4]([OX2H])[CH1X4]([OX2H])',  # Vicinal diol
            'has_catechol': '[c]1[c][c][c][c][c]1-[OX2H]',  # Catechol-like
            'has_nucleoside': '[n]1[c][n][c][c]1',  # Purine/pyrimidine base
            'has_thiol': '[SX2H]',  # Thiol group
            'has_disulfide': '[SX2][SX2]',  # Disulfide bond
        }
        
        functional_groups = {}
        for name, pattern in patterns.items():
            patt = Chem.MolFromSmarts(pattern)
            functional_groups[name] = mol.HasSubstructMatch(patt)
            
            # Count occurrences for some groups
            if functional_groups[name]:
                matches = mol.GetSubstructMatches(patt)
                functional_groups[f'{name}_count'] = len(matches)
        
        return functional_groups
    
    def _detect_metal_binding_sites(self):
        """Detect potential metal binding sites (Rule 7.3 and Rule 10)"""
        mol = self.mol
        
        # Metal binding patterns
        metal_binding_patterns = {
            'carboxylate_binding': '[CX3](=O)[OX1-]',  # Carboxylate anion
            'hydroxamate_binding': '[NX3]C(=O)[OX2H]',  # Hydroxamic acid
            'thiolate_binding': '[SX2H]',  # Thiol
            'imidazole_binding': '[n]1[c][n][c][c]1',  # Imidazole (histidine)
            'phosphate_binding': '[PX4](=O)([OX2H,OX1-])([OX2H,OX1-])[OX2H,OX1-]',  # Phosphate
            'catechol_binding': '[c]1([OX2H])[c]([OX2H])[c][c][c]1',  # Catechol
            'phenolate_binding': '[c][OX1-]',  # Phenolate
        }
        
        metal_sites = {}
        for name, pattern in metal_binding_patterns.items():
            patt = Chem.MolFromSmarts(pattern)
            metal_sites[f'has_{name}'] = mol.HasSubstructMatch(patt)
            
            if metal_sites[f'has_{name}']:
                matches = mol.GetSubstructMatches(patt)
                metal_sites[f'{name}_count'] = len(matches)
        
        # Calculate chelation potential (0-5 scale)
        chelation_score = 0
        if metal_sites.get('has_carboxylate_binding', False):
            chelation_score += 2
        if metal_sites.get('has_thiolate_binding', False):
            chelation_score += 3
        if metal_sites.get('has_imidazole_binding', False):
            chelation_score += 2
        if metal_sites.get('has_catechol_binding', False):
            chelation_score += 4
        if metal_sites.get('has_phosphate_binding', False):
            chelation_score += 3
            
        metal_sites['chelation_potential'] = min(chelation_score, 5)
        
        return metal_sites
    
    def _estimate_uv_absorption(self):
        """Estimate UV absorption characteristics (Rule 2)"""
        mol = self.mol
        
        # Chromophore detection
        chromophores = {
            'has_benzene': '[c]1[c][c][c][c][c]1',  # Benzene ring
            'has_phenol': '[c]1([OX2H])[c][c][c][c]1',  # Phenol
            'has_aniline': '[c]1([NX3])[c][c][c][c]1',  # Aniline
            'has_nitro': '[NX3+](=O)[O-]',  # Nitro group
            'has_nitroso': '[NX2]=O',  # Nitroso group
            'has_azo': '[NX2]=[NX2]',  # Azo group
            'has_carbonyl': '[CX3]=O',  # Carbonyl
            'has_conjugated_diene': '[CX3]=[CX3][CX3]=[CX3]',  # Conjugated diene
            'has_anthracene': '[c]12[c][c][c][c]1[c][c][c][c]2',  # Anthracene
            'has_quinone': '[c]1(=O)[c][c][c](=O)[c][c]1',  # Quinone
        }
        
        uv_properties = {}
        for name, pattern in chromophores.items():
            patt = Chem.MolFromSmarts(pattern)
            uv_properties[name] = mol.HasSubstructMatch(patt)
        
        # Estimate lambda_max based on chromophores
        # Reference: Snyder, L.R., Kirkland, J.J., & Dolan, J.W. (2011). Introduction to Modern Liquid Chromatography (3rd ed.). Chapter 5.
        estimated_lambda_max = 190  # Base UV cutoff
        
        if uv_properties.get('has_benzene', False):
            estimated_lambda_max = max(estimated_lambda_max, 210)
        if uv_properties.get('has_phenol', False):
            estimated_lambda_max = max(estimated_lambda_max, 270)
        if uv_properties.get('has_aniline', False):
            estimated_lambda_max = max(estimated_lambda_max, 280)
        if uv_properties.get('has_nitro', False):
            estimated_lambda_max = max(estimated_lambda_max, 270)
        if uv_properties.get('has_azo', False):
            estimated_lambda_max = max(estimated_lambda_max, 400)
        if uv_properties.get('has_anthracene', False):
            estimated_lambda_max = max(estimated_lambda_max, 360)
        if uv_properties.get('has_quinone', False):
            estimated_lambda_max = max(estimated_lambda_max, 450)
        
        # Conjugation length effect
        conjugation_length = self._estimate_conjugation_length()
        if conjugation_length > 2:
            estimated_lambda_max += (conjugation_length - 2) * 30
        
        uv_properties['estimated_lambda_max_nm'] = estimated_lambda_max
        uv_properties['recommended_detection_wavelength'] = estimated_lambda_max + 20  # Rule: +10-20 nm margin
        
        return uv_properties
    
    def _estimate_conjugation_length(self):
        """Estimate the length of conjugated systems"""
        mol = self.mol
        # Simplified conjugation detection
        # Count alternating double bonds
        conjugated_systems = []
        for atom in mol.GetAtoms():
            if atom.GetIsAromatic():
                # Count ring size as conjugation
                ring_info = mol.GetRingInfo()
                for ring in ring_info.AtomRings():
                    if atom.GetIdx() in ring:
                        conjugated_systems.append(len(ring))
                        break
        
        return max(conjugated_systems) if conjugated_systems else 0
    
    def _predict_charge_states(self):
        """Predict charge states at different pH values (for Rule 1, 4, 8)"""
        mol = self.mol
        
        # Calculate pKa predictions using built-in rules
        # Simplified pKa prediction based on functional groups
        pka_values = []
        
        # Carboxylic acids (pKa ~4-5)
        carboxylic_pattern = Chem.MolFromSmarts('[CX3](=O)[OX2H]')
        if mol.HasSubstructMatch(carboxylic_pattern):
            pka_values.append(4.5)  # Average pKa for carboxylic acids
        
        # Phenols (pKa ~10)
        phenol_pattern = Chem.MolFromSmarts('[c]1([OX2H])[c][c][c][c]1')
        if mol.HasSubstructMatch(phenol_pattern):
            pka_values.append(10.0)
        
        # Primary amines (pKa ~9-10)
        primary_amine_pattern = Chem.MolFromSmarts('[NX3;H2]')
        if mol.HasSubstructMatch(primary_amine_pattern):
            pka_values.append(9.5)
        
        # Secondary amines (pKa ~10-11)
        secondary_amine_pattern = Chem.MolFromSmarts('[NX3;H1]')
        if mol.HasSubstructMatch(secondary_amine_pattern):
            pka_values.append(10.5)
        
        # Tertiary amines (pKa ~9-10)
        tertiary_amine_pattern = Chem.MolFromSmarts('[NX3;H0]')
        if mol.HasSubstructMatch(tertiary_amine_pattern):
            pka_values.append(9.0)
        
        # Phosphates (pKa ~2, 7, 12)
        phosphate_pattern = Chem.MolFromSmarts('[PX4](=O)([OX2H])([OX2H])[OX2H]')
        if mol.HasSubstructMatch(phosphate_pattern):
            pka_values.extend([2.1, 7.2, 12.3])
        
        charge_info = {
            'pka_values': pka_values,
            'is_acidic': any(p < 7 for p in pka_values),
            'is_basic': any(p > 7 for p in pka_values),
            'is_amphoteric': any(p < 7 for p in pka_values) and any(p > 7 for p in pka_values)
        }
        
        return charge_info
    
    def _assess_chemical_stability(self):
        """Assess chemical stability for Rule 9"""
        mol = self.mol
        
        # Hydrolysis-sensitive groups
        hydrolysis_patterns = {
            'ester_hydrolysis': '[#6][CX3](=O)O[#6]',
            'amide_hydrolysis': '[#6][CX3](=O)N[#6]',
            'acetal_hydrolysis': '[#6][OX2][CX4H1]([OX2][#6])[#6]',
            'hemiacetal_hydrolysis': '[#6][OX2][CX4H1]([OX2H])[#6]',
            'lactone_hydrolysis': '[#6]1[CX3](=O)O[#6][#6]1',
        }
        
        # Oxidation-sensitive groups
        oxidation_patterns = {
            'thiol_oxidation': '[SX2H]',
            'phenol_oxidation': '[c]1([OX2H])[c][c][c][c]1',
            'catechol_oxidation': '[c]1([OX2H])[c]([OX2H])[c][c][c]1',
            'aldehyde_oxidation': '[CX3H1](=O)[#6]',
        }
        
        # Photolysis-sensitive groups
        photolysis_patterns = {
            'nitro_photolysis': '[NX3+](=O)[O-]',
            'azo_photolysis': '[NX2]=[NX2]',
            'peroxide_photolysis': '[OX2][OX2]',
        }
        
        stability = {
            'hydrolysis_risk': 0,
            'oxidation_risk': 0,
            'photolysis_risk': 0,
        }
        
        for pattern_name, pattern in hydrolysis_patterns.items():
            patt = Chem.MolFromSmarts(pattern)
            if mol.HasSubstructMatch(patt):
                stability['hydrolysis_risk'] += 2
                stability[f'has_{pattern_name}'] = True
        
        for pattern_name, pattern in oxidation_patterns.items():
            patt = Chem.MolFromSmarts(pattern)
            if mol.HasSubstructMatch(patt):
                stability['oxidation_risk'] += 2
                stability[f'has_{pattern_name}'] = True
        
        for pattern_name, pattern in photolysis_patterns.items():
            patt = Chem.MolFromSmarts(pattern)
            if mol.HasSubstructMatch(patt):
                stability['photolysis_risk'] += 3
                stability[f'has_{pattern_name}'] = True
        
        # Overall stability score (0-10, higher = less stable)
        stability['overall_stability_risk'] = (
            stability['hydrolysis_risk'] + 
            stability['oxidation_risk'] + 
            stability['photolysis_risk']
        )
        
        return stability


class BufferSelectionEngine:
    """
    Main engine for buffer selection based on 10 comprehensive rules
    Each rule is implemented with direct references
    """
    
    def __init__(self):
        self.initialize_buffer_database()
        self.initialize_reference_database()
        # Initialize parameters early so rule weight calculation can reference them safely
        self.params = {}
        self.rule_weights = self.initialize_rule_weights()
        
    def initialize_buffer_database(self):
        """Initialize complete buffer database with all properties from the 10 rules"""
        
        self.buffers = pd.DataFrame([
            # Format: [Name, pKa1, pKa2, pKa3, UV_cutoff_10mM, solubility_ACN, solubility_MeOH, 
            #          volatility_class, temp_coefficient, metal_chelation, storage_stability, ms_compatibility]
            
            # Phosphate buffers
            {'name': 'Phosphate_pKa1', 'base_name': 'Phosphate', 'pka': 2.14, 'pka_temp_coeff': -0.0028,
             'uv_cutoff_10mM': 190, 'uv_cutoff_50mM': 195, 'solubility_ACN_10mM': 65,
             'solubility_MeOH_10mM': 70, 'solubility_IPA_10mM': 60, 'volatility': 'non-volatile',
             'ms_compatible': False, 'max_ms_conc': 0, 'metal_chelation_fe': 3.5,
             'metal_chelation_ca': 2.0, 'storage_stability_days': 180, 'ph_drift_per_day': 0.001,
             'reactivity_risk': 5, 'catalyzes_hydrolysis': True, 'catalyzes_ester_hydrolysis': True,
             'reference': 'Goldberg, R.N. et al. (2002). J. Phys. Chem. Ref. Data, 31(2), 231-370.'},
            
            {'name': 'Phosphate_pKa2', 'base_name': 'Phosphate', 'pka': 7.20, 'pka_temp_coeff': -0.0028,
             'uv_cutoff_10mM': 190, 'uv_cutoff_50mM': 195, 'solubility_ACN_10mM': 65,
             'solubility_MeOH_10mM': 70, 'solubility_IPA_10mM': 60, 'volatility': 'non-volatile',
             'ms_compatible': False, 'max_ms_conc': 0, 'metal_chelation_fe': 3.5,
             'metal_chelation_ca': 2.0, 'storage_stability_days': 180, 'ph_drift_per_day': 0.001,
             'reactivity_risk': 5, 'catalyzes_hydrolysis': True, 'catalyzes_ester_hydrolysis': True,
             'reference': 'Goldberg, R.N. et al. (2002). J. Phys. Chem. Ref. Data, 31(2), 231-370.'},
            
            {'name': 'Phosphate_pKa3', 'base_name': 'Phosphate', 'pka': 12.30, 'pka_temp_coeff': -0.0028,
             'uv_cutoff_10mM': 190, 'uv_cutoff_50mM': 195, 'solubility_ACN_10mM': 65,
             'solubility_MeOH_10mM': 70, 'solubility_IPA_10mM': 60, 'volatility': 'non-volatile',
             'ms_compatible': False, 'max_ms_conc': 0, 'metal_chelation_fe': 3.5,
             'metal_chelation_ca': 2.0, 'storage_stability_days': 180, 'ph_drift_per_day': 0.001,
             'reactivity_risk': 5, 'catalyzes_hydrolysis': True, 'catalyzes_ester_hydrolysis': True,
             'reference': 'Goldberg, R.N. et al. (2002). J. Phys. Chem. Ref. Data, 31(2), 231-370.'},
            
            # Acetate
            {'name': 'Acetate', 'base_name': 'Acetate', 'pka': 4.76, 'pka_temp_coeff': 0.0002,
             'uv_cutoff_10mM': 210, 'uv_cutoff_50mM': 215, 'solubility_ACN_10mM': 85,
             'solubility_MeOH_10mM': 90, 'solubility_IPA_10mM': 80, 'volatility': 'moderate',
             'ms_compatible': True, 'max_ms_conc': 20, 'metal_chelation_fe': 1.5,
             'metal_chelation_ca': 0.5, 'storage_stability_days': 90, 'ph_drift_per_day': 0.003,
             'reactivity_risk': 1, 'catalyzes_hydrolysis': False, 'catalyzes_ester_hydrolysis': False,
             'reference': 'Goldberg, R.N. et al. (2002). J. Phys. Chem. Ref. Data, 31(2), 231-370.'},
            
            # Ammonium Acetate — effective pH range 3.5-5.5 (acetate) and 8.5-10.5 (ammonium)
            # pKa listed here is for the acetate component (4.76); ammonium pKa is 9.25
            # Use Ammonium_Formate or Ammonium_Bicarbonate for near-neutral pH LC-MS work
            {'name': 'Ammonium_Acetate', 'base_name': 'Ammonium Acetate', 'pka': 4.76, 'pka_temp_coeff': -0.031,
             'uv_cutoff_10mM': 210, 'uv_cutoff_50mM': 220, 'solubility_ACN_10mM': 95,
             'solubility_MeOH_10mM': 98, 'solubility_IPA_10mM': 90, 'volatility': 'highly volatile',
             'ms_compatible': True, 'max_ms_conc': 20, 'metal_chelation_fe': 0.5,
             'metal_chelation_ca': 0.1, 'storage_stability_days': 7, 'ph_drift_per_day': 0.05,
             'reactivity_risk': 3, 'catalyzes_hydrolysis': False, 'catalyzes_ester_hydrolysis': False,
             'reacts_with_carbonyls': True,
             'reference': 'Kebarle, P. & Tang, L. (1993). Anal. Chem., 65(22), 972A-986A.'},
            
            # Ammonium Formate
            {'name': 'Ammonium_Formate', 'base_name': 'Ammonium Formate', 'pka': 3.75, 'pka_temp_coeff': 0.000,
             'uv_cutoff_10mM': 210, 'uv_cutoff_50mM': 215, 'solubility_ACN_10mM': 95,
             'solubility_MeOH_10mM': 98, 'solubility_IPA_10mM': 90, 'volatility': 'highly volatile',
             'ms_compatible': True, 'max_ms_conc': 20, 'metal_chelation_fe': 0.5,
             'metal_chelation_ca': 0.1, 'storage_stability_days': 7, 'ph_drift_per_day': 0.05,
             'reactivity_risk': 1, 'catalyzes_hydrolysis': False, 'catalyzes_ester_hydrolysis': False,
             'reacts_with_carbonyls': True,
             'reference': 'Kebarle, P. & Tang, L. (1993). Anal. Chem., 65(22), 972A-986A.'},
            
            # Citrate buffers
            {'name': 'Citrate_pKa1', 'base_name': 'Citrate', 'pka': 3.13, 'pka_temp_coeff': -0.002,
             'uv_cutoff_10mM': 230, 'uv_cutoff_50mM': 240, 'solubility_ACN_10mM': 60,
             'solubility_MeOH_10mM': 65, 'solubility_IPA_10mM': 55, 'volatility': 'non-volatile',
             'ms_compatible': False, 'max_ms_conc': 0, 'metal_chelation_fe': 11.4,
             'metal_chelation_ca': 3.5, 'storage_stability_days': 90, 'ph_drift_per_day': 0.002,
             'reactivity_risk': 5, 'chelates_metals': True,
             'reference': 'Martell, A.E. & Smith, R.M. (2004). Critical Stability Constants.'},
            
            {'name': 'Citrate_pKa2', 'base_name': 'Citrate', 'pka': 4.76, 'pka_temp_coeff': -0.002,
             'uv_cutoff_10mM': 230, 'uv_cutoff_50mM': 240, 'solubility_ACN_10mM': 60,
             'solubility_MeOH_10mM': 65, 'solubility_IPA_10mM': 55, 'volatility': 'non-volatile',
             'ms_compatible': False, 'max_ms_conc': 0, 'metal_chelation_fe': 11.4,
             'metal_chelation_ca': 3.5, 'storage_stability_days': 90, 'ph_drift_per_day': 0.002,
             'reactivity_risk': 5, 'chelates_metals': True,
             'reference': 'Martell, A.E. & Smith, R.M. (2004). Critical Stability Constants.'},
            
            {'name': 'Citrate_pKa3', 'base_name': 'Citrate', 'pka': 6.39, 'pka_temp_coeff': -0.002,
             'uv_cutoff_10mM': 230, 'uv_cutoff_50mM': 240, 'solubility_ACN_10mM': 60,
             'solubility_MeOH_10mM': 65, 'solubility_IPA_10mM': 55, 'volatility': 'non-volatile',
             'ms_compatible': False, 'max_ms_conc': 0, 'metal_chelation_fe': 11.4,
             'metal_chelation_ca': 3.5, 'storage_stability_days': 90, 'ph_drift_per_day': 0.002,
             'reactivity_risk': 5, 'chelates_metals': True,
             'reference': 'Martell, A.E. & Smith, R.M. (2004). Critical Stability Constants.'},
            
            # TRIS
            {'name': 'TRIS', 'base_name': 'TRIS', 'pka': 8.06, 'pka_temp_coeff': -0.031,
             'uv_cutoff_10mM': 220, 'uv_cutoff_50mM': 230, 'solubility_ACN_10mM': 85,
             'solubility_MeOH_10mM': 90, 'solubility_IPA_10mM': 80, 'volatility': 'moderate',
             'ms_compatible': True, 'max_ms_conc': 10, 'metal_chelation_fe': 2.0,
             'metal_chelation_ca': 1.0, 'storage_stability_days': 30, 'ph_drift_per_day': 0.01,
             'reactivity_risk': 5, 'reacts_with_carbonyls': True, 'primary_amine': True,
             'reference': 'Beynon, R.J. & Easterby, J.S. (1996). Buffer Solutions: The Basics.'},
            
            # MES
            {'name': 'MES', 'base_name': 'MES', 'pka': 6.10, 'pka_temp_coeff': -0.011,
             'uv_cutoff_10mM': 230, 'uv_cutoff_50mM': 240, 'solubility_ACN_10mM': 80,
             'solubility_MeOH_10mM': 85, 'solubility_IPA_10mM': 75, 'volatility': 'non-volatile',
             'ms_compatible': False, 'max_ms_conc': 0, 'metal_chelation_fe': 0.5,
             'metal_chelation_ca': 0.1, 'storage_stability_days': 90, 'ph_drift_per_day': 0.002,
             'reactivity_risk': 1,
             'reference': 'Goldberg, R.N. et al. (2002). J. Phys. Chem. Ref. Data, 31(2), 231-370.'},
            
            # HEPES
            {'name': 'HEPES', 'base_name': 'HEPES', 'pka': 7.48, 'pka_temp_coeff': -0.014,
             'uv_cutoff_10mM': 230, 'uv_cutoff_50mM': 240, 'solubility_ACN_10mM': 80,
             'solubility_MeOH_10mM': 85, 'solubility_IPA_10mM': 75, 'volatility': 'non-volatile',
             'ms_compatible': False, 'max_ms_conc': 0, 'metal_chelation_fe': 0.5,
             'metal_chelation_ca': 0.1, 'storage_stability_days': 90, 'ph_drift_per_day': 0.002,
             'reactivity_risk': 1,
             'reference': 'Goldberg, R.N. et al. (2002). J. Phys. Chem. Ref. Data, 31(2), 231-370.'},
            
            # Borate
            {'name': 'Borate', 'base_name': 'Borate', 'pka': 9.23, 'pka_temp_coeff': -0.008,
             'uv_cutoff_10mM': 190, 'uv_cutoff_50mM': 195, 'solubility_ACN_10mM': 40,
             'solubility_MeOH_10mM': 45, 'solubility_IPA_10mM': 35, 'volatility': 'non-volatile',
             'ms_compatible': False, 'max_ms_conc': 0, 'metal_chelation_fe': 1.0,
             'metal_chelation_ca': 0.5, 'storage_stability_days': 180, 'ph_drift_per_day': 0.001,
             'reactivity_risk': 5, 'complexes_diols': True,
             'reference': 'Martell, A.E. & Smith, R.M. (2004). Critical Stability Constants.'},
            
            # Bicarbonate
            {'name': 'Bicarbonate', 'base_name': 'Bicarbonate', 'pka': 6.33, 'pka_temp_coeff': -0.005,
             'uv_cutoff_10mM': 210, 'uv_cutoff_50mM': 220, 'solubility_ACN_10mM': 70,
             'solubility_MeOH_10mM': 75, 'solubility_IPA_10mM': 65, 'volatility': 'moderate',
             'ms_compatible': True, 'max_ms_conc': 20, 'metal_chelation_fe': 1.0,
             'metal_chelation_ca': 0.5, 'storage_stability_days': 1, 'ph_drift_per_day': 0.2,
             'reactivity_risk': 2,
             'reference': 'Sigma-Aldrich. (2020). Buffer Reference Center.'},
            
            # CAPS
            {'name': 'CAPS', 'base_name': 'CAPS', 'pka': 10.33, 'pka_temp_coeff': -0.018,
             'uv_cutoff_10mM': 220, 'uv_cutoff_50mM': 230, 'solubility_ACN_10mM': 80,
             'solubility_MeOH_10mM': 85, 'solubility_IPA_10mM': 75, 'volatility': 'non-volatile',
             'ms_compatible': False, 'max_ms_conc': 0, 'metal_chelation_fe': 0.5,
             'metal_chelation_ca': 0.1, 'storage_stability_days': 90, 'ph_drift_per_day': 0.002,
             'reactivity_risk': 1,
             'reference': 'Goldberg, R.N. et al. (2002). J. Phys. Chem. Ref. Data, 31(2), 231-370.'},
            
            # Formate
            {'name': 'Formate', 'base_name': 'Formate', 'pka': 3.75, 'pka_temp_coeff': 0.000,
             'uv_cutoff_10mM': 210, 'uv_cutoff_50mM': 215, 'solubility_ACN_10mM': 90,
             'solubility_MeOH_10mM': 95, 'solubility_IPA_10mM': 85, 'volatility': 'moderate',
             'ms_compatible': True, 'max_ms_conc': 20, 'metal_chelation_fe': 1.0,
             'metal_chelation_ca': 0.5, 'storage_stability_days': 90, 'ph_drift_per_day': 0.002,
             'reactivity_risk': 1,
             'reference': 'Goldberg, R.N. et al. (2002). J. Phys. Chem. Ref. Data, 31(2), 231-370.'},
        ])
        
        # Convert to DataFrame and set index
        self.buffers = pd.DataFrame(self.buffers)
        self.buffers.set_index('name', inplace=True)
    
    def initialize_reference_database(self):
        """Initialize reference database for tracking rule sources"""
        
        self.references = {
            'rule1': {
                'title': 'Buffer pKa Matching Rule',
                'primary_ref': 'Goldberg, R.N., Kishore, N., & Lennen, R.M. (2002). "Thermodynamic quantities for the ionization reactions of buffers." Journal of Physical and Chemical Reference Data, 31(2), 231-370.',
                'doi': '10.1063/1.1416902',
                'page': '231-370',
                'rule_text': 'Buffer is effective within ±1 pH unit of pKa; optimal within ±0.5 units'
            },
            'rule2': {
                'title': 'UV Transparency Rule',
                'primary_ref': 'Snyder, L.R., Kirkland, J.J., & Dolan, J.W. (2011). Introduction to Modern Liquid Chromatography (3rd ed.). John Wiley & Sons.',
                'chapter': '5',
                'section': '5.2.3',
                'rule_text': 'Buffer UV cutoff must be at least 10 nm below detection wavelength; absorbance <0.05 AU'
            },
            'rule3': {
                'title': 'Organic Solvent Solubility Rule',
                'primary_ref': 'Subirats, X., Bosch, E., & Rosés, M. (2012). "Stripping away the myths surrounding buffer selection in HPLC." LCGC Europe, 25(4), 192-201.',
                'table': '2',
                'rule_text': 'Buffer must remain soluble throughout gradient; HILIC requires ammonium acetate/formate'
            },
            'rule4': {
                'title': 'Buffer Capacity Rule',
                'primary_ref': 'Perrin, D.D., & Dempsey, B. (1974). Buffers for pH and Metal Ion Control. Chapman and Hall.',
                'chapter': '2',
                'table': '2.1',
                'rule_text': 'β ≥ 0.01 for analytical HPLC; β ≥ 0.05 for preparative HPLC; β_max = 0.576 × C'
            },
            'rule5': {
                'title': 'Volatility Rule for Mass Spectrometry',
                'primary_ref': 'Kebarle, P., & Tang, L. (1993). "From ions in solution to ions in the gas phase: The mechanism of electrospray mass spectrometry." Analytical Chemistry, 65(22), 972A-986A.',
                'doi': '10.1021/ac00070a001',
                'rule_text': 'Only volatile buffers (ammonium acetate/formate/bicarbonate) for LC-MS'
            },
            'rule6': {
                'title': 'Temperature Dependence Rule',
                'primary_ref': 'Beynon, R.J., & Easterby, J.S. (1996). Buffer Solutions: The Basics. Oxford University Press.',
                'chapter': '4',
                'table': '4.1',
                'rule_text': 'pKa(T) = pKa(25°C) + dpKa/dT × (T-25); TRIS has high temp coefficient (-0.031/°C)'
            },
            'rule7': {
                'title': 'Chemical Reactivity/Compatibility Rule',
                'primary_ref': "O'Brien, P.J., & Herschlag, D. (2001). 'Functional interrelationships in the alkaline phosphatase superfamily.' Biochemistry, 40(19), 5691-5699.",
                'doi': '10.1021/bi0028892',
                'rule_text': 'Buffers must not react with analyte or catalyze degradation'
            },
            'rule7_1': {
                'title': 'Amine Buffers with Carbonyls',
                'primary_ref': "O'Brien, P.J., & Herschlag, D. (2001). Biochemistry, 40(19), 5691-5699.",
                'rate_constant_TRIS': '2.3 × 10⁻³ M⁻¹s⁻¹',
                'rule_text': 'Avoid TRIS and ammonium with aldehydes/ketones (Schiff base formation)'
            },
            'rule7_2': {
                'title': 'Phosphate with Esters/Amides',
                'primary_ref': "O'Brien, P.J., & Herschlag, D. (2001). Biochemistry, 40(19), 5691-5699.",
                'enhancement_factor': '6.7×',
                'rule_text': 'Avoid phosphate at pH>7 and T>40°C with esters/amides'
            },
            'rule7_3': {
                'title': 'Citrate with Metal Ions',
                'primary_ref': 'Martell, A.E., & Smith, R.M. (2004). Critical Stability Constants. Plenum Press.',
                'logK_Fe3+': '11.4',
                'rule_text': 'Avoid citrate with metal-containing analytes (strong chelation)'
            },
            'rule7_4': {
                'title': 'Borate with Diols',
                'primary_ref': 'Martell, A.E., & Smith, R.M. (2004). Critical Stability Constants. Plenum Press.',
                'logK_glucose': '2.1',
                'rule_text': 'Avoid borate with vicinal diols (carbohydrates, catechols, nucleosides)'
            },
            'rule8': {
                'title': 'Ionic Strength Effect Rule',
                'primary_ref': 'Stahlberg, J. (1999). "Retention models for ions in chromatography." Journal of Chromatography A, 855(1), 3-55.',
                'doi': '10.1016/S0021-9673(99)00176-4',
                'rule_text': 'Optimal ionic strength 20-50 mM for analytical HPLC; I = ½∑cᵢzᵢ²'
            },
            'rule9': {
                'title': 'Buffer Storage Stability Rule',
                'primary_ref': 'Sigma-Aldrich. (2020). "Buffer Reference Center." Sigma-Aldrich Technical Bulletin.',
                'url': 'sigmaaldrich.com/technical-documents/articles/biology/buffer-reference-center.html',
                'rule_text': 'Maximum storage times vary: Phosphate (6mo), Acetate (3mo), Ammonium acetate (1wk)'
            },
            'rule10': {
                'title': 'Metal Complexation Rule',
                'primary_ref': 'Martell, A.E., & Smith, R.M. (2004). Critical Stability Constants (Vols. 1-6). Plenum Press.',
                'volume': '3',
                'rule_text': 'Citrate strongly chelates metals (log K Fe³⁺ = 11.4); phosphate precipitates metals'
            },
            'comprehensive_matrix': {
                'title': 'Comprehensive Buffer Selection Matrix',
                'primary_ref': 'Compilation of all 10 rules above',
                'rule_text': 'Final selection based on weighted scoring of all rules'
            }
        }
    
    def initialize_rule_weights(self):
        """Initialize weights for each rule in the scoring system"""
        
        # Default weights (can be customized by user)
        return {
            'rule1_pka_matching': 0.20,      # Critical for pH control
            'rule2_uv_transparency': 0.15,    # Important for detection
            'rule3_organic_solubility': 0.12,  # Important for gradient methods
            'rule4_buffer_capacity': 0.10,     # Standard requirement
            'rule5_volatility_ms': 0.15 if self.params.get('is_lcms', False) else 0.05,  # Critical for MS
            'rule6_temperature': 0.08,         # Important for high-temp methods
            'rule7_reactivity': 0.10,          # Critical for analyte stability
            'rule8_ionic_strength': 0.03,       # Minor effect
            'rule9_storage_stability': 0.02,    # Practical consideration
            'rule10_metal_complexation': 0.05   # Important for metalloproteins
        }
    
    def set_method_parameters(self, **kwargs):
        """Set method parameters for buffer selection"""
        
        self.params = {
            'target_ph': kwargs.get('target_ph', 7.0),
            'detection_wavelength_nm': kwargs.get('detection_wavelength_nm', 254),
            'is_lcms': kwargs.get('is_lcms', False),
            'organic_modifier': kwargs.get('organic_modifier', 'ACN'),  # ACN, MeOH, IPA
            'max_organic_percent': kwargs.get('max_organic_percent', 80),
            'temperature_c': kwargs.get('temperature_c', 25),
            'buffer_concentration_mM': kwargs.get('buffer_concentration_mM', 50),
            'is_hilic': kwargs.get('is_hilic', False),
            'is_preparative': kwargs.get('is_preparative', False),
            'contains_metals': kwargs.get('contains_metals', False),
            'contains_esters': kwargs.get('contains_esters', False),
            'contains_aldehydes': kwargs.get('contains_aldehydes', False),
            'contains_diols': kwargs.get('contains_diols', False),
            'storage_required_days': kwargs.get('storage_required_days', 7),
            'gradient_elution': kwargs.get('gradient_elution', False),
            'max_allowable_absorbance': kwargs.get('max_allowable_absorbance', 0.05),
        }
        
        # Update rule weights based on method
        self.rule_weights['rule5_volatility_ms'] = 0.15 if self.params['is_lcms'] else 0.05
        
        # Adjust buffer capacity requirement
        if self.params['is_preparative']:
            self.params['required_buffer_capacity'] = 0.05
        elif self.params['gradient_elution']:
            self.params['required_buffer_capacity'] = 0.02
        else:
            self.params['required_buffer_capacity'] = 0.01
    
    def calculate_rule1_score(self, buffer, mol_descriptors):
        """
        Rule 1: Buffer pKa Matching Rule
        Reference: Goldberg, R.N. et al. (2002). J. Phys. Chem. Ref. Data, 31(2), 231-370.
        """
        pka = buffer['pka']
        target_ph = self.params['target_ph']
        temp = self.params['temperature_c']
        
        # Adjust pKa for temperature (Rule 6 integration)
        temp_coeff = buffer.get('pka_temp_coeff', 0)
        pka_adjusted = pka + temp_coeff * (temp - 25)
        
        # Calculate pH difference
        ph_diff = abs(target_ph - pka_adjusted)
        
        # Score calculation (0-1 scale)
        if ph_diff <= 0.5:
            score = 1.0  # Optimal range
        elif ph_diff <= 1.0:
            score = 1.0 - (ph_diff - 0.5) * 2  # Linear decay to 0 at ±1
        else:
            score = 0.0  # Outside effective range
        
        # Reference tracking
        ref_used = self.references['rule1']
        
        return {
            'score': score,
            'pka_adjusted': pka_adjusted,
            'ph_diff': ph_diff,
            'reference': ref_used['primary_ref'],
            'rule_applied': 'Rule 1: Buffer pKa Matching (±0.5 optimal, ±1.0 effective)'
        }
    
    def calculate_rule2_score(self, buffer, mol_descriptors):
        """
        Rule 2: UV Transparency Rule
        Reference: Snyder, L.R., Kirkland, J.J., & Dolan, J.W. (2011). Introduction to Modern LC (3rd ed.). Chapter 5.
        """
        conc = self.params['buffer_concentration_mM']
        detection_wavelength = self.params['detection_wavelength_nm']
        
        # Get UV cutoff based on concentration
        if conc <= 10:
            uv_cutoff = buffer['uv_cutoff_10mM']
        else:
            # Estimate cutoff increase with concentration (5 nm per 50 mM for acetate)
            uv_cutoff = buffer['uv_cutoff_10mM'] + (conc / 50) * 5
        
        # Calculate if detection wavelength is safe
        wavelength_margin = detection_wavelength - uv_cutoff
        
        if wavelength_margin >= 20:
            score = 1.0  # Excellent margin
        elif wavelength_margin >= 10:
            score = 0.8  # Adequate margin
        elif wavelength_margin >= 0:
            score = 0.5  # Marginal
        else:
            score = 0.0  # Unsafe
        
        # Estimate absorbance at detection wavelength
        # Simplified model: absorbance drops exponentially after cutoff
        if detection_wavelength > uv_cutoff:
            estimated_absorbance = 0.1 * np.exp(-0.1 * (detection_wavelength - uv_cutoff))
        else:
            estimated_absorbance = 1.0  # Very high
        
        # Check absorbance limit
        absorbance_ok = estimated_absorbance < self.params['max_allowable_absorbance']
        if not absorbance_ok:
            score = min(score, 0.3)
        
        ref_used = self.references['rule2']
        
        return {
            'score': score,
            'uv_cutoff_nm': uv_cutoff,
            'wavelength_margin': wavelength_margin,
            'estimated_absorbance': estimated_absorbance,
            'absorbance_ok': absorbance_ok,
            'reference': ref_used['primary_ref'],
            'rule_applied': f'Rule 2: UV Transparency (cutoff {uv_cutoff:.0f} nm, detection at {detection_wavelength} nm)'
        }
    
    def calculate_rule3_score(self, buffer, mol_descriptors):
        """
        Rule 3: Organic Solvent Solubility Rule
        Reference: Subirats, X., Bosch, E., & Rosés, M. (2012). LCGC Europe, 25(4), 192-201.
        """
        organic = self.params['organic_modifier']
        max_organic_needed = self.params['max_organic_percent']
        is_hilic = self.params['is_hilic']
        
        # Get solubility limit for the organic modifier
        if organic == 'ACN':
            solubility_limit = buffer['solubility_ACN_10mM']
        elif organic == 'MeOH':
            solubility_limit = buffer['solubility_MeOH_10mM']
        elif organic == 'IPA':
            solubility_limit = buffer['solubility_IPA_10mM']
        else:
            solubility_limit = 50  # Default conservative value
        
        # Adjust for concentration (higher conc = lower solubility)
        conc_factor = 50 / self.params['buffer_concentration_mM']
        solubility_limit_adjusted = solubility_limit * min(conc_factor, 1.0)
        
        # HILIC special rule: >70% organic requires ammonium acetate/formate
        if is_hilic and max_organic_needed > 70:
            if buffer['base_name'] not in ['Ammonium Acetate', 'Ammonium Formate']:
                return {
                    'score': 0.0,
                    'solubility_limit': solubility_limit_adjusted,
                    'organic_needed': max_organic_needed,
                    'hilic_violation': True,
                    'reference': self.references['rule3']['primary_ref'],
                    'rule_applied': 'Rule 3 (HILIC): >70% organic requires ammonium acetate/formate'
                }
        
        # Calculate score based on solubility margin
        solubility_margin = solubility_limit_adjusted - max_organic_needed
        
        if solubility_margin >= 20:
            score = 1.0
        elif solubility_margin >= 10:
            score = 0.8
        elif solubility_margin >= 0:
            score = 0.5
        else:
            score = 0.0  # Will precipitate
        
        ref_used = self.references['rule3']
        
        return {
            'score': score,
            'solubility_limit': solubility_limit_adjusted,
            'solubility_margin': solubility_margin,
            'organic_needed': max_organic_needed,
            'hilic_violation': False,
            'reference': ref_used['primary_ref'],
            'rule_applied': f'Rule 3: Organic Solubility (max {solubility_limit_adjusted:.0f}% {organic})'
        }
    
    def calculate_rule4_score(self, buffer, mol_descriptors):
        """
        Rule 4: Buffer Capacity Rule
        Reference: Perrin, D.D., & Dempsey, B. (1974). Buffers for pH and Metal Ion Control. Chapman and Hall.
        """
        conc = self.params['buffer_concentration_mM'] / 1000  # Convert to M
        pka = buffer['pka']
        target_ph = self.params['target_ph']
        required_beta = self.params['required_buffer_capacity']
        
        # Calculate buffer capacity at target pH
        # β = 2.303 × C × (Ka[H+] / (Ka + [H+])^2)
        Ka = 10 ** -pka
        H = 10 ** -target_ph
        
        beta = 2.303 * conc * (Ka * H) / ((Ka + H) ** 2)
        
        # Maximum buffer capacity (at pH = pKa)
        beta_max = 0.576 * conc
        
        # Score based on meeting requirement
        if beta >= required_beta:
            # Scale score based on how well it meets requirement
            score = min(1.0, beta / (required_beta * 2))
        else:
            score = beta / required_beta * 0.5  # Less than requirement
        
        ref_used = self.references['rule4']
        
        return {
            'score': score,
            'beta': beta,
            'beta_max': beta_max,
            'required_beta': required_beta,
            'meets_requirement': beta >= required_beta,
            'reference': ref_used['primary_ref'],
            'rule_applied': f'Rule 4: Buffer Capacity (β = {beta:.4f}, required {required_beta})'
        }
    
    def calculate_rule5_score(self, buffer, mol_descriptors):
        """
        Rule 5: Volatility Rule for Mass Spectrometry
        Reference: Kebarle, P., & Tang, L. (1993). Anal. Chem., 65(22), 972A-986A.
        """
        is_lcms = self.params['is_lcms']
        
        if not is_lcms:
            # Not using MS, so this rule doesn't apply
            return {
                'score': 1.0,  # No penalty
                'ms_compatible': True,
                'reference': self.references['rule5']['primary_ref'],
                'rule_applied': 'Rule 5: Volatility (not applicable - not LC-MS)'
            }
        
        # Check MS compatibility
        ms_compatible = buffer['ms_compatible']
        max_conc = buffer['max_ms_conc']
        current_conc = self.params['buffer_concentration_mM']
        
        if not ms_compatible:
            score = 0.0
        elif current_conc <= max_conc:
            score = 1.0
        else:
            # Concentration too high, score based on how much over limit
            score = max(0, 1.0 - (current_conc - max_conc) / max_conc)
        
        # ESI signal suppression factors (from reference)
        if buffer['base_name'] == 'Ammonium Acetate':
            if current_conc <= 5:
                signal_factor = 0.95
            elif current_conc <= 20:
                signal_factor = 0.85
            else:
                signal_factor = 0.70
        elif buffer['base_name'] == 'Ammonium Formate':
            if current_conc <= 5:
                signal_factor = 0.96
            elif current_conc <= 20:
                signal_factor = 0.87
            else:
                signal_factor = 0.75
        elif buffer['base_name'] == 'TRIS':
            signal_factor = 0.70 if current_conc <= 5 else 0.40
        else:
            signal_factor = 1.0
        
        ref_used = self.references['rule5']
        
        return {
            'score': score,
            'ms_compatible': ms_compatible,
            'max_ms_conc': max_conc,
            'signal_factor': signal_factor,
            'reference': ref_used['primary_ref'],
            'rule_applied': f'Rule 5: MS Volatility (compatible: {ms_compatible}, signal factor: {signal_factor})'
        }
    
    def calculate_rule6_score(self, buffer, mol_descriptors):
        """
        Rule 6: Temperature Dependence Rule
        Reference: Beynon, R.J., & Easterby, J.S. (1996). Buffer Solutions: The Basics. Oxford University Press.
        """
        temp = self.params['temperature_c']
        pka = buffer['pka']
        temp_coeff = buffer['pka_temp_coeff']
        
        # Calculate pKa change with temperature
        pka_change = temp_coeff * (temp - 25)
        pka_at_temp = pka + pka_change
        
        # Calculate pH shift relative to target
        ph_shift = abs(pka_at_temp - pka)
        
        # Score based on temperature stability
        if temp <= 35:
            score = 1.0  # Ambient temperature, all buffers ok
        elif temp <= 50:
            # Moderate temperature, penalize high temp-coeff buffers
            if abs(temp_coeff) > 0.02:  # TRIS, Ammonium
                score = 0.3
            elif abs(temp_coeff) > 0.01:  # HEPES, MES
                score = 0.7
            else:
                score = 1.0
        elif temp <= 80:
            # High temperature, only phosphate or acetate preferred
            if buffer['base_name'] in ['Phosphate', 'Acetate']:
                score = 0.9
            else:
                score = 0.2
        else:
            # Very high temperature >80°C, phosphate only
            if buffer['base_name'] == 'Phosphate':
                score = 0.8
            else:
                score = 0.0
        
        ref_used = self.references['rule6']
        
        return {
            'score': score,
            'pka_at_temp': pka_at_temp,
            'pka_change': pka_change,
            'temp_coeff': temp_coeff,
            'reference': ref_used['primary_ref'],
            'rule_applied': f'Rule 6: Temperature Dependence (ΔpKa = {pka_change:.3f} at {temp}°C)'
        }
    
    def calculate_rule7_score(self, buffer, mol_descriptors):
        """
        Rule 7: Chemical Reactivity/Compatibility Rule
        Reference: O'Brien, P.J., & Herschlag, D. (2001). Biochemistry, 40(19), 5691-5699.
        Multiple sub-rules for specific incompatibilities
        """
        score = 1.0
        reasons = []
        
        # Rule 7.1: Amine buffers with carbonyls
        if buffer.get('reacts_with_carbonyls', False) and (
            mol_descriptors.get('has_aldehyde', False) or 
            mol_descriptors.get('has_ketone', False)
        ):
            score *= 0.0
            reasons.append('Amine buffer with carbonyl analyte (Schiff base formation)')
            ref_used_7_1 = self.references['rule7_1']
        
        # Rule 7.2: Phosphate with esters/amides at high pH/temp
        if buffer['base_name'] == 'Phosphate' and buffer.get('catalyzes_ester_hydrolysis', False):
            if (self.params['target_ph'] > 7 and self.params['temperature_c'] > 40 and
                (mol_descriptors.get('has_ester', False) or mol_descriptors.get('has_amide', False))):
                score *= 0.1
                reasons.append('Phosphate catalyzes ester/amide hydrolysis at pH>7, T>40°C')
                ref_used_7_2 = self.references['rule7_2']
        
        # Rule 7.3: Citrate with metal ions
        if buffer.get('chelates_metals', False) and (
            self.params['contains_metals'] or mol_descriptors.get('chelation_potential', 0) > 3
        ):
            score *= 0.0
            reasons.append('Citrate strongly chelates metal ions (log K Fe³⁺ = 11.4)')
            ref_used_7_3 = self.references['rule7_3']
        
        # Rule 7.4: Borate with diols
        if buffer.get('complexes_diols', False) and (
            mol_descriptors.get('has_vicinal_diol', False) or
            mol_descriptors.get('has_catechol', False) or
            mol_descriptors.get('has_nucleoside', False)
        ):
            score *= 0.0
            reasons.append('Borate forms complexes with vicinal diols')
            ref_used_7_4 = self.references['rule7_4']
        
        # General reactivity risk
        reactivity_risk = buffer.get('reactivity_risk', 1)
        if reactivity_risk > 3 and score > 0:
            score *= (5 - reactivity_risk) / 2  # Reduce score for high-risk buffers
        
        ref_used = self.references['rule7']
        
        return {
            'score': max(0, score),
            'reasons': reasons,
            'reactivity_risk': reactivity_risk,
            'reference': ref_used['primary_ref'],
            'rule_applied': f'Rule 7: Chemical Reactivity (risk: {reactivity_risk})'
        }
    
    def calculate_rule8_score(self, buffer, mol_descriptors):
        """
        Rule 8: Ionic Strength Effect Rule
        Reference: Stahlberg, J. (1999). J. Chromatogr. A, 855(1), 3-55.
        """
        conc = self.params['buffer_concentration_mM'] / 1000  # Convert to M
        pka = buffer['pka']
        target_ph = self.params['target_ph']
        
        # Estimate ionic strength based on buffer type and pH
        # Simplified model: I ≈ conc * (degree of ionization)
        if target_ph <= pka - 1:
            # Mostly acid form, low charge
            ionic_strength = conc * 0.1
        elif target_ph >= pka + 1:
            # Mostly base form, higher charge
            if buffer['base_name'] in ['Phosphate', 'Citrate']:
                ionic_strength = conc * 2  # Multi-charged
            else:
                ionic_strength = conc * 1  # Single charged
        else:
            # Near pKa, mixture
            if buffer['base_name'] in ['Phosphate', 'Citrate']:
                ionic_strength = conc * 1.5
            else:
                ionic_strength = conc * 0.5
        
        # Optimal ionic strength range
        if self.params['is_preparative']:
            optimal_range = (0.05, 0.15)
        elif self.params['gradient_elution']:
            optimal_range = (0.02, 0.05)
        else:
            optimal_range = (0.02, 0.05)
        
        # Score based on being in optimal range
        if optimal_range[0] <= ionic_strength <= optimal_range[1]:
            score = 1.0
        elif ionic_strength < optimal_range[0]:
            score = ionic_strength / optimal_range[0]  # Below optimal
        else:
            score = optimal_range[1] / ionic_strength  # Above optimal
        
        ref_used = self.references['rule8']
        
        return {
            'score': min(1.0, score),
            'ionic_strength': ionic_strength,
            'optimal_range': optimal_range,
            'reference': ref_used['primary_ref'],
            'rule_applied': f'Rule 8: Ionic Strength (I = {ionic_strength:.3f} M)'
        }
    
    def calculate_rule9_score(self, buffer, mol_descriptors):
        """
        Rule 9: Buffer Storage Stability Rule
        Reference: Sigma-Aldrich. (2020). Buffer Reference Center.
        """
        storage_days_needed = self.params['storage_required_days']
        max_storage_days = buffer['storage_stability_days']
        ph_drift_rate = buffer['ph_drift_per_day']
        
        # Check if buffer can be stored for required time
        if max_storage_days >= storage_days_needed:
            storage_score = 1.0
        else:
            storage_score = max_storage_days / storage_days_needed
        
        # Calculate pH drift over storage period
        ph_drift = ph_drift_rate * storage_days_needed
        
        # pH drift should be < 0.1 pH units for stability
        if ph_drift < 0.1:
            drift_score = 1.0
        elif ph_drift < 0.3:
            drift_score = 0.7
        elif ph_drift < 0.5:
            drift_score = 0.4
        else:
            drift_score = 0.1
        
        # Combined score
        score = (storage_score + drift_score) / 2
        
        ref_used = self.references['rule9']
        
        return {
            'score': score,
            'max_storage_days': max_storage_days,
            'ph_drift_over_period': ph_drift,
            'storage_adequate': max_storage_days >= storage_days_needed,
            'reference': ref_used['primary_ref'],
            'rule_applied': f'Rule 9: Storage Stability (max {max_storage_days} days, drift {ph_drift:.3f} pH)'
        }
    
    def calculate_rule10_score(self, buffer, mol_descriptors):
        """
        Rule 10: Metal Complexation Rule
        Reference: Martell, A.E., & Smith, R.M. (2004). Critical Stability Constants.
        """
        contains_metals = self.params['contains_metals']
        chelation_potential = mol_descriptors.get('chelation_potential', 0)
        
        # Metal chelation constants
        metal_chelation_fe = buffer['metal_chelation_fe']
        metal_chelation_ca = buffer['metal_chelation_ca']
        
        # Calculate risk based on chelation strength
        if contains_metals or chelation_potential > 2:
            # Higher risk for strong chelators
            if metal_chelation_fe > 10:  # Citrate level
                score = 0.0
            elif metal_chelation_fe > 5:
                score = 0.3
            elif metal_chelation_fe > 2:
                score = 0.7
            else:
                score = 1.0
        else:
            score = 1.0  # No metals, no penalty
        
        # Corrosion risk to HPLC system
        corrosion_risk = buffer.get('corrosion_risk', 'LOW')
        if corrosion_risk == 'HIGH':
            score *= 0.5
        
        ref_used = self.references['rule10']
        
        return {
            'score': score,
            'metal_chelation_fe': metal_chelation_fe,
            'metal_chelation_ca': metal_chelation_ca,
            'corrosion_risk': corrosion_risk,
            'reference': ref_used['primary_ref'],
            'rule_applied': f'Rule 10: Metal Complexation (log K Fe³⁺ = {metal_chelation_fe})'
        }
    
    def select_optimal_buffer(self, smiles):
        """
        Main function to select optimal buffer based on all 10 rules
        """
        # Calculate molecular descriptors
        descriptor_calc = MolecularDescriptorCalculator(smiles)
        mol_descriptors = descriptor_calc.calculate_all_descriptors()
        
        # Store results for all buffers
        results = []
        
        for buffer_name, buffer in self.buffers.iterrows():
            buffer_scores = {}
            total_score = 0
            total_weight = 0
            
            # Calculate score for each rule
            rule_functions = [
                ('rule1', self.calculate_rule1_score),
                ('rule2', self.calculate_rule2_score),
                ('rule3', self.calculate_rule3_score),
                ('rule4', self.calculate_rule4_score),
                ('rule5', self.calculate_rule5_score),
                ('rule6', self.calculate_rule6_score),
                ('rule7', self.calculate_rule7_score),
                ('rule8', self.calculate_rule8_score),
                ('rule9', self.calculate_rule9_score),
                ('rule10', self.calculate_rule10_score)
            ]
            
            for rule_name, rule_func in rule_functions:
                rule_result = rule_func(buffer, mol_descriptors)
                buffer_scores[rule_name] = rule_result
                
                # Apply weight — map rule name to weight dict key
                _wkey_map = {
                    'rule1': 'rule1_pka_matching',
                    'rule2': 'rule2_uv_transparency',
                    'rule3': 'rule3_organic_solubility',
                    'rule4': 'rule4_buffer_capacity',
                    'rule5': 'rule5_volatility_ms',
                    'rule6': 'rule6_temperature',
                    'rule7': 'rule7_reactivity',
                    'rule8': 'rule8_ionic_strength',
                    'rule9': 'rule9_storage_stability',
                    'rule10': 'rule10_metal_complexation',
                }
                weight = self.rule_weights.get(_wkey_map.get(rule_name, ''), 0.1)
                
                total_score += rule_result['score'] * weight
                total_weight += weight
            
            # Normalize score
            final_score = total_score / total_weight if total_weight > 0 else 0
            
            results.append({
                'buffer_name': buffer_name,
                'base_name': buffer['base_name'],
                'pka': buffer['pka'],
                'final_score': final_score,
                'scores': buffer_scores,
                'compatibility_notes': self._generate_compatibility_notes(buffer_scores, mol_descriptors)
            })
        
        # Sort by final score
        results.sort(key=lambda x: x['final_score'], reverse=True)
        
        # Compile reference list for Excel export
        reference_list = self._compile_references(results)
        
        return {
            'top_buffers': results[:5],  # Top 5 recommendations
            'all_buffers': results,
            'molecular_descriptors': mol_descriptors,
            'method_parameters': self.params,
            'references': reference_list
        }
    
    def _generate_compatibility_notes(self, scores, mol_descriptors):
        """Generate human-readable compatibility notes"""
        notes = []
        
        # Check critical failures
        for rule_name, rule_result in scores.items():
            if rule_result['score'] == 0:
                if 'reasons' in rule_result and rule_result['reasons']:
                    notes.extend(rule_result['reasons'])
                elif 'rule_applied' in rule_result:
                    notes.append(f"Failed: {rule_result['rule_applied']}")
        
        return notes
    
    def _compile_references(self, results):
        """Compile all references used in the analysis for Excel export"""
        references_used = set()
        
        for result in results:
            for rule_name, rule_result in result['scores'].items():
                if 'reference' in rule_result:
                    references_used.add(rule_result['reference'])
        
        # Convert to list for DataFrame
        ref_list = []
        for ref in references_used:
            ref_list.append({
                'reference': ref,
                'used_in_rules': ', '.join([r for r in self.references.keys() if self.references[r]['primary_ref'] == ref])
            })
        
        return ref_list


class ExcelOutputGenerator:
    """Generate Excel output with all results and references"""
    
    def __init__(self, output_filename=None):
        if output_filename is None:
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            output_filename = f'buffer_selection_{timestamp}.xlsx'
        self.output_filename = output_filename
        
    def generate_excel(self, results):
        """Generate comprehensive Excel report"""
        
        with pd.ExcelWriter(self.output_filename, engine='openpyxl') as writer:
            
            # Sheet 1: Top Recommendations
            self._create_top_recommendations_sheet(writer, results)
            
            # Sheet 2: All Buffers Detailed Scores
            self._create_all_buffers_sheet(writer, results)
            
            # Sheet 3: Molecular Descriptors
            self._create_descriptors_sheet(writer, results)
            
            # Sheet 4: Method Parameters
            self._create_parameters_sheet(writer, results)
            
            # Sheet 5: References
            self._create_references_sheet(writer, results)
            
            # Sheet 6: Rule-by-Rule Analysis
            self._create_rule_analysis_sheet(writer, results)
            
            # Sheet 7: Incompatibility Warnings
            self._create_warnings_sheet(writer, results)
            
        logger.info(f"Excel report generated: {self.output_filename}")
        return self.output_filename
    
    def _create_top_recommendations_sheet(self, writer, results):
        """Create sheet with top 5 recommendations"""
        top_buffers = results['top_buffers']
        
        data = []
        for i, buffer in enumerate(top_buffers, 1):
            data.append({
                'Rank': i,
                'Buffer': buffer['base_name'],
                'pKa': buffer['pka'],
                'Overall Score': f"{buffer['final_score']:.3f}",
                'Compatibility Notes': '; '.join(buffer['compatibility_notes'][:3]) if buffer['compatibility_notes'] else 'Compatible'
            })
        
        df = pd.DataFrame(data)
        df.to_excel(writer, sheet_name='Top Recommendations', index=False)
    
    def _create_all_buffers_sheet(self, writer, results):
        """Create sheet with all buffers and detailed scores"""
        all_buffers = results['all_buffers']
        
        data = []
        for buffer in all_buffers:
            row = {
                'Buffer': buffer['base_name'],
                'pKa': buffer['pka'],
                'Overall Score': f"{buffer['final_score']:.3f}",
            }
            
            # Add individual rule scores
            for rule_name, rule_result in buffer['scores'].items():
                row[f'{rule_name.upper()} Score'] = f"{rule_result['score']:.3f}"
            
            data.append(row)
        
        df = pd.DataFrame(data)
        df.to_excel(writer, sheet_name='All Buffers', index=False)
    
    def _create_descriptors_sheet(self, writer, results):
        """Create sheet with molecular descriptors"""
        descriptors = results['molecular_descriptors']
        
        # Flatten dictionary for DataFrame
        flat_descriptors = {}
        for key, value in descriptors.items():
            if isinstance(value, (int, float, str, bool)):
                flat_descriptors[key] = value
            elif isinstance(value, list):
                flat_descriptors[key] = ', '.join(map(str, value))
            else:
                flat_descriptors[key] = str(value)
        
        df = pd.DataFrame([flat_descriptors])
        df.to_excel(writer, sheet_name='Molecular Descriptors', index=False)
    
    def _create_parameters_sheet(self, writer, results):
        """Create sheet with method parameters"""
        params = results['method_parameters']
        
        df = pd.DataFrame([params])
        df.to_excel(writer, sheet_name='Method Parameters', index=False)
    
    def _create_references_sheet(self, writer, results):
        """Create sheet with all references"""
        references = results['references']
        
        df = pd.DataFrame(references)
        df.to_excel(writer, sheet_name='References', index=False)
    
    def _create_rule_analysis_sheet(self, writer, results):
        """Create sheet with detailed rule-by-rule analysis"""
        top_buffers = results['top_buffers']
        
        # Create rule descriptions
        rule_descriptions = {
            'rule1': 'Buffer pKa Matching (±0.5 optimal)',
            'rule2': 'UV Transparency',
            'rule3': 'Organic Solvent Solubility',
            'rule4': 'Buffer Capacity',
            'rule5': 'MS Volatility',
            'rule6': 'Temperature Dependence',
            'rule7': 'Chemical Reactivity',
            'rule8': 'Ionic Strength',
            'rule9': 'Storage Stability',
            'rule10': 'Metal Complexation'
        }
        
        data = []
        for buffer in top_buffers:
            for rule_name, rule_desc in rule_descriptions.items():
                if rule_name in buffer['scores']:
                    rule_result = buffer['scores'][rule_name]
                    data.append({
                        'Buffer': buffer['base_name'],
                        'Rule': rule_desc,
                        'Score': f"{rule_result['score']:.3f}",
                        'Details': rule_result.get('rule_applied', ''),
                        'Reference': rule_result.get('reference', '')[:100] + '...'  # Truncate for readability
                    })
        
        df = pd.DataFrame(data)
        df.to_excel(writer, sheet_name='Rule Analysis', index=False)
    
    def _create_warnings_sheet(self, writer, results):
        """Create sheet with incompatibility warnings"""
        all_buffers = results['all_buffers']
        
        data = []
        for buffer in all_buffers:
            if buffer['compatibility_notes']:
                for note in buffer['compatibility_notes']:
                    data.append({
                        'Buffer': buffer['base_name'],
                        'Warning': note,
                        'Overall Score': f"{buffer['final_score']:.3f}"
                    })
        
        if data:
            df = pd.DataFrame(data)
            df.to_excel(writer, sheet_name='Warnings', index=False)


# Main execution class
class BufferSelector:
    """
    Main class for buffer selection based on SMILES input
    """
    
    def __init__(self):
        self.engine = BufferSelectionEngine()
        self.excel_generator = ExcelOutputGenerator()
        
    def select_buffer(self, smiles, method_params=None, output_excel=True):
        """
        Select optimal buffer for given SMILES
        
        Parameters:
        -----------
        smiles : str
            SMILES string of the analyte
        method_params : dict
            Method parameters (target_ph, detection_wavelength, etc.)
        output_excel : bool
            Whether to generate Excel output
        
        Returns:
        --------
        dict with results and optionally Excel filename
        """
        
        # Set default method parameters if not provided
        if method_params is None:
            method_params = {}
        
        default_params = {
            'target_ph': 7.0,
            'detection_wavelength_nm': 254,
            'is_lcms': False,
            'organic_modifier': 'ACN',
            'max_organic_percent': 80,
            'temperature_c': 25,
            'buffer_concentration_mM': 50,
            'is_hilic': False,
            'is_preparative': False,
            'contains_metals': False,
            'contains_esters': False,
            'contains_aldehydes': False,
            'contains_diols': False,
            'storage_required_days': 7,
            'gradient_elution': False,
            'max_allowable_absorbance': 0.05
        }
        
        # Update with user parameters
        default_params.update(method_params)
        method_params = default_params
        
        # Set method parameters in engine
        self.engine.set_method_parameters(**method_params)
        
        # Run buffer selection
        results = self.engine.select_optimal_buffer(smiles)
        
        # Generate Excel output if requested
        if output_excel:
            excel_file = self.excel_generator.generate_excel(results)
            results['excel_file'] = excel_file
        
        # Print summary
        self._print_summary(results)
        
        return results
    
    def _print_summary(self, results):
        """Print summary of results to console"""
        print("\n" + "="*80)
        print("BUFFER SELECTION RESULTS")
        print("="*80)
        
        print("\nTop 5 Recommended Buffers:")
        print("-" * 60)
        for i, buffer in enumerate(results['top_buffers'], 1):
            print(f"{i}. {buffer['base_name']} (pKa = {buffer['pka']})")
            print(f"   Score: {buffer['final_score']:.3f}")
            if buffer['compatibility_notes']:
                print(f"   Note: {buffer['compatibility_notes'][0]}")
            print()
        
        print("\nMethod Parameters:")
        print("-" * 60)
        for key, value in results['method_parameters'].items():
            print(f"  {key}: {value}")
        
        print("\n" + "="*80)


# Example usage
if __name__ == "__main__":
    # Initialize selector
    selector = BufferSelector()
    
    # Example 1: Simple carboxylic acid (aspirin)
    print("\n" + "="*80)
    print("EXAMPLE 1: Aspirin (acetylsalicylic acid)")
    print("="*80)
    
    aspirin_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
    
    # Method parameters for aspirin
    aspirin_params = {
        'target_ph': 3.5,  # Acidic pH for carboxylic acid
        'detection_wavelength_nm': 230,  # UV absorption of aspirin
        'is_lcms': False,
        'temperature_c': 30,
        'contains_esters': True,  # Aspirin has ester bond
        'buffer_concentration_mM': 50
    }
    
    results1 = selector.select_buffer(aspirin_smiles, aspirin_params, output_excel=True)
    
    # Example 2: Peptide with metal binding
    print("\n" + "="*80)
    print("EXAMPLE 2: Peptide with histidine (metal binding)")
    print("="*80)
    
    peptide_smiles = "C1=CN=C(N1)CC(C(=O)O)N"  # Simplified histidine
    
    peptide_params = {
        'target_ph': 7.4,  # Physiological pH
        'detection_wavelength_nm': 280,  # Peptide bond
        'is_lcms': True,  # LC-MS analysis
        'contains_metals': True,  # Metalloprotein
        'temperature_c': 37,  # Physiological temp
        'buffer_concentration_mM': 20,  # Lower for MS
        'is_hilic': False
    }
    
    results2 = selector.select_buffer(peptide_smiles, peptide_params, output_excel=True)
    
    # Example 3: Carbohydrate (diol-containing)
    print("\n" + "="*80)
    print("EXAMPLE 3: Glucose (carbohydrate with vicinal diols)")
    print("="*80)
    
    glucose_smiles = "C(C1C(C(C(C(O1)O)O)O)O)O"
    
    glucose_params = {
        'target_ph': 7.0,
        'detection_wavelength_nm': 195,  # Low UV
        'is_lcms': False,
        'contains_diols': True,  # Vicinal diols
        'is_hilic': True,  # HILIC mode for carbohydrates
        'max_organic_percent': 85,  # High organic in HILIC
        'organic_modifier': 'ACN'
    }
    
    results3 = selector.select_buffer(glucose_smiles, glucose_params, output_excel=True)