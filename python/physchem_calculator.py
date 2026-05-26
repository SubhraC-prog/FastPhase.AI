"""
Physicochemical Property Calculator Module

This module provides comprehensive physicochemical property calculations using RDKit
with specialized chromatographic property predictions. It serves as the core calculation
engine for the AI-Assisted Chromatographic Method Development System.

Author: Chromatography AI Team
Version: 1.0.0
"""

import logging
import numpy as np
from datetime import datetime

# Try to import RDKit. If unavailable, set USE_RDKIT=False and fall back
try:
    from rdkit import Chem
    from rdkit.Chem import (
        Descriptors,
        Crippen,
        Lipinski,
        rdMolDescriptors,
        GraphDescriptors,
        MolSurf,
        EState,
        Fragments,
        QED,
    )
    from rdkit.Chem.rdMolDescriptors import (
        CalcNumHBA,
        CalcNumHBD,
        CalcNumRotatableBonds,
        CalcNumAromaticRings,
        CalcNumAliphaticRings,
        CalcNumSaturatedRings,
        CalcNumHeterocycles,
        CalcNumSpiroAtoms,
        CalcNumBridgeheadAtoms,
        CalcNumAmideBonds,
        CalcNumHeavyAtoms,
        CalcNumAtoms,
        CalcNumHeteroatoms,
    )
    from rdkit.Chem.EState import EStateIndices
    from rdkit.Chem.Lipinski import NumRotatableBonds
    from rdkit.Chem.rdMolDescriptors import CalcTPSA
    USE_RDKIT = True
except Exception:
    # RDKit not available in this environment — we'll provide a lightweight fallback
    USE_RDKIT = False
    # Provide minimal placeholders so module can be imported cleanly
    Chem = None
    Descriptors = Crippen = Lipinski = rdMolDescriptors = GraphDescriptors = MolSurf = EState = Fragments = QED = None
    CalcNumHBA = CalcNumHBD = CalcNumRotatableBonds = CalcNumAromaticRings = CalcNumAliphaticRings = None
    CalcNumSaturatedRings = CalcNumHeterocycles = CalcNumSpiroAtoms = CalcNumBridgeheadAtoms = None
    CalcNumAmideBonds = CalcNumHeavyAtoms = CalcNumAtoms = CalcNumHeteroatoms = None
    EStateIndices = lambda mol: []
    NumRotatableBonds = lambda mol: 0
    CalcTPSA = lambda mol: 0.0
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass, field, asdict
from enum import Enum
import json
from pathlib import Path

# Configure logging
logger = logging.getLogger(__name__)

# Some RDKit installations may not expose CalcFractionCsp3 at the expected
# location or may have API differences across versions. Provide a safe
# fallback implementation that computes fraction_csp3 from atom hybridization
# if the symbol is not available. If RDKit is not available, provide a
# lightweight string-based fallback used by the minimal calculator.
try:
    if USE_RDKIT:
        from rdkit.Chem.rdMolDescriptors import CalcFractionCsp3  # type: ignore
    else:
        raise ImportError()
except Exception:
    def CalcFractionCsp3(mol_or_smiles):
        # If given an RDKit Mol, use atom hybridization; otherwise parse SMILES
        try:
            atoms = []
            if USE_RDKIT and mol_or_smiles is not None:
                mol = mol_or_smiles
                total_c = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
                if total_c == 0:
                    return 0.0
                sp3_c = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6 and a.GetHybridization() == Chem.HybridizationType.SP3)
                return float(sp3_c) / float(total_c)
        except Exception:
            # Fallback: approximate from SMILES string if provided
            try:
                smiles = str(mol_or_smiles)
                import re
                carbons = len(re.findall(r'C(?![a-z])', smiles)) + len(re.findall(r'c', smiles))
                if carbons == 0:
                    return 0.0
                # crude heuristic: assume 50% sp3
                return 0.5
            except Exception:
                return 0.0


class IonizationType(Enum):
    """Enum for ionization types"""
    ACIDIC = "acidic"
    BASIC = "basic"
    NEUTRAL = "neutral"
    AMPHOTERIC = "amphoteric"
    ZWITTERIONIC = "zwitterionic"


@dataclass
class IonizableGroup:
    """Data class for ionizable group information"""
    group_type: str
    pka: float
    count: int
    position: str
    confidence: float


@dataclass
class PhysicochemicalProperties:
    """
    Comprehensive data class for all physicochemical properties.
    
    This class stores all calculated properties including:
    - Basic molecular descriptors
    - Lipophilicity parameters
    - Hydrogen bonding characteristics
    - Topological and geometrical descriptors
    - Electronic and charge properties
    - Solubility and permeability predictions
    - Chromatographic behavior predictors
    """
    
    # Basic molecular information
    smiles: str
    inchi: str
    inchikey: str
    molecular_formula: str
    molecular_weight: float
    exact_mass: float
    heavy_atom_count: int
    atom_count: int
    heteroatom_count: int
    
    # Lipophilicity
    logp: float
    logp_crippen: float
    logp_wildman: float
    logd_ph2: float
    logd_ph5: float
    logd_ph74: float
    logd_ph9: float
    logd_ph11: float
    
    # Hydrogen bonding
    hba_lipinski: int
    hbd_lipinski: int
    hba_effective: float
    hbd_effective: float
    hba_estate: float
    hbd_estate: float
    
    # Surface area and volume
    tpsa: float
    labute_asa: float
    molar_volume: float
    molar_refractivity: float
    van_der_waals_volume: float
    solvent_accessible_surface: float
    
    # Topological descriptors
    rotatable_bonds: int
    rotatable_bonds_fraction: float
    aromatic_rings: int
    aliphatic_rings: int
    saturated_rings: int
    heterocycles: int
    spiro_atoms: int
    bridgehead_atoms: int
    balaban_j: float
    bertz_ct: float
    hall_kier_alpha: float
    kappa1: float
    kappa2: float
    kappa3: float
    
    # Electronic properties
    formal_charge: int
    total_charge: float
    max_partial_charge: float
    min_partial_charge: float
    estate_indices: List[float]
    vabc_vol: float
    
    # Constitutional descriptors
    fraction_csp3: float
    num_saturated_carbons: int
    num_unsaturated_carbons: int
    num_amide_bonds: int
    
    # Drug-likeness
    lipinski_violations: int
    ghose_filter: bool
    veber_filter: bool
    muegge_filter: bool
    qed_score: float
    synthetic_accessibility: float
    
    # Ionization profile
    ionization_type: IonizationType
    ionizable_groups: List[IonizableGroup]
    isoelectric_point: float
    pka_acidic_min: float
    pka_acidic_max: float
    pka_basic_min: float
    pka_basic_max: float
    
    # Solubility predictions
    logS: float
    solubility_class: str
    intrinsic_solubility: float  # mg/mL
    
    # Permeability predictions
    logP_app: float  # apparent permeability
    caco2_permeability: float
    mdck_permeability: float
    
    # Chromatographic descriptors
    hydrophobic_index: float
    hydrophilic_index: float
    amphiphilic_moment: float
    chromatographic_hydrophobicity: float
    silanol_interaction_potential: float
    pi_pi_interaction_potential: float
    steric_bulk_parameter: float
    hydrogen_bonding_potential: float
    
    # Fragment-based descriptors
    num_carboxylic_acids: int
    num_primary_amines: int
    num_secondary_amines: int
    num_tertiary_amines: int
    num_amides: int
    num_phenols: int
    num_alcohols: int
    num_thiols: int
    num_halogens: int
    num_nitro_groups: int
    
    # Additional calculated properties
    polarizability: float
    refractivity: float
    parm: float  # Parsimony index
    mr: float  # Molar refractivity
    
    # Metadata
    calculation_timestamp: str
    rdkit_version: str
    confidence_score: float
    
    def to_dict(self) -> Dict:
        """Convert to dictionary for serialization"""
        result = asdict(self)
        # Convert Enum to string
        result['ionization_type'] = self.ionization_type.value
        # Convert ionizable groups to dicts
        result['ionizable_groups'] = [asdict(g) for g in self.ionizable_groups]
        return result
    
    def to_json(self, indent: int = 2) -> str:
        """Convert to JSON string"""
        return json.dumps(self.to_dict(), indent=indent)
    
    def summary_table(self) -> Dict[str, Any]:
        """Return a summary of key properties for reporting"""
        return {
            "Molecular Formula": self.molecular_formula,
            "Molecular Weight": f"{self.molecular_weight:.3f} g/mol",
            "Exact Mass": f"{self.exact_mass:.5f} Da",
            "LogP": f"{self.logp:.3f}",
            "LogD (pH 2.0)": f"{self.logd_ph2:.3f}",
            "LogD (pH 5.0)": f"{self.logd_ph5:.3f}",
            "LogD (pH 7.4)": f"{self.logd_ph74:.3f}",
            "LogD (pH 9.0)": f"{self.logd_ph9:.3f}",
            "LogD (pH 11.0)": f"{self.logd_ph11:.3f}",
            "TPSA": f"{self.tpsa:.2f} Å²",
            "H-Bond Donors": str(self.hbd_lipinski),
            "H-Bond Acceptors": str(self.hba_lipinski),
            "Rotatable Bonds": str(self.rotatable_bonds),
            "Aromatic Rings": str(self.aromatic_rings),
            "Formal Charge": str(self.formal_charge),
            "Ionization Type": self.ionization_type.value.capitalize(),
            "pKa Range": f"{self.pka_acidic_min:.1f}-{self.pka_basic_max:.1f}",
            "LogS": f"{self.logS:.2f}",
            "Solubility Class": self.solubility_class,
            "Lipinski Violations": str(self.lipinski_violations),
            "QED Score": f"{self.qed_score:.3f}",
            "Chromatographic Hydrophobicity": f"{self.chromatographic_hydrophobicity:.3f}",
            "π-π Interaction Potential": f"{self.pi_pi_interaction_potential:.3f}",
            "Silanol Interaction Potential": f"{self.silanol_interaction_potential:.3f}"
        }


class PhysicochemicalCalculator:
    """
    Main calculator class for physicochemical properties.
    
    This class provides comprehensive property calculations including:
    - Basic molecular descriptors
    - Lipophilicity (LogP, LogD at various pH)
    - Ionization (pKa prediction)
    - Solubility prediction
    - Permeability estimation
    - Chromatographic behavior predictors
    
    Attributes:
        mol (Chem.Mol): RDKit molecule object
        smiles (str): Original SMILES string
        properties (PhysicochemicalProperties): Calculated properties
    """
    
    # pKa prediction constants (simplified model - in production use ChemAxon or similar)
    PKA_ACIDIC_GROUPS = {
        'carboxyl': 4.5,
        'phenol': 10.0,
        'thiol': 9.5,
        'sulfonamide': 10.5,
        'sulfonic_acid': -1.0,
        'phosphoric_acid': 2.5,
        'ammonium': 9.5,
    }
    
    PKA_BASIC_GROUPS = {
        'primary_amine': 9.5,
        'secondary_amine': 10.0,
        'tertiary_amine': 10.5,
        'aniline': 4.5,
        'pyridine': 5.5,
        'imidazole': 7.0,
        'guanidine': 13.5,
    }
    
    def __init__(self, smiles: str):
        """
        Initialize the calculator with a SMILES string.
        
        Args:
            smiles: SMILES string of the molecule
            
        Raises:
            ValueError: If SMILES string is invalid
        """
        self.smiles = smiles.strip()
        self.mol = Chem.MolFromSmiles(self.smiles)
        
        if self.mol is None:
            error_msg = f"Invalid SMILES string: {smiles}"
            logger.error(error_msg)
            raise ValueError(error_msg)
        
        # Sanitize molecule
        Chem.SanitizeMol(self.mol)
        
        # Generate 3D coordinates if needed (for volume calculations)
        self.mol_3d = None
        try:
            self.mol_3d = Chem.AddHs(self.mol)
            from rdkit.Chem import AllChem
            AllChem.EmbedMolecule(self.mol_3d, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(self.mol_3d)
        except Exception as e:
            logger.warning(f"Could not generate 3D structure: {e}")
        
        self.properties: Optional[PhysicochemicalProperties] = None
        logger.info(f"Initialized calculator for SMILES: {smiles}")
    
    def calculate_all(self) -> PhysicochemicalProperties:
        """
        Calculate all physicochemical properties.
        
        Returns:
            Complete PhysicochemicalProperties object
        """
        logger.info("Starting comprehensive property calculation")
        
        # Basic molecular information
        basic_info = self._calculate_basic_info()
        
        # Lipophilicity
        lipophilicity = self._calculate_lipophilicity()
        
        # Hydrogen bonding
        hbonding = self._calculate_hydrogen_bonding()
        
        # Surface and volume
        surface_volume = self._calculate_surface_volume()
        
        # Topological descriptors
        topological = self._calculate_topological()
        
        # Electronic properties
        electronic = self._calculate_electronic()
        
        # Constitutional descriptors
        constitutional = self._calculate_constitutional()
        
        # Drug-likeness
        drug_likeness = self._calculate_drug_likeness()
        
        # Ionization profile
        ionization = self._calculate_ionization()
        
        # Solubility predictions
        solubility = self._calculate_solubility()
        
        # Permeability predictions
        permeability = self._calculate_permeability()
        
        # Chromatographic descriptors
        chromatographic = self._calculate_chromatographic()
        
        # Fragment-based descriptors
        fragments = self._calculate_fragments()
        
        # Combine all properties
        self.properties = PhysicochemicalProperties(
            # Basic information
            smiles=self.smiles,
            inchi=basic_info['inchi'],
            inchikey=basic_info['inchikey'],
            molecular_formula=basic_info['formula'],
            molecular_weight=basic_info['mw'],
            exact_mass=basic_info['exact_mass'],
            heavy_atom_count=basic_info['heavy_atoms'],
            atom_count=basic_info['total_atoms'],
            heteroatom_count=basic_info['heteroatoms'],
            
            # Lipophilicity
            logp=lipophilicity['logp'],
            logp_crippen=lipophilicity['logp_crippen'],
            logp_wildman=lipophilicity['logp_wildman'],
            logd_ph2=lipophilicity['logd_ph2'],
            logd_ph5=lipophilicity['logd_ph5'],
            logd_ph74=lipophilicity['logd_ph74'],
            logd_ph9=lipophilicity['logd_ph9'],
            logd_ph11=lipophilicity['logd_ph11'],
            
            # Hydrogen bonding
            hba_lipinski=hbonding['hba_lipinski'],
            hbd_lipinski=hbonding['hbd_lipinski'],
            hba_effective=hbonding['hba_effective'],
            hbd_effective=hbonding['hbd_effective'],
            hba_estate=hbonding['hba_estate'],
            hbd_estate=hbonding['hbd_estate'],
            
            # Surface area and volume
            tpsa=surface_volume['tpsa'],
            labute_asa=surface_volume['labute_asa'],
            molar_volume=surface_volume['molar_volume'],
            molar_refractivity=surface_volume['molar_refractivity'],
            van_der_waals_volume=surface_volume['vdw_volume'],
            solvent_accessible_surface=surface_volume['sas'],
            
            # Topological descriptors
            rotatable_bonds=topological['rotatable_bonds'],
            rotatable_bonds_fraction=topological['rotatable_bonds_fraction'],
            aromatic_rings=topological['aromatic_rings'],
            aliphatic_rings=topological['aliphatic_rings'],
            saturated_rings=topological['saturated_rings'],
            heterocycles=topological['heterocycles'],
            spiro_atoms=topological['spiro_atoms'],
            bridgehead_atoms=topological['bridgehead_atoms'],
            balaban_j=topological['balaban_j'],
            bertz_ct=topological['bertz_ct'],
            hall_kier_alpha=topological['hall_kier_alpha'],
            kappa1=topological['kappa1'],
            kappa2=topological['kappa2'],
            kappa3=topological['kappa3'],
            
            # Electronic properties
            formal_charge=electronic['formal_charge'],
            total_charge=electronic['total_charge'],
            max_partial_charge=electronic['max_partial_charge'],
            min_partial_charge=electronic['min_partial_charge'],
            estate_indices=electronic['estate_indices'],
            vabc_vol=electronic['vabc_vol'],
            
            # Constitutional descriptors
            fraction_csp3=constitutional['fraction_csp3'],
            num_saturated_carbons=constitutional['saturated_carbons'],
            num_unsaturated_carbons=constitutional['unsaturated_carbons'],
            num_amide_bonds=constitutional['amide_bonds'],
            
            # Drug-likeness
            lipinski_violations=drug_likeness['lipinski_violations'],
            ghose_filter=drug_likeness['ghose_filter'],
            veber_filter=drug_likeness['veber_filter'],
            muegge_filter=drug_likeness['muegge_filter'],
            qed_score=drug_likeness['qed_score'],
            synthetic_accessibility=drug_likeness['synthetic_accessibility'],
            
            # Ionization profile
            ionization_type=ionization['ionization_type'],
            ionizable_groups=ionization['ionizable_groups'],
            isoelectric_point=ionization['isoelectric_point'],
            pka_acidic_min=ionization['pka_acidic_min'],
            pka_acidic_max=ionization['pka_acidic_max'],
            pka_basic_min=ionization['pka_basic_min'],
            pka_basic_max=ionization['pka_basic_max'],
            
            # Solubility predictions
            logS=solubility['logS'],
            solubility_class=solubility['solubility_class'],
            intrinsic_solubility=solubility['intrinsic_solubility'],
            
            # Permeability predictions
            logP_app=permeability['logP_app'],
            caco2_permeability=permeability['caco2'],
            mdck_permeability=permeability['mdck'],
            
            # Chromatographic descriptors
            hydrophobic_index=chromatographic['hydrophobic_index'],
            hydrophilic_index=chromatographic['hydrophilic_index'],
            amphiphilic_moment=chromatographic['amphiphilic_moment'],
            chromatographic_hydrophobicity=chromatographic['chrom_hydrophobicity'],
            silanol_interaction_potential=chromatographic['silanol_potential'],
            pi_pi_interaction_potential=chromatographic['pi_pi_potential'],
            steric_bulk_parameter=chromatographic['steric_bulk'],
            hydrogen_bonding_potential=chromatographic['hbond_potential'],
            
            # Fragment-based descriptors
            num_carboxylic_acids=fragments['carboxylic_acids'],
            num_primary_amines=fragments['primary_amines'],
            num_secondary_amines=fragments['secondary_amines'],
            num_tertiary_amines=fragments['tertiary_amines'],
            num_amides=fragments['amides'],
            num_phenols=fragments['phenols'],
            num_alcohols=fragments['alcohols'],
            num_thiols=fragments['thiols'],
            num_halogens=fragments['halogens'],
            num_nitro_groups=fragments['nitro'],
            
            # Additional properties
            polarizability=lipophilicity['polarizability'],
            refractivity=surface_volume['molar_refractivity'],
            parm=topological.get('parm', 0.0),
            mr=surface_volume['molar_refractivity'],
            
            # Metadata
            calculation_timestamp=datetime.now().isoformat(),
            rdkit_version=getattr(Chem, '__version__', getattr(__import__('rdkit'), '__version__', 'unknown')),
            confidence_score=self._calculate_confidence()
        )
        
        logger.info("Property calculation complete")
        return self.properties
    
    def _calculate_basic_info(self) -> Dict:
        """Calculate basic molecular information"""
        return {
            'inchi': Chem.MolToInchi(self.mol),
            'inchikey': Chem.MolToInchiKey(self.mol),
            'formula': rdMolDescriptors.CalcMolFormula(self.mol),
            'mw': Descriptors.MolWt(self.mol),
            'exact_mass': Descriptors.ExactMolWt(self.mol),
            'heavy_atoms': self.mol.GetNumHeavyAtoms(),
            'total_atoms': self.mol.GetNumAtoms(),
            'heteroatoms': CalcNumHeteroatoms(self.mol)
        }
    
    def _calculate_lipophilicity(self) -> Dict:
        """Calculate lipophilicity parameters including LogD at various pH"""
        # LogP calculations
        # RDKit implements the Wildman-Crippen method via Crippen.MolLogP.
        # Both "Crippen" and "Wildman" values refer to the same atom-additive model:
        # Wildman, S. A. & Crippen, G. M. (1999). J. Chem. Inf. Comput. Sci., 39(5), 868-873.
        # DOI: 10.1021/ci990307l
        logp_crippen = Crippen.MolLogP(self.mol)
        logp_wildman = logp_crippen  # Same underlying atom contributions (Wildman & Crippen 1999)
        
        # Simplified LogD calculation based on pKa (in production use proper model)
        ionization = self._calculate_ionization()
        
        # LogD at different pH values (simplified model)
        logd_ph2 = self._calculate_logd_at_ph(2.0, logp_crippen, ionization)
        logd_ph5 = self._calculate_logd_at_ph(5.0, logp_crippen, ionization)
        logd_ph74 = self._calculate_logd_at_ph(7.4, logp_crippen, ionization)
        logd_ph9 = self._calculate_logd_at_ph(9.0, logp_crippen, ionization)
        logd_ph11 = self._calculate_logd_at_ph(11.0, logp_crippen, ionization)
        
        # Polarizability
        polarizability = rdMolDescriptors.CalcNumHeteroatoms(self.mol) * 0.5 + 5
        
        return {
            'logp': round(logp_crippen, 3),
            'logp_crippen': round(logp_crippen, 3),
            'logp_wildman': round(logp_wildman, 3),
            'logd_ph2': round(logd_ph2, 3),
            'logd_ph5': round(logd_ph5, 3),
            'logd_ph74': round(logd_ph74, 3),
            'logd_ph9': round(logd_ph9, 3),
            'logd_ph11': round(logd_ph11, 3),
            'polarizability': round(polarizability, 3)
        }
    
    def _calculate_logd_at_ph(self, ph: float, logp: float, ionization: Dict) -> float:
        """
        Calculate LogD at a specific pH using Henderson-Hasselbalch equation.
        
        This is a simplified model. In production, use proper pKa prediction.
        """
        logd = logp
        
        # Adjust for ionization based on pKa
        if ionization['pka_acidic_min'] < ph < ionization['pka_basic_max']:
            # Near pKa region, reduce LogP
            fraction_ionized = 0.5
            logd -= fraction_ionized * 1.5  # Approximate reduction
        
        # Further adjustments based on pH
        if ph < 3:
            logd = logp + 0.2  # Slightly more hydrophobic in very acid
        elif ph > 10:
            logd = logp - 0.5  # More hydrophilic in basic conditions
        
        return logd
    
    def _calculate_hydrogen_bonding(self) -> Dict:
        """Calculate hydrogen bonding parameters"""
        hba_lipinski = Lipinski.NumHAcceptors(self.mol)
        hbd_lipinski = Lipinski.NumHDonors(self.mol)
        
        # Effective HBA/HBD (weighted by strength)
        hba_effective = hba_lipinski * 0.8
        hbd_effective = hbd_lipinski * 0.9
        
        # EState HBA/HBD
        estate = EStateIndices(self.mol)
        hba_estate = sum(estate) * 0.1  # Simplified
        hbd_estate = sum(estate) * 0.1  # Simplified (same as HBA for now)
        
        return {
            'hba_lipinski': hba_lipinski,
            'hbd_lipinski': hbd_lipinski,
            'hba_effective': round(hba_effective, 2),
            'hbd_effective': round(hbd_effective, 2),
            'hba_estate': round(hba_estate, 2),
            'hbd_estate': round(hbd_estate, 2)
        }
    def _calculate_surface_volume(self) -> Dict:
        """Calculate surface area and volume descriptors"""
        # TPSA
        tpsa = CalcTPSA(self.mol)
        
        # Labute ASA
        labute_asa = MolSurf.LabuteASA(self.mol)
        
        # Molar refractivity
        mr = Crippen.MolMR(self.mol)
        
        # Approximate molar volume (ml/mol)
        molar_volume = Descriptors.MolWt(self.mol) * 1.1
        
        # Van der Waals volume (Å³)
        vdw_volume = sum(atom.GetAtomicNum() * 10 for atom in self.mol.GetAtoms())  # Rough approximation
        
        # Solvent accessible surface (Å²)
        sas = tpsa * 2.5  # Approximation
        
        return {
            'tpsa': round(tpsa, 2),
            'labute_asa': round(labute_asa, 2),
            'molar_volume': round(molar_volume, 2),
            'molar_refractivity': round(mr, 3),
            'vdw_volume': round(vdw_volume, 2),
            'sas': round(sas, 2)
        }
    
    def _calculate_topological(self) -> Dict:
        """Calculate topological descriptors"""
        rotatable_bonds = CalcNumRotatableBonds(self.mol)
        total_bonds = self.mol.GetNumBonds()
        rotatable_fraction = rotatable_bonds / max(total_bonds, 1)
        
        return {
            'rotatable_bonds': rotatable_bonds,
            'rotatable_bonds_fraction': round(rotatable_fraction, 3),
            'aromatic_rings': CalcNumAromaticRings(self.mol),
            'aliphatic_rings': CalcNumAliphaticRings(self.mol),
            'saturated_rings': CalcNumSaturatedRings(self.mol),
            'heterocycles': CalcNumHeterocycles(self.mol),
            'spiro_atoms': CalcNumSpiroAtoms(self.mol),
            'bridgehead_atoms': CalcNumBridgeheadAtoms(self.mol),
            'balaban_j': round(GraphDescriptors.BalabanJ(self.mol), 3),
            'bertz_ct': round(GraphDescriptors.BertzCT(self.mol), 2),
            'hall_kier_alpha': round(Descriptors.HallKierAlpha(self.mol), 3),
            'kappa1': round(Descriptors.Kappa1(self.mol), 2),
            'kappa2': round(Descriptors.Kappa2(self.mol), 2),
            'kappa3': round(Descriptors.Kappa3(self.mol), 2)
        }
    
    def _calculate_electronic(self) -> Dict:
        """Calculate electronic and charge properties"""
        formal_charge = Chem.GetFormalCharge(self.mol)
        
        # Calculate Gasteiger charges
        # RDKit API changes: Try both rdMolDescriptors and AllChem locations
        try:
            from rdkit.Chem.rdMolDescriptors import CalcGasteigerCharges
            CalcGasteigerCharges(self.mol)
        except Exception:
            try:
                from rdkit.Chem.AllChem import ComputeGasteigerCharges
                ComputeGasteigerCharges(self.mol)
            except Exception as e:
                logger.warning(f"Gasteiger charge calculation not available: {e}")
        
        charges = [float(atom.GetProp('_GasteigerCharge')) for atom in self.mol.GetAtoms() 
                   if atom.HasProp('_GasteigerCharge')]
        
        max_charge = max(charges) if charges else 0
        min_charge = min(charges) if charges else 0
        total_charge = sum(charges) if charges else 0
        
        # EState indices
        estate = list(EStateIndices(self.mol))
        
        # VABC volume (not available in all RDKit versions)
        try:
            vabc_vol = rdMolDescriptors.CalcVABC(self.mol)
        except Exception:
            try:
                from rdkit.Chem.AllChem import ComputeMolVolume
                vabc_vol = ComputeMolVolume(self.mol)
            except Exception:
                vabc_vol = 0.0
        
        return {
            'formal_charge': formal_charge,
            'total_charge': round(total_charge, 3),
            'max_partial_charge': round(max_charge, 3),
            'min_partial_charge': round(min_charge, 3),
            'estate_indices': [round(e, 3) for e in estate[:10]],  # First 10
            'vabc_vol': round(vabc_vol, 2)
        }
    
    def _calculate_constitutional(self) -> Dict:
        """Calculate constitutional descriptors"""
        fraction_csp3 = CalcFractionCsp3(self.mol)
        
        # Count carbon types
        total_carbons = sum(1 for atom in self.mol.GetAtoms() if atom.GetAtomicNum() == 6)
        sp3_carbons = sum(1 for atom in self.mol.GetAtoms() 
                         if atom.GetAtomicNum() == 6 and atom.GetHybridization() == Chem.HybridizationType.SP3)
        
        saturated_carbons = sp3_carbons
        unsaturated_carbons = total_carbons - sp3_carbons
        
        # Amide bonds
        amide_bonds = CalcNumAmideBonds(self.mol)
        
        return {
            'fraction_csp3': round(fraction_csp3, 3),
            'saturated_carbons': saturated_carbons,
            'unsaturated_carbons': unsaturated_carbons,
            'amide_bonds': amide_bonds
        }
    
    def _calculate_drug_likeness(self) -> Dict:
        """Calculate drug-likeness parameters"""
        # Lipinski's Rule of Five
        violations = 0
        mw = Descriptors.MolWt(self.mol)
        logp = Crippen.MolLogP(self.mol)
        hbd = Lipinski.NumHDonors(self.mol)
        hba = Lipinski.NumHAcceptors(self.mol)
        
        if mw > 500: violations += 1
        if logp > 5: violations += 1
        if hbd > 5: violations += 1
        if hba > 10: violations += 1
        
        # Ghose filter
        ghose = (160 <= mw <= 480) and (-0.4 <= logp <= 5.6) and (20 <= self._count_atoms() <= 70)
        
        # Veber filter
        rotatable = CalcNumRotatableBonds(self.mol)
        tpsa = CalcTPSA(self.mol)
        veber = (rotatable <= 10) and (tpsa <= 140)
        
        # Muegge filter
        muegge = (200 <= mw <= 600) and (-2 <= logp <= 5) and (hbd <= 5) and (hba <= 10) and (tpsa <= 150)
        
        # QED score
        qed = QED.qed(self.mol)
        
        # Synthetic accessibility (approximation)
        synthetic_accessibility = 5.0 - (qed * 3)  # Rough approximation
        
        return {
            'lipinski_violations': violations,
            'ghose_filter': ghose,
            'veber_filter': veber,
            'muegge_filter': muegge,
            'qed_score': round(qed, 3),
            'synthetic_accessibility': round(synthetic_accessibility, 2)
        }
    
    def _count_atoms(self) -> int:
        """Count atoms excluding hydrogens"""
        return sum(1 for atom in self.mol.GetAtoms() if atom.GetAtomicNum() > 1)
    
    def _calculate_ionization(self) -> Dict:
        """
        Calculate ionization profile including pKa predictions.
        
        This is a simplified model. In production, use ChemAxon or similar.
        """
        ionizable_groups = []
        acidic_pkas = []
        basic_pkas = []
        
        # SMARTS patterns for ionizable groups
        patterns = {
            'carboxyl': ('C(=O)[OH]', 4.5, 'acidic'),
            'phenol': ('c1ccc(O)cc1', 10.0, 'acidic'),
            'primary_amine': ('[NX3;H2]', 9.5, 'basic'),
            'secondary_amine': ('[NX3;H1]', 10.0, 'basic'),
            'tertiary_amine': ('[NX3;H0]', 10.5, 'basic'),
            'aniline': ('[NX3;H2]c1ccccc1', 4.5, 'basic'),
            'pyridine': ('n1ccccc1', 5.5, 'basic'),
            'imidazole': ('c1cnc[nH]1', 7.0, 'amphoteric'),
            'thiol': ('[SH]', 9.5, 'acidic'),
            'sulfonic_acid': ('S(=O)(=O)O', -1.0, 'acidic'),
        }
        
        for name, (smarts, pka, group_type) in patterns.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern:
                matches = self.mol.GetSubstructMatches(pattern)
                if matches:
                    for match in matches:
                        group = IonizableGroup(
                            group_type=name,
                            pka=pka,
                            count=1,
                            position=str(match[0]),
                            confidence=0.8  # Simplified confidence
                        )
                        ionizable_groups.append(group)
                        
                        if group_type == 'acidic':
                            acidic_pkas.append(pka)
                        elif group_type == 'basic':
                            basic_pkas.append(pka)
                        else:  # amphoteric
                            acidic_pkas.append(pka)
                            basic_pkas.append(pka)
        
        # Determine ionization type
        if acidic_pkas and basic_pkas:
            ionization_type = IonizationType.AMPHOTERIC
        elif acidic_pkas:
            ionization_type = IonizationType.ACIDIC
        elif basic_pkas:
            ionization_type = IonizationType.BASIC
        else:
            ionization_type = IonizationType.NEUTRAL
        
        # Calculate isoelectric point (simplified)
        if acidic_pkas and basic_pkas:
            pI = (min(acidic_pkas) + max(basic_pkas)) / 2
        elif acidic_pkas:
            pI = min(acidic_pkas) - 1
        elif basic_pkas:
            pI = max(basic_pkas) + 1
        else:
            pI = 7.0
        
        return {
            'ionization_type': ionization_type,
            'ionizable_groups': ionizable_groups,
            'isoelectric_point': round(pI, 2),
            'pka_acidic_min': round(min(acidic_pkas) if acidic_pkas else 14.0, 2),
            'pka_acidic_max': round(max(acidic_pkas) if acidic_pkas else 14.0, 2),
            'pka_basic_min': round(min(basic_pkas) if basic_pkas else 0.0, 2),
            'pka_basic_max': round(max(basic_pkas) if basic_pkas else 0.0, 2)
        }
    
    def _calculate_solubility(self) -> Dict:
        """
        Calculate solubility predictions.
        
        This is a simplified model. In production, use trained ML models.
        """
        logp = Crippen.MolLogP(self.mol)
        mw = Descriptors.MolWt(self.mol)
        rotatable = CalcNumRotatableBonds(self.mol)
        aromatic = CalcNumAromaticRings(self.mol)
        
        # Simplified logS prediction (ESOL-like model)
        logS = 0.16 - 0.63 * logp - 0.0062 * mw + 0.066 * rotatable - 0.74 * aromatic
        
        # Solubility class
        if logS > -2:
            solubility_class = "High solubility"
            intrinsic_solubility = 10 ** (logS + 1)  # mg/mL approximation
        elif logS > -4:
            solubility_class = "Moderate solubility"
            intrinsic_solubility = 10 ** (logS + 1)
        elif logS > -6:
            solubility_class = "Low solubility"
            intrinsic_solubility = 10 ** (logS + 1)
        else:
            solubility_class = "Very low solubility"
            intrinsic_solubility = 10 ** (logS + 1)
        
        return {
            'logS': round(logS, 2),
            'solubility_class': solubility_class,
            'intrinsic_solubility': round(intrinsic_solubility, 3)
        }
    
    def _calculate_permeability(self) -> Dict:
        """
        Calculate permeability predictions.
        
        This is a simplified model. In production, use trained ML models.
        """
        logp = Crippen.MolLogP(self.mol)
        mw = Descriptors.MolWt(self.mol)
        tpsa = CalcTPSA(self.mol)
        rotatable = CalcNumRotatableBonds(self.mol)
        
        # Apparent permeability (logP_app)
        logP_app = logp * 0.8 - 0.01 * tpsa - 0.1 * (rotatable / 10)
        
        # Caco-2 permeability (10^-6 cm/s) approximation
        if logP_app > 1:
            caco2 = 20 + 10 * logP_app
        elif logP_app > 0:
            caco2 = 10 + 10 * logP_app
        else:
            caco2 = 5 + 5 * logP_app
        
        # MDCK permeability approximation (typically 1.5x Caco-2)
        mdck = caco2 * 1.5
        
        return {
            'logP_app': round(logP_app, 3),
            'caco2': round(caco2, 1),
            'mdck': round(mdck, 1)
        }
    
    def _calculate_chromatographic(self) -> Dict:
        """
        Calculate chromatographic behavior descriptors.
        
        These are specialized descriptors for HPLC method development.
        """
        logp = Crippen.MolLogP(self.mol)
        tpsa = CalcTPSA(self.mol)
        aromatic = CalcNumAromaticRings(self.mol)
        hbd = Lipinski.NumHDonors(self.mol)
        hba = Lipinski.NumHAcceptors(self.mol)
        rotatable = CalcNumRotatableBonds(self.mol)
        
        # Hydrophobic index (for reversed-phase retention)
        hydrophobic_index = 0.7 * logp + 0.3 * (aromatic / 5)
        
        # Hydrophilic index (for HILIC retention)
        hydrophilic_index = 0.5 * (tpsa / 200) + 0.3 * (hba + hbd) / 10
        
        # Amphiphilic moment (simplified)
        amphiphilic_moment = hydrophobic_index - hydrophilic_index
        
        # Chromatographic hydrophobicity (CHI) - approximate
        chrom_hydrophobicity = 10 * hydrophobic_index
        
        # Silanol interaction potential (for basic compounds)
        silanol_potential = (hbd * 0.3) + (tpsa / 500) - (logp * 0.1)
        silanol_potential = max(0, min(1, silanol_potential))  # Clamp to 0-1
        
        # π-π interaction potential
        pi_pi_potential = min(1, aromatic / 3)
        
        # Steric bulk parameter
        steric_bulk = min(1, rotatable / 15)
        
        # Hydrogen bonding potential
        hbond_potential = (hbd * 0.4 + hba * 0.2) / 5
        
        return {
            'hydrophobic_index': round(hydrophobic_index, 3),
            'hydrophilic_index': round(hydrophilic_index, 3),
            'amphiphilic_moment': round(amphiphilic_moment, 3),
            'chrom_hydrophobicity': round(chrom_hydrophobicity, 2),
            'silanol_potential': round(silanol_potential, 3),
            'pi_pi_potential': round(pi_pi_potential, 3),
            'steric_bulk': round(steric_bulk, 3),
            'hbond_potential': round(hbond_potential, 3)
        }
    
    def _calculate_fragments(self) -> Dict:
        """Count specific functional groups"""
        # Define SMARTS patterns
        patterns = {
            'carboxylic_acids': '[CX3](=O)[OX2H]',
            'primary_amines': '[NX3;H2]',
            'secondary_amines': '[NX3;H1]',
            'tertiary_amines': '[NX3;H0]',
            'amides': '[CX3](=O)[NX3]',
            'phenols': '[cX3][OH]',
            'alcohols': '[CX4][OH]',
            'thiols': '[SH]',
            'halogens': '[F,Cl,Br,I]',
            'nitro': '[NX3](=O)=O',
        }
        
        counts = {}
        for name, smarts in patterns.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern:
                counts[name] = len(self.mol.GetSubstructMatches(pattern))
            else:
                counts[name] = 0
        
        return counts
    
    def _calculate_confidence(self) -> float:
        """Calculate confidence score for the predictions"""
        confidence = 0.95  # Base confidence
        
        # Reduce confidence for certain situations
        if self.mol.GetNumAtoms() > 100:
            confidence -= 0.1
        
        if any(atom.GetAtomicNum() > 20 for atom in self.mol.GetAtoms()):
            confidence -= 0.1
        
        if self.mol.GetNumConformers() == 0 and self.mol_3d is None:
            confidence -= 0.05
        
        return round(confidence, 2)
    
    def get_summary(self) -> Dict[str, Any]:
        """Get a summary of key properties"""
        if not self.properties:
            self.calculate_all()
        return self.properties.summary_table()
    
    def export_to_json(self, filename: Optional[str] = None) -> str:
        """Export properties to JSON file"""
        if not self.properties:
            self.calculate_all()
        
        json_str = self.properties.to_json()
        
        if filename:
            with open(filename, 'w') as f:
                f.write(json_str)
            logger.info(f"Properties exported to {filename}")
        
        return json_str
    
    def compare_with(self, other_smiles: str) -> Dict[str, float]:
        """
        Compare this molecule with another SMILES.
        
        Args:
            other_smiles: SMILES string to compare with
            
        Returns:
            Dictionary of similarity metrics
        """
        other_calc = PhysicochemicalCalculator(other_smiles)
        other_props = other_calc.calculate_all()
        
        if not self.properties:
            self.calculate_all()
        
        # Calculate similarities for key properties
        similarities = {
            'mw_similarity': 1 - abs(self.properties.molecular_weight - other_props.molecular_weight) / 500,
            'logp_similarity': 1 - abs(self.properties.logp - other_props.logp) / 10,
            'tpsa_similarity': 1 - abs(self.properties.tpsa - other_props.tpsa) / 200,
            'hba_similarity': 1 - abs(self.properties.hba_lipinski - other_props.hba_lipinski) / 20,
            'hbd_similarity': 1 - abs(self.properties.hbd_lipinski - other_props.hbd_lipinski) / 10,
        }
        
        # Clamp values to 0-1
        for key in similarities:
            similarities[key] = max(0, min(1, similarities[key]))
        
        return similarities


# Fallback implementation for environments without RDKit
if not USE_RDKIT:
    import re
    from types import SimpleNamespace

    ATOMIC_MASSES = {
        'H': 1.0079, 'C': 12.0107, 'N': 14.0067, 'O': 15.999, 'S': 32.065,
        'P': 30.9738, 'Cl': 35.453, 'Br': 79.904, 'F': 18.998, 'I': 126.90
    }

    class PhysicochemicalCalculator:
        """Lightweight fallback for environments without RDKit."""

        def __init__(self, smiles: str):
            self.smiles = smiles or ""

        def _element_counts(self):
            s = self.smiles
            tokens = re.findall(r'Br|Cl|[A-Z][a-z]?|\[.*?\]', s)
            counts = {}
            for t in tokens:
                if t.startswith('[') and t.endswith(']'):
                    t2 = re.sub(r'[^A-Za-z]', '', t)
                    if t2 == '':
                        continue
                    t = t2
                counts[t] = counts.get(t, 0) + 1
            return counts

        def _estimate_mw(self, counts):
            mw = 0.0
            for el, n in counts.items():
                mass = ATOMIC_MASSES.get(el, 12.01)
                mw += mass * n
            return mw

        def calculate_all(self):
            counts = self._element_counts()
            nC = counts.get('C', 0) + counts.get('c', 0)
            nO = counts.get('O', 0)
            nN = counts.get('N', 0)
            mw = self._estimate_mw(counts)

            # Simple heuristic estimates
            logp = round(0.5 * nC - 0.3 * (nO + nN), 3)
            tpsa = round(20.0 * (nO + nN), 2)
            hba = nO + nN
            hbd = max(0, int(nN * 0.5))
            rot = max(0, int(nC * 0.1))
            arom = int(self.smiles.count('c') / 6)

            props = SimpleNamespace()
            props.smiles = self.smiles
            props.molecular_formula = ''.join(f"{k}{v if v>1 else ''}" for k,v in counts.items()) or 'C'
            props.molecular_weight = mw
            props.exact_mass = mw
            props.logp = logp
            props.logp_crippen = logp
            props.logp_wildman = logp
            props.logd_ph2 = logp
            props.logd_ph5 = logp
            props.logd_ph74 = logp
            props.logd_ph9 = logp
            props.logd_ph11 = logp
            props.tpsa = tpsa
            props.hba_lipinski = hba
            props.hbd_lipinski = hbd
            props.rotatable_bonds = rot
            props.aromatic_rings = arom
            props.aliphatic_rings = 0
            props.fraction_csp3 = 0.5
            props.logS = -3.0
            props.solubility_class = 'Moderate' if props.logS > -4 else 'Low'
            props.lipinski_violations = int(props.molecular_weight > 500) + int(props.logp > 5) + int(props.hbd_lipinski > 5) + int(props.hba_lipinski > 10)
            props.qed_score = 0.0
            props.ionization_type = IonizationType.NEUTRAL
            props.isoelectric_point = 7.0
            props.pka_acidic_min = 0.0
            props.pka_acidic_max = 0.0
            props.pka_basic_min = 0.0
            props.pka_basic_max = 0.0

            # Basic chromatographic descriptors for fallback
            props.hydrophobic_index = logp
            props.hydrophilic_index = tpsa / 100.0
            props.chromatographic_hydrophobicity = logp * 10
            props.silanol_interaction_potential = min(1.0, max(0.0, (hbd * 0.2 + tpsa / 500) - logp * 0.1))
            props.pi_pi_interaction_potential = min(1.0, arom / 3)
            props.steric_bulk_parameter = min(1.0, rot / 15)
            props.hydrogen_bonding_potential = min(1.0, (hbd * 0.4 + hba * 0.2) / 5)

            return props


# Additional utility functions for common calculations

def calculate_physchem_properties(smiles: str) -> Dict[str, Any]:
    """
    Convenience function to calculate physicochemical properties.
    
    Args:
        smiles: SMILES string
        
    Returns:
        Dictionary of calculated properties
    """
    calculator = PhysicochemicalCalculator(smiles)
    properties = calculator.calculate_all()
    return properties.summary_table()


def batch_calculate(smiles_list: List[str]) -> List[Dict[str, Any]]:
    """
    Calculate properties for multiple SMILES strings.
    
    Args:
        smiles_list: List of SMILES strings
        
    Returns:
        List of property dictionaries
    """
    results = []
    for smiles in smiles_list:
        try:
            calculator = PhysicochemicalCalculator(smiles)
            props = calculator.calculate_all()
            results.append(props.summary_table())
        except Exception as e:
            logger.error(f"Error processing {smiles}: {e}")
            results.append({"smiles": smiles, "error": str(e)})
    return results


def calculate_chromatographic_descriptors(smiles: str) -> Dict[str, float]:
    """
    Calculate specialized chromatographic descriptors.
    
    Args:
        smiles: SMILES string
        
    Returns:
        Dictionary of chromatographic descriptors
    """
    calculator = PhysicochemicalCalculator(smiles)
    properties = calculator.calculate_all()
    
    return {
        'hydrophobic_index': properties.hydrophobic_index,
        'hydrophilic_index': properties.hydrophilic_index,
        'chrom_hydrophobicity': properties.chromatographic_hydrophobicity,
        'silanol_potential': properties.silanol_interaction_potential,
        'pi_pi_potential': properties.pi_pi_interaction_potential,
        'steric_bulk': properties.steric_bulk_parameter,
        'hbond_potential': properties.hydrogen_bonding_potential
    }


if __name__ == "__main__":
    # Test the calculator with example molecules
    import sys
    from datetime import datetime
    
    test_smiles = [
        "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",  # Ibuprofen
        "CC(=O)OC1=CC=CC=C1C(=O)O",        # Aspirin
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",    # Caffeine
        "C1=CC=C2C(=C1)C(=O)NC2=O",        # Phthalimide
    ]
    
    print("=" * 80)
    print("PHYSICOCHEMICAL PROPERTY CALCULATOR TEST")
    print("=" * 80)
    
    for smiles in test_smiles:
        print(f"\nTesting SMILES: {smiles}")
        try:
            calculator = PhysicochemicalCalculator(smiles)
            props = calculator.calculate_all()
            
            print(f"  Formula: {props.molecular_formula}")
            print(f"  MW: {props.molecular_weight:.2f}")
            print(f"  LogP: {props.logp:.2f}")
            print(f"  LogD (pH 7.4): {props.logd_ph74:.2f}")
            print(f"  TPSA: {props.tpsa:.1f} Å²")
            print(f"  Ionization: {props.ionization_type.value}")
            print(f"  pKa range: {props.pka_acidic_min:.1f}-{props.pka_basic_max:.1f}")
            print(f"  LogS: {props.logS:.2f}")
            print(f"  Solubility: {props.solubility_class}")
            print(f"  CHI: {props.chromatographic_hydrophobicity:.1f}")
            
        except Exception as e:
            print(f"  ERROR: {e}")
    
    print("\n" + "=" * 80)
    print("Test complete")
    print("=" * 80)