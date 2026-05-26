"""
FastPhase.AI - Comprehensive Test Suite
Tests all 5 core modules with real SMILES strings
"""

import sys
import os
import pytest

# Add python directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))

# ============================================================
# Module Import Tests
# ============================================================

def test_import_physchem():
    """Test PhysicochemicalCalculator imports correctly"""
    from physchem_calculator import PhysicochemicalCalculator
    assert PhysicochemicalCalculator is not None

def test_import_hsm():
    """Test HSMEstimator imports correctly"""
    from HSMsolute_check import HSMEstimator
    assert HSMEstimator is not None

def test_import_buffer():
    """Test BufferSelector imports correctly"""
    from buffer_selector import BufferSelector
    assert BufferSelector is not None

def test_import_solvent():
    """Test SolventSelectionSystem imports correctly"""
    from solvent_selector import SolventSelectionSystem
    assert SolventSelectionSystem is not None

def test_import_column():
    """Test HSMColumnSelector imports correctly"""
    from column_selector import HSMColumnSelector
    assert HSMColumnSelector is not None

def test_import_gradient():
    """Test GradientOptimizer imports correctly"""
    from gradient_optimizer import GradientOptimizer
    assert GradientOptimizer is not None

def test_import_reference_manager():
    """Test ReferenceManager imports correctly"""
    from reference_manager import ChromatographyReferenceManager
    assert ChromatographyReferenceManager is not None


# ============================================================
# Physicochemical Calculator Tests
# ============================================================

class TestPhysicochemicalCalculator:
    """Test suite for PhysicochemicalCalculator"""
    
    def test_ibuprofen_properties(self):
        """Test Ibuprofen SMILES processing"""
        from physchem_calculator import PhysicochemicalCalculator
        
        smiles = "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
        calc = PhysicochemicalCalculator(smiles)
        props = calc.calculate_all()
        
        # Assertions
        assert props.molecular_weight > 200
        assert props.molecular_weight < 210
        assert props.logp > 3
        assert props.logp < 5
        assert props.tpsa > 30
        assert props.tpsa < 50
        assert props.hba_lipinski > 1
        assert props.hbd_lipinski >= 1
    
    def test_caffeine_properties(self):
        """Test Caffeine SMILES processing"""
        from physchem_calculator import PhysicochemicalCalculator
        
        smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
        calc = PhysicochemicalCalculator(smiles)
        props = calc.calculate_all()
        
        assert props.molecular_weight > 190
        assert props.molecular_weight < 200
        assert props.logp > -1
        assert props.logp < 1
        assert props.tpsa > 50
        assert props.tpsa < 70
    
    def test_aspirin_properties(self):
        """Test Aspirin SMILES processing"""
        from physchem_calculator import PhysicochemicalCalculator
        
        smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
        calc = PhysicochemicalCalculator(smiles)
        props = calc.calculate_all()
        
        assert props.molecular_weight > 170
        assert props.molecular_weight < 190
        assert props.logp > 1
        assert props.logp < 3


# ============================================================
# HSM Estimator Tests
# ============================================================

class TestHSMEstimator:
    """Test suite for HSMEstimator"""
    
    def test_hsm_ibuprofen(self):
        """Test HSM descriptor calculation for Ibuprofen"""
        from HSMsolute_check import HSMEstimator
        
        estimator = HSMEstimator(pH=7.0)
        smiles = "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
        descriptors = estimator.calculate_from_smiles(smiles)
        
        assert 'eta_prime' in descriptors
        assert 'sigma_prime' in descriptors
        assert 'beta_prime' in descriptors
        assert 'alpha_prime' in descriptors
        assert 'kappa_prime' in descriptors
        
        # Ibuprofen is acidic
        assert descriptors['alpha_prime'] > 0.3
    
    def test_hsm_caffeine(self):
        """Test HSM descriptor calculation for Caffeine"""
        from HSMsolute_check import HSMEstimator
        
        estimator = HSMEstimator(pH=7.0)
        smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
        descriptors = estimator.calculate_from_smiles(smiles)
        
        # Caffeine is neutral/basic
        assert descriptors['kappa_prime'] >= 0


# ============================================================
# Buffer Selector Tests
# ============================================================

class TestBufferSelector:
    """Test suite for BufferSelector"""
    
    def test_buffer_selector_initializes(self):
        """Test BufferSelector initializes correctly"""
        from buffer_selector import BufferSelector
        
        selector = BufferSelector()
        assert selector is not None
    
    def test_buffer_selection_ibuprofen(self):
        """Test buffer selection for Ibuprofen (acidic)"""
        from buffer_selector import BufferSelector
        
        selector = BufferSelector()
        smiles = "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
        
        params = {
            'target_ph': 3.5,
            'detection_wavelength_nm': 230,
            'is_lcms': False
        }
        
        results = selector.select_buffer(smiles, method_params=params, output_excel=False)
        
        assert 'top_buffers' in results
        assert len(results['top_buffers']) > 0


# ============================================================
# Solvent Selector Tests
# ============================================================

class TestSolventSelector:
    """Test suite for SolventSelector"""
    
    def test_solvent_selector_initializes(self):
        """Test SolventSelectionSystem initializes correctly"""
        from solvent_selector import SolventSelectionSystem
        
        selector = SolventSelectionSystem()
        assert selector is not None
    
    def test_solvent_selection_caffeine(self):
        """Test solvent selection for Caffeine"""
        from solvent_selector import SolventSelectionSystem
        
        selector = SolventSelectionSystem()
        smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
        
        results = selector.process_single_smiles(smiles, name="Caffeine")
        
        assert results is not None
        assert 'recommendations' in results
        assert len(results['recommendations']) > 0


# ============================================================
# Column Selector Tests
# ============================================================

class TestColumnSelector:
    """Test suite for ColumnSelector"""
    
    def test_column_selector_initializes(self):
        """Test HSMColumnSelector initializes correctly"""
        from column_selector import HSMColumnSelector
        
        selector = HSMColumnSelector()
        assert selector is not None
    
    def test_column_selection_ibuprofen(self):
        """Test column selection for Ibuprofen"""
        from column_selector import HSMColumnSelector
        
        selector = HSMColumnSelector()
        smiles = "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
        
        results = selector.select_columns_for_smiles(smiles, n_recommendations=5, show_references=False)
        
        assert results is not None
        assert len(results) > 0


# ============================================================
# Gradient Optimizer Tests
# ============================================================

class TestGradientOptimizer:
    """Test suite for GradientOptimizer"""
    
    def test_gradient_optimizer_initializes(self):
        """Test GradientOptimizer initializes correctly"""
        from gradient_optimizer import GradientOptimizer
        
        optimizer = GradientOptimizer()
        assert optimizer is not None
    
    def test_gradient_optimization(self):
        """Test gradient optimization runs"""
        from gradient_optimizer import GradientOptimizer, OptimizationObjective, GradientType
        
        optimizer = GradientOptimizer()
        
        compounds = [{
            'name': 'Test Compound',
            'logp': 3.0,
            'molecular_weight': 200,
            'tpsa': 50
        }]
        
        result = optimizer.optimize_gradient(
            compounds=compounds,
            objective=OptimizationObjective.BALANCED,
            gradient_type=GradientType.LINEAR
        )
        
        assert result is not None
        assert result.gradient_program.total_runtime > 0


# ============================================================
# Integration Tests
# ============================================================

class TestIntegration:
    """Integration tests for the entire system"""
    
    def test_main_controller_import(self):
        """Test main controller imports"""
        # Add parent directory to path for main.py
        sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
        
        from main import ChromatographyAIController
        
        controller = ChromatographyAIController()
        assert controller is not None
    
    def test_full_pipeline_ibuprofen(self):
        """Test full pipeline with Ibuprofen"""
        sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
        
        from main import ChromatographyAIController
        
        controller = ChromatographyAIController()
        smiles = "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
        
        results = controller.process_single_compound(smiles, name="Ibuprofen", project="Test")
        
        assert results is not None
        assert results['status'] == 'Success' or results['status'] == 'Partial'
        assert 'scores' in results


# ============================================================
# Run tests if executed directly
# ============================================================

if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
