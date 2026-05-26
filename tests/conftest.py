"""
Pytest configuration and shared fixtures for FastPhase.AI
"""

import pytest
import sys
import os

# Add paths
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))

@pytest.fixture
def ibuprofen_smiles():
    """Return Ibuprofen SMILES string"""
    return "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"

@pytest.fixture
def caffeine_smiles():
    """Return Caffeine SMILES string"""
    return "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"

@pytest.fixture
def aspirin_smiles():
    """Return Aspirin SMILES string"""
    return "CC(=O)OC1=CC=CC=C1C(=O)O"

@pytest.fixture
def paracetamol_smiles():
    """Return Paracetamol SMILES string"""
    return "CC(=O)NC1=CC=C(C=C1)O"
