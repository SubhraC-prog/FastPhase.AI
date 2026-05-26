# FastPhase.AI  
<div align="center">

# 🧪 FastPhase.AI  
### *AI-Assisted HPLC Method Development System*

<img src="https://img.shields.io/badge/version-3.0-blue?style=for-the-badge" />
<img src="https://img.shields.io/badge/python-3.8+-green?style=for-the-badge" />
<img src="https://img.shields.io/badge/license-MIT-orange?style=for-the-badge" />
<img src="https://img.shields.io/badge/HPLC-AI%20Powered-red?style=for-the-badge" />

<br>

### Intelligent Reversed-Phase HPLC Method Development  
### Built on Peer-Reviewed Science • Powered by AI

</div>

---

# 📌 Overview

**FastPhase.AI** is a comprehensive **AI-driven reversed-phase HPLC method development platform** designed for analytical scientists, pharmaceutical researchers, and chromatography experts.

The system integrates:

- Hydrophobic Subtraction Model (HSM)
- Physicochemical property prediction
- Intelligent buffer & solvent selection
- HSM-based column recommendation
- Gradient optimization
- Automated reference tracking
- Professional Excel & PDF report generation

Using a single **SMILES string**, FastPhase.AI performs multi-stage chromatographic analysis and produces publication-ready outputs aligned with modern regulatory expectations.

---

# ✨ Key Features

## 🔬 AI-Powered Analytical Engine
- Physicochemical descriptor calculation
- HSM descriptor estimation
- Solvent-system optimization
- Buffer recommendation engine
- Column selectivity prediction
- Gradient optimization using LSS theory

## 📚 Reference-Driven Architecture
- 30+ peer-reviewed scientific references
- DOI-linked scientific rules
- Automated citation management
- Vancouver / ACS / APA formatting support

## 📈 Professional Reporting
- Multi-sheet Excel reports
- PDF reports with formatted tables
- Regulatory alignment summaries
- Confidence scoring & recommendations

## ⚙️ Enterprise-Style Workflow
- Excel VBA frontend
- Python backend architecture
- Batch processing support
- Logging & caching system

---

# 🏗️ System Architecture

```text
┌─────────────────────────────────────────────────────────────────────────────────────┐
│                           FASTPHASE.AI SYSTEM ARCHITECTURE                         │
├─────────────────────────────────────────────────────────────────────────────────────┤

    Excel Frontend (VBA)
            │
            ▼
    Python Analytical Engine
            │
 ┌──────────┼──────────────────────────────────────────┐
 │          │          │         │         │           │
 ▼          ▼          ▼         ▼         ▼           ▼

PhysChem   HSM      Buffer    Solvent    Column    Gradient
Engine    Engine    Engine     Engine     Engine    Optimizer

            │
            ▼

     Reference Manager
            │
            ▼

   Excel Report + PDF Report

└─────────────────────────────────────────────────────────────────────────────────────┘
```

---

# 📂 Project Structure

```text
FastPhase.AI/
│
├── FastPhaseAI.xlsm                 # Excel frontend
├── main.py                          # Main controller
│
├── python/
│   ├── physchem_calculator.py
│   ├── HSMsolute_check.py
│   ├── buffer_selector.py
│   ├── solvent_selector.py
│   ├── column_selector.py
│   ├── gradient_optimizer.py
│   ├── report_template.py
│   ├── reference_manager.py
│   └── excel_formatter.py
│
├── output/
│   ├── reports/
│   └── logs/
│
├── rules/
│
├── references/
│   └── master_references.json
│
├── requirements.txt
├── LICENSE
└── README.md
```

---

# 🧠 Core Modules

---

## 1️⃣ Physicochemical Calculator

**File:** `physchem_calculator.py`

Calculates molecular descriptors using **RDKit**.

### Features
- LogP / LogD
- TPSA
- HBA / HBD
- pKa prediction
- QED drug-likeness score
- Silanol interaction potential
- π-π interaction descriptors

---

## 2️⃣ HSM Descriptor Estimator

**File:** `HSMsolute_check.py`

Implements the complete **Hydrophobic Subtraction Model (HSM)**.

### Calculated Parameters

| Descriptor | Meaning |
|---|---|
| η′ | Hydrophobicity |
| σ′ | Steric resistance |
| β′ | H-bond basicity |
| α′ | H-bond acidity |
| κ′ | Cationic interaction |

---

## 3️⃣ Buffer Selection Engine

**File:** `buffer_selector.py`

### 10-Rule Scientific Scoring System

| Rule | Description | Weight |
|---|---|---|
| 1 | pKa Matching | 20% |
| 2 | UV Transparency | 15% |
| 3 | Organic Solubility | 12% |
| 4 | Buffer Capacity | 10% |
| 5 | MS Volatility | 15% |
| 6 | Temperature Dependence | 8% |
| 7 | Chemical Reactivity | 10% |
| 8 | Ionic Strength | 3% |
| 9 | Storage Stability | 2% |
| 10 | Metal Complexation | 5% |

---

## 4️⃣ Solvent Selection Engine

**File:** `solvent_selector.py`

### 7-Rule Intelligent Solvent Ranking

| Rule | Scientific Basis |
|---|---|
| LogP Matching | Solute-solvent compatibility |
| HBD/HBA Complementarity | Hydrogen-bond balancing |
| Polarity Matching | Selectivity alignment |
| Kamlet-Taft Distance | Solvation similarity |
| Viscosity/Pressure | System efficiency |
| UV Transparency | Detector compatibility |
| pH Stability | Column safety |

---

## 5️⃣ Column Recommendation Engine

**File:** `column_selector.py`

Uses **USP PQRI HSM column database** for selectivity-based ranking.

### Features
- Basic compound scoring
- Acidic compound scoring
- Neutral compound scoring
- Fs factor calculation
- Column equivalence assessment

---

## 6️⃣ Gradient Optimizer

**File:** `gradient_optimizer.py`

Implements **Linear Solvent Strength (LSS) Theory**.

### Capabilities
- Gradient prediction
- Peak capacity optimization
- Resolution maximization
- Monte Carlo robustness simulation
- Multi-linear gradient design

---

## 7️⃣ Reference Manager

**File:** `reference_manager.py`

Centralized scientific reference engine.

### Supports
- Vancouver citations
- ACS formatting
- APA formatting
- DOI tracking
- Reference export

---

# 📊 Generated Reports

FastPhase.AI automatically creates comprehensive reports in:

- 📗 Excel (`.xlsx`)
- 📕 PDF (`.pdf`)

---

## Report Sections

| Section | Content |
|---|---|
| 1 | Cover Sheet |
| 2 | Summary & Confidence Scores |
| 3 | Physicochemical Properties |
| 4 | HSM Descriptors |
| 5 | Column Recommendations |
| 6 | Solvent Recommendations |
| 7 | Buffer Recommendations |
| 8 | Gradient Program |
| 9 | References |
| 10 | Regulatory Alignment |

---

# 📖 Rule Repository

The `rules/` directory contains fully documented scientific rules including:

- Mathematical equations
- Threshold criteria
- Optimization ranges
- Literature references
- DOI citations
- Implementation notes

---

# 🚀 Installation

## Prerequisites

- Python 3.8+
- Microsoft Excel (Macros Enabled)
- RDKit
- ReportLab

---

## Clone Repository

```bash
git clone https://github.com/yourusername/FastPhase.AI.git
cd FastPhase.AI
```

---

## Install Dependencies

```bash
pip install -r requirements.txt
```

---

## Install RDKit (Recommended)

```bash
conda install -c conda-forge rdkit
```

---

## Verify Installation

```bash
python main.py --help
```

---

# 💻 Usage

## Method 1 — Excel Frontend (Recommended)

1. Open `FastPhaseAI.xlsm`
2. Enable macros
3. Enter SMILES string
4. Add project information
5. Click **Generate Report**

---

## Method 2 — Command Line

### Process Single Compound

```bash
python main.py \
--smiles "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O" \
--name "Ibuprofen" \
--project "NSAID_Analysis"
```

---

### Generate Excel + PDF

```bash
python main.py \
--smiles "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" \
--name "Caffeine" \
--format both
```

---

### Batch Processing

```bash
python main.py --batch --excel "FastPhaseAI.xlsm"
```

---

### Validate SMILES

```bash
python main.py \
--validate-only \
--smiles "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
```

---

# 🐍 Python API Example

```python
from main import ChromatographyAIController

controller = ChromatographyAIController()

results = controller.process_single_compound(
    smiles="CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    name="Ibuprofen",
    project="Method Development"
)

report_paths = controller.generate_reports(
    results,
    formats=['excel', 'pdf']
)

print(report_paths)
```

---

# 📚 Scientific Foundations

FastPhase.AI is built upon peer-reviewed chromatographic science.

## Key References

| Reference | Application |
|---|---|
| Snyder et al., 2004 | Hydrophobic Subtraction Model |
| Dolan et al., 2004 | Fs Factor |
| Marchand et al., 2005/2008 | HSM Descriptor Estimation |
| Valkó, 2004 | LogP Matching |
| Vitha & Carr, 2006 | HBD/HBA Theory |
| Goldberg et al., 2002 | Buffer pKa Theory |
| Kebarle & Tang, 1993 | MS Volatility |

---

# 📜 Regulatory Alignment

FastPhase.AI incorporates principles from:

- ICH Q2(R2)
- ICH Q8(R2)
- ICH Q9
- ICH Q14
- USP <621>
- USP <1220>

---

# 🛠️ Troubleshooting

## RDKit Not Installed

```bash
conda install -c conda-forge rdkit
```

---

## ReportLab Missing

```bash
pip install reportlab
```

---

## Excel Macro Issues

- Enable macros
- Verify Python path in VBA
- Keep `.xlsm` in project root directory

---

# 📄 License

Distributed under the **MIT License**.

See `LICENSE` for more information.

---

# 📌 Citation

```text
FastPhase.AI: AI-Assisted HPLC Method Development System.

https://github.com/yourusername/FastPhase.AI
```

---

# 🙏 Acknowledgments

- Snyder & Dolan — Hydrophobic Subtraction Model
- USP PQRI — Column Database
- RDKit — Open-source cheminformatics
- Analytical chromatography research community

---

<div align="center">

# 🧪 FastPhase.AI

### Intelligent HPLC Method Development

### Built on Peer-Reviewed Science • Powered by AI

</div>

