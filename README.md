#FastPhase.AI

#Overview
FastPhase.AI is a comprehensive, AI-driven system for reversed-phase HPLC method development. It implements the complete Hydrophobic Subtraction Model (HSM) with full reference tracking, physicochemical property calculations, buffer/solvent/column selection algorithms, and gradient optimization—all backed by peer-reviewed literature.

The system processes SMILES strings through five integrated analytical modules and generates professional, publication-ready reports in Excel and PDF formats.

System Architecture
text
┌─────────────────────────────────────────────────────────────────────────────────────┐
│                           FASTPHASE.AI SYSTEM ARCHITECTURE                           │
├─────────────────────────────────────────────────────────────────────────────────────┤
│                                                                                      │
│  ┌──────────────────────────────────────────────────────────────────────────────┐   │
│  │                         EXCEL FRONTEND (FastPhaseAI.xlsm)                     │   │
│  │  ┌─────────┐ ┌─────────┐ ┌─────────┐ ┌─────────┐ ┌─────────┐ ┌─────────┐    │   │
│  │  │DASHBOARD│ │ History │ │Templates│ │ Settings│ │  Help   │ │  About  │    │   │
│  │  └────┬────┘ └─────────┘ └─────────┘ └─────────┘ └─────────┘ └─────────┘    │   │
│  └───────┼──────────────────────────────────────────────────────────────────────┘   │
│          │ VBA → Command Line                                                       │
│          ▼                                                                           │
│  ┌──────────────────────────────────────────────────────────────────────────────┐   │
│  │                         PYTHON BACKEND (main.py)                              │   │
│  │                                                                               │   │
│  │  ┌─────────────┐ ┌─────────────┐ ┌─────────────┐ ┌─────────────┐ ┌─────────┐│   │
│  │  │ PhysChem    │ │    HSM      │ │   Buffer    │ │   Solvent   │ │ Column  ││   │
│  │  │ Calculator  │ │  Estimator  │ │  Selector   │ │  Selector   │ │Selector ││   │
│  │  └──────┬──────┘ └──────┬──────┘ └──────┬──────┘ └──────┬──────┘ └────┬────┘│   │
│  │         │               │               │               │            │      │   │
│  │         └───────────────┴───────────────┴───────────────┴────────────┘      │   │
│  │                                    │                                         │   │
│  │                          ┌──────────┴──────────┐                             │   │
│  │                          │ Gradient Optimizer  │                             │   │
│  │                          └──────────┬──────────┘                             │   │
│  │                                     │                                        │   │
│  │                          ┌──────────┴──────────┐                             │   │
│  │                          │ Reference Manager   │                             │   │
│  │                          └──────────┬──────────┘                             │   │
│  └─────────────────────────────────────┼────────────────────────────────────────┘   │
│                                        │                                           │
│                          ┌─────────────┴─────────────┐                             │
│                          ▼                           ▼                             │
│              ┌───────────────────┐       ┌───────────────────┐                     │
│              │   Excel Report    │       │    PDF Report     │                     │
│              │   (10 Sections)   │       │  (ReportLab)      │                     │
│              └───────────────────┘       └───────────────────┘                     │
│                                                                                      │
│  ┌──────────────────────────────────────────────────────────────────────────────┐   │
│  │                              OUTPUT DIRECTORY                                 │   │
│  │  ┌─────────────┐ ┌─────────────┐ ┌─────────────┐ ┌─────────────────────────┐│   │
│  │  │  reports/   │ │   logs/     │ │   cache/    │ │  __pycache__/           ││   │
│  │  │  .xlsx .pdf │ │  .log       │ │  .pkl       │ │  compiled .pyc          ││   │
│  │  └─────────────┘ └─────────────┘ └─────────────┘ └─────────────────────────┘│   │
│  └──────────────────────────────────────────────────────────────────────────────┘   │
│                                                                                      │
│  ┌──────────────────────────────────────────────────────────────────────────────┐   │
│  │                          RULE REPOSITORY (rules/)                            │   │
│  │  ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐          │   │
│  │  │Rule01_   │ │Rule02_   │ │Rule03_   │ │Rule04_   │ │Rule05_   │  ...     │   │
│  │  │LogP.txt  │ │HBD_HBA.txt│ │Polarity  │ │Kamlet.txt│ │Viscosity │          │   │
│  │  └──────────┘ └──────────┘ └──────────┘ └──────────┘ └──────────┘          │   │
│  └──────────────────────────────────────────────────────────────────────────────┘   │
│                                                                                      │
└─────────────────────────────────────────────────────────────────────────────────────┘
Directory Structure
text
FastPhase.AI/
│
├── FastPhaseAI.xlsm                 # Excel frontend with VBA macros
├── main.py                          # Main controller (entry point)
│
├── python/                          # Core Python modules
│   ├── __pycache__/                 # Compiled bytecode cache
│   ├── chromatography_ai.log        # System log file
│   ├── physchem_calculator.py       # Physicochemical property calculator
│   ├── HSMsolute_check.py           # HSM descriptor estimator
│   ├── buffer_selector.py           # 10-rule buffer selection engine
│   ├── solvent_selector.py          # 7-rule solvent selection system
│   ├── column_selector.py           # HSM column selector (USP PQRI)
│   ├── gradient_optimizer.py        # LSS gradient optimization
│   ├── report_template.py           # Excel/PDF report generator
│   ├── reference_manager.py         # Reference database manager
│   └── excel_formatter.py           # Excel styling utilities
│
├── output/                          # Generated outputs
│   ├── reports/                     # Excel and PDF reports
│   │   ├── report__YYYYMMDD_HHMMSS.xlsx
│   │   └── report__YYYYMMDD_HHMMSS.pdf
│   └── logs/                        # Additional log files
├── requirements.txt                 # Python dependencies
├── LICENSE                          # MIT License
└── README.md                        # This file

#Report Structure
The system generates comprehensive 10-section reports in both Excel and PDF formats:

Section	Content
1	Cover Sheet
2	Summary & Confidence Scores
3	Physicochemical Properties
4	HSM Descriptors
5	Column Recommendations
6	Solvent Recommendations
7	Buffer Recommendations
8	Gradient Program
9	References (Vancouver style)
10	Regulatory Alignment (ICH/USP)

#Rule Documentation
All 16 scoring rules are documented as individual text files in the rules/ directory, each containing:

Rule description and mathematical formula

Optimal ranges and thresholds

Peer-reviewed references with DOIs


#Installation
Prerequisites
Python 3.8 or higher

Microsoft Excel (with macros enabled for VBA frontend)

RDKit (conda installation recommended)

Step 1: Clone Repository
bash
git clone https://github.com/yourusername/FastPhase.AI.git
cd FastPhase.AI
Step 2: Install Dependencies
bash
pip install -r requirements.txt
For RDKit (recommended via conda):

bash
conda install -c conda-forge rdkit
Step 3: Verify Installation
bash
python main.py --help
Usage
Method 1: Excel Frontend (Recommended)
Open FastPhaseAI.xlsm in Excel

Enable macros when prompted

Enter SMILES in the DASHBOARD sheet (Cell C5)

Enter compound name, project, and notes

Click "Generate Report" button

Method 2: Command Line
Process a single SMILES:

bash
python main.py --smiles "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O" --name "Ibuprofen" --project "NSAID_Analysis"
Process and generate both Excel & PDF:

bash
python main.py --smiles "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" --name "Caffeine" --format both
Batch processing:

bash
python main.py --batch --excel "FastPhaseAI.xlsm"
Validate SMILES only:

bash
python main.py --validate-only --smiles "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
Process from file:

bash
python main.py --smiles-file compound.txt --format both
Method 3: Python API
python
from main import ChromatographyAIController

# Initialize system
controller = ChromatographyAIController()

# Process a compound
results = controller.process_single_compound(
    smiles="CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    name="Ibuprofen",
    project="Method Development"
)

# Generate reports
report_paths = controller.generate_reports(results, formats=['excel', 'pdf'])
print(f"Excel: {report_paths.get('excel')}")
print(f"PDF: {report_paths.get('pdf')}")
Output Files
Reports
Excel: 10-sheet formatted workbook with conditional formatting, color coding, and hyperlinks

PDF: Professional report with cover page, tables, and reference section

Logs
chromatography_ai.log: Full system log with INFO, WARNING, ERROR levels

Cache
__pycache__/: Compiled Python bytecode for faster execution

Key References
The system implements algorithms and scoring rules from the following peer-reviewed sources:

Reference	DOI / URL	Application
Snyder et al., JCA 2004	10.1016/j.chroma.2004.08.121	HSM Foundation
Dolan et al., JCA 2004	10.1016/j.chroma.2004.09.020	Fs Factor, κ′
Marchand et al., JCA 2008	10.1016/j.chroma.2007.11.101	η′ Estimation
Marchand et al., JCA 2005	10.1016/j.chroma.2004.11.014	σ′ Estimation
Abraham, Chem Soc Rev 1993	10.1039/CS9932200073	α′, β′ Scales
Valkó, JCA 2004	10.1016/j.chroma.2004.01.007	LogP Matching
Vitha & Carr, JCA 2006	10.1016/j.chroma.2006.05.023	HBD/HBA
Goldberg et al., JPCRD 2002	10.1063/1.1416902	Buffer pKa
Kebarle & Tang, Anal Chem 1993	10.1021/ac00070a001	MS Volatility
USP PQRI Database	apps.usp.org/app/USPNF/columnsDB.html	Column Data
Complete reference list with all 30+ citations is available in the references/ directory.

Regulatory Alignment
FastPhase.AI implements recommendations from:

ICH Q2(R2) — Analytical Procedure Validation

ICH Q8(R2) — Pharmaceutical Development (Design Space)

ICH Q9 — Quality Risk Management

ICH Q14 — Analytical Procedure Development

USP <621> — Chromatography

USP <1220> — Analytical Procedure Lifecycle

Troubleshooting
"RDKit not available"
bash
conda install -c conda-forge rdkit
"No module named 'reportlab'"
bash
pip install reportlab
Excel macro errors
Ensure FastPhaseAI.xlsm is in the same directory as main.py

Check that Trust Center Settings allow macros

Verify Python path in VBA matches your installation

PDF generation fails
System falls back to text file automatically

Install reportlab for full PDF support: pip install reportlab

License
MIT License — see LICENSE file for details.

Citation
If you use FastPhase.AI in your research, please cite:

text
FastPhase.AI: AI-Assisted HPLC Method Development System.
https://github.com/yourusername/FastPhase.AI
Acknowledgments
Snyder, Dolan, Carr — Hydrophobic Subtraction Model

USP PQRI — Column Equivalence Database

RDKit — Open-source cheminformatics

All peer-reviewed authors whose work made this system possible
