# FastPhase.AI

AI-assisted HPLC method development: Excel front-end + Python engines + VBA bridge.

## Architecture

```
FastPhase.AI.xlsm  ← Excel UI (DASHBOARD, SEARCH HISTORY, SETTINGS …)
        │  VBA bridge (Module1_VBA_Controller.bas)
        ↓
python/main.py                      ← Orchestrator
        ├── physchem_calculator.py  ← LogP, TPSA, drug-likeness, ionisation
        ├── HSMsolute_check.py      ← HSM descriptors (η′, σ′, β′, α′, κ′)
        ├── column_selector.py      ← HSM-based column ranking & Fs factor
        ├── buffer_selector.py      ← 10-rule buffer scoring
        ├── solvent_selector.py     ← Kamlet-Taft / Abraham solvent scoring
        ├── gradient_optimizer.py   ← LSS-theory gradient optimisation
        └── reference_manager.py   ← Unified reference database (Vancouver citations)
```

## Quick Start

### 1. Install dependencies

```bash
# RDKit — must use conda
conda install -c conda-forge rdkit

# Pure-pip packages
pip install -r requirements.txt --break-system-packages
```

### 2. Set Python path (Windows)

In the Excel **SETTINGS** sheet set **Python Executable Path** to the full path of your
Python interpreter, e.g. `C:\ProgramData\miniconda3\envs\chrom\python.exe`.

### 3. Use the system

1. Open `FastPhase.AI.xlsm` and enable macros.
2. Paste a SMILES string into cell **C5** on the **DASHBOARD** sheet.
3. Fill in Compound Name (C6), Project ID (C7), and optional Notes (C8).
4. Click **▶ GENERATE REPORT** — a Python subprocess runs and results are written
   to `output/reports/` and logged in **SEARCH HISTORY**.

### 4. Command-line usage

```bash
cd python/
python main.py --smiles "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O" --name "Ibuprofen" --format both
```

### 5. Batch processing

```bash
python main.py --batch --excel ../FastPhase.AI.xlsm
```

---
## Project Layout

```
FastPhase.AI/
├── README.md
├── requirements.txt
├── Chromatography_AI_System_v3.xlsm   ← Excel workbook (copy here from original)
├── python/
│   ├── main.py
│   ├── physchem_calculator.py
│   ├── HSMsolute_check.py
│   ├── buffer_selector.py
│   ├── column_selector.py
│   ├── gradient_optimizer.py
│   ├── solvent_selector.py
│   └── reference_manager.py
├── vba/
│   └── Module1_VBA_Controller.bas     ← Fixed VBA source (paste into Alt+F11)                      ← Generated Excel + PDF reports land here
```

---

## VBA Installation (after downloading)

1. Open the `.xlsm` file and press **Alt + F11** to open the VBA editor.
2. In the Project Explorer, expand **VBAProject → Modules**.
3. Double-click **Module1** (or insert a new one if absent).
4. Replace all code with the contents of `vba/Module1_VBA_Controller.bas`.
5. Save and close the VBA editor.
6. Back in Excel, assign each macro to its button via **Developer → Assign Macro**.

---

## References

All analytical rules are traceable to peer-reviewed sources. See the
**REFERENCES_MASTER** sheet in the Excel workbook or `reference_manager.py`
for full Vancouver-style citations.

Key references: Snyder et al. (2004) J. Chromatogr. A 1060:77–116 [HSM] ·
Abraham (1993) Chem. Soc. Rev. 22:73–83 [H-bonding] ·
Wildman & Crippen (1999) J. Chem. Inf. Comput. Sci. 39:868–873 [LogP] ·
Goldberg et al. (2002) J. Phys. Chem. Ref. Data 31:231–370 [Buffer pKa] ·
Valkó (2004) J. Chromatogr. A 1037:299–310 [Solvent selection].
