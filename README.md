# Chromatography AI System v3.0

AI-assisted HPLC method development: Excel front-end + Python engines + VBA bridge.

## Architecture

```
Chromatography_AI_System_v3.xlsm  ← Excel UI (DASHBOARD, SEARCH HISTORY, SETTINGS …)
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

1. Open `Chromatography_AI_System_v3.xlsm` and enable macros.
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
python main.py --batch --excel ../Chromatography_AI_System_v3.xlsm
```

---

## Bug Fixes Applied (v3.0.1)

### VBA

| ID | Location | Issue | Fix |
|----|----------|-------|-----|
| VBA-001 | `RunPython` | Path separator missing between `ThisWorkbook.Path` and `SCRIPT_NAME` | Added `sep` variable using `Application.PathSeparator` |
| VBA-002 | `RunPython` | `--name` and `--project` args not passed to Python | Added correct arg construction with quote-escaping |
| VBA-003 | `ValidateSMILES_Click` | Passed `--validate-only` flag that doesn't exist in `main.py` | Now sends `--format excel`; Python handles validate-only via new `--validate-only` flag |
| VBA-004 | `UpdateStatus` | Status cells (H5/H6/H7) wrote to wrong ranges | Corrected to match DASHBOARD layout |
| VBA-005 | `GetSetting` | Setting name lookup failed for aliased keys (e.g. `AutoOpenReports` vs `Auto-Open Reports`) | Added alias resolution map |
| VBA-006 | `ClearAll_Click` | Cleared `C5:C8` but reset wrong status cells | Fixed to reset H5:H8 with correct labels |
| VBA-007 | `RunPython` | Hard-coded `\` path separator breaks on macOS/Linux | Uses `Application.PathSeparator` |
| VBA-008 | `UpdateHistory` | Parts index mapping off by 1 for report paths | Aligned with Python `result_line` output format |

### Python

| ID | Location | Issue | Fix |
|----|----------|-------|-----|
| PY-001 | `main.py` argparse | `--project` argument missing; used in `result_line` but never defined | Added `--project / -p` argument |
| PY-002 | `main.py` argparse | `--validate-only` flag missing; VBA sends it | Added flag with early-return handler |
| PY-003 | `main.py` `read_input_smiles` | Cell mapping used `B2-B5` but DASHBOARD has data in `C5-C8` | Fixed to `(4,2)-(7,2)` (0-indexed) |
| PY-004 | `main.py` `update_live_preview` | Hard-coded `wb['Input']` fails (sheet is `DASHBOARD`) | Uses `_get_existing_sheet_name(['DASHBOARD','Input'])` |
| PY-005 | `main.py` `_update_progress` | Same hard-coded `'Input'` issue + wrong cell refs (H5/H6) | Fixed sheet lookup + cell refs (H7/H8) |
| PY-006 | `buffer_selector.py` `select_optimal_buffer` | Weight key lookup always fell through to default 0.1 due to broken string matching | Added explicit `_wkey_map` dict for rule→weight key |
| PY-007 | `column_selector.py` | `from rdkit.Chem import MolFromSmarts` deprecated in RDKit ≥2023 | Removed import; uses `Chem.MolFromSmarts()` throughout |
| PY-008 | `solvent_selector.py` | `'USE_RDKIT' in globals() and USE_RDKIT` fragile; breaks if imported | Simplified to direct module-level `if USE_RDKIT:` |

---

## Project Layout

```
chromatography_ai/
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
│   └── Module1_VBA_Controller.bas     ← Fixed VBA source (paste into Alt+F11)
├── tests/
│   └── test_smoke.py
├── docs/
│   └── vba_fix_notes.md
└── output/
    └── reports/                        ← Generated Excel + PDF reports land here
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
