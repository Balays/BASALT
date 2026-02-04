# BASALT Pipeline Overview

BASALT is a modular metagenomic binning and refinement pipeline.  
It is organised as a series of numbered steps (S1–S10) plus a few entry-point
and utility modules.

- Main CLI: `BASALT.py`
- CheckM2 branch (default): the recommended modern workflow
- CheckM branch: legacy / compatibility workflow

Each step below lists the main script(s) and core inputs/outputs.

---

## 1. Top-level Entry Point

**File**: `BASALT.py`  
**Role**: Command-line interface; parses arguments and dispatches to appropriate modules.

**Key responsibilities**

- Parse user arguments:
  - Assemblies (`--assemblies`)
  - Short reads (`--shortreads`)
  - Long reads (`--longreads`)
  - HiFi reads (`--hifi`)
  - Threads / RAM (`--threads`, `--ram`)
  - Functional modules (`--module`: `autobinning`, `refinement`, `reassembly`, `all`)
  - Quality check backend (`--quality-check`: `checkm2` or `checkm`)
  - Data feeding / additional binsets / coverage lists
- Normalise inputs into Python structures:
  - `assembly_list`
  - `datasets` (PE reads)
  - `lr_list`, `hifi_list`, `hic_list`
  - `data_feeding_folder`, `binsets_list`, `coverage_list`
- Select branch:
  - **CheckM2**:
    - Data feeding only → `Data_feeding.data_feeding`
    - External binsets derep only → S4 / S5 CheckM2
    - Normal pipeline → `BASALT_main_d.BASALT_main_d`
  - **CheckM**:
    - Data feeding only → `BASALT_main_c_datafeeding.data_feeding_main`
    - External binsets derep only → S4 / S5 CheckM
    - Normal pipeline → modular CheckM entry points:
      - `BASALT_main_c_autobinning`
      - `BASALT_main_c_refinement`
      - `BASALT_main_c_re_assembly`

---

## 2. Core Entry Modules

### 2.1 CheckM2 monolithic main

**File**: `BASALT_main_d.py`  
**Function**: `BASALT_main_d(...)`

- Runs the full pipeline in one function:
  - S1–S3: autobinning + within-group dereplication
  - S4: multi-assembly dereplication
  - S5: DL-based outlier removal
  - S6–S7: contig retrieval + second dereplication
  - S8: OLC-based elongation
  - S9–S10: reassembly + post-reassembly OLC
- Manages:
  - `Basalt_checkpoint.txt` (resume)
  - `Basalt_log.txt`, `BASALT_command.txt`
  - Final bin renaming and final CheckM2 report

### 2.2 CheckM modular mains

**Files**:

- `BASALT_main_c.py` (legacy monolithic CheckM main)
- `BASALT_main_c_autobinning.py` (`BASALT_main_c_autobinning`)
- `BASALT_main_c_refinement.py` (`BASALT_main_c_refinement`)
- `BASALT_main_c_re_assembly.py` (`BASALT_main_c_re_assembly`)
- `BASALT_main_c_datafeeding.py` (`data_feeding_main`)

These provide the same logical pipeline as `BASALT_main_d`, but split by
phase (autobinning / refinement / reassembly / data feeding) and using
CheckM outputs instead of CheckM2.

---

## 3. Step S1 – Autobinning

**Files**:

- `S1_Autobinners_2qc_11152023.py`
- `S1e_extra_binners.py`
- `S1p_Merging_bins_within_group.py`

**Purpose**

- Map reads to assemblies.
- Run multiple binners (MetaBAT2, MaxBin2, CONCOCT; optionally MetaBinner, VAMB, LorBin).
- Merge candidate bin sets within each assembly group.

**Key inputs**

- Assembly FASTA(s)
- PE reads (`datasets`)
- Optional long reads (`lr_list`, `hifi_list`)
- Sensitivity / RAM / threads

**Key outputs**

- Per-assembly bin folders (`*_genomes`)
- Coverage matrices
- PE connection summaries
- Intermediate dictionaries written as text files (connections, depth, bin folders).

---

## 4. Step S2 – Bin abundance & PE connections

**Files**:

- `S2_BinsAbundance_PE_connections_multiple_processes_pool_10032023.py` (CheckM2)
- `S2_BinsAbundance_PE_connections_multiple_processes_pool_checkm.py` (CheckM)

**Purpose**

- Compute coverage matrices and abundance for each bin.
- Integrate PE connections to better characterise bin structure.

**Inputs**

- Bin folders from S1
- Coverage depth files (`*.depth.txt`)
- PE connection tables

**Outputs**

- Combined coverage matrices per bin group.
- PE connection summary files and in-memory dictionaries.

---

## 5. Step S3 – Within-group comparison & selection

**Files**:

- `S3_Bins_comparator_within_group_10042023.py` (CheckM2)
- `S3_Bins_comparator_within_group_checkm.py` (CheckM)

**Purpose**

- Compare bins within each assembly group.
- Use coverage, contig overlap, and quality (CheckM/CheckM2) to pick best bins.
- Filter redundancies **within** each group.

**Outputs**

- `BestBinsSet` folders per assembly.
- Coverage matrix lists and assembly lists for downstream S4.

---

## 6. Step S4 – Multi-assembly comparison & dereplication

**Files**:

- `S4_Multiple_Assembly_Comparitor_multiple_processes_bwt_10242023.py`
- `S4_Multiple_Assembly_Comparitor_multiple_processes_bwt_checkm.py`

**Purpose**

- Compare best binsets across multiple assemblies.
- Use coverage, GC content, contig overlap and CheckM/CheckM2 metrics to:
  - cluster similar bins across assemblies
  - select canonical representatives
- Runs multiple dereplication phases: initial, secondary, final.

**Outputs**

- `BestBinset` and intermediate dereplicated binsets.
- Coverage matrix lists and binset lists used by later steps.

---

## 7. Step S5 – DL-based outlier remover

**Files**:

- `S5_Outlier_remover_DL_11012023.py` (CheckM2)
- `S5_Outlier_remover_DL_checkm.py` (CheckM)
- DL support: `ensemble.py`, `model.py`, `my_dataset.py`, `utils.py`, `BASALT_models_download.py`

**Purpose**

- Compute TNF + coverage features for contigs.
- Run deep-learning models to predict outlier (contaminated) contigs.
- Remove outliers to refine bins.

**DL components**

- `model.py`:
  - `LBR`: Linear–BatchNorm–ReLU residual block.
  - `MLP`: multi-layer perceptron classifier.
- `my_dataset.py`:
  - Dataset loaders for train/val/test contamination feature matrices.
- `ensemble.py`:
  - Loads multiple MLP checkpoints and performs ensemble voting.
- `utils.py`:
  - Model weight download, normalisation, training loops, focal loss, etc.
- `BASALT_models_download.py`:
  - Downloads ensemble CSV + checkpoint files into `BASALT_WEIGHT` cache.

**Outputs**

- Refined binset (`BestBinset_outlier_refined`).
- Updated quality reports.

---

## 8. Step S6 – Contig retrieval from PE contigs

**Files**:

- `S6_retrieve_contigs_from_PE_contigs_10302023.py`
- `S6_retrieve_contigs_from_PE_contigs_checkm.py`
- `S6p_coverage_filtration_mpt_06102022.py`
- `S6p_Run_checkm_taxonomy_wf.py`

**Purpose**

- Use PE connections and coverage patterns to recruit contigs that:
  - were not originally in bins
  - but strongly connect to them
- Optionally filter using CheckM taxonomy.

**Outputs**

- Updated binsets with recruited contigs.
- Coverage- and TNF-based filtration reports.

---

## 9. Step S7 – Within-group contig retrieval & refinement

**Files**:

- `S7_Contigs_retrieve_within_group_10262023.py`
- `S7_Contigs_retrieve_within_group_checkm.py`
- `S7lr_finding_sr_contigs_basing_lr_and_polishing_11022023.py`
- `S7lr_finding_sr_contigs_basing_lr_and_polishing_checkm.py`
- `S7p_Gap_filling2.py`

**Purpose**

- Further within-group contig retrieval based on:
  - coverage consistency
  - PCA / IQR outlier tests
  - PE and long-read connections
- Optional polishing:
  - Short-read polishing (Pilon-style)
  - Long-read polishing
- Gap filling and contig merging for improved N50.

**Outputs**

- Refined binsets (`BestBinset_outlier_refined_filtrated_retrieved*`).
- Polished bins (MAGs) and long-read-merged bins.
- Gap-filled bin FASTAs.

---

## 10. Step S8 – OLC-based reassembly

**Files**:

- `S8_OLC_new_10262023.py`
-, `S8_OLC_new_checkm.py`

**Purpose**

- Use overlap-layout-consensus (OLC) strategy to elongate bin contigs.
- Key operations:
  - Select contigs for elongation (coverage/TNF outliers).
  - BLAST-based similarity grouping (`blast_1`).
  - Sequence merging (`seq_merge`, `elongation_main`).
  - CheckM2/CheckM-informed bin comparison (`merge`, `bin_comparison*`).

**Outputs**

- OLC-elongated bins (`*_OLC` folders).
- Updated CheckM2/CheckM reports after OLC.

---

## 11. Step S9 – Short-read and hybrid reassembly

**Files**:

- `S9_Reassembly_10262023.py`, `S9_Reassembly_checkm.py`
- `S9p_Hybrid_Reassembly_10262023.py`, `S9p_Hybrid_Reassembly_checkm.py`

**Purpose**

- Reassemble each bin using:
  - Short reads (SPAdes, IDBA-UD)
  - Optionally long reads (hybrid SPAdes)
- Evaluate quality (CheckM2 or CheckM) and pick best assembly per bin.

**Outputs**

- `<binset>_re-assembly_binset` with reassembled bins.
- Hybrid assemblies when long reads are available.
- Updated CheckM stats and comparison reports.

---

## 12. Step S10 – Post-reassembly OLC refinement

**Files**:

- `S10_OLC_new_10262023.py`
- `S10_OLC_new_checkm.py`

**Purpose**

- After reassembly, run another round of OLC-based comparison:
  - Compare reassembled bins to original bins.
  - Extend / merge contigs where safe.
  - Apply coverage/TNF-based outlier detection again.
- Final dereplication across reassembly results.

**Outputs**

- Final OLC-refined binset (`*_re-assembly_OLC*`).
- Final quality reports and dereplicated bin folder moved/renamed
  to user-specified `output_folder`.

---

## 13. Data Feeding & Final Dereplication

### 13.1 Data feeding

**File**: `Data_feeding.py`  
**Used via**:

- `BASALT.py` (CheckM2 branch, direct)
- `BASALT_main_c_datafeeding.py` (CheckM branch)

**Purpose**

- Bring external binsets into BASALT:
  - Rename contigs and bin IDs.
  - Map PE reads, recompute coverage matrices.
  - Compute PE connections and CheckM / CheckM2 quality.
- Generate compatible inputs for S4/S5/S6/S7 on mixed binsets.

### 13.2 Final dereplication helper

**File**: `Final_drep.py`  

- PCA-based final dereplication and reporting utilities.
- Used when running extra or standalone dereplication on final binsets.

---

## 14. Cleanup

**File**: `Cleanup.py`  
**Function**: `cleanup(assembly_list)`

- Archive depth, coverage, connection and comparison files
  into compressed archives.
- Remove intermediate mapping/index files, temporary bin folders,
  and other large intermediates.
- Typically called at the end of the full pipeline
  when `functional_module == 'all'`.

