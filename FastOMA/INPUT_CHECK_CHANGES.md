# FastOMA input-check and parameter changes

Date: 5 November 2026

File: `FastOMA/check_input.py` (relative to this FastOMA directory).

## Requested behavior

- Accept a proteome containing at least one protein record. Empty proteomes still fail.
- Remove the amino-acid composition check, so input is no longer rejected based on the proportion of letters distinguishing protein sequences from DNA/RNA.

## Code changes

- Expressed the minimum record count as `if num_prots < 1` and clarified the error message for an empty or unreadable proteome.
- Deleted the residue counting, amino-acid letter selection, fraction calculation, and rejection branch.
- Preserved all original imports, including `itertools` and the BioPython IUPAC alphabet imports, as requested.

```python
# Each proteome must contain at least one protein record.
if num_prots < 1:
    logger.error("The input proteome contains no protein records or Bio.SeqIO couldn't read it, %s/%s.fa", proteome_folder, species_name)
    isOk = False
```

The working copy already contained an uncommitted change from `num_prots <= 2` to `num_prots <= 0`, and from a 20% composition threshold to 5%. This update preserves the intended one-record minimum and removes the composition check entirely. Removing that calculation also prevents division by zero for an empty proteome.

Duplicate protein IDs and IDs longer than 60 characters still fail. The separate requirement for at least three proteome files, species-tree validation, and other input checks are unchanged. Acceptance here does not establish that a sequence is biologically a protein.

## Local validation

Parsed the complete module for Python syntax and executed its extracted `check_proteome` function with synthetic records. All six checks passed:

| Case | Result |
| --- | --- |
| One protein record | Accepted |
| Two protein records | Accepted |
| One DNA-like sequence (`ACGTN`) | Accepted |
| Zero protein records | Rejected without division by zero |
| Duplicate protein IDs | Rejected |
| Protein ID of 61 characters | Rejected |

This was isolated function validation, not a full FastOMA pipeline run or a test of FASTA parsing. No original datasets were modified and no GenomeDK access was used.

## Minimum protein-length parameters — 5 November 2026

Changed the following settings (paths relative to the repository root):

| Configuration file | Parameter | Previous value | New value |
| --- | --- | --- | --- |
| `FastOMA/nextflow.config` | `params.min_sequence_length` (inside the `params` block) | 20 residues | 10 residues |
| `config/cafa_oma.config` | `params.min_sequence_length` | 15 residues | 10 residues |

The existing filter in `FastOMA/FastOMA/_utils_roothog.py`, function `parse_proteomes`, retains records where `len(rec) >= min_sequence_length`. With either updated setting applied, sequences of exactly 10 residues or longer are retained and sequences shorter than 10 residues are excluded. This is a sequence-length threshold, separate from the minimum number of protein records per proteome.

The standalone `fastoma-infer-roothogs` command still defaults to 50 residues unless a value is explicitly supplied. The Nextflow workflow supplies its configured value to that command. Configuration selection and command-line overrides determine the effective threshold for a particular run.

Validation of the parameter edits: reviewed both configuration diffs and ran `git diff --check` for those files successfully. No pipeline was run and no original datasets were modified. This section documents the earlier parameter changes; updating this Markdown did not rerun that validation.

## Species-tree preparation and validation notebooks — 5 November 2026

Changed the CAFA minimum length from 15 to 10 residues in these notebooks (paths relative to the repository root; cell numbers count all cells starting at one):

- `notebooks/Oma_cafa_getting_fastafiles_min15_updated.ipynb`: cell 14 now sets `minimum_length = 10`. Cell 12 uses `m_len = 10` for the matching short-protein statistics. Updated the introductory text to describe the new cutoff. The existing notebook filename is preserved.
- `notebooks/Data_validation.ipynb`: cell 19 now sets `minimum_length = 10`. Updated related descriptions, report columns, and derived report filenames to use 10, including `cafa_min10_length_validation.tsv` and `cafa_proteins_below_10aa.tsv`. Cell 18's heading also reflects 10.

The existing inclusive comparison retains CAFA sequences of exactly 10 residues or longer. CAFA taxa without retained proteins are omitted from the CAFA tree taxa; OMA sequences remain unfiltered by length in this preparation notebook. New tree membership and counts require a future run and have not been established by this edit.

No notebook cells or pipeline were executed. Saved outputs were preserved and must be treated as historical, not results for the new cutoff. Validation was limited to parsing notebook JSON and changed code syntax, confirming saved outputs were unchanged, checking diffs, and an isolated comparison of lengths 9, 10, and 11 against the inclusive cutoff. Original datasets, the older notebook under `scripts/old/`, and the exploratory `protein_count_and_effect.ipynb` were not modified.

## Slurm configuration — 5 November 2026

Added `config/fastoma_genomedk_slurm.config` (repository-relative path) for the `slurm,large` profiles. It sets the task account to `evomics_projects`, queue size to 200, minimum sequence length to 10, and per-task ceilings to 32 CPUs, **256 GB**, and 72 hours. Sascha selected the 256 GB memory cap. The validator starts at 32 GB; other task resource formulas remain inherited from the existing profiles and are subject to the ceilings. Python imports use this repository's modified FastOMA source rather than the old custom installation.

Updated `scripts/job_fastoma_large_dataset.sh` for the confirmed project path `/home/saschabrown/evomics_projects/saschabrown/HOGprop_validation_project`, local repository workflow, and new configuration. The launcher retains 1 CPU and 72 GB, with a 72-hour walltime; this is separate from the task limits. It retains the previously saved Conda environment `fastoma_omamer212` pending confirmation. Outputs and work directories have timestamp/job-ID suffixes to preserve earlier results. The database remains the explicitly supplied prepared `fastoma_input/LUCA.h5` file.

See `config/SLURM_SETUP.md` for settings, remaining environment/input checks, and a configuration-only command for Sascha to run on GenomeDK. Local shell syntax and diff checks passed. Nextflow is not installed locally, so configuration merging, resource-limit enforcement, and task execution remain untested. Nothing was submitted; notebook inputs were not regenerated.

## Species-tree overwrite and combined input selection — 5 November 2026

In `notebooks/Oma_cafa_getting_fastafiles_min15_updated.ipynb`, cell 53 now explicitly opens the generated `irene_species_tree.nwk` with mode `w` before writing the serialized tree. This replaces an existing generated tree rather than appending. The previous `t.write(..., outfile=outSpeciesFile)` already saved to that path; this change makes replacement explicit. The node-name mapping already uses mode `w`. Cell 61 continues to overwrite the derived `fastoma_input/species_tree.nwk` using `shutil.copy2` after checking the selected source files and existing reference copies.

Updated cell 5's project output root to the confirmed `/home/saschabrown/evomics_projects/saschabrown/HOGprop_validation_project`. Original CAFA and OMA source paths were not changed.

The combined dataset in `fastoma_input/proteome/` is a directory of separate selected FASTAs, not one concatenated FASTA:

- CAFA: only taxa remaining after the 10-residue filter and taxon exclusions, with at least one retained protein; their filtered FASTAs are copied and overwritten on refresh.
- OMA: codes in `tax2oma`, built from the included OMA taxonomy metadata after taxon exclusions. This can include OMA taxa absent from CAFA and multiple OMA genomes per taxon; it is not restricted to the CAFA/OMA intersection. OMA sequences are not length-filtered by this preparation cell. Existing matching OMA copies are reused; missing copies are added and differing copies cause an error.
- The saved tree's leaf labels must exactly match the selected CAFA IDs plus OMA codes. Each tree leaf selects its corresponding source FASTA. Unselected extra source FASTAs are not copied. Obsolete numeric CAFA files are removed from the derived combined folder; unexpected OMA files cause an error rather than deletion.

Validation: inspected selection and refresh source cells; parsed edited notebook JSON and changed code syntax, checked whitespace, and confirmed saved outputs were preserved. No notebook cells were executed, no tree or combined dataset regenerated, and no original datasets modified.
