# Stage 01 · FastOMA

The modified FastOMA workflow is at [software/fastoma](../../software/fastoma), with upstream attribution and license alongside it.

Run the [OMA/CAFA preprocessing notebook](notebooks/Oma_cafa_getting_fastafiles_min15_updated.ipynb) to prepare the input for FastOMA, then use the [input validation notebook](notebooks/Data_validation.ipynb). Both notebook sources use an inclusive 10-aa cutoff. Older saved results may reflect 15 aa. The exploratory protein-count and cutoff notebook is also in `notebooks/`. For future notebook runs, set `WP1_ROOT` and `WP1_PREPARED_RUN_ID` before opening Jupyter. Use a new run ID for preprocessing and the same ID when validating its outputs.

The existing FastOMA run is available locally under `results/20261105_101624_7305523/`. The link opens the original report, OrthoXML and FASTA exports in the sibling `HOGprop` workspace.

## Input references for future preprocessing

GenomeDK holds the raw inputs. Under `$WP1_ROOT/input_raw/`, use links with the following layout. The laptop needs the same input files only if running preprocessing there.

| Relative path | Contents |
| --- | --- |
| `cafa6/` | CAFA 6 dataset with `Train/train_sequences.fasta`, `Train/train_terms.tsv`, `Test/testsuperset.fasta` and `Test/testsuperset-taxon-list.tsv` |
| `oma_may2026/` | OMA export with `oma_species.tsv` and per-genome `*.fa` files |
| `reference/` | `oma-species.txt`, `oma-uniprot.May2026.txt.gz`, `taxa.sqlite` and `LUCA.h5` |

Preprocessing copies the source species mapping into its prepared run before appending mappings. For the current project, proceed from the completed FastOMA results in Stage 02.

## FastOMA launch on GenomeDK

From a GenomeDK compute terminal with the intended FastOMA environment activated and `WP1_ROOT` exported, Sascha can run the command below once the launcher has been validated there. Prepare the launcher with the cluster's Slurm resources if submitting it as a batch job.

```bash
bash "$WP1_ROOT/HOGprop_validation_project/analysis/01_fastoma/scripts/run_fastoma.sh" PREPARED_RUN_ID NEW_RUN_ID
```

The prepared run must contain `fastoma_input/proteome/`, `fastoma_input/species_tree.nwk` and `fastoma_input/LUCA.h5` under `data/PREPARED_RUN_ID/`. The launcher retains the existing Slurm profiles and 10-aa setting. It writes FastOMA outputs to `data/NEW_RUN_ID/`, reports to `results/NEW_RUN_ID/` and Nextflow state to `working_directory/01_fastoma/NEW_RUN_ID/`. The shell launcher starts a fresh run. Historical launch and resume scripts in `scripts/legacy/` document the prior setup and reference its paths.
