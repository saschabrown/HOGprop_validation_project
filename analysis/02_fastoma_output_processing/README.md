# Stage 02 · FastOMA output processing

Inputs are the original Stage 01 FastOMA OrthoXML and the OMA-to-NCBI taxon mapping (`irene_cafa_oma_fasta_pairs.tsv` in the repository root). The workflow selects CAFA-training root HOGs, writes an unpruned OrthoXML subset, counts direct same-taxon sisters, prunes selected OMA sisters at speciation-like parents, retains unary nodes and exports uncovered proteins from excluded singleton HOGs and unassigned inputs.

The [sister analysis notebook](notebooks/cafa_train_hog_sisters_reference.ipynb) uses the same paths on both machines. Set `WP1_ROOT` to the absolute workspace location and `WP1_STAGE02_RUN_ID` to a new run name before starting the kernel. Outputs go under `results/<run_id>/`.

| Path | Required input |
| --- | --- |
| `../01_fastoma/results/20261105_101624_7305523/FastOMA_HOGs.orthoxml` | Completed FastOMA OrthoXML |
| `data/reference/all_hog_tree_statistics.tsv` | Saved Step 5 census |
| `data/reference/pruning/summary.json` and `pruning_assignments.tsv` | Saved pruning evidence and assignments |

The laptop reference links point to `Older analysis/`. Set up equivalent references on GenomeDK before running Stage 02 there. This notebook replays those saved inputs and writes new outputs. Census and assignment regeneration is pending. Its final check uses the separate `software/hogprop/` checkout or HOGprop installed in the selected Python environment. See the [extractor notes](scripts/orthoxml_subset_README.md) for the parser dependency.
