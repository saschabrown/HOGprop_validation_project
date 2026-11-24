# Stage 02 · FastOMA output processing

Inputs are the original Stage 01 FastOMA OrthoXML and the OMA-to-NCBI taxon mapping (`irene_cafa_oma_fasta_pairs.tsv` in the repository root). The workflow selects CAFA-training root HOGs, writes an unpruned OrthoXML subset, counts direct same-taxon sisters, prunes selected OMA sisters at speciation-like parents, retains unary nodes and exports uncovered proteins from excluded singleton HOGs and unassigned inputs.

The [sister analysis notebook](notebooks/cafa_train_hog_sisters_reference.ipynb) uses the same paths on both machines. Set `WP1_ROOT` to the absolute workspace location and `WP1_STAGE02_RUN_ID` to a new run name before starting the kernel. Outputs go under `results/<run_id>/`.

| Path | Required input |
| --- | --- |
| `../01_fastoma/results/20261105_101624_7305523/FastOMA_HOGs.orthoxml` | Completed FastOMA OrthoXML |
| `data/reference/all_hog_tree_statistics.tsv` | Saved Step 5 census |
| `data/reference/pruning/summary.json` and `pruning_assignments.tsv` | Saved pruning evidence and assignments |

The laptop reference links point to `Older analysis/`. Set up equivalent references on GenomeDK before running Stage 02 there. This notebook replays those saved inputs and writes new outputs. Census and assignment regeneration is pending. Its final check uses the separate `software/hogprop/` checkout or HOGprop installed in the selected Python environment. See the [extractor notes](scripts/orthoxml_subset_README.md) for the parser dependency.

## Complete-HOG NHX trees

[From complete training HOGs to pruned gene trees](notebooks/cafa_train_hog_sisters_nhx.ipynb) regenerates the root census and sister assignments from Stage 01. It uses Irene’s complete-HOG filter and modified NHX converter, removes the selected same-taxon OMA direct sisters and all CAFA-test tips, collapses all unary groups, and records singleton families separately. It writes a new run under `results/` with exact original root IDs in the filenames. NHX tip labels are XML gene IDs. The leaf audit maps them to protein IDs and source species.

On either machine use a kernel matching [the Stage 02 environment](../../environments/stage02_sister_nhx.yml). In a terminal at the workpackage parent, set `export WP1_ROOT="$(pwd -P)"` before starting Jupyter. Use `WP1_STAGE02_MODE=fixture` for the bundled five-root example or `stage01` (default) for the declared FastOMA run. Each kernel needs the same modified `orthoxml_tools_src` under `$WP1_ROOT/local_config/irene_tree_reference.*/`. The notebook checks the Python-source fingerprint. A full run requires the Stage 01 XML, its sibling `RootHOGsFasta/` directory and the repository-root OMA taxon map. For comparisons, check the input, code and environment entries in `provenance.json` first, then compare table checksums and the per-HOG NHX checksums in `hog_counts.tsv`.

### Converter setup for Filtering_orthoxml.ipynb

Set `WP1_ROOT` to the parent workpackage directory before launching Jupyter. Use the same Stage 02 kernel on both machines (see `environments/stage02_filtering_ete.yml` in the repository). Both NHX conversion and pruned export read Irene’s converter from `$WP1_ROOT/local_config/orthoxml-tools-src`. This directory may be a symlink to a matching orthoxml-tools source checkout containing `cli.py` and the `orthoxml/` package. On GenomeDK, run once in the desktop terminal with WP1_ROOT set:

```bash
mkdir -p "$WP1_ROOT/local_config"
ln -sT /faststorage/project/evomics_projects/Programs/githubscripts/orthoxml-tools/src "$WP1_ROOT/local_config/orthoxml-tools-src"
```

The link command preserves an existing destination. Check the target with `ls -l "$WP1_ROOT/local_config/orthoxml-tools-src/cli.py"`.
