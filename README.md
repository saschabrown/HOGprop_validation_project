# HOGprop validation, workpackage 1

Set `WP1_ROOT` to the absolute outer workspace path in the terminal that starts Jupyter or a stage command. The outer workspace contains this checkout, `input_raw/` for references to versioned originals, and `working_directory/` for temporary execution files. Run directories belong under each stage's Git-excluded `data/` and `results/`. Use unique run IDs and preserve imported historical outputs.

| Stage | Inputs and outputs |
| --- | --- |
| [01 FastOMA](analysis/01_fastoma/README.md) | Versioned OMA/CAFA files and reference databases → prepared inputs, species tree and original FastOMA outputs |
| [02 FastOMA output processing](analysis/02_fastoma_output_processing/README.md) | Original FastOMA XML and taxon mapping → sister analysis, pruned OrthoXML and exclusions |
| [03 Branch-length estimation](analysis/03_branch_length_estimation/README.md) | Pruned trees and matching source protein sequences → alignments and fitted trees |
| [04 HOGprop](analysis/04_hogprop/README.md) | Tree inputs and frozen GO data → masked inputs, predictions and evaluation |

[Environment status](environments/README.md) records known versions and remaining setup work. Legacy notebooks and launchers are source and historical records. See the stage instructions for execution status before using them. Sascha runs GenomeDK steps himself.
