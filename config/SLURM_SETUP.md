# GenomeDK Slurm setup

Prepared: 5 November 2026. Status: local configuration prepared; cluster validation and submission not performed.

## Configuration

Use `config/fastoma_genomedk_slurm.config` with this repository's `FastOMA/FastOMA.nf` and the `slurm,large` profiles. Do not additionally load `config/cafa_oma.config`: its validator override still names the historical `FastOMA_cafa_custom/python` directory.

The new configuration sets:

| Setting | Value |
| --- | --- |
| Task executor | Slurm |
| Account for each task | `evomics_projects` (previously used account; current access unverified) |
| Maximum submitted tasks | 200, matching the saved Slurm profile |
| Minimum protein length | 10 residues |
| Input-check memory | 32 GB initially, multiplied by attempt number |
| Per-task resource ceilings | 32 CPUs, 256 GB, 72 hours; memory cap selected by Sascha; existing CPU/time ceilings |
| Python source | The FastOMA directory containing the selected workflow |
| Partition | Not explicitly selected |

The existing `large` profile supplies the per-label resource settings: single tasks initially request 1 CPU/16 GB/8 hours, medium tasks 6 CPUs/36 GB/12 hours, and high tasks 12 CPUs/128 GB/24 hours. Process-specific rules override these label settings: OMAmer memory depends on database size and rootHOG inference memory depends on the number of proteomes. These are existing formulas, not requirements measured for this dataset. Resource ceilings apply to individual tasks, not total simultaneous resource use.

The account and Slurm executor are applied to the pipeline tasks themselves. If Nextflow is launched inside a batch job, that launcher also needs its own account and resources. The launcher must remain running while it coordinates the task jobs.

## Path and environment checks still needed

- Confirm the current full project path on GenomeDK and the Conda environment containing Nextflow, the FastOMA console commands, OMAmer, and the alignment/tree tools.
- Confirm the appropriate partition, account access, resource ceilings, and task concurrency. The existing large-run limits are a starting point, not verified cluster limits.
- Confirm that the prepared FASTAs and species tree were regenerated using the new 10-residue notebook setting. Editing notebook source did not regenerate inputs.
- Confirm the selected database file and distinct output/work directories for this run. The workflow default currently points to the shared `All.May2026/LUCA.h5` file; an explicit `--omamer_db` selects a different file.

Updated `scripts/job_fastoma_large_dataset.sh` to use `/home/saschabrown/evomics_projects/saschabrown/HOGprop_validation_project`, the repository's `FastOMA/FastOMA.nf`, and the new Slurm configuration. It retains the previously recorded `fastoma_omamer212` Conda environment, which has not been reconfirmed. The launcher requests 1 CPU, 72 GB, and 72 hours; these are separate from the per-task requests. It uses the prepared `fastoma_input/LUCA.h5` explicitly. Output and work directories receive a timestamp/job-ID suffix to preserve previous runs. Queue size comes from the configuration (200); the earlier command-line override of 550 is removed.

Before submission, Sascha must create `logs/slurm` under the project directory: Slurm opens its log files before the script runs. No submission is requested at this stage.

## First check for Sascha to run

After confirming the current project path and activating the chosen environment, run this in a GenomeDK terminal **from the repository root**. This resolves and prints configuration; it does not launch the pipeline or submit Slurm tasks. Nextflow may initialize or download its configured plugins or included configuration during setup.

```bash
nextflow -c ./config/fastoma_genomedk_slurm.config config ./FastOMA -profile slurm,large -flat
```

Share the relevant `process.executor`, `clusterOptions`, `queueSize`, `resourceLimits`, `min_sequence_length`, `omamer_db`, and `beforeScript` settings before preparing the final batch launch command. Do not share credentials. Do not submit until the environment, regenerated inputs, and resolved configuration have been checked.

## Validation and provenance

The configuration was prepared from the local `FastOMA/conf/base.config`, `base_large.config`, `slurm_basic.config`, `FastOMA/nextflow.config`, and existing CAFA configuration and launch script. Local static review checked the task executor, account, source path, cutoff, and resource settings. Nextflow is unavailable in this local environment, so configuration parsing and merging have not been tested. No workflow, notebooks, remote sessions, or Slurm jobs were run.

The task executor and supported Slurm directives follow the [Nextflow Slurm documentation](https://docs.seqera.io/nextflow/reference/executor/slurm). The project account follows [GenomeDK's queue guidance](https://genome.au.dk/docs/interacting-with-the-queue/). Neither source establishes the current limits of Sascha's account.
