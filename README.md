# assemble-cf-isolates
Snakemake pipeline to assemble Oxford Nanopore Technology reads. 
This pipeline was originally applied to P. aeruginosa isolates, 
however the pipeline will also work for other isolates. 
You will need Guppy demultiplexed and basecalled FASTQ files to begin. 

By default, this workflow returns [Homopolished](https://github.com/ythuang0522/homopolish) assemblies. 
Isolates with long and short read data can also be assembled into 
[Polypolished](https://github.com/rrwick/Polypolish) assemblies. 

Prior to polishing with Homopolish or Polypolish, this workflow:
- Filters reads using [Filtlong](https://github.com/rrwick/Filtlong);
- Assembles filtered reads using [Flye](https://github.com/fenderglass/Flye);
- Polishes the Flye assembly using [Medaka](https://github.com/nanoporetech/medaka). 

If using Polypolish, 
short reads are filtered using [Fastp](https://github.com/OpenGene/fastp). 

If desired, the assemblies can be compared to a reference genome using 
[QUAST](https://quast.sourceforge.net/). 

## Setup

### Cloning the repository
In order to start, clone this repository to your local environment 
and navigate into the new directory:

```
git clone https://github.com/ezherman/assemble-cf-isolates
cd assemble-cf-isolates
```

### Creating a Snakemake conda environment
You will need a Snakemake conda environment to begin. To speed up the installation of conda environments,
I recommend including `mamba` in this environment:

```
conda create -c conda-forge -c bioconda -n snakemake snakemake mamba -y
conda activate snakemake
```

### Installing the conda environments used by the workflow
When running the workflow for the first time, Snakemake takes care of
creating the required conda environments (with the exception of one environment, see below). 
You can ask Snakemake
to specifically set up the conda environments with the command below.
This way, the duration of the first execution of the workflow will not
be extended by conda environment installations:

```
snakemake -j1 --use-conda --conda-create-envs-only
```

Regardless of whether you run the command above, you will need to set
up the Homopolish conda environment with the commands below. Homopolish
cannot be installed currently by snakemake due to a dependencies issue
(see [here](https://github.com/ythuang0522/homopolish/issues/57)), 
so this is a temporary workaround:

```
conda config --set channel_priority flexible
conda env create -f workflow/envs/homopolish.yml
conda config --set channel_priority strict
```


### Downloading the homopolish sketch
Homopolish uses a mash sketch, which needs to be downloaded.
This download may take a while, as the gzipped file is 3.9 GB. 
You can ask Snakemake to download the file as follows:

```
snakemake -j1 --use-conda download_homopolish_sketch
```

### Editing the config file

#### Guppy version
You will need to specify the version of Guppy that was used for basecalling.
You can edit the `guppy_version` line in the config file (`config/config.yaml`).

#### Medaka model
Additionally, depending on the version of Guppy that was used, 
as well as the Guppy model used in basecalling, you will need to specify
a particular Medaka model. There are details on selecting the appropriate Medaka
model [here](https://github.com/nanoporetech/medaka#models). And you can
find a full list of the latest and older Medaka consensus models 
[here](https://github.com/nanoporetech/medaka/blob/master/medaka/options.py).
Once you have identified the appropriate Medaka model,
change the `medaka_model` section of the `config.yaml` file as appropriate.

For example, by default the model is set to `r941_min_sup_g507`. This model
is appropriate for MinIon runs, with a R9.4.1 flowcell, basecalled with 
Guppy >= 5.0.7 using the Super-Accuracy model. If your FASTQ files
were generated with Guppy 5.0.0 using the Fast model, then 
you should change the model to `r941_min_fast_g303` in the `config.yaml` file.
Failure to change this setting appropriately may result in Medaka doing
a bad job at polishing your Flye-generated assemblies. 

#### Homopolish model
Finally, you will need to specify the model that is used by Homopolish. 
By default, this is set to `R9.4.pkl`. If your sequencing involved an
`R10` flowcell, set the `homopolish_model` parameter to `R10.3.pkl`. 

### Setting up your input data
Each isolate should have its own directory inside `data/long-read`, named
to identify that isolate (e.g. `barcode02`). This
directory should contain the *gzipped* FASTQ file(s) for that isolate.
If you also have short-read data for this isolate, the files can be
stored in `data/short-read`. For example, if five long-read FASTQ files
were associated with `barcode02`, as well as two short-read FASTQ files,
then the data directory would have the following structure:

```
data
├── long-read
│   └── barcode02
│       ├── fastq_runid_afa204b7accb9021bb58baa5b294dc1a88ee3e9e_110_0.fastq.gz
│       ├── fastq_runid_afa204b7accb9021bb58baa5b294dc1a88ee3e9e_111_0.fastq.gz
│       ├── fastq_runid_afa204b7accb9021bb58baa5b294dc1a88ee3e9e_112_0.fastq.gz
│       ├── fastq_runid_afa204b7accb9021bb58baa5b294dc1a88ee3e9e_114_0.fastq.gz
│       └── fastq_runid_afa204b7accb9021bb58baa5b294dc1a88ee3e9e_117_0.fastq.gz
└── short-read
    ├── barcode02_1.fastq.gz
    └── barcode02_2.fastq.gz
```

## Running the workflow

### Default behaviour: homopolished assemblies
By default, the workflow will create Homopolished assemblies for
all isolate directories found in `data/long-read`. This can be 
achieved with the following command:

```
snakemake -j1 --use-conda
```

Note that you can increase `-j1` to the number of cores available to you,
e.g. `-j8` if you have eight available cores.

You can also explicitly specify that you want Homopolished assemblies:

```
snakemake -j1 --use-conda homopolish_workflow
```

### Polypolished assemblies
If you have short-read data for some or all of your isolates,
you can request polypolished assemblies as follows:

```
snakemake -j1 --use-conda polypolish_workflow
```

### Both homopolished and polypolished assemblies
You can obtain both homopolished and polypolished assemblies
from the same Snakemake run. Note however that this does not
return assemblies that have been both homopolished and polypolished.
Instead, this command runs both of the workflows used above:

```
snakemake -j1 --use-conda homopolish_workflow polypolish_workflow
```

## Running the workflow on a SLURM cluster
This workflow can be run on a SLURM cluster using the scheduler
and profile in the `slurm` directory.

### Setup for running the workflow on a SLURM cluster
The following features need to be specified in the scheduler
script, `slurm/scheduler.sh`:

- `--mail-user`, which will ensure that you are notified by email
when the workflow ends or fails. `--mail-type` and `--mail-user`
can be removed to disable this feature.
- `--account`, for your project account.

The following features need to be specified in the profile config
file, `slurm/slurm_profile/config.yaml`. Note that this config file
functions alongside the workflow config file, `config/config.yaml`.

- `jobs`, which defines the maximum number of cores that snakemake
can request at one time. 
- If you would like to be notified of the failure of individual
jobs (i.e. individual rule executions), you can add the following
to the quoted section of the `cluster` line: 
`--mail-type=FAIL --mail-user youremail@email.com`. 

### SLURM cluster execution
After navigating to the `slurm` directory, run the following
to execute the homopolish workflow:

```
sbatch scheduler_homopolish.sh
```

Snakemake will create a `snakemake_%j.log` log file in the `slurm` 
directory, through which you can monitor the progress of the workflow.
Log files for individual jobs will be written to the `logs_slurm` directory. 

To run the polypolish workflow, run the following:

```
sbatch scheduler_polypolish.sh
```

To run both workflows, run the following:

```
sbatch scheduler_all.sh
```