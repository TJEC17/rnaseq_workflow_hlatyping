# Rnaseqworkflow


## Simple rnaseqworkflow with hlatyping
Simple rnaseqworkflow using trim_galore (with inbuilt fastqc functionality), STAR aligner and Featurecounts. Along with the option of either seq2HLA or Optitype hla typing. Requires both forward and backward fastq reads.


## Requirements 
**Nextflow**
You will need the latest version of nextflow as this pipeline is built in dsl2 and for it to be in your $PATH environment
```
wget -qO- https://get.nextflow.io | bash
```

**seq2HLA**
Seq2HLA was loaded as a module. However if you want this to work you will need to install seq2HLA yourself
and set the path to where that tool is called in the rnaseq2hla.nf script. I recommend either using [bioconda](https://anaconda.org/bioconda/seq2hla) to install seq2HLA, or you can use their [github](https://github.com/TRON-Bioinformatics/seq2HLA).

**Optitype**
Optitype has a variety of prerequisites and packages needed for it to run. I suggest looking at [Optitypes](https://github.com/FRED-2/OptiType/tree/master) github and following the guide on Installation from Source to setup Optitype. Also you can ignore requirement for python2.7, as Optitype works with python3.
Most of the tools can be loaded in via modules. 
You also have to edit the config.ini file, the necessary instructions are found within the file itself

**Important note**
If you have loaded in seq2HLA and wish to use Optitype afterwards, you have to create a new terminal and follow the Optitype steps. You have to set up a new terminal because when you load seq2HLA as a module it sets your python environment to python2 (I know, ancient isn't it), and this will mess with Optitype.


## Installation
Clone the repository:
```
git clone url here
```
Put fastq files (forward and backward reads) into the **fastq_data** folder.
Depending on what hlatyping tool you would like to use, follow the requirements for that tool

## Usage
Before running any of the workflows you will need to generate a genome index and create a fastq_data directory (this is where you will put all your fastq files). To do this run run these commands -
```
mkdir fastq_data
```
Make sure to put all your forwards and backwards read fastq files within this directory. 
Before generating the genome index you need to install [Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa](https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/) and [Homo_sapiens.GRCh38.111.gtf](https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/) files and put them in the **/reference** directory.
Run the command below to generate a genome index that STAR will use to align your fastq files to the human genome.
```
cd reference/
STAR --runMode genomeGenerate --genomeDir ~/Projects/rnaseqworkflow/reference/ --genomeFastaFiles Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa --sjdbGTFfile Homo_sapiens.GRCh38.111.gtf --runThreadN 8
cd ..
```

For **seq2HLA** make sure you've got all the correct requirements for it. Then run -
```
nextflow run rnaseq2hla.nf
```

For **Optitype** make sure you've got all the correct requirements/dependencies. Then run - 
```
nextflow run rnaoptitype.nf
```

You know nextflow script has ran properly when you get 4 output directorys: QC\_trim (contains the trimmed and QC fastq files), seq2hla\_output / Optitype\_output (contains the hla typing output), STAR\_out (contains the STAR align output) and FeatureCounts\_out (contains FeatureCounts output).

You can choose which output directories to output by editing the nextflow file, within each process you will see a line -
```
publishDir('name_of_output', mode : 'copy')
```
If you do not wish for a specific output directory you can delete this line. (I recommend not deleting this line for the hlatyping process and the FeatureCounts process)

## Roadmap
Potentially at the generating genome index as part of the workflow rather than it already being pregenerated in the **reference** folder (depends how long it takes to clone the repo).

## Implimenting Optitype into future pipelines
Optitype is easier to integrate into python base pipelines, however if you clone the optitype repo you will need to change **result.applymap(get_types)** to **result.map(get_types)** 




