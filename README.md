# ConsensuSV-pipeline

Table of Contents
=================

* [What is ConsensuSV?](#what-is-consensusv)
* [Citation](#citation)
* [Requirements](#requirements)
* [Parameters](#parameters)
* [Structure of the data folder](#structure-of-the-data-folder)
* [Comparison to gold-standard set](#comparison-to-gold-standard-set)
* [Examples](#examples)

## What is ConsensuSV?

The tool designed for getting consensus out of multiple SV callers' results.

Important: for the completly automatised fastq-to-vcf (8 SV callers + SNP / Indel calling included) pipeline see: https://github.com/SFGLab/ConsensuSV-pipeline

## Citation

If you use ConsensuSV in your research, we kindly ask you to cite the following publication:

```
@article{Chilinski_ConsensuSVfrom_the_whole-genome_2022,
author = {Chiliński, Mateusz and Plewczynski, Dariusz},
doi = {10.1093/bioinformatics/btac709},
journal = {Bioinformatics},
title = {{ConsensuSV—from the whole-genome sequencing data to the complete variant list}},
year = {2022}
}
```

## Requirements

Requirements:
* bcftools (https://samtools.github.io/bcftools/) in PATH

## Parameters

Options:

Short option | Long option | Description
-------------- | --------------- | ---------------
-f | --sv_folder | older containing folders of samples with raw outputs from SV callers (comma-separated). More information on the structure of the samples folder is shown below.
-mod | --model | Model used for SV discovery (default: pretrained.model).
-o | --output | Output file prefix (default: consensuSV_).
-m | --min_overlap | File with minimum numbers of SVs in the neighbourhood for the SV to be reported (default min_overlaps).
-of | --output_folder | Output folder (default: "output/").
-s | --samples | Samples to include. By default all in the sv_folder. Comma-separated
-c | --callers | Callers to include. By default all in the folders. Comma-separated.
-t | --train | Creates new model. Requires truth.vcf to be present in all the sv folders. VCF file truth.vcf is preprocessed even if flag --no_preprocess is set. If the model is trained, it is required to rerun the program to get the consensus.
-np | --no_preprocess | Flag used for skipping the preprocessing process - all the preprocessed files should be in temp/ folder.

## Structure of the data folder

The samples should follow the rule seen in the following figure:
<p align="center">
<img src="https://github.com/SFGLab/ConsensuSV-core/blob/main/sample_folder_example.png" width="700" height="700" />
</p>

## Implementation details

The workflow of the algorithm is presented in the following figure:

<p align="center">
<img src="https://github.com/SFGLab/ConsensuSV-core/blob/main/workflow.png" width="500"/>
</p>

## Comparison to gold-standard set


## Examples

The example command used for the training of the neural network model:

```shell
python main.py -f /home/ConsensuSV/data/ -t
```
The example command used for getting the consensus SVs (the model included in the package is trained on the 11 SV callers shown on the example sample folder structure):
```shell
python main.py -f /home/ConsensuSV/data/ -o consensuSV
```
