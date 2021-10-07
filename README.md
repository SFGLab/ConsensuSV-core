# ConsensuSV
The tool designed for getting consensus out of multiple SV callers' results.

Important: for the completly automatised fastq-to-vcf (8 SV callers + SNP / Indel calling included) pipeline see: https://github.com/SFGLab/ConsensuSV-pipeline

Requirements:
* bcftools (https://samtools.github.io/bcftools/) in PATH

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

The samples should follow the rule seen in the following figure:

![Sample folder structure](https://github.com/MateuszChilinski/ConsensuSV/blob/master/sample_folder_example.png)

The example command used for the training of the neural network model:
```shell
python main.py -f /home/ConsensuSV/data/ -t
```
The example command used for getting the consensus SVs (the model included in the package is trained on the 11 SV callers shown on the example sample folder structure):
```shell
python main.py -f /home/ConsensuSV/data/ -o consensuSV
```
