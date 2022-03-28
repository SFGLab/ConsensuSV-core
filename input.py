import argparse

def inputHandling():
    """Parsing the arguments.

    Returns:
        argparse.ArgumentParser: Argument parser from argparse.
    """
    parser = argparse.ArgumentParser(description='Gets the SV consensus.')

    parser.add_argument('-f', '--sv_folder', help='Folder containing folders of samples with raw outputs from SV callers (comma-separated). More information on the structure of the samples folder in readme.', required=True)

    parser.add_argument('-mod', '--model', help='Model used for SV discovery (default pretrained.model).', required=False, default="pretrained.model")

    parser.add_argument('-o', '--output', help='Output file prefix.', default="consensuSV_")
    
    parser.add_argument('-of', '--output_folder', help='Output folder. Default in the ConsensuSV folder output/', default="output/")
    
    parser.add_argument('-s', '--samples', help='Samples to include. By default all in the sv_folder. Comma-separated.', required=False, default=None)

    parser.add_argument('-c', '--callers', help='Callers to include. By default all in the folders. Comma-separated.', required=False, default=None)

    parser.add_argument('-m', '--min_overlap', help='File with minimum numbers of SVs in the neighbourhood for the SV to be reported (default min_overlaps).', default="min_overlaps")

    parser.add_argument('-t', '--train', help='Creates new model. Requires truth.vcf to be present in all the sv folders. VCF file truth.vcf is preprocessed even if flag --no_preprocess is set. If you train the model, you need to rerun the program to get the consensus.', action="store_true", required=False)

    parser.add_argument('-np', '--no_preprocess', help='Flag used for skipping the preprocessing process - all the preprocessed files should be in temp/ folder.', action="store_true", required=False)

    args = parser.parse_args()

    return args
