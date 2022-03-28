import os


from subprocess import Popen, PIPE, DEVNULL, STDOUT
from shutil import copyfile
from os import listdir
from os.path import isfile, isdir, join
from SVTools import SVTool
debug = 0
"""Debug level."""

def generate_header(sample_name):
    """Function for generating the header for output VCF file.

    Args:
        sample_name (str): Name of the sample.

    Returns:
        str: Header for the VCF file.
    """
    fin = open("header", "rt")
    data = fin.read()
    data = data.replace('SAMPLENAME', sample_name)
    fin.close()
    return data+"\n"

def execute_command(cmd):
    """Wrapper for executing CLI commands.

    Args:
        cmd (str): Command to execute.
    """
    if(debug):
        print(cmd)
        process = Popen(cmd, shell=True, stdout=PIPE)
    else:
        process = Popen(cmd, shell=True, stdout=DEVNULL, stderr=STDOUT)
    process.communicate()

def reheader_all(dirFrom, dirTo, sv_files, sampleName):
    """Function for changing the header to the desired one for all the input VCF files."""
    # create temp header with sample name
    copyfile("header", "header_temp"+sampleName)
    fin = open("header_temp"+sampleName, "wt")
    fin.write(generate_header(sampleName))
    fin.close()

    # reheader all files
    for file in sv_files:
        cmd = r"bcftools reheader -h header_temp"+sampleName+" -o " + dirTo + file + " " + dirFrom + file
        execute_command(cmd)
    os.remove("header_temp"+sampleName)

def preprocessFile(file, sampleName, header):
    """Function for preprocessing one file.

    Args:
        file (str): Full path to the file.
        sampleName (str): Name of the sample.
        header (str): Header to put into the file."""

    cmd = "sed -i '/:ME:/d' temp/" + file
    execute_command(cmd)

    cmd = "sed -i '/0\/0/d' temp/" + file
    execute_command(cmd)

    cmd = "awk -F " + r"'\t'" + " '{ $4 = ($4 == \"\.\" ? \"N\" : $4) } 1' OFS=" + r"'\t' temp/" + file + " > temp/" + file + "_2"
    execute_command(cmd)
    
    cmd = "cat temp/" + file + r"_2 | awk '$1 ~ /^#/ {print $0;next} {print $0 | "+ "\"sort -k1,1V -k2,2n\"" + r"}' > temp/" + file
    execute_command(cmd)
    
    # remove MEI if there are any

    # ensures there are no . in ref
    additional_filters = r"SVLEN=%SVLEN;SVTYPE=%SVTYPE;CIPOS=%CIPOS;CIEND=%CIEND"

    cmd = r"bcftools query -H -t chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr22,chrX,chrY,chrM -i '(QUAL >= 30 || QUAL = " + "\".\"" + r") && ((SVLEN = " + "\".\"" + r") || (SVLEN < 50000 && SVLEN > 50) || (SVLEN > -50000 && SVLEN < -50))' -f '%CHROM\t%POS\t%ID\t%REF\t%FIRST_ALT\t%QUAL\t%FILTER\tEND=%END;"+additional_filters+r"\tGT\t[%GT]\n' -o temp/"+file+"_2 temp/"+file    
    execute_command(cmd)

    #os.replace("temp/"+sampleName+"/"+file+"_2", "temp/"+sampleName+"/"+file)
    os.replace("temp/"+file+"_2", "temp/"+file)
    
    with open("temp/"+file, 'r') as fin:
        data = fin.read().splitlines(True)
    with open("temp/"+file, 'w') as fout:
        fout.write(header)
        fout.writelines(data[1:])

def preprocessFiles(folder, sampleName, callers):
    """Function for preprocessing all the files in a folder.

    Args:
        folder (str): _description_
        sampleName (str): Name of the sample.
        callers (list of str): List of used SV callers.

    Returns:
        list of SVTool: List with SV tools.
    """
    sv_files = [f for f in listdir(folder) if (isfile(join(folder, f)) and ".vcf" in f and (callers is None or f.split(".vcf")[0] in callers))]

    reheader_all(folder, "temp/"+sampleName+"/", sv_files, sampleName)

    sv_files = [f for f in listdir("temp/"+sampleName+"/") if isfile(join("temp/"+sampleName+"/", f))]

    header = generate_header(sampleName)
    for file in sv_files:
        preprocessFile(sampleName+"/"+file, sampleName, header)
    return loadTempFiles(sampleName)

def loadTempFiles(sampleName):
    """Function for loading the temporary files.

    Args:
        sampleName (str): Name of the sample.

    Returns:
        list of SVTool: List with SV tools.
    """
    sv_tools = list()
    
    sv_files = [f for f in listdir("temp/"+sampleName+"/") if isfile(join("temp/"+sampleName+"/", f))]

    for file in sv_files:
        svtool = SVTool("temp/"+sampleName+"/"+file)
        sv_tools.append(svtool)
    sv_tools.sort(key=lambda x: x.tool)
    return sv_tools

def buildFreqDict(candidates):
    """Function for building the frequency dictionary of the SV candidates.

    Args:
        candidates (list of SVariant): List of structural variants.

    Returns:
        dict of {str : int}: Dictionary containing the frequencies of candidates.
    """
    freqDict = dict()
    for candidate in candidates:
        key = str(candidate.pos)+"-"+str(candidate.end)
        if key not in freqDict:
            freqDict[key] = 1
        else:
            freqDict[key] += 1
    return freqDict

def findMajority(sv, freqDict, candidates):
    """Function for finding the majority of the candidates.

    Args:
        sv (SVariant): Structural variant to find majority for.
        freqDict (str -> int): Frequency dictionary of candidates.
        candidates (list of SVariant): List with all the candidates.

    Returns:
        tuple of (bool, SVariant): Tuple containing information whether the majority was found, and first SV candidate of that majority, which will be used for creation of the consensus SV.
    """
    majorityFound = False
    firstCandidate = None
    winKey = ""

    for key in freqDict:
        if(freqDict[key]/len(candidates) >= 0.7):
            majorityFound = True
            winKey = key
            break
    if(majorityFound):
        for candidate in candidates:
            key = str(candidate.pos)+"-"+str(candidate.end)
            if(key == winKey):
                firstCandidate = candidate
                break
    return (majorityFound, firstCandidate)

def createSVTable():
    """Function for creating a list of SVTools used in all samples.

    Returns:
        list of SVTool: List of SVTools.
    """
    sv_samples = [d.split('/')[-1] for d in listdir("temp/") if isdir(join("temp/", d))]
    
    sv_tools = set()
    for sampleName in sv_samples:
        sv_files = [f for f in listdir("temp/"+sampleName+"/") if isfile(join("temp/"+sampleName+"/", f)) and "vcf" in f]

        for file in sv_files:
            toolname = file.split(".")[0]
            if(toolname == "truth" or toolname == "output" or toolname == "output_sorted" or "_2" in file):
                continue
            sv_tools.add(toolname)
    return sorted(sv_tools)

def preprocess_Y(Y_vector):
    """Preprocessing the Y vector for the ML prediction.

    Args:
        Y_vector (list of SVariant): List with structural variants.

    Returns:
        list of int: List with coordinates of the structural variants (pos and end).
    """
    Y_prepr = list()
    for sv in Y_vector:
        Y_prepr.append(sv.pos)
        Y_prepr.append(sv.end)
    return Y_prepr

def preprocess_X(X_vector):
    """Function for preprocessing the X vector for the ML prediction.

    Args:
        X_vector (list of SVariant): List with structural variants.

    Returns:
        list of int: List with coordinates of the candidate structural variants - one from each tool.
    """
    X_prepr = list()
    sv_all_tools = createSVTable()
    for candidates in X_vector:
        candidatesY_pos = list()
        candidatesY_end = list()
        
        for tool in sv_all_tools:
            found = False
            for sv in candidates:
                if(tool == sv.tool):
                    candidatesY_pos.append(sv.pos)
                    candidatesY_end.append(sv.end)
                    found = True
                    break
            if(found):
                continue
            candidatesY_pos.append(sum(c.pos for c in candidates)/len(candidates)) # tool not present
            candidatesY_end.append(sum(c.end for c in candidates)/len(candidates))
        X_prepr.append(candidatesY_pos)
        X_prepr.append(candidatesY_end)
    return X_prepr

def markUsedCandidates(candidates):
    """Function marks all used Candidates.

    Args:
        candidates (list of SVariant): List with all candidate SV variants.
    """
    for candidate in candidates:
        candidate.used = True

def generateAlgorithmsList(candidates):
    """Creates list of algorithms that support given SV.

    Args:
        candidates (list of SVariant): List with all candidate SV variants.

    Returns:
        str: Comma-separated list of algorithms from all the candidates.
    """
    algorithms = ""
    for candidate in candidates:
        algorithms += candidate.tool+","
    return algorithms[0:-1]

def generateGenotype(candidates):
    """Voting algorithm for generating of the genotype of the SV.

    Args:
        candidates (list of SVariant): List with all candidate SV variants.

    Returns:
        str: Genotype of the variant.
    """
    gt0 = 0
    gt1 = 0
    for candidate in candidates:
        if("0/1" in candidate.gt or "./." in candidate.gt or "." in candidate.gt): # if genotype is ./. it's worth adding as 0/1, because the tool detected there something after all
            gt0 += 1
        elif("1/1" in candidate.gt):
            gt1 += 1
    if(gt0 == 0 and gt1 == 0):
        return "./."
    if(gt1 >= gt0):
        return "1/1"
    else:
        return "0/1"
