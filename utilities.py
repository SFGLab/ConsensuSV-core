import os


from subprocess import Popen, PIPE, DEVNULL, STDOUT
from shutil import copyfile
from os import listdir
from os.path import isfile, isdir, join
from SVTools import SVTool
debug = 0

def generate_header(sample_name):
    fin = open("header", "rt")
    data = fin.read()
    data = data.replace('SAMPLENAME', sample_name)
    fin.close()
    return data+"\n"

def execute_command(cmd):
    if(debug):
        print(cmd)
        process = Popen(cmd, shell=True, stdout=PIPE)
    else:
        process = Popen(cmd, shell=True, stdout=DEVNULL, stderr=STDOUT)
    process.communicate()

def reheader_all(dirFrom, dirTo, sv_files, sampleName):
    # create temp header with sample name
    copyfile("header", "header_temp")
    fin = open("header_temp", "wt")
    fin.write(generate_header(sampleName))
    fin.close()

    # reheader all files
    for file in sv_files:
        cmd = r"bcftools reheader -h header_temp -o " + dirTo + file + " " + dirFrom + file
        execute_command(cmd)
    os.remove("header_temp")

def preprocessFile(file, sampleName, header):

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

def preprocessFiles(folder, sampleName):

    sv_files = [f for f in listdir(folder) if isfile(join(folder, f))]

    reheader_all(folder, "temp/"+sampleName+"/", sv_files, sampleName)

    sv_files = [f for f in listdir("temp/"+sampleName+"/") if isfile(join("temp/"+sampleName+"/", f))]

    header = generate_header(sampleName)
    for file in sv_files:
        preprocessFile(sampleName+"/"+file, sampleName, header)
    return loadTempFiles(sampleName)

def loadTempFiles(sampleName):
    sv_tools = list()
    
    sv_files = [f for f in listdir("temp/"+sampleName+"/") if isfile(join("temp/"+sampleName+"/", f))]

    for file in sv_files:
        svtool = SVTool("temp/"+sampleName+"/"+file)
        sv_tools.append(svtool)
    sv_tools.sort(key=lambda x: x.tool)
    return sv_tools

def buildFreqDict(candidates):
    freqDict = dict()
    for candidate in candidates:
        key = str(candidate.pos)+"-"+str(candidate.end)
        if key not in freqDict:
            freqDict[key] = 1
        else:
            freqDict[key] += 1
    return freqDict

def findMajority(sv, freqDict, candidates):
    majorityFound = False
    firstCandidate = None
    winKey = ""

    for key in freqDict:
        if(freqDict[key]/len(candidates) >= 0.7):
            #print(sv.chrom + " " + sv.svtype + " " + str(sv.pos) + " - " + str(sv.end))
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
    sv_samples = [d.split('/')[-1] for d in listdir("temp/") if isdir(join("temp/", d))]
    
    sv_tools = set()
    for sampleName in sv_samples:
        sv_files = [f for f in listdir("temp/"+sampleName+"/") if isfile(join("temp/"+sampleName+"/", f))]

        for file in sv_files:
            toolname = file.split(".")[0]
            if(toolname == "truth"):
                continue
            sv_tools.add(toolname)
    return sorted(sv_tools)

def preprocess_Y(Y_vector):
    Y_prepr = list()
    for sv in Y_vector:
        Y_prepr.append(sv.pos)
        Y_prepr.append(sv.end)
    return Y_prepr
def preprocess_X(X_vector):
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
    for candidate in candidates:
        candidate.used = True

def generateAlgorithmsList(candidates):
    algorithms = ""
    for candidate in candidates:
        algorithms += candidate.tool+","
    return algorithms[0:-1]

def generateGenotype(candidates):
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
