import shlex
import os
import numpy
import sys
import pickle
import utilities
import shutil
import time

from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPRegressor
from input import inputHandling
from SVTools import SVariant
from shutil import copyfile
from os.path import isdir, join
from os import listdir
import re

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

args = inputHandling()

min_overlap = args.min_overlap
with open(min_overlap) as fh:
    min_overlaps = dict(re.findall(r'(\S+)\s+(.+)', fh.read()))
for key, value in min_overlaps.items():
    min_overlaps[key] = int(value)
# preprocessing of the files
# problems with no svlen?

X_vector = list()
Y_vector = list()

samples_folder = args.sv_folder
output_folder = args.output_folder

folders = [f for f in listdir(samples_folder) if isdir(join(samples_folder, f))]

if(args.samples):
    samples = args.samples.split(",")
else:
    samples = [f.split('/')[-1] for f in folders]

if(args.callers):
    callers = args.callers.split(",")
    callers.append("truth")
else:
    callers = None

if not (args.no_preprocess):
    if os.path.exists("temp") and os.path.isdir("temp"):
        shutil.rmtree("temp")
    os.mkdir("temp")
if not(args.train):
    if not (os.path.exists("output") and os.path.isdir("output")):
        os.mkdir("output")
for sample in samples:
    sample_dir = samples_folder+sample+"/"
    sample_temp_dir = "temp/"+sample+"/"

    if not (args.no_preprocess):
        if os.path.exists(sample_temp_dir) and os.path.isdir(sample_temp_dir):
            pass
        else:
            os.mkdir(sample_temp_dir)

        print("Preprocessing files of "+sample+"...", end='', flush=True)

        sv_tools = utilities.preprocessFiles(sample_dir, sample, callers)

        print(" DONE!", flush=True)
    else:
        if (args.train):
            copyfile(sample_dir+"truth.vcf", sample_temp_dir+"/truth.vcf")
            utilities.preprocessFile(sample+"/truth.vcf", sample, utilities.generate_header(sample))
        sv_tools = utilities.loadTempFiles(sample, callers)
    percDiff = 0.1

    consensusId = 1

    if not(args.train): # load model
        filename = args.model
        loaded_model = pickle.load(open(filename, 'rb'))

    resulting_svs = list()
    for svtool in sv_tools:
        if (args.train):
            if(svtool.tool != "truth"):
                continue
            print("\tProcessing " + sample + "...", end='', flush=True)
        else:
            if(svtool.tool == "truth"):
                continue
            print("\tProcessing tool " + svtool.tool + "...", end='', flush=True)

        for sv in svtool.sv_list:
            if(sv.used): continue
            candidates = list()
            candidates.append(sv)
            for svtool2 in sv_tools:
                if(svtool2.tool == "truth"):
                    continue
                if(svtool.tool == svtool2.tool):
                    continue
                for sv2 in svtool2.sv_list:
                    if(sv.chrom != sv2.chrom):
                        continue
                    if(sv.svtype != sv2.svtype):
                        continue
                    if(sv2.pos > sv.pos+500): # fix later! it should be dependend on ci or % of svlen
                        break
                    if(sv2.used): continue
                    if(sv.checkOverlap(sv2)):
                       candidates.append(sv2)
                       break
            if(len(candidates) < min_overlaps[sv.svtype]): # if fewer than 3 then no point in checking it out
                continue
            if (args.train): # learning phase
                candidates.remove(sv)
                X_vector.append(candidates)
                Y_vector.append(sv)
            else:
                freqDict = utilities.buildFreqDict(candidates)

                # maybe remove all candidates from svtool once consensus was established based on it?
                consensusGT = utilities.generateGenotype(candidates)
                if("./." in consensusGT): # it's not really well-genotyped anyway, so not worth adding
                    continue
                (majorityFound, firstMajor) = utilities.findMajority(sv, freqDict, candidates)
                if(majorityFound):
                    newSv = SVariant("consensus", None, firstMajor.chrom, firstMajor.pos, "consensus_"+str(consensusId), firstMajor.ref[0], firstMajor.end, consensusGT, firstMajor.svlen, firstMajor.svtype, 
                                     -10, 10, -10, 10, utilities.generateAlgorithmsList(candidates))
                    consensusId += 1
                else:
                    result = loaded_model.predict(utilities.preprocess_X([candidates]))
                    pos = result[0]
                    end = result[1]
                    newSv = SVariant("consensus", None, sv.chrom, int(round(pos)), "consensus_"+str(consensusId), sv.ref[0], int(round(end)), consensusGT, int(round(pos-end)), sv.svtype, -10, 10, -10, 10, utilities.generateAlgorithmsList(candidates))
                    consensusId += 1
                resulting_svs.append(newSv)
                utilities.markUsedCandidates(candidates)
        print(" DONE!", flush=True)
    if not(args.train):
        header = utilities.generate_header(sample)
        with open(sample_temp_dir+"output.vcf", 'w') as fout:
            fout.write(header)
            for sv in resulting_svs:
                fout.write(sv.printVcfLine())

        cmd = "cat "+sample_temp_dir+"output.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "+ "\"sort -k1,1V -k2,2n\"" + r"}' > "+sample_temp_dir+"output_sorted.vcf"
        utilities.execute_command(cmd)
        
        shutil.move(sample_temp_dir+"output_sorted.vcf", output_folder+args.output+"_"+sample+".vcf")
    else:
        os.remove("temp/"+sample+"/truth.vcf")

if (args.train): # learning phase
    print("Preparing sets...", end='', flush=True)

    X_preprocessed_vector = utilities.preprocess_X(X_vector)
    Y_preprocessed_vector = utilities.preprocess_Y(Y_vector)

    print(" DONE!", flush=True)

    X_train, X_test, y_train, y_test = train_test_split(X_preprocessed_vector, Y_preprocessed_vector, test_size=0.1, random_state=42, shuffle=True)
    nn = MLPRegressor(hidden_layer_sizes=((len(sv_tools)-1)*2, len(sv_tools)-1), solver='lbfgs', random_state=0)

    print("Creating the model...", end='', flush=True)

    nn.fit(X_train, y_train)


    nn_score = nn.score(X_test, y_test)
    y_pred = nn.predict(X_test)

    print(" DONE!", flush=True)

    print("Score of model: " + str(nn_score))
    error = abs(y_test-y_pred)
    print("Average abs error (testing set of 10%): " + str(numpy.average(error)) + "std: " + str(numpy.std(error)))

    filename = args.model
    pickle.dump(nn, open(filename, 'wb'))
    
    numpy.savetxt("foo.csv", numpy.concatenate((X_test, numpy.vstack((y_test,y_pred)).T), axis=1), delimiter=',', comments="")

# all files are preprocessed now in unified form
