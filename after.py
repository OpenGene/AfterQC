 #!/usr/bin/env python
 
import os,sys
from optparse import OptionParser
import time
import fastq
import preprocesser
from multiprocessing import Process, Queue
import copy
from util import *

def parseCommand():
    usage = "Automatic Filtering, Trimming, and Error Removing for Illumina fastq data(Illumina 1.8 or newer format, see http://support.illumina.com/help/SequencingAnalysisWorkflow/Content/Vault/Informatics/Sequencing_Analysis/CASAVA/swSEQ_mCA_FASTQFiles.htm)\n\nFull command:\n%prog [-d input_dir][-1 read1_file] [-2 read1_file] [-7 index1_file] [-5 index2_file] [-g good_output_folder] [-b bad_output_folder] [-f trim_front] [-t trim_tail] [-q qualified_quality_phred] [-l unqualified_base_limit] [-p poly_size_limit] [-a allow_mismatch_in_poly] [-n n_base_limit] [--debubble=on/off] [--debubble_dir=xxx] [--draw=on/off] [--read1_flag=_R1_] [--read2_flag=_R2_] [--index1_flag=_I1_] [--index2_flag=_I2_] \n\nSimplest usage:\ncd to the folder containing your fastq data, run <python after.py>"
    version = "%prog 0.1.0"
    parser = OptionParser(usage = usage, version = version) 
    parser.add_option("-1", "--read1_file", dest = "read1_file",
        help = "file name of read1, required. If input_dir is specified, then this arg is ignored.")
    parser.add_option("-2", "--read2_file", dest = "read2_file", default = None,
        help = "file name of read2, if paired. If input_dir is specified, then this arg is ignored.")
    parser.add_option("-7", "--index1_file", dest = "index1_file", default = None,
        help = "file name of 7' index. If input_dir is specified, then this arg is ignored.")
    parser.add_option("-5", "--index2_file", dest = "index2_file", default = None,
        help = "file name of 5' index. If input_dir is specified, then this arg is ignored.")
    parser.add_option("-d", "--input_dir", dest = "input_dir", default = None,
        help = "the input dir to process automatically. If read1_file are input_dir are not specified, then current dir (.) is specified to input_dir")
    parser.add_option("-g", "--good_output_folder", dest = "good_output_folder", default = "good",
        help = "the folder to store good reads, by default it is the same folder contains read1")
    parser.add_option("-b", "--bad_output_folder", dest = "bad_output_folder", default = "bad",
        help = "the folder to store bad reads, by default it is same as good_output_folder")
    parser.add_option("", "--read1_flag", dest = "read1_flag", default = "_R1_",
        help = "specify the name flag of read1, default is _R1_, which means a file with name *_R1_* is read1 file")
    parser.add_option("", "--read2_flag", dest = "read2_flag", default = "_R2_",
        help = "specify the name flag of read2, default is _R2_, which means a file with name *_R2_* is read2 file")
    parser.add_option("", "--index1_flag", dest = "index1_flag", default = "_I1_",
        help = "specify the name flag of index1, default is _I1_, which means a file with name *_I1_* is index2 file")
    parser.add_option("", "--index2_flag", dest = "index2_flag", default = "_I2_",
        help = "specify the name flag of index2, default is _I2_, which means a file with name *_I2_* is index2 file")
    parser.add_option("-f", "--trim_front", dest = "trim_front", default = -1, type = "int",
        help = "number of bases to be trimmed in the head of read. -1 means auto detect")
    parser.add_option("-t", "--trim_tail", dest = "trim_tail", default = -1, type = "int",
        help = "number of bases to be trimmed in the tail of read. -1 means auto detect")
    parser.add_option("-q", "--qualified_quality_phred", dest = "qualified_quality_phred", default = 20, type = "int",
        help = "the quality value that a base is qualifyed. Default 20 means phred base quality >=Q20 is qualified.")
    parser.add_option("-u", "--unqualified_base_limit", dest = "unqualified_base_limit", default = 0, type = "int",
        help = "if exists more than unqualified_base_limit bases that quality is lower than qualified quality, then this read/pair is bad. Default 0 means do not filter reads by low quality base count")
    parser.add_option("-p", "--poly_size_limit", dest = "poly_size_limit", default = 40, type = "int",
        help = "if exists one polyX(polyG means GGGGGGGGG...), and its length is >= poly_size_limit, then this read/pair is bad. Default is 40")
    parser.add_option("-a", "--allow_mismatch_in_poly", dest = "allow_mismatch_in_poly", default = 5, type = "int",
        help = "the count of allowed mismatches when evaluating poly_X. Default 5 means disallow any mismatches")
    parser.add_option("-n", "--n_base_limit", dest = "n_base_limit", default = 5, type = "int",
        help = "if exists more than maxn bases have N, then this read/pair is bad. Default is 5")
    parser.add_option("-s", "--seq_len_req", dest = "seq_len_req", default = 35, type = "int",
        help = "if the trimmed read is shorter than seq_len_req, then this read/pair is bad. Default is 35")
    parser.add_option("", "--debubble", dest = "debubble", default = "off",
        help = "specify whether apply debubble algorithm to remove the reads in the bubbles. Default is off")
    parser.add_option("", "--debubble_dir", dest = "debubble_dir", default = "debubble",
        help = "specify the folder to store output of debubble algorithm, default is debubble")
    parser.add_option("", "--draw", dest = "draw", default = "on",
        help = "specify whether draw the pictures or not, when use debubble or QC. Default is on")
    parser.add_option("", "--barcode", dest = "barcode", default = "on",
        help = "specify whether deal with barcode sequencing files, default is on, which means all files with barcode_flag in filename will be treated as barcode sequencing files")
    parser.add_option("", "--barcode_length", dest = "barcode_length", default = 12, type="int",
        help = "specify the designed length of barcode")
    parser.add_option("", "--barcode_flag", dest = "barcode_flag", default = "barcode",
        help = "specify the name flag of a barcoded file, default is barcode, which means a file with name *barcode* is a barcoded file")
    parser.add_option("", "--barcode_verify", dest = "barcode_verify", default = "CAGTA",
        help = "specify the verify sequence of a barcode which is adjunct to the barcode")
    return parser.parse_args()

def processDir(folder, options):
    
    fqext = (".fq", ".fastq", "fq.gz", ".fastq.gz")
    read1name = options.read1_flag
    read2name = options.read2_flag
    index1name = options.index1_flag
    index2name = options.index2_flag
    barcodeflag = options.barcode_flag
    
    #is not a dir
    if not os.path.isdir(folder):
        return
        
    options_list = []
    
    files = os.listdir(folder)
    for f in files:
        path = os.path.join(folder, f)
        if os.path.isdir(path):
            continue
        
        isfq = False
        for ext in fqext:
            if f.endswith(ext):isfq = True
        if isfq == False:
            continue

        # here we skip those files with name starting with Undetermined
        # because these files are usually with unknown barcode and have no need to be processed
        if f.startswith("Undetermined"):
            continue
        
        #find read1 file
        if read1name in f:
            print(f)
            opt = copy.copy(options)
            read1 = path
            opt.read1_file = read1
            read2 = read1.replace(read1name, read2name)
            if os.path.exists(read2):opt.read2_file = read2
            index1 = read1.replace(read1name, index1name)
            if os.path.exists(index1):opt.index1_file = index1
            index2 = read1.replace(read1name, index2name)
            if os.path.exists(index2):opt.index2_file = index2
            if barcodeflag in f and parseBool(options.barcode):
                opt.barcode = True
                #for barcode sequencing we don't trim it at front
                opt.trim_front = 0
                opt.trim_front2 = 0
            else:
                opt.barcode = False
            options_list.append(opt)
    
    if len(options_list) == 0:
        print("no read files to run with, do you call the program correctly?")
        print("see -h for help")
        return
        
    #create a job for each options
    jobs = [Process(target = processOptions, args = (o, )) for o in options_list]
    for job in jobs: job.start()
    #wait for them to complete
    for job in jobs: job.join()
    
def processOptions(options):
    filter = preprocesser.seqFilter(options)
    filter.run()
    
def runDebubble(options):
    #lazy import debubble here because debubble uses PIL, which is not supported by pypy
    #we can run with pypy with debubble off
    import debubble
    if os.path.exists(os.path.join(options.debubble_dir, "circles.csv")):
        return
    print("runDebubble")
    debubble.debubbleDir(options.input_dir, 20, options.debubble_dir, options.draw)
    
def main():
    time1 = time.time()
    if sys.version_info.major >2:
        print('python3 is not supported yet, please use python2')
        sys.exit(1)
    
    (options, args) = parseCommand()
    options.debubble = parseBool(options.debubble)
    options.draw = parseBool(options.draw)
    options.trim_front2 = options.trim_front
    options.trim_tail2 = options.trim_tail
    
    if options.input_dir == None and options.read1_file == None:
        print('specify current dir as input dir')
        options.input_dir="."
    
    if options.input_dir != None:
        if options.debubble:
            runDebubble(options)
        processDir(options.input_dir, options)
    else:    
        processOptions(options)
    
    time2 = time.time()
    print('Time used: ' + str(time2-time1))
    
if __name__  == "__main__":
    main()
