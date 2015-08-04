# AFTER
Automatic Filtering, Trimming, and Error Removing for fastq data
Currently it supports Illumina 1.8 or newer format, see:
http://support.illumina.com/help/SequencingAnalysisWorkflow/Content/Vault/Informatics/Sequencing_Analysis/CASAVA/swSEQ_mCA_FASTQFiles.htm

# Version
1.0

# Features
Do following tasks automatically:
1, Filter PolyA/PolyT/PolyC/PolyG reads
2, Trim reads at front and tail
3, Detect and remove bubble artefacts caused by sequencer
4, Filter low-quality reads

# Simple usage:
1, cd to the folder contains all fastq files
2, run: python after.py

# Debubble
If you want to remove bubble artefacts, run
python after --debubble=on

# Full usage:
python after.py [-d input_dir][-1 read1_file] [-2 read1_file] [-7 index1_file] [-5 index2_file] [-g good_output_folder] [-b bad_output_folder] [-f trim_front] [-t trim_tail] [-m min_quality] [-q qualified_quality] [-l max_low_quality] [-p poly_max] [-a allow_poly_mismatch] [-n max_n_count] [--debubble=on/off] [--debubble_dir=xxx] [--draw=on/off] [--read1_flag=_R1_] [--read2_flag=_R2_] [--index1_flag=_I1_] [--index2_flag=_I2_] 

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -1 READ1_FILE, --read1_file=READ1_FILE
                        file name of read1, required. If input_dir is
                        specified, then this arg is ignored.
  -2 READ2_FILE, --read2_file=READ2_FILE
                        file name of read2, if paired. If input_dir is
                        specified, then this arg is ignored.
  -7 INDEX1_FILE, --index1_file=INDEX1_FILE
                        file name of 7' index. If input_dir is specified, then
                        this arg is ignored.
  -5 INDEX2_FILE, --index2_file=INDEX2_FILE
                        file name of 5' index. If input_dir is specified, then
                        this arg is ignored.
  -g GOOD_OUTPUT_FOLDER, --good_output_folder=GOOD_OUTPUT_FOLDER
                        the folder to store good reads, by default it is the
                        same folder contains read1
  -b BAD_OUTPUT_FOLDER, --bad_output_folder=BAD_OUTPUT_FOLDER
                        the folder to store bad reads, by default it is same
                        as good_output_folder
  -f TRIM_FRONT, --trim_front=TRIM_FRONT
                        number of bases to be trimmed in the head of read. -1
                        means auto detect
  -t TRIM_TAIL, --trim_tail=TRIM_TAIL
                        number of bases to be trimmed in the tail of read. -1
                        means auto detect
  -m MIN_QUALITY, --min_quality=MIN_QUALITY
                        if exists one base has quality < min_quality, then
                        this read/pair will be bad. Default 0 means do not
                        filter reads by the least quality
  -q QUALIFIED_QUALITY, --qualified_quality=QUALIFIED_QUALITY
                        the quality value that a base is qualifyed. Default 20
                        means base quality >=Q20 is qualified.
  -l MAX_LOW_QUALITY, --max_low_quality=MAX_LOW_QUALITY
                        if exists more than maxlq bases that quality is lower
                        than qualified quality, then this read/pair is bad.
                        Default 0 means do not filter reads by low quality
                        base count
  -p POLY_MAX, --poly_max=POLY_MAX
                        if exists one polyX(polyG means GGGGGGGGG...), and its
                        length is >= poly_max, then this read/pair is bad.
                        Default is 40
  -a ALLOW_POLY_MISMATCH, --allow_poly_mismatch=ALLOW_POLY_MISMATCH
                        the count of allowed mismatches when evaluating
                        poly_X. Default 5 means disallow any mismatches
  -n MAX_N_COUNT, --max_n_count=MAX_N_COUNT
                        if exists more than maxn bases have N, then this
                        read/pair is bad. Default is 5
  -s MIN_SEQ_LEN, --min_seq_len=MIN_SEQ_LEN
                        if the trimmed read is shorter than min_seq_len, then
                        this read/pair is bad. Default is 35
  -d INPUT_DIR, --input_dir=INPUT_DIR
                        the input dir to process automatically. If read1_file
                        are input_dir are not specified, then current dir (.)
                        is specified to input_dir
  --debubble=DEBUBBLE   specify whether apply debubble algorithm to remove the
                        reads in the bubbles. Default is off
  --debubble_dir=DEBUBBLE_DIR
                        specify the folder to store output of debubble
                        algorithm, default is debubble
  --draw=DRAW           specify whether draw the pictures or not, when use
                        debubble or QC. Default is on
  --read1_flag=READ1_FLAG
                        specify the name flag of read1, default is _R1_, which
                        means a file with name *_R1_* is read1 file
  --read2_flag=READ2_FLAG
                        specify the name flag of read2, default is _R2_, which
                        means a file with name *_R2_* is read2 file
  --index1_flag=INDEX1_FLAG
                        specify the name flag of index1, default is _I1_,
                        which means a file with name *_I1_* is index2 file
  --index2_flag=INDEX2_FLAG
                        specify the name flag of index2, default is _I2_,
                        which means a file with name *_I2_* is index2 file
