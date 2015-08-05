# AFTER
Automatic Filtering, Trimming, and Error Removing for fastq data＜br /＞
Currently it supports Illumina 1.8 or newer format, see:＜br /＞
http://support.illumina.com/help/SequencingAnalysisWorkflow/Content/Vault/Informatics/Sequencing_Analysis/CASAVA/swSEQ_mCA_FASTQFiles.htm＜br /＞
AFTER can simply go through all fastq files in a folder and then output a <b>good</b> folder and a <b>bad</b> folder, which contains good reads and bad reads of each fastq file＜br /＞

# Version
1.0

# Feedback/contact
infoteam@haplox.com＜br /＞
chen@haplox.com

# Features
AFTER does following tasks automatically:＜br /＞
1, Filter PolyA/PolyT/PolyC/PolyG reads＜br /＞
2, Trim reads at front and tail＜br /＞
3, Detect and remove bubble artefacts caused by sequencer＜br /＞
4, Filter low-quality reads＜br /＞

# Simple usage:
1, cd to the folder contains all fastq files＜br /＞
2, run: python after.py＜br /＞

# Debubble
If you want to remove bubble artefacts, run＜br /＞
python after --debubble=on＜br /＞

# Full usage:
python after.py [-d input_dir][-1 read1_file] [-2 read1_file] [-7 index1_file] [-5 index2_file] [-g good_output_folder] [-b bad_output_folder] [-f trim_front] [-t trim_tail] [-m min_quality] [-q qualified_quality] [-l max_low_quality] [-p poly_max] [-a allow_poly_mismatch] [-n max_n_count] [--debubble=on/off] [--debubble_dir=xxx] [--draw=on/off] [--read1_flag=_R1_] [--read2_flag=_R2_] [--index1_flag=_I1_] [--index2_flag=_I2_] 

Options:
  --version             show program's version number and exit＜br /＞
  -h, --help            show this help message and exit＜br /＞
  -1 READ1_FILE, --read1_file=READ1_FILE＜br /＞
                        file name of read1, required. If input_dir is
                        specified, then this arg is ignored.＜br /＞
  -2 READ2_FILE, --read2_file=READ2_FILE＜br /＞
                        file name of read2, if paired. If input_dir is
                        specified, then this arg is ignored.＜br /＞
  -7 INDEX1_FILE, --index1_file=INDEX1_FILE＜br /＞
                        file name of 7' index. If input_dir is specified, then
                        this arg is ignored.＜br /＞
  -5 INDEX2_FILE, --index2_file=INDEX2_FILE＜br /＞
                        file name of 5' index. If input_dir is specified, then
                        this arg is ignored.＜br /＞
  -g GOOD_OUTPUT_FOLDER, --good_output_folder=GOOD_OUTPUT_FOLDER＜br /＞
                        the folder to store good reads, by default it is the
                        same folder contains read1＜br /＞
  -b BAD_OUTPUT_FOLDER, --bad_output_folder=BAD_OUTPUT_FOLDER＜br /＞
                        the folder to store bad reads, by default it is same
                        as good_output_folder＜br /＞
  -f TRIM_FRONT, --trim_front=TRIM_FRONT＜br /＞
                        number of bases to be trimmed in the head of read. -1
                        means auto detect＜br /＞
  -t TRIM_TAIL, --trim_tail=TRIM_TAIL＜br /＞
                        number of bases to be trimmed in the tail of read. -1
                        means auto detect＜br /＞
  -m MIN_QUALITY, --min_quality=MIN_QUALITY＜br /＞
                        if exists one base has quality < min_quality, then
                        this read/pair will be bad. Default 0 means do not
                        filter reads by the least quality＜br /＞
  -q QUALIFIED_QUALITY, --qualified_quality=QUALIFIED_QUALITY＜br /＞
                        the quality value that a base is qualifyed. Default 20
                        means base quality >=Q20 is qualified.＜br /＞
  -l MAX_LOW_QUALITY, --max_low_quality=MAX_LOW_QUALITY＜br /＞
                        if exists more than maxlq bases that quality is lower
                        than qualified quality, then this read/pair is bad.
                        Default 0 means do not filter reads by low quality
                        base count＜br /＞
  -p POLY_MAX, --poly_max=POLY_MAX＜br /＞
                        if exists one polyX(polyG means GGGGGGGGG...), and its
                        length is >= poly_max, then this read/pair is bad.
                        Default is 40＜br /＞
  -a ALLOW_POLY_MISMATCH, --allow_poly_mismatch=ALLOW_POLY_MISMATCH＜br /＞
                        the count of allowed mismatches when evaluating
                        poly_X. Default 5 means disallow any mismatches＜br /＞
  -n MAX_N_COUNT, --max_n_count=MAX_N_COUNT＜br /＞
                        if exists more than maxn bases have N, then this
                        read/pair is bad. Default is 5＜br /＞
  -s MIN_SEQ_LEN, --min_seq_len=MIN_SEQ_LEN＜br /＞
                        if the trimmed read is shorter than min_seq_len, then
                        this read/pair is bad. Default is 35＜br /＞
  -d INPUT_DIR, --input_dir=INPUT_DIR＜br /＞
                        the input dir to process automatically. If read1_file
                        are input_dir are not specified, then current dir (.)
                        is specified to input_dir＜br /＞
  --debubble=DEBUBBLE   specify whether apply debubble algorithm to remove the
                        reads in the bubbles. Default is off＜br /＞
  --debubble_dir=DEBUBBLE_DIR＜br /＞
                        specify the folder to store output of debubble
                        algorithm, default is debubble＜br /＞
  --draw=DRAW           specify whether draw the pictures or not, when use
                        debubble or QC. Default is on＜br /＞
  --read1_flag=READ1_FLAG＜br /＞
                        specify the name flag of read1, default is _R1_, which
                        means a file with name *_R1_* is read1 file＜br /＞
  --read2_flag=READ2_FLAG＜br /＞
                        specify the name flag of read2, default is _R2_, which
                        means a file with name *_R2_* is read2 file＜br /＞
  --index1_flag=INDEX1_FLAG＜br /＞
                        specify the name flag of index1, default is _I1_,
                        which means a file with name *_I1_* is index2 file＜br /＞
  --index2_flag=INDEX2_FLAG＜br /＞
                        specify the name flag of index2, default is _I2_,
                        which means a file with name *_I2_* is index2 file＜br /＞
