# AfterQC
Automatic Filtering, Trimming, Error Removing and Quality Control for fastq data   
`AfterQC` can simply go through all fastq files in a folder and then output a <b>good</b> folder and a <b>bad</b> folder, which contains good reads and bad reads of each fastq file   
Currently it supports processing data from HiSeq 2000/2500/3000/4000, Nextseq 500/550, MiniSeq...and other [Illumina 1.8 or newer formats](http://support.illumina.com/help/SequencingAnalysisWorkflow/Content/Vault/Informatics/Sequencing_Analysis/CASAVA/swSEQ_mCA_FASTQFiles.htm)   

# Latest release
0.2.0 (Released on 2016-03-28)

# Features:
`AfterQC` does following tasks automatically:  
* Filters reads with too low quality, too short length or too many N
* Filters reads with abnormal PolyA/PolyT/PolyC/PolyG sequences
* Does per-base quality control and plots the figures
* Trims reads at front and tail, according to QC results
* For pair-end sequencing data, `AfterQC` automatically corrects low quality wrong bases in overlapped area of read1/read2
* Detects and eliminates bubble artifact caused by sequencer due to fluid dynamics issues
* Single molecule barcode sequencing support: if all reads have a single molecule barcode (see duplex sequencing), `AfterQC` detects and splits the barcodes from the reads, and put them in the fastq query names


# Dependency:
`AfterQC` uses `editdistance` module, run following before using `AfterQC`:
```shell
pip install editdistance
```
***WARNING: If you haven't installed `editdistance` module, `AfterQC` will use a python implementation of editdistance, but it will be extremely slow.***

# Simple usage:
* Prepare your fastq files in a folder
* For single-end sequencing, the filenames in the folder should be `*R1*`
* For pair-end sequencing, the filenames in the folder should be `*R1*` and `*R2*`
```shell
cd /path/to/fastq/folder
python path/to/AfterQC/after.py
```
* two folders will be automatically generated, a folder 'good' stores the good reads and a folder 'bad' stores the bad reads
* `AfterQC` will print some statistical information after it is done, such how many good reads, how many bad reads, and how many reads are corrected.

# Full options:
```shell
python after.py [-d input_dir][-1 read1_file] [-2 read1_file] [-7 index1_file] [-5 index2_file] [-g good_output_folder] [-b bad_output_folder] [-f trim_front] [-t trim_tail] [-q qualified_quality_phred] [-l unqualified_base_limit] [-p poly_size_limit] [-a allow_mismatch_in_poly] [-n n_base_limit] [--debubble=on/off] [--debubble_dir=xxx] [--draw=on/off] [--read1_flag=_R1_] [--read2_flag=_R2_] [--index1_flag=_I1_] [--index2_flag=_I2_]
```
***Common options***
```shell
  --version             show program's version number and exit
  -h, --help            show this help message and exit
```
***File (name) options***
```

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
  -d INPUT_DIR, --input_dir=INPUT_DIR
                        the input dir to process automatically. If read1_file
                        are input_dir are not specified, then current dir (.)
                        is specified to input_dir
  -g GOOD_OUTPUT_FOLDER, --good_output_folder=GOOD_OUTPUT_FOLDER
                        the folder to store good reads, by default it is the
                        same folder contains read1
  -b BAD_OUTPUT_FOLDER, --bad_output_folder=BAD_OUTPUT_FOLDER
                        the folder to store bad reads, by default it is same
                        as good_output_folder
  --read1_flag=READ1_FLAG
                        specify the name flag of read1, default is R1, which
                        means a file with name *R1* is read1 file
  --read2_flag=READ2_FLAG
                        specify the name flag of read2, default is R2, which
                        means a file with name *R2* is read2 file
  --index1_flag=INDEX1_FLAG
                        specify the name flag of index1, default is I1,
                        which means a file with name *I1* is index2 file
  --index2_flag=INDEX2_FLAG
                        specify the name flag of index2, default is I2,
                        which means a file with name *I2* is index2 file
```
***Filter options***
```
  -f TRIM_FRONT, --trim_front=TRIM_FRONT
                        number of bases to be trimmed in the head of read. -1
                        means auto detect
  -t TRIM_TAIL, --trim_tail=TRIM_TAIL
                        number of bases to be trimmed in the tail of read. -1
                        means auto detect
  --trim_pair_same=TRIM_PAIR_SAME
                        use same trimming configuration for read1 and read2 to
                        keep their sequence length identical, default is true
                        lots of dedup algorithms require this feature
  -q QUALIFIED_QUALITY_PHRED, --qualified_quality_phred=QUALIFIED_QUALITY_PHRED
                        the quality value that a base is qualifyed. Default 20
                        means base quality >=Q20 is qualified.
  -u UNQUALIFIED_BASE_LIMIT, --unqualified_base_limit=UNQUALIFIED_BASE_LIMIT
                        if exists more than unqualified_base_limit bases that
                        quality is lower than qualified quality, then this
                        read/pair is bad. Default 0 means do not filter reads
                        by low quality base count
  -p POLY_SIZE_LIMIT, --poly_size_limit=POLY_SIZE_LIMIT
                        if exists one polyX(polyG means GGGGGGGGG...), and its
                        length is >= poly_size_limit, then this read/pair is
                        bad. Default is 35
  -a ALLOW_MISMATCH_IN_POLY, --allow_mismatch_in_poly=ALLOW_MISMATCH_IN_POLY
                        the count of allowed mismatches when evaluating
                        poly_X. Default 5 means disallow any mismatches
  -n N_BASE_LIMIT, --n_base_limit=N_BASE_LIMIT
                        if exists more than maxn bases have N, then this
                        read/pair is bad. Default is 5
  -s SEQ_LEN_REQ, --seq_len_req=SEQ_LEN_REQ
                        if the trimmed read is shorter than seq_len_req, then
                        this read/pair is bad. Default is 35
```
***Debubble options***   
If you want to eliminate bubble artifact, turn debubble option on (this is slow, usually you don't need to do this): 
```
  --debubble=DEBUBBLE   specify whether apply debubble algorithm to remove the
                        reads in the bubbles. Default is off
  --debubble_dir=DEBUBBLE_DIR
                        specify the folder to store output of debubble
                        algorithm, default is debubble
  --draw=DRAW           specify whether draw the pictures or not, when use
                        debubble or QC. Default is on
```
***Barcoded sequencing options***
```
  --barcode=BARCODE     specify whether deal with barcode sequencing files, default is on
  --barcode_length=BARCODE_LENGTH
                        specify the designed length of barcode
  --barcode_flag=BARCODE_FLAG
                        specify the name flag of a barcoded file, default is
                        barcode, which means a file with name *barcode* is a
                        barcoded file
  --barcode=BARCODE     specify whether deal with barcode sequencing files,
                        default is on, which means all files with barcode_flag
                        in filename will be treated as barcode sequencing
                        files
```
