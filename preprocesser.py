 #!/usr/bin/env python
 
import os,sys
import re
from optparse import OptionParser
import time
import fastq
import util
import barcodeprocesser
import json
from qualitycontrol import *
from qcreporter import QCReporter

def getMainName(filename):
    baseName = os.path.basename(filename)
    mainName = baseName.replace(".fastq", "").replace(".fq", "").replace(".gz", "")
    return mainName

def trim(read, front, tail):
    if tail>0:
        #\n will be trimmed, so add it back
        read[1] = read[1][front:-tail]
        read[3] = read[3][front:-tail]
    else:
        read[1] = read[1][front:]
        read[3] = read[3][front:]
        
    return read
    
def hasPolyX(seq, maxPoly, mismatch):
    if(len(seq)<maxPoly):
        return None
    
    polyCount = {}
    polyArray = ("A", "T", "C", "G", "a", "t", "c", "g", "N")
    for poly in polyArray: polyCount[poly] = 0
    
    for x in xrange(len(seq)):
        frontbase = seq[x]
        
        if not frontbase in polyArray:
            return None
        
        if x >= maxPoly:
            tailbase = seq[x-maxPoly]
            polyCount[tailbase] -= 1
        
        polyCount[frontbase] += 1
        if polyCount[frontbase] >= maxPoly - mismatch:
            return frontbase
    return None
    
def minQuality(read):
    qualStr = read[3]
    minQual = 255
    for q in qualStr:
        if minQual > ord(q):
            minQual = ord(q)
    return minQual - 33
    
def lowQualityNum(read, qual):
    qual += 33
    qualStr = read[3]
    lowQualNum = 0
    for q in qualStr:
        if ord(q) < qual:
            lowQualNum += 1
    return lowQualNum
    
def nNumber(read):
    seqStr = read[1]
    nNum = 0
    for s in seqStr:
        if s == 'N':
            nNum += 1
    return nNum

def getOverlap(r, overlap_len):
    ret = []
    ret.append(r[0])
    ret.append(r[1][len(r[1]) - overlap_len:])
    ret.append(r[2])
    ret.append(r[3][len(r[3]) - overlap_len:])
    return ret

def makeDict(opt):
    d = {
        'index2_flag': opt.index2_flag,
        'draw': opt.draw,
        'barcode':opt.barcode ,
        'index1_flag':opt.index1_flag,
        'seq_len_req': opt.seq_len_req,
        'index1_file': opt.index1_file,
        'overlap_output_folder': opt.overlap_output_folder,
        'trim_tail': opt.trim_tail,
        'trim_pair_same': opt.trim_pair_same,
        'poly_size_limit': opt.poly_size_limit,
        'good_output_folder': opt.good_output_folder,
        'debubble_dir': opt.debubble_dir,
        'index2_file': opt.index2_file,
        'qualified_quality_phred': opt.qualified_quality_phred,
        'barcode_flag': opt.barcode_flag,
        'trim_front': opt.trim_front,
        'barcode_verify': opt.barcode_verify,
        'read2_file': opt.read2_file,
        'n_base_limit': opt.n_base_limit,
        'barcode_length': opt.barcode_length,
        'trim_tail2': opt.trim_tail2,
        'unqualified_base_limit': opt.unqualified_base_limit,
        'allow_mismatch_in_poly': opt.allow_mismatch_in_poly,
        'input_dir': opt.input_dir,
        'read1_file': opt.read1_file,
        'read2_flag': opt.read2_flag,
        'store_overlap': opt.store_overlap,
        'debubble': opt.debubble,
        'read1_flag': opt.read1_flag,
        'trim_front2': opt.trim_front2,
        'bad_output_folder': opt.bad_output_folder,
        'qc_only': opt.qc_only,
        'qc_sample': opt.qc_sample,
        'qc_kmer': opt.qc_kmer
    }
    return d

def init_error_matrix():
    error_matrix = {}
    for correct_base in ALL_BASES:
        error_matrix[correct_base]={}
        for error_base in ALL_BASES:
            if correct_base != error_base:
                error_matrix[correct_base][error_base] = 0
    return error_matrix

def merge_error_matrix(merge_to, merge_from):
    for correct_base in ALL_BASES:
        for error_base in ALL_BASES:
            if correct_base != error_base:
                merge_to[correct_base][error_base] += merge_from[correct_base][error_base]
    
########################### seqFilter
class seqFilter:
    
    #opt is an object contains lots of parameters
    def __init__(self, opt):
        self.options = opt
        self.bubbleCircles = {}
        self.bubbleTiles = []
        
        #detect if the input is paired and if it has index files
        if self.options.read2_file != None:
            self.paired = True
        if self.options.index1_file != None:
            self.hasIndex = True

        self.pattern = re.compile(r'\S+\:\d+\:\S+\:\d+\:\d+\:\d+\:\d+')

    def loadBubbleCircles(self):
        bubbleCircleFile = os.path.join(self.options.debubble_dir, "circles.csv")
        if not os.path.exists(bubbleCircleFile):
            return
        with open(bubbleCircleFile) as f:
            rows = f.readlines()
            for row in rows[1:]:
                r = row.split(",")
                x = float(r[0])
                y = float(r[1])
                radius = float(r[2])
                lane = int(r[3])
                tile = int(r[4])
                circle = (x,y,radius,lane,tile)
                if tile not in self.bubbleTiles:
                    self.bubbleTiles.append(tile)
                    self.bubbleCircles[tile]=[]
                self.bubbleCircles[tile].append(circle)
    
    def isInBubble(self, seqInfo):
        #illumina sequence name line format
        #@<instrument>:<run number>:<flowcell ID>:<lane>:<tile_no>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<index sequence>

        match = self.pattern.search(seqInfo);
        if not match:
            return False

        items = match.group().split(":")
        if len(items) < 7:
            return False
            
        lane = int(items[3])
        tile_no = items[4]
        tile = int(tile_no[1:])
        x = int(items[5])
        y = int(items[6])
        if tile not in self.bubbleTiles:
            return False
        for circle in self.bubbleCircles[tile]:
            cx = circle[0]
            cy = circle[1]
            cr = circle[2]
            clane = circle[3]
            if clane == lane:
                if (cx-x)*(cx-x) + (cy-y)*(cy-y) < cr*cr:
                    return True
                    
        return False

    def writeReads(self, r1, r2, i1, i2, r1_file, r2_file, i1_file, i2_file, flag):
        if self.options.qc_only:
            return

        if r1!=None and r1_file!=None:
            #add flag into the read name
            if flag!=None:
                r1[0] = "@" + flag + r1[0][1:]
            r1_file.writeLines(r1)

        if r2!=None and r2_file!=None:
            #add flag into the read name
            if flag!=None:
                r2[0] = "@" + flag + r2[0][1:]
            r2_file.writeLines(r2)

        if i1!=None and i1_file!=None:
            #add flag into the read name
            if flag!=None:
                i1[0] = "@" + flag + i1[0][1:]
            i1_file.writeLines(i1)

        if i2!=None and i2_file!=None:
            #add flag into the read name
            if flag!=None:
                i2[0] = "@" + flag + i2[0][1:]
            i2_file.writeLines(i2)

    def run(self):
        if self.options.debubble:
            self.loadBubbleCircles()

        #read1_file is required
        read1_file = fastq.Reader(self.options.read1_file)

        #no front trim if sequence is barcoded
        if self.options.barcode:
            self.options.trim_front = 0

        reporter = QCReporter()

        self.r1qc_prefilter = QualityControl(self.options.qc_sample, self.options.qc_kmer)
        self.r2qc_prefilter = QualityControl(self.options.qc_sample, self.options.qc_kmer)
        self.r1qc_prefilter.statFile(self.options.read1_file)
        if self.options.read2_file != None:
            self.r2qc_prefilter.statFile(self.options.read2_file)

        self.r1qc_postfilter = QualityControl(self.options.qc_sample, self.options.qc_kmer)
        self.r2qc_postfilter = QualityControl(self.options.qc_sample, self.options.qc_kmer)

        readLen = self.r1qc_prefilter.readLen
        overlap_histgram = [0 for x in xrange(readLen+1)]
        distance_histgram = [0 for x in xrange(readLen+1)]

        #auto detect trim front and trim tail
        if self.options.trim_front == -1 or self.options.trim_tail == -1:
            #auto trim for read1
            trimFront, trimTail = self.r1qc_prefilter.autoTrim()
            if self.options.trim_front == -1:
                self.options.trim_front = trimFront
            if self.options.trim_tail == -1:
                self.options.trim_tail = trimTail
            #auto trim for read2
            if self.options.read2_file != None:
                # check if we should keep same trimming for read1/read2 to keep their length identical
                # this option is on by default because lots of dedup algorithms require this feature
                if self.options.trim_pair_same:
                    self.options.trim_front2 = self.options.trim_front
                    self.options.trim_tail2 = self.options.trim_tail
                else:
                    trimFront2, trimTail2 = self.r2qc_prefilter.autoTrim()
                    if self.options.trim_front2 == -1:
                        self.options.trim_front2 = trimFront2
                    if self.options.trim_tail2 == -1:
                        self.options.trim_tail2 = trimTail2
                
        print(self.options.read1_file + " options:")
        print(self.options)
        
        #if good output folder not specified, set it as the same folder of read1 file
        good_dir = self.options.good_output_folder
        if good_dir == None:
            good_dir = os.path.dirname(self.options.read1_file)

        #if bad output folder not specified, set it as the same folder of read1 file            
        bad_dir = self.options.bad_output_folder
        if bad_dir == None:
            bad_dir = os.path.join(os.path.dirname(os.path.dirname(good_dir+"/")), "bad")

        #if overlap output folder not specified, set it as the same folder of read1 file
        overlap_dir = self.options.overlap_output_folder
        if overlap_dir == None:
#            overlap_dir = os.path.dirname(self.options.read1_file)
            overlap_dir = os.path.join(os.path.dirname(os.path.dirname(good_dir+"/")), "overlap")

        #save QC results at the same folder of good
        qc_base_folder = self.options.report_output_folder
        if qc_base_folder == None:
            qc_base_folder =  os.path.join(os.path.dirname(os.path.dirname(good_dir+"/")), "QC")
        if not os.path.exists(qc_base_folder):
            os.makedirs(qc_base_folder)
        qc_dir = qc_base_folder
            
        if not os.path.exists(good_dir):
            os.makedirs(good_dir)
            
        if not os.path.exists(bad_dir):
            os.makedirs(bad_dir)

        if self.options.store_overlap and self.options.read2_file != None and (not os.path.exists(overlap_dir)):
            os.makedirs(overlap_dir)
        
        good_read1_file = None
        bad_read1_file = None
        overlap_read1_file = None
        if not self.options.qc_only:
            good_read1_file = fastq.Writer(os.path.join(good_dir, getMainName(self.options.read1_file)+".good.fq"))
            bad_read1_file = fastq.Writer(os.path.join(bad_dir, getMainName(self.options.read1_file)+".bad.fq"))

            overlap_read1_file = None
            if self.options.store_overlap:
                overlap_read1_file = fastq.Writer(os.path.join(overlap_dir, getMainName(self.options.read1_file)+".overlap.fq"))
        
        #other files are optional
        read2_file = None
        good_read2_file = None
        bad_read2_file = None
        overlap_read2_file = None

        index1_file = None
        good_index1_file = None
        bad_index1_file = None
        overlap_index1_file = None

        index2_file = None
        good_index2_file = None
        bad_index2_file = None
        overlap_index2_file = None
        
        #if other files are specified, then read them
        if self.options.read2_file != None:
            read2_file = fastq.Reader(self.options.read2_file)
            if not self.options.qc_only:
                good_read2_file = fastq.Writer(os.path.join(good_dir, getMainName(self.options.read2_file)+".good.fq"))
                bad_read2_file = fastq.Writer(os.path.join(bad_dir, getMainName(self.options.read2_file)+".bad.fq"))
                if self.options.store_overlap and self.options.read2_file != None:
                    overlap_read2_file = fastq.Writer(os.path.join(overlap_dir, getMainName(self.options.read2_file)+".overlap.fq"))
        if self.options.index1_file != None:
            index1_file = fastq.Reader(self.options.index1_file)
            if not self.options.qc_only:
                good_index1_file = fastq.Writer(os.path.join(good_dir, getMainName(self.options.index1_file)+".good.fq"))
                bad_index1_file = fastq.Writer(os.path.join(bad_dir, getMainName(self.options.index1_file)+".bad.fq"))
                if self.options.store_overlap and self.options.read2_file != None:
                    overlap_index1_file = fastq.Writer(os.path.join(overlap_dir, getMainName(self.options.index1_file)+".overlap.fq"))
        if self.options.index2_file != None:
            index2_file = fastq.Reader(self.options.index2_file)
            if not self.options.qc_only:
                good_index2_file = fastq.Writer(os.path.join(good_dir, getMainName(self.options.index2_file)+".good.fq"))
                bad_index2_file = fastq.Writer(os.path.join(bad_dir, getMainName(self.options.index2_file)+".bad.fq"))
                if self.options.store_overlap and self.options.read2_file != None:
                    overlap_index2_file = fastq.Writer(os.path.join(overlap_dir, getMainName(self.options.index2_file)+".overlap.fq"))
            
        r1 = None
        r2 = None
        i1 = None
        i2 = None

        # stat numbers
        TOTAL_BASES = 0
        GOOD_BASES = 0
        TOTAL_READS = 0
        GOOD_READS = 0
        BAD_READS = 0
        BADBCD1 = 0
        BADBCD2 = 0
        BADTRIM1 = 0
        BADTRIM2 = 0
        BADBBL = 0
        BADLEN = 0
        BADPOL = 0
        BADLQC = 0
        BADNCT = 0
        BADINDEL = 0
        BADMISMATCH = 0
        READ_CORRECTED = 0
        BASE_CORRECTED = 0
        BASE_ZERO_QUAL_MASKED = 0
        OVERLAPPED = 0
        OVERLAP_LEN_SUM = 0
        OVERLAP_BASE_SUM = 0
        # error profiling by overlap analysis
        OVERLAP_BASE_ERR = 0
        OVERLAP_ERR_MATRIX = init_error_matrix()

        #adapter trimming by overlap analysis
        TRIMMED_ADAPTER_BASE = 0

        while True:
            r1 = read1_file.nextRead()
            if r1==None:
                break
            else:
                TOTAL_BASES += len(r1[1])
                
            if read2_file != None:
                r2 = read2_file.nextRead()
                if r2==None:
                    break
            if index1_file != None:
                i1 = index1_file.nextRead()
                if i1==None:
                    break
            if index2_file != None:
                i2 = index2_file.nextRead()
                if i2==None:
                    break
                else:
                    TOTAL_BASES += len(r2[1])

            TOTAL_READS += 1
                    
            #barcode processing
            if self.options.barcode:
                barcodeLen1 = barcodeprocesser.detectBarcode(r1[1], self.options.barcode_length, self.options.barcode_verify)
                if barcodeLen1 == 0:
                    self.writeReads(r1, r2, i1, i2, bad_read1_file, bad_read2_file, bad_index1_file, bad_index2_file, "BADBCD1")
                    BADBCD1 += 1
                    continue
                else:
                    if r2 == None:
                        barcodeprocesser.moveBarcodeToName(r1, self.options.barcode_length, self.options.barcode_verify)
                    else:
                        barcodeLen2 = barcodeprocesser.detectBarcode(r2[1], self.options.barcode_length, self.options.barcode_verify)
                        if barcodeLen2 == 0:
                            self.writeReads(r1, r2, i1, i2, bad_read1_file, bad_read2_file, bad_index1_file, bad_index2_file, "BADBCD2")
                            BADBCD2 += 1
                            continue
                        else:
                            barcodeprocesser.moveAndTrimPair(r1, r2, barcodeLen1, barcodeLen2, self.options.barcode_verify)
            
            #trim
            if self.options.trim_front > 0 or self.options.trim_tail > 0:
                r1 = trim(r1, self.options.trim_front, self.options.trim_tail)
                if len(r1[1]) < 5:
                    self.writeReads(r1, r2, i1, i2, bad_read1_file, bad_read2_file, bad_index1_file, bad_index2_file, "BADTRIM1")
                    BADTRIM1 += 1
                    continue
                if r2 != None:
                    r2 = trim(r2, self.options.trim_front2, self.options.trim_tail2)
                    if len(r2[1]) < 5:
                        self.writeReads(r1, r2, i1, i2, bad_read1_file, bad_read2_file, bad_index1_file, bad_index2_file, "BADTRIM2")
                        BADTRIM2 += 1
                        continue

            #filter debubble
            if self.options.debubble:
                if self.isInBubble(r1[0]):
                    self.writeReads(r1, r2, i1, i2, bad_read1_file, bad_read2_file, bad_index1_file, bad_index2_file, "BADBBL")
                    BADBBL += 1
                    continue
            
            #filter sequence length
            if len(r1[1])<self.options.seq_len_req:
                self.writeReads(r1, r2, i1, i2, bad_read1_file, bad_read2_file, bad_index1_file, bad_index2_file, "BADLEN")
                BADLEN += 1
                continue
                    
            #check polyX
            if self.options.poly_size_limit > 0:
                poly1 = hasPolyX(r1[1], self.options.poly_size_limit, self.options.allow_mismatch_in_poly)
                poly2 = None
                if r2!=None:
                    poly2 = hasPolyX(r2[1], self.options.poly_size_limit, self.options.allow_mismatch_in_poly)
                if poly1!=None or poly2!=None:
                    self.writeReads(r1, r2, i1, i2, bad_read1_file, bad_read2_file, bad_index1_file, bad_index2_file, "BADPOL")
                    BADPOL += 1
                    continue
            
            #check low quality count
            if self.options.unqualified_base_limit > 0:
                lowQual1 = lowQualityNum(r1, self.options.qualified_quality_phred)
                lowQual2 = 0
                if r2!=None:
                    lowQual2 = lowQualityNum(r2, self.options.qualified_quality_phred)
                if lowQual1 > self.options.unqualified_base_limit or lowQual1 > self.options.unqualified_base_limit:
                    self.writeReads(r1, r2, i1, i2, bad_read1_file, bad_read2_file, bad_index1_file, bad_index2_file, "BADLQC")
                    BADLQC += 1
                    continue
            
            #check N number
            if self.options.n_base_limit > 0:
                nNum1 = nNumber(r1)
                nNum2 = 0
                if r2!=None:
                    nNum2 = nNumber(r2)
                if nNum1 > self.options.n_base_limit or nNum2 > self.options.n_base_limit:
                    self.writeReads(r1, r2, i1, i2, bad_read1_file, bad_read2_file, bad_index1_file, bad_index2_file, "BADNCT")
                    BADNCT += 1
                    continue

            #check overlap and do error correction
            if r2!=None and (not self.options.no_overlap):
                (offset, overlap_len, distance) = util.overlap(r1[1], r2[1])
                overlap_histgram[overlap_len] += 1
                # deal with the case insert DNA is shorter than read length and cause offset is negative
                # in this case the adapter is sequenced and should be trimmed
                if offset <0 and overlap_len > 30:
                    # shift the junk bases
                    r1[1] = r1[1][0:overlap_len]
                    r1[3] = r1[3][0:overlap_len]
                    r2[1] = r2[1][0:overlap_len]
                    r2[3] = r2[3][0:overlap_len]
                    TRIMMED_ADAPTER_BASE += abs(offset)*2
                    # check the sequence length again after adapter trimmed
                    if len(r1[1])<self.options.seq_len_req:
                        self.writeReads(r1, r2, i1, i2, bad_read1_file, bad_read2_file, bad_index1_file, bad_index2_file, "BADLEN")
                        BADLEN += 1
                        continue
                    # then calc overlap again
                    (offset, overlap_len, distance) = util.overlap(r1[1], r2[1])
                if overlap_len>30:
                    OVERLAPPED += 1
                    distance_histgram[distance] += 1
                    OVERLAP_LEN_SUM += overlap_len
                    # we consider the distance is caused by sequencing error
                    OVERLAP_BASE_SUM += overlap_len * 2
                    OVERLAP_BASE_ERR += distance
                    corrected = 0
                    zero_qual_masked = 0
                    skipped_mismatch = 0
                    if distance>0:
                        #try to fix low quality base
                        hamming = util.hammingDistance(r1[1][len(r1[1]) - overlap_len:], util.reverseComplement(r2[1][len(r2[1]) - overlap_len:]))
                        if hamming != distance:
                            self.writeReads(r1, r2, i1, i2, bad_read1_file, bad_read2_file, bad_index1_file, bad_index2_file, "BADINDEL")
                            BADINDEL += 1
                            continue
                        #print(r1[1][len(r1[1]) - overlap_len:])
                        #print(util.reverseComplement(r2[1][len(r2[1]) - overlap_len:]))
                        #print(r1[3][len(r1[1]) - overlap_len:])
                        #print(util.reverse(r2[3][len(r2[1]) - overlap_len:]))
                        err_mtx = init_error_matrix()
                        for o in xrange(overlap_len):
                            b1 = r1[1][len(r1[1]) - overlap_len + o]
                            b2 = util.complement(r2[1][-o-1])
                            q1 = r1[3][len(r1[3]) - overlap_len + o]
                            q2 = r2[3][-o-1]
                            if b1 != b2:
                                # print(TOTAL_READS, o, b1, b2, q1, q2)
                                this_is_corrected = False
                                if util.qualNum(q1) >= 30 and util.qualNum(q2) <= 14:
                                    if b1!='N' and b2!='N':
                                        err_mtx[util.complement(b1)][util.complement(b2)] += 1
                                    if not self.options.no_correction:
                                        r2[1] = util.changeString(r2[1], -o-1, util.complement(b1))
                                        r2[3] = util.changeString(r2[3], -o-1, q1)
                                        corrected += 1
                                        this_is_corrected = True
                                elif util.qualNum(q2) >= 30 and util.qualNum(q1) <= 14:
                                    if b1!='N' and b2!='N':
                                        err_mtx[b2][b1] += 1
                                    if not self.options.no_correction:
                                        r1[1]= util.changeString(r1[1], len(r1[1]) - overlap_len + o, b2)
                                        r1[3] = util.changeString(r1[3], len(r1[3]) - overlap_len + o, q2)
                                        corrected += 1
                                        this_is_corrected = True
                                if not this_is_corrected:
                                    if self.options.mask_mismatch:
                                        # mask them as zero qual if it is not corrected
                                        zero_qual = '!'
                                        r2[3] = util.changeString(r2[3], -o-1, zero_qual)
                                        r1[3] = util.changeString(r1[3], len(r1[3]) - overlap_len + o, zero_qual)
                                        zero_qual_masked += 1
                                    else:
                                        skipped_mismatch += 1

                                if corrected + zero_qual_masked + skipped_mismatch >= distance:
                                    break
                        #print(r1[1][len(r1[1]) - overlap_len:])
                        #print(util.reverseComplement(r2[1][len(r2[1]) - overlap_len:]))
                        #print(r1[3][len(r1[1]) - overlap_len:])
                        #print(util.reverse(r2[3][len(r2[1]) - overlap_len:]))
                        if corrected + zero_qual_masked + skipped_mismatch == distance:
                            merge_error_matrix(OVERLAP_ERR_MATRIX, err_mtx)
                            if corrected > 0:
                                READ_CORRECTED += 1
                            BASE_CORRECTED += corrected
                            # multiply by 2 since we mask bases by pair
                            BASE_ZERO_QUAL_MASKED += zero_qual_masked * 2
                        else:
                            self.writeReads(r1, r2, i1, i2, bad_read1_file, bad_read2_file, bad_index1_file, bad_index2_file, "BADMISMATCH")
                            BADMISMATCH += 1
                            continue
                    if distance == 0 or distance == corrected:
                        if self.options.store_overlap:
                            self.writeReads(getOverlap(r1, overlap_len), getOverlap(r2, overlap_len), i1, i2, overlap_read1_file, overlap_read2_file, overlap_index1_file, overlap_index2_file, None)

            #write to good       
            self.writeReads(r1, r2, i1, i2, good_read1_file, good_read2_file, good_index1_file, good_index2_file, None)
            GOOD_BASES += len(r1[1])
            if i2 != None:
                GOOD_BASES += len(r2[1])
            if self.options.qc_sample <=0 or TOTAL_READS < self.options.qc_sample:
                self.r1qc_postfilter.statRead(r1)
                if r2 != None:
                    self.r2qc_postfilter.statRead(r2)

            GOOD_READS += 1
            if self.options.qc_only and TOTAL_READS >= self.options.qc_sample:
                break

        self.r1qc_postfilter.qc()
        #self.r1qc_postfilter.plot(qc_dir, "R1-postfilter")
        if self.options.read2_file != None:
            self.r2qc_postfilter.qc()
            #self.r2qc_postfilter.plot(qc_dir, "R2-postfilter")
        
        #close all files
        if not self.options.qc_only:
            good_read1_file.flush()
            bad_read1_file.flush()
            if self.options.read2_file != None:
                good_read2_file.flush()
                bad_read2_file.flush()
            if self.options.index1_file != None:
                good_index1_file.flush()
                bad_index1_file.flush()
            if self.options.index2_file != None:
                good_index2_file.flush()
                bad_index2_file.flush()

        # print stat numbers
        BAD_READS = TOTAL_READS - GOOD_READS
        result = {}
        result['total_bases']=TOTAL_BASES
        result['good_bases']=GOOD_BASES
        result['total_reads']=TOTAL_READS
        result['good_reads']=GOOD_READS
        result['bad_reads']=BAD_READS
        result['bad_reads_with_bad_barcode']= BADBCD1 + BADBCD2
        result['bad_reads_with_reads_in_bubble']= BADBBL
        result['bad_reads_with_bad_read_length']= BADLEN + BADTRIM1 + BADTRIM2
        result['bad_reads_with_polyX']= BADPOL
        result['bad_reads_with_low_quality']=BADLQC
        result['bad_reads_with_too_many_N']= BADNCT
        result['bad_reads_with_bad_overlap']= BADMISMATCH + BADINDEL
        result['readlen'] = readLen

        # plot result bar figure
        labels = ['good reads', 'has_polyX', 'low_quality', 'too_short', 'too_many_N']
        counts = [GOOD_READS, BADPOL, BADLQC, BADLEN + BADTRIM1 + BADTRIM2, BADNCT]
        colors = ['#66BB11', '#FF33AF', '#FFD3F2', '#FFA322', '#FF8899']
        if self.options.read2_file != None:
            labels.append('bad_overlap')
            counts.append(BADMISMATCH + BADINDEL)
            colors.append('#FF6600')
        if self.options.debubble:
            labels.append('in_bubble')
            counts.append(BADBBL)
            colors.append('#EEBB00')
        if self.options.barcode:
            labels.append('bad_barcode')
            counts.append(BADBCD1 + BADBCD2)
            colors.append('#CCDD22')

        for i in xrange(len(counts)):
            type_percent = 0.0
            if TOTAL_READS > 0:
                type_percent = 100.0 * float(counts[i])/TOTAL_READS
            labels[i] = labels[i] + ": " + str(counts[i]) + "(" + str(type_percent) + "%)"

        reporter.addFigure('Good reads and bad reads after filtering', self.r1qc_prefilter.statPlotly(labels, counts, TOTAL_READS, 'filter_stat'), 'filter_stat', "")
        #self.r1qc_prefilter.plotFilterStats(labels, counts, colors, TOTAL_READS, os.path.join(qc_dir, "filter-stat.png"))

        stat={}
        # stat["options"]=self.options
        stat["summary"]=result
        stat["command"]=makeDict(self.options)
        stat["kmer_content"] = {}
        stat["kmer_content"]["read1_prefilter"] = self.r1qc_prefilter.topKmerCount[0:10]
        stat["kmer_content"]["read1_postfilter"] = self.r1qc_postfilter.topKmerCount[0:10]
        if self.options.read2_file != None:
            stat["kmer_content"]["read2_prefilter"] = self.r2qc_prefilter.topKmerCount[0:10]
            stat["kmer_content"]["read2_postfilter"] = self.r2qc_postfilter.topKmerCount[0:10]
            stat["overlap"]={}
            stat["overlap"]['overlapped_pairs']=OVERLAPPED
            if OVERLAPPED > 0:
                stat["overlap"]['average_overlap_length']=float(OVERLAP_LEN_SUM/OVERLAPPED)
            else:
                stat["overlap"]['average_overlap_length']=0.0
            stat["overlap"]['bad_mismatch_reads']=BADMISMATCH
            stat["overlap"]['bad_indel_reads']=BADINDEL
            stat["overlap"]['corrected_reads']=READ_CORRECTED
            stat["overlap"]['corrected_bases']=BASE_CORRECTED
            stat["overlap"]['zero_qual_masked']=BASE_ZERO_QUAL_MASKED
            stat["overlap"]['trimmed_adapter_bases']=TRIMMED_ADAPTER_BASE
            if OVERLAP_BASE_SUM > 0:
                stat["overlap"]['error_rate']=float(OVERLAP_BASE_ERR)/float(OVERLAP_BASE_SUM)
            else:
                stat["overlap"]['error_rate']=0.0
            stat["overlap"]['error_matrix']=OVERLAP_ERR_MATRIX
            stat["overlap"]['edit_distance_histogram']=distance_histgram[0:10]
            reporter.addFigure('Sequence error distribution', self.r1qc_prefilter.errorPlotly(OVERLAP_ERR_MATRIX, 'error_matrix'), 'error_matrix', "")
            reporter.addFigure('Overlap length distribution', self.r1qc_prefilter.overlapPlotly(overlap_histgram, readLen, TOTAL_READS, 'overlap_stat'), 'overlap_stat', "")
            #self.r1qc_prefilter.plotOverlapHistgram(overlap_histgram, readLen, TOTAL_READS, os.path.join(qc_dir, "overlap.png"))

        stat_file = open(os.path.join(qc_dir, os.path.basename(self.options.read1_file) + ".json"), "w")
        stat_json = json.dumps(stat, sort_keys=True,indent=4, separators=(',', ': '))
        stat_file.write(stat_json)
        stat_file.close()

        self.addFiguresToReport(reporter)
        reporter.setStat(stat)
        reporter.setVersion(self.options.version)
        reporter.output(os.path.join(qc_dir, os.path.basename(self.options.read1_file) + ".html"))

    def addFiguresToReport(self, reporter):
        if self.options.read2_file != None:
            reporter.addFigure('Read1 quality curve before filtering', self.r1qc_prefilter.qualityPlotly("r1_pre_quality", 'Read1 quality curve before filtering'), "r1_pre_quality", "")
            reporter.addFigure('Read1 base content distribution before filtering', self.r1qc_prefilter.contentPlotly("r1_pre_content", 'Read1 base content distribution before filtering'), "r1_pre_content", "")
            reporter.addFigure('Read1 GC curve before filtering', self.r1qc_prefilter.gcPlotly("r1_pre_gc", 'Read1 GC curve before filtering'), 'r1_pre_gc', "")
            reporter.addFigure('Read1 per base discontinuity before filtering', self.r1qc_prefilter.discontinuityPlotly("r1_pre_discontinuity", 'Read1 discontinuity curve before filtering'), 'r1_pre_discontinuity', "")
            reporter.addFigure('Read1 kmer strand bias before filtering', self.r1qc_prefilter.strandBiasPlotly("r1_pre_sb", 'Read1 Kmer strand bias before filtering'), 'r1_pre_sb', "")
            reporter.addFigure('Read1 quality curve after filtering', self.r1qc_postfilter.qualityPlotly("r1_post_quality", 'Read1 quality curve after filtering'), "r1_post_quality", "")
            reporter.addFigure('Read1 base content distribution after filtering', self.r1qc_postfilter.contentPlotly("r1_post_content", 'Read1 base content distribution after filtering'), "r1_post_content", "")
            reporter.addFigure('Read1 GC curve after filtering', self.r1qc_postfilter.gcPlotly("r1_post_gc", 'Read1 GC curve after filtering'), 'r1_post_gc', "")
            reporter.addFigure('Read1 per base discontinuity after filtering', self.r1qc_postfilter.discontinuityPlotly("r1_post_discontinuity", 'Read1 discontinuity curve after filtering'), 'r1_post_discontinuity', "")
            reporter.addFigure('Read1 kmer strand bias after filtering', self.r1qc_postfilter.strandBiasPlotly("r1_post_sb", 'Read1 kmer strand bias after filtering'), 'r1_post_sb', "")

            reporter.addFigure('Read2 quality curve before filtering', self.r2qc_prefilter.qualityPlotly("r2_pre_quality", 'Read2 quality curve before filtering'), "r2_pre_quality", "")
            reporter.addFigure('Read2 base content distribution before filtering', self.r2qc_prefilter.contentPlotly("r2_pre_content", 'Read2 base content distribution before filtering'), "r2_pre_content", "")
            reporter.addFigure('Read2 GC curve before filtering', self.r2qc_prefilter.gcPlotly("r2_pre_gc", 'Read2 GC curve before filtering'), 'r2_pre_gc', "")
            reporter.addFigure('Read2 per base discontinuity before filtering', self.r2qc_prefilter.discontinuityPlotly("r2_pre_discontinuity", 'Read2 discontinuity curve before filtering'), 'r2_pre_discontinuity', "")
            reporter.addFigure('Read2 kmer strand bias before filtering', self.r2qc_prefilter.strandBiasPlotly("r2_pre_sb", 'Read2 Kmer strand bias before filtering'), 'r2_pre_sb', "")
            reporter.addFigure('Read2 quality curve after filtering', self.r2qc_postfilter.qualityPlotly("r2_post_quality", 'Read2 quality curve after filtering'), "r2_post_quality", "")
            reporter.addFigure('Read2 base content distribution after filtering', self.r2qc_postfilter.contentPlotly("r2_post_content", 'Read2 base content distribution after filtering'), "r2_post_content", "")
            reporter.addFigure('Read2 GC curve after filtering', self.r2qc_postfilter.gcPlotly("r2_post_gc", 'Read2 GC curve after filtering'), 'r2_post_gc', "")
            reporter.addFigure('Read2 per base discontinuity after filtering', self.r2qc_postfilter.discontinuityPlotly("r2_post_discontinuity", 'Read2 discontinuity curve after filtering'), 'r2_post_discontinuity', "")
            reporter.addFigure('Read2 kmer strand bias after filtering', self.r2qc_postfilter.strandBiasPlotly("r2_post_sb", 'Read2 kmer strand bias after filtering'), 'r2_post_sb', "")

        else:
            reporter.addFigure('Quality curve before filtering', self.r1qc_prefilter.qualityPlotly("r1_pre_quality", 'Quality curve before filtering'), "r1_pre_quality", "")
            reporter.addFigure('Base content distribution before filtering', self.r1qc_prefilter.contentPlotly("r1_pre_content", 'Base content distribution before filtering'), "r1_pre_content", "")
            reporter.addFigure('GC curve before filtering', self.r1qc_prefilter.gcPlotly("r1_pre_gc", 'GC curve before filtering'), 'r1_pre_gc', "")
            reporter.addFigure('Per base discontinuity before filtering', self.r1qc_prefilter.discontinuityPlotly("r1_pre_discontinuity", 'Discontinuity curve before filtering'), 'r1_pre_discontinuity', "")
            reporter.addFigure('Kmer strand bias before filtering', self.r1qc_prefilter.strandBiasPlotly("r1_pre_sb", 'Kmer strand bias before filtering'), 'r1_pre_sb', "")
            reporter.addFigure('Quality curve after filtering', self.r1qc_postfilter.qualityPlotly("r1_post_quality", 'Quality curve after filtering'), "r1_post_quality", "")
            reporter.addFigure('Base content distribution after filtering', self.r1qc_postfilter.contentPlotly("r1_post_content", 'Base content distribution after filtering'), "r1_post_content", "")
            reporter.addFigure('GC curve after filtering', self.r1qc_postfilter.gcPlotly("r1_post_gc", 'GC curve after filtering'), 'r1_post_gc', "")
            reporter.addFigure('Per base discontinuity after filtering', self.r1qc_postfilter.discontinuityPlotly("r1_post_discontinuity", 'Discontinuity curve after filtering'), 'r1_post_discontinuity', "")
            reporter.addFigure('Kmer strand bias after filtering', self.r1qc_postfilter.strandBiasPlotly("r1_post_sb", 'Kmer strand bias after filtering'), 'r1_post_sb', "")


