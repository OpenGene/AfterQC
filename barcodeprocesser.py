 #!/usr/bin/env python
 
import os,sys
from optparse import OptionParser
import time
import util

#how many chars are different of these two strings with same length
def diffNumber(str1, str2):
    diff = 0
    for i in xrange(len(str1)):
        if str1[i]!=str2[i]:
            diff += 1
    return diff

#find the place where barcode actually ends
#barcodeLen is the design length of barcode
#if not found, return 0
def detectBarcode(seq, barcodeLen, verify):
    verifyLen = len(verify)
    if len(seq) <= verifyLen + barcodeLen + 1:
        return 0
    center = seq[barcodeLen:barcodeLen+verifyLen]
    if diffNumber(center, verify) <= 1:
        return barcodeLen
    left = seq[barcodeLen-1:barcodeLen-1+verifyLen]
    if diffNumber(left, verify) == 0:
        return barcodeLen-1
    right = seq[barcodeLen+1:barcodeLen+1+verifyLen]
    if diffNumber(right, verify) == 0:
        return barcodeLen+1
    return 0
        
def moveBarcodeToName(read, barcodeLen, verify):
    verifyLen = len(verify)
    #trick: we stores the barcode in <instrument> to not break the sam/bam format
    #illumina sequence name line format
    #@<instrument>:<run number>:<flowcell ID>:<lane>:<tile_no>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<index sequence>
    colonPos = read[0].find(':')
    barcode = read[1][0:barcodeLen]
    removeLen = verifyLen + barcodeLen
    read[0] = '@' + barcode + read[0][colonPos:]
    read[1] = read[1][removeLen:]
    read[3] = read[3][removeLen:]
    return barcode

def cleanBarcodeTail(read1, read2, readStart1, readStart2):
    reverse1 = util.reverseComplement(readStart1)
    reverse2 = util.reverseComplement(readStart2)
    barcodeStringLen = min(len(readStart1), len(readStart2))
    r1len  = len(read1[1])
    r2len  = len(read2[1])
    compLen = 0
    overlap = False
    for i in xrange(barcodeStringLen):
        compLen = barcodeStringLen - i
        if compLen >= r1len or compLen >= r2len:
            continue
        distance1 = util.editDistance(read1[1][-compLen:], reverse2[i:])
        distance2 = util.editDistance(read2[1][-compLen:], reverse1[i:])
        #if the tail of one end matches the start of the other end
        #we then suspect the template is not shorter than the read length
        #so to trim them on both ends
        threshold  = compLen/5
        if distance1<=threshold  and distance2<=threshold:
            read1[1] = read1[1][:-compLen]
            read1[3] = read1[3][:-compLen]
            read2[1] = read2[1][:-compLen]
            read2[3] = read2[3][:-compLen]
            overlap = True
            break;
    if overlap:
        return compLen
    else:
        return 0

def moveAndTrimPair(read1, read2, barcode1len, barcode2len, verify):
    readStart1  = read1[1][0:barcode1len] + verify
    readStart2  = read2[1][0:barcode2len] + verify
    moveBarcodeToName(read1, barcode1len, verify)
    moveBarcodeToName(read2, barcode2len, verify)
    cleanBarcodeTail(read1, read2, readStart1, readStart2)
