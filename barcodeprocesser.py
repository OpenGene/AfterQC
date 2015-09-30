 #!/usr/bin/env python
 
import os,sys
from optparse import OptionParser
import time

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

#test
if __name__  == "__main__":
    read = ['@NS500713:17:HFG2YBGXX:1:11101:10560:1202 1:N:0:CGAGTA','ATAAAAAAAACACAGTATGGCAAAACCCCATCTCTACTAAAAATACAAAAATTAGCTGGGTGTGGTGGCNNNNNN','+','AAAAAEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEAEEEEEE######']
    verify = 'CAGTA'
    barcodeLen = 12
    moveBarcodeToName(read, barcodeLen, verify)
    print(read)