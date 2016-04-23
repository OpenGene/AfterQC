import os,sys
from optparse import OptionParser
import time
import fastq
import util

maxLen = 1000
allbases = ("A", "T", "C", "G");

########################### QualityControl
class QualityControl:
    readLen = 0
    readCount = 0
    counts = {}
    percents = {}
    totalQual = [0 for x in xrange(maxLen)]
    totalNum = [0 for x in xrange(maxLen)]
    meanQual = [0.0 for x in xrange(maxLen)]

    def __init__(self):
        for base in allbases:
            self.counts[base] = [0 for x in xrange(maxLen)]
            self.percents[base] = [0.0 for x in xrange(maxLen)]

    def statRead(self, read):
        seq = read[1]
        qual = read[3]
        for i in xrange(len(seq)):
            self.totalNum[i] += 1
            self.totalQual[i] += util.qualNum(qual[i])
            b = seq[i]
            if b in allbases:
                self.counts[b][i] += 1

    def calcReadLen(self):
        for pos in xrange(maxLen):
            hasData = False
            for base in allbases:
                if self.counts[base][pos]>0:
                    hasData = True
            if hasData == False:
                self.readLen = pos
                break

    def calcPercents(self):
        #calc percents of each base
        for pos in xrange(self.readLen):
            total = 0
            for base in allbases:
                total += self.counts[base][pos]
            for base in allbases:
                self.percents[base][pos] = float(self.counts[base][pos])/float(total)

    def calcQualities(self):
        for pos in xrange(self.readLen):
            self.meanQual[pos] = float(self.totalQual[pos])/float(self.totalNum[pos])

    def qc(self):    
        self.calcReadLen()
        self.calcPercents()
        self.calcQualities()
        print(self.meanQual)
        
    def statFile(self, filename):
        reader = fastq.Reader(filename)
        #sample up to maxSample reads for stat
        maxSample = 50000000
        while True:
            read = reader.nextRead()
            self.readCount += 1
            if read==None or self.readCount>maxSample:
                break
            self.statRead(read)

        self.qc()

    def autoTrim(self):
        #use (center-5, center+5) as initial good segment        
        center = int(self.readLen/2)
        left = center-5
        right = center+5
        
        threshold = 0.05
        lastStepIsLeft = False
        leftFinished = False
        rightFinished = False
        current = -1
        
        #expand the good segment
        meanPercents = {}
        while not (leftFinished and rightFinished):
            for base in allbases:
                meanPercents[base] = 0.0
                for pos in xrange(left, right):
                    meanPercents[base] += self.percents[base][pos]
                meanPercents[base] /= (right-left)
            
            if leftFinished:
                current = right + 1
                lastStepIsLeft = False
            elif rightFinished:
                current = left - 1
                lastStepIsLeft = True
            elif lastStepIsLeft:
                current = right + 1
                lastStepIsLeft = False
            else:
                current = left - 1
                lastStepIsLeft = True
                                
            percentBias = 0.0
            for base in allbases:
                percentBias += abs(meanPercents[base] - self.percents[base][current])
            
            if percentBias > threshold:
                if lastStepIsLeft:
                    leftFinished = True
                else:
                    rightFinished = True
            else:
                if lastStepIsLeft:
                    left = current
                    if left == 0: leftFinished = True
                else:
                    right = current
                    if right == self.readLen-1: rightFinished = True
                    
        #find the bad segment from front, considering a small window
        #if any in the window is bad, it is bad
        trimFront = left
        window = 3
        for pos in xrange(0, left):
            isGood = True
            for posInWindow in xrange(pos, min(pos+3, self.readLen)):
                percentBias = 0.0
                for base in allbases:
                    percentBias += abs(meanPercents[base] - self.percents[base][posInWindow])
                if percentBias > threshold:
                    isGood = False
            if isGood:    
                trimFront = pos
                break
        #find the bad segment from tail, considering a small window
        #if any in the window is bad, it is bad
        trimTail = right
        for pos in xrange(self.readLen-1, right, -1):
            isGood = True
            for posInWindow in xrange(pos, max(pos-3, 0), -1):
                percentBias = 0.0
                for base in allbases:
                    percentBias += abs(meanPercents[base] - self.percents[base][posInWindow])
                if percentBias > threshold:
                    isGood = False
            if isGood: 
                trimTail = pos
                break
        
        trimFront = min(self.readLen*0.1,trimFront)
        trimTail = min(self.readLen*0.05,self.readLen-1-trimTail)
        # the last base should be definitely trimmed for illumina sequencer output
        trimTail = max(1, trimTail)
        
        return (int(trimFront), int(trimTail))

if __name__  == "__main__":
    qc = QualityControl()
    qc.statFile("R1.fq")
    print(qc.autoTrim())
