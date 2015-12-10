 #!/usr/bin/env python
 
import os,sys
from optparse import OptionParser
import time
import fastq
    
########################### Trimmer
class Trimmer:
        
    def calcTrimLength(self, filename):
        maxLen = 1000
        allbases = ("A", "T", "C", "G");
        counts = {}
        percents = {}
        for base in allbases:
            counts[base] = [0 for x in xrange(maxLen)]
            percents[base] = [0.0 for x in xrange(maxLen)]

        reader = fastq.Reader(filename)
        #sample up to maxSample reads for stat
        maxSample = 5000000
        sampled = 0
        while True:
            read = reader.nextRead()
            sampled += 1
            if read==None or sampled>maxSample:
                break
            seq = read[1]
            for i in xrange(len(seq)):
                b = seq[i]
                if b in allbases:
                    counts[b][i] += 1
                
        #get the length of read
        readLen = 0
        for pos in xrange(maxLen):
            hasData = False
            for base in allbases:
                if counts[base][pos]>0:
                    hasData = True
            if hasData == False:
                readLen = pos
                break
                
        #calc percents of each base
        for pos in xrange(readLen):
            total = 0
            for base in allbases:
                total += counts[base][pos]
            for base in allbases:
                percents[base][pos] = float(counts[base][pos])/float(total)
        
        #use (center-5, center+5) as initial good segment        
        center = int(readLen/2)
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
                    meanPercents[base] += percents[base][pos]
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
                percentBias += abs(meanPercents[base] - percents[base][current])
            
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
                    if right == readLen-1: rightFinished = True
                    
        #find the bad segment from front, considering a small window
        #if any in the window is bad, it is bad
        trimFront = left
        window = 3
        for pos in xrange(0, left):
            isGood = True
            for posInWindow in xrange(pos, min(pos+3, readLen)):
                percentBias = 0.0
                for base in allbases:
                    percentBias += abs(meanPercents[base] - percents[base][posInWindow])
                if percentBias > threshold:
                    isGood = False
            if isGood:    
                trimFront = pos
                break
        #find the bad segment from tail, considering a small window
        #if any in the window is bad, it is bad
        trimTail = right
        for pos in xrange(readLen-1, right, -1):
            isGood = True
            for posInWindow in xrange(pos, max(pos-3, 0), -1):
                percentBias = 0.0
                for base in allbases:
                    percentBias += abs(meanPercents[base] - percents[base][posInWindow])
                if percentBias > threshold:
                    isGood = False
            if isGood: 
                trimTail = pos
                break
        
        trimFront = min(readLen*0.1,trimFront)
        trimTail = min(readLen*0.05,readLen-1-trimTail)
        
        return (int(trimFront), int(trimTail))

