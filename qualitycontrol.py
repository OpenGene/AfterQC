import os,sys
from optparse import OptionParser
import time
import fastq
import util
import matplotlib
# fix matplotlib DISPLAY issue
matplotlib.use('Agg')
import matplotlib.pyplot as plt

MAX_LEN = 200
ALL_BASES = ("A", "T", "C", "G");
KMER_TOP = 10

########################### QualityControl
class QualityControl:

    def __init__(self, qc_sample=1000000, qc_kmer=8):
        self.sampleLimit = qc_sample
        self.kmerLen = qc_kmer
        self.readLen = 0
        self.readCount = 0
        self.baseCounts = {}
        self.percents = {}
        self.baseTotalQual = {}
        self.baseMeanQual = {}
        self.totalQual = [0 for x in xrange(MAX_LEN)]
        self.totalNum = [0 for x in xrange(MAX_LEN)]
        self.meanQual = [0.0 for x in xrange(MAX_LEN)]
        self.gcPercents = [0.0 for x in xrange(MAX_LEN)]
        self.gcHistgram = [0 for x in xrange(100+1)]
        self.kmerCount = {}
        self.topKmerCount = []
        self.totalKmer = 0
        for base in ALL_BASES:
            self.baseCounts[base] = [0 for x in xrange(MAX_LEN)]
            self.percents[base] = [0.0 for x in xrange(MAX_LEN)]
            self.baseMeanQual[base] = [0.0 for x in xrange(MAX_LEN)]
            self.baseTotalQual[base] = [0 for x in xrange(MAX_LEN)]

    def statRead(self, read):
        seq = read[1]
        qual = read[3]
        seqlen = len(seq)
        gc = 0
        for i in xrange(seqlen):
            self.totalNum[i] += 1
            qnum = util.qualNum(qual[i])
            self.totalQual[i] += qnum
            b = seq[i]
            if b=='G' or b=='C':
                gc += 1
            if b in ALL_BASES:
                self.baseCounts[b][i] += 1
                self.baseTotalQual[b][i] += qnum

        gcPer = int(100.0* float(gc)/seqlen)
        self.gcHistgram[gcPer] += 1
        for i in xrange(seqlen - self.kmerLen):
            self.totalKmer += 1
            kmer = seq[i:i+self.kmerLen]
            if kmer in self.kmerCount:
                self.kmerCount[kmer] += 1
            else:
                self.kmerCount[kmer] = 1

    def calcReadLen(self):
        for pos in xrange(MAX_LEN):
            hasData = False
            for base in ALL_BASES:
                if self.baseCounts[base][pos]>0:
                    hasData = True
            if hasData == False:
                self.readLen = pos
                break

    def calcPercents(self):
        #calc percents of each base
        for pos in xrange(self.readLen):
            total = 0
            for base in ALL_BASES:
                total += self.baseCounts[base][pos]
            for base in ALL_BASES:
                self.percents[base][pos] = float(self.baseCounts[base][pos])/float(total)
                self.gcPercents[pos] = float(self.baseCounts['G'][pos] + self.baseCounts['C'][pos])/float(total)

    def calcQualities(self):
        for pos in xrange(self.readLen):
            self.meanQual[pos] = float(self.totalQual[pos])/float(self.totalNum[pos])
            for base in ALL_BASES:
                if self.baseCounts[base][pos] > 0:
                    self.baseMeanQual[base][pos] = float(self.baseTotalQual[base][pos])/float(self.baseCounts[base][pos])

    def sortKmer(self):
        self.topKmerCount = sorted(self.kmerCount.items(), key=lambda x: x[1], reverse=True)

    def plotQuality(self, filename, prefix=""):
        colors = {'A':'red', 'T':'purple', 'C':'blue', 'G':'green'}
        x = range(self.readLen)
        plt.figure(1)
        plt.title(prefix + " base quality" )
        plt.xlim(0, self.readLen)
        plt.ylabel('Quality')
        plt.xlabel('Cycle')
        for base in ALL_BASES:
            plt.plot(x, self.baseMeanQual[base][0:self.readLen], color = colors[base], label=base, alpha=0.3)
        plt.plot(x, self.meanQual[0:self.readLen], color = 'black', label = 'mean')
        plt.legend(loc='upper right', ncol=5)
        plt.savefig(filename)
        plt.close(1)

    def plotContent(self, filename, prefix=""):
        colors = {'A':'red', 'T':'purple', 'C':'blue', 'G':'green'}
        x = range(self.readLen)
        plt.figure(1)
        plt.title(prefix + " base contents" )
        plt.xlim(0, self.readLen)
        plt.ylim(0.0, 0.8)
        plt.ylabel('Percents')
        plt.xlabel('Cycle')
        for base in ALL_BASES:
            plt.plot(x, self.percents[base][0:self.readLen], color = colors[base], label=base, alpha=0.5)
        plt.plot(x, self.gcPercents[0:self.readLen], color = 'black', label='GC')
        plt.legend(loc='upper right', ncol=5)
        plt.savefig(filename)
        plt.close(1)

    def plotGCHistogram(self, filename, prefix=""):
        x = range(100+1)
        plt.figure(1)
        plt.title(prefix + " GC content distribution" )
        plt.ylabel('Count')
        plt.xlabel('GC percentage (%)')
        plt.bar(x, self.gcHistgram, color = 'gray', label='Actual', alpha=0.8)
        # plt.legend(loc='upper right', ncol=5)
        plt.savefig(filename)
        plt.close(1)

    def plot(self, folder=".", prefix=""):
        self.plotQuality(os.path.join(folder, prefix + "-quality.png"), prefix)
        self.plotContent(os.path.join(folder, prefix + "-content.png"), prefix)
        self.plotGCHistogram(os.path.join(folder, prefix + "-gc-curve.png"), prefix)

    def qc(self): 
        self.calcReadLen()
        self.calcPercents()
        self.calcQualities()
        self.sortKmer()
        
    def statFile(self, filename):
        reader = fastq.Reader(filename)
        #sample up to maxSample reads for stat
        while True:
            read = reader.nextRead()
            self.readCount += 1
            if read==None or self.readCount > self.sampleLimit:
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
            for base in ALL_BASES:
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
            for base in ALL_BASES:
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
                for base in ALL_BASES:
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
                for base in ALL_BASES:
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
    qc.plot()
    print(qc.autoTrim())
    print(qc.topKmerCount[0:10])
