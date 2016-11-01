import os,sys
from optparse import OptionParser
import time
import fastq
import util

HAVE_MATPLOTLIB = True
WARNED_PLOT_FAILURE = False

if not sys.executable.endswith("pypy"):
    try:
        import matplotlib
        # fix matplotlib DISPLAY issue
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except Exception:
        HAVE_MATPLOTLIB = False

if sys.executable.endswith("pypy"):
    HAVE_MATPLOTLIB = False

MAX_LEN = 1000
ALL_BASES = ("A", "T", "C", "G");
KMER_TOP = 10

def makeRange(bottom, top):
    return "[" + str(bottom) + "," + str(top) + "]"

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
        self.gcHistogram = [0 for x in xrange(MAX_LEN)]
        self.kmerCount = {}
        self.topKmerCount = []
        self.totalKmer = 0
        self.meanDiscontinuity = [0.0 for x in xrange(MAX_LEN)]
        self.totalDiscontinuity = [0.0 for x in xrange(MAX_LEN)]
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

            # calculate discontinuity
            left = i-2
            right = i+3
            if left<0:
                left = 0
                right = 5
            elif right >= seqlen:
                right = seqlen
                left = seqlen - 5
            discontinuity = 0
            for j in xrange(left, right-1):
                if seq[j] != seq[j+1]:
                    discontinuity += 1
            self.totalDiscontinuity[i] += discontinuity

        #gcPer = int(1000.0* float(gc)/seqlen)
        self.gcHistogram[gc] += 1
        for i in xrange(seqlen - self.kmerLen):
            self.totalKmer += 1
            kmer = seq[i:i+self.kmerLen]
            if kmer in self.kmerCount:
                self.kmerCount[kmer] += 1
            else:
                self.kmerCount[kmer] = 1
                rcKmer = util.reverseComplement(kmer)
                if rcKmer not in self.kmerCount:
                    self.kmerCount[rcKmer] = 0

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

    def calcDiscontinuity(self):
        for pos in xrange(self.readLen):
            self.meanDiscontinuity[pos] = float(self.totalDiscontinuity[pos])/float(self.totalNum[pos])

    def sortKmer(self):
        self.topKmerCount = sorted(self.kmerCount.items(), key=lambda x: x[1], reverse=True)

    def qualityPlotly(self, div, title=""):
        colors = {'A':'rgba(255,0,0,0.5)', 'T':'rgba(128,0,128,0.5)', 'C':'rgba(0,255,0,0.5)', 'G':'rgba(0,0,255,0.5)'}
        json_str = "var data=["
        x = range(self.readLen)
        # four bases
        for base in ALL_BASES:
            json_str += "{"
            json_str += "x:[" + ",".join(map(str, x)) + "],"
            json_str += "y:[" + ",".join(map(str, self.baseMeanQual[base][0:self.readLen])) + "],"
            json_str += "name: '" + base + "',"
            json_str += "mode:'lines',"
            json_str += "line:{color:'" + colors[base] + "', width:1}\n"
            json_str += "},"
        # mean
        json_str += "{"
        json_str += "x:[" + ",".join(map(str, x)) + "],"
        json_str += "y:[" + ",".join(map(str, self.meanQual[0:self.readLen])) + "],"
        json_str += "name: 'mean',"
        json_str += "mode:'lines',"
        json_str += "line:{color:'rgba(20,20,20,255)', width:1}\n"
        json_str += "}\n"
        json_str += "];\n"
        json_str += "var layout={title:'" + title + "', xaxis:{title:'cycles'}, yaxis:{title:'quality'}};\n"
        json_str += "Plotly.newPlot('" + div + "', data, layout);\n"
        return json_str

    def contentPlotly(self, div, title=""):
        colors = {'A':'rgba(255,0,0,0.5)', 'T':'rgba(128,0,128,0.5)', 'C':'rgba(0,255,0,0.5)', 'G':'rgba(0,0,255,0.5)'}
        json_str = "var data=["
        x = range(self.readLen)
        # four bases
        for base in ALL_BASES:
            json_str += "{"
            json_str += "x:[" + ",".join(map(str, x)) + "],"
            json_str += "y:[" + ",".join(map(str, self.percents[base][0:self.readLen])) + "],"
            json_str += "name: '" + base + "',"
            json_str += "mode:'lines',"
            json_str += "line:{color:'" + colors[base] + "', width:1}\n"
            json_str += "},"
        # mean
        json_str += "{"
        json_str += "x:[" + ",".join(map(str, x)) + "],"
        json_str += "y:[" + ",".join(map(str, self.gcPercents[0:self.readLen])) + "],"
        json_str += "name: 'GC',"
        json_str += "mode:'lines',"
        json_str += "line:{color:'rgba(20,20,20,255)', width:1}\n"
        json_str += "}\n"
        json_str += "];\n"
        json_str += "var layout={title:'" + title + "', xaxis:{title:'cycles'}, yaxis:{title:'percents', range:" + makeRange(0.0, 0.8) + "}};\n"
        json_str += "Plotly.newPlot('" + div + "', data, layout);\n"
        return json_str

    def gcPlotly(self, div, title=""):
        if self.readLen == 0:
            return ""
        json_str = "var data=["
        x = range(self.readLen+1)
        xticks = [100.0 * float(t)/self.readLen for t in x]
        json_str += "{"
        json_str += "x:[" + ",".join(map(str, xticks)) + "],"
        json_str += "y:[" + ",".join(map(str, self.gcHistogram[0:self.readLen+1])) + "],"
        json_str += "type:'bar'"
        json_str += "}];"
        json_str += "var layout={title:'" + title + "', xaxis:{title:'percents(%)'}, yaxis:{title:'counts'}};\n"
        json_str += "Plotly.newPlot('" + div + "', data, layout);\n"
        return json_str

    def discontinuityPlotly(self, div, title=""):
        json_str = "var data=["
        x = range(self.readLen)
        json_str += "{"
        json_str += "x:[" + ",".join(map(str, x)) + "],"
        json_str += "y:[" + ",".join(map(str, self.meanDiscontinuity[0:self.readLen])) + "],"
        json_str += "mode:'lines',"
        json_str += "line:{color:'rgba(100,150,0,0.5)', width:2}\n"
        json_str += "}];"
        json_str += "var layout={title:'" + title + "', xaxis:{title:'cycles'}, yaxis:{title:'discontinuity', range:" + makeRange(0.0, max(self.meanDiscontinuity)*1.5) + "}};\n"
        json_str += "Plotly.newPlot('" + div + "', data, layout);\n"
        return json_str

    def strandBiasPlotly(self, div, title=""):
        if self.readLen == 0:
            return ""
        shift = min(50, len(self.topKmerCount)/2)
        # we only sample 1000 points for performance issue
        top = min(len(self.topKmerCount) - shift,1000)
        forward = [0 for i in xrange(top)]
        reverse = [0 for i in xrange(top)]
        step = (len(self.topKmerCount) - shift) / top
        if step == 0:
            step = 1
        maxValue = 0
        for i in xrange(top):
            index = i*step+shift
            if index >= len(self.topKmerCount):
                break
            kmer = self.topKmerCount[i*step+shift][0]
            forward[i] = self.kmerCount[kmer]
            reverse[i] = self.kmerCount[util.reverseComplement(kmer)]
            maxValue = max(max(forward[i], reverse[i]), maxValue)
        json_str = "var data=["
        x = range(self.readLen)
        json_str += "{"
        json_str += "x:[" + ",".join(map(str, forward)) + "],"
        json_str += "y:[" + ",".join(map(str, reverse)) + "],"
        json_str += "mode:'markers',"
        json_str += "type:'scatter',\n"
        json_str += "marker:{size:2, color:'rgba(0,0,50,128)'}\n"
        json_str += "}];"
        json_str += "var layout={title:'" + title + "', xaxis:{title:'relative forward strand KMER count', range:" + makeRange(-10, maxValue) + "}, yaxis:{title:'relative reverse strand KMER count', range:" + makeRange(-10, maxValue) + "}};\n"
        json_str += "Plotly.newPlot('" + div + "', data, layout);\n"
        return json_str


    def overlapPlotly(self, overlap_histgram, readLen, total_reads, div):
        json_str = "var data=["
        x = range(self.readLen+1)
        json_str += "{"
        json_str += "x:[" + ",".join(map(str, x)) + "],"
        json_str += "y:[" + ",".join(map(str, overlap_histgram)) + "],"
        json_str += "type:'bar'"
        json_str += "}];"
        not_overlap_percent = 0
        if total_reads > 0:
            not_overlap_percent = int(overlap_histgram[0]*100.0/total_reads)
        xlabel = 'overlap Length (' + str(not_overlap_percent) + '% not overlapped)'
        json_str += "var layout={title:'Pair overlap Length Histgram', xaxis:{title:'" + xlabel + "', range:" + makeRange(-2, readLen) + "}, yaxis:{title:'counts'}};\n"
        json_str += "Plotly.newPlot('" + div + "', data, layout);\n"
        return json_str

    def errorPlotly(self, error_matrix, div):
        json_str = "var data=["
        names = []
        values = []
        colors = []
        for correct_base in ALL_BASES:
            for error_base in ALL_BASES:
                if correct_base != error_base:
                    name = correct_base + "->" + error_base
                    names.append(name)
                    values.append(error_matrix[correct_base][error_base])
                    if (correct_base=='A' and error_base=='G') or (correct_base=='G' and error_base=='A') or (correct_base=='C' and error_base=='T') or (correct_base=='T' and error_base=='C'):
                        colors.append("'rgba(246, 103, 0,1.0)'")
                    else:
                        colors.append("'rgba(22, 96, 167,1.0)'")
        json_str += "{"
        json_str += "x:['" + "','".join(names) + "'],"
        json_str += "y:[" + ",".join(map(str, values)) + "],"
        json_str += "marker:{color:[" + ",".join(colors) + "]},"
        json_str += "type:'bar'"
        json_str += "}];"
        json_str += "var layout={title:'sequencing error transform distribution', xaxis:{title:'seq error transform'}, yaxis:{title:'counts'}};\n"
        json_str += "Plotly.newPlot('" + div + "', data, layout);\n"
        return json_str

    def statPlotly(self, labels, counts, total_reads, div):
        json_str = "var data=["
        json_str += "{values:[" + ",".join(map(str, counts)) + "],"
        json_str += "labels:['" + "','".join(labels) + "'],"
        json_str += "textinfo: 'none',"
        json_str += "type:'pie'}];\n"
        title = "Filtering statistics of sampled " + str(total_reads) + " reads"
        json_str += "var layout={title:'" + title + "', width:800, height:600};\n"
        json_str += "Plotly.newPlot('" + div + "', data, layout);\n"
        return json_str

    def qc(self): 
        self.calcReadLen()
        self.calcPercents()
        self.calcQualities()
        self.calcDiscontinuity()
        self.sortKmer()
        
    def statFile(self, filename):
        READ_TO_SKIP = 1000
        reader = fastq.Reader(filename)
        stat_reads_num = 0
        skipped_reads = []
        #sample up to maxSample reads for stat
        while True:
            read = reader.nextRead()
            if read==None:
                break
            self.readCount += 1
            # here we skip the first 1000 reads because usually they are usually not stable
            if self.readCount < READ_TO_SKIP:
                skipped_reads.append(read)
                continue
            stat_reads_num += 1
            if stat_reads_num > self.sampleLimit and self.sampleLimit>0:
                break
            self.statRead(read)

        # if the fq file is too small, then we stat the skipped reads again
        if stat_reads_num < READ_TO_SKIP:
            for read in skipped_reads:
                self.statRead(read)

        self.qc()

    def autoTrim(self):
        #use (center-5, center+5) as initial good segment        
        center = int(self.readLen/2)
        front = center
        tail = center
        bad_in_front = False
        bad_in_tail = False

        for front in range(0, center)[::-1]:
            if self.isAbnormalCycle(front, front+1, 0.10):
                bad_in_front = True
                break

        for tail in range(center+1, self.readLen):
            if self.isAbnormalCycle(tail, tail-1, 0.05):
                bad_in_tail = True
                break

        trimFront = 0
        trimTail = 0
        if bad_in_front:
            trimFront = front+1
        if bad_in_tail:
            trimTail = self.readLen-tail
        
        trimFront = min(int(self.readLen*0.1),trimFront)
        trimTail = min(int(self.readLen*0.05),trimTail)
        
        return (trimFront, trimTail)

    def isAbnormalCycle(self, this_cycle, comp_cycle, percent_change_threshold):
        # thresholds
        BASE_TOP = 0.4
        BASE_BOTTOM = 0.15
        GC_TOP = 0.7
        GC_BOTTOM = 0.3
        QUAL_BOTTOM = 20.0

        if self.gcPercents[this_cycle] > GC_TOP or self.gcPercents[this_cycle] < GC_BOTTOM:
            return True

        for base in ALL_BASES:
            if self.percents[base][this_cycle] > BASE_TOP or self.percents[base][this_cycle] < BASE_BOTTOM:
                return True
            if abs(self.percents[base][this_cycle] - self.percents[base][comp_cycle]) > percent_change_threshold:
                return True
            if self.baseMeanQual[base][this_cycle] < QUAL_BOTTOM:
                return True

        return False

if __name__  == "__main__":
    qc = QualityControl()
    qc.statFile("R1.fq")
    qc.plot()
    print(qc.autoTrim())
    print(qc.topKmerCount[0:10])
