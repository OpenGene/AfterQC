 #!/usr/bin/env python
  
import gzip
import os,sys

def isFastq(f):
    fqext = (".fq", ".fastq", "fq.gz", ".fastq.gz")
    for ext in fqext:
        if f.endswith(ext):
            return True
    return False

################################
#fastq.reader

class Reader:
    
    filename = ""
    
    __file = None
    __gz = False
    __eof = False
    
    def __init__(self, fname):
        self.filename = fname
        if self.filename.endswith(".gz"):
            self.__gz = True
            self.__file = gzip.open(self.filename, "r")
        else:
            self.__gz = False
            self.__file = open(self.filename, "r")
        if self.__file == None:
            print("Failed to open file " + self.filename)
            sys.exit(1)
            
    def __del__(self):
        if self.__file != None:
            self.__file.close()
            
    def nextRead(self):
        if self.__eof == True or self.__file == None:
            return None

        lines = []
        #read 4 (lines, name, sequence, strand, quality)
        for i in xrange(0,4):
            line = self.__file.readline()
            if len(line) == 0:
                self.__eof = True
                return None
            lines.append(line)
        return lines

    def isEOF(self):
        return False

################################
#fastq.writer

class Writer:
    
    filename = ""
    
    __file = None
    __gz = False
    
    def __init__(self, fname):
        self.filename = fname
        if self.filename.endswith(".gz"):
            self.__gz = True
            self.__file = gzip.open(self.filename, "w")
        else:
            self.__gz = False
            self.__file = open(self.filename, "w")
        if self.__file == None:
            print("Failed to open file " + self.filename + " to write")
            sys.exit(1)
            
    def __del__(self):
        if self.__file != None:
            self.__file.flush()
            self.__file.close()

    def flush(self):
        if self.__file !=None:
            self.__file.flush()
 
    def writeLines(self, lines):
        if self.__file == None:
            return False
            
        for line in lines:
            self.__file.write(line)
        return True
            
    def writeRead(self, name, seqence, strand, quality):
        if self.__file == None:
            return False
            
        self.__file.write(name)
        self.__file.write(seqence)
        self.__file.write(strand)
        self.__file.write(quality)
        
        return True
