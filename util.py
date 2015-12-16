 #!/usr/bin/env python
 
import os,sys
    
def parseBool(str):
    str = str.lower()
    if str=="true" or str=="yes" or str=="on":
        return True
    else:
        return False

def reverseComplement(origin):
    comp = {"A" : "T", "T" : "A", "C" : "G", "G" : "C", "a" : "t", "t" : "a", "c" : "g", "g" : "c", "N":"N", "\n":"\n"}
    length = len(origin)
    revCompArr = ['' for x in xrange(length)]
    for i in xrange(length):
        orig = origin[length - i -1]
        if orig in comp:
            revCompArr[i] = comp[orig]
        else:
            revCompArr[i] = 'N'
    return ''.join(revCompArr)

#simple edit distance
def editDistance(s1, s2):
    m=len(s1)+1
    n=len(s2)+1

    tbl = {}
    for i in xrange(m): tbl[i,0]=i
    for j in xrange(n): tbl[0,j]=j
    for i in xrange(1, m):
        for j in xrange(1, n):
            cost = 0 if s1[i-1] == s2[j-1] else 1
            tbl[i,j] = min(tbl[i, j-1]+1, tbl[i-1, j]+1, tbl[i-1, j-1]+cost)

    return tbl[i,j]
