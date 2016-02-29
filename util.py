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

    tbl = [([0] * n) for i in xrange(m)]
    for i in xrange(m):tbl[i][0]=i
    for j in xrange(n):tbl[0][j]=j
    for i in xrange(1, m):
        for j in xrange(1, n):
            cost = 0 if s1[i-1] == s2[j-1] else 1
            tbl[i][j] = min(tbl[i][j-1]+1, tbl[i-1][j]+1, tbl[i-1][j-1]+cost)

    return tbl[i][j]

def distance_threshold(overlap_len):
    return min(6, overlap_len/10.0)

def overlap(r1, r2):
    len1 = len(r1)
    len2 = len(r2)
    reverse_r2 = reverseComplement(r2)

    # a match of less than 5 is considered as unconfident
    for offset in range(0,len1-5):
        # the overlap length of r1 & r2 when r2 is move right for offset
        overlap_len = min(len1, (offset + len2 )) - offset

        # remind that Julia is a 1-based coordination system
        distance = editDistance(r1[offset : offset+overlap_len], reverse_r2[0 : overlap_len])
        if distance <= distance_threshold(overlap_len):
            # now we find a good candidate
            # we verify it by moving r2 one more base to see if the distance is getting longer
            # if yes, then current is the best match, otherwise it's not
            next = offset + 1
            next_overlap_len = min(len1, (next + len2 )) - next
            distance2 = editDistance(r1[next : next+next_overlap_len], reverse_r2[0 : next_overlap_len])
            if distance <= distance2:
                return (offset, overlap_len, distance)
    return (0,0,0)
