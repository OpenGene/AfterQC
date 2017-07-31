 #!/usr/bin/env python
 
import os,sys
from ctypes import *

EDIT_DISTANCE_MODULE_EXISTS = True
EDIT_DISTANCE_CTYPES_LOADED = False

# try load editdistance built-in module first
# if editdistance built-in module not exist, try to load editdistance/libed.so, which is made by make in AfterQC folder
try:
    import editdistance
except ImportError:
    EDIT_DISTANCE_MODULE_EXISTS = False
    ed_lib_file  = os.path.join(sys.path[0], "editdistance/libed.so")
    if os.path.exists(ed_lib_file):
        try:
            ed_ctypes = cdll.LoadLibrary(ed_lib_file)
        except Exception:
            EDIT_DISTANCE_CTYPES_LOADED = False
        else:
            EDIT_DISTANCE_CTYPES_LOADED = True
else:
    EDIT_DISTANCE_MODULE_EXISTS = True


COMP = {"A" : "T", "T" : "A", "C" : "G", "G" : "C", "a" : "t", "t" : "a", "c" : "g", "g" : "c", "N":"N", "\n":"\n"}
    
def parseBool(str):
    str = str.lower()
    if str=="true" or str=="yes" or str=="on":
        return True
    else:
        return False

def complement(base):
    return COMP[base]

def qualNum(q):
    return ord(q) - 33

def reverseComplement(origin):
    length = len(origin)
    revCompArr = ['' for x in xrange(length)]
    for i in xrange(length):
        orig = origin[length - i -1]
        if orig in COMP:
            revCompArr[i] = COMP[orig]
        else:
            revCompArr[i] = 'N'
    return ''.join(revCompArr)

def reverse(origin):
    return origin[::-1]

def hammingDistance(s1, s2):
    length = min(len(s1), len(s2))
    d = 0
    for i in xrange(length):
        if s1[i] != s2[i]:
            d += 1
    return d

#simple edit distance
def editDistance(s1, s2):
    # check if editdistance module loaded
    if EDIT_DISTANCE_MODULE_EXISTS:
        return editdistance.eval(s1, s2)
    elif EDIT_DISTANCE_CTYPES_LOADED:
        return ed_ctypes.edit_distance(s1, len(s1), s2, len(s2))

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
    return min(3, overlap_len/10.0)

def overlap(r1, r2):
    return overlap_hm(r1, r2)

def overlap_ed(r1, r2):
    len1 = len(r1)
    len2 = len(r2)
    reverse_r2 = reverseComplement(r2)

    overlapped = False
    overlap_len = 0
    offset = 0
    distance = 0
    offset_0_is_min = True
    # a match of less than 10 is considered as unconfident
    while offset < len1-10 and overlapped==False:
        # the overlap length of r1 & r2 when r2 is move right for offset
        overlap_len = min(len1-offset, len2)

        # remind that Julia is a 1-based coordination system
        distance = editDistance(r1[offset : offset+overlap_len], reverse_r2[0 : overlap_len])
        threshold = distance_threshold(overlap_len)
        if distance <= threshold:
            # now we find a good candidate
            # we verify it by moving r2 one more base to see if the distance is getting longer
            # if yes, then current is the best match, otherwise it's not
            while offset < len1-10:
                next_offset = offset + 1
                next_overlap_len = min(len1-next_offset, len2)
                next_distance = editDistance(r1[next_offset : next_offset+next_overlap_len], reverse_r2[0 : next_overlap_len])
                if distance <= next_distance:
                    overlapped = True
                    break
                else:
                    offset_0_is_min = False
                    offset = next_offset
                    distance = next_distance
                    overlap_len = next_overlap_len
            break
        else:
            offset += max(1, (distance - int(threshold))/2 )
    if offset_0_is_min:
        # in this case, the adapter is sequenced since TEMPLATE_LEN < SEQ_LEN
        # check if distance can get smaller if offset goes negative
        # this only happens when insert DNA is shorter than sequencing read length, and some adapter/primer is sequenced but not trimmed cleanly
        # we go reversely
        offset = 0
        while offset > -(len2-10):
            # the overlap length of r1 & r2 when r2 is move right for offset
            overlap_len = min(len1,  len2- abs(offset))
            distance = editDistance(r1[0:overlap_len], reverse_r2[-offset : -offset + overlap_len])
            threshold = distance_threshold(overlap_len)
            if distance <= threshold:
                while offset > -(len2-10):
                    next_offset = offset - 1
                    next_overlap_len = min(len1,  len2- abs(next_offset))
                    next_distance = editDistance(r1[0:next_overlap_len], reverse_r2[-next_offset : -next_offset + next_overlap_len])
                    if distance <= next_distance:
                        return (offset, overlap_len, distance)
                    else:
                        distance = next_distance
                        overlap_len = next_overlap_len
                        offset = next_offset
            else:
                offset -= max(1, (distance - int(threshold))/2 )
    elif overlapped:
        return (offset, overlap_len, distance)

    return (0,0,0)

# calculate overlap by hamming distance
def overlap_hm(r1, r2):
    len1 = len(r1)
    len2 = len(r2)
    reverse_r2 = reverseComplement(r2)

    limit_distance = 3
    overlap_require = 30
    complete_compare_require = 50

    overlap_len = 0
    offset = 0

    # forward
    # a match of less than 10 is considered as unconfident
    while offset < len1-overlap_require:
        # the overlap length of r1 & r2 when r2 is move right for offset
        overlap_len = min(len1-offset, len2)

        diff = 0
        for i in xrange(overlap_len):
            if r1[offset + i] != reverse_r2[i]:
                diff += 1
                if diff >= limit_distance and i < complete_compare_require:
                    break
        
        if diff < limit_distance or (diff >= limit_distance and i>complete_compare_require):
            return (offset, overlap_len, diff)

        offset += 1


    # reverse
    # in this case, the adapter is sequenced since TEMPLATE_LEN < SEQ_LEN
    # check if distance can get smaller if offset goes negative
    # this only happens when insert DNA is shorter than sequencing read length, and some adapter/primer is sequenced but not trimmed cleanly
    # we go reversely
    offset = 0
    while offset > -(len2-overlap_require):
        # the overlap length of r1 & r2 when r2 is move right for offset
        overlap_len = min(len1,  len2- abs(offset))

        diff = 0
        for i in xrange(overlap_len):
            if r1[i] != reverse_r2[-offset + i]:
                diff += 1
                if diff >= limit_distance and i < complete_compare_require:
                    break
        
        if diff < limit_distance or (diff >= limit_distance and i>complete_compare_require):
            return (offset, overlap_len, diff)

        offset -= 1

    # not matched
    return (0,0,0)

def overlap_hm_cpp(r1, r2):
    len1 = len(r1)
    len2 = len(r2)
    reverse_r2 = reverseComplement(r2)

    limit_distance = 3
    overlap_require = 30
    complete_compare_require = 50

    ret =  ed_ctypes.seek_overlap(r1, len1, reverse_r2, len2, limit_distance, overlap_require, complete_compare_require)
    offset = ret >> 8
    diff = ret & (0xFF)
    if ret == 0x7FFFFFFF:
        return (0,0,0)
    elif offset >=0:
        overlap_len = min(len1-offset, len2)
    else:
        overlap_len = min(len1,  len2- abs(offset))

    return (offset, overlap_len, diff)


def changeString(str, pos, val):
    lst = list(str)
    lst[pos] = val
    return ''.join(lst)

if __name__  == "__main__":
    r1 = "CAGCGCCTACGGGCCCCTTTTTCTGCGCGACCGCGTGGCTGTGGGCGCGGATGCCTTTGAGCGCGGTGACTTCTCACTGCGTATCGAGCCGCTGGAGGTCTCCC"
    r2 = "ACCTCCAGCGGCTCGATACGCAGTGAGAAGTCACCGCGCTCAAAGGCATCCGCGCCCACAGCCACGCGGTCGCGCAGAAAAAGGGGCCCGTAGGCGCGGCTCCC"

    r1 = "CAGCGCCTACGGGCCCCTTTTTCTGCGCGACCGCGTGGCTGTGGGCGCGGATGCCTTTGAGCGCGGTGACTTCTCACTGCGTATCGAGC"
    r2 = "ACCTCCAGCGGCTCGATACGCAGTGAGAAGTCACCGCGCTCAAAGGCATCCGCGCCCACAGCCACGCGGTCGCGCAGAAAAAGGGGTCC"
    print(overlap_ed(r1, r2))
    print(overlap_hm(r1, r2))
    print(overlap_hm_cpp(r1, r2))
    for i in xrange(10000):
        overlap_ed(r1, r2)