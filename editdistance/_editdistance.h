#ifndef ___EDITDISTANCE__H__
#define ___EDITDISTANCE__H__

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// struct PatternMap {
//     uint64_t p_[256][4];
//     unsigned int tmax_;
//     unsigned int tlen_;
// };

unsigned int edit_distance(const char *a, const unsigned int asize, const char *b, const unsigned int bsize);
// void create_patternmap(struct PatternMap *pm, const char *a, const unsigned int size);
// unsigned int edit_distance_by_patternmap(struct PatternMap *mp, const char *b, const unsigned int size);


// return the best offset to align two sequence with hamming distance less than limit_distance
// if distance is found >= limit_distance in <=complete_compare_require first substring, return a flag value 0x7FFFFFFF
int seek_overlap(const char *a, const int asize, const char *b, const int bsize, const int limit_distance, const int complete_compare_require, const int overlap_require);
// void create_patternmap(struct PatternMap *pm, const char *a, const unsigned int size);
// unsigned int edit_distance_by_patternmap(struct PatternMap *mp, const char *b, const unsigned int size);

#ifdef __cplusplus
}
#endif

#endif
