#ifndef SCORINGMATRIX_H_
#define SCORINGMATRIX_H_

#include "BLibDefinitions.h"

inline int32_t ScoringMatrixGetNTScore(char, char, ScoringMatrix*);
inline int32_t ScoringMatrixGetColorScore(char, char, ScoringMatrix*);

int ScoringMatrixRead(char*, ScoringMatrix*, int);
void ScoringMatrixInitialize(ScoringMatrix*);
int32_t ScoringMatrixCheck(ScoringMatrix*, int32_t);

/*JSP: functions for the context sensitive gap alignment */
void ScoringMatrixPrint(ScoringMatrix *sm, int space);
inline double ScoringMatrixGetGapOpenScore(int32_t runLength, ScoringMatrix* sm);
inline double ScoringMatrixGetGapExtensionScore(int32_t runLength, ScoringMatrix* sm );
inline double ScoringMatrixCalculateGapExtensionScore(int32_t runLength, ScoringMatrix* sm);
inline double ScoringMatrixCalculateGapOpenScore(int32_t runLength, ScoringMatrix* sm);

#endif
