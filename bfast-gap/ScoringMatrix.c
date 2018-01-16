#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <math.h>

#include "BLib.h"
#include "BError.h"
#include "BLibDefinitions.h"
#include "ScoringMatrix.h"

/* TODO */
int ScoringMatrixRead(char *scoringMatrixFileName, 
		ScoringMatrix *sm,
		int space)
{
	char *FnName="ScoringMatrixRead";
	FILE *fp;

	/* Open the scoring matrix file */
	if((fp=fopen(scoringMatrixFileName, "r"))==0) {
		PrintError(FnName, scoringMatrixFileName, "Could not open scoringMatrixFileName for reading", Exit, OpenFileError);
	}

	/* Read in the gap open penalty,
	 * gap extension penalty,
	 * nt match score,
	 * nt mismatch score */
//	char com[21];
//	int counter = 0;
//	while(fscanf(fp,"%d",com) == 1)
//	{
//
//	    if(com[0]=='a')
//	        fprintf(dest,"%s 1",com);
//	     if(com[0]=='s')
//	        fprintf(dest,"%s 2",com);
//	}

	/* JSP: Modified to include the homopolymer gap function. */
	if(fscanf(fp, "%d %d %d %d %d %lf %lf %lf %lf %d %lf %lf %lf %lf", &sm->gapOpenPenalty,
				&sm->gapExtensionPenalty,
				&sm->ntMatch,
				&sm->ntMismatch,
				&sm->gapFunction,
				&sm->g1,
				&sm->g2,
				&sm->g3,
				&sm->g4,
				&sm->extFunction,
				&sm->e1,
				&sm->e2,
				&sm->e3,
				&sm->e4) == EOF) {
		PrintError(FnName, scoringMatrixFileName, "Could not read in the gap open penalty, gap extension penalty, nt match score, and nt mismatch score", Exit, OutOfRange);
	}

	if(space == 1) {
		if(fscanf(fp, "%d %d", &sm->colorMatch,
					&sm->colorMismatch) == EOF) {
			PrintError(FnName, scoringMatrixFileName, "Could not read in the color match score and color mismatch score", Exit, OutOfRange);
		}
	}
	//JSP: Initialize the lookup tables for the scoring function.
	int i;
	for (i = 1; i < TABLE_SIZE; i++) {
		sm->gapTable[i] = ScoringMatrixCalculateGapOpenScore(i, sm);
		sm->extTable[i] = ScoringMatrixCalculateGapExtensionScore(i, sm);
	}

	ScoringMatrixCheck(sm, space);

	/* Close the file */
	fclose(fp);
	return 1;
}

/* TODO */
void ScoringMatrixInitialize(ScoringMatrix *sm)
{
	sm->gapOpenPenalty=SCORING_MATRIX_GAP_OPEN;
	sm->gapExtensionPenalty=SCORING_MATRIX_GAP_EXTEND;
	sm->ntMatch=SCORING_MATRIX_NT_MATCH;
	sm->ntMismatch=SCORING_MATRIX_NT_MISMATCH;
	sm->colorMatch=SCORING_MATRIX_COLOR_MATCH;
	sm->colorMismatch=SCORING_MATRIX_COLOR_MISMATCH;
	/*JSP: initialize default function */
	sm->gapFunction=SCORING_FUNCTION;
	sm->g1=G1;
	sm->g2=G2;
	sm->g3=G3;
	sm->g4=G4;
	sm->extFunction=EXT_FUNCTION;
	sm->e1=E1;
	sm->e2=E2;
	sm->e3=E3;
	sm->e4=E4;
	int i;
	for (i = 0; i < TABLE_SIZE; i++) {
		sm->gapTable[i] = SCORING_MATRIX_GAP_OPEN;
		sm->extTable[i] = SCORING_MATRIX_GAP_EXTEND;
	}
}

/* TODO */
/* For color space only */
int32_t ScoringMatrixCheck(ScoringMatrix *sm,
		int space) {
	char *FnName="ScoringMatrixCheck";

	if(0 < sm->gapOpenPenalty) {
		PrintError(FnName, "sm->gapOpenPenalty", "Must be less than or equal to zero", Exit, OutOfRange);
	}
	if(0 < sm->gapExtensionPenalty) {
		PrintError(FnName, "sm->gapExtensionPenalty", "Must be less than or equal to zero", Exit, OutOfRange);
	}
	if(sm->gapExtensionPenalty < sm->gapOpenPenalty) {
		PrintError(FnName, "sm->gapExtensionPenalty < sm->gapOpenPenalty", "Gap extend must be greater than gap open", Exit, OutOfRange);
	}

	if(sm->gapExtensionPenalty < sm->ntMismatch) {
		PrintError(FnName, "sm->gapExtensionPenalty < sm->ntMismatch", "Gap extend must be greater than mismatch", Exit, OutOfRange);
	}

	if(sm->ntMismatch <= sm->gapOpenPenalty) {
		PrintError(FnName, "sm->ntMismatch <= sm->gapOpenPenalty", "Mismatch must be greater than one-base gap", Exit, OutOfRange);
	}

	if(sm->ntMatch < 0) {
		PrintError(FnName, "sm->ntMatch", "Must be greater than or equal to zero", Exit, OutOfRange);
	}
	if(ColorSpace == space && sm->colorMatch < 0) {
		PrintError(FnName, "sm->colorMatch", "Must be greater than or equal to zero", Exit, OutOfRange);
	}
	if(0 < sm->ntMismatch) {
		PrintError(FnName, "sm->ntMismatch", "Must be less than or equal to zero", Exit, OutOfRange);
	}
	if(ColorSpace == space && 0 < sm->colorMismatch) {
		PrintError(FnName, "sm->colorMismatch", "Must be less than or equal to zero", Exit, OutOfRange);
	}
	return 1;
}

inline int32_t ScoringMatrixGetNTScore(char a,
		char b,
		ScoringMatrix *sm)
{
	return (ToUpper(a) == ToUpper(b)) ? sm->ntMatch : sm->ntMismatch;
}

inline int32_t ScoringMatrixGetColorScore(char a, 
		char b, 
		ScoringMatrix *sm) 
{
	return (a == b) ? sm->colorMatch : sm->colorMismatch;
}

/*JSP: Below this comment are functions for the context sensitive gap open strategy.*/

/*
 * Retrieves the gap open score from either the lookup table in ScoringMatrix sm
 * or by calculating the result if the value is not found in the lookup table.
 */
inline double ScoringMatrixGetGapOpenScore(int32_t runLength, ScoringMatrix* sm) {
	if (runLength < TABLE_SIZE && runLength > 0) {
		return sm->gapTable[runLength];
	} else if (runLength > 0) {
		return ScoringMatrixCalculateGapOpenScore(runLength, sm);
	} else {
		return sm->gapOpenPenalty;
	}
}

/*
 * Calculates the gap open score based on the function specified in the ScoringMatrix sm.
 */
inline double ScoringMatrixCalculateGapOpenScore(int32_t runLength, ScoringMatrix* sm) {
	double return_value = (double) sm->gapOpenPenalty;
	if (sm->gapFunction == 1 && (double) runLength >= sm->g1 && (runLength < sm->g3 || sm->g3 <= 0.0 )) { //Piecewise function
		if (sm->g2 >= 0.0) {
			double value = ((double) sm->gapOpenPenalty - (double) sm->ntMismatch) / (double) 2;
			return_value = (double) sm->ntMismatch + value;
		} else {
			return_value = sm->g2;
		}
	} else if (sm->gapFunction == 1 && sm->g3 > 0 && (double) runLength >= sm->g3) { // Piecewise function
		if (sm->g4 >= 0.0) {
			double value = ((double) sm->gapOpenPenalty - (double) sm->ntMismatch) / (double) 2;
			return_value = (double) sm->ntMismatch + value;
		} else {
			return_value = sm->g4;
		}
	} else if (sm->gapFunction == 2) { // Exponential function
		return_value = ((double) sm->gapOpenPenalty - (double) sm->ntMismatch) * pow(M_E, 1 - runLength) + (double) sm->ntMismatch;
	} else if (sm->gapFunction == 3) { // Logistic function
		double numerator = -((double) sm->gapOpenPenalty - (double) sm->ntMismatch);
		double denominator = 1 + pow(M_E, sm->g1 * (runLength - sm->g2));
		return_value = -(numerator / denominator - (double) sm->ntMismatch);
	} else { //Default
		return_value = (double) sm->gapOpenPenalty;
	}
//fprintf(stderr, "GapOpenScore: %d %d %f\n", sm->gapFunction, runLength, return_value);
	return return_value;
}

/*
 * Retrieves the gap extension score from either the lookup table in ScoringMatrix sm
 * or by calculating the result if the value is not found in the lookup table.
 */
inline double ScoringMatrixGetGapExtensionScore(int32_t runLength, ScoringMatrix* sm) {
	if (runLength < TABLE_SIZE && runLength > 0) {
		return sm->extTable[runLength];
	} else if (runLength > 0) {
		return ScoringMatrixCalculateGapExtensionScore(runLength, sm);
	} else {
		return sm->gapExtensionPenalty;
	}
}

/*
 * Calculates the gap extension score based on the function specified in the ScoringMatrix sm.
 */
inline double ScoringMatrixCalculateGapExtensionScore(int32_t runLength, ScoringMatrix* sm) {
	double maximum_value = -1.0;
	double return_value = (double) sm->gapExtensionPenalty;
	if (sm->extFunction == 1 && (double) runLength >= sm->e1 && (runLength < sm->e3 || sm->e3 <= 0.0 )) { //Piecewise function
		if (sm->e2 >= 0.0) {
			double value = ((double) sm->gapExtensionPenalty - (double) maximum_value) / (double) 2;
			return_value =  (double) maximum_value + value;
		} else {
			return_value =  sm->e2;
		}
	} else if (sm->extFunction == 1 && sm->e3 > 0 && (double) runLength >= sm->e3) { // Piecewise function
		if (sm->e4 >= 0.0) {
			double value = ((double) sm->gapExtensionPenalty - (double) maximum_value) / (double) 2;
			return_value =  (double) maximum_value + value;
		} else {
			return_value =  sm->e4;
		}
	} else if (sm->extFunction == 2) { // Exponential function
		return_value =  ((double) sm->gapExtensionPenalty - (double) maximum_value) * pow(M_E, 1 - runLength) + (double) maximum_value;
	} else if (sm->extFunction == 3) { // Logistic function
		double numerator = -((double) sm->gapExtensionPenalty - (double) maximum_value);
		double denominator = 1 + pow(M_E, sm->e1 * (runLength - sm->e2));
		return_value =  -(numerator / denominator - (double) maximum_value);
	} else { //Default
		return_value =  (double) sm->gapExtensionPenalty;
	}
	//fprintf(stderr, "GapExtScore: %d %d %f\n", sm->extFunction, runLength, return_value);
	return return_value;
}


/*
 * Prints out a representation of the scoring matrix.
 */
void ScoringMatrixPrint(ScoringMatrix *sm, int space) {
	if (space != 1) {
		fprintf(stderr, "The following sequence of numbers will be used for the alignment scoring function:\t%d %d %d %d %d %lf %lf %lf %lf %d %lf %lf %lf %lf\n", sm->gapOpenPenalty, sm->gapExtensionPenalty, sm->ntMatch, sm->ntMismatch, sm->gapFunction, sm->g1, sm->g2, sm->g3, sm->g4, sm->extFunction, sm->e1, sm->e2, sm->e3, sm->e4);
	}
}
