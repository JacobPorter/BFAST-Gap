#ifndef BFASTPOSTPROCESS_H_
#define BFASTPOSTPROCESS_H_

#define MIN_FILTER 0
#define MAX_FILTER 4

/* This structure is used by main to communicate with parse_opt. */
struct arguments
{
	char *args[1];							/* No arguments to this function */
	char *fastaFileName;					/* -f */
	char *alignFileName;					/* -i */
	int algorithm;							/* -a */
	int space;								/* -A */
	int strandedness;						/* -S/-Y */
        int positioning;                                                /* -P/-Y */
        int pairing;                                          /* -Y */
	int avgMismatchQuality;					/* -q */
	char *scoringMatrixFileName;			/* -x */
	int randomBest;						    /* -z */
	int minMappingQuality;				    /* -m */
	int minNormalizedScore;				    /* -M */
	int insertSizeSpecified;
	double insertSizeAvg;						/* -v */
	double insertSizeStdDev;					/* -s */
	int numThreads;							/* -n */
	int queueLength;						/* -Q */
	int outputFormat;						/* -O */
	char *outputID;							/* -o */
	char *RGFileName;						/* -r */
	int baseQualityType;						/* -b */
	int timing;                             /* -t */
	int programMode;						/* -h */ 
};

/* Local functions */
int BfastPostProcessValidateInputs(struct arguments*);
void BfastPostProcessAssignDefaultValues(struct arguments*);
void BfastPostProcessPrintProgramParameters(FILE*, struct arguments*);
void BfastPostProcessFreeProgramParameters(struct arguments *args);
void BfastPostProcessPrintGetOptHelp();
void BfastPostProcessGetOptHelp();
void BfastPostProcessPrintGetOptHelp();
struct argp_option {
	char *name; /* Arg name */
	int key;
	char *arg; /* arg symbol */
	int flags;
	char *doc; /* short info about the arg */
	int group;
};
int BfastPostProcessGetOptParse(int, char**, char*, struct arguments*); 
#endif
