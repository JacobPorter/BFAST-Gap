#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <config.h>
#include <zlib.h>
#include <unistd.h>  

#include "../bfast/BLibDefinitions.h"
#include "../bfast/BError.h"
#include "../bfast/AlignedRead.h"
#include "../bfast/RGMatches.h"
#include "../bfast/RGMatch.h"
#include "bevalsim.h"

#define Name "bevalsim"
#define COUNT_ROTATE_NUM 10000

/* Parses a balign file resulting from usin reads generated by
 * bgeneratereads to give accuracy statistics for the mapping.
 * */

int PrintUsage()
{
	fprintf(stderr, "%s %s\n", "bfast", PACKAGE_VERSION);
	fprintf(stderr, "\nUsage:%s [options]\n", Name);
	fprintf(stderr, "\t-i\tFILE\tinput bfast matches file or bfast aligned file\n");
	fprintf(stderr, "\t-r\tFILE\treads file (FASTQ)\n"); 
	fprintf(stderr, "\t-h\t\tprints this help message\n");
	fprintf(stderr, "\nsend bugs to %s\n",
			PACKAGE_BUGREPORT);
	return 1;
}

int main(int argc, char *argv[]) 
{
	char *inputFileName=NULL;
	char *readsFile=NULL;
	int c, type;

	while((c = getopt(argc, argv, "i:r:h")) >= 0) {
		switch(c) {
			case 'h': return PrintUsage();
			case 'i': inputFileName=strdup(optarg); break;
			case 'r': readsFile=strdup(optarg); break;
			default: fprintf(stderr, "Unrecognized option: -%c\n", c); return 1;
		}
	}

	if(1 == argc || argc != optind) {
		return PrintUsage();
	}

	/* Check cmd line options */
	if(NULL == inputFileName) {
		PrintError(Name, "inputFileName", "Command line arguments", Exit, InputArguments);
	}
	if(NULL == readsFile) {
		PrintError(Name, "readsFile", "Command line arguments", Exit, InputArguments);
	}
	if(NULL!=strstr(inputFileName, BFAST_MATCHES_FILE_EXTENSION)) {
		type=BMF;
	}
	else if(NULL!=strstr(inputFileName, BFAST_ALIGNED_FILE_EXTENSION)) {
		type=BAF;
	}
	else {
		type=-1;
		PrintError(Name, "input file name", "Could not recognize file extension", Exit, OutOfRange);
	}

	/* Run program */
	Evaluate(inputFileName, 
			readsFile, 
			type
			);

	/* Terminate */
	fprintf(stderr, "%s", BREAK_LINE);
	fprintf(stderr, "Terminating successfully.\n");
	fprintf(stderr, "%s", BREAK_LINE);
	return 0;
}

void ReadTypeInitialize(ReadType *r)
{
	r->readNum=0;
	r->strand=0;
	r->contig=0;
	r->pos=0;
	r->numEnds=0;
	r->pairedEndLength=0;
	r->readLength=0;
	r->whichReadVariants=0;
	r->startIndel=0;
	r->indelLength=0;
	r->numSNPs=0;
	r->numErrors=0;
	r->deletionLength=0;
	r->insertionLength=0;
	r->aContigOne=NULL;
	r->aPosOne=NULL;
	r->aStrandOne=NULL;
	r->numOne=0;
	r->aContigTwo=NULL;
	r->aPosTwo=NULL;
	r->aStrandTwo=NULL;
	r->numTwo=0;
}

void ReadTypeCopy(ReadType *dest,
		ReadType *src)
{
	/* Only copy meta data */ 
	dest->readNum = src->readNum;
	dest->strand=src->strand;
	dest->contig=src->contig;
	dest->pos=src->pos;
	dest->numEnds=src->numEnds;
	dest->pairedEndLength=src->pairedEndLength;
	dest->readLength=src->readLength;
	dest->whichReadVariants=src->whichReadVariants;
	dest->startIndel=src->startIndel;
	dest->indelLength=src->indelLength;
	dest->numSNPs=src->numSNPs;
	dest->numErrors=src->numErrors;
	dest->deletionLength=src->deletionLength;
	dest->insertionLength=src->insertionLength;
}

void ReadTypePrint(ReadType *r, FILE *fp)
{
	fprintf(fp, "readNum=%d\n", r->readNum);
	fprintf(fp, "strand=%c\n", r->strand);
	fprintf(fp, "contig=%d\n", r->contig);
	fprintf(fp, "pos=%d\n", r->pos);
	fprintf(fp, "numEnds=%d\n", r->numEnds);
	fprintf(fp, "pairedEndLength=%d\n", r->pairedEndLength);
	fprintf(fp, "readLength=%d\n", r->readLength);
	fprintf(fp, "whichReadVariants=%d\n", r->whichReadVariants);
	fprintf(fp, "startIndel=%d\n", r->startIndel);
	fprintf(fp, "indelLength=%d\n", r->indelLength);
	fprintf(fp, "numSNPs=%d\n", r->numSNPs);
	fprintf(fp, "numErrors=%d\n", r->numErrors);
	fprintf(fp, "deletionLength=%d\n", r->deletionLength);
	fprintf(fp, "insertionLength=%d\n", r->insertionLength);
}

int ReadTypeCompare(ReadType *a,
		ReadType *b)
{
	/* Only compare meta data */ 
	/* Nice use of if, else if, and else statements */
	if(a->numEnds != b->numEnds) {
		return (a->numEnds < b->numEnds)?-1:1;
	}
	else if(a->pairedEndLength != b->pairedEndLength) {
		return (a->pairedEndLength < b->pairedEndLength)?-1:1;
	}
	else if(a->readLength != b->readLength) {
		return (a->readLength < b->readLength)?-1:1;
	}
	else if(a->indelLength != b->indelLength) {
		return (a->indelLength < b->indelLength)?-1:1;
	}
	else if(a->numSNPs != b->numSNPs) {
		return (a->numSNPs < b->numSNPs)?-1:1;
	}
	else if(a->numErrors != b->numErrors) {
		return (a->numErrors < b->numErrors)?-1:1;
	}
	else if(a->deletionLength != b->deletionLength) {
		return (a->deletionLength < b->deletionLength)?-1:1;
	}
	else if(a->insertionLength != b->insertionLength) {
		return (a->insertionLength < b->insertionLength)?-1:1;
	}
	else {
		return 0;
	}
}

int ReadTypeRead(ReadType *r, 
		gzFile fp,
		int type)
{
	char *FnName = "../bfast/ReadTypeRead";
	AlignedRead a;
	RGMatches m;
	int i;
	char *readName=NULL;

	/* Initialize */
	AlignedReadInitialize(&a);
	RGMatchesInitialize(&m);

	/* Read in align entries */
	if(BMF==type) {
		if(EOF==RGMatchesRead(fp, &m)) {

			return EOF;
		}
		readName = strdup(m.readName);
		/* Update the number of entries */
		r->numOne = m.ends[0].numEntries;
		r->numTwo = 0;
		if(1 < m.numEnds) {
			r->numTwo = m.ends[1].numEntries;
		}
	}
	else if(BAF==type) {
		if(EOF==AlignedReadRead(&a, fp)) {
			return EOF;
		}
		readName = malloc(sizeof(char)*(2 + a.readNameLength));
		if(NULL == readName) {
			PrintError(FnName, "readName", "Could not allocate memory", Exit, MallocMemory);
		}
		readName[0]='@'; readName[1]='\0';
		strcat(readName, a.readName);
		/* Update the number of entries */
		r->numOne = a.ends[0].numEntries;
		if(1 < a.numEnds) {
			r->numTwo = a.ends[1].numEntries;
		}
	}
	else {
		PrintError(FnName, "type", "Could not recongize file type", Exit, OutOfRange);
	}

	/* Allocate memory */
	if(r->numOne > 0) {
		r->aContigOne = malloc(sizeof(int)*r->numOne);
		if(NULL==r->aContigOne) {
			PrintError(FnName, "r->aContigOne", "Could not allocate memory", Exit, MallocMemory);
		}
		r->aPosOne = malloc(sizeof(int)*r->numOne);
		if(NULL==r->aPosOne) {
			PrintError(FnName, "r->aPosOne", "Could not allocate memory", Exit, MallocMemory);
		}
		r->aStrandOne= malloc(sizeof(char)*r->numOne);
		if(NULL==r->aStrandOne) {
			PrintError(FnName, "r->aStrandOne", "Could not allocate memory", Exit, MallocMemory);
		}
	}
	if(r->numTwo > 0) {
		r->aContigTwo = malloc(sizeof(int)*r->numTwo);
		if(NULL==r->aContigTwo) {
			PrintError(FnName, "r->aContigTwo", "Could not allocate memory", Exit, MallocMemory);
		}
		r->aPosTwo = malloc(sizeof(int)*r->numTwo);
		if(NULL==r->aPosTwo) {
			PrintError(FnName, "r->aPosTwo", "Could not allocate memory", Exit, MallocMemory);
		}
		r->aStrandTwo= malloc(sizeof(char)*r->numTwo);
		if(NULL==r->aStrandTwo) {
			PrintError(FnName, "r->aStrandTwo", "Could not allocate memory", Exit, MallocMemory);
		}
	}
	/* Copy over */
	if(BMF==type) {
		for(i=0;i<r->numOne;i++) {
			r->aContigOne[i] = m.ends[0].contigs[i];
			r->aPosOne[i] = m.ends[0].positions[i];
			r->aStrandOne[i] = m.ends[0].strands[i];
		}
		for(i=0;i<r->numTwo;i++) {
			r->aContigTwo[i] = m.ends[1].contigs[i];
			r->aPosTwo[i] = m.ends[1].positions[i];
			r->aStrandTwo[i] = m.ends[1].strands[i];
		}
	}
	else if(BAF==type) {
		for(i=0;i<r->numOne;i++) {
			r->aContigOne[i] = a.ends[0].entries[i].contig;
			r->aPosOne[i] = a.ends[0].entries[i].position;
			r->aStrandOne[i] = a.ends[0].entries[i].strand;
		}
		for(i=0;i<r->numTwo;i++) {
			r->aContigTwo[i] = a.ends[1].entries[i].contig;
			r->aPosTwo[i] = a.ends[1].entries[i].position;
			r->aStrandTwo[i] = a.ends[1].entries[i].strand;
		}
	}
	else {
		PrintError(FnName, "type", "Could not recongize file type", Exit, OutOfRange);
	}

	/* Convert into read type */
	ReadTypeParseReadName(r, readName);
	free(readName);
	readName=NULL;

	/* Delete align entries */
	AlignedReadFree(&a);

	return 1;
}

void ReadTypeParseReadName(ReadType *r, char *readName)
{
	char *FnName="../bfast/ReadTypeParseR1R2";
	char r1[SEQUENCE_LENGTH]="\0";
	char r2[SEQUENCE_LENGTH]="\0";
	char tempString[SEQUENCE_LENGTH]="\0";
	int i, j;
	char tempChar;
	int state;
	int numRead;

			
	numRead = sscanf(readName, 
				"@readNum=%d_strand=%c_contig=%d_pos=%d_numends=%d_pel=%d_rl=%d_wrv=%d_si=%d_il=%d_r1=%s",
				&r->readNum,
				&r->strand,
				&r->contig,
				&r->pos,
				&r->numEnds,
				&r->pairedEndLength,
				&r->readLength,
				&r->whichReadVariants,
				&r->startIndel,
				&r->indelLength,
				tempString);
	if(11 != numRead) {
		PrintError(FnName, readName, "Could not parse read name (0)", Exit, OutOfRange);
	}
	if(r->numEnds <= 0) {
		fprintf(stderr, "%s\n", readName); // HERE
		fprintf(stderr, "r->readNum=%d\n", r->readNum);
		fprintf(stderr, "r->strand=%c\n", r->strand);
		fprintf(stderr, "r->contig=%d\n", r->contig);
		fprintf(stderr, "r->pos=%d\n", r->pos);
		fprintf(stderr, "r->numEnds=%d\n", r->numEnds);
	}
	assert(0 < r->numEnds);
	if(1 == r->numEnds) {
		strcpy(r1, tempString);
	}
	else if(2 == r->numEnds) {
		/* Parse tempString */
		i=j=0;
		state = 0;
		while(EOF != sscanf(tempString+i, "%c", &tempChar)) {
			switch(tempChar) {
				case '0':
				case '1':
				case '2':
				case '3':
				case '4':
				case '5':
				case '6':
				case '7':
				case '8':
				case '9':
					/*
					   fprintf(stderr, "0-9=%c\n", tempChar);
					   */
					if(0 == state) {
						r1[j] = tempChar;
						j++;
					}
					else if(2 == state) {
						r2[j] = tempChar;
						j++;
					}
					break;
				case '=':
					/*
					   fprintf(stderr, "==%c\n", tempChar);
					   */
					state = 2;
					r1[j] = '\0';
					j=0;
					break;
				default:
					/*
					   fprintf(stderr, "tempChar=%c\n", tempChar);
					   */
					assert(0 <= state && state <= 1);
					state = 1;
					break;
			}
			i++;
		}
		r2[j] = '\0';
		/*
		   fprintf(stderr, "j=%d\n", j);
		   fprintf(stderr, "r->readLength=%d\nr1=%s\nr1(length)=%d\n",
		   r->readLength,
		   r1,
		   (int)strlen(r1));
		   fprintf(stderr, "r->readLength=%d\nr2=%s\nr2(length)=%d\n",
		   r->readLength,
		   r2,
		   (int)strlen(r2));
		   */
		assert(state == 2);
	}
	else {
		PrintError(FnName, "numEnds", "Found more than two ends for a read", Exit, OutOfRange);
	}

	/* Parse r1 and r2 */
	assert(r->readLength == (int)strlen(r1));
	assert(r->numEnds == 1 || r->readLength == (int)strlen(r2));
	r->numSNPs = 0;
	r->numErrors = 0;
	r->deletionLength = 0;
	r->insertionLength = 0;
	for(i=0;i<r->readLength;i++) {
		switch(r1[i]) {
			case '0':
				/* Default */
				break;
			case '1':
				/* Insertion */
				r->insertionLength++;
				break;
			case '2':
				/* SNP */
				r->numSNPs++;
				break;
			case '3':
				/* Error */
				r->numErrors++;
				break;
			case '4':
				/* InsertionAndSNP */
				r->insertionLength++;
				r->numSNPs++;
				break;
			case '5':
				/* InsertionAndError */
				r->insertionLength++;
				r->numErrors++;
				break;
			case '6':
				/* SNPAndError */
				r->numSNPs++;
				r->numErrors++;
				break;
			case '7':
				/* InsertionSNPAndError */
				r->insertionLength++;
				r->numSNPs++;
				r->numErrors++;
				break;
			default:
				PrintError(FnName, "r1[i]", "Could not understand type", Exit, OutOfRange);
		}
		if(2 == r->numEnds) {
			switch(r2[i]) {
				case '0':
					/* Default */
					break;
				case '1':
					/* Insertion */
					r->insertionLength++;
					break;
				case '2':
					/* SNP */
					r->numSNPs++;
					break;
				case '3':
					/* Error */
					r->numErrors++;
					break;
				case '4':
					/* InsertionAndSNP */
					r->insertionLength++;
					r->numSNPs++;
					break;
				case '5':
					/* InsertionAndError */
					r->insertionLength++;
					r->numErrors++;
					break;
				case '6':
					/* SNPAndError */
					r->numSNPs++;
					r->numErrors++;
					break;
				case '7':
					/* InsertionSNPAndError */
					r->insertionLength++;
					r->numSNPs++;
					r->numErrors++;
					break;
				default:
					PrintError(FnName, "r2[i]", "Could not understand type", Exit, OutOfRange);
			}
		}
	}
	if(r->startIndel >= 0 && r->insertionLength == 0) {
		r->deletionLength = r->indelLength;
	}
}

void ReadTypeDelete(ReadType *r)
{
	free(r->aContigOne);
	free(r->aPosOne);
	free(r->aStrandOne);
	free(r->aContigTwo);
	free(r->aPosTwo);
	free(r->aStrandTwo);
	ReadTypeInitialize(r);
}

void StatInitialize(Stat *s, 
		ReadType *r)
{
	s->numCorrectlyAligned[0]=0;
	s->numCorrectlyAligned[1]=0;
	s->numCorrectlyAligned[2]=0;
	s->numCorrectlyAligned[3]=0;
	s->numCorrectlyAligned[4]=0;
	ReadTypeCopy(&s->r, r);
	s->numAligned=0;
	s->numReads=0;
}

void StatPrint(Stat *s, FILE *fp)
{
	fprintf(fp, "%10d %10d %10d %10d %10d %10d %10d ",
			s->numReads,
			s->numAligned,
			s->numCorrectlyAligned[0],
			s->numCorrectlyAligned[1],
			s->numCorrectlyAligned[2],
			s->numCorrectlyAligned[3],
			s->numCorrectlyAligned[4]);
	fprintf(fp, "%d %6d %6d %6d %3d %3d %3d %3d\n",
			s->r.numEnds,
			s->r.pairedEndLength,
			s->r.readLength,
			s->r.indelLength,
			s->r.numSNPs,
			s->r.numErrors,
			s->r.deletionLength,
			s->r.insertionLength);
}

void StatAdd(Stat *s, ReadType *r, int readType)
{
	int i, j;
	int found[5]={0,0,0,0,0};
	int foundAll=0;

	assert(OriginalRead == readType || ReadAligned == readType);

	if(OriginalRead == readType) {
		s->numReads++;
	}
	else {
		if(2 == r->numEnds) {
			for(i=0;i<r->numOne && 0==foundAll;i++) {
				for(j=0;j<r->numTwo && 0==foundAll;j++) {
					foundAll = StatCompare(s, 
							r, 
							r->aContigOne[i],
							r->aPosOne[i],
							r->aStrandOne[i],
							r->aContigTwo[i],
							r->aPosTwo[i],
							r->aStrandTwo[i],
							found);
				}
			}
		}
		else {
			for(i=0;i<r->numOne && 0==foundAll;i++) {
				foundAll = StatCompare(s, 
						r, 
						r->aContigOne[i],
						r->aPosOne[i],
						r->aStrandOne[i],
						0,
						0,
						0,
						found);
			}
		}
		for(i=0;i<5;i++) {
			if(1==found[i]) {
				s->numCorrectlyAligned[i]++;
			}
		}
		if(r->numOne > 0 || (2 == r->numEnds && r->numTwo > 0)) {
			s->numAligned++;
		}
	}
	if(s->numReads < s->numAligned) {
		fprintf(stderr, "\ns->numReads=%d\ns->numAligned=%d\n",
				s->numReads,
				s->numAligned);
	}
	assert(s->numAligned <= s->numReads);
}

int StatCompare(Stat *s, 
		ReadType *r,
		int contigOne,
		int posOne,
		char strandOne,
		int contigTwo,
		int posTwo,
		char strandTwo,
		int *found)
{
	int diffOne, diffTwo;

	/* Must be on the same strand and contig */
	if(r->strand == strandOne && r->contig == contigOne &&
			(1 == r->numEnds || (r->strand == strandTwo && r->contig == contigTwo))) {
		diffOne = (r->pos > posOne)?(r->pos - posOne):(posOne - r->pos);
		diffTwo = 0;
		if(2 == r->numEnds) {
			diffTwo = (r->pos + r->readLength + r->pairedEndLength > posTwo)?(r->pos + r->readLength + r->pairedEndLength - posTwo):(posTwo - (r->pos + r->readLength + r->pairedEndLength)); 
		}

		/* Update */
		if(diffOne <= 10000 && diffTwo <= 10000) {
			found[4]=1;
			if(diffOne <= 1000 && diffTwo <= 1000) {
				found[3]=1;
				if(diffOne <= 100 && diffTwo <= 100) {
					found[2]=1;
					if(diffOne <= 10 && diffTwo <= 10) {
						found[1]=1;
						if(diffOne <= 0 && diffTwo <= 0) {
							found[0]=1;
							return 1;
						}
					}
				}
			}
		}
	}
	return 0;
}

void StatsInitialize(Stats *s) 
{
	s->stats=NULL;
	s->numStats=0;
}

void StatsPrintHeader(FILE *fp)
{
	fprintf(fp, "# COL | Description\n");
	fprintf(fp, "# 0   | number of reads\n");
	fprintf(fp, "# 1   | number of reads aligned\n");
	fprintf(fp, "# 2   | number of correctly aligned within 0 bases\n");
	fprintf(fp, "# 3   | number of correctly aligned within 10 bases\n");
	fprintf(fp, "# 4   | number of correctly aligned within 100 bases\n");
	fprintf(fp, "# 5   | number of correctly aligned within 1000 bases\n");
	fprintf(fp, "# 6   | number of correctly aligned within 10000 bases\n");
	fprintf(fp, "# 7   | paired end\n");
	fprintf(fp, "# 8   | paired end length\n");
	fprintf(fp, "# 9   | read length\n");
	fprintf(fp, "# 10  | indel length\n");
	fprintf(fp, "# 11  | number of snps\n");
	fprintf(fp, "# 12  | number of errors\n");
	fprintf(fp, "# 13  | deletion length\n");
	fprintf(fp, "# 14  | insertion length\n");
}

void StatsPrint(Stats *s, FILE *fp)
{
	int32_t i;
	StatsPrintHeader(fp);
	for(i=0;i<s->numStats;i++) {
		StatPrint(&s->stats[i], fp);
	}
}

void StatsAdd(Stats *s, ReadType *r, int readType)
{
	int32_t i;
	char *FnName="../bfast/StatsAdd";

	assert(OriginalRead == readType || ReadAligned == readType);

	/* Check if it fits somewhere */
	for(i=0;i<s->numStats;i++) {
		if(ReadTypeCompare(r, &s->stats[i].r)==0) {
			/* Add to current */
			StatAdd(&s->stats[i], r, readType);
			return; /* Get out of here */
		}
	}
	if(ReadAligned == readType) {
		ReadTypePrint(r, stderr);
		PrintError(FnName, NULL, "../bfast/Read type was not found in the original reads file", Exit, OutOfRange);
	}
	else {
		/* Otherwise start a new start entry */
		s->numStats++;
		s->stats = realloc(s->stats, sizeof(Stat)*s->numStats);
		if(NULL==s->stats) {
			PrintError(FnName, "s->stats", "Could not allocate memory", Exit, MallocMemory);
		}
		/* Initialize */
		StatInitialize(&s->stats[s->numStats-1], r);
		/* Add */
		StatAdd(&s->stats[s->numStats-1], r, readType);
	}
}

void StatsDelete(Stats *s)
{
	free(s->stats);
	s->stats=NULL;
	s->numStats=0;
}

void Evaluate(char *inputFileName,
		char *readsFile,
		int type)
{
	char *FnName="Evaluate";
	gzFile fpIn;
	FILE *fpOut;
	ReadType r;
	Stats s;
	int32_t count;

	StatsInitialize(&s);

	/* Get the number of reads for each class */
	ReadInReads(readsFile, &s);

	ReadTypeInitialize(&r);

	/* Open the inputFileName file */
	if(!(fpIn=gzopen(inputFileName, "rb"))) {
		PrintError(FnName, inputFileName, "Could not open file for reading", Exit, OpenFileError);
	}

	count = 0;
	fprintf(stderr, "../bfast/Reading in from %s.\nCurrently on:\n%d", inputFileName, 0);
	while(EOF != ReadTypeRead(&r, fpIn, type)) {
		count++;
		if(count % COUNT_ROTATE_NUM == 0) {
			fprintf(stderr, "\r%d", 
					count);
		}

		/* Process the read */
		StatsAdd(&s, &r, ReadAligned);

		/* Reinitialize */
		ReadTypeInitialize(&r);
	}
	fprintf(stderr, "\r%d\n", count);

	/* Open the output file */
	fprintf(stderr, "%s", BREAK_LINE);
	fprintf(stderr, "Outputting...\n");
	if(!(fpOut=fdopen(fileno(stdout), "wb"))) {
		PrintError(FnName, "stdout", "Could not open stdout for writing", Exit, WriteFileError);
	}

	/* Print Stats */
	StatsPrint(&s, fpOut);
	fprintf(stderr, "%s", BREAK_LINE);

	/* Delete Stats */
	StatsDelete(&s);

	/* Close the files */
	gzclose(fpIn);
	fclose(fpOut);
}

void ReadInReads(char *readsFile, Stats *s)
{
	char *FnName="../bfast/ReadInReads";
	FILE *fpIn=NULL;
	int count=0;
	ReadType r;
	char readName[2][SEQUENCE_NAME_LENGTH]={"\0", "\0"};
	char r1[SEQUENCE_LENGTH]="\0";
	char q[SEQUENCE_LENGTH]="\0";
	char r2[SEQUENCE_LENGTH]="\0";

	/* Open the reads file */
	if(!(fpIn=fopen(readsFile, "rb"))) {
		PrintError(FnName, readsFile, "Could not open file for reading", Exit, OpenFileError);
	}

	fprintf(stderr, "../bfast/Reading in original reads from %s.\nCurrently on:\n%d", readsFile, 0);
	count = 0;
	/* Read in read name and read(s) */
	while(EOF != fscanf(fpIn, "%s", readName[0]) &&
			EOF != fscanf(fpIn, "%s", r1) &&
			EOF != fscanf(fpIn, "%s", q) && // comment
			EOF != fscanf(fpIn, "%s", q) && // qual
			(NULL != strpbrk("numends=1", readName[0]) 
			 || (fscanf(fpIn, "%s", readName[1]) &&
				 EOF != fscanf(fpIn, "%s", r2) &&
				 EOF != fscanf(fpIn, "%s", q) && // comment
				 EOF != fscanf(fpIn, "%s", q)))) {
		/* Only accepting 1 or 2 ends currently */
		if(!(NULL != strpbrk("numends=1", readName[0]) ||
					NULL != strpbrk("numends=2", readName[0]))) {
			PrintError(FnName, readName[0], "Only single or paired end reads", Exit, OutOfRange);
		}
		assert(NULL != strpbrk("numends=1", readName[0]) ||
				NULL != strpbrk("numends=2", readName[0]));
		ReadTypeInitialize(&r);

		count++;
		if(count % COUNT_ROTATE_NUM == 0) {
			fprintf(stderr, "\r%d[%d]", 
					count,
					s->numStats
				   );
		}
		if(NULL != strpbrk("numends=2", readName[0]) &&
				0 == strcmp(readName[0], readName[1])) {
		}
		ReadTypeParseReadName(&r, readName[0]);

		/* Add to Stats */
		StatsAdd(s, &r, OriginalRead);

		r1[0]='\0';
		r2[0]='\0';
	}
	fprintf(stderr, "\r%d[%d]\n", 
			count,
			s->numStats
		   );

	/* Close the files */
	fclose(fpIn);
}
