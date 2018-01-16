# BFAST-Gap
BFAST-Gap is a modification of the BFAST short DNA read trimmer.  BFAST-Gap adds functionality for Ion Torrent DNA reads.

BFAST-Gap is backwards compatible with BFAST indexes, and it works similarily to BFAST.  Refer to the BFAST manual and documentation for general usage and installation and compilation instructions.


## Context Sensitive Gap Penalty Function

The one functional addition that BFAST-Gap makes is a context sensitive gap penalty function.  This addresses the overcalling and undercalling error of Ion Torrent read sequencing technology.  During the alignment phase of the program, a penatly for gap open and extension is assigned.  The context sensitive gap penalty function assigns such a penalty based on homopolymer run length of the read in the context where the gap is opened.  In order to use the context sensitive gap penalty function, a scoring file must be used with the '-x' option of the 'localalign' and 'postprocess' commands.

The scoring file consists of 14 numbers each on their own line.  The scoring file has the following format.

```
INT1 -- gap open penalty (negative)
INT2 -- gap extension penalty (negative)
INT3 -- match score (positive)
INT4 -- mismatch penalty (negative)
INT5 -- gap open function
FLOAT1 -- gap open function parameter 1
FLOAT2 -- gap open function parameter 2
FLOAT3 -- gap open function parameter 3
FLOAT4 -- gap open function parameter 4
INT6 -- gap extension function
FLOAT5 -- gap extension function parameter 1
FLOAT6 -- gap extension function parameter 2
FLOAT7 -- gap extension function parameter 3
FLOAT8 -- gap extension function parameter 4
```

There are four gap penalty functions that can be used, and these are given as integer values in the appropriate fields.  These functions consist of the following.

0 - constant (no context sensitive penalty, default BFAST functionality)
1 - piecewise linear
2 - exponential
3 - logistic

Each function has parameters that control the shape of the function.  If a parameter field is unused, it must be set to 0.0.  The four parameters are ordered and come after the integer field that determines the gap penalty function.  The parameters for each function are the following.

0 - constant (no parameters, all must be set to 0.0)
1 - piecewise linear (first parameter, third parameters are cutoffs for the run lenght, and the second and fourth parameters are scores)
2 - exponential -- (no parameters, all must be set to 0.0)
3 - logistic (first parameter -- slope, second parameter -- center, third and fourth parameters are unused and must be set to 0.0)

There is a maximum and a minimum value when using the context sensitive gap functions.  The maximum value of the gap open penalty is INT1, and the minimum value is the mismatch penalty, INT4.  The maximum value of the gap extension penalty is the gap extension penalty, INT2, and the minimum value is -1.0.

The implementation of the scoring functions is in the code file ScoringMatrix.c, and the code file AlignNTSpace.c calls the functions defined in ScoringMatrix.c.

## Example Scoring File

Example scoring files are found in the scoring directory.  The file 'logistic_open_1.0_-15_constant_extension.txt' is recommended as it performed well in tests.  This is that file:

-600
-50
96
-90
3
1.0
-15.0
0.0
0.0
0
0.0
0.0
0.0
0.0

This indicates a logistic gap open penalty with slope 1.0 and center -15.0.  A constant (that is, -50) gap extension penatly is used throughout.  This scoring file was found to perform well in a simulation.

An example of the piecewise function is given in the following.

2
10
-300
25
-150

This indicates that if the homopolymer run length is between 1 and 9 inclusive, then the score is the default score found in INT1 or INT2 depending on the function.  If the run length is between 10 and 24, then the score is -300, and if the run length is larger or equal to 25, then the score is -150.
