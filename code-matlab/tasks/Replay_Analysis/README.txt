REPLAY ANALYSIS FOLDER

Includes:

example workflow script
- ReplayEventDetection

processing functions
- getNAU
- getSWR
- MergeIV
- getCandSeq
- scoreCandSeq

extra utility
- subaxis

----------------------------------------------------------------------
VERSION NOTES:

youkitan 2015-01-03
This is a first attempt at automation of replay/preplay detection. The idea is to have an
easily manipulatable workflow so that various sequence detection and statistical tests can be
run on the same data seamlessly. Currently, only simple sequence detection (multi-unit activity,
SWR occurance etc.) and a rank-order correlation statistic can be run. 

--------------------------------------------------------------------------
ISSUES/AUTHOR NOTES:

- Run the workflow script!!!

- getNAU is not the same as the Dragoi & Tonegawa (2011) method of counting active cells per timebin.
getNAU as a function returns a tsd with NAU. Calling getNAU in getCandSeq restricts the output iv to
those which have at least some NAU. 

- MergeIV is a placeholder for the new InstersectIV. The current implementation of IntersectIV does not
take a set intersect since it keeps the entire interval of the first input iv. In a later update I will
move the current IntersectIV functionality into SelectIV. Merge can do two separate functions: a set
intersection and merging overlapping ivs (specified by input flags) 

- getCandSeq is a bit awkward and might not need to be a function when the helper functions are changed
to work better together independently