Flat profile:

Each sample counts as 0.01 seconds.
  %   cumulative   self              self     total           
 time   seconds   seconds    calls  ms/call  ms/call  name    
 59.53   3531.70  3531.70 19832594     0.18     0.18  PairExactLogLikelihood
 22.31   4855.59  1323.89                             ran2
 17.16   5873.62  1018.03   355000     2.87     2.87  Gaussian
  0.32   5892.55    18.93    74000     0.26    48.23  SimplexMaxLikelihood
  0.25   5907.49    14.94 11660120     0.00     0.00  SimplexSort
  0.21   5919.71    12.22   142000     0.09     0.09  CrossCorr
  0.16   5929.27     9.56    71000     0.13     8.74  PairNGSBData
  0.06   5932.75     3.48 19832594     0.00     0.00  BoundaryCheck
  0.00   5932.85     0.10       37     2.70 111108.72  NonGaussianDistance
  0.00   5932.89     0.04       34     1.18 14646.39  GaussianDistance
  0.00   5932.89     0.00        1     0.00     0.00  SeedRand

 %         the percentage of the total running time of the
time       program used by this function.

cumulative a running sum of the number of seconds accounted
 seconds   for by this function and those listed above it.

 self      the number of seconds accounted for by this
seconds    function alone.  This is the major sort for this
           listing.

calls      the number of times this function was invoked, if
           this function is profiled, else blank.
 
 self      the average number of milliseconds spent in this
ms/call    function per call, if this function is profiled,
	   else blank.

 total     the average number of milliseconds spent in this
ms/call    function and its descendents per call, if this 
	   function is profiled, else blank.

name       the name of the function.  This is the minor sort
           for this listing. The index shows the location of
	   the function in the gprof listing. If the index is
	   in parenthesis it shows where it would appear in
	   the gprof listing if it were to be printed.

		     Call graph (explanation follows)


granularity: each sample hit covers 4 byte(s) for 0.00% of 5932.89 seconds

index % time    self  children    called     name
                                                 <spontaneous>
[1]     77.7    0.00 4609.00                 main [1]
                0.10 4110.92      37/37          NonGaussianDistance [2]
                0.04  497.94      34/34          GaussianDistance [8]
                0.00    0.00       1/1           SeedRand [12]
-----------------------------------------------
                0.10 4110.92      37/37          main [1]
[2]     69.3    0.10 4110.92      37         NonGaussianDistance [2]
               18.93 3550.12   74000/74000       SimplexMaxLikelihood [3]
                4.98  318.31   37000/71000       PairNGSBData [7]
              212.21    0.00   74000/355000      Gaussian [6]
                6.37    0.00   74000/142000      CrossCorr [10]
-----------------------------------------------
               18.93 3550.12   74000/74000       NonGaussianDistance [2]
[3]     60.2   18.93 3550.12   74000         SimplexMaxLikelihood [3]
             3531.70    0.00 19832594/19832594     PairExactLogLikelihood [4]
               14.94    0.00 11660120/11660120     SimplexSort [9]
                3.48    0.00 19832594/19832594     BoundaryCheck [11]
-----------------------------------------------
             3531.70    0.00 19832594/19832594     SimplexMaxLikelihood [3]
[4]     59.5 3531.70    0.00 19832594         PairExactLogLikelihood [4]
-----------------------------------------------
                                                 <spontaneous>
[5]     22.3 1323.89    0.00                 ran2 [5]
-----------------------------------------------
              195.00    0.00   68000/355000      GaussianDistance [8]
              212.21    0.00   74000/355000      NonGaussianDistance [2]
              610.82    0.00  213000/355000      PairNGSBData [7]
[6]     17.2 1018.03    0.00  355000         Gaussian [6]
-----------------------------------------------
                4.58  292.50   34000/71000       GaussianDistance [8]
                4.98  318.31   37000/71000       NonGaussianDistance [2]
[7]     10.5    9.56  610.82   71000         PairNGSBData [7]
              610.82    0.00  213000/355000      Gaussian [6]
-----------------------------------------------
                0.04  497.94      34/34          main [1]
[8]      8.4    0.04  497.94      34         GaussianDistance [8]
                4.58  292.50   34000/71000       PairNGSBData [7]
              195.00    0.00   68000/355000      Gaussian [6]
                5.85    0.00   68000/142000      CrossCorr [10]
-----------------------------------------------
               14.94    0.00 11660120/11660120     SimplexMaxLikelihood [3]
[9]      0.3   14.94    0.00 11660120         SimplexSort [9]
-----------------------------------------------
                5.85    0.00   68000/142000      GaussianDistance [8]
                6.37    0.00   74000/142000      NonGaussianDistance [2]
[10]     0.2   12.22    0.00  142000         CrossCorr [10]
-----------------------------------------------
                3.48    0.00 19832594/19832594     SimplexMaxLikelihood [3]
[11]     0.1    3.48    0.00 19832594         BoundaryCheck [11]
-----------------------------------------------
                0.00    0.00       1/1           main [1]
[12]     0.0    0.00    0.00       1         SeedRand [12]
-----------------------------------------------

 This table describes the call tree of the program, and was sorted by
 the total amount of time spent in each function and its children.

 Each entry in this table consists of several lines.  The line with the
 index number at the left hand margin lists the current function.
 The lines above it list the functions that called this function,
 and the lines below it list the functions this one called.
 This line lists:
     index	A unique number given to each element of the table.
		Index numbers are sorted numerically.
		The index number is printed next to every function name so
		it is easier to look up where the function in the table.

     % time	This is the percentage of the `total' time that was spent
		in this function and its children.  Note that due to
		different viewpoints, functions excluded by options, etc,
		these numbers will NOT add up to 100%.

     self	This is the total amount of time spent in this function.

     children	This is the total amount of time propagated into this
		function by its children.

     called	This is the number of times the function was called.
		If the function called itself recursively, the number
		only includes non-recursive calls, and is followed by
		a `+' and the number of recursive calls.

     name	The name of the current function.  The index number is
		printed after it.  If the function is a member of a
		cycle, the cycle number is printed between the
		function's name and the index number.


 For the function's parents, the fields have the following meanings:

     self	This is the amount of time that was propagated directly
		from the function into this parent.

     children	This is the amount of time that was propagated from
		the function's children into this parent.

     called	This is the number of times this parent called the
		function `/' the total number of times the function
		was called.  Recursive calls to the function are not
		included in the number after the `/'.

     name	This is the name of the parent.  The parent's index
		number is printed after it.  If the parent is a
		member of a cycle, the cycle number is printed between
		the name and the index number.

 If the parents of the function cannot be determined, the word
 `<spontaneous>' is printed in the `name' field, and all the other
 fields are blank.

 For the function's children, the fields have the following meanings:

     self	This is the amount of time that was propagated directly
		from the child into the function.

     children	This is the amount of time that was propagated from the
		child's children to the function.

     called	This is the number of times the function called
		this child `/' the total number of times the child
		was called.  Recursive calls by the child are not
		listed in the number after the `/'.

     name	This is the name of the child.  The child's index
		number is printed after it.  If the child is a
		member of a cycle, the cycle number is printed
		between the name and the index number.

 If there are any cycles (circles) in the call graph, there is an
 entry for the cycle-as-a-whole.  This entry shows who called the
 cycle (as parents) and the members of the cycle (as children.)
 The `+' recursive calls entry shows the number of function calls that
 were internal to the cycle, and the calls entry for each member shows,
 for that member, how many times it was called from other members of
 the cycle.


Index by function name

  [11] BoundaryCheck           [2] NonGaussianDistance     [3] SimplexMaxLikelihood
  [10] CrossCorr               [4] PairExactLogLikelihood  [9] SimplexSort
   [6] Gaussian                [7] PairNGSBData            [5] ran2
   [8] GaussianDistance       [12] SeedRand
