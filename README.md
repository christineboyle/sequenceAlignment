# sequenceAlignment 
find minimum edit distance of two DNA sequences using dynamic programming

Author: Christy Boyle
File: EditDistance.java

Compile:
(command line)
javac EditDistance.java

Execute: 
(command line: any of the following)
java EditDistance (name of input file)

java EditDistance  < (name of input file)

java EditDistance
(then manually enter the input)

Implementation:
In order to implement this algorithm, I used a scanner (java class java.util.Scanner) to read the similarity matrix into a 2D array, the value of delta as an int, and the two sequences as strings. Then, I used a 2D array to contruct the matrix M. In order to construct the alignment sequence and associated costs, I used string buffers (java class java.lang.StringBuffer).
