/**
 * File: EditDistance.java
 * 
 * Finds the best alignment of two DNA sequences that yeilds the minimum edit distance
 *
 * @author Christy Boyle
 * @version 1.0 11/23/14
 */

import java.util.Scanner;
import java.io.*;
import java.lang.StringBuffer;

public class EditDistance{

	static int[][] simMatrix = new int[4][4];
	static int[][] alignMatrix;
	static int gapCost;
	static String allCosts;
	static String endSeq1;
	static String endSeq2;
	static int endCost;

	public static void main(String[] args){
		Scanner sc;
		try{
			if(args.length == 0){
				sc = new Scanner(System.in);
			}
			else{
			  sc = new Scanner(new File(args[0]));
			}

		    String seq1 = "";
		    String seq2 = "";

		    //read in the gap cost
			gapCost = sc.nextInt();

			//read in the similarity matrix
			for (int i = 0; i<4; i++){
				for(int j = 0; j<4; j++){
					simMatrix[i][j] = sc.nextInt();
				}
			}

			sc.nextLine();
			if(sc.hasNext()){
				seq1 = sc.nextLine();
			}
			if(sc.hasNext()){
				seq2 = sc.nextLine();
			}

			//method calls
			findOpt(seq1 , seq2);

			//printMatrix(seq1 , seq2);

			if (seq1.length() > 0 && seq2.length() > 0)
				backtrace(seq1 , seq2);

			printResults(seq1 , seq2);

	    }
		catch(IOException e){ System.out.println(e); }

	}

	/*
	 * This method constructs the dynamic programming array in order to calculate the minimum edit distance
	 * @param seq1 where seq1 is the first DNA sequence
	 * @param seq2 where seq2 is the second DNA sequence
	 */

	public static void findOpt(String seq1 , String seq2){
		int len1 = seq1.length();
		int len2 = seq2.length();

		alignMatrix = new int[len2 + 1][len1 + 1];

		//special case in which only one sequence is given

		if(len1 == 0 || len2 == 0){
			StringBuffer sb1 = new StringBuffer();
		    StringBuffer sb2 = new StringBuffer();
			StringBuffer co = new StringBuffer();
			if(len1 == 0){
				for(int i = 0; i<len2; i++){
					sb2.append(seq2.charAt(i) + " ");
					sb1.append("- ");
					co.append(gapCost + " ");
					}
					endCost = gapCost * len2;
				}
			if(len2 == 0){
				for(int i = 0; i<len1; i++){
					sb1.append(seq1.charAt(i) + " ");
					sb2.append("- ");
					co.append(gapCost + " ");
					}
					endCost = gapCost * len1;
				}

			
			allCosts = co.toString();
			endSeq2 = sb2.toString();
			endSeq1 = sb1.toString();
		}

	else{
		//initialize values in the dynamic programming array
		for(int i = 0; i<=len1; i++)
			alignMatrix[0][i] = i * gapCost;
        for(int j = 0; j<=len2; j++)
			alignMatrix[j][0] = j * gapCost;

		char char1;
		char char2;

		int a;
		int b;
		int c;

		for(int j = 1; j <= len2; j++){
			for(int i = 1; i <= len1; i++){
				char1 = seq1.charAt(i - 1);
				char2 = seq2.charAt(j - 1);

				a = alignMatrix[j - 1][i - 1] + indexSimMatrix(char1 , char2);
				b = alignMatrix[j][i - 1] + gapCost; 
				c = alignMatrix[j - 1][i] + gapCost;

				alignMatrix[j][i] = findMin(a , b , c);

			}
		}

		endCost = alignMatrix[len2][len1];
	}
}

	/*
	 * This method prints out the similarity matrix and the 2D array that calculates the min edit distance
	 * @param seq1 where seq1 is the first DNA sequence
	 * @param seq2 where seq2 is the second DNA sequence
	 */

	public static void printMatrix(String seq1 , String seq2){
		//print the similarity matrix for testing
			System.out.println("Similarity Matrix:");
			for (int k = 0; k<4; k++){
				if (k != 0)
					System.out.println("");
				for(int l = 0; l<4; l++){
					System.out.print(simMatrix[k][l] + " ");
				}
			}
			System.out.println("\n");
		//print the matrix for testing
			System.out.println("Alignment Matrix:");
			for (int k = 0; k<=seq2.length(); k++){
				if (k != 0)
					System.out.println("");
				for(int l = 0; l<=seq1.length(); l++){
					if(alignMatrix[k][l] > 9){
					    System.out.print(alignMatrix[k][l] + " ");
					}
					else{System.out.print(alignMatrix[k][l] + "  ");}
				}
			}
			System.out.println("");


		System.out.println("Minimum edit distance is " + endCost + ".\n");
	}


	/*
	 * This method performs the backtrace in order to find the best alignment
	 * @param seq1 where seq1 is the first DNA sequence
	 * @param seq2 where seq2 is the second DNA sequence
	 */

	public static void backtrace(String seq1 , String seq2){
		StringBuffer alignSeq1 = new StringBuffer();
		StringBuffer alignSeq2 = new StringBuffer();
		StringBuffer costs = new StringBuffer();

		int index1 = seq1.length();
		int index2 = seq2.length();



		while ((index1 > 0) || (index2 > 0)){
			if(index1 == 0){
				for(; index2 > 0; index2--){
					alignSeq1.append("- ");
					alignSeq2.append(seq2.charAt(index2 - 1) + " ");
					costs.append(gapCost + " ");
				}
				break;
			}

			if(index2 == 0){
				for(; index1 > 0; index1--){
					alignSeq2.append("- ");
					alignSeq1.append(seq1.charAt(index1 - 1) + " ");
					costs.append(gapCost + " ");
				}
				break;
			}

			int a = alignMatrix[index2 - 1][index1 - 1] + indexSimMatrix(seq2.charAt(index2 - 1) , seq1.charAt(index1 - 1));
			int b = alignMatrix[index2][index1 - 1] + gapCost;
			int c = alignMatrix[index2 - 1][index1] + gapCost;

			int smallest = findMin(a, b, c);

			if(smallest == a){
				alignSeq1.append(seq1.charAt(index1 - 1) + " ");
				alignSeq2.append(seq2.charAt(index2 - 1) + " ");
				costs.append(indexSimMatrix(seq1.charAt(index1 - 1) , seq2.charAt(index2 - 1)) + " ");
				index1--;
				index2--;
				}

			else if(smallest == c){
				alignSeq1.append("- ");
				alignSeq2.append(seq2.charAt(index2 - 1) + " ");
				costs.append(gapCost + " ");
				index2--;
				}				
				
			else{
				alignSeq1.append(seq1.charAt(index1 - 1) + " ");
				alignSeq2.append("- ");
				costs.append(gapCost + " ");
				index1--;
			}

		}

		allCosts = costs.reverse().toString();
		endSeq2 = alignSeq2.reverse().toString();
		endSeq1 = alignSeq1.reverse().toString();

	}

	/*
	 * This method prints the results of findOpt() and backtrace()
	 * @param seq1 where seq1 is the first DNA sequence
	 * @param seq2 where seq2 is the second DNA sequence
	 */

	public static void printResults(String seq1 , String seq2){
		//print the results
		System.out.println("The best alignment is\n");
		System.out.println(endSeq1);
		System.out.println(endSeq2);
		System.out.println(allCosts + "\n");
		System.out.println("With the minimum edit distance of " + endCost + ".");
	}

	/*
	 * This method constructs the dynamic programming array in order to calculate the minimum edit distance
	 * @param char1 where char1 represents the first index in the 2D array
	 * @param char2 where char2 represents the second index in the 2D array
	 * @return The value of the penalty cost
	 */

	public static int indexSimMatrix(char char1 , char char2){
		int simMatrixIndexS = 0;
		int simMatrixIndexT = 0;
				switch(char1){
					case 'A': simMatrixIndexS = 0;
							  break;
					case 'C': simMatrixIndexS = 1;
							  break;
					case 'G': simMatrixIndexS = 2;
							  break;
					case 'T': simMatrixIndexS = 3;
							  break;
				}
				switch(char2){
					case 'A': simMatrixIndexT = 0;
							  break;
					case 'C': simMatrixIndexT = 1;
							  break;
					case 'G': simMatrixIndexT = 2;
							  break;
					case 'T': simMatrixIndexT = 3;
							  break;
				}

				return simMatrix[simMatrixIndexS][simMatrixIndexT];	
		}


	/*
	 * This method returns the minimum of three integers
	 * @param x where x is the first integer
	 * @param y where y is the second integer
	 * @param z where z is the third integer
	 * @return the minimum of x, y, and z
	 */
	
	public static int findMin(int x, int y, int z){
		int min = z;
		if (min > y)
			min = y;		
		if (min > x)
			min = x;
		return min;
	}

}


