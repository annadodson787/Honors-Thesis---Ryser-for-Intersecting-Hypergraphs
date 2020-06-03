import java.util.HashMap;
import java.util.HashSet;

public class PermutationCycle {
	//public Integer[] seed = {1, 3, 5, 7, 4, 8, 6, 2}; // {1, 5, 7, 3, 11, 15, 13, 9, 16, 4, 10, 6, 12, 8, 14, 2}; //{1, 3, 5, 7, 4, 8, 6, 2}; //
	public Integer[] firstRow;
	
	public HashSet<Integer[]> rows;
	
	public int d;
	public int p;
	public int N;
	
	public HashMap<Integer, Integer[]> permutations;
	
	
	public PermutationCycle(int N, int p) {
		this.N = N;
		this.d = (int) Math.round(logb(N, p));
		this.p = p;
		this.firstRow = new Integer[N];
		this.rows = new HashSet<Integer[]>();
		for (int i = 0; i < N; i++) {
			this.firstRow[i] = i+1;
		}
		//this.generateDivisorSequenceV2(this.N, this.p, this.d);
		Integer[] divisorSequence = this.generateDivisorSequence(this.N, this.p);
		this.printArray1D(divisorSequence);
		System.out.println('\n');
		
		this.permutations=this.generatePermutationSet(this.firstRow, divisorSequence);
		
		//Integer[][] group1 = generateBundle(this.firstRow) //this.seed);

		Integer[][] group1 = generatePermutations(this.firstRow, divisorSequence);
		this.printArray2D(group1);
		//this.addAllRows(group1);
		
		Integer[][][] bundles = new Integer[this.N][this.N][this.N];
		
		Integer[][] currentGroup = group1;
		bundles[0] = group1;
		
		for (int groupnum = 2; groupnum < N; groupnum++) {
			Integer[] newSeed = getCol1(currentGroup);//getCol1(currentGroup); //getDiag(currentGroup);
			//GridPermuter2.printArray1D(diagonal);
			//System.out.println();
			Integer[][] newGroup = generatePermutations(newSeed, divisorSequence); //this.seed);
			this.printArray2D(newGroup);
			
			for (int i = 0; i < N; i++) {
				boolean cycle = this.checkRowIsFirstRow(newGroup[0]);
				if (cycle) {
					System.out.println("Cycle occurred on group number " + groupnum);
					//System.exit(1);
				}
			}
			
			
			for (int i = 0; i < N; i++) {
				boolean works = this.checkRow(newGroup[i]);
				
				if (works == false) {
					boolean identical = this.isSameAsPreviousRow(newGroup[i]);
					
					System.out.print("Row ");
					if (identical == true) System.out.print("identical: ");
					this.printArray1D(newGroup[i]);
					System.out.println(" did not work.");
				}
			}
			this.addAllRows(newGroup);
			
			currentGroup = newGroup;
			bundles[groupnum-1] = currentGroup; 
		}
		
		Integer[][][] MOLS = MOLSto1Intersecting(bundles, N);
		for (int i2 = 0; i2 < N-1; i2++) {
			Integer[][] fixedMOLS = rearrangeSymbols(MOLS[0][0], MOLS[i2]);

			System.out.println("Bundle:");
			PermutationCycle.printArray2D(bundles[i2]);
			//PermutationCycle.printArray2D(MOLS[i2]);
			System.out.println("Reduced MOLS:");
			
			PermutationCycle.printArray2D(fixedMOLS);
		
			
			
			//System.out.println("MOLS:");
			//PermutationCycle.printArray2D(MOLS[i2]);
//			System.out.println("MOLS:");
//			//PermutationCycle.printArray2D(fixedMOLS);
//			PermutationCycle.printArray2D(MOLS[i2]);
//			
//			System.out.println("Corresponding Bundle:");
//			//PermutationCycle.printArray2D(MOLSto1Intersecting(rearrangeSymbols(MOLS[0][0], bundles[i2]), N));
//			PermutationCycle.printArray2D(bundles[i2]);
//			
//			//PermutationCycle.printArray2D(MOLSto1Intersecting(fixedMOLS, N));
		}
		
	}
	
	public Integer[] inverse(Integer[] original) {
		Integer[] inv = new Integer[this.N];
		for (int i = 0; i < this.N; i++) {
			inv[original[i]-1] = i+1;
		}
		return inv;
	}
	
	public Integer[][] rearrangeSymbols(Integer[] I, Integer[][] initialSquare) {
		Integer[][] replace = new Integer[this.N][this.N];
		Integer[] inv = this.inverse(I);
		for (int i = 0; i < this.N; i++) {
			for (int j = 0; j < this.N; j++) {
				int symbol = initialSquare[i][j];
				replace[i][j] = inv[symbol-1];
			}
		}
		return replace;
	}
	
	public Integer[][] MOLSto1Intersecting(Integer[][] MOLS, int q){
		Integer[][] A = new Integer[q][q];
		for (int i = 0; i < q; i++) {
			for (int j = 0; j<q; j++) {
				int k = MOLS[i][j];
				A[k-1][j] = i+1;
			}
		}	
		return A;
	}
	
	
	public Integer[][][] MOLSto1Intersecting (Integer[][][] MOLS, int q){
		Integer[][][] bundles = new Integer[q-1][q][q];
		
		for (int a = 0; a < q-1; a++) {
			Integer[][] M = MOLS[a];
			Integer[][] A = MOLSto1Intersecting(M, q);
			bundles[a] = A;
			printArray2D(A);
		}
		return bundles;
	}
	
	public Integer[][] generateBundle(Integer[] row, Integer[] seed) {
		Integer[][] bundle = new Integer[this.N][this.N];
		
		for (int i = 0; i< N; i++) {
			Integer[] permutation = this.permutations.get(seed[i]);
			
			bundle[i] = this.applyPermutation(row, permutation);
		}
		return bundle;
	}
	
	public Integer[] getCol1(Integer[][] group) {
		Integer[] col = new Integer[this.N];
		for (int i = 0; i < this.N; i++) {
			col[i] = group[i][0];
		}
		return col;
	}
	
	public Integer[] applyPermutation(Integer[] row, Integer[] permutation) {
		Integer[] permutedRow = new Integer[this.N];
		
		for (int i = 0; i < this.N; i++) {
			permutedRow[i] = row[permutation[i]-1];
		}
		return permutedRow;
	}
	
//	public Integer[][] generatePermutationsV2(Integer[] firstRow, Integer[] divisorSequence) {
//		Integer[][] finalArray = new Integer[this.N][this.N];
//		finalArray[0] = firstRow;
//		
//		// [ 1 1 3 1 1 3 1 1]
//				
//		Integer[] currentSeed = firstRow;
//		
//		for (int i = 1; i < N; i++) {
//			int divisor = divisorSequence[i];
//			Integer[] newRow = new Integer[this.N];
//
//			for (int i2 = 0; i2 < N; i2++) {
//	}
	
	public Integer[][] generatePermutations(Integer[] firstRow, Integer[] divisorSequence) {
		Integer[][] finalArray = new Integer[this.N][this.N];
		finalArray[0] = firstRow;
		
		HashMap<Integer, Integer[]> perms = new HashMap<Integer, Integer[]>();
		
		Integer[] currentRow = firstRow;
		
		for (int i = 1; i < N; i++) {
			int divisor = divisorSequence[i];
			Integer[] newRow = new Integer[this.N];

			for (int i2 = 0; i2 < N; i2++) {
				for (int toggle = 0; toggle < this.p; toggle++) {
					int newLocal = (i2 + divisor) % (divisor*this.p);
					int start = i2/(this.p*divisor) * this.p*divisor; 
					//System.out.println("i2 " + i2 + " local " + newLocal + " start " + start);
					int end = newLocal+start;
					//System.out.println("setting element " + end + " from element " + i2);
					newRow[end] = currentRow[i2];
				}
			}
			finalArray[i] = newRow;
			currentRow = newRow;
		}

		return finalArray;
	}
	
	public HashMap<Integer, Integer[]> generatePermutationSet(Integer[] firstRow, Integer[] divisorSequence) {
		Integer[][] finalArray = new Integer[this.N][this.N];
		finalArray[0] = firstRow;
		
		HashMap<Integer, Integer[]> perms = new HashMap<Integer, Integer[]>();
		
		Integer[] currentRow = firstRow;
		
		for (int i = 1; i < N; i++) {
			int divisor = divisorSequence[i];
			Integer[] newRow = new Integer[this.N];

			for (int i2 = 0; i2 < N; i2++) {
				// the "local" index is the place we are in this set of (p) elements - (0 or 1) in 2 case
				int newLocal = (i2 + divisor) % (divisor*this.p);
				// the "start" is where the set of (p) elements begins
				int start = i2/(this.p*divisor) * this.p*divisor; 
				//System.out.println("i2 " + i2 + " local " + newLocal + " start " + start);
				
				// 
				int end = newLocal+start;
				//System.out.println("setting element " + end + " from element " + i2);
				newRow[end] = currentRow[i2];
			}
			finalArray[i] = newRow;
			currentRow = newRow;
		}
		
		for (int i = 0; i < N; i++) {
			perms.put(finalArray[i][0], finalArray[i]);
		}
		return perms;
	}
	
	public static double logb( double a, double b )
	{
	return Math.log(a) / Math.log(b);
	}
	
	public boolean checkRowIsFirstRow(Integer[] newRow) {
		for (int i = 0; i < this.N; i++) {
			if (newRow[i] != i+1) return false;
		}
		return true;
	}
	
	// this method is going to be O(N^2). 
	public boolean checkRow(Integer[] newRow) {
		boolean checksout = true;
		for (Integer[] row: this.rows) {
			boolean good = false;
			for (int i = 0; i < row.length; i++) {
				if (row[i] == newRow[i] && good == false) {
					good = true;
				} else if (row[i] == newRow[i] && good == true) {
					System.out.println("Did not work because 2+ vertices overlap!");
					printArray1D(row);
					printArray1D(newRow);
					System.out.println();
					return false; // 2 edges overlapping
				} else {} //do nothing
			}
			if (good == false) {
				System.out.println("Did not work because no vertices overlap!");
				printArray1D(row);
				printArray1D(newRow);
				System.out.println();
				checksout = false;
			}
		}
		return checksout;	
	}
	
	public boolean isSameAsPreviousRow(Integer[] newRow) {
		boolean checksout = true;
		for (Integer[] row: this.rows) {
			boolean match = false;
			if (newRow[0] == row[0]) match = true;
			for (int i = 0; i < row.length; i++) {
				if (newRow[i] != row[i]) {
					match = false;
				}
			}
			if (match == false) checksout = false;
		}
		return checksout;
	}
	
	public void addAllRows(Integer[][] group) {
		for (int i = 0; i < group.length; i++) {
			this.rows.add(group[i]);
		}
	}
	
	public Integer[] generateDivisorSequenceV2(int N, int p, int d) {
		Integer[] divisorSequence = new Integer[N];
		divisorSequence[0] = 0;
		for (int i = 1; i < N; i++) {
			for (int pow = 0; pow < d; pow++) {
				int power = (int) Math.pow(p, pow);
				//System.out.println(("i is " + i + " power is " + power + " and mod is " + (power%i)));
				if ((i % power) == 0) {
					divisorSequence[i] = power;
					//System.out.println("Setting element " + i + " to " + power );
				}
			}
		}
		return divisorSequence;
	}
	
	public Integer[] generateDivisorSequence(int N, int M) {
		// Generate odd-power sequence
		Integer[] divisorSequence = new Integer[N];
		boolean even = false;//isEven(N);
		
		int spacing = 1;
		int numPlace = N/M;
		while (true) {
			for (int i = 0; i < N; i+= spacing) {
				divisorSequence[i] = numPlace;
			}
			if (numPlace == 1) break;
			spacing = spacing*M;
			numPlace = numPlace/M;
		}
		
		if (even == true) {
			for (int i = 0; i < N; i ++) {
				divisorSequence[i] = divisorSequence[i]/M;
				if (divisorSequence[i] == 0) {
					// we need to treat this as a double local flip; 1 then 8.
					//divisorSequence[i] = N/2;
				}
			}
		}
		
		return divisorSequence;
	}
	
	public static void printArray2D(Integer[][] array2D) {
		for (int i = 0; i < array2D.length; i++) {
			for (int j = 0; j < array2D[0].length; j++) {
				String sp;
				if (array2D[i][j] < 10) { 
					sp = "  ";
				} else {
					sp = " ";
				}
				System.out.print(array2D[i][j] + sp);
			}
			System.out.print('\n');
		}
		System.out.print('\n');
	}
	
	public static void printArray1D(Integer[] array1D) {
		for (int i = 0; i < array1D.length; i++) {
			String sp;
			if (array1D[i] < 10) { 
				sp = "  ";
			} else {
				sp = " ";
			}
			System.out.print(array1D[i] + sp);
		}
		System.out.print('\n');
	}

	
	public static void main(String[] args) {
		
		PermutationCycle p = new PermutationCycle(8, 2);
	}
	
}
