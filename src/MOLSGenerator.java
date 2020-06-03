import java.util.HashSet;

public class MOLSGenerator {
	
	public int d;
	public int p;
	public int q;
	
	public int poly;
	
	public HashSet<Integer[]> rows;
	
	
	public MOLSGenerator(Integer q, Integer p) { //p must be a prime.
		this.q = q;
		this.d = (int) Math.round(logb(q, p));
		
		System.out.println(d);
		this.p = p;
		
		this.rows = new HashSet<Integer[]>();
		
		if (q == 16) {
			this.poly = 0b10011; //x^4 + x + 1
		} else if (q == 8) {
			this.poly = 0b1101; //x^3 + x + 1
		} else if (q == 256) {
			this.poly = 0x11b;
		}
		
		System.out.println(this.poly);
		
		Integer[][] addTable = new Integer[q][q];
		Integer[][] multTable = new Integer[q][q];
		
		for (int i = 0; i < q; i++) {
			for (int j = 0; j < q; j++) {
				addTable[i][j] = GFAdd(i, j, p);
			}
		}
		
		printArray2D(addTable);
		
		for (int i = 0; i < q; i++) {
			for (int j = 0; j < q; j++) {
				multTable[i][j] = GFMult(i, j, p, d);
			}
		}
		printArray2D(multTable);
		
		System.out.println("The MOLS for q = " + q + " are:");
		
		Integer[][][] MOLS = new Integer[q-1][q][q];
		
		// Now we can generate the q-1 MLOS using the definition of product in GF(q).
		for (int a=1; a < q; a++) {
			Integer[][] mols = new Integer[q][q];
			for (int x=0; x<q; x++) {
				for (int y=0; y< q; y++) {
					mols[x][y] = GFAdd(GFMult(a,x, p, d), y, p) + 1;
				}
			}
			printArray2D(mols);
			MOLS[a-1] = mols;
		}
		
		System.out.println("And the corresponding bundles are:");
		Integer[][][] bundles = MOLSto1Intersecting(MOLS, q);
		Integer[][][] bundles2 = MOLSto1Intersecting(bundles, q);
		
		
		for (int a2 = 0; a2 < q-1; a2++) {
			Integer[][] newGroup = bundles[a2];
			Integer[][] newGroup2 = bundles2[a2];
			for (int i = 0; i < q; i++) {
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

		}
		for (int i2 = 0; i2 < q-1; i2++) {
			System.out.println("MOLS:");
			this.printArray2D(bundles2[i2]);
			System.out.println("Corresponding Bundle:");
			this.printArray2D(bundles[i2]);

		}

	}
	
	public Integer[][][] MOLSto1Intersecting (Integer[][][] MOLS, int q){
		Integer[][][] bundles = new Integer[q-1][q][q];
		
		for (int a = 0; a < q-1; a++) {
			Integer[][] M = MOLS[a];
			Integer[][] A = new Integer[q][q];
			for (int i = 0; i < q; i++) {
				for (int j = 0; j<q; j++) {
					int k = M[i][j];
					A[k-1][j] = i+1;
				}
			}
			bundles[a] = A;
			printArray2D(A);
		}
		return bundles;
	}
	
	// as per https://people.cs.clemson.edu/~westall/851/rs-code.pdf
	public int GFAdd(int v1, int v2, int p) {
		if (p==2) {
			return (v1 ^ v2);
		}
		// not implemented: other primes
		return 0;
	}

	public int GFMult(int v1, int v2, int p, int d) {
		if (p==2) {
			int prod = 0;
			int k; 
			int mask;
			
			//bitwise multiplication through the powers of the binary representations
			for (k=0; k < d; k++) {
				if ((v1 & 1) == 1) {
					prod ^= (v2 <<k);
				}
				v1 >>= 1;
				if (v1 == 0) break;
			}
			
			// degree of the polynomial we have is d
			// shift the irreducible polynomial by (o-d) bits where o is the current order of prod
			// as per https://www.doc.ic.ac.uk/~mrh/330tutor/ch04s04.html
			int deg = getDegree(prod);
			while (deg >= d) {
				int numShift = deg - d;
				int tmp = this.poly << numShift;
				prod = prod ^ tmp;
				deg = getDegree(prod);
			}

			return prod;
		}
		return 0;
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
	
	public static double logb( double a, double b )
	{
	return Math.log(a) / Math.log(b);
	}
	
	public static int getDegree(int num) {
		int deg = 0;
	    while (((num >>= 1) != 0))
	        deg += 1;
		return deg;
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
		MOLSGenerator m = new MOLSGenerator(8, 2);
	}

}
