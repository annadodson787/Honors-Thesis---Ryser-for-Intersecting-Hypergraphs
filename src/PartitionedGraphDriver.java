import java.util.ArrayList;
import java.util.Arrays;
import java.util.Dictionary;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.Map.Entry;

import com.opencsv.CSVReader;
import com.opencsv.CSVWriter;

import java.awt.List;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException; 
import java.io.InputStreamReader; 

public class PartitionedGraphDriver {
	
	public int N;
	public int M;
	
	public int maxDepth;
	public int minDepth;
	
	public static void main(String[] args) {
		PartitionedGraphDriver t = new PartitionedGraphDriver(1, 4);
	}
	
	public PartitionedGraphDriver(int M, int N) {
		this.N = N;
		this.M = M;
		
		this.maxDepth = 0;
		this.minDepth = 100;
		
		//playGreedyGame(M, N);
		//playBruteForceGame(M, N);
		boolean greedy = true;
		boolean overlap = false;
		playBruteForceGameRecursive(M, N, greedy, overlap);
		
		try {
			//playCLIGame(M, N);
		} catch (Exception e) {
			System.out.print(e);
		}
		
//		try {
//			System.out.print(formatTreeTable("moves.csv"));
//		} catch (Exception e) {
//			System.out.print(e);
//		}
	}
	
	public void playCLIGame(int M, int N) throws IOException {
		System.out.println("You are playing Bob. Alice will give you a vertex cover, and you must specify some edge which avoids it, "
				+ "and has something in common with all of your other edges. \n");
		
		System.out.println("        Alice                   Bob                                   Counts \n");		
		
		PartitionedGraphRestructured g = new PartitionedGraphRestructured(N);

		Vertex[] aliceMove = new Vertex[(N-1)*M]; // all 0's to start off; it'll actually be 
		
		Edge firstMove = g.generateMove(M, N, aliceMove, true);
		g.addEdge(firstMove, true);
		System.out.println("Bob made a move: " + firstMove);
		
		aliceMove = greedyAlice(M, N, g.vertexToEdgeMap);
		System.out.println("Alice made a move: " + formatLine(g.edges.size(), aliceMove));
		
		int numGames = 0;
		
		BufferedReader reader =  
                new BufferedReader(new InputStreamReader(System.in)); 
      
	    // Reading data using readLine 
	    String edgeString = reader.readLine(); 

		while (edgeString!=null) {
			HashSet<Edge> generateMoves = g.generateAllMoves(M, N, aliceMove);
			
			if (generateMoves.isEmpty()) {
				System.out.println("You lose! No more moves are possible. The transcript for this game is: ");
				System.out.println(g.stringRepresentation());
			}
			
			if (edgeString.equals("")) {
				System.out.print("Your valid moves are: ");
				System.out.println(generateMoves);
				
			} else {
				String[] arr = edgeString.split(" ");
				Edge e = new Edge(N);
				for (String vertexName: arr) {
					if (vertexName.matches("[a-z][1-9]")) {
						if (g.stringToVertexMap.containsKey(vertexName)) {
							e.setVertex(g.stringToVertexMap.get(vertexName), 0);
						} else {
							Vertex newVertex = new Vertex(vertexName);
							// TODO: Fix handling of invalid strings
							e.setVertex(newVertex, 0);
						}	
					}
				}
				boolean flag = false;
				for (Edge ej: generateMoves) {
					if (ej.hashCode() == e.hashCode()) {
						flag = true;
					}
				}
				
				if (flag == true) {
					boolean greedy = g.candidateEdgeIsGreedy(M, N, e, aliceMove);

					g.addEdge(e, greedy);
					System.out.println("You made a move: " + e);
					System.out.println("Your move was greedy: " + greedy);
					aliceMove = greedyAlice(M, N, g.vertexToEdgeMap);
					System.out.println("Alice made a move: " + formatLine(g.edges.size(), aliceMove));
				} 
				else {
					System.out.print("Invalid move. Your valid moves are: ");
					System.out.println(generateMoves);
				}

			}
			edgeString = reader.readLine();
		}
	}
	
	public void playBobOverlapGame(int N, int r) {
		// This kind of game has no vertex cover. 
		// Instead, Bob is constrained in the number of overlaps he can have with his previous edges.
		PartitionedGraphRestructured g = new PartitionedGraphRestructured(N);
		
		Vertex[] aliceMove = new Vertex[g.degree-1]; // all 0's, forever. 
		boolean add_vc = false;
		Edge bobMove = g.generateMove(M, N, aliceMove, add_vc);
		
	}
	
	public void playGreedyGame(int M, int N) {
		System.out.println("        Alice                   Bob                               Counts \n");		
		PartitionedGraphRestructured g = new PartitionedGraphRestructured(N);
		
		int iter = 0;
		Vertex[] aliceMove = new Vertex[g.degree-1]; // all 0's to start off
		System.out.print(formatLine(iter, aliceMove)); 
		
		boolean add_vc = true;
		Edge bobMove = g.generateMove(M, N, aliceMove, add_vc);
		
		while (true) {
			System.out.println(bobMove);
			boolean greedy = g.candidateEdgeIsGreedy(M, N, bobMove, aliceMove);
			g.addEdge(bobMove, greedy);
			
			// Print the updated counts
			System.out.println(g.generateCountString(N));
			iter++;
			
			// Alice reacts
			aliceMove = greedyAlice(M, g.degree, g.vertexToEdgeMap);
			System.out.print(formatLine(iter, aliceMove));
			
			bobMove = g.generateMove(M, N, aliceMove, add_vc);
			if (bobMove == null) {
				break;
			}
			
		}	
	}
	
	public void playBruteForceGameRecursive(int M, int N, boolean greedy, boolean overlap) {
		File file;
		String filename = "greedy.csv";
		file = new File(filename); 
		try { 
			FileWriter outputfile = new FileWriter(file); 
			CSVWriter writer = new CSVWriter(outputfile);
			
			FileWriter logWriter = new FileWriter("gamelog.txt", true);

			PartitionedGraphRestructured g = new PartitionedGraphRestructured(N);

			Vertex[] aliceMove = new Vertex[(N-1)*M]; // all 0's to start off
			Edge firstMove = g.generateMove(M, N, aliceMove, true); // bob's first move is always u1...p_N1

			Counter count = new Counter();
			Counter vcount = new Counter();
			HashSet<Edge> possibleMoves = new HashSet<Edge>();
			possibleMoves.add(firstMove);
			
			if (overlap == false) {
				bruteForceGameRecursiveHelper(M, N, g, possibleMoves, count, vcount, aliceMove, writer, filename, greedy, logWriter);
			} else {
				bruteForceGameOverlapRecursiveHelper(M, N, g, possibleMoves, count, vcount, writer, filename, logWriter);
			}

			logWriter.close();
			writer.close(); 
		}   catch (IOException e) {  
			e.printStackTrace(); 
		}  
	}
		
	public void bruteForceGameRecursiveHelper(int M, int N, PartitionedGraphRestructured g, HashSet<Edge> candidateMoves, Counter count, Counter newVertCount, Vertex[] aliceMove, CSVWriter writer, String filename, boolean greedy, FileWriter logWriter) { //, HashSet<HashSet<Edge> allPreviousMoves) {
		if (candidateMoves.isEmpty()) {
			//End this game and print out the current string.
			g.s = new StringBuilder();
			for (String str: g.stringBuilderList) {
				g.s.append(str);
			}
			//if (g.depth>30 || g.depth<14) {
			if (g.depth >= 13) {//this.maxDepth) {
				count.increment();
				this.maxDepth = g.depth;
				System.out.println();
				System.out.println();
				System.out.print(count.count);
				System.out.print(g.s.toString());
				System.out.println("No possible move.");
			}
//			if (g.depth <= this.minDepth) {
//				this.minDepth = g.depth;
//				System.out.print(g.s.toString());
//				System.out.println("No possible move.");
//			}
			
			//}
			//System.out.print(g.s.toString());
			//System.out.println("No possible move.");
			//System.out.println("Bob just added a new vertex: " + g.newVert);

			if (!filename.equals("dummy.csv")) {
				Object[] array = g.scoreString.toArray();
				String[] sarray = Arrays.copyOf(array, array.length, String[].class);
				writer.writeNext(sarray);
			}
			
			if (g.newVert) {
				newVertCount.increment();
			}
			
			//System.out.println("Total count is now: " + count.count + ". Ending new vertices count: " + newVertCount.count + ". Score for this game: " + g.score + ", measure for this game: " + g.measureScore +  ", depth: " + g.depth +".\n");
			
		} else {
			// Consider a game where we add each edge from the candidate moves.
			for (Edge candidate: candidateMoves) {
				//System.out.println("Depth is: "+g.depth);
				boolean isGreedy = g.candidateEdgeIsGreedy(M, N, candidate, aliceMove);
				int measure = g.scoreCandidateEdge(candidate);
				int hitScore = candidate.hitScore;
				if (measure != candidate.measure) {
					System.err.println("error the measure score was: " + measure + ", while the edge's measure was " + candidate.measure);
				} else {
					//System.out.println("measure score was " + measure + ", edge score was " + candidate.measure);
				}
				
				g.addEdge(candidate, isGreedy);
				Vertex[] newAliceMove = greedyAlice(M, N, g.vertexToEdgeMap);
				String scoreString = " " ; //movesString = candidate.toString();
				HashSet<Edge> newCandidateMoves = g.generateAllMoves(M, g.degree, newAliceMove);
				
				if (filename.equals("gameScores.csv")) {
					g.movesString.add(scoreString);
				}
				else if (filename.equals("hitScores.csv")) {
					g.movesString.add(String.valueOf(hitScore));
				}
				else if (filename.equals("options.csv")) {
					g.movesString.add(newCandidateMoves.toString());
				} else if (filename.equals("depths.csv")) {
					g.movesString.add(Integer.toString(newCandidateMoves.size()));
				} else if (filename.equals("greedy.csv")) {
					g.movesString.add(String.valueOf(isGreedy));
				}
				
				bruteForceGameRecursiveHelper(M, N, g, newCandidateMoves, count, newVertCount, newAliceMove, writer, filename, greedy, logWriter);
				boolean add_vc = true;
				g.removeEdge(candidate, add_vc);
				
				if (!filename.equals("dummy.csv")) {
					int ind  = g.movesString.size() - 1;
					g.movesString.remove(ind);
				}
			}

		}
	}
	
	public void bruteForceGameOverlapRecursiveHelper(int M, int N, PartitionedGraphRestructured g, HashSet<Edge> candidateMoves, Counter count, Counter newVertCount, CSVWriter writer, String filename, FileWriter logWriter) { //, HashSet<HashSet<Edge> allPreviousMoves) {
		if (candidateMoves.isEmpty()) {
			//End this game and print out the current string.
			g.s = new StringBuilder();
			for (String str: g.stringBuilderList) {
				g.s.append(str);
			}
			if (g.depth >= this.maxDepth) {
				this.maxDepth = g.depth;
				//System.out.print(g.s.toString());
				//System.out.println("No possible move.");
			}
			if (g.depth <= this.minDepth) {
				this.minDepth = g.depth;
				//System.out.print(g.s.toString());
				//System.out.println("No possible move.");
			}
			try {
				logWriter.write(count.count);
				logWriter.write('\n');
				logWriter.write(g.s.toString());
				logWriter.write('\n');
			} catch (IOException e) {
				System.out.println("An error occurred.");
				e.printStackTrace();
			}

			if (!filename.equals("dummy.csv")) {
				Object[] array = g.movesString.toArray();
				String[] sarray = Arrays.copyOf(array, array.length, String[].class);
				writer.writeNext(sarray);
			}
			
			if (g.newVert) {
				newVertCount.increment();
			}
			count.increment();
			//System.out.println("Total count is now: " + count.count + ". Ending new vertices count: " + newVertCount.count + ". Score for this game: " + g.score + ", measure for this game: " + g.measureScore +  ", depth: " + g.depth +".\n");
			
		} else {
			// Consider a game where we add each edge from the candidate moves.
			//System.out.println("Candidate moves: " + candidateMoves);
			for (Edge candidate: candidateMoves) {
				//System.out.println("Depth is: "+g.depth);
				boolean flag=true;
				int maxVal = (N-1)*(N-1);
				for (Vertex v: candidate.values) {
					if (v.value > maxVal) {
						flag=false;
					}
				}
				if (flag==true) {
					g.addEdge(candidate, true);
					Vertex[] newAliceMove = greedyAlice(M, N, g.vertexToEdgeMap);
					
					
					boolean add_vc = true;
					HashSet<Edge> newCandidateMoves = g.generateOverlapMoves(N, false, add_vc, newAliceMove);
				
					bruteForceGameOverlapRecursiveHelper(M, N, g, newCandidateMoves, count, newVertCount, writer, filename, logWriter);
					
					g.removeEdge(candidate, add_vc);
					
					//metrics can come later
					//if (!filename.equals("dummy.csv")) {
					//	int ind  = g.movesString.size() - 1;
					//	g.movesString.remove(ind);
					//}
				}
			}
			
		}
	}
	
	// If we want to recursively generate gameplay, we need to build the gameplay string, rather than printing it out.
	// Alternatively, we could associate each "Game" with a new PartitionedGraphRestructured, and build off the previous graphs. 
	public void playBruteForceGame(int M, int N) {
		LinkedList<PartitionedGraphRestructured> graphStack = new LinkedList<PartitionedGraphRestructured>();
		
		PartitionedGraphRestructured g = new PartitionedGraphRestructured(N);
		
		graphStack.addLast(g);
		
		Vertex[] aliceMove = new Vertex[(N-1)*M]; // all 0's to start off; it'll actually be 
		
		Edge firstMove = g.generateMove(M, N, aliceMove, true);
		g.addEdge(firstMove, true);
		
		int numGames = 0;

		while (!graphStack.isEmpty()) {

			g = graphStack.removeLast();
			aliceMove = greedyAlice(M, N, g.vertexToEdgeMap);
			//System.out.println("Alice made a move: " + formatLine(g.edges.size(), aliceMove));
			
			HashSet<Edge> generateMoves = g.generateAllMoves(M, N, aliceMove);
			//Edge greedyMove = g.generateMove(M, N, aliceMove, false);
			
			if (generateMoves.isEmpty()) {
				System.out.println(g.stringRepresentation());
				numGames = numGames + 1;
			}
			
			//System.out.println("The number of possible moves for this move is: " + generateMoves.size());
			//System.out.println("The possible moves are: " + generateMoves);
			
			for (Edge move: generateMoves) {
				//System.out.println("Bob made a move: " + move);
				PartitionedGraphRestructured newGraph = new PartitionedGraphRestructured(g);
				boolean greedy = g.candidateEdgeIsGreedy(M, N, move, aliceMove);
				newGraph.addEdge(move, greedy);
				//newGraph.addEdge(move);
				
				//newGraph.s.append("                                              Greedy Bob: ");
				//newGraph.s.append(greedyMove + "\n");
				
				//System.out.println("Now the graph's game looks like: \n" + newGraph.s.toString());
				//System.out.println("Now the graph's edges are: " + newGraph.edges);
				//System.out.println("Now the graph's counts are: " + newGraph.generateCountString(N));

				graphStack.addLast(newGraph);
				System.out.println("Current stack size: " + graphStack.size() + ", current total games: " + numGames);
			}
			
		}
	}
	
	public Vertex[] aliceAlgorithm(int M, int N, HashMap<Vertex, HashSet<Edge>> vertexToEdgeMap) {
		int vcSize = (N-1)*M;
		Vertex[] vertexCover = new Vertex[vcSize];
		
		HashMap<Vertex, HashSet<Edge>> verticesLeft = createVertMapClone(vertexToEdgeMap);
		
		for (int i = 0; i < vcSize; i++) {
			int max = 0;
			int minpart = -1;
			int minval = -1;
			Vertex bestVertex = null;
			for (Entry<Vertex, HashSet<Edge>> ent: verticesLeft.entrySet()) {
				int size = ent.getValue().size();
				Vertex v = ent.getKey();
				int part = v.part.value;
				int val = v.value;					
				if (size > max) {
					bestVertex = v;
					minpart = part;
					minval = val;
					max = size;
				}
				else if (size == max) {
					if (minpart == -1 || part < minpart) {
						minpart = part;
						bestVertex = v;
					}
					else if (part == minpart) {
						if (minval == -1 || val < minval) {
							minval = val;
							bestVertex = v;
						}
					}
				
				}
			}
			vertexCover[i] = bestVertex;
			verticesLeft.remove(bestVertex);
		}
		
		return vertexCover;
	}
	
	public Vertex[] greedyAlice(int M, int N, HashMap<Vertex, HashSet<Edge>> vertexToEdgeMap) {
		int vcSize = (N-1)*M;
		Vertex[] vertexCover = new Vertex[vcSize];
		
		HashMap<Vertex, HashSet<Edge>> verticesLeft = createVertMapClone(vertexToEdgeMap);

		for (int i = 0; i < vcSize; i++) {
			int max = 0;
			int minpart = -1;
			int minval = -1;
			Vertex bestVertex = null;
			for (Entry<Vertex, HashSet<Edge>> ent: verticesLeft.entrySet()) {
				int size = ent.getValue().size();
				Vertex v = ent.getKey();
				int part = v.part.value;
				int val = v.value;
				
				if (size > max) {
					bestVertex = v;
					minpart = part;
					minval = val;
					max = size;
				}
				else if (size == max) {
					if (minpart == -1 || part < minpart) {
						minpart = part;
						minval = val;
						//System.out.println("Updating bestvertex " + bestVertex + " to " + v + " because part was " + part);
						bestVertex = v;
					}
					else if (part == minpart) {
						if (minval == -1 || val < minval) {
							minval = val;
							//System.out.println("Updating bestvertex " + bestVertex + " to " + v + " because val was " + val);
							bestVertex = v;
							
						}
					}
				
				}
			}
			
			//now, remove the edges in the vertex to edge map which had this vertex in it too.
			// Remove all edges to hit which contain this vertex
			HashSet<Edge> edgesWhichAreNowHit = new HashSet<Edge>();
			HashSet<Edge> ptr = verticesLeft.get(bestVertex);
			
			if (ptr != null) {
				edgesWhichAreNowHit.addAll(ptr);
				for (Edge hit: edgesWhichAreNowHit) {
					Vertex[] containedVertices = hit.values;
					for (Vertex v: containedVertices) {
						HashSet<Edge> tmp = verticesLeft.get(v);
						tmp.remove(hit);
						verticesLeft.put(v, tmp);
					}
				}
			}
			
			vertexCover[i] = bestVertex;
			verticesLeft.remove(bestVertex);
		}
		
		return vertexCover;
	}
	
	
	class PartitionedGraphRestructured {
		public int degree;	
		public HashSet<Edge> edges;
		public HashMap<Vertex, HashSet<Edge>> vertexToEdgeMap;
		public HashMap<Integer, HashSet<Vertex>> partToVertexMap;
		public HashMap<String, Vertex> stringToVertexMap;
		
		public boolean newVert;
		
		public int score;
		public int depth;
		public int measureScore;
		
		public ArrayList<String> stringBuilderList;
		public StringBuilder s;
		public ArrayList<String> scoreString;
		public ArrayList<String> movesString;
		
		public PartitionedGraphRestructured (int N) {
			this.degree = N;
			this.edges = new HashSet<Edge>();
			
			this.score = 0; // increment by x whenever we add an edge which hits x other edges.
			this.depth = 0;
			this.measureScore = 0;
		
			this.vertexToEdgeMap = new HashMap<Vertex, HashSet<Edge>>();
			this.partToVertexMap = new HashMap<Integer, HashSet<Vertex>>();
			this.stringToVertexMap = new HashMap<String, Vertex>();
			
			for (Integer i = 0; i < N; i++) {
				HashSet<Vertex> initialVertices = new HashSet<Vertex>();
				this.partToVertexMap.put(i, initialVertices);
			}
			
			this.stringBuilderList = new ArrayList<String>();
			this.stringBuilderList.add("        Alice                   Bob                       Counts\n");
			
			this.s = new StringBuilder();
			this.s.append("        Alice                   Bob                Counts\n");
			this.scoreString = new ArrayList<String>();
			this.movesString = new ArrayList<String>();
			
			this.newVert = false;
		}
		
		public PartitionedGraphRestructured (PartitionedGraphRestructured g) {
			// Copy over all the information
			this.degree = g.degree;
			this.measureScore = g.measureScore;
			
			this.partToVertexMap = createPartMapClone(g.partToVertexMap); 
			this.vertexToEdgeMap = createVertMapClone(g.vertexToEdgeMap);
			this.stringToVertexMap = createStringMapClone(g.stringToVertexMap);

			this.edges = new HashSet<Edge>();
			this.edges.addAll(g.edges);
			

			this.stringBuilderList = new ArrayList<String>();
			this.stringBuilderList.addAll(g.stringBuilderList);
			
			this.s = new StringBuilder();
			s.append(g.stringRepresentation());
		}
		
		
		public int scoreCandidateEdge(Edge e) {
			
			HashMap<Vertex, HashSet<Edge>> vertMapClone = createVertMapClone(this.vertexToEdgeMap);
			HashSet<Vertex> vertices = new HashSet<Vertex>();
			for (int i = 0; i < e.values.length; i++) {
				vertices.add(e.values[i]);
			}
			
			return scoreCandidateHelper(vertices, vertMapClone, 0);
		}
		
		public int scoreCandidateHelper(HashSet<Vertex> remainingV, HashMap<Vertex, HashSet<Edge>> currMapClone, int currScore) {
			if (remainingV.isEmpty()) return currScore;
			
			int maxAdd = 0;
			
			for (Vertex v: remainingV) {
				HashSet<Vertex> newRemainingV = new HashSet<Vertex>();
				newRemainingV.addAll(remainingV);
				newRemainingV.remove(v);
				
				HashSet<Edge> edgesWhichAreNowHit = new HashSet<Edge>();
				HashSet<Edge> ptr = currMapClone.get(v);
				
				int vScore = 0;

				HashMap<Vertex, HashSet<Edge>> newMapClone = createVertMapClone(currMapClone);
				
				if (ptr != null && ! ptr.isEmpty()) {
					edgesWhichAreNowHit.addAll(ptr);
					for (Edge hit: edgesWhichAreNowHit) {
						vScore = vScore + 1;
						Vertex[] containedVertices = hit.values;
						for (Vertex v2: containedVertices) {
							HashSet<Edge> tmp = newMapClone.get(v2);
							if (tmp != null) {
								tmp.remove(hit);
								newMapClone.put(v2, tmp);
							}
						}
					}
				}
				
				int currAdd = currScore + scoreCandidateHelper(newRemainingV, newMapClone, vScore);
				if (currAdd > maxAdd) {
					maxAdd = currAdd;
				}
			}
			return maxAdd;
		}
		
		public boolean candidateEdgeIsGreedy(int M, int N, Edge e, Vertex[] vc) {
			HashSet<Edge> greedyMoves = generateAllGreedyMoves(M, N, vc);
			//System.out.println(greedyMoves);
			
			for (Edge ej: greedyMoves) {
				//System.out.println(ej.toString() + ", " + e.toString());
				//System.out.println(ej.hashCode() + ", " + e.hashCode());
				if (ej.hashCode() == e.hashCode()) {
					return true;
				}
			}
			return false;
			
			//Greedy chooses the best edge for each iteration.
//			HashSet<Vertex> verticesRemaining = new HashSet<Vertex>(); 
//			for (Vertex edgeValue: e.values) {
//				verticesRemaining.add(edgeValue);
//			}
//			
//			HashMap<Vertex, HashSet<Edge>> remainingVertexToEdgeMap = createVertMapClone(this.vertexToEdgeMap);
//			// Remove the vertex cover from the possible vertex selections
//			for (int i = 0; i < vc.length; i++) {
//				Vertex vertex = vc[i];
//				//decrease the count of vertices in the partition
//				remainingVertexToEdgeMap.remove(vertex);
//			}
//			
//			for (Integer p: this.partToVertexMap.keySet()) {
//				
//				Integer size = this.partToVertexMap.get(p).size();
//				Vertex newVertex = new Vertex(PartName.valueOf(p), size);
//				HashSet<Edge> emptyEdgeSet = new HashSet<Edge>();
//				
//				remainingVertexToEdgeMap.put(newVertex, emptyEdgeSet);
//				
//			}
//			return recursiveGreedyCandidateHelper(verticesRemaining, remainingVertexToEdgeMap);
		}
		
		public boolean recursiveGreedyCandidateHelper(HashSet<Vertex> verticesRemaining, HashMap<Vertex, HashSet<Edge>> remainingVertexToEdgeMap) {
			if (verticesRemaining.isEmpty()) {
				return true;
			} else {
				
				HashSet<Vertex> maxCountVertices = new HashSet<Vertex>();
				
				//first pass get the maxsize
				int maxsize = 0;
				for (HashSet<Edge> edgeSet:remainingVertexToEdgeMap.values()) {
					if (edgeSet.size() > maxsize) {
						maxsize=edgeSet.size();
					}
				}
				
				if (maxsize==0) return true;
				
				//second pass get all vertices corresponding to maxsize
				for (Entry<Vertex, HashSet<Edge>> entry:remainingVertexToEdgeMap.entrySet()) {
					int entrySize = entry.getValue().size();
					if (entrySize==maxsize) {
						maxCountVertices.add(entry.getKey());
					}
				}
				
				boolean correctness = false;
				for (Vertex v: verticesRemaining) {
					if (maxCountVertices.contains(v)) {
						HashSet<Vertex> newVerticesRemaining = new HashSet<Vertex>();
						newVerticesRemaining.addAll(verticesRemaining);
						newVerticesRemaining.remove(v);
												
						// update the remainingVertexToEdgeMap with a new thing for each one 
						HashMap<Vertex, HashSet<Edge>> newRemainingVertexToEdgeMap = createVertMapClone(remainingVertexToEdgeMap);
						
						HashSet<Edge> edgesWhichAreNowHit = new HashSet<Edge>();
						HashSet<Edge> ptr = newRemainingVertexToEdgeMap.get(v);
						
						if (ptr != null) {
							edgesWhichAreNowHit.addAll(ptr);
							for (Edge hit: edgesWhichAreNowHit) {
								Vertex[] containedVertices = hit.values;
								for (Vertex v2: containedVertices) {
									if (remainingVertexToEdgeMap.containsKey(v2)) {
										HashSet<Edge> tmp = newRemainingVertexToEdgeMap.get(v2);
										tmp.remove(hit);
										newRemainingVertexToEdgeMap.put(v2, tmp);
									}
								}
							}
						}
						HashSet<Vertex> verticesInThisPart = this.partToVertexMap.get(v.part.value);
						for (Vertex v3: verticesInThisPart) {
							newRemainingVertexToEdgeMap.remove(v3);
							newVerticesRemaining.remove(v3);
						}
						
						if (recursiveGreedyCandidateHelper(newVerticesRemaining, newRemainingVertexToEdgeMap)) {
							correctness = true;
						}
					}
				}
				return correctness;
			}
		}
		
		public void addEdge(Edge e, boolean isGreedy) {
			this.edges.add(e);
			this.depth = this.depth+1;
			this.measureScore = this.measureScore + e.measure;
			boolean flag = false;
			int localScore = 0;
			
			for (Vertex v:e.values) {
				if (this.vertexToEdgeMap.containsKey(v)) {
					HashSet<Edge> vmap = this.vertexToEdgeMap.get(v);
					this.score = this.score + vmap.size();
					localScore = localScore + vmap.size();
					vmap.add(e);
					this.vertexToEdgeMap.put(v, vmap);
					
					if (this.partToVertexMap.containsKey(v.part.value)) {
						HashSet<Vertex> tmp = this.partToVertexMap.get(v.part.value);
						tmp.add(v);
						this.partToVertexMap.put(v.part.value, tmp);
					} else {
						HashSet<Vertex> newPartSet = this.partToVertexMap.get(v.part.value);
						newPartSet.add(v);
						this.partToVertexMap.put(v.part.value, newPartSet);
					}
					
					StringBuilder vertexName = new StringBuilder();
					PartName name = v.part;
					Integer num = v.value;
					vertexName.append(name);
					vertexName.append(num);

					if (!this.stringToVertexMap.containsKey(vertexName.toString())) {
						this.stringToVertexMap.put(vertexName.toString(), v); // This is ok because we're not gonna get too many, so might as well just add duplicates sometimes
					}
				} else {
					flag = true;
					HashSet<Edge> newMap = new HashSet<Edge>();
					newMap.add(e);
					this.vertexToEdgeMap.put(v, newMap);
					
					if (this.partToVertexMap.containsKey(v.part.value)) {
						HashSet<Vertex> tmp = this.partToVertexMap.get(v.part.value);
						tmp.add(v);
						this.partToVertexMap.put(v.part.value, tmp);
					} else {
						HashSet<Vertex> newPartSet = this.partToVertexMap.get(v.part.value);
						newPartSet.add(v);
						this.partToVertexMap.put(v.part.value, newPartSet);
					}
					
					StringBuilder vertexName = new StringBuilder();
					PartName name = v.part;
					Integer num = v.value;
					vertexName.append(name);
					vertexName.append(num);

					if (!this.stringToVertexMap.containsKey(vertexName.toString())) {
						this.stringToVertexMap.put(vertexName.toString(), v); // This is ok because we're not gonna get too many, so might as well just add duplicates sometimes
					}
				}
			}
			
			if (flag == true) {
				this.newVert = true;
			} else {
				this.newVert = false;
			}
			
			String greedyString;
			if (isGreedy) {
				greedyString = "     Greedy \n";
			} else {
				greedyString = "     Not greedy \n";
			}
			String edgeString = e.toString() + greedyString + this.generateCountString(N) + '\n'; //"   local edge score: " + localScore + ", measure score: " + e.measure + ". Edge was greedy: " + isGreedy + ". Current score: " + this.score + ". Measure score: " + this.measureScore +'\n' + this.generateCountString(e.values.length) + '\n';
			this.stringBuilderList.add(edgeString);
			
			this.scoreString.add(Integer.toString(localScore));
			
			s.append(e.toString() + "   local edge score: " + localScore);
			s.append('\n');
			s.append(this.generateCountString(e.values.length));
			s.append('\n');
		}
		
		public void removeEdge(Edge e, boolean add_vc) {
			int stringListLength = this.stringBuilderList.size();
			this.stringBuilderList.remove(stringListLength-1); // Hopefully this removes the edge string.
			if (add_vc) {
				this.stringBuilderList.remove(stringListLength-2); // And this removes the vertex cover string.
				int scoreListLength = this.scoreString.size();
				this.scoreString.remove(scoreListLength-1);
			}
			this.depth = this.depth - 1;
			this.measureScore = this.measureScore - e.measure;
			
			
			
			for (Vertex v:e.values) {
				if (this.vertexToEdgeMap.containsKey(v)) { // should never get to the else
					HashSet<Edge> vmap = this.vertexToEdgeMap.get(v);
					vmap.remove(e);
					this.score = this.score - vmap.size();
					this.vertexToEdgeMap.put(v, vmap);
					if (vmap.isEmpty()) { //if we removed the last edge, then decrease the vertex counts
						if (this.partToVertexMap.containsKey(v.part.value)) {
							HashSet<Vertex> tmp = this.partToVertexMap.get(v.part.value);
							tmp.remove(v);
							this.partToVertexMap.put(v.part.value, tmp);
						}
						StringBuilder vertexName = new StringBuilder();
						PartName name = v.part;
						Integer num = v.value;
						vertexName.append(name);
						vertexName.append(num);
						this.stringToVertexMap.remove(vertexName.toString());
					}
				} else {
					System.out.println("ERROR should never reach here");
				}
			}
			
			this.edges.remove(e);
		}
		
		public String stringRepresentation() {
			return s.toString();
		}
		
		public HashSet<Edge> generateAllMoves(int M, int N, Vertex[] vc) {
			String vertexString = formatLine(this.edges.size(), vc);
			this.s.append(vertexString); 
			this.stringBuilderList.add(vertexString);

			HashSet<Edge> allPossibleMoves = new HashSet<Edge>();
			
			if (vc[0] == null) {
				Edge newEdge = new Edge(N);
				newEdge.createFirstEdge();
				for (Integer i = 0; i < N; i++) {
					HashSet<Vertex> verts = this.partToVertexMap.get(i);
					Vertex v = newEdge.values[i];
					verts.add(v);
					PartName name = PartName.valueOf(i);
					int num = 1;
					StringBuilder vertexName = new StringBuilder();
					vertexName.append(name);
					vertexName.append(num);
					this.stringToVertexMap.put(vertexName.toString(), v);
					this.partToVertexMap.put(i, verts);
				}
				
				newEdge.hitScore = 0;
				allPossibleMoves.add(newEdge);
				return allPossibleMoves;
			}
			
			// Initialize a for-now-empty edge.
			Edge e = new Edge(N);
			HashMap<Integer, HashSet<Vertex>> initialValidPositions = createPartMapClone(this.partToVertexMap); 
			
			HashMap<Vertex, HashSet<Edge>> currentVertexMap =  createVertMapClone(this.vertexToEdgeMap);
			HashSet<Edge> edgesToHitSoFar = new HashSet<Edge>();
			edgesToHitSoFar.addAll(this.edges);
			
			
			// Remove the vertex cover from the possible vertex selections
			for (int i = 0; i < vc.length; i++) {
				Vertex vertex = vc[i];
				//decrease the count of vertices in the partition
				currentVertexMap.remove(vertex);
				
				// remove it from the set of possible vertices we can pick from
				HashSet<Vertex> amendPossibleVerticesForPosition = initialValidPositions.get(vertex.part.value);
				amendPossibleVerticesForPosition.remove(vertex);
				initialValidPositions.put(vertex.part.value, amendPossibleVerticesForPosition);
			}
			
			// Want to generate all the valid moves; ie ones that differ from no more than a single other move.
			addAllPossibleEdges(currentVertexMap, edgesToHitSoFar, e, initialValidPositions, allPossibleMoves);
			
			//System.out.println("Got all possible edges for the current move: " + allPossibleMoves);
			
			if (allPossibleMoves.isEmpty()) {
				this.s.append("No possible move.\n\n");
				//this.stringBuilderList.add("No possible move.\n\n");
//				this.s.append("                      Bob tried to generate: \n                      " + e.toString() + ".\n");
//				this.s.append("                      However, this misses the edges:\n");
//				this.s.append("                    \n");
//				for (Edge remainingEdge: remainingEdgesToHit) {
//					s.append("  ");
//					s.append(remainingEdge.toString());
//					s.append("\n");
//				}
			}
			//System.out.println("Number of possible moves: " + allPossibleMoves.size());
			return allPossibleMoves;
			
		}
			
		public void addAllPossibleEdges(HashMap<Vertex, HashSet<Edge>> currentVertexMap, HashSet<Edge> edgesToHitSoFar, Edge currentEdge, HashMap<Integer, HashSet<Vertex>> remainingValidPositions, HashSet<Edge> allPossibleEdges) {
			currentEdge.hitScore = currentEdge.hitScore + 1;
			
			// Sets which represent what this edge has left to hit
			HashMap<Vertex, HashSet<Edge>> remainingVertexToEdgeMap = createVertMapClone(currentVertexMap);
			HashSet<Edge> remainingEdgesToHit = new HashSet<Edge>();
			remainingEdgesToHit.addAll(edgesToHitSoFar);
						
			HashMap<Integer, HashSet<Vertex>> remainingPositions = createPartMapClone(remainingValidPositions); 
			
			// Must consider all possible orderings of choices of vertices. So, instead of finding min, we find all possible vertices which we could select this round.
			HashSet<Vertex> possibleVerticesToAdd = new HashSet<Vertex>();
			// get all the candidate new vertices (don't worry about making new ones yet)
			for (HashSet<Vertex> val: remainingValidPositions.values()) {
				possibleVerticesToAdd.addAll(val);
			}
//			System.out.println("Current vertex to edge map: " + currentVertexMap);
//			System.out.println("Current vertex counts: " + currentVertexCounts);
//			System.out.println("Edges to hit so far: " + edgesToHitSoFar);
//			System.out.println("Current edge: " + currentEdge);
//			System.out.println("Remaining valid positions: " + remainingValidPositions);
//			System.out.println("All possible edges: " + allPossibleEdges);
//			System.out.println("Possible vertices to add: " + possibleVerticesToAdd);
			
			Vertex newVertex;
			
			// if empty, we need to consider every choice of new vertex to make.
			if (possibleVerticesToAdd.isEmpty() && remainingEdgesToHit.isEmpty()) {
				Edge e = new Edge(currentEdge); // Make a new edge based on current edge so values don't copy over. 

				// make a new vertex for all remaining valid positions
				for (Integer key: remainingValidPositions.keySet()) {
					PartName name = PartName.valueOf(key);
					int num = this.partToVertexMap.get(key).size() + 1; //[key] + 1;

					newVertex = new Vertex(name, num);
					
					StringBuilder vertexName = new StringBuilder();
					vertexName.append(name);
					vertexName.append(num);
					
//					System.out.println("Setting vertex: " + newVertex);
					e.setVertex(newVertex, 0);
					remainingPositions.remove(newVertex.part.value);
				}
				if (remainingPositions.isEmpty()) {
//					System.out.println("Adding edge: " + e);
					if (!allPossibleEdges.contains(e)) {
						boolean flag = true;
						//System.out.println(e.hashCode());
						for (Edge ej: allPossibleEdges) {
							//System.out.println(ej.hashCode());
							if (ej.hashCode() == e.hashCode()) {
								flag = false;	
							}
						}
						if (flag == true) {
							allPossibleEdges.add(e);
						}
					}
				} else {
					//this isn't a valid edge. 
				}
			} else {
				// Go over all the possible vertices to add and recurse with these.
				for (Vertex v: possibleVerticesToAdd) {
					Edge e = new Edge(currentEdge); // Make a new edge based on current edge so values don't copy over. 

					remainingPositions = createPartMapClone(remainingValidPositions);
					remainingPositions.remove(v.part.value);
					
					remainingVertexToEdgeMap = createVertMapClone(currentVertexMap);
					remainingEdgesToHit = new HashSet<Edge>();
					remainingEdgesToHit.addAll(edgesToHitSoFar);
					
					// Remove all edges to hit which contain this vertex, for this iteration
					HashSet<Edge> edgesWhichAreNowHit = new HashSet<Edge>();
					HashSet<Edge> ptr = remainingVertexToEdgeMap.get(v);
					
					int vScore = 0;
					//int hScore = 0;
					
					if (ptr != null && ! ptr.isEmpty()) {
						//hScore = hScore + 1;
						edgesWhichAreNowHit.addAll(ptr);
						for (Edge hit: edgesWhichAreNowHit) {
							remainingEdgesToHit.remove(hit);
							vScore = vScore + 1;
							Vertex[] containedVertices = hit.values;
							for (Vertex v2: containedVertices) {
								HashSet<Edge> tmp = remainingVertexToEdgeMap.get(v2);
								if (tmp != null) {
									tmp.remove(hit);
									remainingVertexToEdgeMap.put(v2, tmp);
								}
							}
						}
					}
					
					e.setVertex(v, vScore);
					
					if (remainingPositions.isEmpty()) {
						if (remainingEdgesToHit.isEmpty()) {
								allPossibleEdges.add(e);
						} // otherwise, this isn't a valid edge.
					} else {
						addAllPossibleEdges(remainingVertexToEdgeMap, remainingEdgesToHit, e,
							remainingPositions, allPossibleEdges);
					}
				}
				
				
			}	
		}
		
		public HashSet<Edge> generateAllGreedyMoves(int M, int N, Vertex[] vc) {			
			HashSet<Edge> allGreedyMoves = new HashSet<Edge>();
			
			if (vc[0] == null) {
				Edge newEdge = new Edge(N);
				newEdge.createFirstEdge();
				for (Integer i = 0; i < N; i++) {
					HashSet<Vertex> verts = this.partToVertexMap.get(i);
					Vertex v = newEdge.values[i];
					verts.add(v);
					PartName name = PartName.valueOf(i);
					int num = 1;
					StringBuilder vertexName = new StringBuilder();
					vertexName.append(name);
					vertexName.append(num);
					this.stringToVertexMap.put(vertexName.toString(), v);
					this.partToVertexMap.put(i, verts);
				}
				
				allGreedyMoves.add(newEdge);
				return allGreedyMoves;
			}
			
			// Initialize a for-now-empty edge.
			Edge e = new Edge(N);
			HashMap<Integer, HashSet<Vertex>> initialValidPositions = createPartMapClone(this.partToVertexMap); 
			
			HashMap<Vertex, HashSet<Edge>> currentVertexMap =  createVertMapClone(this.vertexToEdgeMap);
			HashSet<Edge> edgesToHitSoFar = new HashSet<Edge>();
			edgesToHitSoFar.addAll(this.edges);
			
			// Remove the vertex cover from the possible vertex selections
			for (int i = 0; i < vc.length; i++) {
				Vertex vertex = vc[i];
				
				// remove it from the set of possible vertices we can pick from
				HashSet<Vertex> amendPossibleVerticesForPosition = initialValidPositions.get(vertex.part.value);
				amendPossibleVerticesForPosition.remove(vertex);
				initialValidPositions.put(vertex.part.value, amendPossibleVerticesForPosition);
				
				//decrease the count of vertices in the partition
				currentVertexMap.remove(vertex);
			}
			//System.out.println(formatLine(this.edges.size(), vc));
			HashMap<Integer, HashSet<Vertex>> validPositions = createPartMapClone(initialValidPositions);
			//System.out.println("Valid positions: " + validPositions);
			//System.out.println(generateCountString(this.degree));

			// Want to generate all the valid moves; ie ones that differ from no more than a single other move.
			addAllGreedyEdges(currentVertexMap, edgesToHitSoFar, e, validPositions, allGreedyMoves);
			
			//System.out.println("Got all possible edges for the current move: " + allPossibleMoves);
			
			if (allGreedyMoves.isEmpty()) {
				this.s.append("No possible move.\n\n");
			}
			//System.out.println("Number of possible moves: " + allGreedyMoves.size());
			return allGreedyMoves;
			
		}
		
		public void addAllGreedyEdges(HashMap<Vertex, HashSet<Edge>> currentVertexMap, HashSet<Edge> edgesToHitSoFar, Edge currentEdge, HashMap<Integer, HashSet<Vertex>> remainingValidPositions, HashSet<Edge> allPossibleEdges) {
			
			// Sets which represent what this edge has left to hit
			HashMap<Vertex, HashSet<Edge>> remainingVertexToEdgeMap = createVertMapClone(currentVertexMap);
			HashSet<Edge> remainingEdgesToHit = new HashSet<Edge>();
			remainingEdgesToHit.addAll(edgesToHitSoFar);
			HashMap<Integer, HashSet<Vertex>> remainingPositions = createPartMapClone(remainingValidPositions); 
			
			// must consider all possible vertices which are greedy
			HashSet<Vertex> possibleVerticesToAdd = new HashSet<Vertex>();
				
			//first pass get the maxsize
			int maxsize = 0;
			for (Entry<Vertex,HashSet<Edge>> entry:remainingVertexToEdgeMap.entrySet()) {
				HashSet<Edge> edgeSet= entry.getValue();
				if (edgeSet.size() > maxsize && remainingPositions.containsKey(entry.getKey().part.value)) {
					maxsize=edgeSet.size();
				}
			}
			
			if (maxsize!=0) {
				for (Entry<Vertex, HashSet<Edge>> entry:remainingVertexToEdgeMap.entrySet()) {
					int entrySize = entry.getValue().size();
					if (entrySize==maxsize && remainingPositions.containsKey(entry.getKey().part.value)) {
						possibleVerticesToAdd.add(entry.getKey());
					}
				}
			} else {
				for (HashSet<Vertex> val: remainingValidPositions.values()) {
					possibleVerticesToAdd.addAll(val);
				}
			}
			
			Vertex newVertex;
			
			// if empty, we need to consider every choice of new vertex to make.
			if (possibleVerticesToAdd.isEmpty() && remainingEdgesToHit.isEmpty()) {
				Edge e = new Edge(currentEdge); // Make a new edge based on current edge so values don't copy over. 

				// make a new vertex for all remaining valid positions
				for (Integer key: remainingPositions.keySet()) {
					PartName name = PartName.valueOf(key);
					int num = this.partToVertexMap.get(key).size() + 1; //[key] + 1;
					System.out.println("NUM: " +num);


					newVertex = new Vertex(name, num);
					
					StringBuilder vertexName = new StringBuilder();
					vertexName.append(name);
					vertexName.append(num);
					
//					System.out.println("Setting vertex: " + newVertex);
					e.setVertex(newVertex, 0);
					remainingPositions.remove(newVertex.part.value);
					for (Vertex v: this.partToVertexMap.get(key)) {
						possibleVerticesToAdd.remove(v);
						//remainingVertexToEdgeMap.remove(v);
					}
				}
				if (remainingPositions.isEmpty()) {
					if (!allPossibleEdges.contains(e)) {
						boolean flag = true;
						//System.out.println(e.hashCode());
						for (Edge ej: allPossibleEdges) {
							//System.out.println(ej.hashCode());
							if (ej.hashCode() == e.hashCode()) {
								flag = false;	
							}
						}
						if (flag == true) {
							allPossibleEdges.add(e);
						}
					}
				} else {
					//this isn't a valid edge. 
				}
			} else {
				// Go over all the possible vertices to add and recurse with these.
				for (Vertex v: possibleVerticesToAdd) {
					Edge e = new Edge(currentEdge); // Make a new edge based on current edge so values don't copy over. 

					remainingPositions = createPartMapClone(remainingValidPositions);
					remainingPositions.remove(v.part.value);
					
					remainingVertexToEdgeMap = createVertMapClone(currentVertexMap);
					
					remainingEdgesToHit = new HashSet<Edge>();
					remainingEdgesToHit.addAll(edgesToHitSoFar);
					
					// Remove all edges to hit which contain this vertex, for this iteration
					HashSet<Edge> edgesWhichAreNowHit = new HashSet<Edge>();
					HashSet<Edge> ptr = remainingVertexToEdgeMap.get(v);
					
					int vScore = 0;
					
					if (ptr != null && ! ptr.isEmpty()) {
						edgesWhichAreNowHit.addAll(ptr);
						for (Edge hit: edgesWhichAreNowHit) {
							remainingEdgesToHit.remove(hit);				
							vScore = vScore + 1;
							Vertex[] containedVertices = hit.values;
							for (Vertex v2: containedVertices) {
								HashSet<Edge> tmp = remainingVertexToEdgeMap.get(v2);
								if (tmp != null) {
									tmp.remove(hit);
									remainingVertexToEdgeMap.put(v2, tmp);
								}
							}
						}
					}
					e.setVertex(v, vScore);
					if (remainingPositions.isEmpty()) {
						if (remainingEdgesToHit.isEmpty()) {
								allPossibleEdges.add(e);
						} // otherwise, this isn't a valid edge.
					} else {
						addAllPossibleEdges(remainingVertexToEdgeMap, remainingEdgesToHit, e,
							remainingPositions, allPossibleEdges);
					}
				}
			}	
		}
		
		public HashSet<Edge> generateOverlapMoves(int N, boolean firstEdge, boolean add_vc, Vertex[] vc) {
			String vertexString = formatLine(this.edges.size(), vc);
			this.s.append(vertexString); 
			this.stringBuilderList.add(vertexString);
			
			HashSet<Edge> allPossibleMoves = new HashSet<Edge>();
			
			if (firstEdge) {
				Edge newEdge = new Edge(N);
				newEdge.createFirstEdge();
				for (Integer i = 0; i < N; i++) {
					HashSet<Vertex> verts = this.partToVertexMap.get(i);
					Vertex v = newEdge.values[i];
					verts.add(v);
					PartName name = PartName.valueOf(i);
					int num = 1;
					StringBuilder vertexName = new StringBuilder();
					vertexName.append(name);
					vertexName.append(num);
					this.stringToVertexMap.put(vertexName.toString(), v);
					this.partToVertexMap.put(i, verts);
				}
				
				allPossibleMoves.add(newEdge);
				return allPossibleMoves;
			}

			// Initialize a for-now-empty edge.
			Edge e = new Edge(N);
			HashMap<Integer, HashSet<Vertex>> initialValidPositions = createPartMapClone(this.partToVertexMap); 
			
			HashMap<Vertex, HashSet<Edge>> currentVertexMap =  createVertMapClone(this.vertexToEdgeMap);
			HashSet<Edge> edgesToHitSoFar = new HashSet<Edge>();
			edgesToHitSoFar.addAll(this.edges);
			
			//If we are playing responsively, remove the vertex cover from the possible vertex selections
			for (int i = 0; i < vc.length; i++) {
				Vertex vertex = vc[i];
				
				// remove it from the set of possible vertices we can pick from
				HashSet<Vertex> amendPossibleVerticesForPosition = initialValidPositions.get(vertex.part.value);
				amendPossibleVerticesForPosition.remove(vertex);
				initialValidPositions.put(vertex.part.value, amendPossibleVerticesForPosition);
				
				//decrease the count of vertices in the partition
				currentVertexMap.remove(vertex);
			}

			// Want to generate all the valid moves; ie ones that differ from no more than a single other move.
			addAllOverlapEdges(N-1, currentVertexMap, edgesToHitSoFar, e, initialValidPositions, allPossibleMoves);

			if (allPossibleMoves.isEmpty()) {
				this.s.append("No possible move.\n\n");
			}
			//System.out.println("Number of possible moves: " + allPossibleMoves.size());
			return allPossibleMoves;
		}
		
		public void addAllOverlapEdges(int maxDegree, HashMap<Vertex, HashSet<Edge>> currentVertexMap, HashSet<Edge> edgesToHitSoFar, Edge currentEdge, HashMap<Integer, HashSet<Vertex>> remainingValidPositions, HashSet<Edge> allPossibleEdges) {
			
			// Sets which represent what this edge has left to hit
			HashMap<Vertex, HashSet<Edge>> remainingVertexToEdgeMap = createVertMapClone(currentVertexMap);
			HashSet<Edge> remainingEdgesToHit = new HashSet<Edge>();
			remainingEdgesToHit.addAll(edgesToHitSoFar);

			HashMap<Integer, HashSet<Vertex>> remainingPositions = createPartMapClone(remainingValidPositions); 

			HashSet<Vertex> possibleVerticesToAdd = new HashSet<Vertex>();
			//modify vertices in the part left to hit
			for (Entry<Integer, HashSet<Vertex>> val: remainingValidPositions.entrySet()) {
				HashSet<Vertex> tmp = new HashSet<Vertex>();
				for (Vertex v: val.getValue()) {
					if (v.value <= maxDegree) {
						possibleVerticesToAdd.add(v);
						tmp.add(v);
					} 
				}
				remainingPositions.put(val.getKey(), tmp);
				// Now remove all the
			}
//									System.out.println("Current vertex to edge map: " + currentVertexMap);
//									//System.out.println("Current vertex counts: " + currentVertexCounts);
//									System.out.println("Edges to hit so far: " + edgesToHitSoFar);
//									System.out.println("Current edge: " + currentEdge);
//									System.out.println("Remaining valid positions: " + remainingValidPositions);
//									System.out.println("All possible edges: " + allPossibleEdges);
//									System.out.println("Possible vertices to add: " + possibleVerticesToAdd);

			Vertex newVertex;

			// if empty, we need to consider every choice of new vertex to make.
			if (possibleVerticesToAdd.isEmpty() && remainingEdgesToHit.isEmpty()) {
				Edge e = new Edge(currentEdge); // Make a new edge based on current edge so values don't copy over. 
				boolean fail = false;
				// make a new vertex for all remaining valid positions
				for (Integer key: remainingValidPositions.keySet()) {
					PartName name = PartName.valueOf(key);
					int num = this.partToVertexMap.get(key).size() + 1; 
					
					//System.out.println("num is " + num + " and max degree is " + maxDegree);
					
					if (num >= maxDegree) {
						fail = true;
					} else {
						newVertex = new Vertex(name, num);
						
						StringBuilder vertexName = new StringBuilder();
						vertexName.append(name);
						vertexName.append(num);
	
						e.setVertex(newVertex, 0);
						remainingPositions.remove(newVertex.part.value);
					}

				}
				if (remainingPositions.isEmpty()) {
					if (!allPossibleEdges.contains(e)) {
						boolean flag = true;
						//System.out.println(e.hashCode());
						for (Edge ej: allPossibleEdges) {
							//System.out.println(ej.hashCode());
							if (ej.hashCode() == e.hashCode()) {
								flag = false;	
							}
						}
						if (flag == true && (fail == false)) {
							//System.out.println("Adding candidate: " + e);
							allPossibleEdges.add(e);
						}
					}
				} else {
					//this isn't a valid edge. 
				}
			} else {
				// Go over all the possible vertices to add and recurse with these.
				for (Vertex v: possibleVerticesToAdd) {
					//System.out.println("Possible vertex to add: " + v);
					Edge e = new Edge(currentEdge); // Make a new edge based on current edge so values don't copy over. 

					remainingPositions = createPartMapClone(remainingValidPositions);
					remainingPositions.remove(v.part.value);
					
					remainingVertexToEdgeMap = createVertMapClone(currentVertexMap);
					remainingEdgesToHit = new HashSet<Edge>();
					remainingEdgesToHit.addAll(edgesToHitSoFar);

					// Remove all edges to hit which contain this vertex, for this iteration
					HashSet<Edge> edgesWhichAreNowHit = new HashSet<Edge>();
					HashSet<Edge> ptr = remainingVertexToEdgeMap.get(v);

					int vScore = 0;

					if (ptr != null && ! ptr.isEmpty()) {
						edgesWhichAreNowHit.addAll(ptr);
						for (Edge hit: edgesWhichAreNowHit) {
							remainingEdgesToHit.remove(hit);
							vScore = vScore + 1;
							Vertex[] containedVertices = hit.values;
							for (Vertex v2: containedVertices) {
								HashSet<Edge> tmp = remainingVertexToEdgeMap.get(v2);
								if (tmp != null) {
									tmp.remove(hit);
									remainingVertexToEdgeMap.put(v2, tmp);
								}
								int vpart = v2.part.value;
								HashSet<Vertex> tmp2 = remainingPositions.get(vpart);
								
								//System.out.println("tmp2 is " + tmp2);
								if (tmp2 != null) {
									tmp2.remove(v2);
									//System.out.println("Removed " + v2 + " because it was in some edge incident on " + v + " so now tmp2 is " + tmp2);
									remainingPositions.put(v2.part.value, tmp2);
								}
							}
						}
					}

					//System.out.println("Setting vertex in edge: " + v);
					e.setVertex(v, vScore);

					if (remainingPositions.isEmpty()) {
						if (remainingEdgesToHit.isEmpty()) {
							//System.out.println("adding candidate: " + e);
							allPossibleEdges.add(e);
						} // otherwise, this isn't a valid edge.
					} else {
						//System.out.println("Remaining positions was: " + remainingPositions);
						addAllPossibleEdges(remainingVertexToEdgeMap, remainingEdgesToHit, e,
								remainingPositions, allPossibleEdges);
					}
				}


			}	
		}

		
		
		public Edge generateMove(int M, int N, Vertex[] vc, boolean add_vc) {
			if (add_vc == true) {
				String vertexCoverString = formatLine(this.edges.size(), vc);
				this.s.append(vertexCoverString); 
				this.stringBuilderList.add(vertexCoverString);
			}
			
			if (vc[0] == null) {
				Edge newEdge = new Edge(N);
				newEdge.createFirstEdge();
				for (Integer i = 0; i < N; i++) {
					HashSet<Vertex> verts = this.partToVertexMap.get(i);
					Vertex v = newEdge.values[i];
					verts.add(v);
					PartName name = PartName.valueOf(i);
					int num = 1;
					StringBuilder vertexName = new StringBuilder();
					vertexName.append(name);
					vertexName.append(num);
					this.stringToVertexMap.put(vertexName.toString(), v);
					this.partToVertexMap.put(i, verts);
				}
				
				s.append(newEdge);
				return newEdge;
			}
			
			// Sets which represent what this edge has left to hit
			HashMap<Vertex, HashSet<Edge>> remainingVertexToEdgeMap = createVertMapClone(this.vertexToEdgeMap);
			HashSet<Edge> remainingEdgesToHit = new HashSet<Edge>();
			remainingEdgesToHit.addAll(this.edges);
			
			// Initialize a for-now-empty edge.
			Edge e = new Edge(N);
			
			HashMap<Integer, HashSet<Vertex>> initialValidPositions = createPartMapClone(this.partToVertexMap); 

			// Remove the vertex cover vertices from possibilities to choose from.
			for (int i = 0; i < vc.length; i++) {
				Vertex vertex = vc[i];
				//decrease the count of vertices in the partition
				remainingVertexToEdgeMap.remove(vertex);
				
				// remove it from the set of possible vertices we can pick from
				HashSet<Vertex> amendPossibleVerticesForPosition = initialValidPositions.get(vertex.part.value);
				amendPossibleVerticesForPosition.remove(vertex);
				initialValidPositions.put(vertex.part.value, amendPossibleVerticesForPosition);
			}
			HashMap<Integer, HashSet<Vertex>> remainingPositions = createPartMapClone(initialValidPositions); 

			// Fill up the edge until there are no positions to fill anymore.
			int numTimesToFill = 0;
			while (remainingPositions.isEmpty() == false) {
				numTimesToFill = numTimesToFill + 1;
				// Get the minimum sized part at this time
				int min = 100000;
				int partChoice = -1;
				
				// Get the sizes 
				for (Entry<Integer, HashSet<Vertex>> ent: remainingPositions.entrySet()) {
					int size = ent.getValue().size();
					if (size > 0 && size < min) {
						min = size;
						partChoice = ent.getKey(); 
					} 
				}
				if (partChoice == -1) {
					int minpart = N;
					for (Entry<Integer, HashSet<Vertex>> ent: remainingPositions.entrySet()) {
						int size = ent.getValue().size();
						int part = ent.getKey();
						if (size == 0 && part < minpart) {
							minpart = part;
							partChoice = part;
						}
					}
				}
				
				//Find the maximum vertex for that part 
				// Choose the vertex from this set that has the maximum size
				int max = 0;
				Vertex bestVertex = null;
				HashSet<Vertex> choices = remainingPositions.get(partChoice);
				// Edge cases: vertex cover covers all the vertices in this part, or there are no more vertices in remainingPositions
				if (choices == null) {
					PartName name = PartName.valueOf(partChoice);
					int num = this.partToVertexMap.get(name.value).size() + 1;
					bestVertex = new Vertex(name, num);
					
					if(add_vc == false) {
						//this.vertexCounts[partChoice] = num;
						StringBuilder vertexName = new StringBuilder();
						vertexName.append(name);
						vertexName.append(num);
						this.stringToVertexMap.put(vertexName.toString(), bestVertex);
						this.vertexToEdgeMap.put(bestVertex, new HashSet<Edge>());
					
						if (this.partToVertexMap.containsKey(partChoice)) {
							HashSet<Vertex> vertexChoices = this.partToVertexMap.get(partChoice);
							vertexChoices.add(bestVertex);
							this.partToVertexMap.put(partChoice, vertexChoices);
						} else {
							HashSet<Vertex> vertexChoices = new HashSet<Vertex>();
							vertexChoices.add(bestVertex);
							this.partToVertexMap.put(partChoice, vertexChoices);
						}
					}
				}
				else if (choices.isEmpty()) {
					if (initialValidPositions.get(partChoice).isEmpty()) {
						//make a new vertex
						PartName name = PartName.valueOf(partChoice);
						int num = this.partToVertexMap.get(partChoice).size() + 1; //vertexCounts[partChoice] + 1;
						bestVertex = new Vertex(name, num);
						
						if (add_vc == false) {
							//this.vertexCounts[partChoice] = num;
							StringBuilder vertexName = new StringBuilder();
							vertexName.append(name);
							vertexName.append(num);
							this.stringToVertexMap.put(vertexName.toString(), bestVertex);
							this.vertexToEdgeMap.put(bestVertex, new HashSet<Edge>());
							
							if (this.partToVertexMap.containsKey(partChoice)) {
								HashSet<Vertex> vertexChoices = this.partToVertexMap.get(partChoice);
								vertexChoices.add(bestVertex);
								this.partToVertexMap.put(partChoice, vertexChoices);
							} else {
								HashSet<Vertex> vertexChoices = new HashSet<Vertex>();
								vertexChoices.add(bestVertex);
								this.partToVertexMap.put(partChoice, vertexChoices);
							}
						}
					} else {
						// OLD CODE: arbitrarily choose a vertex from the initial set
						// bestVertex = initialValidPositions.get(partChoice).iterator().next();
						
						// NEW code: choose 'greedily' (ie the lowest entry available for this part)
						HashSet<Vertex> verticesInPart = initialValidPositions.get(partChoice);
						int lowestVal = 0;
						for (Vertex v: verticesInPart) {
							if (v.value < lowestVal) {
								lowestVal = v.value;
								bestVertex = v;
							}
						}		
					}
				} else {
					int minpart = -1;
					int minval = -1;
					for (Vertex v : choices) {
						HashSet<Edge> edgeCover = remainingVertexToEdgeMap.get(v);
						int size = edgeCover.size();
						int part = v.part.value;
						int val = v.value;
						if (size > max) {
							bestVertex = v;
							minpart = part;
							minval = val;
							max = size;
						} else if (size == max) {
							if (minpart == -1 || part < minpart) {
								minpart = part;
								bestVertex = v;
							} else if (part == minpart) {
								if (minval == -1 || val < minval) {
									minval = val;
									bestVertex = v;
									
								}
							}
						}
					}
				}

				// Fix the vertex in the new edge
				e.setVertex(bestVertex, max);
				
				// Remove this part from the parts remaining and reset its count so it doesn't get chosen again
				remainingPositions.remove(partChoice);	
				//remainingVertexCounts[partChoice] = -1;
				
				// Remove all edges to hit which contain this vertex
				HashSet<Edge> edgesWhichAreNowHit = new HashSet<Edge>();
				HashSet<Edge> ptr = remainingVertexToEdgeMap.get(bestVertex);
				
				int vScore = 0;
				
				if (ptr != null && !ptr.isEmpty()) {
					edgesWhichAreNowHit.addAll(ptr);
					for (Edge hit: edgesWhichAreNowHit) {
						remainingEdgesToHit.remove(hit);
						vScore = vScore + 1;
						Vertex[] containedVertices = hit.values;
						for (Vertex v: containedVertices) {
							HashSet<Edge> tmp = remainingVertexToEdgeMap.get(v);
							if (tmp != null) {
								tmp.remove(hit);
								remainingVertexToEdgeMap.put(v, tmp);
							}
						}
					}
				}
				
				if (vScore != max) System.err.print("Error the vertex score and max vertex score should be same");
			}
			
			if (remainingEdgesToHit.isEmpty()) {
				e.hitScore = numTimesToFill;
				return e;
			} else {
				System.out.println("No possible move.");
				System.out.println("                      Bob tried to generate: \n                      " + e.toString() + ".");
				System.out.println("                      However, this misses the edges:");
				System.out.print("                    ");
				for (Edge remainingEdge: remainingEdgesToHit) {
					System.out.println("  " + remainingEdge.toString());
				}
				// if move not possible, return null
				if (add_vc == true) {
					this.s.append("No possible move.");
					this.s.append("                      Bob tried to generate: \n                      " + e.toString() + ".\n");
					this.s.append("                      However, this misses the edges:\n");
					this.s.append("                    \n");
					for (Edge remainingEdge: remainingEdgesToHit) {
						s.append("  ");
						s.append(remainingEdge.toString());
						s.append("\n");
					}
				}
				return null;
			}
			
		}
		
		
		/// Bob's Algorithm
		// Given N-1 vertices and his own set of old vertices
		// He must determine an edge comprised of N integer selections, specifying the edge across partitions, which can differ from at most M-1 of the prior vertices. 
		// Bob's process:
		// - determine which edge in vc corresponds to the smallest of his own partition sets. 
		// - determine what con
		
		/*public boolean generatePossibleMoves(int M, String[] vertexCoverCandidate) {
			boolean movePossible = true;
			int falsities = 0; // number of edges we miss with the current edge
			while (movePossible = true) {
				Edge candidateEdge = new Edge(this.degree);
				
				//first check: if vc covers all of one of L_i, Bob must generate a new vertex in that position. 
				
				for (String vertex: vertexCoverCandidate) {
					char[] vert = vertex.toCharArray();
					PartName vertLetter = PartName.valueOf(Character.toString(vert[0]));
					Integer vertNumber = Integer.parseInt(Character.toString(vert[1]));
				}	
			}
		
			
		}*/
		
		public String generateCountString(int N) {
			StringBuilder s = new StringBuilder();
			s.append("                                                        {");
			
			// get the count in order ie u's, v's,... 
			for (int i = 0; i < N; i++) {
				int p = this.partToVertexMap.get(i).size();
				for (int j = 1; j <= p; j++) {
					StringBuilder vertexName = new StringBuilder();
					vertexName.append(PartName.valueOf(i));
					vertexName.append(j);
					String lookup = vertexName.toString();
					s.append(lookup);
					s.append(": ");
					Vertex v = this.stringToVertexMap.get(lookup);
					Integer count = this.vertexToEdgeMap.get(v).size();
					s.append(count);
					if (j != (p)) {
						s.append(", ");
					}
				}
				if (i != (N-1)) {
					s.append('\n');
					s.append("                                                          ");
				}
			}
			s.append("}");
			return s.toString();
		}
		
		public String generateCountString(int N, Integer[] vertexCounts, HashMap<Vertex, HashSet<Edge>> vertexToEdgeMap) {
			StringBuilder s = new StringBuilder();
			s.append("                                       {");
			
			// get the count in order ie u's, v's,... 
			for (int i = 0; i < N; i++) {
				for (int j = 1; j <= vertexCounts[i]; j++) {
					StringBuilder vertexName = new StringBuilder();
					vertexName.append(PartName.valueOf(i));
					vertexName.append(j);
					String lookup = vertexName.toString();
					s.append(lookup);
					s.append(": ");
					Vertex v = this.stringToVertexMap.get(lookup);
					Integer count = vertexToEdgeMap.get(v).size();
					s.append(count);
					s.append(", ");
				}
				s.append('\n');
				s.append("                                       ");
			}
			s.append("}");
			return s.toString();
		}
		
	}
	
	enum PartName { 
		v1_(0), v2_(1), v3_(2), v4_(3), v5_(4), v6_(5), v7_(6), v8_(7), v9_(8), v10_(9), v11_(10), v12_(11), v13_(12), v14_(13), v15_(14), v16_(15), v17_(16), v18_(17), v19_(18), v20_(19), v21_(20), v22_(21), v23_(22), v24_(23), v25_(24), v26_(25);
		
		private int value;
		private static HashMap map = new HashMap<>();

		private PartName(int value) {
			this.value = value;
		}

		static {
			for (PartName partNum : PartName.values()) {
				map.put(partNum.value, partNum);
			}
		}

		public static PartName valueOf(int partNum) {
			return (PartName) map.get(partNum);
		}

		public int getValue() {
			return value;
		}
	}
	
	class Vertex {
		public PartName part;
		public int value;
		
		public Vertex(String s) {
			char[] vert = s.toCharArray();
			if (vert.length == 0) return;
			PartName vertLetter = PartName.valueOf(Character.toString(vert[0]));
			Integer vertNumber = Integer.parseInt(Character.toString(vert[1]));
			this.part = vertLetter;
			this.value = vertNumber;
		}
		
		public Vertex(PartName part, Integer value) {
			this.part = part;
			this.value = value;
			
			//System.out.println("Creating new vertex " + this.toString());
		}
		
		
		public String toString() {
			StringBuilder s = new StringBuilder();
			s.append(this.part);
			s.append(this.value);
			return s.toString();
		}
		
		@Override
		public boolean equals(Object v2) {
			if (this == v2)
	            return true;
	        if (v2 == null)
	            return false;
	        if (getClass() != v2.getClass())
	            return false;
	        Vertex vert2 = (Vertex) v2;
	        boolean eq = (this.part == vert2.part && this.value == vert2.value);
			return eq;
		}
		
		@Override
	    public int hashCode() {
	        return this.toString().hashCode();
	    }
	}

	class Edge {
		public int degree;
		public Vertex[] values;
		
		public int numSet;
		public int hitScore;
		public int measure;
		
		public Edge(int N) {
			this.degree = N;
			this.values = new Vertex[N];
			
			this.measure = 0;
			this.hitScore = 0;
		}
		
		public void createFirstEdge() {
			for (int i = 0; i < this.degree; i++) {
				this.values[i] = new Vertex(PartName.valueOf(i), 1);
				this.numSet = this.degree;
			}
			
			this.measure = 0;
			this.hitScore = 0;
		}
		
		// Create a copy of an edge from an existing edge
		public Edge(Edge e) {
			this.degree = e.degree;
			this.values = new Vertex[this.degree];
			this.measure = e.measure;
			this.hitScore = e.hitScore;
		
			for (int i=0; i < e.degree; i++) {
				if (e.values[i] != null) {
					this.values[i] = e.values[i];
				}
			}
		}
		
		@Override
		public boolean equals(Object obj) {
	        Edge other = (Edge) obj;
	        if (this.degree != other.degree) return false;

	        for (int i = 0; i < other.degree; i++) {
	        	if (!this.values[i].equals(other.values[i])) {
	        		return false;
	        	}
	        }
	        if (this.measure != other.measure) {
	        	System.err.print("The edges are the same besides the measures: " + this.measure + ", " + other.measure);
	        }
	        return true;
	    }

		@Override
	    public int hashCode() {
	        return this.toString().hashCode();
	    }
		
		/*public void setVertex(PartName letter, Integer number) {
			Integer tmp = letter.value;
			this.values[tmp] = new Vertex(PartName.valueOf(tmp), number);
			this.numSet = this.numSet+1;
		}*/
		
		public void setVertex(Vertex v, int vScore) {
			this.measure = this.measure + vScore;
			int tmp = v.part.value;
			this.values[tmp] = v;
			this.numSet = this.numSet + 1;
			//System.out.println("Setting vertex " + v);
		}
		
		public boolean isComplete() {
			return (this.numSet == this.degree);
		}
		
		public String toString() {
			StringBuilder s = new StringBuilder();
			
			s.append('(');
			int i = 0;
			for (i = 0; i < this.degree - 1; i++) {
				if (this.values[i] == null) {
					s.append("null");
				} else {
					s.append(this.values[i].toString());
				}
				s.append(", ");
			}
			s.append(this.values[i]);
			s.append(')');
			
			return s.toString();
		}
	}
	
	public class Counter {
		private int count;
		
		public Counter() {
			this.count = 0;
		}
		
		public void increment() {
			this.count = this.count+1;
		}
		
	}
	
	public static HashMap<Vertex, HashSet<Edge>> createVertMapClone(HashMap<Vertex, HashSet<Edge>> map) {
		HashMap<Vertex, HashSet<Edge>> clone = new HashMap<Vertex, HashSet<Edge>>();
		for (Entry<Vertex, HashSet<Edge>> ent: map.entrySet()) {
			HashSet<Edge> copyOfVals = new HashSet<Edge>();
			copyOfVals.addAll(ent.getValue());
			clone.put(ent.getKey(), copyOfVals);
		}
		return clone;
	}
	
	public static HashMap<Integer, HashSet<Vertex>> createPartMapClone(HashMap<Integer, HashSet<Vertex>> map) {
		HashMap<Integer, HashSet<Vertex>> clone = new HashMap<Integer, HashSet<Vertex>>();
		for (Entry<Integer, HashSet<Vertex>> ent: map.entrySet()) {
			HashSet<Vertex> copyOfVals = new HashSet<Vertex>();
			copyOfVals.addAll(ent.getValue());
			clone.put(ent.getKey(), copyOfVals);
		}
		return clone;
	}
	
	public static HashMap<String, Vertex> createStringMapClone(HashMap<String, Vertex> map) {
		HashMap<String, Vertex> clone = new HashMap<String, Vertex>();
		for (Entry<String, Vertex> ent: map.entrySet()) {
			clone.put(ent.getKey(), ent.getValue());
		}
		return clone;
	}
	
	public static String formatLine(int iter, Vertex[] vc) {
		StringBuilder line = new StringBuilder();

		line.append(iter);
		line.append("   [");
		int i = 0;
		if (vc[0] == null) {
			line.append("]                  ");
		} else {
			for (i = 0; i < vc.length - 1; i++) {
				line.append(vc[i].toString());
				line.append(", ");
			}
			line.append(vc[i].toString());
			line.append("]      ");
		}
		
		return line.toString();
	}
	
	public static String formatTreeTable (String inputFileName) {
		StringBuilder sb = new StringBuilder();
		
		ArrayList<String[]> records = new ArrayList<String[]>();
		try (CSVReader csvReader = new CSVReader(new FileReader(inputFileName));) {
		    String[] values = null;
		    while ((values = csvReader.readNext()) != null) {
		        records.add(values);
		    }
		} catch (Exception e) {
			System.out.print(e);
		}
		
		int depth = 0;
		for (int i = 0; i < records.size(); i++) {
			int size = records.get(i).length;
			if (size > depth) {
				depth = size;
			}
		}
		int numLeaves = records.size();
		String[][] array2D = new String[depth][numLeaves];// (String[][]) records.toArray();

		for (int i = 0; i < numLeaves; i++) {
			for (int j = 0; j < depth; j++) {
				String[] instance = records.get(i);
				if (j< instance.length) {
					array2D[j][i] = instance[j];
				} else {
					array2D[j][i] = "";
				}
			}
		}
		
		
		for (int i = 0; i < depth; i ++) {
			StringBuilder line1Builder = new StringBuilder();
			StringBuilder line2Builder = new StringBuilder();
			String[] edgesAtDepth = array2D[i];
			String currentEdge = edgesAtDepth[0];
			HashMap<String,Integer> parentEdges = new HashMap<String,Integer>();
			int edgeCount = 0;
			boolean startFlag = false;
			for (int j = 0; j < numLeaves; j++) {
				String nextEdge = edgesAtDepth[j];
				String parentEdge = array2D[i-1][j];
				if (nextEdge.equals(currentEdge) == false) {
					line1Builder.append(formatEdgeInSpace(edgeCount, currentEdge));
					line2Builder.append(formatDashAndSpace(edgeCount, startFlag, true));
					startFlag = true;
					currentEdge = nextEdge;
					edgeCount = 0;
				}
				edgeCount++;
			}
			line1Builder.append(formatEdgeInSpace(edgeCount, currentEdge));
			line2Builder.append(formatDashAndSpace(edgeCount, startFlag, false));

			sb.append(line2Builder.toString());
			sb.append('\n');
			sb.append(line1Builder.toString());
			sb.append('\n');
		}
		return sb.toString();
	}
	
	public static String formatEdgeInSpace(int space, String edge) {
		String edgeWidthBlank = "                  ";
		int edgeSpace = edgeWidthBlank.toCharArray().length;
		
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < space/2; i++) {
			sb.append(edgeWidthBlank);
		}
		sb.append(edge);
		for (int i = 0; i < space/2; i++) {
			sb.append(edgeWidthBlank);
		}
		return sb.toString();
	}
	
	public static String formatDashAndSpace(int space, boolean leftside, boolean rightside) {
		String edgeWidthBlank = "                  ";
		String dashWidthBlank = "-----------------";
		String delimiter = "|";
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < space/2; i++) {
			if (leftside == true) {
				sb.append(dashWidthBlank);
			} else {
				sb.append(edgeWidthBlank);
			}
		}
		if (leftside == true) {
			sb.append("---------");
		} else {
			sb.append("         ");
		}
		sb.append(delimiter);
		if (rightside == true) {
			sb.append("---------");
		} else {
			sb.append("         ");
		}
		for (int i = 0; i < space/2; i++) {
			if (rightside == true) {
				sb.append(dashWidthBlank);
			} else {
				sb.append(edgeWidthBlank);
			}
		}
		return sb.toString();
	}
	
	
//	public HashSet<Edge> generateAllGreedyMoves(int M, int N, Vertex[] vc) {
//		HashSet<Edge> greedyMoves = new HashSet<Edge>();
//		
//		if (vc[0] == null) {
//			Edge newEdge = new Edge(N);
//			newEdge.createFirstEdge();
//			for (Integer i = 0; i < N; i++) {
//				HashSet<Vertex> verts = this.partToVertexMap.get(i);
//				Vertex v = newEdge.values[i];
//				verts.add(v);
//				PartName name = PartName.valueOf(i);
//				int num = 1;
//				StringBuilder vertexName = new StringBuilder();
//				vertexName.append(name);
//				vertexName.append(num);
//			}
//			
//			greedyMoves.add(newEdge);
//			return greedyMoves;
//		}
//		
//		// Sets which represent what this edge has left to hit
//		HashMap<Vertex, HashSet<Edge>> remainingVertexToEdgeMap = createVertMapClone(this.vertexToEdgeMap);
//		HashSet<Edge> remainingEdgesToHit = new HashSet<Edge>();
//		remainingEdgesToHit.addAll(this.edges);
//		
//		// Initialize a for-now-empty edge.
//		Edge e = new Edge(N);
//		
//		HashMap<Integer, HashSet<Vertex>> initialValidPositions = createPartMapClone(this.partToVertexMap); 
//
//		// Remove the vertex cover vertices from possibilities to choose from.
//		for (int i = 0; i < vc.length; i++) {
//			Vertex vertex = vc[i];
//			//decrease the count of vertices in the partition
//			remainingVertexToEdgeMap.remove(vertex);
//			
//			// remove it from the set of possible vertices we can pick from
//			HashSet<Vertex> amendPossibleVerticesForPosition = initialValidPositions.get(vertex.part.value);
//			amendPossibleVerticesForPosition.remove(vertex);
//			initialValidPositions.put(vertex.part.value, amendPossibleVerticesForPosition);
//		}
//		HashMap<Integer, HashSet<Vertex>> remainingPositions = createPartMapClone(initialValidPositions); 
//		
//		recursiveGreedyHelper(remainingPositions, remainingEdgesToHit, remainingVertexToEdgeMap, e, greedyMoves);
//	}
//	
//	public void recursiveGreedyHelper(HashMap<Integer, HashSet<Vertex>> remainingPositions, 
//			HashSet<Edge> remainingEdgesToHit, HashMap<Vertex, HashSet<Edge>> remainingVertexToEdgeMap, Edge e, HashSet<Edge> greedyMoves) {
//		// Fill up the edge until there are no positions to fill anymore.
//		if (remainingPositions.isEmpty()) {
//			greedyMoves.add(e);
//		} else {
//			// Get the minimum sized part at this time
//			int min = 100000;
//			HashSet<Integer> partChoices = new HashSet<Integer>();
//			
//			for (Entry<Integer, HashSet<Vertex>> ent: remainingPositions.entrySet()) {
//				int size = ent.getValue().size();
//				if (size > 0 && size < min) {
//					min = size;
//					partChoices = new HashSet<Integer>();
//					partChoices.add(ent.getKey());
//				} else if (size > 0 && size == min) {
//					partChoices.add(ent.getKey());
//				}
//			}
//			
//			if (partChoices.isEmpty()) {
//				partChoices.addAll(remainingPositions.keySet());
//			}
//		
//			// Find the maximum vertices for that part 
//			// Choose the vertices from this set that has the maximum size
//			int max = 0;
//			Vertex bestVertex = null;
//			HashSet<Vertex> bestVertices = new HashSet<Vertex>();
//			HashSet<Vertex> vertexChoices = new HashSet<Vertex>();
//			for (Integer partChoice:partChoices) {
//				HashSet<Vertex> choices = remainingPositions.get(partChoice);
//				vertexChoices.addAll(choices);
//			}
//			if (vertexChoices.isEmpty()) {
//				for (Part p: this.partToVertexMap) {
//					// new vertices are the 
//				}
//				PartName name = PartName.valueOf(partChoice);
//				int num = this.partToVertexMap.get(name.value).size() + 1;
//				bestVertex = new Vertex(name, num);
//				
//			} else {
//				for (Vertex vertexChoice: vertexChoices) {
//					HashSet<Edge> edgeCover = remainingVertexToEdgeMap.get(v);
//					int size = edgeCover.size();
//					if (size > max) {
//						min = size;
//						bestVertices = new HashSet<Vertex>();
//						bestVertices.add(vertexChoice);
//					} else if (size == max) {
//						bestVertices.add(vertexChoice);
//					}
//				}
//			}
//			
//				
//				} else {
//					int minpart = -1;
//					int minval = -1;
//					for (Vertex v : choices) {
//						HashSet<Edge> edgeCover = remainingVertexToEdgeMap.get(v);
//						int size = edgeCover.size();
//						int part = v.part.value;
//						int val = v.value;
//						if (size > max) {
//							bestVertex = v;
//							minpart = part;
//							minval = val;
//							max = size;
//						} else if (size == max) {
//							if (minpart == -1 || part < minpart) {
//								minpart = part;
//								bestVertex = v;
//							} else if (part == minpart) {
//								if (minval == -1 || val < minval) {
//									minval = val;
//									bestVertex = v;
//									
//								}
//							}
//						}
//					}
//				}
//				if (choices == null) {
//					PartName name = PartName.valueOf(partChoice);
//					int num = this.partToVertexMap.get(name.value).size() + 1;
//					bestVertex = new Vertex(name, num);
//					
//					if(add_vc == false) {
//						//this.vertexCounts[partChoice] = num;
//						StringBuilder vertexName = new StringBuilder();
//						vertexName.append(name);
//						vertexName.append(num);
//						this.stringToVertexMap.put(vertexName.toString(), bestVertex);
//						this.vertexToEdgeMap.put(bestVertex, new HashSet<Edge>());
//					
//						if (this.partToVertexMap.containsKey(partChoice)) {
//							HashSet<Vertex> vertexChoices = this.partToVertexMap.get(partChoice);
//							vertexChoices.add(bestVertex);
//							this.partToVertexMap.put(partChoice, vertexChoices);
//						} else {
//							HashSet<Vertex> vertexChoices = new HashSet<Vertex>();
//							vertexChoices.add(bestVertex);
//							this.partToVertexMap.put(partChoice, vertexChoices);
//						}
//					}
//				}
//				else if (choices.isEmpty()) {
//					if (initialValidPositions.get(partChoice).isEmpty()) {
//						//make a new vertex
//						PartName name = PartName.valueOf(partChoice);
//						int num = this.partToVertexMap.get(partChoice).size() + 1; //vertexCounts[partChoice] + 1;
//						bestVertex = new Vertex(name, num);
//						
//						if (add_vc == false) {
//							//this.vertexCounts[partChoice] = num;
//							StringBuilder vertexName = new StringBuilder();
//							vertexName.append(name);
//							vertexName.append(num);
//							this.stringToVertexMap.put(vertexName.toString(), bestVertex);
//							this.vertexToEdgeMap.put(bestVertex, new HashSet<Edge>());
//							
//							if (this.partToVertexMap.containsKey(partChoice)) {
//								HashSet<Vertex> vertexChoices = this.partToVertexMap.get(partChoice);
//								vertexChoices.add(bestVertex);
//								this.partToVertexMap.put(partChoice, vertexChoices);
//							} else {
//								HashSet<Vertex> vertexChoices = new HashSet<Vertex>();
//								vertexChoices.add(bestVertex);
//								this.partToVertexMap.put(partChoice, vertexChoices);
//							}
//						}
//					} else {
//						// OLD CODE: arbitrarily choose a vertex from the initial set
//						// bestVertex = initialValidPositions.get(partChoice).iterator().next();
//						
//						// NEW code: choose 'greedily' (ie the lowest entry available for this part)
//						HashSet<Vertex> verticesInPart = initialValidPositions.get(partChoice);
//						int lowestVal = 0;
//						for (Vertex v: verticesInPart) {
//							if (v.value < lowestVal) {
//								lowestVal = v.value;
//								bestVertex = v;
//							}
//						}		
//					} 
//			}
//			
//			for(Vertex bestVertex:bestVertices) {
//				Edge newEdge = new Edge(e);
//				newEdge.setVertex(bestVertex);
//			}
//
//			// Fix the vertex in the new edge
//			e.setVertex(bestVertex);
//			
//			// Remove this part from the parts remaining and reset its count so it doesn't get chosen again
//			remainingPositions.remove(partChoice);	
//			//remainingVertexCounts[partChoice] = -1;
//			
//			// Remove all edges to hit which contain this vertex
//			HashSet<Edge> edgesWhichAreNowHit = new HashSet<Edge>();
//			HashSet<Edge> ptr = remainingVertexToEdgeMap.get(bestVertex);
//			
//			if (ptr != null && !ptr.isEmpty()) {
//				edgesWhichAreNowHit.addAll(ptr);
//				for (Edge hit: edgesWhichAreNowHit) {
//					remainingEdgesToHit.remove(hit);
//					
//					Vertex[] containedVertices = hit.values;
//					for (Vertex v: containedVertices) {
//						HashSet<Edge> tmp = remainingVertexToEdgeMap.get(v);
//						if (tmp != null) {
//							tmp.remove(hit);
//							remainingVertexToEdgeMap.put(v, tmp);
//						}
//					}
//				}
//			}
//		}
//		
//		if (remainingEdgesToHit.isEmpty()) {
//			return e;
//		} else {
//			System.out.println("No possible move.");
//
//			return null;
//		}
//		
//	}

}
