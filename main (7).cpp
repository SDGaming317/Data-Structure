
#include <iostream>
using namespace std;

struct Node {
	int index;
	int weight;

	Node(int i, int w) : index(i), weight(w) {}
};

class Airport {
public:
	string code;    // like "ATL" for the destination
	string city;
	string state;

	Airport() {}
	Airport(string c, string ci, string s)
		: code(c), city(ci), state(s) {}
};

class Edge {
public:
	int destIndex;     // Index of destination airport in airport list
	int cost;
	int distance;

	Edge(int dest, int d, int c) : destIndex(dest), distance(d), cost(c) {}
};

class Graph {
public:
	vector<Airport> airports;
	vector<vector<Edge>> adjList; // Adjacency list

	// Add an airport
	void addAirport(const Airport& a) {
		airports.push_back(a);
		adjList.push_back(std::vector<Edge>());
	}

	// Add a flight (directed edge)
	void addFlight(int originIndex, int destIndex, int distance, int cost) {
		adjList[originIndex].push_back(Edge(destIndex, distance, cost));
	}

	// Helper: Find airport index by code
	int getAirportIndex(const std::string& code) {
		for (int i = 0; i < airports.size(); ++i) {
			if (airports[i].code == code)
				return i;
		}
		return -1; // Not found
	}

	void dijkstraShortestPath(const string& startCode, const string& endCode, bool useCost = false) {
		int n = airports.size();
		int start = getAirportIndex(startCode);
		int end = getAirportIndex(endCode);

		if (start == -1 || end == -1) {
			cout << "Invalid airport code.\n";
			return;
		}

		vector<int> dist(n, INT_MAX);
		vector<int> prev(n, -1);
		vector<bool> visited(n, false);

		MinHeap pq;
		pq.push(Node(start, 0));
		dist[start] = 0;

		while (!pq.isEmpty()) {
			Node current = pq.pop();
			int u = current.index;

			if (visited[u]) continue;
			visited[u] = true;

			for (Edge& e : adjList[u]) {
				int weight = useCost ? e.cost : e.distance;
				int v = e.destIndex;

				if (dist[u] + weight < dist[v]) {
					dist[v] = dist[u] + weight;
					prev[v] = u;
					pq.push(Node(v, dist[v]));
				}
			}
		}

		if (dist[end] == INT_MAX) {
			cout << "No path from " << startCode << " to " << endCode << ".\n";
			return;
		}

		// Reconstruct path
		vector<int> path;
		for (int at = end; at != -1; at = prev[at])
			path.push_back(at);
		reverse(path.begin(), path.end());

		cout << "Shortest route from " << startCode << " to " << endCode << ": ";
		for (int i = 0; i < path.size(); ++i) {
			cout << airports[path[i]].code;
			if (i < path.size() - 1) std::cout << " -> ";
		}
		cout << ". The " << (useCost ? "cost" : "length") << " is " << dist[end] << ".\n";
	}

	vector<string> split(const string& s, char delimiter) {
		vector<string> tokens;
		string token;
		stringstream ss(s);

		while (getline(ss, token, delimiter)) {
			tokens.push_back(token);
		}
		return tokens;
	}

	string extractState(const string& cityState) {
		size_t commaPos = cityState.find(',');
		if (commaPos != string::npos && commaPos + 2 < cityState.size())
			return cityState.substr(commaPos + 2); // skip ", "
		return "";
	}
	void loadFromCSV(const std::string& filename) {
		ifstream file(filename);
		if (!file.is_open()) {
			cout << "Error: Cannot open file.\n";
			return;
		}

		string line;
		map<string, int> airportMap; // airport code -> index

		while (getline(file, line)) {
			vector<string> parts = split(line, ',');

			if (parts.size() < 6) continue;

			string originCode = parts[0];
			string destCode = parts[1];
			string originCity = parts[2];
			string destCity = parts[3];
			int distance = stoi(parts[4]);
			int cost = stoi(parts[5]);

			int originIdx, destIdx;

			// Add origin airport if not already present
			if (airportMap.find(originCode) == airportMap.end()) {
				originIdx = airports.size();
				airportMap[originCode] = originIdx;
				addAirport(Airport(originCode, originCity, extractState(originCity)));
			} else {
				originIdx = airportMap[originCode];
			}

			// Add destination airport
			if (airportMap.find(destCode) == airportMap.end()) {
				destIdx = airports.size();
				airportMap[destCode] = destIdx;
				addAirport(Airport(destCode, destCity, extractState(destCity)));
			} else {
				destIdx = airportMap[destCode];
			}

			addFlight(originIdx, destIdx, distance, cost);
		}

		file.close();
		cout << "Graph loaded with " << airports.size() << " airports.\n";
	}

	void countDirectConnections() {
		int n = airports.size();
		std::vector<int> inbound(n, 0);
		std::vector<int> outbound(n, 0);

		// Count outbound and inbound
		for (int u = 0; u < n; ++u) {
			for (Edge& e : adjList[u]) {
				outbound[u]++;
				inbound[e.destIndex]++;
			}
		}

		// Build a vector of (index, total connections)
		std::vector<std::pair<int, int>> connectionData;
		for (int i = 0; i < n; ++i) {
			int total = inbound[i] + outbound[i];
			connectionData.push_back(std::make_pair(i, total));
		}

		// Sort manually (selection sort)
		for (int i = 0; i < connectionData.size(); ++i) {
			int maxIdx = i;
			for (int j = i + 1; j < connectionData.size(); ++j) {
				if (connectionData[j].second > connectionData[maxIdx].second)
					maxIdx = j;
			}
			std::swap(connectionData[i], connectionData[maxIdx]);
		}

		// Output
		std::cout << "\nAirport\tConnections\n";
		for (auto& pair : connectionData) {
			int idx = pair.first;
			std::cout << airports[idx].code << "\t" << pair.second << "\n";
		}
	}
	Graph createUndirectedGraph();
	// Functions to implement:
	// - DFS/BFS (Task 4)
	// - Count connections (Task 5)
	// - Create G_u (Task 6)
	// - Prim's (Task 7), Kruskal's (Task 8)
};


Graph createUndirectedGraph() {
	Graph G_u;

	int n = airports.size();
	for (int i = 0; i < n; ++i)
		G_u.addAirport(airports[i]);

	// Matrix to track if edge has been handled
	vector<vector<bool>> processed(n, vector<bool>(n, false));

	for (int u = 0; u < n; ++u) {
		for (Edge& e : adjList[u]) {
			int v = e.destIndex;
			if (processed[u][v] || processed[v][u]) continue;

			int costUV = -1, costVU = -1;

			// Find (u b v)
			for (Edge& edge : adjList[u]) {
				if (edge.destIndex == v) {
					costUV = edge.cost;
					break;
				}
			}

			// Find (v b u)
			for (Edge& edge : adjList[v]) {
				if (edge.destIndex == u) {
					costVU = edge.cost;
					break;
				}
			}

			if (costUV != -1 && costVU != -1) {
				// Both directions exist b keep cheaper one
				int minCost = std::min(costUV, costVU);
				G_u.addFlight(u, v, 0, minCost);
				G_u.addFlight(v, u, 0, minCost);
			} else if (costUV != -1) {
				// Only u b
				G_u.addFlight(u, v, 0, costUV);
				G_u.addFlight(v, u, 0, costUV);
			} else if (costVU != -1) {
				// Only v b u
				G_u.addFlight(u, v, 0, costVU);
				G_u.addFlight(v, u, 0, costVU);
			}

			processed[u][v] = processed[v][u] = true;
		}
	}

	cout << "Undirected graph G_u created with " << G_u.airports.size() << " airports.\n";
	return G_u;
}


class MinHeap {
	vector<Node> heap;

	void heapifyUp(int i) {
		while (i > 0 && heap[i].weight < heap[(i - 1) / 2].weight) {
			swap(heap[i], heap[(i - 1) / 2]);
			i = (i - 1) / 2;
		}
	}

	void heapifyDown(int i) {
		int left = 2 * i + 1, right = 2 * i + 2, smallest = i;
		if (left < heap.size() && heap[left].weight < heap[smallest].weight)
			smallest = left;
		if (right < heap.size() && heap[right].weight < heap[smallest].weight)
			smallest = right;
		if (smallest != i) {
			swap(heap[i], heap[smallest]);
			heapifyDown(smallest);
		}
	}

public:
	void push(Node n) {
		heap.push_back(n);
		heapifyUp(heap.size() - 1);
	}

	Node pop() {
		Node minNode = heap[0];
		heap[0] = heap.back();
		heap.pop_back();
		heapifyDown(0);
		return minNode;
	}

	bool isEmpty() const {
		return heap.empty();
	}
};

void Graph::primMST() {
	int n = airports.size();
	vector<int> key(n, INT_MAX);  // edge weights
	vector<int> parent(n, -1);    // store MST structure
	vector<bool> inMST(n, false); // track visited nodes

	key[0] = 0;  // start from node 0

	for (int count = 0; count < n - 1; ++count) {
		// Find the minimum key vertex not yet included
		int minKey = INT_MAX, u = -1;
		for (int v = 0; v < n; ++v) {
			if (!inMST[v] && key[v] < minKey) {
				minKey = key[v];
				u = v;
			}
		}

		if (u == -1) {
			std::cout << "Graph is disconnected. MST cannot be formed.\n";
			return;
		}

		inMST[u] = true;

		// Update keys for neighbors
		for (Edge& e : adjList[u]) {
			int v = e.destIndex;
			int weight = e.cost;
			if (!inMST[v] && weight < key[v]) {
				key[v] = weight;
				parent[v] = u;
			}
		}
	}

	// Print the MST
	int totalCost = 0;
	cout << "\nMinimal Spanning Tree (Primbs):\n";
	cout << "Edge\tWeight\n";

	for (int i = 1; i < n; ++i) {
		if (parent[i] == -1) continue; // skip unconnected nodes
		std::cout << airports[parent[i]].code << " - " << airports[i].code
		          << "\t" << key[i] << "\n";
		totalCost += key[i];
	}

	cout << "Total Cost of MST: " << totalCost << "\n";
}


void Graph::kruskalMST() {
	struct FullEdge {
		int u, v, cost;
		FullEdge(int _u, int _v, int _c) : u(_u), v(_v), cost(_c) {}
	};

	// Step 1: Extract all undirected edges
	vector<FullEdge> edges;
	vector<vector<bool>> added(airports.size(), vector<bool>(airports.size(), false));

	for (int u = 0; u < airports.size(); ++u) {
		for (Edge& e : adjList[u]) {
			int v = e.destIndex;
			if (!added[u][v] && !added[v][u]) {
				edges.push_back(FullEdge(u, v, e.cost));
				added[u][v] = added[v][u] = true;
			}
		}
	}

	// Step 2: Sort edges by cost using selection sort (no std::sort)
	for (int i = 0; i < edges.size(); ++i) {
		int minIdx = i;
		for (int j = i + 1; j < edges.size(); ++j) {
			if (edges[j].cost < edges[minIdx].cost)
				minIdx = j;
		}
		swap(edges[i], edges[minIdx]);
	}

	// Step 3: Union-Find setup
	int n = airports.size();
	vector<int> parent(n), rank(n, 0);
	for (int i = 0; i < n; ++i) parent[i] = i;

	// Helper: find root
	auto find = [&](int i) -> int {
		while (i != parent[i]) i = parent[i];
		return i;
	};

	// Helper: union by rank
	auto unionSet = [&](int x, int y) {
		int rootX = find(x), rootY = find(y);
		if (rootX == rootY) return false;
		if (rank[rootX] < rank[rootY]) parent[rootX] = rootY;
		else if (rank[rootX] > rank[rootY]) parent[rootY] = rootX;
		else {
			parent[rootY] = rootX;
			rank[rootX]++;
		}
		return true;
	};

	// Step 4: Build MST
	vector<FullEdge> mst;
	int totalCost = 0;

	for (FullEdge& e : edges) {
		if (unionSet(e.u, e.v)) {
			mst.push_back(e);
			totalCost += e.cost;
		}
	}

	// Output
	cout << "\nMinimal Spanning Tree (Kruskalbs):\n";
	cout << "Edge\t\tWeight\n";
	for (FullEdge& e : mst) {
		cout << airports[e.u].code << " - " << airports[e.v].code
		     << "\t" << e.cost << "\n";
	}
	cout << "Total Cost of MST: " << totalCost << "\n";

	if (mst.size() < n



	        int main()
{


	return 0;
}