#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <climits>
using namespace std;

//class node indexes for diktra algirith
struct Node {
    int index;
    int weight;

    Node(int i, int w) : index(i), weight(w) {}
};
//airport class
class Airport {
public:
    string code;    // like "ATL" for the destination
    string city;
    string state;

    Airport() {}
    Airport(string c, string ci, string s)
        : code(c), city(ci), state(s) {}
};
//edge class
class Edge {
public:
    int destIndex;
    int cost;
    int distance;

    Edge(int dest, int d, int c) : destIndex(dest), distance(d), cost(c) {}
};
//minhead class with all functions
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

class Graph {
public:
    vector<Airport> airports;
    vector<vector<Edge>> adjList;

    void addAirport(const Airport& a) {
        airports.push_back(a);
        adjList.push_back(vector<Edge>());
    }

    void addFlight(int originIndex, int destIndex, int distance, int cost) {
        adjList[originIndex].push_back(Edge(destIndex, distance, cost));
    }

    int getAirportIndex(const string& code) {
        for (int i = 0; i < airports.size(); ++i) {
            if (airports[i].code == code)
                return i;
        }
        return -1;
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

        vector<int> path;
        for (int at = end; at != -1; at = prev[at])
            path.push_back(at);
        reverse(path.begin(), path.end());

        cout << "Shortest route from " << startCode << " to " << endCode << ": ";
        for (int i = 0; i < path.size(); ++i) {
            cout << airports[path[i]].code;
            if (i < path.size() - 1) cout << " -> ";
        }
        cout << ". The " << (useCost ? "cost" : "length") << " is " << dist[end] << ".\n";
    }
    vector<string> split(const string& line) {
        vector<string> tokens;
        string token;
        bool inQuotes = false;

        for (char c : line) {
            if (c == '"') {
                inQuotes = !inQuotes;  // toggle quotes
            } else if (c == ',' && !inQuotes) {
                tokens.push_back(token);
                token.clear();
            } else {
                token += c;
            }
        }
        tokens.push_back(token); // add last token
        return tokens;
    }


    string extractState(const string& cityState) {
        size_t commaPos = cityState.find(',');
        if (commaPos != string::npos && commaPos + 2 < cityState.size())
            return cityState.substr(commaPos + 2);
        return "";
    }

    void loadFromCSV(const string& filename) {
        ifstream file(filename);
        if (!file.is_open()) {
            cout << "Error: Cannot open file.\n";
            return;
        }

        string line;
        map<string, int> airportMap;
        int lineNum = 0;

        while (getline(file, line)) {
            lineNum++;
            vector<string> parts = splitCSVLine(line);  // Use the new function

            if (parts.size() < 6) {
                cout << "Skipping line " << lineNum << ": not enough parts.\n";
                continue;
            }

            string originCode = parts[0];
            string destCode = parts[1];
            string originCity = parts[2];
            string destCity = parts[3];

            int distance, cost;
            try {
                distance = stoi(parts[4]);
                cost = stoi(parts[5]);
            } catch (const std::exception& e) {
                cout << "Skipping line " << lineNum << ": invalid number. " << e.what() << "\n";
                continue;
            }

            int originIdx, destIdx;

            if (airportMap.find(originCode) == airportMap.end()) {
                originIdx = airports.size();
                airportMap[originCode] = originIdx;
                addAirport(Airport(originCode, originCity, extractState(originCity)));
            } else {
                originIdx = airportMap[originCode];
            }

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
//helper function for leading file
    vector<string> splitCSVLine(const string& line) {
        vector<string> tokens;
        string token;
        bool inQuotes = false;

        for (char c : line) {
            if (c == '"') {
                inQuotes = !inQuotes;
            } else if (c == ',' && !inQuotes) {
                tokens.push_back(token);
                token.clear();
            } else {
                token += c;
            }
        }
        tokens.push_back(token);
        return tokens;
    }



    void countDirectConnections() {
        int n = airports.size();
        vector<int> inbound(n, 0);
        vector<int> outbound(n, 0);

        for (int u = 0; u < n; ++u) {
            for (Edge& e : adjList[u]) {
                outbound[u]++;
                inbound[e.destIndex]++;
            }
        }

        vector<pair<int, int>> connectionData;
        for (int i = 0; i < n; ++i) {
            int total = inbound[i] + outbound[i];
            connectionData.push_back({i, total});
        }

        for (int i = 0; i < connectionData.size(); ++i) {
            int maxIdx = i;
            for (int j = i + 1; j < connectionData.size(); ++j) {
                if (connectionData[j].second > connectionData[maxIdx].second)
                    maxIdx = j;
            }
            swap(connectionData[i], connectionData[maxIdx]);
        }

        cout << "\nAirport\tConnections\n";
        for (auto& pair : connectionData) {
            cout << airports[pair.first].code << "\t" << pair.second << "\n";
        }
    }

    Graph createUndirectedGraph();
    void primMST();
    void kruskalMST();
};

Graph Graph::createUndirectedGraph() {
    Graph G_u;
    int n = airports.size();
    for (int i = 0; i < n; ++i)
        G_u.addAirport(airports[i]);

    vector<vector<bool>> processed(n, vector<bool>(n, false));

    for (int u = 0; u < n; ++u) {
        for (Edge& e : adjList[u]) {
            int v = e.destIndex;
            if (processed[u][v] || processed[v][u]) continue;

            int costUV = -1, costVU = -1;

            for (Edge& edge : adjList[u]) if (edge.destIndex == v) costUV = edge.cost;
            for (Edge& edge : adjList[v]) if (edge.destIndex == u) costVU = edge.cost;

            int minCost = min((costUV != -1 ? costUV : INT_MAX), (costVU != -1 ? costVU : INT_MAX));
            G_u.addFlight(u, v, 0, minCost);
            G_u.addFlight(v, u, 0, minCost);

            processed[u][v] = processed[v][u] = true;
        }
    }

    cout << "Undirected graph G_u created with " << G_u.airports.size() << " airports.\n";
    return G_u;
}

//prism algorithm
void Graph::primMST() {
    int n = airports.size();
    vector<int> key(n, INT_MAX), parent(n, -1);
    vector<bool> inMST(n, false);

    key[0] = 0;
    for (int count = 0; count < n - 1; ++count) {
        int minKey = INT_MAX, u = -1;
        for (int v = 0; v < n; ++v)
            if (!inMST[v] && key[v] < minKey) {
                minKey = key[v];
                u = v;
            }

        if (u == -1) {
            cout << "Graph is disconnected. MST cannot be formed.\n";
            return;
        }

        inMST[u] = true;
        for (Edge& e : adjList[u]) {
            int v = e.destIndex;
            if (!inMST[v] && e.cost < key[v]) {
                key[v] = e.cost;
                parent[v] = u;
            }
        }
    }

    int totalCost = 0;
    cout << "\nMinimal Spanning Tree (Prim's):\nEdge\tWeight\n";
    for (int i = 1; i < n; ++i) {
        if (parent[i] != -1) {
            cout << airports[parent[i]].code << " - " << airports[i].code << "\t" << key[i] << "\n";
            totalCost += key[i];
        }
    }
    cout << "Total Cost of MST: " << totalCost << "\n";
}

//kruskals algorithm
void Graph::kruskalMST() {
    struct FullEdge {
        int u, v, cost;
        FullEdge(int _u, int _v, int _c) : u(_u), v(_v), cost(_c) {}
    };

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

    for (int i = 0; i < edges.size(); ++i) {
        int minIdx = i;
        for (int j = i + 1; j < edges.size(); ++j)
            if (edges[j].cost < edges[minIdx].cost) minIdx = j;
        swap(edges[i], edges[minIdx]);
    }

    int n = airports.size();
    vector<int> parent(n), rank(n, 0);
    for (int i = 0; i < n; ++i) parent[i] = i;

    function<int(int)> find = [&](int i) {
        while (i != parent[i]) i = parent[i];
        return i;
    };

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

    vector<FullEdge> mst;
    int totalCost = 0;

    for (FullEdge& e : edges)
        if (unionSet(e.u, e.v)) {
            mst.push_back(e);
            totalCost += e.cost;
        }

    cout << "\nMinimal Spanning Tree (Kruskal's):\nEdge\t\tWeight\n";
    for (FullEdge& e : mst)
        cout << airports[e.u].code << " - " << airports[e.v].code << "\t" << e.cost << "\n";
    cout << "Total Cost of MST: " << totalCost << "\n";
}

//main function calling loading file and displaying information
int main() {
    Graph g;
    g.loadFromCSV("airports.csv");
    g.countDirectConnections();
    Graph g_u = g.createUndirectedGraph();
    g_u.primMST();
    g_u.kruskalMST();
    return 0;
}
