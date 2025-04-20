class Edge {
public:
    std::string destination;
    int distance;
    int cost;
    Edge* next; // For adjacency list

    Edge(std::string dest, int dist, int c);
};