#include <iostream>
#include <vector>
#include <algorithm>

const int MAX_V = 100;
const int INF = 1000000000;

class Edge {
private:
    int source;
    int destination;
    int weight;

public:
    Edge(int src, int dest, int w) : source(src), destination(dest), weight(w) {}

    int getSource() const { return source; }
    int getDestination() const { return destination; }
    int getWeight() const { return weight; }
};

class Graph {
private:
    int adjacencyMatrix[MAX_V][MAX_V];
    int numVertices;

public:
    Graph(int vertices) : numVertices(vertices) {
        for (int i = 0; i < MAX_V; i++) {
            for (int j = 0; j < MAX_V; j++) {
                adjacencyMatrix[i][j] = 0;
            }
        }
    }

    void addEdge(int u, int v, int weight) {
        adjacencyMatrix[u][v] = weight;
        adjacencyMatrix[v][u] = weight;  // For undirected graph
    }

    std::vector<int> getNeighbors(int vertex) {
        std::vector<int> neighbors;
        for (int i = 0; i < numVertices; i++) {
            if (adjacencyMatrix[vertex][i] != 0) {
                neighbors.push_back(i);
            }
        }
        return neighbors;
    }

    int getWeight(int u, int v) {
        return adjacencyMatrix[u][v];
    }

    std::vector<int> getAllVertices() {
        std::vector<int> vertices;
        for (int i = 0; i < numVertices; i++) {
            vertices.push_back(i);
        }
        return vertices;
    }

    bool isConnected(int u, int v) {
        return adjacencyMatrix[u][v] != 0;
    }

    int getNumVertices() const { return numVertices; }
};

class PriorityQueue {
private:
    int distances[MAX_V];
    int nodes[MAX_V];
    int size;

    void heapify(int i) {
        int left = 2 * i + 1;
        int right = 2 * i + 2;
        int min = i;

        if (left < size && distances[left] < distances[min]) {
            min = left;
        }
        if (right < size && distances[right] < distances[min]) {
            min = right;
        }
        if (min != i) {
            std::swap(distances[i], distances[min]);
            std::swap(nodes[i], nodes[min]);
            heapify(min);
        }
    }

public:
    PriorityQueue() : size(0) {}

    bool isEmpty() const { return size == 0; }

    void add(int distance, int node) {
        distances[size] = distance;
        nodes[size] = node;
        size++;
        for (int i = (size - 1) / 2; i >= 0; i--) {
            heapify(i);
        }
    }

    void pop() {
        if (!isEmpty()) {
            distances[0] = distances[size - 1];
            nodes[0] = nodes[size - 1];
            size--;
            heapify(0);
        }
    }

    int topNode() const { return nodes[0]; }
    int topDistance() const { return distances[0]; }
};

class AlgorithmResult {
private:
    enum class ResultType { SHORTEST_PATH, MST };
    ResultType resultType;
    std::vector<int> path;
    int distance;
    std::vector<Edge> mstEdges;
    int totalWeight;

public:
    void setPathResult(const std::vector<int>& p, int d) {
        resultType = ResultType::SHORTEST_PATH;
        path = p;
        distance = d;
    }

    void setMSTResult(const std::vector<Edge>& edges, int weight) {
        resultType = ResultType::MST;
        mstEdges = edges;
        totalWeight = weight;
    }
};

class DijkstraAlgorithm {
private:
    Graph& graph;
    PriorityQueue priorityQueue;
    std::vector<int> path;
    int shortestDistance;

    void calculateDistances(int start) {
        std::vector<int> distances(graph.getNumVertices(), INF);
        std::vector<bool> visited(graph.getNumVertices(), false);
        std::vector<int> previous(graph.getNumVertices(), -1);

        distances[start] = 0;
        priorityQueue.add(0, start);

        while (!priorityQueue.isEmpty()) {
            int current = priorityQueue.topNode();
            int currentDist = priorityQueue.topDistance();
            priorityQueue.pop();

            if (visited[current]) continue;
            visited[current] = true;

            for (int neighbor : graph.getNeighbors(current)) {
                if (!visited[neighbor]) {
                    int newDist = currentDist + graph.getWeight(current, neighbor);
                    if (newDist < distances[neighbor]) {
                        distances[neighbor] = newDist;
                        previous[neighbor] = current;
                        priorityQueue.add(newDist, neighbor);
                    }
                }
            }
        }

        shortestDistance = distances[graph.getNumVertices() - 1];
        reconstructPath(previous, start, graph.getNumVertices() - 1);
    }

    void reconstructPath(const std::vector<int>& previous, int start, int end) {
        path.clear();
        for (int at = end; at != -1; at = previous[at]) {
            path.push_back(at);
        }
        std::reverse(path.begin(), path.end());
    }

public:
    DijkstraAlgorithm(Graph& g) : graph(g) {}

    void findShortestPath(int start, int end) {
        calculateDistances(start);
    }

    std::vector<int> getPath() const { return path; }
    int getDistance() const { return shortestDistance; }
};

class PrimAlgorithm {
private:
    Graph& graph;
    PriorityQueue priorityQueue;
    std::vector<Edge> mstEdges;
    int totalWeight;

public:
    PrimAlgorithm(Graph& g) : graph(g), totalWeight(0) {}

    void findMST() {
        std::vector<bool> inMST(graph.getNumVertices(), false);
        std::vector<int> key(graph.getNumVertices(), INF);
        std::vector<int> parent(graph.getNumVertices(), -1);

        priorityQueue.add(0, 0);
        key[0] = 0;

        while (!priorityQueue.isEmpty()) {
            int u = priorityQueue.topNode();
            priorityQueue.pop();

            if (inMST[u]) continue;
            inMST[u] = true;

            for (int v : graph.getNeighbors(u)) {
                int weight = graph.getWeight(u, v);
                if (!inMST[v] && weight < key[v]) {
                    key[v] = weight;
                    parent[v] = u;
                    priorityQueue.add(key[v], v);
                }
            }
        }

        // Construct MST edges
        mstEdges.clear();
        totalWeight = 0;
        for (int i = 1; i < graph.getNumVertices(); i++) {
            if (parent[i] != -1) {
                int weight = graph.getWeight(parent[i], i);
                mstEdges.emplace_back(parent[i], i, weight);
                totalWeight += weight;
            }
        }
    }

    std::vector<Edge> getMSTEdges() const { return mstEdges; }
    int getTotalWeight() const { return totalWeight; }

};


int main() {
    // Create a graph with 9 vertices (0 to 8)
    Graph graph(9);
    
    // Add edges with weights
    graph.addEdge(0, 1, 4);
    graph.addEdge(0, 7, 8);
    graph.addEdge(1, 2, 8);
    graph.addEdge(1, 7, 11);
    graph.addEdge(2, 3, 7);
    graph.addEdge(2, 8, 2);
    graph.addEdge(2, 5, 4);
    graph.addEdge(3, 4, 9);
    graph.addEdge(3, 5, 14);
    graph.addEdge(4, 5, 10);
    graph.addEdge(5, 6, 2);
    graph.addEdge(6, 7, 1);
    graph.addEdge(6, 8, 6);
    graph.addEdge(7, 8, 7);

    // Test Dijkstra's algorithm
    DijkstraAlgorithm dijkstra(graph);
    dijkstra.findShortestPath(0, 8);
    
    // Print shortest path and distance
    std::cout << "Shortest path from 0 to 8:" << std::endl;
    for (int vertex : dijkstra.getPath()) {
        std::cout << vertex << " ";
    }
    std::cout << "\nDistance: " << dijkstra.getDistance() << std::endl;
    
    // Test Prim's algorithm
    PrimAlgorithm prim(graph);
    prim.findMST();
    
    // Print MST edges and total weight
    std::cout << "\nMinimum Spanning Tree edges:" << std::endl;
    for (const Edge& edge : prim.getMSTEdges()) {
        std::cout << edge.getSource() << " - " << edge.getDestination() 
                 << " (weight: " << edge.getWeight() << ")" << std::endl;
    }
    std::cout << "Total MST weight: " << prim.getTotalWeight() << std::endl;

    return 0;
}
