#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <chrono>
#include <set>
#include <functional>

using namespace std;
// Struct for storing point data
struct Point {
    int id;
    double x, y;
    
    Point(int id, double x, double y) : id(id), x(x), y(y) {}
};

// Struct for a dataset
struct Dataset {
    string name;
    vector<Point> points;
    vector<vector<double>> distanceMatrix;
    Dataset(string name, vector<Point> &&points, vector<vector<double>> &&distanceMatrix) : name(name), points(points), distanceMatrix(distanceMatrix) {}
};


// Function to calculate Euclidean distance between two points
double calculateDistance(const Point& a, const Point& b) 
{
    return round(sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y)));
}

Dataset parseTSPFile(const std::string& filename) {
    std::ifstream file(filename);
    std::string line;
    std::string name;
    std::vector<Point> points;

    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return Dataset("", {}, {});
    }

    while (std::getline(file, line)) {
        // Parse the name
        if (line.rfind("NAME", 0) == 0) {
            std::istringstream iss(line);
            std::string discard;
            iss >> discard >> name;
        }
        
        // Start reading points after the NODE_COORD_SECTION line
        else if (line == "NODE_COORD_SECTION") {
            while (std::getline(file, line) && line != "EOF") {
                std::istringstream iss(line);
                int id;
                double x, y;
                if (iss >> id >> x >> y) {
                    points.emplace_back(id, x, y);
                }
            }
        }
    }
    vector<std::vector<double>> distanceMatrix(points.size(), vector<double>(points.size(), 0.0));
    for (int i = 0; i < points.size(); ++i) {
        for (int j = i + 1; j < points.size(); ++j) {
            double dx = points[i].x - points[j].x;
            double dy = points[i].y - points[j].y;
            double distance = std::round(std::sqrt(dx * dx + dy * dy));
            distanceMatrix[i][j] = distance;
            distanceMatrix[j][i] = distance; // Symmetric matrix
        }
    }


    file.close();
    return Dataset(name, std::move(points), std::move(distanceMatrix));
}
double bruteForceRecursive(size_t start, size_t end, vector<size_t> &sequence, const vector<vector<double>> &distanceMatrix)
{
    if (start == end)
    {
        // Calculate accumulated distance
        // Store the answer
    }
    for (int i = start; i <= end; i++)
    {
        swap(sequence[i], start);
        bruteForceRecursive(start + 1, end, sequence, distanceMatrix);
        swap(sequence[i], start);
    }
    return 0.0;
}
// Placeholder for exact algorithm (brute-force)
vector<int> exactAlgorithm(const Dataset& dataset, int timeLimit) {
    // Implement brute-force algorithm here
    vector<size_t> sequence;
    for (int i = 0; i < dataset.points.size(); i++)
    {
        sequence.push_back(i);
    }
    bruteForceRecursive(0, dataset.points.size() - 1, sequence, dataset.distanceMatrix);
    return {0}; // Replace with computed tour
}

// Placeholder for approximation algorithm (MST-based 2-approximation)
// vector<int> approximateAlgorithm(const vector<Point>& points) {
vector<int> approximateAlgorithm(const Dataset& dataset) {
    // Implement 2-approximation algorithm here
    // Placeholder for result
    int n = dataset.points.size();
    const auto& distanceMatrix = dataset.distanceMatrix; // 使用已建立的距離矩陣

    // Step 1: Build MST using Prim's algorithm
    vector<double> minEdge(n, numeric_limits<double>::max());
    vector<int> parent(n, -1);
    vector<bool> inMST(n, false);
    minEdge[0] = 0.0;

    for (int i = 0; i < n; ++i) {
        int u = -1;

        // Find the vertex with the minimum edge that is not yet included in MST
        for (int v = 0; v < n; ++v) {
            if (!inMST[v] && (u == -1 || minEdge[v] < minEdge[u])) {
                u = v;
            }
        }

        inMST[u] = true;

        // Update the edges for the MST
        for (int v = 0; v < n; ++v) {
            if (distanceMatrix[u][v] && !inMST[v] && distanceMatrix[u][v] < minEdge[v]) {
                minEdge[v] = distanceMatrix[u][v];
                parent[v] = u;
            }
        }
    }

    // Step 2: Perform a pre-order traversal of the MST to create a TSP tour
    vector<int> tour;
    set<int> visited;

    function<void(int)> preorderTraversal = [&](int node) {
        if (visited.count(node)) return;
        visited.insert(node);
        tour.push_back(node);

        for (int i = 0; i < n; ++i) {
            if (parent[i] == node || parent[node] == i) {
                preorderTraversal(i);
            }
        }
    };

    preorderTraversal(0); // Start the traversal from the first node

    // Step 3: Calculate the cost of the approximate tour
    double tourCost = 0.0;
    for (int i = 0; i < tour.size() - 1; ++i) {
        tourCost += distanceMatrix[tour[i]][tour[i + 1]];
    }
    // Add the cost to return to the starting point
    tourCost += distanceMatrix[tour.back()][tour[0]];

    // Print the approximate tour cost
    cout << "Approximate Tour Cost: " << tourCost << endl;

    return tour;
    // return {0}; // Replace with computed tour
}

// Placeholder for local search algorithm (e.g., Simulated Annealing)
vector<int> localSearchAlgorithm(const vector<Point>& points, int timeLimit, int seed) {
    // Implement local search algorithm here
    // Placeholder for result
    return {0}; // Replace with computed tour
}

// Function to save solution to file
void saveSolution(const string& instance, const string& method, int timeLimit, int seed, double quality, const vector<int>& tour) {
    string filename = instance + " " + method + " " + to_string(timeLimit) + (method == "LS" ? " " + to_string(seed) : "") + ".sol";
    ofstream outFile(filename);
    
    if (outFile.is_open()) {
        outFile << quality << "\n";
        for (size_t i = 0; i < tour.size(); ++i) {
            outFile << tour[i] << (i < tour.size() - 1 ? "," : "");
        }
        outFile.close();
        cout << "Solution saved to " << filename << endl;
    } else {
        cerr << "Error: Unable to open file " << filename << " for writing." << endl;
    }
}

// Main function to handle input and select algorithm
int main(int argc, char* argv[]) {
    if (argc < 4) {
        cerr << "Usage: " << argv[0] << " <filename> <method> <timeLimit> [seed]" << endl;
        return 1;
    }
    
    string filename = argv[1];
    // string filename = "/Users/alyssatsai/Documents/GitHub/CSE-6140-Project-YT/DATA/DATA/Atlanta.tsp";
    string method = argv[2];
    // string method = "Approx";
    int timeLimit = atoi(argv[3]);
    // int timeLimit = atoi("60");
    int seed = argc == 5 ? atoi(argv[4]) : 0;
    // int seed = argc == 5 ? atoi("42") : 0;
    
    // Load dataset
    Dataset dataset = parseTSPFile(filename);

    // 輸出距離矩陣
/*     cout << "Distance Matrix:" << endl;
    for (const auto& row : dataset.distanceMatrix) {
        for (double dist : row) {
            cout << dist << " ";
        }
        cout << endl;
    } */
    
    // Select and run the algorithm
    vector<int> tour;
    double quality = numeric_limits<double>::max();
    
    cout << "Received method: " << method << endl;

    if (method == "BF") 
    {
        tour = exactAlgorithm(dataset, timeLimit);
        quality = 0; // Set this to the computed solution quality
    } 
    else if (method == "Approx") 
    {   
        cerr << "Method Approx" << endl;
        tour = approximateAlgorithm(dataset);
        quality = 0; // Set this to the computed solution quality
    } 
    else if (method == "LS") 
    {
        tour = localSearchAlgorithm(dataset.points, timeLimit, seed);
        quality = 0; // Set this to the computed solution quality
    } 
    else 
    {
        cerr << "Error: Unknown method " << method << endl;
        return 1;
    }
    
    // Save the result
    saveSolution(filename, method, timeLimit, seed, quality, tour);
    
    return 0;
}