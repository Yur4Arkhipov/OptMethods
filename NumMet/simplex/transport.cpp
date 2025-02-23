#include <fstream>
#include <iostream>
#include <string>
#include <vector>

struct TransportProblem {
    std::vector<double> supply;      // Supply quantities
    std::vector<double> demand;      // Demand quantities  
    std::vector<std::vector<double>> costs; // Transportation costs matrix
};

void transportToLinear(
    const TransportProblem& tp,
    std::vector<std::vector<double>>& coefficients,
    std::vector<double>& b,
    std::vector<int>& signs,
    std::vector<double>& objective
) {
    int m = tp.supply.size();    // Number of suppliers
    int n = tp.demand.size();    // Number of consumers
    int total_vars = m * n;      // Total number of variables xij
    int total_constraints = m + n; // Supply + demand constraints

    // Resize output vectors
    coefficients.resize(total_constraints, std::vector<double>(total_vars, 0.0));
    b.resize(total_constraints);
    signs.resize(total_constraints, 0); // All equations are '='
    objective.resize(total_vars);

    // Fill objective function coefficients (costs)
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            objective[i * n + j] = -tp.costs[i][j];
        }
    }

    // Supply constraints: sum(xij) = ai for each i
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            coefficients[i][i * n + j] = 1.0;
        }
        b[i] = tp.supply[i];
    }

    // Demand constraints: sum(xij) = bj for each j
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < m; i++) {
            coefficients[m + j][i * n + j] = 1.0;
        }
        b[m + j] = tp.demand[j];
    }
}

bool isBalanced(const TransportProblem& tp) {
    double totalSupply = 0.0;
    double totalDemand = 0.0;
    
    for (double s : tp.supply) totalSupply += s;
    for (double d : tp.demand) totalDemand += d;
    
    return std::abs(totalSupply - totalDemand) < 1e-10;
}

extern void open_files(std::ifstream &input_file, std::ofstream &output_file, 
                      const std::string &input_fileName, const std::string &output_fileName);
extern int simplex_main();

void solveTransportProblem(const TransportProblem& tp) {
    if (!isBalanced(tp)) {
        std::cout << "Error: Transportation problem is not balanced\n";
        return;
    }

    std::vector<std::vector<double>> coefficients;
    std::vector<double> b;
    std::vector<int> signs;
    std::vector<double> objective;

    // Convert transport problem to linear form
    transportToLinear(tp, coefficients, b, signs, objective);

    // Write the linear programming problem to input.txt
    std::ofstream input("input.txt");
    if (!input.is_open()) {
        std::cerr << "Error: Cannot create input file\n";
        return;
    }

    int m = tp.supply.size();
    int n = tp.demand.size();
    input << m + n << " " << m * n << "\n";  // num_equations num_variables

    // Write constraints
    for (int i = 0; i < coefficients.size(); i++) {
        for (int j = 0; j < coefficients[i].size(); j++) {
            input << coefficients[i][j] << " ";
        }
        input << b[i] << " " << signs[i] << "\n";
    }

    // Write objective function coefficients
    for (double c : objective) {
        input << c << " ";
    }
    input.close();

    // Call simplex solver
    simplex_main();

    // Read and interpret results
    std::ifstream output("output.txt");
    if (!output.is_open()) {
        std::cout << "Error: Cannot read results\n";
        return;
    }

    // Print transportation solution
    std::cout << "\nTransportation Problem Solution:\n";
    std::cout << "================================\n";

    std::string line;
    bool readingSolution = false;
    double objectiveValue = 0.0;
    while (std::getline(output, line)) {
        if (line.find("Objective Function Value:") != std::string::npos) {
            sscanf(line.c_str(), "Objective Function Value: %lf", &objectiveValue);
            std::cout << "Total transportation cost: " << -objectiveValue << std::endl; // Negate the value
            readingSolution = true;
            continue;
        }
        if (readingSolution) {
            if (line.find("==========") != std::string::npos) break;
            
            int var_idx;
            double value;
            if (sscanf(line.c_str(), "X%d = %lf", &var_idx, &value) == 2 && value > 0) {
                var_idx--; // Convert from 1-based to 0-based indexing
                int i = var_idx / n;  // supplier index
                int j = var_idx % n;  // consumer index
                std::cout << "From Supplier " << (i + 1) << " to Consumer " << (j + 1) 
                         << ": " << value << " units (cost: " << tp.costs[i][j] << ")\n";
            }
        }
    }
    output.close();
}

int main() {
    TransportProblem tp;
    
    tp.supply = {10, 20, 30};
    tp.demand = {15, 20, 25};
    // tp.demand = {25, 30, 25};
    
    tp.costs = {
        {5, 3, 1},
        {3, 2, 4},
        {4, 1, 2}  
    };

    solveTransportProblem(tp);
    return 0;
}