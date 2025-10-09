#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>

using std::cout; using std::string; using std::to_string; using std::setw; using std::left;
class Configuration {
public:
    Configuration() = default; 

    int dimension = 10;
    string representation = "binary";
    int population_size = 100 ;
    bool uniform_crossover = false;
    string crossover_method = "2-point";
    float cross_prob = 0.9;
    float mut_prob = 0.1;
    int num_generations = 500 ; 
    bool debug = false;

    int tournament_size = 2 ; 
    double alpha = 1.0 ;
    double beta = 1.0 ;
    double gamma = 1.0 ;
    void print_configuration() const {
        auto configurations = get_configurations();
        const int column_width = 20;
        const int num_columns = 2;
        const int frame_width = num_columns * column_width + 3;
        
        cout << string(frame_width, '-') << "\n";
        cout << '|' << left << setw(column_width) << "Parameter"
             << '|' << left << setw(column_width) << "Value" << '|' << "\n";
        cout << string(frame_width, '-') << "\n";
        
        for (const auto& [key, value] : configurations) {
            cout << '|' << left << setw(column_width) << key 
                 << '|' << left << setw(column_width) << value << "|\n";
        }
        cout << string(frame_width, '-') << "\n";
    }

private:
    std::vector<std::pair<string, string>> get_configurations() const {
        string cross_prob_num = to_string(cross_prob);
        string cross_prob_rounded = cross_prob_num.substr(0, cross_prob_num.find(".")+2);

        string mut_prob_num = to_string(mut_prob);
        string mut_prob_rounded = mut_prob_num.substr(0, mut_prob_num.find(".")+2);
      
        return {
            {"dimension", to_string(dimension)},
            {"representation", representation},
            {"population_size", to_string(population_size)},
            {"uniform_crossover", uniform_crossover ? "true" : "false"},
            {"crossover_method", crossover_method},
            {"cross_prob", cross_prob_rounded},
            {"mut_prob", mut_prob_rounded},
            {"num_generations", to_string(num_generations)},
            // {"debug", (debug) ? "true" : "false"},
        };
    }
};
