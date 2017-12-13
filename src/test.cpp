// #include <cstdlib>
// #include <cstdio>
// #include <ctime>
// #include <vector>

// #include "evolutionary.hpp"

// /* test problem - find maximum of a parabola  */

// struct generator {
//     using solution_type = double;
//     solution_type operator()() { return (rand() % 1000) / 100.0 - 5.0; }
// };

// struct fitness_evaluator {
//     using fitness_type = double;
//     double operator()(double x) { return 10 + 7*x + 4*x*x - x*x*x - 0.5*x*x*x*x; }
// };

// struct mutator {
//     double operator()(double x) { return x + ((rand() % 10000000) / 100000000.0 - 0.05); }
// };

// struct crosser {
//     auto operator()(double x, double y) { return std::vector<double> {(x + y) / 2}; }
// };

// int main() {
//     srand(time(NULL));
//     auto best = ea::genetic_algorithm(ea::generator_initializer<generator>(generator(), 10),
//                                       ea::cross_mutation_reproducer<ea::rank_parent_selector, crosser, mutator>(ea::rank_parent_selector(16), crosser(), mutator(), 2),
//                                       //ea::mutation_reproducer<mutator>(mutator(), 2),
//                                       ea::best_n_selector(10),
//                                       ea::const_num_iter_terminator(10000),
//                                       fitness_evaluator());

//     printf("Best solution: %.8lf <%.8lf>\n", best.second, best.first);
// }