#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <iterator>
#include <limits>
#include <random>
#include <set>
#include <utility>
#include <vector>

#include "ctpl.h"

namespace {
    template <typename random_engine = std::default_random_engine>
    struct randomizer {
        randomizer()
            : seed(time(NULL)), engine(seed) { srand(seed); }
        int operator()(int min, int max)
            { return rand() % (max - min + 1) + min; }

        unsigned seed;
        random_engine engine;
    };

    randomizer<> rng;
}

namespace ea {
    template <typename T>
    using container = std::vector<T>;

    template <typename fitness_type, typename solution_type>
    using entity = std::pair<fitness_type, solution_type>;

    /*** Utility functions ***/

#ifdef MULTITHREADED
    template <typename solution_container, typename fitness_policy>
    inline auto evaluate(const solution_container& sc, fitness_policy& evaluator) {
        static ctpl::thread_pool pool(std::thread::hardware_concurrency());


    }
#else
    template <typename solution_container, typename fitness_policy>
    inline auto evaluate(const solution_container& sc, fitness_policy& evaluator) {
        using solution_entity = entity<typename fitness_policy::fitness_type, typename solution_container::value_type>;

        container<solution_entity> solutions;

        for (const auto& sol : sc)
            solutions.emplace_back(evaluator(sol), sol);

        return solutions;
    }
#endif

    template <typename T>
    std::ostream& operator<<(std::ostream& os, std::vector<T> v) {
        os << "< ";
        for (const auto& e : v) os << e << " ";
        os << ">";
        return os;
    }

    template <typename solution_container>
    inline void show(const solution_container& sc) {
        for (auto& x : sc) std::cout << "[" << x.first << ", " << x.second << "] ";
        std::cout << "\n";
    }

    /*** Genetic algorithm itself ***/

    template <typename initializer_policy, typename reproduction_policy,
              typename selection_policy, typename termination_policy, typename fitness_policy>
    auto genetic_algorithm(initializer_policy initializer, reproduction_policy reproducer,
                           selection_policy selector, termination_policy terminator, fitness_policy evaluator) {
        using solution_entity = entity<typename fitness_policy::fitness_type, typename initializer_policy::solution_type>;

        container<solution_entity> solutions = evaluate(initializer(), evaluator);

        while (!terminator(solutions)) {
            solutions = selector(evaluate(reproducer(solutions), evaluator));
            //show(solutions);
        }

        return *(std::max_element(solutions.begin(), solutions.end()));
    }

    /*** Common population initializers ***/

    template <typename solution_type>
    struct const_initializer {
        const_initializer(std::initializer_list<solution_type> cnsts)
            : consts(cnsts.begin(), cnsts.end()) {}
        const_initializer(const const_initializer& ci)
            : consts(ci.consts) {}
        const_initializer(const_initializer&& ci)
            : consts(std::move(ci.consts)) {}
        const std::vector<solution_type>& operator()() const { return consts; }
    private:
        const std::vector<solution_type> consts;
    };

    template <typename generator>
    struct generator_initializer {
        using solution_type = typename generator::solution_type;

        generator_initializer(generator g, int n)
            { for (int i = 0; i < n; ++i) consts.push_back(g()); }
        generator_initializer(const generator_initializer& ci)
            : consts(ci.consts) {}
        generator_initializer(generator_initializer&& ci)
            : consts(std::move(ci.consts)) {}
        const std::vector<solution_type>& operator()() const { return consts; }
    private:
        std::vector<solution_type> consts;
    };

    /*** Common reproducers ***/

    template <typename mutation_policy>
    struct mutation_reproducer {
        mutation_reproducer(int en = 0)
            : elite_num(en) {}
        mutation_reproducer(mutation_policy mp, int en = 0)
            : elite_num(en), mutator(mp) {}
        template <typename solution_container>
        const std::vector<typename solution_container::value_type::second_type> operator()(solution_container& sc) {
            std::vector<typename solution_container::value_type::second_type> offspring(elite_num);

            std::partial_sort(sc.begin(), sc.end(), sc.begin() + elite_num,
                              std::greater<typename solution_container::value_type>());

            std::transform(sc.begin(), sc.begin() + elite_num, offspring.begin(),
                           [](auto x){ return x.second; });

            for (const auto& sol : sc)
                offspring.push_back(mutator(sol.second));
            return offspring;
        }
    private:
        int elite_num;
        mutation_policy mutator;
    };

    template <typename parent_select_policy, typename crossover_policy, typename mutation_policy>
    struct cross_mutation_reproducer {
        cross_mutation_reproducer(int en = 0)
            : elite_num(en) {}
        cross_mutation_reproducer(parent_select_policy psp, crossover_policy cp, mutation_policy mp, int en = 0)
            : elite_num(en), parent_selector(psp), crosser(cp), mutator(mp) {}
        template <typename solution_container>
        auto operator()(solution_container& sc) {
            std::vector<typename solution_container::value_type::second_type> offspring(elite_num);

            std::stable_sort(sc.begin(), sc.end(), std::greater<typename solution_container::value_type>());

            std::transform(sc.begin(), sc.begin() + elite_num, offspring.begin(),
                           [](auto x){ return x.second; });

            for (const auto parents : parent_selector(sc))
                for (auto child : crosser(parents.first, parents.second))
                    offspring.push_back(mutator(child));

            return offspring;
        }
    private:
        int elite_num;
        parent_select_policy parent_selector;
        crossover_policy crosser;
        mutation_policy mutator;
    };

    /*** Common crossover parent selectors ***/

    struct rank_parent_selector {
        rank_parent_selector(unsigned num)
            : n(num) {}
        template <typename solution_container>
        auto operator()(const solution_container& sc) {
            using solution_type = typename solution_container::value_type::second_type;
            std::vector<std::pair<solution_type, solution_type>> parents;

            auto sz = sc.size();

            for (auto i = 0u; i < n; ++i) {
                int x = select(sz), y;
                while ((y = select(sz)) == x) {}
                parents.emplace_back(sc[x].second, sc[y].second);
            }

            return parents;
        }
    private:
        int select(int max)
            { return max - static_cast<int>((sqrt(8 * rng(0, max*(max + 1) / 2 - 1) + 1) - 1) / 2.0) - 1; }
        unsigned n;
    };

    /*** Common population selectors ***/

    struct best_n_selector {
        best_n_selector(unsigned num)
            : n(num) {}
        template <typename solution_container>
        solution_container operator()(solution_container sc) {
            if (n < sc.size()) {
                std::partial_sort(sc.begin(), sc.begin() + n, sc.end(),
                                  std::greater<typename solution_container::value_type>());
                sc.resize(n);
            }
            return sc;
        }
    private:
        unsigned n;
    };

    /*** Common algorithm terminators ***/

    struct const_num_iter_terminator {
        const_num_iter_terminator(unsigned num)
            : num_iters(num), curr_iter(0) {}
        template <typename solution_container>
        bool operator()(const solution_container&)
            { return ++curr_iter >= num_iters; }
    private:
        unsigned num_iters, curr_iter;
    };

    struct finite_time_terminator {
        finite_time_terminator(unsigned dur)
            : duration(dur) { start_time = std::time(nullptr); }
        template <typename solution_container>
        bool operator()(const solution_container&)
            { return std::time(nullptr) - start_time >= duration; }
    private:
        unsigned duration;
        std::time_t start_time;
    };

    template <typename fitness_t>
    struct const_max_fitness_terminator {
        const_max_fitness_terminator(unsigned mui)
            : max_unchanged_iters(mui), curr_unchanged_iters(0), max_fitness(std::numeric_limits<fitness_t>::min()) {}
        template <typename solution_container>
        bool operator()(const solution_container& sc) {
            fitness_t new_max = std::max_element(sc.begin(), sc.end())->first;
            if (new_max > max_fitness) {
                max_fitness = new_max;
                curr_unchanged_iters = 0;
            } else {
                ++curr_unchanged_iters;
            }

            return curr_unchanged_iters >= max_unchanged_iters;
        }
    private:
        unsigned max_unchanged_iters, curr_unchanged_iters;
        fitness_t max_fitness;
    };

    /*** Common initial solution generator functors ***/

    template <typename permutation_solution>
    struct permutation_shuffle_generator {
        using solution_type = permutation_solution;

        permutation_shuffle_generator(const permutation_solution& ps)
            : base(ps) {}
        permutation_solution operator()()
            { permutation_solution ps(base); std::shuffle(ps.begin(), ps.end(), rng.engine); return ps; }
        private:
            permutation_solution base;
    };

    /*** Common crossover operators ***/

    struct ordered_crossover {
        template <typename solution_type>
        std::vector<solution_type> operator()(const solution_type& x, const solution_type& y) {
            auto n = x.size();
            auto max = *std::max_element(x.begin(), x.end());

            int c1 = rng(0, n-2);
            int c2 = rng(c1 + 1, n-1);

            solution_type z(n, 0), w(n, 0);
            std::vector<bool> z_occ(max + 1), w_occ(max + 1);

            for (auto i = c1; i < c2; ++i) {
                z[i] = x[i];
                w[i] = y[i];
                z_occ[x[i]] = w_occ[y[i]] = true;
            }

            int zpos = c2, wpos = c2;
            for (auto i = 0u; i < n; ++i) {
                if (!z_occ[y[i]]) {
                    z[zpos] = y[i];
                    zpos = (zpos + 1) % n;
                }
                if (!w_occ[x[i]]) {
                    w[wpos] = x[i];
                    wpos = (wpos + 1) % n;
                }
            }

            return std::vector<solution_type> {z, w};
        }
    };

    /*** Common mutation operators ***/

    struct swap_mutator {
        template <typename solution_type>
        solution_type operator()(solution_type& s)
            { std::iter_swap(s.begin() + rng(0, s.size() - 1), s.begin() + rng(0, s.size() - 1)); return s; }
    };
}
