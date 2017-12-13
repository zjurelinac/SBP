#include "evolutionary.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>

using dist_t = double;
using dist_matrix = std::vector<std::vector<dist_t>>;

#define POPULATION_SIZE 20
#define NUM_ITERATIONS 1000

#define SCHOOL 0

template <typename T>
std::ostream& operator<<(std::ostream& os, std::vector<T> v) {
    os << "< ";
    for (const auto& e : v) os << e << " ";
    os << ">";
    return os;
}

struct point {
    dist_t x, y;
    point() {}
    point(dist_t x, dist_t y)
        : x(x), y(y) {}
    bool operator<(const point& pt) const
        { return x == pt.x ? y < pt.y : x < pt.x; }
};

inline dist_t pt_dist(const point& a, const point& b) {
    return sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y));
}

struct distance_fitness {
    using fitness_type = double;

    distance_fitness(int c, std::vector<int>& stc, dist_matrix& ssd)
        : C(c), STC(stc), SSD(ssd) {}
    dist_t operator()(std::vector<int> route, bool show_extra_cap = false) {
        dist_t dist = 0.0;
        int last_stop = SCHOOL;
        int studs = 0;

        for (auto stop : route) {
            if (studs + STC[stop] > C) {
                dist += SSD[last_stop][SCHOOL] + SSD[SCHOOL][stop];
                studs = STC[stop];

                if (show_extra_cap)
                    printf("Move : %d -> %d, returning to SCHOOL with extra capacity: %d\n", last_stop, stop, C - studs);
            } else {
                dist += SSD[last_stop][stop];
                studs += STC[stop];
            }

            last_stop = stop;
        }

        dist += SSD[last_stop][SCHOOL];

        return -dist;   // lesser distance == greater fitness
    }

    private:
        int C;
        std::vector<int>& STC;
        dist_matrix& SSD;
};

auto best_city_route(int C, std::vector<int>& STC, dist_matrix& SSD) {
    std::vector<int> base;
    for (auto i = 0u; i < STC.size(); ++i)
        if (STC[i] > 0)
            base.push_back(i);

    for (auto& x : base) printf("%d ", x); puts("");

    auto best = ea::genetic_algorithm(ea::generator_initializer<ea::permutation_shuffle_generator<std::vector<int>>>(ea::permutation_shuffle_generator<std::vector<int>>(base), POPULATION_SIZE),
                                      ea::cross_mutation_reproducer<ea::rank_parent_selector, ea::ordered_crossover, ea::swap_mutator>(ea::rank_parent_selector(POPULATION_SIZE), ea::ordered_crossover(), ea::swap_mutator(), 2),
                                      ea::best_n_selector(POPULATION_SIZE),
                                      ea::const_num_iter_terminator(NUM_ITERATIONS),
                                      distance_fitness(C, STC, SSD));

    return best;
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        puts("Usage: sbr infile");
        return 1;
    }

    freopen(argv[1], "r", stdin);

    int M, N, C;
    dist_t MAX_WALK, x, y;

    // Instance data loading

    scanf("%d stops, %d students, %lf maximum walk, %d capacity\n", &M, &N, &MAX_WALK, &C);

    std::vector<point> stops, students;

    for (int i = 0; i < M; ++i) {
        scanf("%*d %lf %lf", &x, &y);
        stops.emplace_back(x, y);
    }

    for (int i = 0; i < N; ++i) {
        scanf("%*d %lf %lf", &x, &y);
        students.emplace_back(x, y);        // IMPORTANT: original student indexes 1-based
    }

    // Distance preprocessing

    dist_matrix SSD(M);   // Stop-Stop-Distance

    for (int i = 0; i < M; ++i)
        SSD[i].resize(M, 0);

    for (int i = 0; i < M; ++i)
        for (int j = i + 1; j < M; ++j)
            SSD[i][j] = SSD[j][i] = pt_dist(stops[i], stops[j]);

    std::vector<int> SS(M-1);               // Sorted-Stops

    for (int i = 1; i < M; ++i)
        SS[i-1] = i;

    std::stable_sort(SS.begin(), SS.end(), [SSD](int x, int y){ return SSD[0][x] < SSD[0][y]; });

    std::vector<std::vector<int>> TNS(N);   // sTudent-Nearby-Stops
    std::vector<std::vector<int>> SNT(M);   // Stop-Nearby-sTudents

    for (int i = 0; i < N; ++i)
        for (auto j : SS)
            if (pt_dist(students[i], stops[j]) < MAX_WALK) {
                TNS[i].push_back(j);
                SNT[j].push_back(i);
            }

    // Initial stop assignment - greedy

    std::vector<int> STC(M, 0);             // Stop-sTudent-Count
    std::vector<int> TAS(N, -1);            // sTudent-Assigned-Stop

    for (auto s : SS) printf("%d, ", s); puts("");

    for (auto j : SS)
        for (auto i : SNT[j]) {
            if (STC[j] >= C)
                break;
            else if (TAS[i] == -1) {
                TAS[i] = j;
                ++STC[j];
            }
        }

    for (int i = 0; i < N; ++i) printf("%d %d\n", i + 1, TAS[i]); puts("");

    for (int j = 0; j < M; ++j) printf("%d %d\n", j, STC[j]); puts("");

    // Actual problem solving
    auto best = best_city_route(C, STC, SSD);
    printf("Best solution: %.8lf\n", best.first);
    std::cout << best.second << "\n";

    auto df = distance_fitness(C, STC, SSD);
    df(best.second, true);

    return 0;
}