#include "evolutionary.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <stack>
#include <vector>

#define MAX_FIT_CHANGE_ITERS 50000

using dist_t = float;
using dist_matrix = std::vector<std::vector<dist_t>>;

#define SCHOOL 0
#define INF 1e18f

static unsigned int FITNESS_EVALS = 0;

template <typename T>
std::ostream& operator<<(std::ostream& os, std::vector<T> v) {
    os << "< ";
    for (const auto& e : v) os << (int) e << " ";
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
    using fitness_type = dist_t;

    distance_fitness(int c, std::vector<int>& stc, dist_matrix& ssd)
        : C(c), stop_students_count(stc), stop_stop_dist(ssd) {}
    dist_t operator()(std::vector<int> route) {
        ++FITNESS_EVALS;
        dist_t dist = 0.0;
        int last_stop = SCHOOL;
        int studs = 0;

        for (auto stop : route) {
            if (studs + stop_students_count[stop] > C) {
                dist += stop_stop_dist[last_stop][SCHOOL] + stop_stop_dist[SCHOOL][stop];
                studs = stop_students_count[stop];
            } else {
                dist += stop_stop_dist[last_stop][stop];
                studs += stop_students_count[stop];
            }

            last_stop = stop;
        }

        dist += stop_stop_dist[last_stop][SCHOOL];

        return -dist;   // lesser distance == greater fitness
    }

    void reconstruct(std::vector<int> route) {
        int studs = 0;
        for (auto stop : route) {
            if (studs + stop_students_count[stop] > C) {
                printf("\n%d ", stop);
                studs = stop_students_count[stop];
            } else {
                printf("%d ", stop);
                studs += stop_students_count[stop];
            }
        }
        puts("\n");
    }

    private:
        int C;
        std::vector<int>& stop_students_count;
        dist_matrix& stop_stop_dist;
};

struct dynamic_distance_fitness {
    using fitness_type = dist_t;

    dynamic_distance_fitness(int c, std::vector<int>& stc, dist_matrix& ssd)
        : C(c), stop_students_count(stc), stop_stop_dist(ssd) {}
    dist_t operator()(std::vector<int> route) {
        ++FITNESS_EVALS;
        this->route = route;
        this->N = route.size();

        DP.clear(); DP.resize(N*N, INF);
        auto sol = stop_stop_dist[SCHOOL][route[0]] + solve(0, 0, 0);

        return -sol;
    }

    dist_t solve(int i, int l, int L) {
        if (DP[i*N + l] != INF) return DP[i*N + l];  // If this state was already calculated
        if (i == N-1) return DP[i*N + l] = stop_stop_dist[route[i]][SCHOOL];    // If reached the last stop on the route
        L += stop_students_count[route[i]];          // Pick all students at the stop

        dist_t sol_direct = (L + stop_students_count[route[i+1]] <= C)
            ? solve(i+1, l+1, L) + stop_stop_dist[route[i]][route[i+1]] : INF;  // Go directly to the next stop if possible
        dist_t sol_school = solve(i+1, 0, 0) + stop_stop_dist[route[i]][SCHOOL] + stop_stop_dist[SCHOOL][route[i+1]];  // Return to school and then continue
        return DP[i*N + l] = std::min(sol_direct, sol_school);
    }

    void reconstruct(std::vector<int> route) {
        (*this)(route);

        int l = 0, L = 0;
        for (int i = 0; i < N; ++i) {
            L += stop_students_count[route[i]];
            printf("%d ", route[i]);
            dist_t sol_direct = (L + stop_students_count[route[i+1]] <= C)
                ? DP[(i+1)*N + l+1] + stop_stop_dist[route[i]][route[i+1]] : INF;
            dist_t sol_school = DP[(i+1)*N] + stop_stop_dist[route[i]][SCHOOL] + stop_stop_dist[SCHOOL][route[i+1]];

            if (sol_direct <= sol_school) {
                ++l;
            } else {
                L = l = 0;
                putchar('\n');
            }
        }

        putchar('\n');
    }

    private:
        int C;
        std::vector<int>& stop_students_count;
        dist_matrix& stop_stop_dist;

        int N;
        bool *R;
        std::vector<dist_t> DP;
        std::vector<int> route;
};

auto best_city_route(int C, std::vector<int>& stop_students_count, dist_matrix& stop_stop_dist) {
    std::vector<int> base;
    for (auto i = 0u; i < stop_students_count.size(); ++i)
        if (stop_students_count[i] > 0)
            base.push_back(i);

    return ea::genetic_algorithm(ea::generator_initializer<ea::permutation_shuffle_generator<std::vector<int>>>(ea::permutation_shuffle_generator<std::vector<int>>(base), POPULATION),
                                 ea::cross_mutation_reproducer<ea::rank_parent_selector, ea::ordered_crossover, ea::swap_mutator>(ea::rank_parent_selector(POPULATION), ea::ordered_crossover(), ea::swap_mutator(), ELITENESS),
                                 ea::best_n_selector(POPULATION),

#if defined(TERMINATE_ITERS)
                                 ea::const_num_iter_terminator(TERMINATE_ITERS),
#elif defined(TERMINATE_TIME)
                                 ea::finite_time_terminator(TERMINATE_TIME),
#else /* Terminate on max fitness reached */
                                 ea::const_max_fitness_terminator<dist_t>(MAX_FIT_CHANGE_ITERS),
#endif

#ifdef DYNAMIC_FITNESS
                                 dynamic_distance_fitness(C, stop_students_count, stop_stop_dist)
#else
                                 distance_fitness(C, stop_students_count, stop_stop_dist)
#endif
                                );
}

template <typename solution_type>
void reconstruct(solution_type best, int N, int C, std::vector<int>& stop_students_count, std::vector<int>& student_chosen_stop, dist_matrix& stop_stop_dist) {
#ifdef DYNAMIC_FITNESS
    dynamic_distance_fitness ddf(C, stop_students_count, stop_stop_dist);
    ddf.reconstruct(best.second);
#else
    distance_fitness df(C, stop_students_count, stop_stop_dist);
    df.reconstruct(best.second);
#endif
    for (int i = 0; i < N; ++i) printf("%d %d\n", i + 1, student_chosen_stop[i]); puts("");
}

struct nearest_students_greedy {
    nearest_students_greedy(int M, int N, int C, std::vector<int>& sorted_stops, std::vector<std::vector<int>>& stop_nearby_students, dist_matrix& stop_stop_dist)
        : M(M), N(N), C(C), sorted_stops(sorted_stops), stop_nearby_students(stop_nearby_students), stop_stop_dist(stop_stop_dist), stop_students_count(M, 0), student_chosen_stop(N, -1) {}

    void assign() {
        for (auto s : sorted_stops) {
            int observed = 0, occupied = 0;
            for (auto t : stop_nearby_students[s]) {
                ++observed;
                if (stop_students_count[s] >= C) break;
                if (student_chosen_stop[t] == -1) {
                    student_chosen_stop[t] = s;
                    ++stop_students_count[s];
                } else {
                    ++occupied;
                }
            }
        }
    }

    auto calculate_route()
        { return best_city_route(C, stop_students_count, stop_stop_dist); }

    int M, N, C;
    std::vector<int>& sorted_stops;
    std::vector<std::vector<int>>& stop_nearby_students;
    dist_matrix& stop_stop_dist;

    std::vector<int> stop_students_count;
    std::vector<int> student_chosen_stop;
};

struct least_stops_greedy {
    least_stops_greedy(int M, int N, int C, std::vector<int>& sorted_stops, std::vector<std::vector<int>>& stop_nearby_students, dist_matrix& stop_stop_dist)
        : M(M), N(N), C(C), sorted_stops(sorted_stops), stop_nearby_students(stop_nearby_students), stop_stop_dist(stop_stop_dist), stop_students_count(M, 0), student_chosen_stop(N, -1) {}

    void assign() {
        int unassigned = N;
        while (unassigned > 0) {
            int best_stop = -1, best_potential = 0;
            for (auto s : sorted_stops) {
                if (stop_students_count[s] > 0) continue;       // If some students already assigned to this stop, skip it

                int potential = 0;
                bool filled = false;
                for (auto t : stop_nearby_students[s]) {
                    if (student_chosen_stop[t] != -1) continue; // If student already assigned somewhere, skip him
                    if (++potential >= C) {                       // If stop filled to its capacity, break
                        filled = true;
                        break;
                    }
                }

                if (filled) {                                   // If some stop can be filled, do it
                    best_stop = s;
                    best_potential = C;
                    break;
                } else if (potential > best_potential) {        // Else, choose the one with most possible students
                    best_stop = s;
                    best_potential = potential;
                }
            }

            int assigned = 0;
            for (auto t : stop_nearby_students[best_stop]) {    // Record the assignments
                if (student_chosen_stop[t] != -1) continue;
                student_chosen_stop[t] = best_stop;
                if (++assigned >= C) break;
            }
            stop_students_count[best_stop] = assigned;

            unassigned -= assigned;
        }
    }

    auto calculate_route() { return best_city_route(C, stop_students_count, stop_stop_dist); }

    int M, N, C;
    std::vector<int>& sorted_stops;
    std::vector<std::vector<int>>& stop_nearby_students;
    dist_matrix& stop_stop_dist;

    std::vector<int> stop_students_count;
    std::vector<int> student_chosen_stop;
};

int main(int argc, char* argv[]) {
    if (argc != 2) {
        puts("Usage: sbr infile");
        return 1;
    }

    freopen(argv[1], "r", stdin);

    int M, N, C;
    dist_t MAX_WALK, x, y;


    /*** Instance data loading ***/

    scanf("%d stops, %d students, %f maximum walk, %d capacity\n", &M, &N, &MAX_WALK, &C);

    std::vector<point> stops, students;

    for (int i = 0; i < M; ++i) {
        scanf("%*d %f %f", &x, &y);
        stops.emplace_back(x, y);
    }

    for (int i = 0; i < N; ++i) {
        scanf("%*d %f %f", &x, &y);
        students.emplace_back(x, y);    // IMPORTANT: original student indexes 1-based
    }


    /*** Distance preprocessing ***/

    dist_matrix stop_stop_dist(M);      // Distances between two stops
    dist_matrix stop_student_dist(M);   // Distances between stops and students

    for (int s = 0; s < M; ++s) {
        stop_stop_dist[s].resize(M, 0);
        stop_student_dist[s].resize(N, 0);
    }

    // Stop-stop distance calculation
    for (int s1 = 0; s1 < M; ++s1)
        for (int s2 = s1 + 1; s2 < M; ++s2)
            stop_stop_dist[s1][s2] = stop_stop_dist[s2][s1] = pt_dist(stops[s1], stops[s2]);

    // Stop-student distance calculation
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < N; ++j)
            stop_student_dist[i][j] = pt_dist(stops[i], students[j]);


    /*** Stops sorting to determine their priority (closest to the school comes first) ***/

    std::vector<int> sorted_stops(M-1);        // Sorted stop indices

    for (int i = 1; i < M; ++i)
        sorted_stops[i-1] = i;

    std::stable_sort(sorted_stops.begin(), sorted_stops.end(),
        [stop_stop_dist](auto x, auto y){ return stop_stop_dist[0][x] < stop_stop_dist[0][y]; });


    /*** Matching stops to nearby students and the opposite***/

    std::vector<std::vector<int>> stop_nearby_students(M);   // Students near the stop

    for (auto s : sorted_stops)
        for (int t = 0; t < N; ++t)
            if (stop_student_dist[s][t] <= MAX_WALK)
                stop_nearby_students[s].push_back(t);

    for (int s = 1; s < M; ++s)
        std::stable_sort(stop_nearby_students[s].begin(), stop_nearby_students[s].end(),
            [s, stop_student_dist](auto x, auto y){ return stop_student_dist[s][x] < stop_student_dist[s][y]; });


    /*** Initial stop assignment - try two greedys ***/

    std::cerr << "Done preprocessing\n";

    auto nearest = nearest_students_greedy(M, N, C, sorted_stops, stop_nearby_students, stop_stop_dist);
    nearest.assign();
    auto best_nearest = nearest.calculate_route();
    std::cerr << "Done greedy nearest :: " << best_nearest.first << "\n";

    auto least = least_stops_greedy(M, N, C, sorted_stops, stop_nearby_students, stop_stop_dist);
    least.assign();
    auto best_least = least.calculate_route();
    std::cerr << "Done greedy least :: " << best_least.first << "\n";

    /*** Reconstruct the best route and output it ***/

    if (best_nearest.first > best_least.first) {
        reconstruct(best_nearest, N, C, nearest.stop_students_count, nearest.student_chosen_stop, stop_stop_dist);
        std::cerr << -best_nearest.first << "\n";
    } else {
        reconstruct(best_least, N, C, least.stop_students_count, least.student_chosen_stop, stop_stop_dist);
        std::cerr << -best_least.first << "\n";
    }
    std::cerr << "Fitness evals = " << FITNESS_EVALS << "\n";

    return 0;
}