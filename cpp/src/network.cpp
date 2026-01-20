#include <algorithm>
#include <functional>
#include <limits>
#include <queue>
#include <utility>

#include "network.hpp"

namespace rcmc {

void dfs2(const UndirectedNetwork &N, int v, std::vector<int> &component_eqs,
          std::vector<bool> &visited) {
    visited[v] = true;
    component_eqs.push_back(v);
    for (const auto &[u, _] : N.edges[v]) {
        if (!visited[u]) {
            dfs2(N, u, component_eqs, visited);
        }
    }
}

std::vector<Component> UndirectedNetwork::decompose() const {
    const int n = num_vertices();
    std::vector<Component> components;
    std::vector<bool> visited(n, false);

    std::cout << n << " vertices in total" << std::endl;
    for (int i = 0; i < n; ++i) {
        if (!visited[i]) {
            std::vector<int> global_to_local;
            dfs2(*this, i, global_to_local, visited);
            std::sort(global_to_local.begin(), global_to_local.end());

            std::unordered_map<int, int> local_to_global;
            for (std::size_t j = 0; j < global_to_local.size(); ++j) {
                local_to_global[global_to_local[j]] = j;
            }

            UndirectedNetwork component(global_to_local.size());
            for (std::size_t j = 0; j < global_to_local.size(); ++j) {
                component.vertices[j] = vertices[global_to_local[j]];
            }
            for (std::size_t j = 0; j < global_to_local.size(); ++j) {
                for (const auto &[k, value] : edges[global_to_local[j]]) {
                    component.set_value(j, local_to_global[k], value);
                }
            }

            components.push_back({component, global_to_local, local_to_global});
        }
    }

    return components;
}

struct Edge {
    int src;
    int dst;
    Real value;
};

std::pair<Real, std::vector<std::pair<int, Real>>>
dijkstra_max(const DirectedNetwork &N, const int src, const int dst) {
    if (src == dst) {
        return {0.0, {}};
    }

    const int n = N.num_vertices();

    std::vector<Edge> highest(
        n, Edge{-1, -1, std::numeric_limits<Real>::infinity()});

    std::priority_queue<std::pair<Real, int>, std::vector<std::pair<Real, int>>,
                        std::greater<std::pair<Real, int>>>
        Q;
    Q.emplace(0.0, src);
    highest[src] = Edge{-1, -1, 0.0};

    while (!Q.empty()) {
        const auto [value, i] = Q.top();
        Q.pop();

        if (value > highest[i].value) {
            continue;
        }
        if (i == dst) {
            break;
        }

        for (const auto &[j, v] : N.edges[i]) {
            const Edge next = value > v ? highest[i] : Edge{i, j, v};
            if (next.value < highest[j].value) {
                highest[j] = next;
                Q.emplace(next.value, j);
            }
        }
    }

    if (highest[dst].src == -1) {
        return {std::numeric_limits<Real>::infinity(), {}};
    }

    auto [tmp1, path1] = dijkstra_max(N, src, highest[dst].src);
    const auto [tmp2, path2] = dijkstra_max(N, highest[dst].dst, dst);

    path1.emplace_back(highest[dst].dst, highest[dst].value);
    path1.insert(path1.end(), path2.begin(), path2.end());

    return {highest[dst].value, path1};
}

std::pair<Real, std::vector<std::pair<int, Real>>>
dijkstra_min(const DirectedNetwork &N, const int src, const int dst) {
    if (src == dst) {
        return {std::numeric_limits<Real>::infinity(), {}};
    }

    const int n = N.num_vertices();
    std::vector<Edge> lowest(n, Edge{-1, -1, 0.0});

    std::priority_queue<std::pair<Real, int>> Q;
    Q.emplace(std::numeric_limits<Real>::infinity(), src);
    lowest[src] = Edge{-1, -1, std::numeric_limits<Real>::infinity()};

    while (!Q.empty()) {
        const auto [value, i] = Q.top();
        Q.pop();

        if (value > lowest[i].value) {
            continue;
        }
        if (i == dst) {
            break;
        }

        for (const auto &[j, v] : N.edges[i]) {
            const Edge next = value < v ? lowest[i] : Edge{i, j, v};
            if (next.value > lowest[j].value) {
                lowest[j] = next;
                Q.emplace(next.value, j);
            }
        }
    }

    if (lowest[dst].src == -1) {
        return {0.0, {}};
    }

    auto [tmp1, path1] = dijkstra_min(N, src, lowest[dst].src);
    const auto [tmp2, path2] = dijkstra_min(N, lowest[dst].dst, dst);

    path1.emplace_back(lowest[dst].dst, lowest[dst].value);
    path1.insert(path1.end(), path2.begin(), path2.end());

    return {lowest[dst].value, path1};
}
} // namespace rcmc
