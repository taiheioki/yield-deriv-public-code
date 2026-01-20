#ifndef RCMC_NETWORK_HPP
#define RCMC_NETWORK_HPP

#include <cassert>
#include <map>
#include <numeric>
#include <utility>
#include <vector>

#include "numdef.hpp"

namespace rcmc {

struct DirectedNetwork {
    std::vector<Real> vertices;
    std::vector<std::map<int, Real>> edges;

    explicit DirectedNetwork(const int n) : vertices(n), edges(n) {}

    int num_vertices() const { return vertices.size(); }

    int num_edges() const {
        return std::accumulate(
            edges.begin(), edges.end(), 0,
            [](const int acc, const auto &es) { return acc + int(es.size()); });
    }

    // Parallel edges are not allowed.
    void set_value(const int i, const int j, const Real value) {
        assert(0 <= i && i < num_vertices());
        assert(0 <= j && j < num_vertices());

        edges[i][j] = value;
    }

    Real get_value(const int i, const int j) const {
        assert(0 <= i && i < num_vertices());
        assert(0 <= j && j < num_vertices());
        const auto itr = edges[i].find(j);
        return itr == edges[i].end() ? 0.0 : itr->second;
    }
};

struct Component;

struct UndirectedNetwork : public DirectedNetwork {
    using DirectedNetwork::DirectedNetwork;
    // using DirectedNetwork::get_value;

    int num_edges() const { return DirectedNetwork::num_edges() / 2; }

    void set_value(const int i, const int j, const Real value) {
        DirectedNetwork::set_value(i, j, value);
        DirectedNetwork::set_value(j, i, value);
    }

  private:
    void dfs(int v, std::vector<bool> &visited) const {
        visited[v] = true;
        for (const auto &[u, _] : edges[v]) {
            if (!visited[u]) {
                dfs(u, visited);
            }
        }
    }

  public:
    int num_components() const {
        std::vector<bool> visited(num_vertices(), false);
        int count = 0;

        for (int i = 0; i < num_vertices(); ++i) {
            if (!visited[i]) {
                ++count;
                dfs(i, visited);
            }
        }
        return count;
    }

  public:
    std::vector<Component> decompose() const;
};

struct Component {
    UndirectedNetwork network;
    std::vector<int> global_to_local;
    std::unordered_map<int, int> local_to_global;
};

std::pair<Real, std::vector<std::pair<int, Real>>>
dijkstra_max(const DirectedNetwork &N, const int src, const int dst);

std::pair<Real, std::vector<std::pair<int, Real>>>
dijkstra_min(const DirectedNetwork &N, const int src, const int dst);
} // namespace rcmc

#endif
