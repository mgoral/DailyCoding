/**
 * Copyright (C) Michal Goral, 2014
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef ALGO_GRAPHSEARCH_HPP_
#define ALGO_GRAPHSEARCH_HPP_

#include <vector>
#include <functional>

#include "Graph.hpp"
#include "detail/GraphSearchHelpers.hpp"
#include "detail/Node.hpp"

namespace algo
{

/**
 * Just a helper class containing useful predicates.
 * Available publicly so clients might use it too!
 */
template <typename GraphType>
struct GraphPred
{
    static bool nodeTrue(const typename GraphType::NodeType&) { return true; }
    static bool nodeFalse(const typename GraphType::NodeType&) { return false; }
};

/**
 * The same as GraphPred but Graph Key and Value can be used instead of full Graph type.
 * In some cases it might be easier to type.
 */
template <typename Key, typename Value>
struct GraphPredV : public GraphPred<algo::Graph<Key, Value>> { };

/**
 * findShortestPathBfs
 *
 * Finds the shortest path from the node indexed with a start key until a given predicate returns
 * true. In BFS the shortest path is counted by the number of edges (i.e. paths between edges are
 * assumed to be the same length). If it's not the case, use different algorithm, like Dijkstra's or
 * A*.
 *
 * @param graph: graph that should be searched
 * @param startPoint: index to the node stored in graph from which search is started
 * @param endNodeMatch: condition which, when met, stops the search
 * @param frontierMatch: when this condition is met, node is not added to frontier, i.e. search
 *        won't go through it (creates a "border" in a graph)
 * @param stopMatch: stop graph search when given condition is met. Note however that BFS doesn't
 *        stop immediately. It will first evaluate all frontiers on depth where stopMatch returned
 *        true.
 * @return vector of references to data stored in nodes.
 */
template <typename GraphType, typename EndNodePred, typename FrontierPred, typename StopPred>
std::vector<std::reference_wrapper<typename GraphType::ValueType>> findShortestPathBfs(
    GraphType& graph,
    const typename GraphType::KeyType& startPoint,
    EndNodePred endNodeMatch,
    FrontierPred frontierMatch,
    StopPred stopMatch)
{
    std::vector<std::reference_wrapper<typename GraphType::ValueType>> ret;
    auto endNodeCallback = [&ret](typename GraphType::NodeType* foundNode)
                    { ret = detail::reconstructPath(foundNode); return true; };
    detail::bfsImpl(graph, startPoint, endNodeMatch, frontierMatch, stopMatch, endNodeCallback);
    return ret;
}

template <typename GraphType, typename EndNodePred>
std::vector<std::reference_wrapper<typename GraphType::ValueType>> findShortestPathBfs(
    GraphType& graph,
    const typename GraphType::KeyType& startPoint,
    EndNodePred endNodeMatch)
{
    auto neverStop = GraphPred<GraphType>::nodeFalse;
    auto alwaysFrontier = GraphPred<GraphType>::nodeTrue;
    return findShortestPathBfs(graph, startPoint, endNodeMatch, alwaysFrontier, neverStop);
}

/**
 * findMatchingNodesBfs
 *
 * Finds all nodes matching a predicate in a graph. It's particulary useful with stopMatch because
 * knowing a simple stop search you avoid searching all graph elements.
 *
 * Params are the same as in findShortestPathBfs
 */
template <typename GraphType, typename EndNodePred, typename FrontierPred, typename StopPred>
std::vector<std::reference_wrapper<typename GraphType::ValueType>> findMatchingNodesBfs(
    GraphType& graph,
    const typename GraphType::KeyType& startPoint,
    EndNodePred endNodeMatch,
    FrontierPred frontierMatch,
    StopPred stopMatch)
{
    std::vector<std::reference_wrapper<typename GraphType::ValueType>> ret;
    auto endNodeCallback = [&ret](typename GraphType::NodeType* foundNode)
                    { ret.push_back(foundNode->data()); return false; };
    detail::bfsImpl(graph, startPoint, endNodeMatch, frontierMatch, stopMatch, endNodeCallback);
    return ret;
}

template <typename GraphType, typename Predicate>
std::vector<std::reference_wrapper<typename GraphType::ValueType>> findMatchingNodesBfs(
    GraphType& graph,
    const typename GraphType::KeyType& startPoint,
    Predicate endNodeMatch)
{
    auto neverStop = GraphPred<GraphType>::nodeFalse;
    auto alwaysFrontier = GraphPred<GraphType>::nodeTrue;
    return findMatchingNodesBfs(graph, startPoint, endNodeMatch, alwaysFrontier, neverStop);
}

} // namespace algo

#endif // ALGO_GRAPHSEARCH_HPP_
