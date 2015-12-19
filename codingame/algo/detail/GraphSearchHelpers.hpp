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

#ifndef ALGO_DETAIL_GRAPHSEARCHHELPERS_HPP_
#define ALGO_DETAIL_GRAPHSEARCHHELPERS_HPP_

#include <functional>
#include <algorithm>
#include <limits>
#include <queue>
#include <stdexcept>

namespace algo
{
namespace detail
{

struct DefaultDistanceCalculator
{
    constexpr DefaultDistanceCalculator() = default;

    template<typename T>
    constexpr unsigned operator()(const detail::Node<T>&, const detail::Node<T>&) const
    {
        return 1;
    }
};

template<typename NodeType>
std::vector<std::reference_wrapper<typename NodeType::DataType>> reconstructPath(NodeType* node)
{
    if (node == nullptr)
        throw std::invalid_argument("reconstructPath: got nullptr");

    std::vector<std::reference_wrapper<typename NodeType::DataType>> ret;

    auto* currentNode = node;
    while (currentNode->parent != currentNode)
    {
        ret.push_back(std::ref(currentNode->data()));
        currentNode = currentNode->parent;
    }
    ret.push_back(std::ref(currentNode->data()));
    std::reverse(ret.begin(), ret.end());
    return ret;
}

// EndNodeMatchCallBack: bool endNodeCallback(GraphType::NodeType*)
// returns true when bfs should be stopped after it's called
template <typename GraphType,
          typename EndNodePred,
          typename FrontierPred,
          typename StopPred,
          typename EndNodeMatchCallBack>
inline void bfsImpl(
    GraphType& graph,
    const typename GraphType::KeyType& startPoint,
    EndNodePred endNodeMatch,
    FrontierPred frontierMatch,
    StopPred stopMatch,
    EndNodeMatchCallBack endNodeCallback)
{
    graph.reset();

    auto* startNode = graph.findNode(startPoint);
    startNode->state = detail::NodeState::Open;
    startNode->distance = 0;

    std::queue<decltype(startNode)> frontier;
    frontier.push(startNode);

    bool stop = false;
    unsigned stopDepth = std::numeric_limits<unsigned>::max();

    while (frontier.size() > 0)
    {
        auto* currentNode = frontier.front();
        frontier.pop();

        // On some occasions deeper nodes might be added to frontier before the stopMatch matches.
        // Disregard these nodes.
        if (currentNode->distance <= stopDepth && endNodeMatch(*currentNode))
        {
            if (endNodeCallback(currentNode))
                return;
        }

        if (!stop && stopMatch(*currentNode))
        {
            stop = true;
            stopDepth = currentNode->distance;
        }

        if (stop)
            continue;

        for (auto* neighbour : currentNode->neighbours())
        {
            if (!frontierMatch(*neighbour))
                continue;

            if (neighbour->state == detail::NodeState::NotVisited)
            {
                neighbour->state = detail::NodeState::Open;
                neighbour->distance = currentNode->distance + 1;
                neighbour->parent = currentNode;
                frontier.push(neighbour);
            }
        }
        currentNode->state = detail::NodeState::Closed;
    }
}

} // namespace detail
} // namespace algo

#endif // ALGO_DETAIL_GRAPHSEARCHHELPERS_HPP_
