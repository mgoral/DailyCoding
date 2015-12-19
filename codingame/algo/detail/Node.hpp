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

#ifndef ALGO_NODE_HPP_
#define ALGO_NODE_HPP_

#include <vector>
#include <functional>
#include <limits>
#include <algorithm>

namespace algo
{
namespace detail
{

enum class NodeState
{
    NotVisited = 1,
    Open = 2,
    Closed = 3
};

template <typename T>
class Node
{
public:
    using DataType = T;
    using NodeType = Node<DataType>;
    using NeighboursContainer = std::vector<NodeType*>;

    // constructors
    Node(const DataType& data, const NeighboursContainer&  neighbours) :
        data_(data), neighbours_(neighbours) { }
    explicit Node(const DataType& data) : Node(data, NeighboursContainer()) { }

    const DataType& data() const
    {
        return data_;
    }

    DataType& data()
    {
        return data_;
    }

    void addNeighbour(NodeType* n)
    {
        // Filter out incorrect neighbours.
        if (n == nullptr)
            return;

        if (this == n)
            return;

        if (std::find(neighbours_.begin(), neighbours_.end(), n) != neighbours_.end())
            return;

        neighbours_.push_back(n);
    }

    NeighboursContainer neighbours() const
    {
        return neighbours_;
    }

    void resetInfo()
    {
        parent = this;
        state = NodeState::NotVisited;
        distance = std::numeric_limits<unsigned>::max();
    }

    /*
    template<typename DistanceCalculator>
    auto distance(const NodeType& other, DistanceCalculator calc) -> decltype(calc(*this, other))
    {
        return calc(*this, other);
    }
    */

public:
    NodeType* parent;
    NodeState state;
    unsigned distance;

private:
    NeighboursContainer neighbours_;
    DataType data_;
};

} // namespace detail
} // namespace algo

#endif // ALGO_NODE_HPP_
