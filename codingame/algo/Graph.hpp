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

#ifndef ALGO_GRAPH_HPP_
#define ALGO_GRAPH_HPP_

#include <unordered_map>
#include <stdexcept>
#include "detail/Node.hpp"

namespace algo
{

/**
 * @param Key: meets requirements of unordered_map::key_type (which means that it must be hashable
 *        with std::hash<Key> and hashes shouldn't collide for optimal complexity). ints are usually
 *        a good choice ;)
 * @param Value: stored data type which meets requirements of unordered_map::mapped_type. 
 *        Must be copy-construible.
 */
template <typename Key, typename Value>
class Graph
{
public:
    using KeyType = Key;
    using ValueType = Value;
    using NodeType = detail::Node<ValueType>;
    using NodesContainer = std::unordered_map<KeyType, detail::Node<ValueType>>;

    NodeType* addNode(const KeyType& key, const ValueType& data)
    {
        auto insertResult = nodes_.insert(std::make_pair(key, NodeType(data)));
        return &(insertResult.first->second);
    }

    NodeType* findNode(const KeyType& key)
    {
        auto it = nodes_.find(key);
        if (it == nodes_.end())
            throw std::logic_error("Node with a given key not found in graph!");
        return &(it->second);
    }

    const NodeType* findNode(const KeyType& key) const
    {
        auto it = nodes_.find(key);
        if (it == nodes_.end())
            throw std::logic_error("Node with a given key not found in graph!");
        return &(it->second);
    }

    // Adds two way connection
    void addConnection(const KeyType& lhs, const KeyType& rhs)
    {
        auto lhsIt = nodes_.find(lhs);
        auto rhsIt = nodes_.find(rhs);

        if (lhsIt == nodes_.end() || rhsIt == nodes_.end())
            return;

        NodeType& lhsNode = lhsIt->second;
        NodeType& rhsNode = rhsIt->second;

        lhsNode.addNeighbour(&rhsNode);
        rhsNode.addNeighbour(&lhsNode);
    }

    void addConnection(
        const KeyType& lhsKey, const ValueType& lhsData,
        const KeyType& rhsKey, const ValueType& rhsData)
    {
        NodeType* lhsNode = addNode(lhsKey, lhsData);
        NodeType* rhsNode = addNode(rhsKey, rhsData);

        lhsNode->addNeighbour(rhsNode);
        rhsNode->addNeighbour(lhsNode);
    }

    // Adds one way connection
    void addDirectedConnection(const KeyType& from, const KeyType& to)
    {
        auto fromIt = nodes_.find(from);
        auto toIt = nodes_.find(to);

        if (fromIt == nodes_.end() || toIt == nodes_.end())
            return;

        NodeType& fromNode = fromIt->second;
        NodeType& toNode = toIt->second;

        fromNode.addNeighbour(&toNode);
    }

    void addDirectedConnection(
        const KeyType& fromKey, const ValueType& fromData,
        const KeyType& toKey, const ValueType& toData)
    {
        NodeType* fromNode = addNode(fromKey, fromData);
        NodeType*& toNode = addNode(toKey, toData);

        fromNode->addNeighbour(toNode);
    }

    void reset()
    {
        for (auto& node : nodes_)
        {
            node.second.resetInfo();
        }
    }

    typename NodesContainer::iterator begin()
    {
        return nodes_.begin();
    }

    typename NodesContainer::iterator end()
    {
        return nodes_.end();
    }

    // Accessors
    ValueType& data(const KeyType& key) { return nodes_.at(key).data(); }
    const ValueType& data(const KeyType& key) const { return nodes_.at(key).data(); }

    NodeType& node(const KeyType& key) { return nodes_.at(key); }
    const NodeType& node(const KeyType& key) const { return nodes_.at(key); }

private:
    NodesContainer nodes_;
};

} // namespace algo

#endif // ALGO_GRAPH_HPP_
