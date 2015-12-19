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

#include <iostream>
#include <array>
#include <vector>
#include <algorithm>

#include "algo/Graph.hpp"
#include "algo/GraphSearch.hpp"

const constexpr unsigned POD_COST = 20;

struct Assignment
{
    int priority = -1;
    std::vector<unsigned> path;
    unsigned neededPods = 0;

    bool operator==(const Assignment& other) const
    {
        if (!path.empty() && !other.path.empty())
            return path.back() == other.path.back();
        return true;
    }

    bool operator<(const Assignment& other) const { return priority < other.priority; }
    bool operator>(const Assignment& other) const { return priority > other.priority; }
};

class AssignmentsManager
{
public:
    bool hasAssignmentFor(unsigned zoneId) const
    {
        return moves_.find(zoneId) != moves_.end();
    }

    void addMoveCommand(const Assignment& assignment)
    {
        if (assignment.path.size() < 1)
            return;

        unsigned srcZone = assignment.path[0];
        auto& assignments = moves_[srcZone];
        auto it = std::find(assignments.begin(), assignments.end(), assignment);
        if (it != assignments.end())
        {
            Assignment concated = concatAssignments(*it, assignment);
            assignments.erase(it);
            assignments.push_back(concated);
        }
        else
        {
            assignments.push_back(assignment);
        }
    }

    template<typename GraphType>
    void performMoves(const GraphType& graph, unsigned playerNo)
    {
        if (moves_.empty())
        {
            std::cerr << "No moves planned" << std::endl;
            std::cout << "WAIT" << std::endl;
            return;
        }

        bool commandGiven = false;

        for (auto& val : moves_)
        {
            auto zoneId = val.first;
            auto& assignments = val.second;
            std::sort(assignments.begin(), assignments.end(), std::greater<Assignment>());
            auto playerPods = graph.data(zoneId).pods[playerNo];

            auto assignmentId = 0;

            while (playerPods > 0 && assignments.size() > assignmentId)
            {
                auto& a = assignments[assignmentId];
                unsigned usedPods = std::min(playerPods, a.neededPods);
                // if there's a path - it's a move command
                // otherwise it's a defence command
                if (a.path.size() > 1)
                {
                    commandGiven = true;
                    std::cout << usedPods << " " << a.path.at(0) << " " << a.path.at(1) << " ";
                }
                playerPods -= usedPods;
                ++assignmentId;
            }

        }

        if (commandGiven)
        {
            std::cout << std::endl;
        }
        else
        {
            std::cerr << "Only defence commands this turn." << std::endl;
            std::cout << "WAIT" << std::endl;
        }

        moves_.clear();
    }

    void addBuyCommand(const Assignment& assignment)
    {
        if (assignment.path.empty())
            return;

        unsigned srcZone = assignment.path[0];
        auto it = std::find_if(buys_.begin(), buys_.end(),
            [srcZone](const Assignment& a){ return a.path[0] == srcZone; });
        if (it != buys_.end() && it->neededPods >= assignment.neededPods)
            return;
        buys_.push_back(assignment);
    }

    void performBuys(unsigned availablePlatinum)
    {
        if (buys_.empty())
        {
            std::cerr << "No buy commands." << std::endl;
            std::cout << "WAIT" << std::endl;;
            return;
        }

        std::sort(buys_.begin(), buys_.end(), std::greater<Assignment>());

        const unsigned podCost = 20;
        auto buyId = 0;
        bool boughtSth = false;

        while (availablePlatinum >= podCost && buys_.size() > buyId)
        {
            auto& buy = buys_[buyId];
            unsigned neededPlat = buy.neededPods * podCost;
            if (neededPlat <= availablePlatinum)
            {
                std::cout << buy.neededPods << " " << buy.path[0] << " ";
                availablePlatinum -= (buy.neededPods * podCost);
                boughtSth = true;
            }
            ++buyId;
        }

        if (boughtSth) // any buy command was given
        {
            std::cout << std::endl;
        }
        else
        {
            std::cerr << "Not enough platinum to perform any buy." << std::endl;
            std::cout << "WAIT" << std::endl;;
        }

        buys_.clear();
    }

private:
    Assignment concatAssignments(const Assignment& lhs, const Assignment& rhs)
    {
        Assignment ret;
        ret.neededPods = lhs.neededPods + rhs.neededPods;
        if (lhs < rhs)
        {
            ret.priority = rhs.priority;
            ret.path = rhs.path;
        }
        else
        {
            ret.priority = lhs.priority;
            ret.path = lhs.path;
        }
        return ret;
    }

private:
    std::unordered_map<int, std::vector<Assignment>> moves_;
    std::vector<Assignment> buys_;
};

struct Zone
{
    unsigned id;
    int owner = -1;
    int platinum = 0;
    int neighbourhoodValue = 0;
    std::array<unsigned, 4> pods = {{0, 0, 0, 0}}; // single brace works in clang, not in g++...
    std::vector<Assignment> assignments;

    unsigned maxEnemyCount(unsigned playerNo) const
    {
        auto podsCopy = pods;
        podsCopy[playerNo] = 0;
        return *std::max_element(podsCopy.begin(), podsCopy.end());
    }

    bool hasPlayerPods(unsigned playerNo) const
    {
        return pods[playerNo] > 0;
    }
};

struct Player
{
    unsigned id = 0;
    unsigned platinum = 200;
};

class AssignmentGenerator
{
    using GraphType = algo::Graph<unsigned, Zone>;
    using NodeType = GraphType::NodeType;

public:
    AssignmentGenerator(unsigned playerNo, unsigned totalPlayers, GraphType& graph) :
        playerNo_(playerNo), totalPlayers_(totalPlayers), graph_(graph)
    {
        for (auto& kv : graph_)
            addPlatinumNode(&kv.second);
    }

    void generateAssignments(unsigned availablePlatinum)
    {
        auto zonesWithMyPods = getZonesWithMyPods();

        conquerZones();
        defendZones();
        freeRoamAll(zonesWithMyPods);
        buyPods();

        assMan.performMoves(graph_, playerNo_);
        assMan.performBuys(availablePlatinum);
    }

private:
    void addPlatinumNode(NodeType* node)
    {
        // fill neighbourhoodValue for all nodes
        auto neighbours = node->neighbours();
        auto neighboursCount = neighbours.size();
        unsigned totalPlatinumInNeighbourhood = std::accumulate(
            neighbours.begin(), neighbours.end(), 0,
            [neighboursCount](unsigned acc, const NodeType* n){ return acc + n->data().platinum; } );

        node->data().neighbourhoodValue = node->data().platinum;// +  totalPlatinumInNeighbourhood / neighboursCount;

        if (node->data().platinum > 0)
        {
            platinumNodes_.push_back(node);
            ++platinumNodesCount_[node->data().platinum - 1];
        }
    }

    bool tryNinjaAttack(NodeType* platinumNode)
    {
        unsigned maxEnemiesOnPlat = platinumNode->data().maxEnemyCount(playerNo_);
        unsigned platinumAvailable = platinumNode->data().neighbourhoodValue;
        unsigned pno = playerNo_;

        // Ninja attack
        if (maxEnemiesOnPlat <= 3)
        {
            auto neighbours = platinumNode->neighbours();
            auto neededPodsToAttack = maxEnemiesOnPlat + 1;

            std::sort(neighbours.begin(), neighbours.end(),
                [pno](NodeType* lhs, NodeType* rhs)
                { return lhs->data().pods[pno] > rhs->data().pods[pno]; });

            std::vector<Assignment> tempAsses;

            for (auto* neighbour : neighbours)
            {
                Assignment a;
                a.priority = platinumAvailable;
                a.neededPods = std::min(neededPodsToAttack, neighbour->data().pods[pno]);
                a.path = {neighbour->data().id, platinumNode->data().id};

                neededPodsToAttack -= a.neededPods;

                if (a.neededPods > 0)
                    tempAsses.push_back(a);

                if (neededPodsToAttack == 0)
                    break;
            }

            if (neededPodsToAttack == 0 && tempAsses.size() > 1)
            {
                for (auto& a : tempAsses)
                {
                    assMan.addMoveCommand(a);
                    std::cerr << "Becoming a ninja: " << a.path[0]  << "->" << a.path[1] << std::endl;
                }
                return true;
            }
        }
        return false;
    }

    void conquerZones()
    {
        for (auto* platinumNode : platinumNodes_)
        {
            if (platinumNode->data().owner == playerNo_)
                continue;

            if (tryNinjaAttack(platinumNode))
                continue;

            unsigned maxEnemiesOnPlat = platinumNode->data().maxEnemyCount(playerNo_);
            unsigned platinumAvailable = platinumNode->data().neighbourhoodValue;
            unsigned pno = playerNo_;

            auto pathToClosestZoneWithEnoughPods  = algo::findShortestPathBfs(
                graph_, platinumNode->data().id,
                [pno, maxEnemiesOnPlat](const NodeType& n)
                    { return n.data().pods[pno] > maxEnemiesOnPlat;  },
                [pno, maxEnemiesOnPlat](const NodeType& frontier)
                    {
                        auto mec = frontier.data().maxEnemyCount(pno);
                        return (mec == 0 || mec < maxEnemiesOnPlat);
                    },
                algo::GraphPred<GraphType>::nodeFalse);

            if (!pathToClosestZoneWithEnoughPods.empty())
            {
                unsigned pathSize = pathToClosestZoneWithEnoughPods.size();
                unsigned distance = pathSize - std::min(2u, pathSize);

                Assignment a;
                a.priority = platinumAvailable - std::min(platinumAvailable, distance);
                a.neededPods = maxEnemiesOnPlat + 1;  // let's hope it won't overflow ;)

                // we searched starting from platinum node, but assignment is set starting from node
                // with pod (from found path end)
                std::vector<unsigned> path;
                std::transform(
                    pathToClosestZoneWithEnoughPods.rbegin(),
                    pathToClosestZoneWithEnoughPods.rend(),
                    std::back_inserter(a.path), [](Zone& z) { return z.id; });
                assMan.addMoveCommand(a);
            }
        }
    }

    void defendZones()
    {
        for (auto* platinumNode : platinumNodes_)
        {
            if (platinumNode->data().owner != playerNo_)
                continue;

            std::array<unsigned, 4> enemyCount = {{0, 0, 0, 0}};

            for (auto* neighbour : platinumNode->neighbours())
            {
                enemyCount[0] += neighbour->data().pods[0];
                enemyCount[1] += neighbour->data().pods[1];
                enemyCount[2] += neighbour->data().pods[2];
                enemyCount[3] += neighbour->data().pods[3];
            }

            enemyCount[playerNo_] = 0;

            unsigned maxEnemies = *std::max_element(enemyCount.begin(), enemyCount.end());

            if (maxEnemies > 0)
            {
                bool enoughExistingPods = false;
                Assignment a;
                a.priority = platinumNode->data().neighbourhoodValue + 1; // defensive assignments have higher priority
                a.neededPods = maxEnemies;

                for (auto* neighbour : platinumNode->neighbours())
                {
                    if (neighbour->data().pods[playerNo_] >= maxEnemies)
                    {
                        a.path = {neighbour->data().id, platinumNode->data().id};
                        assMan.addMoveCommand(a);
                        enoughExistingPods = true;
                        break;
                    }
                }

                if (!enoughExistingPods)
                {
                    a.path.push_back(platinumNode->data().id);
                    assMan.addBuyCommand(a);
                }

            }
        }
    }

    void freeRoamAll(const std::vector<std::reference_wrapper<Zone>>& zonesWithMyPods)
    {

        for (Zone& zone : zonesWithMyPods)
        {
            std::vector<unsigned> takenDestinations;
            std::vector<unsigned> firstFound;
            for (unsigned i = 0; i < zone.pods[playerNo_]; ++i)
            {
                auto closestUntakenPath = algo::findShortestPathBfs(graph_, zone.id,
                    [&takenDestinations, this](const NodeType& n)
                        { return n.data().owner != playerNo_ && !zoneFound(n.data().id, takenDestinations); },
                    algo::GraphPred<GraphType>::nodeTrue,
                    algo::GraphPred<GraphType>::nodeFalse);

                Assignment a;
                a.priority = 0;

                if (!closestUntakenPath.empty())
                {
                    std::vector<unsigned> path;
                    std::transform(closestUntakenPath.begin(), closestUntakenPath.end(),
                        std::back_inserter(a.path), [](Zone& z) { return z.id; });
                    if (firstFound.empty())
                        firstFound = a.path;

                    a.neededPods = 1;

                    takenDestinations.push_back(a.path.back());
                }
                else if (!firstFound.empty())
                {
                    a.path = firstFound;
                    a.neededPods = zone.pods[playerNo_]; // assignment manager will correctly move remaining pods
                    continue;
                }

                assMan.addMoveCommand(a);
            }
        }
    }

    std::vector<std::reference_wrapper<Zone>> getZonesWithMyPods()
    {
        std::vector<std::reference_wrapper<Zone>> ret;
        for (auto& kv : graph_)
        {
            Zone& z = kv.second.data();
            if (z.hasPlayerPods(playerNo_))
                ret.push_back(z);
        }
        return ret;
    }

    void buyPods()
    {
        for (auto* pn : platinumNodes_)
        {
            if (pn->data().owner != playerNo_)
            {
                unsigned neighbourhoodValue = pn->data().neighbourhoodValue;
                unsigned platinum = pn->data().platinum;

                // usually optimal placement during the first turn when there are more than 2
                // players ends with disaster...
                if (firstBuyEver_ && totalPlayers_ > 2 && platinum >= 5)
                    continue;

                auto nearbyPods = findMatchingNodesBfs(graph_, pn->data().id,
                    [this](const NodeType& n) { return n.data().pods[playerNo_] > 0; },
                    algo::GraphPred<GraphType>::nodeTrue,
                    [](const NodeType& n) { return n.distance > 1; }); // FIXME: empirical test what radius of check is the best

                decltype(nearbyPods) path;

                if (!nearbyPods.empty())
                {
                    auto it = std::max_element(nearbyPods.begin(), nearbyPods.end(),
                        [this](const Zone& lhs, const Zone& rhs)
                            { return lhs.pods[playerNo_] < rhs.pods[playerNo_]; });
                    path = findShortestPathBfs(graph_, pn->data().id,
                        [&it](const NodeType& n) { return n.data().id == it->get().id; },
                        [this](const NodeType& n) { return n.data().maxEnemyCount(playerNo_) == 0; },
                        algo::GraphPred<GraphType>::nodeFalse);
                }

                if (path.empty())
                {
                    path = findShortestPathBfs(graph_, pn->data().id,
                        [this](const NodeType& n)
                            { auto owner = n.data().owner; return (owner == -1 || owner == playerNo_); },
                        algo::GraphPred<GraphType>::nodeTrue,
                        algo::GraphPred<GraphType>::nodeFalse);
                }

                if (path.empty())
                    continue;

                unsigned podsToBuy = pn->data().maxEnemyCount(playerNo_) + 1;
                unsigned alreadyAvailablePods = path.back().get().pods[playerNo_];
                if (alreadyAvailablePods > podsToBuy)
                    continue;
                podsToBuy -= alreadyAvailablePods;

                Assignment a;
                a.priority = static_cast<int>(neighbourhoodValue);
                a.path.push_back(path.back().get().id);
                a.neededPods = podsToBuy;

                if (path.size() > 1)
                    a.priority = 1;

                if (firstBuyEver_ && neighbourhoodValue > 5 && platinum > 3)
                    a.neededPods += 1;

                assMan.addBuyCommand(a);
            }
        }

        for (auto kv : graph_)
        {
            const auto& node = kv.second;
            if (node.data().owner == -1)
            {
                Assignment a;
                a.priority = 0;
                a.neededPods = 1;
                a.path.push_back(node.data().id);
                assMan.addBuyCommand(a);
            }
        }

        firstBuyEver_ = false;
    }

    bool zoneFound(unsigned zoneId, const std::vector<unsigned>& foundZones) const
    {
        return std::find(foundZones.begin(), foundZones.end(), zoneId) != foundZones.end();
    }

private:
    unsigned playerNo_;
    unsigned totalPlayers_;
    GraphType& graph_;
    std::vector<NodeType*> platinumNodes_;
    AssignmentsManager assMan;
    bool firstBuyEver_ = true;
    std::array<int, 6> platinumNodesCount_ = {{0,0,0,0,0,0}};
};

class Game
{
    using GraphType = algo::Graph<unsigned, Zone>;

public:
    void run()
    {
        init();
        AssignmentGenerator assGen(me_.id, playersNo_,  graph_);

        while (true)
        {
            parseRoundInput();
            assGen.generateAssignments(me_.platinum);
        }
    }

private:
    void init()
    {
        unsigned linksNo;

        std::cin >> playersNo_ >> me_.id >> zoneNo_ >> linksNo;

        for (unsigned i = 0; i < zoneNo_; ++i)
        {
            Zone z;
            std::cin >> z.id >> z.platinum;
            auto* node = graph_.addNode(z.id, z);
        }

        for (unsigned i = 0; i < linksNo; ++i)
        {
            unsigned zid1, zid2;
            std::cin >> zid1 >> zid2;
            graph_.addConnection(zid1, zid2);
        }
    }

    void parseRoundInput()
    {
        std::cin >> me_.platinum;
        for (unsigned i = 0; i < zoneNo_; ++i)
        {
            unsigned zid;
            std::cin >> zid;
            auto& z = graph_.data(zid);
            std::cin >> z.owner >> z.pods[0] >> z.pods[1] >> z.pods[2] >> z.pods[3];
        }
    }

private:
    unsigned playersNo_;
    unsigned zoneNo_;
    Player me_;

    GraphType graph_;
};

int main()
{
    try
    {
        Game g;
        g.run();
    }
    catch (const std::exception& e)
    {
        std::cerr << "Uncaught exception. Program will now terminate." << std::endl;
        std::cerr << "what(): " << e.what() << std::endl;
        return 1;
    }
}
