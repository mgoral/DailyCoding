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
 */

#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <memory>
#include <limits>
#include <cmath>

struct Coordinates;
class Zone;
struct Drone;
struct Team;

static int TURNS_LEFT = 200;

using CoordinatesPtr= std::shared_ptr<Coordinates>;
using ZonePtr = std::shared_ptr<Zone>;
using DronePtr = std::shared_ptr<Drone>;
using TeamPtr = std::shared_ptr<Team>;

struct Coordinates
{
    Coordinates(int x, int y) : x(x), y(y)
    {
    }

    Coordinates() : Coordinates(0, 0)
    {
    }

    double distance(const Coordinates& p) const
    {
        unsigned xSquare = std::pow(x - p.x, 2);
        unsigned ySquare = std::pow(y - p.y, 2);

        return std::sqrt(xSquare + ySquare);
    }

    int x;
    int y;

    friend std::ostream& operator<<(std::ostream& stream, const Coordinates& p)
    {
        return stream << p.x << " " << p.y;
    }

    friend bool operator==(const Coordinates& rhs, const Coordinates& lhs)
    {
        return (rhs.x == lhs.x && rhs.y == lhs.y);
    }

    friend bool operator!=(const Coordinates& rhs, const Coordinates& lhs)
    {
        return !(rhs == lhs);
    }
};

struct Objective
{
    Objective(
        const Coordinates& coord,
        int score = 0,
        int neededDrones = std::numeric_limits<int>::max()) :
            coord(coord),
            commandedAtTurn(200 - TURNS_LEFT),
            possibleScore(score),
            neededDrones(neededDrones)
    {
    }

    bool operator<(const Objective& other)
    {
        return possibleScore < other.possibleScore;
    }

    Coordinates coord = Coordinates();
    int commandedAtTurn = 0;
    int maxTurn = 200;
    int possibleScore = 0;
    int neededDrones = 1;

};

class Drone
{
public:
    Drone(unsigned tid, unsigned id, const Coordinates& coord) : tid_(tid), id_(id), position_(coord)
    {
    }

    void shout()
    {
        if (!objectives_.empty())
        {
            auto obj = std::max_element(objectives_.begin(), objectives_.end());
            std::cerr << "obj: " << obj->possibleScore << std::endl;
            std::cout << obj->coord << std::endl;
        }
        else
        {
            std::cerr << id_ << ": no objective!!!" << std::endl;
            std::cout << position_ << std::endl;;
        }
    }

    bool hasObjective()
    {
        return !objectives_.empty();
    }

    void changePosition(const Coordinates& coord)
    {
        position_ = coord;
        unsigned currTurn = 200 - TURNS_LEFT;
        objectives_.erase(
            std::remove_if(objectives_.begin(), objectives_.end(),
                [&coord, currTurn](const Objective& obj)
                {return obj.coord == coord || obj.maxTurn <= currTurn;}),
            objectives_.end());
    }

    const Coordinates& position() const
    {
        return position_;
    }

    void addObjective(const Objective& obj)
    {
        auto sameObj = std::find_if(objectives_.begin(), objectives_.end(),
            [&obj](const Objective& lhs) { return obj.coord == lhs.coord; });
        if (sameObj != objectives_.end())
            *sameObj = obj;
        else
            objectives_.push_back(obj);
    }

    unsigned owner() const
    {
        return tid_;
    }

    unsigned id() const
    {
        return id_;
    }

    bool operator<(const Drone& other)
    {
        return id_ < other.id_;
    }

private:
    unsigned tid_;
    unsigned id_;
    Coordinates position_;
    std::vector<Objective> objectives_;
};

class Zone
{
public:
    static const unsigned RADIUS = 100;
    static const unsigned CLOSE_DISTANCE = 201;

    Zone(unsigned id, const Coordinates& c, unsigned teamNo) :
        zid_(id), center_(c), drones_(teamNo)
    {
    }

    void addDrone(const DronePtr& drone)
    {
        unsigned teamId = drone->owner();
        drones_.at(teamId).push_back(drone);
    }

    Coordinates center() const
    {
        return center_;
    }

    void changeOwner(int owner)
    {
        owner_ = owner;
    }

    int owner() const
    {
        return owner_;
    }

    bool isOwner(int playerId)
    {
        return owner_ == playerId;
    }

    std::vector<DronePtr> allDrones(unsigned playerId) const
    {
        return drones_[playerId];
    }

    unsigned maxEnemiesCountInRadius(unsigned playerId, unsigned radius = Zone::RADIUS) const
    {
        unsigned ret = 0;
        for (int i = 0; i < drones_.size(); ++i)
        {
            if (i != playerId)
                ret = std::max(dronesCountInRadius(i, radius), ret);
        }

        return ret;
    }

    unsigned dronesCountInRadius(unsigned playerId, unsigned radius = Zone::RADIUS) const
    {
        return std::count_if(drones_[playerId].begin(), drones_[playerId].end(),
            [this, radius](const DronePtr& lhs) { return distance(*lhs) <= radius; });
    }

    std::vector<DronePtr> dronesInsideRadius(unsigned playerId, unsigned radius = Zone::RADIUS) const
    {
        std::vector<DronePtr> ret;
        for (const auto& drone : drones_[playerId])
        {
            if (distance(*drone) <= radius)
                ret.push_back(drone);
        }
        return ret;
    }

    std::vector<DronePtr> dronesOutsideRadius(unsigned playerId, unsigned radius = Zone::RADIUS) const
    {
        std::vector<DronePtr> ret;
        for (const auto& drone : drones_[playerId])
        {
            if (distance(*drone) > radius)
                ret.push_back(drone);
        }
        return ret;
    }

    double distance(const Drone& drone) const
    {
        return center_.distance(drone.position());
    }

    unsigned id() const
    {
        return zid_;
    }

    // Zone should be notified each turn that drones have finished their movements.
    void allDronesUpdated()
    {
        for (auto& droneVec : drones_)
        {
            std::stable_sort(droneVec.begin(), droneVec.end(),
                [this](const DronePtr& lhs, const DronePtr& rhs)
                { return distance(*lhs) < distance(*rhs); });
        }
    }

    friend std::ostream& operator<<(std::ostream& stream, const Zone& z)
    {
        return stream << "zid[" << z.zid_ << "]:";
    }

private:
    Coordinates center_;
    int owner_ = -1;
    unsigned zid_ = 0; // debug purposes

    std::vector<std::vector<DronePtr>> drones_;
};

struct Team
{
    Team(unsigned tid, std::vector<DronePtr>::size_type numberOfDrones)
    {
        for (unsigned id = 0; id < numberOfDrones; ++id)
            drones.push_back(std::make_shared<Drone>(tid, id, Coordinates()));
    }

    void addObjective(const Objective& obj)
    {
        for (auto& drone : drones)
            drone->addObjective(obj);
    }

    std::vector<DronePtr> drones;
};

class Game
{
public:
    Game()
    {
        unsigned playerNo, zoneNo;
        std::cin >> playerNo >> playerId_ >> dronesPerTeam_ >> zoneNo;

        initTeams(playerNo);
        initZones(zoneNo, playerNo);
    }

    void run()
    {
        while(true)
        {
            --TURNS_LEFT;
            unsigned roundNo = 200 - TURNS_LEFT;

            std::cerr << "Starting round << " << roundNo << std::endl;;

            update();
            play(roundNo);
        }
    }

private:
    void initTeams(unsigned playerNo)
    {
        for (unsigned teamId = 0; teamId < playerNo; ++teamId)
            teams_.push_back(std::make_shared<Team>(teamId, dronesPerTeam_));
    }

    void initZones(unsigned zoneNo, unsigned playerNo)
    {
        for (unsigned zoneId = 0; zoneId < zoneNo; ++zoneId)
        {
            int zoneX, zoneY;
            std::cin >> zoneX >> zoneY;
            ZonePtr zone = std::make_shared<Zone>(zoneId, Coordinates{zoneX, zoneY}, playerNo);

            for (auto& team : teams_)
                for (auto drone : team->drones)
                    zone->addDrone(drone);

            zones_.push_back(zone);
        }
    }

    void update()
    {
        for (auto& zone : zones_)
        {
            int owner;
            std::cin >> owner;
            zone->changeOwner(owner);
        }

        for (unsigned teamNo = 0; teamNo < teams_.size(); ++teamNo)
        {
            for (auto& drone : teams_[teamNo]->drones)
            {
                Coordinates p;
                std::cin >> p.x >> p.y;
                drone->changePosition(Coordinates(p.x, p.y));
            }
        }

        for (auto& zone : zones_)
            zone->allDronesUpdated();

    }

    int enemiesToMyDronesDifference(const ZonePtr& zone, unsigned radius = Zone::RADIUS)
    {
        auto myDrones = zone->dronesCountInRadius(playerId_, radius);
        auto enemies = zone->maxEnemiesCountInRadius(playerId_, radius);
        return enemies - myDrones;
    }

    void play(unsigned roundNo)
    {
        // logic here
        auto myTeam = teams_[playerId_];

        for (auto& zone : zones_)
        {
            Objective obj(zone->center(), 0, std::numeric_limits<int>::max());

            int dronesToOwnZone = 1 + zone->maxEnemiesCountInRadius(playerId_);

            if (!zone->isOwner(playerId_) && dronesToOwnZone > dronesPerTeam_)
            {
                myTeam->addObjective(Objective(zone->center()));
                continue;
            }

            unsigned modifier = 10;

            if (zone->isOwner(playerId_))
            {
                // defence parameters
                std::cerr << *zone << " defending" << std::endl;

                unsigned bigRadius = 510;
                unsigned defenders = zone->dronesCountInRadius(playerId_);
                unsigned potentialDefenders = zone->dronesCountInRadius(playerId_, bigRadius);
                unsigned attackers = zone->maxEnemiesCountInRadius(playerId_, bigRadius);
                unsigned diff = enemiesToMyDronesDifference(zone, bigRadius);

                if (defenders < attackers  && potentialDefenders >= attackers)
                {
                    modifier = 10;
                    dronesToOwnZone += diff;
                }
                else if (defenders >= attackers && attackers > 0)
                {
                    modifier = 10;
                    dronesToOwnZone = attackers;
                }
                else
                {
                    // nothing nearby or not worth defending, may leave zone
                    modifier = 0;
                }
            }
            else if (zone->isOwner(-1))
            {
                // zone not taken
                std::cerr << *zone << " capturing" << std::endl;
                modifier = 8;
            }
            else
            {
                // attack
                std::cerr << *zone << " attacking" << std::endl;
                modifier = 6;
            }

            std::vector<DronePtr> dronesInside = zone->dronesInsideRadius(playerId_);
            std::vector<DronePtr> dronesOutside = zone->dronesOutsideRadius(playerId_);
            std::vector<DronePtr> allDrones = zone->allDrones(playerId_);

            int maximumIndex = static_cast<int>(dronesOutside.size()) - 1;
            DronePtr furthestNeededDrone = allDrones[std::max(0, dronesToOwnZone - 1)];
            //    dronesOutside[std::min(maximumIndex, std::max(0, dronesToOwnZone - 1))];
            double maxDistance = zone->distance(*furthestNeededDrone);
            unsigned maxTurnsToReach = static_cast<int>((maxDistance - 1.0) / 100.0);

            int calculatedScore = (TURNS_LEFT - dronesToOwnZone * maxTurnsToReach) * modifier;

            obj.possibleScore = std::max(0, calculatedScore);
            obj.maxTurn = roundNo + maxTurnsToReach;
            obj.neededDrones = dronesToOwnZone;

            std::cerr << *zone << " score=" << obj.possibleScore << ", turns=" << obj.maxTurn << std::endl;

            for (int i = 0; i < dronesToOwnZone && i < allDrones.size(); ++i)
                allDrones[i]->addObjective(obj);

            if (!zone->isOwner(playerId_))
            {
                Objective roamingObj = obj;
                roamingObj.possibleScore /= 2;
                roamingObj.maxTurn = roundNo + 1;
                for (auto& drone : allDrones)
                {
                    if (!drone->hasObjective())
                        drone->addObjective(roamingObj);
                }
            }
        }

        for (auto& drone : myTeam->drones)
        {
            drone->shout();
        }
    }

private:
    unsigned playerId_;
    std::vector<DronePtr>::size_type dronesPerTeam_;
    std::vector<ZonePtr> zones_;
    std::vector<TeamPtr> teams_;
};

int main()
{
    Game g;
    g.run();

    return 0;
}
