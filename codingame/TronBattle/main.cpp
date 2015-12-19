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
#include <stdexcept>
#include <memory>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>

using Direction = std::string;
using DirectionVector = std::vector<Direction>;

static const Direction DU = "UP";
static const Direction DD = "DOWN";
static const Direction DL = "LEFT";
static const Direction DR = "RIGHT";
static const DirectionVector DIRECTIONS{DU, DD, DL, DR};

struct Coordinates
{
    int x = 0;
    int y = 0;

    Coordinates() : Coordinates(-1, -1)
    {
    }

    Coordinates(int x, int y) : x(x), y(y)
    {
    }

    bool areCorrect() const
    {
        return x >= 0 && y >= 0;
    }

    double distance(const Coordinates& coord) const
    {
        unsigned xSquare = std::pow(x - coord.x, 2);
        unsigned ySquare = std::pow(y - coord.y, 2);

        return std::sqrt(xSquare + ySquare);
    }

    std::string directionTo(const Coordinates& lhs) const
    {
        if (*this == lhs)
            throw std::invalid_argument("directionTo received the same coordinates as in this");

        if (areCorrect() && lhs.areCorrect())
        {
            Direction xDir;
            Direction yDir;

            if (lhs.x < x)
                xDir = DL;
            else if (lhs.x > x)
                xDir = DR;

            if (lhs.y < y)
                yDir = DU;
            else if (lhs.y > y)
                yDir = DD;

            std::cerr << "directionTo" << *this << "->" << lhs << "xDir=" << xDir << ", yDir=" << yDir << std::endl;

            if (xDir.empty())
                return yDir;
            else if (yDir.empty())
                return xDir;

            if (std::abs(lhs.x - x) < std::abs(lhs.y - y))
                return xDir;
            return yDir;
        }

        throw std::runtime_error("incorrect coordinates!");
    }

    friend bool operator==(const Coordinates& lhs, const Coordinates& rhs)
    {
        return (lhs.x == rhs.x && lhs.y == rhs.y);
    }

    friend bool operator!=(const Coordinates& lhs, const Coordinates& rhs)
    {
        return !(lhs == rhs);
    }

    friend std::ostream& operator<<(std::ostream& stream, const Coordinates& c)
    {
        return stream << "(" << c.x << "," << c.y << ")";
    }
};

class Field
{
    using Grid = std::vector<std::vector<bool>>;

public:
    Field() : grid_(30, std::vector<bool>(20, false))
    {
    }

    void take(const Coordinates& coord)
    {
        grid_[coord.x][coord.y] = true;
    }

    void clear(const Coordinates& coord)
    {
        grid_[coord.x][coord.y] = false;
    }

    bool at(const Coordinates& coord) const
    {
        if (coord.x < grid_.size() && coord.y < grid_[0].size())
            return grid_[coord.x][coord.y];
        return true;  // imaginary border around a grid
    }

private:
    Grid grid_;
};

class DirectionFinder
{
public:
    DirectionFinder(const DirectionFinder&) = delete;
    DirectionFinder operator=(const DirectionFinder) = delete;

    explicit DirectionFinder(const Field& field) : field_(field)
    {
    }

    std::vector<std::string> freeMoves(const Coordinates& coord) const
    {
        std::vector<std::string> ret;

        if (!field_.at(Coordinates(coord.x, coord.y - 1)))
            ret.push_back(DU);
        if (!field_.at(Coordinates(coord.x, coord.y + 1)))
            ret.push_back(DD);
        if (!field_.at(Coordinates(coord.x - 1, coord.y)))
            ret.push_back(DL);
        if (!field_.at(Coordinates(coord.x + 1, coord.y)))
            ret.push_back(DR);

        return ret;
    }

    std::string findDirection(const Coordinates& from, const Coordinates& to,
        DirectionVector forbidden = DirectionVector()) const
    {
        return findDirection(from, to, "", forbidden);
    }

    std::string findDirection(const Coordinates& from, const Coordinates& to,
        const Direction& preference, DirectionVector forbidden = DirectionVector()) const
    {
        std::vector<std::string> availableDirs = freeMoves(from);
        std::string neededDir = from.directionTo(to);

        std::sort(availableDirs.begin(), availableDirs.end());
        std::sort(forbidden.begin(), forbidden.end());
        std::vector<std::string> allowedDirs;
        std::set_difference(
            availableDirs.begin(), availableDirs.end(),
            forbidden.begin(), forbidden.end(), std::back_inserter(allowedDirs));

        if (std::find(allowedDirs.begin(), allowedDirs.end(), neededDir) == allowedDirs.end())
        {
            // There is no straight path to the given point
            // When using find_first_of, the order of elements is more than significant.
            std::vector<std::string> dirs;
            if (neededDir == DR)
                dirs = {DU, DD, DL};
            else if (neededDir == DL)
                dirs = {DU, DL, DR};
            else if (neededDir == DU)
                dirs = {DL, DR, DD};
            else if (neededDir == DD)
                dirs = {DL, DR, DU};

            if (!preference.empty()) 
            {
                auto found = std::find(dirs.begin(), dirs.end(), preference);
                if (found != dirs.end())
                {
                    if (dirs[0] == DU && (preference == DL || preference == DR) ||
                        (dirs[0] == DL && (preference == DU || preference == DD)))
                    {
                        std::cerr << "Moving preference to the begin: " << preference << std::endl;
                        dirs.erase(found);
                        dirs.insert(dirs.begin(), preference);
                    }
                }
            }

            auto foundDir = std::find_first_of(
                allowedDirs.begin(), allowedDirs.end(), dirs.begin(), dirs.end());
            if (foundDir != allowedDirs.end())
                neededDir = *foundDir;
            else
                neededDir = "";
        }

        return neededDir;
    }

    Coordinates directionToCoord(const Coordinates& from, const std::string& dir) const
    {
        if (dir == DR)
            return Coordinates(from.x + 1, from.y);
        else if(dir == DL)
            return Coordinates(from.x - 1, from.y);
        else if(dir == DU)
            return Coordinates(from.x, from.y - 1);
        else if(dir == DD)
            return Coordinates(from.x, from.y + 1);

        throw std::invalid_argument("incorrect direction: " + dir);
    }

private:
    const Field& field_;
};

class Player
{
public:
    using Ribbon = std::vector<Coordinates>;

    explicit Player(unsigned id, const Coordinates& initialPos = Coordinates()) : id_(id)
    {
        moveTo(initialPos);
    }

    unsigned id() const
    {
        return id_;
    }

    void moveTo(const Coordinates& coord)
    {
        if (!coord.areCorrect())
            return;
        if (!ribbon_.empty() && ribbon_.back() == coord)
            return;
        ribbon_.push_back(coord);
    }

    const Ribbon& ribbon() const
    {
        return ribbon_;
    }

    void clearRibbon()
    {
        ribbon_.clear();
    }

    Coordinates head() const
    {
        if (!ribbon_.empty())
            return ribbon_.back();
        return Coordinates(-1, -1);
    }

private:
    unsigned id_;
    Ribbon ribbon_;
};

class Strategy
{
public:
    virtual std::string evaluate(const Field&, const Player&) const = 0;

protected:
    DirectionVector unsafeDirections(Coordinates& c, const Field& field) const
    {
        std::vector<std::string> ret;
        for (auto dir : DIRECTIONS)
        {
            if (!directionSafe(c, dir, field))
                ret.push_back(dir);
        }
        return ret;
    }

    bool directionSafe(const Coordinates& c, std::string dir, const Field& field) const
    {
        Coordinates moveTo;
        Coordinates moveToNext;
        Coordinates moveToPrev;

        if (dir == DU)
        {
            moveTo = Coordinates(c.x, c.y - 1);
            moveToNext = Coordinates(c.x - 1, c.y - 1);
            moveToPrev = Coordinates(c.x + 1, c.y - 1);
        }
        else if (dir == DD)
        {
            moveTo = Coordinates(c.x, c.y + 1);
            moveToNext = Coordinates(c.x - 1, c.y + 1);
            moveToPrev = Coordinates(c.x + 1, c.y + 1);
        }
        else if (dir == DL)
        {
            moveTo = Coordinates(c.x - 1, c.y);
            moveToNext = Coordinates(c.x - 1, c.y - 1);
            moveToPrev = Coordinates(c.x - 1, c.y + 1);
        }
        else if (dir == DR)
        {
            moveTo = Coordinates(c.x + 1, c.y);
            moveToNext = Coordinates(c.x + 1, c.y - 1);
            moveToPrev = Coordinates(c.x + 1, c.y + 1);
        }
        else
        {
            throw std::invalid_argument("incorrect direction: " + dir);
        }

        return (!field.at(moveTo)) && (!field.at(moveToNext) || !field.at(moveToPrev));
    }
};

class FreeRoam : public Strategy
{
public:
    std::string evaluate(const Field& field, const Player& myPlayer) const
    {
        std::cerr << "FreeRoam::evaluate" << std::endl;

        Coordinates coord = myPlayer.ribbon().back();

        DirectionFinder finder(field);
        for (auto direction : DIRECTIONS)
        {
            Coordinates dirCoord = finder.directionToCoord(coord, direction);
            if (!(field.at(dirCoord)))
                return direction;
        }

        return "I want to break free!";
    }
};

class ToCoordinates : public Strategy
{
public:
    explicit ToCoordinates(const Coordinates& coord) : requiredCoord_(coord)
    {
    }

    std::string evaluate(const Field& field, const Player& myPlayer) const
    {
        std::cerr << "ToCoordinates::evaluate: " << requiredCoord_ << std::endl;

        Coordinates coord = myPlayer.ribbon().back();

        std::cerr << "coords: " << coord << std::endl;

        DirectionVector unsafe = unsafeDirections(coord, field);

        std::cerr << "unsafe:";
        for (auto us : unsafe)
            std::cerr << " " << us;
        std::cerr << std::endl;

        DirectionFinder finder(field);
        std::string direction = finder.findDirection(coord, requiredCoord_, unsafe);
        if (direction.empty())
            direction = finder.findDirection(coord, requiredCoord_);

        return direction;
    }

private:
    Coordinates requiredCoord_;
};

/*
class HugTheWall : public Strategy
{
public:
    HugTheWall(unsigned myId) : myId_(myId) {}

    std::string evaluate(const Field& field, const std::vector<Player>& players) const
    {
        const Player& myPlayer = players[myId_];

        Coordinates coord = myPlayer.ribbon().back();
        auto ribbonSize = myPlayer.ribbon().size();

        Coordinates previousCoord;
        if (ribbonSize >= 2)
        {
            previousCoord = myPlayer.ribbon()[ribbonSize - 2];
        }

        FreeDirectionFinder finder(field);

        if (previousCoord.areCorrect() && isNextToWall(coord))
        {
            if (!inCorner(coord))
            {
                std::string dir = previousCoord.directionTo(coord);
                Coordinates next = finder.directionToCoord(coord, dir);
                return finder.findDirection(coord, next);
            }
        }

        std::vector<Coordinates> edges;
        if (directionSafe(coord, DU))
            edges.push_back(Coordinates(coord.x, 0));
        if (directionSafe(coord, DD))
            edges.push_back(Coordinates(coord.x, 19));
        if (directionSafe(coord, DL))
            edges.push_back(Coordinates(0, coord.y));
        if (directionSafe(coord, DR))
            edges.push_back(Coordinates(29, coord.y));

        if (edges.empty())
        {
            // there's nothing more to do but to enter the unknown tunnel...
            Coordinates up = finder.directionToCoord(coord, DU);
            Coordinates down = finder.directionToCoord(coord, DD);
            Coordinates left = finder.directionToCoord(coord, DL);
            Coordinates right = finder.directionToCoord(coord, DR);

            if (up.areCorrect() && !field.at(up))
                edges.push_back(up);
            if (down.areCorrect() && !field.at(down))
                edges.push_back(down);
            if (left.areCorrect() && !field.at(left))
                edges.push_back(left);
            if (right.areCorrect() && !field.at(right))
                edges.push_back(right);
        }

        if (edges.empty())
        {
            // we lost so just continue the move just to return something
            return previousCoord.directionTo(coord);
        }

        auto closestEdge = std::min_element(edges.begin(), edges.end(),
            [&coord](const Coordinates& lhs, const Coordinates& rhs)
            { return coord.distance(lhs) < coord.distance(rhs); });

        return finder.findDirection(coord, *closestEdge);
    }

private:
    bool isNextToWall(const Coordinates& c) const
    {
        return (c.x == 0 || c.x == 29 || c.y == 0 || c.y == 19);
    }

    bool inCorner(const Coordinates& c) const
    {
        return (c == Coordinates(0, 0) || c == Coordinates(0, 19)
                || c == Coordinates(29, 0) || c == Coordinates(29, 19));
    }

private:
    unsigned myId_;
};
*/

class Game
{
public:
    Game()
    {
        req_ = {
            Coordinates(14, 9),
            Coordinates(14, 0),
            Coordinates(15, 19),
            Coordinates(0, 19),
            Coordinates(0, 0),
            Coordinates(0, 19),
            Coordinates(29, 19)
        };
    }

    void run()
    {
        while (true)
        {
            parseInput();
            play();
        }
    }

private:
    void parseInput()
    {
        std::cin >> totalPlayers_ >> myPlayerId_;

        for (unsigned id = 0; id < totalPlayers_; ++id)
        {
            int tailX, tailY, headX, headY;
            std::cin >> tailX >> tailY >> headX >> headY;
            Coordinates movedTo(headX, headY);

            allocatePlayer(id);
            Player& player = players_[id];

            if (movedTo.areCorrect())
            {
                field_.take(movedTo);
                player.moveTo(movedTo);
            }
            else
            {
                const Player::Ribbon& ribbon = player.ribbon();
                for (const Coordinates& ribbonCoord : ribbon)
                {
                    field_.clear(ribbonCoord);
                }
                player.clearRibbon();
            }

            if (id == myPlayerId_ && !movedTo.areCorrect())
            {
                std::cerr << "I lost! :(" << std::endl;
                return;
            }
        }
    }

    void allocatePlayer(unsigned id)
    {
        if (players_.size() <= id)
        {
            for (unsigned newId = players_.size(); newId <= id; ++newId)
                players_.emplace_back(newId);
        }
    }

    void play()
    {
        const Player& myPlayer = players_[myPlayerId_];
        std::unique_ptr<Strategy> strategy = chooseStrategy();
        std::string direction = strategy->evaluate(field_, myPlayer);
        std::cout << direction << std::endl;
    }

    std::unique_ptr<Strategy> chooseStrategy() const
    {
        std::unique_ptr<Strategy> strategy(new FreeRoam());

        auto freeCoord = std::find_if(req_.begin(), req_.end(),
            [this](const Coordinates& c) { return !field_.at(c); });
        if (freeCoord != req_.end())
        {
            strategy.reset(new ToCoordinates(*freeCoord));
        }

        return strategy;
    }

private:
    unsigned myPlayerId_ = std::numeric_limits<unsigned>::max();
    unsigned totalPlayers_ = 0;

    Field field_;
    std::vector<Player> players_;

    std::vector<Coordinates> req_;
};

int main()
{
    Game g;
    g.run();
    return 0;
}
