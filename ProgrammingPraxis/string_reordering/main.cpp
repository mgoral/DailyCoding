/**
 * String Re-ordering. March 27, 2015
 * http://programmingpraxis.com/2015/03/27/string-re-ordering/
 *
 * You are given a string O that specifies the desired ordering of letters in a target string T. 
 * For example, given string O = “eloh” the target string T = “hello” would be re-ordered “elloh” 
 * and the target string T = “help” would be re-ordered “pelh” (letters not in the order string are
 * moved to the beginning of the output in some unspecified order).
 *
 * Your task is to write a program that produces the requested string re-ordering.
 */

#include <iostream>
#include <string>
#include <unordered_map>
#include <algorithm>

using OrderT = std::unordered_map<char, int>;

OrderT makeCmpMap(const std::string& orderStr)
{
    OrderT order;

    int i = 0;
    for (const auto& ch : orderStr)
        order[ch] = ++i;

    return order;
}

std::string reorder(std::string target, const OrderT& order)
{
    std::sort(target.begin(), target.end(), 
        [&order](const char& lhs, const char& rhs) 
        {
            auto lhsIt = order.find(lhs);
            auto rhsIt = order.find(rhs);
            int lhsOrder = lhsIt == order.end() ? 0 : lhsIt->second;
            int rhsOrder = rhsIt == order.end() ? 0 : rhsIt->second;
            return lhsOrder < rhsOrder;
        });
    return target;
}

int main()
{
    std::string order, target;

    std::cout << "Order letters: ";
    std::cin >> order;
    std::cout << "Target letters: ";
    std::cin >> target;

    auto orderMap =  makeCmpMap(order);
    auto reordered = reorder(target, orderMap);

    std::cout << "Reordered target string: \"" + reordered + "\"" << std::endl;
}
