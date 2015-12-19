#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

bool cmp (int lhs, int rhs)
{
    auto alhs = std::abs(lhs);
    auto arhs = std::abs(rhs);

    if (alhs == arhs)
        return lhs >= 0;
    return alhs < arhs;
}

int main()
{
    std::vector<int> t;

	int n;
	std::cin >> n;
	for (int i = 0; i < n; ++i) {
	    int temp;
	    std::cin >> temp;
	    t.push_back(temp);
	}

	auto found = std::min_element(t.begin(), t.end(), cmp);

	if (found != t.end())
	    std::cout << *found << std::endl;
    else
    	std::cout << "0" << std::endl;

	return 0;
}
