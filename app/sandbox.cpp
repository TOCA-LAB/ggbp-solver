#include <iostream>
#include <map>
#include <ranges>
#include <set>
#include <unordered_set>
#include <vector>

#include "utils/fix_c++_print.hpp"
#include "utils/contract.hpp"

#include "utils/combinations.hpp"

using namespace std;

int main() {

    vector<int> vec;

    for (int i = 0; i < 112; i++)
        vec.push_back(i);

    ggbp::utils::Combinations comb(vec, 2);

    for (auto c : comb) {
        std::cout << c << std::endl;
        int i;
        cin >> i;
    }


    return 0;
}

#if 0
int main() {

    std::vector<int> S = {1,2,3};
    std::set<std::string> names = {"Ana", "Bob"};
    std::unordered_set<int> T = {2, 6, 1, 3};
    std::map<int, std::string> dict;
    dict[1] = "bola";
    dict[2] = "ovo";
    dict[3] = "natal";

    std::cout << S << names <<  T << dict << std::endl;
}
#endif

#if 0
int main (int argc, char *argv[]) {
    

    std::vector<int> v = {1, 2, 3, 4, 5};
    PRECONDITION(all_in_interval(v, 2, 6));

    auto doubled_evens =
        v 
        | std::views::filter([](int x){ return x % 2 == 0; })
        | std::views::transform([](int x){ return x * 2; });

     std::vector<int> result(doubled_evens.begin(), doubled_evens.end());
     cout << result << endl;
        return 0;
}
#endif
