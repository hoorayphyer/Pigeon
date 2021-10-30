#include "particle/map.hpp"
#include "testfw/testfw.hpp"

using namespace particle;

void check_consistency(const map<int>& m, const std::vector<species>& sps,
                       const std::vector<int>& values) {
  {  // loop order. NOTE this also checks that all keys are distinct
    int sp_prev = -1;
    bool is_ordered = true;
    for (auto sp : m) {
      is_ordered = sp_prev < static_cast<int>(sp);
      sp_prev = static_cast<int>(sp);
    }
    REQUIRE(is_ordered);
  }

  {  // has
    for (int s = 0; s < NUM_SPECIES; ++s) {
      auto sp = static_cast<species>(s);
      REQUIRE(m.has(sp) ==
              (std::find(sps.begin(), sps.end(), sp) != sps.end()));
    }
  }

  {  // values
    for (auto sp : m) {
      int i = 0;
      for (; i < sps.size(); ++i) {
        if (sp == sps[i]) break;
      }
      REQUIRE(m[sp] == values[i]);
    }
  }
}

SCENARIO("Test Map", "[particle]") {
  map<int> m;

  check_consistency(m, {}, {});

  m.insert(PO, 5);
  check_consistency(m, {PO}, {5});

  m.insert(EL, 3);
  check_consistency(m, {PO, EL}, {5, 3});

  m.insert(PO, 9);
  check_consistency(m, {PO, EL}, {9, 3});

  m.insert(PH, 8);
  check_consistency(m, {PH, PO, EL}, {8, 9, 3});

  m.erase(IO);
  check_consistency(m, {PH, PO, EL}, {8, 9, 3});

  m.erase(PO);
  check_consistency(m, {PH, EL}, {8, 3});

  m[EL] = 7;
  check_consistency(m, {PH, EL}, {8, 7});
}
