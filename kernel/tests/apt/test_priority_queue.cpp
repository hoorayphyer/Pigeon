#include "all_in_one.hpp"
#include "apt/priority_queue.hpp"

using namespace apt;

SCENARIO("Test apt::priority_queue", "[apt]") {
  priority_queue<long double> pq;
  pq.push(0, 1.00);
  pq.push(1, 2.00);
  {
    auto [i , val] = pq.top();
    REQUIRE( i == 1 );
    REQUIRE( val == 2.00 );
  }
  pq.pop();
  {
    auto [i , val] = pq.top();
    REQUIRE( i == 0 );
    REQUIRE( val == 1.00 );
  }
  pq.push(1, 0.5);
  {
    auto [i , val] = pq.top();
    REQUIRE( i == 0 );
    REQUIRE( val == 1.00 );
  }
  pq.push(2, 2.5);
  {
    auto [i , val] = pq.top();
    REQUIRE( i == 2 );
    REQUIRE( val == 2.5 );
  }
}
