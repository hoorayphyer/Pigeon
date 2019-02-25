#include "apt/handle.hpp"
#include "catch2/catch.hpp"

bool Foo_destructor_called;
bool FooDeleter_called;

void reset_foo_status() {
  Foo_destructor_called = false;
  FooDeleter_called = false;
};

struct Foo {
  ~Foo() { Foo_destructor_called = true; }
};

void FooDeleter ( Foo* foo) {
  FooDeleter_called = true;
  delete foo;
}

SCENARIO("handle", "[apt]") {
  using FooHdl = apt::Handle<Foo, FooDeleter>;

  FooHdl h1;
  REQUIRE( h1.use_count() == 0 );
  auto* foo = new Foo;

  reset_foo_status();
  h1.reset(foo);
  REQUIRE( h1.use_count() == 1 );

  WHEN("being reset from holding nothing") {
    REQUIRE_FALSE(FooDeleter_called);
    REQUIRE_FALSE(Foo_destructor_called);
  }

  FooHdl h2 = h1;
  REQUIRE( h1.use_count() == 2 );
  reset_foo_status();
  h2.reset();
  REQUIRE( h1.use_count() == 1 );
  WHEN("a copy resets, nothing happens because h1 still references the resource") {
    REQUIRE_FALSE(FooDeleter_called);
    REQUIRE_FALSE(Foo_destructor_called);
  }

  reset_foo_status();
  h1.reset();
  REQUIRE( h1.use_count() == 0 );
  THEN("original copy frees") {
    REQUIRE(FooDeleter_called);
    REQUIRE(Foo_destructor_called);
  }
}

