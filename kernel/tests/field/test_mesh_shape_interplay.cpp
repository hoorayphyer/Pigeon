#include "field/mesh_shape_interplay.cpp"
#include "kernel/grid.hpp"
#include "kernel/shapef.hpp"
#include "apt/print.hpp"
#include "apt/pair.hpp"
#include "catch2/catch.hpp"
#include <iostream>

using namespace field;
using apt::array;

constexpr int DPtc = 3;
