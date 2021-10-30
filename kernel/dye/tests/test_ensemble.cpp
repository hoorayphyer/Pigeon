#include "dye/ensemble_impl.hpp"
#include "testfw/testfw.hpp"

using namespace dye;

TEMPLATE_TEST_CASE("Create Trivial Ensemble From Cartesian Topology",
                   "[dye][mpi]"
                   // NOTE Notation: XxYxZ is the cartesian partition. The
                   // cartesian topology is periodic in all directions
                   ,
                   (aio::IndexType<-1, -1>), (aio::IndexType<-2, -1>),
                   (aio::IndexType<-1, -2>), (aio::IndexType<-2, -2>),
                   (aio::IndexType<-4, -4>), (aio::IndexType<-8, -8>),
                   (aio::IndexType<-16, -16>), (aio::IndexType<-32, -32>)) {
  constexpr int DGrid = 2;

  auto cart_opt = aio::make_cart(TestType::get(), mpi::world);
  auto ens_opt = create_ensemble<DGrid>(cart_opt);
  REQUIRE(bool(cart_opt) == bool(ens_opt));
  if (ens_opt) {
    const auto& cart = *cart_opt;
    const auto& ens = *ens_opt;
    REQUIRE(ens.intra.size() == 1);
    for (int i = 0; i < DGrid; ++i) {
      auto [src, dest] = cart.shift(i);
      if (TestType::get()[i] == -1) {
        REQUIRE_FALSE(ens.inter[i][LFT]);
        REQUIRE_FALSE(ens.inter[i][RGT]);
      } else {
        REQUIRE(bool(src) == bool(ens.inter[i][LFT]));
        REQUIRE(bool(dest) == bool(ens.inter[i][RGT]));
      }
    }

    REQUIRE(ens.chief == 0);
    REQUIRE(ens.chief_cart_rank == cart.rank());
    auto [c, d, p] = cart.coords_dims_periodic();
    for (int i = 0; i < DGrid; ++i) {
      REQUIRE(ens.cart_coords[i] == c[i]);
      REQUIRE(ens.cart_dims[i] == d[i]);
      REQUIRE(ens.is_periodic[i] == p[i]);
    }
  }
  mpi::world.barrier();
}

SCENARIO("Test intercommunicator created in the ensemble in 2x1 topology",
         "[dye][mpi]") {
  constexpr int DGrid = 2;

  auto cart_opt = aio::make_cart(apt::array<int, DGrid>{2, 1}, mpi::world);
  auto ens_opt = create_ensemble<DGrid>(cart_opt);
  if (ens_opt) {
    const auto& cart = *cart_opt;
    auto& ens = *ens_opt;
    if (cart.rank() == 0) {
      REQUIRE(bool(ens.inter[0][1]));
      const auto& comm = *(ens.inter[0][1]);
      int X = -565;
      auto req = comm.Isend(0, 147, &X, 1);
      mpi::wait(req);
    } else {
      REQUIRE(bool(ens.inter[0][0]));
      const auto& comm = *(ens.inter[0][0]);
      int X = 0;
      auto req = comm.Irecv(0, 147, &X, 1);
      mpi::wait(req);
      REQUIRE(X == -565);
    }
  }
  mpi::world.barrier();
}

TEMPLATE_TEST_CASE("Create one ensemble with nonzero replicas", "[dye][mpi]",
                   (std::integral_constant<int, 2>),
                   (std::integral_constant<int, 10>),
                   (std::integral_constant<int, 100>),
                   (std::integral_constant<int, 1000>)) {
  constexpr auto EnsSize = TestType::value;
  if (mpi::world.size() >= EnsSize) {
    constexpr int DGrid = 2;
    std::vector<int> cart_dims = {1, 1};
    std::vector<bool> periodic = {false, false};

    // first create a trivial ensemble
    auto cart_opt = aio::make_cart(cart_dims, periodic, mpi::world);
    auto intra = mpi::world.split((mpi::world.rank() < EnsSize));

    auto ens_opt = create_ensemble<DGrid>(cart_opt, intra);
    if (mpi::world.rank() < EnsSize) {
      const auto& cart = *cart_opt;
      const auto& ens = *ens_opt;

      REQUIRE(ens.intra.size() == EnsSize);
      REQUIRE(ens.chief == 0);
      REQUIRE(ens.chief_cart_rank == 0);
      for (int i = 0; i < DGrid; ++i) {
        REQUIRE(ens.cart_coords[i] == 0);
        REQUIRE(ens.cart_dims[i] == 1);
        REQUIRE(ens.is_periodic[i] == false);
      }
    }
  }
}
