#include "mpipp/mpi++.hpp"
#include "testfw/testfw.hpp"

using namespace mpi;

SCENARIO("World", "[parallel][mpi]") {
  if (world.size() >= 2) {
    int myrank = world.rank();
    // std::cout << "world size = " << world.size() << " ," << myrank <<
    // std::endl;
    if (0 == myrank) {
      int msg = 147;
      world.send(1, 22, &msg);
    } else if (1 == myrank) {
      int msg;
      world.recv(0, 22, &msg, 1);
      REQUIRE(msg == 147);
    }
  }
  world.barrier();
}

// SCENARIO("Group", "[parallel][mpi]") {
//   auto grp_world = world.group();
//   SECTION("Constructors") {
//     WHEN("Default") {
//       Group grp;
//       REQUIRE( MPI_UNDEFINED == grp.rank() );
//       REQUIRE( 0 == grp.size() );
//     }

//     WHEN("From MPI_WORLD_GROUP") {
//       AND_WHEN("no duplicates") {
//       }
//       AND_WHEN("there are duplicates") {
//       }
//     }
//   }

//   SECTION("Operations") {

//   }

// }

SCENARIO("Comm", "[parallel][mpi]") {
  auto comm_opt = world.split(true);
  REQUIRE(comm_opt);
  auto& comm = *comm_opt;
  if (world.size() >= 2) {
    int myrank = comm.rank();
    if (0 == myrank) {
      int msg = 147;
      comm.send(1, 22, &msg);
    } else if (1 == myrank) {
      int msg;
      comm.recv(0, 22, &msg, 1);
      REQUIRE(147 == msg);
    }
  }
  world.barrier();
}

// TODOL restore after typelist is installed
// SCENARIO("Datatype", "[parallel][mpi]") {
// #define TestType(_T_, _MPIT_)                        \
//   REQUIRE( datatype((_T_*)0) == MPI_##_MPIT_ );      \
//   REQUIRE( datatype((const _T_*)0) == MPI_##_MPIT_); \
//   REQUIRE( datatype((_T_)0) == MPI_##_MPIT_);        \

//   TestType(int, INT);
//   TestType(char, CHAR);
//   TestType(short, SHORT);
//   TestType(long, LONG);

//   TestType(unsigned char, UNSIGNED_CHAR);
//   TestType(unsigned short, UNSIGNED_SHORT);
//   TestType(unsigned int, UNSIGNED);
//   TestType(unsigned long, UNSIGNED_LONG);

//   TestType(float, FLOAT);
//   TestType(double, DOUBLE);
//   TestType(bool, CXX_BOOL);
// #undef TestType
// }

SCENARIO("intra P2P communications", "[parallel][mpi]") {
  SECTION("Nonblocking") {
    if (world.size() >= 2) {
      SECTION("one sided communication, test on wait") {
        if (world.rank() == 0) {
          int send_msg = 351;
          Request req = world.Isend(1, 147, &send_msg, 1);
          mpi::wait(req);
        } else if (world.rank() == 1) {
          int recv_msg;
          Request req = world.Irecv(0, 147, &recv_msg, 1);
          mpi::wait(req);
          REQUIRE(recv_msg == 351);
        }
        world.barrier();
      }

      SECTION("double sided communication, test on waitall") {
        if (world.rank() < 2) {
          std::vector<Request> reqs(2);
          int myrank = world.rank();
          int send_msg = myrank + 147;
          reqs[0] = world.Isend(1 - myrank, 147, &send_msg, 1);
          int recv_msg = -1;
          reqs[1] = world.Irecv(1 - myrank, 147, &recv_msg, 1);
          waitall(reqs);
          REQUIRE(send_msg == myrank + 147);
          REQUIRE(recv_msg == 1 - myrank + 147);
        }
        world.barrier();
      }

      SECTION("test wait on null requests") {
        if (world.rank() == 0) {
          Request req;
          mpi::wait(req);

          std::vector<Request> reqs(2);
          mpi::waitall(reqs);
        }
        world.barrier();
      }
    }
  }
}

SCENARIO("intra collective communications", "[parallel][mpi]") {
  SECTION("allreduce") {
    auto rw_opt = aio::reduced_world(4, world);
    if (rw_opt) {
      auto& rw = *rw_opt;
      WHEN("out of place") {
        std::vector<int> buf(100, rw.rank());
        auto res = rw.allreduce(mpi::by::SUM, buf.data(), buf.size());
        REQUIRE(bool(res));
        for (auto x : buf) REQUIRE(x == rw.rank());
        for (auto x : *res) REQUIRE(x == rw.size() * (rw.size() - 1) / 2);
      }
      WHEN("in place") {
        std::vector<int> buf(100, rw.rank());
        auto res =
            rw.allreduce<mpi::IN_PLACE>(mpi::by::SUM, buf.data(), buf.size());
        REQUIRE(!res);
        for (auto x : buf) REQUIRE(x == rw.size() * (rw.size() - 1) / 2);
      }
    }
    world.barrier();
  }

  SECTION("broadcast") {
    auto rw_opt = aio::reduced_world(4, world);
    if (rw_opt) {
      auto& rw = *rw_opt;
      std::vector<int> buf(100, rw.rank());
      const int root = 2;
      rw.broadcast(root, buf.data(), buf.size());
      for (auto x : buf) REQUIRE(x == root);
    }
    world.barrier();
  }

  SECTION("Ibroadcast") {
    auto rw_opt = aio::reduced_world(4, world);
    if (rw_opt) {
      auto& rw = *rw_opt;
      std::vector<int> buf(100, rw.rank());
      const int root = 2;
      auto req = rw.Ibroadcast(root, buf.data(), buf.size());
      wait(req);
      for (auto x : buf) REQUIRE(x == root);
    }
    world.barrier();
  }
}

SCENARIO("intercomm P2P communications", "[parallel][mpi]") {
  SECTION("each comm has one member") {
    if (world.size() >= 2 && world.rank() < 2) {
      int myrank = world.rank();
      InterComm inter(self, 0, {world}, 1 - myrank, 147);
      if (0 == myrank) {
        int msg = 147;
        inter.send(0, 22, &msg, 1);
      } else if (1 == myrank) {
        int msg;
        inter.recv(0, 22, &msg, 1);
        REQUIRE(msg == 147);
      }
    }
    world.barrier();
  }
}

// TODO
// SCENARIO("inter collective communications", "[parallel][mpi]") {
//   SECTION("allreduce") {
//     auto rw_opt = aio::reduced_world(4,world);
//     if(rw_opt) {
//       auto& rw = *rw_opt;
//       std::vector<int> buf(100, rw.rank());
//       auto res = rw.allreduce(mpi::by::SUM, buf.data(), buf.size());
//       REQUIRE(bool(res));
//       for ( auto x : buf )
//         REQUIRE(x == rw.rank());
//       for ( auto x : *res )
//         REQUIRE(x == rw.size() * (rw.size() - 1) / 2);
//     }
//     world.barrier();
//   }
// }

SCENARIO("cartesian communicator", "[parallel][mpi]") {
  WHEN("1x1 and periodic") {
    if (world.rank() == 0) {
      CartComm cart(self, {1, 1}, {true, true});
      for (int i = 0; i < 2; ++i) {
        auto [src_rank, dest_rank] = cart.shift(i, 1);
        REQUIRE(src_rank);
        REQUIRE(dest_rank);
        REQUIRE(*src_rank == 0);
        REQUIRE(*dest_rank == 0);
      }
    }
    world.barrier();
  }
}

SCENARIO("cartesian sub", "[parallel][mpi]") {
  constexpr int Nr = 2;
  constexpr int Nth = 4;
  auto rw_opt = aio::reduced_world(Nr * Nth, world);
  if (rw_opt) {
    auto& rw = *rw_opt;

    CartComm cart(rw, {Nr, Nth}, {true, true});

    std::vector<int> coords_bottom_left{0, 0};
    std::vector<int> coords_top_right{1, 3};

    auto sub_cart = cart.sub(coords_bottom_left, coords_top_right);
    std::vector<int> msg{-1, -1, -1, -1, -1, -1, -1};
    msg[0] = bool(sub_cart);
    {
      auto c = cart.coords();
      msg[1] = c[0];
      msg[2] = c[1];
    }
    if (sub_cart) {
      std::vector<int> crds;
      std::vector<Topo> topos;
      std::tie(crds, topos) = sub_cart->coords_topos();
      msg[3] = crds[0];
      msg[4] = crds[1];
      msg[5] = topos[0].signed_dim();
      msg[6] = topos[1].signed_dim();
    }

    auto buf_opt = cart.gather(0, msg.data(), msg.size());
    REQUIRE(bool(buf_opt) == bool(cart.rank() == 0));
    if (buf_opt) {
      const auto& buf = *buf_opt;
      for (int i = 0; i < cart.size(); ++i) {
        const auto* p = buf.data() + i * msg.size();
        const auto is_in = p[0];
        const auto cr = p[1];
        const auto cth = p[2];
        const auto sub_cr = p[3];
        const auto sub_cth = p[4];
        const auto sub_topo_r = p[5];
        const auto sub_topo_th = p[6];
        REQUIRE(i / Nth == cr);
        REQUIRE(i % Nth == cth);
        REQUIRE(is_in ==
                (coords_bottom_left[0] <= cr and cr < coords_top_right[0] and
                 coords_bottom_left[1] <= cth and cth < coords_top_right[1]));
        if (not is_in) {
          REQUIRE(sub_cr == -1);
          REQUIRE(sub_cth == -1);
          REQUIRE(sub_topo_r == -1);
          REQUIRE(sub_topo_th == -1);
        } else {
          REQUIRE(sub_cr == cr - coords_bottom_left[0]);
          REQUIRE(sub_cth == cth - coords_bottom_left[1]);
          REQUIRE(sub_topo_r == coords_top_right[0] - coords_bottom_left[0]);
          REQUIRE(sub_topo_th == coords_top_right[1] - coords_bottom_left[1]);
        }
      }
    }
  }
  world.barrier();
}
