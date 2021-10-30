#ifndef _SILOPP_PMPIO_HPP_
#define _SILOPP_PMPIO_HPP_

#include "silopp/silo++.hpp"

namespace mpi {
struct Comm;
}

namespace silo {
struct Pmpio {
  std::optional<mpi::Comm> comm;
  std::string filename;
  std::string dirname;
  Mode mode = Mode::Write;

  template <typename F>
  void operator()(const F& f) const {
    if (comm) {
      const auto myrank = comm->rank();
      if (myrank != 0) {
        int baton = 0;
        comm->recv(myrank - 1, 147, &baton, 1);
      }

      {  // force cleanup dbfile before possible mpi::send
        auto dbfile = open(filename, mode);
        if (Mode::Write == mode && !dbfile.exists(dirname)) {
          dbfile.mkdir(dirname);
        }
        // set current directory within the silo file
        dbfile.cd(dirname);

        f(dbfile);
      }

      if (myrank != comm->size() - 1) {
        int baton = 0;
        comm->send(myrank + 1, 147, &baton, 1);
      }

      comm->barrier();  // NOTE this barrier is necessary to prevent early
                        // finished process from reentering from a later while
                        // other processes haven't finished the present one.
    }
  }
};

}  // namespace silo

#endif
