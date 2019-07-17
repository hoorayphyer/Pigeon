#include "dye/dynamic_balance.hpp"
#include "mpipp/mpi++.hpp"
#include "particle/particle.hpp"
#include "apt/priority_queue.hpp" // used in calc_new_nprocs

// bifurcate
namespace dye::impl {
  auto bifurcate( const mpi::Comm& parent, bool color ) {
    mpi::Comm child ( *(parent.split(color)) );
    std::optional<mpi::InterComm> itc;
    if ( child.size() == parent.size() ) return std::make_tuple( child, itc );

    auto is_found =
      []( const auto& ranks, int r ) {
        for ( const auto elm : ranks )
          if ( elm == r ) return true;
        return false;
      };

    // NOTE assume ranks_in_child.size() < ranks_in_parent.size()
    auto find_first_in_other =
      [is_found] ( const auto& ranks_in_child, auto parent_size ) {
        for ( int i = 0; i < parent_size; ++i )
          if ( !is_found(ranks_in_child, i) ) return i;
      };

    // RATIONALE the one with smallest rank in parent in each branch becomes the leader.
    int local_leader = 0; // NOTE local_leader is the rank in the child comm], using parent.rank() in calling parent.split guarantees the correctness
    int remote_leader = find_first_in_other(child.group().translate_all(parent.group()), parent.size());
    itc.emplace(mpi::InterComm(child, local_leader, parent, remote_leader, 80));

    return std::make_tuple( child, itc );
  }
}

// calc_nprocs_new
namespace dye::impl {
  // loads = [ ens_tot_load_0, ens_size_0, ens_tot_load_1, ens_size_1... ]
  std::vector<int> calc_new_nprocs( const std::vector<particle::load_t>& loads_and_nprocs, particle::load_t target_load, const unsigned int max_nprocs ) {
    const auto nens = loads_and_nprocs.size() / 2;
    std::vector<int> nproc;
    nproc.reserve(nens);
    nproc.resize(nens, 1);

    apt::priority_queue<long double> pq; // store average load in pq
    // start with one process per ensemble
    for ( int i = 0; i < nens; ++i )
      pq.push( i, static_cast<long double>(loads_and_nprocs[2*i]) );

    int total_nprocs = 0;
    for ( int i = 0; i < nens; ++i )
      total_nprocs += loads_and_nprocs[2*i+1];
    total_nprocs = std::max<int>(total_nprocs, max_nprocs);

    int nfree = std::max<int>( total_nprocs - nens, 0 );

    while ( nfree ) {
      auto[i, l] = pq.top();
      pq.pop();
      if ( l <= target_load ) break;
      ++nproc[i];
      --nfree;
      pq.push( i, static_cast<long double>(loads_and_nprocs[2*i]) / nproc[i] );
    }

    return nproc;
  }
}

// relinguish
namespace dye::impl {
  template < typename T, template < typename > class PtcSpecs >
  void relinguish_data( particle::array<T, PtcSpecs>& ptc_array,
                        const mpi::InterComm& itc, int root ) {
    std::vector<particle::Particle<T,PtcSpecs>> buf;
    int count = 0;
    if ( root == MPI_ROOT ) {
      const auto size = ptc_array.size();
      count = ( size / itc.remote_size() ) + ( size % itc.remote_size() != 0 );
      itc.broadcast( root, &count, 1 );
      buf.reserve( count * itc.remote_size() );
    } else {
      itc.broadcast( root, &count, 1 );
      if ( root != MPI_PROC_NULL )
        buf.reserve( count );
    }

    buf.resize( buf.capacity() );

    if ( root == MPI_ROOT ) {
      for ( int i = 0; i < ptc_array.size(); ++i )
        buf[i] = std::move(ptc_array[i]);
      ptc_array.resize(0);
    }

    itc.scatter(root, buf.data(), count );

    if ( root != MPI_ROOT && root != MPI_PROC_NULL ) {
      auto old_size = ptc_array.size();
      ptc_array.resize( old_size + count );
      for ( int i = 0; i < count; ++i )
        ptc_array[ old_size + i ] = buf[i];
    }
  }

}

// assign_labels
namespace dye::impl{
  std::optional<int> assign_labels( const std::optional<mpi::InterComm>& job_market_opt,
                                    const std::vector<int>& deficits,
                                    std::optional<int> cur_label ) {
    // RATIONALE job_market contains all primaries and all idles including those just fired from shrinking ensembles. The goal here is for all primaries to coordinate whom to send their labels to.
    // NOTE `deficits` should be valid already. This means
    //   - deficits[i] >= 0
    //   - the sum of deficits <= num_idles. In the case of inequality, idles will be picked by their rank in the job_market.
    // Edge cases:
    //   - no idles
    std::optional<int> new_label;
    if ( cur_label ) // primary
      new_label = cur_label;

    if ( !job_market_opt ) return new_label;
    const auto& job_market = *job_market_opt;
    const bool is_prmy {cur_label};

    int num_offers = 0;
    if ( is_prmy ) {
      for ( int i = 0; i < deficits.size(); ++i ) num_offers += deficits[i];

      int root = ( 0 == job_market.rank() ) ? MPI_ROOT : MPI_PROC_NULL;
      job_market.broadcast( root, &num_offers, 1 );
    } else {
      job_market.broadcast(0, &num_offers, 1 );
    }


    if ( cur_label ) { // primary

      int myrank = job_market.rank();
      int num_offers_ahead = 0;
      for ( int i = 0; i < myrank; ++i )
        num_offers_ahead += deficits[i];

      int label = *cur_label;
      std::vector<mpi::Request> reqs(deficits[myrank]);
      for ( int i = 0; i < deficits[myrank]; ++i )
        reqs[i] = job_market.Isend(num_offers_ahead + i, 835, &label, 1 );

      mpi::waitall(reqs);

    } else {
      if ( job_market.rank() < num_offers ) {
        int label = 0;
        job_market.recv(MPI_ANY_SOURCE, 835, &label, 1);
        new_label.emplace(label);
      }
    }

    return new_label;
  };
}

namespace dye::impl {
  struct balance_code {
    static constexpr int send = 1;
    static constexpr int none = 0;
    static constexpr int recv = -1;
  };

  std::vector<int> get_ptc_num_surplus ( const std::vector<particle::load_t>& nums_ptc ) {
    // NOTE: surplus = actual number - expected number
    std::vector<int> spls( nums_ptc.size() );
    particle::load_t num_tot = 0;
    for ( auto x : nums_ptc ) num_tot += x;

    particle::load_t average = num_tot / nums_ptc.size();
    int remain = num_tot % nums_ptc.size();
    for ( unsigned int rank = 0; rank < nums_ptc.size(); ++rank ) {
      particle::load_t num_exp = average + ( rank < remain );
      spls[rank] = nums_ptc[rank] - num_exp;
    }

    return spls;
  }

  // The ith element of nums_adjust, which should be ordered by ensemble ranks, specifies the number of particles ( could be positive, zero, or negative ) to be adjusted for the process i. Postive indicates sending. The constraint is that the sum of nums_adjust elements must be zero in order to conserve total particle number. The return value contains the actual sendrecv scheme for the corresponding process. Each scheme is of the format, action, member1, number1, member2, number2, etc, where numbers here are all positive.
  std::vector<int> get_instr( const std::vector<int>& nums_surplus, int my_rank ) {
    std::vector<std::vector<int>> instructions ( nums_surplus.size() );
    // split the input array into positive(send), zero, and negative(recv)
    std::vector<int> procs_send;
    std::vector<int> procs_recv;
    for ( unsigned int i = 0; i < nums_surplus.size(); ++i ) {
      if ( 0 == nums_surplus[i] )
        instructions[i] = { balance_code::none };
      else if ( nums_surplus[i] > 0 ) {
        instructions[i] = { balance_code::send };
        procs_send.push_back(i);
      }
      else {
        instructions[i] = { balance_code::recv };
        procs_recv.push_back(i);
      }
    }

    // match postives with negatives
    auto it_s = procs_send.begin();
    auto it_r = procs_recv.begin();
    int num_s = it_s != procs_send.end() ? std::abs( nums_surplus[*it_s] ) : 0;
    int num_r = it_r != procs_recv.end() ? std::abs( nums_surplus[*it_r] ) : 0;

    while ( it_s != procs_send.end() && it_r != procs_recv.end() ) {
      auto& instr_s = instructions[ *it_s ];
      auto& instr_r = instructions[ *it_r ];
      int num_transfer = std::min( num_s, num_r );
      instr_s.push_back( *it_r );
      instr_s.push_back(num_transfer);
      instr_r.push_back( *it_s );
      instr_r.push_back(num_transfer);

      num_s -= num_transfer;
      num_r -= num_transfer;
      if ( 0 == num_s && (++it_s) != procs_send.end() )
        num_s = std::abs( nums_surplus[*it_s] );
      if ( 0 == num_r && (++it_r) != procs_recv.end() )
        num_r = std::abs( nums_surplus[*it_r] );
    }

    return instructions[my_rank];

  }

}

// detailed balance
namespace dye {
  template < typename T, template < typename > class PtcSpecs >
  void detailed_balance ( particle::array<T, PtcSpecs>& ptcs,
                          const mpi::Comm& intra ) {
    particle::load_t my_load = ptcs.size();
    auto loads = intra.allgather(&my_load, 1);
    auto instr = impl::get_instr( impl::get_ptc_num_surplus( loads ), intra.rank() );

    if ( impl::balance_code::none == instr[0] ) return;

    // Start from the back of the particle array for both send and recv actions
    int position = ptcs.size();
    const int num_comms = ( instr.size() - 1 ) / 2; // number of communications to make

    std::vector<mpi::Request> reqs(num_comms);

    std::vector<particle::Particle<T,PtcSpecs>> buffer;

    for ( int i = 0; i < num_comms; ++i ) {
      int other_rank = instr[ 2 * i + 1 ];
      int num = instr[ 2 * i + 2 ];
      buffer.resize( num );
      int tag = 147;
      if ( impl::balance_code::send == instr[0] ) {
        position -= num;
        for ( int j = 0; j < num; ++j ) buffer[j] = std::move(ptcs[position + j]);
        reqs[i] = intra.Isend( other_rank, tag, buffer.data(), num );
      } else {
        // TODOL double check if using Irecv/Isend would cause race condition?
        // NOTE: using Irecv here causes hanging behavior on some platforms. So we use recv.
        intra.recv( other_rank, tag, buffer.data(), num );
        ptcs.resize(ptcs.size() + num);
        for ( int j = 0; j < num; ++j ) ptcs[position + j] = std::move(buffer[j]);
        position += num;
      }
    }
    mpi::waitall(reqs);
    if ( impl::balance_code::send == instr[0] ) ptcs.resize( position );
  }
}

// dynamic_load_balance
namespace dye {
  template < int DGrid >
  std::optional<Ensemble<DGrid>>
  deploy ( particle::load_t my_tot_load,
           const std::optional<Ensemble<DGrid>>& ens_opt_old,
           const std::optional<mpi::CartComm>& cart_opt,
           unsigned int target_load ) {
    std::optional<Ensemble<DGrid>> ens_opt_new;

    // NOTE deficit = desired number - current number
    std::vector<int> nproc_deficit; // significant only at primaries
    bool is_active = ens_opt_old;

    { // Step 1. based on particle load, primaries figure out the surpluses. If an ensemble has surplus > 0, the primary will flag the last few processes to be leaving the ensemble.
      if ( ens_opt_old ) {
        const auto& ens = *ens_opt_old;
        ens.intra.template reduce< mpi::IN_PLACE >( mpi::by::SUM, ens.chief, &my_tot_load, 1 );
        int new_ens_size = 0;
        if ( cart_opt ) {
          particle::load_t my_load[2] = { my_tot_load, ens.intra.size() };
          auto loads_and_nprocs = cart_opt->allgather(my_load, 2);
          nproc_deficit = impl::calc_new_nprocs( loads_and_nprocs, target_load, mpi::world.size() );
          new_ens_size = nproc_deficit[cart_opt->rank()];
          for ( int i = 0; i < nproc_deficit.size(); ++i )
            nproc_deficit[i] -= loads_and_nprocs[2*i+1];
        }
        ens.intra.broadcast(ens.chief, &new_ens_size, 1);
        if ( ens.intra.rank() >= new_ens_size ) is_active = false;
        if ( cart_opt ) { // clear negative deficits
          for ( auto& x : nproc_deficit ) x = ( x > 0 ) ? x : 0;
        }
      }
    }

    { // Step 2. Assign ensemble labels to all processes and create new ensemble. Need an intercommunicator between all primaries and all idles
      std::optional<unsigned int> new_label;
      bool is_replica = is_active && !cart_opt;
      auto prmy_idle_comm = std::get<0>( impl::bifurcate( mpi::world, is_replica ) );
      if ( is_replica ) new_label.emplace(ens_opt_old->label());
      else {
        // treat primaries and idles here
        std::optional<unsigned int> cur_label;
        if ( cart_opt ) cur_label.emplace( ens_opt_old->label() );

        auto[comm, job_market] = impl::bifurcate( prmy_idle_comm, is_active );
        new_label = impl::assign_labels( job_market, nproc_deficit, cur_label );
      }
      auto new_ens_intra_opt = mpi::world.split( new_label );

      ens_opt_new = create_ensemble<DGrid>( cart_opt, new_ens_intra_opt );
    }

    return ens_opt_new;
  }

  template < typename T, template < typename > class PtcSpecs, int DGrid >
  void dynamic_load_balance( particle::map<particle::array<T, PtcSpecs>>& particles,
                             std::optional<Ensemble<DGrid>>& ens_opt,
                             const std::optional<mpi::CartComm>& cart_opt,
                             unsigned int target_load ) {
    std::optional<Ensemble<DGrid>> ens_opt_new;
    { // Deploy
      particle::load_t my_tot_load = 0;
      if ( ens_opt ) {
        for ( auto&[sp,ptcs] : particles ) {
          if ( sp == particle::species::photon )
            my_tot_load += ptcs.size() / 3; // empirically photon performance is about 3 times faster than that of a particle
          else
            my_tot_load += ptcs.size();
        }
      }
      ens_opt_new = deploy( my_tot_load, ens_opt, cart_opt, target_load );
    }

    { // Leaving processes clear up
      bool is_leaving = ens_opt && ( !ens_opt_new || ( ens_opt->label() != ens_opt_new->label() ) );
      if ( ens_opt ) {
        auto[comm, itc_opt] = impl::bifurcate( ens_opt->intra, is_leaving );
        if ( itc_opt ) {
          const auto& itc = *itc_opt;

          for ( auto&[sp, ptcs] : particles ) {
            if (is_leaving) {
              for ( int i = 0; i < itc.size(); ++i ) {
                int root = (itc.rank() == i) ? MPI_ROOT : MPI_PROC_NULL;
                impl::relinguish_data( ptcs, itc, root );
              }
            } else {
              for ( int i = 0; i < itc.remote_size(); ++i )
                impl::relinguish_data( ptcs, itc, i );
            }
          }
        }
      }
    }
    { // clear ens_opt properly, which needs 1. we do the following outside the above scope so as to make sure [comm, itc_opt] are deallocated, 2. freeing communicator is collective though implemented as local operation, so everyone needs to call.
      ens_opt.reset();
      ens_opt.swap(ens_opt_new);
    }

    { // Perform a complete rebalance on all ensembles together. NOTE TODOL touchcreate is the preferred way to initialize particle array of a species. But in communication, touchcreate is not possible unless how others are sending you can be known a priori. While one could use a general buffer for all species, here we assume all species are touch created before dynamic_load_balance.
      if ( ens_opt ) {
        for ( auto&[sp, ptcs] : particles ) {
          detailed_balance(ptcs, ens_opt->intra);
        }
      }
    }

  }
}
