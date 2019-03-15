#include "parallel/dynamic_balance.hpp"
#include "parallel/mpi++.hpp"

// need intercomm between employers and employees.
// when one party is done, the leader will notify the other leader, who then broadcast this to the members
int assign_ensemble_label() {
  std::optional<mpi::Comm> employers;
  std::optional<mpi::Comm> employees;
  std::optional<mpi::InterComm> job_market; // connecting these two
  // TODO check optionals

  int target_load;

  const int offer_tag = 835;

  auto f_employer =
    [&](){
      // Find out number of employees needed.
      int my_load [2] = {0,1}; // TODO need local ensemble. ave_load, ens_size
      // First, all employers all_gather the load.
      auto loads = employers->allgather(my_load, 2);

      unsigned long total_load = 0;
      for ( int i = 0; i < loads.size() / 2; ++i )
        total_load += static_cast<unsigned long>(loads[2*i]) * loads[2*i+1];
      int ave_load = total_load / mpi::world.size() + 1; // +1 to enforce an attainable ave_load
      ave_load = ( ave_load > target_load ) ? ave_load : target_load;

      for ( int i = 0; i < loads.size() / 2; ++i ) {
        // store in loads[2*i] the deficit in number of nodes. If a negative value, the ensemble will fire an employee
        loads[2*i] = ( static_cast<unsigned long>(loads[2*i]) * loads[2*i+1] ) / ave_load;
        loads[2*i] -= loads[2*i+1];
        // NOTE taking the floor in 1st line ensures that this is an attainable deployment
        // TODO there maybe some residual processes left unused. Fine a new scheme to include them in
      }

      int num_offers = 0;
      for ( int i = 0; i < loads.size() / 2; ++i ) num_offers += loads[2*i];

      // a leader in the employers is selected to notify all idle employees how many of them will be recruited
      if ( job_markdet.remote_size() > 0 ) {
        int root = ( 0 == employers->rank() ) ? MPI_ROOT : MPI_PROC_NULL;
        job_markdet.broadcast( &num_offers, 1, root );
      }

      // TODO here prmies need to do something

    };

  auto f_employee =
    [&]() {
      int num_offers = 0;
      job_market.broadcast(&num_offers, 1, 0);
      if ( employees.rank() >= num_offers ) return;

      int offer = 0;
      job_market.recv(MPI_ANY_SOURCE, offer_tag, &offer, 1);
    };

}

namespace parallel::impl {
  struct Rebalancer {
  private:
    struct rebalance_code {
      static constexpr int send = 1;
      static constexpr int none = 0;
      static constexpr int recv = -1;
    };

    // NOTE: number to be adjusted = actual number - expected number
    std::vector<int> AverageAll ( const std::vector<unsigned long>& nums_ptc ) {
      std::vector<int> nums_adjust( nums_ptc.size() );
      unsigned long num_tot = 0;
      for ( auto x : nums_ptc ) num_tot += x;

      unsigned long average = num_tot / nums_ptc.size();
      int remain = num_tot % nums_ptc.size();
      for ( unsigned int rank = 0; rank < nums_ptc.size(); ++rank ) {
        unsigned long num_exp = average + ( rank < remain );
        nums_adjust[rank] = nums_ptc[rank] - num_exp;
      }

      return nums_adjust;
    }

    // The ith element of nums_adjust, which should be ordered by ensemble ranks, specifies the number of particles ( could be positive, zero, or negative ) to be adjusted for the process i. Postive indicates sending. The constraint is that the sum of nums_adjust elements must be zero in order to conserve total particle number. The return value contains the actual sendrecv scheme for the corresponding process. Each scheme is of the format, action, member1, number1, member2, number2, etc, where numbers here are all positive.
    std::vector<int> get_instr( const std::vector<int>& nums_adjust, int my_rank ) {
      std::vector<std::vector<int>> instructions ( nums_adjust.size() );
      // split the input array into positive(send), zero, and negative(recv)
      std::vector<int> procs_send;
      std::vector<int> procs_recv;
      for ( unsigned int i = 0; i < nums_adjust.size(); ++i ) {
        if ( 0 == nums_adjust[i] )
          instructions[i] = { rebalance_code::none };
        else if ( nums_adjust[i] > 0 ) {
          instructions[i] = { rebalance_code::send };
          procs_send.push_back(i);
        }
        else {
          instructions[i] = { rebalance_code::recv };
          procs_recv.push_back(i);
        }
      }

      // match postives with negatives
      auto it_s = procs_send.begin();
      auto it_r = procs_recv.begin();
      int num_s = it_s != procs_send.end() ? std::abs( nums_adjust[*it_s] ) : 0;
      int num_r = it_r != procs_recv.end() ? std::abs( nums_adjust[*it_r] ) : 0;

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
          num_s = std::abs( nums_adjust[*it_s] );
        if ( 0 == num_r && (++it_r) != procs_recv.end() )
          num_r = std::abs( nums_adjust[*it_r] );
      }

      return instructions[my_rank];

    }

    template < typename Ptc_t >
    void exec_instr ( std::vector<Ptc_t>& ptcs, const mpi::Comm& ensemble, const std::vector<int>& instr ) {
      if ( rebalance_code::none == scheme[0] ) return;

      // Start from the back of the particle array for both send and recv actions
      int position = ptcs.size();
      const int num_comms = ( instr.size() - 1 ) / 2; // number of communications to make

      std::vector<mpi::Request> reqs(num_comms);

      for ( int i = 0; i < num_comms; ++i ) {
        int other_rank = instr[ 2 * i + 1 ];
        int num = instr[ 2 * i + 2 ];
        int tag = 147;
        if ( rebalance_code::send == scheme[0] ) {
          position -= num;
          reqs[i] = ensemble.Isend( other_rank, tag, ptcs.data() + position, num );
        } else {
          // TODO double check if using Irecv/Isend would cause race condition?
          // NOTE: using Irecv here causes hanging behavior on some platforms. So we use recv.
          // ensemble.Irecv( other_rank, tag, ptcs.PtcData().data() + position, num, requests.data() + i );
          ensemble.recv( other_rank, tag, ptcs.data() + position, num );
          position += num;
        }
      }
      mpi::waitall(reqs);

      if ( rebalance_code::send == scheme[0] )
        ptcs.erase( position, ptcs.size() - position );

      ptcs.resize( position );
    }

  public:
    template < typename Ptc_t >
    void balance_out ( std::vector<Ptc_t>& ptcs, const mpi::Comm& ensemble ) {
      unsigned long my_load = ptcs.size();
      auto loads = ensemble.allgather(&my_load, 1);
      auto instr = reb.get_instr( reb.AverageAll( loads ), ensemble.rank() );
      reb.exec_instr( ptcs, ensemble, instr );
    }
  };
}

namespace parallel {
  void dynamic_adjust() {
    std::optional<mpi::Comm> cart_opt;
    std::optional<mpi::Comm> ens_opt;
    std::optional<Locale<DGrid>> loc_opt; // really locale and ensemble are same concept

    int new_label = 0;
    { // Step 1. based on particle load, assign ensemble labels to each process. // TODO label for the idle group?
      unsigned long my_load = 0; // a process within an ensemble
      // // get particle load TODO
      // std::unordered_map<ParticleType, unsigned long> N_ptcs;
      // for ( const auto& elm : _data.particles )
      //   load += elm.second.Number();
      // load += _data.photons.Number() / 3;// empirically photon performance is about 3 times faster than that of a particle

      if ( ens_opt ) {
        const auto& locale;
        ens_opt -> reduce<mpi::by::SUM, mpi::IN_PLACE>(&my_load, 1, locale.chief);
      }
      // RATIONALE
      // - leaving processes start from bottom
      new_label = assign_ensemble_label(); // everyone participates
    }

    { // Step 2. Leaving processes clear up
      if ( ens_opt ) {
        bool is_leaving = ( new_label != locale.lable );
        // split the ensemble into staying and leaving
        auto split_opt = ens_opt -> split( is_leaving, ens_opt->rank() );
        if ( split_opt->size() != ens_opt->size() ) {
          int local_leader = is_leaving ? ens_opt->size() - 1 : 0;
          int remote_leader = ens_opt_size() - local_leader - 1;
          auto inter_comm = mpi::InterComm(*split_opt, local_leader, ens_opt, remote_leader, 315);

          // auto f_reb_empty_leavers; // TODO use alltoallv to scatter
          // for ( auto& elm : data.particles ) {
          //   auto& ptc = elm.second;
          //   f_reb_empty_leavers( ptc );
          // }
          // f_reb_empty_leavers( data.photons );

          // reduce fields calculated on local particle information accummulated over timesteps onto the locale chief.
          // Basically its pair creation events. Don't forget to zero out fields on leaving processes
          // TODO: Can we do DA after data export to avoid this? I think so.
        }
      }

    }

    { // Step 3. Create new ensemble communicator and update domain members.
      auto new_ensemble = mpi::world.split( new_label, mpi::world.rank() );
      if ( new_label != IDLE_LABEL ) {
        // update locale
        // update local grid
        // update neighbors intercomm
        // update ensemble
      } else {
        // clear the above
        // clear all data
        // must be fresh clean
      }
    }

    { // Step 4. Perform a complete rebalance on all ensembles together
      // TODO transform to std::vector<CParticle>
      // if ( new_ensemble_exists ) {
      //   impl::Rebalancer reb;
      //   for ( auto& elm : data.particles ) {
      //     auto& ptc = elm.second;
      //     reb.balance_out( ptc, ensemble );
      //   }
      //   reb.balance_out( photons, ensemble );
      // }
    }

    { // Step 5. After locale is updated, sort particles and photons in case of actions such as annihilation FIXME is this needed?
    }

  }
}
