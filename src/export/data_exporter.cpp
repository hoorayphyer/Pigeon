#include "export/io.hpp"
#include "export/silo++.hpp"
#include "parallel/communicator.hpp"
#include "utility/filesys.hpp"

#include "dynamic_variables.hpp"
#include "parameters.hpp"

// #include <chrono>

// using namespace std::chrono;

// namespace {
//   std::array<int, 3> lowOffset;
//   std::array<int, 3> hiOffset;
// }

// template < typename T, typename U >
// void SampleData (const MultiArray<T>& field, MultiArray<U>& sample, bool is_do_average) {
//   // assert(n >= 0 && n < VECTOR_DIM);
//   // int factor = _grid.reducedDim(n) / _dataGrid.reducedDim(n);
//   // need to assign zeros in all guard cells
//   const auto& grid = _domain -> grid();
//   const auto& dataGrid = _domain -> dataGrid();
//   sample.assign(static_cast<U>(0.0));
//   Vec3<int> factor;
//   for (int comp = 0; comp < VECTOR_DIM; comp++)
//     factor[comp] = grid.reducedDim(comp) / dataGrid.reducedDim(comp);

//   for (int k = 0; k < dataGrid.reducedDim(2); k++) {
//     int field_start2 = k * factor[2];
//     int field_end2 = std::min( (k+1) * factor[2], grid.reducedDim(2) );
//     for (int j = 0; j < dataGrid.reducedDim(1); j++) {
//       int field_start1 = j * factor[1];
//       int field_end1 = std::min( (j+1) * factor[1], grid.reducedDim(1) );
//       for (int i = 0; i < dataGrid.reducedDim(0); i++) {
//         int field_start0 = i * factor[0];
//         int field_end0 = std::min( (i+1) * factor[0], grid.reducedDim(0) );
//         // average over values in the subfield given by field modulo dataGrid
//         auto sum = static_cast<T>(0) + static_cast<U>(0); // force type promotion
//         for ( int k_f = field_start2; k_f < field_end2; ++k_f ) {
//           for ( int j_f = field_start1; j_f < field_end1; ++j_f ) {
//             for ( int i_f = field_start0; i_f < field_end0; ++i_f ) {
//               sum += field( i_f + grid.guard[0], j_f + grid.guard[1], k_f + grid.guard[2] );
//             }
//           }
//         }

//         sample(i + dataGrid.guard[0], j + dataGrid.guard[1], k + dataGrid.guard[2])
//           = static_cast<U>( is_do_average ? sum / ( (field_end2 - field_start2) * (field_end1 - field_start1) * (field_end0 - field_start0) ) : sum );

//         // sample(i + dataGrid.guard[0], j + dataGrid.guard[1], k + dataGrid.guard[2])
//         //   = static_cast<U>(
//         //                    field((i + 1) * factor[0] - 1 + grid.guard[0],
//         //                          (j + 1) * factor[1] - 1 + grid.guard[1],
//         //                          (k + 1) * factor[2] - 1 + grid.guard[2])
//         //                    );
//       }
//     }
//   }

//   // due to existence of replicas, doing SendGuardCells here is not always correct. So it is deferred.
//   //   _domain -> SendGuardCells(sample, dataGrid);

// }

// template < typename T, typename U >
// void SampleData (const ScalarField<T>& field, ScalarField<U>& sample, bool is_do_average=true) {
//   SampleData(field.data(), sample.data(), is_do_average);
// }

// template < typename T, typename U >
// void SampleData (const VectorField<T>& field, VectorField<U>& sample, bool is_do_average=true) {
//   for (int n = 0; n < VECTOR_DIM; n++) {
//     SampleData(field.data(n), sample.data(n), is_do_average);
//   }
// }

// template < typename FiniteDiff >
// void GetDiv( const AperData &data, PICOutputData &data_out, FiniteDiff& fd, std::array<bool,6> is_bdry ) {
//   data_out.DivE.assign(0.0);
//   data_out.DivB.assign(0.0);
//   const VectorField<Scalar>& Efield = data.Efield;
//   const VectorField<Scalar>& Bfield = data.Bfield;
//   fd.ComputeDivergence( Efield, _scalar_tmp, FieldType::ETYPE, is_bdry.data() );
//   SampleData( _scalar_tmp, data_out.DivE );
//   fd.ComputeDivergence( Bfield, _scalar_tmp, FieldType::BTYPE, is_bdry.data() );
//   SampleData( _scalar_tmp, data_out.DivB );
// }

// void AverageMomAndGamma( const AperData &data, PICOutputData &data_out ) {
//   // rationale. The total P and gamma are first obtained and sampled before doing the average. This saves some operation and is currently the only way in the case of ensemble.

//   VectorField<Scalar>& Pfield = _vector_tmp;
//   ScalarField<Scalar>& Gamma = _scalar_tmp;
//   for ( auto it = data.particles.cbegin(); it != data.particles.cend(); ++it ) {
//     auto ptcType = it->first;
//     Pfield.assign( 0.0 );
//     Gamma.assign( 0.0 );
//     // sum over all particles
//     const auto& particles = data.particles.at(ptcType);
//     for ( int idx = 0; idx < particles.Number(); ++idx ) {
//       if ( particles.IsEmpty(idx) ) continue;

//       const auto& p = particles.PtcData()[idx].p;
//       const auto& cell = particles.PtcData()[idx].cell;
//       Pfield.data(0) [cell] += p.x;
//       Pfield.data(1) [cell] += p.y;
//       Pfield.data(2) [cell] += p.z;
//       Gamma.data() [cell] += std::sqrt( 1.0 + p.dot(p) );
//     }

//     // FIXME: workaround for retaining the double precision of Gamma and Pfield. The most efficient way is to gather just the sampled fields. But that causes some loss of precision and hampers the comparison between results with different replica maps.
//     _domain -> EnsReduceFields( Gamma, Pfield );

//     // sample total P and gamma. Note that we sum the values when sampling. This is because later on we will divide P and gamma by ptcNumbers which is obtained by summation during sampling.
//     auto& Pfield_out = data_out.Pfields.at(ptcType);
//     auto& Gamma_out = data_out.Gammas.at(ptcType);
//     SampleData( Gamma,  Gamma_out, false );
//     SampleData( Pfield, Pfield_out, false );

//     // if ( _domain -> inEnsemble() ) {
//     //   _domain -> EnsGatherFields( Gamma_out, Pfield_out );
//     // }

//     // average. Only the root here is significant
//     // NOTE data_out.ptcNumbers on primary is assumed to have been reduced outside
//     const auto& ptcNum = data_out.ptcNumbers.at(ptcType);
//     for ( int i = 0; i < Pfield_out.gridSize(); ++i ) {
//       int count = static_cast<int>( ptcNum.data()[i] + 0.5 ); // ptcNum, which is a float, should be very close to whole integers
//       if ( 0 == count ) {
//         Pfield_out.data(0) [i] = 0.0;
//         Pfield_out.data(1) [i] = 0.0;
//         Pfield_out.data(2) [i] = 0.0;
//         Gamma_out.data() [i] = 0.0;
//       } else {
//         using T_out = typename std::remove_reference<decltype(Pfield_out)>::type::data_type;

//         Pfield_out.data(0) [i] /= static_cast<T_out>( count );
//         Pfield_out.data(1) [i] /= static_cast<T_out>( count );
//         Pfield_out.data(2) [i] /= static_cast<T_out>( count );
//         Gamma_out.data() [i] /= static_cast<T_out>( count );
//       }
//     }

//   }

// }

// void DataExporter::GetFlux( const AperData &data, PICOutputData &data_out ) {
//   ScalarField<Scalar>& flux = _scalar_tmp;
//   Domain& domain = *_domain;
//   const Grid& grid = _domain -> grid();
//   Index B1stagger = GetStagProperty( FieldType::BTYPE, 0 );
//   Scalar* sendbuffer = new Scalar [ grid.reducedDim(0) * grid.reducedDim(2) ];
//   Scalar* recvbuffer = new Scalar [ grid.reducedDim(0) * grid.reducedDim(2) ];
//   // NOTE the following implementation is for 2D LogSpherical. It is second order.
//   // integrate over theta
//   for ( int k = grid.guard[2]; k < grid.dims[2] - grid.guard[2]; ++k ) {
//     int buffer_index_phi = ( k - grid.guard[2] ) * grid.reducedDim(0);
//     for ( int i = grid.guard[0]; i < grid.dims[0] - grid.guard[0]; ++i ) {
//       int buffer_index = ( i - grid.guard[0] ) + buffer_index_phi;
//       Scalar r_sp = std::exp( grid.pos( 0, i, B1stagger.x ) );
//       // set guard cell on the left which is next to bulk to zero. NOTE this is enough for second order accuracy
//       flux( i, grid.guard[1] - 1, k ) = 0;
//       // scan over the cells in the bulk
//       for ( int j = grid.guard[1]; j < grid.dims[1] - grid.guard[1]; ++j ) {
//         Scalar theta = grid.pos( 1, j , B1stagger.y );
//         flux( i, j , k ) = flux( i, j-1, k ) + grid.delta[1] * r_sp * r_sp * std::sin(theta) * data.Bfield(0, i, j-1, k);
//       }

//       // register the last flux value in the bulk in the sendbuffer
//       sendbuffer[buffer_index] = flux( i, grid.dims[1] - grid.guard[1] - 1, k );
//     }
//   }

//   // MPI_ExScan
//   domain.comm().cartesian().scan( sendbuffer, recvbuffer, grid.reducedDim(0) * grid.reducedDim(2), 1, true );
//   // set the 0th element of recvbuffer to 0 in case it is intialized to some random value
//   recvbuffer[0] = 0.0;

//   // add recvbuffers to local values
//   for ( int k = grid.guard[2]; k < grid.dims[2] - grid.guard[2]; ++k ) {
//     int buffer_index_phi = ( k - grid.guard[2] ) * grid.reducedDim(0);
//     for ( int i = grid.guard[0]; i < grid.dims[0] - grid.guard[0]; ++i ) {
//       int buffer_index = ( i - grid.guard[0] ) + buffer_index_phi;
//       for ( int j = grid.guard[1]; j < grid.dims[1] - grid.guard[1]; ++j ) {
//         flux( i, j, k ) += recvbuffer[buffer_index];
//       }
//     }
//   }

//   // sample the data. NOTE that sampleData performs sendGuardCells in the end.
//   SampleData( flux, data_out.magneticFlux );

//   delete [] sendbuffer;
//   delete [] recvbuffer;

// }

// template < typename P >
// void CountNumbers( const Particles<P>& particles, ScalarField<unsigned int>& ptcNum ) {
//   ptcNum.assign(0);
//   for ( unsigned int i = 0 ; i < particles.Number(); ++i ) {
//     if ( particles.IsEmpty(i) ) continue;
//     ptcNum.data() [ particles.PtcData()[i].cell ] += 1;
//   }

// }


namespace io {

  // FIXME TODO
  // void ExportImpl( DBfile* dbfile, const AperData& data, const AperParams& params, Domain& domain ) {
  //   if ( _comm->is_idle() ) return;

  //   const Grid& grid = params.grid();

  //   auto& contents = InfoCollector::Instance().contents;
  //   contents.t_export_ens = 0.0;
  //   auto t0 = high_resolution_clock::now();

  //   ScalarField<Scalar> _scalar_tmp(grid);
  //   VectorField<Scalar> _vector_tmp(grid);

  //   // Fields export
  //   // NOTE sendguardcells for E, B and J should be done outside.
  //   if ( _comm -> is_primary() ) {
  //     put_var(dbfile, "E1", data.Efield.ptr(0));
  //     put_var(dbfile, "E2", data.Efield.ptr(1));
  //     put_var(dbfile, "E3", data.Efield.ptr(2));

  //     put_var(dbfile, "B1", data.Bfield.ptr(0));
  //     put_var(dbfile, "B2", data.Bfield.ptr(1));
  //     put_var(dbfile, "B3", data.Bfield.ptr(2));

  //     put_var(dbfile, "J1", data.Jfield.ptr(0));
  //     put_var(dbfile, "J2", data.Jfield.ptr(1));
  //     put_var(dbfile, "J3", data.Jfield.ptr(2));

  //     fd.ComputeDivergence( Efield, _scalar_tmp, FieldType::ETYPE, params.ens_specs.is_bdry.data() );
  //     domain.SendGuardCells( _scalar_tmp );
  //     put_var(dbfile, "divE", _scalar_tmp.ptr());

  //     fd.ComputeDivergence( Bfield, _scalar_tmp, FieldType::BTYPE, params.ens_specs.is_bdry.data() );
  //     domain.SendGuardCells( _scalar_tmp );
  //     put_var(dbfile, "divB", _scalar_tmp.ptr());
  //   }

  //   // FIXME TODO save gamma P density, which is total gamma / physical cell volume. Later this divided by number density gives us gamma P per unit particle. TODO: think this over.
  //   // Particles export
  //   // NOTE deposited fields, need ens reduce followed by sendGuardCells
  //   for ( auto& elm : data.particles ) {
  //     auto ptcType = elm.first;
  //     auto& ptcs = elm.second;

  //     auto ptcStr = PtcType2Str(ptcType);

  //     // FIXME TODO save number density instead of just number
  //     { auto& num_field = _scalar_tmp;
  //       num_field.assign(0.0);

  //       ParticleDepositer pd;
  //       std::visit( [&]( auto& x ) {
  //                     pd.Deposit(num_field, x, params, comm, domain,
  //                                [](auto&& a) {return 1.0;});
  //                   } , ptcs );

  //       domain.EnsReduceFields(num_field);
  //       if ( _comm -> is_primary() ) {
  //         DBoptlist* optlist = DBMakeOptlist(1);
  //         auto charge = 1.0; // FIXME TODO write a function to return physical charge for species
  //         DBAddOption(optlist, DBOPT_DTIME, &charge); // used to store charge of the species
  //         domain.SendGuardCells(num_field);
  //         put_var( dbfile, "n_"+ptcStr, num_field.ptr(), optlist );
  //         DBFreeOptlist(optlist);
  //       }
  //     }


  //     { auto& E_field = _scalar_tmp;
  //       E_field.assign(0.0);
  //       ParticleDepositer pd;
  //       std::visit( [&]( auto& x ) {
  //                     pd.Deposit(E_field, x, params, comm, domain,
  //                                [](const auto& ptc) {
  //                                  if constexpr
  //                                    (true)
  //                                      return sqrt(1.0 + ptc.p.dot(ptc.p) );
  //                                  else return ptc.E;
  //                                } ) } , ptcs );

  //       domain.EnsReduceFields(E_field);
  //       if ( _comm -> is_primary() ) {
  //         domain.SendGuardCells(E_field);
  //         put_var( dbfile, "E_" + ptcStr, E_field.ptr() );
  //       }
  //     }


  //     { auto& P = _vector_tmp;
  //       P.assign(0.0);

  //       ParticleDepositer pd;
  //       std::visit( [&]( auto& x ) {
  //                     pd.Deposit(P, x, params, comm, domain,
  //                                [](const auto& ptc) { return ptc.p;});
  //                   } , ptcs );
  //       domain.EnsReduceFields(P);
  //       if ( _comm -> is_primary() ) {
  //         domain.SendGuardCells(P);
  //         for ( int i = 0; i < 3; ++i )
  //           put_var( dbfile, "P"+std::to_string(i+1)+"_"+ptcStr, P.ptr(i), optlist );
  //       }
  //     }

  //     // charged species do the following. Check charge
  //     // FIXME in the presence of this, may we just drop exporting J_total?
  //     if ( is_charged(ptcType) ) {
  //       auto& J_sp = _vector_tmp;
  //       J_sp.assign(0.0);
  //       CurrentDepositer j_dep;
  //       std::visit([&]( auto& x ) {
  //                    j_dep( J_sp, dt, x, params, comm, domain );
  //                  }, ptcs );
  //       if ( _comm -> is_primary() ) {
  //         put_var( dbfile, "J1_"+ptcStr, J_sp.ptr(0) );
  //         put_var( dbfile, "J2_"+ptcStr, J_sp.ptr(1) );
  //         put_var( dbfile, "J3_"+ptcStr, J_sp.ptr(2) );
  //       }

  //     }

  //   }

  //   { auto& pc_rate = _scalar_tmp;
  //     pc_rate.copyFrom(data.pairCreationEvents);
  //     domain.EnsReduceFields(pc_rate);
  //     if ( _comm -> is_primary() ) {
  //       // convert to rate
  //       pc_rate.multiplyBy( 1.0 / ( _pane.interval * params.dt ) );
  //       domain.SendGuardCells( pc_rate );
  //       put_var( dbfile, "Pair_Creation_Rate", pc_rate.ptr() );
  //     }
  //   }

  //   auto t1 = high_resolution_clock::now();
  //   auto dur = duration_cast<clock_cycle>(t1 - t0);
  //   contents.t_export = dur.count();
  // }


  void export_data( int timestep, const DynamicVars& dvars, const Params& params, const std::optional<mpi::Comm>& primary, const mpi::Comm& ensemble ) {
    constexpr int Padding = 6;
    char str_ts [Padding + 1];
    sprintf(str_ts, "%06d", timestep); // TODOL use Padding instead of 6

    std::string dirname = filesys::append_slash(params.this_run_dir) + "data/timestep" + str_ts + "/";
    filesys::create_directories(dirname);

    if ( !primary ) {
      // ExportImpl(dbfile, dvars, params, domain); //TODO
      return;
    }

    auto time = timestep * params.dt;
    // DBoptlist* optlist = DBMakeOptlist(8);
    // DBAddOption(optlist, DBOPT_TIME, &time);
    // DBAddOption(optlist, DBOPT_CYCLE, &timeStep);
    // if (displayGuard) {
    //   DBAddOption(optlist, DBOPT_LO_OFFSET, lowOffset.data());
    //   DBAddOption(optlist, DBOPT_HI_OFFSET, hiOffset.data());
    // }
    // DBAddOption(optlist, DBOPT_BASEINDEX, pos.data());

    // DBSetFriendlyHDF5Names(1);
    // if (_pane.useCompression) DBSetCompression("METHOD=GZIP");

    auto dbfile = silo::pmpio::open<silo::Mode::Write>( dirname, *primary );

    // put_mesh(dbfile, grid, optlist);
    // ExportImpl(dbfile, dvars, params, domain);

    // DBFreeOptlist(optlist);
    silo::close( dbfile );

    if ( primary->rank() != 0 ) return;

    // generate master file
    std::string masterfile = filesys::append_slash(params.this_run_dir) + "data/timestep" + str_ts + ".silo";

    dbfile = silo::open<silo::Mode::Write>( masterfile ); // TODOL double check if recycling dbfile is OK

    // put_multimesh ( dbfielMaster, timestep );

    // for (unsigned int i = 0; i < quadVars.size(); i++) {
    //   put_multivar( dbfileMaster, timestep, varname, vartype );
    // }

    silo::close( dbfile );

  }
}

