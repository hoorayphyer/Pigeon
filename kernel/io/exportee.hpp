#ifndef _IO_EXPORTEE_HPP_
#define _IO_EXPORTEE_HPP_

namespace io {
  template < int DGrid,
             typename Real,
             template < typename > class S,
             typename ShapeF,
             typename RealJ,
             typename Metric >
  struct Exportee {
    void (*action)( const mani::Grid<Real,DGrid>& grid; // local grid
                    field::Field<Real, 3, DGrid>& E;
                    field::Field<Real, 3, DGrid>& B;
                    field::Field<RealJ, 3, DGrid>& J;// J is Jmesh on a replica
                    particle::map<particle::array<Real,S>>& particles;
                   ) = nullptr;

    template < typename Texp >
    static void save_to( dbfile*);
  };

}

#endif
