#ifndef  _PARAMETERS_H_
#define  _PARAMETERS_H_

// TODOL traits and params are basically the same thing. Just compile time known or not
template < typename Real >
struct Params {
  Real dt;
  int total_timesteps;

  Real e; // electric charge
};

#endif
