#include "silopp/silo_toc.hpp"

#include <silo.h>

namespace silo {
template <typename file_t>
std::vector<std::string> Toc<file_t>::toc_impl(int num, char** names) {
  std::vector<std::string> res;
  res.reserve(num);
  for (int i = 0; i < num; ++i) res.emplace_back(names[i]);
  return res;
}

template <typename file_t>
std::vector<std::string> Toc<file_t>::toc_dir() {
  DBtoc* toc = DBGetToc(_dbfile());
  return toc_impl(toc->ndir, toc->dir_names);
}

// template < typename file_t >
// std::vector<std::string> Toc<file_t>::toc_array() {
//   DBtoc* toc = DBGetToc( _dbfile() );
//   return toc_impl(toc->narray, toc->array_names);
// }

template <typename file_t>
std::vector<std::string> Toc<file_t>::toc_var() {
  DBtoc* toc = DBGetToc(_dbfile());
  return toc_impl(toc->nvar, toc->var_names);
}

template <typename file_t>
std::vector<std::string> Toc<file_t>::toc_multimesh() {
  DBtoc* toc = DBGetToc(_dbfile());
  return toc_impl(toc->nmultimesh, toc->multimesh_names);
}

template <typename file_t>
std::vector<std::string> Toc<file_t>::toc_multivar() {
  DBtoc* toc = DBGetToc(_dbfile());
  return toc_impl(toc->nmultivar, toc->multivar_names);
}

}  // namespace silo

#include "silopp/silo++.hpp"
namespace silo {
template struct Toc<file_t>;
}
