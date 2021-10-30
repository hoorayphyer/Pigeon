#include "filesys.hpp"

#include <filesystem>
#include <iostream>

namespace filesystem = std::filesystem;
using std::error_code;

namespace fs {
std::string absolute(std::string dir) {
  filesystem::path p = dir;
  return filesystem::absolute(p);
}

// std::string canonical( std::string dir ) {
//   filesystem::path p = dir;
//   return filesystem::canonical(p);
// }

std::string& append_slash(std::string& dir) {
  if (dir != "" && dir.back() != '/') dir.append("/");
  return dir;
}

std::string& remove_slash(std::string& dir) {
  if (dir.back() == '/') dir.pop_back();
  return dir;
}

void create_directories(std::string dir) {
  // NOTE on some platforms path ending with a slash is not created as a
  // directory
  remove_slash(dir);
  filesystem::path p(dir);
  error_code returnedError;
  if (!filesystem::create_directories(p, returnedError) &&
      !filesystem::exists(p)) {
    std::cout << "Directory neither existed nor created: " << dir << std::endl;
  }
  if (returnedError) {
    std::cerr << "ERROR creating directories " << dir
              << "\nError code = " << returnedError.value() << std::endl;
  }
}

bool exists(std::string dir) {
  filesystem::path p(dir);
  return filesystem::exists(p);
}

void remove(std::string file) {
  filesystem::path p(file);
  filesystem::remove(p);
}

void remove_all(std::string dir) {
  filesystem::path p(dir);
  filesystem::remove_all(p);
}

void create_directory_symlink(std::string to, std::string new_link) {
  // NOTE to and new_link must not end with '/'
  remove_slash(to);
  remove_slash(new_link);

  filesystem::path p_to(to);
  filesystem::path p_new_link(new_link);
  filesystem::create_directory_symlink(filesystem::absolute(p_to),
                                       filesystem::absolute(p_new_link));
}

void copy_file(std::string target, std::string dest) {
  filesystem::path target_path(target);
  filesystem::path dest_path(dest);

  error_code returnedError;
  filesystem::copy_file(target_path, dest_path, returnedError);
  if (returnedError) {
    std::cerr << "ERROR copying " << target << " to " << dest
              << "\nError code = " << returnedError.value() << std::endl;
  }
}

void rename(std::string old_path, std::string new_path) {
  filesystem::rename({old_path}, {new_path});
}

bool equivalent(std::string p1, std::string p2) {
  return filesystem::equivalent({p1}, {p2});
}
}  // namespace fs

namespace fs {
bool is_directory(std::string dir) { return filesystem::is_directory({dir}); }

bool is_regular_file(std::string dir) {
  return filesystem::is_regular_file({dir});
}
}  // namespace fs

namespace fs {
using std_dir_itr = filesystem::directory_iterator;

directory_iterator::directory_iterator() : _itr(std_dir_itr()) {}

directory_iterator::directory_iterator(std::string dir)
    : _itr(std_dir_itr(filesystem::path(dir))) {}

bool directory_iterator::operator!=(const directory_iterator& other) const {
  return std::any_cast<const std_dir_itr&>(_itr) !=
         std::any_cast<const std_dir_itr&>(other._itr);
}

directory_iterator::value_type directory_iterator::operator*() const {
  return (*(std::any_cast<const std_dir_itr&>(_itr))).path();
}

directory_iterator& directory_iterator::operator++() {
  ++(std::any_cast<std_dir_itr&>(_itr));
  return *this;
}

directory_iterator directory_iterator::end() const noexcept { return {}; }
}  // namespace fs
