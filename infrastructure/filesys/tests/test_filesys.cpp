#include "testfw/testfw.hpp"
#include "filesys/filesys.hpp"
#include <filesystem>
#include <fstream>

using namespace fs;

SCENARIO("Test directory_iterator", "[filesys]") {
  const std::string dir = "TestFS";
  {
    create_directories(dir);
    std::ofstream of;
    of.open(dir + "/abc.txt");
    of.close();
    of.open(dir + "/xyz.silo");
    of.close();
  }
  {
    auto std_dir_itr = std::filesystem::directory_iterator(dir);
    auto std_dir_itr_end = std::filesystem::end(std_dir_itr);
    auto my_dir_itr = directory_iterator(dir);
    auto my_dir_itr_end = my_dir_itr.end();

    while ( std_dir_itr != std_dir_itr_end ) {
      std::string my_entry = *my_dir_itr;
      std::string std_entry = (*std_dir_itr).path();
      REQUIRE(my_entry == std_entry);
      ++std_dir_itr;
      ++my_dir_itr;
    }
    REQUIRE_FALSE( my_dir_itr != my_dir_itr_end );
  }

  remove_all(dir);
}
