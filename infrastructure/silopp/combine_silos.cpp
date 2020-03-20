#include "filesys/filesys.hpp"
#include <iostream>
#include <unistd.h>
#include <hdf5.h>
#include <memory>
#include <unordered_set>

int main(int argc, char** argv) {
  if ( argc != 2 ) {
    std::cout << "Incorrect number of arguments" << std::endl;
    return 0;
  }

  chdir(argv[1]);

  constexpr int Ts_start = -11;
  constexpr int Ts_len = 6;

  std::unordered_set<int> timesteps_treated;

  for ( const auto dir : fs::directory_iterator(".") ) {
    if (! fs::is_directory(dir) ) continue;
    for ( const auto fl : fs::directory_iterator(dir) ) {
      if ( fl.find(".silo") == std::string::npos ) continue;
      auto ts_str = fl.substr( fl.size() + Ts_start, Ts_len );
      int ts = std::stoi(ts_str);
      if ( timesteps_treated.find(ts) != timesteps_treated.end() ) {
        // if found ts, skip
        continue;
      } else {
        timesteps_treated.insert(ts);
        std::string target_file = "timestep"+ts_str+".silo";
        if ( fs::exists(target_file) ) fs::remove(target_file);
        fs::copy_file( fl, target_file );

        auto file_id = H5Fopen(target_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        auto group_id = H5Gopen(file_id, "/.silo", H5P_DEFAULT);

        hsize_t  num_obj;
        H5Gget_num_objs(group_id, &num_obj);
        for ( int i = 1; i <= num_obj; i += 2 ) { // only odd i need be modified
          char name[20];
          snprintf(name,10,"#%06d",i); // null terminator is automatically appended

          std::string file_ns;
          hid_t datatype;
          { // retrieve stored file_ns
            auto did = H5Dopen(group_id, name, H5P_DEFAULT);
            std::unique_ptr<char[]> tmp;
            auto size = H5Dget_storage_size(did);
            datatype = H5Dget_type(did);
            tmp.reset(new char[size]);
            H5Dread(did, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void*)(tmp.get())); // it is null terminated

            H5Dclose(did);
            file_ns = {tmp.get()};
          }
          H5Ldelete(group_id, name, H5P_DEFAULT); // deletes the dataset

          { // first create a new dataset
            std::string file_ns_new = file_ns[0] + dir.substr(2) + "/" + file_ns.substr(1);
            hsize_t len = file_ns_new.size() + 1; // +1 for null terminator
            hid_t dataspace = H5Screate_simple(1, &len, NULL); // second arg is the current dimensions

            hid_t did = H5Dcreate(group_id, name, datatype, dataspace,
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

            // then write the new file_ns
            H5Dwrite(did,H5T_STD_U8LE,H5S_ALL,H5S_ALL,H5P_DEFAULT,(void*)(file_ns_new.c_str()));

            H5Dclose(did);
            H5Sclose(dataspace);
          }
          H5Tclose(datatype);

        }
        H5Gclose(group_id);
        H5Fclose(file_id);


      }
    }
  }

  chdir("..");
  return 0;
}
