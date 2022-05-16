
#ifdef HAS_HDF5
#include <hdf5.h>
#include <ctime>

#include "io.h"

void save_simulation(const Bbh_system& bhs, const Stellar_system& stars, std::string filename, std::string name, std::string description) {
  hid_t file_id, group_id, dataspace_id, attribute_id, dataset_id, strtype;
  herr_t status;
  hsize_t dims[2];
  time_t rawtime;

  if (filename.empty()) filename = "Untitled.h5";
  if (name.empty()) name = "Untitled";
  if (description.empty()) description = "No description";

  file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  dims[0] = 1;
  dataspace_id = H5Screate_simple(1, dims, NULL);
  attribute_id = H5Acreate(file_id, "Time", H5T_STD_U64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT); // H5T_UNIX_D64LE
  time(&rawtime);
  H5Awrite(attribute_id, H5T_NATIVE_ULONG, &rawtime);
  H5Aclose(attribute_id);

  strtype = H5Tcopy(H5T_C_S1);
  H5Tset_size(strtype, 30);
  attribute_id = H5Acreate(file_id, "Name", strtype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute_id, strtype, name.c_str());
  H5Aclose(attribute_id);
  H5Tset_size(strtype, 512);
  attribute_id = H5Acreate(file_id, "Description", strtype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute_id, strtype, description.c_str());
  H5Aclose(attribute_id);
  H5Sclose(dataspace_id);

  group_id = H5Gcreate(file_id, "/Stars", H5P_DEFAULT, H5P_DEFAULT,
                       H5P_DEFAULT); // Create the "Stars" group.
  /* Attributes: m, sigma */
  dims[0] = 1;
  dataspace_id = H5Screate_simple(
      1, dims, NULL); // Dataspace for a single scalar attribute
  attribute_id = H5Acreate(group_id, "m", H5T_IEEE_F64LE, dataspace_id,
                           H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &stars.m);
  H5Aclose(attribute_id);
  H5Sclose(dataspace_id);

  /* Stars' coordinates */
  dims[0] = stars.coords.size();
  dims[1] = 6;
  dataspace_id = H5Screate_simple(2, dims, NULL);
  dataset_id = H5Dcreate(group_id, "Data", H5T_IEEE_F64LE, dataspace_id,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (double*)stars.coords.data());
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
  /* Status */
  dataspace_id = H5Screate_simple(1, dims, NULL);
  dataset_id = H5Dcreate(group_id, "Status", H5T_STD_U8LE, dataspace_id,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, stars.stat.data());
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
  /* E & L */
  dataspace_id = H5Screate_simple(1, dims, NULL);
  dataset_id = H5Dcreate(group_id, "E", H5T_IEEE_F64LE, dataspace_id,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, stars.E.data());
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
  dataspace_id = H5Screate_simple(1, dims, NULL);
  dataset_id = H5Dcreate(group_id, "L", H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, stars.L.data());
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);

  H5Gclose(group_id); // Close "Stars" group.

  group_id = H5Gcreate(file_id, "/BHs", H5P_DEFAULT, H5P_DEFAULT,
                       H5P_DEFAULT); // Create the "BHs" group.
  /* Attributes: q */
  dims[0] = 1;
  dataspace_id = H5Screate_simple(
      1, dims, NULL); // Dataspace for a single scalar attribute
  attribute_id = H5Acreate(group_id, "q", H5T_IEEE_F64LE, dataspace_id,
                           H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &bhs.q);
  H5Aclose(attribute_id);
  H5Sclose(dataspace_id);

  /* BHs' coordinates */
  dims[0] = 8;
  dataspace_id = H5Screate_simple(1, dims, NULL);
  dataset_id = H5Dcreate(group_id, "Data", H5T_IEEE_F64LE, dataspace_id,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, bhs.coords.data);
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
  /* Forces */
  dims[0] = 4;
  dataspace_id = H5Screate_simple(1, dims, NULL);
  dataset_id = H5Dcreate(group_id, "Forces", H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  double friction[4];
  friction[0] = bhs.f1;
  friction[1] = bhs.ff1;
  friction[2] = bhs.f2;
  friction[3] = bhs.ff2;
  H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, friction);
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
  H5Gclose(group_id); // Close "Stars" group.

  status = H5Fclose(file_id);
}
#endif