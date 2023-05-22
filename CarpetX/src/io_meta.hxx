#ifndef CARPETX_CARPETX_IO_META_HXX
#define CARPETX_CARPETX_IO_META_HXX

#include <cctk.h>

#include <yaml-cpp/yaml.h>

#include <array>
#include <string>
#include <vector>

namespace CarpetX {

enum class group_type_t {
  none,
  gf,
  array,
};

enum class reduction_t {
  volume,
  sum,
  sum_abs,
  sum_squared,
  average,
  standard_deviation,
  minimum,
  maximum,
  maximum_abs,
  norm1,
  norm2,
  norm_inf,
  minimum_location,
  maximum_location,
};

struct output_file_description_t {
  std::string filename;                //
  std::string description;             // human readable
  std::string writer_thorn;            //
  std::vector<std::string> variables;  //
  std::vector<int> iterations;         //
  std::vector<int> output_directions;  //
  std::vector<reduction_t> reductions; //
  std::string format_name;             //
  std::array<int, 3> format_version;   //
};

void OutputMeta_RegisterOutputFile(
    output_file_description_t output_file_description);

void OutputMeta(const cGH *);

} // namespace CarpetX

#endif // #ifndef CARPETX_CARPETX_IO_META_HXX
