#ifndef CARPETX_CARPETX_IO_SILO_HXX
#define CARPETX_CARPETX_IO_SILO_HXX

#include <cctk.h>

#ifdef HAVE_CAPABILITY_Silo

#include <string>
#include <vector>

namespace CarpetX {

int InputSiloParameters(const std::string &input_dir,
                        const std::string &input_file);
void InputSiloGridStructure(cGH *cctkGH, const std::string &input_dir,
                            const std::string &input_file, int input_iteration);
void InputSilo(const cGH *cctkGH, const std::vector<bool> &input_group,
               const std::string &input_dir, const std::string &input_file);

void OutputSilo(const cGH *cctkGH, const std::vector<bool> &output_group,
                const std::string &output_dir, const std::string &output_file);

} // namespace CarpetX

#endif

#endif // #ifndef CARPETX_CARPETX_IO_SILO_HXX
