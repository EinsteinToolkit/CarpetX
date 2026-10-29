#ifndef CARPETX_CARPETX_IO_CONDUIT_HXX
#define CARPETX_CARPETX_IO_CONDUIT_HXX

#include <cctk.h>

#ifdef HAVE_CAPABILITY_Conduit

#include <string>
#include <vector>

namespace CarpetX {

int InputConduitParameters(const std::string &input_dir,
                           const std::string &input_file);
void InputConduitGridStructure(cGH *cctkGH, const std::string &input_dir,
                               const std::string &input_file,
                               int input_iteration);
void InputConduit(const cGH *cctkGH, const std::vector<bool> &input_group,
                  const std::string &input_dir, const std::string &input_file);

void OutputConduit(const cGH *cctkGH, const std::vector<bool> &output_group,
                   const std::string &output_dir,
                   const std::string &output_file);

} // namespace CarpetX

#endif

#endif // #ifndef CARPETX_CARPETX_IO_CONDUIT_HXX
