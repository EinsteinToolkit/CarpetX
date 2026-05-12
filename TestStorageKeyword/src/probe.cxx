#include <cctk.h>
#include <cctk_Arguments.h>

extern "C" void TestStorageKeyword_Probe(CCTK_ARGUMENTS) {
  struct expectation {
    const char *name;
    int active_tl;
  };
  const expectation specs[] = {
      {"TestStorageKeyword::storaged_gf", 2},
      {"TestStorageKeyword::unstoraged_gf", 0},
      {"TestStorageKeyword::storaged_arr", 2},
      {"TestStorageKeyword::unstoraged_arr", 0},
  };
  for (const auto &spec : specs) {
    const int gi = CCTK_GroupIndex(spec.name);
    if (gi < 0)
      CCTK_VERROR("Group '%s' not found", spec.name);
    const int has = CCTK_QueryGroupStorage(cctkGH, spec.name);
    const int got = CCTK_ActiveTimeLevelsGI(cctkGH, gi);
    const int want_has = spec.active_tl > 0 ? 1 : 0;
    if (has != want_has)
      CCTK_VERROR("Group '%s': CCTK_QueryGroupStorage=%d, expected %d",
                  spec.name, has, want_has);
    if (got != spec.active_tl)
      CCTK_VERROR("Group '%s': CCTK_ActiveTimeLevelsGI=%d, expected %d",
                  spec.name, got, spec.active_tl);
  }
}
