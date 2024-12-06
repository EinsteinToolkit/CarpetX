#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "loop.hxx"

using Loop::dim;

extern "C" void TestArrayGroup_DynamicData(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS_TestArrayGroup_DynamicData;

  struct errorcount_t {
    int lsh, ash, gsh, lbnd, ubnd, bbox, nghostzones;
  };

  // Validate grid function dynamic data
  {
    const int gf_gi = CCTK_GroupIndex("TestArrayGroup::test_gf");
    cGroupDynamicData gf_data;
    CCTK_GroupDynamicData(cctkGH, gf_gi, &gf_data);

    cGroup gf_group;
    const int ierr = CCTK_GroupData(gf_gi, &gf_group);
    if (ierr)
      CCTK_ERROR("error in GroupData for grid functions");

    if (gf_data.dim != gf_group.dim || gf_data.dim != 3)
      CCTK_ERROR("incorrect dimension in grid function dynamic data");
    if (gf_data.activetimelevels != CCTK_ActiveTimeLevelsGI(cctkGH, gf_gi))
      CCTK_ERROR("incorrect activetimelevels in grid function dynamic data");

    errorcount_t error{0, 0, 0, 0, 0, 0, 0};

    for (int i = 0; i < gf_data.dim; i++) {
      if (gf_data.lsh[i] != cctkGH->cctk_lsh[i])
        error.lsh += 1;
      if (gf_data.ash[i] != cctkGH->cctk_ash[i])
        error.ash += 1;
      if (gf_data.gsh[i] != cctkGH->cctk_gsh[i])
        error.gsh += 1;
      if (gf_data.lbnd[i] != cctkGH->cctk_lbnd[i])
        error.lbnd += 1;
      if (gf_data.ubnd[i] != cctkGH->cctk_ubnd[i])
        error.lbnd += 1;
      if (gf_data.bbox[2 * i] != cctkGH->cctk_bbox[2 * i])
        error.bbox += 1;
      if (gf_data.bbox[2 * i + 1] != cctkGH->cctk_bbox[2 * i + 1])
        error.bbox += 1;
      if (gf_data.nghostzones[i] != cctkGH->cctk_nghostzones[i])
        error.nghostzones += 1;
    }
    if (error.lsh)
      CCTK_ERROR("incorrect lsh data in grid function dynamic data");
    if (error.ash)
      CCTK_ERROR("incorrect ash data in grid function dynamic data");
    if (error.gsh)
      CCTK_ERROR("incorrect gsh data in grid function dynamic data");
    if (error.lbnd)
      CCTK_ERROR("incorrect lbnd data in grid function dynamic data");
    if (error.ubnd)
      CCTK_ERROR("incorrect ubnd data in grid function dynamic data");
    if (error.bbox)
      CCTK_ERROR("incorrect bbox data in grid function dynamic data");
    if (error.nghostzones)
      CCTK_ERROR("incorrect nghostzones data in grid function dynamic data");
  }

  // Validate grid scalar dynamic data
  char const *const scalars[] = {
      "TestArrayGroup::test_scalar",
      "TestArrayGroup::test_scalar_int",
      "TestArrayGroup::test_scalar_complex",
  };
  for (auto scalar : scalars) {
    const int scalar_gi = CCTK_GroupIndex(scalar);
    cGroupDynamicData scalar_data;
    CCTK_GroupDynamicData(cctkGH, scalar_gi, &scalar_data);

    cGroup scalar_group;
    const int ierr = CCTK_GroupData(scalar_gi, &scalar_group);
    if (ierr)
      CCTK_ERROR("error in GroupData for scalars");

    if (scalar_data.dim != scalar_group.dim || scalar_data.dim != 0)
      CCTK_ERROR("incorrect dimension in grid scalar dynamic data");
    if (scalar_data.activetimelevels !=
        CCTK_ActiveTimeLevelsGI(cctkGH, scalar_gi))
      CCTK_ERROR("incorrect activetimelevels in grid scalar dynamic data");

    errorcount_t error{0, 0, 0, 0, 0, 0, 0};

    for (int i = 0; i < scalar_data.dim; i++) {
      if (scalar_data.lsh[i] != -1)
        error.lsh += 1;
      if (scalar_data.ash[i] != -1)
        error.ash += 1;
      if (scalar_data.gsh[i] != -1)
        error.gsh += 1;
      if (scalar_data.lbnd[i] != -1)
        error.lbnd += 1;
      if (scalar_data.ubnd[i] != -1)
        error.lbnd += 1;
      if (scalar_data.bbox[2 * i] != -1)
        error.bbox += 1;
      if (scalar_data.bbox[2 * i + 1] != -1)
        error.bbox += 1;
      if (scalar_data.nghostzones[i] != -1)
        error.nghostzones += 1;
    }
    // data is padded to Loop:dim with "neutral"
    for (int i = scalar_data.dim; i < dim; i++) {
      if (scalar_data.lsh[i] != 1)
        error.lsh += 1;
      if (scalar_data.ash[i] != 1)
        error.ash += 1;
      if (scalar_data.gsh[i] != 1)
        error.gsh += 1;
      if (scalar_data.lbnd[i] != 0)
        error.lbnd += 1;
      if (scalar_data.ubnd[i] != 0)
        error.lbnd += 1;
      if (scalar_data.bbox[2 * i] != 1)
        error.bbox += 1;
      if (scalar_data.bbox[2 * i + 1] != 1)
        error.bbox += 1;
      if (scalar_data.nghostzones[i] != 0)
        error.nghostzones += 1;
    }
    if (error.lsh)
      CCTK_ERROR("incorrect lsh data in scalar dynamic data");
    if (error.ash)
      CCTK_ERROR("incorrect ash data in scalar dynamic data");
    if (error.gsh)
      CCTK_ERROR("incorrect gsh data in scalar dynamic data");
    if (error.lbnd)
      CCTK_ERROR("incorrect lbnd data in scalar dynamic data");
    if (error.ubnd)
      CCTK_ERROR("incorrect ubnd data in scalar dynamic data");
    if (error.bbox)
      CCTK_ERROR("incorrect bbox data in scalar dynamic data");
    if (error.nghostzones)
      CCTK_ERROR("incorrect nghostzones data in scalar dynamic data");
  }

  // Validate grid array dynamic data
  char const *const arrays[] = {
      "TestArrayGroup::test_array",
      "TestArrayGroup::test_array_int",
      "TestArrayGroup::test_array_complex",
  };
  for (auto array : arrays) {
    const int array_gi = CCTK_GroupIndex(array);
    cGroupDynamicData array_data;
    CCTK_GroupDynamicData(cctkGH, array_gi, &array_data);

    cGroup array_group;
    const int ierr = CCTK_GroupData(array_gi, &array_group);
    if (ierr)
      CCTK_ERROR("error in GroupData for arrays");

    const int sz[2] = {5, 6};

    if (array_data.dim != array_group.dim || array_data.dim != 2)
      CCTK_ERROR("incorrect dimension in array dynamic data");
    if (array_data.activetimelevels != 1)
      CCTK_ERROR("incorrect activetimelevels in array dynamic data");

    errorcount_t error{0, 0, 0, 0, 0, 0, 0};

    for (int i = 0; i < array_data.dim; i++) {
      if (array_data.lsh[i] != sz[i])
        error.lsh += 1;
      if (array_data.ash[i] != sz[i])
        error.ash += 1;
      if (array_data.gsh[i] != sz[i])
        error.gsh += 1;
      if (array_data.lbnd[i] != 0)
        error.lbnd += 1;
      if (array_data.ubnd[i] != sz[i] - 1)
        error.lbnd += 1;
      if (array_data.bbox[2 * i] != 1)
        error.bbox += 1;
      if (array_data.bbox[2 * i + 1] != 1)
        error.bbox += 1;
      if (array_data.nghostzones[i] != 0)
        error.nghostzones += 1;
    }
    // data is padded to Loop:dim with "neutral"
    for (int i = array_data.dim; i < dim; i++) {
      if (array_data.lsh[i] != 1)
        error.lsh += 1;
      if (array_data.ash[i] != 1)
        error.ash += 1;
      if (array_data.gsh[i] != 1)
        error.gsh += 1;
      if (array_data.lbnd[i] != 0)
        error.lbnd += 1;
      if (array_data.ubnd[i] != 0)
        error.lbnd += 1;
      if (array_data.bbox[2 * i] != 1)
        error.bbox += 1;
      if (array_data.bbox[2 * i + 1] != 1)
        error.bbox += 1;
      if (array_data.nghostzones[i] != 0)
        error.nghostzones += 1;
    }
    if (error.lsh)
      CCTK_ERROR("incorrect lsh data in array dynamic data");
    if (error.ash)
      CCTK_ERROR("incorrect ash data in array dynamic data");
    if (error.gsh)
      CCTK_ERROR("incorrect gsh data in array dynamic data");
    if (error.lbnd)
      CCTK_ERROR("incorrect lbnd data in array dynamic data");
    if (error.ubnd)
      CCTK_ERROR("incorrect ubnd data in array dynamic data");
    if (error.bbox)
      CCTK_ERROR("incorrect bbox data in array dynamic data");
    if (error.nghostzones)
      CCTK_ERROR("incorrect nghostzones data in array dynamic data");
  }
}
