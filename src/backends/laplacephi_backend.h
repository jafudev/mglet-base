#pragma once

#include <cstddef>

#include "f_arr_view.h"
#include "mglet_precision.h"

namespace mglet::backend
{

void laplacephi_backend(
    FArrView<mgletreal> res,
    FArrView<const mgletreal> phi,
    FArrView<const mgletreal> gsaw,
    FArrView<const mgletreal> gsae,
    FArrView<const mgletreal> gsas,
    FArrView<const mgletreal> gsan,
    FArrView<const mgletreal> gsab,
    FArrView<const mgletreal> gsat,
    FArrView<const mgletreal> gsap,
    FArrView<const mgletreal> bp,
    FArrView<const mgletint> mygrids,
    FArrView<const mgletint> kkk,
    FArrView<const mgletint> jjj,
    FArrView<const mgletint> iii,
    FArrView<const mgletint> ip3d,
    FArrView<const mgletint> ip1dx,
    FArrView<const mgletint> ip1dy,
    FArrView<const mgletint> ip1dz);

void laplacephi_level_backend(
    FArrView<mgletreal> res,
    FArrView<const mgletreal> phi,
    FArrView<const mgletreal> gsaw,
    FArrView<const mgletreal> gsae,
    FArrView<const mgletreal> gsas,
    FArrView<const mgletreal> gsan,
    FArrView<const mgletreal> gsab,
    FArrView<const mgletreal> gsat,
    FArrView<const mgletreal> gsap,
    FArrView<const mgletreal> bp,
    FArrView<const mgletint> mygridsonlvl,
    FArrView<const mgletint> kkk,
    FArrView<const mgletint> jjj,
    FArrView<const mgletint> iii,
    FArrView<const mgletint> ip3d,
    FArrView<const mgletint> ip1dx,
    FArrView<const mgletint> ip1dy,
    FArrView<const mgletint> ip1dz);

} // namespace mglet::backend
