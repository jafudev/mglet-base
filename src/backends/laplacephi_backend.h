#pragma once

#include <cstddef>

#include "f_arr_view.h"
#include "mglet_precision.h"

namespace mglet::backend
{

void laplacephi(
    FArrView<mgletreal> res,
    FArrView<mgletreal> phi,
    FArrView<mgletreal> gsaw,
    FArrView<mgletreal> gsae,
    FArrView<mgletreal> gsas,
    FArrView<mgletreal> gsan,
    FArrView<mgletreal> gsab,
    FArrView<mgletreal> gsat,
    FArrView<mgletreal> gsap,
    FArrView<mgletreal> bp,
    FArrView<mgletint> mygrids,
    FArrView<mgletint> kkk,
    FArrView<mgletint> jjj,
    FArrView<mgletint> iii,
    FArrView<mgletint> ip3d,
    FArrView<mgletint> ip1dx,
    FArrView<mgletint> ip1dy,
    FArrView<mgletint> ip1dz);

} // namespace mglet::backend
