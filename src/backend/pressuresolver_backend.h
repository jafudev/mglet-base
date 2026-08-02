#pragma once

#include <cstddef>

#include "f_arr_view.h"
#include "mglet_precision.h"

namespace mglet::backend
{

void maxabscal_backend(
    FArrView<mgletreal> maxabsgrid,
    FArrView<const mgletreal> phi,
    FArrView<const mgletint> mygrids,
    FArrView<const mgletint> kkk,
    FArrView<const mgletint> jjj,
    FArrView<const mgletint> iii,
    FArrView<const mgletint> ip3d);

void rescal_backend(
    FArrView<mgletreal> rhs_view,
    FArrView<const mgletreal> res_view,
    FArrView<const mgletint> mygrids,
    FArrView<const mgletint> kkk,
    FArrView<const mgletint> jjj,
    FArrView<const mgletint> iii,
    FArrView<const mgletint> ip3d);

void sipiter1_hyperplane_level_backend(
    FArrView<mgletreal> res,
    FArrView<const mgletreal> rhs,
    FArrView<const mgletreal> siplw,
    FArrView<const mgletreal> sipls,
    FArrView<const mgletreal> siplb,
    FArrView<const mgletreal> siplpr,
    FArrView<const mgletifk> miphp,
    FArrView<const mgletifk> idxhp,
    FArrView<const mgletint> mygridsonlvl,
    FArrView<const mgletint> kkk,
    FArrView<const mgletint> jjj,
    FArrView<const mgletint> iii,
    FArrView<const mgletint> ip3d);

void sipiter2_hyperplane_level_backend(
    FArrView<mgletreal> dp,
    FArrView<mgletreal> res,
    FArrView<const mgletreal> sipue,
    FArrView<const mgletreal> sipun,
    FArrView<const mgletreal> siput,
    FArrView<const mgletifk> miphp,
    FArrView<const mgletifk> idxhp,
    FArrView<const mgletint> mygridsonlvl,
    FArrView<const mgletint> kkk,
    FArrView<const mgletint> jjj,
    FArrView<const mgletint> iii,
    FArrView<const mgletint> ip3d);

} // namespace mglet::backend
