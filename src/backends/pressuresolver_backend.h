#pragma once

#include <cstddef>

#include "f_arr_view.h"
#include "mglet_precision.h"

namespace mglet::backend
{

void maxabscal_backend(
    const FArrView<mgletreal> maxabsgrid,
    FArrView<mgletreal> phi,
    const FArrView<mgletint> mygrids,
    mgletint nmygrids,
    const FArrView<mgletint> kkk,
    const FArrView<mgletint> jjj,
    const FArrView<mgletint> iii,
    const FArrView<mgletint> ip3d);

void accumulate_pcorr_backend(
    FArrView<mgletreal> dp_view,
    const FArrView<mgletreal> hilf_view);

void rescal_backend(
    FArrView<mgletreal> rhs_view,
    const FArrView<mgletreal> res_view,
    mgletint nmygrids,
    const FArrView<mgletint> mygrids,
    const FArrView<mgletint> kkk,
    const FArrView<mgletint> jjj,
    const FArrView<mgletint> iii,
    const FArrView<mgletint> ip3d);

void sipiter1_hyperplane_level_backend(
    FArrView<mgletreal> res,
    const FArrView<mgletreal> rhs,
    const FArrView<mgletreal> siplw,
    const FArrView<mgletreal> sipls,
    const FArrView<mgletreal> siplb,
    const FArrView<mgletreal> siplpr,
    const FArrView<mgletifk> miphp,
    const FArrView<mgletifk> idxhp,
    mgletint nmygridsonlvl,
    const FArrView<mgletint> mygridsonlvl,
    const FArrView<mgletint> kkk,
    const FArrView<mgletint> jjj,
    const FArrView<mgletint> iii,
    const FArrView<mgletint> ip3d);

void sipiter2_hyperplane_level_backend(
    FArrView<mgletreal> dp,
    FArrView<mgletreal> res,
    const FArrView<mgletreal> sipue,
    const FArrView<mgletreal> sipun,
    const FArrView<mgletreal> siput,
    const FArrView<mgletifk> miphp,
    const FArrView<mgletifk> idxhp,
    mgletint nmygridsonlvl,
    const FArrView<mgletint> mygridsonlvl,
    const FArrView<mgletint> kkk,
    const FArrView<mgletint> jjj,
    const FArrView<mgletint> iii,
    const FArrView<mgletint> ip3d);

} // namespace mglet::backend
