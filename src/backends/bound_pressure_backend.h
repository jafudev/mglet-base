#pragma once

#include "f_arr_view.h"
#include "mglet_precision.h"

namespace mglet::backend
{

void bound_pressure_bp_backend(
    FArrView<mgletreal> p,
    FArrView<mgletreal> pbuffer,
    const FArrView<mgletreal> bp,
    const FArrView<mgletreal> dx,
    const FArrView<mgletreal> dy,
    const FArrView<mgletreal> dz,
    const FArrView<mgletreal> ddx,
    const FArrView<mgletreal> ddy,
    const FArrView<mgletreal> ddz,
    mgletint nboundtasks,
    const FArrView<mgletint> boundtasks_lvl,
    const FArrView<mgletint> kkk,
    const FArrView<mgletint> jjj,
    const FArrView<mgletint> iii,
    const FArrView<mgletint> ip3d,
    const FArrView<mgletint> ip1dx,
    const FArrView<mgletint> ip1dy,
    const FArrView<mgletint> ip1dz,
    const FArrView<mgletint> ipbb);

void bound_pressure_nobp_backend(
    FArrView<mgletreal> p,
    FArrView<mgletreal> pbuffer,
    const FArrView<mgletreal> dx,
    const FArrView<mgletreal> dy,
    const FArrView<mgletreal> dz,
    const FArrView<mgletreal> ddx,
    const FArrView<mgletreal> ddy,
    const FArrView<mgletreal> ddz,
    mgletint nboundtasks,
    const FArrView<mgletint> boundtasks_lvl,
    const FArrView<mgletint> kkk,
    const FArrView<mgletint> jjj,
    const FArrView<mgletint> iii,
    const FArrView<mgletint> ip3d,
    const FArrView<mgletint> ip1dx,
    const FArrView<mgletint> ip1dy,
    const FArrView<mgletint> ip1dz,
    const FArrView<mgletint> ipbb);

} // namespace mglet::backend
