#pragma once

#include "f_arr_view.h"
#include "mglet_precision.h"

namespace mglet::backend
{

void bound_pressure_bp_backend(
    FArrView<mgletreal> p,
    FArrView<const mgletreal> pbuffer,
    FArrView<const mgletreal> bp,
    FArrView<const mgletreal> dx,
    FArrView<const mgletreal> dy,
    FArrView<const mgletreal> dz,
    FArrView<const mgletreal> ddx,
    FArrView<const mgletreal> ddy,
    FArrView<const mgletreal> ddz,
    mgletint nboundtasks,
    FArrView<const mgletint> boundtasks_lvl,
    FArrView<const mgletint> kkk,
    FArrView<const mgletint> jjj,
    FArrView<const mgletint> iii,
    FArrView<const mgletint> ip3d,
    FArrView<const mgletint> ip1dx,
    FArrView<const mgletint> ip1dy,
    FArrView<const mgletint> ip1dz,
    FArrView<const mgletint> ipbb);

void bound_pressure_nobp_backend(
    FArrView<mgletreal> p,
    FArrView<const mgletreal> pbuffer,
    FArrView<const mgletreal> dx,
    FArrView<const mgletreal> dy,
    FArrView<const mgletreal> dz,
    FArrView<const mgletreal> ddx,
    FArrView<const mgletreal> ddy,
    FArrView<const mgletreal> ddz,
    mgletint nboundtasks,
    FArrView<const mgletint> boundtasks_lvl,
    FArrView<const mgletint> kkk,
    FArrView<const mgletint> jjj,
    FArrView<const mgletint> iii,
    FArrView<const mgletint> ip3d,
    FArrView<const mgletint> ip1dx,
    FArrView<const mgletint> ip1dy,
    FArrView<const mgletint> ip1dz,
    FArrView<const mgletint> ipbb);

} // namespace mglet::backend
