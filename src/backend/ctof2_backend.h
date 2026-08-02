#pragma once

#include "f_arr_view.h"
#include "mglet_precision.h"

namespace mglet::backend
{

void process_selftasks_ctof2_backend(
    FArrView<const mgletreal> fc,
    FArrView<mgletreal> ff,
    FArrView<const mgletint> selftasks,
    FArrView<const mgletint> kkk,
    FArrView<const mgletint> jjj,
    FArrView<const mgletint> iii,
    FArrView<const mgletint> ip3d);

} // namespace mglet::backend
