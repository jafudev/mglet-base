#pragma once

#include "f_arr_view.h"
#include "mglet_precision.h"

namespace mglet::backend
{

void process_selftasks_ctof2_backend(
    const FArrView<mgletreal> fc,
    FArrView<mgletreal> ff,
    mgletint nselftasks,
    const FArrView<mgletint> selftasks,
    const FArrView<mgletint> kkk,
    const FArrView<mgletint> jjj,
    const FArrView<mgletint> iii,
    const FArrView<mgletint> ip3d);

} // namespace mglet::backend
