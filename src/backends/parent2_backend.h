#pragma once

#include <cstddef>

#include "f_arr_view.h"
#include "mglet_precision.h"

namespace mglet::backend
{

void process_selftasks_backend(
    FArrView<mgletreal> a1,
    FArrView<mgletreal> a2,
    FArrView<mgletreal> a3,
    FArrView<mgletreal> a4,
    FArrView<mgletreal> a5,
    FArrView<mgletreal> a6,
    FArrView<mgletreal> b1,
    FArrView<mgletreal> b2,
    FArrView<mgletreal> b3,
    FArrView<mgletreal> b4,
    FArrView<mgletreal> b5,
    FArrView<mgletreal> b6,
    FArrView<mgletint> selftasks,
    FArrView<mgletint> ip3d,
    FArrView<mgletint> kkk,
    FArrView<mgletint> jjj,
    FArrView<mgletint> iii);

} // namespace mglet::backend
