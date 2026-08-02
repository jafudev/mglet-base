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
    FArrView<const mgletint> selftasks,
    FArrView<const mgletint> ip3d,
    FArrView<const mgletint> kkk,
    FArrView<const mgletint> jjj,
    FArrView<const mgletint> iii);

} // namespace mglet::backend
