#pragma once

#include <cstddef>

#include "f_arr_view.h"
#include "mglet_precision.h"

namespace mglet::backend
{

void process_selftasks_conn2_backend(
    FArrView<mgletreal> a1,
    FArrView<mgletreal> a2,
    FArrView<mgletreal> a3,
    FArrView<mgletreal> a4,
    FArrView<mgletreal> a5,
    FArrView<mgletreal> a6,
    FArrView<mgletint> selftasks,
    FArrView<mgletint> kkk,
    FArrView<mgletint> jjj,
    FArrView<mgletint> iii,
    FArrView<mgletint> ip3d);

} // namespace mglet::backend
