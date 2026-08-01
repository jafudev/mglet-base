#pragma once

#include <cstddef>

#include "f_arr_view.h"
#include "mglet_precision.h"

namespace mglet::backend
{

void accumulate_pcorr_backend(FArrView<mgletreal> dp_view, const FArrView<mgletreal> hilf_view);
void rescal_backend(
    FArrView<mgletreal> rhs_view,
    const FArrView<mgletreal> res_view,
    mgletint nmygrids,
    const FArrView<mgletint> mygrids,
    const FArrView<mgletint> kkk,
    const FArrView<mgletint> jjj,
    const FArrView<mgletint> iii,
    const FArrView<mgletint> ip3d
);

} // namespace mglet::backend
