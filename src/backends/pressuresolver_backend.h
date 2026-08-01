#pragma once

#include <cstddef>

#include "f_arr_view.h"
#include "mglet_precision.h"

namespace mglet::backend
{

void accumulate_pcorr_backend(FArrView<mgletreal> dp_arr, const FArrView<mgletreal> hilf_arr);

}
