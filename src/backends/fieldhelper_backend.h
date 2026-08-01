#pragma once

#include <cstddef>

#include "f_arr_view.h"
#include "mglet_precision.h"

namespace mglet::backend
{

void set_field_arr_realk_backend(FArrView<mgletreal> arr_view, mgletreal val);

}
