#pragma once

#include "f_arr_view.h"
#include "mglet_precision.h"

namespace mglet::gpu
{

void set_farr_realk(FArrView<mgletreal> farr, mgletreal val);

void set_farr_ifk(FArrView<mgletifk> farr, mgletifk val);

void add_farr_realk(FArrView<mgletreal> lhs, FArrView<const mgletreal> rhs);

} // namespace mglet::gpu
