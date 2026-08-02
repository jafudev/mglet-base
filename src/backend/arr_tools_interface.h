#pragma once

#include "mapped_arr_view.h"
#include "mglet_precision.h"

namespace mglet::gpu
{

void set_farr_realk(MappedArrView<mgletreal> farr, mgletreal val);

void set_farr_ifk(MappedArrView<mgletifk> farr, mgletifk val);

void add_farr_realk(MappedArrView<mgletreal> lhs, MappedArrView<const mgletreal> rhs);

} // namespace mglet::gpu
