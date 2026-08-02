#ifdef _MGLET_USE_BACKEND_

extern "C"
{
#include <ISO_Fortran_binding.h>
}
#include <cstddef>

#include <omp.h>

#include "f_arr_view.h"
#include "mglet_precision.h"
#include "tooling_interface.h"

extern "C" void set_field_arr_realk_c(CFI_cdesc_t* farr, mgletreal val)
{
    mglet::backend::set_farr_realk(mglet::backend::FArrView<mgletreal>(farr), val);
}

extern "C" void set_field_arr_ifk_c(CFI_cdesc_t* farr, mgletifk val)
{
    mglet::backend::set_farr_ifk(mglet::backend::FArrView<mgletifk>(farr), val);
}

#endif // _MGLET_USE_BACKEND_
