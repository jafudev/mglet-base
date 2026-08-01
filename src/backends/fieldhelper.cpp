#ifdef _MGLET_USE_BACKEND_

extern "C"
{
#include <ISO_Fortran_binding.h>
}
#include <cstddef>

#include <omp.h>

#include "errr.h"
#include "f_arr_view.h"
#include "fieldhelper_backend.h"
#include "mglet_precision.h"

extern "C" void set_field_arr_realk_c(CFI_cdesc_t* field_arr, mgletreal val)
{
    const auto arr = mglet::backend::FArrView<mgletreal>(field_arr);

    mglet::backend::set_field_arr_realk_backend(arr, val);
}

#endif // _MGLET_USE_BACKEND_
