
#ifdef _MGLET_USE_BACKEND_

extern "C"
{
#include <ISO_Fortran_binding.h>
}
#include <cstddef>

#include <omp.h>

#include "errr.h"
#include "pressuresolver_backend.h"
#include "mglet_precision.h"

extern "C" void accumulate_pcorr_c(CFI_cdesc_t* dp, CFI_cdesc_t* hilf)
{
    const auto dp_arr = mglet::backend::FArrView<mgletreal>(dp);
    const auto hilf_arr = mglet::backend::FArrView<mgletreal>(hilf);

    mglet::backend::accumulate_pcorr_backend(dp_arr, hilf_arr);
}

#endif // _MGLET_USE_BACKEND_
