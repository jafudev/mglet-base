#ifdef _MGLET_USE_BACKEND_

extern "C"
{
#include <ISO_Fortran_binding.h>
}
#include <cstddef>

#include <omp.h>

#include "ctof2_backend.h"
#include "errr.h"
#include "mglet_precision.h"

extern "C" void process_selftasks_ctof2_c(
    CFI_cdesc_t* fc,
    CFI_cdesc_t* ff,
    CFI_cdesc_t* selftasks,
    CFI_cdesc_t* kkk,
    CFI_cdesc_t* jjj,
    CFI_cdesc_t* iii,
    CFI_cdesc_t* ip3d)
{
    mglet::backend::process_selftasks_ctof2_backend(
        mglet::backend::FArrView<const mgletreal>(fc),
        mglet::backend::FArrView<mgletreal>(ff),
        mglet::backend::FArrView<const mgletint>(selftasks),
        mglet::backend::FArrView<const mgletint>(kkk),
        mglet::backend::FArrView<const mgletint>(jjj),
        mglet::backend::FArrView<const mgletint>(iii),
        mglet::backend::FArrView<const mgletint>(ip3d));
}

#endif // _MGLET_USE_BACKEND_