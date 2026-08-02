
#ifdef _MGLET_USE_BACKEND_

extern "C"
{
#include <ISO_Fortran_binding.h>
}
#include <cstddef>

#include <omp.h>

#include "errr.h"
#include "mglet_precision.h"
#include "parent2_backend.h"

extern "C" void process_selftasks_c(
    CFI_cdesc_t* a1,
    CFI_cdesc_t* a2,
    CFI_cdesc_t* a3,
    CFI_cdesc_t* a4,
    CFI_cdesc_t* a5,
    CFI_cdesc_t* a6,
    CFI_cdesc_t* b1,
    CFI_cdesc_t* b2,
    CFI_cdesc_t* b3,
    CFI_cdesc_t* b4,
    CFI_cdesc_t* b5,
    CFI_cdesc_t* b6,
    CFI_cdesc_t* selftasks,
    CFI_cdesc_t* ip3d,
    CFI_cdesc_t* kkk,
    CFI_cdesc_t* jjj,
    CFI_cdesc_t* iii)
{
    mglet::backend::process_selftasks_backend(
        mglet::backend::FArrView<mgletreal>(a1),
        mglet::backend::FArrView<mgletreal>(a2),
        mglet::backend::FArrView<mgletreal>(a3),
        mglet::backend::FArrView<mgletreal>(a4),
        mglet::backend::FArrView<mgletreal>(a5),
        mglet::backend::FArrView<mgletreal>(a6),
        mglet::backend::FArrView<mgletreal>(b1),
        mglet::backend::FArrView<mgletreal>(b2),
        mglet::backend::FArrView<mgletreal>(b3),
        mglet::backend::FArrView<mgletreal>(b4),
        mglet::backend::FArrView<mgletreal>(b5),
        mglet::backend::FArrView<mgletreal>(b6),
        mglet::backend::FArrView<const mgletint>(selftasks),
        mglet::backend::FArrView<const mgletint>(ip3d),
        mglet::backend::FArrView<const mgletint>(kkk),
        mglet::backend::FArrView<const mgletint>(jjj),
        mglet::backend::FArrView<const mgletint>(iii));
}

#endif // _MGLET_USE_BACKEND_
