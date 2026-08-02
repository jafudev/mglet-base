#ifdef _MGLET_USE_BACKEND_

extern "C"
{
#include <ISO_Fortran_binding.h>
}
#include <cstddef>

#include <omp.h>

#include "conn2_backend.h"
#include "errr.h"
#include "mglet_precision.h"

extern "C"
{

void process_selftasks_conn2_c(
    CFI_cdesc_t* a1,
    CFI_cdesc_t* a2,
    CFI_cdesc_t* a3,
    CFI_cdesc_t* a4,
    CFI_cdesc_t* a5,
    CFI_cdesc_t* a6,
    CFI_cdesc_t* selftasks,
    CFI_cdesc_t* kkk,
    CFI_cdesc_t* jjj,
    CFI_cdesc_t* iii,
    CFI_cdesc_t* ip3d)
{
    mglet::backend::process_selftasks_conn2_backend(
        mglet::backend::FArrView<mgletreal>(a1),
        mglet::backend::FArrView<mgletreal>(a2),
        mglet::backend::FArrView<mgletreal>(a3),
        mglet::backend::FArrView<mgletreal>(a4),
        mglet::backend::FArrView<mgletreal>(a5),
        mglet::backend::FArrView<mgletreal>(a6),
        mglet::backend::FArrView<mgletint>(selftasks),
        mglet::backend::FArrView<mgletint>(kkk),
        mglet::backend::FArrView<mgletint>(jjj),
        mglet::backend::FArrView<mgletint>(iii),
        mglet::backend::FArrView<mgletint>(ip3d));
}
}

#endif // _MGLET_USE_BACKEND_
