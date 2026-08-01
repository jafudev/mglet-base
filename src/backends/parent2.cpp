
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
    const auto a1_view = mglet::backend::FArrView<mgletreal>(a1);
    const auto a2_view = mglet::backend::FArrView<mgletreal>(a2);
    const auto a3_view = mglet::backend::FArrView<mgletreal>(a3);
    const auto a4_view = mglet::backend::FArrView<mgletreal>(a4);
    const auto a5_view = mglet::backend::FArrView<mgletreal>(a5);
    const auto a6_view = mglet::backend::FArrView<mgletreal>(a6);
    const auto b1_view = mglet::backend::FArrView<mgletreal>(b1);
    const auto b2_view = mglet::backend::FArrView<mgletreal>(b2);
    const auto b3_view = mglet::backend::FArrView<mgletreal>(b3);
    const auto b4_view = mglet::backend::FArrView<mgletreal>(b4);
    const auto b5_view = mglet::backend::FArrView<mgletreal>(b5);
    const auto b6_view = mglet::backend::FArrView<mgletreal>(b6);
    const auto selftasks_view = mglet::backend::FArrView<mgletint>(selftasks);
    const auto ip3d_view = mglet::backend::FArrView<mgletint>(ip3d);
    const auto kkk_view = mglet::backend::FArrView<mgletint>(kkk);
    const auto jjj_view = mglet::backend::FArrView<mgletint>(jjj);
    const auto iii_view = mglet::backend::FArrView<mgletint>(iii);

    mglet::backend::process_selftasks_backend(
        a1_view,
        a2_view,
        a3_view,
        a4_view,
        a5_view,
        a6_view,
        b1_view,
        b2_view,
        b3_view,
        b4_view,
        b5_view,
        b6_view,
        selftasks_view,
        ip3d_view,
        kkk_view,
        jjj_view,
        iii_view);
}

#endif // _MGLET_USE_BACKEND_
