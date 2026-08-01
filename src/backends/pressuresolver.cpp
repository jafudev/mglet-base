
#ifdef _MGLET_USE_BACKEND_

extern "C"
{
#include <ISO_Fortran_binding.h>
}
#include <cstddef>

#include <omp.h>

#include "errr.h"
#include "mglet_precision.h"
#include "pressuresolver_backend.h"

extern "C" void accumulate_pcorr_c(CFI_cdesc_t* dp, CFI_cdesc_t* hilf)
{
    const auto dp_view = mglet::backend::FArrView<mgletreal>(dp);
    const auto hilf_view = mglet::backend::FArrView<mgletreal>(hilf);

    mglet::backend::accumulate_pcorr_backend(dp_view, hilf_view);
}

extern "C" void rescal_c(
    CFI_cdesc_t* rhs,
    CFI_cdesc_t* res,
    mgletint nmygrids,
    CFI_cdesc_t* mygrids,
    CFI_cdesc_t* kkk,
    CFI_cdesc_t* jjj,
    CFI_cdesc_t* iii,
    CFI_cdesc_t* idim3d)
{
    const auto rhs_view = mglet::backend::FArrView<mgletreal>(rhs);
    const auto res_view = mglet::backend::FArrView<mgletreal>(res);
    const auto mygrids_view = mglet::backend::FArrView<mgletint>(mygrids);
    const auto kkk_view = mglet::backend::FArrView<mgletint>(kkk);
    const auto jjj_view = mglet::backend::FArrView<mgletint>(jjj);
    const auto iii_view = mglet::backend::FArrView<mgletint>(iii);
    const auto idim3d_view = mglet::backend::FArrView<mgletint>(idim3d);

    mglet::backend::rescal_backend(
        rhs_view, res_view, nmygrids, mygrids_view, kkk_view, jjj_view, iii_view, idim3d_view);
}

#endif // _MGLET_USE_BACKEND_
