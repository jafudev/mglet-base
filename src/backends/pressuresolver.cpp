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

extern "C"
{

void maxabscal_c(
    CFI_cdesc_t* maxabsgrid,
    CFI_cdesc_t* phi,
    CFI_cdesc_t* mygrids,
    mgletint nmygrids,
    CFI_cdesc_t* kkk,
    CFI_cdesc_t* jjj,
    CFI_cdesc_t* iii,
    CFI_cdesc_t* idim3d)
{
    mglet::backend::maxabscal_backend(
        mglet::backend::FArrView<mgletreal>(maxabsgrid),
        mglet::backend::FArrView<mgletreal>(phi),
        mglet::backend::FArrView<mgletint>(mygrids),
        nmygrids,
        mglet::backend::FArrView<mgletint>(kkk),
        mglet::backend::FArrView<mgletint>(jjj),
        mglet::backend::FArrView<mgletint>(iii),
        mglet::backend::FArrView<mgletint>(idim3d));
}

void accumulate_pcorr_c(CFI_cdesc_t* dp, CFI_cdesc_t* hilf)
{
    const auto dp_view = mglet::backend::FArrView<mgletreal>(dp);
    const auto hilf_view = mglet::backend::FArrView<mgletreal>(hilf);

    mglet::backend::accumulate_pcorr_backend(dp_view, hilf_view);
}

void rescal_c(
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

void sipiter1_hyperplane_level_c(
    CFI_cdesc_t* res,
    CFI_cdesc_t* rhs,
    CFI_cdesc_t* siplw,
    CFI_cdesc_t* sipls,
    CFI_cdesc_t* siplb,
    CFI_cdesc_t* siplpr,
    CFI_cdesc_t* miphp,
    CFI_cdesc_t* idxhp,
    mgletint nmygridsonlvl,
    CFI_cdesc_t* mygridsonlvl,
    CFI_cdesc_t* kkk,
    CFI_cdesc_t* jjj,
    CFI_cdesc_t* iii,
    CFI_cdesc_t* ip3d)
{
    mglet::backend::sipiter1_hyperplane_level_backend(
        mglet::backend::FArrView<mgletreal>(res),
        mglet::backend::FArrView<mgletreal>(rhs),
        mglet::backend::FArrView<mgletreal>(siplw),
        mglet::backend::FArrView<mgletreal>(sipls),
        mglet::backend::FArrView<mgletreal>(siplb),
        mglet::backend::FArrView<mgletreal>(siplpr),
        mglet::backend::FArrView<mgletifk>(miphp),
        mglet::backend::FArrView<mgletifk>(idxhp),
        nmygridsonlvl,
        mglet::backend::FArrView<mgletint>(mygridsonlvl),
        mglet::backend::FArrView<mgletint>(kkk),
        mglet::backend::FArrView<mgletint>(jjj),
        mglet::backend::FArrView<mgletint>(iii),
        mglet::backend::FArrView<mgletint>(ip3d));
}

void sipiter2_hyperplane_level_c(
    CFI_cdesc_t* dp,
    CFI_cdesc_t* res,
    CFI_cdesc_t* sipue,
    CFI_cdesc_t* sipun,
    CFI_cdesc_t* siput,
    CFI_cdesc_t* miphp,
    CFI_cdesc_t* idxhp,
    mgletint nmygridsonlvl,
    CFI_cdesc_t* mygridsonlvl,
    CFI_cdesc_t* kkk,
    CFI_cdesc_t* jjj,
    CFI_cdesc_t* iii,
    CFI_cdesc_t* ip3d)
{
    mglet::backend::sipiter2_hyperplane_level_backend(
        mglet::backend::FArrView<mgletreal>(dp),
        mglet::backend::FArrView<mgletreal>(res),
        mglet::backend::FArrView<mgletreal>(sipue),
        mglet::backend::FArrView<mgletreal>(sipun),
        mglet::backend::FArrView<mgletreal>(siput),
        mglet::backend::FArrView<mgletifk>(miphp),
        mglet::backend::FArrView<mgletifk>(idxhp),
        nmygridsonlvl,
        mglet::backend::FArrView<mgletint>(mygridsonlvl),
        mglet::backend::FArrView<mgletint>(kkk),
        mglet::backend::FArrView<mgletint>(jjj),
        mglet::backend::FArrView<mgletint>(iii),
        mglet::backend::FArrView<mgletint>(ip3d));
}

}

#endif // _MGLET_USE_BACKEND_
