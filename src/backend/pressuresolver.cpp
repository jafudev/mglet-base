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
#include "tooling_interface.h"

extern "C"
{

void maxabscal_c(
    CFI_cdesc_t* maxabsgrid,
    CFI_cdesc_t* phi,
    CFI_cdesc_t* mygrids,
    CFI_cdesc_t* kkk,
    CFI_cdesc_t* jjj,
    CFI_cdesc_t* iii,
    CFI_cdesc_t* idim3d)
{
    mglet::backend::maxabscal_backend(
        mglet::backend::FArrView<mgletreal>(maxabsgrid),
        mglet::backend::FArrView<const mgletreal>(phi),
        mglet::backend::FArrView<const mgletint>(mygrids),
        mglet::backend::FArrView<const mgletint>(kkk),
        mglet::backend::FArrView<const mgletint>(jjj),
        mglet::backend::FArrView<const mgletint>(iii),
        mglet::backend::FArrView<const mgletint>(idim3d));
}

void accumulate_pcorr_c(CFI_cdesc_t* dp, CFI_cdesc_t* hilf)
{
    mglet::backend::add_farr_realk(
        mglet::backend::FArrView<mgletreal>(dp), mglet::backend::FArrView<const mgletreal>(hilf));
}

void rescal_c(
    CFI_cdesc_t* rhs,
    CFI_cdesc_t* res,
    CFI_cdesc_t* mygrids,
    CFI_cdesc_t* kkk,
    CFI_cdesc_t* jjj,
    CFI_cdesc_t* iii,
    CFI_cdesc_t* idim3d)
{
    mglet::backend::rescal_backend(
        mglet::backend::FArrView<mgletreal>(rhs),
        mglet::backend::FArrView<const mgletreal>(res),
        mglet::backend::FArrView<const mgletint>(mygrids),
        mglet::backend::FArrView<const mgletint>(kkk),
        mglet::backend::FArrView<const mgletint>(jjj),
        mglet::backend::FArrView<const mgletint>(iii),
        mglet::backend::FArrView<const mgletint>(idim3d));
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
    CFI_cdesc_t* mygridsonlvl,
    CFI_cdesc_t* kkk,
    CFI_cdesc_t* jjj,
    CFI_cdesc_t* iii,
    CFI_cdesc_t* ip3d)
{
    mglet::backend::sipiter1_hyperplane_level_backend(
        mglet::backend::FArrView<mgletreal>(res),
        mglet::backend::FArrView<const mgletreal>(rhs),
        mglet::backend::FArrView<const mgletreal>(siplw),
        mglet::backend::FArrView<const mgletreal>(sipls),
        mglet::backend::FArrView<const mgletreal>(siplb),
        mglet::backend::FArrView<const mgletreal>(siplpr),
        mglet::backend::FArrView<const mgletifk>(miphp),
        mglet::backend::FArrView<const mgletifk>(idxhp),
        mglet::backend::FArrView<const mgletint>(mygridsonlvl),
        mglet::backend::FArrView<const mgletint>(kkk),
        mglet::backend::FArrView<const mgletint>(jjj),
        mglet::backend::FArrView<const mgletint>(iii),
        mglet::backend::FArrView<const mgletint>(ip3d));
}

void sipiter2_hyperplane_level_c(
    CFI_cdesc_t* dp,
    CFI_cdesc_t* res,
    CFI_cdesc_t* sipue,
    CFI_cdesc_t* sipun,
    CFI_cdesc_t* siput,
    CFI_cdesc_t* miphp,
    CFI_cdesc_t* idxhp,
    CFI_cdesc_t* mygridsonlvl,
    CFI_cdesc_t* kkk,
    CFI_cdesc_t* jjj,
    CFI_cdesc_t* iii,
    CFI_cdesc_t* ip3d)
{
    mglet::backend::sipiter2_hyperplane_level_backend(
        mglet::backend::FArrView<mgletreal>(dp),
        mglet::backend::FArrView<mgletreal>(res),
        mglet::backend::FArrView<const mgletreal>(sipue),
        mglet::backend::FArrView<const mgletreal>(sipun),
        mglet::backend::FArrView<const mgletreal>(siput),
        mglet::backend::FArrView<const mgletifk>(miphp),
        mglet::backend::FArrView<const mgletifk>(idxhp),
        mglet::backend::FArrView<const mgletint>(mygridsonlvl),
        mglet::backend::FArrView<const mgletint>(kkk),
        mglet::backend::FArrView<const mgletint>(jjj),
        mglet::backend::FArrView<const mgletint>(iii),
        mglet::backend::FArrView<const mgletint>(ip3d));
}
}

#endif // _MGLET_USE_BACKEND_
