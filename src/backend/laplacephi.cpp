#ifdef _MGLET_USE_BACKEND_

extern "C"
{
#include <ISO_Fortran_binding.h>
}
#include <cstddef>

#include <omp.h>

#include "errr.h"
#include "laplacephi_backend.h"
#include "mglet_precision.h"

extern "C"
{

void laplacephi_c(
    CFI_cdesc_t* res,
    CFI_cdesc_t* phi,
    CFI_cdesc_t* gsaw,
    CFI_cdesc_t* gsae,
    CFI_cdesc_t* gsas,
    CFI_cdesc_t* gsan,
    CFI_cdesc_t* gsab,
    CFI_cdesc_t* gsat,
    CFI_cdesc_t* gsap,
    CFI_cdesc_t* bp,
    CFI_cdesc_t* mygrids,
    CFI_cdesc_t* kkk,
    CFI_cdesc_t* jjj,
    CFI_cdesc_t* iii,
    CFI_cdesc_t* ip3d,
    CFI_cdesc_t* ip1dx,
    CFI_cdesc_t* ip1dy,
    CFI_cdesc_t* ip1dz)
{
    mglet::backend::laplacephi_backend(
        mglet::backend::FArrView<mgletreal>(res),
        mglet::backend::FArrView<const mgletreal>(phi),
        mglet::backend::FArrView<const mgletreal>(gsaw),
        mglet::backend::FArrView<const mgletreal>(gsae),
        mglet::backend::FArrView<const mgletreal>(gsas),
        mglet::backend::FArrView<const mgletreal>(gsan),
        mglet::backend::FArrView<const mgletreal>(gsab),
        mglet::backend::FArrView<const mgletreal>(gsat),
        mglet::backend::FArrView<const mgletreal>(gsap),
        mglet::backend::FArrView<const mgletreal>(bp),
        mglet::backend::FArrView<const mgletint>(mygrids),
        mglet::backend::FArrView<const mgletint>(kkk),
        mglet::backend::FArrView<const mgletint>(jjj),
        mglet::backend::FArrView<const mgletint>(iii),
        mglet::backend::FArrView<const mgletint>(ip3d),
        mglet::backend::FArrView<const mgletint>(ip1dx),
        mglet::backend::FArrView<const mgletint>(ip1dy),
        mglet::backend::FArrView<const mgletint>(ip1dz));
}


void laplacephi_level_c(
    CFI_cdesc_t* res,
    CFI_cdesc_t* phi,
    CFI_cdesc_t* gsaw,
    CFI_cdesc_t* gsae,
    CFI_cdesc_t* gsas,
    CFI_cdesc_t* gsan,
    CFI_cdesc_t* gsab,
    CFI_cdesc_t* gsat,
    CFI_cdesc_t* gsap,
    CFI_cdesc_t* bp,
    mgletint nmygridsonlvl,
    CFI_cdesc_t* mygridsonlvl,
    CFI_cdesc_t* kkk,
    CFI_cdesc_t* jjj,
    CFI_cdesc_t* iii,
    CFI_cdesc_t* ip3d,
    CFI_cdesc_t* ip1dx,
    CFI_cdesc_t* ip1dy,
    CFI_cdesc_t* ip1dz)
{
    mglet::backend::laplacephi_level_backend(
        mglet::backend::FArrView<mgletreal>(res),
        mglet::backend::FArrView<const mgletreal>(phi),
        mglet::backend::FArrView<const mgletreal>(gsaw),
        mglet::backend::FArrView<const mgletreal>(gsae),
        mglet::backend::FArrView<const mgletreal>(gsas),
        mglet::backend::FArrView<const mgletreal>(gsan),
        mglet::backend::FArrView<const mgletreal>(gsab),
        mglet::backend::FArrView<const mgletreal>(gsat),
        mglet::backend::FArrView<const mgletreal>(gsap),
        mglet::backend::FArrView<const mgletreal>(bp),
        mglet::backend::FArrView<const mgletint>(mygridsonlvl),
        mglet::backend::FArrView<const mgletint>(kkk),
        mglet::backend::FArrView<const mgletint>(jjj),
        mglet::backend::FArrView<const mgletint>(iii),
        mglet::backend::FArrView<const mgletint>(ip3d),
        mglet::backend::FArrView<const mgletint>(ip1dx),
        mglet::backend::FArrView<const mgletint>(ip1dy),
        mglet::backend::FArrView<const mgletint>(ip1dz));
}

}

#endif // _MGLET_USE_BACKEND_
