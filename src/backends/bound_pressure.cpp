#ifdef _MGLET_USE_BACKEND_

extern "C"
{
#include <ISO_Fortran_binding.h>
}
#include <cstddef>

#include <omp.h>

#include "bound_pressure_backend.h"
#include "errr.h"
#include "mglet_precision.h"

extern "C"
{

void bound_pressure_bp_c(
    CFI_cdesc_t* p,
    CFI_cdesc_t* pbuffer,
    CFI_cdesc_t* bp,
    CFI_cdesc_t* dx,
    CFI_cdesc_t* dy,
    CFI_cdesc_t* dz,
    CFI_cdesc_t* ddx,
    CFI_cdesc_t* ddy,
    CFI_cdesc_t* ddz,
    mgletint nboundtasks,
    CFI_cdesc_t* boundtasks_lvl,
    CFI_cdesc_t* kkk,
    CFI_cdesc_t* jjj,
    CFI_cdesc_t* iii,
    CFI_cdesc_t* ip3d,
    CFI_cdesc_t* ip1dx,
    CFI_cdesc_t* ip1dy,
    CFI_cdesc_t* ip1dz,
    CFI_cdesc_t* ipbb)
{
    mglet::backend::bound_pressure_bp_backend(
        mglet::backend::FArrView<mgletreal>(p),
        mglet::backend::FArrView<mgletreal>(pbuffer),
        mglet::backend::FArrView<mgletreal>(bp),
        mglet::backend::FArrView<mgletreal>(dx),
        mglet::backend::FArrView<mgletreal>(dy),
        mglet::backend::FArrView<mgletreal>(dz),
        mglet::backend::FArrView<mgletreal>(ddx),
        mglet::backend::FArrView<mgletreal>(ddy),
        mglet::backend::FArrView<mgletreal>(ddz),
        nboundtasks,
        mglet::backend::FArrView<mgletint>(boundtasks_lvl),
        mglet::backend::FArrView<mgletint>(kkk),
        mglet::backend::FArrView<mgletint>(jjj),
        mglet::backend::FArrView<mgletint>(iii),
        mglet::backend::FArrView<mgletint>(ip3d),
        mglet::backend::FArrView<mgletint>(ip1dx),
        mglet::backend::FArrView<mgletint>(ip1dy),
        mglet::backend::FArrView<mgletint>(ip1dz),
        mglet::backend::FArrView<mgletint>(ipbb));
}

void bound_pressure_nobp_c(
    CFI_cdesc_t* p,
    CFI_cdesc_t* pbuffer,
    CFI_cdesc_t* dx,
    CFI_cdesc_t* dy,
    CFI_cdesc_t* dz,
    CFI_cdesc_t* ddx,
    CFI_cdesc_t* ddy,
    CFI_cdesc_t* ddz,
    mgletint nboundtasks,
    CFI_cdesc_t* boundtasks_lvl,
    CFI_cdesc_t* kkk,
    CFI_cdesc_t* jjj,
    CFI_cdesc_t* iii,
    CFI_cdesc_t* ip3d,
    CFI_cdesc_t* ip1dx,
    CFI_cdesc_t* ip1dy,
    CFI_cdesc_t* ip1dz,
    CFI_cdesc_t* ipbb)
{
    mglet::backend::bound_pressure_nobp_backend(
        mglet::backend::FArrView<mgletreal>(p),
        mglet::backend::FArrView<mgletreal>(pbuffer),
        mglet::backend::FArrView<mgletreal>(dx),
        mglet::backend::FArrView<mgletreal>(dy),
        mglet::backend::FArrView<mgletreal>(dz),
        mglet::backend::FArrView<mgletreal>(ddx),
        mglet::backend::FArrView<mgletreal>(ddy),
        mglet::backend::FArrView<mgletreal>(ddz),
        nboundtasks,
        mglet::backend::FArrView<mgletint>(boundtasks_lvl),
        mglet::backend::FArrView<mgletint>(kkk),
        mglet::backend::FArrView<mgletint>(jjj),
        mglet::backend::FArrView<mgletint>(iii),
        mglet::backend::FArrView<mgletint>(ip3d),
        mglet::backend::FArrView<mgletint>(ip1dx),
        mglet::backend::FArrView<mgletint>(ip1dy),
        mglet::backend::FArrView<mgletint>(ip1dz),
        mglet::backend::FArrView<mgletint>(ipbb));
}
}

#endif // _MGLET_USE_BACKEND_
