#include "ctof2_backend.h"

#include <cstddef>
#include <cstdint>

#include <cuda_runtime.h>

#include "cutools.h"
#include "errr.h"
#include "f_arr_view.h"

namespace mglet::backend
{

namespace
{

static constexpr mgletint selftask_size = 8;

__global__ void process_selftasks_ctof2_kernel(
    const mgletreal* __restrict__ fc,
    mgletreal* __restrict__ ff,
    const mgletint* __restrict__ selftasks,
    mgletint nselftasks,
    const mgletint* __restrict__ kkk,
    const mgletint* __restrict__ jjj,
    const mgletint* __restrict__ iii,
    const mgletint* __restrict__ ip3d)
{
    const auto itask = blockIdx.x;
    if (itask >= nselftasks)
    {
        return;
    }

    const mgletint* task = selftasks + itask * selftask_size;

    const auto igridf = task[0] - 1; // 1-based Fortran grid id -> 0-based
    const auto igridc = task[1] - 1;
    const auto ista = task[2];
    const auto jsta = task[3];
    const auto ksta = task[4];
    const auto isto = task[5];
    const auto jsto = task[6];
    const auto ksto = task[7];

    const auto kkf = kkk[igridf];
    const auto jjf = jjj[igridf];
    const auto iif = iii[igridf];
    const auto kkc = kkk[igridc];
    const auto jjc = jjj[igridc];

    const auto ip3f = ip3d[igridf] - 1;
    const auto ip3c = ip3d[igridc] - 1;

    const auto n = kkf * jjf * iif;

    for (mgletint lin = threadIdx.x; lin < n; lin += blockDim.x)
    {
        const auto k = lin % kkf + 1;
        const auto j = (lin / kkf) % jjf + 1;
        const auto i = lin / (kkf * jjf) + 1;

        const auto ic = (i - 1) / 2 + ista;
        const auto jc = (j - 1) / 2 + jsta;
        const auto kc = (k - 1) / 2 + ksta;

        if (ic > isto || jc > jsto || kc > ksto)
        {
            continue;
        }

        const auto ff_idx = ip3f + (k - 1) + (j - 1) * kkf + (i - 1) * kkf * jjf;
        const auto fc_idx = ip3c + (kc - 1) + (jc - 1) * kkc + (ic - 1) * kkc * jjc;

        ff[ff_idx] = fc[fc_idx];
    }
}

} // namespace

void process_selftasks_ctof2_backend(
    FArrView<const mgletreal> fc,
    FArrView<mgletreal> ff,
    FArrView<const mgletint> selftasks,
    FArrView<const mgletint> kkk,
    FArrView<const mgletint> jjj,
    FArrView<const mgletint> iii,
    FArrView<const mgletint> ip3d)
{
    const auto ntasks = static_cast<int>((selftasks.flat_size() / selftask_size)) - 1;

    if (ntasks < 0)
    {
        MGLET_ERRR();
    }

    if (ntasks == 0)
    {
        return;
    }

    const auto threads = ::dim3{256};
    const auto blocks = ::dim3{static_cast<unsigned>(ntasks)};

    process_selftasks_ctof2_kernel<<<blocks, threads>>>(
        fc.device_ptr(),
        ff.device_ptr(),
        selftasks.device_ptr(),
        ntasks,
        kkk.device_ptr(),
        jjj.device_ptr(),
        iii.device_ptr(),
        ip3d.device_ptr());

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
}

} // namespace mglet::backend
