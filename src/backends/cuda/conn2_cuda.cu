#include "parent2_backend.h"

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

static constexpr mgletint num_selftasks = 15;

__device__ void arr_to_arr_task(
    mgletreal* __restrict__ dst,
    const mgletreal* __restrict__ src,
    mgletint kk,
    mgletint jj,
    mgletint istart,
    mgletint jstart,
    mgletint kstart,
    mgletint istart_d,
    mgletint istop_d,
    mgletint jstart_d,
    mgletint jstop_d,
    mgletint kstart_d,
    mgletint kstop_d)
{
    const auto koff = kstart - kstart_d;
    const auto joff = jstart - jstart_d;
    const auto ioff = istart - istart_d;

    const auto ni = istop_d - istart_d + 1;
    const auto nj = jstop_d - jstart_d + 1;
    const auto nk = kstop_d - kstart_d + 1;
    if (ni <= 0 || nj <= 0 || nk <= 0)
    {
        return;
    }

    const auto n = ni * nj * nk;

    for (mgletint lin = threadIdx.x; lin < n; lin += blockDim.x)
    {
        const auto dk = lin % nk;
        const auto dj = (lin / nk) % nj;
        const auto di = lin / ((std::int64_t)nk * nj);

        const auto i = istart_d + di;
        const auto j = jstart_d + dj;
        const auto k = kstart_d + dk;

        const auto dst_idx = (k - 1) + (j - 1) * kk + (i - 1) * kk * jj;
        const auto src_idx = (k + koff - 1) + (j + joff - 1) * kk + (i + ioff - 1) * kk * jj;

        dst[dst_idx] = src[src_idx];
    }
}

__global__ void process_selftasks_conn2_kernel(
    mgletreal* __restrict__ a1,
    mgletreal* __restrict__ a2,
    mgletreal* __restrict__ a3,
    mgletreal* __restrict__ a4,
    mgletreal* __restrict__ a5,
    mgletreal* __restrict__ a6,
    mgletint nstasks,
    const mgletint* __restrict__ stasks,
    const mgletint* __restrict__ kkk,
    const mgletint* __restrict__ jjj,
    const mgletint* __restrict__ ip3d)
{
    const int block_idx = blockIdx.x;
    if (block_idx >= nstasks)
    {
        return;
    }

    const mgletint* task = stasks + block_idx * num_selftasks;

    const auto fieldid = task[0];
    const auto igrid = task[1] - 1;
    const auto igrid_d = task[2] - 1;
    const auto istart = task[3];
    // const auto istop = task[4];
    const auto jstart = task[5];
    // const auto jstop = task[6];
    const auto kstart = task[7];
    // const auto kstop = task[8];
    const auto istart_d = task[9];
    const auto istop_d = task[10];
    const auto jstart_d = task[11];
    const auto jstop_d = task[12];
    const auto kstart_d = task[13];
    const auto kstop_d = task[14];

    const auto kk = kkk[igrid];
    const auto jj = jjj[igrid];
    const auto ip3 = ip3d[igrid] - 1;
    const auto ip3_d = ip3d[igrid_d] - 1;

    mgletreal* arr;
    switch (fieldid)
    {
        case 1: arr = a1; break;
        case 2: arr = a2; break;
        case 3: arr = a3; break;
        case 4: arr = a4; break;
        case 5: arr = a5; break;
        case 6: arr = a6; break;
        default: return;
    }

    arr_to_arr_task(
        arr + ip3_d,
        arr + ip3,
        kk,
        jj,
        istart,
        jstart,
        kstart,
        istart_d,
        istop_d,
        jstart_d,
        jstop_d,
        kstart_d,
        kstop_d);
}

} // namespace

void process_selftasks_conn2_backend(
    FArrView<mgletreal> a1,
    FArrView<mgletreal> a2,
    FArrView<mgletreal> a3,
    FArrView<mgletreal> a4,
    FArrView<mgletreal> a5,
    FArrView<mgletreal> a6,
    FArrView<const mgletint> selftasks,
    FArrView<const mgletint> kkk,
    FArrView<const mgletint> jjj,
    FArrView<const mgletint> iii,
    FArrView<const mgletint> ip3d)
{
    const auto ntasks = static_cast<int>(selftasks.flat_size() / num_selftasks);

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

    process_selftasks_conn2_kernel<<<blocks, threads>>>(
        a1.device_ptr(),
        a2.device_ptr(),
        a3.device_ptr(),
        a4.device_ptr(),
        a5.device_ptr(),
        a6.device_ptr(),
        ntasks,
        selftasks.device_ptr(),
        kkk.device_ptr(),
        jjj.device_ptr(),
        ip3d.device_ptr());

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
}

} // namespace mglet::backend
