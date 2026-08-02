#include <cstdint>

#include "errr.h"
#include "gpu_check.h"
#include "gpu_include.h"
#include "mapped_arr_view.h"

namespace mglet::gpu
{

namespace
{

constexpr int num_selftasks = 15;

__device__ __forceinline__ mgletreal get_parent_face_value(
    const mgletreal* __restrict__ arr, int kk, int jj,
    int istart, int istop, int jstart, int jstop, int kstart, int kstop,
    int jc, int ic)
{
    // Same three-way branch as the Fortran: constant start/stop pair
    // identifies which direction is normal to the face.
    int k, j, i;
    if (istart == istop) {
        k = kstart + jc - 1; j = jstart + ic - 1; i = istart;
    } else if (jstart == jstop) {
        k = kstart + jc - 1; j = jstart; i = istart + ic - 1;
    } else {
        k = kstart; j = jstart + jc - 1; i = istart + ic - 1;
    }
    return arr[(k - 1) + (long)(j - 1) * kk + (long)(i - 1) * kk * jj];
}

// Both prolong_face_first_direction and prolong_face_second_direction in
// the Fortran have identical bodies -- one shared device function covers both.
__device__ __forceinline__ mgletreal prolong_face(
    mgletreal val_c, mgletreal val_m1, int stag, bool odd)
{
    if (stag == 0) return val_c;
    return odd ? mgletreal(0.5) * (val_c + val_m1) : val_c;
}

__device__ void arr_to_buffers_task(
    const mgletreal* __restrict__ arr, mgletreal* __restrict__ buffer,
    int kk, int jj,
    int ibb, int istart, int istop, int jstart, int jstop,
    int kstart, int kstop, int jj2d, int ii2d, int stag1, int stag2)
{
    // ii2d/jj2d are per-task and arbitrary in size, so this grid-strides
    // over both -- same pattern as the "support any kk/jj/ii" kernels above.
    for (int i = 1 + threadIdx.y; i <= ii2d; i += blockDim.y) {
        for (int j = 1 + threadIdx.x; j <= jj2d; j += blockDim.x) {

            const bool odd_j = (j % 2) == 1;
            const bool odd_i = (i % 2) == 1;

            const int jc = 2 + (j - 1) / 2;
            const int ic = 2 + (i - 1) / 2;

            mgletreal val_c = get_parent_face_value(
                arr, kk, jj, istart, istop, jstart, jstop, kstart, kstop, jc, ic);
            mgletreal val_jm1 = val_c;
            if (stag1 == 1 && odd_j) {
                val_jm1 = get_parent_face_value(
                    arr, kk, jj, istart, istop, jstart, jstop, kstart, kstop, jc - 1, ic);
            }
            mgletreal val_i = prolong_face(val_c, val_jm1, stag1, odd_j);

            mgletreal val_im1 = val_i;
            if (stag2 == 1 && odd_i) {
                val_c = get_parent_face_value(
                    arr, kk, jj, istart, istop, jstart, jstop, kstart, kstop, jc, ic - 1);
                val_jm1 = val_c;
                if (stag1 == 1 && odd_j) {
                    val_jm1 = get_parent_face_value(
                        arr, kk, jj, istart, istop, jstart, jstop, kstart, kstop, jc - 1, ic - 1);
                }
                val_im1 = prolong_face(val_c, val_jm1, stag1, odd_j);
            }

            const mgletreal val_out = prolong_face(val_i, val_im1, stag2, odd_i);

            const long idst = ibb + (j - 1) + (long)(i - 1) * jj2d;
            buffer[idst - 1] = val_out;   // buffer is Fortran 1-based; -1 for device offset
        }
    }
}

__global__ void process_selftasks_parent2_kernel(
    const mgletint* __restrict__ stasks, int selftasksize, int ntasks,
    const mgletint* __restrict__ ip3d,
    const mgletint* __restrict__ kkgrid,
    const mgletint* __restrict__ jjgrid,
    mgletreal* __restrict__ a1, mgletreal* __restrict__ a2, mgletreal* __restrict__ a3,
    mgletreal* __restrict__ a4, mgletreal* __restrict__ a5, mgletreal* __restrict__ a6,
    mgletreal* __restrict__ b1, mgletreal* __restrict__ b2, mgletreal* __restrict__ b3,
    mgletreal* __restrict__ b4, mgletreal* __restrict__ b5, mgletreal* __restrict__ b6)
{
    const int itask = blockIdx.x;
    if (itask >= ntasks) return;

    const mgletint* task = stasks + (long)itask * selftasksize;

    const int fieldid = task[0];
    const int igridc  = task[1] - 1;    // grid ids are 1-based in Fortran
    const int ibb     = task[2];
    const int istart  = task[3];
    const int istop   = task[4];
    const int jstart  = task[5];
    const int jstop   = task[6];
    const int kstart  = task[7];
    const int kstop   = task[8];
    // task[9], task[10] unused, matching the Fortran (indices 10, 11 never read)
    const int jj2d  = task[11];
    const int ii2d  = task[12];
    const int stag1 = task[13];
    const int stag2 = task[14];

    const int kk = kkgrid[igridc];
    const int jj = jjgrid[igridc];
    const long ip3 = ip3d[igridc] - 1;

    const mgletreal* arr;
    mgletreal* buf;

    switch (fieldid) {
        case 1: arr = a1 + ip3; buf = b1; break;
        case 2: arr = a2 + ip3; buf = b2; break;
        case 3: arr = a3 + ip3; buf = b3; break;
        case 4: arr = a4 + ip3; buf = b4; break;
        case 5: arr = a5 + ip3; buf = b5; break;
        case 6: arr = a6 + ip3; buf = b6; break;
        default: return;   // release-mode equivalent of the Fortran's #ifdef _MGLET_DEBUG_ errr()
    }

    arr_to_buffers_task(arr, buf, kk, jj, ibb, istart, istop, jstart, jstop,
        kstart, kstop, jj2d, ii2d, stag1, stag2);
}

} // namespace

void process_selftasks_parent2(
    MappedArrView<mgletreal> a1,
    MappedArrView<mgletreal> a2,
    MappedArrView<mgletreal> a3,
    MappedArrView<mgletreal> a4,
    MappedArrView<mgletreal> a5,
    MappedArrView<mgletreal> a6,
    MappedArrView<mgletreal> b1,
    MappedArrView<mgletreal> b2,
    MappedArrView<mgletreal> b3,
    MappedArrView<mgletreal> b4,
    MappedArrView<mgletreal> b5,
    MappedArrView<mgletreal> b6,
    MappedArrView<const mgletint> selftasks,
    MappedArrView<const mgletint> ip3d,
    MappedArrView<const mgletint> kkk,
    MappedArrView<const mgletint> jjj,
    MappedArrView<const mgletint> iii)
{
    const auto ntasks = static_cast<int>((selftasks.flat_size() / num_selftasks)) - 1;

    if (ntasks < 0)
    {
        MGLET_ERRR();
    }

    if (ntasks == 0) 
    {
        return;
    }

    const auto threads = ::dim3{32, 8};
    const auto blocks = ::dim3{static_cast<unsigned>(ntasks)};

    process_selftasks_parent2_kernel<<<blocks, threads>>>(
        selftasks.device_ptr(), num_selftasks, ntasks,
        ip3d.device_ptr(), kkk.device_ptr(), jjj.device_ptr(),
        a1.device_ptr(), a2.device_ptr(), a3.device_ptr(),
        a4.device_ptr(), a5.device_ptr(), a6.device_ptr(),
        b1.device_ptr(), b2.device_ptr(), b3.device_ptr(),
        b4.device_ptr(), b5.device_ptr(), b6.device_ptr());

    GPU_CHECK(gpuGetLastError());
    GPU_CHECK(gpuDeviceSynchronize());
}

} // namespace mglet::gpu

#ifdef _MGLET_USE_BACKEND_

extern "C" void process_selftasks_parent2_c(
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
    mglet::gpu::process_selftasks_parent2(
        mglet::gpu::MappedArrView<mgletreal>(a1),
        mglet::gpu::MappedArrView<mgletreal>(a2),
        mglet::gpu::MappedArrView<mgletreal>(a3),
        mglet::gpu::MappedArrView<mgletreal>(a4),
        mglet::gpu::MappedArrView<mgletreal>(a5),
        mglet::gpu::MappedArrView<mgletreal>(a6),
        mglet::gpu::MappedArrView<mgletreal>(b1),
        mglet::gpu::MappedArrView<mgletreal>(b2),
        mglet::gpu::MappedArrView<mgletreal>(b3),
        mglet::gpu::MappedArrView<mgletreal>(b4),
        mglet::gpu::MappedArrView<mgletreal>(b5),
        mglet::gpu::MappedArrView<mgletreal>(b6),
        mglet::gpu::MappedArrView<const mgletint>(selftasks),
        mglet::gpu::MappedArrView<const mgletint>(ip3d),
        mglet::gpu::MappedArrView<const mgletint>(kkk),
        mglet::gpu::MappedArrView<const mgletint>(jjj),
        mglet::gpu::MappedArrView<const mgletint>(iii));
}

#endif // _MGLET_USE_BACKEND_
