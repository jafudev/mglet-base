#include <iostream>

extern "C" {
#include <ISO_Fortran_binding.h>
}
#include <omp.h>
#include <stdint.h>
#include <stdio.h>

#include "mglet_precision.h"  // declares launch_set_field_arr_realk()
#include "errr.h"

#ifdef __cplusplus
extern "C" {
#endif

void set_field_arr_realk_backend(CFI_cdesc_t *arr, mgletreal val)
{
    if (arr->rank != 1) MGLET_ERRR();

    mgletreal* host_ptr = (mgletreal*)arr->base_addr;
    int device_num = omp_get_default_device();

    mgletreal* dev_ptr = (mgletreal*)omp_get_mapped_ptr(host_ptr, device_num);
    if (dev_ptr == NULL) MGLET_ERRR();


}

#ifdef __cplusplus
}
#endif
