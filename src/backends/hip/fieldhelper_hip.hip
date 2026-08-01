#include "fieldhelper_backend.h"

#include <cstddef>

#include <hip/hip_runtime.h>

#include "errr.h"

namespace mglet::backend {


void set_field_arr_realk_backend(mgletreal* device_ptr, std::size_t n, mgletreal val)
{
    MGLET_ERRR();
}

}  // namespace mglet::backend
