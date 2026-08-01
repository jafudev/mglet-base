#pragma once

#include <ISO_Fortran_binding.h>
#include <omp.h>

#include "errr.h"
#include "mglet_precision.h"

namespace mglet::backend
{

template <class T>
class FArrView
{
  public:
    FArrView(CFI_cdesc_t* arr)
    {
        if (arr == nullptr)
        {
            MGLET_ERRR();
        }
        if (arr->base_addr == nullptr)
        {
            MGLET_ERRR();
        }

        host_ptr_ = static_cast<mgletreal*>(arr->base_addr);

        const auto device_num_ = omp_get_default_device();
        device_ptr_ = static_cast<mgletreal*>(omp_get_mapped_ptr(host_ptr_, device_num_));
        if (device_ptr_ == nullptr)
        {
            MGLET_ERRR();
        }

        if (arr->rank <= 0)
        {
            MGLET_ERRR();
        }
        flat_size_ = 1;
        for (std::size_t i = 0; i < arr->rank; ++i)
        {
            flat_size_ *= static_cast<std::size_t>(arr->dim[i].extent);
        }
    }

    T* host_ptr() const
    {
        return host_ptr_;
    }
    T* device_ptr() const
    {
        return device_ptr_;
    }
    std::size_t flat_size() const
    {
        return flat_size_;
    }

  private:
    T* host_ptr_;
    T* device_ptr_;
    std::size_t flat_size_;
};

} // mglet::backend
