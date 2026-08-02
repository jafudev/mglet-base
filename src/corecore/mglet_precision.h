#ifndef __MGLET_PRECISION_H__
#define __MGLET_PRECISION_H__

#include <cstdint>

#ifdef __cplusplus
extern "C" {
#endif
#include <ISO_Fortran_binding.h>
#ifdef __cplusplus
}
#endif

#ifdef _MGLET_DOUBLE_PRECISION_
typedef double mgletreal;
#define CFI_type_mgletreal CFI_type_double
#else
typedef float mgletreal;
#define CFI_type_mgletreal CFI_type_float
#endif

#ifdef _MGLET_INT64_
typedef std::int64_t mgletint;

// If clang is used as a "companion compiler" to NAG nagfor, this will fail or
// not work. I am not sure if that is supported, though. Usually we use GCC
// for this.
// Cray cc/CC also set __clang__ as a macro - but does not behave like it. So
// We need to check for __cray__ first
#if __cray__
#define CFI_type_mgletint CFI_type_long_long
#elif __NVCOMPILER || __clang__
#define CFI_type_mgletint CFI_type_int64_t
#else
#define CFI_type_mgletint CFI_type_long_long
#endif

#else
typedef std::int32_t mgletint;

#ifdef _MGLET_IFK64_
typedef std::int64_t mgletifk;
#else
typedef std::int32_t mgletifk;
#endif

// If clang is used as a "companion compiler" to NAG nagfor, this will fail or
// not work. I am not sure if that is supported, though. Usually we use GCC
// for this.
// Cray cc/CC also set __clang__ as a macro - but does not behave like it. So
// We need to check for __cray__ first
#if __cray__
#define CFI_type_mgletint CFI_type_int
#elif __NVCOMPILER || __clang__
#define CFI_type_mgletint CFI_type_int32_t
#else
#define CFI_type_mgletint CFI_type_int
#endif

#endif

#endif /* __MGLET_PRECISION_H__ */
