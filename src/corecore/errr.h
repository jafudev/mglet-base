#pragma once

#include <cstdint>

#include "mglet_precision.h"

#define MGLET_ERRR() errr_c(__FILE__, static_cast<mgletint>(__LINE__))

extern "C"
{
    [[noreturn]] void errr_c(const char *fname, mgletint line);
}

[[noreturn]] inline void errr(const char *fname, mgletint line)
{
    errr_c(fname, line);
}
