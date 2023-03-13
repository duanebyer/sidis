#ifndef SIDIS_UTILITY_HPP
#define SIDIS_UTILITY_HPP

namespace sidis {

#if defined(__GNUC__)
#define SIDIS_UNREACHABLE() __builtin_unreachable()
#elif defined(_MSC_VER)
#define SIDIS_UNREACHABLE() __assume(0)
#else
#define SIDIS_UNREACHABLE() exit(1)
#endif

}

#endif

