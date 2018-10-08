/* Taken and modified from
 * https://github.com/JGU-HPC/parallelprogrammingbook/blob/master/chapter4/all_pairs_distance_matrix/data/mnist_exporter.py
 */

#ifndef HPC_HELPERS_HPP
#define HPC_HELPERS_HPP

#include <iostream>
#include <cstdint>
#include <cmath>
#include <limits>
#include "cuda_runtime.h"

#define TIMERSTART(label)                                                  \
    cudaEvent_t start##label, stop##label;                                 \
    float time##label;                                                     \
    cudaEventCreate(&start##label);                                        \
    cudaEventCreate(&stop##label);                                         \
    cudaEventRecord(start##label, 0);

#define TIMERSTOP(label)                                                   \
    cudaEventRecord(stop##label, 0);                                       \
    cudaEventSynchronize(stop##label);                                     \
    cudaEventElapsedTime(&time##label, start##label, stop##label);         \
    std::cout << "TIMING: " << time##label << " ms (" << #label << ")"     \
                << std::endl;

// See http://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
template<typename T>
bool almost_equal(T x, T y, int ulp=2) {
    // the machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    T epsilon = std::numeric_limits<T>::epsilon();
    // Wee need to adjust epsilon because floats have a
    // precision gap of about 7 digits
    bool a = std::fabs(x-y) < epsilon * std::fabs(x+y) * ulp;
    // Is the difference smaller than 1.17549e-38? Not a normal result actually.
    bool b = std::fabs(x-y) < std::numeric_limits<T>::min();
    return a || b;
}

#define CUERR {                                                            \
    cudaError_t err;                                                       \
    if ((err = cudaGetLastError()) != cudaSuccess) {                       \
        std::cout << "CUDA error: " << cudaGetErrorString(err) << " : "    \
                    << __FILE__ << ", line " << __LINE__ << std::endl;       \
        exit(1);                                                           \
    }                                                                      \
}

// safe division
#define SDIV(x,y)(((x)+(y)-1)/(y))

// no_init_t
#include <type_traits>

template<class T>
class no_init_t {
public:

    static_assert(std::is_fundamental<T>::value &&
                  std::is_arithmetic<T>::value, 
                  "wrapped type must be a fundamental, numeric type");

    //do nothing
    constexpr no_init_t() noexcept {}

    //convertible from a T
    constexpr no_init_t(T value) noexcept: v_(value) {}

    //act as a T in all conversion contexts
    constexpr operator T () const noexcept { return v_; }

    // negation on value and bit level
    constexpr no_init_t& operator - () noexcept { v_ = -v_; return *this; }
    constexpr no_init_t& operator ~ () noexcept { v_ = ~v_; return *this; }

    // prefix increment/decrement operators
    constexpr no_init_t& operator ++ ()    noexcept { v_++; return *this; }
    constexpr no_init_t& operator -- ()    noexcept { v_--; return *this; }

    // postfix increment/decrement operators
    constexpr no_init_t operator ++ (int) noexcept {
       auto old(*this);
       v_++; 
       return old; 
    }
    constexpr no_init_t operator -- (int) noexcept {
       auto old(*this);
       v_--; 
       return old; 
    }

    // assignment operators
    constexpr no_init_t& operator  += (T v) noexcept { v_  += v; return *this; }
    constexpr no_init_t& operator  -= (T v) noexcept { v_  -= v; return *this; }
    constexpr no_init_t& operator  *= (T v) noexcept { v_  *= v; return *this; }
    constexpr no_init_t& operator  /= (T v) noexcept { v_  /= v; return *this; }

    // bit-wise operators
    constexpr no_init_t& operator  &= (T v) noexcept { v_  &= v; return *this; }
    constexpr no_init_t& operator  |= (T v) noexcept { v_  |= v; return *this; }
    constexpr no_init_t& operator  ^= (T v) noexcept { v_  ^= v; return *this; }
    constexpr no_init_t& operator >>= (T v) noexcept { v_ >>= v; return *this; }
    constexpr no_init_t& operator <<= (T v) noexcept { v_ <<= v; return *this; }

private:
   T v_;
};

// transfer constants
#define H2D (cudaMemcpyHostToDevice)
#define D2H (cudaMemcpyDeviceToHost)
#define H2H (cudaMemcpyHostToHost)
#define D2D (cudaMemcpyDeviceToDevice)

#endif