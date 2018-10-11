/* Taken and modified from
 * https://github.com/JGU-HPC/parallelprogrammingbook/blob/master/chapter4/all_pairs_distance_matrix/data/mnist_exporter.py
 */

#ifndef HPC_HELPERS_HPP
#define HPC_HELPERS_HPP

#include <iostream>
#include <cstdint>
#include <cmath>
#include <limits>
#include <chrono>

void start(const char *label_);
void stop();

// See http://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
template<typename T>
bool almost_equal(T x, T y, int ulp=2);

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

// template<class T>
// class no_init_t {
// public:

//     static_assert(std::is_fundamental<T>::value &&
//                   std::is_arithmetic<T>::value, 
//                   "wrapped type must be a fundamental, numeric type");

//     //do nothing
//     constexpr no_init_t() noexcept {}

//     //convertible from a T
//     constexpr no_init_t(T value) noexcept: v_(value) {}

//     //act as a T in all conversion contexts
//     constexpr operator T () const noexcept { return v_; }

//     // negation on value and bit level
//     constexpr no_init_t& operator - () noexcept { v_ = -v_; return *this; }
//     constexpr no_init_t& operator ~ () noexcept { v_ = ~v_; return *this; }

//     // prefix increment/decrement operators
//     constexpr no_init_t& operator ++ ()    noexcept { v_++; return *this; }
//     constexpr no_init_t& operator -- ()    noexcept { v_--; return *this; }

//     // postfix increment/decrement operators
//     constexpr no_init_t operator ++ (int) noexcept {
//        auto old(*this);
//        v_++; 
//        return old; 
//     }
//     constexpr no_init_t operator -- (int) noexcept {
//        auto old(*this);
//        v_--; 
//        return old; 
//     }

//     // assignment operators
//     constexpr no_init_t& operator  += (T v) noexcept { v_  += v; return *this; }
//     constexpr no_init_t& operator  -= (T v) noexcept { v_  -= v; return *this; }
//     constexpr no_init_t& operator  *= (T v) noexcept { v_  *= v; return *this; }
//     constexpr no_init_t& operator  /= (T v) noexcept { v_  /= v; return *this; }

//     // bit-wise operators
//     constexpr no_init_t& operator  &= (T v) noexcept { v_  &= v; return *this; }
//     constexpr no_init_t& operator  |= (T v) noexcept { v_  |= v; return *this; }
//     constexpr no_init_t& operator  ^= (T v) noexcept { v_  ^= v; return *this; }
//     constexpr no_init_t& operator >>= (T v) noexcept { v_ >>= v; return *this; }
//     constexpr no_init_t& operator <<= (T v) noexcept { v_ <<= v; return *this; }

// private:
//    T v_;
// };

// transfer constants
#define H2D (cudaMemcpyHostToDevice)
#define D2H (cudaMemcpyDeviceToHost)
#define H2H (cudaMemcpyHostToHost)
#define D2D (cudaMemcpyDeviceToDevice)

#endif