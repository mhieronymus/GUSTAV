#include "helper/timing.hpp"


std::chrono::time_point<std::chrono::system_clock> start_c, stop_c;
std::string label;

void start(const char *label_) {
    start_c = std::chrono::system_clock::now();
    label = label_;
}

void stop() {
    stop_c = std::chrono::system_clock::now();
    std::chrono::duration<double> delta = stop_c-start_c;
    std::cout << "# elapsed time ("<< label <<"): "                       \
              << delta.count()  << "s" << std::endl;
}

// See http://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
template<typename T>
bool almost_equal(T x, T y, int ulp) {
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