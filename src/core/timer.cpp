#include "core/timer.h"
#include "core/logger.h"
#include <iostream>
#include <iomanip>

namespace fem {

void Timer::start() {
    start_   = Clock::now();
    running_ = true;
}

void Timer::stop() {
    end_     = Clock::now();
    running_ = false;
}

Real Timer::elapsed_s() const {
    auto end = running_ ? Clock::now() : end_;
    return std::chrono::duration<Real>(end - start_).count();
}

Real Timer::elapsed_ms() const {
    return elapsed_s() * 1000.0;
}

void Timer::print(std::string_view label, Real seconds) {
    std::cout << std::fixed << std::setprecision(4);
    FEM_INFO(std::string(label) + ": " + std::to_string(seconds) + " s");
}

// ── ScopedTimer ──
ScopedTimer::ScopedTimer(std::string_view label)
    : label_(label)
{
    timer_.start();
}

ScopedTimer::~ScopedTimer() {
    timer_.stop();
    Timer::print(label_, timer_.elapsed_s());
}

}  // namespace fem
