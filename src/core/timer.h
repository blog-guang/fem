#pragma once

#include "core/types.h"
#include <chrono>
#include <string_view>

namespace fem {

// ── 高精度计时器 ──
class Timer {
public:
    void start();
    void stop();

    Real elapsed_s()  const;   // 秒
    Real elapsed_ms() const;   // 毫秒

    // 打印辅助: "label: X.XXX s"
    static void print(std::string_view label, Real seconds);

private:
    using Clock = std::chrono::high_resolution_clock;
    Clock::time_point start_{};
    Clock::time_point end_{};
    bool running_{false};
};

// ── RAII 计时器: 作用域退出自动 stop + print ──
class ScopedTimer {
public:
    explicit ScopedTimer(std::string_view label);
    ~ScopedTimer();

private:
    std::string label_;
    Timer timer_;
};

#define FEM_TIMER(label) fem::ScopedTimer _timer_##__LINE__(label)

}  // namespace fem
