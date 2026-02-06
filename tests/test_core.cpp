#include <gtest/gtest.h>
#include "core/types.h"
#include "core/logger.h"
#include "core/timer.h"
#include <thread>

using namespace fem;

// ── types ──
TEST(TypesTest, BasicSizes) {
    EXPECT_EQ(sizeof(Real), 8);
    EXPECT_EQ(MAX_NODES, 8);
}

TEST(TypesTest, Vec3) {
    Vec3 v = {1.0, 2.0, 3.0};
    EXPECT_DOUBLE_EQ(v[0], 1.0);
    EXPECT_DOUBLE_EQ(v[1], 2.0);
    EXPECT_DOUBLE_EQ(v[2], 3.0);
}

TEST(TypesTest, FmtSci) {
    std::string s = fmt_sci(1.23e-4, 2);
    EXPECT_NE(s.find("1.23"), std::string::npos);
}

// ── logger ──
TEST(LoggerTest, LevelFilter) {
    auto& log = Logger::instance();
    LogLevel orig = log.level();

    log.set_level(LogLevel::WARN);
    // DEBUG/INFO 不会输出 (无 crash 即可)
    log.debug("should be filtered");
    log.info("should be filtered");
    log.warn("this should appear");

    log.set_level(orig);  // 恢复
}

// ── timer ──
TEST(TimerTest, ElapsedPositive) {
    Timer t;
    t.start();
    // 短暂忙等
    volatile double sum = 0;
    for (int i = 0; i < 100000; ++i) sum += i;
    t.stop();

    EXPECT_GT(t.elapsed_ms(), 0.0);
    EXPECT_GT(t.elapsed_s(),  0.0);
}

TEST(TimerTest, RunningElapsed) {
    Timer t;
    t.start();
    Real mid = t.elapsed_s();  // 仍在 running 中调用
    t.stop();
    Real final_val = t.elapsed_s();

    EXPECT_GE(final_val, mid);
}
