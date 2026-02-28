#pragma once

#include <cstdint>
#include <string_view>

namespace fem {

enum class LogLevel : uint8_t {
    DEBUG = 0,
    INFO,
    WARN,
    ERROR
};

class Logger {
public:
    static Logger& instance();

    void set_level(LogLevel level) { level_ = level; }
    LogLevel level() const         { return level_; }

    void debug(std::string_view msg);
    void info(std::string_view msg);
    void warn(std::string_view msg);
    void error(std::string_view msg);

private:
    Logger() = default;
    void log(LogLevel lvl, const char* prefix, std::string_view msg);

    LogLevel level_{LogLevel::INFO};
};

// ── 宏接口 ──
#define FEM_DEBUG(msg) fem::Logger::instance().debug(msg)
#define FEM_INFO(msg)  fem::Logger::instance().info(msg)
#define FEM_WARN(msg)  fem::Logger::instance().warn(msg)
#define FEM_ERROR(msg) fem::Logger::instance().error(msg)

}  // namespace fem
