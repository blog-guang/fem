#include "core/logger.h"
#include <iostream>
#include <chrono>
#include <iomanip>
#include <ctime>

namespace fem {

Logger& Logger::instance() {
    static Logger logger;
    return logger;
}

void Logger::debug(std::string_view msg) { log(LogLevel::DEBUG, "[DEBUG]", msg); }
void Logger::info(std::string_view msg)  { log(LogLevel::INFO,  "[INFO] ", msg); }
void Logger::warn(std::string_view msg)  { log(LogLevel::WARN,  "[WARN] ", msg); }
void Logger::error(std::string_view msg) { log(LogLevel::ERROR, "[ERROR]", msg); }

void Logger::log(LogLevel lvl, const char* prefix, std::string_view msg) {
    if (lvl < level_) return;

    auto now  = std::chrono::system_clock::now();
    auto time = std::chrono::system_clock::to_time_t(now);
    std::tm* tm = std::localtime(&time);

    std::cout << std::put_time(tm, "%H:%M:%S") << " " << prefix << " " << msg << "\n";
}

}  // namespace fem
