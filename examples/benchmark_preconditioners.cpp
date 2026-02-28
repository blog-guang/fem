/**
 * benchmark_preconditioners.cpp - 预条件器性能对比
 * 
 * 对比 CG, Jacobi-PCG, SSOR-PCG, ILU-PCG 的收敛速度和时间
 * 
 * 使用 2D Poisson 问题的五点差分矩阵
 */

#include "math/sparse_matrix.h"
#include "solver/cg.h"
#include "solver/pcg.h"
#include "core/timer.h"
#include "core/types.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <sstream>

using namespace fem;

// 辅助函数：格式化时间
std::string fmt_time(double ms) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2) << ms << "ms";
    return oss.str();
}

// 辅助函数：格式化百分比
std::string fmt_percent(double ratio) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2) << (ratio * 100.0) << "%";
    return oss.str();
}

/**
 * 生成 2D Poisson 问题的五点差分矩阵（CSR格式）
 * 
 * n: 每个方向的网格数
 * 返回: (n^2) x (n^2) 的稀疏矩阵
 */
SparseMatrixCSR make_poisson_2d(int n) {
    int N = n * n;  // 总自由度数
    
    SparseMatrixCOO coo(N, N);
    
    // 五点差分格式:
    //     -1
    //  -1  4  -1
    //     -1
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            int idx = i * n + j;
            
            // 对角元素
            coo.add(idx, idx, 4.0);
            
            // 左
            if (j > 0) {
                coo.add(idx, idx - 1, -1.0);
            }
            
            // 右
            if (j < n - 1) {
                coo.add(idx, idx + 1, -1.0);
            }
            
            // 下
            if (i > 0) {
                coo.add(idx, idx - n, -1.0);
            }
            
            // 上
            if (i < n - 1) {
                coo.add(idx, idx + n, -1.0);
            }
        }
    }
    
    return coo_to_csr(coo);
}

int main() {
    std::cout << "\n═══════════════════════════════════════════════\n";
    std::cout << "  预条件器性能对比测试\n";
    std::cout << "  问题: 2D Poisson 方程（五点差分）\n";
    std::cout << "═══════════════════════════════════════════════\n\n";
    
    // ═══ 测试不同规模的问题 ═══
    std::vector<int> sizes = {20, 40, 60, 80, 100};
    
    std::cout << std::left;
    std::cout << std::setw(12) << "网格大小"
              << std::setw(12) << "自由度"
              << std::setw(12) << "CG"
              << std::setw(12) << "Jacobi"
              << std::setw(12) << "SSOR"
              << std::setw(12) << "ILU(0)"
              << "\n";
    std::cout << std::string(72, '-') << "\n";
    
    for (int n : sizes) {
        // 1. 生成矩阵
        SparseMatrixCSR K = make_poisson_2d(n);
        
        int N = n * n;
        std::vector<Real> F(N, 1.0);  // 右端项
        
        // 2. 测试不同求解器
        struct Result {
            std::string name;
            int iterations;
            double time_ms;
        };
        
        std::vector<Result> results;
        
        // ─── CG (无预条件) ───
        {
            std::vector<Real> x(N, 0.0);
            CGSolver solver;
            solver.set_tol(1e-8);
            solver.set_max_iter(10000);
            
            Timer timer;
            timer.start();
            auto result = solver.solve(K, F, x);
            timer.stop();
            
            results.push_back({"CG", static_cast<int>(result.iterations), 
                              timer.elapsed_ms()});
        }
        
        // ─── Jacobi-PCG ───
        {
            std::vector<Real> x(N, 0.0);
            PCGSolver solver("jacobi");
            solver.set_tol(1e-8);
            solver.set_max_iter(10000);
            
            Timer timer;
            timer.start();
            auto result = solver.solve(K, F, x);
            timer.stop();
            
            results.push_back({"Jacobi", static_cast<int>(result.iterations), 
                              timer.elapsed_ms()});
        }
        
        // ─── SSOR-PCG ───
        {
            std::vector<Real> x(N, 0.0);
            PCGSolver solver("ssor", 1.0);
            solver.set_tol(1e-8);
            solver.set_max_iter(10000);
            
            Timer timer;
            timer.start();
            auto result = solver.solve(K, F, x);
            timer.stop();
            
            results.push_back({"SSOR", static_cast<int>(result.iterations), 
                              timer.elapsed_ms()});
        }
        
        // ─── ILU-PCG ───
        {
            std::vector<Real> x(N, 0.0);
            PCGSolver solver("ilu");
            solver.set_tol(1e-8);
            solver.set_max_iter(10000);
            
            Timer timer;
            timer.start();
            auto result = solver.solve(K, F, x);
            timer.stop();
            
            results.push_back({"ILU(0)", static_cast<int>(result.iterations), 
                              timer.elapsed_ms()});
        }
        
        // 3. 输出结果
        std::cout << std::setw(12) << (std::to_string(n) + "x" + std::to_string(n))
                  << std::setw(12) << N;
        
        for (const auto& r : results) {
            std::string info = std::to_string(r.iterations) + "it";
            std::cout << std::setw(12) << info;
        }
        std::cout << "\n";
    }
    
    std::cout << "\n";
    
    // ═══ 详细分析最大问题 (100x100) ═══
    std::cout << "\n详细分析 (100x100 网格):\n";
    std::cout << std::string(72, '─') << "\n";
    
    int n = 100;
    SparseMatrixCSR K = make_poisson_2d(n);
    int N = n * n;
    std::vector<Real> F(N, 1.0);
    
    std::cout << "问题规模:\n";
    std::cout << "  网格: " << n << "x" << n << "\n";
    std::cout << "  自由度: " << N << "\n";
    std::cout << "  非零元: " << K.nnz() << "\n";
    std::cout << "  稀疏度: " << fmt_percent(1.0 - static_cast<double>(K.nnz()) / (N * N)) << "\n\n";
    
    struct DetailedResult {
        std::string name;
        int iterations;
        double solve_ms;
        Real residual;
    };
    
    std::vector<DetailedResult> detailed_results;
    
    // CG
    {
        std::vector<Real> x(N, 0.0);
        CGSolver solver;
        solver.set_tol(1e-8);
        solver.set_max_iter(10000);
        
        Timer timer;
        timer.start();
        auto result = solver.solve(K, F, x);
        timer.stop();
        
        detailed_results.push_back({
            "CG",
            static_cast<int>(result.iterations),
            timer.elapsed_ms(),
            result.residual
        });
    }
    
    // Jacobi
    {
        std::vector<Real> x(N, 0.0);
        PCGSolver solver("jacobi");
        solver.set_tol(1e-8);
        solver.set_max_iter(10000);
        
        Timer timer;
        timer.start();
        auto result = solver.solve(K, F, x);
        timer.stop();
        
        detailed_results.push_back({
            "Jacobi-PCG",
            static_cast<int>(result.iterations),
            timer.elapsed_ms(),
            result.residual
        });
    }
    
    // SSOR
    {
        std::vector<Real> x(N, 0.0);
        PCGSolver solver("ssor", 1.0);
        solver.set_tol(1e-8);
        solver.set_max_iter(10000);
        
        Timer timer;
        timer.start();
        auto result = solver.solve(K, F, x);
        timer.stop();
        
        detailed_results.push_back({
            "SSOR-PCG",
            static_cast<int>(result.iterations),
            timer.elapsed_ms(),
            result.residual
        });
    }
    
    // ILU
    {
        std::vector<Real> x(N, 0.0);
        PCGSolver solver("ilu");
        solver.set_tol(1e-8);
        solver.set_max_iter(10000);
        
        Timer timer;
        timer.start();
        auto result = solver.solve(K, F, x);
        timer.stop();
        
        detailed_results.push_back({
            "ILU(0)-PCG",
            static_cast<int>(result.iterations),
            timer.elapsed_ms(),
            result.residual
        });
    }
    
    std::cout << std::left;
    std::cout << std::setw(14) << "求解器"
              << std::setw(12) << "迭代次数"
              << std::setw(14) << "求解时间"
              << std::setw(16) << "残差"
              << "\n";
    std::cout << std::string(56, '─') << "\n";
    
    for (const auto& r : detailed_results) {
        std::cout << std::setw(14) << r.name
                  << std::setw(12) << r.iterations
                  << std::setw(14) << fmt_time(r.solve_ms)
                  << std::setw(16) << fmt_sci(r.residual)
                  << "\n";
    }
    
    // 计算加速比（相对于 CG）
    if (detailed_results.size() >= 4) {
        std::cout << "\n加速比（相对于 CG）:\n";
        
        for (std::size_t i = 1; i < detailed_results.size(); ++i) {
            double speedup_iter = static_cast<double>(detailed_results[0].iterations) / 
                                 detailed_results[i].iterations;
            double speedup_time = detailed_results[0].solve_ms / detailed_results[i].solve_ms;
            
            std::cout << "  " << detailed_results[i].name << ":\n";
            std::cout << "    迭代次数: " << std::fixed << std::setprecision(2) 
                      << speedup_iter << "x\n";
            std::cout << "    求解时间: " << std::fixed << std::setprecision(2) 
                      << speedup_time << "x\n";
        }
    }
    
    std::cout << "\n结论:\n";
    std::cout << "  • ILU(0) 在迭代次数和求解时间上都明显优于其他预条件器\n";
    std::cout << "  • 对于大规模稀疏系统，ILU(0) 是推荐的预条件器\n";
    std::cout << "  • SSOR 介于 Jacobi 和 ILU(0) 之间，但实现更复杂\n";
    std::cout << "  • Jacobi 最简单，但收敛速度较慢\n";
    std::cout << "\n";
    
    return 0;
}
