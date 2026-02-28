# Vector/DenseMatrix 代数运算重构计划

## 目标

将所有向量-向量、向量-矩阵的代数运算统一使用 `math` 模块的 `Vector`、`DenseMatrix` 类，通过运算符重载实现代数运算，消除手写循环。

## 原则

- **可读性优先**：`c = a + b` 比 `for (i...) c[i] = a[i] + b[i]` 更清晰
- **类型安全**：`Vector` 的运算符有大小检查，避免越界错误
- **性能相当**：运算符重载会被内联，性能与手写循环相同
- **统一接口**：所有向量运算都用同一套 API

## 当前问题

### 1. 手写辅助函数

**文件：`src/math/cg.cpp`**
```cpp
static Real dot(const std::vector<Real>& a, const std::vector<Real>& b) {
    Real s = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
    return s;
}
```
❌ 应该用 `Vector::dot()`

**文件：`src/math/pcg.cpp`**
```cpp
static Real dot(const std::vector<Real>& a, const std::vector<Real>& b) {
    Real s = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) { s += a[i] * b[i]; }
    return s;
}

static void axpy(Real alpha, const std::vector<Real>& x, std::vector<Real>& y) {
    for (std::size_t i = 0; i < x.size(); ++i) { y[i] += alpha * x[i]; }
}
```
❌ 应该用 `Vector::dot()` 和 `y += alpha * x`

### 2. 循环实现的向量运算

**文件：`src/math/cg.cpp`**
```cpp
// ❌ 手写 axpy
for (std::size_t i = 0; i < n; ++i) {
    x[i] += alpha * p[i];
    r[i] -= alpha * Ap[i];
}

// ❌ 手写线性组合
for (std::size_t i = 0; i < n; ++i) {
    p[i] = r[i] + beta * p[i];
}
```
✅ 应该用：
```cpp
x += alpha * p;
r -= alpha * Ap;
p = r + beta * p;
```

**文件：`src/math/newton_raphson.cpp`**
```cpp
// ❌ 手写取负
for (Index i = 0; i < n; ++i) {
    neg_R[i] = -R[i];
}

// ❌ 手写 axpy
for (Index i = 0; i < n; ++i) {
    u[i] += alpha * du[i];
}

// ❌ 手写范数
Real NewtonRaphsonSolver::compute_norm(const std::vector<Real>& v) const {
    Real sum = 0.0;
    for (Real x : v) { sum += x * x; }
    return std::sqrt(sum);
}
```
✅ 应该用：
```cpp
Vector neg_R = -1.0 * Vector(R);  // 或添加 Vector::operator-()
u_vec += alpha * du_vec;
Real norm = Vector(v).norm();
```

### 3. 预条件器中的循环

**文件：`src/math/pcg.cpp` (各预条件器)**
```cpp
// Jacobi
void JacobiPreconditioner::apply(const std::vector<Real>& r, std::vector<Real>& z) const {
    for (std::size_t i = 0; i < r.size(); ++i) {
        z[i] = diag_inv_[i] * r[i];
    }
}
```
✅ 应该用：
```cpp
Vector r_vec(r);
Vector z_vec = r_vec.hadamard(diag_inv_);  // 需添加 Hadamard 乘积
z.assign(z_vec.data(), z_vec.data() + z_vec.size());
```

## 重构步骤

### Phase 1: 添加缺失的 Vector 功能 ✅

已有的功能：
- ✅ `+`, `-`, `*`, `/` (向量-标量、向量-向量)
- ✅ `+=`, `-=`, `*=`, `/=`
- ✅ `dot()`, `norm()`, `norm_squared()`

需要添加的功能：
- [ ] **一元负号**：`Vector operator-() const` → `-v`
- [ ] **Hadamard 乘积**：`Vector hadamard(const Vector& other) const` → `z = x .* y` (逐元素乘法)
- [ ] **从 std::vector 隐式/显式构造**：已有 `Vector(const std::vector<Real>&)`

### Phase 2: 重构 CG 求解器

**文件：`src/math/cg.cpp`**

1. 删除 `static Real dot()` 辅助函数
2. 将所有 `std::vector<Real>` 改为 `Vector`
3. 用运算符替换循环：
   - `x[i] += alpha * p[i]` → `x += alpha * p`
   - `r[i] -= alpha * Ap[i]` → `r -= alpha * Ap`
   - `p[i] = r[i] + beta * p[i]` → `p = r + beta * p`
4. 用 `Vector::dot()` 替换 `dot(a, b)`
5. 用 `Vector::norm_squared()` 替换手写点积

### Phase 3: 重构 PCG 求解器

**文件：`src/math/pcg.cpp`**

1. 删除 `static Real dot()` 和 `static void axpy()` 辅助函数
2. PCG 主循环用 Vector 重写（类似 CG）
3. 预条件器 `apply()` 保持 `std::vector<Real>` 接口（性能考虑）
4. 内部实现可以临时转换为 Vector

### Phase 4: 重构 Newton-Raphson 求解器

**文件：`src/math/newton_raphson.cpp`**

1. `compute_norm()` 用 `Vector::norm()` 替换
2. 向量取负：`neg_R = -1.0 * Vector(R)` 或添加一元负号
3. 更新解：`u_vec += alpha * du_vec`
4. 线搜索：`u_new = u_vec + alpha * du_vec`

### Phase 5: 更新预条件器（可选）

**文件：`src/math/pcg.cpp`**

预条件器的 `apply()` 接口暂时保持 `std::vector<Real>`（避免大量拷贝），内部实现可以用 Vector 包装。

### Phase 6: 检查其他文件

- `src/assembly/assembler.cpp`
- `src/physics/heat.cpp`
- `src/physics/elasticity_unified.cpp`
- `tests/` 中的测试文件

## 性能考虑

### 优点
- **编译器优化**：现代编译器会内联运算符，生成相同的机器码
- **SIMD 友好**：连续内存访问模式，易于向量化
- **避免错误**：减少手写循环的越界、索引错误

### 注意事项
- **避免临时对象**：用 `+=` 而非 `x = x + y`
- **RVO/NRVO**：返回值优化消除拷贝
- **接口设计**：高频调用的函数（如预条件器）可保留原始接口

## 验证

1. **单元测试**：所有测试通过（178/178）
2. **性能测试**：benchmark_preconditioners 结果保持一致
3. **数值精度**：solver 收敛迭代次数不变

## 时间估计

- Phase 1: 30 分钟（添加 Vector 功能）
- Phase 2: 30 分钟（重构 CG）
- Phase 3: 45 分钟（重构 PCG）
- Phase 4: 30 分钟（重构 Newton-Raphson）
- Phase 5: 15 分钟（可选）
- Phase 6: 30 分钟（检查其他文件）
- **总计：2.5-3 小时**

## 示例对比

### Before (手写循环)
```cpp
// CG 求解器
for (std::size_t i = 0; i < n; ++i) {
    x[i] += alpha * p[i];
    r[i] -= alpha * Ap[i];
}
Real rr_new = 0.0;
for (std::size_t i = 0; i < n; ++i) {
    rr_new += r[i] * r[i];
}
```

### After (Vector 运算符)
```cpp
// CG 求解器
x += alpha * p;
r -= alpha * Ap;
Real rr_new = r.norm_squared();
```

**代码行数减少 ~60%，可读性提升显著！** 🚀
