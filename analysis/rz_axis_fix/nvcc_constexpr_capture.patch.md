# nvcc constexpr-if 捕获修复(需回补到本地工作树)

**症状**(只有 nvcc 报,CPU 编译不报):
```
error: An extended __device__ lambda cannot first-capture variable
       in constexpr-if context
```
**根因**:变量在 device lambda **外**声明,却**首次**在 `if constexpr (...)` 分支里使用。
**修法**:在 lambda 顶部(任何 `if constexpr` 之前)用 `amrex::ignore_unused(...)` 强制预捕获。

---

## 1) `src/rhs/Weno.h`
锚点(lambda 内已有的一行):
```cpp
            amrex::ignore_unused(ibm_markers, face_flux);
```
在其**后**插入:
```cpp
            // nvcc: force-capture before any if-constexpr branch
            amrex::ignore_unused(pressure_reconstruction_states,
                                 pressure_reconstruction_component,
                                 radial_domain_high_index);
```

## 2) `src/rhs/Afd.h`
锚点(lambda 开头):
```cpp
            const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
            const bool trace_this_face =
```
改成:
```cpp
            const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
            // nvcc: force-capture before any if-constexpr branch
            amrex::ignore_unused(radial_domain_high_index,
                                 radial_pressure_face_flux);
            const bool trace_this_face =
```

---

**注**:0807c 那轮我只补在 AWS 副本上,本地未回补,所以 0807d 上传后同样的错误又出现了一次。
建议把这两处直接改进本地 `src/`,否则每次上 GPU 都要重补。
另一种等价写法:把相关 `constexpr` 常量**声明在 lambda 内部**。
