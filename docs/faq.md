# 物理问题 FAQ

[`../README.md`](../README.md) 讲的是"命令怎么敲、YAML 怎么写"。这份讲的是
"为什么这么写"，也就是用 vc-NEB 算相变能垒时绕不开的那几个物理问题。

目录

1. [vc-NEB 的原理](#1-vc-neb-的原理)
2. [超胞与公度](#2-超胞与公度)
3. [k 点怎么选](#3-k-点怎么选)
4. [压力怎么加](#4-压力怎么加)
5. [DFT+U、UKS、多重度、展宽](#5-dftuuks多重度展宽)
6. [磁性在路径上变了怎么办](#6-磁性在路径上变了怎么办)
7. [初末态已经优化过了，还要再优化吗](#7-初末态已经优化过了还要再优化吗)
8. [共格能垒 ≠ 成核能垒](#8-共格能垒--成核能垒)
9. [运行环境与 ASE](#9-运行环境与-ase)

---

## 1. vc-NEB 的原理

### 1.1 为什么不能只动原子

结构相变里晶胞也要变。如果做固定晶胞的 NEB，你必须先决定"用哪个晶胞"：
用初态的？末态的？中间某个？无论选哪个，另一端的结构就被压进了错误的晶胞，
能量里混进一大块弹性应变能，路径也会被拉歪。所以变胞 NEB 把**晶胞和原子
一起**当自由度来优化，两端点各自在自己该有的体积上弛豫（有外压时按焓）。

### 1.2 联合坐标

每张图的状态写成一个长向量

```
y = [ w·cell(9 个分量) | pos(3N 个分量) ]
```

`cell` 是 3×3 晶胞矩阵按行拉直，`pos` 是笛卡尔坐标。权重 `w` 来自
`neb.cell_weight`，默认取 `max(1, max_step / cell_max_step)`（默认
`max_step: 0.2 Å`、`cell_max_step: 0.02 Å` → `w = 10`）：意思是
"1 Å 的晶胞变化，在路径度量里等价于 10 Å 的原子位移"。

这个 `w` 的唯一作用是让两把尺子对齐：优化器在 `y` 空间里每步最多走
`max_step`，于是原子最多动 `max_step` Å，晶胞分量最多动
`max_step / w = cell_max_step` Å。想让晶胞动得更自由，就增大 `cell_max_step`
（`w` 会自动变小）。

### 1.3 NEB 力怎么分解

对中间的图 `i`，先取路径切向（也在 `y` 空间，所以含晶胞分量）

```
τ   = (y[i+1] − y[i-1]) / |y[i+1] − y[i-1]|
F_i = −∇E(y[i])               （真实力，晶胞部分的共轭力是 F_cell/w）
f⊥  = F_i − (F_i·τ) τ          （垂直路径，负责往路径上"贴")
f∥  = (F_i·τ) τ
F_neb = f⊥ + k·(|y[i+1]−y[i]| − |y[i]−y[i-1]|)·τ
```

`k = neb.spring_k`。弹簧力只改**路径上图的分布**和收敛速度，不改收敛后的
势垒（收敛时弹簧力沿路径的分量被平衡掉）。

### 1.4 攀爬像（climbing image）

* 选哪张图爬：取中间图里能量最高的那张；之后若别的图比它高出一个
  滞后阈值 `tol = max(1e-5, 1e-3·|E_cand − E_0|)` 才切换（避免每步跳来跳去）；
* 攀爬像的力：`F = f⊥ − f∥`（把沿路径的分量**反号**），于是它沿路径往
  上爬到一阶鞍点，同时仍被约束在路径附近；
* 日志每步都会写 `CI=n`；`summary.json` 的 `image_max` 和
  `history.jsonl` 的 `climbing_image` 记的就是它。

### 1.5 晶胞上的力

```
F_cell = −V · σ · cell^(−T)
```

`σ` 是引擎给的应力（有外压时已加 `P·I`）。这是"固定分数坐标"约定下
能量对晶胞矩阵的导数：晶胞随应力松弛，原子同时用自己的笛卡尔力松弛。

**不要**再补一项 `−fracᵀ·F_atom`。在自由原子的情况下那一项的量级和原子力
一样大，会让晶胞去跟随原子而不是跟随应力；`vcneb/band.py` 的注释里记了实测：
在 EMT 的 fcc→hcp 测试上加了这一项，晶胞被剪到 ~10°，势垒高了 4 倍。
`tests/test_cellforce.py` 用有限差分（整体缩放晶胞）验证的就是这条公式。

### 1.6 什么时候算收敛

```
max|F_atom| < neb.f_th        并且        max|F_cell| / n_atoms < neb.f_th
```

晶胞力按原子数归一化，否则它和原子力不是一个量纲（`F_cell` 是"每个 Å 的
能量"，`F_atom` 是"每个原子的力"）。

---

## 2. 超胞与公度

### 2.1 为什么必须同一个超胞

插值时晶胞是**线性插值**的（`(1−t)·cell_A + t·cell_B`），原子也是按序对应
线性插值。这就要求两件事：

1. 两个端点的晶胞**形状相近**（否则中间图是畸变的晶胞，能量没有意义）；
2. 两个端点的**原子数相同、顺序对应**（否则你会把 A 相的 Bi1 和 B 相的 Bi3
   连成一条线，路径穿过原子）。

所以：把两个相各自扩成同一个超胞，是这套方法的**前提**，不是可选项。

### 2.2 怎么选超胞矩阵

设两个相的常规胞（conventional cell，不是原胞）分别有 `Z_A`、`Z_B` 个原子，
选整数矩阵 `M_A`、`M_B`，要求

```
det(M_A) · Z_A  =  det(M_B) · Z_B        （原子数相同）
M_B⁻¹ · M_A    ≈  I                       （晶胞形状相近，应变小）
```

第二个条件决定了你能把应变压到多小。Bi 的例子：

| 相 | 空间群 | 常规胞原子数 | V (Å³) | V/atom |
|----|--------|--------------|--------|--------|
| Bi A7 | R-3m | 6 | 219.52 | 36.59 |
| Bi-II | C2/m | 4 | 139.23 | 34.81 |

`6·d₀ = 4·d₁` 的最小解是 `d₀ = 2, d₁ = 3` → 12 原子。于是
`M_A`（det 2，从 A7 常规胞）和 `M_B`（det 3，从 Bi-II 常规胞）都给出 12 原子。
这就是 `docs/example-bi.md` 里那对矩阵的来源。

### 2.3 为什么还要 Niggli 约化 + 刚体对齐

整数超胞矩阵几乎总会给出一个"合法的丑基"：三个矢量很长而且夹角接近平行
（实测有的组合夹角只有 19°）。平面波/GPW 代码在这种基上要么极慢要么直接
报警。所以 `tools/build_common_cell.py` 做四件事：

1. **spglib 标准化**：把输入（可能是原胞、可能是 2×1×1）变成常规胞；
2. **整数超胞 + Niggli 约化**：保持晶格和原子，换一个规整的基；
3. **刚体旋转**（SVD/Kabsch）：把 B 的晶胞整体转到 A 的取向上，这样应变
   矩阵 `F = cell_B·cell_A⁻¹` 才接近单位矩阵；
4. **原子匹配**：在最近像（MIC）距离下用 `linear_sum_assignment` 求
   一一对应，并对原点平移做搜索，取总位移最小的解。

脚本会打印 `residual strain`（剩下来的、必须由相变本身吸收的应变）和
每个原子的位移：

```
residual strain B->A: 4.76 %   (isotropic part +1.67 %)
matching: total |dr| = 16.615 A  (mean 1.385 A)
  worst pair: 2.40 A
```

### 2.4 代价

原子数是最小公倍数，代价通常是"原子数 × 图数"次 SCF。Bi 这个例子里
12 原子 × 9 图。如果应变要求你上 24 原子，成本就翻倍——先用
`residual strain` 比较几个候选矩阵，再决定。

### 2.5 实在做不到公度怎么办

两个相没法用同一个超胞描述（原子数怎么都凑不上），那 vc-NEB 给不出严格
的变胞能垒。可选的路：

* 找一个"母相超胞 + 切变"描述（很多马氏体相变可以）；
* 退一步做固定晶胞的 NEB（同一晶胞、只是原子重排），明确说明这是近似，
  能垒里不含晶胞形变的贡献；
* 用别的办法（变胞的 CI-NEB 加约束、或者直接算两相的焓差 + 相界面能）。

---

## 3. k 点怎么选

### 3.1 为什么整条路径只用一个网格

变胞时每个图的晶胞体积不同。如果每张图按自己的晶胞重新生成 k 点，那么
某张图可能 8×5×3、下一张 9×5×3——**网格一变，能量就有台阶**，这个台阶会被
当成力，最后混进势垒里。所以 vcneb 在开始之前就把网格定死，全程共用。

代价：体积变化大时，某些图上的网格偏粗或偏密。所以网格要按"路径上最需要
密网格的那个晶胞"来选。

### 3.2 具体是怎么算出来的

```
grid_i = max(1, int(|b_i| / spacing) + 1)         # |b_i| = 第 i 个倒格矢长度（不含 2π）
```

正交晶胞里 `|b_i| = 1/a_i`，于是就是"`1/(a_i·spacing)` 向上取整"；单斜/三斜
里要用真正的倒格矢长度，心算别硬套 `1/a`。

`spacing` 就是 `calc.kpts: {density: X}` 里的 `X`，单位 ≈ 1/Å，含义是倒空间
网格间隔的上限。两个端点的晶胞逐轴取**最大的那个 grid**，然后冻结。

Bi 例子（`density: 0.03`，a,b,c = 4.538 / 7.860 / 12.308 Å）：

```
1/(0.03×4.538) + 1 = 8.35 → 8
1/(0.03×7.860) + 1 = 5.24 → 5
1/(0.03×12.308)+ 1 = 3.71 → 3

k-points: [8, 5, 3] (frozen from kspacing 0.03)   (same grid for all images)
```

这条日志每次 `--dry-run` 都会打印，**换胞以后先看这一行**。

### 3.3 大晶胞 → 小晶胞，k 点怎么换

k 点数与晶格长度成反比（`k_i · a_i ≈ 常数`）。所以：

| 你做的事 | k 点怎么变 |
|----------|-----------|
| 单胞 → 2×2×2 超胞 | 每个方向大约减半 |
| 12 原子超胞 → 96 原子超胞 | 8×5×3 量级 → 4×2×1 量级，甚至只留 Γ |
| 原胞 → 常规胞 | 按新的三个晶格矢量长度重算，不是简单乘 |

不想手算就用 `{density: X}`：它会按两端点的实际晶格自动算，你只管挑 `X`。
经验值：金属 `0.02–0.03`，半导体/绝缘体 `0.04–0.06`（配合更严格的 `scf_thr`）。

### 3.4 变胞/加压时要注意什么

* 变胞路径的网格由**两个端点**决定（逐轴取大）。路径中间如果体积掉到比
  两端都小很多，网格相对变密（安全方向），日志会 WARN
  `frozen k-mesh ... is too coarse for this cell`——那说明晶胞跑出了端点
  覆盖的范围，要么加密网格，要么检查路径是不是失控了；
* **加压**让晶胞变小，所以用 `{density}` 算出来的网格在压力下天然偏密一点点，
  不会偏粗。反过来说：0 GPa 下定的网格不能直接拿去算 20 GPa（晶胞小了 10 %，
  网格该密 10 %），所以**不要**把 0 压的 `config.resolved.yaml` 直接改成
  `pressure: 20` 复用，重新 `--dry-run` 看网格；
* 一旦开始跑就不能换网格（换了就是另一条曲线），要换就 `--force`。

### 3.5 网格够不够

做一次 k 收敛测试：`density` 从 0.03 缩到 0.02（网格更密），势垒变化
< 10 meV 就可以停了。对金属小胞，`8×8×8` 起步是常识；对 12 原子以上的
超胞，`8×5×3` 这种量级通常够。

---

## 4. 压力怎么加

### 4.1 工具做的是什么

`neb.pressure: 2.55`（GPa）之后，整条带的能量换成焓：

```
H = E + P·V
dH/dh = −V·(σ + P·I)·h^(−T)          h = 晶胞矩阵
```

具体表现：

* 端点弛豫是**按焓**弛豫的（等于把晶胞压到该压力下的平衡体积）；
* 每张图的能量、晶胞力、`summary.json`、`barrier.png` 全是焓；
* 势垒 = 焓差，天然含 `P·ΔV`。

### 4.2 三条必须遵守的

1. **DFT 输入里的压力留 0**：QE 的 `press`、CP2K 的 `PRESSURE`、VASP 的
   `PSTRESS`、ABACUS 的 `press*` 全都别设。工具自己写的是 0 压输入，它
   在能量层加 `P·V`。两边都加就是加了两次。
2. `pressure` 只对 `variable_cell: true` 有意义。固定晶胞时 `P·V` 是常数，
   能垒一点都不变（只是能量零点平移）。
3. 只支持**静水压**（标量 `P`）。各向异性应力要自己改
   `EnthalpyCalculator`（`vcneb/band.py`，约 30 行）。

### 4.3 "文献说压力改变了相变"，该怎么算

正确做法是**在目标压力下算共格能垒**，而不是"0 压路径 + 事后加个 PV"：

* 0 压下的最小能量路径，在 2.55 GPa 下一般不是最小焓路径（原子配位、
  晶胞形状都会变）；
* 工具会把端点弛豫到该压力的平衡体积，路径由两个"该压力下的端点"插出来，
  这才对得上文献里"相变压力 P"的物理含义。

建议的做法（也是审稿人会问的）：

1. 先确认两相在目标压力下都能稳定存在（各自做变胞弛豫，看焓-体积曲线）；
2. 在该压力下做 vc-NEB，报 `ΔH‡ = ΔE + P·ΔV`，并**分开**报 `ΔE` 和 `P·ΔV`
   两个分量（相变体积差大时，`P·ΔV` 往往是主角）；
3. 扫一遍压力（例如 1.5 / 2.0 / 2.55 / 3.0 GPa），看能垒随压力的变化：
   能垒趋近 0 的那个压力就是这套共格路径的力学失稳点，通常和实验相变压力
   在同一量级；
4. 检查路径上的体积是单调/连续变化的。若某张图体积跳变，说明路径不物理。

### 4.4 一个可以复现的数字例子

内置的 LJ Ar 例子（fcc→hcp，4 原子，`examples/lj_ar.yaml`，改
`neb.pressure` 就能复现）：

| `neb.pressure` | 端点体积 | 势垒 |
|----------------|----------|------|
| 0 GPa | 145.99 Å³ | 0.0142 eV |
| 5 GPa | 100.91 Å³ | 0.0700 eV |

注意这是**验证"压力确实进了焓"的玩具计算**：LJ 氩的体模量很小，5 GPa 压掉了
30 % 体积，所以势垒涨了 5 倍，不要当物理结果引用。

---

## 5. DFT+U、UKS、多重度、展宽

### 5.1 CP2K

```yaml
calc:
  engine: cp2k
  params:
    uks: true                  # 自旋极化；磁性 + DFT+U 必须开
    multiplicity: 3            # 2S+1；开 smearing 时只是初猜
    hubbard_u: {Ni: 4.0}       # U_eff = U − J，单位 eV，按元素给
    hubbard_l: 2               # 2 = d，3 = f
    plus_u_method: mulliken    # mulliken / lowdin / mulliken_charges
    mulliken: true             # 打印 Mulliken 布居（默认就开）
    smearing: true
    smearing_method: marzari_vanderbilt
    smearing_sigma: 0.05       # eV（写进 SIGMA [eV]）
```

生成的输入是每个 `&KIND` 一个 `&DFT_PLUS_U` 块（Dudarev 形式，`U_MINUS_J`
+ `L`）。`doctor`/`keys` 会替你检查几件容易忘的事：

* `multiplicity` 没开 `uks` → 直接被忽略（变成闭壳层）；
* 开了 `uks` 但没 `added_mos` 就用 `relax_multiplicity` → 无效；
* 闭壳层 + DFT+U → 有能量修正但没有自旋，分不出 FM/AFM；
* DFT+U 不开展宽 → SCF 容易掉到亚稳占据。

### 5.2 QE

```yaml
calc:
  engine: qe
  params:
    nspin: 2
    starting_magnetization: {Ni: 0.5}
    lda_plus_u: true
    Hubbard_U: {Ni: 4.0}          # eV
    U_projection_type: ortho-atomic
    occupations: smearing
    smearing: marzari-vanderbilt
    degauss: 0.02                 # Ry
```

QE 的 `params` 只认 ASE 登记过的关键词，**拼错会静默丢失**：写之前
`vcneb keys qe | grep -i -E 'mag|hubbard|smear'` 对一遍，写之后打开
`run/images/02/espresso.pwi` 看实际内容。

### 5.3 ABACUS

```yaml
calc:
  engine: abacus
  params:
    nspin: 2
    dft_plus_u: 1
    orbital_corr: [2, -1]        # 每个元素一项：2 = d，3 = f，-1 = 不加 U
    hubbard_u: [4.0, 0.0]        # eV，顺序与 STRU 里的元素顺序一致
    smearing_method: mp          # mp / mv / fd / fixed
    smearing_sigma: 0.01         # Ry
    out_mul: 1                   # 打印 Mulliken（可选）
```

ABACUS 的 `params` 基本原样进 `INPUT`（ASE 的 calculator 负责写），
所以 ABACUS 手册里的关键词都能用；vcneb 只额外代写赝势/轨道/`KPT` 和几个
安全默认值（`basis_type: lcao`、`ks_solver: genelpa`、`symmetry: 0`、
`cal_force: 1`、`cal_stress: 1`、`scf_thr: 1e-7`、`scf_nmax: 100`）。

### 5.4 smearing：fermi 还是高斯？

先说结论（金属相变能垒的常见选择）：

| 体系 | 推荐 | 说明 |
|------|------|------|
| 金属（Bi、Sn、过渡金属…） | MP（`methfessel_paxton`）或 MV（`marzari_vanderbilt`） | 对总能"不敏感"，力平滑，适合跑 NEB |
| 半金属/小带隙 | Fermi-Dirac（`fermi`） | 物理上就是电子温度；展宽当温度用 |
| 绝缘体/分子 | GAUSSIAN（小 `sigma`）或 `occupations: fixed` | 只为帮 SCF 收敛；`sigma` 大了会抹平能垒 |
| 收敛困难 | 先用大 `sigma` 跑到收敛，再逐步减小 | 直接上小 `sigma` 常常抖 |

三个引擎的写法差异：

| | 方法名 | 宽度关键字 | 单位 |
|--|--------|-----------|------|
| QE | `smearing: mp / mv / fd / gaussian` | `degauss` | Ry |
| CP2K | `smearing_method: mp / mv / fermi / gaussian / energy_window` | `smearing_sigma`（MP/MV/GAUSSIAN）或 `electronic_temperature`（FERMI_DIRAC） | **eV** / K |
| ABACUS | `smearing_method: mp / mv / fd / fixed` | `smearing_sigma` | Ry |

CP2K 的 `SIGMA` 在 vcneb 里写的是**eV**（写手会自动带 `[eV]` 单位），
而 CP2K 原生默认是 Hartree——这就是为什么同一套体系在 QE/CP2K 之间搬参数时
数值看起来差很多。`gauss`/`gaussian` 从 CP2K 2026.2 开始支持（`METHOD
GAUSSIAN`），vcneb 的映射表里已经有了；旧版本会直接报"未知方法"。

选多大的宽度：金属 QE `degauss` 0.01–0.03 Ry / CP2K `smearing_sigma`
0.05–0.1 eV；并且**至少检查一次**：`sigma` 加倍，势垒变化应该在 10 meV 量级，
如果变化很大，说明这个量你要在文章里报出来。

### 5.5 为什么 vcneb 只解析 CP2K 的 Mulliken

`barrier.csv` 里的 `total_spin` / `net_charge` 两列只有在引擎输出里能读到
布居时才有。CP2K 的写手会自动加 `&DFT/&PRINT/&MULLIKEN`，解析器读那张表；
QE 和 ABACUS 的 ASE calculator 不解析布居，所以那两列不会出现。要看自旋，
去 `run/images/NN/` 里自己 grep：

```bash
grep -i "total magnetization" run/images/*/espresso.pwo      # QE
grep -i "total magnetism"     run/images/*/OUT.*.ABACUS/running_scf.log  # ABACUS
```

---

## 6. 磁性在路径上变了怎么办

这是相变能垒里最容易出错的地方：**如果两个端点的磁态不一样（FM → 非磁、
FM → AFM、高自旋 → 低自旋），那么"能垒"里可能有一大块是"占据约束"造成的，
不是几何造成的。**

### 6.1 先判断有没有发生

* `barrier.csv` 有 `total_spin` 列时，看它沿路径怎么变；
* 日志里会有一条逐图自旋，以及一条 WARNING（变化 > 0.05 μB 时）：

```
WARNING the total spin moment changes along the path. With a magnetization
forced by the input (CP2K MULTIPLICITY without smearing, QE tot_magnetization,
ABACUS nupdown) part of the barrier can be an occupation constraint rather
than the geometry. Compare the Mulliken table of the endpoint and the climbing
image.
```

### 6.2 三种做法（按推荐顺序）

**A. 让自旋自己浮动（最常用）**

打开展宽，别钉多重度：

```yaml
    uks: true
    smearing: true
    smearing_method: marzari_vanderbilt
    smearing_sigma: 0.05
    added_mos: 12          # 空带，展宽才转得动
```

CP2K 里 `&SMEAR` 之后 α/β 共用一个费米能级，最终自旋由能量决定（实测：
`multiplicity: 3` + 展宽的 Ni 体系最后得了 M = 0.81，而不是 2）。这是最接近
"体系自己决定磁态"的做法。代价是自旋可能落在两个极小之间来回跳，看
`history.jsonl` 的能量曲线有没有台阶。

**B. 只给初猜，让 CP2K 沿路径放开多重度**

```yaml
    uks: true
    multiplicity: 3
    relax_multiplicity: 0.00001   # 需要 added_mos > 0
```

它让 α/β 电子数在 SCF 内部浮动。比 A 更"约束在磁态附近"，但每一步都会多花
SCF 时间。

**C. 分开算两条路径（磁态差别很大时）**

FM 和 AFM 用不同的 `starting_magnetization`（QE）或不同的初始 `MAGMOM`
（CP2K 用 `&KIND &BS`，走 `params.extra`），各跑一条 vc-NEB，比较两条能垒。
这是最老实、也最容易被审稿人接受的做法：报"从 FM 态到非磁态，能垒 X eV；
从 AFM 态到非磁态，能垒 Y eV"。

### 6.3 DFT+U 时的额外注意

* 必须 `uks: true`（否则 +U 只改能量，没有自旋）；
* `U` 值必须在整条路径上固定（工具本来就固定，别在续算时改，见 README §11.4）；
* `plus_u_method` 固定用 `mulliken`（默认），换 LOWDIN 会得到不同的 U 投影；
* +U 的 SCF 更硬：`cutoff`/`rel_cutoff`（CP2K）或 `ecutwfc`（QE/ABACUS）
  要重新收敛一次，`eps_scf` 收紧到 `1e-7`。

---

## 7. 初末态已经优化过了，还要再优化吗

`neb.relax_endpoints: true` 是默认：两条端点会各自再弛豫到
`neb.endpoint_fmax`（默认 `0.01 eV/Å`）。**大多数情况下应该保留**，因为：

* 端点如果不是极小点，它自己的能量就偏高/偏低，能垒直接带误差；
* 换成公度超胞之后，端点的应变状态已经变了（Bi 例子里残差应变 4.76 %），
  必须重新弛豫；
* 有 `pressure` 时，端点要按**焓**弛豫（0 压下的优化结构在 2.55 GPa 下不是
  极小点）。

什么情况下可以关掉：

```yaml
neb:
  relax_endpoints: false
```

* 端点确实已经是你用的这套参数（同一赝势、同一截断、同一 `k` 网格、同一展宽、
  同一 DFT+U）下的极小点，而且是**同一个超胞**；
* 你想要"固定端点"的可复现性（比如对比不同路径时端点必须一模一样）。

判断端点干不干净就看日志：

```
WARNING endpoint 'initial' did not reach 0.05 eV/A in 300 steps;
the barrier will be less accurate.  Increase neb.endpoint_steps.
```

出现这条就说明 `relax_endpoints: false` 会更糟。

---

## 8. 共格能垒 ≠ 成核能垒

vc-NEB 的两个端点必须能塞进同一个晶胞，所以它算的是**共格（coherent）
转变**：整个晶体同时变形，没有界面。真实的相变一般从界面/缺陷/表面开始，
先形成临界核，成核能垒里包含界面能，量级可以差一个数量级（共格的通常更大，
因为要把整块晶体的应变能都付掉）。

写文章时的三句话：

1. 报的数是"coherent, uniaxial/stress-free …"的**共格转变能垒**，不是成核能垒；
2. 如果实验/文献报的是成核能垒或相变动力学温度，比较时要说清差异；
3. vc-NEB 给的是 0 K 的势能（或焓）差，**不含零点能和熵**。要自由能垒，
   拿鞍点图（`band.traj` 里 `image_max` 那张）做频率计算再加 ZPE/熵修正。

---

## 9. 运行环境与 ASE

### 9.1 底层是 ASE 吗

是。vcneb 用 ASE 读结构（POSCAR/cif/xyz/traj…）、写 `Trajectory`
（`band.traj`）、做 IDPP 插值（`ase.mep.idpp_interpolate`）和最优化
（FIRE）。所以 `ase` 是硬依赖，装 `vcneb` 时会自动带上。

### 9.2 什么环境能跑

| 要求 | 说明 |
|------|------|
| Python ≥ 3.9 | 只用标准库 + numpy + ase + pyyaml（画图另加 matplotlib） |
| 操作系统 | Linux / WSL / macOS / Windows 都行（引擎可执行文件能在 PATH 里就行） |
| 并行 | 只通过 `calc.command` 调外部程序（`mpirun`、`srun`、`sbatch` 包装脚本、容器都行）；`neb.parallel_images` 可以多图并行 |
| 队列系统 | 把 `command` 写成 `srun ...`，或写一个提交脚本；登录节点没有可执行文件时用 `vcneb doctor` 会报错，正式跑加 `--skip-engine-check` |
| GPU | 不支持（也不需要）：一切靠文件接口 |

### 9.3 一个典型的集群用法

```yaml
calc:
  engine: abacus
  command: srun -n 64 abacus
  mpi: 64
```

```bash
# 交互节点上检查
vcneb doctor neb.yaml
vcneb run neb.yaml --dry-run
# 提交
sbatch -N 1 -n 64 -t 24:00:00 job.sh     # job.sh 里就一行 vcneb run neb.yaml --restart
```

每次排队时间到、任务被杀，就再交一次 `--restart`（进度在 `state.json` 里，
见 README §11）。
