# vcneb — 傻瓜式变胞 NEB（QE / CP2K / ABACUS）

> GitHub 仓库名 **VCTST**；命令行和包名仍是 `vcneb`（`pip install -e .` 之后敲 `vcneb`）。

一个 YAML 文件进去，一条命令出来：**端点弛豫 → 插值（IDPP）→ 变胞 NEB
（含攀爬像）→ 势垒图**。不用写 Python，所有开关都在同一个配置文件里。

```bash
vcneb init -e cp2k -o neb.yaml     # 1. 生成带注释的模板
vcneb doctor neb.yaml              # 2. 检查结构/路径/赝势/程序能不能跑（不花机时）
vcneb run neb.yaml --dry-run       # 3. 只建路径，导出两个端点，先看一眼
vcneb run neb.yaml                 # 4. 正式跑：summary.json / barrier.csv / barrier.png
vcneb status ./run                 # 5. 随时看进度、看要不要续算
```

| engine   | 说明 | 需要自己准备 |
|----------|------|--------------|
| `qe`     | Quantum ESPRESSO（pw.x），平面波 | UPF 伪势 + `pw.x` 可执行文件 |
| `cp2k`   | CP2K（cp2k.popt/psmp），GPW | CP2K `data` 目录（BASIS_MOLOPT、GTH_POTENTIALS） |
| `abacus` | ABACUS，LCAO 或 PW | UPF 伪势（LCAO 还要 `.orb` 轨道文件） |
| `emt`    | ASE 内置 EMT（只支持少数金属，秒级） | 什么都不用（自检/教学） |
| `lj`     | Lennard-Jones（Ar 的 fcc→hcp 经典例子） | 什么都不用（秒级） |

目录

1. [它做什么、不做什么](#1-它做什么不做什么)
2. [安装](#2-安装)
3. [30 秒自检（不需要任何 DFT 程序）](#3-30-秒自检不需要任何-dft-程序)
4. [命令一览](#4-命令一览)
5. [逐条命令详解](#5-逐条命令详解)
6. [YAML 字段全表](#6-yaml-字段全表)
7. [三种 DFT 引擎怎么配](#7-三种-dft-引擎怎么配)
8. [路径是怎么插出来的（插值四步）](#8-路径是怎么插出来的插值四步)
9. [删点：min_distance / prune](#9-删点min_distance--prune)
10. [加压力：neb.pressure](#10-加压力nebpressure)
11. [续算与重来](#11-续算与重来)
12. [输出文件](#12-输出文件)
13. [常见错误](#13-常见错误)
14. [限制与注意](#14-限制与注意)
15. [例子与测试](#15-例子与测试)

原理、超胞/公度、k 点、DFT+U、磁性、smearing 这些"物理问题"另见
[`docs/faq.md`](docs/faq.md)；一个完整的真实算例（Bi 的 A7→Bi-II，2.55 GPa）
见 [`docs/example-bi.md`](docs/example-bi.md)；想知道每条命令最后落到哪个函数、
插值和删点具体改了哪几行，见 [`docs/internals.md`](docs/internals.md)。

---

## 1. 它做什么、不做什么

**做**（一条 `vcneb run` 全包）：

* 读两个端点结构（POSCAR/CONTCAR/cif/xyz/traj/…，ASE 能读的都行），自动把
  末态的原子顺序对齐到初态；
* 按给定精度弛豫两个端点（变胞时连晶胞一起弛豫；有外压时按焓弛豫）；
* 生成整条初始路径（线性 + IDPP 预优化 + 周期性展开 + 坏点剔除）；
* 跑变胞 NEB：原子和晶胞一起动，攀爬像（CI）打开以拿到真鞍点；
* 全程只用一个固定的 k 点网格（变胞时不重新生成 k 点，避免能垒里混进
  网格切换造成的假台阶）；
* 输出势垒、正/反向势垒、逐图能量/受力/体积/自旋，画图，写轨迹，断点续算。

**不做**：

* 不替你判断过渡态/路径在物理上是否合理——请用 `band.traj` 自己看；
* 不做截断能/k 点/展宽/赝势的收敛测试（模板里的值只保证"能跑通"）；
* 不改变"共格（coherent）路径"这个前提：给出的能垒是两相共格转变的能垒，
  不是成核能垒（见 `docs/faq.md` 第 8 节）。

---

## 2. 安装

```bash
git clone <repo> vcneb && cd vcneb
python -m pip install -e .            # 依赖 numpy / ase / pyyaml
python -m pip install -e ".[plot]"    # 想要 barrier.png 再加 matplotlib
vcneb --version                       # vcneb 0.1.0
```

不想装也行，直接把目录塞进 `PYTHONPATH`：

```bash
export PYTHONPATH=/path/to/vcneb:$PYTHONPATH   # PowerShell: $env:PYTHONPATH="E:\path\to\vcneb"
python -m vcneb --version
```

要求 Python ≥ 3.9。`spglib` / `scipy` 只有 `tools/` 里的辅助脚本用到，
主程序不需要。

---

## 3. 30 秒自检（不需要任何 DFT 程序）

```bash
vcneb init -e lj -o ar.yaml
vcneb doctor ar.yaml
vcneb run ar.yaml          # 2 秒，Ar 的 fcc→hcp，势垒 ~0.014 eV
vcneb status ./run_demo
```

真实输出（本机实测）：

```
barrier (forward) = 0.0142 eV at image 3
barrier (reverse) = 0.0143 eV
plot written: /tmp/doc_demo/run_demo/barrier.png

barrier   : 0.0142 eV
reverse   : 0.0143 eV
steps     : 27  converged: True
results   : /tmp/doc_demo/run_demo
```

`examples/emt_cu111.yaml`（Cu(111) 滑移，固定晶胞）也是几十秒就跑完。

---

## 4. 命令一览

```
vcneb [-h] [--version] {init,doctor,keys,run,plot,status}

  init   [-o OUTPUT] [-e {qe,cp2k,abacus,emt,lj}]   写模板（不花机时）
  doctor CONFIG                                     体检（不花机时）
  keys   [project|structure|calc|neb|output|qe|cp2k|abacus]
  run    CONFIG [--restart] [--force] [--dry-run]
                 [--max-steps N] [--skip-engine-check] [-v]
  status RUN_DIR                                    看进度/势垒
  plot   RUN_DIR [-o OUT]                           只重画图
```

| 命令 | 花 DFT 机时? | 作用 |
|------|--------------|------|
| `init` | 否 | 生成带注释的模板 YAML（按引擎不同） |
| `keys` | 否 | 打印所有字段/默认值，或某引擎支持的 `calc.params` 键 |
| `doctor` | 否 | 依赖、配置、结构、引擎（可执行文件、赝势、轨道、data 目录）全查一遍 |
| `run --dry-run` | 否 | 建路径、打印间距、导出 `endpoint_*.vasp`，不调用任何电子结构程序 |
| `run` | 是 | 端点弛豫 + 变胞 NEB |
| `status` | 否 | 读 `summary.json` / `state.json` 报进度和势垒 |
| `plot` | 否 | 用 `summary.json` 重画 `barrier.png` |

---

## 5. 逐条命令详解

### 5.1 `vcneb init`

```bash
vcneb init -e qe     -o neb_qe.yaml
vcneb init -e cp2k   -o neb_cp2k.yaml
vcneb init -e abacus -o neb_abacus.yaml
vcneb init                    # 默认 -e lj -o neb.yaml
```

* `-e/--engine`：`qe` / `cp2k` / `abacus` / `emt` / `lj`，注释跟着引擎变
  （CP2K 模板自带 `uks`/`hubbard_u`/`mulliken` 说明，QE 模板自带 smearing
  取舍说明）。
* `-o/--output`：目标文件。**文件已存在会直接报错**
  （`... already exists; delete it first or pick another output file`），
  不会覆盖你改过的配置。
* 模板开头有中文提示：实际只需改 `structure`、`calc.command`、赝势路径、
  `neb.images/f_th/max_steps`。

### 5.2 `vcneb keys`

```bash
vcneb keys            # 五段（project/structure/calc/neb/output）全部字段+默认值
vcneb keys neb        # 只看一段
vcneb keys qe         # QE: ASE 认识的 pw.x 全部关键词，按 namelist 分组，
                      #     并标出哪些由工具代写（<- default set by vcneb）
vcneb keys cp2k       # CP2K: 内置写手真正会渲染的 calc.params 键
vcneb keys abacus     # ABACUS: 无白名单，提示哪些键由工具代写
```

`vcneb keys neb` 的真实输出（可以拿来核对 yaml 有没有拼错）：

```yaml
[neb]
images: 5
interpolation: idpp
align_cell: true
tangent: improved
min_distance: 2.0
prune: auto
variable_cell: true
pressure: 0.0
spring_k: 0.1
f_th: 0.03
max_steps: 300
dt: 0.1
dt_max: 0.3
max_step: 0.2
cell_max_step: 0.02
cell_weight: null
climbing_image: true
relax_endpoints: true
endpoint_fmax: 0.02
endpoint_steps: 200
parallel_images: 1
keep_intermediate_images: true
```

两个进阶键不在这张表里（代码里按需读取，不给就用默认值）：
`neb.idpp_steps`（200）、`neb.idpp_fmax`（0.02 eV/Å）。

> 字段拼错时 YAML 不会报错、只会被忽略。不确定就跑一次 `vcneb keys`，
> 或者跑完看 `<workdir>/config.resolved.yaml`——那是**实际生效**的完整配置。

### 5.3 `vcneb doctor`

```bash
vcneb doctor neb.yaml
```

四段检查，任何一段 FAIL 都返回非零。

```
vcneb 0.1.0 doctor: ar.yaml

[dependencies]
  OK    numpy      2.5.3
  OK    ase        3.26.0b1
  OK    yaml       6.0.3
  OK    matplotlib 3.11.2

[config]
  OK    config looks consistent

[structures]
  OK    initial: 4 atoms, V=145.5316 A^3, formula=Ar4
  OK    final  : 4 atoms, V=145.7304 A^3, formula=Ar4
        initial: built from config ({'formula': 'Ar', 'crystal': 'fcc', 'a': 5.26, 'cubic': True})
        final: built from config ({'formula': 'Ar', 'crystal': 'hcp', 'a': 3.72, 'c': 6.08, 'orthorhombic': True})

[engine]
  Lennard-Jones (built-in, for workflow tests only)
  OK    engine configuration looks usable

all checks passed -> vcneb run ar.yaml
```

`[engine]` 段的 `WARN` 是**建议**不是错误，例如：

* CP2K 的 `calc.params` 里有写手不认识的键（会被静默忽略，提示改用 `params.extra`）；
* `multiplicity` 没配 `uks: true`（等于闭壳层，磁矩没有意义）；
* `basis_set`/`potential` 的名字在你的 `BASIS_MOLOPT`/`GTH_POTENTIALS` 里
  对该元素不存在（会顺便列出该元素可用的 PBE 势）；
* `smearing_method: mv` 但没给 `smearing_sigma`（CP2K 会用默认 0.002 Ha ≈ 0.054 eV）。

`doctor` 查不到的东西：运行时才缺的共享库（`libxxx.so`，用
`calc.env.LD_LIBRARY_PATH` 解决）、许可证、队列限额。

### 5.4 `vcneb run`

```bash
vcneb run neb.yaml                    # 正式跑
vcneb run neb.yaml --dry-run          # 只建路径（第一次一定先跑这个）
vcneb run neb.yaml --restart          # 断点续算
vcneb run neb.yaml --force            # 清掉本工具自己的输出后重跑
vcneb run neb.yaml --max-steps 800    # 临时改 neb.max_steps
vcneb run neb.yaml --skip-engine-check
vcneb run neb.yaml -v                 # 或 VCNEB_VERBOSE=1，出错时打印 traceback
```

五个开关的准确语义：

| 开关 | 做什么 | 不做什么 |
|------|--------|----------|
| `--dry-run` | 建路径、打印间距、写 `endpoint_initial.vasp`/`endpoint_final.vasp`；**端点不弛豫**（日志写明 "endpoints are used exactly as given"）；引擎没装好也只给 warning | 不写 `state.json`、不算能量，所以 dry-run 之后可以直接 `vcneb run` 正式跑 |
| `--restart` | 从 `<workdir>/state.json` 载入整条带（每图晶胞/坐标/能量、FIRE 步长、`pruned_images`）接着算；张数由存档决定 | 不会重新插值；`state.json` 不存在直接报错 |
| `--force` | 先删 `images/`、`logs/`、`state.json`、`summary.json`、`barrier.*`、`band.traj`、`history.jsonl`、`endpoint_*.vasp`、`config.resolved.yaml` 再重跑 | 只删本工具自己的输出，不删你放进 workdir 的其它文件（所以 `workdir` 别写 `.`） |
| `--max-steps N` | 覆盖 `neb.max_steps`（本次有效，不写回 yaml） | 不影响 `endpoint_steps` |
| `--skip-engine-check` | 引擎报问题时也照跑（例如集群上 `mpirun` 只在计算节点可见） | 不会让配置错误消失 |

**执行流程**（`cli.cmd_run`）：

1. `validate()`：引擎名、端点文件、`images/f_th/max_steps/interpolation`、
   `min_distance/prune`、`pressure`、`kpts` 合法性；有问题直接列出退出。
2. `engine.plan_kpoints([initial.cell, final.cell])`：**一次性**定下 k 网格，
   之后每张图都用它（日志 `k-points: [8, 5, 3]   (same grid for all images)`）。
3. `engine.check()`：可执行文件、赝势、轨道、data 目录。
4. `Band(...)` → `relax_endpoints()`（`neb.relax_endpoints: false` 可关）→
   `interpolate()` → `save_state()`。
5. NEB 主循环：每步 `evaluate_all()`（可 `parallel_images` 并发）→ NEB 力 →
   FIRE 走一步 → `_guard_geometry()` → `save_state()`。
6. 收敛后写 `summary.json` / `barrier.csv` / `band.traj` / `barrier.png` /
   `config.resolved.yaml`。

**返回码**：`0` 收敛；`2` 跑完但没收敛；`1` 配置/引擎/计算错误；
`130` 被 Ctrl-C 打断（进度已在 `state.json`）。

**日志逐行**（真实输出）：

```
vcneb 0.1.0  config: /tmp/p5/ar.yaml
engine: Lennard-Jones (built-in, for workflow tests only)
  project      : ar_fcc_to_hcp_demo
  workdir      : /tmp/p5/run_demo
  engine       : lj   command: -
  images       : 5 movable (+2 endpoints), interpolation=idpp
  variable cell: True   spring_k=0.1   f_th=0.002 eV/A
  max_steps    : 400   dt=0.1   climbing_image=True
  pressure     : 5.0 GPa (every energy below is the enthalpy H = E + P*V)
k-points: Gamma   (same grid for all images)
endpoints: 4 atoms, V0=145.5316 -> V1=145.7304 A^3
relaxing initial endpoint (fmax=0.002 eV/A, cell=free) ...
  initial: E=3.230105 eV  max|F|=0.0000 eV/A  converged  (2s)
relaxing final endpoint (fmax=0.002 eV/A, cell=free) ...
  final: E=3.230234 eV  max|F|=0.0000 eV/A  converged  (3s)
  IDPP pre-relaxation of the movable images done
initial path built: 7 images (5 movable), interpolation=idpp
smallest interatomic distance on the path: 3.229 A
path spacing (RMS atom displacement between neighbouring images): 0.095, 0.098, 0.101, 0.104, 0.107, 0.111 A
------------------------------------------------------------------------------
starting vc-NEB: 5 movable images, f_th=0.002 eV/A, max_steps=400
------------------------------------------------------------------------------
...
step   52  E[max]=3.30006  barrier=0.0700  max|F_atom|=0.0019  max|F_cell|/atom=0.0013  dt=0.150  (0s)  CI=3
converged: max|F_perp| = 0.0019 eV/A < 0.002
```

`pressure` 那行**只在 `neb.pressure ≠ 0` 时出现**（这里把 `vcneb init -e lj`
的产物加了 `neb.pressure: 5.0`，所以是 52 步 / 0.0700 eV；同一份配置
0 GPa 是 27 步 / 0.0142 eV，见 §3、§10）。

`step` 行含义：`E[max]` 是最高那张图的能量，`barrier` 是当前势垒，
`max|F_atom|` 是各图垂直力（去掉平行分量）的原子部分最大值，
`max|F_cell|/atom` 是晶胞部分的垂直力除以原子数，`dt` 是 FIRE 当前步长，
`CI` 是当前攀爬像编号。**两个判据同时满足才算收敛**：
`max|F_atom| < f_th` 且 `max|F_cell|/atom < f_th`。

### 5.5 `vcneb status`

```bash
vcneb status ./run
```

跑过的话读 `summary.json`：

```
run     : run_demo
engine  : lj
status  : converged after 27 steps
barrier : 0.0142 eV (reverse 0.0143 eV)

 image      E (eV)      dE (eV)    max|F| (eV/A)     V (A^3)
     0     -0.33022      0.00000         0.0000    145.9930
     1     -0.32720      0.00302         0.0006    147.7561
     2     -0.32128      0.00894         0.0016    148.9084
     3     -0.31601      0.01420         0.0016    149.3339
     4     -0.31726      0.01296         0.0001    148.9618
     5     -0.32459      0.00563         0.0001    147.7834
     6     -0.33035     -0.00013         0.0000    145.9052
```

没收敛时额外给一行可以直接复制的命令：

```
not converged yet; continue with
  vcneb run <config> --restart --max-steps 506
```

只有 `state.json`（被 Ctrl-C 打断、还没写 summary 时）：

```
run     : ./run
progress: step 12 of an unfinished run
resume  : vcneb run <config> --restart
```

两个都没有：`... does not look like a vcneb run directory`。

### 5.6 `vcneb plot`

```bash
vcneb plot ./run                 # 重画 ./run/barrier.png
vcneb plot ./run -o fig/band.pdf # 换输出路径/格式
```

只读 `summary.json`，不会再碰 DFT。左边是能量曲线（标出势垒和最高点），
有 Mulliken 布居时右边多一个自旋面板。

---

## 6. YAML 字段全表

### 6.1 `project`

| 键 | 默认 | 说明 |
|----|------|------|
| `name` | `vcneb` | 只用于日志/`summary.json` 显示 |
| `workdir` | `./run` | 所有输出放这里；相对路径按 **YAML 文件所在目录**解析。别写 `.`（`--force` 会在这下面删东西） |

### 6.2 `structure`

| 键 | 默认 | 说明 |
|----|------|------|
| `initial` / `final` | 无 | 端点结构文件（ASE 能读的任意格式）。原子数必须相同；同种元素内部的顺序会自动对齐（`linear_sum_assignment` + 最小像距离） |
| `build_initial` / `build_final` | 无 | 不用文件，直接给 ASE `bulk()` 的参数：`{formula, crystal, a, c, cubic, orthorhombic, supercell, vacuum, ...}` |
| `vacuum` | 无 | （仅 `build_*` 用）沿 z 加真空层 |

### 6.3 `calc`

| 键 | 默认 | 说明 |
|----|------|------|
| `engine` | 无 | `qe` / `cp2k` / `abacus` / `emt` / `lj` |
| `command` | 无 | 完整命令行，如 `mpirun -np 8 pw.x`；每张图都跑一次 |
| `mpi` | 1 | 只用于记账（`mpi × parallel_images` 超核数时提醒），不会自动加进 `command` |
| `env` | `{}` | 程序需要的环境变量（`LD_LIBRARY_PATH`、`OMP_NUM_THREADS`…） |
| `pseudo_dir` | 无 | UPF 所在目录（QE / ABACUS） |
| `data_dir` | 无 | CP2K 的 `data` 目录；ABACUS 里放 `.orb` |
| `timeout` | 3600 | 单次调用超时（秒） |
| `kpts` | 无 | `[4,4,4]` / `{gamma: true}` / `{density: 0.03}`（见 §6.6） |
| `pseudopotentials` | `{}` | 元素→文件名，如 `Cu: Cu_ONCV_PBE-1.0.upf` |
| `orbitals` | `{}` | 元素→`.orb` 文件名（ABACUS LCAO） |
| `params` | `{}` | 引擎关键字，透传（规则见 §7.4） |

### 6.4 `neb`

| 键 | 默认 | 说明 |
|----|------|------|
| `images` | 5 | **可动** image 数（不含两个端点，总图数 = images + 2） |
| `interpolation` | `idpp` | `idpp`（推荐）或 `linear` |
| `align_cell` | `true` | 插值前把末态晶胞换个写法对齐初态（§8.1） |
| `tangent` | `improved` | 路径切向：`improved`（Henkelman–Jónsson 能量加权）或 `central`（中心差分，v0.1.0 的旧行为）。对称/直线带上两者等价 |
| `min_distance` | 2.0 | 插值后某张图最近原子间距小于它就删掉（Å）；`0` = 关 |
| `prune` | `auto` | `auto` 删点继续 / `error` 直接报错 / `false` 只留警告 |
| `variable_cell` | `true` | 变胞 NEB；`false` 固定晶胞（slab、团簇建议关） |
| `pressure` | 0.0 | 外压（GPa），整条带用焓 `H = E + P·V`（§10） |
| `spring_k` | 0.1 | 弹簧常数（eV/Å²） |
| `f_th` | 0.03 | 收敛标准：垂直力 max\|F⊥\|（eV/Å），原子和晶胞两部分都要满足 |
| `max_steps` | 300 | 最大步数（`--max-steps` 可临时覆盖） |
| `dt` | 0.1 | FIRE 初始步长 |
| `dt_max` | 0.3 | FIRE 步长上限 |
| `max_step` | 0.2 | 单步原子最大位移（Å） |
| `cell_max_step` | 0.02 | 单步晶胞最大变化（Å） |
| `cell_weight` | `null` | 晶胞在联合坐标里的权重；`null` = `max_step / cell_max_step` |
| `climbing_image` | `true` | 打开才有真鞍点；用 1e-3·ΔE 的滞后避免最高点每步来回跳 |
| `relax_endpoints` | `true` | 先弛豫两个端点。端点已充分弛豫可以关，但**必须**保证它们在该精度下是极小点 |
| `endpoint_fmax` | 0.02 | 端点弛豫收敛标准（eV/Å） |
| `endpoint_steps` | 200 | 端点弛豫最大步数 |
| `parallel_images` | 1 | 同时跑几张图（线程级）。`calc.mpi × parallel_images` 别超核数 |
| `keep_intermediate_images` | `true` | `false` 时跑完删中间图目录（能量/受力仍在 `summary.json`/`band.traj`） |
| `idpp_steps` | 200 | （进阶，不在 `vcneb keys` 表里）IDPP 预优化每图 MDMin 步数上限 |
| `idpp_fmax` | 0.02 | （进阶）IDPP 预优化收敛标准（eV/Å） |

### 6.5 `output`

| 键 | 默认 | 说明 |
|----|------|------|
| `plot` | `true` | 结束时画 `barrier.png` |
| `save_trajectory` | `true` | 每步写 `band.traj`（`ase gui band.traj` 直接看） |

### 6.6 k 点怎么写（`calc.kpts`）

| 写法 | 含义 |
|------|------|
| `kpts: [4, 4, 4]` | 固定 Monkhorst-Pack 网格，原样用 |
| `kpts: {gamma: true}` | 只算 Gamma 点（分子/大真空层） |
| `kpts: {density: 0.03}` | 按倒空间间距 0.03 Å⁻¹ 反算：`n_i = int(\|b_i\|/density) + 1`，再对**两个端点晶胞**逐轴取最大，全程冻结 |

变胞时必须"一个网格到底"：每张图重新生成 k 点会让采样数在路径上跳整数格，
能量出现小台阶 → 假力、假能垒。`density` 模式还会检查路径上是否出现比两端
更"吃" k 点的晶胞，出现就 WARN 让你加密。

---

## 7. 三种 DFT 引擎怎么配

### 7.1 QE

```yaml
calc:
  engine: qe
  command: mpirun -np 8 pw.x
  mpi: 8
  pseudo_dir: /home/me/pseudos       # 放 *.upf
  env:
    OMP_NUM_THREADS: "1"
    LD_LIBRARY_PATH: /opt/qe/lib
  pseudopotentials:
    Cu: Cu.pbe-dn-rrkjus_psl.1.0.0.UPF
  kpts: [4, 4, 4]
  params:
    ecutwfc: 60
    ecutrho: 480
    occupations: smearing
    smearing: marzari-vanderbilt      # 变胞 NEB 推荐 mp/mv，不建议 gauss/fermi-dirac
    degauss: 0.02                     # Ry
    conv_thr: 1.0e-8
```

工具代写：`calculation=scf`、`tprnfor=.true.`、`tstress=.true.`、
`disk_io=minimal`、`verbosity=high`、`outdir`、`prefix`、`pseudo_dir`、`kpts`。

**坑**：QE 的 `params` 只认 ASE `ase/io/espresso_namelist/keys.py` 登记过的
关键词，**拼错会被静默丢弃**、不报错。改完一定看一眼
`images/02/espresso.pwi` 到底写了什么（`vcneb keys qe` 能列全）。

### 7.2 CP2K

```yaml
calc:
  engine: cp2k
  command: mpirun -np 8 cp2k.popt
  mpi: 8
  data_dir: /opt/cp2k/data
  kpts: [4, 4, 4]
  params:
    basis_set: DZVP-MOLOPT-SR-GTH
    potential: GTH-PBE-q11
    xc: PBE
    cutoff: 600              # Ry
    rel_cutoff: 60           # Ry
    eps_scf: 1.0e-6
    smearing: true
    smearing_method: marzari_vanderbilt
    smearing_sigma: 0.05     # eV
    added_mos: 12
    # ---- 自旋 / DFT+U / 布居（可选）----
    uks: true                # 磁性体系（含 +U）必须打开
    multiplicity: 3          # 2S+1，只有 uks 时生效
    hubbard_u: {Ni: 4.0}     # U_eff = U-J [eV]，按元素给
    hubbard_l: 2             # 2 = d，3 = f
    plus_u_method: mulliken  # mulliken（默认）/ lowdin / mulliken_charges
    mulliken: true           # 打印 Mulliken（每原子 α/β 布居、净电荷、自旋）
```

写手支持的全部键（`vcneb keys cp2k` 原样列出）：
`basis_set, potential, xc, cutoff, rel_cutoff, ngrids, eps_scf, max_scf,
eps_default, smearing, smearing_method, smearing_sigma, electronic_temperature,
window_size, added_mos, poisson_solver, vdw, extra, print_level, basis_set_file,
potential_file, uks, multiplicity, charge, relax_multiplicity, plus_u_method,
hubbard_u, hubbard_l, mulliken, mulliken_orbitals`。
**其它键一律被忽略**（`doctor` 会 WARN）。要写没被支持的东西，用
`params.extra` 原样插进 `&DFT`：

```yaml
  params:
    extra: |
      &SMEAR ON
        METHOD FERMI_DIRAC
        ELECTRONIC_TEMPERATURE [K] 300
      &END SMEAR
```

工具自己写、**不要**在 `params` 里重复的：`&GLOBAL`（`RUN_TYPE
ENERGY_FORCE`、`PRINT_LEVEL`）、`STRESS_TENSOR ANALYTICAL`、`&PRINT
&FORCES/&STRESS_TENSOR`、`BASIS_SET_FILE_NAME`/`POTENTIAL_FILE_NAME`、
`&MGRID`、`&QS`、`&KPOINTS`、`&SCF`（`SCF_GUESS`、`EPS_SCF`、`MAX_SCF`、
smearing、OT/DIAGONALIZATION、MIXING）、`&PRINT &MULLIKEN`、`&XC`、
`&SUBSYS`（`&CELL`/`&COORD`/`&KIND`，含 `&DFT_PLUS_U`）。

两个容易踩的点：

* 展宽分两类：`gaussian / methfessel_paxton / marzari_vanderbilt` 用
  `SIGMA [eV]`（写 `smearing_sigma`）；`fermi_dirac` 用
  `ELECTRONIC_TEMPERATURE [K]`（写 `electronic_temperature`）。写错 `doctor`
  会提示"该参数被忽略"。
* `multiplicity` 在开了 `smearing` 时只是初猜（最终自旋由费米能级决定）；
  不开 smearing 时 CP2K 把 α/β 电子数钉死成 `(N ± (M-1))/2`，等于强制磁态，
  要放开写 `relax_multiplicity: 0.00001`（需要 `added_mos > 0`）。

### 7.3 ABACUS

```yaml
calc:
  engine: abacus
  command: mpirun -np 8 abacus
  mpi: 8
  pseudo_dir: /opt/PP_ORB     # *.upf
  data_dir: /opt/PP_ORB       # *.orb
  pseudopotentials:
    Cu: Cu_ONCV_PBE-1.0.upf
  orbitals:
    Cu: Cu_gga_9au_100Ry_4s2p2d1f.orb
  kpts: {density: 0.03}
  params:
    basis_type: lcao
    ks_solver: genelpa
    ecutwfc: 100               # Ry
    dft_functional: pbe
    smearing_method: mp        # gauss / fd / mp / mv / fixed
    smearing_sigma: 0.01       # Ry
    scf_thr: 1.0e-7
    symmetry: 0
```

ABACUS **没有白名单**：`params` 里任何键都会原样写进 `INPUT`（代价是拼错也
照写、`doctor` 查不出来，只能自己看 `images/02/INPUT`）。工具代写的是
`pp, basis, kpts, basis_type, ks_solver, symmetry, scf_thr, scf_nmax,
cal_force, cal_stress, device, ntype`。

多值关键词要写成**一个空格分隔的字符串**（写成 YAML 列表会被 `str()` 成
`[4.0, 4.0]`，ABACUS 读不了）：

```yaml
    orbital_corr: "2 2 -1"
    hubbard_u: "4.0 4.0 0.0"
    nupdown: 2.0
```

找全字段：在自己 ABACUS 源码树里
`grep -oP '^### \K.*' docs/advanced/input_files/input-main.md`。

### 7.4 三个引擎对 `calc.params` 的差别（别想当然）

| | 未知键会怎样 | 怎么查 | 补救 |
|---|---|---|---|
| QE | **静默丢弃**（ASE 只认登记过的 namelist 键） | `vcneb keys qe`；看 `images/*/espresso.pwi` | 换等效写法，或 `params.input_data` 里按 namelist 分组写 |
| CP2K | **静默忽略**（内置写手只渲染上表 30 个键） | `vcneb keys cp2k`；`doctor` 会 WARN | `params.extra` 原样插 `&DFT` |
| ABACUS | **照写**（无白名单） | `vcneb keys abacus` 只列工具代写项 | 自己看 `images/*/INPUT` |

三种情形的共同核对方法：`<workdir>/config.resolved.yaml`（实际生效的配置）
+ `<workdir>/images/<idx>/`（实际生成的输入文件）。

---

## 8. 路径是怎么插出来的（插值四步）

`neb.interpolation` 只决定第 3 步走不走 IDPP，整条流程固定是：

```
align_cell → 线性插值（分数坐标 + 晶胞）→ IDPP 预优化 → 周期性展开 → 删点检查
```

### 8.1 第一步：`align_cell`（只换"写法"，不动结构）

同一个晶体可以写成不同的晶格矢量顺序/方向（`a1↔a2`、`a1→-a1`…）。插值在
**分数坐标**里做，两个物理上等价但索引不同的描述混用，路径会把原子互相塞进
对方身上。

代码（`band._align_final_cell`）枚举末态的 **6 种轴置换 × 8 种符号翻转 = 48
种**等价写法，按

```
score = ‖sym(F) − I‖ + RMS(原子位移) / 平均晶格长度,     F = c0⁻¹ · c1'
```

打分取最小；**只改末态的晶胞描述**（`states[-1].cell`），原子和晶胞本身不动。
用了非平凡变换时写日志（实测）：

```
  final cell re-indexed to line up with the initial one: (a,b,c) -> (b,c,a), signs +++
  final cell re-indexed to line up with the initial one: (a,b,c) -> (a,b,c), signs +--
```

`signs` 是三个轴的符号，`+++` 表示只换了顺序、`+--` 表示还翻了两个轴。

如果 48 种里最好的仍然不行（`strain > 15 %` 或 `disp > 2.5 Å`），说明两个端点
根本不是同一个超胞；只 WARN 不报错（可能真的就是要跑一个大变形路径），原文：

```
WARNING: the two endpoint cells do not line up (strain 40 %, atom RMS displacement 0.62 A)
and no re-indexing can fix that. The two structures are probably not the same supercell
(e.g. 2x1x1 vs 3x1x1) or the atom order does not match, so the interpolated images may
contain overlapping atoms; see neb.min_distance / neb.prune.
```

看到这句就该回去把两个端点做成公度超胞（`docs/faq.md` 第 2 节），而不是靠
删点硬撑。

### 8.2 第二步：线性插值

分数坐标做**最小像**差分（`dfrac -= round(dfrac)`），再按 `t = k/(n_free+1)`
线性混合；晶胞用同一个 `t` 线性插值（所以体积沿路径线性变化）：

```
cell_k = (1-t)·c0 + t·c1
frac_k = frac0 + t·dfrac
pos_k  = frac_k · cell_k
```

### 8.3 第三步：IDPP 预优化（`interpolation: idpp`，默认开）

每张可动图单独用 ASE 的 `IDPP` 计算器（目标距离矩阵 = 两端距离矩阵的线性
插值 `d1 + i·(d2−d1)/n`）跑 `MDMin`，`fmax = neb.idpp_fmax`（默认 0.02）、
最多 `neb.idpp_steps`（默认 200）步；**晶胞固定**，只动原子。

为什么不直接用 ASE 的 `idpp_interpolate`：它拒绝"晶胞沿路径变化"的带，而
变胞 NEB 恰恰就是这种。IDPP 的目标函数只用原子间距离，所以逐图单独优化
完全等价，又允许晶胞不同。失败时只 WARN 并保留线性路径：

```
WARNING: IDPP interpolation failed (...); keeping the linear path
```

`interpolation: linear` 就是跳过这一步；单原子体系会提示
`only one atom: IDPP has nothing to do, keeping the linear path`。

### 8.4 第四步：周期性展开

IDPP 可能把原子推过晶胞边界。路径的切向/弹簧/间距都用**直角坐标差**，一个被
wrap 的原子看起来就像跳了一整个晶格长度（~5 Å），整条带会炸。所以沿带做
一次最小像展开 `x_i = x_{i−1} + MIC(x_i − x_{i−1})`——描述同一组结构，但坐标
连续了（真正算能量时 `atoms.wrap()`，物理量不变）。有动作时日志：

```
  unwrapped 3 image(s) so the path is continuous across the periodic boundary
```

### 8.5 建完路径后打印什么

```
initial path built: 9 images (7 movable), interpolation=idpp
smallest interatomic distance on the path: 2.605 A
path spacing (RMS atom displacement between neighbouring images): 0.143, 0.133, 0.140, 0.158, 0.158, 0.123, 0.140, 0.139 A
```

（这六行是 `docs/example-bi.md` §4 里那个真跑出来的 dry-run 输出，可复现。）

* `smallest interatomic distance`：全路径最短原子间距；`< 1.2 Å` 会额外 WARN
  （原子顺序或超胞不对）；
* `path spacing`：相邻图的 RMS 原子位移，正常应大致均匀。最大/最小超过 3 倍
  会 WARN（`the path spacing is very uneven…`，通常就是端点原子顺序不对）；
  如果是因为删点才不均匀，会自动改口：
  `(uneven on purpose: images 2, 3, 4 were dropped, see neb.min_distance)`。

---

## 9. 删点：`min_distance` / `prune`

IDPP 修不了"晶胞本身被插值成错误形状"的图（比如相变中间某张图晶胞被压到
接近零体积），这种图里两个原子会重合到 0.1 Å，任何 DFT 都算不出来。做法和
USPEX 丢掉坏候选一样：直接丢掉这张图，剩下的继续跑 NEB。

规则（`band._prune_close_images`）：

* 只检查**可动图**，端点永远不动；
* 判据：该图最小像最近原子间距 `d_min < neb.min_distance`（默认 2.0 Å）；
* `min_distance: 0` 或 `prune: false` → 完全跳过（只保留 §8.5 的 WARN）；
* `prune: error` → 报错停下，列出图号和间距；
* `prune: auto` → 删点继续；`pruned` 记原图号，`n_free` 变成剩余可动图数；
  全删光 → 报错；剩不到 3 张再 WARN。

下面三段日志都能用这个命令复现（不用任何 DFT 程序，全程 `--dry-run`）：

```bash
# 一条命令跑完下面四段（自动生成/删除临时 yaml，几秒钟）
bash examples/prune_demo.sh

# 或者手动来一遍。注意复制到 examples/ 里，别复制到 /tmp：
# 配置里的 ./structures/... 是相对 YAML 文件所在目录解析的
cp examples/lj_ar.yaml examples/pr.yaml
#   把 neb.images 改成 8，再加一行 neb.min_distance: 3.1（或下面的 3.3 / 3.5），
#   prune 依次取 auto / auto / error / auto
vcneb run examples/pr.yaml --dry-run
```

**auto，删掉一部分**（`min_distance: 3.1`）：Ar 的两个端点最近邻是 3.72 Å，
但线性插出来的中间图会被压到 3.08 Å（这就是"插值本身不物理"的典型），
于是正好中间两张被丢掉：

```
  unwrapped 1 image(s) so the path is continuous across the periodic boundary
WARNING: dropped images 4, 5 (closest contact 3.08, 3.08 A < neb.min_distance=3.1 A); band now has 6 movable images
initial path built: 8 images (6 movable), interpolation=idpp
smallest interatomic distance on the path: 3.120 A
path spacing (RMS atom displacement between neighbouring images): 0.112, 0.112, 0.112, 0.337, 0.112, 0.112, 0.112 A
  (uneven on purpose: images 4, 5 were dropped, see neb.min_distance)
```

注意 `band now has 6 movable images`：**删点是在插值之后、开跑之前做的**，
所以图总数从 10（8 可动 + 2 端点）变成 8，相邻图距离从 0.112 变成
`0.112, ..., 0.337, ...`（被删处变宽），日志最后一行专门说明这是故意的。

**auto，删太多（会额外警告）**（`min_distance: 3.3`）：

```
WARNING: dropped images 2, 3, 4, 5, 6, 7 (closest contact 3.26, 3.12, 3.08, 3.08, 3.12, 3.26 A < neb.min_distance=3.3 A); band now has 2 movable images
WARNING: only 2 movable image(s) left: a band that coarse usually overestimates the barrier badly.  The real fix is to make the endpoints commensurate (same supercell, same atom order), not to prune images.
initial path built: 4 images (2 movable), interpolation=idpp
smallest interatomic distance on the path: 3.466 A
path spacing (RMS atom displacement between neighbouring images): 0.112, 0.787, 0.112 A
  (uneven on purpose: images 2, 3, 4, 5, 6, 7 were dropped, see neb.min_distance)
```

只剩 2 张可动图还能跑，但这条带粗到几乎没有意义了 —— 该回去改超胞/原子顺序。

**prune: error，直接停**（`min_distance: 3.5`）：

```
ERROR: 8 interpolated image(s) put atoms closer than neb.min_distance=3.5 A (images 1, 2, 3, 4, 5, 6, 7, 8: 3.47, 3.26, 3.12, 3.08, 3.08, 3.12, 3.26, 3.47 A).
  Fix the endpoints (atom order / supercell must match) or set neb.prune: auto to drop those images and continue.
```

**auto，全删光**（`min_distance: 3.5` + `prune: auto`）：

```
ERROR: every interpolated image is unphysical (closest contact 3.47, 3.26, 3.12, 3.08, 3.08, 3.12, 3.26, 3.47 A < neb.min_distance=3.5 A), so there is no band left.
  The two endpoint structures almost certainly do not correspond to each other (different supercell, atom order or species order).  Increase neb.images only helps if the endpoints are right.
```

两种 `min_distance` 写法都合法（`3.3` 和 `3.30` 等价），日志里按原样打印
`neb.min_distance=3.3 A`。

删点**不是**修复手段，只是"让这条带还能跑完"的兜底：被删的图不会回来，带变
粗（相邻图距离变大），能垒精度下降。删点结果记在 `state.json` 的
`pruned_images` 里，`--restart` 原样沿用（不重新插值、不重新删点）。

---

## 10. 加压力：`neb.pressure`

相变能垒一般要在**给定压力**下算。写 `neb.pressure: 2.55`（GPa），程序把整条
带的能量换成焓、晶胞力换成焓的导数：

```
H = E + P·V
dH/dh = −V·(σ + P·I)·h⁻ᵀ
```

实现（`band.EnthalpyCalculator`）：包住引擎的 calculator，能量加 `P·V`、应力
的三个正应力分量加 `P`。于是

* 端点弛豫会弛豫到"焓的极小点"（该压力下的平衡体积），不是能量的极小点；
* 每张图的晶胞力、日志能量、`summary.json`、`barrier.png` 全是**焓**；
* 势垒 = 焓差，天然包含 `P·ΔV`。

两条必须记住的：

1. **DFT 输入里的压力要留 0**（QE `press`、CP2K `PRESSURE`、VASP `PSTRESS`、
   ABACUS 的 `press*`）。工具自己只写 0 压的输入，你别在 `params` 里再加。
2. `pressure` 只对 `variable_cell: true` 有意义。固定晶胞时 `P·V` 是常数，
   能垒不变（能量零点平移而已）。

实测（内置 LJ Ar 的 fcc→hcp，4 原子；**只用来验证"压力确实进了焓"，不要当
物理结果**——LJ 氩体模量很小，5 GPa 会压掉约 30 % 体积）：

| `neb.pressure` | 端点体积 | 势垒 |
|----------------|----------|------|
| 0 GPa | 145.99 Å³ | 0.0142 eV |
| 5 GPa | 100.91 Å³ | 0.0700 eV |

日志会多一行：

```
  pressure     : 5.0 GPa (every energy below is the enthalpy H = E + P*V)
```

---

## 11. 续算与重来

DFT 的 NEB 一定会被打断（排队时间到、机器重启、参数写错）。`vcneb run` 的
五个开关里有三个是管这件事的：

| 开关 | 做的事 | 什么时候用 |
|------|--------|------------|
| `--restart` | 从 `state.json` 读回上次停下的那一步（每张图的晶胞、坐标、能量、优化器步长），从下一步接着跑 | 中断后续算、步数不够加步数 |
| `--force` | 删掉上一次的产物，从第一步重新开始 | 换了结构/引擎/展宽，或想把路径重插 |
| `--max-steps N` | 只覆盖这一次的 `neb.max_steps`，不写回 YAML | 续算时加步数 |
| `--dry-run` | 只弛豫端点 + 插值 + 打印，不算能量 | 第一次配 YAML，先看路径对不对 |
| `--skip-engine-check` | 跳过"可执行文件/赝势在不在"的检查 | 登录节点上没有 `pw.x`、要靠提交脚本跑时 |

### 11.1 `state.json` 里到底存了什么

```json
{
  "step": 137,                          // 已完成的 NEB 步数
  "variable_cell": true,
  "n_atoms": 12,
  "n_movable": 7,                       // 中间可动图数（删点后 < images-2）
  "pruned_images": [2, 3],              // 被删掉的图
  "optimizer_dt": {"1": 0.14, "2": 0.1},// 每张图各自的 FIRE 时间步
  "images": [                           // 所有图，含两个端点
    {"cell": [[...]], "pos": [[...]], "energy": -1234.5678}
  ]
}
```

存的是**几何 + 能量 + 优化器状态**。不存的是 SCF 波函数、`images/NN/` 里的
中间文件、任何 Hamiltonian。所以：

* `--restart` 不是"接着上一步的 SCF 收敛"，而是"接着第 137 步的构型继续降力"，
  每一步的每张图仍然从头做一次 SCF；
* 端点弛豫**不会重做**——端点也在 `images` 里，已经带着弛豫后的晶胞和坐标；
* 续算的开销 = 剩下的步数 × 每步的图数，不会因为"续"而变便宜。

### 11.2 会拦住你的几条报错（原文）

**① 目录里已经有东西，但没带开关**

```
ERROR: ./run already contains a NEB run.
  use `--restart` to continue it, or `--force` to start over
```

想接着跑就加 `--restart`；想从头来就加 `--force`（或自己删掉 `run/state.json`）。

**② 步数用完了**

```
ERROR: this band is already at step 400 while max_steps=400: there is nothing
left to do.
  Raise neb.max_steps in the config (or pass e.g. `--max-steps 900`) and run
  with --restart again.
```

`vcneb run neb.yaml --restart --max-steps 900`。注意 `--max-steps` 是**总步数
上限**，不是"再加 900 步"。

**③ 结构变了**

```
ERROR: restart state in ./run/state.json has 12 atoms, the current setup has 24
ERROR: restart state ./run/state.json is inconsistent: 9 images but
n_movable=7.  Delete it and rebuild the band.
```

端点原子数/图数变了就没法续，只能 `--force`。

**④ 只改了 `neb.images`**

```
WARNING: the saved band has 5 movable images, the config asks for 7; keeping
the saved band.  The image count cannot change in the middle of a run -- if you
really want a different path, delete ./run/state.json and start again.
```

续算时图数按 `state.json` 里那条带走，YAML 的 `images` 被忽略。真要改就 `--force`。

### 11.3 一次典型的中断 → 续算

```bash
vcneb doctor neb.yaml                # 先查配置（不花机时）
vcneb run neb.yaml --dry-run         # 再查路径（不花机时）

vcneb run neb.yaml --max-steps 60    # 正式跑；或跑到一半 Ctrl-C
# Ctrl-C 时日志最后一行：
#   interrupted; progress is in state.json, continue with --restart
# 退出码 130

vcneb status ./run                   # 第几步、每张图的能量/力/体积
# 没收敛时它直接把下一条命令打给你：
#   not converged yet; continue with
#     vcneb run <config> --restart --max-steps {steps+500}

vcneb run neb.yaml --restart --max-steps 560
```

续算这次启动时**不会重做端点弛豫**，日志里直接是：

```
k-points: Gamma   (same grid for all images)
endpoints: 6 atoms, V0=218.2974 -> V1=218.2974 A^3
restarting from step 20 (7 images)
------------------------------------------------------------------------------
starting vc-NEB: 5 movable images, f_th=0.01 eV/A, max_steps=500
------------------------------------------------------------------------------
```

（`restarting from step N (M images)` 里的 M 是含端点的总图数，N 是
`state.json` 里记的已完成步数，所以下一次迭代的步号是 N+1。）

`vcneb status` 对**中断后还没写 `summary.json`** 的运行也能用：

```
run     : ./run
progress: step 137 of an unfinished run
resume  : vcneb run <config> --restart
```

### 11.4 续算时能改什么、不能改什么

| 改这个 | 结果 |
|--------|------|
| `neb.max_steps` / `--max-steps` | 随便改，这就是续算的主要用途 |
| `neb.f_th` | 可以，判据立刻按新值执行 |
| `neb.dt` / `max_step` / `cell_max_step` / `spring_k` | 可以，下一步生效 |
| `calc.params`（截断、展宽、k 点、DFT+U …） | **不要**。前后两段不是同一套参数，势垒是两段数拼起来的 |
| `neb.pressure` / `neb.variable_cell` / `calc.engine` | **不要**。能量零点（焓/能量）和自由度全变了 |
| `structure.initial` / `final`（换超胞） | 不行，会报 atoms 不一致 |

要改"不要"那几行，就 `--force` 从头跑——那不是续算，是另一次计算。

---

## 12. 输出文件

跑完之后 `project.workdir`（例子里是 `./run`）长这样：

```
run/
├── summary.json            # 全部结果（收敛后写一次）
├── barrier.csv             # 每张图的 E / dE / max|F| / V（+ 自旋、电荷）
├── barrier.png             # 势垒图（output.plot: true 且装了 matplotlib）
├── band.traj               # 最近一步所有图的 ASE 轨迹（save_trajectory）
├── state.json              # 续算状态，每一步覆盖写
├── history.jsonl           # 每一步一行
├── config.resolved.yaml    # 实际生效的配置（含默认值），复现用
├── endpoint_initial.vasp   # 只有 --dry-run 才写：弛豫+重排后的两个端点
├── endpoint_final.vasp
├── images/
│   ├── end_initial/        # 端点 0 的引擎输入/输出
│   ├── 01/ 02/ ... 07/     # 中间图；单图失败重试时多出 01_retry1/ 01_retry2/
│   └── end_final/
└── logs/
    ├── vcneb.log           # 主日志（终端上的一切 + 时间戳 + 累计耗时）
    ├── relax_initial.log   # 端点弛豫日志
    └── relax_final.log
```

每个 `images/NN/` 里是该引擎自己写的文件：

| 引擎 | 目录里的文件 |
|------|--------------|
| `qe` | `espresso.pwi`、`espresso.pwo`、`vcneb.save/`（`outdir`） |
| `cp2k` | `cp2k.inp`、`cp2k.out`、`stdout.log` |
| `abacus` | `INPUT`、`STRU`、`KPT`、拷进来的赝势/轨道文件 |
| `emt` / `lj` | 空目录（没有文件级输入） |

`neb.keep_intermediate_images: false` 时，跑完会删掉中间图的目录（两个端点的
留着），能量/力仍在 `summary.json` 和 `band.traj` 里。DFT 跑十几张图时这些
目录能占好几个 GB。

### 12.1 `summary.json`

```json
{
  "engine": "abacus",
  "tangent": "improved",
  "converged": true,
  "steps": 214,
  "kpts": [8, 5, 3],
  "barrier_eV": 0.0431,
  "reverse_barrier_eV": 0.0512,
  "energy_initial_eV": -12345.6789,
  "energy_final_eV": -12345.6712,
  "energy_max_eV": -12345.6355,
  "image_max": 4,
  "f_th_eV_per_A": 0.05,
  "variable_cell": true,
  "n_atoms": 12,
  "n_movable_images": 7,
  "seconds": 18342.5,
  "final_max_force_atom_eV_per_A": 0.0418,
  "final_max_force_cell_per_atom_eV_per_A": 0.0061,
  "images": [
    {"image": 0, "energy_eV": -12345.6789, "delta_eV": 0.0,
     "max_force_eV_per_A": 0.0, "volume_A3": 439.032,
     "total_spin": null, "net_charge": null, "mulliken": null}
  ]
}
```

（上面的**字段名和嵌套结构**与真实文件一致；里面的能量/步数/耗时是占位值，
别拿去对账。`kpts`、`converged`、`image_max` 这些才是你要读的。）

* `barrier_eV = max(E) − E[0]`，`reverse_barrier_eV = max(E) − E[-1]`；
  设了 `neb.pressure` 时这里的 `E` 全是焓 `H = E + P·V`；
* `kpts` 是那条"全程冻结"的网格（`calc.gamma: true` 时是字符串 `"gamma"`）；
* `tangent` 是这次实际用的路径切向（`improved` / `central`），`state.json` 里也存了一份，
  `vcneb status` 会打印出来 —— 续算换了算法时能一眼看出跑的是哪一个；
* 两个 `final_max_force_*` 只有收敛时才有；
* `image_max` 是能量最高的图号。它等于 0 或最后一张，说明**这条路径上没有
  鞍点**（端点自己最高），此时 `barrier_eV` 约等于 0。

### 12.2 `barrier.csv`

```
image,energy_eV,delta_eV,max_force_eV_per_A,volume_A3[,total_spin,net_charge]
0,-12345.67890000,0.00000000,0.000000,439.032000
1,-12345.66120000,0.01770000,0.213456,436.477000
```

（列名和行的顺序是真的；数值是占位，`volume_A3` 那一列取自真跑过的
`docs/example-bi.md` 算例。）

后面两列只有在**至少一张图给出了 Mulliken 布居**时才追加：目前只有 `cp2k`
引擎会自动在输入里写 `&MULLIKEN ON` 并解析它；`qe`/`abacus` 引擎不解析布居，
所以那两列不会出现。有布居时日志会多一行逐图自旋：

```
Mulliken total spin moment along the path (images 0..8): 0.00  0.12  0.45 ...
```

整条路径的自旋总量变化超过 0.05 μB 时还会给你一条 WARNING（见 §14）。

### 12.3 `history.jsonl`

每一步一行，用来画收敛曲线：

```json
{"step": 1, "energies": [-0.4953, -0.4785, ...], "max_force_atom": 0.3576,
 "max_force_cell_per_atom": 0.0012, "climbing_image": 3}
```

`energies` 是**全部图**（含端点）的能量数组，长度 = 图数；`climbing_image`
是这一步攀爬像的编号，`null` 表示这一步没开 CI。

### 12.4 `config.resolved.yaml`

`vcneb` 补全默认值之后的配置原样写出来。**复现请存这个文件**，而不是你手写
的那份——你没写的键（`spring_k`、`cell_max_step`、`f_th` …）都在这份里定死了。

---

## 13. 常见错误

所有消息都在 `run/logs/vcneb.log` 里，`WARNING` 行不影响继续跑，`ERROR` 行
会终止（退出码 1）。

### 13.1 配置/环境（`vcneb doctor` 就能查出来）

| 日志里看到 | 意思 | 怎么办 |
|------------|------|--------|
| `engine is not ready:` + `executable 'pw.x' not found in PATH` | `calc.command` 的第一个词找不到 | 写全路径，或换 `mpirun -np 8 /opt/qe/bin/pw.x` |
| `pseudopotential for Bi not found: ./Bi.upf` | `calc.pseudo_dir` + `calc.pseudopotentials` 对不上 | 改文件名（区分大小写），或改 `pseudo_dir` |
| `calc.pseudo_dir does not exist` | 路径相对**配置文件所在目录**（不是你的当前目录） | 用绝对路径最省事 |
| `calc.params.ecutwfc is not set` | QE 必须给截断 | 加 `ecutwfc: 60` |
| `no k-points given: add calc.kpts: [4, 4, 4]` | 三种 k 点写法一个都没给 | 见 §6.6 |
| `set either calc.kpts or kspacing, not both` | 两种写法同时给了 | 删掉一个；变胞路径建议只用 `kpts` 整数网格或 `{density: X}` |
| `calc.params.multiplicity is ignored without uks: true` | CP2K 没开自旋极化却给了多重度 | 加 `uks: true`（磁性体系必须） |
| `ERROR: ./run already contains a NEB run.` | 目录里有旧结果 | `--restart` 或 `--force`（§11） |
| `ERROR: this band is already at step 400 while max_steps=400` | 步数到顶 | `--restart --max-steps 900` |
| `ERROR: ./run does not look like a vcneb run directory` | `status` 指错了目录 | 指向含 `state.json`/`summary.json` 的那层 |
| `ERROR: no summary.json in ./run; run the NEB first` | 还没跑完就想 `plot` | 先 `run`；只想看图可以先 `--dry-run` 看结构 |

### 13.2 结构/插值

| 日志里看到 | 意思 | 怎么办 |
|------------|------|--------|
| `the two endpoint cells do not line up (strain 95 %, atom RMS displacement 4.27 A) and no re-indexing can fix that` | 两个端点根本不是同一个超胞（典型：六方 2×1×1 vs 单斜 3×1×1） | 用 `tools/build_common_cell.py` 做成公度超胞，见 `docs/example-bi.md` §1（这条 95 % 是把 `Bi.vasp` 和 `Bi2.vasp` 直接喂给 vcneb 实测出来的） |
| `the path spacing is very uneven (see the list above)` | 原子顺序没对齐（除了删点造成的，删点会另说 "uneven on purpose"） | 看 `run/endpoint_initial.vasp` 和 `endpoint_final.vasp`，重排末态原子顺序 |
| `atoms are very close on the initial path (< 1.2 A)` | 同上，或者端点本身有问题 | 同上；顺便查最近邻距离 |
| `dropped images 4, 5 (closest contact 3.08, 3.08 A < neb.min_distance=3.1 A); band now has 6 movable images` | `min_distance` 太大，把合法图也删了 | 把 `min_distance` 降到"真实最近邻 − 0.6 Å"左右（见 §9，`min_distance: 3.1` 那个例子就是这条） |
| `every interpolated image is unphysical (closest contact 0.42 A < neb.min_distance=2.0 A), so there is no band left` | 整条路径塌了 | 端点不公度或原子顺序错，先修这个 |
| `only 1 movable image(s) left: a band that coarse usually overestimates the barrier badly` | 删点删太狠 | 提高 `neb.images`，并降低 `min_distance` |
| `IDPP interpolation failed (...); keeping the linear path` | IDPP 不收敛（极少见） | 可忽略；或在 YAML 里写 `interpolation: linear` |
| `unwrapped 1 image(s) so the path is continuous across the periodic boundary` | 正常行为（不是错误）：某张图的某个原子被折叠到另一侧 | 不用管 |

### 13.3 跑的过程中

| 日志里看到 | 意思 | 怎么办 |
|------------|------|--------|
| `image 3: attempt 1 failed (CalledProcessError: ...), retrying in .../03_retry1` | 单张图第一次失败，自动重试（最多 3 次） | 一般不用管；连续失败会变成下面那条 |
| `ERROR: image 3: SCF/calculator failed after 3 attempts.` | 三次都失败 | 去 `run/images/03_retry2/` 看 `cp2k.out` / `espresso.pwo` / ABACUS 输出 |
| `ERROR: two atoms are only 0.55 A apart in the band` | 优化步太大，原子撞在一起 | 降 `neb.max_step`（原子）/ `neb.cell_max_step`（晶胞） |
| `ERROR: image 4 cell collapsed or became left-handed (volume -0.42 A^3)` | 晶胞步太大，把晶胞推翻了 | 降 `neb.cell_max_step`；变胞问题严重时也降 `neb.dt` |
| `WARNING every image returned the same energy (spread < 1e-9 eV)` | 能量完全不动：要么势能面真的平，要么引擎根本没算（原子被固定、calculator 没生效） | 看 `images/NN/` 里的输入文件 |
| `WARNING frozen k-mesh [8, 5, 3] is too coarse for this cell ...` | 路径上晶胞跑出了两端点覆盖的范围，冻结的网格偏粗 | 换成 `calc.kpts: {density: 0.03}` 让它按整条路径定网格 |
| `WARNING endpoint 'initial' did not reach 0.01 eV/A in 200 steps` | 端点没弛豫干净 | 加 `neb.endpoint_steps`；或该压力/该精度下端点不是极小点 |
| `WARNING reached max_steps=400 without converging (last max|F_perp|=0.0831 eV/A, target 0.05)` | 步数不够 | `--restart --max-steps 900` |
| `WARNING the total spin moment changes along the path` + Mulliken 一行 | 磁性在路径上变了（见 §14） | 看 `barrier.csv` 的 `total_spin`；必要时改用展宽让自旋浮动 |
| `interrupted; progress is in state.json, continue with --restart`（退出码 130） | 你按了 Ctrl-C | `--restart` 继续（§11） |

---

## 14. 限制与注意

* **切向**默认用 Henkelman–Jónsson 能量加权差分（`neb.tangent: improved`）：
  上坡图用前向差、下坡图用后向差，平台区按能量差加权，避免带在弯折处"切角"。
  旧的中心差分仍可选（`neb.tangent: central`），在对称/近似直线的带上两者
  给出**完全相同**的结果（实测：LJ 34 步 / 0.0383 eV、删点后的不均匀 LJ 28 步 /
  0.0343 eV 都一模一样；EMT Cu(111) 41 步 → 39 步，能垒 0.4435 eV 不变），
  但在弯折带上中心差分的切向会偏，垂直力方向偏差可达 20° 以上。
这些不是 bug，是方法的边界，写文章时要用得上。

**1. 共格（coherent）能垒 ≠ 成核能垒。**
两个端点必须能塞进同一个晶胞、原子一一对应，这条路径描述的是"整块晶体同时
变形"的能垒。真实相变通常从界面/缺陷处开始，成核能垒一般比这个小得多。
报告里要写清是 coherent 的。

**2. 整条路径只有一个 k 网格。**
这是有意的：变胞时若每张图按自己的晶胞重算 k 点，网格切换会给能量带来
台阶，混进能垒里。代价是体积变化大时这个网格对某几张图偏粗/偏密，
所以 `calc.kpts` 要按"两端点里最大的那个"来选（`{density: X}` 就是干这个的）。

**3. 晶胞力。**
内部用的是 `F_cell = −V·σ·cell⁻ᵀ`（`σ` 已含外压 `P·I`），直接作用在晶胞
矩阵的 9 个分量上。**不要**再额外减 `−fracᵀ·F_atom` ——原子的贡献由原子
自由度自己走。手写 ASE filter 时这一项最容易重复计算。

**4. `f_th` 的下限由 DFT 噪声定。**
SCF 只收到 `1e-6`、k 点又稀的时候，力的噪声可能有 0.03 eV/Å，`f_th: 0.05`
就只是在拟合噪声。经验做法：先跑一张图，连续两次 SCF 的力差 < `f_th`/3，
再决定 `f_th` 能压到多少（`tests/test_cellforce.py` 里有类似的判据）。

**5. 攀爬像（CI）给的是近似鞍点。**
CI 把能量最高那张图推到一阶鞍点上，但有限图数下它只是近似；要频率/要
确认鞍点，拿 `band.traj` 里那张图去做 Hessian（`ase.vibrations`）。

**6. 展宽会改变能垒。**
金属体系里 `smearing_sigma` 从 0.01 Ry 加到 0.05 Ry，势垒动 10 meV 是正常的。
报告里给出你用的展宽，并至少检查一次收敛。

**7. 压力只支持静水压（标量）。**
`neb.pressure` 把三个正应力分量同时加 `P`。要各向异性应力（例如外延应变
+ 单轴压），得自己写一个 `EnthalpyCalculator`（`vcneb/band.py` 里 30 行）。

**8. 端点必须同赝势、同截断、同展宽、同 k 网格。**
这不是"建议"，是这套方法的前提：能垒是两个态的差，两个态用两套参数就没有
可比性。工具强制了同一个 k 网格，其余要你自己保证。

**9. `spring_k` 不影响收敛后的势垒。**
它只影响收敛速度和路径的形状（NEB 的弹簧力在收敛时不做功）。调它是因为
"图挤在一起/散开"，不是因为它改变物理。

**10. 一个机制 ≠ 所有机制。**
你给两个端点，它就给你这一条路的能垒。换原子顺序/换超胞/换界面取向就是
另一条路径，能垒会不一样。挑哪条是物理判断，工具不替你做。

---

## 15. 例子与测试

### 15.1 `examples/`

| 文件 | 引擎 | 算什么 | 花多久 |
|------|------|--------|--------|
| `lj_ar.yaml` | `lj` | Ar 的 fcc(111) ABC→ABA 堆垛滑移，变胞 | 2 秒 |
| `emt_cu111.yaml` | `emt` | Cu(111) 表面层滑移，固定晶胞 | 几十秒 |
| `qe_li.yaml` | `qe` | Li 的 bcc→fcc，变胞 | 看机器 |
| `cp2k_li.yaml` | `cp2k` | 同上（GPW） | 看机器 |
| `abacus_cu.yaml` | `abacus` | Cu 的层错路径（LCAO） | 看机器 |
| `prune_demo.sh` | `lj` | 一次跑出 §9 里四段删点日志（临时 yaml 自动删） | 10 秒 |

没有 DFT 也能先把流程走通：

```bash
vcneb run examples/lj_ar.yaml
vcneb status examples/runs/lj_ar
vcneb plot examples/runs/lj_ar -o /tmp/ar.png
```

三个 DFT 例子在跑之前一定先 `vcneb doctor`，它会告诉你缺什么：

```bash
vcneb doctor examples/qe_li.yaml
vcneb doctor examples/cp2k_li.yaml
vcneb doctor examples/abacus_cu.yaml
```

`examples/structures/` 是这些例子的端点结构；`examples/make_111_stack.py`
是生成 fcc/hcp 堆垛结构的脚本（`python3 examples/make_111_stack.py`）。

### 15.2 `tools/`

| 脚本 | 干什么 | 额外依赖 |
|------|--------|----------|
| `tools/build_common_cell.py` | 把两个相放进同一个公度超胞：常规胞（spglib）→ 整数超胞矩阵 → Niggli 约化 → 刚体对齐（SVD）→ 原子匹配（`linear_sum_assignment`），最后写两个原子数相同、顺序对应的 POSCAR | `spglib`、`scipy` |
| `tools/setup_bi_case.sh` | 用上面两个 POSCAR 生成 Bi 的算例目录和 YAML | 无 |

用法（这就是 `docs/example-bi.md` 里的 Bi 例子）：

```bash
python3 tools/build_common_cell.py \
  --a Bi.vasp --b Bi2.vasp --out ./build \
  --ma "-2 0 1, 0 2 0, -1 2 0" \
  --mb "-2 -2 1, -1 -1 -1, -1 -2 -1"
```

输出里要盯三行：两个超胞的原子数必须相等、`residual strain B->A`（越小越好，
Bi 这个例子里是 4.76 %）、`matching: ... (mean X A)` 和 `worst pair`（原子
位移的平均值/最大值，Bi 里是 1.385 / 2.40 Å）。strain 大过 15 % 脚本会警告。

### 15.3 测试

```bash
python -m pytest tests -q                        # 9 passed（约 9 秒）

# 另外两个是独立脚本（有 main()，pytest 不收集，用 python 直接跑）：
PYTHONPATH=. python3 tests/test_cellforce.py         # 晶胞力 vs 有限差分（fractional 固定）
PYTHONPATH=. python3 tests/test_cellforce_basis.py   # 同上，Cartesian 基组那一路
PYTHONPATH=. python3 tests/t_plumbing_pop.py         # 布居 → barrier.csv → summary.json → 图
```

`PYTHONPATH=.` 只在没 `pip install -e .` 时需要（脚本从仓库根目录 import `vcneb`）。
`pytest tests -q` 的 9 个用例 = `test_interpolation_guards.py` 7 个 +
`test_pressure.py` 2 个。

* `test_cellforce.py`：把 `F_cell = −V·σ·cell⁻ᵀ` 与有限差分（整体缩放晶胞）
  对拍，确认没有重复计算 `−fracᵀ·F_atom`；
* `test_pressure.py`：`H = E + P·V` 的端点和力与 ASE 的 `ExpCellFilter`
  对拍；
* `test_interpolation_guards.py`：7 个用例，删点/间距/几何守卫
  （`min_distance` 太大、只剩一张图、晶胞塌缩），每个用例真的建一遍带、
  然后检查日志/异常文本，也就是 §9 那几段
* `tests/cp2k_u.yaml`、`tests/cp2k_u_mini.yaml`：不是测试，是 `vcneb doctor`
  的夹具（UKS + DFT+U + MULTIPLICITY + Mulliken），用来看 CP2K 关键字怎么
  被解析/校验，不跑 CP2K：

```bash
vcneb doctor tests/cp2k_u.yaml
#   CP2K  [...]  cutoff=600 Ry  basis={'Ni': 'DZVP-MOLOPT-SR-GTH-q18'}  kpts=[6, 6, 6]
#               UKS mult=3  DFT+U(Ni U-J=4eV L=2)  mulliken=on
```

### 15.4 看完了下一步

* 物理问题（超胞、k 点、压力、DFT+U、磁性、smearing、续算策略）：
  [`docs/faq.md`](docs/faq.md)
* 一个完整的真实算例（Bi A7 → Bi-II，2.55 GPa，12 原子，ABACUS）：
  [`docs/example-bi.md`](docs/example-bi.md)
* 代码地图（命令 → 函数 → 行号，插值五步、删点分支表、`state.json` 字段）：
  [`docs/internals.md`](docs/internals.md)
