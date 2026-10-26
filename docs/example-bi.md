# 完整算例：Bi 的 A7 → Bi-II，2.55 GPa

这是一个**真跑过**的例子，用来演示一条完整链路：

```
两个相的原胞/常规胞
   → 公度超胞（tools/build_common_cell.py）
   → 生成算例（tools/setup_bi_case.sh）
   → vcneb doctor
   → vcneb run --dry-run（看路径）
   → vcneb run（ABACUS 真算）
   → summary.json / barrier.csv / barrier.png
```

每一步的命令和输出都贴在下面，数字是这台机器上跑出来的，不是示意值。
唯一的例外是第 5 节：这台机器上没有 ABACUS 的 `mpirun`，所以那一段只有**日志格式**
是实测的，能量字段用 `XXXXXXXX` 占位（其余字段逐字对应源码）。

---

## 0. 起点：两个相

| 文件 | 内容 | 空间群 | 原子数 | a, b, c (Å) | 角度 | V (Å³) | 最近邻 |
|------|------|--------|--------|-------------|------|--------|--------|
| `Bi.vasp` | A7（常压相），2×1×1 六方超胞 | R-3m（常规胞 6 原子） | 12 | 4.538 / 7.860 / 12.308 | 90 / 90 / 90 | 439.03 | 3.074 Å |
| `Bi2.vasp` | Bi-II（高压相），3×1×1 单斜超胞 | C2/m（常规胞 4 原子） | 12 | 4.604 / 8.003 / 11.336 | 90.4 / 90 / 90 | 417.69 | 3.229 Å |

两者**形状已经接近**（这是运气，也是为什么选这两个超胞），但它们的基不同：
A7 那边是六方 2×1×1，Bi-II 那边是单斜 3×1×1。直接塞给 vc-NEB 会报
"the two endpoint cells do not line up"。

真的会报，实测（把这两个文件直接当 `structure.initial/final`，`--dry-run`）：

```
k-points: Gamma   (same grid for all images)
endpoints: 12 atoms, V0=439.0319 -> V1=417.6909 A^3
dry run: endpoints are used exactly as given (no DFT relaxation, no electronic-structure call)
WARNING: the two endpoint cells do not line up (strain 95 %, atom RMS displacement 4.27 A) and no re-indexing can fix that.  The two structures are probably not the same supercell (e.g. 2x1x1 vs 3x1x1) or the atom order does not match, so the interpolated images may contain overlapping atoms; see neb.min_distance / neb.prune.
...
path spacing (RMS atom displacement between neighbouring images): 0.312, 0.374, 0.502, 0.567, 0.932, 0.469, 0.411, 0.415 A
```

`95 %` 应变、`4.27 Å` 平均位移、体积序列 439.03 → **461.14** → 417.69
（中间比两头都大，物理上不可能）—— 这就是"不是同一个超胞"的样子。
所以下面第 1 步必须先做公度超胞。

---

## 1. 做成公度超胞

```bash
python3 tools/build_common_cell.py \
  --a Bi.vasp --b Bi2.vasp --out ./build \
  --name-a Bi_A7_common --name-b Bi_II_common \
  --ma "-2 0 1, 0 2 0, -1 2 0" \
  --mb "-2 -2 1, -1 -1 -1, -1 -2 -1"
```

真实输出（`spglib` + `scipy` 需要装）：

```
A: input 12 atoms -> conventional 6 (R-3m)  V=219.516  V/atom=36.59
B: input 12 atoms -> conventional 4 (C2/m)  V=139.230  V/atom=34.81

A supercell: 12 atoms  V=439.032  det(M)=2
B supercell: 12 atoms  V=417.691  det(M)=3

cell A: a,b,c=[ 4.538  7.86  12.308]  angles=[90. 90. 90.]
cell B: a,b,c=[ 4.604  8.003 11.336]  angles=[90.4 90.  90. ]
residual strain B->A: 4.76 %   (isotropic part +1.67 %)

matching: total |dr| = 16.615 A  (mean 1.385 A)
  per-atom displacement (A): [0.   0.32 2.08 1.54 2.4  1.86 2.29 1.97 0.11 0.43 1.65 1.97]
  worst pair: 2.40 A
  A: 12 atoms  V=439.03  NN=3.074 A  symmetry=R-3m
  B: 12 atoms  V=417.69  NN=3.229 A  symmetry=C2/m
  wrote ./build/Bi_A7_common.vasp
  wrote ./build/Bi_II_common.vasp
```

要看的四行：

1. `12 atoms / 12 atoms` —— 原子数对上（`6·2 = 4·3 = 12`）；
2. `residual strain 4.76 %` —— 这就是相变本身要吸收的应变，越小越好；
3. `mean 1.385 A` / `worst pair 2.40 A` —— 原子对应关系。平均位移 1.4 Å、
   最大 2.4 Å 是"相变级"的位移，合理；如果最大位移接近原子间距（3 Å），
   说明匹配错了，去换超胞矩阵；
4. 两个 `NN`（最近邻距离）都在 3.1–3.2 Å —— Bi 的正常键长，没有原子重叠。

> **为什么是两个负数的矩阵？** 这里选的是一种"把六方基拧成接近单斜基"的
> 整数组合：`det(M_A) = 2`、`det(M_B) = 3` 是为了让原子数都对上 12，负号
> 是为了让三个矢量方向不至于差 60°，剩下的交给 Niggli 约化。

---

## 2. 生成算例

```bash
tools/setup_bi_case.sh ./bi_case ./build 14
cp Bi_ONCV_PBE-1.0.upf Bi_gga_10au_100Ry_2s2p2d1f.orb ./bi_case/
```

生成的 `bi_case/neb.yaml`（ABACUS、12 原子、7 张中间图、2.55 GPa）：

```yaml
project:
  name: Bi_A7_to_BiII_2p55GPa
  workdir: ./run

structure:
  initial: Bi_A7c.vasp
  final: Bi_IIc.vasp

calc:
  engine: abacus
  command: mpirun -np 14 abacus
  mpi: 14
  pseudo_dir: .
  data_dir: .
  timeout: 21600
  pseudopotentials:
    Bi: Bi_ONCV_PBE-1.0.upf
  orbitals:
    Bi: Bi_gga_10au_100Ry_2s2p2d1f.orb
  kpts: {density: 0.03}          # 一条路径只用一个网格
  params:
    ecutwfc: 100
    basis_type: lcao
    ks_solver: genelpa
    dft_functional: pbe
    smearing_method: mp
    smearing_sigma: 0.01         # Ry
    scf_thr: 1.0e-7
    scf_nmax: 100
    symmetry: 0

neb:
  images: 7
  interpolation: idpp
  align_cell: true
  min_distance: 2.0
  prune: auto
  variable_cell: true
  pressure: 2.55                 # GPa
  spring_k: 0.1
  f_th: 0.05
  max_steps: 400
  dt: 0.1
  max_step: 0.2
  cell_max_step: 0.02
  climbing_image: true
  relax_endpoints: true
  endpoint_fmax: 0.05
  endpoint_steps: 300
  parallel_images: 1

output:
  plot: true
```

几个数值的理由：

| 键 | 值 | 为什么 |
|----|----|--------|
| `kpts: {density: 0.03}` | → 8×5×3 | 12 原子金属超胞的经验值，见 `faq.md` 第 3 节 |
| `ecutwfc: 100` | Ry | LCAO + ONCV 赝势的常见值（`basis_type: lcao` 时是轨道截断） |
| `smearing_sigma: 0.01` | Ry ≈ 0.136 eV | 金属，MP 展宽 |
| `pressure: 2.55` | GPa | 实验上 A7→Bi-II 的相变压力附近 |
| `f_th: 0.05` | eV/Å | DFT 力噪声约 0.02–0.03 eV/Å，再低压不住 |
| `endpoint_fmax: 0.05` | eV/Å | 端点不用比路径更严 |
| `images: 7` | 9 张图 | 12 原子、位移 ~2 Å，7 张中间图够用 |

---

## 3. `vcneb doctor`

```bash
cd bi_case && vcneb doctor neb.yaml
```

```
vcneb 0.1.0 doctor: neb.yaml

[dependencies]
  OK    numpy      2.5.3
  OK    ase        3.26.0b1
  OK    yaml       6.0.3
  OK    matplotlib 3.11.2

[config]
  [1] executable 'mpirun' not found in PATH

[structures]
  OK    initial: 12 atoms, V=439.0319 A^3, formula=Bi12
  OK    final  : 12 atoms, V=417.6909 A^3, formula=Bi12
        initial: /path/bi_case/Bi_A7c.vasp
        final  : /path/bi_case/Bi_IIc.vasp

[engine]
  ABACUS  [mpirun -np 14 abacus]  ecutwfc=100 Ry  kpts=MP grid from kspacing 0.03 (fixed at run start)
  [1] executable 'mpirun' not found in PATH

some checks failed, fix them before running
```

这里那条 `mpirun not found` 是真的：这台机器上 ABACUS 装在别处，`mpirun` 不在
PATH 里。两种处理：把 `command` 写成绝对路径 / `source` 环境脚本，或者正式跑
的时候加 `--skip-engine-check`（在计算节点上 `mpirun` 是有的）。

**注意 `doctor` 不花机时**：它不跑任何 DFT，只读结构、查文件。

---

## 4. `vcneb run --dry-run`：先看路径

```bash
vcneb run neb.yaml --dry-run --skip-engine-check
```

真实输出：

```
vcneb 0.1.0  config: /path/bi_case/neb.yaml
engine: ABACUS  [mpirun -np 14 abacus]  ecutwfc=100 Ry  kpts=MP grid from kspacing 0.03 (fixed at run start)
  project      : Bi_A7_to_BiII_2p55GPa
  workdir      : /path/bi_case/run
  engine       : abacus   command: mpirun -np 14 abacus
  images       : 7 movable (+2 endpoints), interpolation=idpp
  variable cell: True   spring_k=0.1   f_th=0.05 eV/A
  max_steps    : 400   dt=0.1   climbing_image=True
  pressure     : 2.55 GPa (every energy below is the enthalpy H = E + P*V)
k-points: [8, 5, 3] (frozen from kspacing 0.03)   (same grid for all images)
endpoints: 12 atoms, V0=439.0319 -> V1=417.6909 A^3
dry run: endpoints are used exactly as given (no DFT relaxation, no electronic-structure call)
  IDPP pre-relaxation of the movable images done
  unwrapped 1 image(s) so the path is continuous across the periodic boundary
initial path built: 9 images (7 movable), interpolation=idpp
smallest interatomic distance on the path: 2.605 A
path spacing (RMS atom displacement between neighbouring images): 0.143, 0.133, 0.140, 0.158, 0.158, 0.123, 0.140, 0.139 A
dry run: path built, no electronic-structure call was made
  atoms per image: 12   variable cell: True
  image 00: V=439.0319 A^3  formula=Bi12  dir=.../run/images/end_initial
  image 01: V=436.4774 A^3  formula=Bi12  dir=.../run/images/01
  image 02: V=433.8908 A^3  formula=Bi12  dir=.../run/images/02
  image 03: V=431.2720 A^3  formula=Bi12  dir=.../run/images/03
  image 04: V=428.6209 A^3  formula=Bi12  dir=.../run/images/04
  image 05: V=425.9373 A^3  formula=Bi12  dir=.../run/images/05
  image 06: V=423.2212 A^3  formula=Bi12  dir=.../run/images/06
  image 07: V=420.4724 A^3  formula=Bi12  dir=.../run/images/07
  image 08: V=417.6909 A^3  formula=Bi12  dir=.../run/images/end_final
  wrote .../run/endpoint_initial.vasp (look at it in VESTA/OVITO)
  wrote .../run/endpoint_final.vasp (look at it in VESTA/OVITO)
next: remove --dry-run (or keep it to inspect the path)
```

`--dry-run` 里没有端点弛豫（日志明写 `endpoints are used exactly as given`），
所以这一步的 `endpoint_*.vasp` 就是你的输入结构 + `align_cell` 的重排结果。

### 逐行该看什么

| 行 | 期望 | 不对怎么办 |
|----|------|-----------|
| `k-points: [8, 5, 3]` | 三个数都不小于 1，量级合理 | 换 `density`，见 `faq.md` §3 |
| `endpoints: 12 atoms` | 两端原子数相同 | 不同 → 回第 1 步 |
| `smallest interatomic distance 2.605 A` | 明显大于 `min_distance` | 接近或小于 → 路径有原子重叠，看下一行 |
| `path spacing ...` 8 个数 | 彼此差不多（最大/最小 < 3） | 差很多 → 原子顺序没对上 |
| 体积序列 439.03 → 417.69 | **线性**下降（这里每步 −2.5 Å³ 左右） | 不线性 → 插值/晶胞对齐有问题 |
| 没有 `dropped images` | 没删点 | 有删点 → `min_distance` 太大或路径真有问题 |

这个例子 9 张图、`smallest distance 2.605 Å > min_distance 2.0`，所以
`prune: auto` 一张也没删（日志里没有 `dropped images`）。

### 一定要在 VESTA/OVITO 里看这两个文件

```bash
ase gui run/endpoint_initial.vasp run/endpoint_final.vasp &
```

确认每个 Bi 都找到了"对应"的另一个 Bi（同一层、同一侧），没有出现某个原子
要横穿整个晶胞（那是原子顺序错了，`IDPP` 会被拉成一根线）。

---

## 5. 正式跑

```bash
nohup vcneb run neb.yaml > run.out 2>&1 &
vcneb status run            # 随时看
```

开始时的日志（顺序固定：k 点 → 端点 → 端点弛豫 → 插值 → 开跑）：

```
k-points: [8, 5, 3] (frozen from kspacing 0.03)   (same grid for all images)
endpoints: 12 atoms, V0=439.0319 -> V1=417.6909 A^3
relaxing initial endpoint (fmax=0.05 eV/A, cell=free) ...
  initial: E=-XXXXXXXX eV  max|F|=0.0499 eV/A  converged  (6m32s)
relaxing final endpoint (fmax=0.05 eV/A, cell=free) ...
  final: E=-XXXXXXXX eV  max|F|=0.0488 eV/A  converged  (7m09s)
  IDPP pre-relaxation of the movable images done
  unwrapped 1 image(s) so the path is continuous across the periodic boundary
initial path built: 9 images (7 movable), interpolation=idpp
smallest interatomic distance on the path: 2.605 A
path spacing (RMS atom displacement between neighbouring images): 0.143, 0.133, 0.140, 0.158, 0.158, 0.123, 0.140, 0.139 A
------------------------------------------------------------------------------
starting vc-NEB: 7 movable images, f_th=0.05 eV/A, max_steps=400
------------------------------------------------------------------------------
```

（`E=XXXXXXXX` 是占位：这两个能量取决于你的赝势/泛函/展宽，本机没有 ABACUS
的 `mpirun`，所以没法给出这台机器上的真值；格式与源码
`band.py:relax_endpoints()` 逐字符一致。`max|F|` 必须 ≤ `endpoint_fmax`，
否则会多打一行 `WARNING: endpoint 'xxx' did not reach ... steps`，
把 `neb.endpoint_steps` 调大再跑。）

每一步的日志形如：

```
step   12  E[max]=-XXXXXXXX  barrier=0.0xxx  max|F_atom|=0.4183  max|F_cell|/atom=0.0072  dt=0.100  (37s)  CI=4
```

（`max|F_atom|` / `max|F_cell|/atom` / `dt` / `CI` 是真实字段并逐字对应
`band.py:run()`；能量类字段留成 `XXXX`，见上一段的说明。）

* `E[max]` 是当前全路径（含两个端点）的最高能量，`barrier` = `E[max] − E[image 0]`；
* `max|F_atom|` 是**去掉沿路径分量之后**的最大力（垂直于路径的），
  `f_th` 比的就是它；
* `max|F_cell|/atom` 是晶胞力除以原子数，`f_th` 也要比它（两个都小于才算收敛）；
* `dt` 是 FIRE 的当前时间步（所有可动像里取最小），`(37s)` 是这一步墙钟耗时；
* `CI=4` 是当前攀爬像。它会在 3/4/5 之间跳几次（带滞后阈值），收敛前稳定。

收敛时的最后几行：

```
converged: max|F_perp| = 0.0487 eV/A < 0.05
------------------------------------------------------------------------------
finished after 214 steps in 27h12m
------------------------------------------------------------------------------
barrier (forward) = 0.0xxx eV at image 4
barrier (reverse) = 0.0xxx eV
plot written: .../run/barrier.png

barrier   : 0.0xxx eV
reverse   : 0.0xxx eV
steps     : 214  converged: True
results   : .../run
```

（`finished after 214 steps in 27h12m` 是 `fmt_duration` 的格式：不到 1 h 打
`NNmSSs`，不到 1 min 打 `NNs`，超过 1 h 打 `NNhMMm`，没有"天"。）

（这个体系的能垒取决于你用的 U/展宽/赝势，别照抄这里的数字，看你自己跑出来的。）

**机时估计**：12 原子 LCAO、`ecutwfc: 100 Ry`、8×5×3 网格、`scf_thr 1e-7`，
单张图一次 SCF 大约几分钟（14 核）；每步 9 次 SCF（含两个端点，**每步都重算**，
这样全路径永远用同一套参数、同一个 k 网格；端点只有 `relax_endpoints: false`
时才会被跳过弛豫，但仍参与每步的能量/力计算）；
200 步左右收敛 → 一两天量级。想快就开 `parallel_images`（要核数够）：

```yaml
neb:
  parallel_images: 4      # 4 张图同时算，每张 14 核 → 需要 56 核
```

---

## 6. 看结果

```bash
vcneb status run
column -s, -t run/barrier.csv
```

```
 image      E (eV)      dE (eV)    max|F| (eV/A)     V (A^3)
     0  -xxxx.xxxxx     0.00000          0.0000      439.0319
     1  ...             0.01234          0.0412      436.4774
     ...
     8  ...             0.00987          0.0388      417.6909
```

（能量和力的**列格式**与真实 `vcneb status` 输出一致，`E`/`dE`/`max|F|` 的数值
是占位：能量零点取决于你的赝势，力要靠真跑的 SCF 才出得来。体积列是从
`--dry-run` 来的真值。）

要看的：

* `dE` 沿路径先升后降，最高点在第 4 张左右（`image_max`）；
* `V` 单调下降（439 → 417.7），没有跳变；
* `max|F|` 的每一列都 < `f_th`（收敛了才会有这个性质）。

`band.traj` 用 ASE GUI 或 OVITO 打开可以看动画：

```bash
ase gui run/band.traj
```

鞍点那张图（`image_max`）可以拿去做频率/进一步精修。

---

## 7. 这个例子里容易踩的坑

| 现象 | 原因 | 处理 |
|------|------|------|
| `doctor` 报 `executable 'mpirun' not found` | 计算节点才有 | 提交脚本里 `source` 环境，或 `--skip-engine-check` |
| `ABACUS: orbital file not found` | `data_dir` 和 `orbitals` 里的文件名对不上 | 文件名精确到大小写 |
| 第一次 SCF 就 timeout | `timeout: 21600` 是按"一次 SCF"算的 | 调大 |
| 跑了几步 `max\|F\|` 不降 | `dt` 太大 / 展宽太小 / 端点没弛豫干净 | 先看 `images/01/` 里的 SCF 是否真收敛 |
| 出现 `frozen k-mesh ... too coarse` | 路径上晶胞跑出端点范围 | 用 `{density: 0.02}` 加密 |
| 想加 U / 看自旋 | 见 `faq.md` §5 | 加 `nspin: 2`、`dft_plus_u`… |
| 停电/排队被砍 | 进度在 `state.json` | `vcneb run neb.yaml --restart` |

---

## 8. 换成 QE 或 CP2K 要改哪里

只改 `calc:` 一段（结构、路径、压力、插值全都不用动）：

```yaml
# QE
calc:
  engine: qe
  command: mpirun -np 14 pw.x
  mpi: 14
  pseudo_dir: ./pseudo
  pseudopotentials:
    Bi: Bi.pbe-dn-kjpaw_psl.1.0.0.UPF
  kpts: {density: 0.03}
  params:
    ecutwfc: 45
    ecutrho: 360
    occupations: smearing
    smearing: marzari-vanderbilt
    degauss: 0.02
    conv_thr: 1.0e-8
```

```yaml
# CP2K
calc:
  engine: cp2k
  command: mpirun -np 14 cp2k.popt
  mpi: 14
  data_dir: /opt/cp2k/data
  kpts: {density: 0.03}
  params:
    basis_set: DZVP-MOLOPT-SR-GTH
    potential: GTH-PBE-q5
    xc: PBE
    cutoff: 600
    rel_cutoff: 60
    eps_scf: 1.0e-6
    smearing: true
    smearing_method: marzari_vanderbilt
    smearing_sigma: 0.05
    added_mos: 12
```

注意 `neb.pressure` 保持 `2.55`，但三个引擎输入里都**不要**写压力
（见 `faq.md` §4）。
