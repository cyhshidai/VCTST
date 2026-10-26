# 代码地图：每个命令、插值、删点到底动了哪几行

本文只回答一件事：**你敲的那条命令，最后落到了哪个函数、改了什么、写了哪些文件。**
物理背景在 [`faq.md`](faq.md)，按步操作的算例在 [`example-bi.md`](example-bi.md)。

> 行号对应本文写作时的代码（`vcneb/` 共 3802 行）。改过代码行号就会漂，
> 文末「§9 行号漂了怎么找」给了重新定位用的 `grep`。

## 0. 一张图

```
neb.yaml ──load_config──> Config ──validate──> [有问题 → 打印，退出 1]
                                     │
                     get_endpoints ──┴──> 两个 Atoms（末态原子顺序已对到初态）
                                     │
                      make_engine + plan_kpoints  → k 点定死，全带共用
                                     │
                         Band(initial, final, cfg, engine, log)
                                     │
         ┌───────────────────────────┴────────────────────────────┐
    relax_endpoints                                     interpolate（4 步 + 删点）
         └───────────────────────────┬────────────────────────────┘
                              run()：NEB 迭代 → write_outputs()
```

## 1. 命令 → 代码

| 命令 | 入口 | 干什么 | 写文件 |
|------|------|--------|--------|
| `vcneb init` | `cli.py:41 cmd_init` | 拷一份带注释的模板（模板文本来自 `config.py:580 default_template`） | `<out>.yaml` |
| `vcneb keys` | `cli.py:188 cmd_keys` | 把 `DEFAULTS`（`config.py:24`）渲染出来 | 无 |
| `vcneb doctor` | `cli.py:56 cmd_doctor` | `validate()` + `engine.check()`，不花机时 | 无 |
| `vcneb run` | `cli.py:229 cmd_run` | 见 §2 | `run/` 全套 |
| `vcneb status` | `cli.py:397 cmd_status` | 只读 `state.json` / `summary.json` / `barrier.csv` | 无 |
| `vcneb plot` | `cli.py:388 cmd_plot` → `band.py:1040 plot_band` | 画带 | `barrier.png` |

所有命令行开关（`--restart --force --dry-run --max-steps --skip-engine-check -v`）
都在 `cli.py:430 build_parser` 里定义，没有藏到别处的"隐式开关"。

## 2. `vcneb run` 的时序（`cli.py:229`）

1. `load_config`（`config.py:278`）：YAML 覆盖 `DEFAULTS`（`config.py:24`），相对路径转绝对（`Config._resolve_paths`，`config.py:145`）。
2. `cfg.validate()`（`config.py:159`）：有问题就打印清单并 `return 1`，一个 SCF 都不跑。
3. 建 `run/`、`run/logs/`、`run/images/`；`Logger`（`util.py:127`）同时写 `run/logs/vcneb.log`。
4. state 检查：给了 `--restart` 但没有 `state.json` → 报错；有 `state.json` 却没给 `--restart` → 报错（给 `--force` 则先 `_clean_workdir`，`cli.py:213`）。
5. `get_endpoints`（`structs.py:124`）→ `_match_endpoints`（`structs.py:60`）把末态原子顺序对到初态。
6. `make_engine`（`engines.py:1281`）+ `engine.plan_kpoints([c0, c1])`（`engines.py:97`）→ **k 点在这里定死**（见 §6）。
7. `engine.check()`：查可执行文件 / 赝势 / data 目录（`engines.py:41 BaseEngine`、`165 QEEngine`、`940 CP2KEngine`、`1181 AbacusEngine`）；`--dry-run` 时只 warning 不退出。
8. `Band(...)`（`band.py:223`）。
9. 没有 `--restart` 且不是 dry-run：`relax_endpoints`（`band.py:371`）→ `interpolate`（`band.py:487`）→ `save_state`（`band.py:887`）。
   * `--restart`：只 `load_state`（`band.py:902`），**不重新插值、不重新删点**；
   * `--dry-run`：`relax_endpoints` + `interpolate` + `_dry_run_report`（`cli.py:328`）后 `return 0`，**不写 `state.json`**。
10. `band.run()`（`band.py:753`）：每步 `evaluate_all`（`band.py:352`）→ `neb_forces`（`band.py:675`）→ `_optimizer_step`（`band.py:745`）→ 判 `max|F_atom| < f_th` 且 `max|F_cell|/atom < f_th`。
11. `write_outputs`（`band.py:941`）→ `summary.json` / `barrier.csv` / `barrier.png` / `band.traj`。
12. `_write_resolved_config`（`cli.py:375`）→ `config.resolved.yaml`；`_maybe_clean_images`（`cli.py:350`）按 `neb.keep_intermediate_images` 清理 `images/`。
13. 返回码：`0` 收敛 / `2` 没收敛（可 `--restart`）/ `1` 配置或引擎问题 / `130` Ctrl-C。

## 3. 插值：`interpolate()`（`band.py:487`）到底改了什么

`neb.interpolation` 只决定第 3 步做不做，五步顺序固定：

| # | 代码 | 改的东西 | 物理结构会变吗 |
|---|------|----------|----------------|
| 1 | `_align_final_cell` `band.py:433` | 末态**晶胞写法**：枚举 6 种轴置换 × 8 种符号翻 = 48 种，按 `score = ‖sym(F)−I‖ + RMS(原子位移)/L`（`F = c0⁻¹·c1'`）取最小，只改 `states[-1].cell` | 不变（同一组原子、同一个晶胞，只换索引/符号） |
| 2 | 线性插值 `band.py:496–506` | 分数坐标做最小像差分 `dfrac −= round(dfrac)`，按 `t = k/(n_free+1)` 混合 `frac` 与 `cell` | 生成 `n_free` 张新图；这一步本身会插出不物理的结构 |
| 3 | `_idpp_relax` `band.py:611` | 逐图 IDPP + `MDMin`，目标距离矩阵 `d1 + i·(d2−d1)/n`，`fmax=neb.idpp_fmax`（0.02）、`steps=neb.idpp_steps`（200），**晶胞固定** | 只动原子；失败只 WARN 并保留线性路径 |
| 4 | `_unwrap_images` `band.py:636` | 沿带 `x_i = x_{i−1} + MIC(x_i − x_{i−1})`，把跨边界的原子接回来 | 不变（同一组周期像，只是坐标连续了） |
| 5 | `_prune_close_images` `band.py:557` | 可能**删图**（见 §4） | 图数变少 |

第 2–4 步做完后统一打印（`band.py:515–550`）：

```
initial path built: N images (M movable), interpolation=idpp
smallest interatomic distance on the path: 2.605 A
path spacing (RMS atom displacement between neighbouring images): 0.143, 0.133, ... A
```

这三行就是 README §8.5 那段。`path spacing` 的均匀性检查也在这里：最大值超过
最小值 3 倍会 WARN，若是因为删点导致的，会自动改口成 `(uneven on purpose: ...)`。

为什么不用 `ase.neb.idpp_interpolate`：它拒绝"晶胞沿路径变化"的带，而变胞 NEB
全是这种。IDPP 的目标函数只用距离矩阵，所以逐图单独做完全等价，又允许晶胞不同
（`band.py:611` 的 docstring 写了这段话）。

## 4. 删点：`_prune_close_images()`（`band.py:557`）

前提：`neb.min_distance > 0`；只检查**可动图**（`range(1, n_images−1)`），端点不动。
判据：该图最小像最近间距 `d_min < neb.min_distance`（`_image_min_distance`，`band.py:552`）。

| 条件 | 代码行为 | 你会看到 |
|------|----------|----------|
| `min_distance: 0`（或没写）、或 `n_free < 1` | `return`（`band.py:569`） | 什么都不删，只有 §8.5 那条 `< 1.2 Å` 的 WARN |
| `prune: false` | `return` | 同上 |
| 没有坏图 | `return` | 什么都不做 |
| `prune: error` + 有坏图 | `raise VcnebError`（`band.py:579`） | 列出图号和间距，`vcneb run` 返回 1 |
| `prune: auto` + 有坏图 + 至少剩 1 张 | 删图（`band.py:591`）、`self.pruned=[...]`、`n_free=len(keep)−2` | `WARNING: dropped images ...`；剩不到 3 张再补一条 WARN |
| `prune: auto` + 全删光 | `raise VcnebError`（`band.py:584`） | `ERROR: every interpolated image is unphysical ...` |

要点：

* 删点发生在 `interpolate()` 内部、`band.run()` 之前（`band.py:521`），所以**被删的图不会再回来**，带会变粗（`path spacing` 在被删处变大，日志标 `uneven on purpose`）；
* 图号与 `pruned` 存进 `state.json.pruned_images`（`band.py:893`），`--restart` 原样恢复（`band.py:929`），不会重新插值/重新删点；
* 想把删掉的图要回来：改 `neb.min_distance`，删掉 `run/state.json`，重跑。没有别的开关。

四段日志（删一部分 / 删太多 / `error` / 全删光）可以用
`bash examples/prune_demo.sh` 一次复现，不需要任何 DFT 程序。

## 5. `state.json` ↔ 代码

`save_state`（`band.py:887`）/ `load_state`（`band.py:902`）：

| 字段 | 含义 | 续算时用在哪 |
|------|------|--------------|
| `step` | 已经走完的步数 | 日志 `restarting from step N` |
| `images[].cell` / `.pos` | 每张图的晶胞与坐标 | 直接接上，不重新插值 |
| `images[].energy` | 上次的能量 | 续算第一屏就能画带 |
| `n_movable` | 可动图数 | 和当前配置不一致时 WARN `keeping the saved band` |
| `pruned_images` | 被删的图号 | 恢复 `self.pruned` |
| `optimizer_dt` | 每张图的步长 | 恢复 FIRE 的自适应步长 |
| `variable_cell` / `n_atoms` | 守门 | 原子数不一致直接报错，不猜 |

## 6. k 点为什么"冻结"

`BaseEngine.plan_kpoints`（`engines.py:97`）：`calc.kpts: {density: X}` 时，对初末
两个晶胞各算一套 MP 网格（`grid_for_cell`，`engines.py:91`），取**逐轴最大**，
整个 run 共用这一套。

理由（源码注释 `engines.py:73–81`）：沿路径晶胞在变，逐图重算网格会让 k 点数
按整数跳变，每次跳变都在总能上留一个小台阶 = 假力 = 假能垒。冻结网格让残差
沿路径连续变化（残差多大要自己做 k 点收敛测试）。

某个中间图其实需要更密网格时，`kpoints_too_coarse`（`engines.py:117`）会把轴报出来。

## 7. 引擎层：`calc.params` 落到哪

| 引擎 | 类 | 输入怎么来的 | `calc.params` 有白名单吗 |
|------|----|--------------|--------------------------|
| QE | `QEEngine` `engines.py:165` | ASE `espresso` 模板 + `input_data` | 直接并进 `input_data` |
| CP2K | `CP2KEngine` `engines.py:940` | `_render_cp2k_input` `engines.py:745` | 有：`_CP2K_PARAM_KEYS`（`engines.py:588`），表外的键 `vcneb doctor` 会指出来 |
| ABACUS | `AbacusEngine` `engines.py:1181` | 生成 `INPUT` / `STRU` / `KPT` | 见 README §7.3 |

CP2K 的结果回读：`_parse_cp2k_output`(`252`)、`_parse_cp2k_forces`(`297`)、
`_parse_cp2k_stress`(`361`)、`_parse_cp2k_mulliken`(`405`) —— `barrier.csv` 里的
`total_spin` / `net_charge` 就来自最后一个。

焓：`EnthalpyCalculator`（`band.py:40`）+ `cell_force_from_stress`（`band.py:89`），
只有 `neb.pressure` 非零时才包在真正计算器外面。

## 8. 测试在守什么

| 文件 | 守的 | 怎么跑 |
|------|------|--------|
| `tests/test_interpolation_guards.py` | 7 个：`align_cell` 开/关、删点三态（partial / error / all）、`min_distance: 0` 关掉、`state.json` 往返后带大小不变 | `pytest` |
| `tests/test_pressure.py` | 2 个：`H = E + P·V` 与 ASE `ExpCellFilter` 对拍、加压后晶胞确实变小 | `pytest` |
| `tests/test_cellforce.py` | `F_cell = −V·σ·cell⁻ᵀ` 与有限差分对拍（fractional 固定） | `python3 tests/...`（有 `main()`） |
| `tests/test_cellforce_basis.py` | 同上，Cartesian 基组那一路 | 同上 |
| `tests/t_plumbing_pop.py` | 布居 → `barrier.csv` → `summary.json` → 图 | 同上 |
| `tests/cp2k_u.yaml`、`cp2k_u_mini.yaml` | 不是测试，是 `vcneb doctor` 的夹具（UKS + DFT+U + MULTIPLICITY + Mulliken） | `vcneb doctor tests/cp2k_u.yaml` |

## 9. 行号漂了怎么找

```bash
grep -n 'def interpolate\|def _align_final_cell\|def _idpp_relax\|def _unwrap_images\|def _prune_close_images\|def _image_min_distance' vcneb/band.py
grep -n 'def cmd_' vcneb/cli.py
grep -n 'def plan_kpoints\|def grid_for_cell\|class .*Engine\|def make_engine' vcneb/engines.py
grep -n '^DEFAULTS\|def validate\|def load_config\|def default_template' vcneb/config.py
```

想确认**这一版**的默认值，别读文档，读 `vcneb keys` —— 它就是 `DEFAULTS` 的渲染。
