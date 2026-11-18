# wave-ice SPH 程序（二维水-冰耦合）

二维弱可压缩 SPH 程序，用于模拟波浪与浮冰（及冰的弯曲、压缩破坏）的相互作用。
源码基于 HEU 早期溃坝 SPH 程序改造（`dambreaking_sph.f` 第 3 行注释：Improved by Zheng Xing in HEU 2013），
后续加入了冰的弹塑性破坏与本构、冰-水双向耦合、造波/消波与 OpenMP 并行。

## 目录结构

```
gutishuntaili.sln                 Visual Studio 解决方案
gutishuntaili/
  gutishuntaili.vfproj            工程文件（Intel Fortran）
  dambreaking_sph.f               主程序：几何/粒子初始化、参数、主循环入口
  menu.inc                        全局公共块与物理参数（E、rho0、dt、bs、s_rad 等）
  timeintegration.for             时间推进（预测-校正、冰/水双时间步、输出控制）
  singlestep.for                  冰相单步：状态方程 -> sp（邻域）-> interforce
  singlestep1.for                 水相单步：造波活塞 -> 虚粒子 -> sp1 -> interforce1
  interforce.for                  冰相内力：SPH_SFDI 修正、D-P 塑性、软化、断裂屏蔽、冰-水作用
  interforce1.for                 水相内力：压力梯度、人工粘性、自由面/边界压力外推、冰反力累积
  sp.for / sp1.for                链式网格邻域搜索（冰 / 水）
  kernel.f                        核函数：三次样条、五次样条、C4 型（kernelw）及导数版本
  wavemaker.f                     Goring 推板造波（当前未启用，见 singlestep1.for）
  s_grid.f / rn_grid.f90          邻域网格模块
  time_print.f90 / time_elapsed.f90  时间戳与计时
```

## 编译

需要 Intel Fortran（ifort/ifx）+ Visual Studio：

- 打开 `gutishuntaili.sln`，选 `x64` + `Release`，生成解决方案；
- 命令行示例（在 `gutishuntaili\` 下）：
  `ifx /Qopenmp /O2 /fpp dambreaking_sph.f kernel.f interforce.for interforce1.for singlestep.for singlestep1.for sp.for sp1.for timeintegration.for wavemaker.f s_grid.f rn_grid.f90 time_print.f90 time_elapsed.f90`

## 运行

1. 双击/命令行运行 `gutishuntaili.exe`；
2. 按提示输入总时间步数 `tsn`（`dt = 5e-6 s`，`tsn = 200000` 约对应 1 s 物理时间）；
3. 结束后可选择继续计算（输入 1）或退出（输入 0）。

## 输出（均为 Tecplot/ASCII，已在 .gitignore 中忽略）

| 文件 | 内容 |
|---|---|
| `animate.dat` | 水粒子场动画（最大，约几十 MB） |
| `ice.dat` | 冰粒子场（坐标、速度、应力、塑性应变、破碎标记） |
| `freesurface.dat` | 自由面粒子 |
| `force.dat` | 指定冰粒子的位移与应变率时程（列名沿用旧名，并非力） |
| `wave1/2/3.dat` | x=9 m / 14 m / 22 m 处自由面高程时程（x=22 超出当前 16 m 计算域） |
| `eng.dat` | 冰相势能/动能/总能量 |

## 已知硬编码（改造前需先处理）

- `menu.inc`：全部物性与数值参数（`E=1.3865e8`、`rho0=900`、`rho01=1000`、`dt=5e-6`、`bs=45`、`s_rad=3.1`、`deltf=0.1` 等）；
- `dambreaking_sph.f`：计算域 `x0=16, y0=0.6`、水深 `ysl=0.4`、冰板 `3 m x 0.02 m`、粒子数 `600x4`（冰）与 `3180x80`（水）；
- `sp.for` 第 34/76 行：邻域搜索被限制在 `x in [5,11], y >= 0.3` 的冰区；
- `singlestep1.for`：造波参数 `wh=0.025 m`、`period=1.2 s`；消波区见 `timeintegration.for`（`lf=16, ld=3`）；
- `interforce.for` 第 963 行、`interforce1.for` 第 1179/1236 行：δ-SPH 密度耗散项（`deltf`）目前为注释状态。

## 备注

- 仓库只跟踪源码、工程文件与说明文档；`*.dat`、编译中间件、可执行文件不入库。
- 数据文件可由程序重新生成，需要保存某次算例结果时请单独打包留存，不要直接提交。
