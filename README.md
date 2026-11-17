# Chapter 4: Medium Scale DSGE Model

## 简介

《高级宏观经济学》第四章配套 MATLAB/Dynare 代码，复现中等规模DSGE模型的贝叶斯估计、七类冲击的IRF、货币政策传导机制对比、方差分解、历史分解及模型矩与数据矩比较。

完整模型设定、参数估计和定量分析内容参见 [lecture note](medium_dsge_lecture_note.pdf)。

## 目录结构

- `codes/`：主程序、Dynare模型、估计、分析及绘图函数。
- `data/`：美国原始数据、观测序列和数据处理程序，来源与单位见Excel。
- `estimates/`：课件使用的独立估计快照、后验众数、保留的后验抽样和诊断结果。
- `reference_results/`：课件定量结果，供快速绘图和复现核对。
- `medium_dsge_lecture_note.pdf`：本章讲义。
- `output/`：运行时生成的图表、数值结果及Dynare工作文件。

本包不依赖其他章节或父目录文件。`estimates/`与`reference_results/`仅供程序读取，重新运行分析或估计均不会覆盖。发布结果与自行重新估计的输出分别保存在`output/published/`和`output/reestimated/`。

## 运行环境

- MATLAB R2023a、Dynare 7.1。
- Statistics and Machine Learning Toolbox用于先验和后验密度绘图。Parallel Computing Toolbox可选，用于并行估计。
- 中文绘图使用Microsoft YaHei，其他系统可在绘图函数中替换字体。
- 已包含所需数据和估计结果，运行无需联网。

## 使用方法

将MATLAB当前文件夹切换到`codes/`，添加本机Dynare路径，例如：

```matlab
addpath('C:/Program Files/dynare/7.1/matlab');
DSGE_main;                       % 基于随包估计结果重新计算全部定量分析并绘图
```

首次使用推荐直接读取随包估计结果，无需重跑MCMC。IRF报告2000组后验抽样下的逐期均值与90% HPD区间，机制比较使用同一批抽样和紧缩性货币政策冲击。历史分解、方差分解和模型矩在后验均值参数处计算，参数表同时保存后验众数。

其他入口：

```matlab
DSGE_main('figures');                    % 仅用随包数值结果重绘课件图像
DSGE_main('check');                      % 核验稳态与BK条件
DSGE_main('data');                       % 从Excel重建观测序列至output/data/
DSGE_main('estimate');                   % 从头估计，耗时较长
DSGE_main('analyze','reestimated');       % 使用自己的估计结果进行分析
DSGE_main('extend',[],100000);            % 自己的每条MCMC链增加100000次抽样
```

图像保存为`output/published/figures/`或`output/reestimated/figures/`中的PNG，数值结果保存于相应`results/`目录。重新估计默认4条链、每链100000次正式抽样，并包含独立的众数搜索及试运行阶段。已有自行估计结果时，`estimate`直接重建分析；重新从头估计可先将`output/reestimated/`另存备份。
