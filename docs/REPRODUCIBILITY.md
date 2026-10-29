# 复现与核验

## 环境与入口

MATLAB R2024a，Optimization Toolbox 24.1。基础输入数据内嵌于run_motor_modeling.m，不依赖私有数据或联网服务。

在仓库根目录执行verify_motor_results可核验基础快照；执行run_motor_modeling及paper_figures重算基础结果；执行run_revision_experiments重算补充结果。补充入口读取results/results.mat中的参数和已知可行初值，不能只下载单独的.m文件。

## 结果对应关系

- results/summary.csv对应基础情景，论文表8至表12、附录A使用该快照。
- revision_results/revision_summary.csv对应12套补充计划，论文表13至表16使用该快照。
- revision_results/margin88_schedule.csv为名义5000件、88℃控制上限计划。
- results/robust_worst_schedule.csv为5500件联合最坏工况计划，不应与88℃计划混用。
- 各目录的MAT文件保存计划、参数和求解器输出，CSV便于无需MATLAB直接查看。

## 核验标准

基础实验有10套有效计划通过独立积分与交付核验，名义方案直接在最坏工况回放是明确保留的不可行反例。27组鲁棒参数回放均合格。

12套补充计划均完成5000件，并核验分阶段交付、电机最短10分钟启停、风机最短5分钟持续时间、费用及温度。补充实验ODE积分与解析递推最大差约6.65×10⁻⁹℃。这些检查均在原计算阶段执行；上传仓库时未重新运行全部优化。

HiGHS相对间隙输出为比值，汇总的Gap_percent已乘100。退出标志2表示求解被时间限制中止但取得了整数可行解；必须结合独立可行性检查与上下界理解结果。

## 时间网格注意事项

改变控制间隔时，最短持续时间按分钟重新换算，不能一直写成两个区间。代码在2.5分钟网格下分别要求电机至少4段、风机至少2段。解析热递推消除了常值区间的积分离散误差，但没有消除控制网格对最优值的影响。

日志保留为batch.log、verification.log、revision_batch.log，其中本机用户名和工作目录已作替换。仓库不包含机器配置、凭据、工具缓存或内部排版预览。
