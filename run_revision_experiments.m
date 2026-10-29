function experiments=run_revision_experiments()
%RUN_REVISION_EXPERIMENTS 三项补强实验；保持原始论文结果不变。
% 入口：experiments=run_revision_experiments。需要 Optimization Toolbox。
root=fileparts(mfilename('fullpath')); out=fullfile(root,'revision_results');
if ~isfolder(out), mkdir(out); end
loaded=load(fullfile(root,'results','results.mat')); r=loaded.results;
th=r.thermal; beta=r.fit.beta; c=r.cfg;
c.speeds=(500:50:1500)'; c.maxTime=60; c.gapTolerance=1e-4;
names={'joint90','margin89','margin88','margin87','margin85',...
    'fixed_speed_fan_optimized','speed_optimized_fan_on',...
    'time10','time2p5','speed100','speed25','joint90_long'};
experiments=cell(numel(names),1);
for j=1:numel(names)
    cfg=c;
    switch names{j}
        case 'margin89', cfg.tempLimit=89;
        case 'margin88', cfg.tempLimit=88;
        case 'margin87', cfg.tempLimit=87;
        case 'margin85', cfg.tempLimit=85;
        case 'fixed_speed_fan_optimized'
            cfg.fixedSpeed=[1250*ones(400/c.dt,1);zeros(80/c.dt,1)];
        case 'speed_optimized_fan_on', cfg.alwaysFan=true;
        case 'time10', cfg.dt=10;
        case 'time2p5', cfg.dt=2.5;
        case 'speed100', cfg.speeds=unique([(500:100:1500)';1250]);
        case 'speed25', cfg.speeds=(500:25:1500)';
        case 'joint90_long', cfg.maxTime=180;
    end
    % 从已通过原始独立核验的计划构造整数可行初值。
    if strcmp(names{j},'time10') || isfield(cfg,'fixedSpeed')
        n=[1250*ones(400/cfg.dt,1);zeros(80/cfg.dt,1)]; fan=ones(size(n));
    elseif strcmp(names{j},'speed100')
        n=r.nominal.n; fan=r.nominal.fan;
    elseif strcmp(names{j},'joint90_long')
        n=experiments{1}.n; fan=experiments{1}.fan;
    else
        rep=5/cfg.dt; n=repelem(r.fine.n,rep); fan=repelem(r.fine.fan,rep);
    end
    seed=run_motor_modeling('simulate',n,fan,cfg,th,beta);
    if ~seed.feasible || (isfield(cfg,'alwaysFan') && cfg.alwaysFan)
        fan=ones(size(n)); seed=run_motor_modeling('simulate',n,fan,cfg,th,beta);
    end
    assert(seed.feasible,'Warm start must be feasible');
    cfg.seed=struct('n',n,'fan',fan,'temp',seed.temp);
    s=run_motor_modeling('solve',cfg,th,beta,names{j});
    s.cfg=rmfield(s.cfg,'seed');
    s=independentCheck(s,th,beta);
    experiments{j}=s;
    run_motor_modeling('export',s,out,names{j});
    save(fullfile(out,'revision_checkpoint.mat'),'experiments');
    publish(experiments(1:j),out,r);
end
save(fullfile(out,'revision_results.mat'),'experiments');
makeRevisionFigures(experiments,out);
fprintf('Revision experiments complete: %d verified plans\n',numel(experiments));
end

function s=independentCheck(s,th,beta)
% 独立用ode45积分，用数值积分累计电费以外的能量代数重算。
c=s.cfg; N=numel(s.n); T=zeros(N+1,1); T(1)=c.initialTemp;
P=zeros(N,1); q=zeros(N,1);
for k=1:N
    n=s.n(k); torque=c.loadFactor*(30+.02*n); x=n/1000; y=torque/60;
    L=c.lossFactor*([1 x y x*x x*y y*y]*beta);
    if n==0, L=0; end
    R=2*th.R; if s.fan(k), R=th.R*c.coolFactor; end
    [~,traj]=ode45(@(~,v)(L-(v-c.ambient)/R)/th.C,[0 c.dt*60],T(k),...
        odeset('RelTol',1e-9,'AbsTol',1e-10));
    T(k+1)=traj(end); P(k)=2*pi*n*torque/60+L+c.fanPower*s.fan(k);
    q(k)=.01*n*c.dt;
end
on=s.n>0; starts=on & ~[false;on(1:end-1)];
time=(0:N-1)'*c.dt; price=.9*ones(N,1);
price(time<120)=.65; price(time>=120 & time<240)=1.1; price(time>=240 & time<360)=.45;
energy=P*c.dt/60000+.03*starts;
assert(abs(sum(energy)-s.energy)<1e-8);
assert(abs(sum(price.*energy)-s.cost)<1e-8);
s.odeMaxError=max(abs(T-s.temp)); assert(s.odeMaxError<1e-6);
assert(max(T)<=c.tempLimit+1e-6);
Q=[0;cumsum(q)]; idx=[120 240 360 480]/c.dt+1;
assert(all(Q(idx)>=[1000;2300;3600;5000]-1e-6));
edge=[1;find(diff(on)~=0)+1;N+1]; len=diff(edge)*c.dt; states=on(edge(1:end-1));
assert(all(len(states)>=10-1e-8));
for k=2:numel(len)-1, assert(states(k) || len(k)>=10-1e-8); end
edge=[1;find(diff(s.fan)~=0)+1;N+1]; assert(all(diff(edge)*c.dt>=5-1e-8));
assert(all(s.fan==0 | s.fan==1));
assert(all(s.n==0 | ismember(s.n,c.speeds)));
if isfield(c,'fixedSpeed'), assert(all(s.n==c.fixedSpeed)); end
if isfield(c,'alwaysFan'), assert(all(s.fan==1)); end
s.fanSwitches=sum(diff(s.fan)~=0); s.speedChanges=sum(diff(s.n)~=0);
s.fanMinutes=sum(s.fan)*c.dt;
% HiGHS absolute gap is U-B. Retain the solver-reported certificate.
assert(isfield(s.solverOutput,'absolutegap'));
s.lowerBound=s.solverObjective-s.solverOutput.absolutegap;
s.verified=true;
end

function publish(sets,out,r)
m=zeros(numel(sets),14); labels=strings(numel(sets),1);
for k=1:numel(sets)
    s=sets{k}; labels(k)=string(s.label);
    m(k,:)=[s.cfg.tempLimit s.cfg.dt numel(s.cfg.speeds) s.cost s.lowerBound ...
        100*s.relativeGap s.energy s.maxTemp s.fanMinutes s.fanSwitches ...
        s.speedChanges s.solveSeconds s.odeMaxError s.exitflag];
end
t=array2table(m,'VariableNames',{'Limit_C','Step_min','SpeedCount','Cost_yuan',...
    'LowerBound_yuan','Gap_percent','Energy_kWh','MaxTemp_C','FanMinutes',...
    'FanSwitches','SpeedChanges','SolveSeconds','ODEError_C','Exitflag'});
t=addvars(t,labels,'Before',1,'NewVariableNames','Experiment');
writetable(t,fullfile(out,'revision_summary.csv')); disp(t);
fid=fopen(fullfile(out,'revision_summary.json'),'w','n','UTF-8');
fprintf(fid,'%s',jsonencode(struct('baselineCost',r.baseline.cost,...
    'baselineEnergy',r.baseline.energy,'experiments',table2struct(t)),'PrettyPrint',true)); fclose(fid);
end

function makeRevisionFigures(sets,out)
font='Microsoft YaHei';
f=figure('Color','w','Position',[80 80 1200 560]); t=tiledlayout(f,1,2,'TileSpacing','compact');
ids=[5 4 3 2 1]; x=cellfun(@(s)s.cfg.tempLimit,sets(ids));
ub=cellfun(@(s)s.cost,sets(ids)); lb=cellfun(@(s)s.lowerBound,sets(ids));
ax=nexttile(t); errorbar(ax,x,ub,ub-lb,zeros(size(ub)),'o-','LineWidth',1.6);
xlabel(ax,'控制温度上限 / ℃'); ylabel(ax,'班次电费 / 元'); grid(ax,'on');
title(ax,'点为可行费用 下端为求解下界');
ax=nexttile(t); plot(ax,x,cellfun(@(s)s.fanMinutes,sets(ids)),'s-','LineWidth',1.6);
xlabel(ax,'控制温度上限 / ℃'); ylabel(ax,'风机开启时间 / min'); grid(ax,'on');
title(ax,'控制上限变化的冷却需求');
set(findall(f,'Type','axes'),'FontName',font,'FontSize',14);
exportgraphics(f,fullfile(out,'fig5_temperature_margin.png'),'Resolution',200); savefig(f,fullfile(out,'fig5_temperature_margin.fig'));
f=figure('Color','w','Position',[80 80 1100 540]); ax=axes(f);
ids=[8 1 9 10 11]; labels={'10min 50转','5min 50转','2.5min 50转','5min 100转','5min 25转'};
ub=cellfun(@(s)s.cost,sets(ids)); lb=cellfun(@(s)s.lowerBound,sets(ids));
errorbar(ax,1:5,ub,ub-lb,zeros(size(ub)),'o','LineWidth',1.8);
xticks(ax,1:5); xticklabels(ax,labels); xlim(ax,[.5 5.5]); grid(ax,'on');
ylabel(ax,'班次电费 / 元'); title(ax,'不同控制网格的费用上下界 均限时60秒');
set(ax,'FontName',font,'FontSize',14);
exportgraphics(f,fullfile(out,'fig6_grid_bounds.png'),'Resolution',200); savefig(f,fullfile(out,'fig6_grid_bounds.fig'));
end
