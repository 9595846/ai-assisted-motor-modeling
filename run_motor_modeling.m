function results = run_motor_modeling(action,varargin)
%RUN_MOTOR_MODELING 电机数学建模完整求解；MATLAB R2024a + Optimization Toolbox。
% 在本目录执行 results = run_motor_modeling。全部数据内嵌，无外部数据依赖。
% 研究扩展接口复用同一模型，避免不同实验使用不同物理定义。
if nargin>0
    switch action
        case 'solve', results=optimizeSchedule(varargin{:});
        case 'simulate', results=simulate(varargin{:});
        case 'export', writeSchedule(varargin{:}); results=[];
        otherwise, error('Unknown action');
    end
    return
end
root = fileparts(mfilename('fullpath'));
out = fullfile(root,'results');
if ~isfolder(out), mkdir(out); end
diary(fullfile(out,'matlab_run.log'));
cleanup = onCleanup(@() diary('off'));
fprintf('MATLAB %s\n',version);
assert(license('test','Optimization_Toolbox'),'需要 Optimization Toolbox');

%% 问题一：有量纲归一化的二次损耗拟合与留一验证
nData = [500;750;1000;1250;1500];
tData = [30 60 90];
lossData = [300 470 760;350 530 840;420 620 950;510 730 1090;620 860 1250];
[nn,tt] = ndgrid(nData,tData);
X = design(nn(:),tt(:)); y = lossData(:);
beta = X\y; residual = X*beta-y;
loo = zeros(size(y));
for k=1:numel(y)
    keep = true(size(y)); keep(k)=false;
    loo(k)=X(k,:)*(X(keep,:)\y(keep))-y(k);
end
fit.beta = beta;
fit.rmse = sqrt(mean(residual.^2));
fit.r2 = 1-sum(residual.^2)/sum((y-mean(y)).^2);
fit.maxRelativeError = max(abs(residual)./y);
fit.looRMSE = sqrt(mean(loo.^2));
lossfun = @(n,T) reshape(design(n(:),T(:))*beta,size(n));
ng = (500:1:1500)'; tg=30+0.02*ng;
lg=lossfun(ng,tg); pg=2*pi*ng.*tg/60;
eg=(pg+lg)./(600*ng); % kWh/件；q=0.01*n 件/min
[emin,imin]=min(eg); fit.bestSpeed=ng(imin); fit.minEnergy=emin;
curve=table(ng,tg,lg,pg,(pg./(pg+lg)),eg,...
    'VariableNames',{'RPM','Torque_Nm','Loss_W','Output_W','Efficiency','kWh_per_item'});
writetable(curve,fullfile(out,'speed_curve.csv'));
writetable(table(nn(:),tt(:),y,X*beta,residual,loo,...
    'VariableNames',{'RPM','Torque_Nm','MeasuredLoss_W','FittedLoss_W','Residual_W','LOOResidual_W'}),...
    fullfile(out,'loss_fit.csv'));
fprintf('Loss beta = '); fprintf('%.9g ',beta); fprintf('\n');
disp(fit);

%% 问题二：同时利用加热、冷却数据识别 R、C
tsec=[0;5;10;15;20;30]*60;
heat=[25;45.5;58.9;67.8;73.7;80.1];
cool=[80;61.3;48.9;40.8;35.4;29.5];
resfun=@(v) [25+1000*v(1)*(1-exp(-tsec/(v(1)*v(2))))-heat;...
    25+55*exp(-tsec/(v(1)*v(2)))-cool];
opts=optimoptions('lsqnonlin','Display','off');
[v,~,rr,ef,oo]=lsqnonlin(resfun,[0.06 12000],[0.001 100],[1 1e6],opts);
assert(ef>0,'热参数识别没有收敛：%s',oo.message);
thermal.R=v(1); thermal.C=v(2); thermal.tau=v(1)*v(2);
thermal.rmse=sqrt(mean(rr.^2));
thermal.Roff=2*v(1);
fprintf('Thermal R=%.9f K/W C=%.4f J/K tau=%.5f min RMSE=%.6f C\n',...
    thermal.R,thermal.C,thermal.tau/60,thermal.rmse);
eqOn=25+thermal.R*lg; eqOff=25+thermal.Roff*lg;
temp8=eqOn+(30-eqOn).*exp(-480*60/thermal.tau);
hit=inf(size(ng)); h=eqOff>90;
hit(h)=-thermal.Roff*thermal.C/60*log((eqOff(h)-90)./(eqOff(h)-30));
continuous=table(ng,eqOn,eqOff,temp8,hit,eqOn<=90,eqOff<=90,...
    'VariableNames',{'RPM','EquilibriumOn_C','EquilibriumOff_C','Temperature8hOn_C',...
    'TimeTo90Off_min','IndefiniteOn','IndefiniteOff'});
writetable(continuous,fullfile(out,'continuous_operation.csv'));
writetable(table(tsec/60,heat,heat+rr(1:6),cool,cool+rr(7:12),...
    'VariableNames',{'Time_min','HeatingMeasured_C','HeatingFit_C','CoolingMeasured_C','CoolingFit_C'}),...
    fullfile(out,'thermal_fit.csv'));

%% 问题三、四：精确离散热递推 + 混合整数线性规划
cfg.dt=5; cfg.horizon=480; cfg.speeds=(500:100:1500)';
cfg.speeds=unique([cfg.speeds;1250]);
cfg.ambient=25; cfg.loadFactor=1; cfg.demand=5000;
cfg.fanPower=120; cfg.coolFactor=1; cfg.lossFactor=1;
cfg.initialTemp=30; cfg.tempLimit=90; cfg.maxTime=90;
cfg.deliveryEnd=480;
base=baseline(cfg,thermal,beta);
writeSchedule(base,out,'baseline');
nominal=optimizeSchedule(cfg,thermal,beta,'nominal');
writeSchedule(nominal,out,'nominal');
save(fullfile(out,'checkpoint.mat'),'fit','thermal','cfg','base','nominal');

% 更细的50 r/min速度网格：仅作为离散精度对比，不是连续最优证书。
finecfg=cfg; finecfg.speeds=(500:50:1500)';
fine=optimizeSchedule(finecfg,thermal,beta,'fine_speed_grid');
writeSchedule(fine,out,'fine_speed_grid');

% 原方案在不确定条件下的直接回放
stressCfg=cfg; stressCfg.ambient=35; stressCfg.loadFactor=1.10;
stressCfg.demand=5500; stressCfg.lossFactor=1.03;
replay=simulate(nominal.n,nominal.fan,stressCfg,thermal,beta);
replay.label='nominal_replayed_worst';
writeSchedule(replay,out,'nominal_replayed_worst');

% 单因素、订单扰动与联合最坏情景重优化
names={'hot_35C','load_plus10','order_5500','robust_worst','robust_extended','robust_upgraded','robust_both'};
scenarios=cell(size(names));
for k=1:numel(names)
    c=cfg;
    switch names{k}
        case 'hot_35C', c.ambient=35;
        case 'load_plus10', c.loadFactor=1.10;
        case 'order_5500', c.demand=5500;
        otherwise
            c=stressCfg;
            if ismember(names{k},{'robust_extended','robust_both'})
                c.horizon=510; c.deliveryEnd=510;
            end
            if ismember(names{k},{'robust_upgraded','robust_both'})
                c.fanPower=200; c.coolFactor=0.8;
            end
    end
    scenarios{k}=optimizeSchedule(c,thermal,beta,names{k});
    writeSchedule(scenarios{k},out,names{k});
    save(fullfile(out,'scenario_checkpoint.mat'),'scenarios','names');
end

% 鲁棒方案在27个环境/负载/损耗场景下验证温度和交付
robust=scenarios{4}; checks=zeros(27,6); row=0;
for ta=[25 30 35]
    for lf=[0.9 1 1.1]
        for pf=[0.97 1 1.03]
            c=stressCfg; c.ambient=ta; c.loadFactor=lf; c.lossFactor=pf;
            s=simulate(robust.n,robust.fan,c,thermal,beta); row=row+1;
            checks(row,:)=[ta lf pf s.maxTemp s.cost s.feasible];
        end
    end
end
writetable(array2table(checks,'VariableNames',{'Ambient_C','LoadFactor','LossFactor','MaxTemp_C','Cost_yuan','Feasible'}),...
    fullfile(out,'robust_validation_27.csv'));
assert(all(checks(:,6)==1),'鲁棒回放存在不可行情景');

allSchedules=[{base,nominal,fine,replay},scenarios];
labels=[{'baseline','nominal','fine_speed_grid','nominal_replayed_worst'},names];
metrics=zeros(numel(allSchedules),9);
for k=1:numel(allSchedules)
    s=allSchedules{k};
    metrics(k,:)=[s.energy s.cost s.maxTemp s.production(end) s.starts s.feasible s.exitflag s.relativeGap s.maxViolation];
end
summary=array2table(metrics,'VariableNames',{'Energy_kWh','Cost_yuan','MaxTemp_C','Production','Starts','Feasible',...
    'Exitflag','RelativeGap','MaxConstraintViolation'});
summary=addvars(summary,string(labels(:)),'Before',1,'NewVariableNames','Scenario');
writetable(summary,fullfile(out,'summary.csv'));
disp(summary);
results=struct('fit',fit,'thermal',thermal,'cfg',cfg,'baseline',base,'nominal',nominal,...
    'fine',fine,'replay',replay,'scenarios',{scenarios},'summary',summary,'robustValidation',checks);
save(fullfile(out,'results.mat'),'results');
makeFigures(results,ng,lg,pg,eg,tsec,heat,cool,rr,out);
fid=fopen(fullfile(out,'numerical_summary.json'),'w','n','UTF-8');
fprintf(fid,'%s',jsonencode(struct('fit',fit,'thermal',thermal,'summary',table2struct(summary)),'PrettyPrint',true));
fclose(fid);
verify_motor_results();
write_motor_report(results,root);
fprintf('Finished. Results: %s\n',out);
end

function X=design(n,T)
% x=n/1000, y=T/60；完整二次响应面
x=n/1000; y=T/60;
X=[ones(size(x)),x,y,x.^2,x.*y,y.^2];
end

function s=optimizeSchedule(c,th,beta,label)
fprintf('\nSolving %s ...\n',label);
N=c.horizon/c.dt; ns=numel(c.speeds); n=c.speeds;
load=c.loadFactor*(30+0.02*n);
L=c.lossFactor*(design(n,load)*beta);
P=2*pi*n.*load/60+L;
assert(all(L>=0)&all(2*pi*n.*load/60<=15000),'工况超出额定功率或损耗异常');
price=tariff((0:N-1)'*c.dt);
y=optimvar('mode',N,ns,'Type','integer','LowerBound',0,'UpperBound',1);
u=optimvar('on',N,1,'Type','integer','LowerBound',0,'UpperBound',1);
fan=optimvar('fan',N,1,'Type','integer','LowerBound',0,'UpperBound',1);
start=optimvar('start',N,1,'Type','integer','LowerBound',0,'UpperBound',1);
stop=optimvar('stop',N,1,'Type','integer','LowerBound',0,'UpperBound',1);
T=optimvar('temp',N+1,1,'LowerBound',min(c.ambient,c.initialTemp),'UpperBound',c.tempLimit);
prob=optimproblem('Objective',sum(price.*((y*P+c.fanPower*fan)*c.dt/60000+0.03*start)));
prob.Constraints.mode=sum(y,2)==u;
previous=[0;u(1:end-1)];
prob.Constraints.transition=u-previous==start-stop;
prob.Constraints.noSimultaneous=start+stop<=1;
% 一般时间网格：持续时间按分钟定义，不能固定为两个区间。
runSlots=ceil(10/c.dt); fanSlots=ceil(5/c.dt);
minRun=optimconstr(N,1); minStop=optimconstr(N,1);
for k=1:N
    last=min(N,k+runSlots-1); len=last-k+1;
    if len==runSlots
        minRun(k)=sum(u(k:last))>=runSlots*start(k);
    else
        minRun(k)=-start(k)>=0;
    end
    minStop(k)=sum(1-u(k:last))>=len*stop(k);
end
prob.Constraints.minRun=minRun; prob.Constraints.minStop=minStop;
if fanSlots>1
    % 初始风机可自由选状态，但班次内新状态须保持至少5min。
    fanHoldUp=optimconstr(N-1,1); fanHoldDown=optimconstr(N-1,1);
    for k=2:N
        last=min(N,k+fanSlots-1); len=last-k+1;
        if len<fanSlots
            fanHoldUp(k-1)=fan(k)-fan(k-1)>=0;
            fanHoldDown(k-1)=fan(k-1)-fan(k)>=0;
        else
            fanHoldUp(k-1)=sum(fan(k:last))>=fanSlots*(fan(k)-fan(k-1));
            fanHoldDown(k-1)=sum(1-fan(k:last))>=fanSlots*(fan(k-1)-fan(k));
        end
    end
    prob.Constraints.fanHoldUp=fanHoldUp; prob.Constraints.fanHoldDown=fanHoldDown;
    prob.Constraints.fanInitial=fan(1:fanSlots)==fan(1);
end
if isfield(c,'fixedSpeed')
    prob.Constraints.fixedSpeed=y*n==c.fixedSpeed;
end
if isfield(c,'alwaysFan') && c.alwaysFan
    prob.Constraints.alwaysFan=fan==1;
end
prod=cumsum(y*(0.01*n*c.dt));
idx=[120 240 360 c.deliveryEnd]/c.dt;
prob.Constraints.delivery=prod(idx)>=[1000;2300;3600;c.demand];
if c.deliveryEnd>480
    prob.Constraints.originalDeadline=prod(480/c.dt)>=5000;
end
prob.Constraints.initial=T(1)==c.initialTemp;
Ron=th.R*c.coolFactor; Roff=2*th.R;
aon=exp(-c.dt*60/(Ron*th.C)); aoff=exp(-c.dt*60/(Roff*th.C));
heat=y*L;
predOn=aon*T(1:end-1)+(1-aon)*(c.ambient+Ron*heat);
predOff=aoff*T(1:end-1)+(1-aoff)*(c.ambient+Roff*heat);
% M=200℃覆盖温度范围及全部可选工况递推残差；风机状态选取准确递推式。
M=200;
prob.Constraints.onHi=T(2:end)-predOn<=M*(1-fan);
prob.Constraints.onLo=T(2:end)-predOn>=-M*(1-fan);
prob.Constraints.offHi=T(2:end)-predOff<=M*fan;
prob.Constraints.offLo=T(2:end)-predOff>=-M*fan;
opts=optimoptions('intlinprog','Display','final','MaxTime',c.maxTime);
if isfield(c,'gapTolerance'), opts.RelativeGapTolerance=c.gapTolerance; end
timer=tic;
if isfield(c,'seed')
    seed=c.seed; x0.mode=double(seed.n==n'); x0.on=double(seed.n>0);
    x0.fan=seed.fan; x0.start=max(0,[x0.on(1);diff(x0.on)]);
    x0.stop=max(0,-[x0.on(1);diff(x0.on)]); x0.temp=seed.temp;
    assert(isfinite(evaluate(prob.Objective,x0)));
    for k=1:numel(fieldnames(prob.Constraints))
        fn=fieldnames(prob.Constraints);
        assert(max(infeasibility(prob.Constraints.(fn{k}),x0),[],'all')<1e-4,'初始计划不满足约束');
    end
    [sol,fval,ef,output]=solve(prob,x0,'Options',opts);
else
    [sol,fval,ef,output]=solve(prob,'Options',opts);
end
solveSeconds=toc(timer);
if isempty(fieldnames(sol)) || isempty(sol.mode)
    error('场景%s未找到可行解：%s',label,output.message);
end
% 允许时间上限返回已验证可行解，但必须保留gap而不能声称精确最优。
fields=fieldnames(prob.Constraints); maxViol=0;
for k=1:numel(fields)
    vv=infeasibility(prob.Constraints.(fields{k}),sol);
    maxViol=max(maxViol,max(vv,[],'all'));
end
assert(maxViol<1e-4,'求解约束违反量过大：%g',maxViol);
yy=round(sol.mode); chosen=yy*n;
s=simulate(chosen,round(sol.fan),c,th,beta);
s.label=label; s.exitflag=double(ef); s.maxViolation=maxViol;
s.relativeGap=NaN;
if isfield(output,'relativegap'), s.relativeGap=output.relativegap; end
s.solverOutput=output; s.solverObjective=fval;
s.solveSeconds=solveSeconds;
assert(s.feasible,'连续热模型回放未通过');
assert(abs(s.cost-fval)<1e-3,'目标函数和回放电费不一致');
assert(max(abs(s.temp-sol.temp))<1e-3,'优化热递推和独立回放不一致');
fprintf('%s: cost %.6f energy %.6f Tmax %.6f starts %d flag %g gap %g\n',...
    label,s.cost,s.energy,s.maxTemp,s.starts,s.exitflag,s.relativeGap);
end

function s=simulate(n,fan,c,th,beta)
N=numel(n); n=n(:); fan=fan(:); on=n>0;
load=c.loadFactor*(30+0.02*n);
L=c.lossFactor*(design(n,load)*beta); L(~on)=0;
P=2*pi*n.*load/60+L+c.fanPower*fan;
starts=on&~[false;on(1:end-1)];
energySlot=P*c.dt/60000+0.03*starts;
price=tariff((0:N-1)'*c.dt);
T=zeros(N+1,1); T(1)=c.initialTemp;
for k=1:N
    R=2*th.R;
    if fan(k), R=th.R*c.coolFactor; end
    a=exp(-c.dt*60/(R*th.C));
    T(k+1)=c.ambient+R*L(k)+(T(k)-c.ambient-R*L(k))*a;
end
production=[0;cumsum(0.01*n*c.dt)];
req=[1000;2300;3600;c.demand];
delivery=production([120;240;360;c.deliveryEnd]/c.dt+1);
originalDeadlineOK=c.deliveryEnd<=480 || production(480/c.dt+1)>=5000-1e-5;
% 一阶系统在常值输入区间内单调，端点界即整个区间的精确温度界。
tempok=max(T)<=c.tempLimit+1e-5;
edges=[1;find(diff(on)~=0)+1;N+1]; lengths=diff(edges)*c.dt;
states=on(edges(1:end-1));
durationok=all(lengths(states)>=10-1e-8);
for j=1:numel(lengths)
    if ~states(j) && j~=1 && j~=numel(lengths)
        durationok=durationok && lengths(j)>=10-1e-8;
    end
end
fanEdges=[1;find(diff(fan)~=0)+1;N+1];
fanDurationOK=all(diff(fanEdges)*c.dt>=5-1e-8);
s=struct('n',n,'fan',fan,'temp',T,'production',production,'power',P,...
    'energySlot',energySlot,'energy',sum(energySlot),'cost',sum(price.*energySlot),...
    'maxTemp',max(T),'starts',sum(starts),'delivery',delivery,...
    'feasible',tempok && durationok && fanDurationOK && originalDeadlineOK && all(delivery>=req-1e-5),...
    'exitflag',NaN,'relativeGap',NaN,'maxViolation',NaN,'cfg',c);
end

function s=baseline(c,th,beta)
% 基准在400min恰好完成5000件，无需截断到时间网格之外。
N=c.horizon/c.dt; n=zeros(N,1); n(1:400/c.dt)=1250;
s=simulate(n,ones(N,1),c,th,beta); s.label='baseline';
end

function p=tariff(t)
p=0.90*ones(size(t));
p(t<120)=0.65; p(t>=120 & t<240)=1.10; p(t>=240 & t<360)=0.45;
end

function writeSchedule(s,out,label)
N=numel(s.n); dt=s.cfg.dt;
clock=@(m) string(datetime(2026,1,1,8,0,0)+minutes(m),'HH:mm');
edges=[1;find(diff(s.n)~=0 | diff(s.fan)~=0)+1;N+1];
a=edges(1:end-1); b=edges(2:end)-1;
startTime=clock((a-1)*dt); endTime=clock(b*dt);
tab=table(startTime,endTime,s.n(a),s.fan(a),s.temp(a),s.temp(b+1),s.production(b+1),...
    'VariableNames',{'Start','End','RPM','FanOn','TempStart_C','TempEnd_C','CumulativeItems'});
writetable(tab,fullfile(out,[label '_schedule.csv']));
slot=table((0:N-1)'*dt,s.n,s.fan,s.power,s.temp(1:end-1),s.temp(2:end),s.production(2:end),...
    'VariableNames',{'Elapsed_min','RPM','FanOn','Power_W','TempStart_C','TempEnd_C','CumulativeItems'});
writetable(slot,fullfile(out,[label '_timeseries.csv']));
end

function makeFigures(r,n,L,P,e,t,heat,cool,res,out)
% 多幅静态研究图；均明确坐标、单位并保存可编辑FIG。
f=figure('Color','w','Position',[80 80 1150 760],'Name','Motor parameter identification');
tl=tiledlayout(f,2,2,'TileSpacing','compact');
ax=nexttile(tl); plot(ax,n,L,'LineWidth',1.8); grid(ax,'on');
xlabel(ax,'Speed (r/min)'); ylabel(ax,'Loss (W)'); title(ax,'Load-following motor loss');
ax=nexttile(tl); plot(ax,n,P./(P+L)*100,'LineWidth',1.8); grid(ax,'on');
xlabel(ax,'Speed (r/min)'); ylabel(ax,'Efficiency (%)');
ax=nexttile(tl); plot(ax,n,e*1000,'LineWidth',1.8); grid(ax,'on');
xlabel(ax,'Speed (r/min)'); ylabel(ax,'Energy (Wh/item)');
ax=nexttile(tl); plot(ax,t/60,[heat cool],'o','LineWidth',1.3); hold(ax,'on');
plot(ax,t/60,[heat+res(1:6),cool+res(7:12)],'-','LineWidth',1.8);
grid(ax,'on'); xlabel(ax,'Time (min)'); ylabel(ax,'Temperature (C)');
legend(ax,{'Heating data','Cooling data','Heating model','Cooling model'},'Location','best');
exportgraphics(f,fullfile(out,'01_identification.png'),'Resolution',180);
savefig(f,fullfile(out,'01_identification.fig'));

f=figure('Color','w','Position',[100 60 1250 900],'Name','Production scheduling');
tl=tiledlayout(f,4,1,'TileSpacing','compact');
ss={r.baseline,r.nominal}; colors=lines(2);
for j=1:4
    ax=nexttile(tl); hold(ax,'on');
    for k=1:2
        s=ss{k}; time=(0:numel(s.n))'*s.cfg.dt/60+8;
        switch j
            case 1, stairs(ax,time,[s.n;s.n(end)],'Color',colors(k,:),'LineWidth',1.4);
            case 2, plot(ax,time,s.production,'Color',colors(k,:),'LineWidth',1.5);
            case 3, plot(ax,time,s.temp,'Color',colors(k,:),'LineWidth',1.5);
            case 4, stairs(ax,time,[s.power;s.power(end)]/1000,'Color',colors(k,:),'LineWidth',1.4);
        end
    end
    if j==2, scatter(ax,[10 12 14 16],[1000 2300 3600 5000],35,'k','filled'); end
    if j==3, yline(ax,90,'r--'); end
    labels={'Speed (r/min)','Cumulative items','Temperature (C)','Power (kW)'};
    ylabel(ax,labels{j}); grid(ax,'on'); xlim(ax,[8 16]);
    if j==1, legend(ax,{'Baseline','Optimized'},'Location','best'); end
end
xlabel(ax,'Clock time (h)');
exportgraphics(f,fullfile(out,'02_schedule.png'),'Resolution',180);
savefig(f,fullfile(out,'02_schedule.fig'));

f=figure('Color','w','Position',[120 100 1150 660],'Name','Robustness and investment');
tl=tiledlayout(f,2,1,'TileSpacing','compact');
ax=nexttile(tl); hold(ax,'on');
sets={r.nominal,r.replay,r.scenarios{4}};
for k=1:3
    s=sets{k}; plot(ax,(0:numel(s.n))*s.cfg.dt/60+8,s.temp,'LineWidth',1.5);
end
yline(ax,90,'r--'); grid(ax,'on'); ylabel(ax,'Temperature (C)'); xlabel(ax,'Clock time (h)');
legend(ax,{'Nominal optimized','Nominal replayed worst case','Robust optimized'},'Location','best');
ax=nexttile(tl); ids=8:11;
bar(ax,categorical({'Robust 8h','Extend 30min','Upgrade fan','Both'}),r.summary.Cost_yuan(ids));
ylabel(ax,'Electricity cost (yuan)'); grid(ax,'on');
exportgraphics(f,fullfile(out,'03_robustness.png'),'Resolution',180);
savefig(f,fullfile(out,'03_robustness.fig'));
end
