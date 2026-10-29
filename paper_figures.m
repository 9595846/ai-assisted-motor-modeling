function paper_figures()
%PAPER_FIGURES 论文所需补充统计和中文图形，不重新运行整数优化。
root=fileparts(mfilename('fullpath')); d=load(fullfile(root,'results','results.mat')); r=d.results;
out=fullfile(root,'paper','figures'); if ~isfolder(out), mkdir(out); end
data=readtable(fullfile(root,'results','loss_fit.csv'));
x=data.RPM/1000; y=data.Torque_Nm/60; loss=data.MeasuredLoss_W;
XX={[ones(size(x)) x y], [ones(size(x)) x y x.^2 y.^2], ...
    [ones(size(x)) x y x.^2 x.*y y.^2]};
stats=zeros(3,4);
for j=1:3
    X=XX{j}; b=X\loss; res=X*b-loss; loo=zeros(size(loss));
    for k=1:numel(loss)
        keep=true(size(loss)); keep(k)=false;
        loo(k)=X(k,:)*(X(keep,:)\loss(keep))-loss(k);
    end
    stats(j,:)=[size(X,2),sqrt(mean(res.^2)),sqrt(mean(loo.^2)),max(abs(res)./loss)*100];
end
tab=array2table(stats,'VariableNames',{'Parameters','RMSE_W','LOORMSE_W','MaxRelativeError_percent'});
tab=addvars(tab,["Linear";"Quadratic_without_interaction";"Full_quadratic"],'Before',1,'NewVariableNames','Model');
writetable(tab,fullfile(root,'paper','model_comparison.csv')); disp(tab);
plans={r.baseline,r.fine}; metrics=zeros(2,4); phase=zeros(8,6);
for j=1:2
    s=plans{j}; n=s.n; T=30+.02*n; xx=n/1000; yy=T/60;
    L=[ones(size(n)),xx,yy,xx.^2,xx.*yy,yy.^2]*r.fit.beta; L(n==0)=0;
    mech=2*pi*n.*T/60;
    metrics(j,:)=[sum(mech)*5/60000,sum(L)*5/60000,sum(s.fan)*120*5/60000,.03*s.starts];
    for k=1:4
        idx=(k-1)*24+(1:24); prices=[.65 1.1 .45 .9];
        e=sum(s.energySlot(idx));
        phase((j-1)*4+k,:)=[j,k,mean(n(idx)),sum(n(idx))*.01*5,e,e*prices(k)];
    end
end
writetable(array2table(metrics,'VariableNames',{'Mechanical_kWh','MotorLoss_kWh','Fan_kWh','Startup_kWh'}),...
    fullfile(root,'paper','energy_components.csv'));
writetable(array2table(phase,'VariableNames',{'Plan','Phase','AverageRPM','Production','Energy_kWh','Cost_yuan'}),...
    fullfile(root,'paper','phase_statistics.csv'));
font='SimSun'; color=lines(4);
f=figure('Color','w','Position',[50 50 1040 400]); tl=tiledlayout(f,1,2,'TileSpacing','compact');
ax=nexttile(tl); stem(ax,1:15,data.Residual_W,'filled','Color',color(1,:),'LineWidth',1.2);
yline(ax,0,'k-'); xlabel(ax,'观测序号'); ylabel(ax,'拟合残差 / W'); title(ax,'二次损耗模型残差'); grid(ax,'on');
ax=nexttile(tl); curve=readtable(fullfile(root,'results','speed_curve.csv'));
yyaxis(ax,'left'); plot(ax,curve.RPM,100*curve.Efficiency,'LineWidth',1.6); ylabel(ax,'电机效率 / %');
yyaxis(ax,'right'); plot(ax,curve.RPM,1000*curve.kWh_per_item,'LineWidth',1.6); ylabel(ax,'单件电耗 / Wh');
xlabel(ax,'转速 / (r/min)'); title(ax,'效率与单件电耗'); grid(ax,'on'); style(f,font);
exportgraphics(f,fullfile(out,'fig1_loss_efficiency.png'),'Resolution',220);

f=figure('Color','w','Position',[50 50 1040 400]); tl=tiledlayout(f,1,2,'TileSpacing','compact');
ax=nexttile(tl); t=(0:1:1800)'; R=r.thermal.R; C=r.thermal.C;
plot(ax,t/60,25+1000*R*(1-exp(-t/(R*C))),'LineWidth',1.6); hold(ax,'on');
plot(ax,t/60,25+55*exp(-t/(R*C)),'LineWidth',1.6);
th=readtable(fullfile(root,'results','thermal_fit.csv'));
scatter(ax,th.Time_min,th.HeatingMeasured_C,30,color(1,:),'o');
scatter(ax,th.Time_min,th.CoolingMeasured_C,30,color(2,:),'s');
xlabel(ax,'时间 / min'); ylabel(ax,'绕组温度 / ℃'); grid(ax,'on');
legend(ax,{'加热模型','冷却模型','加热观测','冷却观测'},'Location','best');
ax=nexttile(tl); cont=readtable(fullfile(root,'results','continuous_operation.csv'));
plot(ax,cont.RPM,cont.EquilibriumOn_C,'LineWidth',1.6); hold(ax,'on');
plot(ax,cont.RPM,cont.EquilibriumOff_C,'LineWidth',1.6); yline(ax,90,'r--');
xlabel(ax,'转速 / (r/min)'); ylabel(ax,'平衡温度 / ℃'); grid(ax,'on');
legend(ax,{'风机开启','风机关闭','温度上限'},'Location','best'); style(f,font);
exportgraphics(f,fullfile(out,'fig2_thermal.png'),'Resolution',220);

f=figure('Color','w','Position',[50 50 1060 700]); tl=tiledlayout(f,3,1,'TileSpacing','compact');
for j=1:3
    ax=nexttile(tl); hold(ax,'on');
    for k=1:2
        s=plans{k}; time=(0:96)'/12+8;
        switch j
            case 1, stairs(ax,time,[s.n;s.n(end)],'LineWidth',1.5);
            case 2, plot(ax,time,s.production,'LineWidth',1.5);
            case 3, plot(ax,time,s.temp,'LineWidth',1.5);
        end
    end
    if j==2, scatter(ax,[10 12 14 16],[1000 2300 3600 5000],25,'k','filled'); end
    if j==3, yline(ax,90,'r--'); end
    yl={'转速 / (r/min)','累计产量 / 件','绕组温度 / ℃'}; ylabel(ax,yl{j});
    grid(ax,'on'); xlim(ax,[8 16]);
    if j==1, legend(ax,{'基准方案','50转网格方案'},'Location','best'); end
end
xlabel(ax,'时刻 / h'); style(f,font);
exportgraphics(f,fullfile(out,'fig3_nominal_schedule.png'),'Resolution',220);

f=figure('Color','w','Position',[50 50 1060 450]); ax=axes(f); hold(ax,'on');
ss={r.nominal,r.replay,r.scenarios{4}};
for j=1:3
    s=ss{j}; plot(ax,(0:96)'/12+8,s.temp,'LineWidth',1.5);
end
yline(ax,90,'r--'); xlabel(ax,'时刻 / h'); ylabel(ax,'绕组温度 / ℃');
grid(ax,'on'); xlim(ax,[8 16]); legend(ax,{'名义条件原计划','最坏条件直接回放','最坏条件鲁棒计划'},'Location','best');
style(f,font); exportgraphics(f,fullfile(out,'fig4_robustness.png'),'Resolution',220);
save(fullfile(root,'paper','supplementary_statistics.mat'),'stats','metrics','phase');
fprintf('Paper figures and supplementary statistics completed.\n');
end

function style(f,font)
ax=findall(f,'Type','axes'); set(ax,'FontName',font,'FontSize',14);
lg=findall(f,'Type','legend'); set(lg,'FontName',font,'FontSize',13);
end
