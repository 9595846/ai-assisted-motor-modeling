function verification = verify_motor_results()
%VERIFY_MOTOR_RESULTS 独立复核已保存的生产计划及热解析递推。
root=fileparts(mfilename('fullpath'));
data=load(fullfile(root,'results','results.mat'),'results'); r=data.results;
plans=[{r.baseline,r.nominal,r.fine},r.scenarios];
N=numel(plans); labels=strings(N,1); maxODEError=zeros(N,1);
originalDeadline=zeros(N,1); deliveryPass=false(N,1); durationPass=false(N,1);
for j=1:N
    s=plans{j}; c=s.cfg; labels(j)=string(s.label);
    temp=c.initialTemp; temps=zeros(numel(s.n)+1,1); temps(1)=temp;
    for k=1:numel(s.n)
        n=s.n(k); torque=c.loadFactor*(30+0.02*n);
        x=n/1000; y=torque/60;
        loss=[1 x y x*x x*y y*y]*r.fit.beta*c.lossFactor;
        if n==0, loss=0; end
        R=2*r.thermal.R;
        if s.fan(k), R=r.thermal.R*c.coolFactor; end
        rhs=@(~,T) (loss-(T-c.ambient)/R)/r.thermal.C;
        [~,vals]=ode45(rhs,[0 c.dt*60],temp,odeset('RelTol',1e-10,'AbsTol',1e-10));
        temp=vals(end); temps(k+1)=temp;
    end
    maxODEError(j)=max(abs(temps-s.temp));
    assert(maxODEError(j)<1e-6,'%s温度ODE复核失败',s.label);
    assert(max(temps)<=c.tempLimit+1e-5,'安全温度约束失败');
    q=[0;cumsum(s.n*0.01*c.dt)];
    originalDeadline(j)=q(480/c.dt+1);
    req=[1000;2300;3600;c.demand];
    deliveryPass(j)=all(q([120;240;360;c.deliveryEnd]/c.dt+1)>=req-1e-5);
    if c.deliveryEnd>480
        deliveryPass(j)=deliveryPass(j) && originalDeadline(j)>=5000-1e-5;
    end
    assert(deliveryPass(j),'交付数量复核失败');
    on=s.n>0; edges=[1;find(diff(on)~=0)+1;numel(on)+1];
    len=diff(edges)*c.dt; state=on(edges(1:end-1));
    durationPass(j)=all(len(state)>=10);
    for k=2:numel(len)-1
        if ~state(k), durationPass(j)=durationPass(j)&&len(k)>=10; end
    end
    assert(durationPass(j),'最短启停时间复核失败');
    assert(c.dt>=5,'风机时间网格不足以自动保证最短开关时长');
    assert(abs(q(end)-s.production(end))<1e-6,'产量记录不一致');
    assert(s.starts==sum(on&~[false;on(1:end-1)]),'启动次数不一致');
    if isfield(s,'solverOutput')
        assert(s.maxViolation<1e-4,'整数规划约束违反量过大');
    end
end
assert(all(r.robustValidation(:,6)==1),'27情景鲁棒验证失败');
verification=table(labels,maxODEError,originalDeadline,deliveryPass,durationPass);
writetable(verification,fullfile(root,'results','independent_verification.csv'));
disp(verification);
fprintf('Independent verification passed for %d schedules and 27 robust scenarios.\n',N);
end
