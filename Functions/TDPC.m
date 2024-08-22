function [t_TDPC,mean_TDPC,var_TDPC,skew_TDPC,kurt_TDPC,th_sample_TDPC,...
    dth_sample_TDPC,Nupdate] = TDPC(endtime,stopcrt_t,tspan,E,E1,E2,E3,E4,E_L2...
    ,basis_sample,xi_sample,theta,h,o,Lb,Ld,g,m,TD_sample)
%% property
opts_PC = odeset('RelTol',1e-5,'AbsTol',1e-8,'Events',@stopopt);
DOF = 3;
P = length(E);
t = 0;
iter = ceil(endtime/stopcrt_t);
tspan_TDPC = tspan;
E_new = E; E1_new = E1; E2_new = E2; E3_new = E3; E_L2_new = E_L2;
TDbasis = basis_sample;
Nupdate = 0;

%% initial set 
x_i = [ E4*cos(theta); E4*sin(theta); theta; zeros(P-1,1)];
xD_i = zeros(DOF*P,1);
init_TDPC = [x_i;xD_i];
t_TDPC = []; mean_TDPC = []; var_TDPC = []; skew_TDPC = []; kurt_TDPC = [];
th_sample_TDPC = []; dth_sample_TDPC = [];
ctrv_old = zeros(TD_sample,2);

for i = 1:iter
    
    [te,ye,~,~,ie] = ode45(@(t,x) SP_ode_TDPC(t,x, g, m, E_new, E1_new, E2_new, E3_new, E_L2_new, DOF)...
        ,tspan_TDPC,init_TDPC,opts_PC);
    if length(te) == 1 || abs(te(1)-te(end)) < 1E-2
        break;
    end
    
    if i ~= iter
        P = length(E_new);
        t = te(end-1,:);
        lent = length(te)-2;
        [mean_temp,var_temp,skew_temp,kurt_temp,rvth,rvdth] = statistics(TDbasis,DOF,ye,P,lent);
        
        t_TDPC = cat(1,t_TDPC,te(1:lent));
        mean_TDPC = cat(2,mean_TDPC,mean_temp);
        var_TDPC = cat(2,var_TDPC,var_temp);
        th_sample_TDPC = cat(2,th_sample_TDPC,rvth);
        dth_sample_TDPC = cat(2,dth_sample_TDPC,rvdth);
        skew_TDPC = cat(2,skew_TDPC,skew_temp);
        kurt_TDPC = cat(2,kurt_TDPC,kurt_temp);
        
        % random variable
        Xe = ye(end-1,:);
        th  =  Xe(2*P+1:3*P);
        dth  =  Xe(5*P+1:6*P);
        
        % Correlation
        rve = [ xi_sample,TDbasis*th',TDbasis*dth'];       
        [rve,co] = correlation(rve);
        
        % select random variable via correlation
        trig = 0.80;
        [update_flag, rv] = select_rv(rve,co,ctrv_old,trig,i);
        ctrv_old = rve;
        if update_flag ==1
            [TDbasis,E_new,E1_new, E2_new, E3_new, E_L2_new,init_TDPC]...
                = updategPC_GramSchmidt(t,rv,Xe,xi_sample,o,Lb,Ld,TDbasis);
            Nupdate = Nupdate + 1;
        else
            init_TDPC = Xe';
        end
        tspan_TDPC = round(t,5) : h : endtime;
        fprintf("Event time : %.2d , Nupdate : %d, size : %d \n",t,Nupdate,length(E_new))
    else
        P = length(E_new);
        t = te(end,:);
        lent = length(te);
        [mean_temp,var_temp,skew_temp,kurt_temp,rvth,rvdth] = statistics(TDbasis,DOF,ye,P,lent);
        t_TDPC = cat(1,t_TDPC,te(1:end));
        mean_TDPC = cat(2,mean_TDPC,mean_temp);
        var_TDPC = cat(2,var_TDPC,var_temp);
        skew_TDPC = cat(2,skew_TDPC,skew_temp);
        kurt_TDPC = cat(2,kurt_TDPC,kurt_temp);
        th_sample_TDPC = cat(2,th_sample_TDPC,rvth);
        dth_sample_TDPC = cat(2,dth_sample_TDPC,rvdth);
        fprintf("Finish time : %d , Nupdate : %d \n",t,Nupdate)
    end
end
th_sample_TDPC = th_sample_TDPC(:,1:1/h:end);
dth_sample_TDPC = dth_sample_TDPC(:,1:1/h:end);


    function [Condition, isterminal, direction] = stopopt(te,y)
    Condition = te > t+stopcrt_t;
    isterminal = 1; 
    direction = 0;
    end
end
