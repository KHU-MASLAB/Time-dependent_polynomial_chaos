function PC = SolveStochasticDynamics_TDPC(Prop,init,Galerkin,XiSample,BasisSample,opts)

%% initial config. (theta1 = 0)
P = length(Galerkin.E2);
DOF = length(init)/2;
init_PC = zeros(DOF*P*2,1);
init_PC(1:P) = 0.5*Galerkin.E1L*cos(Prop.th1_0);
init_PC(P+1:2*P) = 0.5*Galerkin.E1L*sin(Prop.th1_0);
init_PC(2*P+1:3*P) = Galerkin.E1*Prop.th1_0;
tspan_TDPC = Prop.tspan;
%% stochastic solver
stopcrt_t = opts.stopcrt_t;
t = 0; 
Nupdate = 0;
h_digit = ceil(log10(1/Prop.h));
iter = ceil(Prop.endTime/stopcrt_t);
ode_opt = odeset('RelTol',1e-6,'AbsTol',1e-8,'Events',@stopopt);
t_TDPC = []; mean_TDPC = []; var_TDPC = []; skew_TDPC = []; kurt_TDPC = [];
th1Sample_TDPC = []; th1DotSample_TDPC = [];
ctrv_old = [];
for i = 1:iter
    
    [te,ye] = ode45(@(t,x) SolveDAE_stochastic(t,x,Prop,Galerkin),tspan_TDPC,init_PC,ode_opt);
    
    P = length(Galerkin.E2);
    
    if  round(te(end),h_digit) < round(t+stopcrt_t,h_digit) % If the solver diverges, terminate the process.
        lent = length(te);
        [mean_temp,var_temp,skew_temp,kurt_temp,rvth,rvdth] = statistics(BasisSample,DOF,ye,P,lent);
        t_TDPC = cat(1,t_TDPC,te(1:end));
        mean_TDPC = cat(2,mean_TDPC,mean_temp);
        var_TDPC = cat(2,var_TDPC,var_temp);
        skew_TDPC = cat(2,skew_TDPC,skew_temp);
        kurt_TDPC = cat(2,kurt_TDPC,kurt_temp);
        th1Sample_TDPC = cat(2,th1Sample_TDPC,rvth);
        th1DotSample_TDPC = cat(2,th1DotSample_TDPC,rvdth);
        break;
    end
    
    if i ~= iter
        t = te(end-1,:);
        lent = length(te)-2;
        [mean_temp,var_temp,skew_temp,kurt_temp,rvth,rvdth] = statistics(BasisSample,DOF,ye,P,lent);
        t_TDPC = cat(1,t_TDPC,te(1:lent));
        mean_TDPC = cat(2,mean_TDPC,mean_temp);
        var_TDPC = cat(2,var_TDPC,var_temp);
        skew_TDPC = cat(2,skew_TDPC,skew_temp);
        kurt_TDPC = cat(2,kurt_TDPC,kurt_temp);
        th1Sample_TDPC = cat(2,th1Sample_TDPC,rvth);
        th1DotSample_TDPC = cat(2,th1DotSample_TDPC,rvdth);
        
        %% random variable
        Xe = ye(end-1,:);
        th  =  Xe(2*P+1:3*P);
        dth  =  Xe((DOF+2)*P+1:(DOF+3)*P);
        candidate = [th',dth'];
        rve = [ XiSample, BasisSample*candidate];
%         Xi_RvPlot(rve,1)        
 
        %% Correlation
        cor = correlation(rve);
        
        %% select random variable via correlation
        trig = 0.80;
        [update_flag, rv, rv_list] = select_rv(rve,cor,ctrv_old,trig,i);
        ctrv_old = rve;
        fprintf(['Update rv : ' repmat('%d\t',1,length(rv_list)) '\n'],rv_list)
        %% update basis and initial condtion
        if update_flag == 1
            [Galerkin,BasisSample,init_PC]...
                    = updategPC_GramSchmidt_sampling(opts,rv,Xe,XiSample,BasisSample);
            Nupdate = Nupdate + 1;   
        else
            init_PC = Xe';
        end
        tspan_TDPC = round(t,5) : Prop.h : Prop.endTime;
        P = length(Galerkin.E2);
        fprintf("Event time : %.2d , Nupdate : %d, size : %d \n",t,Nupdate,P)
    else
        t = te(end,:);
        lent = length(te);
        [mean_temp,var_temp,skew_temp,kurt_temp,rvth,rvdth] = statistics(BasisSample,DOF,ye,P,lent);
        t_TDPC = cat(1,t_TDPC,te(1:end));
        mean_TDPC = cat(2,mean_TDPC,mean_temp);
        var_TDPC = cat(2,var_TDPC,var_temp);
        skew_TDPC = cat(2,skew_TDPC,skew_temp);
        kurt_TDPC = cat(2,kurt_TDPC,kurt_temp);
        th1Sample_TDPC = cat(2,th1Sample_TDPC,rvth);
        th1DotSample_TDPC = cat(2,th1DotSample_TDPC,rvdth);
        fprintf("Finish time : %d , Nupdate : %d \n",t,Nupdate)
    end

end

PC.th1Sample_TDPC = th1Sample_TDPC';
PC.th1DotSample_TDPC = th1DotSample_TDPC';
PC.mean_TDPC = mean_TDPC';
PC.var_TDPC = var_TDPC';
PC.skew_TDPC = skew_TDPC';
PC.kurt_TDPC = kurt_TDPC';
PC.t_TDPC = t_TDPC;

    function [Condition, isterminal, direction] = stopopt(te,y)
        Condition = te > t+stopcrt_t;
        isterminal = 1; 
        direction = 0;
    end
end
