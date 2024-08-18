function [t_TDPC,mean_TDPC,var_TDPC,skew_TDPC,kurt_TDPC,th_sample_TDPC,X_sample_TDPC,...
    Y_sample_TDPC,dth_sample_TDPC,Nupdate] = TDPC(endtime,stopcrt_t,tspan,Ephi2,E,E1,E2,E3,E4,E_L2...
    ,basis_sample,xi_sample,theta,h,o,Lb,Ld,g,m,TD_sample)
%% property
opts_PC = odeset('RelTol',1e-5,'AbsTol',1e-8,'Events',@stopopt);
DOF = 3;
P = length(E);
t = 0;
iter = ceil(endtime/stopcrt_t);
tspan_TDPC = tspan;
Ephi2_new = Ephi2;
E_new = E;
E1_new = E1;
E2_new = E2;
E3_new = E3;
E_L2_new = E_L2;
TDbasis = basis_sample;
Nupdate = 0;    % number of update basis
%% initial set 
x_i = [ E4*cos(theta); E4*sin(theta); theta; zeros(P-1,1)];
xD_i = zeros(DOF*P,1);
xDD_i = zeros(1,DOF*P);
xDD_TDPC = xDD_i;
init_TDPC = [x_i;xD_i];

t_TDPC = [];
mean_TDPC = [];
var_TDPC = [];
mean_TDPC2 = [];
var_TDPC2 = [];
skew_TDPC = [];
kurt_TDPC = [];
th_sample_TDPC = [];
X_sample_TDPC = [];
Y_sample_TDPC = [];
dth_sample_TDPC = [];
% Cov = zeros(iter,9);
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
        mu_basis = sum(TDbasis)/TD_sample;
        mean_temp = zeros(2*DOF,length(te)-2);
        var_temp  = zeros(2*DOF,length(te)-2);
        skew_temp = zeros(2*DOF,length(te)-2);
        kurt_temp = zeros(2*DOF,length(te)-2);
        for j = 1:2*DOF
            mean_temp(j,:) =  mu_basis * ye(1:end-2,1+P*(j-1):j*P)';
            rv_temp = TDbasis * ye(1:end-2,1+P*(j-1):j*P)';
            mu2 = sum(rv_temp.^2)/TD_sample;
            mu3 = sum(rv_temp.^3)/TD_sample;
            mu4 = sum(rv_temp.^4)/TD_sample;
            var_temp(j,:) = mu2 - mean_temp(j,:).^2;
            skew_temp(j,:) = (mu3 - 3*mean_temp(j,:).*var_temp(j,:) - mean_temp(j,:).^3)./(sqrt(var_temp(j,:)).^3);
            kurt_temp(j,:) = (mu4-4*mean_temp(j,:).*mu3+6*mean_temp(j,:).^2.*mu2-3*mean_temp(j,:).^4)./var_temp(j,:).^2;
        end
        rvth  =  TDbasis * ye(1:end-2,1+2*P:3*P)';
        rvX  =  TDbasis * ye(1:end-2,1:P)';
        rvY  =  TDbasis * ye(1:end-2,1+P:2*P)';
        rvdth  =  TDbasis * ye(1:end-2,1+5*P:6*P)';

        t_TDPC = cat(1,t_TDPC,te(1:end-2));
        mean_TDPC = cat(2,mean_TDPC,mean_temp);
        var_TDPC = cat(2,var_TDPC,var_temp);
        th_sample_TDPC = cat(2,th_sample_TDPC,rvth);
%         X_sample_TDPC = cat(2,X_sample_TDPC,rvX);
%         Y_sample_TDPC = cat(2,Y_sample_TDPC,rvY);
        dth_sample_TDPC = cat(2,dth_sample_TDPC,rvdth);
        skew_TDPC = cat(2,skew_TDPC,skew_temp);
        kurt_TDPC = cat(2,kurt_TDPC,kurt_temp);
        % Calculate Variance 
        Xe = ye(end-1,:);
        id = 1;
        % random variable
        Nsample = length(TDbasis);
        th  =  Xe(2*P+1:3*P);
        dth  =  Xe(5*P+1:6*P);
        % Correlation
        rve = [ xi_sample,TDbasis*th',TDbasis*dth'];       
        rve = rve - sum(rve)/Nsample;        
        std_rv2 = sqrt(sum(rve.^2)/Nsample);
        co = zeros(size(rve,2),size(rve,2));
        for j = 1:size(rve,2)
            for k = 1:size(rve,2)
                co(j,k) = (sum( rve(:,j) .* rve(:,k) )/Nsample)/std_rv2(j)/std_rv2(k);
            end
        end
        trig = 0.80;
        if i ==1
            ctrv_old = rve(:,2:3);
        end
        % select random variable via correlation
        rv = [ xi_sample,TDbasis*th',TDbasis*dth'];
        rv_list = 1;
        for j = 2:size(rv,2)
            cr = 1;
            for k = 1:j-1
                if abs(co(k,j)) >= trig
                   cr = 0; 
                end
            end
            if cr == 1
                rv_list = [rv_list,j];
            end
        end
        if i ~= 1
            co_rv = zeros(1,size(rv,2)-1);
            for j = 2:size(rv,2)
                std_n = sqrt(sum(rve(:,j).^2)/Nsample);
                std_o = sqrt(sum(ctrv_old(:,j-1).^2)/Nsample);
                co_rv(j-1) = sum( rve(:,j) .* ctrv_old(:,j-1)  )/Nsample/ std_n/ std_o;
            end
            if abs(co_rv) >= 0.97
                    rv_list = 1;
            end
        end
        ctrv_old = rve(:,2:3);
        rv = rv(:,rv_list);
        if length(rv_list) == 1
            id = 0;
        end
%         fprintf("update random varable : %d \n",rv_list)
        if id ==1
            [TDbasis,Ephi2_new,E_new,E1_new, E2_new, E3_new, E_L2_new,init_TDPC]...
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
        mu_basis = sum(TDbasis)/TD_sample;
        mean_temp = zeros(2*DOF,length(te));
        var_temp  = zeros(2*DOF,length(te));
        skew_temp = zeros(2*DOF,length(te));
        kurt_temp = zeros(2*DOF,length(te));
        for j = 1:2*DOF
            mean_temp(j,:) =  mu_basis * ye(1:end,1+P*(j-1):j*P)';
            rv_temp = TDbasis * ye(1:end,1+P*(j-1):j*P)';
            mu2 = sum(rv_temp.^2)/TD_sample;
            mu3 = sum(rv_temp.^3)/TD_sample;
            mu4 = sum(rv_temp.^4)/TD_sample;
            var_temp(j,:) = mu2 - mean_temp(j,:).^2;
            skew_temp(j,:) = (mu3 - 3*mean_temp(j,:).*var_temp(j,:) - mean_temp(j,:).^3)./(sqrt(var_temp(j,:)).^3);
            kurt_temp(j,:) = (mu4-4*mean_temp(j,:).*mu3+6*mean_temp(j,:).^2.*mu2-3*mean_temp(j,:).^4)./var_temp(j,:).^2;
        end
        rvth  =  TDbasis * ye(1:end,1+2*P:3*P)';
        rvX  =  TDbasis * ye(1:end,1:P)';
        rvY  =  TDbasis * ye(1:end,1+P:2*P)';
        rvdth  =  TDbasis * ye(1:end,1+5*P:6*P)';

        t_TDPC = cat(1,t_TDPC,te(1:end));
        mean_TDPC = cat(2,mean_TDPC,mean_temp);
        var_TDPC = cat(2,var_TDPC,var_temp);
        th_sample_TDPC = cat(2,th_sample_TDPC,rvth);
%         X_sample_TDPC = cat(2,X_sample_TDPC,rvX);
%         Y_sample_TDPC = cat(2,Y_sample_TDPC,rvY);
        dth_sample_TDPC = cat(2,dth_sample_TDPC,rvdth);
        skew_TDPC = cat(2,skew_TDPC,skew_temp);
        kurt_TDPC = cat(2,kurt_TDPC,kurt_temp);
%         fprintf("Finish time : %d , Nupdate : %d \n",t,Nupdate)
    end
end
th_sample_TDPC = th_sample_TDPC(:,1:1/h:end);
dth_sample_TDPC = dth_sample_TDPC(:,1:1/h:end);
% X_sample_TDPC = X_sample_TDPC(:,1:1/h:end);
% Y_sample_TDPC = Y_sample_TDPC(:,1:1/h:end);

    function [Condition, isterminal, direction] = stopopt(te,y)
    Condition = te > t+stopcrt_t;
    isterminal = 1; 
    direction = 0;
    end
end
