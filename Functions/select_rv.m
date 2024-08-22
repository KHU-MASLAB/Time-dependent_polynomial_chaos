function [update_flag, rve] = select_rv(rve,co,ctrv_old,trig,i)
    rv_list = 1;
    [Nsample,Nvar] = size(rve);
    for j = 2:Nvar
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
        co_rv = zeros(1,Nvar-2);
        for j = 2:Nvar
            std_n = sqrt(sum(rve(:,j).^2)/Nsample);
            std_o = sqrt(sum(ctrv_old(:,j-1).^2)/Nsample);
            co_rv(j-1) = sum( rve(:,j) .* ctrv_old(:,j-1)  )/Nsample/ std_n/ std_o;
        end
        if abs(co_rv) >= 0.97
                rv_list = 1;
        end
    end
    rve = rve(:,rv_list);
    if length(rv_list) == 1
        update_flag = 0;
    else
        update_flag = 1;
    end