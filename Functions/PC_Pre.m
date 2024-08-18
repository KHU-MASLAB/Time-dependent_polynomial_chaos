function [E, Ephi2, E_L2, E1, E2, E3, E4] = PC_Pre(dist,basis,loc,wts,Lb,Ld)
Npoint = length(loc);
P = size(basis,2);
%% distribution
if dist == "Uniform"
    W = ones(length(loc),1)*0.5;
elseif dist == "Normal"
    W = 1./(sqrt(2*pi)*exp(-loc.^2/2));
    % W = exp(-0.5*loc.^2)/sqrt(2*pi);
end
%% calculate Galerkin Projetion
for i = 1:P
    basis(:,i) = basis(:,i)/sqrt(sum(basis(:,i).^2 .* W .* wts));
end
E = zeros(P,P);
for i = 1:P
    for j = i:P
        Pe = sum(basis(:,i) .* basis(:,j) .* W .* wts);
        E(i,j)= Pe;
        E(j,i)= Pe;
    end
end
Ephi2 = diag(E);

E_L2 = zeros(P,P);
for i=1:P
    for j=1:Npoint
        PCB_loc = basis(j,:);
        Pe = PCB_loc(i) .* PCB_loc * (Lb + Ld * loc(j))^2 * wts(j) * W(j);        
        E_L2(i,:)= E_L2(i,:) + Pe;
    end
end

E1 = zeros(P,P*P);
n=1;
for i=1:P
    for j=1:P
        for k=1:Npoint
            PCB_loc = basis(k,:);
            Pe = PCB_loc(i) * PCB_loc(j) * PCB_loc * wts(k) * W(k);
            E1(:,n)= E1(:,n) + Pe';
        end
        n=n+1;  
    end
end

E2 = zeros(P,P*P*P);
n = 1;
for i = 1:P    
    for j = 1:P
        for k = 1:P
            for l = 1:Npoint
                PCB_loc = basis(l,:);
                Pe = PCB_loc(i) * PCB_loc(j) * PCB_loc(k) * PCB_loc * wts(l) * W(l);
                E2(:,n)= E2(:,n) + Pe';    
            end
             n=n+1;
        end
    end
end

E3 = zeros(P,1);
for i = 1:Npoint
    PCB_loc = basis(i,:);
    Pe = PCB_loc * wts(i) * W(i);
    E3 = E3 + Pe';    
end

E4 = zeros(P,1);
for i=1:Npoint
    PCB_loc = basis(i,:);
    Pe = (PCB_loc./Ephi2') * (Lb + Ld * loc(i)) * wts(i) * W(i);
    E4 = E4 + Pe';
end



end