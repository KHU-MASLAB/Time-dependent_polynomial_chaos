function midx = multi_idx(d,n)

P = partitions(d,1:d,n);
N = size(P,1);
num = 1;
for i = 1:N
    a = sum(P(num,:));
    if a > n
        P(num,:) = [];
        num = num-1;
    end
    num = num+1;
end
N = size(P,1);
C = cell(1,N);
for k = 1:N
    tmp = repelem(1:d,P(k,:));
    tmp(end+1:n) = 0;
    C{k} = flip(unique(perms(tmp),'rows'));
end
midx = vertcat(C{:});

end