function [ W ] = norm_adj( adj )

degree = sum(adj');
W = zeros(length(adj(:,1)),length(adj(1,:)));
for i = 1:length(adj(:,1)),
    for j = 1:length(adj(1,:)),
        if adj(i,j) ~= 0,
            W(i,j) = 1/degree(i);
        end
    end
end

end
