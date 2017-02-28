function [ F ] = insulated_diff( W, b )

    temp = eye(length(W)) - (1-b)*W;
    F = b*inv(temp);

end
