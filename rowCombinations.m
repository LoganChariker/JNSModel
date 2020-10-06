function [ XY ] = rowCombinations( X, Y )
%rowCombinations Gives all combinations of rows of X with rows of Y
%   This function takes an m-by-n matrix X and an r-by-s matrix Y, and
%   outputs an mr-by-(n+s) matrix XY, where the rows of XY give all the
%   combinations of rows of X with rows of Y
[m, n] = size(X);
[r, s] = size(Y);

XY = zeros(m*r, n+s);

for i = 1 : r
    for j = 1 : m
        XY((i-1)*m + j,:) = [X(j,:), Y(i,:)];
    end
end

end

