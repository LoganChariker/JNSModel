function [ combinations ] = rowCombineMany( cells )
%rowCombineMany Same as rowCombinations, except takes in a cell of matrices
%   Detailed explanation goes here
if isempty(cells)
    combinations = [];
elseif length(cells) == 1
    combinations = cells{1};
elseif length(cells) == 2
    combinations = rowCombinations(cells{1},cells{2});
else
    combinations = rowCombineMany( {rowCombinations(cells{1}, cells{2}), cells{3:end}} );
end

end

