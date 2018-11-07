function [min, index] = getMinimum(m, totalRows)
min(1,1) = m(1,totalRows);
index = 1;
for i = 2:length(m)-1
    if m(i,totalRows) < min(1,1)
        min(1,1) = m(i,totalRows);
        index = i;
    end
end
end