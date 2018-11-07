function gp = getGridPosition (index, width)
rowIndex = ceil(index / width);
colIndex = mod(index, width);
gp(1,1) = rowIndex;
gp(1,2) = colIndex;
end