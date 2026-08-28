function xSnv = applySnv(x)
rowMean = mean(x, 2);
rowStd = std(x, 0, 2);
rowStd(rowStd == 0) = 1;
xSnv = (x - rowMean) ./ rowStd;
end
