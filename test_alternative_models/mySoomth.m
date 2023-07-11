function Z = mySoomth(X, Bin)
X = [ones(Bin,1)*X(1); X; X(end)*ones(Bin, 1)];
X = conv(X, ones(1,Bin)/Bin, 'same');
X(1:Bin) = [];
X(end-Bin+1:end) = [];
Z = X;
end