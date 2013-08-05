VGS = (0 : 0.001 : 1);
I = zeros(size(VGS));

for i = 1 : numel(VGS)
   I(i) = current(4e-5, 1e-5, VGS(i), 1);
end