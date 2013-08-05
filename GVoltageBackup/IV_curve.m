VGS = (0.001 : 0.001 : 3);
I = zeros(size(VGS));

for i = 1 : numel(VGS)
    I(i) = current_b(5e-6, 3e-5, VGS(i), C_top_1nm_fit);
    disp(i);
end