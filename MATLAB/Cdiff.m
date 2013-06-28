dV1 = 0.001;
dV2 = 0.020;

pb_Cdiff_C_0 = zeros(size(pb_Qt_C_0));
mpb_Cdiff_C_0 = zeros(size(mpb_Qt_C_0));

for i = 1 : 5
    for j = 1 : 249
        pb_Cdiff_C_0(i, j + 1) = (pb_Qt_C_0(i, j + 1) - pb_Qt_C_0(i, j)) / dV1;
        mpb_Cdiff_C_0(i, j + 1) = (mpb_Qt_C_0(i, j + 1) - mpb_Qt_C_0(i, j)) / dV1;
    end
    for j = 250 : 286
        pb_Cdiff_C_0(i, j + 1) = (pb_Qt_C_0(i, j + 1) - pb_Qt_C_0(i, j)) / dV2;
        mpb_Cdiff_C_0(i, j + 1) = (mpb_Qt_C_0(i, j + 1) - mpb_Qt_C_0(i, j)) / dV2;
    end
end

pb_Cdiff_E_R = zeros(size(pb_Qt_E_R));
mpb_Cdiff_E_R = zeros(size(mpb_Qt_E_R));

for i = 1 : 5
    for j = 1 : 249
        pb_Cdiff_E_R(i, j + 1) = (pb_Qt_E_R(i, j + 1) - pb_Qt_E_R(i, j)) / dV1;
        mpb_Cdiff_E_R(i, j + 1) = (mpb_Qt_E_R(i, j + 1) - mpb_Qt_E_R(i, j)) / dV1;
    end
    for j = 250 : 286
        pb_Cdiff_E_R(i, j + 1) = (pb_Qt_E_R(i, j + 1) - pb_Qt_E_R(i, j)) / dV2;
        mpb_Cdiff_E_R(i, j + 1) = (mpb_Qt_E_R(i, j + 1) - mpb_Qt_E_R(i, j)) / dV2;
    end
end

pb_Cdiff_EFF = zeros(size(pb_Qt_EFF));
mpb_Cdiff_EFF = zeros(size(mpb_Qt_EFF));

for j = 1 : 249
    pb_Cdiff_EFF(1, j + 1) = (pb_Qt_EFF(1, j + 1) - pb_Qt_EFF(1, j)) / dV1;
end
for j = 250 : 286
    pb_Cdiff_EFF(1, j + 1) = (pb_Qt_EFF(1, j + 1) - pb_Qt_EFF(1, j)) / dV2;
end

for i = 1 : 5
    for j = 1 : 249
        mpb_Cdiff_EFF(i, j + 1) = (mpb_Qt_EFF(i, j + 1) - mpb_Qt_EFF(i, j)) / dV1;
    end
    for j = 250 : 286
        mpb_Cdiff_EFF(i, j + 1) = (mpb_Qt_EFF(i, j + 1) - mpb_Qt_EFF(i, j)) / dV2;
    end
end

pb_Cdiff_Z = zeros(size(pb_Qt_Z));
mpb_Cdiff_Z = zeros(size(mpb_Qt_Z));

for i = 1 : 2
    for j = 1 : 249
        pb_Cdiff_Z(i, j + 1) = (pb_Qt_Z(i, j + 1) - pb_Qt_Z(i, j)) / dV1;
        mpb_Cdiff_Z(i, j + 1) = (mpb_Qt_Z(i, j + 1) - mpb_Qt_Z(i, j)) / dV1;
    end
    for j = 250 : 286
        pb_Cdiff_Z(i, j + 1) = (pb_Qt_Z(i, j + 1) - pb_Qt_Z(i, j)) / dV2;
        mpb_Cdiff_Z(i, j + 1) = (mpb_Qt_Z(i, j + 1) - mpb_Qt_Z(i, j)) / dV2;
    end
end