# CalcRhoTP     

> code

```cpp
double CalcRhoTP(double TK, double TC, double PBar, double Patm) {
    // 计算由于压力引起的过量摩尔体积变化，计算由于压力 1 atm 到 1001 atm 而引起的活度系数的导数
    // 按设定的TK值计算。如果不是STP条件（Patm !=1），则计算Δ压力=1e-6 bar时的过量性能

    double AphiP, AphiPPlus, X14, gX14, gpX14, X20, gX20, gpX20, X12, gX12, gpX12;
    double mt, Av, Fv, Ex_Pitzer, V_ex, V_ion, dens, VperKgWater, MassperKgwater;
    int m, a, c, n, iden;

    C2_Pitzer2019();
    AphiP = fAphicalc();

    // 计算各种X值和对应的gX、gpX
    X14 = 1.4 * sqrt(Ist); // For 2:(-2) pairs or ions
    gX14 = 2 * (1 - (1 + X14) * exp(-X14)) / (X14 * X14);
    gpX14 = -2 * (1 - (1 + X14 + 0.5 * X14 * X14) * exp(-X14)) / (X14 * X14);

    X20 = 2 * sqrt(Ist); // For 1:(-2), (-1):2, or 1:(-1) pairs
    gX20 = 2 * (1 - (1 + X20) * exp(-X20)) / (X20 * X20);
    gpX20 = -2 * (1 - (1 + X20 + 0.5 * X20 * X20) * exp(-X20)) / (X20 * X20);

    X12 = 12 * sqrt(Ist);
    gX12 = 2 * (1 - (1 + X12) * exp(-X12)) / (X12 * X12);
    gpX12 = -2 * (1 - (1 + X12 + 0.5 * X12 * X12) * exp(-X12)) / (X12 * X12);

    mt = fBtermcalc();

    // 保存当前bterm和CPhi值
    for (m = 0; m < NumCat; m++) {
        for (a = 0; a < NumAn; a++) {
            btermP[m][a] = bterm[m][a];
            CtermP[m][a] = CPhi[m][a];
        }
    }

    // 增加压力并重新计算
    PBar = PBar + 0.000001;
    Patm = PBar / 1.013254;
    // Ppsia = PBar * 14.503774; // 注释掉，如果不需要可以删除

    C2_Pitzer2019();

    AphiPPlus = fAphicalc();
    mt = fBtermcalc();

    // 保存增加压力后的值
    for (m = 0; m < NumCat; m++) {
        for (a = 0; a < NumAn; a++) {
            btermPPlus[m][a] = bterm[m][a];
            CtermPPlus[m][a] = CPhi[m][a];
        }
    }

    // 恢复原始压力
    PBar = PBar - 0.000001;
    Patm = PBar / 1.013254;
    // Ppsia = PBar * 14.503774;

    V0TP(); // 在T,P下重新计算V0

    // 计算体积相关项
    Av = -4 * RBar * TK * (AphiPPlus - AphiP) / 0.000001;
    Fv = Av / RBar / TK * Ist / 1.2 * log(1 + 1.2 * sqrt(Ist)); // Unit mol/Kg/bar

    // 计算bVterm和cVterm
    for (m = 0; m < NumCat; m++) {
        for (a = 0; a < NumAn; a++) {
            bVterm[m][a] = (btermPPlus[m][a] - btermP[m][a]) / 0.000001; // Unit Kg/mol/bar
            cVterm[m][a] = (CtermPPlus[m][a] - CtermP[m][a]) / 0.000001; // Unit (Kg/mol)^2/bar
        }
    }

    // 计算Pitzer超额项
    Ex_Pitzer = 0;
    for (m = 0; m < NumCat; m++) {
        for (a = 0; a < NumAn; a++) {
            Ex_Pitzer = Ex_Pitzer + 2 * mc[m] * ma[a] *
                (bVterm[m][a] + mc[m] * ChCat[m] * cVterm[m][a] /
                    (2 * sqrt(fabs(ChCat[m] * ChAn[a])))); // unit mol/Kg/bar
        }
    }

    // 计算超额体积和离子体积
    V_ex = (Fv + Ex_Pitzer) * RBar * TK; // L/Kg
    V_ion = V_ex * 1000; // cm3/Kg

    // 添加各离子的偏摩尔体积贡献
    for (a = 0; a < NumAn; a++) {
        V_ion = V_ion + ma[a] * V0_a[a]; // Note Unit of V0 is cm3/Kg
    }

    for (c = 0; c < NumCat; c++) {
        V_ion = V_ion + mc[c] * V0_c[c];
    }

    for (n = 0; n < NumNeut; n++) {
        V_ion = V_ion + mn[n] * V0_n[n];
    }

    // 计算密度
    // If useEOS = 1 Then
    //     dens = density(3) * 1000 'unit g/L
    // Else
    dens = fH2ODensity(TK, PBar); // g/L density of pure water at T, P
    // End If

    // 计算总体积和质量
    VperKgWater = (1.0 / dens * 1000000.0) + V_ion; // volume in cm3 of 1 kg water at T
    MassperKgwater = 1000; // g of water of 1 kg water

    // 添加各离子的质量贡献
    for (iden = 0; iden < NumCat; iden++) {
        MassperKgwater = MassperKgwater + mc[iden] * MWCat[iden];
    }

    for (iden = 0; iden < NumAn; iden++) {
        MassperKgwater = MassperKgwater + ma[iden] * MWAn[iden];
    }

    for (iden = 0; iden < NumNeut; iden++) {
        MassperKgwater = MassperKgwater + mn[iden] * MWNeut[iden]; // g
    }

    // 返回密度 (g/cm3)
    return MassperKgwater / VperKgWater;
}

```