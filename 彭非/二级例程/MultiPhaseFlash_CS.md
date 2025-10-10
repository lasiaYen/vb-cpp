# MultiPhaseFlash_CS

```
注：其中调用了OrderPhases函数，
源代码：Call OrderPhases(T, P, MW, nCompA, nPhaseMax, idx_Water, beta_mol, Zk, MWk, vol_cm3mol, dens_gcm3, yik, logPHI)

但OrderPhases的定义是：Sub OrderPhases(compositions, phi, beta, MW_Phase, mass_phase, density, Compr, nComp, nPhase)

compositions、Phi是二维数组，而传入的T和P均为标量，存在传入参数与函数参数列表不匹配的问题，仍未解决。
```


> code

```cpp
/**
 * @brief Wrapper: 输出新的 ShellFlash 到旧 SSP 格式
 *
 * @param ioSheet (void*)  输入：用于查找 EOS 参数的表格指针（在C里只是占位，具体实现交由调用方）
 * @param ioCol   (int)    输入：ioSheet中用于查找EOS参数的列号
 * @param eosProps(double**) 输入/输出：EOS参数矩阵 [nComp,6]
 * @param kij     (double**) 输入/输出：二元相互作用参数 [nComp,nComp]
 * @param TK      (double) 输入：温度 (K)
 * @param PBar    (double) 输入：压力 (bar)
 * @param total_moles (double) 输入：混合物总摩尔数
 * @param zInput  (double*) 输入：整体摩尔组成（旧顺序，长度 nCompMax）
 * @param gNeut   (double*) 输入：中性组分活度系数，gNeut[0]=CO2, gNeut[1]=H2S
 * @param aH2O    (double) 输入：水活度
 * @param density (double*) 输出：各相的密度 (g/cm^3)，长度 nPhaseMax
 * @param compositions (double**) 输出：组分组成矩阵 [nCompMax, nPhaseMax+1]（第一列为整体组成）
 * @param phi     (double**) 输出：逸度系数矩阵 [nCompMax, nPhaseMax]
 * @param Compr   (double*) 输出：各相压缩因子 Z
 * @param beta    (double*) 输出：摩尔相分率
 * @param zOutput (double*) 输出：整体摩尔组成
 * @param mass_phase (double*) 输出：各相质量 (g)
 * @param MW_Phase (double*) 输出：各相平均分子量 (g/mol)
 * @param No_Phases (int*) 输出：非零相的数量
 */
void MultiPhaseFlash_CS(int* ioSheet, int ioCol, double TK, double PBar,
    double* zInput, double aH2O, double* zOutput, double* mass_phase, double* MW_Phase, int* No_Phases)
{
    int i, j, k;
    int FlashType = 1;  /* modified Wilson PT-flash */
    int EOS = 2;        /* Peng-Robinson EOS */
    int SSPactive = 1;

    int nCompMax = 15;  /* 最大组分数 */

    /* ========= 输入预处理：调用外部函数 ========== */
    double* zi = NULL;
    int* idx_CompA = NULL, nCompA = 0;
    int* idx_Henry = NULL, idx_Water = -1;
    int* phaseID = NULL, nPhaseMax = 0;

    Flash_Input_Processing(ioSheet, ioCol, eosProps, kij, zInput, nCompMax,
        &zi, &idx_CompA, &nCompA, &idx_Henry, &idx_Water,
        &phaseID, &nPhaseMax);

    /* ========= 将输入转为C数组 ========= */
    double T = TK;
    double P = PBar;
    double* MW = (double*)malloc(nCompA * sizeof(double));
    double gamma[2];
    gamma[0] = gNeut[0];
    gamma[1] = gNeut[1];
    double aW = aH2O;

    for (i = 0; i < nCompA; i++) {
        MW[i] = eosProps[i][0]; /* VB的 eosProps(i,1) → C的 eosProps[i][0] */
    }

    /* ========= 调用 ShellFlash ========= */
    double* beta_mol = (double*)malloc(nPhaseMax * sizeof(double));
    double* Zk = (double*)malloc(nPhaseMax * sizeof(double));
    double** yik = (double**)malloc(nCompA * sizeof(double*));
    double** logPHI = (double**)malloc(nCompA * sizeof(double*));
    for (i = 0; i < nCompA; i++) {
        yik[i] = (double*)malloc(nPhaseMax * sizeof(double));
        logPHI[i] = (double*)malloc(nPhaseMax * sizeof(double));
    }

    ShellFlash(FlashType, EOS, T, P, zi, eosProps, kij, nCompA, nPhaseMax,
        phaseID, beta_mol, Zk, yik, logPHI,
        SSPactive, idx_Water, idx_Henry, aW, gamma);

    /* ========= 排序相并计算性质 ========= */
    int nPhase = 0;
    double MWk[3], vol_cm3mol[3], dens_gcm3[3];

    //为什么这里的参数和OrderPhases定义的参数不同？
    //C语言：OrderPhases(double** compositions, double** phi, double* beta, double* MW_Phase,
    //    double* mass_phase, double* density, double* Compr, int* nComp, int* nPhase)

    //vb原定义：Sub OrderPhases(compositions, phi, beta, MW_Phase, mass_phase, density, Compr, nComp, nPhase)
    //vb调用：Call OrderPhases(T, P, MW, nCompA, nPhaseMax, idx_Water, beta_mol, Zk, MWk, vol_cm3mol, dens_gcm3, yik, logPHI)
    OrderPhases(T, P, MW, nCompA, nPhaseMax, idx_Water, beta_mol, Zk,
        MWk, vol_cm3mol, dens_gcm3, yik, logPHI);

    /* ========= 输出到调用方变量 ========= */
    for (i = 0; i < nCompMax; i++) {
        compositions[i][0] = zInput[i];
        zOutput[i] = zInput[i];
    }

    for (k = 0; k < 3; k++) {
        if (beta_mol[k] > 0.0) {
            nPhase++;
            density[k] = dens_gcm3[k];
            Compr[k] = Zk[k];
            beta[k] = beta_mol[k];
            mass_phase[k] = beta_mol[k] * MWk[k] * total_moles;
            MW_Phase[k] = MWk[k];

            for (i = 0; i < nCompA; i++) {
                j = idx_CompA[i];
                compositions[j][k + 1] = yik[i][k];
                phi[j][k] = exp(logPHI[i][k]);
            }
        }
    }
    *No_Phases = nPhase;

    /* 计算气相体积 (m^3) —— 如果需要返回，需加到参数列表 */
    double mf_Vg = Compr[0] * total_moles * beta[0] * 0.00008314461 * TK / PBar;
    (void)mf_Vg; /* 防止未使用警告 */

    /* ========= 内存释放 ========= */
    free(MW);
    free(beta_mol);
    free(Zk);
    for (i = 0; i < nCompA; i++) {
        free(yik[i]);
        free(logPHI[i]);
    }
    free(yik);
    free(logPHI);
}

```

