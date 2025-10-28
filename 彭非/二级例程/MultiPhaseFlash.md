# MultiPhaseFlash

> code

注：里面调用了TureFlash，参数列表如何传入等待决定

```cpp
/**
 * @brief  多相闪蒸计算
 *
 *  多相闪蒸计算的总控制函数
 *
 * @param   mf_ParametersWereRead
 * @param   TCr
 * @param   PCr
 * @param   Omega
 * @param   MWgas
 * @param   kPr
 * @param   c0
 * @param   c1
 * @param   TK                  输入 - double, 温度，[开尔文]
 * @param   PBar                输入 - double, 压力，[巴]
 * @param   total_moles         输入 - double, 包括水在内的储层流体每天的总摩尔数，[摩尔]
 * @param   z                   输入 - double数组, 包括水在内的储层流体的摩尔组成，[无量纲]
 * @param   gNeut               输入 - double数组, 水相中 CH4、CO2、H2S 和 H2O 的活度系数
 * @param   aH2O                输入 - double, 水相中水的活度； mf_Vg 气体体积 (m^3)
 * @param   density             输出 - double数组, 各相的密度（气相、有机相和水相），[g/cm3]
 * @param   compositions        输出 - double**, 整个体系的摩尔组成，行：组分，列：相，[无量纲]
 * @param   phi                 输出 - double**, 整个体系的逸度系数，行：组分，列：相，[无量纲]
 * @param   Compr               输出 - double数组, 各相的压缩性（气相、有机相和水相），[无量纲]
 * @param   beta                输出 - double数组, 各相的摩尔相分数（气相、有机相和水相），[无量纲]
 * @param   zOutput             输出 - double数组, 整体体系的摩尔相分数构图，[无量纲]
 * @param   mass_phase          输出 - double数组, 每日各相的质量（气相、有机相和水相），[kg]
 * @param   MW_Phase            输出 - 各相的平均分子量
 * @param   No_Phases           输出 - 实际存在的相的个数
 */
void MultiPhaseFlash(bool* mf_ParametersWereRead, double* TCr, double* PCr, double* Omega, double* MWgas, double** kPr, double* c0, double* c1,
    double TK, double PBar, double total_moles, double* z, double* gNeut, double aH2O, double* density, double** compositions, double** phi,
    double* Compr, double* beta, double* zOutput, double** mass_phase, double** MW_Phase, int* No_Phases)
{

    char EOS[] = "PR"; //建立了状态方程为 Peng-Robinson: PR，或 Soave-RK: SRK
    int max_NumGases = 15; //这确定了组件的数量和可能的最大阶段数

    //变量定义
    int NonZeroNumGases = 0, mf_NumGases = 6, MaxBeta, i, j, k, L, m;
    double mf_Vg, mf_mass_w_GO;
    int counterEquilibrium_final, counter_final;
    bool eqAqueous, eqVapor;
    double mf_nTCO2_GO, mf_nTH2S_GO;
    int iPure;
    int Numbeta;
    int pure;
    double Sz;
    double* logphipure;
    double PsatPure;
    char* phase;
    char** phaseName;
    double logphipureL, ComprL, logphipureV, ComprV;


    //这表明，如果没有水，就不应该存在水相
    z[max_NumGases - 1] > 0 ? MaxBeta = 3 : MaxBeta = 2;
    //ReDim mf_gNeut(UBound(gNeut)), zInput(max_NumGases)
    double* mf_gNeut = (double*)malloc(max_NumGases * sizeof(double));
    double* zInput = (double*)malloc(max_NumGases * sizeof(double));

    //将值从外部变体类型变量输出传递到双精度类型变量
    double mf_TK = TK, mf_PBar = PBar, mf_aH2O = aH2O;

    for (i = 0; i < max_NumGases; i++)
    {
        zInput[i] = z[i];
        mf_gNeut[i] = gNeut[i];
    }

    for (i = 0; i < max_NumGases - 2; i++)
    {
        if (zInput[i] > 0)
        {
            NonZeroNumGases++;
            pure = i;
            if (i > 3)
                mf_NumGases++;
        }
    }
    if (zInput[max_NumGases - 2] > 0)
    {
        NonZeroNumGases++;//Counting Nitrogen, N2
        pure = max_NumGases - 1;
    }
    if (zInput[max_NumGases - 1] > 0)
    {
        NonZeroNumGases++;//Counting Water
        pure = max_NumGases;
    }
    if (zInput[max_NumGases - 1] > 0 && NonZeroNumGases == 2 && (zInput[0] > 0 || zInput[1] > 0 || zInput[2] > 0 ||
        zInput[3] > 0 || zInput[max_NumGases - 2] > 0))
    {
        MaxBeta = 2;
    }
    if (NonZeroNumGases == 1)MaxBeta = 3;
    //赋值后，就可以开辟一些需要的数组了。

    double* zGlobal = (double*)malloc(mf_NumGases * sizeof(double));
    int* iFlash = (int*)malloc(max_NumGases * sizeof(int));
    double* zOut = (double*)malloc(mf_NumGases * sizeof(double));

    double* mf_beta = (double*)malloc(MaxBeta * sizeof(double));
    double* mf_Compr = (double*)malloc(MaxBeta * sizeof(double));
    double* mf_density = (double*)malloc(MaxBeta * sizeof(double));

    double** mf_phi = (double**)malloc(mf_NumGases * sizeof(double*));
    double** mf_compositions = (double**)malloc(mf_NumGases * sizeof(double*));
    for (i = 0; i < mf_NumGases; i++) {
        mf_phi[i] = (double*)malloc(MaxBeta * sizeof(double));
        mf_compositions[i] = (double*)malloc((MaxBeta + 1) * sizeof(double));
    }
    //ReDim mf_Compr(MaxBeta), mf_density(MaxBeta), mass_phase(MaxBeta)，注，参数列表参数mass_phase的ReDim行为
    free(*mass_phase);
    *mass_phase = (double*)malloc(MaxBeta * sizeof(double));
    double* x = (double*)malloc(mf_NumGases * sizeof(double));

    //获取参数值、缩小组成以适应现有化合物并使其标准化的子程序
    *mf_ParametersWereRead = false;
    //Call InitialPreparationSSP(mf_ParametersWereRead, EOS, zInput, iFlash, zGlobal, TCr, PCr, Omega, MWgas, kPr, c0, c1)
    InitialPreparationSSP(mf_ParametersWereRead, EOS, zInput, iFlash, zGlobal, &TCr, &PCr, &Omega,
        &MWgas, &kPr, &c0, &c1, max_NumGases, 15);

    int counterEquilibrium_final = 0, iter_final = 0, counter_final = 0;//初始化防止传入空值

    //计算多次闪蒸平衡的主要子程序
    if (NonZeroNumGases > 1)
    {
        Numbeta = 0;//原代码无这个，这里赋值防止传入未初始化的值。
        phaseName = NULL;//原代码无这个，这里赋值防止传入未初始化的值。

        //怎么调用待定
        TrueFlash(EOS, mf_TK, mf_PBar, NonZeroNumGases, zGlobal, mf_gNeut, mf_aH2O, TCr, PCr, Omega, MWgas, kPr, c0, c1, &Numbeta, phaseName, mf_beta, mf_compositions, mf_phi, mf_Compr, mf_density,
            &counterEquilibrium_final, &iter_final, &counter_final, zOut);
    }
    else if (NonZeroNumGases == 1)//纯组分情况下的性质计算
    {
        Numbeta = 1;
        for (i = 0; i < mf_NumGases; i++)
            if (iFlash[i] == pure)iPure = i;
        if (pure < 3 || (pure == max_NumGases))
        {
            phase = (char*)malloc(6); // 分配足够容纳"vapor"及结束符'\0'的空间
            strcpy(phase, "vapor");
            k = 1;
            if (TK < TCr[iPure])
            {
                if (pure == 1)PsatPure = 1;//CH4的位置
                if (pure == 2)PsatPure = PsatCO2(mf_TK);//CO2的位置
                if (pure == max_NumGases)PsatPure = PsatH2O(mf_TK);//水的位置，因为默认情况下水是最后一个

                if (PBar > PsatPure && PsatPure != 0)
                {
                    free(phase);
                    phase = (char*)malloc(7);
                    strcpy(phase, "liquid");
                    if (pure = max_NumGases)k = 3;
                }
            }
        }
        phi_calc(false, false, EOS, phase, mf_TK, mf_PBar, zGlobal, zGlobal, mf_gNeut, mf_aH2O, TCr, PCr,
            Omega, c0, c1, kPr, logphipure, &mf_Compr[k], NonZeroNumGases);//NumGases
        if (strcmp(phase, "vapor") == 0)
            Compr[k] = ComprV;
        else
        {
            phi_calc(false, false, EOS, "vapor", mf_TK, mf_PBar, zGlobal, zGlobal, mf_gNeut, mf_aH2O, TCr, PCr,
                Omega, c0, c1, kPr, logphipure, &ComprV, NonZeroNumGases);//NumGases
            logphipureV = logphipure[iPure];
            phi_calc(false, false, EOS, "liquid", mf_TK, mf_PBar, zGlobal, zGlobal, mf_gNeut, mf_aH2O, TCr, PCr,
                Omega, c0, c1, kPr, logphipure, &ComprL, NonZeroNumGases);//NumGases
            logphipureL = logphipure[iPure];

            if (logphipureV > logphipureL)
            {
                k = 2;
                mf_Compr[k] = ComprL;
                logphipure[iPure] = logphipureL;
            }
            else
            {
                k = 1;
                mf_Compr[k] = ComprV;
                logphipure[iPure] = logphipureV;
            }
        }
        mf_phi[iPure][k - 1] = exp(logphipure[iPure]);//Exp
        mf_compositions[iPure][k] = 1.0;

        //相密度计算
        if (mf_Compr[k - 1] > 0)
        {
            mf_density[k - 1] = MWgas[iPure] * PBar / RBar / mf_Compr[k - 1] / TK;
            mf_beta[k - 1] = 1;
        }
        else
        {
            mf_density[k - 1] = 0.0;
            mf_beta[k - 1] = 0.0;
        }
    }

    if (MaxBeta < 3)
    {
        if (zInput[max_NumGases - 1] > 0 && NonZeroNumGases == 2 && (zInput[0] > 0 || zInput[1] > 0 || zInput[2] > 0 || zInput[3] > 0 ||
            zInput[max_NumGases - 2] > 0))
        {
            //ReDim Preserve mf_beta(MaxBeta), mf_density(MaxBeta), mf_Compr(MaxBeta)
            // Dim mf_compositions_temp() As Double, mf_phi_temp() As Double
            // ReDim mf_compositions_temp(mf_NumGases, MaxBeta + 1), mf_phi_temp(mf_NumGases, MaxBeta)        
            //在MaxBeta获得新的值之前，先保存原来的数组
            int temMaxBeta = MaxBeta;
            MaxBeta = 3;
            //创建新空间
            double* temp_mf_beta = (double*)malloc(temMaxBeta * sizeof(double));
            double* temp_mf_density = (double*)malloc(temMaxBeta * sizeof(double));
            double* temp_mf_Compr = (double*)malloc(temMaxBeta * sizeof(double));
            //搬运值过去
            for (i = 0; i < temMaxBeta; i++)
            {
                temp_mf_beta[i] = mf_beta[i];
                temp_mf_density[i] = mf_density[i];
                temp_mf_Compr[i] = mf_Compr[i];
            }
            //销毁原空间，指针指向新的空间
            free(mf_beta); free(mf_density); free(mf_Compr);
            mf_beta = (double*)malloc(MaxBeta * sizeof(double));
            mf_density = (double*)malloc(MaxBeta * sizeof(double));
            mf_Compr = (double*)malloc(MaxBeta * sizeof(double));

            for (i = 0; i < temMaxBeta; i++)
            {
                mf_beta[i] = temp_mf_beta[i];
                mf_density[i] = temp_mf_density[i];
                mf_Compr[i] = temp_mf_Compr[i];
            }
            free(temp_mf_beta); free(temp_mf_density); free(temp_mf_Compr);

            double** mf_compositions_temp = (double**)malloc(mf_NumGases * sizeof(double*));
            double** mf_phi_temp = (double**)malloc(mf_NumGases * sizeof(double*));
            for (i = 0; i < mf_NumGases; i++)
            {
                mf_compositions_temp[i] = (double*)malloc((MaxBeta + 1) * sizeof(double));
                mf_phi_temp[i] = (double*)malloc((MaxBeta) * sizeof(double));
            }

            if (zInput[max_NumGases] > 0)
            {
                mf_beta[2] = mf_beta[1];
                mf_beta[1] = 0;
                mf_density[2] = mf_density[1];
                mf_density[1] = 0;
                mf_Compr[2] = mf_Compr[1];
                mf_Compr[1] = 0;
            }
            for (i = 0; i < mf_NumGases; i++)
            {
                for (j = 0; j < MaxBeta - 1; j++)
                {
                    mf_compositions_temp[i][j + 1] = mf_compositions[i][j + 1];
                    mf_phi_temp[i][j] = mf_phi[i][j];
                }
            }

            //ReDim mf_compositions(mf_NumGases, MaxBeta + 1), mf_phi(mf_NumGases, MaxBeta)
            // 重新分配主数组
            for (i = 0; i < mf_NumGases; i++)
            {
                free(mf_compositions[i]);
                free(mf_phi[i]);
            }
            free(mf_compositions);
            free(mf_phi);
            mf_compositions = (double**)malloc(mf_NumGases * sizeof(double*));
            mf_phi = (double**)malloc(mf_NumGases * sizeof(double*));

            for (i = 0; i < mf_NumGases; i++)
            {
                mf_compositions[i] = (double*)malloc((MaxBeta + 1) * sizeof(double));
                mf_phi[i] = (double*)malloc((MaxBeta) * sizeof(double));
            }

            // 从临时数组复制回主数组
            for (i = 0; i < mf_NumGases; i++)
            {
                if (zInput[max_NumGases] > 0)
                {
                    mf_compositions[i][3] = mf_compositions_temp[i][2];
                    mf_compositions[i][1] = mf_compositions_temp[i][1];
                    mf_phi[i][0] = mf_phi_temp[i][0];
                    mf_phi[i][2] = mf_phi_temp[i][1];
                }
                else {
                    for (j = 0; j < MaxBeta - 1; j++)
                    {
                        mf_compositions[i][j] = mf_compositions_temp[i][j];
                        mf_phi[i][j] = mf_phi_temp[i][j];
                    }
                }
            }
            // 释放临时数组
            for (i = 0; i < mf_NumGases; i++)
            {
                free(mf_compositions_temp[i]);
                free(mf_phi_temp[i]);
            }
            free(mf_compositions_temp);
            free(mf_phi_temp);
        }
    }
    //ReDim MW_Phase(MaxBeta), mass_phase(MaxBeta)
    *MW_Phase = (double*)realloc(*MW_Phase, MaxBeta * sizeof(double));
    *mass_phase = (double*)realloc(*mass_phase, MaxBeta * sizeof(double));

    //计算分子量、相质量
    //将双精度型变量输出的值传递到外部变量类型变量
    //并放大到组件的原始尺寸
    for (j = 0; j < MaxBeta; j++)
    {
        for (i = 0; i < mf_NumGases; i++)
        {
            x[i] = mf_compositions[i][j + 1];
            compositions[iFlash[i]][j + 1] = mf_compositions[i][j + 1];
            phi[iFlash[i]][j] = mf_phi[i][j];
            zOutput[iFlash[i]] = zOut[i];
        }

        beta[j] = mf_beta[j];
        (*MW_Phase)[j] = weighted_mean(x, MWgas, mf_NumGases);
        (*mass_phase)[j] = total_moles * mf_beta[j] * (*MW_Phase)[j] / 1000;
        density[j] = mf_density[j];
        Compr[j] = mf_Compr[j];
    }
    // 重新排序输出（CALEB）
    OrderPhases(&compositions, &phi, &beta, MW_Phase, mass_phase, &density, &Compr, &max_NumGases, &MaxBeta);
    *No_Phases = Numbeta;
    mf_Vg = Compr[0] * total_moles * beta[0] * 0.00008314461 * mf_TK / PBar; // 气体体积 m^3
}
```

> MultiPhaseFlash调用的函数

1、TrueFlash(闫晓晓师兄实现)
2、PsatCO2
3、PsatH2O
4、weighted_mean
5、OrderPhases （一级例程-彭非）

```cpp
double PsatCO2(double TK)
{
    const double TC = 304.1282;   // 临界温度 [K]
    const double Pc = 73.773;     // 临界压力 [bar]
    const double a1 = -7.0602087;
    const double a2 = 1.9391218;
    const double a3 = -1.6463597;
    const double a4 = -3.2995634;
    const double t1 = 1.0;
    const double t2 = 1.5;
    const double t3 = 2.0;
    const double t4 = 4.0;

    double Theta = 1.0 - TK / TC;
    double exponent = (TC / TK) * (a1 * pow(Theta, t1) +
        a2 * pow(Theta, t2) +
        a3 * pow(Theta, t3) +
        a4 * pow(Theta, t1));  // 注意：VB 代码写的是 t1，这里保持一致
    return Pc * exp(exponent);
}


double PsatH2O(double TK)
{
    const double TC = 647.096; /* 临界温度 K */
    const double Pc = 220.64;  /* 临界压力 Bar */
    const double a1 = -7.85951783;
    const double a2 = 1.84408259;
    const double a3 = -11.7866497;
    const double a4 = 22.6807411;
    const double a5 = -15.9618719;
    const double a6 = 1.80122502;

    double T = 1.0 - TK / TC;

    double Psat = Pc * exp((TC / TK) * (a1 * T + a2 * pow(T, 1.5) + a3 * pow(T, 3.0) +
                                        a4 * pow(T, 3.5) + a5 * pow(T, 4.0) + a6 * pow(T, 7.5)));
    return Psat;
}


double weighted_mean(double* weight, double* property, int n)
{
    double sum = 0.0;
    double sumW = 0.0;
    for (int i = 0; i < n; i++) {
        sum += weight[i] * property[i];
        sumW += weight[i];
    }
    if (sumW > 0.0)
        return sum / sumW;
    else
        return 0.0;
}

```

