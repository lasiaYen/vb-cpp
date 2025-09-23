# TrueFlash

## PstaH2O

> 代码部分

```cpp
/**
 * @brief 计算水的饱和蒸气压 (IAPWS-95)
 *
 * 根据 IAPWS 1995 标准，返回给定温度下水的饱和蒸气压 (单位: Bar)。
 *
 * @param TK 温度 [K]
 * @return double 水的饱和蒸气压 [Bar]
 */
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
```

> 测试代码

```cpp
void testPsatH2O() {
    double temps[] = { 300.0, 373.15, 500.0, 600.0 }; // 几个测试温度(K)
    int n = sizeof(temps) / sizeof(temps[0]);
    int i;

    for (i = 0; i < n; i++)
    {
        double TK = temps[i];
        double Psat = PsatH2O(TK);
        printf("T = %.2f K, Psat = %.6f Bar\n", TK, Psat);
    }

}
```

> 预期结果

```text
T = 300.00 K, Psat = 0.035367 Bar
T = 373.15 K, Psat = 1.014180 Bar
T = 500.00 K, Psat = 26.392227 Bar
T = 600.00 K, Psat = 123.448374 Bar
```

## InitialBeta

> 代码部分

```cpp
/**
 * @brief 初始化相分数 beta
 *
 * 此函数为多相平衡计算提供初始相分数的猜测值。
 *
 * 处理流程：
 * 1. 默认将前两个相的分数设为 0.5。
 * 2. 检查最后一个组分（通常为水）是否存在（z[NumGases-1] > 0）。
 * 3. 若存在水，判断温度与压力条件：
 *    - 当温度低于临界温度且压力高于水的饱和蒸气压时，水可能形成独立相。
 *    - 根据相数 (MaxBeta) 调整水的分配。
 * 4. 特殊情况：当仅有两个非零气体组分时，强制设定两个相各占 0.5。
 *
 * @param eqVapor   是否计算气相平衡 (未在此函数中使用)
 * @param eqAqueous 是否计算水相平衡 (未在此函数中使用)
 * @param TK        体系温度 (K)
 * @param PBar      体系压力 (bar)
 * @param NonZeroNumGases 非零组分数目
 * @param z         [输入] 组分摩尔分数数组，长度为 NumGases
 * @param TCr       [输入] 组分临界温度数组，长度为 NumGases
 * @param PCr       [输入] 组分临界压力数组，长度为 NumGases
 * @param MWgas     [输入] 组分摩尔质量数组，长度为 NumGases
 * @param beta      [输出] 相分数数组，长度为 MaxBeta
 * @param betaexists [输出] 相存在标记数组
 */
void InitialBeta(bool eqVapor, bool eqAqueous,
                 double TK, double PBar,
                 int NonZeroNumGases,
                 double *z, double *TCr, double *PCr, double *MWgas,
                 double *beta, int NumGases, int MaxBeta,
                 bool *betaexists)
{
    int m, j; // 源代码定义了，但未被调用
    double sumBeta;

    beta[0] = 0.5;
    beta[1] = 0.5;

    if (z[NumGases - 1] > 0)
    {
        if (TK < TCr[NumGases - 1] && (PBar > PsatH2O(TK)))
        {
            if (MaxBeta > 2)
            {
                beta[0] = (1 - z[NumGases - 1]) / 2.0;
                beta[1] = (1 - z[NumGases - 1]) / 2.0;
                beta[MaxBeta - 1] = z[NumGases - 1];
            }
            else
            {
                beta[0] = 1 - z[NumGases - 1];
                beta[MaxBeta - 1] = z[NumGases - 1];
            }
        }
        else
        {
            if (MaxBeta > 2)
            {
                beta[0] = (1 - z[NumGases - 1]) / 2.0 + z[NumGases - 1];
                beta[1] = (1 - z[NumGases - 1]) / 2;
                beta[MaxBeta - 1] = 0;
            }
            else
            {
                beta[0] = 1;
                beta[MaxBeta - 1] = 0;
            }
        }

        if (NonZeroNumGases == 2)
        {
            /* VB 中 ReDim beta(MaxBeta) 等价于重置，此处仅重新赋值 */
            beta[0] = 0.5;
            beta[MaxBeta - 1] = 0.5;
        }
    }

    sumBeta = 0;

    for (int j = 0; j < MaxBeta; j++)
    {
        sumBeta = beta[j] + sumBeta;
    }

    for (int j = 0; j < MaxBeta; j++)
    {
        beta[j] = beta[j] / sumBeta;
    }

    for (int j = 0; j < MaxBeta; j++)
    {
        if (beta[j] > 0)
        {
            betaexists[j] = true;
        }
    }
}
```

> 测试代码

```cpp
void testInitialBeta() {
    const int NumGases = 3;  // 例如：CO2, CH4, H2O
    const int MaxBeta = 3;

    double TCr[3] = { 304.2, 190.6, 647.1 };   // 临界温度
    double PCr[3] = { 73.8, 46.0, 220.6 };     // 临界压力
    double MWgas[3] = { 44.01, 16.04, 18.02 }; // 摩尔质量
    double z[3], beta[MaxBeta];
    bool betaexists[MaxBeta];

    // 测试 1：无水
    memset(beta, 0, sizeof(beta));
    memset(betaexists, 0, sizeof(betaexists));
    double z1[3] = { 0.5, 0.5, 0.0 };
    InitialBeta(true, true, 300, 50, 2, z1, TCr, PCr, MWgas, beta, NumGases, MaxBeta, betaexists);
    printf("Test1 (no water):beta = %.2f %.2f %.2f,betaexists=%d,%d,%d\n", beta[0], beta[1], beta[2], betaexists[0], betaexists[1], betaexists[2]);

    // 测试 2：有水，低温高压，MaxBeta=3 → 独立水相
    memset(beta, 0, sizeof(beta));
    memset(betaexists, 0, sizeof(betaexists));
    double z2[3] = { 0.3, 0.3, 0.4 };
    InitialBeta(true, true, 300, 50, 3, z2, TCr, PCr, MWgas, beta, NumGases, MaxBeta, betaexists);
    printf("Test2 (water, condense):beta = %.2f %.2f %.2f,betaexists=%d,%d,%d\n", beta[0], beta[1], beta[2], betaexists[0], betaexists[1], betaexists[2]);

    // 测试 3：有水，但条件不满足 (高温)
    memset(beta, 0, sizeof(beta));
    memset(betaexists, 0, sizeof(betaexists));
    double z3[3] = { 0.4, 0.4, 0.2 };
    InitialBeta(true, true, 700, 50, 3, z3, TCr, PCr, MWgas, beta, NumGases, MaxBeta, betaexists);
    printf("Test3 (water, mix):beta = %.2f %.2f %.2f,betaexists=%d,%d,%d\n", beta[0], beta[1], beta[2], betaexists[0], betaexists[1], betaexists[2]);

    // 测试 4：仅两个非零组分 → 强制 0.5, 0.5
    memset(beta, 0, sizeof(beta));
    memset(betaexists, 0, sizeof(betaexists));
    double z4[3] = { 0.5, 0.5, 0.0 };
    InitialBeta(true, true, 300, 50, 2, z4, TCr, PCr, MWgas, beta, NumGases, 2, betaexists);
    printf("Test4 (2 nonzero comps): beta = %.2f %.2f %.2f,betaexists=%d,%d,%d\n", beta[0], beta[1], beta[2], betaexists[0], betaexists[1], betaexists[2]);

    // 测试 5：极端情况，水占 100%
    memset(beta, 0, sizeof(beta));
    memset(betaexists, 0, sizeof(betaexists));
    double z5[3] = { 0.0, 0.0, 1.0 };
    InitialBeta(true, true, 300, 50, 1, z5, TCr, PCr, MWgas, beta, NumGases, MaxBeta, betaexists);
    printf("Test5 (pure water): beta = %.2f %.2f %.2f,betaexists=%d,%d,%d\n", beta[0], beta[1], beta[2], betaexists[0], betaexists[1], betaexists[2]);
}
```

> 预期结果

```cpp
Test1 (no water):           beta = 0.50 0.50 0.00,  betaexists=1,1,0
Test2 (water, condense):    beta = 0.30 0.30 0.40,  betaexists=1,1,1
Test3 (water, mix):         beta = 0.60 0.40 0.00,  betaexists=1,1,0
Test4 (2 nonzero comps):    beta = 0.50 0.50 0.00,  betaexists=1,1,0
Test5 (pure water):         beta = 0.00 0.00 1.00,  betaexists=0,0,1
```
