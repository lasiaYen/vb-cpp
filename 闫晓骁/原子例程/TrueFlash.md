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
void test_PsatH2O() {
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
void test_InitialBeta() {
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

## myLimitIndex

> 代码部分

```cpp
/**
 * @brief 根据组分存在性筛选属性的最大或最小值索引
 *
 * 此函数在给定组分摩尔分数和属性值数组中，按照要求查找满足条件的
 * 最大值或最小值对应的索引。若未找到有效元素则返回 -1。
 *
 * 处理流程：
 * 1. 判断参数 max_or_min：
 *    - "max"：寻找属性值最大的组分索引；
 *    - "min"：寻找属性值最小的组分索引。
 * 2. 遍历 property 数组，对应的 composition[i] 必须大于阈值 tol (1e-6)，
 *    才认为该组分有效。
 * 3. 在有效组分中，记录当前最优值（最大或最小）及其索引。
 * 4. 返回找到的索引；若没有符合条件的组分，则返回 -1。
 *
 * @param max_or_min   字符串参数，取值 "max" 或 "min"，分别表示寻找最大值或最小值
 * @param composition  [输入] 组分摩尔分数数组，长度为 k
 * @param property     [输入] 对应的属性数组（如临界温度、压力等），长度为 k
 * @param k            数组长度，即组分总数
 *
 * @return long        满足条件的属性值对应索引；
 *                     若未找到有效元素，返回 -1
 */
long myLimitIndex(const char* max_or_min, const double* composition, const double* property, long k)
{
    long i;
    double lastProperty;
    double tol = 1e-6;
    long index = -1; /* 初始化返回值，未找到有效元素时返回 -1 */

    if (strcmp(max_or_min, "max") == 0)
    {
        lastProperty = -1.79769313486231E+307; /* 非常小的初始值 */
        for (i = 0; i < k; ++i)
        { /* C 数组从 0 开始 */
            if ((property[i] > lastProperty) && (composition[i] > tol))
            {
                lastProperty = property[i];
                index = i;
            }
        }
    }
    else if (strcmp(max_or_min, "min") == 0)
    {
        lastProperty = 1.79769313486231E+307; /* 非常大的初始值 */
        for (i = 0; i < k; ++i)
        {
            if ((property[i] < lastProperty) && (composition[i] > tol))
            {
                lastProperty = property[i];
                index = i;
            }
        }
    }

    return index;
}
```

> 测试代码

```cpp
void test_myLimitIndex()
{
    printf("==== Testing myLimitIndex ====\n");

    /* 1. 基础功能测试（max） */
    {
        const double comp[] = { 0.5, 0.3, 0.2 };
        const double prop[] = { 10.0, 20.0, 15.0 };
        long idx = myLimitIndex("max", comp, prop, 3);
        printf("Test 1 (max): result = %ld (expect 1)\n", idx);
    }

    /* 2. 基础功能测试（min） */
    {
        const double comp[] = { 0.5, 0.3, 0.2 };
        const double prop[] = { 10.0, 20.0, 15.0 };
        long idx = myLimitIndex("min", comp, prop, 3);
        printf("Test 2 (min): result = %ld (expect 0)\n", idx);
    }

    /* 3. 存在无效组分 */
    {
        const double comp[] = { 0.0, 1e-8, 0.5, 0.4 };
        const double prop[] = { 100.0, 200.0, 300.0, 250.0 };
        long idx = myLimitIndex("max", comp, prop, 4);
        printf("Test 3 (skip invalid): result = %ld (expect 2)\n", idx);
    }

    /* 4. 全部无效 */
    {
        const double comp[] = { 0.0, 1e-10, 5e-7 };
        const double prop[] = { 10.0, 20.0, 30.0 };
        long idx = myLimitIndex("min", comp, prop, 3);
        printf("Test 4 (all invalid): result = %ld (expect -1)\n", idx);
    }

    /* 5. 属性相等 */
    {
        const double comp[] = { 0.4, 0.4, 0.4 };
        const double prop[] = { 50.0, 50.0, 50.0 };
        long idx = myLimitIndex("max", comp, prop, 3);
        printf("Test 5 (equal values): result = %ld (expect 0)\n", idx);
    }

    /* 6. 单元素 */
    {
        const double comp[] = { 0.9 };
        const double prop[] = { 42.0 };
        long idx = myLimitIndex("min", comp, prop, 1);
        printf("Test 6 (single element): result = %ld (expect 0)\n", idx);
    }

    /* 7. 无效的参数 */
    {
        const double comp[] = { 0.5, 0.5 };
        const double prop[] = { 5.0, 10.0 };
        long idx = myLimitIndex("other", comp, prop, 2);
        printf("Test 7 (invalid param): result = %ld (expect -1)\n", idx);
    }

    /* 8. 负数属性 */
    {
        const double comp[] = { 0.3, 0.7 };
        const double prop[] = { -5.0, -10.0 };
        long idx = myLimitIndex("min", comp, prop, 2);
        printf("Test 8 (negative values): result = %ld (expect 1)\n", idx);
    }

    printf("==== Testing Finished ====\n");
}
```

> 预期结果

```text
Test 1 (max):               result = 1 (expect 1)
Test 2 (min):               result = 0 (expect 0)
Test 3 (skip invalid):      result = 2 (expect 2)
Test 4 (all invalid):       result = -1 (expect -1)
Test 5 (equal values):      result = 0 (expect 0)
Test 6 (single element):    result = 0 (expect 0)
Test 7 (invalid param):     result = -1 (expect -1)
Test 8 (negative values):   result = 1 (expect 1)
```

## Initialization

> 代码部分

Initialization主要用来修改logphi3phase的值，这里通过函数参数传入，对其进行原地修改，且logphi3phase为二维数组。

```cpp
/**
 * @brief 初始化 logphi3phase 数组（fugacity coefficient 对数）
 *
 * 为多相平衡计算提供初始对数 fugacity coefficient。
 * 使用 Wilson 公式计算初始 Ki，并根据水和烃类分子量做启发式修正。
 *
 * @param zGlobalWater 水的摩尔分数
 * @param NonZeroNumGases 非零气体组分数量
 * @param MaxBeta 最大相数
 * @param EOS 状态方程类型 ("SRK" 或 "PR")
 * @param MWgas[] 摩尔质量数组，长度 NumGases
 * @param TCr[] 临界温度数组，长度 NumGases
 * @param PCr[] 临界压力数组，长度 NumGases
 * @param Omega[] 非球形因子数组，长度 NumGases
 * @param TK 温度 [K]
 * @param PBar 压力 [Bar]
 * @param logphi3phase[][] 输出数组，大小 [NumGases][MaxBeta]
 * @param NumGases 气体组分数量
 */
void Initialization(double zGlobalWater, int NonZeroNumGases, int MaxBeta, const char* EOS,
    const double* MWgas, const double* TCr, const double* PCr, const double* Omega,
    double TK, double PBar, double** logphi3phase, int NumGases)
{
    int i, j;
    double correction;
    double* lnKi = (double*)malloc(NumGases * sizeof(double));

    /* Wilson approximation for phase equilibrium constant */
    for (i = 0; i < NumGases; ++i)
    {
        lnKi[i] = log(PCr[i] / PBar) + 5.373 * (1.0 + Omega[i]) * (1.0 - TCr[i] / TK);
    }

    /* Initialize logphi3phase */
    for (i = 0; i < NumGases; ++i)
    {
        logphi3phase[i][0] = 0.0; /* phase 1 = vapor (ideal gas) */

        for (j = 1; j < MaxBeta; ++j)
        {
            logphi3phase[i][j] = lnKi[i]; /* phase j = hydrocarbon liquids */

            if (zGlobalWater > 0.0 && NonZeroNumGases == 2 && MWgas[i] >= 210.0)
            {
                logphi3phase[i][j] = lnKi[i] + 10.0;
            }
        }
    }

    /* Heuristic correction if water is present */
    if (zGlobalWater > 0.0)
    {

        if (MaxBeta > 2 && NonZeroNumGases > 2)
        {
            for (j = 1; j < MaxBeta - 1; ++j)
            {
                logphi3phase[NumGases - 1][j] = lnKi[NumGases - 1] + 5.0; /* Water in organic liquids */
            }
        }

        for (i = 0; i < NumGases - 1; ++i)
        {
            if (MWgas[i] < 45.0)
            {
                logphi3phase[i][MaxBeta - 1] = lnKi[i] + 5.0; /* Light species in aqueous phase */
            }
            else
            {
                correction = MWgas[i] / 3.8261;
                logphi3phase[i][MaxBeta - 1] = lnKi[i] + correction; /* Heavy species in aqueous phase */
            }
        }
    }

    /* Minimum value threshold to avoid numerical instability */
    for (j = 0; j < MaxBeta; ++j)
    {
        for (i = 0; i < NumGases; ++i)
        {
            if (logphi3phase[i][j] < -20.0)
            {
                logphi3phase[i][j] = -20.0;
            }
        }
    }

    free(lnKi);
}
```

> 测试代码

```cpp
void test_Initialization()
{
    printf("==== Testing Initialization ====\n");

    /* 测试参数设置 */
    int NumGases = 3;
    int MaxBeta = 3;
    double zGlobalWater = 0.1;       /* 水存在 */
    int NonZeroNumGases = 3;         /* 非零气体组分数 */
    const char* EOS = "SRK";         /* 状态方程 */
    double MWgas[3] = { 16.0, 44.0, 210.0 };  /* 气体摩尔质量 */
    double TCr[3] = { 190.0, 305.0, 647.0 }; /* 临界温度 */
    double PCr[3] = { 46.0, 73.0, 220.0 };   /* 临界压力 */
    double Omega[3] = { 0.011, 0.099, 0.344 }; /* 非球形因子 */
    double TK = 298.15;  /* 温度 */
    double PBar = 50.0;  /* 压力 */

    /* 动态分配 logphi3phase 二维数组 */
    double** logphi3phase = (double**)malloc(NumGases * sizeof(double*));
    for (int i = 0; i < NumGases; ++i)
    {
        logphi3phase[i] = (double*)malloc(MaxBeta * sizeof(double));
        for (int j = 0; j < MaxBeta; ++j)
            logphi3phase[i][j] = 0.0; /* 初始化为 0 */
    }

    /* 调用被测函数 */
    Initialization(zGlobalWater, NonZeroNumGases, MaxBeta, EOS, MWgas, TCr, PCr, Omega, TK, PBar, logphi3phase, NumGases);

    /* 打印结果 */
    for (int i = 0; i < NumGases; ++i)
    {
        printf("Component %d: ", i);
        for (int j = 0; j < MaxBeta; ++j)
        {
            printf("%8.4f ", logphi3phase[i][j]);
        }
        printf("\n");
    }

    /* 释放内存 */
    for (int i = 0; i < NumGases; ++i)
        free(logphi3phase[i]);
    free(logphi3phase);

    printf("==== Testing Finished ====\n");
```

> 预期结果

```text
Component 0:   0.0000   1.8870   6.8870
Component 1:   0.0000   0.2428   5.2428
Component 2:   0.0000  -1.9677  -6.9677
```
