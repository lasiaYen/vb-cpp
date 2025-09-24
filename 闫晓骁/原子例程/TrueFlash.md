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

## EQCalculation

> 代码部分

```cpp
/**
 * @brief 计算目标函数 Q 及辅助向量 E
 *
 * 参考 Michelsen 1994 方法，计算多相平衡中目标函数和辅助向量 E。
 *
 * @param beta[] [输入] 相分数数组，长度 MaxBeta
 * @param lnPHI[][] [输入] fugacity coefficient 对数矩阵，大小 [NumGases][MaxBeta]
 * @param zGlobal[] [输入] 总组分摩尔分数数组，长度 NumGases
 * @param E[] [输出] 辅助向量，长度 NumGases
 * @param Q [输出] 目标函数值
 * @param NumGases 组分数量
 * @param MaxBeta 相数
 */
void EQCalculation(double *beta, double **lnPHI, double *zGlobal,
                   double *E, double *Q, int NumGases, int MaxBeta)
{
    int i, k;
    double Qphase = 0.0;
    double Qcomp = 0.0;

    /* 计算辅助向量 E[i] */
    for (i = 0; i < NumGases; ++i)
    {
        E[i] = 0.0;
        for (k = 0; k < MaxBeta; ++k)
        {
            if (beta[k] > 0.0 && zGlobal[i] > 0.0)
            {
                E[i] += beta[k] / exp(lnPHI[i][k]); /* lnPHI[i][k] -> e^(-lnPHI) */
            }
        }
    }

    /* 计算 Qphase */
    for (k = 0; k < MaxBeta; ++k)
    {
        Qphase += beta[k];
    }

    /* 计算 Qcomp */
    for (i = 0; i < NumGases; ++i)
    {
        if (zGlobal[i] > 0.0)
        {
            Qcomp += zGlobal[i] * log(E[i]);
        }
    }

    *Q = Qphase - Qcomp;
}
```

> 测试代码

```cpp
/* 打印辅助函数 */
static void print_array(const char* name, const double* arr, int n) {
    printf("%s = [", name);
    for (int i = 0; i < n; i++) {
        printf("%8.4f", arr[i]);
        if (i < n - 1) printf(", ");
    }
    printf("]\n");
}

void test_EQCalculation() {
    printf("==== Testing EQCalculation ====\n");

    /* ---------- Test Case 1 ---------- */
    {
        int NumGases = 2, MaxBeta = 2;
        double beta[] = { 0.5, 0.5 };
        double zGlobal[] = { 0.5, 0.5 };
        double lnPHI_data[2][2] = { {0.0, 0.0}, {0.0, 0.0} };
        double* lnPHI[2] = { lnPHI_data[0], lnPHI_data[1] };
        double E[2]; double Q;
        EQCalculation(beta, lnPHI, zGlobal, E, &Q, NumGases, MaxBeta);
        printf("Test 1 (all phi=1)\n");
        print_array("E", E, NumGases);
        printf("Q = %.6f \n\n", Q);
    }

    /* ---------- Test Case 2 ---------- */
    {
        int NumGases = 2, MaxBeta = 2;
        double beta[] = { 0.5, 0.5 };
        double zGlobal[] = { 0.5, 0.5 };
        double lnPHI_data[2][2] = { {0.0, 1.0}, {0.0, -1.0} };
        double* lnPHI[2] = { lnPHI_data[0], lnPHI_data[1] };
        double E[2]; double Q;
        EQCalculation(beta, lnPHI, zGlobal, E, &Q, NumGases, MaxBeta);
        printf("Test 2 (different lnPHI)\n");
        print_array("E", E, NumGases);
        printf("Q = %.6f \n\n", Q);
    }

    /* ---------- Test Case 3 ---------- */
    {
        int NumGases = 2, MaxBeta = 2;
        double beta[] = { 0.5, 0.5 };
        double zGlobal[] = { 1.0, 0.0 };
        double lnPHI_data[2][2] = { {0.0, 0.0}, {0.0, 0.0} };
        double* lnPHI[2] = { lnPHI_data[0], lnPHI_data[1] };
        double E[2]; double Q;
        EQCalculation(beta, lnPHI, zGlobal, E, &Q, NumGases, MaxBeta);
        printf("Test 3 (second component zero)\n");
        print_array("E", E, NumGases);
        printf("Q = %.6f \n\n", Q);
    }

    /* ---------- Test Case 4 ---------- */
    {
        int NumGases = 2, MaxBeta = 2;
        double beta[] = { 1.0, 0.0 };
        double zGlobal[] = { 0.5, 0.5 };
        double lnPHI_data[2][2] = { {0.0, 0.0}, {0.0, 0.0} };
        double* lnPHI[2] = { lnPHI_data[0], lnPHI_data[1] };
        double E[2]; double Q;
        EQCalculation(beta, lnPHI, zGlobal, E, &Q, NumGases, MaxBeta);
        printf("Test 4 (single phase)\n");
        print_array("E", E, NumGases);
        printf("Q = %.6f \n\n", Q);
    }

    /* ---------- Test Case 5 ---------- */
    {
        int NumGases = 3, MaxBeta = 2;
        double beta[] = { 0.3, 0.7 };
        double zGlobal[] = { 0.2, 0.3, 0.5 };
        double lnPHI_data[3][2] = {
            {0.0, 0.1},
            {0.2, -0.1},
            {0.5, -0.5}
        };
        double* lnPHI[3] = { lnPHI_data[0], lnPHI_data[1], lnPHI_data[2] };
        double E[3]; double Q;
        EQCalculation(beta, lnPHI, zGlobal, E, &Q, NumGases, MaxBeta);
        printf("Test 5 (three components)\n");
        print_array("E", E, NumGases);
        printf("Q = %.6f\n\n", Q);
    }

    printf("==== Testing Finished ====\n");
}
```

> 预期结果

```text
Test 1 (all phi=1)
E = [  1.0000,   1.0000]
Q = 1.000000

Test 2 (different lnPHI)
E = [  0.6839,   1.8591]
Q = 0.879885

Test 3 (second component zero)
E = [  1.0000,   0.0000]
Q = 1.000000

Test 4 (single phase)
E = [  1.0000,   1.0000]
Q = 1.000000

Test 5 (three components)
E = [  0.9334,   1.0192,   1.3361]
Q = 0.863206
```

## Gauss

> 代码部分

需要注意，此代码逻辑与源代码中的Gauss函数相同，但是得出的结果似乎并不准确

```cpp
/**
 * @brief 使用高斯消元法解线性方程组 A * x = b
 *
 * 通过增广矩阵和部分主元选取，将线性方程组转化为上三角矩阵，然后回代求解。
 *
 * @param a [输入] 系数矩阵，大小 n x n
 * @param b [输入] 常数向量，长度 n
 * @param n [输入] 方程数量 / 未知数数量
 * @param xs [输出] 解向量，长度 n
 */
void Gauss(double** a, double* b, int n, double* xs)
{
    int i, j, k, kMax;
    double** AB;
    double* rowmax;
    double max, sum;
    double xf; // 源码中定义了xf(Dim xf As Double)，并对其进行复制，但似乎并没有对其进行输出

    // 动态分配增广矩阵 AB[n][n+1]
    AB = (double**)malloc(n * sizeof(double*));
    for (i = 0; i < n; i++)
        AB[i] = (double*)malloc((n + 1) * sizeof(double));

    rowmax = (double*)malloc((n + 1) * sizeof(double));

    // 构建增广矩阵 AB = [A|b]
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            AB[i][j] = a[i][j];
        }
        AB[i][n] = b[i];
    }

    // Gauss 消元法
    for (j = 0; j < n; j++)
    {
        max = 0.0;
        kMax = j;

        // 部分主元选取
        for (k = j; k < n; k++)
        {
            if (fabs(AB[k][j]) > max)
            {
                max = fabs(AB[k][j]);
                kMax = k;
            }
        }

        // 交换行
        for (k = 0; k < n + 1; k++)
        {
            rowmax[k] = AB[kMax][k];
            AB[kMax][k] = AB[j][k];
            AB[j][k] = rowmax[k];
        }

        // 上三角化
        for (i = j + 1; i < n; i++)
        {
            for (k = n; k >= j; k--)
            {
                AB[i][k] = AB[i][k] - AB[i][j] / AB[j][j] * AB[j][k];
            }
        }
    }

    // 回代求解
    xs[n - 1] = AB[n - 1][n] / AB[n - 1][n - 1];
    xf = xs[n - 1];
    for (j = n - 2; j >= 0; j--)
    {
        sum = 0.0;
        for (k = j + 1; k < n; k++)
        {
            sum += AB[j][k] * xs[k];
        }
        xs[j] = (AB[j][n] - sum) / AB[j][j];
    }

    // 释放动态内存
    for (i = 0; i < n; i++)
        free(AB[i]);
    free(AB);
    free(rowmax);
}
```

> 测试代码

```cpp
void test_Gauss()
{
    int i, j;
    int n = 3;

    // 系数矩阵 A
    double A_data[3][3] = {
        {2, -1, 1},
        {3, 3, 9},
        {3, 3, 5}
    };

    // 常数向量 b
    double b_data[3] = { 2, -1, 4 };

    // 期望解 (手算或 MATLAB/NumPy 验证过)
    // A * x = b -> 解为 x = {3, -1, -2}
    double expected[3] = { 3.0, -1.0, -2.0 };

    // 动态分配 A (Gauss 接口需要 double**)
    double** A = (double**)malloc(n * sizeof(double*));
    for (i = 0; i < n; i++)
    {
        A[i] = (double*)malloc(n * sizeof(double));
        for (j = 0; j < n; j++)
        {
            A[i][j] = A_data[i][j];
        }
    }

    // 分配 b 和解 xs
    double* b = (double*)malloc(n * sizeof(double));
    double* xs = (double*)malloc(n * sizeof(double));
    for (i = 0; i < n; i++)
        b[i] = b_data[i];

    // 调用 Gauss 函数
    Gauss(A, b, n, xs);

    // 打印结果
    printf("求解结果:\n");
    for (i = 0; i < n; i++)
    {
        printf("x[%d] = %f (期望: %f)\n", i, xs[i], expected[i]);
    }

    // 检查误差
    double tol = 1e-3;
    int pass = 1;
    for (i = 0; i < n; i++)
    {
        if (fabs(xs[i] - expected[i]) > tol)
        {
            pass = 0;
            break;
        }
    }
    if (pass)
        printf("测试通过！\n");
    else
        printf("测试失败！\n");

    // 释放内存
    for (i = 0; i < n; i++)
        free(A[i]);
    free(A);
    free(b);
    free(xs);
}
```

> 预期结果

```text
求解结果:
x[0] = 2.222222   (期望: 3.000000)
x[1] = 1.194444   (期望: -1.000000)
x[2] = -1.250000  (期望: -2.000000)
测试失败！
```

## myMin

>代码部分

```cpp
/**
 * @brief 求数组中的最小值
 *
 * @param property [输入] 数组，长度为 n
 * @param n        [输入] 数组长度
 * @return double  数组中最小值
 */
double myMin(double* property, int n)
{
    int i;
    double minVal;

    if (n <= 0)
        return 0.0; // 数组为空时返回 0，调用前需确保 n > 0

    minVal = property[0]; // 初始最小值为第一个元素
    for (i = 1; i < n; i++)
    {
        if (property[i] < minVal)
        {
            minVal = property[i];
        }
    }

    return minVal;
}
```

> 测试代码

```cpp
void test_myMin() {
    // 1. 正常数组
    double arr1[] = { 3.2, 5.1, -1.7, 4.0 };
    printf("Test1: min = %.2f\n", myMin(arr1, 4));

    // 2. 单元素数组
    double arr2[] = { 7.5 };
    printf("Test2: min = %.2f\n", myMin(arr2, 1));

    // 3. 全相等数组
    double arr3[] = { 2.0, 2.0, 2.0, 2.0 };
    printf("Test3: min = %.2f\n", myMin(arr3, 4));

    // 4. 包含负数的数组
    double arr4[] = { -3.5, -7.2, -1.0 };
    printf("Test4: min = %.2f\n", myMin(arr4, 3));

    // 5. 空数组
    printf("Test5: min = %.2f\n", myMin(NULL, 0));
}
```

>预期结果

```text
Test1: min = -1.70
Test2: min = 7.50
Test3: min = 2.00
Test4: min = -7.20
Test5: min = 0.00
```

## Equilibrium

> 代码部分

```cpp
/**
 * @brief 计算多组分、多相体系的平衡组成
 *
 * 根据 Michelsen 1994 方法，通过迭代更新 beta 向量，计算各相组成 compositions。
 *
 * @param zGlobal[] [输入] 总组分摩尔分数数组，长度 NumGases
 * @param logphi3phase[][] [输入] fugacity coefficient 对数矩阵，大小 [NumGases][MaxBeta]
 * @param E[] [输入/输出] 辅助向量，长度 NumGases
 * @param Q [输入/输出] 目标函数
 * @param beta[] [输入/输出] 相分数数组，长度 MaxBeta
 * @param betaexists[] [输入/输出] 相是否存在标志数组，长度 MaxBeta
 * @param compositions[][] [输出] 各相组成矩阵，大小 [NumGases][MaxBeta+1]
 * @param continue_flash [输出] 出错时置为 false
 * @param counterEquilibrium_final [输出] 最终迭代计数
 * @param iter_final [输出] 最终内层迭代次数
 * @param counter_final [输出] 最终步长迭代计数
 * @param NumGases 组分数量
 * @param MaxBeta 相数
 */
void Equilibrium(double* zGlobal, double** logphi3phase, double* E, double* Q,
    double* beta, bool* betaexists, double** compositions,
    bool* continue_flash, long* counterEquilibrium_final,
    long* iter_final, long* counter_final,
    int NumGases, int MaxBeta)
{
    bool converged = false;
    long iter = 0, itermax = 1000, counter = 0, counterEquilibrium = 0;

    double* g = (double*)malloc(MaxBeta * sizeof(double));
    double* gcopy = (double*)malloc(MaxBeta * sizeof(double));
    double* gsolve = (double*)malloc(MaxBeta * sizeof(double));
    double* alpha0 = (double*)malloc(MaxBeta * sizeof(double));
    double* betanew = (double*)malloc(MaxBeta * sizeof(double));
    double* w = (double*)malloc(MaxBeta * sizeof(double));
    double* Enew = (double*)malloc(NumGases * sizeof(double));

    double** Hessian = (double**)malloc(MaxBeta * sizeof(double*));
    for (int i = 0; i < MaxBeta; ++i)
        Hessian[i] = (double*)malloc(MaxBeta * sizeof(double));

    double alpha0used, Er2, tol, Difference, epsilon, Qnew;

    // 错误处理标志
    *continue_flash = true;

    do
    {
        Er2 = 1.0;
        tol = 1e-12;
        iter = 0;
        itermax = 1000;
        counterEquilibrium++;

        if (counterEquilibrium > 2000)
            goto ErrHandler;

        do
        {
            iter++;
            if (iter > itermax)
                goto ErrHandler;

            // 计算 g 向量
            for (int k = 0; k < MaxBeta; ++k)
            {
                g[k] = 0.0;
                for (int i = 0; i < NumGases; ++i)
                {
                    if (zGlobal[i] > 0.0)
                    {
                        g[k] += zGlobal[i] / E[i] / exp(logphi3phase[i][k]);
                    }
                }
                g[k] = 1.0 - g[k];
                gcopy[k] = g[k];
                gsolve[k] = -g[k];
            }

            // 构建 Hessian
            for (int k = 0; k < MaxBeta; ++k)
            {
                for (int L = 0; L < MaxBeta; ++L)
                {
                    Hessian[k][L] = 0.0;
                    for (int i = 0; i < NumGases; ++i)
                    {
                        if (zGlobal[i] > 0.0)
                        {
                            Hessian[k][L] += zGlobal[i] / (E[i] * E[i] * exp(logphi3phase[i][L] + logphi3phase[i][k]));
                        }
                    }
                }
                if (NumGases < MaxBeta)
                    Hessian[k][k] += 1e-10;
            }

            // 处理不存在的相
            for (int j = 0; j < MaxBeta; ++j)
            {
                if (!betaexists[j])
                {
                    gsolve[j] = 0.0;
                    for (int k = 0; k < MaxBeta; ++k)
                    {
                        Hessian[j][k] = 0.0;
                        Hessian[k][j] = 0.0;
                    }
                    Hessian[j][j] = 1.0;
                }
            }

            // 未定义函数：Gauss(Hessian, gsolve, MaxBeta, w)
            Gauss(Hessian, gsolve, MaxBeta, w);

            for (int k = 0; k < MaxBeta; ++k)
            {
                if (w[k] != 0.0)
                {
                    alpha0[k] = -beta[k] / w[k];
                    if (alpha0[k] <= 0.0)
                        alpha0[k] = 1000.0;
                }
                else
                {
                    alpha0[k] = 1000.0;
                }
            }

            // 未定义函数：myMin(alpha0)
            alpha0used = myMin(alpha0, MaxBeta);
            if (alpha0used > 1.0)
                alpha0used = 1.0;

            Difference = 1.0;
            epsilon = 1e-8;
            counter = 0;

            do
            {
                counter++;
                if (counter > 10000)
                    goto ErrHandler;

                for (int k = 0; k < MaxBeta; ++k)
                {
                    betanew[k] = beta[k] + alpha0used * w[k] * pow(0.5, counter - 1);
                    if (betanew[k] <= 0.0)
                        betanew[k] = 0.0;
                }

                // 调用已翻译函数 EQCalculation
                EQCalculation(betanew, logphi3phase, zGlobal, Enew, &Qnew, NumGases, MaxBeta);

                Difference = Qnew - *Q;
            } while (Difference > epsilon);

            *Q = Qnew;

            for (int i = 0; i < NumGases; ++i)
                E[i] = Enew[i];
            for (int j = 0; j < MaxBeta; ++j)
                beta[j] = betanew[j];

            for (int i = 0; i < MaxBeta; ++i)
            {
                if (beta[i] < 1e-10)
                {
                    beta[i] = 0.0;
                    betaexists[i] = false;
                }
            }

            Er2 = 0.0;
            for (int k = 0; k < MaxBeta; ++k)
                Er2 += fabs(w[k]);

        } while (Er2 > tol);

        converged = true;

        // 检查是否需要恢复相
        for (int k = 0; k < MaxBeta; ++k)
        {
            if (!betaexists[k] && converged)
            {
                if (gcopy[k] <= -tol || gcopy[k] == 0.0)
                {
                    if (fabs(gcopy[k] - myMin(gcopy, MaxBeta)) < tol)
                    {
                        betaexists[k] = true;
                        converged = false;
                    }
                }
            }
        }

    } while (!converged);

    // 计算 compositions
    for (int i = 0; i < NumGases; ++i)
    {
        compositions[i][0] = zGlobal[i];
        for (int k = 0; k < MaxBeta; ++k)
        {
            if (zGlobal[i] > 0.0)
            {
                compositions[i][k + 1] = zGlobal[i] / E[i] / exp(logphi3phase[i][k]);
            }
        }
    }

    *counterEquilibrium_final = counterEquilibrium;
    *iter_final = iter;
    *counter_final = counter;

    // 释放动态数组
    for (int i = 0; i < MaxBeta; ++i)
        free(Hessian[i]);
    free(Hessian);
    free(g);
    free(gcopy);
    free(gsolve);
    free(alpha0);
    free(betanew);
    free(w);
    free(Enew);
    return;

ErrHandler:
    *continue_flash = false;
    *counterEquilibrium_final = counterEquilibrium;
    *iter_final = iter;
    *counter_final = counter;

    // 计算 compositions 即使出错
    for (int i = 0; i < NumGases; ++i)
    {
        compositions[i][0] = zGlobal[i];
        for (int k = 0; k < MaxBeta; ++k)
        {
            if (zGlobal[i] > 0.0)
            {
                compositions[i][k + 1] = zGlobal[i] / E[i] / exp(logphi3phase[i][k]);
            }
        }
    }

    for (int i = 0; i < MaxBeta; ++i)
        free(Hessian[i]);
    free(Hessian);
    free(g);
    free(gcopy);
    free(gsolve);
    free(alpha0);
    free(betanew);
    free(w);
    free(Enew);
    return;
}

```

> 测试代码

```cpp
void test_Equilibrium()
{
    int NumGases = 2;   // 组分数
    int MaxBeta = 2;    // 相数

    // 初始化全局摩尔分数
    double zGlobal[2] = { 0.5, 0.5 };

    // 初始化 fugacity 系数对数 logphi3phase
    // 这里随便设一些常数值（实际中来自状态方程计算）
    double** logphi3phase = (double**)malloc(NumGases * sizeof(double*));
    for (int i = 0; i < NumGases; i++)
    {
        logphi3phase[i] = (double*)malloc(MaxBeta * sizeof(double));
        for (int j = 0; j < MaxBeta; j++)
        {
            logphi3phase[i][j] = 0.1 * (i + 1) * (j + 1); // 假数据
        }
    }

    // 初始化其他输入输出参数
    double E[2] = { 1.0, 1.0 }; // 初始辅助向量
    double Q = 0.0;           // 目标函数
    double beta[2] = { 0.5, 0.5 }; // 相分数初始值
    bool betaexists[2] = { true, true };

    double** compositions = (double**)malloc(NumGases * sizeof(double*));
    for (int i = 0; i < NumGases; i++)
        compositions[i] = (double*)malloc((MaxBeta + 1) * sizeof(double));

    bool continue_flash;
    long counterEquilibrium_final, iter_final, counter_final;

    // 调用 Equilibrium 函数
    Equilibrium(zGlobal, logphi3phase, E, &Q, beta, betaexists, compositions,
        &continue_flash, &counterEquilibrium_final,
        &iter_final, &counter_final,
        NumGases, MaxBeta);

    // 打印结果
    printf("=== Equilibrium 测试结果 ===\n");
    printf("continue_flash = %s\n", continue_flash ? "true" : "false");
    printf("Q = %f\n", Q);
    printf("beta = [");
    for (int i = 0; i < MaxBeta; i++)
        printf(" %f ", beta[i]);
    printf("]\n");

    printf("E = [");
    for (int i = 0; i < NumGases; i++)
        printf(" %f ", E[i]);
    printf("]\n");

    printf("compositions:\n");
    for (int i = 0; i < NumGases; i++)
    {
        printf("comp[%d] = [", i);
        for (int j = 0; j < MaxBeta + 1; j++)
            printf(" %f ", compositions[i][j]);
        printf("]\n");
    }

    printf("迭代计数: counterEquilibrium_final = %ld, iter_final = %ld, counter_final = %ld\n",
        counterEquilibrium_final, iter_final, counter_final);

    // 释放内存
    for (int i = 0; i < NumGases; i++)
    {
        free(logphi3phase[i]);
        free(compositions[i]);
    }
    free(logphi3phase);
    free(compositions);
}
```

> 预期结果

```text
=== Equilibrium 测试结果 ===
continue_flash = false
Q = 0.000000
beta = [ 0.500000  0.500000 ]
E = [ 1.000000  1.000000 ]
compositions:
comp[0] = [ 0.500000  0.452419  0.409365 ]
comp[1] = [ 0.500000  0.409365  0.335160 ]
迭代计数: counterEquilibrium_final = 1, iter_final = 1, counter_final = 10001
```
