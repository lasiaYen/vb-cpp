# 原子例程

## V0TP

> code

```cpp

/**
 * @brief   计算 fV0 函数值
 *
 * 对应 VB 中的 Function fV0(...)，用于根据温度 TK、压力 PBar 以及
 * 一组经验参数 (q1~q18) 计算摩尔体积的经验公式。
 *
 * @param TK    温度 (K)
 * @param PBar  压力 (bar)
 *
 * @return double 计算得到的摩尔体积
 */
double fV0(double TK, double PBar, double q1, double q2, double q3, double q4, double q5, double q6, double q7, double q8,
	double q9, double q10, double q11, double q12, double q13, double q14, double q15, double q16, double q17, double q18)
{
    double result = 0.0;

    result = q1
        + q2 / TK
        + q3 * TK
        + q4 * TK * TK
        + q5 / (TK - 227.0)
        + q6 / (647.0 - TK)
        + PBar * (q7 + q8 / TK + q9 * TK + q10 * TK * TK
            + q11 / (TK - 227.0) + q12 / (647.0 - TK))
        + PBar * PBar * (q13 + q14 / TK + q15 * TK + q16 * TK * TK
            + q17 / (TK - 227.0) + q18 / (647.0 - TK));

    return result;
}

/**
 * @brief   计算多种离子的 V0_c 和 V0_a 数组
 *
 * 对应 VB 中的 Sub V0TP()，根据给定的温度 TK 和压力 PBar，
 * 计算 Na、K、Mg、Ca、Ba、Sr、Cl、SO4 等离子的摩尔体积，
 * 并写入 V0_c[] 和 V0_a[] 数组。
 *
 * @param TK      温度 (K)
 * @param PBar    压力 (bar)
 * @param V0_c    存放阳离子结果的数组 (输出)
 * @param V0_a    存放阴离子结果的数组 (输出)
 * @return void   结果直接写入 V0_c[] 与 V0_a[]
 */
void V0TP(double TK, double PBar, double* V0_c, double* V0_a)
{
    // Na
    V0_c[iNa] = fV0(TK, PBar,
        8.76686173829976, 10.7463747460684, -3.94438220704875E-02,
        6.24254747432051E-05, -270.67565216552, 3.71906197249154,
        -0.548144793641968, 60.4681423698375, 1.6487628506002E-03,
        -1.78543820686743E-06, -7.22279554933479E-02, 6.1016223460645,
        5.37401040542647E-04, -5.78337684998825E-02, -1.64057441413787E-06,
        1.85284769422731E-09, 2.33307684437836E-04, -8.30091414725064E-03);

    // K
    V0_c[iK] = fV0(TK, PBar,
        1365.58178, -146187.60179, -4.00314, 0.004292463,
        0.0, -18894.00317, -6.66675, 715.61054, 0.019742057,
        -0.000020810979, 0.0, 81.91098, 0.00534941, -0.573121,
        -0.0000158576885, 0.0000000166987, 0.0, -0.0649312);

    // Mg
    V0_c[iMg] = fV0(TK, PBar,
        15.5470319999999, 18848.33832, -0.380047233, 0.00100500148,
        -761.133887, -22952.60934, -1.27205782, 142.833027,
        0.00343937874, -0.00000368366162, 0.0, 36.3788742,
        0.0, 0.0, 0.000000400785822, 0.0, 0.0, -0.0429183805);

    // Ca
    V0_c[iCa] = fV0(TK, PBar,
        136.567817, 3135.53072, -0.68264817, 0.0012917585,
        -761.133887, -22952.60934, -1.6506531, 162.991634,
        0.0048786976, -0.00000616046222, 0.0, 71.877495,
        0.0, 0.0, 0.000000400785822, 0.0, 0.0, -0.0429183805);

    // Ba
    V0_c[iBa] = fV0(TK, PBar,
        131.77473, 3135.53072, -0.65074076, 0.0012917585,
        -761.133887, -22952.60934, -1.6506531, 162.991634,
        0.0048786976, -0.00000616046222, 0.0, 71.877495,
        0.0, 0.0, 0.000000400785822, 0.0, 0.0, -0.0429183805);

    // Sr
    V0_c[iSr] = fV0(TK, PBar,
        131.475189, 3135.53072, -0.66249868, 0.0012917585,
        -761.133887, -22952.60934, -1.6506531, 167.104986,
        0.0048786976, -0.00000616046222, 0.0, 69.810131,
        0.0, 0.0, 0.000000400785822, 0.0, 0.0, -0.0429183805);

    // Cl
    V0_a[iCl] = fV0(TK, PBar,
        195.93822, -23046.39821, -0.29604, 0.000299867,
        0.0, -13674.59683, -0.20212, 19.60946, 0.000482443,
        -0.000000766921, 0.0, 21.30102, 0.0, 0.0,
        0.0000000714885, 0.0, 0.0, -0.00727);

    // SO4 (取后面 Dai 20151027 的参数)
    V0_a[iSO4] = fV0(TK, PBar,
        20.1255488592824, 8382.68656863382, 8.63029649907132E-02,
        3.86933893465996E-05, -262.073120048569, -20014.7787488445,
        1.93246311860811E-02, 1.16742679461389, 9.95461079547389E-06,
        7.58647008800297E-08, -1.30203576761739, -2.31700638251545,
        -3.61401771714286E-05, -4.71960303681654E-03, -7.8135942334434E-09,
        4.3129299430231E-11, 1.77316241419105E-03, 7.24281664341625E-03);
}
```

> Test

```cpp
void test_V0TP() {
	//2 3 4 5 6 7以及 2 6
	double V0_c[8] = { 0 };//下标0-7
	double V0_a[8] = { 0 };//0-7
	double TK = 273, PBar = 0.83;
	V0TP(TK, PBar, V0_c, V0_a);
	printf("预期输出，V0_c的0 1位为0，V0_a除了2 6 其它位为0\n");
	printf("V0_c:");
	for (int i = 0; i < 8; i++)
		printf("%2.f, ", V0_c[i]);
	printf("\n");
	printf("V0_a:");
	for (int i = 0; i < 8; i++)
		printf("%2.f, ", V0_a[i]);
	printf("\n");
}
```

> Result

```text
预期输出，V0_c的0 1位为0，V0_a除了2 6 其它位为0
V0_c: 0,  0, -3,  7, -22, -20, -20, -16,
V0_a: 0,  0, 16,  0,  0,  0, 18,  0,
```


## Log10

> code

```cpp
static double Log10(double x)
{
    return log10(x);
}
```


## CubicRoots

> code

```cpp
/**
 * @brief 求解三次方程的实根
 *
 * 三次方程形式：coef1 * x^3 + coef2 * x^2 + coef3 * x + coef4 = 0
 * 按照 Numerical Recipes 的方法实现，可能有一个或三个实根。
 *
 * @param coef1(double) 三次项系数（一般为1）
 * @param coef2(double) 二次项系数
 * @param coef3(double) 一次项系数
 * @param coef4(double) 常数项
 * @param root1(double*) 输出根1
 * @param root2(double*) 输出根2
 * @param root3(double*) 输出根3
 */
 void CubicRoots(double coef1, double coef2, double coef3, double coef4, double* root1, double* root2, double* root3)
{
    double QCubic, Rcubic, Xcubic, Theta;

    *root1 = 0.0;
    *root2 = 0.0;
    *root3 = 0.0;

    QCubic = (pow(coef2, 2.0) - 3.0 * coef3) / 9.0;
    Rcubic = (2.0 * pow(coef2, 3.0) - 9.0 * coef2 * coef3 + 27.0 * coef4) / 54.0;

    if (pow(Rcubic, 2.0) - pow(QCubic, 3.0) > 0.0) {
        /* 只有一个实根 */
        double A = pow(fabs(Rcubic) + sqrt(pow(Rcubic, 2.0) - pow(QCubic, 3.0)), 1.0 / 3.0);
        double B = QCubic / A;
        double sigh = (Rcubic > 0) ? 1.0 : ((Rcubic < 0) ? -1.0 : 0.0);
        *root1 = -sigh * (A + B) - coef2 / 3.0;
    }
    else {
        /* 三个实根 */
        Xcubic = Rcubic / sqrt(pow(QCubic, 3.0));
        /* 使用反三角恒等式 arccos(x) = atan(-x/sqrt(-x^2+1)) + pi/2 */
        Theta = atan(-Xcubic / sqrt(1.0 - Xcubic * Xcubic)) + pi / 2.0;
        *root1 = -2.0 * sqrt(QCubic) * cos(Theta / 3.0) - coef2 / 3.0;
        *root2 = -2.0 * sqrt(QCubic) * cos((Theta + 2.0 * pi) / 3.0) - coef2 / 3.0;
        *root3 = -2.0 * sqrt(QCubic) * cos((Theta + 4.0 * pi) / 3.0) - coef2 / 3.0;
    }
}
```

> Test

```cpp
void Test_CubicRoots() 
{
    double coef1 = 1,  coef2 = -6,  coef3 = 9,  coef4 = -2,  root1=-100,  root2=-100,  root3=-100;
    //这个三次函数的极大值是f(1)=3,极小值是f(3)=-2,必有3个实数根
    CubicRoots(coef1, coef2, coef3, coef4, &root1, &root2, &root3);
    printf("root1:%.2f, root2:%.2f, root3:%.2f\n", root1, root2, root3);
}
```

> Result

```text
root1:0.27, root2:3.73, root3:2.00
```


## QuadraticRoots

> code

```cpp
/**
 * @brief 求解二次方程的实根
 *
 * 二次方程形式：coef1 * x^2 + coef2 * x + coef3 = 0
 * 使用 Numerical Recipes 推荐的高精度算法
 *
 * @param coef1(double) 二次项系数
 * @param coef2(double) 一次项系数
 * @param coef3(double) 常数项
 * @param root1(double*) 输出根1
 * @param root2(double*) 输出根2
 */
 void QuadraticRoots(double coef1, double coef2, double coef3, double* root1, double* root2)
{
    double qroot;
    double discriminant = pow(coef2, 2.0) - 4.0 * coef1 * coef3;
    if (discriminant < 0.0) {
        /* 无实根 */
        *root1 = NAN;
        *root2 = NAN;
        return;
    }
    qroot = -0.5 * (coef2 - sqrt(discriminant));
    *root1 = coef3 / qroot;
    *root2 = qroot / coef1;
}
```

> Test
```cpp
void Test_QuadraticRoots()
{
    //对应上面例子， 求极值点，此处root1和root2结果必为1 3
    double coef1 = 3, coef2 = -12, coef3 = 9, root1 = -100, root2 = -100;
    QuadraticRoots(coef1, coef2, coef3, &root1, &root2);
    printf("root1:%.2f, root2:%.2f", root1, root2);
}
```

> Result

```text
root1:1.00, root2:3.00
```


## OrderPhases

> code

```cpp
/**
 * @brief  按密度排序相和相属性， 保证只有水相类型可以输出到第三列
 *
 * 因为会对参数列表修改，因此采用指针传递
 *
 * @param   compositions    各相中各组分的组成
 * @param   phi
 * @param   beta  表示系统中每个相所占的比例
 * @param   MW_Phase  各相的分子量
 * @param   mass_phase    各相的质量
 * @param   density   各相的密度
 * @param   Compr 各相的压缩因子
 * @param   nComp
 * @param   nPhase
 */
 void OrderPhases(double*** compositions, double*** phi, double** beta, double** MW_Phase,
    double** mass_phase, double** density, double** Compr, int* nComp, int* nPhase)
{
    int nC = *nComp;
    int nP = *nPhase;
    int i, k, L;

    double** compositions0 = (double**)malloc(nC * sizeof(double*));
    double** phi0 = (double**)malloc(nC * sizeof(double*));
    for (i = 0; i < nC; i++)
    {
        compositions0[i] = (double*)malloc((nP + 1) * sizeof(double));
        phi0[i] = (double*)malloc(nP * sizeof(double));
    }
    double* BETA0 = (double*)malloc(nP * sizeof(double));
    double* MW_Phase0 = (double*)malloc(nP * sizeof(double));
    double* mass_phase0 = (double*)malloc(nP * sizeof(double));
    double* density0 = (double*)malloc(nP * sizeof(double));
    double* Compr0 = (double*)malloc(nP * sizeof(double));

    //按密度重新排序相
    int ascending = 1;
    int* idx_dens = (int*)malloc(nP * sizeof(int));
    for (k = 0; k < nP; k++)
    {
        density0[k] = (*density)[k];
        idx_dens[k] = k;
    }
    BubbleSort(ascending, density0, idx_dens, nP);

    //将所有变量分配给临时数组，按密度排序
    for (k = 0; k < nP; k++)
    {
        L = idx_dens[k];
        for (i = 0; i < nC; i++)
        {
            /*
            原代码：
               compositions0(i, k + 1) = compositions(i, L + 1)
               phi0(i, k) = phi(i, L)
            */
            compositions0[i][k + 1] = (*compositions)[i][L];
            phi0[i][k] = (*phi)[i][L - 1];
        }
        BETA0[k] = (*beta)[L - 1];
        MW_Phase0[k] = (*MW_Phase)[L - 1];
        mass_phase0[k] = (*mass_phase)[L - 1];
        density0[k] = (*density)[L - 1];
        Compr0[k] = (*Compr)[L - 1];
    }

    //确定平衡状态下存在的相类型
    //这里将其全部初始化为-1(在vb中是0，即数组最小索引的前一位)
    int gas = -1, oil = -1, aqueous = -1, zeroPhase = -1, HC1 = -1, HC2 = -1, HCPhases = 0;
    for (k = 0; k < nP; k++) {
        if (density0[k] >= 0.3 && compositions0[nC - 1][k + 1] > 0.5)
            aqueous = k;
        else if (density0[k] == 0)  // 判断double类型是否为0一般不这样写，此处对应原Vb代码。
            zeroPhase = k;
        else {
            if (HC1 == 0) HC1 = k;
            else HC2 = k;
        }
    }
    //特别的，HC1和HC2的赋值与之前的循环量有关，而循环量已严格匹配，因此不需要偏移
    if (HC2 > 0) HCPhases = 2; else HCPhases = 1;
    if (HCPhases == 1)
    {
        if (density0[HC1] > 1e-20 && density0[HC1] < 0.3) gas = HC1; // 蒸汽状 vapor-like
        else if (density0[HC1] > 1e-20) oil = HC1;  // 富含水分的液体状密度 any non-vapor, non-aqueous phase
    }
    if (HCPhases == 2)
    {
        if (density0[HC1] < density0[HC2])
        {
            gas = HC1;
            oil = HC2;
        }
        else {
            oil = HC1;
            gas = HC2;
        }
    }

    //如果存在水相，则输出三相
    if (aqueous != -1 && nP != 3) {
        /*
        vb原码：nPhase = 3 ,redim compositions(nComp, nPhase + 1)等直接修改了参数列表。
        推测: 使用引用传递
        */
        nP = 3;
        *nPhase = 3;

        for (i = 0; i < nC; i++) {
            free((*compositions)[i]);
            free((*phi)[i]);
        }
        free(*compositions);
        free(*phi);
        free(*beta);
        free(*MW_Phase);
        free(*mass_phase);
        free(*density);
        free(*Compr);

        *compositions = (double**)malloc(nC * sizeof(double*));
        *phi = (double**)malloc(nC * sizeof(double*));
        for (i = 0; i < nC; i++) {
            (*compositions)[i] = (double*)malloc((nP + 1) * sizeof(double));
            (*phi)[i] = (double*)malloc(nP * sizeof(double));
        }
        *beta = (double*)malloc(nP * sizeof(double));
        *MW_Phase = (double*)malloc(nP * sizeof(double));
        *mass_phase = (double*)malloc(nP * sizeof(double));
        *density = (double*)malloc(nP * sizeof(double));
        *Compr = (double*)malloc(nP * sizeof(double));
    }

    // output phases
    L = 0;
    for (k = 0; k < nP; k++)
    {
        if (k == gas) L = 1;//gas由HC赋值，HC由k赋值，因此直接照搬原代码
        else if (k == oil) L = 2;
        else if (k == aqueous) L = 3;

        if (L != 0)
        {
            if (BETA0[k] > 0)
            {
                //注：原代码中是beta[L] = BETA0[k]，那么存在beta的第1、2、3个元素存在被修改的情况，C中的第1、2、3个元素，对应的索引为0、1、2，此处注意偏移量。
                (*beta)[L - 1] = BETA0[k];
                (*MW_Phase)[L - 1] = MW_Phase0[k];
                (*mass_phase)[L - 1] = mass_phase0[k];
                (*density)[L - 1] = density0[k];
                (*Compr)[L - 1] = Compr0[k];
                for (i = 0; i < nC; i++)
                {
                    (*compositions)[i][L] = compositions0[i][k + 1];
                    (*phi)[i][L - 1] = phi0[i][k];
                }
            }
            else {
                (*beta)[L - 1] = 0;
                (*MW_Phase)[L - 1] = 0;
                (*mass_phase)[L - 1] = 0;
                (*density)[L - 1] = 0;
                (*Compr)[L - 1] = 0;
                for (i = 0; i < nC; i++)
                {
                    (*compositions)[i][L] = 0;
                    (*phi)[i][L - 1] = 0;
                }
            }
        }
    }

    // 从输出中删除任何不存在的阶段。 注释：本质上第一个修改beat[1],对应C语言的就是beta[0].，可以省略L，直接写数字即可，未修改是为了对应vb代码。
    L = 1;
    if (gas == 0)
    {
        (*beta)[L - 1] = 0;
        (*MW_Phase)[L - 1] = 0;
        (*mass_phase)[L - 1] = 0;
        (*density)[L - 1] = 0;
        (*Compr)[L - 1] = 0;
        for (i = 0; i < nC; i++)
        {
            (*compositions)[i][L] = 0;
            (*phi)[i][L - 1] = 0;
        }
    }
    L = 2;
    if (oil == 0)
    {
        (*beta)[L - 1] = 0;
        (*MW_Phase)[L - 1] = 0;
        (*mass_phase)[L - 1] = 0;
        (*density)[L - 1] = 0;
        (*Compr)[L - 1] = 0;
        for (i = 0; i < nC; i++)
        {
            (*compositions)[i][L] = 0;
            (*phi)[i][L - 1] = 0;
        }
    }
    L = 3;
    if (aqueous == 0 && nP == 3)
    {
        (*beta)[L - 1] = 0;
        (*MW_Phase)[L - 1] = 0;
        (*mass_phase)[L - 1] = 0;
        (*density)[L - 1] = 0;
        (*Compr)[L - 1] = 0;
        for (i = 0; i < nC; i++)
        {
            (*compositions)[i][L] = 0;
            (*phi)[i][L - 1] = 0;
        }
    }

    // 释放临时开辟的空间
    for (i = 0; i < nC; i++)
    {
        free(compositions0[i]);
        free(phi0[i]);
    }
    free(compositions0);
    free(phi0);
    free(BETA0);
    free(MW_Phase0);
    free(mass_phase0);
    free(density0);
    free(Compr0);
    free(idx_dens);
}
```

> Test
```
已确认能运行
```


