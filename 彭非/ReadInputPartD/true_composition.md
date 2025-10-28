# true_composition

> code

```cpp
void true_composition(double TK, double PBar, double mol_HC, double mol_W, double aH2O, double* gNeut, double nTCO2, double nTH2S,
    int useEOS, double* reservoir_Composition,double* feed_Composition, double* total_moles, 
    double** MWgas, double** TCr, double** PCr, double** Omega, double** mf_c0, double** mf_c1, double*** kPr,
    double** composition_G, double** lnphi_Gas ,double** composition_Aq, double** lnphi_Aq, double* Compr_composition_Aq, double** mf_reservoir_Composition, 
    double** mf_feed_Composition, int* NumGases

) 
{

    /*
     * 此子程序计算总摩尔组成，用于输入
     * 多闪蒸子程序，包括油、气和水
     *
     * ======= INPUT =======
     * TK, Temperature, [Kelvin]
     * PBar, Pressure, [Bar]
     * mol_HC: Moles of hydrocarbon
     * mol_W: Moles of water
     * aH2O: Water activity
     * gNeut: Array of neutral species parameters
     * nTCO2, nTH2S: Total moles of CO2 and H2S
     * useEOS: EOS option flag
     * reservoir_Composition: The reservoir composition (all components but water)
     *
     * ======= OUTPUT =======
     * total_moles: Total number of moles
     * feed_composition: Molar composition of 15 components to make it compatible
     *                  with other routines (Multiflash)
     */

    int i;
    // double mol_O, mol_G; // 注释掉未使用的变量

    // 动态数组声明（在C中需要手动管理内存）

    //double* composition_G = NULL;
    //double* lnphi_Gas = NULL;
    //double* composition_Aq = NULL;
    //double* lnphi_Aq = NULL;
    //double Compr_composition_Aq;

    //double* mf_reservoir_Composition = NULL;
    //double* mf_feed_Composition = NULL;
    double mf_gNeut[2];


    char EOS[3] = "PR"; // EOS类型

    // 分配内存
    *MWgas = (double*)malloc(*NumGases * sizeof(double));
    *TCr = (double*)malloc(*NumGases * sizeof(double));
    *PCr = (double*)malloc(*NumGases * sizeof(double));
    *Omega = (double*)malloc(*NumGases * sizeof(double));
    *mf_c0 = (double*)malloc(*NumGases * sizeof(double));
    *mf_c1 = (double*)malloc(*NumGases * sizeof(double));

    // 分配二维数组kPr
    *kPr = (double**)malloc(*NumGases * sizeof(double*));
    for (i = 0; i < *NumGases; i++) {
        (*kPr)[i] = (double*)malloc(*NumGases * sizeof(double));
    }

    *lnphi_Gas = (double*)malloc(5 * sizeof(double));
    *composition_Aq = (double*)malloc(*NumGases * sizeof(double));
    *lnphi_Aq = (double*)malloc(*NumGases * sizeof(double));

    *mf_reservoir_Composition = (double*)malloc((*NumGases - 1) * sizeof(double));
    *mf_feed_Composition = (double*)malloc((*NumGases) * sizeof(double));


    // 初始化mf_gNeut数组
    mf_gNeut[0] = gNeut[0]; // 注意：VB数组从1开始，C从0开始
    mf_gNeut[1] = gNeut[1];

    // 获取EOS参数
    Get_EOS_Parameters(*NumGases, EOS, MWgas, TCr, PCr, Omega, mf_c0, mf_c1, kPr);

    // 初始化feed_Composition数组
    for (i = 0; i < *NumGases; i++) {
        feed_Composition[i] = 0.0;
    }

    // 计算各组分的摩尔数（不包括水）
    // 注意：VB中reservoir_Composition索引从1到NumGases-1，对应C中的0到NumGases-2
    for (i = 0; i < (*NumGases - 1); i++) {
        feed_Composition[i] = reservoir_Composition[i] * mol_HC; // in moles
    }

    // 如果useEOS = 2，则使用指定的CO2和H2S摩尔数
    if (useEOS == 2) {
        feed_Composition[1] = nTCO2; // 原feed_Composition(2)对应feed_Composition[1]
        feed_Composition[2] = nTH2S; // 原feed_Composition(3)对应feed_Composition[2]
    }

    // 设置水的摩尔数（第15个组分，对应索引14）
    feed_Composition[14] = mol_W; // mole of water in the feed

    // 计算总摩尔数
    *total_moles = 0.0;
    for (i = 0; i < *NumGases; i++) {
        // if UseEOS=2 then total_moles are different from that calc in partD
        *total_moles = *total_moles + feed_Composition[i];
    }

    // 计算摩尔分数
    for (i = 0; i < *NumGases; i++) {
        feed_Composition[i] = feed_Composition[i] / (*total_moles);
    }

}

```