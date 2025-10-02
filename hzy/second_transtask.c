#include <stdio.h>
#include <math.h>
#include "publicpara.h"

int i;            // 当前样本索引
double TK;        // 温度 [K]，用于密度/溶液计算
double TC;        // 温度 [°C]，用于 CalcRhoTP
double PBar;      // 压力 [bar]，用于密度计算
double Patm;      // 大气压 [bar]，用于密度计算


// 假设的离子总数量
const int NUM_SAMPLES = 100;
const int NUM_CAT = 100;
const int NUM_AN = 100;
const int NUM_NEUT = 100;
const int MaxGases=100;
const int MaxComponents=100;

// 初始混合物浓度
double AlkMix[NUM_SAMPLES];
double TAcMix[NUM_SAMPLES];
double TH2SaqMix[NUM_SAMPLES];
double TH4SiO4Mix[NUM_SAMPLES];
double TH3BO3Mix[NUM_SAMPLES];
double TNH4Mix[NUM_SAMPLES];
double NaMix[NUM_SAMPLES];
double KMix[NUM_SAMPLES];
double MgMix[NUM_SAMPLES];
double CaMix[NUM_SAMPLES];
double SrMix[NUM_SAMPLES];
double BaMix[NUM_SAMPLES];
double FeMix[NUM_SAMPLES];
double FMix[NUM_SAMPLES];
double ZnMix[NUM_SAMPLES];
double PbMix[NUM_SAMPLES];
double ClMix[NUM_SAMPLES];
double SO4Mix[NUM_SAMPLES];
double BrMix[NUM_SAMPLES];

// 暂存浓度
double HCO3stpMix[NUM_SAMPLES];
double ACstpMix[NUM_SAMPLES];
double HstpMix[NUM_SAMPLES];
double OHstpMix[NUM_SAMPLES];
double CO3stpMix[NUM_SAMPLES];
double HSstpMix[NUM_SAMPLES];
double NH4STPMix[NUM_SAMPLES];
double H2BO3stpMix[NUM_SAMPLES];

// pH 测量值
double pHMeterStpMix[NUM_SAMPLES];

double TDSMix[NUM_SAMPLES];  // 每个样本的总溶解固体浓度 [mg/L]
double MWCat[NUM_CAT];   // 阳离子分子量 [g/mol]
double MWAn[NUM_AN];     // 阴离子分子量 [g/mol]
double MWNeut[NUM_NEUT]; // 中性组分分子量 [g/mol]


double rho25c;          // 溶液参考密度 [kg/L]
int UseMolal;           // 是否使用摩尔浓度模式 (0=否, 1=是)
int useEOSmix[NUM_SAMPLES]; // EOS 使用标志
int kk;                 // 当前混合物索引
double IStCosolvent;    // 溶液离子强度
double Ist;             // 当前溶液离子强度
double xMeOH, xMEG;     // 共溶剂摩尔分数
double DpHj;            // pH 校正值

void D2_CalcDensitypH(int i); 
void CalcIonicStrength(void);
double fTPFunc(int iTP); 
double CalcRhoTP(double TK, double TC, double PBar, double Patm);
void C1_ThermodynamicEquilConsts(void);
void C2_PitzerActCoefs_T_P_ISt(int gNeut, double aH2O, double TK, double TC, double PBar, double Patm);
void PengRobinson3(void);
void C5_CalcpHPCO2PH2SSTP(void);

double mc[NUM_CAT];
double ma[NUM_AN];
double mn[NUM_NEUT];

// 不知道是不是离子
int H4SiO4, TH3BO3;

// C2_PitzerActCoefs_T_P_ISt函数参数
double gNeut[NUM_NEUT]; // 中性组分活度系数或状态
double aH2O;            // 水活度，0~1
double rho_Mix[NUM_SAMPLES]; // 每个样品的混合物密度

// mc 阳离子浓度
// ma 阴离子浓度
// mn 中性物质浓度
void D1_CalcDensity(double* mc, double* ma, double* mn, int i)
{
    double Alk, TAc, TCO2, TNH4, TH3BO3, TH2Saq, TH4SiO4, TFe;
    double TCa;

    //************************************
    // 初始化混合物浓度
    HCO3stpMix[i] = AlkMix[i];
    ACstpMix[i] = TAcMix[i];
    HstpMix[i] = 0.000001;
    OHstpMix[i] = 0.0000001;
    CO3stpMix[i] = 0.0;
    HSstpMix[i] = 0.0;
    NH4STPMix[i] = TNH4Mix[i];
    H2BO3stpMix[i] = 0.0;

    int Iteration2 = 0;

    mc[iH] = HstpMix[i];
    mc[iNa] = NaMix[i];
    mc[iK]  = KMix[i];
    mc[iMg] = MgMix[i];
    mc[iCa] = CaMix[i];
    double TCa = mc[iCa];
    mc[iSr] = SrMix[i];
    mc[iBa] = BaMix[i];
    mc[iFe] = FeMix[i];
    mc[iZn] = ZnMix[i];
    mc[iPb] = PbMix[i];

    ma[iOH] = OHstpMix[i];
    ma[iCl] = ClMix[i];
    ma[iAc] = ACstpMix[i];
    mc[iNH4] = NH4STPMix[i];
    ma[iH2BO3] = H2BO3stpMix[i];
    ma[iHCO3] = HCO3stpMix[i];
    ma[iCO3] = CO3stpMix[i];
    ma[iH3SiO4] = 0.0;
    ma[iH2SiO4] = 0.0;
    ma[iSO4] = SO4Mix[i];
    ma[iHS] = HSstpMix[i];
    ma[intF] = FMix[i];
    ma[iBr] = BrMix[i];

    double Alk = AlkMix[i];
    double TAc = TAcMix[i];
    double TH2Saq = TH2SaqMix[i];
    double TH4SiO4 = TH4SiO4Mix[i];
    double TH3BO3 = TH3BO3Mix[i];
    double TNH4 = TNH4Mix[i];

    mn[iNH3] = 0.0;
    mn[iH3BO3] = TH3BO3;
    mn[iH4SiO4aq] = TH4SiO4;
    double TFe = mc[iFe];

    double TDSOld = TDSMix[i];
    double rhoOld = rho_Mix[i];
    double TDS = TDSMix[i];
    double TDSSSE = 1.0;

    if (UseMolal == 0) {
        while (TDSSSE > 1e-8 && Iteration2 < 20) {
            D2_CalcDensitypH(i); // 计算离子强度、密度及HCO3, AC, HS speciation

            double CalculateTDSDen = 0.0;
            TDS = 0.0;

            for (int iden = 1; iden < NUM_CAT; iden++) {
                TDS += 1000.0 * rho25c * mc[iden] * MWCat[iden];
                CalculateTDSDen += 0.001 * mc[iden] * MWCat[iden];
            }

            for (int iden = 1; iden < NUM_AN; iden++) {
                TDS += 1000.0 * rho25c * ma[iden] * MWAn[iden];
                CalculateTDSDen += 0.001 * ma[iden] * MWAn[iden];
            }

            for (int iden = 1; iden < NUM_NEUT; iden++) {
                TDS += 1000.0 * rho25c * mn[iden] * MWNeut[iden];
                CalculateTDSDen += 0.001 * mn[iden] * MWNeut[iden];
            }

            TDS = TDS / (1.0 + CalculateTDSDen);

            // 根据新TDS更新摩尔浓度
            for (int iden = 1; iden < NUM_CAT; iden++)
                mc[iden] *= (rho25c - TDSOld / 1e6) / (rho25c - TDS / 1e6);

            for (int iden = 1; iden < NUM_AN; iden++)
                ma[iden] *= (rho25c - TDSOld / 1e6) / (rho25c - TDS / 1e6);

            for (int iden = 1; iden < NUM_NEUT; iden++)
                mn[iden] *= (rho25c - TDSOld / 1e6) / (rho25c - TDS / 1e6);

            Alk *= (rho25c - TDSOld / 1e6) / (rho25c - TDS / 1e6);
            TAc *= (rho25c - TDSOld / 1e6) / (rho25c - TDS / 1e6);
            TCO2 *= (rho25c - TDSOld / 1e6) / (rho25c - TDS / 1e6);
            TNH4 *= (rho25c - TDSOld / 1e6) / (rho25c - TDS / 1e6);
            TH3BO3 *= (rho25c - TDSOld / 1e6) / (rho25c - TDS / 1e6);
            TH2Saq *= (rho25c - TDSOld / 1e6) / (rho25c - TDS / 1e6);
            TH4SiO4 *= (rho25c - TDSOld / 1e6) / (rho25c - TDS / 1e6);
            TFe *= (rho25c - TDSOld / 1e6) / (rho25c - TDS / 1e6);

            D2_CalcDensitypH(i); // 再次计算密度和分配

            if (TDSOld == 0.0)
                break;

            TDSSSE = pow(TDS / TDSOld - 1.0, 2.0);
            TDSOld = TDS;
            Iteration2++;
        }
    } else {
        // UseMolal != 0 时的处理
        CalcIonicStrength();
        xMeOH = 0.0;
        xMEG = 0.0;
        IStCosolvent = Ist;
        // mt = fTPFunc(0);
        rho25c = CalcRhoTP(TK, TC, PBar, Patm);
        // pH = pHMeterStpMix[i] + DpHj;

        C1_ThermodynamicEquilConsts();
        C2_PitzerActCoefs_T_P_ISt(gNeut, aH2O, TK, TC, PBar, Patm);
        PengRobinson3();
        C5_CalcpHPCO2PH2SSTP();
        
        // 源代码中前面都没有i
        mc[iH] = iH;
        ma[iOH] = iOH;
        ma[iAc] = iAc;
        mn[iNH3] = iNH3;
        ma[iH2BO3] = iH2BO3;
        ma[iHCO3] = iHCO3;
        ma[iCO3] = iCO3;
        ma[iHS] = iHS;
        ma[iH3SiO4] = iH3SiO4;
        ma[iH2SiO4] = iH2SiO4;
        mn[iH4SiO4aq] = H4SiO4;
        mn[iH3BO3] = TH3BO3 - iH2BO3;
        mn[iCO2aq] = iCO2aq;
        mn[iH2Saq] = iH2Saq;
        mn[iHAcaq] = iHAcaq;

        if (useEOSmix[kk] == 1) {
            mc[iH] = 1e-7;
            ma[iOH] = 1e-7;
            ma[iAc] = TAc;
            ma[iHCO3] = Alk;
            ma[iHS] = 0.0;
            mn[iH4SiO4aq] = TH4SiO4;
            mn[iH3BO3] = TH3BO3;
            mc[iNH4] = TNH4;
        }

        // mt = fTPFunc(0);
        CalcIonicStrength();
        rho25c = CalcRhoTP(TK, TC, PBar, Patm);
    }
}

// 全局变量声明（假设这些是全局的，根据上下文）
double Ppsia;  // 压力（psia）
double Patm;   // 大气压力（atm）
double PBar;   // 压力（bar）
double TF;     // 华氏温度
double TC;     // 摄氏温度
double TK;     // 绝对温度（K）

// 其他输入变量（假设为全局，根据上下文）
double Pvol;   // 体积压力（psia）
double TVol;   // 体积温度（F）
double PpH;    // pH压力（psia）
double TpH;    // pH温度（F）

/**
 * @brief 设置温度和压力参数（fTPFunc）
 *
 * 此函数根据iTP值设置温度和压力参数。
 * - iTP=0: 标准条件 (77F, 14.696 psia)
 * - iTP=1: 使用Pvol和TVol
 * - iTP=2: 使用PpH和TpH
 * - iTP=3: API标准条件 (60F, 14.696 psia)
 *
 * @param iTP 温度压力选项
 * @return TK (绝对温度，K) 或 -1（错误）
 */
double fTPFunc(int iTP) {
    switch (iTP) {
        case 0: // 标准条件
            Ppsia = 14.696;
            Patm = Ppsia / 14.696;
            PBar = Ppsia / 14.503774;
            TF = 77;
            TC = (TF - 32) * 5.0 / 9.0;
            TK = TC + 273.15;
            break;
        case 1: // 使用 Pvol 和 TVol
            Ppsia = Pvol;
            Patm = Ppsia / 14.696;
            PBar = Ppsia / 14.503774;
            TF = TVol;
            TC = (TF - 32) * 5.0 / 9.0;
            TK = TC + 273.15;
            break;
        case 2: // 使用 PpH 和 TpH
            Ppsia = PpH;
            Patm = Ppsia / 14.696;
            PBar = Ppsia / 14.503774;
            TF = TpH;
            TC = (TF - 32) * 5.0 / 9.0;
            TK = TC + 273.15;
            break;
        case 3: // API 标准条件
            Ppsia = 14.696;
            Patm = Ppsia / 14.696;
            PBar = Ppsia / 14.503774;
            TF = 60;
            TC = (TF - 32) * 5.0 / 9.0;
            TK = TC + 273.15;
            break;
        default:
            // 处理无效输入
            fprintf(stderr, "Invalid iTP value: %d\n", iTP);
            return -1.0; // 返回错误码
    }
    return TK; // 返回绝对温度（K）
}

/**
 * @brief 读取输入部分C（ReadInputPartC）
 *
 * 此函数读取并设置混合物kk的输入参数，包括标准条件下的浓度、pH选项、气体组成等。
 * 计算密度和TDS，更新各种浓度数组，并处理摩尔浓度到TDS的转换。
 * 适用于电解质溶液混合物的输入处理，支持多种运行模式（如测试案例、海水混合等）。
 *
 * @param kk 当前混合物索引
 * @param mt [输出] 温度压力函数值
 * @param UseTPpHMix [数组] 使用TpH混合标志 [MaxMix]
 * @param Run10TestCases 运行10测试案例标志
 * @param Loop10 循环10计数
 * @param Run_Seawater_Mixing 海水混合运行标志
 * @param LoopMixing 混合循环计数
 * @param Run_MixingTwoWells 运行两个井混合标志
 * @param RunMultiMix 多混合运行标志
 * @param LoopResChem 残余化学循环计数
 * @param RunStatMix 统计混合运行标志
 * @param HCO3stpMix [数组] 标准条件下HCO3-浓度 [MaxMix]
 * @param AlkMix [数组] 碱度混合 [MaxMix]
 * @param ACstpMix [数组] 标准条件下Ac-浓度 [MaxMix]
 * @param TAcMix [数组] 总Ac混合 [MaxMix]
 * @param HstpMix [数组] 标准条件下H+浓度 [MaxMix]
 * @param OHstpMix [数组] 标准条件下OH-浓度 [MaxMix]
 * @param CO3stpMix [数组] 标准条件下CO3^2-浓度 [MaxMix]
 * @param HSstpMix [数组] 标准条件下HS-浓度 [MaxMix]
 * @param NH4STPMix [数组] 标准条件下NH4+浓度 [MaxMix]
 * @param TNH4Mix [数组] 总NH4混合 [MaxMix]
 * @param H2BO3stpMix [数组] 标准条件下H2BO3-浓度 [MaxMix]
 * @param TDS [输出] 总溶解固体（mg/L）
 * @param yH2S [输入/输出] H2S气体摩尔分数
 * @param yCO2 [输入/输出] CO2气体摩尔分数
 * @param Iteration2 迭代2计数
 * @param mc 阳离子摩尔浓度数组 [MaxSpecies]
 * @param ma 阴离子摩尔浓度数组 [MaxSpecies]
 * @param mn 中性物种摩尔浓度数组 [MaxSpecies]
 * @param iH H+ 索引
 * @param iNa Na+ 索引
 * @param iK K+ 索引
 * @param iMg Mg^2+ 索引
 * @param iCa Ca^2+ 索引
 * @param iSr Sr^2+ 索引
 * @param iBa Ba^2+ 索引
 * @param iFe Fe^2+ 索引
 * @param iZn Zn^2+ 索引
 * @param iPb Pb^2+ 索引
 * @param iRa Ra^2+ 索引
 * @param iOH OH- 索引
 * @param iCl Cl- 索引
 * @param iAc Ac- 索引
 * @param iNH4 NH4+ 索引
 * @param iH2BO3 H2BO3- 索引
 * @param iHCO3 HCO3- 索引
 * @param iCO3 CO3^2- 索引
 * @param iH3SiO4 H3SiO4- 索引
 * @param iH2SiO4 H2SiO4^2- 索引
 * @param iSO4 SO4^2- 索引
 * @param iHS HS- 索引
 * @param intF F- 索引
 * @param iBr Br- 索引
 * @param Alk [输入/输出] 碱度
 * @param TAc [输入/输出] 总醋酸浓度
 * @param TH2Saq [输入/输出] 总H2S(aq)浓度
 * @param TH4SiO4 [输入/输出] 总H4SiO4浓度
 * @param TH3BO3 [输入/输出] 总H3BO3浓度
 * @param TNH4 [输入/输出] 总NH4浓度
 * @param TFe [输入/输出] 总Fe浓度
 * @param TPb [输入/输出] 总Pb浓度
 * @param TZn [输入/输出] 总Zn浓度
 * @param iNH3 NH3 索引
 * @param iH3BO3 H3BO3 索引
 * @param iH4SiO4aq H4SiO4(aq) 索引
 * @param use_pH [输入/输出] pH使用选项
 * @param usepHmix [数组] pH混合使用 [MaxMix]
 * @param UseH2Sgas [输入/输出] 使用H2S气体标志
 * @param UseH2SgasMix [数组] H2S气体混合使用 [MaxMix]
 * @param TCO2 [输入/输出] 总CO2浓度
 * @param TCO2Mix [数组] 总CO2混合 [MaxMix]
 * @param yCO2Mix [数组] CO2气体摩尔分数混合 [MaxMix]
 * @param yH2SMix [数组] H2S气体摩尔分数混合 [MaxMix]
 * @param useEOSmix [数组] 使用EOS混合 [MaxMix]
 * @param SumofZMix [数组] Z混合总和 [MaxMix]
 * @param zMix [数组] Z混合 [MaxMix][MaxComponents]
 * @param CalculatedTDSMix [数组] 计算TDS混合 [MaxMix]
 * @param NaMix [数组] Na混合 [MaxMix]
 * @param KMix [数组] K混合 [MaxMix]
 * @param MgMix [数组] Mg混合 [MaxMix]
 * @param CaMix [数组] Ca混合 [MaxMix]
 * @param TCa [输出] 总Ca
 * @param SrMix [数组] Sr混合 [MaxMix]
 * @param BaMix [数组] Ba混合 [MaxMix]
 * @param FeMix [数组] Fe混合 [MaxMix]
 * @param ZnMix [数组] Zn混合 [MaxMix]
 * @param PbMix [数组] Pb混合 [MaxMix]
 * @param RaMix [数组] Ra混合 [MaxMix]
 * @param ClMix [数组] Cl混合 [MaxMix]
 * @param SO4Mix [数组] SO4混合 [MaxMix]
 * @param FMix [数组] F混合 [MaxMix]
 * @param BrMix [数组] Br混合 [MaxMix]
 * @param rho25CMix [数组] 25°C密度混合 [MaxMix]
 * @param H3SiO4Mix [数组] H3SiO4混合 [MaxMix]
 * @param H2SiO4Mix [数组] H2SiO4混合 [MaxMix]
 * @param NH3Mix [数组] NH3混合 [MaxMix]
 * @param H4SiO4Mix [数组] H4SiO4混合 [MaxMix]
 * @param H3BO3Mix [数组] H3BO3混合 [MaxMix]
 * @param CO2aqMix [数组] CO2(aq)混合 [MaxMix]
 * @param H2SaqMix [数组] H2S(aq)混合 [MaxMix]
 * @param HACaqMix [数组] HAc(aq)混合 [MaxMix]
 * @param AlkMix [数组] 碱度混合 [MaxMix]
 * @param TAcMix [数组] 总Ac混合 [MaxMix]
 * @param TH2SaqMix [数组] 总H2S(aq)混合 [MaxMix]
 * @param TH4SiO4Mix [数组] 总H4SiO4混合 [MaxMix]
 * @param TNH4Mix [数组] 总NH4混合 [MaxMix]
 * @param TH3BO3Mix [数组] 总H3BO3混合 [MaxMix]
 * @param yCH4Mix [数组] CH4气体摩尔分数混合 [MaxMix]
 * @param UseTPVolMix [数组] 使用TpV混合 [MaxMix]
 * @param WaterDensityMix [数组] 水密度混合 [MaxMix]
 * @param rho25c [输入] 25°C密度
 * @param UseMolal 使用摩尔浓度标志
 * @param CalculateTDSDen [输出] 计算TDS密度
 * @param NumCat 阳离子数量
 * @param NumAn 阴离子数量
 * @param NumNeut 中性物种数量
 * @param MWCat 阳离子分子量数组 [MaxCat]
 * @param MWAn 阴离子分子量数组 [MaxAnion]
 * @param MWNeut 中性物种分子量数组 [MaxNeut]
 * @param iCO2aq CO2(aq) 索引
 * @param iH2Saq H2S(aq) 索引
 * @param iHAcaq HAc(aq) 索引
 */
void ReadInputPartC(int kk, double *mt, int UseTPpHMix[], int Run10TestCases, int Loop10, int Run_Seawater_Mixing, int LoopMixing,
                    int Run_MixingTwoWells, int RunMultiMix, int LoopResChem, int RunStatMix,
                    double HCO3stpMix[], double AlkMix[], double ACstpMix[], double TAcMix[], double HstpMix[], double OHstpMix[],
                    double CO3stpMix[], double HSstpMix[], double NH4STPMix[], double TNH4Mix[], double H2BO3stpMix[], double *TDS,
                    double *yH2S, double *yCO2, int *Iteration2,
                    double mc[], double ma[], double mn[],
                    /* int iH, int iNa, int iK, int iMg, int iCa, int iSr, int iBa, int iFe,
                    int iZn, int iPb, int iRa, int iOH, int iCl, int iAc, int iNH4, int iH2BO3, int iHCO3, int iCO3, int iH3SiO4,
                    int iH2SiO4, int iSO4, int iHS, int intF, int iBr,
                    */
                    double *Alk, double *TAc, double *TH2Saq, double *TH4SiO4, double *TH3BO3, double *TNH4, double *TFe, double *TPb,
                    double *TZn, 
                    /*int iNH3, int iH3BO3, int iH4SiO4aq,*/
                    int *use_pH, int usepHmix[], int *UseH2Sgas, int UseH2SgasMix[],
                    double *TCO2, double TCO2Mix[], double yCO2Mix[], double yH2SMix[], int useEOSmix[], double SumofZMix[],
                    double zMix[][MaxComponents], double CalculatedTDSMix[], double NaMix[], double KMix[], double MgMix[], double CaMix[],
                    double *TCa, double SrMix[], double BaMix[], double FeMix[], double ZnMix[], double PbMix[], double RaMix[],
                    double ClMix[], double SO4Mix[], double FMix[], double BrMix[], double rho25CMix[], double H3SiO4Mix[], double H2SiO4Mix[],
                    double NH3Mix[], double H4SiO4Mix[], double H3BO3Mix[], double CO2aqMix[], double H2SaqMix[], double HACaqMix[],
                    double AlkMix2[], double TAcMix2[], double TH2SaqMix[], double TH4SiO4Mix2[], double TNH4Mix2[], double TH3BO3Mix2[],
                    double yCH4Mix[], int UseTPVolMix[], double WaterDensityMix[], double rho25c, int UseMolal, double *CalculateTDSDen,
                    int NumCat, int NumAn, int NumNeut, double MWCat[], double MWAn[], double MWNeut[]
                    /*int iCO2aq, int iH2Saq, int iHAcaq*/
                ) {
    *mt = fTPFunc(0);  // Densitym TDS, and m calculated at STP condition
    if (UseTPpHMix[kk] == 1) *mt = fTPFunc(2);

    if (Run10TestCases == 1 && Loop10 > 1) goto label100;
    if (Run_Seawater_Mixing == 1 && LoopMixing > 1) goto label100;
    if (Run_MixingTwoWells == 1 && LoopMixing > 1) goto label100;
    if (RunMultiMix == 1 && LoopResChem > 1) goto label100;
    if (RunStatMix == 1 && LoopMixing > 1) goto label100;

    HCO3stpMix[kk] = AlkMix[kk];
    ACstpMix[kk] = TAcMix[kk];
    HstpMix[kk] = 0.000001;
    OHstpMix[kk] = 0.0000001;
    CO3stpMix[kk] = 0;
    HSstpMix[kk] = 0;
    NH4STPMix[kk] = TNH4Mix[kk];
    H2BO3stpMix[kk] = 0;
    *TDS = 0;
    *yH2S = 0;
    *yCO2 = 0;

    *Iteration2 = 0;

    mc[iH] = HstpMix[kk];
    mc[iNa] = NaMix[kk];
    mc[iK] = KMix[kk];
    mc[iMg] = MgMix[kk];
    mc[iCa] = CaMix[kk];
    *TCa = mc[iCa];
    mc[iSr] = SrMix[kk];
    mc[iBa] = BaMix[kk];
    mc[iFe] = FeMix[kk];
    mc[iZn] = ZnMix[kk];
    mc[iPb] = PbMix[kk];
    mc[iRa] = RaMix[kk];

    ma[iOH] = OHstpMix[kk];
    ma[iCl] = ClMix[kk];
    ma[iAc] = ACstpMix[kk];
    mc[iNH4] = NH4STPMix[kk];
    ma[iH2BO3] = H2BO3stpMix[kk];
    ma[iHCO3] = HCO3stpMix[kk];
    ma[iCO3] = CO3stpMix[kk];

    ma[iH3SiO4] = 0;
    ma[iH2SiO4] = 0;

    ma[iSO4] = SO4Mix[kk];
    ma[iHS] = HSstpMix[kk];
    ma[intF] = FMix[kk];
    ma[iBr] = BrMix[kk];

    *Alk = AlkMix[kk];
    *TAc = TAcMix[kk];
    *TH2Saq = TH2SaqMix[kk];
    *TH4SiO4 = TH4SiO4Mix[kk];
    *TH3BO3 = TH3BO3Mix[kk];
    *TNH4 = TNH4Mix[kk];

    *TFe = FeMix[kk];
    *TPb = PbMix[kk];
    *TZn = ZnMix[kk];

    mn[iNH3] = 0;
    mn[iH3BO3] = *TH3BO3;
    mn[iH4SiO4aq] = *TH4SiO4;

    *use_pH = usepHmix[kk];
    *UseH2Sgas = UseH2SgasMix[kk];

    if (*use_pH == 3) {
        *TCO2 = TCO2Mix[kk];
    } else {
        *TCO2 = 0;
    }

    *yCO2 = yCO2Mix[kk];
    *yH2S = yH2SMix[kk];  // yCH4 = 1 - yCO2 - yH2S 'assume dry gas

    // If yCH4 < 0 Then yCH4 = 0  '?????????????????????????????????????????????????????????????

    if (useEOSmix[kk] == 1 && *yCO2 == 0 && SumofZMix[kk] > 0) *yCO2 = zMix[kk][2];  // if UseEOS=1 then set YCO2 and YH2S to reservoir condition to calculate density and TDS only if the reservoir fluid comp is given
    if (useEOSmix[kk] == 1 && *yH2S == 0 && SumofZMix[kk] > 0) *yH2S = zMix[kk][3];

    // Call C1_ThermodynamicEquilConsts  'Only function of T,P, does not recalculate in D1_CalcDensity
    // Call PengRobinson3

    // Call D1_CalcDensity(kk) 'Calculate TDS, density and fix molality based on the predicted density and TDS
    D1_CalcDensity(mc, ma, mn, kk);

    CalculatedTDSMix[kk] = *TDS;

    HstpMix[kk] = mc[iH];
    NaMix[kk] = mc[iNa];
    KMix[kk] = mc[iK];
    MgMix[kk] = mc[iMg];
    CaMix[kk] = mc[iCa];
    *TCa = mc[iCa];

    SrMix[kk] = mc[iSr];
    BaMix[kk] = mc[iBa];
    FeMix[kk] = *TFe;
    ZnMix[kk] = mc[iZn];
    PbMix[kk] = mc[iPb];
    RaMix[kk] = mc[iRa];

    OHstpMix[kk] = ma[iOH];
    ClMix[kk] = ma[iCl];
    ACstpMix[kk] = ma[iAc];
    NH4STPMix[kk] = mc[iNH4];
    H2BO3stpMix[kk] = ma[iH2BO3];
    HCO3stpMix[kk] = ma[iHCO3];
    CO3stpMix[kk] = ma[iCO3];

    SO4Mix[kk] = ma[iSO4];
    HSstpMix[kk] = ma[iHS];
    FMix[kk] = ma[intF];
    BrMix[kk] = ma[iBr];

    rho25CMix[kk] = rho25c;

    H3SiO4Mix[kk] = ma[iH3SiO4];
    H2SiO4Mix[kk] = ma[iH2SiO4];

    NH3Mix[kk] = mn[iNH3];
    H4SiO4Mix[kk] = mn[iH4SiO4aq];
    H3BO3Mix[kk] = mn[iH3BO3];

    CO2aqMix[kk] = mn[iCO2aq];
    H2SaqMix[kk] = mn[iH2Saq];
    HACaqMix[kk] = mn[iHAcaq];

    AlkMix[kk] = *Alk;
    TAcMix[kk] = *TAc;
    TH2SaqMix[kk] = *TH2Saq;
    TH4SiO4Mix[kk] = *TH4SiO4;
    TNH4Mix[kk] = *TNH4;
    TH3BO3Mix[kk] = *TH3BO3;

    yCO2Mix[kk] = *yCO2;
    yH2SMix[kk] = *yH2S;
    yCH4Mix[kk] = 1 - *yCO2 - *yH2S;  // Set YCO2 and YH2S to the calculated value if pH option is used.

    if (UseTPVolMix[kk] == 0) WaterDensityMix[kk] = rho25c;

    if (UseMolal == 1) {
        *TDS = 0;
        *CalculateTDSDen = 0;  // Calculate TDS from density

        int iden;
        for (iden = 2; iden <= NumCat; iden++) {
            *CalculateTDSDen += 0.001 * mc[iden] * MWCat[iden];  // =Sum of g salt/Kg H2O
        }
        for (iden = 2; iden <= NumAn; iden++) {
            *CalculateTDSDen += 0.001 * ma[iden] * MWAn[iden];
        }
        for (iden = 2; iden <= NumNeut; iden++) {
            *CalculateTDSDen += 0.001 * mn[iden] * MWNeut[iden];
        }

        *TDS = *CalculateTDSDen / (1 + *CalculateTDSDen) * rho25c * 1000000.0;  // TDS in unit of mg/L,  numerator=(Kgsalt/KgH2O), denometer=(1+Kgsalt/KgH2O)=(Kgsoln/KgH2O);density Kgsoln/Lsoln
        CalculatedTDSMix[kk] = *TDS;
    }
label100: ;
}

/**
 * @brief 获取EOS参数（Get_EOS_Parameters）
 *
 * 此函数从输入源读取EOS（状态方程）参数，包括气体分子量、临界温度、临界压力、Omega参数、mf_c0、mf_c1以及kPr交互参数矩阵。
 * 支持PR（Peng-Robinson）和SRK（Soave-Redlich-Kwong）EOS类型。
 * 后续MySQL读取将替换当前输入源。
 *
 * @param NumGases 气体数量
 * @param EOS EOS类型字符串 ("PR" 或 "SRK")
 * @param MWgas 气体分子量数组 [NumGases+1]
 * @param TCr 临界温度数组 [NumGases+1]
 * @param PCr 临界压力数组 [NumGases+1]
 * @param Omega Omega参数数组 [NumGases+1]
 * @param mf_c0 mf_c0参数数组 [NumGases+1]
 * @param mf_c1 mf_c1参数数组 [NumGases+1]
 * @param kPr kPr交互参数矩阵 [NumGases+1][MaxGases+1]
 */
void Get_EOS_Parameters(int NumGases, char* EOS, double MWgas[], double TCr[], double PCr[], double Omega[], double mf_c0[], double mf_c1[], double kPr[][MaxGases]) {
    int i, j;

    for (i = 1; i <= NumGases; i++) {
        /*MySQL 读入*/ MWgas[i] = 0; // Worksheets("Input").Cells(4 + i, 24)
        /*MySQL 读入*/ TCr[i] = 0; // Worksheets("Input").Cells(4 + i, 26)
        /*MySQL 读入*/ PCr[i] = 0; // Worksheets("Input").Cells(4 + i, 27)
        /*MySQL 读入*/ Omega[i] = 0; // Worksheets("Input").Cells(4 + i, 28)
        
        if (strcmp(EOS, "PR") == 0) {
            /*MySQL 读入*/ mf_c0[i] = 0; // Worksheets("Input").Cells(4 + i, 29)
            /*MySQL 读入*/ mf_c1[i] = 0; // Worksheets("Input").Cells(4 + i, 30)
        }
        
        if (strcmp(EOS, "SRK") == 0) {
            /*MySQL 读入*/ mf_c0[i] = 0; // Worksheets("Input").Cells(4 + i, 31)
            /*MySQL 读入*/ mf_c1[i] = 0; // Worksheets("Input").Cells(4 + i, 32)
        }
        
        for (j = 1; j <= i; j++) {
            /*MySQL 读入*/ kPr[i][j] = 0; // Worksheets("Input").Cells(4 + i, 34 + j)
        }
    }
}

/**
 * @brief 计算伪组成（pseudo_composition）
 * 
 * 此函数基于API重力和气体比重计算油气水的伪组成，用于EOS计算。
 * 输出15个组分的摩尔组成以兼容其他例程（Multiflash）。
 * 条件假设在给定TK、PBar下。
 * 
 * @param API 伪API重力（给定TP下的）
 * @param SGG 气体比重（给定TP下气体密度/空气密度）
 * @param VgTP 气体体积（MSCFPD，给定TP下）
 * @param mol_opd 油摩尔数（per day）
 * @param mol_wpd 水摩尔数（per day）
 * @param TK 绝对温度（K）
 * @param PBar 压力（bar）
 * @param aH2O 水活度
 * @param gNeut 中性活度系数数组 [2]
 * @param nTCO2 CO2摩尔数
 * @param nTH2S H2S摩尔数
 * @param yCO2 CO2摩尔分数
 * @param yH2S H2S摩尔分数
 * @param YH2O 水蒸气摩尔分数
 * @param total_moles [输出] 总摩尔数（per day）
 * @param feed_Composition [输出] 进料组成数组（15个组分摩尔分数）[15]
 * @param mol_HC [输出] 烃类总摩尔数（per day）
 */
void pseudo_composition(double API, double SGG, double VgTP, double mol_opd, double mol_wpd, double TK, double PBar, double aH2O, double gNeut[], double nTCO2, double nTH2S, double yCO2, double yH2S, double YH2O, double *total_moles, double feed_Composition[], double *mol_HC) {
    // ======= INPUT =======
    // API: Pseudo API gravity at given TP
    // SGG: Specific gravity of gas (density of gas / density of air at given TP)
    // VgTP: Volume of Gas in MSCFPD (thousands of cubic feet per day at given TP)
    // mol_opd: Moles of oil per day
    // mol_wpd: Moles of Water per day
    // TK, Temperature, [Kelvin]
    // PBar,  Pressure, [Bar]

    // Conditions are assumed at a given TK, PBar

    // ======= OUTPUT =======
    // feed_Composition: Molar composition of 15 components to make it compatible with other routines (Multiflash)

    double SGO, sgLightO, sgHeavyO, mwLightO, mwHeavyO, mwLightG, mwHeavyG, mwAir;
    double xLightO, mwO, mwG, yLightG, yHeavyG;
    double mol_gpd, mol_LightOpd, mol_HeavyOpd, mol_LightGpd;
    double mol_HeavyGpd;
    double composition_G[15];
    double MWgas[15], TCr[15], PCr[15], Omega[15], mf_c0[15], mf_c1[15], kPr[15][15];
    double lnphi_Gas[15];
    double lnphiH2O_Aq, composition_Aq[15], lnphi_Water[15], Compr_composition_Aq;
    double Compr_G;
    int i, j;
    char* EOS;
    double mf_gNeut[2];

    // const double RBar = 0.000083144621;  // m3*Bar/mol*K

    EOS = "PR";
    int NumGases = 15;

    // ReDim MWgas(NumGases), TCr(NumGases), PCr(NumGases), Omega(NumGases), mf_c0(NumGases), mf_c1(NumGases), kPr(NumGases, NumGases)
    // ReDim composition_Aq(NumGases), lnphi_Water(NumGases), composition_G(NumGases), lnphi_Gas(NumGases)

    SGO = 141.5 / (API + 131.5);  // specific gravity of the mixture (i.e. oil)

    sgLightO = 0.6548;  // specific gravity of hexane
    sgHeavyO = 0.92;  // specific gravity of C7+ equivalent

    mwLightO = 86.18;  // molecular weight of hexane
    mwHeavyO = 321.26;  // previously 450#, now molecular weight of C26-C80 fraction for consistency with critical properties

    mwLightG = 16.04;  // molecular weight of methane
    mwHeavyG = 58.12;  // molecular weight of butane
    mwAir = 28.97;  // molecular weight of air

    // MIN API = 22
    // MAX API = 84

    // min SGG = 0.55
    // max SGG = 2.00

    // calculation of mole fractions of hexane and toluene equivalent in oil phase (from API gravity)

    xLightO = mwHeavyO * (sgHeavyO - SGO) / (mwLightO * (sgHeavyO / sgLightO) * (SGO - sgLightO) + mwHeavyO * (sgHeavyO - SGO));  // mole fraction of hexane in the oil

    mwO = xLightO * mwLightO + (1 - xLightO) * mwHeavyO;  // calculated molecular weight of oil

    // calculation of mole fractions of methane and butane equivalent in the gas phase (from SG of gas)

    mwG = SGG * mwAir;  // calculated molar mass of gas

    // Subroutine that gets critical properties and parameters for the components
    Get_EOS_Parameters(NumGases, EOS, MWgas, TCr, PCr, Omega, mf_c0, mf_c1, kPr);

    mf_gNeut[0] = gNeut[0];
    mf_gNeut[1] = gNeut[1];
    // Sub phi_calc(eqVapor,eqAqueous, EOS As String, phase As String, TK As Double, PBar As Double, x() As Double, xGlobal() As Double, gNeut() As Double, aH2O As Double, TCr() As Double, PCr() As Double, Omega() As Double, c0() As Double, c1() As Double, kPr() As Double, lnphi() As Double, z As Double)

    // Call phi_calc(False, False, EOS, "liquid", CDbl(TK), CDbl(PBar), composition_Aq, composition_Aq, mf_gNeut, CDbl(aH2O), TCr, PCr, Omega, mf_c0, mf_c1, kPr, lnphi_Water, Compr_composition_Aq)  'commented out by Amy ???
    // lnphiH2O_Aq = lnphi_Water(NumGases)

    yLightG = (mwHeavyG - mwG - mwHeavyG * (yCO2 + yH2S + YH2O) + yCO2 * 44.01 + yH2S * 34.08 + YH2O * 18.01528) / (mwHeavyG - mwLightG);  // mole fraction of methane in the gas phase
    if (yLightG < 0) yLightG = 0;

    yHeavyG = 1 - (yLightG + yCO2 + yH2S + YH2O);
    if (yHeavyG < 0) yHeavyG = 0;
    // This calculates the compressibility factor of the vapor phase
    for (i = 3; i < NumGases - 1; i++) {  // 0-based: i=3 to 12 for components 4 to 13
        composition_G[i] = 0.0;
    }
    double denom = yLightG + yCO2 + yH2S + yHeavyG + YH2O;
    composition_G[0] = yLightG / denom;  // methane
    composition_G[1] = yCO2 / denom;  // CO2
    composition_G[2] = yH2S / denom;  // H2S
    composition_G[9] = yHeavyG / denom;  // n-butane (index 10 in 1-based)
    composition_G[14] = YH2O;  // water (index 15 in 1-based)

    phi_calc(0, 0, EOS, "vapor", TK, PBar, composition_G, composition_G, mf_gNeut, aH2O, TCr, PCr, Omega, mf_c0, mf_c1, kPr, lnphi_Gas, &Compr_G);
    mol_gpd = VgTP * PBar / (Compr_G * RBar * TK);  // moles of gas per day -- R is in m3/bar/K/g-mol

    mol_LightOpd = xLightO * mol_opd;  // moles of hexane equivalent per day

    mol_HeavyOpd = (1 - xLightO) * mol_opd;  // moles of toluene equivalent per day

    mol_LightGpd = yLightG * mol_gpd;  // moles of methane equivalent per day

    mol_HeavyGpd = (1 - yLightG) * mol_gpd;  // moles of butane equivalent per day
    // total_moles = mol_gpd + mol_wpd + mol_opd  ' TOTAL MOLES PER DAY (GAS + OIL + WATER)

    for (i = 0; i < NumGases; i++) {
        feed_Composition[i] = 0.0;
    }

    feed_Composition[0] = mol_LightGpd;  // mole fraction of methane in the feed
    feed_Composition[1] = nTCO2;  // mole fraction of CO2 in the feed
    feed_Composition[2] = nTH2S;  // mole fraction of H2S in the feed
    feed_Composition[6] = mol_HeavyGpd;  // mole fraction of n-butane in the feed (index 7 in 1-based)
    feed_Composition[9] = mol_LightOpd;  // mole fraction of n-hexane in the feed (index 10)
    feed_Composition[12] = mol_HeavyOpd;  // mole fraction of toluene in the feed (index 13)
    feed_Composition[14] = mol_wpd;  // mole fraction of water in the feed (index 15)

    *total_moles = 0;
    for (i = 0; i < NumGases; i++) {
        *total_moles += feed_Composition[i];
    }
    *mol_HC = *total_moles - mol_wpd;
    for (i = 0; i < NumGases; i++) {
        feed_Composition[i] /= *total_moles;
    }
}

int main() {
    return 0;
}