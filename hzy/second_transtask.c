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

int main() {
    return 0;
}