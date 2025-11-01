#include <stdio.h>
#include<mysql.h>//头文件
#include<stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#pragma warning(disable:4996)


double* NaMix;
double* MgMix;
double* CaMix;
double* SrMix;
double* BaMix;
double* FeMix;
double* ZnMix;
double* ClMix;
double* PbMix;
double* BrMix;
double* RaMix;

double* NH3Mix;
double* H3SiO4Mix;
double* H2SiO4Mix;
double* H4SiO4Mix;
double* H3BO3Mix;
double* CO2aqMix;
double* H2SaqMix;
double* HACaqMix;

double* UseH2SgasMix;
double* SO4Mix;
double* FMix;
double* TDSMix;
double* AlkMix;
double* TAcMix;
double* KMix;
double* MixFrac;

double* rho_Mix;
double* TH2SaqMix;
double* pHMeterStpMix;
double* TH4SiO4Mix;
double* TNH4Mix;
double* TH3BO3Mix;

double* SampleIDMix;
char* SampleDateMix;
double* OperatorMix;
char* WellNameMix;
double* FieldMix;
char* StateMix;

double* VgTPMix;
double* VoMix;
double* VwMix;
double* VMeOHMix;
double* VMEGMix;
double* oilAPIgravMix;
double* gasSpGravMix;
double* MixFracGas;

double* nTCO2Mix;
double* nTCH4Mix;
double* nTH2SMix;
double* mass_w_Mix;
double* mass_o_Mix;
double* MixFracOil;

double* mass_MeOH_mix;
double* mass_MEG_mix;
double* Qheat;
double* yCO2Mix;
double* yH2SMix;
double* yCH4Mix;

double* YCH4stpmix;
double* RatioOilBPointsmix;
double* CalculatedTDSMix;
double* rho25CMix;

double* HstpMix;
double* OHstpMix;
double* HCO3stpMix;
double* CO3stpMix;
double* ACstpMix;
double* HSstpMix;
double* NH4STPMix;
double* H2BO3stpMix;

double* HCO3AlkMix;
double* CO3AlkMix;
double* HAlkMix;
double* OHAlkMix;
double* ConcFactor;

double* TCO2Mix;
double* TofpH;
double* PofpH;
double* TofVol;
double* PofVol;

double* OilDensityMix;
double* GasDensityMix;
double* WaterDensityMix;
double* UseTPpHMix;
double* UseTPVolMix;
double* useEOSmix;

double* molAlk;
double* molTAC;
double* molTNH4;
double* molTH3BO3;
double* molTH2Saq;
double* molTH4SiO4;

double* mol_g_origMix;
double* mol_o_OrigMix;
double* mol_w_OrigMix;
double* mol_g_finalMix;
double* mol_o_finalMix;
double* mol_w_finalMix;
double* mol_w_evapMix;

double* Total_molesMix;
double* SumofZMix;
double* nTCO2MixEOS;
double* nTH2SMixEOS;

double* usepHmix;

int nob = 1;
int nob_Input = 2;
int nob_InputII = 0;
int Read_InputII = 0;
int Run1000Cases;
int CaseCount[2] = {1,2};

double LoopTP1000Cases;

int j;
int kk;

double** zMix;

// ????无法确定数据类型

int RunStat;

int RunH2SGUI;

int Run_CalcConcFactor;

double rho25c;

double RunQualityControlChecks;

double TBH;
int RunMultiMix;

double RunWhatIf;
double TWIInit;
double UseSI;
double TWH;
double PBH;
double PWIInit;

double PWH;

double tInh;
double SelectInh;
double InhNo;
double InhNo1;
double FracInhNo1;
double InhNo2;

char InhName[21][50]; // 最长为50，C语言不支持字符串，这一点很要命

double MaxInh;
double NoRiskcalc;
double ConcInhCelRisk[10];
double ConcInhCalRisk[10];
double ConcInhBarRisk[10];
double ConcInhAnRisk[10];
double ConcInhGypRisk[10];

MYSQL ConnectMysql(void) {
    MYSQL mysql;//定义一个MySQL类型的变量

    if (mysql_init(&mysql) == NULL)//mysql变量初始化
    {
        printf("error\n");
        exit(1);//如果发生错误就退出程序
    }
    if (mysql_library_init(0, NULL, NULL) != 0)//数据库初始化
    {
        printf("error\n");
        exit(1);
    }
    if (mysql_real_connect(&mysql, "localhost", "root", "19742017Ab", "test_new_new", 3306, 0, 0) == NULL)//连接
    {
        printf("%s\n", mysql_error(&mysql));//输出错误原因(字符串)
        exit(1);
    }
    if (mysql_set_character_set(&mysql, "utf8mb4") != 0)//设置字符集为utf8
    {
        printf("erroe\n");
        exit(1);
    }


    return mysql;
}



MYSQL_RES* search(char* sql) {

    MYSQL mysql;
    MYSQL_RES* res;
    mysql = ConnectMysql();
    int flag = mysql_real_query(&mysql, sql, (unsigned int)strlen(sql));

    /*mysql_store_result讲全部的查询结果读取到客户端*/
    res = mysql_store_result(&mysql);
    return res;
}

typedef struct
{

    /* ---------- sample_basic_information ---------- */
    char SampleID[256]; // 样品ID
    char Date[256];     // 采样日期
    char Operator[256]; // 采样人员
    char WellName[256]; // 油井名称
    char Location[256]; // 所在位置
    char Field[256];    // 油田名称

    /* ---------- sample_composition_information ---------- */
    char SampleInfo[256]; // 样品信息（用户自由输入）
    double Na_aq;
    double K_aq;
    double Mg_aq;
    double Ca_aq;
    double Sr_aq;
    double Ba_aq;
    double FeII_aq;
    double Zn_aq;
    double Pb_aq;
    double Cl_aq;
    double SO4_aq;
    double F_aq;
    double Br_aq;
    double Si_aq;
    double FeIII_aq;
    double Li_aq;
    double Be_aq;
    double Ra_aq;
    double Mn_aq;
    double Cu_aq;
    double Al_aq;
    double P_aq;
    double I_aq;
    double U_aq;
    double Alk_Bicarbonate_aq;
    double Alk_Carbonate_aq;
    double OrgAcid_Acetate_aq;
    double Ammonia_aq;
    double B_aq;
    double TDS_aq;
    double Density_STP;
    double CO2_pct_g;
    int Option_Use_H2Sg;
    double H2S_pct_g;
    double H2S_aq;
    double pH_STP;
    double Q_Gas;
    double Q_Oil;
    double Q_Water;

    /* ---------- sample_temperature_pressure_information ---------- */
    double T_initial;
    double T_final;
    double P_initial;
    double P_final;
    double API;
    double SG_g;
    double Q_MeOH;
    double Q_MEG;
    double StrongAcid_aq;
    double StrongBase_aq;
    double Conc_Multiplier;
    double T_pH;
    double P_pH;
    double T_Q;
    double P_Q;

    /* ---------- sample_oil_phase_information ---------- */
    double C1_o;
    double CO2_o;
    double H2S_o;
    double C2_o;
    double C3_o;
    double iC4_o;
    double nC4_o;
    double iC5_o;
    double nC5_o;
    double C6_o;
    double C7_C12_o;
    double C13_C25_o;
    double C26_C80_o;
    double N2_o;

    /* ---------- config_options_information ---------- */
    int Option_Alk;
    int Option_Defined_TP;
    int Option_TP_for_pH;
    int Option_TP_for_Q;
    int Option_EoS;
    int Option_Water_HC;

} SampleData;


void writeSampleBasicInformation(MYSQL_RES* res, SampleData* data)
{
    MYSQL_ROW row;
    while ((row = mysql_fetch_row(res)))
    {
        strcpy(data->SampleID, row[0] ? row[0] : "");
        strcpy(data->Date, row[1] ? row[1] : "");
        strcpy(data->Operator, row[2] ? row[2] : "");
        strcpy(data->WellName, row[3] ? row[3] : "");
        strcpy(data->Location, row[4] ? row[4] : "");
        strcpy(data->Field, row[5] ? row[5] : "");
        printf("test %s \n", row[2]);
    }
    printf("Write Sample Basic Information Success\n");
}

// ---------- 写入样品成分信息 ----------
void writeSampleCompositionInformation(MYSQL_RES* res, SampleData* data)
{
    MYSQL_ROW row;
    while ((row = mysql_fetch_row(res)))
    {
        strcpy(data->SampleID, row[0] ? row[0] : "");
        strcpy(data->SampleInfo, row[1] ? row[1] : "");


        data->Na_aq = atof(row[2]);
        data->K_aq = atof(row[3]);
        data->Mg_aq = atof(row[4]);
        data->Ca_aq = atof(row[5]);
        data->Sr_aq = atof(row[6]);
        data->Ba_aq = atof(row[7]);
        data->FeII_aq = atof(row[8]);
        data->Zn_aq = atof(row[9]);
        data->Pb_aq = atof(row[10]);
        data->Cl_aq = atof(row[11]);
        data->SO4_aq = atof(row[12]);
        data->F_aq = atof(row[13]);
        data->Br_aq = atof(row[14]);
        data->Si_aq = atof(row[15]);
        data->FeIII_aq = atof(row[16]);
        data->Li_aq = atof(row[17]);
        data->Be_aq = atof(row[18]);
        data->Ra_aq = atof(row[19]);
        data->Mn_aq = atof(row[20]);
        data->Cu_aq = atof(row[21]);
        data->Al_aq = atof(row[22]);
        data->P_aq = atof(row[23]);
        data->I_aq = atof(row[24]);
        data->U_aq = atof(row[25]);
        data->Alk_Bicarbonate_aq = atof(row[26]);
        data->Alk_Carbonate_aq = atof(row[27]);
        data->OrgAcid_Acetate_aq = atof(row[28]);
        data->Ammonia_aq = atof(row[29]);
        data->B_aq = atof(row[30]);
        data->TDS_aq = atof(row[31]);
        data->Density_STP = atof(row[32]);
        data->CO2_pct_g = atof(row[33]);
        data->Option_Use_H2Sg = atoi(row[34]);
        data->H2S_pct_g = atof(row[35]);
        data->H2S_aq = atof(row[36]);
        data->pH_STP = atof(row[37]);
        data->Q_Gas = atof(row[38]);
        data->Q_Oil = atof(row[39]);
        data->Q_Water = atof(row[40]);
    }
    printf("Write Sample Composition Information Success\n");
}

// ---------- 写入温压条件信息 ----------
void writeSampleTemperaturePressureInformation(MYSQL_RES* res, SampleData* data)
{
    MYSQL_ROW row;
    while ((row = mysql_fetch_row(res)))
    {
        strcpy(data->SampleID, row[0] ? row[0] : "");

        data->T_initial = atof(row[1]);
        data->T_final = atof(row[2]);
        data->P_initial = atof(row[3]);
        data->P_final = atof(row[4]);
        data->API = atof(row[5]);
        data->SG_g = atof(row[6]);
        data->Q_MeOH = atof(row[7]);
        data->Q_MEG = atof(row[8]);
        data->StrongAcid_aq = atof(row[9]);
        data->StrongBase_aq = atof(row[10]);
        data->Conc_Multiplier = atof(row[11]);
        data->T_pH = atof(row[12]);
        data->P_pH = atof(row[13]);
        data->T_Q = atof(row[14]);
        data->P_Q = atof(row[15]);
    }
    printf("Write Sample Temperature & Pressure Information Success\n");
}

// ---------- 写入油相信息 ----------
void writeSampleOilPhaseInformation(MYSQL_RES* res, SampleData* data)
{
    MYSQL_ROW row;
    while ((row = mysql_fetch_row(res)))
    {
        strcpy(data->SampleID, row[0] ? row[0] : "");

        data->C1_o = atof(row[1]);
        data->CO2_o = atof(row[2]);
        data->H2S_o = atof(row[3]);
        data->C2_o = atof(row[4]);
        data->C3_o = atof(row[5]);
        data->iC4_o = atof(row[6]);
        data->nC4_o = atof(row[7]);
        data->iC5_o = atof(row[8]);
        data->nC5_o = atof(row[9]);
        data->C6_o = atof(row[10]);
        data->C7_C12_o = atof(row[11]);
        data->C13_C25_o = atof(row[12]);
        data->C26_C80_o = atof(row[13]);
        data->N2_o = atof(row[14]);
    }
    printf("Write Sample Oil Phase Information Success\n");
}

// ---------- 写入配置参数信息 ----------
void writeConfigOptionsInformation(MYSQL_RES* res, SampleData* data)
{
    MYSQL_ROW row;
    while ((row = mysql_fetch_row(res)))
    {
        data->Option_Alk = atoi(row[0]);
        data->Option_Defined_TP = atoi(row[1]);
        data->Option_TP_for_pH = atoi(row[2]);
        data->Option_TP_for_Q = atoi(row[3]);
        data->Option_EoS = atoi(row[4]);
        data->Option_Water_HC = atoi(row[5]);
    }
    printf("Write Config Options Information Success\n");
}

void printSampleData(const SampleData* data)
{
    printf("\n================= Sample Basic Information =================\n");
    printf("SampleID: %s\n", data->SampleID);
    printf("Date: %s\n", data->Date);
    printf("Operator: %s\n", data->Operator);
    printf("WellName: %s\n", data->WellName);
    printf("Location: %s\n", data->Location);
    printf("Field: %s\n", data->Field);

    printf("\n================= Sample Composition Information =================\n");
    printf("SampleInfo: %s\n", data->SampleInfo);
    printf("Na_aq: %.4f\n", data->Na_aq);
    printf("K_aq: %.4f\n", data->K_aq);
    printf("Mg_aq: %.4f\n", data->Mg_aq);
    printf("Ca_aq: %.4f\n", data->Ca_aq);
    printf("Sr_aq: %.4f\n", data->Sr_aq);
    printf("Ba_aq: %.4f\n", data->Ba_aq);
    printf("FeII_aq: %.4f\n", data->FeII_aq);
    printf("Zn_aq: %.4f\n", data->Zn_aq);
    printf("Pb_aq: %.4f\n", data->Pb_aq);
    printf("Cl_aq: %.4f\n", data->Cl_aq);
    printf("SO4_aq: %.4f\n", data->SO4_aq);
    printf("F_aq: %.4f\n", data->F_aq);
    printf("Br_aq: %.4f\n", data->Br_aq);
    printf("Si_aq: %.4f\n", data->Si_aq);
    printf("FeIII_aq: %.4f\n", data->FeIII_aq);
    printf("Li_aq: %.4f\n", data->Li_aq);
    printf("Be_aq: %.4f\n", data->Be_aq);
    printf("Ra_aq: %.4f\n", data->Ra_aq);
    printf("Mn_aq: %.4f\n", data->Mn_aq);
    printf("Cu_aq: %.4f\n", data->Cu_aq);
    printf("Al_aq: %.4f\n", data->Al_aq);
    printf("P_aq: %.4f\n", data->P_aq);
    printf("I_aq: %.4f\n", data->I_aq);
    printf("U_aq: %.4f\n", data->U_aq);
    printf("Alk_Bicarbonate_aq: %.4f\n", data->Alk_Bicarbonate_aq);
    printf("Alk_Carbonate_aq: %.4f\n", data->Alk_Carbonate_aq);
    printf("OrgAcid_Acetate_aq: %.4f\n", data->OrgAcid_Acetate_aq);
    printf("Ammonia_aq: %.4f\n", data->Ammonia_aq);
    printf("B_aq: %.4f\n", data->B_aq);
    printf("TDS_aq: %.4f\n", data->TDS_aq);
    printf("Density_STP: %.4f\n", data->Density_STP);
    printf("CO2_pct_g: %.4f\n", data->CO2_pct_g);
    printf("Option_Use_H2Sg: %d\n", data->Option_Use_H2Sg);
    printf("H2S_pct_g: %.4f\n", data->H2S_pct_g);
    printf("H2S_aq: %.4f\n", data->H2S_aq);
    printf("pH_STP: %.4f\n", data->pH_STP);
    printf("Q_Gas: %.4f\n", data->Q_Gas);
    printf("Q_Oil: %.4f\n", data->Q_Oil);
    printf("Q_Water: %.4f\n", data->Q_Water);

    printf("\n================= Temperature & Pressure Information =================\n");
    printf("T_initial: %.4f\n", data->T_initial);
    printf("T_final: %.4f\n", data->T_final);
    printf("P_initial: %.4f\n", data->P_initial);
    printf("P_final: %.4f\n", data->P_final);
    printf("API: %.4f\n", data->API);
    printf("SG_g: %.4f\n", data->SG_g);
    printf("Q_MeOH: %.4f\n", data->Q_MeOH);
    printf("Q_MEG: %.4f\n", data->Q_MEG);
    printf("StrongAcid_aq: %.4f\n", data->StrongAcid_aq);
    printf("StrongBase_aq: %.4f\n", data->StrongBase_aq);
    printf("Conc_Multiplier: %.4f\n", data->Conc_Multiplier);
    printf("T_pH: %.4f\n", data->T_pH);
    printf("P_pH: %.4f\n", data->P_pH);
    printf("T_Q: %.4f\n", data->T_Q);
    printf("P_Q: %.4f\n", data->P_Q);

    printf("\n================= Oil Phase Information =================\n");
    printf("C1_o: %.4f\n", data->C1_o);
    printf("CO2_o: %.4f\n", data->CO2_o);
    printf("H2S_o: %.4f\n", data->H2S_o);
    printf("C2_o: %.4f\n", data->C2_o);
    printf("C3_o: %.4f\n", data->C3_o);
    printf("iC4_o: %.4f\n", data->iC4_o);
    printf("nC4_o: %.4f\n", data->nC4_o);
    printf("iC5_o: %.4f\n", data->iC5_o);
    printf("nC5_o: %.4f\n", data->nC5_o);
    printf("C6_o: %.4f\n", data->C6_o);
    printf("C7_C12_o: %.4f\n", data->C7_C12_o);
    printf("C13_C25_o: %.4f\n", data->C13_C25_o);
    printf("C26_C80_o: %.4f\n", data->C26_C80_o);
    printf("N2_o: %.4f\n", data->N2_o);

    printf("\n================= Configuration Options =================\n");
    printf("Option_Alk: %d\n", data->Option_Alk);
    printf("Option_Defined_TP: %d\n", data->Option_Defined_TP);
    printf("Option_TP_for_pH: %d\n", data->Option_TP_for_pH);
    printf("Option_TP_for_Q: %d\n", data->Option_TP_for_Q);
    printf("Option_EoS: %d\n", data->Option_EoS);
    printf("Option_Water_HC: %d\n", data->Option_Water_HC);

    printf("===========================================================\n");
}


void fetchData(SampleData *data) {
    char sql[1024];
    snprintf(sql, sizeof(sql), "SELECT * FROM sample_basic_information Where SampleId ='SMP01'");
    MYSQL_RES* res = search(sql);
    writeSampleBasicInformation(res, data);

    snprintf(sql, sizeof(sql), "SELECT * FROM sample_composition_information Where SampleId ='SMP01'");
    res = search(sql);
    writeSampleCompositionInformation(res, data);

    snprintf(sql, sizeof(sql), "SELECT * FROM sample_temperature_pressure_information  Where SampleId ='SMP01'");
    res = search(sql);
    writeSampleTemperaturePressureInformation(res, data);

    snprintf(sql, sizeof(sql), "SELECT * FROM sample_oil_phase_information  Where SampleId ='SMP01'");
    res = search(sql);
    writeSampleOilPhaseInformation(res, data);

    snprintf(sql, sizeof(sql), "SELECT * FROM config_options_information");
    res = search(sql);
    writeConfigOptionsInformation(res, data);

}

void initData()
{
    int len = nob_Input + nob_InputII;
    NaMix = (double*)malloc(len * sizeof(double));
    MgMix = (double*)malloc(len * sizeof(double));
    CaMix = (double*)malloc(len * sizeof(double));
    SrMix = (double*)malloc(len * sizeof(double));
    BaMix = (double*)malloc(len * sizeof(double));
    FeMix = (double*)malloc(len * sizeof(double));
    ZnMix = (double*)malloc(len * sizeof(double));
    ClMix = (double*)malloc(len * sizeof(double));
    PbMix = (double*)malloc(len * sizeof(double));
    BrMix = (double*)malloc(len * sizeof(double));
    RaMix = (double*)malloc(len * sizeof(double));

    NH3Mix = (double*)malloc(len * sizeof(double));
    H3SiO4Mix = (double*)malloc(len * sizeof(double));
    H2SiO4Mix = (double*)malloc(len * sizeof(double));
    H4SiO4Mix = (double*)malloc(len * sizeof(double));
    H3BO3Mix = (double*)malloc(len * sizeof(double));
    CO2aqMix = (double*)malloc(len * sizeof(double));
    H2SaqMix = (double*)malloc(len * sizeof(double));
    HACaqMix = (double*)malloc(len * sizeof(double));

    UseH2SgasMix = (double*)malloc(len * sizeof(double));
    SO4Mix = (double*)malloc(len * sizeof(double));
    FMix = (double*)malloc(len * sizeof(double));
    TDSMix = (double*)malloc(len * sizeof(double));
    AlkMix = (double*)malloc(len * sizeof(double));
    TAcMix = (double*)malloc(len * sizeof(double));
    KMix = (double*)malloc(len * sizeof(double));
    MixFrac = (double*)malloc(len * sizeof(double));

    rho_Mix = (double*)malloc(len * sizeof(double));
    TH2SaqMix = (double*)malloc(len * sizeof(double));
    pHMeterStpMix = (double*)malloc(len * sizeof(double));
    TH4SiO4Mix = (double*)malloc(len * sizeof(double));
    TNH4Mix = (double*)malloc(len * sizeof(double));
    TH3BO3Mix = (double*)malloc(len * sizeof(double));
    SampleIDMix = (double*)malloc(len * sizeof(double));
    SampleDateMix = (char*)malloc(len * sizeof(char));
    OperatorMix = (double*)malloc(len * sizeof(double));
    WellNameMix = (char*)malloc(len * sizeof(char));
    FieldMix = (double*)malloc(len * sizeof(double));
    StateMix = (char*)malloc(len * sizeof(char));
    VgTPMix = (double*)malloc(len * sizeof(double));
    VoMix = (double*)malloc(len * sizeof(double));
    VwMix = (double*)malloc(len * sizeof(double));
    VMeOHMix = (double*)malloc(len * sizeof(double));
    VMEGMix = (double*)malloc(len * sizeof(double));
    oilAPIgravMix = (double*)malloc(len * sizeof(double));
    gasSpGravMix = (double*)malloc(len * sizeof(double));
    MixFracGas = (double*)malloc(len * sizeof(double));
    nTCO2Mix = (double*)malloc(len * sizeof(double));
    nTCH4Mix = (double*)malloc(len * sizeof(double));
    nTH2SMix = (double*)malloc(len * sizeof(double));
    mass_w_Mix = (double*)malloc(len * sizeof(double));
    mass_o_Mix = (double*)malloc(len * sizeof(double));
    MixFracOil = (double*)malloc(len * sizeof(double));
    mass_MeOH_mix = (double*)malloc(len * sizeof(double));
    mass_MEG_mix = (double*)malloc(len * sizeof(double));
    Qheat = (double*)malloc(len * sizeof(double));
    yCO2Mix = (double*)malloc(len * sizeof(double));
    yH2SMix = (double*)malloc(len * sizeof(double));
    yCH4Mix = (double*)malloc(len * sizeof(double));
    YCH4stpmix = (double*)malloc(len * sizeof(double));
    RatioOilBPointsmix = (double*)malloc(len * sizeof(double));
    CalculatedTDSMix = (double*)malloc(len * sizeof(double));
    rho25CMix = (double*)malloc(len * sizeof(double));
    HstpMix = (double*)malloc(len * sizeof(double));
    OHstpMix = (double*)malloc(len * sizeof(double));
    HCO3stpMix = (double*)malloc(len * sizeof(double));
    CO3stpMix = (double*)malloc(len * sizeof(double));
    ACstpMix = (double*)malloc(len * sizeof(double));
    HSstpMix = (double*)malloc(len * sizeof(double));
    NH4STPMix = (double*)malloc(len * sizeof(double));
    H2BO3stpMix = (double*)malloc(len * sizeof(double));
    HCO3AlkMix = (double*)malloc(len * sizeof(double));
    CO3AlkMix = (double*)malloc(len * sizeof(double));
    HAlkMix = (double*)malloc(len * sizeof(double));
    OHAlkMix = (double*)malloc(len * sizeof(double));
    ConcFactor = (double*)malloc(len * sizeof(double));
    TCO2Mix = (double*)malloc(len * sizeof(double));
    TofpH = (double*)malloc(len * sizeof(double));
    PofpH = (double*)malloc(len * sizeof(double));
    TofVol = (double*)malloc(len * sizeof(double));
    PofVol = (double*)malloc(len * sizeof(double));
    OilDensityMix = (double*)malloc(len * sizeof(double));
    GasDensityMix = (double*)malloc(len * sizeof(double));
    WaterDensityMix = (double*)malloc(len * sizeof(double));
    UseTPpHMix = (double*)malloc(len * sizeof(double));
    UseTPVolMix = (double*)malloc(len * sizeof(double));
    useEOSmix = (double*)malloc(len * sizeof(double));
    molAlk = (double*)malloc(len * sizeof(double));
    molTAC = (double*)malloc(len * sizeof(double));
    molTNH4 = (double*)malloc(len * sizeof(double));
    molTH3BO3 = (double*)malloc(len * sizeof(double));
    molTH2Saq = (double*)malloc(len * sizeof(double));
    molTH4SiO4 = (double*)malloc(len * sizeof(double));
    mol_g_origMix = (double*)malloc(len * sizeof(double));
    mol_o_OrigMix = (double*)malloc(len * sizeof(double));
    mol_w_OrigMix = (double*)malloc(len * sizeof(double));
    mol_g_finalMix = (double*)malloc(len * sizeof(double));
    mol_o_finalMix = (double*)malloc(len * sizeof(double));
    mol_w_finalMix = (double*)malloc(len * sizeof(double));
    mol_w_evapMix = (double*)malloc(len * sizeof(double));
    Total_molesMix = (double*)malloc(len * sizeof(double));
    SumofZMix = (double*)malloc(len * sizeof(double));
    nTCO2MixEOS = (double*)malloc(len * sizeof(double));
    nTH2SMixEOS = (double*)malloc(len * sizeof(double));
    usepHmix = (double*)malloc(len * sizeof(double));

    /*zMix = malloc((nob_Input + nob_InputII) * sizeof(double*));
    for (int i = 0; i < (nob_Input + nob_InputII); i++)
    {
        zMix[i] = malloc(15 * sizeof(double));
    }*/
}



int Run10TestCases;
int Loop10;
int Run_Seawater_Mixing;
int LoopMixing;
int LoopResChem;
int Run_MixingTwoWells;

// 源代码中未说明类型，按逻辑推测为int型
int RunNORM;
//

double UseTPCalciteSheet;
int useTPVol;
int ING;


int CaseCountWI[40];
int Loop1WI;
int Loop2WI;

double NaWI[40][10];
double CaWI[40][10];
double SrWI[40][10];
double BaWI[40][10];
double ClWI[40][10];
double SO4WI[40][10];
double HCO3AlkWI[40][10];
double CO3AlkWI[40][10];
double TACWI[40][10];
double YCO2WI[40][10];
double YH2SWI[40][10];
double TH2SaqWI[40][10];
double pHMeterSTPWI[40][10];
double ConcFactorWI[40][10];

int UseMolal;

int errmsg[20];

double yCO2;
double CO2aq;
double HCO3;
double CO3;
double yH2S;
double HS;
double H2Saq;

void ReadInputPartA(int kk,SampleData* data)
{
    if (Run10TestCases == 1 && Loop10 > 1)
    {
        return;
    }
    if (Run_Seawater_Mixing == 1 && LoopMixing > 1)
    {
        return;
    }
    if (Run_MixingTwoWells == 1 && LoopMixing > 1)
    {
        return;
    }
    if (RunMultiMix == 1 && LoopResChem > 1)
    {
        return;
    }

    // 源代码中的这个逻辑为判断 Gas/Day 的单位，但是我们在数据库中没有设置这个单位的字段
    // If Worksheets(mySheet).Cells(35, 2).Value = "(kSm^3/D)" Then UseSI = 1

    // 源代码中的这个逻辑为判断 Na+ 的单位，但是我们在数据库中没有设置这个单位的字段
    // If Worksheets(mySheet).Cells(10, 2).Value = "(m)" Then UseMolal = 1

    //SampleDateMix[kk] = Date;
    //WellNameMix[kk] = WellName;
    NaMix[kk] = data->Na_aq;
    KMix[kk] = data->K_aq;
    MgMix[kk] = data->Mg_aq;
    CaMix[kk] = data->Ca_aq;
    SrMix[kk] = data->Sr_aq;
    BaMix[kk] = data->Ba_aq;
    FeMix[kk] = data->FeII_aq;
    ZnMix[kk] = data->Zn_aq;
    PbMix[kk] = data->Pb_aq;
    ClMix[kk] = data->Cl_aq;
    SO4Mix[kk] = data->SO4_aq;
    FMix[kk] = data->F_aq;
    BrMix[kk] = data->Br_aq;
    TH4SiO4Mix[kk] = data->Si_aq;              // TH4SiO4Mix(kk) = Worksheets(mySheet).Cells(23, j + 2).Value 其中23行对应Silica
    HCO3AlkMix[kk] = data->Alk_Bicarbonate_aq; // HCO3AlkMix(kk) = Worksheets(mySheet).Cells(24, j + 2).Value 其中24行对应Total Alkalinity
    CO3AlkMix[kk] = data->Alk_Carbonate_aq;    // CO3AlkMix(kk) = Worksheets(mySheet).Cells(25, j + 2).Value 其中25行对应 CO3 Alkalinity
    TAcMix[kk] = data->OrgAcid_Acetate_aq;     // TAcMix(kk) = Worksheets(mySheet).Cells(26, j + 2).Value 其中26行对应Carboxylates
    TNH4Mix[kk] = data->Ammonia_aq;            // TNH4Mix(kk) = Worksheets(mySheet).Cells(27, j + 2).Value其中27行对应 Ammonia
    TH3BO3Mix[kk] = data->B_aq;                // TH3BO3Mix(kk) = Worksheets(mySheet).Cells(28, j + 2).Value其中28行对应Borate;
    yCO2Mix[kk] = data->CO2_pct_g / 100;       // yCO2Mix(kk) = Worksheets(mySheet).Cells(31, j + 2).Value / 100 其中31行对应 Co2 Gas Analysis
    UseH2SgasMix[kk] = data->Option_Use_H2Sg;  // UseH2SgasMix(kk) = Worksheets(mySheet).Cells(32, j + 2).Value 其中32行对应 Use H2S Gas Analysis
    if (UseH2SgasMix[kk] == 1)
    {
        yH2SMix[kk] = data->H2S_pct_g / 100; // yH2SMix(kk) = Worksheets(mySheet).Cells(33, j + 2).Value / 100 其中33行对应Gas H2S% or H2Saq
        TH2SaqMix[kk] = 0;
    }
    else
    {
        TH2SaqMix[kk] = data->H2S_pct_g; // TH2SaqMix(kk) = Worksheets(mySheet).Cells(33, j + 2).Value
        yH2SMix[kk] = 0;
    }
    pHMeterStpMix[kk] = data->pH_STP;       // pHMeterStpMix(kk) = Worksheets(mySheet).Cells(34, j + 2).Value其中34行对应pH2 measured;
    HAlkMix[kk] = data->StrongAcid_aq;      // HAlkMix(kk) = Worksheets(mySheet).Cells(47, j + 2).Value其中47行对应H+(Strong acid);
    OHAlkMix[kk] = data->StrongBase_aq;     // OHAlkMix(kk) = Worksheets(mySheet).Cells(48, j + 2).Value其中48行对应OH-(Strong base);
    ConcFactor[kk] = data->Conc_Multiplier; // ConcFactor(kk) = Worksheets(mySheet).Cells(49, j + 2).Value其中49行对应Conc. multiplier;
    usepHmix[kk] = data->Option_Alk;        // usepHmix(kk) = Worksheets(mySheet).Cells(51, j + 2).Value其中51行对应Four Cal options;
    if (usepHmix[kk] == 0)            // If usepHmix(kk) = Empty Then usepHmix(kk) = 0这里的Empty设为0
    {
        usepHmix[kk] = 0;
    }
    UseTPpHMix[kk] = data->Option_TP_for_pH; // UseTPpHMix(kk) = Worksheets(mySheet).Cells(53, j + 2).Value其中53行对应T,P for pH;
    UseTPVolMix[kk] = data->Option_TP_for_Q; // UseTPVolMix(kk) = Worksheets(mySheet).Cells(54, j + 2).Value其中54行对应 T,P for G/O/W;
    if (nob_Input + nob_InputII == 1)
    {
        UseTPCalciteSheet = data->Option_Defined_TP; // UseTPCalciteSheet = Worksheets(mySheet).Cells(52, j + 2).Value其中52行对应 Use TP on Calcite sheet?;
        useEOSmix[kk] = data->Option_EoS;            // useEOSmix(kk) = Worksheets(mySheet).Cells(55, j + 2).Value其中55行对应Use Flash Calculator;
    }
    else
    {
        // UseTPCalciteSheet = Worksheets(MySheetMix).Cells(52, 8).Value 这里找不到第8列对应着什么
        // useEOSmix(kk) = Worksheets(MySheetMix).Cells(55, 8).Value 同上
    }

    // If UseTPCalciteSheet = "" Then UseTPCalciteSheet = 0
    // If UseTPpHMix(kk) = "" Then UseTPpHMix(kk) = 0
    // If UseTPVolMix(kk) = "" Then UseTPVolMix(kk) = 0
    // If useEOSmix(kk) = "" Then useEOSmix(kk) = 0 这里的逻辑判断似乎有问题，因为这些几乎全为double数组，而在判断中判断是否为空字符串;

    TofpH[kk] = data->T_pH; // TofpH(kk) = Worksheets(mySheet).Cells(58, j + 2).Value其中58行对应Temp. for pH meas.;
    PofpH[kk] = data->P_pH; // PofpH(kk) = Worksheets(mySheet).Cells(59, j + 2).Value其中59行对应Pres. for pH meas.;
    TofVol[kk] = data->T_Q; // TofVol(kk) = Worksheets(mySheet).Cells(60, j + 2).Value其中60行对应 T for fluids meas.;
    PofVol[kk] = data->P_Q; // PofVol(kk) = Worksheets(mySheet).Cells(61, j + 2).Value其中61行对P for for fluids meas.;

    if (RunNORM == 1)
    {
        // RaMix(kk) = Worksheets(mySheet).Cells(62, j + 2).Value / 1000000000000#表中62行没有字段
    }
    else
    {
        RaMix[kk] = 0;
    }
    useTPVol = UseTPVolMix[kk];
    SumofZMix[kk] = 0;
    // For iNG = 1 To 14
    // zMix(kk, iNG) = Worksheets(mySheet).Cells(65 + iNG, j + 2) / 100:
    // SumofZMix(kk) = SumofZMix(kk) + zMix(kk, iNG) 这里是对原表中66-79行进行读取，对应为 sample_oil_phase_information数据表，这里不采用循环
    zMix[kk][0] = data->C1_o / 100;
    zMix[kk][1] = data->CO2_o / 100;
    zMix[kk][2] = data->H2S_o / 100;
    zMix[kk][3] = data->C2_o / 100;
    zMix[kk][4] = data->C3_o / 100;
    zMix[kk][5] = data->iC4_o / 100;
    zMix[kk][6] = data->nC4_o / 100;
    zMix[kk][7] = data->iC5_o / 100;
    zMix[kk][8] = data->nC5_o / 100;
    zMix[kk][9] = data->C6_o / 100;
    zMix[kk][10] = data->C7_C12_o / 100;
    zMix[kk][11] = data->C13_C25_o / 100;
    zMix[kk][12] = data->C26_C80_o / 100;
    zMix[kk][13] = data->N2_o / 100;

    for (int i = 0; i < 14; i++)
    {
        SumofZMix[kk] = SumofZMix[kk] + zMix[kk][i];
    }

    if (SumofZMix[kk] > 0)
    {
        for (int i = 0; i < 14; i++)
        {
            zMix[kk][i] = zMix[kk][i] / SumofZMix[kk];
        }
    }
    zMix[kk][14] = 0.0;
    if (RunH2SGUI != 1)
    {
        // MultiplePpt = Worksheets("Input").Range("S11").Value
    }

    if (RunWhatIf == 1)
    {
        if (CaseCountWI[Loop1WI])
        {
            NaMix[kk] = NaWI[Loop1WI][Loop2WI];
        }
        if (CaseCountWI[Loop1WI] == 4)
        {
            CaMix[kk] = CaWI[Loop1WI][Loop2WI];
        }
        if (CaseCountWI[Loop1WI] == 5)
        {
            SrMix[kk] = SrWI[Loop1WI][Loop2WI];
        }

        if (CaseCountWI[Loop1WI] == 6)
        {
            BaMix[kk] = BaWI[Loop1WI][Loop2WI];
        }
        if (CaseCountWI[Loop1WI] == 10)
        {
            ClMix[kk] = ClWI[Loop1WI][Loop2WI];
        }
        if (CaseCountWI[Loop1WI] == 11)
        {
            SO4Mix[kk] = SO4WI[Loop1WI][Loop2WI];
        }
        if (CaseCountWI[Loop1WI] == 15)
        {
            HCO3AlkMix[kk] = HCO3AlkWI[Loop1WI][Loop2WI];
        }

        if (CaseCountWI[Loop1WI] == 16)
        {
            CO3AlkMix[kk] = CO3AlkWI[Loop1WI][Loop2WI];
        }
        if (CaseCountWI[Loop1WI] == 17)
        {
            TAcMix[kk] = TACWI[Loop1WI][Loop2WI];
        }
        if (CaseCountWI[Loop1WI] == 19)
        {
            usepHmix[kk] = 0;
            yCO2Mix[kk] = YCO2WI[Loop1WI][Loop2WI] / 100;
        }

        if (CaseCountWI[Loop1WI] == 20)
        {
            UseH2SgasMix[kk] = 1;
            yH2SMix[kk] = YH2SWI[Loop1WI][Loop2WI] / 100;
        }

        if (CaseCountWI[Loop1WI] == 21)
        {
            UseH2SgasMix[kk] = 1;
            TH2SaqMix[kk] = TH2SaqWI[Loop1WI][Loop2WI] / 100;
        }

        if (CaseCountWI[Loop1WI] == 22)
        {
            usepHmix[kk] = 1;
            pHMeterStpMix[kk] = pHMeterSTPWI[Loop1WI][Loop2WI] / 100;
        }

        if (CaseCountWI[Loop1WI] == 32)
        {
            ConcFactor[kk] = ConcFactorWI[Loop1WI][Loop2WI] / 100;
        }
    }

    if (TDSMix[kk] = 0) // If TDSMix(kk) = 0 Or TDSMix(kk) = "" Then 这里不对空字符串做判断
    {
        if (UseMolal == 0)
        {
            NaMix[kk] = NaMix[kk] / (22990.0 * (0.9991 - 0.0000003612 * TDSMix[kk])) * ConcFactor[kk];
            KMix[kk] = KMix[kk] / (39098.0 * (0.9991 - 0.0000003612 * TDSMix[kk])) * ConcFactor[kk];
            MgMix[kk] = MgMix[kk] / (24305.0 * (0.9991 - 0.0000003612 * TDSMix[kk])) * ConcFactor[kk];
            CaMix[kk] = CaMix[kk] / (40080.0 * (0.9991 - 0.0000003612 * TDSMix[kk])) * ConcFactor[kk];
            SrMix[kk] = SrMix[kk] / (87620.0 * (0.9991 - 0.0000003612 * TDSMix[kk])) * ConcFactor[kk];
            BaMix[kk] = BaMix[kk] / (137330.0 * (0.9991 - 0.0000003612 * TDSMix[kk])) * ConcFactor[kk];
            FeMix[kk] = FeMix[kk] / (55847.0 * (0.9991 - 0.0000003612 * TDSMix[kk])) * ConcFactor[kk];
            ZnMix[kk] = ZnMix[kk] / (65380.0 * (0.9991 - 0.0000003612 * TDSMix[kk])) * ConcFactor[kk];
            PbMix[kk] = PbMix[kk] / (207200.0 * (0.9991 - 0.0000003612 * TDSMix[kk])) * ConcFactor[kk];
            ClMix[kk] = ClMix[kk] / (35450.0 * (0.9991 - 0.0000003612 * TDSMix[kk])) * ConcFactor[kk];
            SO4Mix[kk] = SO4Mix[kk] / (96064.0 * (0.9991 - 0.0000003612 * TDSMix[kk])) * ConcFactor[kk];
            FMix[kk] = FMix[kk] / (18998.0 * (0.9991 - 0.0000003612 * TDSMix[kk])) * ConcFactor[kk];
            BrMix[kk] = BrMix[kk] / (79904.0 * (0.9991 - 0.0000003612 * TDSMix[kk])) * ConcFactor[kk];
            TH4SiO4Mix[kk] = TH4SiO4Mix[kk] / (28085.0 * (0.9991 - 0.0000003612 * TDSMix[kk])) * ConcFactor[kk];
            HCO3AlkMix[kk] = HCO3AlkMix[kk] / (61019.0 * (0.9991 - 0.0000003612 * TDSMix[kk])) * ConcFactor[kk];
            CO3AlkMix[kk] = CO3AlkMix[kk] / (60019.0 * (0.9991 - 0.0000003612 * TDSMix[kk])) * ConcFactor[kk];
            TAcMix[kk] = TAcMix[kk] / (59046.0 * (0.9991 - 0.0000003612 * TDSMix[kk])) * ConcFactor[kk];
            TNH4Mix[kk] = TNH4Mix[kk] / (17031.0 * (0.9991 - 0.0000003612 * TDSMix[kk])) * ConcFactor[kk];
            TH3BO3Mix[kk] = TH3BO3Mix[kk] / (10811.0 * (0.9991 - 0.0000003612 * TDSMix[kk])) * ConcFactor[kk];

            if (UseH2SgasMix[kk] == 0)
            {
                TH2SaqMix[kk] = TH2SaqMix[kk] / (34080.0 * (0.9991 - 0.0000003612 * TDSMix[kk])) * ConcFactor[kk];
            }

            HAlkMix[kk] = HAlkMix[kk] / (0.9991 - 0.0000003612 * TDSMix[kk]);
            OHAlkMix[kk] = OHAlkMix[kk] / (0.9991 - 0.0000003612 * TDSMix[kk]);
            RaMix[kk] = RaMix[kk] / (226.0 * (0.9991 - 0.0000003612 * TDSMix[kk])) * ConcFactor[kk];
        }
        else if (UseMolal == 1)
        {
            NaMix[kk] = NaMix[kk] * ConcFactor[kk]; // Convert mg/L to molality
            KMix[kk] = KMix[kk] * ConcFactor[kk];
            MgMix[kk] = MgMix[kk] * ConcFactor[kk];
            CaMix[kk] = CaMix[kk] * ConcFactor[kk];
            SrMix[kk] = SrMix[kk] * ConcFactor[kk];
            BaMix[kk] = BaMix[kk] * ConcFactor[kk];
            FeMix[kk] = FeMix[kk] * ConcFactor[kk];
            ZnMix[kk] = ZnMix[kk] * ConcFactor[kk];
            PbMix[kk] = PbMix[kk] * ConcFactor[kk]; // Pb added
            ClMix[kk] = ClMix[kk] * ConcFactor[kk];
            SO4Mix[kk] = SO4Mix[kk] * ConcFactor[kk];
            FMix[kk] = FMix[kk] * ConcFactor[kk];
            BrMix[kk] = BrMix[kk] * ConcFactor[kk];           // Br added
            TH4SiO4Mix[kk] = TH4SiO4Mix[kk] * ConcFactor[kk]; // Input silica as SiO2
            HCO3AlkMix[kk] = HCO3AlkMix[kk] * ConcFactor[kk];
            CO3AlkMix[kk] = CO3AlkMix[kk] * ConcFactor[kk];
            TAcMix[kk] = TAcMix[kk] * ConcFactor[kk];
            TNH4Mix[kk] = TNH4Mix[kk] * ConcFactor[kk];
            TH3BO3Mix[kk] = TH3BO3Mix[kk] * ConcFactor[kk];

            if (UseH2SgasMix[kk] == 0)
            {
                TH2SaqMix[kk] = TH2SaqMix[kk] * ConcFactor[kk]; // Used to calculate yH2Sstp
            }

            // HAlkMix 和 OHAlkMix 被注释掉，不处理
            RaMix[kk] = RaMix[kk] * ConcFactor[kk];
        }
        AlkMix[kk] = HCO3AlkMix[kk] + 2 * CO3AlkMix[kk] - HAlkMix[kk] + OHAlkMix[kk];
        TCO2Mix[kk] = HCO3AlkMix[kk] + CO3AlkMix[kk];
        if (UseTPCalciteSheet != 1)
        {
            UseTPCalciteSheet = 0;
        }

        if (Run_Seawater_Mixing == 1 && j == 2)
        {
            yCO2Mix[kk] = pow(10.0, -3.5);
            yH2SMix[kk] = 0.0;
            TH2SaqMix[kk] = 0.0;
            UseH2SgasMix[kk] = 0; // 假设是 int 类型
        }
        if (UseTPpHMix[kk] == 1)
        {
            if (TofpH[kk] == 0.0 && UseSI == 0)
                TofpH[kk] = 77.0;
            if (PofpH[kk] == 0.0 && UseSI == 0)
                PofpH[kk] = 14.696;
            if (TofpH[kk] == 0.0 && UseSI == 1)
                TofpH[kk] = 25.0;
            if (PofpH[kk] == 0.0 && UseSI == 1)
                PofpH[kk] = 1.0;
        }
        else
        {
            if (UseSI == 0)
            {
                TofpH[kk] = 77.0;
                PofpH[kk] = 14.696; // set to 77 F and 14.696 psia as default for pH calculation
            }
            else
            {
                TofpH[kk] = 25.0;
                PofpH[kk] = 1.0;
            }
        }
    }

    if (UseTPVolMix[kk] == 1)
    {
        if (TofVol[kk] == 0.0 && UseSI == 0)
            TofVol[kk] = 77.0; // set default to 77 F
        if (PofVol[kk] == 0.0 && UseSI == 0)
            PofVol[kk] = 14.696; // set default to 1 atm
        if (TofVol[kk] == 0.0 && UseSI == 1)
            TofVol[kk] = 25.0;
        if (PofVol[kk] == 0.0 && UseSI == 1)
            PofVol[kk] = 1.013254;
    }
    else
    {
        if (UseSI == 0)
        {
            TofVol[kk] = 77.0;
            PofVol[kk] = 14.696;
        }
        else
        { // set to 25 C and 1.013254 atm as default for EOS calculation
            TofVol[kk] = 25.0;
            PofVol[kk] = 1.013254;
        }
    }

    yCH4Mix[kk] = 1 - (yCO2Mix[kk] + yH2SMix[kk]);
    if (yCH4Mix[kk] < 0)
    {
        yCH4Mix[kk] = 0;
    }

    if (UseH2SgasMix[kk] != 0 && UseH2SgasMix[kk] != 1)
    {
        // MsgBox("Row 32 is expecting a value of 1 or 0, please enter a vlue of 1 or 0 and rerun")
        //     End
    }

    if (yCO2Mix[kk] > 1.0)
    {
        errmsg[1] = 1;
        yCO2Mix[kk] = 1.0;
    }

    if (yCO2Mix[kk] < 0.0)
    {
        errmsg[2] = 2;
        yCO2Mix[kk] = 0.0;
        yCO2 = 0.0;
        CO2aq = 0.0;
        HCO3 = 0.0;
        CO3 = 0.0;
    }

    if (yH2SMix[kk] > 1.0)
    {
        errmsg[3] = 3;
        yH2SMix[kk] = 1.0;
        yH2S = 1.0;
        TH2SaqMix[kk] = 0.0; // This will cause the program to use yH2Sstp as the calculation for TH2Saq instead of the input sheet value
    }

    if (yH2SMix[kk] < 0.0)
    {
        errmsg[4] = 4;
        yH2SMix[kk] = 0.0;
        yH2S = 0.0;
        HS = 0.0;
        H2Saq = 0.0;
        TH2SaqMix[kk] = 0.0;
    }
}


double VgTPWI[40][10];
double VoWI[40][10];
double VwWI[40][10];
double VwSW1;
double VoSW1;
double VgSW1;
double VMeOHSW1;
double VMEGSW1;
int RunMultiMixSlb;
double TpH;
double PpH;
double TVol;
double Pvol;
double CO2aqOld;
double TCO2SSE;
double TCO2Old;
// string mt;
double Patm;
double PBar;
double TK;

double SArea;
double QBrineFlow;
double radiusC[11];
double radiusA[12];
int iCa, iSr, iBa, iFe, iHCO3, iHS, iCO3, iSO4;




void ReadInputPartB(int kk,SampleData *data)
{
    VgTPMix[kk] = data->Q_Gas;
    VoMix[kk] = data->Q_Oil;
    VwMix[kk] = data->Q_Water;
    oilAPIgravMix[kk] = data->API;
    gasSpGravMix[kk] = data->SG_g;
    VMeOHMix[kk] = data->Q_MeOH;
    VMEGMix[kk] = data->Q_MEG;

    if (RunWhatIf == 1)
    {
        if (CaseCountWI[Loop1WI] == 25)
        {
            VgTPMix[kk] = VgTPWI[Loop1WI][Loop2WI];
        }
        if (CaseCountWI[Loop1WI] == 26)
        {
            VoMix[kk] = VoWI[Loop1WI][Loop2WI];
        }
        if (CaseCountWI[Loop1WI] == 27)
        {
            VwMix[kk] = VwWI[Loop1WI][Loop2WI];
        }
    }

    if (Run_Seawater_Mixing == 1 && LoopMixing == 1 && kk == 1)
    {
        VwSW1 = VwMix[1];
        VgSW1 = VgTPMix[1];
        VoSW1 = VoMix[1];
        VMeOHSW1 = VMeOHMix[1];
        VMEGSW1 = VMEGMix[1];
    }

    // 舍弃此字段
    if (Run10TestCases == 1)
    {
        usepHmix[kk] = 0;
        // yCO2Mix(kk) = Worksheets("Calcite").Cells(15 + Loop10, 3) / 100
    }

    if (RunMultiMix == 1)
    {
        // VgTPMix(kk) = Worksheets("MultiMix").Cells(2 + LoopResChem, 14 + j).Value
        // VoMix(kk) = Worksheets("MultiMix").Cells(2 + LoopResChem, 9 + j).Value
        // VwMix(kk) = Worksheets("MultiMix").Cells(2 + LoopResChem, 4 + j).Value
    }
    if (RunMultiMixSlb == 1)
    {
        // VgTPMix[kk] = Worksheets("MultiMix_Slb").Cells(2 + LoopResChem, 14 + j).Value;
        // VoMix[kk] = Worksheets("MultiMix_Slb").Cells(2 + LoopResChem, 9 + j).Value;
        // VwMix[kk] = Worksheets("MultiMix_Slb").Cells(2 + LoopResChem, 4 + j).Value
    }
    if (Run_MixingTwoWells == 1)
    {
        VgTPMix[kk] = VgTPMix[kk] * MixFrac[kk];
        VoMix[kk] = VoMix[kk] * MixFrac[kk];
        VwMix[kk] = VwMix[kk] * MixFrac[kk];
        VMeOHMix[kk] = VMeOHMix[kk] * MixFrac[kk];
        VMEGMix[kk] = VMEGMix[kk] * MixFrac[kk];
    }
    if (Run_Seawater_Mixing == 1)
    {
        VwMix[kk] = VwSW1 * MixFrac[kk];

        if (kk == 1)
        {
            VoMix[kk] = VoSW1 * MixFrac[kk];
            VgTPMix[kk] = VgSW1 * MixFrac[kk];
            VMeOHMix[kk] = VMeOHSW1 * MixFrac[kk];
            VMEGMix[kk] = VMEGSW1 * MixFrac[kk];
        }
        else
        {
            VoMix[kk] = 0;
            VgTPMix[kk] = 0;
            VMeOHMix[kk] = 0;
            VMEGMix[kk] = 0;
        }
    }

    if (UseSI == 0)
    {
        // convert 1000 ft^3 to m^3
        VgTPMix[kk] = VgTPMix[kk] * 28.31685;
    }
    else
    {
        // convert 1000 m^3 to m^3
        VgTPMix[kk] = VgTPMix[kk] * 1000;
    }

    if (UseSI == 1)
    {
        // Convert m^3 to barrels
        VoMix[kk] = VoMix[kk] / 0.159;
        VwMix[kk] = VwMix[kk] / 0.159;
        VMeOHMix[kk] = VMeOHMix[kk] / 0.159;
        VMEGMix[kk] = VMEGMix[kk] / 0.159;

        // Convert temperature from Celsius to Fahrenheit
        TofVol[kk] = TofVol[kk] * 9.0 / 5.0 + 32.0;

        // Convert pressure from Bar to psia
        PofVol[kk] = PofVol[kk] * 14.503774;

        // Same conversions for pH-related data
        TofpH[kk] = TofpH[kk] * 9.0 / 5.0 + 32.0;
        PofpH[kk] = PofpH[kk] * 14.503774;
    }

    TpH = TofpH[kk];
    PpH = PofpH[kk];
    TVol = TofVol[kk];
    Pvol = PofVol[kk];

    // 若为 0 或空字符串则赋默认值
    // 在 C 中不能判断 ""，只能判断是否为 0
    if (oilAPIgravMix[kk] == 0)
    {
        oilAPIgravMix[kk] = 30.0;
    }
    if (gasSpGravMix[kk] == 0)
    {
        gasSpGravMix[kk] = 0.6;
    }

    if (UseTPVolMix[kk] == 0)
    {
        // mt = fTPFunc(0); // STP condition

        GasDensityMix[kk] = gasSpGravMix[kk] *
            (Patm * 28.97 / (0.08206 * TK)); // kg/m^3

        // 调用外部函数，暂时忽略fH2ODensity
        //OilDensityMix[kk] = (141.5 / (oilAPIgravMix[kk] + 131.5)) *
        //    (fH2ODensity(TK, PBar) / 1000.0);
    }
    else
    { // UseTPVolMix == 1

        // mt = fTPFunc(1); // Use actual T, P

        GasDensityMix[kk] = gasSpGravMix[kk] *
            (Patm * 28.97 / (0.08206 * TK)); // kg/m^3
        // 调用外部函数，暂时忽略fH2ODensity
        //OilDensityMix[kk] = (141.5 / (oilAPIgravMix[kk] + 131.5)) *
        //    (fH2ODensity(TK, PBar) / 1000.0);
    }

    // If volume of gas is equal to zero, then add one mL of gas
    if (VgTPMix[kk] == 0.0)
    {
        VgTPMix[kk] = 1.0 / 1000000.0;
    }

    // If no oil is reported, add one ml of oil per day, to avoid singularities
    if (VoMix[kk] == 0.0)
    {
        VoMix[kk] = 1.0 / 159.0 / 1000.0;
    }

    // If no water is reported, add one ml of water per day, to avoid singularities
    if (VwMix[kk] == 0.0)
    {
        VwMix[kk] = 1.0 / 159.0 / 1000.0;
    }

    // MeOH and MEG cannot both be > 0
    if (VMeOHMix[kk] > 0.0 && VMEGMix[kk] > 0.0)
    {
        // printf("MeOH and MEG cannot both be > 0. Please check the Input Sheet and try it again.\n");
        // exit(1); // equivalent to VB's 'End'
        return;
    }

    // Mass of methanol (Kg)
    mass_MeOH_mix[kk] = 159.0 * VMeOHMix[kk] * 0.7914;

    // Mass of MEG (Kg)
    mass_MEG_mix[kk] = 159.0 * VMEGMix[kk] * 1.1098;

    // Brine flow rate, unit: cm3/sec
    QBrineFlow = VwMix[kk] * 159.0 * 1000.0 / 86400.0;

    // Calculate Surface Area contacted by brine, cm²
    // SArea = M_PI * PipeID * PipeL;  无M_PI参数

    // Define hydrated cationic radii (Angstroms)
    radiusC[iCa] = 4.12;
    radiusC[iFe] = 4.28;
    radiusC[iBa] = 4.04;
    radiusC[iSr] = 4.12;

    // Define hydrated anion radii (Angstroms)
    radiusA[iHCO3] = 3.04;
    radiusA[iHS] = 2.07;
    radiusA[iCO3] = 3.94;
    radiusA[iSO4] = 3.79;
}



double gNeut[16]; double zOutput[15]; double z[20]; double density[3]; double mc[15]; double ma[15];
double ChCat[15]; double ChAn[15];

//全局变量
/*
   iH = 1: iNa = 2: iK = 3: iMg = 4: iCa = 5: iSr = 6: iBa = 7: iFe = 8: iZn = 9: iPb = 10: iNH4 = 11: iRa = 12  注：改成偏移

   iOH = 0, iCl = 1, iAc = 2, iHCO3 = 3, iCO3 = 4, iSO4 = 5, iHS = 6, intF = 7, iBr = 8,  iH2BO3 = 9, iH3SiO4 = 10, iH2SiO4 = 11 'Anion indexes

   iCH4aq = 1: iCO2aq = 2: iH2Saq = 3: iHAcaq = 4: iH4SiO4aq = 5: iNH3 = 6: iH3BO3 = 7: iFeSaq = 8 'Neutral aquatic indexes

*/
int iH = 0, iNa = 1, iK = 2, iMg = 3, iCa = 4, iSr = 5, iBa = 6, iFe = 7, iZn = 8, iPb = 9, iNH4 = 10, iRa = 11;
int iOH = 0, iCl = 1, iAc = 2, iHCO3 = 3, iCO3 = 4, iSO4 = 5, iHS = 6, intF = 7, iBr = 8, iH2BO3 = 9, iH3SiO4 = 10, iH2SiO4 = 11;
int iCH4aq = 0, iCO2aq = 1, iH2Saq = 2, iHAcaq = 3, iH4SiO4aq = 4, iNH3 = 5, iH3BO3 = 6, iFeSaq = 7;


double Alk, TAc, TNH4, TH3BO3, TH2Saq, TCO2, VW, VgTP, VO, VMeOH, VMEG, mass_MeOH, mass_MEG;
double yCO2, yH2S, yCH4;
double TFe; double SGG; double API;
double Ppsia, TF, TC;

double mtotal, MoleCharge; // 仅在C2_PitzerActCoefs_T_P_ISt、CalcIonicStrength使用

double SumOfCations, SumOfAnions; // 仅在CalcIonicStrength、QualityControlCalculations使用

double RatioOilBPoints;

double Ist; double DpHj;
double pH; 

int NumCat = 12, NumAn = 13;//final不变量

// 函数定义
void fTPFunc(int iTP) {
    if (iTP == 0) {
        Ppsia = 14.696;
        Patm = Ppsia / 14.696;
        PBar = Ppsia / 14.503774;
        TF = 77;
        TC = (TF - 32) * 5.0 / 9.0;
        TK = TC + 273.15;
    }

    if (iTP == 1) {
        Ppsia = Pvol;
        Patm = Ppsia / 14.696;
        PBar = Ppsia / 14.503774;
        TF = TVol;
        TC = (TF - 32) * 5.0 / 9.0;
        TK = TC + 273.15;
    }

    if (iTP == 2) {
        Ppsia = PpH;
        Patm = Ppsia / 14.696;
        PBar = Ppsia / 14.503774;
        TF = TpH;
        TC = (TF - 32) * 5.0 / 9.0;
        TK = TC + 273.15;
    }

    if (iTP == 3) { // for API
        Ppsia = 14.696;
        Patm = Ppsia / 14.696;
        PBar = Ppsia / 14.503774;
        TF = 60;
        TC = (TF - 32) * 5.0 / 9.0;
        TK = TC + 273.15;
    }
}

void CalcIonicStrength(void)
{
    // Calcule common system terms and functions for Pitzer theory: 
    // ISt, Z, m-total gX, gpX, JX, and JpX.

    mtotal = 0;
    MoleCharge = 0.0;
    Ist = 0.0;
    SumOfCations = 0.0;

    for (int c = 0; c < NumCat; c++) {
        MoleCharge = MoleCharge + ChCat[c] * mc[c];
        SumOfCations = SumOfCations + ChCat[c] * mc[c];
        Ist = Ist + pow(ChCat[c], 2.0) * mc[c];
        mtotal = mtotal + mc[c];
    }

    // SumOfAnions = -0.0000001
    SumOfAnions = 0.0;

    for (int a = 0; a < NumAn; a++) {
        MoleCharge = MoleCharge + fabs(ChAn[a]) * ma[a];
        SumOfAnions = SumOfAnions + ChAn[a] * ma[a];
        Ist = Ist + pow(ChAn[a], 2.0) * ma[a];
        mtotal = mtotal + ma[a];
    }

    if (Ist <= 0.0) {
        // When only negative alkalinity is entered in the input sheet,
        // set ionic molality = (H+OH) of pure water
        Ist = 2.0 * 0.0000001;
    }

    Ist = Ist / 2.0;
    DpHj = 0.129 * pow(Ist, 0.5);
}

double fRatioOilBPoints() {
    double base = SGG / 0.6;
    double rankine = TF + 460.0; // 转换为兰金温度
    double api_diff = API - 30.0;

    if (API >= 30.0) {
        return pow(base * exp(23.93 * api_diff / rankine), 1.0 / 1.187);
    }
    else {
        return pow(base * exp(25.724 * api_diff / rankine), 1.0 / 1.0937);
    }
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


double KgwCO2, K1H2CO3, KHAc, K2HCO3, KH2O, K1H2S, K2HS, KNH4, KH3BO3;
double KgoCO2, KgoH2S;  //A13_CalcH2SFugacity、C1_ThermodynamicEquilConsts、C4_SSPEquilCalcs、fTotalCO2H2Smoles
double KgwCH4, KgoCH4; // C1_ThermodynamicEquilConsts、C4_SSPEquilCalcs、fTotalCO2H2Smoles
double KgwH2S;         //C1_ThermodynamicEquilConsts、C4_SSPEquilCalcs、C4_EOS_TCO2_SSPEquilCalcs、C5_CalcpHPCO2PH2SSTP、QualityControlCalculations、fTotalCO2H2Smoles
double RBar = 0.083144; double R = 83.144;

// ==================== C1_ThermodynamicEquilConsts ====================
//  This subroutine is called by subroutines D, E, and G above.
// =====================================================================
void C1_ThermodynamicEquilConsts(void)
{
    TC = TK - 273.15;
    double Psat = PsatH2O(TK);

    KgwCO2 = 55.508 / exp(
        log(Psat)
        - 9.14122 / (TK / 647.096)
        + 2.8192 * pow(1 - TK / 647.096, 0.355) / (TK / 647.096)
        + 11.28516 * pow(TK / 647.096, -0.41) * exp(1 - TK / 647.096)
        - 0.8066
    ); // IUPAC 2012 KH, unit m/bar

    double dV = (37.88 - 0.14442 * TC + 0.001243 * pow(TC, 2)
        - 0.0000044738 * pow(TC, 3) + 0.000000005726 * pow(TC, 4)) * 0.001;

    if (TK <= 373.15)
        KgwCO2 = KgwCO2 * pow(exp(dV * (PBar - 1) / RBar / TK
            + 0.0000067 / (2 * RBar * TK) * pow(PBar - Psat, 2)), -1) / 14.503774;

    if (TK > 373.15)
        KgwCO2 = KgwCO2 * pow(exp(dV * (PBar - 1) / RBar / TK
            + 0.0000067 / (2 * RBar * TK) * pow(PBar - Psat, 2)), -1) / 14.503774;

    if (TC <= 100)
        Psat = 1.013254; // set 1 atm when < 100°C

    // ================= Carbonic acid system (IUPAC 2012) =================
    double q1 = -441.490479, q2 = 26901.0527, q3 = 157.2016907, q4 = -0.07219967, q5 = -2003878.4, q6 = -19.57801521, 
        q7 = 925.6200149, q8 = 6.714256299, q9 = 0.003645431058, q10 = -0.1743884044, q11 = -0.00124018735;

    K1H2CO3 = pow(10,
        q1 + q2 / TK + q3 * log10(TK) + q4 * TK + q5 / pow(TK, 2)
        + (PBar - Psat) * (q6 / TK + q7 / pow(TK, 2) + q8 * log10(TK) / TK)
        + pow(PBar - Psat, 2) * (q9 / TK + q10 / pow(TK, 2) + q11 * log10(TK) / TK)
    );

    q1 = -332.5306; q2 = 17540.07; q3 = 120.13393;
    q4 = -0.06545969; q5 = -1277752.3;
    q6 = -12.81797624; q7 = 603.2417035; q8 = 4.419625804;
    q9 = 0.00139842542; q10 = -0.07141847943; q11 = -0.0004736672395;

    K2HCO3 = pow(10,
        q1 + q2 / TK + q3 * log10(TK) + q4 * TK + q5 / pow(TK, 2)
        + (PBar - Psat) * (q6 / TK + q7 / pow(TK, 2) + q8 * log10(TK) / TK)
        + pow(PBar - Psat, 2) * (q9 / TK + q10 / pow(TK, 2) + q11 * log10(TK) / TK)
    );

    // Acetic acid
    KHAc = pow(10, -1 * (-66.227 + 3216.269 / TK + 10.566 * log(TK)))
        * exp(-(-15.82 - 0.0219 * TC) * (Patm - 1) / (R * TK));

    // Water density (Wagner & Kruse, 1998)
    double RhoH2OTP = (999.83952 + 16.945176 * TC
        - 0.0079870401 * pow(TC, 2)
        - 0.000046170461 * pow(TC, 3)
        + 0.00000010556302 * pow(TC, 4)
        - 2.8054253E-10 * pow(TC, 5))
        / (1 + 0.01687985 * TC);

    RhoH2OTP += (0.043922 - 0.000076 * TC + 0.000000126 * pow(TC, 2)
        + 0.00000000319 * pow(TC, 3)) * (Patm - 1);

    KH2O = pow(10,
        -4.098 - 3245.2 / TK + 223620 / pow(TK, 2) - 39840000 / pow(TK, 3)
        + (13.957 - 1262.3 / TK + 856410 / pow(TK, 2)) * log10(RhoH2OTP / 1000)
    );

    // ================= H2S system =================
    q1 = 782.43945; q2 = 0.361261; q3 = -0.00016722;
    q4 = -20565.7315; q5 = -142.741722;

    K1H2S = pow(10,
        q1 + q2 * TK + q3 * pow(TK, 2) + q4 / TK + q5 * log(TK))
        * exp(-(-14.8 + 0.002 * TC - 0.0004 * pow(TC, 2)) * (Patm - Psat) / (R * TK));
    // ================= H2S Second dissociation =================
    q1 = 137.2755; q2 = -8235.4184; q3 = -22.5503;
    q4 = 0.0162; q6 = -12.81797624;
    q7 = 603.2417035; q8 = 4.419625804;
    q9 = 0.00139842542; q10 = -0.07141847943; q11 = -0.0004736672395;

    K2HS = pow(10,
        q1 + q2 / TK + q3 * log(TK) + q4 * TK
        + (PBar - Psat) * (q6 / TK + q7 / pow(TK, 2) + q8 * log10(TK) / TK)
        + pow(PBar - Psat, 2) * (q9 / TK + q10 / pow(TK, 2) + q11 * log10(TK) / TK)
    );

    // ================= Ammonium ion =================
    KNH4 = pow(10, -16.79744216 - 1893.658478 / TK + 5.614302869 * log10(TK))
        * exp(((-26.43 + 0.0889 * TC - 0.000905 * pow(TC, 2)) / (R * TK)) * (PBar - Psat))
        * exp(0.5 * 0.001 * ((-5.03 + 0.0814 * TC) / (R * TK)) * pow(PBar - Psat, 2));

    // ================= Boric acid =================
    KH3BO3 = pow(10,
        50.0370439 - 3283.059242 / TK - 19.50251187 * log10(TK))
        * exp(((-29.48 - 0.1622 * TC + 0.002608 * pow(TC, 2)) / (R * TK)) * (PBar - Psat))
        * exp(0.5 * 0.001 * ((-2.84) / (R * TK)) * pow(PBar - Psat, 2));

    // ================= Henry's law constants (molality/psia) =================
    KgoCO2 = pow(10, -1 * (2.198 + 0.00253 * TF - 0.000001681 * pow(TF, 2) + 0.000038 * Ppsia));
    KgwCH4 = pow(10, -1 * (3.798 + 0.004136 * TF - 0.000009486 * pow(TF, 2) + 0.000038 * Ppsia));
    KgoCH4 = pow(10, -1 * (2.728 + 0.002546 * TF - 0.000004265 * pow(TF, 2) + 0.000038 * Ppsia));
    KgoH2S = pow(10, -1 * (1.57 + 0.005043 * TF - 0.000005311 * pow(TF, 2) + 0.000038 * Ppsia));

    // ================= Henry's constant for H2S (Kassa) =================
    double ln_H_T_H2S = (13788 / TK - 185.19 + 29.087 * log(TK)
        - 0.027637 * TK - 1445200 / pow(TK, 2));

    dV = (33.18 + 0.092661 * TC - 0.00054853 * pow(TC, 2)
        + 0.0000015354 * pow(TC, 3) - 0.0000000015459 * pow(TC, 4)) * 0.001;

    KgwH2S = pow(exp(ln_H_T_H2S) * exp(dV * (Patm * 1.01353 - Psat) / RBar / TK), -1)
        / 14.503774;

    // ================= Carbonate minerals =================
    if ((PBar - Psat) <= 500)
        KspCalcite = pow(10, -171.9065 - 0.077993 * TK + 2839.319 / TK
            + 71.595 * log10(TK))
        * pow(10, (0.514 - 0.000197 * TC + 1.7096E-15 * pow(TC, 6))
            * (PBar - Psat) / 500);

    if ((PBar - Psat) > 500)
        KspCalcite = pow(10, -171.9065 - 0.077993 * TK + 2839.319 / TK
            + 71.595 * log10(TK))
        * pow(10, 0.514 - 0.000197 * TC + 1.7096E-15 * pow(TC, 6))
        * pow(10, ((0.928 - 0.000079 * TC + 2.1485E-15 * pow(TC, 6))
            - (0.514 - 0.000197 * TC + 1.7096E-15 * pow(TC, 6)))
            * ((PBar - Psat) - 500) / 500);

    if ((PBar - Psat) <= 500)
        KspBarite = pow(10, 136.035 - 7680.41 / TK - 48.595 * log10(TK))
        * pow(10, (0.394 - 0.0001119 * TC + 1.5305E-15 * pow(TC, 6))
            * (PBar - Psat) / 500);

    if ((PBar - Psat) > 500)
        KspBarite = pow(10, 136.035 - 7680.41 / TK - 48.595 * log10(TK))
        * pow(10, 0.394 - 0.0001119 * TC + 1.5305E-15 * pow(TC, 6))
        * pow(10, ((0.674 + 0.0001229 * TC + 1.9202E-15 * pow(TC, 6))
            - (0.394 - 0.0001119 * TC + 1.5305E-15 * pow(TC, 6)))
            * ((PBar - Psat) - 500) / 500);

    // ================= Celestite =================
    KspCelestite = exp(-4762.71 - 0.878035 * TK + 0.000184788 * pow(TK, 2)
        - 320587.4 / TK + 731.756 * log(TK)
        + 99430.6 * log(TK) / TK
        - 2.9811 * (TK - 263) / TK * log(TK - 263));

    dV = -(337.2 + -1.754 * TK + 0.002658 * pow(TK, 2)) * 0.001;
    dk = (-388.1 + 2.26 * TK + -0.00338 * pow(TK, 2)) * 0.000001;

    if (TK <= 373.15)
        KspCelestite *= exp(-dV * (PBar - 1) / RBar / TK
            + dk / (2 * RBar * TK) * pow(PBar - 1, 2));
    else
        KspCelestite *= exp(-dV * (PBar - Psat) / RBar / TK
            + dk / (2 * RBar * TK) * pow(PBar - Psat, 2));

    // ================= Hemihydrate (CaSO4·0.5H2O) =================
    if ((PBar - Psat) <= 500)
        KspHemihydrate = pow(10, -1 * (-183.603 + 8033.771 / TK + 28.173 * log(TK)))
        * pow(10, (0.423 - 0.0001 * TC + 1.6176E-15 * pow(TC, 6))
            * (PBar - Psat) / 500);

    if ((PBar - Psat) > 500)
        KspHemihydrate = pow(10, -1 * (-183.603 + 8033.771 / TK + 28.173 * log(TK)))
        * pow(10, 0.423 - 0.0001 * TC + 1.6176E-15 * pow(TC, 6))
        * pow(10, ((0.729 + 0.0001576 * TC + 2.0302E-15 * pow(TC, 6))
            - (0.423 - 0.0001 * TC + 1.6176E-15 * pow(TC, 6)))
            * ((PBar - Psat) - 500) / 500);

    // ================= Gypsum (CaSO4·2H2O) =================
    if ((PBar - Psat) <= 500)
        KspGypsum = pow(10, -128.9622 + 1.761129 * TK - 0.01484172 * pow(TK, 2)
            + 0.000837554 * pow(TK, 2.5) - 0.00001384613 * pow(TK, 3))
        * pow(10, (0.328 + 0.0000589 * TC + 1.612E-15 * pow(TC, 6))
            * (PBar - Psat) / 500);

    if ((PBar - Psat) > 500)
        KspGypsum = pow(10, -128.9622 + 1.761129 * TK - 0.01484172 * pow(TK, 2)
            + 0.000837554 * pow(TK, 2.5) - 0.00001384613 * pow(TK, 3))
        * pow(10, (0.328 + 0.0000589 * TC + 1.612E-15 * pow(TC, 6)))
        * pow(10, ((0.553 + 0.0004785 * TC + 2.018E-15 * pow(TC, 6))
            - (0.328 + 0.0000589 * TC + 1.612E-15 * pow(TC, 6)))
            * ((PBar - Psat) - 500) / 500);

    // ================= Anhydrite =================
    KspAnhydrite = pow(10,
        (-44.9380794269878 + 1.08048295151073E-02 * PBar) * log10(TK)
        + (1.818855325 * PBar - 4943.71187) / TK
        - 3.19664522256897E-02 * PBar + 123.46568);

    // ================= Siderite (FeCO3) =================
    KspSiderite = exp(129.97 / TK - 50.205 + 7.3143 * log(TK) - 0.052913 * TK)
        * exp(-(-48.76 - 0.5304 * TC) * (Patm - Psat) / (R * TK));

    // ================= RaSO4 (same pressure dependence as Barite) =================
    if ((PBar - Psat) <= 500)
        KspRaSO4 = pow(10, 137.98 - 8346 / TK - 48.595 * log10(TK))
        * pow(10, (0.394 - 0.0001119 * TC + 1.5305E-15 * pow(TC, 6))
            * (PBar - Psat) / 500);

    if ((PBar - Psat) > 500)
        KspRaSO4 = pow(10, 137.98 - 8346 / TK - 48.595 * log10(TK))
        * pow(10, 0.394 - 0.0001119 * TC + 1.5305E-15 * pow(TC, 6))
        * pow(10, ((0.674 + 0.0001229 * TC + 1.9202E-15 * pow(TC, 6))
            - (0.394 - 0.0001119 * TC + 1.5305E-15 * pow(TC, 6)))
            * ((PBar - Psat) - 500) / 500);

    // ================= Halite (NaCl) =================
    KspHalite = pow(10, 462.55683140897 + 0.225843153104659 * TK
        - 11062.1698808439 / TK
        - 84.5268944091097 * log(TK)
        - 1.08055601606581E-04 * pow(TK, 2));

    dV = (-14.156 + 0.22255 * TC - 0.0031441 * pow(TC, 2)
        + 0.00001744 * pow(TC, 3) - 0.000000045213 * pow(TC, 4)) * 0.001;

    dk = -1.183398E-14 * pow(TK, 4)
        + 2.019932E-11 * pow(TK, 3)
        - 0.00000001332746 * pow(TK, 2)
        + 0.000003882627 * TK - 0.0004193603;

    KspHalite *= exp(-dV * (PBar - Psat) / RBar / TK
        + dk / (2 * RBar * TK) * pow(PBar - Psat, 2));

    // ================= FeS (mackinawite / troilite series) =================
    KspFeS = pow(10, -1 * (-129.964067035986 + 6151.50521920783 / TK
        + 19.7824266766714 * log(TK)))
        * pow(10, (0.2 - 0.00005 * TC + 1.044E-15 * pow(TC, 6)) * PBar / 500);

    KspZnS = pow(10, -1 * (-138.1059 + 8665.5519 / TK + 21.0391 * log(TK)))
        * pow(10, (0.198 + 0.0000094 * TC + 1.049E-15 * pow(TC, 6)) * PBar / 500);

    KspPbS = pow(10, -1 * (-113.4984 + 9436.5178 / TK + 16.9658 * log(TK)))
        * pow(10, (0.202 - 0.0000663 * TC + 9.142E-16 * pow(TC, 6)) * PBar / 500);

    KspFeSAm = pow(10, -1 * (-14.12 + 2731.3 / TK + 0.02654 * TK))
        * pow(10, (0.2 - 0.00005 * TC + 1.044E-15 * pow(TC, 6)) * PBar / 500);

    KspTrot = pow(10, -101.92925047219 + 2310.1324958895 / TK
        + 17.89882023677 * log(TK)
        - 0.039206108448706 * TK)
        * pow(10, (0.2 - 0.00005 * TC + 1.044E-15 * pow(TC, 6)) * PBar / 500);

    KspPyrr = pow(10, -108.04725226365 + 2391.1732648763 / TK
        + 19.064810519867 * log(TK)
        - 0.0415135620025018 * TK)
        * pow(10, (0.2 - 0.00005 * TC + 1.044E-15 * pow(TC, 6)) * PBar / 500);

    // ================= Fluorite (CaF2) =================
    KspCaF2 = pow(10, 66.348 - 4298.2 / TK - 25.271 * log10(TK))
        * pow(10, (0.399 - 0.0000047 * TC) * Patm / 500);

    // ================= Dolomite (CaMg(CO3)2) =================
    KspDol = pow(10, -1 * (257.2181 + 0.1438 * TK - 3197.3234 / TK - 109.6227 * log10(TK)))
        * pow(10, (0.982 - 0.0003101 * TC) * Patm / 500);

    KspDol = pow(10, -1 * (-175.9 + 6835.4 / TK + 68.727 * log10(TK)))
        * pow(10, (0.982 - 0.0003101 * TC) * Patm / 500);

    // ================= Cosolvent dielectric constant =================
    DielecConst = 1 / (0.0068 * pow(xMeOH, 2) + 0.0115 * xMeOH + 0.0128);

    KstCaHCO3 = 0.034877
        * exp(139190.0 / DielecConst / TK)
        * exp(-5.1328E-38
            * pow((2.528E+23 * IStCosolvent / DielecConst / TK), 0.5)
            / (1.5365E-33 * DielecConst * TK
                * (1 + 0.00000000024
                    * pow((2.528E+23 * IStCosolvent / DielecConst / TK), 0.5))));

    // ================= Smithsonite (ZnCO3) =================
    KspZnCO3 = pow(10, -1 * (8.7334 - 2173.9249 / TK + 1.4880105 * log(TK)))
        * pow(10, (0.468 - 0.000176 * TC) * Patm / 500);

    // ================= Ca(OH)2 & Mg(OH)2 =================
    KspCaOH2 = pow(10, 110.3408 - 4516.4997 / TK - 40.6657 * log10(TK))
        * pow(10, (0.528 - 0.0002789 * TC) * Patm / 500);

    KspMgOH2 = pow(10, 92.6256 - 4665.57436 / TK - 15.52458654 * log(TK))
        * pow(10, (0.484 - 0.0001633 * TC) * Patm / 500);

    // ================= SrCO3 and BaCO3 =================
    if ((PBar - Psat) <= 500)
        KspSrCO3 = pow(10, 155.6841 - 7272.6012 / TK - 56.8052 * log10(TK))
        * pow(10, (0.528 - 0.000259 * TC + 1.682E-15 * pow(TC, 6))
            * (PBar - Psat) / 500);

    if ((PBar - Psat) > 500)
        KspSrCO3 = pow(10, 155.6841 - 7272.6012 / TK - 56.8052 * log10(TK))
        * pow(10, (0.528 - 0.000259 * TC + 1.682E-15 * pow(TC, 6))
            * pow(10, ((0.951 - 0.0001671 * TC + 2.114E-15 * pow(TC, 6))
                - (0.528 - 0.000259 * TC + 1.682E-15 * pow(TC, 6)))
                * ((PBar - Psat) - 500) / 500));

    if ((PBar - Psat) <= 500)
        KspBaCO3 = pow(10, 244.3819 - 11526.0874 / TK - 86.6577 * log10(TK))
        * pow(10, (0.523 - 0.0003039 * TC + 1.631E-15 * pow(TC, 6))
            * (PBar - Psat) / 500);

    if ((PBar - Psat) > 500)
        KspBaCO3 = pow(10, 244.3819 - 11526.0874 / TK - 86.6577 * log10(TK))
        * pow(10, (0.523 - 0.0003039 * TC + 1.631E-15 * pow(TC, 6)))
        * pow(10, ((0.936 - 0.0002343 * TC + 2.05E-15 * pow(TC, 6))
            - (0.523 - 0.0003039 * TC + 1.631E-15 * pow(TC, 6)))
            * ((PBar - Psat) - 500) / 500);

    // ================= Silica equilibria (KH4SiO4, KH3SiO3) =================
    q1 = -69.27744384; q2 = -1.893100838; q3 = 11.76344126;
    q4 = -0.025416705; q5 = 0.0000102738; q6 = -19.57801521;
    q7 = 925.6200149; q8 = 6.714256299; q9 = 0.003645431058;
    q10 = -0.1743884044; q11 = -0.00124018735;

    if (TK <= 373.15)
        KH4SiO4 = pow(10, q1 + q2 / TK + q3 * log(TK) + q4 * TK + q5 / pow(TK, 2)
            + (PBar - 1) * (q6 / TK + q7 / pow(TK, 2) + q8 * log10(TK) / TK)
            + pow(PBar - 1, 2) * (q9 / TK + q10 / pow(TK, 2)
                + q11 * log10(TK) / TK));

    else
        KH4SiO4 = pow(10, q1 + q2 / TK + q3 * log(TK) + q4 * TK + q5 / pow(TK, 2)
            + (PBar - Psat) * (q6 / TK + q7 / pow(TK, 2)
                + q8 * log10(TK) / TK)
            + pow(PBar - Psat, 2) * (q9 / TK + q10 / pow(TK, 2)
                + q11 * log10(TK) / TK));

    q1 = -104.8656231; q2 = 1.44107352; q3 = 18.76152428;
    q4 = -0.046106584; q5 = -0.0000299188; q6 = -12.81797624;
    q7 = 603.2417035; q8 = 4.419625804; q9 = 0.00139842542;
    q10 = -0.07141847943; q11 = -0.0004736672395;

    if (TK <= 373.15)
        KH3SiO3 = pow(10, q1 + q2 / TK + q3 * log(TK) + q4 * TK + q5 / pow(TK, 2)
            + (PBar - 1) * (q6 / TK + q7 / pow(TK, 2)
                + q8 * log10(TK) / TK)
            + pow(PBar - 1, 2) * (q9 / TK + q10 / pow(TK, 2)
                + q11 * log10(TK) / TK));

    else
        KH3SiO3 = pow(10, q1 + q2 / TK + q3 * log(TK) + q4 * TK + q5 / pow(TK, 2)
            + (PBar - Psat) * (q6 / TK + q7 / pow(TK, 2)
                + q8 * log10(TK) / TK)
            + pow(PBar - Psat, 2) * (q9 / TK + q10 / pow(TK, 2)
                + q11 * log10(TK) / TK));

    // ================= Chrysotile (Mg3Si2O5(OH)4) =================
    q1 = 417.2882131; q2 = -26905.36959; q3 = -36.96412567;
    q4 = 0.030495614; q5 = 0.0000073516; q6 = -19.57801521;
    q7 = 925.6200149; q8 = 6.714256299; q9 = 0.003645431058;
    q10 = -0.1743884044; q11 = -0.00124018735;

    if (TK <= 373.15)
        KspChrysotile = pow(10, q1 + q2 / TK + q3 * log(TK)
            + q4 * TK + q5 / pow(TK, 2)
            + (PBar - 1) * (q6 / TK + q7 / pow(TK, 2)
                + q8 * log10(TK) / TK)
            + pow(PBar - 1, 2) * (q9 / TK + q10 / pow(TK, 2)
                + q11 * log10(TK) / TK));
    else
        KspChrysotile = pow(10, q1 + q2 / TK + q3 * log(TK)
            + q4 * TK + q5 / pow(TK, 2)
            + (PBar - Psat) * (q6 / TK + q7 / pow(TK, 2)
                + q8 * log10(TK) / TK)
            + pow(PBar - Psat, 2) * (q9 / TK + q10 / pow(TK, 2)
                + q11 * log10(TK) / TK));

    // ================= Diopside (CaMgSi2O6) =================
    q1 = 395.6671951; q2 = -25105.31181; q3 = -33.08370906;
    q4 = 0.02733812; q5 = 0.0000065452; q6 = -19.57801521;
    q7 = 925.6200149; q8 = 6.714256299; q9 = 0.003645431058;
    q10 = -0.1743884044; q11 = -0.00124018735;

    if (TK <= 373.15)
        KspDiopside = pow(10, q1 + q2 / TK + q3 * log(TK)
            + q4 * TK + q5 / pow(TK, 2)
            + (PBar - 1) * (q6 / TK + q7 / pow(TK, 2)
                + q8 * log10(TK) / TK)
            + pow(PBar - 1, 2) * (q9 / TK + q10 / pow(TK, 2)
                + q11 * log10(TK) / TK));
    else
        KspDiopside = pow(10, q1 + q2 / TK + q3 * log(TK)
            + q4 * TK + q5 / pow(TK, 2)
            + (PBar - Psat) * (q6 / TK + q7 / pow(TK, 2)
                + q8 * log10(TK) / TK)
            + pow(PBar - Psat, 2) * (q9 / TK + q10 / pow(TK, 2)
                + q11 * log10(TK) / TK));

    // ================= Greenalite (Fe3Si2O5(OH)4) =================
    q1 = 374.226; q2 = -25313.92; q3 = -30.77; q4 = 0.0225; q5 = 0.0000062;
    q6 = -19.57801521; q7 = 925.6200149; q8 = 6.714256299;
    q9 = 0.003645431058; q10 = -0.1743884044; q11 = -0.00124018735;

    if (TK <= 373.15)
        KspGreenalite = pow(10, q1 + q2 / TK + q3 * log(TK)
            + q4 * TK + q5 / pow(TK, 2)
            + (PBar - 1) * (q6 / TK + q7 / pow(TK, 2)
                + q8 * log10(TK) / TK)
            + pow(PBar - 1, 2) * (q9 / TK + q10 / pow(TK, 2)
                + q11 * log10(TK) / TK));
    else
        KspGreenalite = pow(10, q1 + q2 / TK + q3 * log(TK)
            + q4 * TK + q5 / pow(TK, 2)
            + (PBar - Psat) * (q6 / TK + q7 / pow(TK, 2)
                + q8 * log10(TK) / TK)
            + pow(PBar - Psat, 2) * (q9 / TK + q10 / pow(TK, 2)
                + q11 * log10(TK) / TK));

    // ================= Quartz =================
    KspQuartz = pow(10, -1 * (170.73 - 0.0174 * TK + 10.0 * log10(TK)))
        * pow(10, (0.254 - 0.00017 * TC) * (PBar - Psat) / 500);

    // ================= Amorphous Silica =================
    KspAmSilica = pow(10, -1 * (130.3 - 0.0004 * TK + 9.0 * log10(TK)))
        * pow(10, (0.230 - 0.00016 * TC) * (PBar - Psat) / 500);

    // ================= Fe(OH)3(am) =================
    KspFeOH3Am = pow(10, -1 * (-84.555 + 7804.43 / TK + 22.87 * log(TK)))
        * pow(10, (0.4 - 0.0001 * TC + 1.0E-15 * pow(TC, 6))
            * (PBar - Psat) / 500);

    // ================= Zn(OH)2(am) =================
    KspZnOH2Am = pow(10, -1 * (-96.157 + 9246.8 / TK + 25.74 * log(TK)))
        * pow(10, (0.4 - 0.0001 * TC + 1.0E-15 * pow(TC, 6))
            * (PBar - Psat) / 500);

    // ================= Metal chloride complexes =================
    BetaDot_ZnCl = pow(10, -110.54 + 0.0359 * TK + 3611.3 / TK + 16.132 * log10(TK));
    BetaDot_FeCl = pow(10, -86.84 + 0.0281 * TK + 3307.6 / TK + 13.56 * log10(TK));
    BetaDot_FeCl2 = pow(10, -92.41 + 0.0313 * TK + 3422.8 / TK + 14.67 * log10(TK));
    BetaDot_FeCl3 = pow(10, -101.09 + 0.0362 * TK + 3549.9 / TK + 15.88 * log10(TK));

    // ================= Miscellaneous =================
    KspAlOH3 = pow(10, -1 * (-115.116 + 9827.943 / TK + 30.343 * log(TK)))
        * pow(10, (0.4 - 0.0001 * TC + 1.0E-15 * pow(TC, 6))
            * (PBar - Psat) / 500);

    KspMgCO3 = pow(10, -1 * (145.291 - 7153.2 / TK - 52.657 * log10(TK)))
        * pow(10, (0.51 - 0.0002 * TC + 1.6E-15 * pow(TC, 6))
            * (PBar - Psat) / 500);

    // 结束函数
}

}



void ReadInputPartD(int kk, int j, SampleData* data) 
{

    int mol_w3;

    double feed_Composition[15];
    for (int i = 0; i < 15; i++) {
        feed_Composition[i] = 0;
    }

    double mw_oil = 0;

    double tempgNeut[2] = { gNeut[iCO2aq], gNeut[iH2Saq] };

    //-----------------------------------------------------
    int iNG = 0;
    for (; iNG < 15; iNG++) {
        zOutput[iNG] = 0;
        z[iNG] = 0;
    }
    for (; iNG < 2; iNG++)
        density[iNG] = 0;

    if (RunStat == 0) {
        //将 A2 - A4 和 A6 - A8 的闪光灯设置为关闭
        if (Run10TestCases == 1 || Run_MassTransfer == 1 || Run_Seawater_Mixing == 1 || Run_MixingTwoWells == 1 || RunWhatIf == 1 || RunMultiMix == 1) {
            if (nob_Input + nob_InputII == 1) {
                //待定：Worksheets(mySheet).Cells(55, j + 2).Value = 0
            }
            else {
                //待定：Worksheets(MySheetMix).Cells(55, 8).Value = Empty
            }
        }
        if (nob_Input + nob_InputII == 1) {
            //待定：useEOSmix(kk) = Worksheets(mySheet).Cells(55, j + 2).Value
        }
        else {
            //待定：useEOSmix(kk) = Worksheets(MySheetMix).Cells(55, 8).Value 'UseEOSMix  values for mixed solution IS set on Column H
        }

        //待定：usedryHC = Worksheets(mySheet).Cells(56, j + 2).Value

        if (RunShellMultiflash == 1) useEOSmix[kk] = 0;

        /*
        * 待定，注：useEOSmix是整形的，为什么会用""来判断？？？？？
        if (useEOSmix[kk] == "") useEOSmix(kk) = 0
        */
        if (RunMultiMix == 1) useEOSmix[kk] = 0;

        //注意，由于 ZMix (kk,i) 在下面重新分配，因此值将从此 Sub 开头的输入表中重新读取。
        for (iNG = 0; iNG < 14; iNG++) {
            //待定：zMix(kk, iNG) = Worksheets(mySheet).Cells(65 + iNG, j + 2) / 100:
        }
        zMix[kk][14] = 0.0;
    }

    SumofZMix[kk] = 0.0;
    for (iNG = 0; iNG < 14; iNG++) {
        if (zMix[kk][iNG] < 1e-7) zMix[kk][iNG] = 0;
        SumofZMix[kk] += zMix[kk][iNG];
    }
    if (SumofZMix[kk] > 0) {
        for (iNG = 0; iNG < 14; iNG++)//Normalized z[0] to z[13]
            zMix[kk][iNG] /= SumofZMix[kk];
    }

    //为该子程序分配局部变量
    mc[iH] = HstpMix[kk]; mc[iNa] = NaMix[kk]; mc[iK] = KMix[kk]; mc[iMg] = MgMix[kk]; mc[iCa] = CaMix[kk];
    int TCa = mc[iCa];
    mc[iSr] = SrMix[kk]; mc[iBa] = BaMix[kk]; mc[iFe] = FeMix[kk]; mc[iZn] = ZnMix[kk]; mc[iPb] = PbMix[kk]; mc[iRa] = RaMix[kk];

    ma[iOH] = OHstpMix[kk]; ma[iCl] = ClMix[kk]; ma[iAc] = ACstpMix[kk]; mc[iNH4] = NH4STPMix[kk]; ma[iH2BO3] = H2BO3stpMix[kk];
    ma[iHCO3] = HCO3stpMix[kk]; ma[iCO3] = CO3stpMix[kk];
    ma[iSO4] = SO4Mix[kk]; ma[iHS] = HSstpMix[kk]; ma[intF] = FMix[kk]; ma[iBr] = BrMix[kk];

    Alk = AlkMix[kk]; TAc = TAcMix[kk]; TNH4 = TNH4Mix[kk]; TH3BO3 = TH3BO3Mix[kk]; TH2Saq = TH2SaqMix[kk]; TCO2 = TCO2Mix[kk];

    VW = VwMix[kk]; VgTP = VgTPMix[kk]; VO = VoMix[kk]; VMeOH = VMeOHMix[kk]; VMEG = VMEGMix[kk]; mass_MeOH = mass_MeOH_mix[kk]; mass_MEG = mass_MEG_mix[kk];

    yCO2 = yCO2Mix[kk], yH2S = yH2SMix[kk], yCH4 = yCH4Mix[kk];   // Local variable values; in this loop only.

    int useEOS = useEOSmix[kk]; int use_pH = usepHmix[kk]; int UseH2Sgas = UseH2SgasMix[kk];

    void* mt;


    TFe = mc[iFe];

    for (iNG = 0; iNG < 14; iNG++) {
        z[iNG] = zMix[kk][iNG];
        if (useEOS == 3 || useEOS == 0)
            z[iNG] = 0;
    }

    SGG = gasSpGravMix[kk]; API = oilAPIgravMix[kk];

    //根据Cragoe, C.S.1929《石油产品的热力学性质》（美国商务部标准局杂项出版物第97号）计算石油的平均分子量。
    if (API > 20 && API < 80) mw_oil = 6084.0 / (API - 5.9);
    if (API <= 20) mw_oil = 6084 / (20 - 5.9);
    if (API >= 80) mw_oil = 6084 / (80 - 5.9);

    //根据流体温度和压力计算水、油的质量以及水、油和气体的摩尔数

    if (UseTPVolMix[kk] == 1) fTPFunc(1);//检查 T、P 以进行 mass_W、mass_O、mol_g 计算
    if (UseTPVolMix[kk] == 0) fTPFunc(0);

    CalcIonicStrength();
    pH = pHMeterStpMix[kk] + DpHj;
    RatioOilBPoints = fRatioOilBPoints();
    C1_ThermodynamicEquilConsts();
    C2_PitzerActCoefs_T_P_ISt(gNeut, aH2O, TK, TC, PBar, Patm); //计算在标准温度和压力下水合物抑制剂存在下的作用系数
    PengRobinson3();
    //注意，如果 useEOSmix(kk)<>0，则省略此 pH 计算步骤。换句话说，如果此处 useEOS<>0，则不会运行 pH 和形态分析。形态分析已在 ReadInputPartC 中完成
    C5_CalcpHPCO2PH2SSTP();
    //重新分配 CO2aq、HCO3、CO3、H2Saq、HS 并重新计算离子强度
    //mt = fmn();
    fmn();

    CalcIonicStrength();
    C2_PitzerActCoefs_T_P_ISt(gNeut, aH2O, TK, TC, PBar, Patm);
    C5_CalcpHPCO2PH2SSTP();

    //使用温度 TVol 和 PVol 下计算的水密度
    if (UseTPVolMix[kk] == 1) WaterDensityMix[kk] = CalcRhoTP(TK, TC, PBar, Patm);

    mass_o_Mix[kk] = 159 * VoMix[kk] * OilDensityMix[kk]; // 这些用于在 useEOS=0 和 useEOS=3 中计算 nTCO2 和 nTH2S。这些值在 useEOS=1 或 2 中重新计算。
    mass_w_Mix[kk] = 159 * VwMix[kk] * WaterDensityMix[kk]; // 换算成公斤盐水
    mass_w_Mix[kk] = mass_w_Mix[kk] * (1 - CalculatedTDSMix[kk] / rho25CMix[kk] * 0.000001);

    Mass_o = mass_o_Mix[kk];
    mass_w = mass_w_Mix[kk];

    mt = fTotalCO2H2Smoles(*yCO2, RatioOilBPoints, Znew); //计算每种气体成分的总摩尔数（无论是在气体中还是在油或水中）
    // 请注意，VgTP 的单位是 m^3，当 Vg 的单位是 MMCF 时，829 是从 He 的系数 78740 转换而来的。

    yCO2Mix[kk] = *yCO2;
    yCH4Mix[kk] = *yCH4;
    yH2SMix[kk] = *yH2S;

    nTCO2Mix[kk] = nTCO2;
    nTCH4Mix[kk] = nTCH4;
    nTH2SMix[kk] = nTH2S;

    *YH2O = PsatH2O(TK) / PBar;

    *mol_w_Orig = 1000 * mass_w / 18.01528; // moles of water per day

    if (*usedryHC == 1)
        *mol_W = *mol_w_Orig;
    else
        *mol_W = 1000 * mass_w / 18.01528 + PsatH2O(TK) * VgTP * 1000 / Znew / RBar / TK; //包括气相中的 H2O，仅用作使用 EOS > 0 时的初始猜测



    //***************步骤 1  碳氢化合物调节**************

    //If useEOS <> 0 Then
    if (useEOS != 0)
    {
        if (useEOS == 3) {
            *mol_o = 1000 * Mass_o / (*mw_oil); // moles of oil per day
            *mol_g = VgTP * PBar * 1000 / (Znew * RBar * TK); // moles of gas per day
            *mol_HC = (*mol_o) + (*mol_g);

            if ((VoMix[kk] == (1.0 / 159.0 / 1000.0)) && (VgTP == 0.000001)) { // 当气体和油都为零并混合水时
                if (nob == 1) {
                    errmsg[13] = 14;
                    useEOS = 0;
                    goto label_3003;
                }
                else {
                    z_before_precipitation[1] = nTCO2Mix[kk] / (*mol_W);
                    z_before_precipitation[2] = nTH2SMix[kk] / (*mol_W);
                    z_before_precipitation[14] = 1 - z_before_precipitation[1] - z_before_precipitation[2];
                    Total_molesMix[kk] = (*mol_W);
                    nTCO2MixEOS[kk] = nTCO2;
                    nTH2SMixEOS[kk] = nTH2S;
                    mt = fTPFunc(0);

                    if (usepHmix[kk] == 1) *pH = pHMeterStpMix[kk] + DpHj;

                    CalcIonicStrength();
                    *RatioOilBPoints = fRatioOilBPoints();
                    C1_ThermodynamicEquilConsts();
                    C2_PitzerActCoefs_T_P_ISt(gNeut, aH2O, TK, TC, PBar, Patm); // 计算在标准温度和压力下水合物抑制剂存在下的作用系数
                    PengRobinson3();

                    // ***注意，如果 useEOSmix[kk]!=0，则省略此 pH 计算步骤。换句话说，如果此处 useEOS!=0，则不会运行 pH 和形态分析。形态分析已在 ReadInputPartC 中完成。
                    C5_CalcpHPCO2PH2SSTP();

                    // ******重新分配 CO2aq、HCO3、CO3、H2Saq、HS 并重新计算离子强度
                    mt = fmn;
                    CalcIonicStrength();
                    C2_PitzerActCoefs_T_P_ISt(gNeut, aH2O, TK, TC, PBar, Patm);
                    C5_CalcpHPCO2PH2SSTP();
                    *pHMeterReading = *pH - DpHj;
                    goto label_3003;
                }
            }
        }
        if (useEOS == 1 || useEOS == 2) {
            if (SumofZMix[kk] == 0) {       // 仅当给出碳氢化合物成分时，才运行 HC 调节
                if (Run_Seawater_Mixing == 1 || Run_MixingTwoWells == 1) {
                    goto label_3002;
                }
                else {
                    if (nob == 1) {
                        errmsg[13] = 14;  // 数组下标从0开始，14对应索引13
                        useEOS = 0;
                        goto label_3003;
                    }
                    else {
                        if (mol_W > 0) {
                            z_before_precipitation[1] = nTCO2Mix[kk] / (*mol_W);  // 原下标2对应索引1
                            z_before_precipitation[2] = nTH2SMix[kk] / (*mol_W);  // 原下标3对应索引2
                            z_before_precipitation[14] = 1 - z_before_precipitation[1] - z_before_precipitation[2];  // 原下标15对应索引14
                            Total_molesMix[kk] = (*mol_W);
                            nTCO2MixEOS[kk] = nTCO2;
                            nTH2SMixEOS[kk] = nTH2S;
                            goto label_3003;
                        }
                        else {
                            errmsg[13] = 14;  // 数组下标从0开始，14对应索引13
                            useEOS = 0;
                            goto label_3003;
                        }
                    }
                }
            }
            else {
                total_moles = 1;

                MultiPhaseFlash(mf_ParametersWereRead, mf_TCr, mf_PCr, mf_Omega, mf_MWgas, mf_kPr, mf_c0, mf_c1, TK, PBar, total_moles, z, tempgNeut,
                    aH2O, density, compositions, phi, Compr, beta, zOutput, mass_phase, MW_Phase, No_Phases);

                if (beta[0] > 0 && beta[1] > 0) {
                    if (VoMix[kk] > 1.0 / 159.0 / 1000.0) { // key to oil if vo is greater than 0
                        mass_o_Mix[kk] = 159 * VoMix[kk] * density[1]; // Kg
                        *mol_o = mass_o_Mix[kk] * 1000 / (*MW_Phase)[1];
                        *mol_g = (*mol_o) / beta[1] * beta[0];
                    }
                    else if (VgTP > 0.000001) {    // 仅当 Vg 大于 0 且 Vo=0 时才允许
                        *mol_g = VgTP * PBar * 1000 / (Compr[0] * RBar * TK);
                        *mol_o = (*mol_g) / beta[0] * beta[1];
                    }
                    else { // 当没有给出油气体积时，只允许计算混合
                        if (Run_Seawater_Mixing == 1 || Run_MixingTwoWells == 1) {
                            goto label_3002;
                        }
                        else {
                            if (nob == 1) {
                                mol_g = 0;
                                mol_o = 0;
                                errmsg[15] = 16;  // 数组下标从0开始，16对应索引15
                                useEOS = 0;
                                goto label_3003;
                            }
                            else {
                                if (mol_W > 0) {
                                    z_before_precipitation[1] = nTCO2Mix[kk] / (*mol_W);
                                    z_before_precipitation[2] = nTH2SMix[kk] / (*mol_W);
                                    z_before_precipitation[14] = 1 - z_before_precipitation[1] - z_before_precipitation[2];
                                    Total_molesMix[kk] = (*mol_W);
                                    nTCO2MixEOS[kk] = nTCO2;
                                    nTH2SMixEOS[kk] = nTH2S;
                                    goto label_3003;
                                }
                                else {
                                    errmsg[15] = 16;  // 数组下标从0开始，16对应索引15
                                    useEOS = 0;
                                    goto label_3003;
                                }
                            }
                        }
                    }
                }

                if (beta[0] > 0 && beta[1] == 0) { // 仅存在气相
                    if (VoMix[kk] > 1.0 / 159.0 / 1000.0) { // 如果存在石油量，则石油的关键
                        if ((z[7] + z[8] + z[9] + z[10] + z[11] + z[12]) > 0) { // 只有当 HC 大于 C4 时，才是石油的关键
                            mass_o_Mix[kk] = 159 * VoMix[kk] * density[0]; // Kg
                            (*mol_o) = mass_o_Mix[kk] * 1000 / (*MW_Phase)[0];
                            (*mol_g) = 0;
                        }
                        else {
                            (*mol_g) = VgTP * PBar * 1000 / (Compr[0] * RBar * TK); // 如果不存在大于 C4 的 HC，则为气体键
                            (*mol_o) = 0;
                        }
                    }
                    else if (VgTP > 0.000001) { // key to gas if vol of oil=0
                        (*mol_g) = VgTP * PBar * 1000 / (Compr[0] * RBar * TK);
                        (*mol_o) = 0;
                    }
                    else { // 如果没有给出石油或天然气产量
                        if (Run_Seawater_Mixing == 1 || Run_MixingTwoWells == 1) {
                            goto label_3002;
                        }
                        else {
                            if (nob == 1) {
                                mol_g = 0;
                                mol_o = 0;
                                errmsg[15] = 16;  // 数组下标从0开始，16对应索引15
                                useEOS = 0;
                                goto label_3003;
                            }
                            else {
                                if ((*mol_W) > 0) {
                                    z_before_precipitation[1] = nTCO2Mix[kk] / (*mol_W);
                                    z_before_precipitation[2] = nTH2SMix[kk] / (*mol_W);
                                    z_before_precipitation[14] = 1 - z_before_precipitation[1] - z_before_precipitation[2];
                                    Total_molesMix[kk] = (*mol_W);
                                    nTCO2MixEOS[kk] = nTCO2;
                                    nTH2SMixEOS[kk] = nTH2S;
                                    goto label_3003;
                                }
                                else {
                                    errmsg[15] = 16;  // 数组下标从0开始，16对应索引15
                                    useEOS = 0;
                                    goto label_3003;
                                }
                            }
                        }
                    }
                }

                if (beta[0] == 0 && beta[1] > 0) {
                    if (VoMix[kk] > 1.0 / 159.0 / 1000.0) { // key to oil if vol of oil is present
                        mass_o_Mix[kk] = 159 * VoMix[kk] * density[1]; // Kg
                        (*mol_o) = mass_o_Mix[kk] * 1000 / (*MW_Phase)[1];
                        (*mol_g) = 0;
                    }
                    else if (VgTP > 0.000001) { // key to gas if vol of oil=0
                        (*mol_g) = VgTP * PBar * 1000 / (Compr[1] * RBar * TK);
                        (*mol_o) = 0;
                    }
                    else { // if neither oil or gas vol is given
                        if (Run_Seawater_Mixing == 1 || Run_MixingTwoWells == 1) {
                            goto label_3002;
                        }
                        else {
                            if (nob == 1) {
                                errmsg[15] = 16;  // 数组下标从0开始，16对应索引15
                                useEOS = 0;
                                goto label_3003;
                            }
                            else {
                                if ((*mol_W) > 0) {
                                    z_before_precipitation[1] = nTCO2Mix[kk] / (*mol_W);
                                    z_before_precipitation[2] = nTH2SMix[kk] / (*mol_W);
                                    z_before_precipitation[14] = 1 - z_before_precipitation[1] - z_before_precipitation[2];
                                    Total_molesMix[kk] = (*mol_W);
                                    nTCO2MixEOS[kk] = nTCO2;
                                    nTH2SMixEOS[kk] = nTH2S;
                                    goto label_3003;
                                }
                                else {
                                    errmsg[15] = 16;  // 数组下标从0开始，16对应索引15
                                    useEOS = 0;
                                    goto label_3003;
                                }
                            }
                        }
                    }
                }

                (*mol_HC) = (*mol_o) + (*mol_g);

                if (useEOS == 1) {
                    nTCO2EOS = z[1] * (*mol_HC);  // 原z(2)对应z[1]
                    nTH2sEOS = z[2] * (*mol_HC);  // 原z(3)对应z[2]
                    nTCO2 = nTCO2EOS;
                    nTH2S = nTH2sEOS;
                }

                if (useEOS == 2) {
                    nTCO2EOS = nTCO2;
                    nTH2sEOS = nTH2S;
                    double zHC = z[0]; // 通过在 STP 条件下分解出（HCO3、CO3 和 HS）的总摩尔数来重新调整 nTCO2EOS 和 nTH2SEOS

                    for (iNG = 3; iNG < 14; iNG++) zHC = zHC + z[iNG]; // 总量中 z(碳氢化合物-CO2-H2s) 的摩尔分数


                    (*mol_HC) = (*mol_HC) * zHC + nTCO2EOS + nTH2sEOS;  // 如果 useEOS=2 nTCO2EOS，则修改 HC 的总摩尔数，计算气体、油和 CO2aq 中的 CO2 摩尔数（不包括 HCO3 和 CO3）

                    z[0] = z[0] * (*mol_HC);
                    z[1] = nTCO2EOS;
                    z[2] = nTH2sEOS;

                    for (iNG = 3; iNG < 14; iNG++) z[iNG] = z[iNG] * (*mol_HC);

                    for (iNG = 0; iNG < 14; iNG++) z[iNG] = z[iNG] / (*mol_HC);
                }
            } // end SumofZ>0
        }

        if (VW == 1.0 / 159.0 / 1000.0) {               // if water volume=0, skip Flash
            if (Run_Seawater_Mixing == 1 || Run_MixingTwoWells == 1) {
                goto label_3002;
            }
            else {
                errmsg[13] = 14;  // 数组下标从0开始，14对应索引13
                useEOS = 0;
                goto label_3003;
            }
        }
        //'*****使用EOS选项1或2
        //仅当给定实际的G/O/W体积时才允许进行Flash计算，例外情况是混合两口井和海水
        if (useEOS == 1 || useEOS == 2) {
            // 首先利用气相中的纯水蒸气计算储层成分
            if ((*usedryHC) == 0) { // 湿烃情况，总水量=输入水量+HC中的水量
                mol_w3 = 0;
                while (mol_w3 <= 0) { // 确保添加足够的水以开始迭代
                    (*mol_W) = (*mol_W) * 1.1;

                    true_composition(TK, PBar, *mol_HC, *mol_W, aH2O, tempgNeut, nTCO2EOS, nTH2sEOS, useEOS, z, feed_Composition, &total_moles,
                        MWgas, TCr, PCr, Omega, &mf_c0, &mf_c1, kPr,
                        composition_G, lnphi_Gas, composition_Aq, lnphi_Aq, Compr_composition_Aq, mf_reservoir_Composition, mf_feed_Composition, NumGases);

                    if (feed_Composition[0] < 0.0000001) {  // 原feed_Composition(1)对应feed_Composition[0]
                        feed_Composition[0] = 0;
                    }

                    MultiPhaseFlash(mf_ParametersWereRead, mf_TCr, mf_PCr, mf_Omega, mf_MWgas, mf_kPr, mf_c0, mf_c1, TK, PBar, total_moles, feed_Composition,
                        tempgNeut, aH2O, density, compositions, phi, Compr, beta, zOutput, mass_phase, MW_Phase, No_Phases);

                    if (compositions[14][3] > 0.5) {  // 原compositions(15,4)对应compositions[14][3]
                        mol_w3 = total_moles * beta[2] * compositions[14][3]; // H2O moles in aqueous phase
                    }
                    else {
                        mol_w3 = 0;
                    }
                }

                if (mol_w3 > 0) {
                    if (compositions[14][3] > 0.5) {
                        while (((mol_w3 - (*mol_w_Orig)) / (*mol_w_Orig)) * ((mol_w3 - (*mol_w_Orig)) / (*mol_w_Orig)) > 0.0001) {
                            (*mol_W) = (*mol_W) - mol_w3 + (*mol_w_Orig);
                            total_moles = (*mol_HC) + (*mol_W);
                            true_composition(TK, PBar, *mol_HC, *mol_W, aH2O, tempgNeut, nTCO2EOS, nTH2sEOS, useEOS, z, feed_Composition, &total_moles,
                                MWgas, TCr, PCr, Omega, &mf_c0, &mf_c1, kPr,
                                composition_G, lnphi_Gas, composition_Aq, lnphi_Aq, Compr_composition_Aq, mf_reservoir_Composition, mf_feed_Composition, NumGases);
                            MultiPhaseFlash(0, mf_TCr, mf_PCr, mf_Omega, mf_MWgas, mf_kPr, mf_c0, mf_c1, TK, PBar, total_moles, feed_Composition,
                                tempgNeut, aH2O, density, compositions, phi, Compr, beta, zOutput, mass_phase, MW_Phase, No_Phases);
                            mol_w3 = total_moles * beta[2] * compositions[14][3];
                        }
                    }
                    else {
                        errmsg[13] = 14;  // 数组下标从0开始，14对应索引13
                        useEOS = 0;
                        goto label_3003;
                    }
                }
                else {
                    errmsg[13] = 14;  // 数组下标从0开始，14对应索引13
                    useEOS = 0;
                    goto label_3003;
                }
            }
            else if (*usedryHC == 1) { // dry hydrocarbon, only water from Input
                true_composition(TK, PBar, *mol_HC, *mol_W, aH2O, tempgNeut, nTCO2EOS, nTH2sEOS, useEOS, z, feed_Composition, &total_moles,
                    MWgas, TCr, PCr, Omega, &mf_c0, &mf_c1, kPr,
                    composition_G, lnphi_Gas, composition_Aq, lnphi_Aq, Compr_composition_Aq, mf_reservoir_Composition, mf_feed_Composition, NumGases);

                if (feed_Composition[0] < 0.0000001) {
                    feed_Composition[0] = 0;
                }

                MultiPhaseFlash(mf_ParametersWereRead, mf_TCr, mf_PCr, mf_Omega, mf_MWgas, mf_kPr, mf_c0, mf_c1, TK, PBar, total_moles, feed_Composition,
                    tempgNeut, aH2O, density, compositions, phi, Compr, beta, zOutput, mass_phase, MW_Phase, No_Phases);
            }

            nTCO2MixEOS[kk] = total_moles * zOutput[1];  // 原zOutput(2)对应zOutput[1]
            nTH2SMixEOS[kk] = total_moles * zOutput[2];  // 原zOutput(3)对应zOutput[2]
            mass_w = total_moles * beta[2] * compositions[14][3] * 0.01801528; // mass_w is only the aqueous phase H2O
            Total_molesMix[kk] = total_moles;
            mass_w_Mix[kk] = mass_w;

            if (*No_Phases == 3) { // Output oil and gas density from flash calculation if useEOS=1 or useEOS=2
                GasDensityMix[kk] = density[0] * 1000;  // 原density(1)对应density[0]
                OilDensityMix[kk] = density[1];         // 原density(2)对应density[1]
            }

            if (*No_Phases == 2) { // Only g/o, G/W or O/W phases exist
                if (beta[2] == 0) {
                    //本例中不存在水相。增加水量并重新运行。程序中止。
                    printf("Aqueous phase does not existed in this case. Increase water volume and run again. Program abort.\n"); // Aqueous phase does not exist
                    exit(1);  // 替代VB的End
                }
                else {

                    if (density[0] > 0 && density[0] < 0.3) GasDensityMix[kk] = density[0] * 1000;// phase 1 is gas

                    else if (density[0] > 0 && density[0] > 0.3) OilDensityMix[kk] = density[0]; // phase 1 is oil

                    else if (density[1] > 0 && density[1] < 0.3)  GasDensityMix[kk] = density[1] * 1000;// Phase 2 is gas

                    else if (density[1] > 0 && density[1] > 0.3) OilDensityMix[kk] = density[1];

                }
            }

            if (*No_Phases == 1) {
                if (density[0] > 0 && compositions[14][1] < 0.8) { // hydrocarbon phase
                    printf("Aqueous phase does not existed in this case. Increase water volume and run again. Program abort.\n");
                    exit(1);
                }
                else if (density[1] > 0 && compositions[14][2] < 0.8) {
                    printf("Aqueous phase does not existed in this case. Increase water volume and run again. Program abort.\n");
                    exit(1);
                }
                else if (density[2] > 0 && compositions[14][3] < 0.8) {
                    printf("Aqueous phase does not existed in this case. Increase water volume and run again. Program abort.\n");
                    exit(1);
                }
            }
        }

        if (useEOS == 3) {
            // 仅在给定天然气和石油产量时运行 useEOS；第一次迭代涵盖 usedryHC=0 或 1
            pseudo_composition(API, SGG, VgTP, *mol_o, *mol_W, TK, PBar, aH2O, tempgNeut, nTCO2, nTH2S, *yCO2, *yH2S, *YH2O, &total_moles, feed_Composition, mol_HC,
                MWgas, TCr, PCr, Omega, &mf_c0, &mf_c1, kPr,
                composition_G, lnphi_Gas, composition_Aq, lnphi_Aq, Compr_composition_Aq, mf_reservoir_Composition, mf_feed_Composition, NumGases);

            if (feed_Composition[0] < 0.0000001) {
                feed_Composition[0] = 0;
            }

            SumofZ = feed_Composition[0];
            for (iNG = 1; iNG < 14; iNG++) {
                if (feed_Composition[iNG] < 0.0000001) feed_Composition[iNG] = 0;
                SumofZ = SumofZ + feed_Composition[iNG];
            }

            if (SumofZ < 0.000001) {
                errmsg[13] = 14;
                useEOS = 0;
            }

            MultiPhaseFlash(0, mf_TCr, mf_PCr, mf_Omega, mf_MWgas, mf_kPr, mf_c0, mf_c1, TK, PBar, total_moles, feed_Composition, tempgNeut,
                aH2O, density, compositions, phi, Compr, beta, zOutput, mass_phase, MW_Phase, No_Phases);

            if (*usedryHC == 0) {
                mol_w3 = total_moles * beta[2] * compositions[14][3];

                while (mol_w3 <= 0) { // 确保添加足够的水以开始迭代
                    (*mol_W) = (*mol_W) * 1.1;
                    pseudo_composition(API, SGG, VgTP, *mol_o, *mol_W, TK, PBar, aH2O, tempgNeut, nTCO2, nTH2S, *yCO2, *yH2S, *YH2O, &total_moles, feed_Composition, mol_HC,
                        MWgas, TCr, PCr, Omega, &mf_c0, &mf_c1, kPr,
                        composition_G, lnphi_Gas, composition_Aq, lnphi_Aq, Compr_composition_Aq, mf_reservoir_Composition, mf_feed_Composition, NumGases);
                    MultiPhaseFlash(mf_ParametersWereRead, mf_TCr, mf_PCr, mf_Omega, mf_MWgas, mf_kPr, mf_c0, mf_c1, TK, PBar, total_moles, feed_Composition,
                        tempgNeut, aH2O, density, compositions, phi, Compr, beta, zOutput, mass_phase, MW_Phase, No_Phases);

                    if (compositions[14][3] > 0.5)  mol_w3 = total_moles * beta[2] * compositions[14][3]; // 水相中的 H2O 摩尔数
                    else mol_w3 = 0;

                }

                if (mol_w3 > 0) {
                    while (((mol_w3 - (*mol_w_Orig)) / (*mol_w_Orig)) * ((mol_w3 - (*mol_w_Orig)) / (*mol_w_Orig)) > 0.0001) {
                        (*mol_W) = (*mol_W) - mol_w3 + (*mol_w_Orig);
                        total_moles = total_moles - mol_w3 + (*mol_w_Orig);

                        // recalculate z(i) into reservoir composition to use true_composition sub
                        for (iNG = 0; iNG < 14; iNG++) {
                            z[iNG] = zOutput[iNG] / (1 - zOutput[14]);
                        }

                        pseudo_composition(API, SGG, VgTP, *mol_o, *mol_W, TK, PBar, aH2O, tempgNeut,
                            nTCO2, nTH2S, *yCO2, *yH2S, *YH2O, &total_moles, feed_Composition, mol_HC,
                            MWgas, TCr, PCr, Omega, &mf_c0, &mf_c1, kPr,
                            composition_G, lnphi_Gas, composition_Aq, lnphi_Aq, Compr_composition_Aq, mf_reservoir_Composition, mf_feed_Composition, NumGases);

                        MultiPhaseFlash(0, mf_TCr, mf_PCr, mf_Omega, mf_MWgas, mf_kPr, mf_c0, mf_c1, TK, PBar, total_moles, feed_Composition, tempgNeut,
                            aH2O, density, compositions, phi, Compr, beta, zOutput, mass_phase, MW_Phase, No_Phases);
                        mol_w3 = total_moles * beta[2] * compositions[14][3];
                    }
                }

                nTCO2MixEOS[kk] = total_moles * zOutput[1];
                nTH2SMixEOS[kk] = total_moles * zOutput[2];
                mass_w = total_moles * beta[2] * compositions[14][3] * 0.01801528; // mass_w 仅为水相 H2O
                Total_molesMix[kk] = total_moles;
                mass_w_Mix[kk] = mass_w;
            }
            else {
                Total_molesMix[kk] = total_moles;
            }
        } // Correspond to useEOS=3



        // ***** Run speciation and pH calculation for useEOS>0
        nTCO2_before_precipitation = nTCO2EOS;
        nTH2S_before_precipitation = nTH2sEOS;
        Total_moles_before_precipitation = total_moles;

        for (int iz = 0; iz < 15; iz++) {
            z_before_precipitation[iz] = zOutput[iz];
            zMix[kk][iz] = zOutput[iz];
        }

        if (*usedryHC == 1) {  // 如果假设碳氢化合物为干碳氢化合物，则重新计算 T、P 处的真实水浓度
            for (int c = 0; c < NumCat; c++) mc[c] = mc[c] * (*mol_w_Orig) * 0.01801528 / mass_w;

            for (int a = 0; a < NumAn; a++)  ma[a] = ma[a] * (*mol_w_Orig) * 0.01801528 / mass_w;


            for (int n = 0; n < NumNeut; n++) mn[n] = mn[n] * (*mol_w_Orig) * 0.01801528 / mass_w;


            Alk = Alk * (*mol_w_Orig) * 0.01801528 / mass_w;
            TAc = TAc * (*mol_w_Orig) * 0.01801528 / mass_w;
            TNH4 = TNH4 * (*mol_w_Orig) * 0.01801528 / mass_w;
            TH3BO3 = TH3BO3 * (*mol_w_Orig) * 0.01801528 / mass_w;
            TH2Saq = TH2Saq * (*mol_w_Orig) * 0.01801528 / mass_w;
            TH4SiO4 = TH4SiO4 * (*mol_w_Orig) * 0.01801528 / mass_w;

            HstpMix[kk] = mc[iH];  // 假设iH等是1-based索引
            NaMix[kk] = mc[iNa];
            KMix[kk] = mc[iK];
            MgMix[kk] = mc[iMg];
            CaMix[kk] = mc[iCa];
            TCa = mc[iCa];
            SrMix[kk] = mc[iSr];
            BaMix[kk] = mc[iBa];
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

            AlkMix[kk] = Alk;
            TAcMix[kk] = TAc;
            TH2SaqMix[kk] = TH2Saq;
            TH4SiO4Mix[kk] = TH4SiO4;
            TNH4Mix[kk] = TNH4;
            TH3BO3Mix[kk] = TH3BO3;
        }

        // **** flash and convert to STP condition
        // Call C4_EOS_TCO2_SSPEquilCalcs(0, 2, 2, KspCalcite) 'calculate Separator or STP pH if useEOS>0
        // *********
        // mt = fmn  'Assign neutral species and HCO3, CO3, HS and recalculate ionic strength, Pitzer Coeff and pH
        // Call CalcIonicStrength
        // Call C2_PitzerActCoefs_T_P_ISt(gNeut, aH2O, TK, TC, PBar, Patm)
        // Call C4_EOS_TCO2_SSPEquilCalcs(0, 5, 2, KspCalcite) '(ppt_or_not, iMetals, iGas, Ksp)

        mt = fTPFunc(0);
        C4_SSPEquilCalcs(0, 5, 2, KspCalcite);  // (ppt_or_not, iMetals, iGas, Ksp)
        mt = fmn;        // CO2aq, HCO3, CO3, H2Saq, HS
        C2_PitzerActCoefs_T_P_ISt(gNeut, aH2O, TK, TC, PBar, Patm);
        C4_SSPEquilCalcs(0, 5, 2, KspCalcite);  // (ppt_or_not, iMetals, iGas, Ksp)

        *pHMeterReading = *pH - DpHj; // at separator T&P or STP
        *rhoTP = CalcRhoTP(TK, TC, PBar, Patm); // at separator T&P or STP

        for (int c = 0; c < NumCat; c++)
            molc[c][kk] = mc[c] * mass_w;

        for (int a = 0; a < NumAn; a++)
            mola[a][kk] = ma[a] * mass_w;

        for (int n = 0; n < NumNeut; n++)
            moln[n][kk] = mn[n] * mass_w;

        molAlk[kk] = Alk * mass_w;
        molTAC[kk] = TAc * mass_w;
        molTNH4[kk] = TNH4 * mass_w;
        molTH3BO3[kk] = TH3BO3 * mass_w;
        molTH2Saq[kk] = TH2Saq * mass_w;
        molTH4SiO4[kk] = TH4SiO4 * mass_w;

    }

label_3002:
    if (Run_Seawater_Mixing == 1 || Run_MixingTwoWells == 1) {
        if (Run_Seawater_Mixing == 1) {
            if (LoopMixing == 1 && kk == 1) { // when Fr of brine#1=0
                for (iNG = 0; iNG < 15; iNG++)
                    z_before_precipitation[iNG] = 0;

                Total_molesMix[kk] = 0;
                nTCO2MixEOS[kk] = 0;
                nTH2SMixEOS[kk] = 0;
                mass_w_Mix[kk] = 0;
                useEOS = 0;
                goto label_3003;
            }

            if (LoopMixing == 11 && kk == 2) { // when Fr of seawater=0
                for (iNG = 0; iNG < 15; iNG++)
                    z_before_precipitation[iNG] = 0;

                Total_molesMix[kk] = 0;
                nTCO2MixEOS[kk] = 0;
                nTH2SMixEOS[kk] = 0;
                mass_w_Mix[kk] = 0;
                goto label_3003;
            }

            if (kk == 2) { // Set EOS parameter for seawater
                Total_molesMix[kk] = (*mol_W);
                z_before_precipitation[1] = nTCO2Mix[kk] / Total_molesMix[kk];
                z_before_precipitation[2] = nTH2SMix[kk] / Total_molesMix[kk];
                z_before_precipitation[14] = 1 - z_before_precipitation[1] - z_before_precipitation[2];
                goto label_3003;
            }
        }

        if (Run_MixingTwoWells == 1) {
            if (LoopMixing == 1 && kk == 1) { // when Fr Brine#1 =0
                for (iNG = 0; iNG < 15; iNG++)
                    z_before_precipitation[iNG] = 0;

                Total_molesMix[kk] = 0;
                nTCO2MixEOS[kk] = 0;
                nTH2SMixEOS[kk] = 0;
                mass_w_Mix[kk] = 0;
                goto label_3003;
            }

            if (LoopMixing == 11 && kk == 2) { // When Fr brine#2 =0
                for (iNG = 0; iNG < 15; iNG++)
                    z_before_precipitation[iNG] = 0;

                Total_molesMix[kk] = 0;
                nTCO2MixEOS[kk] = 0;
                nTH2SMixEOS[kk] = 0;
                mass_w_Mix[kk] = 0;
                goto label_3003;
            }

            if (VgTP == 0.000001 && VO == 1.0 / 159.0 / 1000.0 && VW > 1.0 / 159.0 / 1000.0) { // when only water present for either brine 1 or 2
                Total_molesMix[kk] = (*mol_W);
                z_before_precipitation[1] = nTCO2Mix[kk] / Total_molesMix[kk];
                z_before_precipitation[2] = nTH2SMix[kk] / Total_molesMix[kk];
                z_before_precipitation[14] = 1 - z_before_precipitation[0] - z_before_precipitation[1];
                goto label_3003;
            }

            if (VW == 1.0 / 159.0 / 1000.0) { // When water volume is not given for either Brine#1 or Brine#2
                useEOS = 0;
                goto label_3003;
            }
        }
    }

label_3003:
    // 将z_before_precipitation复制到zMix
    for (iNG = 0; iNG < 15; iNG++)
        zMix[kk][iNG] = z_before_precipitation[iNG];


    // ***** QC calculation
    if (RunStat == 1) {
        if (RunQualityControlChecks == 1) {
            // Worksheets(mySheet).Activate 
            //int kk, int j,
            //    int* UseH2Sgas, double* pH, double* aH, double* yH2S, int useEOS, double K1H2S, double KgwH2S,
            //    double Ppsia, double KstFeSaq, double TZn, double TPb, double hydHS, double ppt, double* BetaDot, double* gDot,
            //    double KH2O, double KHAc, double KH3BO3, double KNH4, double KH4SiO4, double KH3SiO3, double K1H2CO3, double KgwCO2,
            //    double K2HCO3,
            //    double* yCO2,
            //    double* pHMeterReading
            QualityControlCalculations(kk, j, &UseH2Sgas, pH, aH, yH2S, useEOS, K1H2S, KgwH2S, Ppsia, KstFeSaq, TZn, TPb, hydHS, ppt, BetaDot, gDot,
                KH2O, KHAc, KH3BO3, KNH4, KH4SiO4, KH3SiO3, K1H2CO3, KgwCO2, K2HCO3, yCO2, pHMeterReading);
        }
    }
    else {
        if (RunQualityControlChecks == 1 && kk <= nob_Input) { // only run QC if requested from Input Sheet.
            // Worksheets(mySheet).Activate
            QualityControlCalculations(kk, j, &UseH2Sgas, pH, aH, yH2S, useEOS, K1H2S, KgwH2S, Ppsia, KstFeSaq, TZn, TPb, hydHS, ppt, BetaDot, gDot,
                KH2O, KHAc, KH3BO3, KNH4, KH4SiO4, KH3SiO3, K1H2CO3, KgwCO2, K2HCO3, yCO2, pHMeterReading);
            if (kk == nob) {
                goto label_123;   // Last Brine QC has been printed to Input Sheet; exit calculations.
            }
        }
        if (RunQualityControlChecks_II == 1 && kk > nob_Input) { // only run QC if requested from InputII Sheet.
            // Worksheets(mySheet).Activate
            QualityControlCalculations(kk, j, &UseH2Sgas, pH, aH, yH2S, useEOS, K1H2S, KgwH2S, Ppsia, KstFeSaq, TZn, TPb, hydHS, ppt, BetaDot, gDot,
                KH2O, KHAc, KH3BO3, KNH4, KH4SiO4, KH3SiO3, K1H2CO3, KgwCO2, K2HCO3, yCO2, pHMeterReading);
            if (kk == nob_Input + nob_InputII) {
                exit(0);//vb源码 ： End 
            }
        }
    }

    // *****End of QC calculation
    // ***** 将 pH、PCO2、PH2S、TH2Saq 打印到输入表
    // If RunShellMultiflash <> 1 Then
    if (use_pH == 0 && useEOS == 0) {        // 然后，必须指定yCO2stp。需要计算pH值和yH2Sstp。
        if (Run_Seawater_Mixing == 0 && Run_MixingTwoWells == 0 && RunMultiMix == 0 && Run10TestCases == 0 && RunWhatIf == 0 && RunStat == 0) {
            if (Run_CalcConcFactor == 1) {
                // 在C中，可能需要调用特定的输出函数
                //待定：Worksheets(MySheetMix).Cells(34, 8) = pHMeterReading
            }
            else {
                //待定： Worksheets(mySheet).Cells(34, j + 2) = pHMeterReading
            }
        }
    }

    if (use_pH == 1 && useEOS == 0) {                            // Use the pH measured at STP.
        if (RunWhatIf != 1 && RunStat == 0) {
            // 待定： Worksheets(mySheet).Cells(31, j + 2) = yCO2 * 100

        }
    }

    if (use_pH == 2 && useEOS == 0) {                            // Use the pH measured at STP.
        if (RunWhatIf != 1 && RunStat == 0) {
            // 待定： Worksheets(mySheet).Cells(24, j + 2) = Alk * (61019 * (rho25c - TDS / 1000000.0))
        }
        AlkMix[kk] = Alk;
    }

    if (use_pH == 3 && useEOS == 0) {         // use TCO2 to calculate pH
        if (Run_Seawater_Mixing == 0 && Run_MixingTwoWells == 0 && RunMultiMix == 0 && Run10TestCases == 0 && RunWhatIf == 0 && RunStat == 0) {
            *pHMeterReading = *pH - DpHj;
            if (Run_CalcConcFactor == 1) {
                // 待定： Worksheets(MySheetMix).Cells(34, 8) = pHMeterReading
            }
            else {
                //待定： Worksheets(mySheet).Cells(34, j + 2) = pHMeterReading
            }
        }
        if (RunWhatIf != 1) {
            //待定：  Worksheets(mySheet).Cells(31, j + 2) = yCO2 * 100
        }
    }

    if (useEOS != 0) {
        if (Run_Seawater_Mixing == 0 && Run_MixingTwoWells == 0 && RunMultiMix == 0 && Run10TestCases == 0 && RunWhatIf == 0 && RunStat == 0) {
            if (Run_CalcConcFactor == 1) {
                //待定：  Worksheets(MySheetMix).Cells(34, 8) = pHMeterReading
            }
            else {
                //待定：  Worksheets(mySheet).Cells(34, j + 2) = pHMeterReading
            }
        }
    }

    if (UseMolal == 1) {
        //待定：  Worksheets(mySheet).Cells(29, j + 2) = CalculatedTDSMix(kk)
    }

    if (mass_MeOH > 0) {
        xMeOH = (mass_MeOH / 32.0) / (mass_MeOH / 32.0 + mass_w / 18.0); // ?????????????????????????????
    }

    if (mass_MEG > 0) {
        xMEG = (mass_MEG / 62.07) / (mass_MEG / 62.07 + mass_w / 18.0);
    }

    if (mass_MeOH > 0) {
        IStCosolvent = Ist * mass_w / (mass_w + mass_MeOH);
    }

    if (mass_MEG > 0) {
        IStCosolvent = Ist * mass_w / (mass_w + mass_MEG);
    }

label_123:
    // End Sub - 在C中，这可能是函数的结束
    // 假设这是一个函数的结尾
    return;


}










void B2_ReadinAllData(SampleData* data)
{
    if (Read_InputII == 1)
    {
        nob = nob_Input;
    }

    if (Run1000Cases == 1)
    {
        nob = nob_Input;
    }
    for (int iRead = 0; iRead < nob; iRead++)
    {
        printf("%d", iRead);
        if (Run1000Cases == 1)
        {
            CaseCount[iRead] = LoopTP1000Cases;
        }
        j = CaseCount[iRead];
        kk = iRead;
        ReadInputPartA(kk, data);
    }

    for (int iRead = 0; iRead < nob; iRead++)
    {
        j = CaseCount[iRead];
        kk = iRead;
        ReadInputPartB(kk,data);
    }

    for (int iRead = 0; iRead < nob; iRead++)
    {
        j = CaseCount[iRead];
        kk = iRead;
        //ReadInputPartC(kk);
        if (RunStat == 0)
        {
            if (RunH2SGUI != 1)
            {
                if (Run_CalcConcFactor == 1)
                {
                    // Worksheets(MySheetMix).Cells(30, 8).Value = rho25c
                }
                else
                {
                    // Worksheets(mySheet).Cells(30, j + 2).Value = rho25c
                }
            }
        }
    }
    for (int iRead = 0; iRead < nob; iRead++)
    {
        j = CaseCount[iRead];
        kk = iRead;
        //ReadInputPartD(kk, j);
    }

    if (nob > 0)
    {
        if (nob + nob_InputII > 1 && RunQualityControlChecks == 0)
        {
            j = 6;
        }
        // TBH = Worksheets(mySheet).Cells(39, j + 2)
        if (RunMultiMix == 1)
        {
            // TBH = Worksheets("MultiMix").Cells(2 + LoopResChem, 3).Value
        }
        if (RunWhatIf == 1)
        {
            TBH = TWIInit;
        }

        if (UseSI == 1)
        {
            TBH = TBH * 9 / 5 + 32;
        }

        if (TBH < 24.8 && Run1000Cases != 1)
        {
            // MsgBox("Initial temperature is below 24.8F or -4C. Initial temperature will be set to 77 F or 25 C.")
        }

        if (TBH < 24.8)
        {
            if (UseSI == 1)
            {
                // Worksheets(mySheet).Cells(39, j + 2) = 250;
            }
            else
            {
                // Worksheets(mySheet).Cells(39, j + 2) = 482;
            }
        }

        if (TBH < 24.8)
        {
            TBH = 77;
        }

        if (TBH > 482 && Run1000Cases != 1)
        {
            // MsgBox ("Initial temperature is above 482 F. Initail temperature will be set to 482 F or 250 C.")
        }
        if (TBH > 482)
        {
            if (UseSI == 1)
            {
                // Worksheets(mySheet).Cells(39, j + 2) = 250
            }
            else
            {
                // Worksheets(mySheet).Cells(39, j + 2) = 482
            }
        }
        if (TBH > 482)
        {
            TBH = 482;
        }

        // TWH = Worksheets(mySheet).Cells(40, j + 2)
        if (RunMultiMix == 1)
        {
            // TWH = Worksheets("MultiMix").Cells(2 + LoopResChem, 3).Value
        }
        if (RunWhatIf == 1)
        {
            TWH = TWIInit;
        }
        if (UseSI == 1)
        {
            TWH = TWH * 9 / 5 + 32;
        }

        if (TWH < 24.8 && Run1000Cases != 1)
        {
            // MsgBox("Final temperature is below 24.8F or -4C. Final temperature will be set to 77 F or 25 C.")
        }

        if (TWH < 24.8)
        {
            if (UseSI == 1)
            {
                // Worksheets(mySheet).Cells(40, j + 2) = 25
            }
            else
            {
                // Worksheets(mySheet).Cells(40, j + 2) = 77
            }
        }
        if (TWH < 24.8)
        {
            TWH = 77;
        }
        if (TWH > 482 && Run1000Cases != 1)
        {
            // MsgBox ("Final temperature is above 482 F. Final temperature will be set to 482 F or 250 C.")
        }
        if (TWH > 482)
        {
            if (UseSI == 1)
            {
                // Worksheets(mySheet).Cells(40, j + 2) = 250
            }
            else
            {
                // Worksheets(mySheet).Cells(40, j + 2) = 482
            }
        }
        if (TWH > 482)
        {
            TWH = 482;
        }

        // PBH = Worksheets(mySheet).Cells(41, j + 2);

        if (RunMultiMix == 1)
        {
            // PBH = Worksheets("MultiMix").Cells(2 + LoopResChem, 4).Value;
        }

        if (RunWhatIf == 1)
        {
            PBH = PWIInit;
        }

        if (UseSI == 1)
        {
            PBH = PBH * 14.503774; // note that TBH, TWH, PBH, PWH is defaulted to F, and psia
        }
        if (PBH < 12 && Run1000Cases != 1)
        {
            // MsgBox("Initial pressure is below 12 psia or 0.827 bar. Initial pressure will be set to 14.7 psi or 1.01325 bar.");
        }
        if (PBH < 12)
        {
            if (UseSI == 1)
            {
                // Worksheets(mySheet).Cells(41, j + 2) = 1.01325;
            }
            else
            {
                // Worksheets(mySheet).Cells(41, j + 2) = 14.7;
            }
        }

        if (PBH < 12)
        {
            PBH = 14.7; // If P is zero, probably it is gauge and therefore add 1 atm
        }
        if (PBH > 30000 && Run1000Cases != 1)
        {
            // MsgBox("Initial pressure is above 30,000 psi. Initial pressure will be set to 30,000 psia or 2,068.4 bar.");
        }
        if (PBH > 30000)
        {
            if (UseSI == 1)
            {
                // Worksheets(mySheet).Cells(41, j + 2) = 30000 / 14.503774;
            }
            else
            {
                // Worksheets(mySheet).Cells(41, j + 2) = 30000;
            }
        }

        if (PBH > 30000)
        {
            PBH = 30000;
        }

        // PWH = Worksheets(mySheet).Cells(42, j + 2);

        if (RunMultiMix == 1)
        {
            // PWH = Worksheets("MultiMix").Cells(2 + LoopResChem, 4).Value;
        }

        if (RunWhatIf == 1)
        {
            PWH = PWIInit;
        }

        if (UseSI == 1)
        {
            PWH = PWH * 14.503774;
        }

        if (PWH < 12 && Run1000Cases != 1)
        {
            // MsgBox("Final pressure is below 12 psia or 0.827 bar. Final pressure will be set to 14.7 psia or 1.01325 bar.");
        }

        if (PWH < 12)
        {
            if (UseSI == 1)
            {
                // Worksheets(mySheet).Cells(42, j + 2) = 1.01325;
            }
            else
            {
                // Worksheets(mySheet).Cells(42, j + 2) = 14.7;
            }
        }

        if (PWH < 12)
        {
            PWH = 14.7; // If P is zero, probably it is gauge and therefore add 1 atm
        }

        if (PWH > 30000 && Run1000Cases != 1)
        {
            // MsgBox("Final pressure is above 30,000 psi. Final pressure will be set to 30000 psia or 2068.4 bar.");
        }

        if (PWH > 30000)
        {
            if (UseSI == 1)
            {
                // Worksheets(mySheet).Cells(42, j + 2) = 30000 / 14.503774;
            }
            else
            {
                // Worksheets(mySheet).Cells(42, j + 2) = 30000;
            }
        }

        if (PWH > 30000)
        {
            PWH = 30000;
        }
    }

    if (tInh == 0)
    {
        tInh = 1;
    }

    for (int i = 0; i < 10; i++)
    {
        // InhName[i] = Worksheets("Input").Cells(1 + i, 18);
    }
    for (int i = 10; i < 20; i++)
    {
        // InhName(i) = Worksheets("Input").Cells(1 + i - 10, 20)
    }

    NoRiskcalc = 3;
    for (int i = 1; i < (NoRiskcalc - 1); i++)
    {
        ConcInhBarRisk[i] = MaxInh / (NoRiskcalc - 1) * i;
        ConcInhCalRisk[i] = MaxInh / (NoRiskcalc - 1) * i;
        ConcInhAnRisk[i] = MaxInh / (NoRiskcalc - 1) * i;
        ConcInhGypRisk[i] = MaxInh / (NoRiskcalc - 1) * i;

        ConcInhCelRisk[i] = MaxInh / (NoRiskcalc - 1) * i;
    };
}

int main()
{
    SampleData data;
    fetchData(&data);
    initData();
    B2_ReadinAllData(&data);
}


