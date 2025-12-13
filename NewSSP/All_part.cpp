#include <stdio.h>
#include <math.h>
#include<stdlib.h>
#include <stdio.h>
#include <string.h>
#include<float.h>
#include <vector>
#include <iostream>
#include<fstream>
#include<string>
#include <map>
#pragma warning(disable:4996)

template<typename T>
using Vec = std::vector<T>;
using DoubleVec = Vec<double>;
using DoubleMatrix = Vec<Vec<double>>;
using IntVec = Vec<int>;
using IntMatrix = Vec<Vec<int>>;
using BoolVec = Vec<bool>;
using StrVec = Vec<std::string>;


class SimulationContext {
public:
    SimulationContext() {
        int len = nob_Input + nob_InputII;
        NaMix.resize(len, 0);
        MgMix.resize(len, 0);
        CaMix.resize(len, 0);
        SrMix.resize(len, 0);
        BaMix.resize(len, 0);
        FeMix.resize(len, 0);
        ZnMix.resize(len, 0);
        ClMix.resize(len, 0);
        PbMix.resize(len, 0);
        BrMix.resize(len, 0);
        RaMix.resize(len, 0);

        NH3Mix.resize(len, 0);
        H3SiO4Mix.resize(len, 0);
        H2SiO4Mix.resize(len, 0);
        H4SiO4Mix.resize(len, 0);
        H3BO3Mix.resize(len, 0);
        CO2aqMix.resize(len, 0);
        H2SaqMix.resize(len, 0);
        HACaqMix.resize(len, 0);

        UseH2SgasMix.resize(len, 0);
        SO4Mix.resize(len, 0);
        FMix.resize(len, 0);
        TDSMix.resize(len, 0);
        AlkMix.resize(len, 0);
        TAcMix.resize(len, 0);
        KMix.resize(len, 0);
        MixFrac.resize(len, 0);

        rho_Mix.resize(len, 0);
        TH2SaqMix.resize(len, 0);
        pHMeterStpMix.resize(len, 0);
        TH4SiO4Mix.resize(len, 0);
        TNH4Mix.resize(len, 0);
        TH3BO3Mix.resize(len, 0);

        SampleIDMix.resize(len, 0);
        SampleDateMix = "";
        OperatorMix.resize(len, 0);
        WellNameMix = "";
        FieldMix.resize(len, 0);
        StateMix = "";

        VgTPMix.resize(len, 0);
        VoMix.resize(len, 0);
        VwMix.resize(len, 0);
        VMeOHMix.resize(len, 0);
        VMEGMix.resize(len, 0);
        oilAPIgravMix.resize(len, 0);
        gasSpGravMix.resize(len, 0);
        MixFracGas.resize(len, 0);

        nTCO2Mix.resize(len, 0);
        nTCH4Mix.resize(len, 0);
        nTH2SMix.resize(len, 0);
        mass_w_Mix.resize(len, 0);
        mass_o_Mix.resize(len, 0);
        MixFracOil.resize(len, 0);

        mass_MeOH_mix.resize(len, 0);
        mass_MEG_mix.resize(len, 0);
        Qheat.resize(len, 0);
        yCO2Mix.resize(len, 0);
        yH2SMix.resize(len, 0);
        yCH4Mix.resize(len, 0);

        YCH4stpmix.resize(len, 0);
        RatioOilBPointsmix.resize(len, 0);
        CalculatedTDSMix.resize(len, 0);
        rho25CMix.resize(len, 0);

        HstpMix.resize(len, 0);
        OHstpMix.resize(len, 0);
        HCO3stpMix.resize(len, 0);
        CO3stpMix.resize(len, 0);
        ACstpMix.resize(len, 0);
        HSstpMix.resize(len, 0);
        NH4STPMix.resize(len, 0);
        H2BO3stpMix.resize(len, 0);

        HCO3AlkMix.resize(len, 0);
        CO3AlkMix.resize(len, 0);
        HAlkMix.resize(len, 0);
        OHAlkMix.resize(len, 0);
        ConcFactor.resize(len, 0);

        TCO2Mix.resize(len, 0);
        TofpH.resize(len, 0);
        PofpH.resize(len, 0);
        TofVol.resize(len, 0);
        PofVol.resize(len, 0);

        OilDensityMix.resize(len, 0);
        GasDensityMix.resize(len, 0);
        WaterDensityMix.resize(len, 0);
        UseTPpHMix.resize(len, 0);
        UseTPVolMix.resize(len, 0);
        useEOSmix.resize(len, 0);

        molAlk.resize(len, 0);
        molTAC.resize(len, 0);
        molTNH4.resize(len, 0);
        molTH3BO3.resize(len, 0);
        molTH2Saq.resize(len, 0);
        molTH4SiO4.resize(len, 0);

        mol_g_origMix.resize(len, 0);
        mol_o_OrigMix.resize(len, 0);
        mol_w_OrigMix.resize(len, 0);
        mol_g_finalMix.resize(len, 0);
        mol_o_finalMix.resize(len, 0);
        mol_w_finalMix.resize(len, 0);
        mol_w_evapMix.resize(len, 0);

        Total_molesMix.resize(len, 0);
        SumofZMix.resize(len, 0);
        nTCO2MixEOS.resize(len, 0);
        nTH2SMixEOS.resize(len, 0);

        usepHmix.resize(len, 0);


        zMix.resize(len, DoubleVec(15, 0));


        //---------------------------------
        CaseCount.resize(5); CaseCount[0] = 1; CaseCount[1] = 2;

        InhName.resize(21, "");

        ConcInhCelRisk.resize(10, 0);
        ConcInhCalRisk.resize(10, 0);
        ConcInhBarRisk.resize(10, 0);
        ConcInhAnRisk.resize(10, 0);
        ConcInhGypRisk.resize(10, 0);

        CaseCountWI.resize(40, 0);

        NaWI.resize(40, DoubleVec(10, 0));
        CaWI.resize(40, DoubleVec(10, 0));
        SrWI.resize(40, DoubleVec(10, 0));
        BaWI.resize(40, DoubleVec(10, 0));
        ClWI.resize(40, DoubleVec(10, 0));
        SO4WI.resize(40, DoubleVec(10, 0));
        HCO3AlkWI.resize(40, DoubleVec(10, 0));
        CO3AlkWI.resize(40, DoubleVec(10, 0));
        TACWI.resize(40, DoubleVec(10, 0));
        YCO2WI.resize(40, DoubleVec(10, 0));
        YH2SWI.resize(40, DoubleVec(10, 0));
        TH2SaqWI.resize(40, DoubleVec(10, 0));
        pHMeterSTPWI.resize(40, DoubleVec(10, 0));
        ConcFactorWI.resize(40, DoubleVec(10, 0));

        nob = 1;

        errmsg.resize(20, 0);
    }

public:
    DoubleVec NaMix;
    DoubleVec MgMix;
    DoubleVec CaMix;
    DoubleVec SrMix;
    DoubleVec BaMix;
    DoubleVec FeMix;
    DoubleVec ZnMix;
    DoubleVec ClMix;
    DoubleVec PbMix;
    DoubleVec BrMix;
    DoubleVec RaMix;

    DoubleVec NH3Mix;
    DoubleVec H3SiO4Mix;
    DoubleVec H2SiO4Mix;
    DoubleVec H4SiO4Mix;
    DoubleVec H3BO3Mix;
    DoubleVec CO2aqMix;
    DoubleVec H2SaqMix;
    DoubleVec HACaqMix;

    IntVec UseH2SgasMix;
    DoubleVec SO4Mix;
    DoubleVec FMix;
    DoubleVec TDSMix;
    DoubleVec AlkMix;
    DoubleVec TAcMix;
    DoubleVec KMix;
    DoubleVec MixFrac;

    DoubleVec rho_Mix;
    DoubleVec TH2SaqMix;
    DoubleVec pHMeterStpMix;
    DoubleVec TH4SiO4Mix;
    DoubleVec TNH4Mix;
    DoubleVec TH3BO3Mix;

    DoubleVec SampleIDMix;
    std::string SampleDateMix;
    DoubleVec OperatorMix;
    std::string WellNameMix;
    DoubleVec FieldMix;
    std::string StateMix;

    DoubleVec VgTPMix;
    DoubleVec VoMix;
    DoubleVec VwMix;
    DoubleVec VMeOHMix;
    DoubleVec VMEGMix;
    DoubleVec oilAPIgravMix;
    DoubleVec gasSpGravMix;
    DoubleVec MixFracGas;

    DoubleVec nTCO2Mix;
    DoubleVec nTCH4Mix;
    DoubleVec nTH2SMix;
    DoubleVec mass_w_Mix;
    DoubleVec mass_o_Mix;
    DoubleVec MixFracOil;

    DoubleVec mass_MeOH_mix;
    DoubleVec mass_MEG_mix;
    DoubleVec Qheat;
    DoubleVec yCO2Mix;
    DoubleVec yH2SMix;
    DoubleVec yCH4Mix;

    DoubleVec YCH4stpmix;
    DoubleVec RatioOilBPointsmix;
    DoubleVec CalculatedTDSMix;
    DoubleVec rho25CMix;

    DoubleVec HstpMix;
    DoubleVec OHstpMix;
    DoubleVec HCO3stpMix;
    DoubleVec CO3stpMix;
    DoubleVec ACstpMix;
    DoubleVec HSstpMix;
    DoubleVec NH4STPMix;
    DoubleVec H2BO3stpMix;

    DoubleVec HCO3AlkMix;
    DoubleVec CO3AlkMix;
    DoubleVec HAlkMix;
    DoubleVec OHAlkMix;
    DoubleVec ConcFactor;

    DoubleVec TCO2Mix;
    DoubleVec TofpH;
    DoubleVec PofpH;
    DoubleVec TofVol;
    DoubleVec PofVol;

    DoubleVec OilDensityMix;
    DoubleVec GasDensityMix;
    DoubleVec WaterDensityMix;
    DoubleVec UseTPpHMix;
    DoubleVec UseTPVolMix;
    DoubleVec useEOSmix;

    DoubleVec molAlk;
    DoubleVec molTAC;
    DoubleVec molTNH4;
    DoubleVec molTH3BO3;
    DoubleVec molTH2Saq;
    DoubleVec molTH4SiO4;
    DoubleVec mol_g_origMix;
    DoubleVec mol_o_OrigMix;
    DoubleVec mol_w_OrigMix;
    DoubleVec mol_g_finalMix;
    DoubleVec mol_o_finalMix;
    DoubleVec mol_w_finalMix;
    DoubleVec mol_w_evapMix;

    DoubleVec Total_molesMix;
    DoubleVec SumofZMix;
    DoubleVec nTCO2MixEOS;
    DoubleVec nTH2SMixEOS;

    // edit by hzy - 改为int类型
    IntVec usepHmix;

    DoubleMatrix zMix;

    double MaxInh;
    double NoRiskcalc;
    DoubleVec ConcInhCelRisk;
    DoubleVec ConcInhCalRisk;
    DoubleVec ConcInhBarRisk;
    DoubleVec ConcInhAnRisk;
    DoubleVec ConcInhGypRisk;

    DoubleMatrix NaWI;
    DoubleMatrix CaWI;
    DoubleMatrix SrWI;
    DoubleMatrix BaWI;
    DoubleMatrix ClWI;
    DoubleMatrix SO4WI;
    DoubleMatrix HCO3AlkWI;
    DoubleMatrix CO3AlkWI;
    DoubleMatrix TACWI;
    DoubleMatrix YCO2WI;
    DoubleMatrix YH2SWI;
    DoubleMatrix TH2SaqWI;
    DoubleMatrix pHMeterSTPWI;
    DoubleMatrix ConcFactorWI;

public:
    int nob = 1;
    /*下面这三个变量原先是const  我先去掉const   -彭非*/
    int nob_Input = 1;
    int nob_InputII = 0;
    int Read_InputII = 0;
    int Run1000Cases;
    IntVec CaseCount;

    int LoopTP1000Cases;

    int j;
    int kk;

    StrVec InhName;
    int Run_CalcConcFactor;
    double TBH;
    double TWH;
    double PBH;
    double PWH;
    int RunMultiMix;
    double rho25c = 1.0;       // 25°C密度

    double tInh;

    // const 类型，控制参数，可以考虑放入控制类中   （在A1种会被赋值）
    int RunStat = 0;
    int RunH2SGUI = 0;
    int RunQualityControlChecks = 0;

    const int RunStatMix = 0;
    const double RunWhatIf = -1; // vb 中没有赋初值
    const double TWIInit = -1;   // vb 中没有赋初值
    const double UseSI = 0;    // vb 中没有赋初值，但是判断UseSI = 0会进入，暂时 表示0, 在InputPartB中会影响simContext.VgTPMix
    const double PWIInit = -1;       // vb 中没有赋初值
    double SIPyrr = 0; // vb 中没有赋初值

    int Run10TestCases = 0;
    int Loop10;
    int Run_Seawater_Mixing;
    int LoopMixing;
    int LoopResChem;
    int Run_MixingTwoWells;
    int Run_MassTransfer_WhatIf = 0;

    int UseTPCalciteSheet;
    double useTPVol;//   - 修改  by彭非 //int simContext.useTPVol;
    IntVec CaseCountWI;
    int Loop1WI;
    int Loop2WI;

    int OutPutActCoefs = 1;

    int UseMolal;

    IntVec errmsg;

    double yCO2;
    double CO2aq;
    double HCO3;
    double CO3;
    double yH2S;
    double HS;
    double H2Saq;

    int useSR = 0;

    // 源代码中未说明类型，按逻辑推测为int
    int RunNORM;

    // 在单组份计算中还未使用的变量
    double SelectInh;
    double InhNo;
    double InhNo1;
    double FracInhNo1;
    double InhNo2;
    int ING;

public:

};
SimulationContext simContext;

// ==================== Excel操作接口类（抽象层）====================
class ExcelWriter {
public:
    virtual ~ExcelWriter() = default;

    virtual void writeCell(const std::string& sheetName, int row, int col, double value) = 0;
    virtual void writeCell(const std::string& sheetName, int row, int col, const std::string& value) = 0;
    virtual void writeCell(const std::string& sheetName, int row, int col, int value) = 0;
    virtual void writeCell(const std::string& sheetName, int row, int col, bool value) = 0;
    virtual void clearCell(const std::string& sheetName, int row, int col) = 0;
    virtual void setRange(const std::string& sheetName, const std::string& range, const std::string& value) = 0;

protected:
    // 行/列计算辅助方法
    int calculateRow(int baseRow, int loopTP, int offset = 0) {
        return baseRow - loopTP + offset;
    }

    int calculateRowWithStat(int baseRow, int statTPnob, int loopTP, int offset = 0) {
        return baseRow + statTPnob - loopTP + offset;
    }

    int calculateRowWithIter(int baseRow, int loopTP, int iter, int offset = 0) {
        return baseRow - loopTP + (iter - 1) * 32 + offset;
    }
};

// ==================== 数据存储结构 ====================
struct OutputData {
    // 按工作表组织的数据结构
    struct WorksheetData {
        std::string name;
        std::map<std::pair<int, int>, double> numericData;      // (row, col) -> value
        std::map<std::pair<int, int>, std::string> textData;    // (row, col) -> text
        std::map<std::pair<int, int>, bool> booleanData;        // (row, col) -> boolean
    };

    std::map<std::string, WorksheetData> worksheets;

    void addData(const std::string& sheet, int row, int col, double value) {
        worksheets[sheet].numericData[{row, col}] = value;
    }

    void addData(const std::string& sheet, int row, int col, const std::string& value) {
        worksheets[sheet].textData[{row, col}] = value;
    }

    void addData(const std::string& sheet, int row, int col, bool value) {
        worksheets[sheet].booleanData[{row, col}] = value;
    }
};



// ==================== 辅助函数实现 ====================
int calculateRowWithStat(int baseRow, int statTPnob, int loopTP, int offset = 0) {
    return baseRow + statTPnob - loopTP + offset;
}

// ==================== ExcelWriter的具体实现示例 ====================
// 注意：以下是几种可能的实现方案

// 方案1：使用CSV文件输出（跨平台）
class CSVExcelWriter : public ExcelWriter {
private:
    std::string basePath;
    std::map<std::string, std::ofstream> fileStreams;

public:
    CSVExcelWriter(const std::string& path) : basePath(path) {}

    ~CSVExcelWriter() {
        for (auto& stream : fileStreams) {
            if (stream.second.is_open()) {
                stream.second.close();
            }
        }
    }

    void writeCell(const std::string& sheetName, int row, int col, double value) override {
        std::string filename = basePath + "/" + sheetName + ".csv";
        ensureFileOpen(filename);

        // CSV格式：行, 列, 值
        fileStreams[filename] << row << "," << col << "," << value << std::endl;
    }

    void writeCell(const std::string& sheetName, int row, int col, const std::string& value) override {
        std::string filename = basePath + "/" + sheetName + ".csv";
        ensureFileOpen(filename);

        // 处理可能包含逗号的字符串
        std::string escapedValue = value;
        escapedValue.replace(escapedValue.begin(), escapedValue.end(), ',', ';'); // 替换引号
        fileStreams[filename] << row << "," << col << ",\"" << escapedValue << "\"" << std::endl;
    }

    void writeCell(const std::string& sheetName, int row, int col, int value) override {
        writeCell(sheetName, row, col, static_cast<double>(value));
    }

    void writeCell(const std::string& sheetName, int row, int col, bool value) override {
        writeCell(sheetName, row, col, value ? "TRUE" : "FALSE");
    }

    void clearCell(const std::string& sheetName, int row, int col) override {
        // CSV中清空单元格可以写空值
        writeCell(sheetName, row, col, "");
    }

    void setRange(const std::string& sheetName, const std::string& range, const std::string& value) override {
        // 简化实现：将范围视为单个单元格
        // 实际应用中需要解析范围字符串（如"A1:B10"）
        std::cout << "CSVWriter: setRange not fully implemented for range: " << range << std::endl;
    }

private:
    void ensureFileOpen(const std::string& filename) {
        if (fileStreams.find(filename) == fileStreams.end() || !fileStreams[filename].is_open()) {
            fileStreams[filename].open(filename, std::ios::app);
            if (!fileStreams[filename].is_open()) {
                std::cerr << "无法打开文件: " << filename << std::endl;
            }
        }
    }
};

const double pi = 3.14159265358979;    // π (圆周率)

// 气体常数（不同单位）
const double RBar = 0.083144;           // Gas constant (L·bar/(K·mol))// 气体常数，单位：bar·m³/(kmol·K)
const double R = 83.144;                // Gas constant (cm³·bar/(K·mol))
const double RAtm = 0.082057;            // Gas constant (L·atm/(K·mol))

//-------------------宏定义  以下------
//下面的define索引进行偏移
// Cation indexes: 阳离子索引（宏定义，编译时替换为对应整数）
enum class CationIndex : size_t {
    iH = 0,
    iNa,
    iK,
    iMg,
    iCa,
    iSr,
    iBa,
    iFe,
    iZn,
    iPb,
    iNH4,
    iRa,
    COUNT
};
#define iH    0
#define iNa   1
#define iK    2
#define iMg   3
#define iCa   4
#define iSr   5
#define iBa   6
#define iFe   7
#define iZn   8
#define iPb   9
#define iNH4  10
#define iRa   11

// Anion indexes: 阴离子索引（宏定义，对应数组下标）
enum class AnionIndex : size_t {
    iOH = 0,
    iCl,
    iAc,
    iHCO3,
    iCO3,
    iSO4,
    iHS,
    iF,      // 将 intF 改为更有意义的名称
    iBr,
    iH2BO3,
    iH3SiO4,
    iH2SiO4,
    iSion,
    COUNT
};
#define iOH     0
#define iCl     1
#define iAc     2
#define iHCO3   3
#define iCO3    4
#define iSO4    5
#define iHS     6
#define intF    7
#define iBr     8
#define iH2BO3  9
#define iH3SiO4 10
#define iH2SiO4 11
#define iSion   12



// 'Neutral aquatic indexes
enum class NeutralAquaticIndex : size_t {
    iCH4aq = 0,
    iCO2aq,
    iH2Saq,
    iHAcaq,
    iH4SiO4aq,
    iNH3,
    iH3BO3,
    iFeSaq,
    COUNT
};
#define iCH4aq     0
#define iCO2aq     1
#define iH2Saq     2
#define iHAcaq     3 
#define iH4SiO4aq  4
#define iNH3       5 
#define iH3BO3     6
#define iFeSaq     7

//// Oil phase indexes
enum class OilPhaseIndex : size_t {
    iCH4o = 0,
    iCO2o,
    iH2So,
    COUNT
};
#define iCH4o  0 
#define iCO2o  1
#define iH2So  2

/* bDot物质种类索引定义（用于标识不同络合物种或离子的编号） */
enum class BDotSpeciesIndex : size_t {
    iZnDot = 0,
    iPbDot,
    iHSDot,
    iClDot,
    iZnCl,
    iZnCl2,
    iZnCl3,
    iZnCl4,
    iZnHS2,
    iZnHS3,
    iPbCl,
    iPbCl2,
    iPbCl3,
    iPbCl4,
    iPbHS2,
    iPbHS3,
    iOHDot,
    iZnOH,
    iZnOH2,
    COUNT
};
#define iZnDot   0	// Zn²⁺离子（锌离子）
#define iPbDot   1	// Pb²⁺离子（铅离子）
#define iHSDot   2	// HS⁻离子（硫氢根离子）
#define iClDot   3	// Cl⁻离子（氯离子）
#define iZnCl    4	// ZnCl⁺络离子（一氯合锌离子）
#define iZnCl2   5	// ZnCl₂中性分子（二氯合锌）
#define iZnCl3   6	// ZnCl₃⁻络离子（三氯合锌离子）
#define iZnCl4   7	// ZnCl₄²⁻络离子（四氯合锌离子）
#define iZnHS2   8	// Zn(simContext.HS)₂中性分子（二硫氢合锌）
#define iZnHS3   9	// Zn(simContext.HS)₃⁻络离子（三硫氢合锌离子）
#define iPbCl    10  // PbCl⁺络离子（一氯合铅离子）
#define iPbCl2   11	 // PbCl₂中性分子（二氯合铅）
#define iPbCl3   12  // PbCl₃⁻络离子（三氯合铅离子）
#define iPbCl4   13	 // PbCl₄²⁻络离子（四氯合铅离子）
#define iPbHS2   14  // Pb(simContext.HS)₂中性分子（二硫氢合铅）
#define iPbHS3   15  // Pb(simContext.HS)₃⁻络离子（三硫氢合铅离子）
#define iOHDot   16
#define iZnOH    17
#define iZnOH2   18

// 气体组分索引定义
enum class GasComponentIndex : size_t {
    iCH4g = 0,
    iCO2g,
    iH2Sg,
    iC2g,
    iC3g,
    iC4ig,
    iC4ng,
    iC5ig,
    iC5ng,
    iC6g,
    iC7_12g,
    iC13_25g,
    iC26_80g,
    iN2g,
    iH2Og,
    COUNT
};
#define iCH4g    0
#define iCO2g    1
#define iH2Sg    2
#define iC2g     3
#define iC3g     4
#define iC4ig    5
#define iC4ng    6
#define iC5ig    7
#define iC5ng    8
#define iC6g     9
#define iC7_12g  10
#define iC13_25g 11
#define iC26_80g 12
#define iN2g     13
#define iH2Og    14

// 矿物和沉淀物索引定义
enum class MineralIndex : size_t {
    iBaSO4 = 0,
    iCaSO42H2O,
    iSrSO4,
    ihemiCaSO4,
    iCaSO4,
    iNaCl,
    iCaHCO32,
    iFeHCO32,
    iFeHS2,
    COUNT
};
#define iBaSO4      0
#define iCaSO42H2O  1
#define iSrSO4      2
#define ihemiCaSO4  3
#define iCaSO4      4
#define iNaCl       5
#define iCaHCO32    6
#define iFeHCO32    7
#define iFeHS2      8
;
//------------------------    defind   分割线 -------------------

double VMeOH, VMEG;// 仅在PartD赋值，但从未见使用过，不明所以的变量

double VW, VgTP, VO, mass_MeOH, mass_MEG;
double yCH4;
double SGG;
double Ppsia, TF, TC;

double mtotal;// - 局部变量，在CalcIonicStrength计算，然后C2_PitzerActCoefs_T_P_ISt会调用计算好的它    可以改成传参

double MoleCharge; // - 局部变量，在CalcIonicStrength计算，然后C2_PitzerActCoefs_T_P_ISt会调用计算好的它  

double SumOfCations, SumOfAnions; // 仅在CalcIonicStrength、QualityControlCalculations使用

double RatioOilBPoints;

double ISt; double DpHj;

int NumCat = 12, NumAn = 13, NumNeut = 8;//final不变量

double KgwCO2 = 0, K1H2CO3 = 0, KHAc = 0, K2HCO3 = 0, KH2O = 0, K1H2S = 0, K2HS = 0, KNH4 = 0, KH3BO3 = 0;
double KgoCO2 = 0, KgoH2S = 0;  //A13_CalcH2SFugacity、C1_ThermodynamicEquilConsts、C4_SSPEquilCalcs、fTotalCO2H2Smoles
double KgwCH4 = 0, KgoCH4 = 0; // C1_ThermodynamicEquilConsts、C4_SSPEquilCalcs、fTotalCO2H2Smoles
double KgwH2S = 0;         //C1_ThermodynamicEquilConsts、C4_SSPEquilCalcs、C4_EOS_TCO2_SSPEquilCalcs、C5_CalcpHPCO2PH2SSTP、QualityControlCalculations、fTotalCO2H2Smoles

// 平衡常数
double KspCalcite = 0, KspBarite = 0, KspCelestite = 0;
double KspHemihydrate = 0, KspGypsum = 0, KspAnhydrite = 0, KspSiderite = 0, KspHalite = 0;

double IStCosolvent = 0;

//------------------未记录生命周期----------------
double KspFeS = 0;//确认仅在C1被赋值，因此在C1首次调用
double KspZnS = 0;//确认仅在C1被赋值，因此在C1首次调用
double KspPbS = 0;//确认仅在C1被赋值，因此在C1首次调用
double KspFeSAm = 0, KspTrot = 0, KspPyrr = 0, KspCaF2 = 0, KspDol = 0;//确认仅在C1被赋值，因此在C1首次调用
double xMeOH = 0; // 多处被赋值多处被调用
double KspZnCO3 = 0, KspCaOH2 = 0, KspMgOH2 = 0, KspSrCO3 = 0, KspBaCO3 = 0, KH4SiO4 = 0, KH3SiO3 = 0;//确认仅在C1被赋值，因此在C1首次调用
double KspChrysotile = 0, KspDiopside = 0, KspGreenalite = 0, KspQuartz = 0, KspAmSilica = 0;//确认仅在C1被赋值，因此在C1首次调用

// 镭共沉淀参数
double RaBaSO4ppt = 0.0;
double RaCaSO4ppt = 0.0;
double RaSrSO4ppt = 0.0;
double DRaBarite = 0.0, DRaCelestite = 0.0, DRaAnhydrite = 0.0, GammaSolidRaBarite = 0.0, GammaSolidRaCelestite = 0.0, GammaSolidRaAnhydrite = 0.0;//确认仅在C1被赋值，因此在C1首次调用，生命周期似乎很简单
double KstFeSaq = 0.0;//确认仅在C1被赋值，因此在C1首次调用

double BetaDot[20] = { 0 };// 

double xMEG = 0;
double aH2O = 0;

double zOutput[15] = { 0 }; DoubleVec z(20, 0);
DoubleVec density(3, 0);
double ChCat[15] = { 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 1, 2, 0, 0, 0 }; double ChAn[15] = { -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -1, -2, -2, 0 ,0 };

double b0[15][15] = { 0 }; double b1[15][15] = { 0 }; double b2[15][15] = { 0 };
double CPhi[15][15] = { 0 }; double Lnc[15][15] = { 0 }; double Lna[15][15] = { 0 };
double zeta[15][15][15] = { 0 }; double Taap[15][15] = { 0 };
double Yaapc[15][15][15] = { 0 }; double Yccpa[15][15][15] = { 0 };
double Tccp[15][15] = { 0 };

double V0_c[15] = { 0 }; double V0_a[15] = { 0 }; //B1_InitializeIndices赋值、V0TP赋值、CalcRhoTP调用

double V0_n[10] = { 0 };//V0_n在B1中赋值，仅在CalcRhoTP中使用


//---------------以下4个数组fgammaN调用的东西在B1中赋值
double gNNeut[10] = { 1,1,1,1,1,1,1,1,0,0 };
double gNAn[15] = { 1,1,1,1,1,1,1,1,1,1,1,1,1,0,0 };
double gNCat[15] = { 1,1,1,1,1,1,1,1,1,1,1,1,0,0,0 };
double gNMean[10] = { 1,1,1,1,1,1,0,0,0,0 };

//---------------
double aNH2O = 0;
//double gCat[15] = { 0 };
//double gAn[15] = { 0 };
double gDot[20] = { 0 };
//double gGas[15] = { 0 };

double MWCat[15] = { 0 }, MWAn[15] = { 0 }, MWNeut[10] = { 0 };

double z_before_precipitation[15] = { 0 };


double Mass_o = 0, mass_w = 0;
double nTCO2 = 0, Znew = 0;
double nTCO2EOS = 0, nTCH4 = 0, nTH2S = 0, nTH2sEOS = 0, aH = 0;
double pHMeterReading = 0;
double total_moles;

//simContext.Run10TestCases = 1, 
int Run_MassTransfer = 0, RunShellMultiflash = 0;



bool mf_ParametersWereRead;
//   以下几个数组，只会被用在MultiPhaseFlash，并且MultiPhaseFlash调用InitialPreparationSSP对它们进行初始化
//  下面这几个量管理起来非常困难， 其定义、利用、释放对MultiPhaseFlash函数造成重大干扰
DoubleVec mf_TCr;
DoubleVec mf_PCr;
DoubleVec mf_Omega;
DoubleVec mf_MWgas;
DoubleVec mf_c0;
DoubleVec mf_c1;
DoubleMatrix mf_kPr;

// 下面这几个变量长度应该是固定的，但是需要弄成指针形式ReDim density(3), compositions(15, 4), phi(15, 3), Compr(3), beta(3), zOutput(15)
DoubleVec mass_phase; //3
DoubleMatrix compositions;  //15*4
DoubleMatrix phi;  //15*3
DoubleVec beta; //3
DoubleVec MW_Phase; //3
DoubleVec Compr;    //3 class.beta

int No_Phases = 0;
double SumofZ = 0;

double nTCO2_before_precipitation, nTH2S_before_precipitation, Total_moles_before_precipitation;

DoubleVec lnphi_Gas; //这个和pseudo_composition、phi_calc、true_composition有关

//----C4
double ppt;
double Vg_ZRT; //仅在C4和B5_CalculateSIvalues使用过
double pHHigh, pHLow;
double H, OH;
int Iteration;
int useEOS;
double t1, t2, t3;
double root1, root2, root3;
double a, b, cc;
double total_moles_Temp;  //用在C4_SSPEquilCalcs、C4_EOS_TCO2_SSPEquilCalcs、fEOS_Speciation
DoubleVec zTemp(15, 0);//用在C4_SSPEquilCalcs、C4_EOS_TCO2_SSPEquilCalcs、fEOS_Speciation   // 长度15

char* myiosheet, * myiocol;  // C4 -> MultiPhaseFlash_CS ->  Flash_Input_Processing
double eosProps[15][6];//ReDim eosProps(nComponents, 1 To 6)      - B1显示，nComponents=15

double kij[15][15];//kij(nComponents, nComponents)
double ypH[30], xpH[30], xppt[30];//ypH(30), xpH(30), xppt(30)

double PCO2, PH2S, PCH4;
double TH4SiO4_Greenalite, TH4SiO4_Dipside, TH4SiO4_Chrysotile; // 出现在C4_SSPEquilCalcs、C4_EOS_TCO2_SSPEquilCalcs  两个函数
double hydHS;
double S, hydAc, AC, HAcaq, hydH2BO3, H2BO3, hydNH3, NH3, hydH2SiO4, H2SiO4, H3SiO4, H4SiO4;
double faH;

//TrueFlash新加入的全局变量
DoubleVec E;

//----QualityControlCalculations
int UseH2Sgas;

//----PartD
double rhoTP;
int RunQualityControlChecks_II;
//molc(NumCat, nob_Input + nob_InputII),mola(NumAn, nob_Input + nob_InputII), moln(NumNeut, nob_Input + nob_InputII)
DoubleMatrix molc; DoubleMatrix mola; DoubleMatrix moln;

//全局变量 -- by 黄志缘
/********************************************************************************************************/

// 物理常数
#define NAv         6.0221367E+23        // 阿伏伽德罗常数 (mol⁻¹)
#define eElec       1.60217733E-19       // 电子电荷 (C)
#define eps0        8.854187818E-12      // 真空介电常数 (F/m)
#define kBoltz      1.380658E-23         // 玻尔兹曼常数 (J/K)

double Lnn[16][16], a0[16];
double Alk = 0.0;          // 碱度
double TAc = 0.0;          // 总醋酸浓度  
double TH2Saq = 0.0;       // 总H2S(aq)浓度
double TH4SiO4 = 0.0;      // 总H4SiO4浓度
double TH3BO3 = 0.0;       // 总H3BO3浓度
double TNH4 = 0.0;         // 总NH4浓度
double TFe = 0.0;          // 总Fe浓度
double TPb = 0.0;          // 总Pb浓度
double TZn = 0.0;          // 总Zn浓度
double TCO2 = 0.0;         // 总CO2浓度
double TDS = 0.0;          // 总溶解固体
double TCa = 0.0;          // 总Ca
double CalculateTDSDen = 0.0; // 计算TDS密度

// C5_CalcpHPCO2PH2SSTP 变量
// double KgwCO2;
// double KgwH2S;
int use_pH;
//int UseH2Sgas,useEOS;    -已定义
// double *pHMeterReading;

// D2_CalcDensitypH
double rhoOld, rhoSSE;
//int Iteration;    -已定义

// D1_CalcDensity
//double  H, OH, AC, simContext.HS, NH3, H2BO3, simContext.HCO3, simContext.CO3, H3SiO4, H2SiO4, H4SiO4;  -全部已定义
//double simContext.CO2aq, simContext.H2Saq, HAcaq;  -全部已定义
double TDSSSE, TDSOld;

const int MaxMix = 100, MaxSpecies = 13, MaxComponents = 100, MaxCat = 100, MaxNeut = 100, MaxAnion = 100;

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

double SArea = 0.0;        // 表面积
double QBrineFlow = 0.0;   // 盐水流量
DoubleVec radiusC(11, 0.0);
DoubleVec radiusA(11, 0.0);

//PartB新加入的两个全局变量，目前单组份用不到，先保留
double PipeID = 0; // 管道内径
double PipeL = 0; // 管道长度


// -------------------------------------------------------------
// 12月11日后新增的全局变量 
// 新增的变量
double dSIMeOHBar;
double dSIMEGBar;
double dSIMeOHcal;
double dSIMEGcal;
double dSIMeOHHal;
double dSIMEGHal;

double rgas_water;
double tc_water;
double pc_water;
double dc_water;

int UseSR;
int H2Oevap;
int RunGoalSeek, RunStatGoalSeek;
int RunSMT, RunH2SPartition;
int usePTB;

int Ncount, NCount_II;
int RunStatReservoirCalc = 0;
int RunStatSICalcSSP = 0;

DoubleVec InhSolLimit(15, 0);
double QTotal;
//仅在B3和B5使用，我们可以设为局部的，然后通过参数列表传给B3和B5
DoubleVec mcInit(15);
DoubleVec maInit(15);
double AlkInit;
double TH4SiO4Init;
double TNH4Init;
double TH3BO3Init;
double TAcInit;
double TFeInit;
double mass_w_0;

int LoopTP;

// 饱和指数 (SI)
// double SIBar = 0.0;
// double SIGyp = 0.0;
// double SIHemi = 0.0;
// double SIAn = 0.0;
// double SICel = 0.0;
// double SIHal = 0.0;
// double SICaF2 = 0.0;
// double SICal = 0.0;
// double SISid = 0.0;
// double SIZnCO3 = 0.0;
// double SIFeSAm = 0.0;
// double SIFeS = 0.0;
// double SITrot = 0.0;
// double SIPyrr = 0.0;
// double SIZnS = 0.0;
// double SIPbS = 0.0;
// double SICaOH2 = 0.0;
// double SIMgOH2 = 0.0;
// double SIAmSilica = 0.0;
// double SIQuartz = 0.0;
// double SIGreenalite = 0.0;
// double SIDiopside = 0.0;
// double SIChrysotile = 0.0;
// double SISrCO3 = 0.0;
// double SIBaCO3 = 0.0;
double SICal, SIBar, SIGyp, SIHemi, SIAn, SICel, SIHal, SICaF2, SISid, SIFeS, SIZnS, SIBaCO3, SISrCO3;
double SIDol, SIFeSAm, SITrot, SIPbS, SIZnCO3, SIAmSilica, SIQuartz, SIDiopside;
double SIChrysotile, SIGreenalite;

//左边的可能是局部的，待优化
double SICalBh;
double SIBarBH;
double SIGypBH;
double SIHemiBH;
double SIAnBH;
double SICelBH;
double SIHalBH;
double SICaF2BH;
double SISidBH;
double SIFeSBH;
double SIZnSBH;
double SIBaCO3BH;
double SISrCO3BH;

double ConcInhCal, ConcInhBar;

double WHConcInhCal, WHConcInhBar;

double dSICal, dSIBar, dSIHal, dSISid, dSIFeS;

int RunSSP = 1;

// 沉淀量变量
double pptCal=0.0, pptBar=0.0, pptSid=0.0, pptGyp=0.0, pptHemi=0.0, pptAn=0.0, pptCel=0.0, pptHal=0.0;
double pptFeS=0.0, pptZnS=0.0, pptZnCO3=0.0, pptPbS=0.0, pptCaF2=0.0;
double pptMgOH2=0.0, pptCaOH2=0.0, pptFeSAm=0.0, pptTrot=0.0;
double pptQuartz=0.0, pptGreenalite=0.0, pptDiopside=0.0, pptChrysotile=0.0, pptAmSilica=0.0;
double pptSrCO3 = 0.0;
double pptBaCO3 = 0.0;

double pptFeCO3_NoMassTransfer, pptFeCO3_MassTransfer = 0;
double pptCalcite_MassTransfer, pptCalcite_NoMassTransfer = 0;
double pptBarite_MassTransfer, pptBarite_NoMassTransfer = 0;
double pptCelestite_MassTransfer, pptCelestite_NoMassTransfer = 0;
double pptFeS_MassTransfer, pptFeS_NoMassTransfer = 0;

// Dai 2020 deposition 沉积速率
double pptCalcite_MassTransfer_V = 0.0;
double pptBarite_MassTransfer_V = 0.0;
double pptFeCO3_MassTransfer_V = 0.0;

DoubleVec mcMT(11, 0); // 传质后的阳离子浓度
DoubleVec maMT(11, 0);  // 传质后的阴离子浓度
DoubleVec mcEq(11, 0);   // 平衡浓度
DoubleVec mcKin(11, 0);  // 动力学浓度

double TFe_MT; // 传质后的总铁
double AlkMT; // 传质后的碱度

double PBubblePt; double VgTP_MT, VO_MT, VW_MT;
double a1;

int QPBubblePt;

double SIMgOH2, SICaOH2;

DoubleVec CaseCount_II;
int Iter_MT_WI = 0;

DoubleMatrix bi(20, DoubleVec(20, 0));
DoubleMatrix ci(20, DoubleVec(20, 0));
DoubleMatrix gi(20, DoubleVec(20, 0));
DoubleMatrix ai(20, DoubleVec(20, 0));
DoubleMatrix celi(20, DoubleVec(20, 0));

DoubleVec bInhBar(20, 0);
DoubleVec EaRInh(20, 0);
DoubleVec LnAInh(20, 0);

DoubleVec bInhCal(20);
DoubleVec bInhGyp(20, 0);
DoubleVec bInhAn(20, 0);
DoubleVec bInhCel(20, 0);
double InhNoBar;    
double tInh;
int InhNoCal;
double t0Cal;
double t0Bar;
double t0An;
double t0Cel;
double t0Gyp;
int InhNoGyp;
int ConcInhAn;
int InhNoAn;
int InhNoCel;
double ConcInhCel;
double dens;
std::string Flow_Pattern;
std::string Flow_Regime;
double Eff_D, Eff_V;
DoubleVec DiffC(11, 0);
DoubleVec DiffA(12, 0);
DoubleVec Diff(10, 0);

// SIRisk[抑制剂编号][LoopTP][参数索引][物质类型]
// std::vector<std::vector<std::vector<DoubleVec>>> SIRisk;
std::vector<std::vector<std::vector<DoubleVec>>> SIRisk(
    20,  // 第一维大小
    std::vector<std::vector<DoubleVec>>(
        10,  // 第二维大小
        std::vector<DoubleVec>(
            10,  // 第三维大小
            DoubleVec(4)  // DoubleVec初始化为大小为4
        )
    )
);
// 井名数组
StrVec WellNameMix;

// 工作表名称
std::string myname;

// 传质系数数组
DoubleVec km;

double fSafetyBar;
double fSafetyCal;
double fsafetyAn;
double fSafetyCel;
double fSafetyGyp;

// 化学物质浓度数组
DoubleVec mc(15, 0);  // 阳离子浓度
DoubleVec ma(15, 0);  // 阴离子浓度
DoubleVec mn(15, 0);  // 中性离子浓度

// 活度系数数组
DoubleVec gCat(15, 0);   // 阳离子活度系数
DoubleVec gAn(15, 0);    // 阴离子活度系数
DoubleVec gNeut(15, 0); // 中性分子活度系数
DoubleVec gGas(15, 0);   // 气体活度系数
DoubleVec gL(15, 0);     // 油相活度系数


double ConcInhGyp;
double m_CR_selected = 0.0;
// double PBubblePt = 0.0;
int Use_Corr_in_Deposition = 0;


// pH变量
double pH=0.0;
double pH_before_precipitation=0.0;
double pHaftercalciteppt = 0.0;
double pHafterAmsilicappt = 0.0;
double pHafterQuartzppt = 0.0;
double pHafterGreenaliteppt = 0.0;
double pHafterDiopsideppt = 0.0;
double pHafterChrysotileppt = 0.0;
double pHafterMgOH2ppt = 0.0;
double pHafterCaOH2ppt = 0.0;
double pHafterFeCO3ppt = 0.0;
   

// 条件标志
int MultiplePpt = 0;     

// 最大沉积速率
double pptMax_MassTransfer_V = 0.0;
double pptTotal_MassTransfer_V = 0.0;

int StatTPnob = 0;
double location = 0.0;
double depth = 0.0;
// double Time_lsl_burst = 0.0;
double lgt_calcite_MT = 0.0;
double lgt_barite_MT = 0.0;
// double ViscWatIst = 0.0;
// double CpPerMl = 0.0;
// double pptCalcite_NoMassTransfer = 0.0;
// double pptBarite_NoMassTransfer = 0.0;
double TDSHalite = 0.0;
double Rho25cHalite = 0.0;
// double pptSrCO3 = 0.0;
// double pptBaCO3 = 0.0;

int i_Inhibitor_MT=0;
double c_Inhibitor_MT=0.0;
double Time_lsl_burst=0.0;
int Use_lsl_stability=0;
double PipeRough=0.0;
double ReNO = 0.0;
double Hold_l=0.0;
/********************************************************************************************************/

// 针对Shell Input Excel表格的样品数据结构体
struct SampleData
{

    /* ---------- sample_basic_information ---------- */
    char SampleID[256]; // 样品ID Excel row 4
    char Date[256];     // 采样日期 Excel row 5
    char Operator[256]; // 采样人员 Excel row 6
    char WellName[256]; // 油井名称 Excel row 7
    char Location[256]; // 所在位置 Excel row 8
    char Field[256];    // 油田名称 Excel row 9

    /* ---------- sample_composition_information ---------- */
    char SampleInfo[256]; // 样品信息（用户自由输入）
    double Na_aq;   // Execl row 10
    double K_aq;    // Excel row 11
    double Mg_aq;   // Excel row 12
    double Ca_aq;   // Excel row 13
    double Sr_aq;   // Excel row 14
    double Ba_aq;   // Excel row 15
    double FeII_aq; // Excel row 16
    double Zn_aq;   // Excel row 17
    double Pb_aq;   // Excel row 18
    double Cl_aq;   // Excel row 19
    double SO4_aq;  // Excel row 20
    double F_aq;    // Excel row 21
    double Br_aq;   // Excel row 22
    double Si_aq;   // Excel row 23
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
    double Alk_Bicarbonate_aq;  // Excel row 24
    double Alk_Carbonate_aq;    // Excel row 25
    double OrgAcid_Acetate_aq;  // Excel row 26
    double Ammonia_aq;          // Excel row 27
    double B_aq;                // Excel row 28
    double TDS_aq;              // Excel row 29
    double Density_STP;         // Excel row 30
    double CO2_pct_g;           // Excel row 31
    int Option_Use_H2Sg;        // Excel row 32
    double H2S_pct_g;           // Excel row 33
    double H2S_aq;              // Excel row
    double pH_STP;              // Excel row 34
    double Q_Gas;               // Excel row 35
    double Q_Oil;               // Excel row 36
    double Q_Water;             // Excel row 37

    /* ---------- sample_temperature_pressure_information ---------- */
    double T_initial;   // Excel row 39
    double T_final;     // Excel row 40
    double P_initial;   // Excel row 41
    double P_final;     // Excel row 42
    double API;         // Excel row 43
    double SG_g;        // Excel row 44
    double Q_MeOH;      // Excel row 45
    double Q_MEG;       // Excel row 46
    double StrongAcid_aq;   // Excel row 47
    double StrongBase_aq;   // Excel row 48
    double Conc_Multiplier; // Excel row 49
    double T_pH;            // Excel row 58
    double P_pH;            // Excel row 59
    double T_Q;             // Excel row 60
    double P_Q;             // Excel row 61

    /* ---------- sample_oil_phase_information ---------- */
    double C1_o;    // Excel row 66
    double CO2_o;   // Excel row 67
    double H2S_o;   // Excel row 68
    double C2_o;    // Excel row 69
    double C3_o;    // Excel row 70
    double iC4_o;   // Excel row 71
    double nC4_o;   // Excel row 72
    double iC5_o;   // Excel row 73
    double nC5_o;   // Excel row 74
    double C6_o;    // Excel row 75
    double C7_C12_o;// Excel row 76
    double C13_C25_o;// Excel row 77
    double C26_C80_o;// Excel row 78
    double N2_o;    // Excel row 79

    /* ---------- config_options_information ---------- */
    int Option_Alk;         // Excel row 51
    int Option_Defined_TP;  // Excel row 52
    int Option_TP_for_pH;   // Excel row 53
    int Option_TP_for_Q;    // Excel row 54
    int Option_EoS;         // Excel row 55
    int Option_Water_HC;    // Excel row 56

};

// 初始化测试数据函数
void mockData_sheetInput(SampleData* data)
{

    /* ---------- sample_basic_information ---------- */
    strcpy(data->SampleID, "smp01");
    strcpy(data->Date, "2025-01-01");
    strcpy(data->Operator, "张三");
    strcpy(data->WellName, "井-01");
    strcpy(data->Location, "大庆油田一区");
    strcpy(data->Field, "大庆油田");

    /* ---------- sample_composition_information ---------- */
    strcpy(data->SampleInfo, "常规采出水样");

    data->Na_aq = 19872.0;
    data->K_aq = 500.0;
    data->Mg_aq = 54.0;
    data->Ca_aq = 6500.0;
    data->Sr_aq = 700.0;
    data->Ba_aq = 550.0;
    data->FeII_aq = 12.0;
    data->Zn_aq = 10.0;
    data->Pb_aq = 1.0;
    data->Cl_aq = 43000.0;
    data->SO4_aq = 10.0;
    data->F_aq = 1.0;
    data->Br_aq = 10.0;
    data->Si_aq = 10.0;
    data->FeIII_aq = 0.1;
    data->Li_aq = 0.03;
    data->Be_aq = 0.01;
    data->Ra_aq = 0.002;
    data->Mn_aq = 0.4;
    data->Cu_aq = 0.06;
    data->Al_aq = 0.7;
    data->P_aq = 0.9;
    data->I_aq = 0.04;
    data->U_aq = 0.005;
    data->Alk_Bicarbonate_aq = 281.0;
    data->Alk_Carbonate_aq = 0.0;
    data->OrgAcid_Acetate_aq = 0.0;
    data->Ammonia_aq = 0.0;
    data->B_aq = 0.0;
    data->TDS_aq = 70000.0;
    data->Density_STP = 1.05;
    data->CO2_pct_g = 0.28;
    data->Option_Use_H2Sg = 0;
    data->H2S_pct_g = 4.32;
    data->H2S_aq = 1.2;
    data->pH_STP = 7.61618544407364;
    data->Q_Gas = 8500.0;
    data->Q_Oil = 1000.0;
    data->Q_Water = 100.0;

    /* ---------- sample_temperature_pressure_information ---------- */
    data->T_initial = 340;
    data->T_final = 140.0;
    data->P_initial = 7000.0;
    data->P_final = 400.0;
    data->API = 40.0;
    data->SG_g = 0.80;
    data->Q_MeOH = 0.0;
    data->Q_MEG = 0.0;
    data->StrongAcid_aq = 0.0;
    data->StrongBase_aq = 0.0;
    data->Conc_Multiplier = 1.0;
    data->T_pH = 0.0;
    data->P_pH = 0.0;
    data->T_Q = 340.0;
    data->P_Q = 7000.0;

    /* ---------- sample_oil_phase_information ---------- */
    data->C1_o = 74.16;
    data->CO2_o = 0.28;
    data->H2S_o = 0.013;
    data->C2_o = 7.9;
    data->C3_o = 4.15;
    data->iC4_o = 0.71;
    data->nC4_o = 1.44;
    data->iC5_o = 0.53;
    data->nC5_o = 0.66;
    data->C6_o = 0.81;
    data->C7_C12_o = 4.04;
    data->C13_C25_o = 1.45897808911913;
    data->C26_C80_o = 0.189774787571497;
    data->N2_o = 0.6;

    /* ---------- config_options_information ---------- */
    data->Option_Alk = 0;
    data->Option_Defined_TP = 0;
    data->Option_TP_for_pH = 0;
    data->Option_TP_for_Q = 0;
    data->Option_EoS = 0;
    data->Option_Water_HC = 0;
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


void initData()
{
    //int len = simContext.nob_Input + simContext.nob_InputII;

    //zMix = (double**)malloc((simContext.nob_Input + simContext.nob_InputII) * sizeof(double*));
    //for (int i = 0; i < (simContext.nob_Input + simContext.nob_InputII); i++)
    //{
    //    zMix[i] = (double*)malloc(15 * sizeof(double));
    //}
}

void pointerInit_pf() {
    mass_phase = DoubleVec(3, 0);
    beta = DoubleVec(3, 0);
    MW_Phase = DoubleVec(3, 0);
    Compr = DoubleVec(3, 0);

    compositions = DoubleMatrix(15, DoubleVec(4, 0));
    phi = DoubleMatrix(15, DoubleVec(3, 0));
    //----PartD
    molc.resize(NumCat, DoubleVec(simContext.nob_Input + simContext.nob_InputII, 0));
    mola.resize(NumAn, DoubleVec(simContext.nob_Input + simContext.nob_InputII, 0));
    moln.resize(NumNeut, DoubleVec(simContext.nob_Input + simContext.nob_InputII, 0));
}


void ReadInputPartA(int kk, SampleData* data)
{
    if (simContext.Run10TestCases == 1 && simContext.Loop10 > 1)
    {
        return;
    }
    if (simContext.Run_Seawater_Mixing == 1 && simContext.LoopMixing > 1)
    {
        return;
    }
    if (simContext.Run_MixingTwoWells == 1 && simContext.LoopMixing > 1)
    {
        return;
    }
    if (simContext.RunMultiMix == 1 && simContext.LoopResChem > 1)
    {
        return;
    }

    // 源代码中的这个逻辑为判断 Gas/Day 的单位，但是我们在数据库中没有设置这个单位的字段
    // If Worksheets(mySheet).Cells(35, 2).Value = "(kSm^3/D)" Then simContext.UseSI = 1

    // 源代码中的这个逻辑为判断 Na+ 的单位，但是我们在数据库中没有设置这个单位的字段
    // If Worksheets(mySheet).Cells(10, 2).Value = "(m)" Then simContext.UseMolal = 1

    //SampleDateMix[kk] = Date;
    //WellNameMix[kk] = WellName;
    simContext.NaMix[kk] = data->Na_aq;
    simContext.KMix[kk] = data->K_aq;
    simContext.MgMix[kk] = data->Mg_aq;
    simContext.CaMix[kk] = data->Ca_aq;
    simContext.SrMix[kk] = data->Sr_aq;
    simContext.BaMix[kk] = data->Ba_aq;
    simContext.FeMix[kk] = data->FeII_aq;
    simContext.ZnMix[kk] = data->Zn_aq;
    simContext.PbMix[kk] = data->Pb_aq;
    simContext.ClMix[kk] = data->Cl_aq;
    simContext.SO4Mix[kk] = data->SO4_aq;
    simContext.FMix[kk] = data->F_aq;
    simContext.BrMix[kk] = data->Br_aq;
    simContext.TH4SiO4Mix[kk] = data->Si_aq;              // TH4SiO4Mix(kk) = Worksheets(mySheet).Cells(23, j + 2).Value 其中23行对应Silica
    simContext.HCO3AlkMix[kk] = data->Alk_Bicarbonate_aq; // HCO3AlkMix(kk) = Worksheets(mySheet).Cells(24, j + 2).Value 其中24行对应Total Alkalinity
    simContext.CO3AlkMix[kk] = data->Alk_Carbonate_aq;    // CO3AlkMix(kk) = Worksheets(mySheet).Cells(25, j + 2).Value 其中25行对应 simContext.CO3 Alkalinity
    simContext.TAcMix[kk] = data->OrgAcid_Acetate_aq;     // simContext.TAcMix(kk) = Worksheets(mySheet).Cells(26, j + 2).Value 其中26行对应Carboxylates
    simContext.TNH4Mix[kk] = data->Ammonia_aq;            // TNH4Mix(kk) = Worksheets(mySheet).Cells(27, j + 2).Value其中27行对应 Ammonia
    simContext.TH3BO3Mix[kk] = data->B_aq;                // TH3BO3Mix(kk) = Worksheets(mySheet).Cells(28, j + 2).Value其中28行对应Borate;
    simContext.yCO2Mix[kk] = data->CO2_pct_g / 100;       // yCO2Mix(kk) = Worksheets(mySheet).Cells(31, j + 2).Value / 100 其中31行对应 Co2 Gas Analysis
    simContext.UseH2SgasMix[kk] = data->Option_Use_H2Sg;  // simContext.UseH2SgasMix(kk) = Worksheets(mySheet).Cells(32, j + 2).Value 其中32行对应 Use H2S Gas Analysis
    if (simContext.UseH2SgasMix[kk] == 1)
    {
        simContext.yH2SMix[kk] = data->H2S_pct_g / 100; // yH2SMix(kk) = Worksheets(mySheet).Cells(33, j + 2).Value / 100 其中33行对应Gas H2S% or simContext.H2Saq
        simContext.TH2SaqMix[kk] = 0;
    }
    else
    {
        simContext.TH2SaqMix[kk] = data->H2S_pct_g; // TH2SaqMix(kk) = Worksheets(mySheet).Cells(33, j + 2).Value
        simContext.yH2SMix[kk] = 0;
    }
    simContext.pHMeterStpMix[kk] = data->pH_STP;       // pHMeterStpMix(kk) = Worksheets(mySheet).Cells(34, j + 2).Value其中34行对应pH2 measured;
    simContext.HAlkMix[kk] = data->StrongAcid_aq;      // HAlkMix(kk) = Worksheets(mySheet).Cells(47, j + 2).Value其中47行对应H+(Strong acid);
    simContext.OHAlkMix[kk] = data->StrongBase_aq;     // OHAlkMix(kk) = Worksheets(mySheet).Cells(48, j + 2).Value其中48行对应OH-(Strong base);
    simContext.ConcFactor[kk] = data->Conc_Multiplier; // ConcFactor(kk) = Worksheets(mySheet).Cells(49, j + 2).Value其中49行对应Conc. multiplier;
    simContext.usepHmix[kk] = data->Option_Alk;        // usepHmix(kk) = Worksheets(mySheet).Cells(51, j + 2).Value其中51行对应Four Cal options;
    if (simContext.usepHmix[kk] == 0)            // If usepHmix(kk) = Empty Then usepHmix(kk) = 0这里的Empty设为0
    {
        simContext.usepHmix[kk] = 0;
    }
    simContext.UseTPpHMix[kk] = data->Option_TP_for_pH; // UseTPpHMix(kk) = Worksheets(mySheet).Cells(53, j + 2).Value其中53行对应T,P for pH;
    simContext.UseTPVolMix[kk] = data->Option_TP_for_Q; // UseTPVolMix(kk) = Worksheets(mySheet).Cells(54, j + 2).Value其中54行对应 T,P for G/O/W;
    if (simContext.nob_Input + simContext.nob_InputII == 1)
    {
        simContext.UseTPCalciteSheet = data->Option_Defined_TP; // simContext.UseTPCalciteSheet = Worksheets(mySheet).Cells(52, j + 2).Value其中52行对应 Use TP on Calcite sheet?;
        simContext.useEOSmix[kk] = data->Option_EoS;            // useEOSmix(kk) = Worksheets(mySheet).Cells(55, j + 2).Value其中55行对应Use Flash Calculator;
    }
    else
    {
        // simContext.UseTPCalciteSheet = Worksheets(MySheetMix).Cells(52, 8).Value 这里找不到第8列对应着什么
        // useEOSmix(kk) = Worksheets(MySheetMix).Cells(55, 8).Value 同上
    }

    // If simContext.UseTPCalciteSheet = "" Then simContext.UseTPCalciteSheet = 0
    // If UseTPpHMix(kk) = "" Then UseTPpHMix(kk) = 0
    // If UseTPVolMix(kk) = "" Then UseTPVolMix(kk) = 0
    // If useEOSmix(kk) = "" Then useEOSmix(kk) = 0 这里的逻辑判断似乎有问题，因为这些几乎全为double数组，而在判断中判断是否为空字符串;

    simContext.TofpH[kk] = data->T_pH; // TofpH(kk) = Worksheets(mySheet).Cells(58, j + 2).Value其中58行对应Temp. for pH meas.;
    simContext.PofpH[kk] = data->P_pH; // PofpH(kk) = Worksheets(mySheet).Cells(59, j + 2).Value其中59行对应Pres. for pH meas.;
    simContext.TofVol[kk] = data->T_Q; // TofVol(kk) = Worksheets(mySheet).Cells(60, j + 2).Value其中60行对应 T for fluids meas.;
    simContext.PofVol[kk] = data->P_Q; // PofVol(kk) = Worksheets(mySheet).Cells(61, j + 2).Value其中61行对P for for fluids meas.;

    if (simContext.RunNORM == 1)
    {
        // simContext.RaMix(kk) = Worksheets(mySheet).Cells(62, j + 2).Value / 1000000000000#表中62行没有字段
    }
    else
    {
        simContext.RaMix[kk] = 0;
    }
    simContext.useTPVol = simContext.UseTPVolMix[kk];
    simContext.SumofZMix[kk] = 0;
    // For iNG = 1 To 14
    // zMix(kk, iNG) = Worksheets(mySheet).Cells(65 + iNG, j + 2) / 100:
    // SumofZMix(kk) = SumofZMix(kk) + zMix(kk, iNG) 这里是对原表中66-79行进行读取，对应为 sample_oil_phase_information数据表，这里不采用循环
    simContext.zMix[kk][0] = data->C1_o / 100;
    simContext.zMix[kk][1] = data->CO2_o / 100;
    simContext.zMix[kk][2] = data->H2S_o / 100;
    simContext.zMix[kk][3] = data->C2_o / 100;
    simContext.zMix[kk][4] = data->C3_o / 100;
    simContext.zMix[kk][5] = data->iC4_o / 100;
    simContext.zMix[kk][6] = data->nC4_o / 100;
    simContext.zMix[kk][7] = data->iC5_o / 100;
    simContext.zMix[kk][8] = data->nC5_o / 100;
    simContext.zMix[kk][9] = data->C6_o / 100;
    simContext.zMix[kk][10] = data->C7_C12_o / 100;
    simContext.zMix[kk][11] = data->C13_C25_o / 100;
    simContext.zMix[kk][12] = data->C26_C80_o / 100;
    simContext.zMix[kk][13] = data->N2_o / 100;

    for (int i = 0; i < 14; i++)
    {
        simContext.SumofZMix[kk] = simContext.SumofZMix[kk] + simContext.zMix[kk][i];
    }

    if (simContext.SumofZMix[kk] > 0)
    {
        for (int i = 0; i < 14; i++)
        {
            simContext.zMix[kk][i] = simContext.zMix[kk][i] / simContext.SumofZMix[kk];
        }
    }
    simContext.zMix[kk][14] = 0.0;
    if (simContext.RunH2SGUI != 1)
    {
        // MultiplePpt = Worksheets("Input").Range("S11").Value
    }

    if (simContext.RunWhatIf == 1)
    {
        if (simContext.CaseCountWI[simContext.Loop1WI])
        {
            simContext.NaMix[kk] = simContext.NaWI[simContext.Loop1WI][simContext.Loop2WI];
        }
        if (simContext.CaseCountWI[simContext.Loop1WI] == 4)
        {
            simContext.CaMix[kk] = simContext.CaWI[simContext.Loop1WI][simContext.Loop2WI];
        }
        if (simContext.CaseCountWI[simContext.Loop1WI] == 5)
        {
            simContext.SrMix[kk] = simContext.SrWI[simContext.Loop1WI][simContext.Loop2WI];
        }

        if (simContext.CaseCountWI[simContext.Loop1WI] == 6)
        {
            simContext.BaMix[kk] = simContext.BaWI[simContext.Loop1WI][simContext.Loop2WI];
        }
        if (simContext.CaseCountWI[simContext.Loop1WI] == 10)
        {
            simContext.ClMix[kk] = simContext.ClWI[simContext.Loop1WI][simContext.Loop2WI];
        }
        if (simContext.CaseCountWI[simContext.Loop1WI] == 11)
        {
            simContext.SO4Mix[kk] = simContext.SO4WI[simContext.Loop1WI][simContext.Loop2WI];
        }
        if (simContext.CaseCountWI[simContext.Loop1WI] == 15)
        {
            simContext.HCO3AlkMix[kk] = simContext.HCO3AlkWI[simContext.Loop1WI][simContext.Loop2WI];
        }

        if (simContext.CaseCountWI[simContext.Loop1WI] == 16)
        {
            simContext.CO3AlkMix[kk] = simContext.CO3AlkWI[simContext.Loop1WI][simContext.Loop2WI];
        }
        if (simContext.CaseCountWI[simContext.Loop1WI] == 17)
        {
            simContext.TAcMix[kk] = simContext.TACWI[simContext.Loop1WI][simContext.Loop2WI];
        }
        if (simContext.CaseCountWI[simContext.Loop1WI] == 19)
        {
            simContext.usepHmix[kk] = 0;
            simContext.yCO2Mix[kk] = simContext.YCO2WI[simContext.Loop1WI][simContext.Loop2WI] / 100;
        }

        if (simContext.CaseCountWI[simContext.Loop1WI] == 20)
        {
            simContext.UseH2SgasMix[kk] = 1;
            simContext.yH2SMix[kk] = simContext.YH2SWI[simContext.Loop1WI][simContext.Loop2WI] / 100;
        }

        if (simContext.CaseCountWI[simContext.Loop1WI] == 21)
        {
            simContext.UseH2SgasMix[kk] = 1;
            simContext.TH2SaqMix[kk] = simContext.TH2SaqWI[simContext.Loop1WI][simContext.Loop2WI] / 100;
        }

        if (simContext.CaseCountWI[simContext.Loop1WI] == 22)
        {
            simContext.usepHmix[kk] = 1;
            simContext.pHMeterStpMix[kk] = simContext.pHMeterSTPWI[simContext.Loop1WI][simContext.Loop2WI] / 100;
        }

        if (simContext.CaseCountWI[simContext.Loop1WI] == 32)
        {
            simContext.ConcFactor[kk] = simContext.ConcFactorWI[simContext.Loop1WI][simContext.Loop2WI] / 100;
        }
    }

    simContext.TDSMix[kk] = 0.0000000000;
    if (simContext.TDSMix[kk] == 0) // If simContext.TDSMix(kk) = 0 Or simContext.TDSMix(kk) = "" Then 这里不对空字符串做判断
    {
        //修改   by   -彭非
        if (simContext.UseMolal == 0) {
            //for (int iTDS = 0; iTDS < 19; iTDS++);
                /*simContext.TDSMix(kk) = simContext.TDSMix(kk) + Worksheets(mySheet).Cells(9 + iTDS, j + 2).Value*/
            simContext.TDSMix[kk] += data->Na_aq;
            simContext.TDSMix[kk] += data->K_aq;
            simContext.TDSMix[kk] += data->Mg_aq;
            simContext.TDSMix[kk] += data->Ca_aq;
            simContext.TDSMix[kk] += data->Sr_aq;
            simContext.TDSMix[kk] += data->Ba_aq;
            simContext.TDSMix[kk] += data->FeII_aq;
            simContext.TDSMix[kk] += data->Zn_aq;
            simContext.TDSMix[kk] += data->Pb_aq;
            simContext.TDSMix[kk] += data->Cl_aq;
            simContext.TDSMix[kk] += data->SO4_aq;
            simContext.TDSMix[kk] += data->F_aq;
            simContext.TDSMix[kk] += data->Br_aq;
            simContext.TDSMix[kk] += data->Si_aq;
            simContext.TDSMix[kk] += data->Alk_Bicarbonate_aq;
            simContext.TDSMix[kk] += data->Alk_Carbonate_aq;
            simContext.TDSMix[kk] += data->OrgAcid_Acetate_aq;
            simContext.TDSMix[kk] += data->Ammonia_aq;
            simContext.TDSMix[kk] += data->B_aq;
        }
        else {
            //Use corelation of NaCl molal to TDS, initial guess for molality input
            //for (int iTDS = 0; iTDS < 19; iTDS++);
                //simContext.TDSMix(kk) = simContext.TDSMix(kk) + Worksheets(mySheet).Cells(9 + iTDS, j + 2).Value * 53459 / 2 
            simContext.TDSMix[kk] += data->Na_aq * 53459 / 2;
            simContext.TDSMix[kk] += data->K_aq * 53459 / 2;
            simContext.TDSMix[kk] += data->Mg_aq * 53459 / 2;
            simContext.TDSMix[kk] += data->Ca_aq * 53459 / 2;
            simContext.TDSMix[kk] += data->Sr_aq * 53459 / 2;
            simContext.TDSMix[kk] += data->Ba_aq * 53459 / 2;
            simContext.TDSMix[kk] += data->FeII_aq * 53459 / 2;
            simContext.TDSMix[kk] += data->Zn_aq * 53459 / 2;
            simContext.TDSMix[kk] += data->Pb_aq * 53459 / 2;
            simContext.TDSMix[kk] += data->Cl_aq * 53459 / 2;
            simContext.TDSMix[kk] += data->SO4_aq * 53459 / 2;
            simContext.TDSMix[kk] += data->F_aq * 53459 / 2;
            simContext.TDSMix[kk] += data->Br_aq * 53459 / 2;
            simContext.TDSMix[kk] += data->Si_aq * 53459 / 2;
            simContext.TDSMix[kk] += data->Alk_Bicarbonate_aq * 53459 / 2;
            simContext.TDSMix[kk] += data->Alk_Carbonate_aq * 53459 / 2;
            simContext.TDSMix[kk] += data->OrgAcid_Acetate_aq * 53459 / 2;
            simContext.TDSMix[kk] += data->Ammonia_aq * 53459 / 2;
            simContext.TDSMix[kk] += data->B_aq * 53459 / 2;
        }
    }

    //


    //If rho_Mix(kk) = 0 Or rho_Mix(kk) = "" Then rho_Mix(kk) = 0.9991 + 0.0000006398 * simContext.TDSMix(kk)
    simContext.rho_Mix[kk] = 0.00000000;
    if (simContext.rho_Mix[kk] == 0) simContext.rho_Mix[kk] = 0.9991 + 0.0000006398 * simContext.TDSMix[kk];

    //If the cell is leave empty assume a concentration factor=1
    if (simContext.ConcFactor[kk] == 0 || simContext.ConcFactor[kk] < 0) simContext.ConcFactor[kk] = 1;

    if (simContext.ConcFactor[kk] != 1) {
        simContext.Run_CalcConcFactor = 1;
        simContext.VwMix[kk] /= simContext.ConcFactor[kk];
    }

    if (simContext.UseMolal == 0)
    {
        simContext.NaMix[kk] = simContext.NaMix[kk] / (22990.0 * (0.9991 - 0.0000003612 * simContext.TDSMix[kk])) * simContext.ConcFactor[kk];
        simContext.KMix[kk] = simContext.KMix[kk] / (39098.0 * (0.9991 - 0.0000003612 * simContext.TDSMix[kk])) * simContext.ConcFactor[kk];
        simContext.MgMix[kk] = simContext.MgMix[kk] / (24305.0 * (0.9991 - 0.0000003612 * simContext.TDSMix[kk])) * simContext.ConcFactor[kk];
        simContext.CaMix[kk] = simContext.CaMix[kk] / (40080.0 * (0.9991 - 0.0000003612 * simContext.TDSMix[kk])) * simContext.ConcFactor[kk];
        simContext.SrMix[kk] = simContext.SrMix[kk] / (87620.0 * (0.9991 - 0.0000003612 * simContext.TDSMix[kk])) * simContext.ConcFactor[kk];
        simContext.BaMix[kk] = simContext.BaMix[kk] / (137330.0 * (0.9991 - 0.0000003612 * simContext.TDSMix[kk])) * simContext.ConcFactor[kk];
        simContext.FeMix[kk] = simContext.FeMix[kk] / (55847.0 * (0.9991 - 0.0000003612 * simContext.TDSMix[kk])) * simContext.ConcFactor[kk];
        simContext.ZnMix[kk] = simContext.ZnMix[kk] / (65380.0 * (0.9991 - 0.0000003612 * simContext.TDSMix[kk])) * simContext.ConcFactor[kk];
        simContext.PbMix[kk] = simContext.PbMix[kk] / (207200.0 * (0.9991 - 0.0000003612 * simContext.TDSMix[kk])) * simContext.ConcFactor[kk];
        simContext.ClMix[kk] = simContext.ClMix[kk] / (35450.0 * (0.9991 - 0.0000003612 * simContext.TDSMix[kk])) * simContext.ConcFactor[kk];
        simContext.SO4Mix[kk] = simContext.SO4Mix[kk] / (96064.0 * (0.9991 - 0.0000003612 * simContext.TDSMix[kk])) * simContext.ConcFactor[kk];
        simContext.FMix[kk] = simContext.FMix[kk] / (18998.0 * (0.9991 - 0.0000003612 * simContext.TDSMix[kk])) * simContext.ConcFactor[kk];
        simContext.BrMix[kk] = simContext.BrMix[kk] / (79904.0 * (0.9991 - 0.0000003612 * simContext.TDSMix[kk])) * simContext.ConcFactor[kk];
        simContext.TH4SiO4Mix[kk] = simContext.TH4SiO4Mix[kk] / (28085.0 * (0.9991 - 0.0000003612 * simContext.TDSMix[kk])) * simContext.ConcFactor[kk];
        simContext.HCO3AlkMix[kk] = simContext.HCO3AlkMix[kk] / (61019.0 * (0.9991 - 0.0000003612 * simContext.TDSMix[kk])) * simContext.ConcFactor[kk];
        simContext.CO3AlkMix[kk] = simContext.CO3AlkMix[kk] / (60019.0 * (0.9991 - 0.0000003612 * simContext.TDSMix[kk])) * simContext.ConcFactor[kk];
        simContext.TAcMix[kk] = simContext.TAcMix[kk] / (59046.0 * (0.9991 - 0.0000003612 * simContext.TDSMix[kk])) * simContext.ConcFactor[kk];
        simContext.TNH4Mix[kk] = simContext.TNH4Mix[kk] / (17031.0 * (0.9991 - 0.0000003612 * simContext.TDSMix[kk])) * simContext.ConcFactor[kk];
        simContext.TH3BO3Mix[kk] = simContext.TH3BO3Mix[kk] / (10811.0 * (0.9991 - 0.0000003612 * simContext.TDSMix[kk])) * simContext.ConcFactor[kk];

        if (simContext.UseH2SgasMix[kk] == 0)
        {
            simContext.TH2SaqMix[kk] = simContext.TH2SaqMix[kk] / (34080.0 * (0.9991 - 0.0000003612 * simContext.TDSMix[kk])) * simContext.ConcFactor[kk];
        }

        simContext.HAlkMix[kk] = simContext.HAlkMix[kk] / (0.9991 - 0.0000003612 * simContext.TDSMix[kk]);
        simContext.OHAlkMix[kk] = simContext.OHAlkMix[kk] / (0.9991 - 0.0000003612 * simContext.TDSMix[kk]);
        simContext.RaMix[kk] = simContext.RaMix[kk] / (226.0 * (0.9991 - 0.0000003612 * simContext.TDSMix[kk])) * simContext.ConcFactor[kk];
    }
    else if (simContext.UseMolal == 1)
    {
        simContext.NaMix[kk] = simContext.NaMix[kk] * simContext.ConcFactor[kk]; // Convert mg/L to molality
        simContext.KMix[kk] = simContext.KMix[kk] * simContext.ConcFactor[kk];
        simContext.MgMix[kk] = simContext.MgMix[kk] * simContext.ConcFactor[kk];
        simContext.CaMix[kk] = simContext.CaMix[kk] * simContext.ConcFactor[kk];
        simContext.SrMix[kk] = simContext.SrMix[kk] * simContext.ConcFactor[kk];
        simContext.BaMix[kk] = simContext.BaMix[kk] * simContext.ConcFactor[kk];
        simContext.FeMix[kk] = simContext.FeMix[kk] * simContext.ConcFactor[kk];
        simContext.ZnMix[kk] = simContext.ZnMix[kk] * simContext.ConcFactor[kk];
        simContext.PbMix[kk] = simContext.PbMix[kk] * simContext.ConcFactor[kk]; // Pb added
        simContext.ClMix[kk] = simContext.ClMix[kk] * simContext.ConcFactor[kk];
        simContext.SO4Mix[kk] = simContext.SO4Mix[kk] * simContext.ConcFactor[kk];
        simContext.FMix[kk] = simContext.FMix[kk] * simContext.ConcFactor[kk];
        simContext.BrMix[kk] = simContext.BrMix[kk] * simContext.ConcFactor[kk];           // Br added
        simContext.TH4SiO4Mix[kk] = simContext.TH4SiO4Mix[kk] * simContext.ConcFactor[kk]; // Input silica as SiO2
        simContext.HCO3AlkMix[kk] = simContext.HCO3AlkMix[kk] * simContext.ConcFactor[kk];
        simContext.CO3AlkMix[kk] = simContext.CO3AlkMix[kk] * simContext.ConcFactor[kk];
        simContext.TAcMix[kk] = simContext.TAcMix[kk] * simContext.ConcFactor[kk];
        simContext.TNH4Mix[kk] = simContext.TNH4Mix[kk] * simContext.ConcFactor[kk];
        simContext.TH3BO3Mix[kk] = simContext.TH3BO3Mix[kk] * simContext.ConcFactor[kk];

        if (simContext.UseH2SgasMix[kk] == 0)
        {
            simContext.TH2SaqMix[kk] = simContext.TH2SaqMix[kk] * simContext.ConcFactor[kk]; // Used to calculate yH2Sstp
        }

        // HAlkMix 和 OHAlkMix 被注释掉，不处理
        simContext.RaMix[kk] = simContext.RaMix[kk] * simContext.ConcFactor[kk];
    }
    simContext.AlkMix[kk] = simContext.HCO3AlkMix[kk] + 2 * simContext.CO3AlkMix[kk] - simContext.HAlkMix[kk] + simContext.OHAlkMix[kk];
    simContext.TCO2Mix[kk] = simContext.HCO3AlkMix[kk] + simContext.CO3AlkMix[kk];
    if (simContext.UseTPCalciteSheet != 1)
    {
        simContext.UseTPCalciteSheet = 0;
    }

    if (simContext.Run_Seawater_Mixing == 1 && simContext.j == 2)
    {
        simContext.yCO2Mix[kk] = pow(10.0, -3.5);
        simContext.yH2SMix[kk] = 0.0;
        simContext.TH2SaqMix[kk] = 0.0;
        simContext.UseH2SgasMix[kk] = 0; // 假设是 int 类型
    }
    if (simContext.UseTPpHMix[kk] == 1)
    {
        if (simContext.TofpH[kk] == 0.0 && simContext.UseSI == 0)
            simContext.TofpH[kk] = 77.0;
        if (simContext.PofpH[kk] == 0.0 && simContext.UseSI == 0)
            simContext.PofpH[kk] = 14.696;
        if (simContext.TofpH[kk] == 0.0 && simContext.UseSI == 1)
            simContext.TofpH[kk] = 25.0;
        if (simContext.PofpH[kk] == 0.0 && simContext.UseSI == 1)
            simContext.PofpH[kk] = 1.0;
    }
    else
    {
        if (simContext.UseSI == 0)
        {
            simContext.TofpH[kk] = 77.0;
            simContext.PofpH[kk] = 14.696; // set to 77 F and 14.696 psia as default for pH calculation
        }
        else
        {
            simContext.TofpH[kk] = 25.0;
            simContext.PofpH[kk] = 1.0;
        }
    }

    if (simContext.UseTPVolMix[kk] == 1)
    {
        if (simContext.TofVol[kk] == 0.0 && simContext.UseSI == 0)
            simContext.TofVol[kk] = 77.0; // set default to 77 F
        if (simContext.PofVol[kk] == 0.0 && simContext.UseSI == 0)
            simContext.PofVol[kk] = 14.696; // set default to 1 atm
        if (simContext.TofVol[kk] == 0.0 && simContext.UseSI == 1)
            simContext.TofVol[kk] = 25.0;
        if (simContext.PofVol[kk] == 0.0 && simContext.UseSI == 1)
            simContext.PofVol[kk] = 1.013254;
    }
    else
    {
        if (simContext.UseSI == 0)
        {
            simContext.TofVol[kk] = 77.0;
            simContext.PofVol[kk] = 14.696;
        }
        else
        { // set to 25 C and 1.013254 atm as default for EOS calculation
            simContext.TofVol[kk] = 25.0;
            simContext.PofVol[kk] = 1.013254;
        }
    }

    simContext.yCH4Mix[kk] = 1 - (simContext.yCO2Mix[kk] + simContext.yH2SMix[kk]);
    if (simContext.yCH4Mix[kk] < 0)
    {
        simContext.yCH4Mix[kk] = 0;
    }

    if (simContext.UseH2SgasMix[kk] != 0 && simContext.UseH2SgasMix[kk] != 1)
    {
        // MsgBox("Row 32 is expecting a value of 1 or 0, please enter a vlue of 1 or 0 and rerun")
        //     End
    }

    if (simContext.yCO2Mix[kk] > 1.0)
    {
        simContext.errmsg[0] = 1;  // 源代码：errmsg(1) = 1
        simContext.yCO2Mix[kk] = 1.0;
    }

    if (simContext.yCO2Mix[kk] < 0.0)
    {
        simContext.errmsg[1] = 2;
        simContext.yCO2Mix[kk] = 0.0;
        simContext.yCO2 = 0.0;
        simContext.CO2aq = 0.0;
        simContext.HCO3 = 0.0;
        simContext.CO3 = 0.0;
    }

    if (simContext.yH2SMix[kk] > 1.0)
    {
        simContext.errmsg[2] = 3;
        simContext.yH2SMix[kk] = 1.0;
        simContext.yH2S = 1.0;
        simContext.TH2SaqMix[kk] = 0.0; // This will cause the program to use yH2Sstp as the calculation for TH2Saq instead of the input sheet value
    }

    if (simContext.yH2SMix[kk] < 0.0)
    {
        simContext.errmsg[3] = 4;   //
        simContext.yH2SMix[kk] = 0.0;
        simContext.yH2S = 0.0;
        simContext.HS = 0.0;
        simContext.H2Saq = 0.0;
        simContext.TH2SaqMix[kk] = 0.0;
    }
}


double fTPFunc(int iTP);   //PartB要使用，在这声明了
double fH2ODensity(double TK, double PBar);

void ReadInputPartB(int kk, SampleData* data)
{
    simContext.VgTPMix[kk] = data->Q_Gas;
    simContext.VoMix[kk] = data->Q_Oil;
    simContext.VwMix[kk] = data->Q_Water;
    simContext.oilAPIgravMix[kk] = data->API;
    simContext.gasSpGravMix[kk] = data->SG_g;
    simContext.VMeOHMix[kk] = data->Q_MeOH;
    simContext.VMEGMix[kk] = data->Q_MEG;

    if (simContext.RunWhatIf == 1)
    {
        if (simContext.CaseCountWI[simContext.Loop1WI] == 25)
        {
            simContext.VgTPMix[kk] = VgTPWI[simContext.Loop1WI][simContext.Loop2WI];
        }
        if (simContext.CaseCountWI[simContext.Loop1WI] == 26)
        {
            simContext.VoMix[kk] = VoWI[simContext.Loop1WI][simContext.Loop2WI];
        }
        if (simContext.CaseCountWI[simContext.Loop1WI] == 27)
        {
            simContext.VwMix[kk] = VwWI[simContext.Loop1WI][simContext.Loop2WI];
        }
    }

    if (simContext.Run_Seawater_Mixing == 1 && simContext.LoopMixing == 1 && kk == 0)
    {
        //offset:VwSW1 = VwMix(1): VgSW1 = VgTPMix(1): VoSW1 = VoMix(1): VMeOHSW1 = VMeOHMix(1): VMEGSW1 = VMEGMix(1)
        VwSW1 = simContext.VwMix[0];
        VgSW1 = simContext.VgTPMix[0];
        VoSW1 = simContext.VoMix[0];
        VMeOHSW1 = simContext.VMeOHMix[0];
        VMEGSW1 = simContext.VMEGMix[0];
    }

    // 舍弃此字段
    if (simContext.Run10TestCases == 1)
    {
        simContext.usepHmix[kk] = 0;
        // yCO2Mix(kk) = Worksheets("Calcite").Cells(15 + simContext.Loop10, 3) / 100
    }

    if (simContext.RunMultiMix == 1)
    {
        // VgTPMix(kk) = Worksheets("MultiMix").Cells(2 + simContext.LoopResChem, 14 + j).Value
        // VoMix(kk) = Worksheets("MultiMix").Cells(2 + simContext.LoopResChem, 9 + j).Value
        // VwMix(kk) = Worksheets("MultiMix").Cells(2 + simContext.LoopResChem, 4 + j).Value
    }
    if (RunMultiMixSlb == 1)
    {
        // VgTPMix[kk] = Worksheets("MultiMix_Slb").Cells(2 + simContext.LoopResChem, 14 + j).Value;
        // VoMix[kk] = Worksheets("MultiMix_Slb").Cells(2 + simContext.LoopResChem, 9 + j).Value;
        // VwMix[kk] = Worksheets("MultiMix_Slb").Cells(2 + simContext.LoopResChem, 4 + j).Value
    }
    if (simContext.Run_MixingTwoWells == 1)
    {
        simContext.VgTPMix[kk] = simContext.VgTPMix[kk] * simContext.MixFrac[kk];
        simContext.VoMix[kk] = simContext.VoMix[kk] * simContext.MixFrac[kk];
        simContext.VwMix[kk] = simContext.VwMix[kk] * simContext.MixFrac[kk];
        simContext.VMeOHMix[kk] = simContext.VMeOHMix[kk] * simContext.MixFrac[kk];
        simContext.VMEGMix[kk] = simContext.VMEGMix[kk] * simContext.MixFrac[kk];
    }
    if (simContext.Run_Seawater_Mixing == 1)
    {
        simContext.VwMix[kk] = VwSW1 * simContext.MixFrac[kk];

        if (kk == 0)
        {
            simContext.VoMix[kk] = VoSW1 * simContext.MixFrac[kk];
            simContext.VgTPMix[kk] = VgSW1 * simContext.MixFrac[kk];
            simContext.VMeOHMix[kk] = VMeOHSW1 * simContext.MixFrac[kk];
            simContext.VMEGMix[kk] = VMEGSW1 * simContext.MixFrac[kk];
        }
        else
        {
            simContext.VoMix[kk] = 0;
            simContext.VgTPMix[kk] = 0;
            simContext.VMeOHMix[kk] = 0;
            simContext.VMEGMix[kk] = 0;
        }
    }

    if (simContext.UseSI == 0)
    {
        // convert 1000 ft^3 to m^3
        simContext.VgTPMix[kk] = simContext.VgTPMix[kk] * 28.31685;
    }
    else
    {
        // convert 1000 m^3 to m^3
        simContext.VgTPMix[kk] = simContext.VgTPMix[kk] * 1000;
    }

    if (simContext.UseSI == 1)
    {
        // Convert m^3 to barrels
        simContext.VoMix[kk] = simContext.VoMix[kk] / 0.159;
        simContext.VwMix[kk] = simContext.VwMix[kk] / 0.159;
        simContext.VMeOHMix[kk] = simContext.VMeOHMix[kk] / 0.159;
        simContext.VMEGMix[kk] = simContext.VMEGMix[kk] / 0.159;

        // Convert temperature from Celsius to Fahrenheit
        simContext.TofVol[kk] = simContext.TofVol[kk] * 9.0 / 5.0 + 32.0;

        // Convert pressure from Bar to psia
        simContext.PofVol[kk] = simContext.PofVol[kk] * 14.503774;

        // Same conversions for pH-related data
        simContext.TofpH[kk] = simContext.TofpH[kk] * 9.0 / 5.0 + 32.0;
        simContext.PofpH[kk] = simContext.PofpH[kk] * 14.503774;
    }

    TpH = simContext.TofpH[kk];
    PpH = simContext.PofpH[kk];
    TVol = simContext.TofVol[kk];
    Pvol = simContext.PofVol[kk];

    // 若为 0 或空字符串则赋默认值
    // 在 C 中不能判断 ""，只能判断是否为 0
    if (simContext.oilAPIgravMix[kk] == 0)
    {
        simContext.oilAPIgravMix[kk] = 30.0;
    }
    if (simContext.gasSpGravMix[kk] == 0)
    {
        simContext.gasSpGravMix[kk] = 0.6;
    }

    if (simContext.UseTPVolMix[kk] == 0)
    {
        fTPFunc(0); // STP condition  mt = fTPFunc(0)

        simContext.GasDensityMix[kk] = simContext.gasSpGravMix[kk] *
            (Patm * 28.97 / (0.08206 * TK)); // kg/m^3

        // OilDensityMix(kk) = (141.5 / (oilAPIgravMix(kk) + 131.5)) * (fH2ODensity(TK, PBar) / 1000) 'conver density from API gravity
        simContext.OilDensityMix[kk] = (141.5 / (simContext.oilAPIgravMix[kk] + 131.5)) *
            (fH2ODensity(TK, PBar) / 1000.0);
    }
    else
    { // UseTPVolMix == 1

        fTPFunc(1); // Use actual T, P

        simContext.GasDensityMix[kk] = simContext.gasSpGravMix[kk] *
            (Patm * 28.97 / (0.08206 * TK)); // kg/m^3
        // 调用外部函数，暂时忽略fH2ODensity
        simContext.OilDensityMix[kk] = (141.5 / (simContext.oilAPIgravMix[kk] + 131.5)) *
            (fH2ODensity(TK, PBar) / 1000.0);
    }

    // If volume of gas is equal to zero, then add one mL of gas
    if (simContext.VgTPMix[kk] == 0.0)
    {
        simContext.VgTPMix[kk] = 1.0 / 1000000.0;
    }

    // If no oil is reported, add one ml of oil per day, to avoid singularities
    if (simContext.VoMix[kk] == 0.0)
    {
        simContext.VoMix[kk] = 1.0 / 159.0 / 1000.0;
    }

    // If no water is reported, add one ml of water per day, to avoid singularities
    if (simContext.VwMix[kk] == 0.0)
    {
        simContext.VwMix[kk] = 1.0 / 159.0 / 1000.0;
    }

    // MeOH and MEG cannot both be > 0
    if (simContext.VMeOHMix[kk] > 0.0 && simContext.VMEGMix[kk] > 0.0)
    {
        // printf("MeOH and MEG cannot both be > 0. Please check the Input Sheet and try it again.\n");
        // exit(1); // equivalent to VB's 'End'
        return;
    }

    // Mass of methanol (Kg)
    simContext.mass_MeOH_mix[kk] = 159.0 * simContext.VMeOHMix[kk] * 0.7914;

    // Mass of MEG (Kg)
    simContext.mass_MEG_mix[kk] = 159.0 * simContext.VMEGMix[kk] * 1.1098;

    // Brine flow rate, unit: cm3/sec
    QBrineFlow = simContext.VwMix[kk] * 159.0 * 1000.0 / 86400.0;

    // Calculate Surface Area contacted by brine, cm²
    SArea = pi * PipeID * PipeL;

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

void CalcIonicStrength()
{
    // Calcule common system terms and functions for Pitzer theory: 
    // ISt, Z, m-total gX, gpX, JX, and JpX.

    mtotal = 0;
    MoleCharge = 0.0;
    ISt = 0.0;
    SumOfCations = 0.0;

    for (int c = 0; c < NumCat; c++) {
        MoleCharge = MoleCharge + ChCat[c] * mc[c];
        SumOfCations = SumOfCations + ChCat[c] * mc[c];
        ISt = ISt + pow(ChCat[c], 2.0) * mc[c];
        mtotal = mtotal + mc[c];
    }

    // SumOfAnions = -0.0000001
    SumOfAnions = 0.0;

    for (int a = 0; a < NumAn; a++) {
        MoleCharge = MoleCharge + fabs(ChAn[a]) * ma[a];
        SumOfAnions = SumOfAnions + ChAn[a] * ma[a];
        ISt = ISt + pow(ChAn[a], 2.0) * ma[a];
        mtotal = mtotal + ma[a];
    }

    if (ISt <= 0.0) {
        // When only negative alkalinity is entered in the input sheet,
        // set ionic molality = (H+OH) of pure water
        ISt = 2.0 * 0.0000001;
    }

    ISt = ISt / 2.0;
    DpHj = 0.129 * pow(ISt, 0.5);
}

double fRatioOilBPoints(double API) {
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

double PsatH2O(double tk)
{
    const double TC = 647.096; /* 临界温度 K */
    const double Pc = 220.64;  /* 临界压力 Bar */
    const double a1 = -7.85951783;
    const double a2 = 1.84408259;
    const double a3 = -11.7866497;
    const double a4 = 22.6807411;
    const double a5 = -15.9618719;
    const double a6 = 1.80122502;

    double T = 1.0 - tk / TC;

    double Psat = Pc * exp((TC / tk) * (a1 * T + a2 * pow(T, 1.5) + a3 * pow(T, 3.0) +
        a4 * pow(T, 3.5) + a5 * pow(T, 4.0) + a6 * pow(T, 7.5)));
    return Psat;
}


/**
 * @brief 计算热力学平衡常数
 *
 * 此函数计算多种化学物质在水溶液中的热力学平衡常数，包括：
 * - 气体溶解度常数（CO2、CH4、H2S等）
 * - 酸解离常数（碳酸、硫化氢、乙酸、硅酸等）
 * - 矿物溶度积常数（方解石、石膏、重晶石、闪锌矿等）
 * - 配合物稳定常数（Zn、Pb的氯配合物和硫配合物等）
 *
 * 计算考虑了温度、压力和离子强度的影响，适用于地球化学模拟。
 *
 */
double DRaCalcite; //--注：C1特有的变量，在C1里被赋值但从未被任何地方调用
void C1_ThermodynamicEquilConsts()
{
    double dk;
    double KspRaSO4;//仅在C1使用的局部变量


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
    else
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
    double ln_H_T_H2S = (13788.0 / TK - 185.19 + 29.087 * log(TK)
        - 0.027637 * TK - 1445200.0 / pow(TK, 2));

    dV = (33.18 + 0.092661 * TC - 0.00054853 * pow(TC, 2)
        + 0.0000015354 * pow(TC, 3) - 0.0000000015459 * pow(TC, 4)) * 0.001;

    KgwH2S = pow(exp(ln_H_T_H2S) * exp(dV * (Patm * 1.01353 - Psat) / RBar / TK), -1)
        / 14.503774;

    // ================= Carbonate minerals =================
    if ((PBar - Psat) <= 500) {
        // -1
        KspCalcite = pow(10, -171.9065 - 0.077993 * TK + 2839.319 / TK
            + 71.595 * log10(TK))
            * pow(10, (0.514 - 0.000197 * TC + 1.7096E-15 * pow(TC, 6))
                * (PBar - Psat) / 500);

        // -2
        KspBarite = pow(10, 136.035 - 7680.41 / TK - 48.595 * log10(TK))
            * pow(10, (0.394 - 0.0001119 * TC + 1.5305E-15 * pow(TC, 6))
                * (PBar - Psat) / 500);
    }
    else {
        //-1
        KspCalcite = pow(10, -171.9065 - 0.077993 * TK + 2839.319 / TK
            + 71.595 * log10(TK))
            * pow(10, 0.514 - 0.000197 * TC + 1.7096E-15 * pow(TC, 6))
            * pow(10, ((0.928 - 0.000079 * TC + 2.1485E-15 * pow(TC, 6))
                - (0.514 - 0.000197 * TC + 1.7096E-15 * pow(TC, 6)))
                * ((PBar - Psat) - 500) / 500);
        // -2
        KspBarite = pow(10, 136.035 - 7680.41 / TK - 48.595 * log10(TK))
            * pow(10, 0.394 - 0.0001119 * TC + 1.5305E-15 * pow(TC, 6))
            * pow(10, ((0.674 + 0.0001229 * TC + 1.9202E-15 * pow(TC, 6))
                - (0.394 - 0.0001119 * TC + 1.5305E-15 * pow(TC, 6)))
                * ((PBar - Psat) - 500) / 500);
    }

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

    // ================= -1  Hemihydrate (CaSO4·0.5H2O) =================
    // ================= -2  Gypsum (CaSO4·2H2O) =================
    if ((PBar - Psat) <= 500) {
        // -1 
        KspHemihydrate = pow(10, -1 * (-183.603 + 8033.771 / TK + 28.173 * log(TK)))
            * pow(10, (0.423 - 0.0001 * TC + 1.6176E-15 * pow(TC, 6))
                * (PBar - Psat) / 500);

        // -2
        KspGypsum = pow(10, -128.9622 + 1.761129 * TK - 0.01484172 * pow(TK, 2)
            + 0.000837554 * pow(TK, 2.5) - 0.00001384613 * pow(TK, 3))
            * pow(10, (0.328 + 0.0000589 * TC + 1.612E-15 * pow(TC, 6))
                * (PBar - Psat) / 500);
    }

    else {
        // - 1
        KspHemihydrate = pow(10, -1 * (-183.603 + 8033.771 / TK + 28.173 * log(TK)))
            * pow(10, 0.423 - 0.0001 * TC + 1.6176E-15 * pow(TC, 6))
            * pow(10, ((0.729 + 0.0001576 * TC + 2.0302E-15 * pow(TC, 6))
                - (0.423 - 0.0001 * TC + 1.6176E-15 * pow(TC, 6)))
                * ((PBar - Psat) - 500) / 500);

        // -2
        KspGypsum = pow(10, -128.9622 + 1.761129 * TK - 0.01484172 * pow(TK, 2)
            + 0.000837554 * pow(TK, 2.5) - 0.00001384613 * pow(TK, 3))
            * pow(10, (0.328 + 0.0000589 * TC + 1.612E-15 * pow(TC, 6)))
            * pow(10, ((0.553 + 0.0004785 * TC + 2.018E-15 * pow(TC, 6))
                - (0.328 + 0.0000589 * TC + 1.612E-15 * pow(TC, 6)))
                * ((PBar - Psat) - 500) / 500);
    }


    // ================= Anhydrite =================
    KspAnhydrite = pow(10,
        (-44.9380794269878 + 1.08048295151073E-02 * PBar) * log10(TK)
        + (1.818855325 * PBar - 4943.71187) / TK
        - 3.19664522256897E-02 * PBar + 123.46568);

    // ================= Siderite (FeCO3) =================
    KspSiderite = exp(129.97 / TK - 50.205 + 7.3143 * log(TK) - 0.052913 * TK)
        * exp(-(-48.76 - 0.5304 * TC) * (Patm - Psat) / (R * TK));

    // ================= RaSO4 (same pressure dependence as Barite) =================
    if ((PBar - Psat) <= 500) {
        KspRaSO4 = pow(10, 137.98 - 8346 / TK - 48.595 * log10(TK))
            * pow(10, (0.394 - 0.0001119 * TC + 1.5305E-15 * pow(TC, 6))
                * (PBar - Psat) / 500);

    }
    else {
        KspRaSO4 = pow(10, 137.98 - 8346 / TK - 48.595 * log10(TK))
            * pow(10, 0.394 - 0.0001119 * TC + 1.5305E-15 * pow(TC, 6))
            * pow(10, ((0.674 + 0.0001229 * TC + 1.9202E-15 * pow(TC, 6))
                - (0.394 - 0.0001119 * TC + 1.5305E-15 * pow(TC, 6)))
                * ((PBar - Psat) - 500) / 500);
    }

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

    // ================= Dolomite (CaMg(simContext.CO3)2) =================
    KspDol = pow(10, -1 * (257.2181 + 0.1438 * TK - 3197.3234 / TK - 109.6227 * log10(TK)))
        * pow(10, (0.982 - 0.0003101 * TC) * Patm / 500);

    KspDol = pow(10, -1 * (-175.9 + 6835.4 / TK + 68.727 * log10(TK)))
        * pow(10, (0.982 - 0.0003101 * TC) * Patm / 500);

    // ================= Cosolvent dielectric constant =================
    double DielecConst = 1 / (0.0068 * pow(xMeOH, 2) + 0.0115 * xMeOH + 0.0128);

    double KstCaHCO3 = 0.034877
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
    if ((PBar - Psat) <= 500) {
        // -1 
        KspSrCO3 = pow(10, 155.6841 - 7272.6012 / TK - 56.8052 * log10(TK))
            * pow(10, (0.528 - 0.000259 * TC + 1.682E-15 * pow(TC, 6))
                * (PBar - Psat) / 500);

        // -2
        KspBaCO3 = pow(10, 244.3819 - 11526.0874 / TK - 86.6577 * log10(TK))
            * pow(10, (0.523 - 0.0003039 * TC + 1.631E-15 * pow(TC, 6))
                * (PBar - Psat) / 500);
    }
    else {
        // - 1
        KspSrCO3 = pow(10, 155.6841 - 7272.6012 / TK - 56.8052 * log10(TK))
            * pow(10, (0.528 - 0.000259 * TC + 1.682E-15 * pow(TC, 6))
                * pow(10, ((0.951 - 0.0001671 * TC + 2.114E-15 * pow(TC, 6))
                    - (0.528 - 0.000259 * TC + 1.682E-15 * pow(TC, 6)))
                    * ((PBar - Psat) - 500) / 500));

        // -2
        KspBaCO3 = pow(10, 244.3819 - 11526.0874 / TK - 86.6577 * log10(TK))
            * pow(10, (0.523 - 0.0003039 * TC + 1.631E-15 * pow(TC, 6)))
            * pow(10, ((0.936 - 0.0002343 * TC + 2.05E-15 * pow(TC, 6))
                - (0.523 - 0.0003039 * TC + 1.631E-15 * pow(TC, 6)))
                * ((PBar - Psat) - 500) / 500);
    }

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

    // ---------------------- Silicates (from SOLMINEQ.88) ----------------------

    q1 = 417.2882131; q2 = 1.474994811; q3 = -73.40388764; q4 = 0.110687481; q5 = 0.0000797508;
    q6 = 0.432; q7 = 0.0001186; q8 = 1.7613E-15; q9 = 0.661; q10 = 0.0007287; q11 = 2.1954E-15; // Ksp Chrysotile, Mg3Si2O5(OH)4

    if ((PBar - Psat) <= 500)
        KspChrysotile = pow(10, q1 + q2 / TK + q3 * log(TK) + q4 * TK + q5 / pow(TK, 2))
        * pow(10, (q6 + q7 * TC + q8 * pow(TC, 6)) * (PBar - Psat) / 500);
    else
        KspChrysotile = pow(10, q1 + q2 / TK + q3 * log(TK) + q4 * TK + q5 / pow(TK, 2))
        * pow(10, q6 + q7 * TC + q8 * pow(TC, 6))
        * pow(10, ((q9 + q10 * TC + q11 * pow(TC, 6))
            - (q6 + q7 * TC + q8 * pow(TC, 6)))
            * ((PBar - Psat) - 500) / 500);

    // -------------------------------------------------------------------------

    q1 = 210.5283437; q2 = 1.482041773; q3 = -35.80594147; q4 = 0.043892489; q5 = -0.0000593566;
    q6 = 0.422; q7 = -0.0000558; q8 = 2.113E-15; q9 = 0.638; q10 = 0.0001832; q11 = 1.37E-15; // Ksp Diopside, CaMgSi2O6

    if ((PBar - Psat) <= 500)
        KspDiopside = pow(10, q1 + q2 / TK + q3 * log(TK) + q4 * TK + q5 / pow(TK, 2))
        * pow(10, (q6 + q7 * TC + q8 * pow(TC, 6)) * (PBar - Psat) / 500);
    else
        KspDiopside = pow(10, q1 + q2 / TK + q3 * log(TK) + q4 * TK + q5 / pow(TK, 2))
        * pow(10, q6 + q7 * TC + q8 * pow(TC, 6))
        * pow(10, ((q9 + q10 * TC + q11 * pow(TC, 6))
            - (q6 + q7 * TC + q8 * pow(TC, 6)))
            * ((PBar - Psat) - 500) / 500);

    // -------------------------------------------------------------------------

    q1 = 201.8967952; q2 = 6.283898117; q3 = -32.54753804; q4 = 0.020272784; q5 = 0.044972535;
    q6 = 0.645; q7 = -0.000798; q8 = 1.4829E-15; q9 = 1.156; q10 = -0.0007986; q11 = 2.0008E-15; // Ksp Greenalite, Fe3Si2O5(OH)4

    if ((PBar - Psat) <= 500)
        KspGreenalite = pow(10, q1 + q2 / TK + q3 * log(TK) + q4 * TK + q5 / pow(TK, 2))
        * pow(10, (q6 + q7 * TC + q8 * pow(TC, 6)) * (PBar - Psat) / 500);
    else
        KspGreenalite = pow(10, q1 + q2 / TK + q3 * log(TK) + q4 * TK + q5 / pow(TK, 2))
        * pow(10, q6 + q7 * TC + q8 * pow(TC, 6))
        * pow(10, ((q9 + q10 * TC + q11 * pow(TC, 6))
            - (q6 + q7 * TC + q8 * pow(TC, 6)))
            * ((PBar - Psat) - 500) / 500);

    // -------------------------------------------------------------------------

    q1 = -38.98708976; q2 = 1.465048127; q3 = 6.604416649; q4 = -0.008608808; q5 = -0.0000366598;
    q6 = 0.111; q7 = -0.0000936; q8 = 4.53E-17; q9 = 0.208; q10 = -0.000153; q11 = 5.59E-17; // Ksp Quartz, SiO2

    if ((PBar - Psat) <= 500)
        KspQuartz = pow(10, q1 + q2 / TK + q3 * log(TK) + q4 * TK + q5 / pow(TK, 2))
        * pow(10, (q6 + q7 * TC + q8 * pow(TC, 6)) * (PBar - Psat) / 500);
    else
        KspQuartz = pow(10, q1 + q2 / TK + q3 * log(TK) + q4 * TK + q5 / pow(TK, 2))
        * pow(10, q6 + q7 * TC + q8 * pow(TC, 6))
        * pow(10, ((q9 + q10 * TC + q11 * pow(TC, 6))
            - (q6 + q7 * TC + q8 * pow(TC, 6)))
            * ((PBar - Psat) - 500) / 500);

    // -------------------------------------------------------------------------

    q1 = -22.01528297; q2 = 1.353509335; q3 = 3.633407549; q4 = -0.004631088; q5 = -0.0000366598;
    q6 = 0.118; q7 = -0.0000936; q8 = 4.53E-17; q9 = 0.216; q10 = -0.000153; q11 = 5.59E-17; // Ksp Amorphous Silica

    if ((PBar - Psat) <= 500)
        KspAmSilica = pow(10, q1 + q2 / TK + q3 * log(TK) + q4 * TK + q5 / pow(TK, 2))
        * pow(10, (q6 + q7 * TC + q8 * pow(TC, 6)) * (PBar - Psat) / 500);
    else
        KspAmSilica = pow(10, q1 + q2 / TK + q3 * log(TK) + q4 * TK + q5 / pow(TK, 2))
        * pow(10, q6 + q7 * TC + q8 * pow(TC, 6))
        * pow(10, ((q9 + q10 * TC + q11 * pow(TC, 6))
            - (q6 + q7 * TC + q8 * pow(TC, 6)))
            * ((PBar - Psat) - 500) / 500);

    // -------------------------------------------------------------------------

    Psat = PsatH2O(TK); // Reset vapor pressure of water

    // ---------------------- Ra partition coefficients ----------------------
    DRaBarite = pow(10, (428.2 / TK - 1.181));
    DRaCelestite = pow(10, (1438 / TK - 2.375));
    DRaAnhydrite = pow(10, (1930 / TK - 3.57));
    DRaCalcite = pow(10, (40.53 / TK - 0.2221));         //--注：DRaCalcite被赋值但从未被调用

    GammaSolidRaBarite = 266.6258 / 1.987 / TK;
    GammaSolidRaCelestite = 1950.176 / 1.987 / TK;
    GammaSolidRaAnhydrite = GammaSolidRaCelestite;

    // ---------------------- FeS aqueous complex ----------------------
    KstFeSaq = pow(10, -1 * (-4.57800335606811 + 1772.69547293922 / TK
        + 0.386641958182017 * log10(TK)));

    // ---------------------- Zn and Pb chloride complexes ----------------------
    BetaDot[iZnCl] = pow(10, 659.6174 + 0.3209708 * TK - 0.000115975 * pow(TK, 2)
        - 15493.94 / TK - 121.5629 * log(TK));
    BetaDot[iZnCl2] = pow(10, -4283.384 - 2.106234 * TK + 0.0009509436 * pow(TK, 2)
        + 97639.55 / TK + 789.8018 * log(TK));
    BetaDot[iZnCl3] = pow(10, 679.9657 + 0.3263276 * TK - 0.0001204947 * pow(TK, 2)
        - 15801.06 / TK - 125.1792 * log(TK));
    BetaDot[iZnCl4] = pow(10, -972.732 + 12.21668 * TK - 0.06111154 * pow(TK, 2)
        + 0.0001521141 * pow(TK, 3)
        - 0.0000001881333 * pow(TK, 4)
        + 9.252311E-11 * pow(TK, 5));

    BetaDot[iPbCl] = pow(10, -52.97153 - 0.02182295 * TK + 0.00002386252 * pow(TK, 2)
        + 1994.203 / TK + 9.141143 * log(TK));
    BetaDot[iPbCl2] = pow(10, 673.5621 + 0.2532472 * TK - 0.00006306824 * pow(TK, 2)
        - 17402.88 / TK - 119.8971 * log(TK));
    BetaDot[iPbCl3] = pow(10, -1819.472 - 0.8596028 * TK + 0.0003943687 * pow(TK, 2)
        + 43775.55 / TK + 332.6902 * log(TK));
    BetaDot[iPbCl4] = pow(10, -509.5157 + 6.837851 * TK - 0.03624178 * pow(TK, 2)
        + 0.00009494584 * pow(TK, 3)
        - 0.000000122756 * pow(TK, 4)
        + 6.262392E-11 * pow(TK, 5));

    // ---------------------- Zn & Pb simContext.HS complexes ----------------------
    q1 = 848.820711081324; q2 = 8523.81138931672; q3 = -200.318510598844; q4 = 1.12515623321474; q5 = -6.87407323192827e-04;
    BetaDot[iZnHS2] = pow(10.0, q1 + q2 / TK + q3 * log(TK) + q4 * TK + q5 * (TK * TK));

    q1 = 14.7708017586138; q2 = -19.937638657592; q3 = 0.877588609670279; q4 = -2.29953487957037e-02; q5 = 7.70611430037756e-06;
    BetaDot[iZnHS3] = pow(10.0, q1 + q2 / TK + q3 * log(TK) + q4 * TK + q5 * (TK * TK));

    q1 = 24.303146180657; q2 = -7.41314320085126; q3 = -1.11914072688766; q4 = -2.55711457549767e-02; q5 = 4.86459072972289e-05;
    BetaDot[iZnOH2] = pow(10.0, q1 + q2 / TK + q3 * log(TK) + q4 * TK + q5 * (TK * TK));

    q1 = -1622.2843; q2 = 42144.5235; q3 = 297.1852; q4 = -0.7658; q5 = 0.000345;
    BetaDot[iPbHS2] = pow(10, q1 + q2 / TK + q3 * log(TK) + q4 * TK + q5 * pow(TK, 2));

    q1 = -832.0208; q2 = 25256.693; q3 = 150.534; q4 = -0.3694; q5 = 0.000177;
    BetaDot[iPbHS3] = pow(10, q1 + q2 / TK + q3 * log(TK) + q4 * TK + q5 * pow(TK, 2));
}



double fV0(double tk, double pBar, double q1, double q2, double q3, double q4, double q5, double q6, double q7, double q8,
    double q9, double q10, double q11, double q12, double q13, double q14, double q15, double q16, double q17, double q18)
{
    double result = 0.0;

    result = q1
        + q2 / tk
        + q3 * tk
        + q4 * tk * tk
        + q5 / (tk - 227.0)
        + q6 / (647.0 - tk)
        + pBar * (q7 + q8 / tk + q9 * tk + q10 * tk * tk
            + q11 / (tk - 227.0) + q12 / (647.0 - tk))
        + pBar * pBar * (q13 + q14 / tk + q15 * tk + q16 * tk * tk
            + q17 / (tk - 227.0) + q18 / (647.0 - tk));

    return result;
}


void V0TP(double tk, double pBar)
{
    // Na
    V0_c[iNa] = fV0(tk, pBar,
        8.76686173829976, 10.7463747460684, -3.94438220704875E-02,
        6.24254747432051E-05, -270.67565216552, 3.71906197249154,
        -0.548144793641968, 60.4681423698375, 1.6487628506002E-03,
        -1.78543820686743E-06, -7.22279554933479E-02, 6.1016223460645,
        5.37401040542647E-04, -5.78337684998825E-02, -1.64057441413787E-06,
        1.85284769422731E-09, 2.33307684437836E-04, -8.30091414725064E-03);

    // K
    V0_c[iK] = fV0(tk, pBar,
        1365.58178, -146187.60179, -4.00314, 0.004292463,
        0.0, -18894.00317, -6.66675, 715.61054, 0.019742057,
        -0.000020810979, 0.0, 81.91098, 0.00534941, -0.573121,
        -0.0000158576885, 0.0000000166987, 0.0, -0.0649312);

    // Mg
    V0_c[iMg] = fV0(tk, pBar,
        15.5470319999999, 18848.33832, -0.380047233, 0.00100500148,
        -761.133887, -22952.60934, -1.27205782, 142.833027,
        0.00343937874, -0.00000368366162, 0.0, 36.3788742,
        0.0, 0.0, 0.000000400785822, 0.0, 0.0, -0.0429183805);

    // Ca
    V0_c[iCa] = fV0(tk, pBar,
        136.567817, 3135.53072, -0.68264817, 0.0012917585,
        -761.133887, -22952.60934, -1.6506531, 162.991634,
        0.0048786976, -0.00000616046222, 0.0, 71.877495,
        0.0, 0.0, 0.000000400785822, 0.0, 0.0, -0.0429183805);

    // Ba
    V0_c[iBa] = fV0(tk, pBar,
        131.77473, 3135.53072, -0.65074076, 0.0012917585,
        -761.133887, -22952.60934, -1.6506531, 162.991634,
        0.0048786976, -0.00000616046222, 0.0, 71.877495,
        0.0, 0.0, 0.000000400785822, 0.0, 0.0, -0.0429183805);

    // Sr
    V0_c[iSr] = fV0(tk, pBar,
        131.475189, 3135.53072, -0.66249868, 0.0012917585,
        -761.133887, -22952.60934, -1.6506531, 167.104986,
        0.0048786976, -0.00000616046222, 0.0, 69.810131,
        0.0, 0.0, 0.000000400785822, 0.0, 0.0, -0.0429183805);

    // Cl
    V0_a[iCl] = fV0(tk, pBar,
        195.93822, -23046.39821, -0.29604, 0.000299867,
        0.0, -13674.59683, -0.20212, 19.60946, 0.000482443,
        -0.000000766921, 0.0, 21.30102, 0.0, 0.0,
        0.0000000714885, 0.0, 0.0, -0.00727);

    // SO4 (取后面 Dai 20151027 的参数)
    double q1 = 20.1255488592824;
    double q2 = 8382.68656863382;
    double q3 = 8.63029649907132E-02;
    double q4 = 3.86933893465996E-05;
    double q5 = -262.073120048569;
    double q6 = -20014.7787488445;
    double q7 = 1.93246311860811E-02;
    double q8 = 1.16742679461389;
    double q9 = 9.95461079547389E-06;
    double q10 = 7.58647008800297E-08;
    double q11 = -1.30203576761739;
    double q12 = -2.31700638251545;
    double q13 = -3.61401771714286E-05;
    double q14 = -4.71960303681654E-03;
    double q15 = -7.8135942334434E-09;
    double q16 = 4.3129299430231E-11;
    double q17 = 1.77316241419105E-03;
    double q18 = 7.24281664341625E-03;

    V0_a[iSO4] =
        q1 + q2 / tk + q3 * tk + q4 * pow(tk, 2)
        + q5 / (tk - 227.0)
        + q6 / (647.0 - tk)
        + pBar * (q7 + q8 / tk + q9 * tk + q10 * pow(tk, 2)
            + q11 / (tk - 227.0)
            + q12 / (647.0 - tk))
        + pow(pBar, 2) * (q13 + q14 / tk + q15 * tk + q16 * pow(tk, 2)
            + q17 / (tk - 227.0)
            + q18 / (647.0 - tk));
}



double fH2ODensity(double tk, double pBar) {
    double nreg1[34] = { 0.14632971213167,-0.84548187169114,-3.756360367204,3.3855169168385,-0.95791963387872,
        0.15772038513228,-0.016616417199501,8.1214629983568E-04,2.8319080123804E-04,-6.0706301565874E-04,-0.018990068218419,-0.032529748770505,-0.021841717175414,-5.283835796993E-05,
        -4.7184321073267E-04, -3.0001780793026E-04, 4.7661393906987E-05,  -4.4141845330846E-06,  -7.2694996297594E-16, -3.1679644845054E-05,
        -2.8270797985312E-06,-8.5205128120103E-10, -2.2425281908E-06, -6.5171222895601E-07, -1.4341729937924E-13,-4.0516996860117E-07, -1.2734301741641E-09, -1.7424871230634E-10,
        -6.8762131295531E-19,1.4478307828521E-20, 2.6335781662795E-23,  -1.1947622640071E-23,1.8228094581404E-24, -9.3537087292458E-26
    };
    int ireg1[34] = { 0,0,0,0,0,0,0,0,1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, 29, 30, 31, 32 };
    int jreg1[34] = { -2, -1, 0,  1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0,  1, 3, 17, -4, 0, 6, -5, -2, 10, -8, -11,-6, -29, -31, -38, -39, -40, -41 };

    double tau, piDen, gammapireg1, volreg1;
    int igammapireg;
    double rgas_water = 461.526;    //'gas constant in J/(kg K)

    tau = 1386.0 / tk;
    piDen = 0.1 * pBar / 16.53;
    gammapireg1 = 0.0;

    for (igammapireg = 0; igammapireg < 34; igammapireg++) {
        gammapireg1 = gammapireg1 - nreg1[igammapireg] * ireg1[igammapireg] * pow(7.1 - piDen, ireg1[igammapireg] - 1.0) * pow(tau - 1.222, jreg1[igammapireg]);
    }

    volreg1 = rgas_water * tk * piDen * gammapireg1 / (pBar * 100000.0);

    // 函数通过 return 语句返回值
    return 1.0 / volreg1;
}

double fPP(double tk, double p0, double p1, double p2, double p3, double p4, double p5) {
    double Tr = 298.15;

    // 计算fPP
    double temp = p0 + p1 * (tk - Tr) + p2 * (pow(tk, 2) - pow(Tr, 2)) +
        p3 * (1.0 / tk - 1.0 / Tr) + p4 * log(tk / Tr) +
        p5 * (1.0 / pow(tk, 2) - 1.0 / pow(Tr, 2));

    return temp;
}

double fPZ6(double tk, double z1, double z2, double z3, double z4, double z5, double z6) {
    //fPZ6 = z1 + z2 * 100#  / tk + z3 * 0.01 * tk + z4 * 0.0001 * tk ^ 2 + z5 * 10# / (tk - 227) + z6 * Log(tk)
    double temp = z1 + z2 * 100.0 / tk + z3 * 0.01 * tk + z4 * 0.0001 * pow(tk, 2) + z5 * 10.0 / (tk - 227.0) + z6 * log(tk);
    return temp;
}

double fHolmes(double tk, double pBar, double q1, double q2, double q3, double q4, double q5, double q6, double q7,
    double q8, double q9, double q10, double q11, double q12, double q13, double q14, double q15, double q16, double q17)
{
    double temp = q1 + 1.0 / 2.0 * q2 * tk + 1.0 / 6.0 * q3 * pow(tk, 2) +
        1.0 / 12.0 * q4 * pow(tk, 3) +
        1.0 / 6.0 * q5 * pow(tk, 2) * (log(tk) - 5.0 / 6.0) +
        q6 * (tk / 2.0 + 3.0 * pow(227.0, 2) / 2.0 / tk +
            227.0 * (tk - 227.0) / tk * log(tk - 227.0)) +
        q7 * (2.0 * (647.0 - tk) / tk + 1.0) * log(647.0 - tk) +
        pBar * (q8 + q9 / tk + q10 * tk + q11 * pow(tk, 2) +
            q12 / (tk - 227.0) + q13 / (647.0 - tk)) +
        pow(pBar, 2) * (q14 + q15 / tk + q16 * tk + q17 * pow(tk, 2));

    return temp;
}

double fDuan1(double tk, double pBar, double q1, double q2, double q3, double q4, double q5, double q6, double q7, double q8) {
    double temp = q1 + q2 * tk + q3 / tk + q4 / (tk - 210.0) +
        q5 / (647.0 - tk) +
        q6 * pow(tk - 443.0, 3) / 3.0 +
        q7 * (pBar - 1.0) +
        q8 * pow(pBar - 1.0, 2) / 2.0;

    return temp;
}

double  fDuan2(double tk, double pBar, double q1, double q2, double q3, double q4, double q5, double q6, double q7,
    double q8, double q9, double q10, double q11) {
    double temp = q1 + q2 * tk + q3 / tk + q4 * pow(tk, 2) +
        q5 / (630.0 - tk) + q6 * pBar +
        q7 * pBar * log(tk) + q8 * pBar / tk +
        q9 / (630.0 - tk) +
        q10 * pow(pBar, 2) / pow(630.0 - tk, 2) +
        q11 * tk * log(pBar);
    return temp;
}

//'from Duan and Li 2008
double  fPZ_DL_6(double tk, double patm, double z1, double z2, double z3, double z4, double z5, double z6) {
    double temp = z1 + z2 * tk + z3 / (tk - 210.0) +
        z4 / (patm * 1.013254 + 100.0) +
        z5 * tk * log(patm * 1.013254) +
        z6 * tk * patm * 1.013254 * log(patm * 1.013254);
    return temp;
}


void C2_Pitzer2019(double tk, double tc, double pBar, double patm) {
    V0TP(tk, pBar);

    // b0(iH, iCl) = 0.1769 + -0.0914 * Log(fH2ODensity(tk, pBar) / 997.048) + 0 * (fH2ODensity(tk, pBar) - 997.048) / 1 + -0.0004034 * (tc - 25) / 1 + 0.000062 * (pBar - 1) / 10 'Holmes et al. 1987_Model I_BP
    b0[iH][iCl] = 0.1769 - 0.0914 * log(fH2ODensity(tk, pBar) / 997.048) + 0 * (fH2ODensity(tk, pBar) - 997.048) / 1 - 0.0004034 * (tc - 25) / 1 + 0.000062 * (pBar - 1) / 10;

    // b0(iH, iSO4) = fPP(tk, 0.08198, -0.17932, 0.000106, 4655, 49.798, 0)
    b0[iH][iSO4] = fPP(tk, 0.08198, -0.17932, 0.000106, 4655, 49.798, 0);

    // b0(iNa, iOH) = 276.33247 + -7310.7724 / (tk) + -49.35997 * Log(tk) + 0.11070737 * (tk) + -0.000041248335 * (tk) ^ 2 + 11.931122 / (tk - 227) + 1.6386916 / (647 - tk) 'from Pabalan and Pitzer 1987
    b0[iNa][iOH] = 276.33247 - 7310.7724 / tk - 49.35997 * log(tk) + 0.11070737 * tk - 0.000041248335 * tk * tk + 11.931122 / (tk - 227.0) + 1.6386916 / (647.0 - tk);

    // q1 = -656.81518: q2 = 24.8691295: q3 = 0.000053812752667: q4 = -5.58874699E-08: q5 = 6.589326333E-12
    double q1 = -656.81518;
    double q2 = 24.8691295;
    double q3 = 0.000053812752667;
    double q4 = -5.58874699e-08;
    double q5 = 6.589326333e-12;
    // q6 = -4.4640952: q7 = 0.01110991383: q8 = -2.657339906E-07: q9 = 1.746006963E-10: q10 = 1.0462619E-14
    double q6 = -4.4640952;
    double q7 = 0.01110991383;
    double q8 = -2.657339906e-07;
    double q9 = 1.746006963e-10;
    double q10 = 1.0462619e-14;
    // q11 = -0.000005307012889: q12 = 8.634023325E-10: q13 = -4.1785962E-13
    double q11 = -0.000005307012889;
    double q12 = 8.634023325e-10;
    double q13 = -4.1785962e-13;
    // q14 = -1.579365943: q15 = 0.002202282079: q16 = -1.310550324E-07: q17 = -6.381368333E-11
    double q14 = -1.579365943;
    double q15 = 0.002202282079;
    double q16 = -1.310550324e-07;
    double q17 = -6.381368333e-11;
    // q18 = 9.706578079: q19 = -0.02686039622: q20 = 0.00001534474401: q21 = -3.215398267E-09
    double q18 = 9.706578079;
    double q19 = -0.02686039622;
    double q20 = 0.00001534474401;
    double q21 = -3.215398267e-09;
    b0[iNa][iCl] = q1 / tk + q2 + q3 * pBar + q4 * pow(pBar, 2) + q5 * pow(pBar, 3) + q6 * log(tk) + (q7 + q8 * pBar + q9 * pow(pBar, 2) + q10 * pow(pBar, 3)) * tk
        + (q11 + q12 * pBar + q13 * pow(pBar, 2)) * pow(tk, 2) + (q14 + q15 * pBar + q16 * pow(pBar, 2) + q17 * pow(pBar, 3)) / (tk - 227)
        + (q18 + q19 * pBar + q20 * pow(pBar, 2) + q21 * pow(pBar, 3)) / (680 - tk);

    // b0(iNa, iHCO3) = fPP(tk, 0.02837, -0.0103, 0.00000782, -579.52, 0, 0)
    b0[iNa][iHCO3] = fPP(tk, 0.02837, -0.0103, 0.00000782, -579.52, 0, 0);

    // b0(iNa, iCO3) = fPP(tk, 0.03856, 0.00128, 0.0000056, -1986.1, -7.5408, 0)
    b0[iNa][iCO3] = fPP(tk, 0.03856, 0.00128, 0.0000056, -1986.1, -7.5408, 0);

    // %based on Pabalan & Pitzer 1988 and Dai 20151112   - 注意涉及到除法必须给数字加小数
    q1 = 0.60955633; q2 = -16.090797; q3 = -0.10932828; q4 = 0.00025321479; q5 = -0.000000099384034; q6 = 0.040107638; q7 = 0.021711348; q8 = 0.001738512; q9 = 0.001722469; q10 = -0.01255087; q11 = -0.0104936; q12 = 92308.895357; q13 = 963.974106;
    b0[iNa][iSO4] = q1 * pow(tk, 2) / 6.0 + q2 * tk / 2.0 + q3 * pow(tk, 2) * (log(tk) / 2.0 - 5.0 / 12.0) / 3.0 + q4 * pow(tk, 3) / 12.0 + q5 * pow(tk, 4) / 20.0 +
        q6 * (tk / 2.0 + 3 * pow(227, 2) / 2.0 / tk + 227 * (tk - 227) * log(tk - 227) / tk) - q7 * (tk / 2.0 + 3 * pow(647, 2) / 2.0 / tk - 647 * (647 - tk) * log(647 - tk) / tk) - q12 / tk - q9 * (pow(298.15, 2) / tk) + q13 + q11
        + 0.001 * (pBar - 200) * fPZ6(tk, 0.047242702840047, -5.38099320830228e-02, 1.18935539713076e-02, -1.47824630406919e-02, -0.733032327494331, 3.62293033499918e-02) + 0.000001 * pow((pBar - 200), 2) * fPZ6(tk, -5.32348171981031e-02, -0.101775451027319, 2.21588585653303e-02, 4.05106192359346e-03, 0.661622313615445, -2.06389886993178e-02);

    // b0(iK, iOH) = fPP(tk, 0.1298, 0.00721, -0.00000184, -863.72, -4.5653, 0)    'Kaasa's book... added by Dai
    b0[iK][iOH] = fPP(tk, 0.1298, 0.00721, -0.00000184, -863.72, -4.5653, 0);

    q1 = 0.04808; q2 = -758.48; q3 = -4.7062; q4 = 0.010072;  q5 = -0.0000037599; q6 = 0;
    b0[iK][iCl] = q1 + q2 * (1.0 / tk - 1.0 / 298.15) + q3 * log(tk / 298.15) + q4 * (tk - 298.15) + q5 * (pow(tk, 2) - pow(298.15, 2)) + q6 * log(tk - 260);

    // b0(iK, iHCO3) = fPP(tk, -0.01344, -0.21739, 0.0000921, 10020, 82.417, 0)    'Kaasa's book... added by Dai
    b0[iK][iHCO3] = fPP(tk, -0.01344, -0.21739, 0.0000921, 10020, 82.417, 0);

    // b0(iK, iCO3) = fPP(tk, 0.12765, 0.0151, -0.0000148, 454.53, 0, 0)    'Kaasa's book... added by Dai
    b0[iK][iCO3] = fPP(tk, 0.12765, 0.0151, -0.0000148, 454.53, 0, 0);

    // b0(iK, iSO4) = fPP(tk, 0, 0.00424, 0.00000000977, -606.85, -3.0789, 0) + 0.001 * (pBar - 200) * fPZ6(tk, -3.78606959461161E-02, -0.117878400340705, 0.274800180783657, -6.43463283494378E-02, -1.96395999732603, 2.16674763010344E-02) + 0.000001 * (pBar - 200) ^ 2 * fPZ6(tk, -0.285218812395032, -0.853413242667859, -0.228114332820375, 9.55341201750416E-02, 7.18488790131111, -0.11317948093169)
    b0[iK][iSO4] = fPP(tk, 0, 0.00424, 0.00000000977, -606.85, -3.0789, 0) + 0.001 * (pBar - 200) * fPZ6(tk, -3.78606959461161e-02, -0.117878400340705, 0.274800180783657, -6.43463283494378e-02, -1.96395999732603, 2.16674763010344e-02) + 0.000001 * pow((pBar - 200), 2) * fPZ6(tk, -0.285218812395032, -0.853413242667859, -0.228114332820375, 9.55341201750416e-02, 7.18488790131111, -0.11317948093169);

    // b0(iK, iHS) = fPP(tk, -0.337, 0, 0, 0, 0, 0)    'Kaasa's book... added by Dai
    b0[iK][iHS] = fPP(tk, -0.337, 0, 0, 0, 0, 0);

    q1 = 0; q2 = 0.00414544383; q3 = -0.0000276747461; q4 = 3.37946704e-08; q5 = 0;    // Holmes Ca-Cl term
    q6 = 0; q7 = 0.00118276629; q8 = 0.00126084149; q9 = -0.158424548; q10 = -0.0000032972643;
    q11 = 3.37768212e-09; q12 = 0.00241466763; q13 = -0.0229175172; q14 = 0; q15 = 0;
    q16 = -1.2497591e-10; q17 = 3.54502058e-13;
    b0[iCa][iCl] = fHolmes(tk, pBar, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, q14, q15, q16, q17);

    // b0(iCa, iHCO3) = fPP(tk, 0.25604, 1.9273, -0.000813, -89647, -733.92, 0)    'Kaasa's book... added by Dai
    b0[iCa][iHCO3] = fPP(tk, 0.25604, 1.9273, -0.000813, -89647, -733.92, 0);

    // added by GD 20191021
    b0[iCa][iSO4] = fHolmes(tk, pBar, 20.4885394057851, 0.182733456879443, -1.35163233606994e-03, 5.14105639711679e-07, 8.70101445510704e-05, -2.52953035845223e-02, -1.00755838424923, -0.057005580645621, 8.93006760077719, 1.34345737705482e-04, -1.1135110963711e-07, -0.211564774545699, 0.153436546399729, -8.89933782533691e-05, 1.13985897736532e-02, 2.25220460716217e-07, -1.85334832426362e-10);

    q1 = 0.405500216; q2 = 0.004145444; q3 = -0.000228457; q4 = -0.0000000633123; q5 = 0.0000401087;
    q6 = 0; q7 = -0.001712441; q8 = 0.001260841; q9 = -0.152128885; q10 = -0.00000346379;
    q11 = 0.00000000370249; q12 = 0.002414668; q13 = -0.022917517; q14 = 0; q15 = 0;
    q16 = -0.000000000124976; q17 = 3.05038e-13;
    b0[iMg][iCl] = fHolmes(tk, pBar, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, q14, q15, q16, q17);

    // b0(iMg, iHCO3) = fPP(tk, 0.0385, 0.8738, -0.00037, -40167, -330.82, 0)    'Kaasa's book... added by Dai
    b0[iMg][iHCO3] = fPP(tk, 0.0385, 0.8738, -0.00037, -40167, -330.82, 0);

    // ''From SSP2013
    q1 = -1.0282; q2 = 0.008479; q3 = -0.0000233667; q4 = 0.000000021575; q5 = 0.00068402; q6 = 0.21499;
    b0[iMg][iSO4] = q1 * (tk / 2 + pow(298, 2) / (2 * tk) - 298) + q2 * (pow(tk, 2) / 6
        + pow(298, 3) / (3 * tk) - pow(298, 2) / 2) + q3 * (pow(tk, 3) / 12 + pow(298, 4) / (4 * tk)
            - pow(298, 3) / 3) + q4 * (pow(tk, 4) / 20 + pow(298, 5) / (5 * tk) - pow(298, 4) / 4) + q5 * (298 - pow(298, 2) / tk) + q6
        + 0.001 * (pBar - 20) * fPZ6(tk, 4.69338657997024e-02, 0.1885487503942, -2.57681401400115e-02, -2.92174438263678e-03, -5.28701612182822e-02, 0.010175024490052) + 0.000001 * pow((pBar - 20), 2) * fPZ6(tk, 4.21601882371307e-02, 0.470528627913928, -1.38985610560422e-02, -8.14864903217854e-03, -0.63735169569643, 3.21374826812395e-03);

    q1 = -2.84735725; q2 = 0.061583644; q3 = -0.001004361; q4 = 0.0000000337947; q5 = 0.000140902;
    q6 = -0.002234499; q7 = 0.001182766; q8 = 0.001235248; q9 = -0.158424548; q10 = -0.00000329726;
    q11 = 0.00000000337768; q12 = 0.002414668; q13 = -0.022917517; q14 = 0; q15 = 0;
    q16 = -0.000000000124976; q17 = 3.54502e-13;
    b0[iBa][iCl] = fHolmes(tk, pBar, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, q14, q15, q16, q17);

    // b0(iBa, iSO4) = 0.536160193861812 * b0(iCa, iSO4) 'based on Dai 20151112
    b0[iBa][iSO4] = 0.536160193861812 * b0[iCa][iSO4];

    q1 = -137.411207; q2 = 2.111859321; q3 = -0.049522204; q4 = -0.00000719953; q5 = 0.007892885;
    q6 = -0.027503059; q7 = 1.086188061; q8 = 1.26948963407883e-03; q9 = -0.158424548; q10 = -0.00000329726;
    q11 = 0.00000000337768; q12 = 0.002414668; q13 = -0.022917517; q14 = 0; q15 = 0;
    q16 = -0.000000000124976; q17 = 3.54502e-13;
    b0[iSr][iCl] = fHolmes(tk, pBar, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, q14, q15, q16, q17);

    // b0(iSr, iSO4) = 0.78474727943517 * b0(iCa, iSO4)  'based on Dai 20151112
    b0[iSr][iSO4] = 0.78474727943517 * b0[iCa][iSO4];

    // b1(iH, iCl) = 0.2973 + 16.147 * Log(fH2ODensity(tk, pBar) / 997.048) + -0.0176 * (fH2ODensity(tk, pBar) - 997.048) / 1 + 0 * (tc - 25) / 1 + 0.00072 * (pBar - 1) / 10 'Holmes et al. 1987_Model I_BP
    b1[iH][iCl] = 0.2973 + 16.147 * log(fH2ODensity(tk, pBar) / 997.048) - 0.0176 * (fH2ODensity(tk, pBar) - 997.048) / 1 + 0 * (tc - 25) / 1 + 0.00072 * (pBar - 1) / 10;

    // b1(iNa, iOH) = 462.86977 + -10294.181 / (tk) + -85.960581 * Log(tk) + 0.23905969 * (tk) + -0.00010795894 * (tk) ^ 2 'from Pabalan and Pitzer 1987
    b1[iNa][iOH] = 462.86977 - 10294.181 / tk - 85.960581 * log(tk) + 0.23905969 * tk - 0.00010795894 * pow(tk, 2);

    // b1(iNa, iCl) = 119.31966 / tk - 0.48309327 + 0.0014068095 * tk - 4.2345814 / (tk - 227)
    b1[iNa][iCl] = 119.31966 / tk - 0.48309327 + 0.0014068095 * tk - 4.2345814 / (tk - 227);

    // b1(iNa, iHCO3) = fPP(tk, 0.0475, 0.14998, -0.0000574, -8814.8, -63.826, 0)
    b1[iNa][iHCO3] = fPP(tk, 0.0475, 0.14998, -0.0000574, -8814.8, -63.826, 0);

    // b1(iNa, iCO3) = fPP(tk, 1.5208, -0.66971, 0.000333, 19064, 204.49, 0)
    b1[iNa][iCO3] = fPP(tk, 1.5208, -0.66971, 0.000333, 19064, 204.49, 0);

    // '%based on Pabalan & Pitzer 1988
    q1 = 1.1040235; q2 = -25.758534; q3 = -0.20290775; q4 = 0.00053309441; q5 = -0.00000023576724; q6 = 0; q7 = 0.14455381; q8 = 0.005820066; q9 = 0.005512612; q10 = 0.703766; q11 = 0.690077; q12 = 363078.71668; q13 = 1926.602872;
    b1[iNa][iSO4] = q1 * pow(tk, 2) / 6.0 + q2 * tk / 2 + q3 * pow(tk, 2) * (log(tk) / 2.0 - 5.0 / 12.0) / 3.0 + q4 * pow(tk, 3) / 12.0 + q5 * pow(tk, 4) / 20.0 + q6 * (tk / 2.0 + 3 * pow(227, 2) / 2.0 / tk + 227 * (tk - 227) * log(tk - 227) / tk) - q7 * (tk / 2.0 + 3 * pow(647, 2) / 2.0 / tk - 647 * (647 - tk) * log(647 - tk) / tk) - q12 / tk - q9 * (pow(298.15, 2) / tk) + q13 + q11;

    q1 = 0.0476; q2 = 303.9; q3 = 1.066; q4 = 0; q5 = 0; q6 = 0.047;
    b1[iK][iOH] = fPP(tk, 0.32, 0.0634, -0.0000276, -3232.7, -24.573, 0);

    b1[iK][iCl] = q1 + q2 * (1 / tk - 1 / 298.15) + q3 * log(tk / 298.15) + q4
        * (tk - 298.15) + q5 * (pow(tk, 2) - pow(298.15, 2)) + q6 * log(tk - 260);

    // b1(iK, iHCO3) = fPP(tk, 0.0401, -0.31301, 0.000135, 13868, 116.4, 0)    'Kaasa's book... added by Dai
    b1[iK][iHCO3] = fPP(tk, 0.0401, -0.31301, 0.000135, 13868, 116.4, 0);

    // b1(iK, iCO3) = fPP(tk, 1.4248, 0.0977, -0.0000968, 3164, 0, 0)    'Kaasa's book... added by Dai
    b1[iK][iCO3] = fPP(tk, 1.4248, 0.0977, -0.0000968, 3164, 0, 0);

    // b1(iK, iSO4) = fPP(tk, 0.6179, 0.042, 0.0000169, -8714, -42.582, 0)    'Kaasa's book... added by Dai
    b1[iK][iSO4] = fPP(tk, 0.6179, 0.042, 0.0000169, -8714, -42.582, 0);

    // b1(iK, iHS) = fPP(tk, 0.884, 0, 0, 0, 0, 0)    'Kaasa's book... added by Dai
    b1[iK][iHS] = fPP(tk, 0.884, 0, 0, 0, 0, 0);

    q1 = 0; q2 = -0.16737337; q3 = 0.0195851174; q4 = 0.00000751975973; q5 = -0.00367501519; // Holmes Ca-Cl interaction
    q6 = -0.0239198164; q7 = 0; q8 = 0; q9 = 0; q10 = 0.00000107765583;
    q11 = -3.96914481e-09; q12 = 0; q13 = 0; q14 = 0; q15 = 0;
    q16 = 0; q17 = 0;
    b1[iCa][iCl] = fHolmes(tk, pBar, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, q14, q15, q16, q17);

    // b1(iCa, iHCO3) = fPP(tk, 0.30575, -0.28251, 0.000118, 14501, 110.86, 0)    'Kaasa's book... added by Dai
    b1[iCa][iHCO3] = fPP(tk, 0.30575, -0.28251, 0.000118, 14501, 110.86, 0);

    // Added by GD 20191021
    b1[iCa][iSO4] = fPZ6(tk, -58.3369742068748, 240.112816949155, 14.6657417194044, 0.874458506940693, 2.86581921261989, -11.9026193709332);

    q1 = 0; q2 = -0.16737337; q3 = 0.019728358; q4 = 0.00000753744; q5 = -0.003696071;
    q6 = -0.025038112; q7 = 0; q8 = 0; q9 = 0; q10 = 0.00000107766;
    q11 = -0.00000000396914; q12 = 0; q13 = 0; q14 = 0; q15 = 0; q16 = 0; q17 = 0;
    b1[iMg][iCl] = fHolmes(tk, pBar, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, q14, q15, q16, q17);

    // b1(iMg, iHCO3) = fPP(tk, 0.87005, -8.1001, 0.00344, 373490, 3064.9, 0)    'Kaasa's book... added by Dai
    b1[iMg][iHCO3] = fPP(tk, 0.87005, -8.1001, 0.00344, 373490, 3064.9, 0);

    // based on Phutela & Pitzer
    b1[iMg][iSO4] = -0.29596 * (tk / 2 + pow(298, 2) / (2 * tk) - 298) + 0.00094564 * (pow(tk, 2) / 6 + pow(298, 3) / (3 * tk) - pow(298, 2) / 2) + 0.01028 * (298 - pow(298, 2) / tk) + 3.3646;

    q1 = 0; q2 = -0.583253079; q3 = 0.062151491; q4 = 0.0000222759; q5 = -0.011539945;
    q6 = -0.078415969; q7 = 0; q8 = 0.000197778; q9 = 0; q10 = 0.00000107766;
    q11 = -0.00000000396914; q12 = 0; q13 = 0; q14 = 0; q15 = 0; q16 = 0; q17 = 0;
    b1[iBa][iCl] = fHolmes(tk, pBar, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, q14, q15, q16, q17);

    // b1(iBa, iSO4) = 0.536160193861812 * b1(iCa, iSO4) 'based on Dai 20151112
    b1[iBa][iSO4] = 0.536160193861812 * b1[iCa][iSO4];

    q1 = 0; q2 = -0.164948621; q3 = 0.01951972; q4 = 0.00000772188; q5 = -0.003675015;
    q6 = -0.023919816; q7 = 0; q8 = -0.000144388; q9 = 0; q10 = 0.00000107766;
    q11 = -0.00000000396914; q12 = 0; q13 = 0; q14 = 0; q15 = 0; q16 = 0; q17 = 0;
    b1[iSr][iCl] = fHolmes(tk, pBar, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, q14, q15, q16, q17);

    // b1(iSr, iSO4) = 0.78474727943517 * b1(iCa, iSO4)  'based on Dai 20151112
    b1[iSr][iSO4] = 0.78474727943517 * b1[iCa][iSO4];

    // b2(iCa, iCl) = -0.5 / Exp(-16.5 + 7150 / tk) 'Holmes
    b2[iCa][iCl] = -0.5 / exp(-16.5 + 7150 / tk);

    // b2(iMg, iCl) = b2(iCa, iCl)
    b2[iMg][iCl] = b2[iCa][iCl];

    // b2(iBa, iCl) = b2(iCa, iCl)
    b2[iBa][iCl] = b2[iCa][iCl];

    // b2(iSr, iCl) = b2(iCa, iCl)
    b2[iSr][iCl] = b2[iCa][iCl];

    // added by GD 20191021
    b2[iCa][iSO4] = fHolmes(tk, pBar, -3451.92417674836, -25.0510497042672, 0.52463858314361, -6.03901870423439e-04, -5.41906113936402e-02, 4.17388513213358, 96.8711207569075, -1.38634677934667, 452.409611072822, -3.37062199562341e-04, 4.42711788329082e-06, -15.5294298988312, -58.8816456334969, -4.86122268847504e-03, 0.458609072994063, 1.52062601592044e-05, -1.41724741011899e-08);

    // based on Phutela & Pitzer and P dependence from Dai 20151112
    b2[iMg][iSO4] = -13.764 * (tk / 2.0 + pow(298, 2) / (2 * tk) - 298) + 0.12121 * (pow(tk, 2) / 6.0 + pow(298, 3) / (3 * tk) - pow(298, 2) / 2) + -0.00027642 * (pow(tk, 3) / 12.0 + pow(298, 4) / (4 * tk) - pow(298, 3) / 3.0) + -0.21515 * (298 - pow(298, 2) / tk) + -32.743 +
        0.1 * (pBar - 20) * fPZ6(tk, 0.412262082080429, 0.400058114021816, 3.00888166907448e-02, -2.94302940870418e-02, -4.42981067782703, 5.66583267829718e-02) +
        0.0001 * pow((pBar - 20), 2) * fPZ6(tk, 0.234558885627083, -0.245867648958111, 6.55732421096113e-02, -2.74923796315413e-02, -4.12894058959444, 0.011915963054195);

    // b2(iBa, iSO4) = 0.536160193861812 * b2(iCa, iSO4) 'based on Dai 20151112
    b2[iBa][iSO4] = 0.536160193861812 * b2[iCa][iSO4];

    // b2(iSr, iSO4) = 0.78474727943517 * b2(iCa, iSO4)  'based on Dai 20151112
    b2[iSr][iSO4] = 0.78474727943517 * b2[iCa][iSO4];

    // CPhi(iH, iCl) = 2 * (0.000362 + 0 * Log(fH2ODensity(tk, pBar) / 997.048) + 0 * (fH2ODensity(tk, pBar) - 997.048) / 1 + -0.00003036 * (tc - 25) / 1 + 0 * (pBar - 1) / 10) 'Holmes et al. 1987_Model I_BP
    CPhi[iH][iCl] = 2 * (0.000362 + 0 * log(fH2ODensity(tk, pBar) / 997.048) + 0 * (fH2ODensity(tk, pBar) - 997.048) / 1 - 0.00003036 * (tc - 25) / 1 + 0 * (pBar - 1) / 10);

    // CPhi(iH, iSO4) = fPP(tk, 0.06375, -0.11362, 0.0000582, 3218.7, 34.424, 0)
    CPhi[iH][iSO4] = fPP(tk, 0.06375, -0.11362, 0.0000582, 3218.7, 34.424, 0);

    // CPhi(iNa, iOH) = -16.615961 + 444.59966 / (tk) + 2.9680772 * Log(tk) + -0.0067008449 * (tk) + 0.000002533892 * (tk) ^ 2 + -0.68923899 / (tk - 227) + -0.081156286 / (647 - tk) 'from Pabalan and Pitzer 1987
    CPhi[iNa][iOH] = -16.615961 + 444.59966 / tk + 2.9680772 * log(tk) - 0.0067008449 * tk + 0.000002533892 * pow(tk, 2) - 0.68923899 / (tk - 227) - 0.081156286 / (647 - tk);

    q1 = -6.1084589; q2 = 0.40217793; q3 = 0.000022902837; q4 = -0.075354649;
    q5 = 0.0001531767295; q6 = -0.000000090550901; q7 = -1.53860082e-08; q8 = 8.69266e-11;
    q9 = 0.353104136; q10 = -0.00043314252; q11 = -0.09187145529; q12 = 0.00051904777;
    CPhi[iNa][iCl] = q1 / tk + q2 + q3 * pBar + q4 * log(tk) + (q5 + q6 * pBar) * tk
        + (q7 + q8 * pBar) * pow(tk, 2) + (q9 + q10 * pBar) / (tk - 227) + (q11 + q12 * pBar) / (680 - tk);

    // based on Pabalan & Pitzer 1988
    q1 = -0.291330454580456; q2 = 7.47067054163403; q3 = 5.25526087271343e-02; q4 = -1.25440539335741e-04; q5 = 5.05817041380687e-08; q6 = -1.46616028524056e-02; q7 = -1.09759151043306e-02; q8 = -3.16065983167313e-04; q9 = -3.26114253370686e-04; q10 = 1.07722061259521e-02; q11 = 1.05906493888146e-02; q12 = -43691.7324087285; q13 = -441.068472680535;
    CPhi[iNa][iSO4] = q1 * pow(tk, 2) / 6.0 + q2 * tk / 2.0 + q3 * pow(tk, 2) * (log(tk) / 2.0 - 5.0 / 12.0) / 3.0 + q4 * pow(tk, 3) / 12.0 + q5 * pow(tk, 4) / 20.0 + q6 * (tk / 2 + 3 * pow(227, 2) / 2.0 / tk + 227 * (tk - 227) * log(tk - 227) / tk) - q7 * (tk / 2 + 3 * pow(647, 2) / 2.0 / tk - 647 * (647 - tk) * log(647 - tk) / tk) - q12 / tk - q9 * (pow(298.15, 2) / tk) + q13 + q11;

    // CPhi(iK, iOH) = fPP(tk, 0.0041, -0.00274, 0.00000104, 195.95, 1.2493, 0)    'Kaasa's book... added by Dai
    CPhi[iK][iOH] = fPP(tk, 0.0041, -0.00274, 0.00000104, 195.95, 1.2493, 0);

    q1 = -0.000788; q2 = 91.27; q3 = 0.58643; q4 = -0.001298; q5 = 0.00000049567; q6 = 0;
    CPhi[iK][iCl] = q1 + q2 * (1 / tk - 1 / 298.15) + q3 * log(tk / 298.15) + q4
        * (tk - 298.15) + q5 * (pow(tk, 2) - pow(298.15, 2)) + q6 * log(tk - 260);

    // CPhi(iK, iSO4) = fPP(tk, 0.00915, -0.000181, 0, -16.092, 0, 0)    'Kaasa's book... added by Dai
    CPhi[iK][iSO4] = fPP(tk, 0.00915, -0.000181, 0, -16.092, 0, 0);

    q1 = -0.131583284; q2 = 0; q3 = 0.000289257572; q4 = 0.000000128494802; q5 = -0.000056273068; // Holmes
    q6 = -0.000594574164; q7 = 0; q8 = -0.000000958297102; q9 = 0; q10 = 0;
    q11 = 6.34029223e-12; q12 = 0; q13 = 0; q14 = -5.60197799e-09; q15 = 0.0000017747878;
    q16 = 0; q17 = 0;
    CPhi[iCa][iCl] = fHolmes(tk, pBar, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, q14, q15, q16, q17);

    // added by GD 20191021
    CPhi[iCa][iSO4] = fHolmes(tk, pBar, 5.77985723123435, -0.121409849839165, -1.87628978264948e-04, -3.77175580664723e-07, 1.40342341987879e-04, 3.77829703206553e-03, 0.20380678760812, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0); // No pressure dependence, by GD

    q1 = -0.131583284; q2 = -0.000958991; q3 = 0.000341089; q4 = 0.000000128495; q5 = -0.0000644255;
    q6 = -0.00067376; q7 = 0.00079875; q8 = -0.00000118509; q9 = 0; q10 = 0; q11 = 6.34029e-12;
    q12 = 0; q13 = 0; q14 = -0.00000000560198; q15 = 0.00000177479; q16 = 0; q17 = 1.31968e-14;
    CPhi[iMg][iCl] = fHolmes(tk, pBar, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, q14, q15, q16, q17);

    // based on Phutela & Pitzer
    CPhi[iMg][iSO4] = 4 * (0.10541 * (tk / 2 + pow(298, 2) / (2 * tk) - 298) + -0.00089316 * (pow(tk, 2) / 6 + pow(298, 3) / (3 * tk) - pow(298, 2) / 2) + 0.00000251 * (pow(tk, 3) / 12 + pow(298, 4) / (4 * tk)
        - pow(298, 3) / 3) + -0.0000000023436 * (pow(tk, 4) / 20 + pow(298, 5) / (5 * tk) - pow(298, 4) / 4) + -0.000087899 * (298 - pow(298, 2) / tk) + 0.006993);

    q1 = -0.131583284; q2 = 0; q3 = 0.000247683; q4 = 0.0000000895556; q5 = -0.0000465877;
    q6 = -0.000630297; q7 = 0; q8 = 0.00000625234; q9 = 0; q10 = 0;
    q11 = 6.34029e-12; q12 = 0; q13 = 0; q14 = -0.00000000560198; q15 = 0.00000177479; q16 = 0; q17 = 0;
    CPhi[iBa][iCl] = fHolmes(tk, pBar, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, q14, q15, q16, q17);

    // CPhi(iBa, iSO4) = 0.536160193861812 * CPhi(iCa, iSO4) 'based on Dai 20151112
    CPhi[iBa][iSO4] = 0.536160193861812 * CPhi[iCa][iSO4];

    q1 = -0.131583284; q2 = 0.002686696; q3 = 0.0000837031; q4 = 0.0000000361626; q5 = -0.0000167901;
    q6 = -0.000594574; q7 = 0; q8 = -0.00000601422; q9 = 0; q10 = 0;
    q11 = 6.34029e-12; q12 = 0; q13 = 0; q14 = -0.00000000560198; q15 = 0.00000177479; q16 = 0; q17 = 0;
    CPhi[iSr][iCl] = fHolmes(tk, pBar, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, q14, q15, q16, q17);

    // CPhi(iSr, iSO4) = 0.78474727943517 * CPhi(iCa, iSO4)  'based on Dai 20151112
    CPhi[iSr][iSO4] = 0.78474727943517 * CPhi[iCa][iSO4];

    // Lnc(iCO2aq, iNa) = fPP(tk, 0.0803, 0.0955, -0.0000382, -4937.1, -38.454, 0)  ' Below LN values are from Kaasa
    Lnc[iCO2aq][iNa] = fPP(tk, 0.0803, 0.0955, -0.0000382, -4937.1, -38.454, 0);

    // Lnc(iCO2aq, iK) = fPP(tk, 0.04856, 0.0955, -0.0000382, -4937.1, -38.454, 0)
    Lnc[iCO2aq][iK] = fPP(tk, 0.04856, 0.0955, -0.0000382, -4937.1, -38.454, 0);

    // Lnc(iCO2aq, iMg) = fPP(tk, 0.14998, 0.00237, -0.0000008, 377.48, 0, 0)
    Lnc[iCO2aq][iMg] = fPP(tk, 0.14998, 0.00237, -0.0000008, 377.48, 0, 0);

    // Lnc(iCO2aq, iCa) = fPP(tk, 0.15267, 0.00203, -0.000000465, 370.68, 0, 0)
    Lnc[iCO2aq][iCa] = fPP(tk, 0.15267, 0.00203, -0.000000465, 370.68, 0, 0);

    // Lnc(iCO2aq, iBa) = fPP(tk, 0.15267, 0.00203, -0.000000465, 370.68, 0, 0)
    Lnc[iCO2aq][iBa] = fPP(tk, 0.15267, 0.00203, -0.000000465, 370.68, 0, 0);

    // Lnc(iCO2aq, iSr) = fPP(tk, 0.15267, 0.00203, -0.000000465, 370.68, 0, 0)
    Lnc[iCO2aq][iSr] = fPP(tk, 0.15267, 0.00203, -0.000000465, 370.68, 0, 0);

    // Lnc(iCO2aq, iFe) = fPP(tk, 0.15267, 0.00203, -0.000000465, 370.68, 0, 0)
    Lnc[iCO2aq][iFe] = fPP(tk, 0.15267, 0.00203, -0.000000465, 370.68, 0, 0);

    // Lna(iCO2aq, iCl) = fPP(tk, 0.01919, -0.00527, 0.00000164, 492.38, 2.7967, 0)
    Lna[iCO2aq][iCl] = fPP(tk, 0.01919, -0.00527, 0.00000164, 492.38, 2.7967, 0);

    // Lna(iCO2aq, iBr) = fPP(tk, 0.01919, -0.00527, 0.00000164, 492.38, 2.7967, 0)
    Lna[iCO2aq][iBr] = fPP(tk, 0.01919, -0.00527, 0.00000164, 492.38, 2.7967, 0);

    // Lna(iCO2aq, iSO4) = fPP(tk, 0.19003, 0.0295, -0.0000215, 1970.8, 0, 0)
    Lna[iCO2aq][iSO4] = fPP(tk, 0.19003, 0.0295, -0.0000215, 1970.8, 0, 0);

    // *********************************** Holmes iK, iCl term
    q1 = -0.0210289; q2 = 0.603967; q3 = 0.00367768; q4 = -0.00000705537; q5 = 0.00000000197968;
    q6 = -0.00247588; q7 = 0.14416; q8 = 0.000677136; q9 = 0.000656838; q10 = 0.04808;
    q11 = 0.050038; q12 = -2931.268116; q13 = -33.953143; q14 = 0; q15 = 0; q16 = 0; q17 = 0; q18 = 0;
    b0[iK][iCl] = q1 * pow(tk, 2) / 6 + q2 * tk / 2 + q3 * pow(tk, 2) * (log(tk) / 2 - 5.0 / 12.0) / 3 + q4 * pow(tk, 3) / 12 + q5 * pow(tk, 4) / 20 + q6 * (tk / 2 + 3 * pow(227, 2) / 2 / tk + 227 * (tk - 227) * log(tk - 227) / tk) - q7 * (2 * (647 - tk) * log(647 - tk) / tk + log(647 - tk)) - q12 / tk - q9 * (pow(298.15, 2) / tk) + q13 + q11;

    q1 = 0; q2 = 0; q3 = 0.0000000945016; q4 = -0.000000000290741; q5 = 0;
    q6 = 0.00326205; q7 = 0.000000839662; q8 = 0; q9 = -0.00000000441638; q10 = 6.71235e-12;
    q11 = 0; q12 = -0.0000442327; q13 = -0.000000000797437; q14 = 0; q15 = 4.12771e-12;
    q16 = -6.24996e-15; q17 = 0; q18 = 0.0000000416221;
    b0[iK][iCl] = b0[iK][iCl] + (q1 + q2 / tk + q3 * tk + q4 * pow(tk, 2) + q6 / (647 - tk)) * (pBar - 179) + (pow(pBar, 2) - pow(179, 2)) / 2 * (q7 + q8 / tk + q9 * tk + q10 * pow(tk, 2) + q12 / (647 - tk)) + (pow(pBar, 3) - pow(179, 3)) / 3 * (q13 + q14 / tk + q15 * tk + q16 * pow(tk, 2) + q18 / (647 - tk));

    q1 = 0.220813; q2 = -4.61849; q3 = -0.0410116; q4 = 0.000110445; q5 = -0.0000000473196;
    q6 = -0.027412; q7 = 0.332883; q8 = 0.000967854; q9 = 0.000967854; q10 = 0.218752;
    q11 = 0.218752; q12 = 6353.355434; q13 = 193.004059; q14 = 0; q15 = 0; q16 = 0; q17 = 0; q18 = 0;
    b1[iK][iCl] = q1 * pow(tk, 2) / 6 + q2 * tk / 2 + q3 * pow(tk, 2) * (log(tk) / 2 - 5.0 / 12.0) / 3 + q4 * pow(tk, 3) / 12 + q5 * pow(tk, 4) / 20 + q6 * (tk / 2 + 3 * pow(227, 2) / 2 / tk + 227 * (tk - 227) * log(tk - 227) / tk) - q7 * (2 * (647 - tk) * log(647 - tk) / tk + log(647 - tk)) - q12 / tk - q9 * (pow(298.15, 2) / tk) + q13 + q11;

    q1 = 0; q2 = 0.000764891; q3 = 0; q4 = -0.0000000112131; q5 = 1.72256e-11;
    q6 = 0; q7 = -0.00571188; q8 = -0.0000412364; q9 = -0.0000412364; q10 = -0.000394;
    q11 = -0.000394; q12 = 28.17218; q13 = -0.125567; q14 = 0; q15 = 0; q16 = 0; q17 = 0; q18 = 0;
    CPhi[iK][iCl] = 2 * (q1 * pow(tk, 2) / 6 + q2 * tk / 2 + q3 * pow(tk, 2) * (log(tk) / 2 - 5.0 / 12.0) / 3 + q4 * pow(tk, 3) / 12 + q5 * pow(tk, 4) / 20 + q6 * (tk / 2 + 3 * pow(227, 2) / 2 / tk + 227 * (tk - 227) * log(tk - 227) / tk) - q7 * (2 * (647 - tk) * log(647 - tk) / tk + log(647 - tk)) - q12 / tk - q9 * (pow(298.15, 2) / tk) + q13 + q11);

    // '%Pabalan and Pitzer
    Taap[iCl][iSO4] = 0.03; Taap[iSO4][iCl] = 0.03;

    // '%based on Pabalan & Pitzer 1988
    Yaapc[iCl][iSO4][iNa] = -0.016958 + 3.13544 / tk + 0.0000216352 * tk + -131254.0 / pow((647 - tk), 4);
    Yaapc[iSO4][iCl][iNa] = -0.016958 + 3.13544 / tk + 0.0000216352 * tk + -131254.0 / pow((647 - tk), 4);

    // 'Duan and Li 2008 Carbonate
    q1 = 0.0661; q2 = 0; q3 = 0; q4 = 0; q5 = 0; q6 = 0.0000000375951; q7 = 0; q8 = 0;
    b0[iNa][iHCO3] = fDuan1(tk, pBar, q1, q2, q3, q4, q5, q6, q7, q8);

    q1 = 0.5153; q2 = -0.0005991; q3 = 0; q4 = -25.81; q5 = -2.659; q6 = 0; q7 = 0.0000875; q8 = -0.0000000266;
    b0[iNa][iCO3] = fDuan1(tk, pBar, q1, q2, q3, q4, q5, q6, q7, q8);

    q1 = -4.116; q2 = 0.006309; q3 = 924; q4 = -52.02; q5 = -80.26; q6 = 0; q7 = 0.0001634; q8 = -0.000000139;
    b1[iNa][iHCO3] = fDuan1(tk, pBar, q1, q2, q3, q4, q5, q6, q7, q8);

    q1 = 2.044; q2 = -0.004303; q3 = 0; q4 = -25.45; q5 = 361.8; q6 = 0; q7 = 0; q8 = 0;
    b1[iNa][iCO3] = fDuan1(tk, pBar, q1, q2, q3, q4, q5, q6, q7, q8);

    CPhi[iNa][iHCO3] = 0;

    q1 = -0.0914; q2 = 0; q3 = 0; q4 = 6.482; q5 = 8.048; q6 = 0; q7 = -0.0000289; q8 = 0;
    CPhi[iNa][iCO3] = fDuan1(tk, pBar, q1, q2, q3, q4, q5, q6, q7, q8);

    q1 = -0.2739092216; q2 = 0.0007399855859; q3 = 55.5213285; q4 = 0; q5 = 0; q6 = 0; q7 = 0;
    q8 = 0.005683638727; q9 = -0.0008009093476; q10 = 0; q11 = -0.0000174562027;
    Lnc[iCO2aq][iNa] = fDuan2(tk, pBar, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11);
    Lna[iCO2aq][iCl] = 0;

    q1 = -0.01665719188; q2 = 0.0000013916186; q3 = 0; q4 = 0; q5 = 0;
    q6 = 0; q7 = 0; q8 = -0.001873812115; q9 = -0.001577400757; q10 = 0; q11 = 0;
    zeta[iCO2aq][iNa][iCl] = fDuan2(tk, pBar, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11);

    // based on Dai 20151112
    Tccp[iNa][iCa] = fPZ_DL_6(tk, patm, -0.147510859131801, 9.76229611297803e-04, -4.82248417271452, -2.06519624782849, -1.51280704258336e-05, -1.04898649291008e-08);
    Tccp[iCa][iNa] = Tccp[iNa][iCa];

    Tccp[iNa][iBa] = 1.14420890430833 * Tccp[iNa][iCa];
    Tccp[iBa][iNa] = Tccp[iNa][iBa];
    Tccp[iNa][iSr] = 1.89691943139874 * Tccp[iNa][iCa];
    Tccp[iSr][iNa] = Tccp[iNa][iSr];

    Yccpa[iNa][iCa][iCl] = fPZ_DL_6(tk, patm, 2.26910136864354e-02, -2.43069462162211e-04, 2.40163930955923, 0.900545162398897, 9.17559621916219e-06, 1.11650476459493e-09);
    Yccpa[iCa][iNa][iCl] = Yccpa[iNa][iCa][iCl];

    Yccpa[iNa][iBa][iCl] = 1.29191707034912 * Yccpa[iNa][iCa][iCl];
    Yccpa[iBa][iNa][iCl] = Yccpa[iNa][iBa][iCl];
    Yccpa[iNa][iSr][iCl] = 2.32054171945995 * Yccpa[iNa][iCa][iCl];
    Yccpa[iSr][iNa][iCl] = Yccpa[iNa][iSr][iCl];

    // based on Dai 20151112
    Yaapc[iCl][iHCO3][iNa] = fPZ_DL_6(tk, patm, 3.29114286722895e-02, -9.905848162123e-05, -1.1319579697159, -6.66253293573831e-02, -1.209935715211e-05, 5.05498885272366e-09);
    Yaapc[iHCO3][iCl][iNa] = Yaapc[iCl][iHCO3][iNa];

    Yaapc[iSO4][iHCO3][iNa] = fPZ6(tk, -2.30508813504961, 10.9078063774135, 0.173488671930141, 2.38184731898223e-02, -10.9432999062871, 6.12652697241691e-03);
    Yaapc[iHCO3][iSO4][iNa] = Yaapc[iSO4][iHCO3][iNa];

    // '%based on IUPAC
    Lnc[iCO2aq][iCa] = 0.2045 - 0.000284 * (tk - 298.15);
    b0[iCa][iHCO3] = -3.7313 + 1371.42 / tk - 57330.0 / pow(tk, 2);
    b1[iCa][iHCO3] = 4.3005 - 2819.46 / tk + 483720.0 / pow(tk, 2);
    b0[iMg][iHCO3] = -1.9113 + 769.53 / tk - 57330.0 / pow(tk, 2);
    b1[iMg][iHCO3] = 14.3043 - 5590.6 / tk + 483720.0 / pow(tk, 2);

    Lna[iCO2aq][iHCO3] = 0;
    Lna[iCO2aq][iCO3] = 0;
    Lnn[iCO2aq][iCO2aq] = 0;

    Tccp[iNa][iCa] = fPZ_DL_6(tk, patm, -9.00319213981419e-02, 3.32942474038081e-04, -1.71678275211677, 6.96629148452292, 3.53927512842587e-05, -1.39564921696874e-09);
    Tccp[iCa][iNa] = Tccp[iNa][iCa];

    Yccpa[iNa][iCa][iCl] = fPZ_DL_6(tk, patm, -2.03123697088405e-02, 6.79617962918556e-05, 1.46202927499665, -2.31987379795446, -1.23017318725963e-05, -2.5328414139378e-09);

    Tccp[iNa][iCa] = fPZ_DL_6(tk, PsatH2O(tk), -9.00319213981419e-02, 3.32942474038081e-04, -1.71678275211677, 6.96629148452292, 3.53927512842587e-05, -1.39564921696874e-09);

    Tccp[iNa][iCa] = Tccp[iNa][iCa] + -0.00000025 * tk * (patm - PsatH2O(tk));

    Tccp[iCa][iNa] = Tccp[iNa][iCa];
    Yccpa[iNa][iCa][iCl] = fPZ_DL_6(tk, PsatH2O(tk), -2.03123697088405e-02, 6.79617962918556e-05, 1.46202927499665, -2.31987379795446, -1.23017318725963e-05, -2.5328414139378e-09);

    Yccpa[iNa][iCa][iCl] = Yccpa[iNa][iCa][iCl] + 0.00000001 * tk * (patm - PsatH2O(tk));

    Yccpa[iCa][iNa][iCl] = Yccpa[iNa][iCa][iCl];

    // ''===========================simContext.H2Saq coefficients
    q1 = -0.10242; q2 = 0.000322; q3 = 27.88934; q4 = 0; q5 = 0; q6 = 0; q7 = 0;            // Noe that q1-q3 are fitted by Dai based on Barrett's data, q7-q11 are set to equal to simContext.CO2aq parameters
    q8 = 0.005683638727; q9 = -0.0008009093476; q10 = 0; q11 = -0.0000174562027;
    Lnc[iH2Saq][iNa] = fDuan2(tk, pBar, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11);

    Lna[iH2Saq][iCl] = 0;
    q1 = 0.06745; q2 = -0.000104; q3 = -13.71127; q4 = 0; q5 = 0;
    q6 = 0; q7 = 0; q8 = -0.001873812115; q9 = -0.001577400757; q10 = 0; q11 = 0;
    zeta[iH2Saq][iNa][iCl] = fDuan2(tk, pBar, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11);
    // For m = 1 To NumCat  'fbtermma calculation in density calculation, set parameter for sulfide = parameters for carbonate
    int m;
    for (m = 0; m < NumCat; m++) {
        b0[m][iSion] = b0[m][iCO3];
        b1[m][iSion] = b1[m][iCO3];
        b2[m][iSion] = b2[m][iCO3];
        CPhi[m][iSion] = CPhi[m][iCO3];
    }

    // Dai fit some Fe related paramters 20161116
    // Xin 12/01/2020, ensures gc(iFe) is reasonable at all 10 m
    b0[iFe][iCl] = 0.05 + 0.001 * (tk - 298.15);
    b1[iFe][iCl] = 5.639834738 - 7.21091902166862e-03 * tk;
    b2[iFe][iCl] = -36.51133921 - 6.59720850886405e-03 * tk;
    CPhi[iFe][iCl] = 0.05;
    b0[iFe][iHS] = 1.999536302;
    b1[iFe][iHS] = -12.2860517351668;
    b2[iFe][iHS] = -4.94686944425638e-05;
    CPhi[iFe][iHS] = 0.5;
    b0[iFe][iSO4] = -0.999999943766306 - 1.99999911511231e-03 * (tk - 298.15);
    b1[iFe][iSO4] = -36.1108185261268 + 0.22019249845693 * tk - 2.23069409321093e-04 * pow(tk, 2);
    b2[iFe][iSO4] = -0.5 * pow(10, (-0.8398 - 0.895 * log10(tk) + 0.012 * tk));
    CPhi[iFe][iSO4] = -0.05;
    b0[iFe][iAc] = -0.5; b1[iFe][iAc] = -5; b2[iFe][iAc] = 0; CPhi[iFe][iAc] = -0.05;

}



void fgammaN() {
    double dSIMeOHcal, dSIMeOHBar, dSIMeOHHal;
    double dSIMEGcal, dSIMEGBar, dSIMEGHal;

    if (xMeOH > 0) {
        // Methanol calculations
        gNNeut[iCO2aq] = pow(10.0, ((-6.029 + 1444.9 / TK) * xMeOH - 1.17 * pow(xMeOH, 2)));

        dSIMeOHcal = (108.0317 - 4898.9661 / TK - 15.0753 * log(TK) - 0.8632 * ISt / (1.0 + sqrt(ISt))) * xMeOH - 2.3459 * pow(xMeOH, 2);

        gNAn[iHCO3] = pow(pow(10.0, dSIMeOHcal), (1.0 / 6.0));
        gNCat[iCa] = pow(pow(10.0, dSIMeOHcal), (2.0 / 3.0));
        gNAn[iCO3] = gNCat[iCa];
        gNAn[iAc] = pow(10.0, ((2.8478 - 0.2753 * ISt) * xMeOH));
        gNNeut[iHAcaq] = 1.0;
        gNNeut[iH2Saq] = gNNeut[iCO2aq];
        gNAn[iHS] = gNAn[iHCO3];

        dSIMeOHBar = (112.1 - 2094.5 / TK - 16.47 * log(TK)) * xMeOH + (8.437 - 5721.4 / TK) * pow(xMeOH, 2);
        gNMean[iBaSO4] = pow(pow(10.0, dSIMeOHBar), 0.5);

        gNMean[iCaSO42H2O] = pow(10.0, ((-5.6666 + 30.0535 / (1.0 + sqrt(ISt))) * xMeOH - (47.31 / pow(1.0 + sqrt(ISt), 2)) * pow(xMeOH, 2)));
        gNMean[iCaSO42H2O] = sqrt(gNMean[iCaSO42H2O]); // ^0.5 equivalent to sqrt()

        gNMean[iSrSO4] = pow(10.0, ((4.963 - 0.206 * ISt) * xMeOH - 4.017 * pow(xMeOH, 2)));

        dSIMeOHHal = (-107.31 + 6661.08 / TK + 15.3994 * log(TK)) * xMeOH + (2.6705 - 1392.3 / TK) * pow(xMeOH, 2);
        gNMean[iNaCl] = pow(pow(10.0, dSIMeOHHal), 0.5);

        aNH2O = 1.0 - 0.5403 * xMeOH - 0.4339 * pow(xMeOH, 2);
        gNCat[iH] = pow(10.0, ((-1.0667 - 0.304 * ISt) * xMeOH));
        gNAn[iOH] = gNCat[iH];
    }

    if (xMEG > 0) {
        // MEG calculations
        gNNeut[iCO2aq] = pow(10.0, ((-2.954 + 691.6 / TK - 0.382 * ISt / (1.0 + sqrt(ISt))) * xMEG - (151.9 / TK) * pow(xMEG, 2) - 0.67 * pow(xMEG, 4)));

        dSIMEGcal = (180.9605 - 9190.3924 / TK - 25.7158 * log(TK) - 1.7492 * ISt / pow(1.0 + sqrt(ISt), 0.5)) * xMEG + (-5.4646 + 3.8981 * ISt / (1.0 + sqrt(ISt))) * pow(xMEG, 2);

        gNAn[iHCO3] = pow(pow(10.0, dSIMEGcal), (1.0 / 6.0));
        gNCat[iCa] = pow(pow(10.0, dSIMEGcal), (2.0 / 3.0));
        gNAn[iCO3] = gNCat[iCa];
        gNAn[iAc] = pow(10.0, ((1.9008 - 0.2753 * ISt) * xMEG));
        gNNeut[iHAcaq] = 1.0;
        gNNeut[iH2Saq] = gNNeut[iCO2aq];
        gNAn[iHS] = gNAn[iHCO3];

        dSIMEGBar = (172.18 - 8125.4 / TK - 24.527 * log(TK)) * xMEG - 5.8668 * pow(xMEG, 2);
        gNMean[iBaSO4] = pow(pow(10.0, dSIMEGBar), 0.5);

        gNMean[iCaSO42H2O] = pow(10.0, ((10.5486 - 1767.92 / TK - 3.4044 * ISt / pow(1.0 + sqrt(ISt), 0.5) + 1.1963 * ISt / (1.0 + sqrt(ISt))) * xMEG + (-5.2799 + 6.5046 * ISt / (1.0 + sqrt(ISt))) * pow(xMEG, 2)));
        gNMean[iCaSO42H2O] = sqrt(gNMean[iCaSO42H2O]);

        gNMean[iCaSO4] = pow(10.0, ((7.585 - 1277.76 / TK + 3.1211 * ISt / pow(1.0 + sqrt(ISt), 2)) * xMEG + (-3.543 - 1.6694 * ISt / (1.0 + sqrt(ISt))) * pow(xMEG, 2)));
        gNMean[iCaSO4] = sqrt(gNMean[iCaSO4]);

        gNMean[ihemiCaSO4] = gNMean[iCaSO4];
        gNMean[iSrSO4] = gNMean[iBaSO4];

        dSIMEGHal = (-1.8742 + 911.0668 / TK) * xMEG + (5.3434 - 2939.27 / TK) * pow(xMEG, 2);
        gNMean[iNaCl] = pow(pow(10.0, dSIMEGHal), 0.5);

        aNH2O = 1.0 - 1.32 * xMEG + 0.33 * pow(xMEG, 2);
        gNCat[iH] = pow(10.0, ((0.1806 - 0.2456 * ISt) * xMEG + 1.6212 * pow(xMEG, 2.3622)));
        gNAn[iOH] = gNCat[iH];
    }
}



/**
 * @brief 计算Pitzer活性系数（温度、压力、离子强度版本）
 *
 * 此函数使用温度(TK)、温度(摄氏度 TC)、压力(PBar)和大气压(Patm)计算阴离子、阳离子和中性化合物的活性系数。
 * 基于Pitzer离子相互作用模型，包含二元、三元和中性项的计算。还计算水的渗透系数和活性，以及特定络合物的Bdot参数。
 * 注意：函数依赖外部数组和全局变量（如mc, ma, b0, b1等），这些需在调用前定义。
 * 特殊处理：对于某些离子对（如Na-SO4, K-SO4, Ca-SO4）使用自定义alpha值。
 */
void C2_PitzerActCoefs_T_P_ISt(DoubleVec& gNeut, double& aH2O, double tk, double tc, double pBar, double patm) {

    if (xMeOH > 0 || xMEG > 0) {
        //mt = fgammaN();//在该函数中只调用了一次
        fgammaN();
    }

    int a, c, m;

    C2_Pitzer2019(tk, tc, pBar, patm);  //此函数不修改这四个值

    double U1 = 342.79;
    double U2 = -0.0050866;
    double U3 = 0.0000009469;
    double U4 = -2.0525;
    double U5 = 3115.9;
    double U6 = -182.89;
    double U7 = -8032.5;
    double U8 = 4214200;
    double U9 = 2.1417;
    double D1000 = U1 * exp(U2 * tk + U3 * pow(tk, 2));
    double cc = U4 + U5 / (U6 + tk);
    double b = U7 + U8 / tk + U9 * tk;
    double Dielec = D1000 + cc * log((b + pBar) / (b + 1000));

    a1 = 16.945176;

    // dens0 等系数用于密度计算，但后续直接调用 fH2ODensity(tk, pBar)，因此此处计算被覆盖。
    dens = fH2ODensity(tk, pBar);

    // APhi 计算：Debye-Hückel 参数
    double APhi = (1.0 / 3.0) * pow((2 * pi * NAv * dens), 0.5) * pow((eElec * eElec / (4 * pi * eps0 * Dielec * kBoltz * tk)), 1.5);

    // 计算 gX 和 gpX 函数（依赖离子强度 ISt）
    double X14 = 1.4 * sqrt(ISt);  // 对于 2:(-2) 对或离子
    double gX14 = 2 * (1 - (1 + X14) * exp(-X14)) / (X14 * X14);
    double gpX14 = -2 * (1 - (1 + X14 + 0.5 * X14 * X14) * exp(-X14)) / (X14 * X14);

    double X20 = 2 * sqrt(ISt);  // 对于 1:(-2), (-1):2 或 1:(-1) 对
    double gX20 = 2 * (1 - (1 + X20) * exp(-X20)) / (X20 * X20);
    double gpX20 = -2 * (1 - (1 + X20 + 0.5 * X20 * X20) * exp(-X20)) / (X20 * X20);

    double X12 = 12 * sqrt(ISt);
    double gX12 = 2 * (1 - (1 + X12) * exp(-X12)) / (X12 * X12);
    double gpX12 = -2 * (1 - (1 + X12 + 0.5 * X12 * X12) * exp(-X12)) / (X12 * X12);

    // 添加 by GD 20191021：针对 Ca-SO4
    double xCaSO4 = 32 * APhi * sqrt(ISt);
    double gXCaSO4 = 2 * (1 - (1 + xCaSO4) * exp(-xCaSO4)) / (xCaSO4 * xCaSO4);
    double gpXCaSO4 = -2 * (1 - (1 + xCaSO4 + 0.5 * xCaSO4 * xCaSO4) * exp(-xCaSO4)) / (xCaSO4 * xCaSO4);

    // 计算 JX 和 JpX 函数（同电荷但不同离子，如 1:2 或 -1:-2）
    X12 = 6 * 1 * 2 * APhi * sqrt(ISt);
    double JX12 = X12 / (4 + 4.581 * pow(X12, -0.7237) * exp(-0.012 * pow(X12, 0.528)));
    X12 = 1.001 * X12;
    double JX12delta = X12 / (4 + 4.581 * pow(X12, -0.7237) * exp(-0.012 * pow(X12, 0.528)));
    X12 = X12 / 1.001;
    double JpX12 = (JX12delta - JX12) / (0.001 * X12);

    double X11 = 6 * 1 * 1 * APhi * sqrt(ISt);
    double JX11 = X11 / (4 + 4.581 * pow(X11, -0.7237) * exp(-0.012 * pow(X11, 0.528)));
    X11 = 1.001 * X11;
    double JX11delta = X11 / (4 + 4.581 * pow(X11, -0.7237) * exp(-0.012 * pow(X11, 0.528)));
    X11 = X11 / 1.001;
    double JpX11 = (JX11delta - JX11) / (0.001 * X11);

    double X22 = 6 * 2 * 2 * APhi * sqrt(ISt);
    double JX22 = X22 / (4 + 4.581 * pow(X22, -0.7237) * exp(-0.012 * pow(X22, 0.528)));
    X22 = 1.001 * X22;
    double JX22delta = X22 / (4 + 4.581 * pow(X22, -0.7237) * exp(-0.012 * pow(X22, 0.528)));
    X22 = X22 / 1.001;
    double JpX22 = (JX22delta - JX22) / (0.001 * X22);

    double ETh = (0.5 / ISt) * (JX12 - 0.5 * (JX11 + JX22));
    double EThp = (0.25 / (ISt * ISt)) * (X12 * JpX12 - 0.5 * (X11 * JpX11 + X22 * JpX22)) - ETh / ISt;
    double Phip = EThp;

    // 计算 f_gamma（渗透系数部分）
    double f_gamma = -APhi * (sqrt(ISt) / (1 + 1.2 * sqrt(ISt)) + (2 / 1.2) * log(1 + 1.2 * sqrt(ISt)));
    for (c = 0; c < NumCat; c++) {
        for (a = 0; a < NumAn; a++) {
            double Bpca = b1[c][a] * gpX14 / ISt + b2[c][a] * gpX12 / ISt;
            // Holmes and Dai 特殊处理
            if (ChCat[c] == 1) {
                X20 = 2 * sqrt(ISt);
                if (c == 1 && a == 5) X20 = 1.4 * sqrt(ISt);  // Na-SO4
                if (c == 2 && a == 5) X20 = 1.4 * sqrt(ISt);  // K-SO4
                gpX20 = -2 * (1 - (1 + X20 + 0.5 * X20 * X20) * exp(-X20)) / (X20 * X20);
                Bpca = b1[c][a] * gpX20 / ISt + b2[c][a] * gpX12 / ISt;
            }
            if (ChCat[c] == 2 && ChAn[a] == -1) {
                X20 = (2 - 0.00181 * (tk - 298.15)) * sqrt(ISt);
                gpX20 = -2 * (1 - (1 + X20 + 0.5 * X20 * X20) * exp(-X20)) / (X20 * X20);
                Bpca = b1[c][a] * gpX20 / ISt + b2[c][a] * gpX12 / ISt;
            }
            // 添加 by GD 20191021：针对 Ca-SO4
            if (c == 4 && a == 5) {
                Bpca = b1[c][a] * gpX14 / ISt + b2[c][a] * gpXCaSO4 / ISt;
            }
            f_gamma += mc[c] * ma[a] * Bpca;
        }
    }

    for (c = 0; c < NumCat - 1; c++) {
        for (int cp = c + 1; cp < NumCat; cp++) {
            if (ChCat[c] != ChCat[cp]) {
                f_gamma += mc[c] * mc[cp] * Phip;
            }
        }
    }

    for (a = 0; a < NumAn - 1; a++) {
        for (int ap = a + 1; ap < NumAn; ap++) {
            if (ChAn[a] != ChAn[ap]) {
                f_gamma += ma[a] * ma[ap] * Phip;
            }
        }
    }

    // 阳离子活性系数循环
    for (m = 0; m < NumCat; m++) {
        double term1 = pow(ChCat[m], 2) * f_gamma;  // D-H 项

        double term2 = 0;
        for (a = 0; a < NumAn; a++) {
            double BMa = b0[m][a] + b1[m][a] * gX14 + b2[m][a] * gX12;
            // Holmes and Dai 特殊处理（注意：此处 c 应为 m）
            if (ChCat[m] == 1) {
                X20 = 2 * sqrt(ISt);
                if (c == 1 && a == 5) X20 = 1.4 * sqrt(ISt);  // Na-SO4
                if (c == 2 && a == 5) X20 = 1.4 * sqrt(ISt);  // K-SO4
                gX20 = 2 * (1 - (1 + X20) * exp(-X20)) / (X20 * X20);
                BMa = b0[m][a] + b1[m][a] * gX20 + b2[m][a] * gX12;
            }
            if (ChCat[m] == 2 && ChAn[a] == -1) {
                X20 = (2 - 0.00181 * (tk - 298.15)) * sqrt(ISt);
                gX20 = 2 * (1 - (1 + X20) * exp(-X20)) / (X20 * X20);
                BMa = b0[m][a] + b1[m][a] * gX20 + b2[m][a] * gX12;
            }
            // 添加 by GD 20191021：针对 Ca-SO4
            if (m == 4 && a == 5) {
                BMa = b0[m][a] + b1[m][a] * gX14 + b2[m][a] * gXCaSO4;
            }
            double CMa = CPhi[m][a] / (2 * sqrt(fabs(ChCat[m] * ChAn[a])));
            term2 += ma[a] * (2 * BMa + MoleCharge * CMa);
        }

        double term3 = 0;
        for (c = 0; c < NumCat; c++) {
            double PhiMc = Tccp[m][c];
            if (ChCat[m] != ChCat[c]) PhiMc = Tccp[m][c] + ETh;
            double SumYMca = 0;
            for (a = 0; a < NumAn; a++) {
                SumYMca += ma[a] * Yccpa[m][c][a];
            }
            term3 += mc[c] * (2 * PhiMc + SumYMca);
        }

        double term4 = 0;
        for (a = 0; a < NumAn - 1; a++) {
            for (int ap = a + 1; ap < NumAn; ap++) {
                term4 += ma[a] * ma[ap] * Yaapc[a][ap][m];
            }
        }

        double term5 = 0;
        for (c = 0; c < NumCat; c++) {
            for (a = 0; a < NumAn; a++) {
                double Cca = CPhi[c][a] / (2 * sqrt(fabs(ChCat[c] * ChAn[a])));
                term5 += mc[c] * ma[a] * Cca;
            }
        }
        term5 = fabs(ChCat[m]) * term5;

        double term6 = 0;
        for (int n = 0; n < NumNeut; n++) {
            term6 += mn[n] * 2 * Lnc[n][m];
        }

        double term7 = 0;
        for (int n = 0; n < NumNeut; n++) {
            for (a = 0; a < NumAn; a++) {
                // pH = pH;  // VB调试语句，忽略
                term7 += 6 * mn[n] * ma[a] * zeta[n][m][a];
            }
        }

        double termsum = term1 + term2 + term3 + term4 + term5 + term6 + term7;
        double ActCoefM = exp(termsum);
        gCat[m] = ActCoefM;
    }

    // 特殊设置：Pb, NH4, Ra 的活性系数
    gCat[iPb] = gCat[iZn];
    gCat[iNH4] = gCat[iK];
    gCat[iRa] = gCat[iCa];

    // 阴离子活性系数循环
    for (int iPz = 0; iPz < NumAn; iPz++) {
        double term1 = pow(ChAn[iPz], 2) * f_gamma;  // D-H 项

        double term2 = 0;
        for (c = 0; c < NumCat; c++) {
            double BcX = b0[c][iPz] + b1[c][iPz] * gX14 + b2[c][iPz] * gX12;
            // Holmes and Dai 特殊处理
            if (ChCat[c] == 1) {
                X20 = 2 * sqrt(ISt);
                if (c == 1 && iPz == 5) X20 = 1.4 * sqrt(ISt);  // Na-SO4
                if (c == 2 && iPz == 5) X20 = 1.4 * sqrt(ISt);  // K-SO4
                gX20 = 2 * (1 - (1 + X20) * exp(-X20)) / (X20 * X20);
                BcX = b0[c][iPz] + b1[c][iPz] * gX20 + b2[c][iPz] * gX12;
            }
            if (ChCat[c] == 2 && ChAn[iPz] == -1) {
                X20 = (2 - 0.00181 * (tk - 298.15)) * sqrt(ISt);
                gX20 = 2 * (1 - (1 + X20) * exp(-X20)) / (X20 * X20);
                BcX = b0[c][iPz] + b1[c][iPz] * gX20 + b2[c][iPz] * gX12;
            }
            // 添加 by GD 20191021：针对 Ca-SO4
            if (c == 4 && iPz == 5) {
                BcX = b0[c][iPz] + b1[c][iPz] * gX14 + b2[c][iPz] * gXCaSO4;
            }
            double CcX = CPhi[c][iPz] / (2 * sqrt(fabs(ChAn[iPz] * ChCat[c])));
            term2 += mc[c] * (2 * BcX + MoleCharge * CcX);
        }

        double term3 = 0;
        for (a = 0; a < NumAn; a++) {
            double PhiXa = Taap[iPz][a];
            if (ChAn[iPz] != ChAn[a]) PhiXa = Taap[iPz][a] + ETh;
            double SumYXac = 0;
            for (c = 0; c < NumCat; c++) {
                SumYXac += mc[c] * Yaapc[iPz][a][c];
            }
            term3 += ma[a] * (2 * PhiXa + SumYXac);
        }

        double term4 = 0;
        for (c = 0; c < NumCat - 1; c++) {
            for (int cp = c + 1; cp < NumCat; cp++) {
                term4 += mc[c] * mc[cp] * Yccpa[c][cp][iPz];
            }
        }

        double term5 = 0;
        for (c = 0; c < NumCat; c++) {
            for (a = 0; a < NumAn; a++) {
                double Cca = CPhi[c][a] / (2 * sqrt(fabs(ChCat[c] * ChAn[a])));
                term5 += mc[c] * ma[a] * Cca;
            }
        }
        term5 = fabs(ChAn[iPz]) * term5;

        double term6 = 0;
        for (int n = 0; n < NumNeut; n++) {
            term6 += mn[n] * 2 * Lna[n][iPz];
        }

        double term7 = 0;
        for (int n = 0; n < NumNeut; n++) {
            for (c = 0; c < NumCat; c++) {
                // pH = pH;  // VB调试语句，忽略
                term7 += 6 * mn[n] * mc[c] * zeta[n][c][iPz];
            }
        }

        double termsum = term1 + term2 + term3 + term4 + term5 + term6 + term7;
        double ActCoefX = exp(termsum);
        gAn[iPz] = ActCoefX;
    }
    gAn[iSion] = gAn[iCO3];  // 设置硫化物活性系数等于碳酸盐

    // 中性化合物活性系数
    for (int n = 0; n < NumNeut; n++) {
        double termn = 0;
        for (c = 0; c < NumCat; c++) {
            termn += mc[c] * 2 * Lnc[n][c];
        }
        for (a = 0; a < NumAn; a++) {
            termn += ma[a] * 2 * Lna[n][a];
        }
        //if (n == 2)
        //    pH = pH;

        //Lnn(15, 15)
        termn += 2 * mn[n] * Lnn[n][n];

        for (c = 0; c < NumCat; c++) {
            for (a = 0; a < NumAn; a++) {
                if (c == 1 && a == 1) pH = pH;
                // pH = pH;  // VB调试语句，忽略
                termn += mc[c] * ma[a] * zeta[n][c][a];
            }
        }

        double ActCoefn = exp(termn);
        gNeut[n] = ActCoefn;
    }
    // 特殊设置：CH4, HAc, NH3, FeS 的活性系数
    gNeut[iCH4aq] = gNeut[iCO2aq];
    gNeut[iHAcaq] = gNeut[iCO2aq];
    gNeut[iNH3] = gNeut[iCO2aq];
    gNeut[iFeSaq] = gNeut[iH2Saq];
    gNeut[iH4SiO4aq] = pow(10, (0.00978 * pow(10, 280 / tk) * ISt));  // Solmineq 88 page 46

    // 水的渗透系数 PhiH2O 和活性 aH2O
    double term1 = -APhi * pow(ISt, 1.5) / (1 + 1.2 * sqrt(ISt));
    double term2 = 0;
    for (c = 0; c < NumCat; c++) {
        for (a = 0; a < NumAn; a++) {
            double BPhica = b0[c][a] + b1[c][a] * exp(-1.4 * sqrt(ISt)) + b2[c][a] * exp(-12 * sqrt(ISt));
            // 修正 by GD 20191021
            if (ChCat[c] == 1) {
                X20 = 2 * sqrt(ISt);
                if (c == 1 && a == 5) X20 = 1.4 * sqrt(ISt);  // Na-SO4
                if (c == 2 && a == 5) X20 = 1.4 * sqrt(ISt);  // K-SO4
                BPhica = b0[c][a] + b1[c][a] * exp(-X20) + b2[c][a] * exp(-12 * sqrt(ISt));
            }
            if (ChCat[c] == 2 && ChAn[a] == -1) {
                X20 = (2 - 0.00181 * (tk - 298.15)) * sqrt(ISt);
                BPhica = b0[c][a] + b1[c][a] * exp(-X20) + b2[c][a] * exp(-12 * sqrt(ISt));
            }
            // 添加 by GD 20191021：针对 Ca-SO4
            if (c == 4 && a == 5) {
                BPhica = b0[c][a] + b1[c][a] * exp(-1.4 * sqrt(ISt)) + b2[c][a] * exp(-xCaSO4);
            }
            double Cca = CPhi[c][a] / (2 * sqrt(fabs(ChCat[c] * ChAn[a])));
            term2 += mc[c] * ma[a] * (BPhica + MoleCharge * Cca);
        }
    }

    double term3 = 0;
    for (c = 0; c < NumCat - 1; c++) {
        for (int cp = c + 1; cp < NumCat; cp++) {
            double PhiPhiccp = Tccp[c][cp];
            if (ChCat[c] != ChCat[cp]) PhiPhiccp = Tccp[c][cp] + ETh + ISt * EThp;
            double Sumccpa = 0;
            for (int a = 0; a < NumAn; a++) {
                Sumccpa += ma[a] * Yccpa[c][cp][a];
            }
            term3 += mc[c] * mc[cp] * (PhiPhiccp + Sumccpa);
        }
    }

    double term4 = 0;
    for (a = 0; a < NumAn - 1; a++) {
        for (int ap = a + 1; ap < NumAn; ap++) {
            double PhiPhiaap = Taap[a][ap];
            if (ChAn[a] != ChAn[ap]) PhiPhiaap = Taap[a][ap] + ETh + ISt * EThp;
            double Sumaapc = 0;
            for (c = 0; c < NumCat; c++) {
                Sumaapc += mc[c] * Yaapc[a][ap][c];
            }
            term4 += ma[a] * ma[ap] * (PhiPhiaap + Sumaapc);
        }
    }

    double term5 = 0;
    for (int n = 0; n < NumNeut; n++) {
        for (c = 0; c < NumCat; c++) {
            term5 += mn[n] * mc[c] * Lnc[n][c];
        }
    }

    double term6 = 0;
    for (int n = 0; n < NumNeut; n++) {
        for (a = 0; a < NumAn; a++) {
            term6 += mn[n] * ma[a] * Lna[n][a];
        }
    }

    double term7 = 0;
    for (int n = 0; n < NumNeut; n++) {
        for (c = 0; c < NumCat; c++) {
            for (a = 0; a < NumAn; a++) {
                // pH = pH;  // VB调试语句，忽略
                term7 += mn[n] * mc[c] * ma[a] * zeta[n][c][a];
            }
        }
    }

    double termsum = term1 + term2 + term3 + term4 + term5 + term6 + term7;
    double PhiH2O = 2 * termsum / mtotal + 1;  // 水的渗透系数
    aH2O = 1;
    if (fabs(18 * mtotal * PhiH2O / 1000) < 600) {
        aH2O = exp(-18 * mtotal * PhiH2O / 1000);  // 水的活性
    }

    // Bdot 参数用于 Zn 和 Pb 络合物（基于 Solmineq）
    double B_dot = 0.03804695 + 0.0001031039 * tc + 0.0000007119498 * pow(tc, 2) - 0.00000001968215 * pow(tc, 3) + 1.276773E-10 * pow(tc, 4) - 2.71893E-13 * pow(tc, 5);
    double B_gamma = (50.29158649 * sqrt(0.001 * dens)) / sqrt(Dielec * tk);
    double Lambda_gamma = -log10(1 + 0.0180153 * mtotal);

    // a0 参数设置
    double a0[20] = { 0 };
    a0[iClDot] = 3; a0[iZnDot] = 6; a0[iPbDot] = 4.5; a0[iHSDot] = 3; a0[iZnCl] = 4; a0[iZnCl2] = 0; a0[iZnCl3] = 4; a0[iZnCl4] = 5;
    a0[iPbCl] = 4; a0[iPbCl2] = 0; a0[iPbCl3] = 4; a0[iPbCl4] = 5; a0[iZnHS2] = 0; a0[iZnHS3] = 4; a0[iPbHS2] = 0; a0[iPbHS3] = 4;
    a0[iZnOH] = 4; a0[iZnOH2] = 0; a0[iOHDot] = 3.5;

    /* 优化：提前计算，防止重复调用log和sqrt */
    double sqrtI = sqrt(ISt);
    double LOG10 = log(10.0);

    gDot[iClDot] = pow(10.0,
        -3.0 / LOG10 * APhi * 1.0 * sqrtI / (1.0 + a0[iClDot] * B_gamma * sqrtI) + Lambda_gamma + B_dot * ISt);
    gDot[iZnDot] = pow(10.0,
        -3.0 / LOG10 * APhi * 4.0 * sqrtI / (1.0 + a0[iZnDot] * B_gamma * sqrtI) + Lambda_gamma + B_dot * ISt);
    gDot[iPbDot] = pow(10.0,
        -3.0 / LOG10 * APhi * 4.0 * sqrtI / (1.0 + a0[iPbDot] * B_gamma * sqrtI) + Lambda_gamma + B_dot * ISt);
    gDot[iHSDot] = pow(10.0,
        -3.0 / LOG10 * APhi * 4.0 * sqrtI / (1.0 + a0[iHSDot] * B_gamma * sqrtI) + Lambda_gamma + B_dot * ISt);
    gDot[iOHDot] = pow(10.0,
        -3.0 / LOG10 * APhi * 4.0 * sqrtI / (1.0 + a0[iOHDot] * B_gamma * sqrtI) + Lambda_gamma + B_dot * ISt);

    /* Zn-Cl complexes */
    gDot[iZnCl] = pow(10.0,
        -3.0 / LOG10 * APhi * 1.0 * sqrtI / (1.0 + a0[iZnCl] * B_gamma * sqrtI) + Lambda_gamma + B_dot * ISt);
    gDot[iZnCl2] = 1.0;
    gDot[iZnCl3] = pow(10.0,
        -3.0 / LOG10 * APhi * 1.0 * sqrtI / (1.0 + a0[iZnCl3] * B_gamma * sqrtI) + Lambda_gamma + B_dot * ISt);

    gDot[iZnCl4] = pow(10.0,
        -3.0 / LOG10 * APhi * 4.0 * sqrtI / (1.0 + a0[iZnCl4] * B_gamma * sqrtI) + Lambda_gamma + B_dot * ISt);

    /* Zn–HS complexes */
    gDot[iZnHS2] = 1.0;
    gDot[iZnHS3] = pow(10.0,
        -3.0 / LOG10 * APhi * 1.0 * sqrtI / (1.0 + a0[iZnHS3] * B_gamma * sqrtI) + Lambda_gamma + B_dot * ISt);

    /* Pb–Cl complexes */
    gDot[iPbCl] = pow(10.0,
        -3.0 / LOG10 * APhi * 1.0 * sqrtI / (1.0 + a0[iPbCl] * B_gamma * sqrtI) + Lambda_gamma + B_dot * ISt);
    gDot[iPbCl2] = 1.0;
    gDot[iPbCl3] = pow(10.0,
        -3.0 / LOG10 * APhi * 1.0 * sqrtI / (1.0 + a0[iPbCl3] * B_gamma * sqrtI) + Lambda_gamma + B_dot * ISt);
    gDot[iPbCl4] = pow(10.0,
        -3.0 / LOG10 * APhi * 4.0 * sqrtI / (1.0 + a0[iPbCl4] * B_gamma * sqrtI) + Lambda_gamma + B_dot * ISt);

    /* Pb–HS complexes */
    gDot[iPbHS2] = 1.0;
    gDot[iPbHS3] = pow(10.0,
        -3.0 / LOG10 * APhi * 1.0 * sqrtI / (1.0 + a0[iPbHS3] * B_gamma * sqrtI) + Lambda_gamma + B_dot * ISt);

    /* Zn-OH complexes */
    gDot[iZnOH] = pow(10.0,
        -3.0 / LOG10 * APhi * 1.0 * sqrtI / (1.0 + a0[iZnOH] * B_gamma * sqrtI) + Lambda_gamma + B_dot * ISt);
    gDot[iZnOH2] = 1.0;
}


void CubicRoots(double coef1, double coef2, double coef3, double coef4, double& root1, double& root2, double& root3)
{
    double QCubic, Rcubic, Xcubic, Theta;
    double pi = 3.14159265358979;
    root1 = 0.0;
    root2 = 0.0;
    root3 = 0.0;

    QCubic = (pow(coef2, 2.0) - 3.0 * coef3) / 9.0;
    Rcubic = (2.0 * pow(coef2, 3.0) - 9.0 * coef2 * coef3 + 27.0 * coef4) / 54.0;

    if (pow(Rcubic, 2.0) - pow(QCubic, 3.0) > 0.0) {
        /* 只有一个实根 */
        double A = pow(fabs(Rcubic) + sqrt(pow(Rcubic, 2.0) - pow(QCubic, 3.0)), 1.0 / 3.0);
        double B = QCubic / A;
        double sigh = (Rcubic > 0) ? 1.0 : ((Rcubic < 0) ? -1.0 : 0.0);
        root1 = -sigh * (A + B) - coef2 / 3.0;
    }
    else {
        /* 三个实根 */
        Xcubic = Rcubic / sqrt(pow(QCubic, 3.0));
        /* 使用反三角恒等式 arccos(x) = atan(-x/sqrt(-x^2+1)) + pi/2 */
        Theta = atan(-Xcubic / sqrt(1.0 - Xcubic * Xcubic)) + pi / 2.0;
        root1 = -2.0 * sqrt(QCubic) * cos(Theta / 3.0) - coef2 / 3.0;
        root2 = -2.0 * sqrt(QCubic) * cos((Theta + 2.0 * pi) / 3.0) - coef2 / 3.0;
        root3 = -2.0 * sqrt(QCubic) * cos((Theta + 4.0 * pi) / 3.0) - coef2 / 3.0;
    }
}

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
void QuadraticRoots(double coef1, double coef2, double coef3, double& root1, double& root2)
{
    double qroot;
    double discriminant = pow(coef2, 2.0) - 4.0 * coef1 * coef3;
    if (discriminant < 0.0) {
        /* 无实根 */
        root1 = NAN;
        root2 = NAN;
        return;
    }
    qroot = -0.5 * (coef2 - sqrt(discriminant));
    root1 = coef3 / qroot;
    root2 = qroot / coef1;
}


void PengRobinson3() {
    int iNG, jNG;
    double aPrmix, bPRmix, AStarPR, BStarPr;
    double coef1, coef2, coef3, coef4;
    double root1, root2, root3, VolGasMolar;
    double Sum_aijPR[3];
    double term1, term2, term3, term4;

    // 压力单位转换
    PBar = Patm * 1.013254;

    // 如果没有气体，设置默认值避免除零错误
    if (simContext.yCO2 == 0.0 && simContext.yH2S == 0.0) {
        yCH4 = 1.0;
    }

    double yGas[20] = { 0 };
    // 设置气体摩尔分数
    yGas[iCH4g] = yCH4;
    yGas[iCO2g] = simContext.yCO2;
    yGas[iH2Sg] = simContext.yH2S;
    double TCr[20] = { 0 };
    // 设置临界温度 (K)
    TCr[iCH4g] = 190.4;
    TCr[iCO2g] = 304.1;
    TCr[iH2Sg] = 373.2;
    double Pc[20] = { 0 };
    // 设置临界压力 (bar)
    Pc[iCH4g] = 46.0;
    Pc[iCO2g] = 73.8;
    Pc[iH2Sg] = 89.4;
    double Omega[20] = { 0 };// 局部数组
    // 设置偏心因子
    Omega[iCH4g] = 0.011;
    Omega[iCO2g] = 0.239;
    Omega[iH2Sg] = 0.081;

    double kPr[20][20] = { 0 };// 局部数组

    // 初始化二元交互参数矩阵
    for (iNG = 0; iNG < 3; iNG++) {
        for (jNG = 0; jNG < 3; jNG++) {
            kPr[iNG][jNG] = 0.0;
        }
    }

    // 设置二元交互参数
    kPr[iCH4g][iCO2g] = 0.0919;
    kPr[iCO2g][iCH4g] = 0.0919;
    kPr[iH2Sg][iCO2g] = 0.0974;
    kPr[iCO2g][iH2Sg] = 0.0974;
    kPr[iH2Sg][iCH4g] = 0.084;
    kPr[iCH4g][iH2Sg] = 0.084;
    double F_Omega[20] = { 0 };
    double aPR[20] = { 0 };
    double bPR[20] = { 0 };
    // 计算Peng-Robinson参数
    for (iNG = 0; iNG < 3; iNG++) {
        F_Omega[iNG] = 0.37464 + 1.54226 * Omega[iNG] - 0.26992 * pow(Omega[iNG], 2);

        double temp = 1.0 + F_Omega[iNG] * (1.0 - sqrt(TK / TCr[iNG]));
        aPR[iNG] = (0.45724 * pow(RBar, 2) * pow(TCr[iNG], 2) / Pc[iNG]) * pow(temp, 2);

        bPR[iNG] = 0.0778 * RBar * TCr[iNG] / Pc[iNG];
    }

    // 混合规则
    aPrmix = 0.0;
    bPRmix = 0.0;

    for (iNG = 0; iNG < 3; iNG++) {
        bPRmix += yGas[iNG] * bPR[iNG];

        for (jNG = 0; jNG < 3; jNG++) {
            aPrmix += yGas[iNG] * yGas[jNG] * sqrt(aPR[iNG] * aPR[jNG]) * (1.0 - kPr[iNG][jNG]);
        }
    }

    // 计算无量纲参数
    AStarPR = aPrmix * PBar / pow(RBar * TK, 2);
    BStarPr = bPRmix * PBar / (RBar * TK);

    // 设置立方状态方程的系数
    coef1 = 1.0;
    coef2 = -(1.0 - BStarPr);
    coef3 = AStarPR - 2.0 * BStarPr - 3.0 * pow(BStarPr, 2);
    coef4 = -AStarPR * BStarPr + pow(BStarPr, 2) + pow(BStarPr, 3);

    // 求解立方方程
    CubicRoots(coef1, coef2, coef3, coef4, root1, root2, root3);

    // 选择最大的根作为气体压缩因子
    Znew = root1;
    if (root2 > Znew) Znew = root2;
    if (root3 > Znew) Znew = root3;

    // 计算摩尔体积
    VolGasMolar = Znew * RBar * TK / PBar;

    // 计算Sum_aijPR
    for (iNG = 0; iNG < 3; iNG++) {
        Sum_aijPR[iNG] = 0.0;
        for (jNG = 0; jNG < 3; jNG++) {
            Sum_aijPR[iNG] += 2.0 * yGas[jNG] * sqrt(aPR[iNG] * aPR[jNG]) * (1.0 - kPr[iNG][jNG]);
        }
    }

    // 计算逸度系数
    for (iNG = 0; iNG < 3; iNG++) {
        term1 = (aPrmix / (2.8284 * bPRmix * RBar * TK)) *
            (Sum_aijPR[iNG] / aPrmix - bPR[iNG] / bPRmix);

        term2 = log((VolGasMolar + 2.4142 * bPRmix) / (VolGasMolar - 0.4142 * bPRmix));

        term3 = (bPR[iNG] / bPRmix) * (Znew - 1.0);

        term4 = log(Znew - BStarPr);

        gGas[iNG] = exp(term3 - term4 - term1 * term2);
    }
}

/**
 * @brief 计算STP条件下的pH、PCO2和PH2S
 *
 * 此函数用于在标准温度压力条件下计算水溶液的pH值、CO2分压和H2S分压。
 * 根据不同的输入条件组合（pH测量值、碱度、气体分压等），采用二分法求解pH值，
 * 并计算相关化学物种的浓度分布。
 *
 * @param use_pH pH使用标志: 0-计算pH, 1-使用pH计算PCO2, 2-使用pH计算碱度, 3-使用TCO2计算pH
 * @param UseH2Sgas H2S气体使用标志: 0-使用TH2Saq, 1-使用yH2S
 * @param useEOS 状态方程使用标志
 * @param TK 温度 (K)
 * @param Ppsia 压力 (psia)
 * @param simContext.yCO2 CO2摩尔分数
 * @param simContext.yH2S H2S摩尔分数
 * @param Alk 碱度 (mol/kg)
 * @param TAc 总乙酸浓度 (mol/kg)
 * @param TH2Saq 总溶解H2S浓度 (mol/kg)
 * @param TFe 总铁浓度 (mol/kg)
 * @param TCO2 总无机碳浓度 (mol/kg)
 * @param TNH4 总铵浓度 (mol/kg)
 * @param TH3BO3 总硼酸浓度 (mol/kg)
 * @param TH4SiO4 总硅酸浓度 (mol/kg)
 * @param pH [输入/输出] pH值
 * @param pHMeterReading [输出] pH计读数
 * @param simContext.errmsg [输出] 错误信息数组
 */
double fMeSSpeciation(int im, int igas);  //C5需要用的函数，声明在这里
void C5_CalcpHPCO2PH2SSTP(int use_pH, int UseH2Sgas, int useEOS) {

    // 局部变量声明
    double tHCO3, tCO3;
    int k;

    // 情况1: 使用P-CO2和碱度计算pH (use_pH = 0, UseH2Sgas = 0)
    if (use_pH == 0 && UseH2Sgas == 0 && useEOS == 0) {
        pHHigh = 14.0;
        pHLow = 0.0;

        for (k = 0; k < 30; k++) {
            pH = (pHHigh + pHLow) / 2.0;
            aH = pow(10.0, -(pH));

            // 计算H+和OH-浓度
            H = aH / (gCat[iH] * gNCat[iH]);
            OH = KH2O / (aH * gAn[iOH] * gNAn[iOH]);

            // 计算碳酸系统物种 gGas 小数点后第五位不一致
            simContext.CO2aq = KgwCO2 * Ppsia * (simContext.yCO2) * gGas[iCO2g] / (gNeut[iCO2aq] * gNNeut[iCO2aq]);
            simContext.HCO3 = (K1H2CO3 * aH2O) * simContext.CO2aq * gNeut[iCO2aq] * gNNeut[iCO2aq] /
                (aH * gAn[iHCO3] * gNAn[iHCO3]);
            simContext.CO3 = K2HCO3 * simContext.HCO3 * gAn[iHCO3] * gNAn[iHCO3] /
                (aH * gAn[iCO3] * gNAn[iCO3]);

            // 计算硫系统物种
            hydHS = aH * gAn[iHS] * gNAn[iHS] / (K1H2S * gNeut[iH2Saq] * gNNeut[iH2Saq]) +
                1.0 + (K2HS * gAn[iHS] * gNAn[iHS]) / (aH * gAn[iSion] * gNAn[iSion]);
            simContext.HS = TH2Saq / hydHS;

            // // 金属硫化物形态计算
            if (TH2Saq > 0) fMeSSpeciation(1, 2);  // 源代码是(2, 2)，因为第一个参数im涉及到索引，做一个偏移

            simContext.H2Saq = aH * simContext.HS * gAn[iHS] * gNAn[iHS] / (K1H2S * gNeut[iH2Saq] * gNNeut[iH2Saq]);
            S = K2HS * simContext.HS * gAn[iHS] * gNAn[iHS] / (aH * gAn[iSion] * gNAn[iSion]);
            simContext.yH2S = simContext.H2Saq * gNeut[iH2Saq] * gNNeut[iH2Saq] / (KgwH2S * Ppsia * gGas[iH2Sg]);

            // 计算其他弱酸物种
            hydAc = aH * gAn[iAc] * gNAn[iAc] / (KHAc * gNeut[iHAcaq] * gNNeut[iHAcaq]) + 1.0;
            AC = TAc / hydAc;
            HAcaq = TAc - AC;

            hydH2BO3 = aH * gAn[iH2BO3] * gNAn[iH2BO3] / (KH3BO3 * gNeut[iH3BO3] * gNNeut[iH3BO3]) + 1.0;
            H2BO3 = TH3BO3 / hydH2BO3;

            hydNH3 = aH * gNeut[iNH3] * gNNeut[iNH3] / (KNH4 * gCat[iNH4] * gNCat[iNH4]) + 1.0;
            NH3 = TNH4 / hydNH3;

            hydH2SiO4 = aH * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]) +
                aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]) + 1.0;

            H2SiO4 = TH4SiO4 / hydH2SiO4;
            H3SiO4 = H2SiO4 * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]);
            H4SiO4 = H2SiO4 * aH * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]);

            // 计算碱度残差
            faH = Alk - (simContext.HCO3 + 2.0 * simContext.CO3 + simContext.HS + 2.0 * S + AC + NH3 +
                H2BO3 + H3SiO4 + 2.0 * H2SiO4 + OH - H);

            // 二分法更新pH范围
            if (faH > 0) pHLow = pH;
            else pHHigh = pH;
        }

        pHMeterReading = pH - DpHj;
        mn[iFeSaq] = KstFeSaq * mc[iFe] * simContext.HS * gAn[iHS] * gNAn[iHS] *
            gCat[iFe] * gNCat[iFe] / (gNeut[iFeSaq] * gNNeut[iFeSaq] * aH);

        if (simContext.yH2S > 1.0) {
            simContext.errmsg[2] = 3;
            simContext.yH2S = 1.0;
        }
    }

    if (use_pH == 0 && UseH2Sgas == 1 && useEOS == 0) {
        pHHigh = 14.0;
        pHLow = 0.0;

        for (k = 0; k < 30; k++) {
            pH = (pHHigh + pHLow) / 2.0;
            aH = pow(10.0, -(pH));

            // 计算H+和OH-浓度
            H = aH / (gCat[iH] * gNCat[iH]);
            OH = KH2O / (aH * gAn[iOH] * gNAn[iOH]);

            // 计算碳酸系统物种
            simContext.CO2aq = KgwCO2 * Ppsia * (simContext.yCO2) * gGas[iCO2g] / (gNeut[iCO2aq] * gNNeut[iCO2aq]);
            simContext.HCO3 = (K1H2CO3 * aH2O) * simContext.CO2aq * gNeut[iCO2aq] * gNNeut[iCO2aq] /
                (aH * gAn[iHCO3] * gNAn[iHCO3]);
            simContext.CO3 = K2HCO3 * simContext.HCO3 * gAn[iHCO3] * gNAn[iHCO3] /
                (aH * gAn[iCO3] * gNAn[iCO3]);

            // 计算硫系统物种
            simContext.H2Saq = KgwH2S * Ppsia * (simContext.yH2S) * gGas[iH2Sg] / (gNeut[iH2Saq] * gNNeut[iH2Saq]);
            simContext.HS = K1H2S * simContext.H2Saq * gNeut[iH2Saq] * gNNeut[iH2Saq] /
                (aH * gAn[iHS] * gNAn[iHS]);
            S = K2HS * simContext.HS * gAn[iHS] * gNAn[iHS] /
                (aH * gAn[iSion] * gNAn[iSion]);

            // 计算其他弱酸物种
            hydAc = aH * gAn[iAc] * gNAn[iAc] / (KHAc * gNeut[iHAcaq] * gNNeut[iHAcaq]) + 1.0;
            AC = TAc / hydAc;
            HAcaq = TAc - AC;

            hydH2BO3 = aH * gAn[iH2BO3] * gNAn[iH2BO3] / (KH3BO3 * gNeut[iH3BO3] * gNNeut[iH3BO3]) + 1.0;
            H2BO3 = TH3BO3 / hydH2BO3;

            hydNH3 = aH * gNeut[iNH3] * gNNeut[iNH3] / (KNH4 * gCat[iNH4] * gNCat[iNH4]) + 1.0;
            NH3 = TNH4 / hydNH3;

            hydH2SiO4 = aH * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]) +
                aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]) + 1.0;

            H2SiO4 = TH4SiO4 / hydH2SiO4;
            H3SiO4 = H2SiO4 * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]);
            H4SiO4 = H2SiO4 * aH * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]);

            // 计算碱度残差
            faH = Alk - (simContext.HCO3 + 2.0 * simContext.CO3 + simContext.HS + 2.0 * S + AC + NH3 +
                H2BO3 + H3SiO4 + 2.0 * H2SiO4 + OH - H);

            // 二分法更新pH范围
            if (faH > 0) pHLow = pH;
            else pHHigh = pH;
        }

        pHMeterReading = pH - DpHj;
        mc[iFe] = TFe / (1.0 + KstFeSaq * simContext.HS * gAn[iHS] * gNAn[iHS] *
            gCat[iFe] * gNCat[iFe] / (gNeut[iFeSaq] * gNNeut[iFeSaq] * aH));
        mn[iFeSaq] = KstFeSaq * mc[iFe] * simContext.HS * gAn[iHS] * gNAn[iHS] *
            gCat[iFe] * gNCat[iFe] / (gNeut[iFeSaq] * gNNeut[iFeSaq] * aH);
    }
    // 情况2: 使用pH和碱度计算P-CO2 (use_pH = 1)
    if (use_pH == 1) {
        aH = pow(10.0, -(pH));

        // H2S系统计算
        if (UseH2Sgas == 0) {
            hydHS = aH * gAn[iHS] * gNAn[iHS] / (K1H2S * gNeut[iH2Saq] * gNNeut[iH2Saq]) +
                1.0 + (K2HS * gAn[iHS] * gNAn[iHS]) / (aH * gAn[iSion] * gNAn[iSion]);
            simContext.HS = TH2Saq / hydHS;

            if (TH2Saq > 0) {
                fMeSSpeciation(1, 2);
            }

            simContext.H2Saq = aH * simContext.HS * gAn[iHS] * gNAn[iHS] / (K1H2S * gNeut[iH2Saq] * gNNeut[iH2Saq]);
            S = K2HS * simContext.HS * gAn[iHS] * gNAn[iHS] / (aH * gAn[iSion] * gNAn[iSion]);
            simContext.yH2S = simContext.H2Saq * gNeut[iH2Saq] * gNNeut[iH2Saq] / (KgwH2S * Ppsia * gGas[iH2Sg]);

            if (simContext.yH2S > 1.0) {
                simContext.errmsg[2] = 3;
                simContext.yH2S = 1.0;
            }
        }
        else {
            simContext.H2Saq = KgwH2S * Ppsia * (simContext.yH2S) * gGas[iH2Sg] / gNeut[iH2Saq] / gNNeut[iH2Saq];
            simContext.HS = K1H2S * simContext.H2Saq * gNeut[iH2Saq] * gNNeut[iH2Saq] / (aH * gAn[iHS] * gNAn[iHS]);
            S = K2HS * simContext.HS * gAn[iHS] * gNAn[iHS] / (aH * gAn[iSion] * gNAn[iSion]);
            mc[iFe] = TFe / (1.0 + KstFeSaq * simContext.HS * gAn[iHS] * gNAn[iHS] *
                gCat[iFe] * gNCat[iFe] / (gNeut[iFeSaq] * gNNeut[iFeSaq] * aH));
            TH2Saq = simContext.H2Saq + simContext.HS + S + KstFeSaq * mc[iFe] * simContext.HS * gAn[iHS] * gNAn[iHS] *
                gCat[iFe] * gNCat[iFe] / (gNeut[iFeSaq] * gNNeut[iFeSaq] * aH);
        }

        if (TH2Saq == 0.0 && simContext.yH2S == 0.0) {
            simContext.H2Saq = 0.0; simContext.HS = 0.0; TH2Saq = 0.0; simContext.yH2S = 0.0;
        }

        // 计算其他物种
        H = aH / gCat[iH] / gNCat[iH];
        OH = KH2O / (aH * gAn[iOH] * gNAn[iOH]);
        hydAc = aH * gAn[iAc] * gNAn[iAc] / (KHAc * gNeut[iHAcaq] * gNNeut[iHAcaq]) + 1.0;
        AC = TAc / hydAc;
        HAcaq = TAc - AC;

        hydH2BO3 = aH * gAn[iH2BO3] * gNAn[iH2BO3] / (KH3BO3 * gNeut[iH3BO3] * gNNeut[iH3BO3]) + 1.0;
        H2BO3 = TH3BO3 / hydH2BO3;

        hydNH3 = aH * gNeut[iNH3] * gNNeut[iNH3] / (KNH4 * gCat[iNH4] * gNCat[iNH4]) + 1.0;
        NH3 = TNH4 / hydNH3;

        hydH2SiO4 = aH * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]) +
            aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]) + 1.0;

        H2SiO4 = TH4SiO4 / hydH2SiO4;
        H3SiO4 = H2SiO4 * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]);
        H4SiO4 = H2SiO4 * aH * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]);

        // 计算yCO2
        tHCO3 = (K1H2CO3 * aH2O) * KgwCO2 * Ppsia * gGas[iCO2g] /
            (aH * gAn[iHCO3] * gNAn[iHCO3]);
        tCO3 = (K1H2CO3 * aH2O) * K2HCO3 * KgwCO2 * Ppsia * gGas[iCO2g] /
            (aH * aH * gAn[iCO3] * gNAn[iCO3]);

        simContext.yCO2 = (Alk + H - AC - simContext.HS - OH - NH3 - H2BO3 - H3SiO4 - 2.0 * H2SiO4) /
            (tHCO3 + 2.0 * tCO3);

        // 计算碳酸物种浓度
        simContext.HCO3 = (K1H2CO3 * aH2O) * KgwCO2 * Ppsia * (simContext.yCO2) * gGas[iCO2g] /
            (aH * gAn[iHCO3] * gNAn[iHCO3]);
        simContext.CO3 = (K1H2CO3 * aH2O) * K2HCO3 * KgwCO2 * Ppsia * (simContext.yCO2) * gGas[iCO2g] /
            (aH * aH * gAn[iCO3] * gNAn[iCO3]);
        simContext.CO2aq = KgwCO2 * Ppsia * (simContext.yCO2) * gGas[iCO2g] / gNeut[iCO2aq] / gNNeut[iCO2aq];

        mn[iFeSaq] = KstFeSaq * mc[iFe] * simContext.HS * gAn[iHS] * gNAn[iHS] *
            gCat[iFe] * gNCat[iFe] / (gNeut[iFeSaq] * gNNeut[iFeSaq] * aH);

        // 错误检查
        if (simContext.yCO2 > 1.0) {
            simContext.errmsg[0] = 1;
            simContext.yCO2 = 1.0;
        }
        if (simContext.yCO2 < 0.0) {
            simContext.errmsg[1] = 2;
            simContext.yCO2 = 0.0;
            simContext.HCO3 = 0.0; simContext.CO3 = 0.0; simContext.CO2aq = 0.0;
        }
    }
    // 情况3: 使用pH和yCO2计算Alk (use_pH = 2)
    if (use_pH == 2 && useEOS == 0) {
        aH = pow(10.0, -(pH));

        // H2S系统计算
        if (UseH2Sgas == 0) {
            hydHS = aH * gAn[iHS] * gNAn[iHS] / (K1H2S * gNeut[iH2Saq] * gNNeut[iH2Saq]) +
                1.0 + (K2HS * gAn[iHS] * gNAn[iHS]) / (aH * gAn[iSion] * gNAn[iSion]);
            simContext.HS = TH2Saq / hydHS;

            if (TH2Saq > 0) fMeSSpeciation(1, 2);

            simContext.H2Saq = aH * simContext.HS * gAn[iHS] * gNAn[iHS] / (K1H2S * gNeut[iH2Saq] * gNNeut[iH2Saq]);
            S = K2HS * simContext.HS * gAn[iHS] * gNAn[iHS] / (aH * gAn[iSion] * gNAn[iSion]);
            simContext.yH2S = simContext.H2Saq * gNeut[iH2Saq] * gNNeut[iH2Saq] / (KgwH2S * Ppsia * gGas[iH2Sg]);

            if (simContext.yH2S > 1.0) {
                simContext.errmsg[2] = 3;
                simContext.yH2S = 1.0;
            }
        }
        else {
            simContext.H2Saq = KgwH2S * Ppsia * (simContext.yH2S) * gGas[iH2Sg] / gNeut[iH2Saq] / gNNeut[iH2Saq];
            simContext.HS = K1H2S * simContext.H2Saq * gNeut[iH2Saq] * gNNeut[iH2Saq] / (aH * gAn[iHS] * gNAn[iHS]);
            S = K2HS * simContext.HS * gAn[iHS] * gNAn[iHS] / (aH * gAn[iSion] * gNAn[iSion]);
            mc[iFe] = TFe / (1.0 + KstFeSaq * simContext.HS * gAn[iHS] * gNAn[iHS] *
                gCat[iFe] * gNCat[iFe] / (gNeut[iFeSaq] * gNNeut[iFeSaq] * aH));
            TH2Saq = simContext.H2Saq + simContext.HS + S + KstFeSaq * mc[iFe] * simContext.HS * gAn[iHS] * gNAn[iHS] *
                gCat[iFe] * gNCat[iFe] / (gNeut[iFeSaq] * gNNeut[iFeSaq] * aH);
        }

        if (TH2Saq == 0.0 && simContext.yH2S == 0.0) {
            simContext.H2Saq = 0.0; simContext.HS = 0.0; TH2Saq = 0.0; simContext.yH2S = 0.0;
            simContext.TH2SaqMix[simContext.kk] = 0.0; simContext.yH2SMix[simContext.kk] = 0.0;  // Assuming i is defined elsewhere
        }

        // 计算其他弱酸物种
        hydH2SiO4 = aH * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]) +
            aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]) + 1.0;

        H2SiO4 = TH4SiO4 / hydH2SiO4;
        H3SiO4 = H2SiO4 * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]);
        H4SiO4 = H2SiO4 * aH * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]);

        H = aH / gCat[iH] / gNCat[iH];
        OH = KH2O / (aH * gAn[iOH] * gNAn[iOH]);
        hydAc = aH * gAn[iAc] * gNAn[iAc] / (KHAc * gNeut[iHAcaq] * gNNeut[iHAcaq]) + 1.0;
        AC = TAc / hydAc;
        HAcaq = TAc - AC;

        hydH2BO3 = aH * gAn[iH2BO3] * gNAn[iH2BO3] / (KH3BO3 * gNeut[iH3BO3] * gNNeut[iH3BO3]) + 1.0;
        H2BO3 = TH3BO3 / hydH2BO3;

        hydNH3 = aH * gNeut[iNH3] * gNNeut[iNH3] / (KNH4 * gCat[iNH4] * gNCat[iNH4]) + 1.0;
        NH3 = TNH4 / hydNH3;

        // 计算碳酸物种浓度
        simContext.HCO3 = (K1H2CO3 * aH2O) * KgwCO2 * Ppsia * (simContext.yCO2) * gGas[iCO2g] /
            (aH * gAn[iHCO3] * gNAn[iHCO3]);
        simContext.CO3 = (K1H2CO3 * aH2O) * K2HCO3 * KgwCO2 * Ppsia * (simContext.yCO2) * gGas[iCO2g] /
            (aH * aH * gAn[iCO3] * gNAn[iCO3]);

        // 计算Alk
        Alk = simContext.HCO3 + 2.0 * simContext.CO3 + AC + NH3 + H2BO3 + simContext.HS + 2.0 * S + H3SiO4 + 2.0 * H2SiO4 + OH - H;

        simContext.CO2aq = KgwCO2 * Ppsia * (simContext.yCO2) * gGas[iCO2g] / gNeut[iCO2aq] / gNNeut[iCO2aq];

        mn[iFeSaq] = KstFeSaq * mc[iFe] * simContext.HS * gAn[iHS] * gNAn[iHS] *
            gCat[iFe] * gNCat[iFe] / (gNeut[iFeSaq] * gNNeut[iFeSaq] * aH);
    }
    // 情况4: 使用TCO2和碱度计算pH (use_pH = 3, UseH2Sgas = 0)
    if (use_pH == 3 && UseH2Sgas == 0 && useEOS == 0) {
        pHHigh = 14.0;
        pHLow = 0.0;

        for (k = 0; k < 30; k++) {
            pH = (pHHigh + pHLow) / 2.0;
            aH = pow(10.0, -(pH));

            // 计算H+和OH-浓度
            H = aH / (gCat[iH] * gNCat[iH]);
            OH = KH2O / (aH * gAn[iOH] * gNAn[iOH]);

            // 计算碳酸系统物种 (从TCO2)
            simContext.HCO3 = TCO2 / (aH / (K1H2CO3 * aH2O) / gNeut[iCO2aq] / gNNeut[iCO2aq] * gAn[iHCO3] * gNAn[iHCO3] +
                1.0 + K2HCO3 / aH / gAn[iCO3] / gNAn[iCO3] * gAn[iHCO3] * gNAn[iHCO3]);
            simContext.CO2aq = aH * simContext.HCO3 * gAn[iHCO3] * gNAn[iHCO3] / ((K1H2CO3 * aH2O) * gNeut[iCO2aq] * gNNeut[iCO2aq]);
            simContext.CO3 = K2HCO3 * simContext.HCO3 * gAn[iHCO3] * gNAn[iHCO3] / (aH * gAn[iCO3] * gNAn[iCO3]);
            simContext.yCO2 = simContext.CO2aq * gNeut[iCO2aq] * gNNeut[iCO2aq] / (KgwCO2 * Ppsia * gGas[iCO2g]);

            // 计算硫系统物种
            hydHS = aH * gAn[iHS] * gNAn[iHS] / (K1H2S * gNeut[iH2Saq] * gNNeut[iH2Saq]) +
                1.0 + (K2HS * gAn[iHS] * gNAn[iHS]) / (aH * gAn[iSion] * gNAn[iSion]);
            simContext.HS = TH2Saq / hydHS;

            // // 金属硫化物形态计算
            if (TH2Saq > 0) fMeSSpeciation(1, 2);

            simContext.H2Saq = aH * simContext.HS * gAn[iHS] * gNAn[iHS] / (K1H2S * gNeut[iH2Saq] * gNNeut[iH2Saq]);
            S = K2HS * simContext.HS * gAn[iHS] * gNAn[iHS] / (aH * gAn[iSion] * gNAn[iSion]);
            simContext.yH2S = simContext.H2Saq * gNeut[iH2Saq] * gNNeut[iH2Saq] / (KgwH2S * Ppsia * gGas[iH2Sg]);

            // 计算其他弱酸物种
            hydAc = aH * gAn[iAc] * gNAn[iAc] / (KHAc * gNeut[iHAcaq] * gNNeut[iHAcaq]) + 1.0;
            AC = TAc / hydAc;
            HAcaq = TAc - AC;

            hydH2BO3 = aH * gAn[iH2BO3] * gNAn[iH2BO3] / (KH3BO3 * gNeut[iH3BO3] * gNNeut[iH3BO3]) + 1.0;
            H2BO3 = TH3BO3 / hydH2BO3;

            hydNH3 = aH * gNeut[iNH3] * gNNeut[iNH3] / (KNH4 * gCat[iNH4] * gNCat[iNH4]) + 1.0;
            NH3 = TNH4 / hydNH3;

            hydH2SiO4 = aH * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]) +
                aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]) + 1.0;

            H2SiO4 = TH4SiO4 / hydH2SiO4;
            H3SiO4 = H2SiO4 * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]);
            H4SiO4 = H2SiO4 * aH * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]);

            // 计算碱度残差
            faH = Alk - (simContext.HCO3 + 2.0 * simContext.CO3 + simContext.HS + 2.0 * S + AC + NH3 +
                H2BO3 + H3SiO4 + 2.0 * H2SiO4 + OH - H);

            // 二分法更新pH范围
            if (faH > 0) pHLow = pH;
            else pHHigh = pH;
        }

        pHMeterReading = pH - DpHj;
        mn[iFeSaq] = KstFeSaq * mc[iFe] * simContext.HS * gAn[iHS] * gNAn[iHS] *
            gCat[iFe] * gNCat[iFe] / (gNeut[iFeSaq] * gNNeut[iFeSaq] * aH);

        if (simContext.yH2S > 1.0) {
            simContext.errmsg[2] = 3;
            simContext.yH2S = 1.0;
        }
    }

    if (use_pH == 3 && UseH2Sgas == 1 && useEOS == 0) {
        pHHigh = 14.0;
        pHLow = 0.0;

        for (k = 0; k < 30; k++) {
            pH = (pHHigh + pHLow) / 2.0;
            aH = pow(10.0, -(pH));

            // 计算H+和OH-浓度
            H = aH / (gCat[iH] * gNCat[iH]);
            OH = KH2O / (aH * gAn[iOH] * gNAn[iOH]);

            // 计算碳酸系统物种 (从TCO2)
            simContext.HCO3 = TCO2 / (aH / (K1H2CO3 * aH2O) / gNeut[iCO2aq] / gNNeut[iCO2aq] * gAn[iHCO3] * gNAn[iHCO3] +
                1.0 + K2HCO3 / aH / gAn[iCO3] / gNAn[iCO3] * gAn[iHCO3] * gNAn[iHCO3]);
            simContext.CO2aq = aH * simContext.HCO3 * gAn[iHCO3] * gNAn[iHCO3] / ((K1H2CO3 * aH2O) * gNeut[iCO2aq] * gNNeut[iCO2aq]);
            simContext.CO3 = K2HCO3 * simContext.HCO3 * gAn[iHCO3] * gNAn[iHCO3] / (aH * gAn[iCO3] * gNAn[iCO3]);

            // 计算硫系统物种
            simContext.H2Saq = KgwH2S * Ppsia * (simContext.yH2S) * gGas[iH2Sg] / (gNeut[iH2Saq] * gNNeut[iH2Saq]);
            simContext.HS = K1H2S * simContext.H2Saq * gNeut[iH2Saq] * gNNeut[iH2Saq] /
                (aH * gAn[iHS] * gNAn[iHS]);
            S = K2HS * simContext.HS * gAn[iHS] * gNAn[iHS] /
                (aH * gAn[iSion] * gNAn[iSion]);

            // 计算其他弱酸物种
            hydAc = aH * gAn[iAc] * gNAn[iAc] / (KHAc * gNeut[iHAcaq] * gNNeut[iHAcaq]) + 1.0;
            AC = TAc / hydAc;
            HAcaq = TAc - AC;

            hydH2BO3 = aH * gAn[iH2BO3] * gNAn[iH2BO3] / (KH3BO3 * gNeut[iH3BO3] * gNNeut[iH3BO3]) + 1.0;
            H2BO3 = TH3BO3 / hydH2BO3;

            hydNH3 = aH * gNeut[iNH3] * gNNeut[iNH3] / (KNH4 * gCat[iNH4] * gNCat[iNH4]) + 1.0;
            NH3 = TNH4 / hydNH3;

            hydH2SiO4 = aH * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]) +
                aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]) + 1.0;

            H2SiO4 = TH4SiO4 / hydH2SiO4;
            H3SiO4 = H2SiO4 * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]);
            H4SiO4 = H2SiO4 * aH * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]);

            // 计算碱度残差
            faH = Alk - (simContext.HCO3 + 2.0 * simContext.CO3 + simContext.HS + 2.0 * S + AC + NH3 +
                H2BO3 + H3SiO4 + 2.0 * H2SiO4 + OH - H);

            // 二分法更新pH范围
            if (faH > 0) pHLow = pH;
            else pHHigh = pH;
        }

        simContext.yCO2 = simContext.CO2aq * gNeut[iCO2aq] * gNNeut[iCO2aq] / (KgwCO2 * Ppsia * gGas[iCO2g]);
        mc[iFe] = TFe / (1.0 + KstFeSaq * simContext.HS * gAn[iHS] * gNAn[iHS] *
            gCat[iFe] * gNCat[iFe] / (gNeut[iFeSaq] * gNNeut[iFeSaq] * aH));
        mn[iFeSaq] = KstFeSaq * mc[iFe] * simContext.HS * gAn[iHS] * gNAn[iHS] *
            gCat[iFe] * gNCat[iFe] / (gNeut[iFeSaq] * gNNeut[iFeSaq] * aH);
        pHMeterReading = pH - DpHj;
    }
}



void fmn() {
    // 初始化中性物种浓度
    mn[0] = 0;           // 原mn(1)对应mn[0]
    mn[1] = simContext.CO2aq;       // 原mn(2)对应mn[1]
    mn[2] = simContext.H2Saq;       // 原mn(3)对应mn[2]
    mn[3] = HAcaq;       // 原mn(4)对应mn[3]
    mn[4] = H4SiO4;      // 原mn(5)对应mn[4]
    mn[5] = NH3;         // 原mn(6)对应mn[5]
    mn[6] = 0;           // 原mn(7)对应mn[6]

    // 计算FeSaq浓度（基于平衡常数）
    mn[7] = KstFeSaq * mc[iFe] * simContext.HS * gAn[iHS] * gNAn[iHS] * gCat[iFe] * gNCat[iFe] /
        (gNeut[iFeSaq] * gNNeut[iFeSaq] * aH);  // 原mn(8)对应mn[7]

    // 设置阴离子浓度
    ma[0] = OH;          // 原ma(1)对应ma[0]
    ma[2] = AC;          // 原ma(3)对应ma[2]
    ma[3] = simContext.HCO3;        // 原ma(4)对应ma[3]
    ma[4] = simContext.CO3;         // 原ma(5)对应ma[4]
    ma[6] = simContext.HS;          // 原ma(7)对应ma[6]
    ma[9] = H2BO3;       // 原ma(10)对应ma[9]
    ma[10] = H3SiO4;     // 原ma(11)对应ma[10]
    ma[11] = H2SiO4;     // 原ma(12)对应ma[11]

    // 设置阳离子浓度
    mc[0] = H;           // 原mc(1)对应mc[0]
}


double fAphicalc(double tk, double pBar)
{
    double U1 = 342.79, U2 = -0.0050866, U3 = 0.0000009469;
    double U4 = -2.0525, U5 = 3115.9, U6 = -182.89;
    double U7 = -8032.5, U8 = 4214200.0, U9 = 2.1417;

    double D1000, cc, b, Dielec;
    double term1, term2, result;

    D1000 = U1 * exp(U2 * tk + U3 * tk * tk);
    cc = U4 + U5 / (U6 + tk);
    b = U7 + U8 / tk + U9 * tk;

    Dielec = D1000 + cc * log((b + pBar) / (b + 1000.0));

    dens = fH2ODensity(tk, pBar);  // 单位：kg/m³

    term1 = sqrt(2.0 * pi * NAv * dens);
    term2 = pow(eElec * eElec / (4.0 * pi * eps0 * Dielec * kBoltz * tk), 1.5);

    result = (1.0 / 3.0) * term1 * term2;

    return result;
}


// b0,b1,b2 在C2_Pitzer2019里面计算
void fBtermcalc(double** bterm, double gX14, double gX12)
{
    int m, a;
    double X20, gX20;
    int c = 13;  //用于占位

    for (m = 0; m < NumCat; ++m)
    {
        for (a = 0; a < NumAn; ++a)
        {
            // 基本 2:(-2) 型离子计算 bterm(15, 15)
            bterm[m][a] = b0[m][a] + b1[m][a] * gX14 + b2[m][a] * gX12;

            // 1:(-1) 或 1:(-2) 阳离子处理
            if (ChCat[m] == 1.0)
            {
                X20 = 2.0 * sqrt(ISt);

                // 特例 Na–SO4, K–SO4
                // 注意：VB 原代码使用 c=2,3，这里推测应为 m==1-based索引，对应 m=1,2。  
                // c在原代码里用了，但是不知道是什么
                if ((c == 1 && a == 5) || (c == 2 && a == 5))  // 转成 0-based 索引
                    X20 = 1.4 * sqrt(ISt);

                gX20 = 2.0 * (1.0 - (1.0 + X20) * exp(-X20)) / (X20 * X20);
                //b0b1b2在更高一层的B1、C2_PitzerActCoefsConstants被赋值过，若要使bterm保持一致，必须引入！！
                bterm[m][a] = b0[m][a] + b1[m][a] * gX20 + b2[m][a] * gX12;
            }

            // 2:(-1) 型电解质 (如 CaCl2)
            if (ChCat[m] == 2.0 && ChAn[a] == -1.0)
            {
                X20 = (2.0 - 0.00181 * (TK - 298.15)) * sqrt(ISt);
                gX20 = 2.0 * (1.0 - (1.0 + X20) * exp(-X20)) / (X20 * X20);
                bterm[m][a] = b0[m][a] + b1[m][a] * gX20 + b2[m][a] * gX12;
            }
        }
    }
    //printf("\nbterm:\n");
    //for (m = 0; m < NumCat; ++m)
    //{
    //    printf("[");
    //    for (a = 0; a < NumAn; ++a)
    //    {
    //        printf("%lf,", bterm[m][a]);
    //    }
    //    printf("]\n");
    //}
}




double CalcRhoTP(double tk, double tc, double pBar, double patm) {
    // 计算由于压力引起的过量摩尔体积变化，计算由于压力 1 atm 到 1001 atm 而引起的活度系数的导数
    // 按设定的tk值计算。如果不是STP条件（patm !=1），则计算Δ压力=1e-6 bar时的过量性能

    double AphiP, AphiPPlus, X14, gX14, gpX14, X20, gX20, gpX20, X12, gX12, gpX12;
    double mt, Av, Fv, Ex_Pitzer, V_ex, V_ion, VperKgWater, MassperKgwater;
    int m, a, c, n, iden;

    C2_Pitzer2019(tk, tc, pBar, patm);
    AphiP = fAphicalc(tk, pBar);

    // 计算各种X值和对应的gX、gpX
    X14 = 1.4 * sqrt(ISt); // For 2:(-2) pairs or ions
    gX14 = 2 * (1 - (1 + X14) * exp(-X14)) / (X14 * X14);
    gpX14 = -2 * (1 - (1 + X14 + 0.5 * X14 * X14) * exp(-X14)) / (X14 * X14);

    X20 = 2 * sqrt(ISt); // For 1:(-2), (-1):2, or 1:(-1) pairs
    gX20 = 2 * (1 - (1 + X20) * exp(-X20)) / (X20 * X20);
    gpX20 = -2 * (1 - (1 + X20 + 0.5 * X20 * X20) * exp(-X20)) / (X20 * X20);

    X12 = 12 * sqrt(ISt);
    gX12 = 2 * (1 - (1 + X12) * exp(-X12)) / (X12 * X12);
    gpX12 = -2 * (1 - (1 + X12 + 0.5 * X12 * X12) * exp(-X12)) / (X12 * X12);

    //这六个数组都是局部的, 均为15*15
    double** btermP = (double**)malloc(15 * sizeof(double*));
    double** CtermP = (double**)malloc(15 * sizeof(double*));
    double** btermPPlus = (double**)malloc(15 * sizeof(double*));
    double** CtermPPlus = (double**)malloc(15 * sizeof(double*));
    double** bVterm = (double**)malloc(15 * sizeof(double*));
    double** cVterm = (double**)malloc(15 * sizeof(double*));

    //bterm（15*15）仅在fBtermcalc赋值，并且仅在CalcRhoTP使用，因此局部化
    double** bterm = (double**)malloc(15 * sizeof(double*));

    for (int i = 0; i < 15; i++) {
        btermP[i] = (double*)malloc(15 * sizeof(double));
        CtermP[i] = (double*)malloc(15 * sizeof(double));
        btermPPlus[i] = (double*)malloc(15 * sizeof(double));
        CtermPPlus[i] = (double*)malloc(15 * sizeof(double));
        bVterm[i] = (double*)malloc(15 * sizeof(double));
        cVterm[i] = (double*)malloc(15 * sizeof(double));

        bterm[i] = (double*)malloc(15 * sizeof(double));
    }

    fBtermcalc(bterm, gX14, gX12);

    // 保存当前bterm和CPhi值
    for (m = 0; m < NumCat; m++) {
        for (a = 0; a < NumAn; a++) {
            btermP[m][a] = bterm[m][a];
            CtermP[m][a] = CPhi[m][a];
        }
    }

    // 增加压力并重新计算
    pBar = pBar + 0.000001;
    patm = pBar / 1.013254;
    Ppsia = pBar * 14.503774;

    C2_Pitzer2019(tk, tc, pBar, patm);

    AphiPPlus = fAphicalc(tk, pBar);
    fBtermcalc(bterm, gX14, gX12);

    // 保存增加压力后的值
    for (m = 0; m < NumCat; m++) {
        for (a = 0; a < NumAn; a++) {
            btermPPlus[m][a] = bterm[m][a];
            CtermPPlus[m][a] = CPhi[m][a];
        }
    }

    // 恢复原始压力
    pBar = pBar - 0.0000001;
    patm = pBar / 1.013254;
    Ppsia = pBar * 14.503774;

    V0TP(tk, pBar); // 在T,P下重新计算V0

    // 计算体积相关项
    Av = -4.0 * RBar * tk * (AphiPPlus - AphiP) / 0.000001;
    Fv = Av / RBar / tk * ISt / 1.2 * log(1 + 1.2 * sqrt(ISt)); // Unit mol/Kg/bar

    // 计算bVterm和cVterm
    for (m = 0; m < NumCat; m++) {
        for (a = 0; a < NumAn; a++) {
            bVterm[m][a] = (btermPPlus[m][a] - btermP[m][a]) / 0.000001; // Unit Kg/mol/bar
            cVterm[m][a] = (CtermPPlus[m][a] - CtermP[m][a]) / 0.000001; // Unit (Kg/mol)^2/bar
        }
    }

    // 计算Pitzer超额项
    Ex_Pitzer = 0;
    for (m = 0; m < NumCat; m++) {
        for (a = 0; a < NumAn; a++) {
            Ex_Pitzer = Ex_Pitzer + 2 * mc[m] * ma[a] *
                (bVterm[m][a] + mc[m] * ChCat[m] * cVterm[m][a] /
                    (2 * sqrt(fabs(ChCat[m] * ChAn[a])))); // unit mol/Kg/bar
        }
    }

    // 计算超额体积和离子体积
    V_ex = (Fv + Ex_Pitzer) * RBar * tk; // L/Kg
    V_ion = V_ex * 1000; // cm3/Kg

    // 添加各离子的偏摩尔体积贡献
    for (a = 0; a < NumAn; a++) {
        V_ion = V_ion + ma[a] * V0_a[a]; // Note Unit of V0 is cm3/Kg
    }

    for (c = 0; c < NumCat; c++) {
        V_ion = V_ion + mc[c] * V0_c[c];
    }

    for (n = 0; n < NumNeut; n++) {
        V_ion = V_ion + mn[n] * V0_n[n];
    }

    dens = fH2ODensity(tk, pBar); // g/L density of pure water at T, P

    // 计算总体积和质量
    VperKgWater = (1.0 / dens * 1000000.0) + V_ion; // volume in cm3 of 1 kg water at T
    MassperKgwater = 1000; // g of water of 1 kg water

    // 添加各离子的质量贡献
    for (iden = 0; iden < NumCat; iden++) {
        MassperKgwater = MassperKgwater + mc[iden] * MWCat[iden];
    }

    for (iden = 0; iden < NumAn; iden++) {
        MassperKgwater = MassperKgwater + ma[iden] * MWAn[iden];
    }

    for (iden = 0; iden < NumNeut; iden++) {
        MassperKgwater = MassperKgwater + mn[iden] * MWNeut[iden]; // g
    }

    //释放内存
    int i;
    if (btermP) {
        for (i = 0; i < 15; i++) free(btermP[i]);
        free(btermP);
    }
    if (CtermP) {
        for (i = 0; i < 15; i++) free(CtermP[i]);
        free(CtermP);
    }
    if (btermPPlus) {
        for (i = 0; i < 15; i++) free(btermPPlus[i]);
        free(btermPPlus);
    }
    if (CtermPPlus) {
        for (i = 0; i < 15; i++) free(CtermPPlus[i]);
        free(CtermPPlus);
    }
    if (bVterm) {
        for (i = 0; i < 15; i++) free(bVterm[i]);
        free(bVterm);
    }
    if (cVterm) {
        for (i = 0; i < 15; i++) free(cVterm[i]);
        free(cVterm);
    }
    if (bterm) {
        for (i = 0; i < 15; i++) free(bterm[i]);
        free(bterm);
    }

    //修改全局的变量
    TK = tk;
    PBar = pBar;
    Patm = patm;

    // 返回密度 (g/cm3)
    return MassperKgwater / VperKgWater;
}


void fTotalCO2H2Smoles() {

    // 计算总CO2摩尔数
    nTCO2 = Patm * 14.696 * simContext.yCO2 * (829.0 * VgTP / (Znew * TK) +
        gGas[iCO2g] * (KgwCO2 * mass_w / gNeut[iCO2aq] / gNNeut[iCO2aq] +
            RatioOilBPoints * KgoCO2 * Mass_o / gL[iCO2o])) +
        (simContext.HCO3 + simContext.CO3) * mass_w;

    // nTCO2EOS 等于气体中的CO2加上油中的CO2加上水中的CO2，不包括HCO3和CO3
    nTCO2EOS = nTCO2;

    // 计算总CH4摩尔数
    nTCH4 = Patm * 14.696 * yCH4 * (829.0 * VgTP / (Znew * TK) +
        gGas[iCH4g] * (KgwCH4 * mass_w / gNeut[iCH4aq] +
            RatioOilBPoints * KgoCH4 * Mass_o / gL[iCH4g]));

    // 计算总H2S摩尔数（包括HS-和FeS配合物）
    // nTH2S vb 536.894123326537 cpp 537.8953812802223
    // vb simContext.yH2S 5.41442967642879E-05 cpp 5.4245428347734764e-05
    // mc[iFe] vb 0.000219991780725818 cpp 0.00011386369818347661 QualityControlCalculations D1_CalcDensity C5_CalcpHPCO2PH2SSTP
    // gAn[iHS] vb 0.645200696723431 cpp 0.64535500358708076
    // gCat[iFe] vb 0.111366817290854 cpp 0.1113832023689522
    nTH2S = Patm * 14.696 * simContext.yH2S * (829.0 * VgTP / (Znew * TK) +
        gGas[iH2Sg] * (KgwH2S * mass_w / gNeut[iH2Saq] / gNNeut[iH2Saq] +
            RatioOilBPoints * KgoH2S * Mass_o / gL[iH2So])) +
        simContext.HS * (1.0 + KstFeSaq * mc[iFe] * simContext.HS * gAn[iHS] * gNAn[iHS] *
            gCat[iFe] * gNCat[iFe] / (gNeut[iFeSaq] * gNNeut[iFeSaq] * aH)) * mass_w;

    // nTH2sEOS 等于气体、油和水中的H2S，不包括HS-
    nTH2sEOS = nTH2S;
}


/*      存在 内存泄露 的可能性  ， 详细查看代码：if (!(*mf_ParametersWereRead))    - 彭非  */
void InitialPreparationSSP(
    bool& mf_ParametersWereRead,
    const char* EOS,
    DoubleVec& zInput,    // 输入组分
    IntVec& iFlash,      // 闪蒸索引
    DoubleVec& zGlobal,   // 全局组分
    DoubleVec& mf_TCr,
    DoubleVec& mf_PCr,
    DoubleVec& mf_Omega,
    DoubleVec& mf_MWgas,
    DoubleMatrix& mf_kPr,
    DoubleVec& mf_c0,
    DoubleVec& mf_c1,
    int max_NumGases,
    int NumGases
)
{
    int i, j, k;
    double zSum = 0.0;

    for (i = 0; i < 4; i++) {
        iFlash[i] = i;
        zGlobal[i] = zInput[i];
    }
    k = 4;
    for (i = 4; i < max_NumGases; i++) {
        if (zInput[i] > 0) {
            iFlash[k] = i;
            zGlobal[k] = zInput[i];
            k++;
        }
        else if (i >= max_NumGases - 2) {
            iFlash[k] = i;
            zGlobal[k] = zInput[i];
            k++;
        }
    }
    // -------------------------------------------
    //考虑=》mf_相关变量是从MultiPhaseFlash传过来的，此处又要清理空间，重新New出空间，且请确保-mf_通过指针传入
    // 
    // 举例子：MultiPhaseFlash的数组int* mf_a，若要在InitialPreparationSSP中ReDim，请传入int** mf_a，即数组mf_a的地址，请务必仔细考虑
    // 
    // by   - 彭非
    //--------------------------------------------
    if (!mf_ParametersWereRead) {

        //mf_kPr比较特殊，我不知道它的行列是否一直是同一个数组，如果上次调用MultiPhaseFlash得到的NumGases是6，而这次的NumGases是5，绝对会出现内存泄露。

        //源代码：ReDim mf_MWgas(NumGases), mf_TCr(NumGases), mf_PCr(NumGases), mf_Omega(NumGases), mf_c0(NumGases), mf_c1(NumGases), mf_kPr(NumGases, NumGases)
        mf_MWgas.resize(NumGases, 0);
        mf_TCr.resize(NumGases, 0);
        mf_PCr.resize(NumGases, 0);
        mf_Omega.resize(NumGases, 0);
        mf_c0.resize(NumGases, 0);
        mf_c1.resize(NumGases, 0);
        mf_kPr = DoubleMatrix(NumGases, DoubleVec(NumGases, 0));
        /*
        for (i = 0; i < NumGases; i++) {
            // mf_MWgas[i] = Worksheets("Input").Cells[4 + iFlash[i]][24];
            // mf_TCr[i]   = Worksheets("Input").Cells[4 + iFlash[i]][26];
            // mf_PCr[i]   = Worksheets("Input").Cells[4 + iFlash[i]][27];
            // mf_Omega[i] = Worksheets("Input").Cells[4 + iFlash[i]][28];

            if (strcmp(EOS, "PR") == 0) {
                // mf_c0[i] = Worksheets("Input").Cells[4 + iFlash[i]][29];
                // mf_c1[i] = Worksheets("Input").Cells[4 + iFlash[i]][30];
            }

            if (strcmp(EOS, "SRK") == 0) {
                // mf_c0[i] = Worksheets("Input").Cells[4 + iFlash[i]][31];
                // mf_c1[i] = Worksheets("Input").Cells[4 + iFlash[i]][32];
            }

            for (j = 0; j < i; j++) {
                // mf_kPr[i][j] = Worksheets("Input").Cells[4 + iFlash[i]][34 + iFlash[j]];
            }
        }
        */
        //手动模拟上面的输入
        mf_MWgas[0] = 16.043;
        mf_TCr[0] = 190.6;
        mf_PCr[0] = 46;
        mf_Omega[0] = 0.008;

        mf_MWgas[1] = 44.01;
        mf_TCr[1] = 304.2;
        mf_PCr[1] = 73.76;
        mf_Omega[1] = 0.225;

        mf_MWgas[2] = 34.08;
        mf_TCr[2] = 373.5;
        mf_PCr[2] = 89.63;
        mf_Omega[2] = 0.094;

        mf_MWgas[3] = 30.07;
        mf_TCr[3] = 305.5;
        mf_PCr[3] = 48.84;
        mf_Omega[3] = 0.098;

        mf_MWgas[4] = 44.094;
        mf_TCr[4] = 369.9;
        mf_PCr[4] = 42.46;
        mf_Omega[4] = 0.152;

        mf_MWgas[5] = 58.124;
        mf_TCr[5] = 408.2;
        mf_PCr[5] = 36.48;
        mf_Omega[5] = 0.176;

        mf_MWgas[6] = 58.124;
        mf_TCr[6] = 425.3;
        mf_PCr[6] = 38;
        mf_Omega[6] = 0.193;

        mf_MWgas[7] = 72.151;
        mf_TCr[7] = 460.45;
        mf_PCr[7] = 33.84;
        mf_Omega[7] = 0.227;

        mf_MWgas[8] = 72.151;
        mf_TCr[8] = 469.6;
        mf_PCr[8] = 33.74;
        mf_Omega[8] = 0.251;

        mf_MWgas[9] = 86.178;
        mf_TCr[9] = 507.5;
        mf_PCr[9] = 29.69;
        mf_Omega[9] = 0.296;

        mf_MWgas[10] = 108.46355910141;
        mf_TCr[10] = 575.5;
        mf_PCr[10] = 26.26;
        mf_Omega[10] = 0.545;

        mf_MWgas[11] = 203.768522560204;
        mf_TCr[11] = 708.6;
        mf_PCr[11] = 16.76;
        mf_Omega[11] = 0.857;

        mf_MWgas[12] = 321.359072235909;
        mf_TCr[12] = 945.6;
        mf_PCr[12] = 13.36;
        mf_Omega[12] = 1.248;

        mf_MWgas[13] = 28.013;
        mf_TCr[13] = 126.2;
        mf_PCr[13] = 33.94;
        mf_Omega[13] = 0.040;

        mf_MWgas[14] = 18;
        mf_TCr[14] = 647.1;
        mf_PCr[14] = 220.55;
        mf_Omega[14] = 0.345;

        if (strcmp(EOS, "PR") == 0) {
            mf_c0[0] = -5.19998589998324;
            mf_c1[0] = 0.0;

            mf_c0[1] = -1.90887591231191;
            mf_c1[1] = 0.0;

            mf_c0[2] = -3.92147265302369;
            mf_c0[3] = -5.7950356085208;
            mf_c0[4] = -6.35369532428396;
            mf_c0[5] = -7.18062223841032;
            mf_c0[6] = -6.48762839907065;
            mf_c0[7] = -6.19849948277172;
            mf_c0[8] = -5.12105162468565;
            mf_c0[9] = -3.4814313848538;
            mf_c0[10] = -14.8622158330974;
            mf_c0[11] = -29.8260575190994;
            mf_c0[12] = 5.77517376259963;
            mf_c0[13] = -4.23208315053643;
            mf_c0[14] = 3.00773525559068;

            mf_c1[2] = 0.0;
            mf_c1[3] = 0.0;
            mf_c1[4] = 0.0;
            mf_c1[5] = 0.0;
            mf_c1[6] = 0.0;
            mf_c1[7] = 0.0;
            mf_c1[8] = 0.0;
            mf_c1[9] = 0.0;
            mf_c1[10] = -0.00560031997876096;
            mf_c1[11] = 0.101896336756916;
            mf_c1[12] = 0.279177061769215;
            mf_c1[13] = 0.0;
            mf_c1[14] = 0.00873186552687633;
        }
        if (strcmp(EOS, "SRK") == 0) {
            // ??
        }

        mf_kPr[0][0] = 0;

        mf_kPr[1][0] = 0.12;

        mf_kPr[2][0] = 0.08;
        mf_kPr[2][1] = 0.12;

        mf_kPr[3][1] = 0.15;
        mf_kPr[3][2] = 0.07;

        mf_kPr[4][1] = 0.15;
        mf_kPr[4][2] = 0.07;

        mf_kPr[5][1] = 0.15;
        mf_kPr[5][2] = 0.06;

        mf_kPr[6][1] = 0.15;
        mf_kPr[6][2] = 0.06;

        mf_kPr[7][1] = 0.15;
        mf_kPr[7][2] = 0.06;

        mf_kPr[8][1] = 0.15;
        mf_kPr[8][2] = 0.06;

        mf_kPr[9][1] = 0.15;
        mf_kPr[9][2] = 0.05;

        mf_kPr[10][1] = 0.15;
        mf_kPr[10][2] = 0.05;

        mf_kPr[11][1] = 0.15;
        mf_kPr[11][2] = 0.05;

        mf_kPr[12][1] = 0.15;
        mf_kPr[12][2] = 0.05;

        mf_kPr[13][0] = 0.02;
        mf_kPr[13][3] = 0.06;
        mf_kPr[13][4] = 0.08;
        mf_kPr[13][5] = 0.08;
        mf_kPr[13][6] = 0.08;
        mf_kPr[13][7] = 0.08;
        mf_kPr[13][8] = 0.08;
        mf_kPr[13][9] = 0.08;
        mf_kPr[13][10] = 0.08;
        mf_kPr[13][11] = 0.08;
        mf_kPr[13][12] = 0.08;

        mf_kPr[14][0] = 0.45;
        mf_kPr[14][3] = 0.45;
        mf_kPr[14][4] = 0.53;
        mf_kPr[14][5] = 0.52;
        mf_kPr[14][6] = 0.52;
        mf_kPr[14][7] = 0.5;
        mf_kPr[14][8] = 0.5;
        mf_kPr[14][9] = 0.5;
        mf_kPr[14][10] = 0.5;
        mf_kPr[14][11] = 0.5;
        mf_kPr[14][12] = 0.5;
        mf_ParametersWereRead = true;
    }

    // --- 4. 归一化全局摩尔分数 ---
    zSum = 0.0;
    for (i = 0; i < NumGases; i++)
        zSum += zGlobal[i];

    for (j = 0; j < NumGases; j++)
        zGlobal[j] /= zSum;
}


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



//---------------------------------------------闫师兄的：

bool isAqueous_2016(bool eqAqueous, DoubleVec& composition, DoubleVec& globalComposition, int NumGases)
{
    double composition_ratio = 1.001; // 用于比较的比例因子
    double tol = composition_ratio - 1.0;
    double total_dry = 0.0, total_comp_dry = 0.0;
    double global_comp_Water = 0.0, composition_Water = 0.0;
    int NonWaterComponents = 0;
    int i, nComponents = NumGases; //sizeof(composition) / sizeof(composition[0]); // 数组长度由调用方控制

    if (globalComposition[nComponents - 1] > 0)
    { // 如果全局水组分大于0
        for (i = 0; i < nComponents - 1; i++)
        {
            total_dry += globalComposition[i];
            total_comp_dry += composition[i];
            if (composition[i] > 0.0)
                NonWaterComponents++;
        }

        global_comp_Water = globalComposition[nComponents - 1] / (total_dry + globalComposition[nComponents - 1]);
        composition_Water = composition[nComponents - 1] / (total_comp_dry + composition[nComponents - 1]);
        total_comp_dry = total_comp_dry / (total_comp_dry + composition[nComponents - 1]);

        if (total_comp_dry > 0.0)
        {
            if (composition_Water > total_comp_dry)
            {
                if ((composition_Water - composition_ratio * global_comp_Water) > 0.0)
                {
                    return true;
                }
                else if (composition_Water > 0.5 && NonWaterComponents == 1)
                {
                    return true;
                }
                else if (eqAqueous && (composition_Water - global_comp_Water < tol * global_comp_Water))
                {
                    return true;
                }
                else if (eqAqueous && composition_Water > 0.99 && (composition_Water - global_comp_Water < 0.001))
                {
                    return true;
                }
            }
        }
        else if (composition[nComponents - 1] > 0.0)
        {
            return true;
        }
    }

    return false;
}



static double arcCos(double x)
{
    double tol = 1.0 - 1e-14;  // 容差 0.00000000000001
    double result;

    if (fabs(x) >= tol)
        result = 0.0;
    else
        result = atan(-x / sqrt(-x * x + 1.0)) + (2.0 * atan(1.0));

    return result;
}

static double QBRT(double x)
{
    return (x >= 0.0) ? pow(x, 1.0 / 3.0) : -pow(-x, 1.0 / 3.0);
}

void QUBIC(double P, double Q, double R, double* root, int* root_num)
{
    const double DEG120 = 2.09439510239319;  // 120度 (弧度)
    const double tolerance = 0.00001;
    const double Tol2 = 1E-20;

    double z[3] = { 0 };
    double p2, RMS, a, b, DISCR, t1, t2, Ratio;
    double sum, DIF, AD3, E0, CPhi, PhiD3, PD3;
    int i, nRoots;

    // --- 1. 转换为标准形式 Z^3 + aZ + b = 0 ---
    p2 = P * P;
    a = Q - p2 / 3.0;
    b = P * (2.0 * p2 - 9.0 * Q) / 27.0 + R;
    RMS = sqrt(a * a + b * b);

    // --- 2. 三重根情形 ---
    if (RMS < Tol2)
    {
        nRoots = 3;
        for (i = 0; i < 3; ++i)
            root[i] = -P / 3.0;
        *root_num = nRoots;
        return;
    }

    // --- 3. 判别式 ---
    DISCR = pow(a / 3.0, 3.0) + pow(b / 2.0, 2.0);

    if (DISCR > 0.0)
    {
        // 判别式 > 0 → 一个实根，两个复根
        t1 = -b / 2.0;
        t2 = sqrt(DISCR);

        Ratio = (t1 == 0.0) ? 1.0 : t2 / t1;

        if (fabs(Ratio) < tolerance)
        {
            // 三个实根，其中两个相等
            nRoots = 3;
            z[0] = 2.0 * QBRT(t1);
            z[1] = QBRT(-t1);
            z[2] = z[1];
        }
        else
        {
            nRoots = 1;
            sum = t1 + t2;
            DIF = t1 - t2;
            z[0] = QBRT(sum) + QBRT(DIF);
        }
    }
    else
    {
        // 判别式 <= 0 → 三个不同实根（使用三角法）
        nRoots = 3;
        AD3 = a / 3.0;
        E0 = 2.0 * sqrt(-AD3);
        CPhi = -b / (2.0 * sqrt(-pow(AD3, 3.0)));
        PhiD3 = arcCos(CPhi) / 3.0;

        z[0] = E0 * cos(PhiD3);
        z[1] = E0 * cos(PhiD3 + DEG120);
        z[2] = E0 * cos(PhiD3 - DEG120);
    }

    // --- 4. 还原到原始方程的根 ---
    PD3 = P / 3.0;
    for (i = 0; i < nRoots; ++i)
        root[i] = z[i] - PD3;

    *root_num = nRoots;
}


void cubic(double a3, double a2, double a1, double a0, const char* min_or_max,
    double* z, double bm, double PBar, double TK)
{

    const double RBar = 83.144621; // cm3*bar/(mol*K)
    double P, Q, R;
    double root[3]; // 最多三根
    int numRoots, i;

    P = a2 / a3;
    Q = a1 / a3;
    R = a0 / a3;

    // 调用外部函数求解三次方程
    QUBIC(P, Q, R, root, &numRoots);

    if (strcmp(min_or_max, "max") == 0)
    {
        *z = root[0];
        for (i = 1; i < numRoots; i++)
        {
            if (root[i] > *z)
            {
                *z = root[i];
            }
        }
    }
    else
    { // "min"
        *z = root[0];
        for (i = 1; i < numRoots; i++)
        {
            if ((root[i] < *z) && ((root[i] * TK * RBar / PBar) > bm))
            {
                *z = root[i];
            }
        }
    }
}


double logKHenry_CH4(double TK)
{
    const double TC = 647.096;   // Tc for water Kelvin
    const double a = -10.44708;
    const double b = 4.66491;
    const double c = 12.12986;

    double Tr, tau;

    Tr = TK / TC;
    tau = 1.0 - Tr;

    return log(PsatH2O(TK)) +
        (a / Tr +
            b * pow(tau, 0.355) / Tr +
            c * pow(Tr, -0.41) * exp(tau));
}


double KH_CO2(double TK, double PBar)
{
    double dV, RBar;
    RBar = 0.083144;

    double KH_CO2_val;

    KH_CO2_val = 55.508 / exp(
        log(PsatH2O(TK))
        - 9.14122 / (TK / 647.096)
        + 2.8192 * pow(1.0 - TK / 647.096, 0.355) / (TK / 647.096)
        + 11.28516 * pow(TK / 647.096, -0.41) * exp(1.0 - TK / 647.096)
        - 0.8066
    );

    dV =
        (37.88
            - 0.14442 * (TK - 273.15)
            + 0.001243 * pow(TK - 273.15, 2)
            - 0.0000044738 * pow(TK - 273.15, 3)
            + 0.000000005726 * pow(TK - 273.15, 4)) * 0.001;

    if (TK <= 373.15)
    {
        KH_CO2_val = KH_CO2_val * pow(
            exp(dV * (PBar - 1.0) / RBar / TK
                + 0.0000067 / (2.0 * RBar * TK) * pow(PBar - 1.0, 2)),
            -1.0);
    }
    else
    {
        KH_CO2_val = KH_CO2_val * pow(
            exp(dV * (PBar - PsatH2O(TK)) / RBar / TK
                + 0.0000067 / (2.0 * RBar * TK) * pow(PBar - PsatH2O(TK), 2)),
            -1.0);
    }

    return KH_CO2_val;
}


double KH_H2S(double TK, double PBar)
{
    double dV, ln_H_T_H2S, RBar;
    RBar = 0.083144;

    ln_H_T_H2S =
        (13788.0 / TK
            - 185.19
            + 29.087 * log(TK)
            - 0.027637 * TK
            - 1445200.0 / pow(TK, 2));

    dV =
        (33.18
            + 0.092661 * (TK - 273.15)
            - 0.00054853 * pow(TK - 273.15, 2)
            + 0.0000015354 * pow(TK - 273.15, 3)
            - 0.0000000015459 * pow(TK - 273.15, 4)) * 0.001;

    if (TK <= 373.15)
    {
        return pow(exp(ln_H_T_H2S) * exp(dV * (PBar - 1.0) / RBar / TK), -1.0);
    }
    else
    {
        return pow(exp(ln_H_T_H2S) * exp(dV * (PBar - PsatH2O(TK)) / RBar / TK), -1.0);
    }
}


double logKHenry_C2H6(double TK)
{
    const double TC = 647.096;
    const double a = -19.67563;
    const double b = 4.51222;
    const double c = 20.62567;

    double Tr = TK / TC;
    double tau = 1.0 - Tr;

    return log(PsatH2O(TK)) +
        (a / Tr
            + b * pow(tau, 0.355) / Tr
            + c * pow(Tr, -0.41) * exp(tau));
}

double KH_CH4(double TK, double PBar)
{
    double TF;
    TF = ((9.0 / 5.0) * (TK - 273.15)) + 32.0;

    return 14.504 * pow(10.0,
        -1.0 * (3.798
            + 0.004136 * TF
            - 0.000009486 * pow(TF, 2)
            - 0.000038 * PBar * 14.504));
}

double logKHenry_N2(double TK)
{
    const double TC = 647.096;
    const double a = -9.67578;
    const double b = 4.72162;
    const double c = 11.70585;

    double Tr = TK / TC;
    double tau = 1.0 - Tr;

    return log(PsatH2O(TK)) +
        (a / Tr
            + b * pow(tau, 0.355) / Tr
            + c * pow(Tr, -0.41) * exp(tau));
}


void phi_calc(
    bool eqVapor, bool eqAqueous, const std::string& EOS, const std::string& phase,
    double TK, double PBar, DoubleVec& x, DoubleVec& xGlobal,
    DoubleVec& gNeut, double aH2O, DoubleVec& TCr, DoubleVec& PCr,
    DoubleVec& Omega, DoubleVec& c0, DoubleVec& c1, DoubleMatrix& kPr,
    DoubleVec& lnPHI, double* z, int NumGases)
{
    //注意：源代码中是Const RBar As Double = 83.14 ，这个RBar和全局的那个重名但数值不一样，这里改成temp_RBar
    const double RBar = 83.14; // cm3·bar/(mol·K)
    DoubleVec Tr(NumGases, 0);

    DoubleVec aPR(NumGases, 0);

    DoubleVec bPR(NumGases, 0);

    DoubleVec lnPHI_Peneloux(NumGases, 0);

    double aPrmix = 0.0, bm = 0.0, cm = 0.0, F_Omega = 0, BStarPr, QQ;
    double sigma, epsilon, sigmaP, epsilonP, S;
    double a0, a1, a2, a3, asum, ai_bar, qi_bar, ii, Z_Peneloux;
    int i, j, k, comp;
    double sum_X = 0;
    double Peneloux_yes_no = 1;
    // 初始化混合物 a 和 b
    for (i = 0; i < NumGases; i++)
    {
        aPR[i] = 0.0;
        bPR[i] = 0.0;
        lnPHI[i] = 0.0;
        lnPHI_Peneloux[i] = 0.0;
    }

    for (i = 0; i < NumGases; i++)
        sum_X += x[i];

    if (sum_X > 1e-20)
    {
        // 计算各组分 a、b 以及混合物 bm、cm
        for (k = 0; k < NumGases; k++)
        {
            if (xGlobal[k] > 0.0)
            {
                Tr[k] = TK / TCr[k];
                if (EOS == "SRK")
                {
                    aPR[k] = 0.42748 * RBar * RBar * TCr[k] * TCr[k] / PCr[k];
                    F_Omega = 0.48 + 1.574 * Omega[k] - 0.176 * Omega[k] * Omega[k];
                    bPR[k] = 0.08664 * RBar * TCr[k] / PCr[k];
                }
                else if (EOS == "PR")
                {
                    aPR[k] = 0.45724 * RBar * RBar * TCr[k] * TCr[k] / PCr[k];
                    F_Omega = 0.37464 + 1.54226 * Omega[k] - 0.26992 * Omega[k] * Omega[k];
                    bPR[k] = 0.07779 * RBar * TCr[k] / PCr[k];
                }

                aPR[k] = aPR[k] * pow(1.0 + F_Omega * (1.0 - sqrt(Tr[k])), 2.0);
                bm += x[k] * bPR[k];
                cm += x[k] * (c0[k] + c1[k] * (TK - 288.15));
            }
        }
        // 计算混合物 aPrmix
        aPrmix = 0.0;
        for (i = 0; i < NumGases; i++)
        {
            for (j = 0; j < i; j++)
            {
                aPrmix += 2.0 * x[i] * x[j] * sqrt(aPR[i] * aPR[j]) * (1.0 - kPr[i][j]);
            }
            aPrmix += x[i] * x[i] * aPR[i];
        }

        // 计算 BStarPr 和 QQ
        BStarPr = bm * PBar / (RBar * TK);
        QQ = aPrmix / (bm * RBar * TK);
        // EOS 系数 sigma, epsilon
        if (EOS == "SRK")
        {
            sigma = 1.0;
            epsilon = 0.0;
        }
        else if (EOS == "PR")
        {
            sigma = 1.0 + sqrt(2.0);
            epsilon = 1.0 - sqrt(2.0);
        }

        epsilonP = epsilon;
        sigmaP = sigma;
        S = 1.0;

        // 三次方程求 z
        a3 = 1.0;
        a2 = (epsilonP + sigmaP) * BStarPr - 1.0 - S * BStarPr;
        a1 = (epsilonP * sigmaP) * BStarPr * BStarPr - (epsilonP + sigmaP) * BStarPr - (epsilonP + sigmaP) * S * BStarPr * BStarPr + QQ * BStarPr;
        a0 = -sigmaP * epsilonP * BStarPr * BStarPr - sigmaP * epsilonP * BStarPr * BStarPr * BStarPr * S - QQ * S * BStarPr * BStarPr;

        if (phase == "vapor")
        {
            cubic(a3, a2, a1, a0, "max", z, bm, PBar, TK);
        }
        else if (phase == "liquid")
        {
            cubic(a3, a2, a1, a0, "min", z, bm, PBar, TK);
        }

        Z_Peneloux = *z - Peneloux_yes_no * cm * PBar / (RBar * TK); // Peneloux 修正

        // 计算各组分 ln(phi)
        for (comp = 0; comp < NumGases; comp++)
        {
            if (xGlobal[comp] > 0)
            {
                asum = 0.0;
                for (k = 0; k < NumGases; k++)
                {
                    if (k > comp)
                    {
                        asum += x[k] * sqrt(aPR[comp] * aPR[k]) * (1.0 - kPr[k][comp]);
                    }
                    else
                    {
                        asum += x[k] * sqrt(aPR[comp] * aPR[k]) * (1.0 - kPr[comp][k]);
                    }
                }

                ai_bar = 2.0 * asum - aPrmix;
                qi_bar = QQ * (1.0 + ai_bar / aPrmix - bPR[comp] / bm);
                ii = 1.0 / (sigma - epsilon) * log((*z + sigma * BStarPr) / (*z + epsilon * BStarPr));

                lnPHI[comp] = bPR[comp] / bm * (*z - 1.0) - log(*z - BStarPr) - qi_bar * ii;
                lnPHI_Peneloux[comp] = lnPHI[comp] - Peneloux_yes_no * (c0[comp] + c1[comp] * (TK - 288.15)) * PBar / RBar / TK;
                lnPHI[comp] = lnPHI_Peneloux[comp];
            }
        }
        // 修正水相逸度系数   -这里暂时没有进去
        if (phase == "liquid" && isAqueous_2016(eqAqueous, x, xGlobal, NumGases))
        {
            if (xGlobal[0] > 0)
                lnPHI[0] = logKHenry_CH4(TK) - log(PBar) + (34.5) * (PBar - 1) / (RBar * TK);

            if (xGlobal[1] > 0) {
                lnPHI[1] = log(gNeut[0] * (1000.0 / (18.015 * x[NumGases - 1])) / (PBar * KH_CO2(TK, PBar)));
            }

            if (xGlobal[2] > 0) {
                lnPHI[2] = log(gNeut[1] * (1000.0 / (18.015 * x[NumGases - 1])) / (PBar * KH_H2S(TK, PBar)));
            }

            if (xGlobal[3] > 0) {
                lnPHI[3] = logKHenry_C2H6(TK) - log(PBar) + (52.9) * (PBar - 1) / (RBar * TK);
            }

            if (xGlobal[NumGases - 1] > 0) {
                lnPHI[NumGases - 1] = logKHenry_N2(TK) - log(PBar) + (35.7) * (PBar - 1) / (RBar * TK);
            }

            double Psat = PsatH2O(TK);
            double qW = aPR[NumGases - 1] / (bPR[NumGases - 1] * RBar * TK);
            double BStarW = bPR[NumGases - 1] * Psat / (RBar * TK);
            double ZW, iW, lnPHIWsat;

            a3 = 1.0;
            a2 = (epsilonP + sigmaP) * BStarW - 1.0 - S * BStarW;
            a1 = (epsilonP * sigmaP) * BStarW * BStarW - (epsilonP + sigmaP) * BStarW - (epsilonP + sigmaP) * S * BStarW * BStarW + qW * BStarW;
            a0 = -sigmaP * epsilonP * BStarW * BStarW - sigmaP * epsilonP * BStarW * BStarW * BStarW * S - qW * S * BStarW * BStarW;

            cubic(a3, a2, a1, a0, "min", &ZW, bPR[NumGases - 1], Psat, TK);
            iW = log((ZW + sigmaP * BStarW) / (ZW + epsilonP * BStarW)) / (sigmaP - epsilonP);
            lnPHIWsat = ZW - 1.0 - log(ZW - BStarW) - qW * iW;

            lnPHI[NumGases - 1] = log(aH2O * Psat / (x[NumGases - 1] * PBar)) + lnPHIWsat + ZW * (1.0 - Psat / PBar);
        }
        *z = Z_Peneloux;
    }
    else
    {
        for (i = 0; i < NumGases; i++)
            lnPHI[i] = 0;
    }
}



IntVec BubbleSort(int Ascend, DoubleVec& ArrayIn, IntVec& index, int n)
{
    double SrtTemp;
    int indSrtTemp;
    int i, j;

    if (Ascend == 1)
    {
        for (i = 0; i < n; i++)
        {
            for (j = i + 1; j < n; j++)
            {
                if (ArrayIn[i] > ArrayIn[j])
                {
                    SrtTemp = ArrayIn[j];
                    ArrayIn[j] = ArrayIn[i];
                    ArrayIn[i] = SrtTemp;

                    indSrtTemp = index[j];
                    index[j] = index[i];
                    index[i] = indSrtTemp;
                }
            }
        }
    }
    else
    {
        for (i = 0; i < n; i++)
        {
            for (j = i + 1; j < n; j++)
            {
                if (ArrayIn[i] < ArrayIn[j])
                {
                    SrtTemp = ArrayIn[j];
                    ArrayIn[j] = ArrayIn[i];
                    ArrayIn[i] = SrtTemp;

                    indSrtTemp = index[j];
                    index[j] = index[i];
                    index[i] = indSrtTemp;
                }
            }
        }
    }
    return index;
}

//-------------------------------------------------


double weighted_mean(DoubleVec& weight, DoubleVec& property, int n)
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


void OrderPhases(DoubleMatrix& compositions, DoubleMatrix& phi,
    DoubleVec& beta, DoubleVec& MW_Phase, DoubleVec& mass_phase,
    DoubleVec& density, DoubleVec& Compr, int& nComp, int& nPhase)
{
    //int nC = *nComp;
    //int nP = *nPhase;
    int i, k, L;

    DoubleMatrix compositions0(nComp, DoubleVec(nPhase + 1, 0));
    DoubleMatrix phi0(nComp, DoubleVec(nPhase, 0));

    DoubleVec BETA0(nPhase, 0);
    DoubleVec MW_Phase0(nPhase, 0);
    DoubleVec mass_phase0(nPhase, 0);
    DoubleVec density0(nPhase, 0);
    DoubleVec Compr0(nPhase, 0);

    //按密度重新排序相
    int ascending = 1;
    IntVec idx_dens(nPhase, 0);
    for (k = 0; k < nPhase; k++)
    {
        density0[k] = density[k];
        idx_dens[k] = k;
    }
    BubbleSort(ascending, density0, idx_dens, nPhase);

    //将所有变量分配给临时数组，按密度排序
    for (k = 0; k < nPhase; k++)
    {
        L = idx_dens[k];
        for (i = 0; i < nComp; i++)
        {
            /*
            原代码：
               compositions0(i, k + 1) = compositions(i, L + 1)
               phi0(i, k) = phi(i, L)
            */
            compositions0[i][k + 1] = compositions[i][L + 1];
            phi0[i][k] = phi[i][L];
        }
        BETA0[k] = beta[L];
        MW_Phase0[k] = MW_Phase[L];
        mass_phase0[k] = mass_phase[L];
        density0[k] = density[L];
        Compr0[k] = Compr[L];
    }

    //确定平衡状态下存在的相类型
    //这里将其全部初始化为-1(在vb中是0，即数组最小索引的前一位)
    int gas = -1, oil = -1, aqueous = -1, zeroPhase = -1, HC1 = -1, HC2 = -1, HCPhases = 0;
    for (k = 0; k < nPhase; k++) {
        if (density0[k] >= 0.3 && compositions0[nComp - 1][k + 1] > 0.5)
            aqueous = k;
        else if (density0[k] == 0)  // 判断double类型是否为0一般不这样写，此处对应原Vb代码。
            zeroPhase = k;
        else {
            if (HC1 == -1) HC1 = k;
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
    if (aqueous != -1 && nPhase != 3) {
        /*
        vb原码：nPhase = 3 ,redim compositions(nComp, nPhase + 1)等直接修改了参数列表。
        推测: 使用引用传递
        */
        nPhase = 3;
        for (auto& row : compositions)
            row.resize(nPhase);
        for (auto& row : phi)
            row.resize(nPhase);
        beta.resize(nPhase, 0);
        MW_Phase.resize(nPhase, 0);
        mass_phase.resize(nPhase, 0);
        density.resize(nPhase, 0);
        Compr.resize(nPhase, 0);
    }

    // output phases
    L = 0;
    for (k = 0; k < nPhase; k++)
    {
        if (k == gas) L = 1;//gas由HC赋值，HC由k赋值，因此直接照搬原代码
        else if (k == oil) L = 2;
        else if (k == aqueous) L = 3;

        if (L != 0)
        {
            if (BETA0[k] > 0)
            {
                //注：原代码中是beta[L] = BETA0[k]，那么存在beta的第1、2、3个元素存在被修改的情况，C中的第1、2、3个元素，对应的索引为0、1、2，此处注意偏移量。
                beta[L - 1] = BETA0[k];
                MW_Phase[L - 1] = MW_Phase0[k];
                mass_phase[L - 1] = mass_phase0[k];
                density[L - 1] = density0[k];
                Compr[L - 1] = Compr0[k];
                for (i = 0; i < nComp; i++)
                {
                    compositions[i][L] = compositions0[i][k + 1];
                    phi[i][L - 1] = phi0[i][k];
                }
            }
            else {
                beta[L - 1] = 0;
                MW_Phase[L - 1] = 0;
                mass_phase[L - 1] = 0;
                density[L - 1] = 0;
                Compr[L - 1] = 0;
                for (i = 0; i < nComp; i++)
                {
                    compositions[i][L] = 0;
                    phi[i][L - 1] = 0;
                }
            }
        }
    }

    // 从输出中删除任何不存在的阶段。 注释：本质上第一个修改beat[1],对应C语言的就是beta[0].，可以省略L，直接写数字即可，未修改是为了对应vb代码。
    L = 1;
    if (gas == -1)
    {
        beta[L - 1] = 0;
        MW_Phase[L - 1] = 0;
        mass_phase[L - 1] = 0;
        density[L - 1] = 0;
        Compr[L - 1] = 0;
        for (i = 0; i < nComp; i++)
        {
            compositions[i][L] = 0;
            phi[i][L - 1] = 0;
        }
    }
    L = 2;
    if (oil == -1)
    {
        beta[L - 1] = 0;
        MW_Phase[L - 1] = 0;
        mass_phase[L - 1] = 0;
        density[L - 1] = 0;
        Compr[L - 1] = 0;
        for (i = 0; i < nComp; i++)
        {
            compositions[i][L] = 0;
            phi[i][L - 1] = 0;
        }
    }
    L = 3;
    if (aqueous == -1 && nPhase == 3)
    {
        beta[L - 1] = 0;
        MW_Phase[L - 1] = 0;
        mass_phase[L - 1] = 0;
        density[L - 1] = 0;
        Compr[L - 1] = 0;
        for (i = 0; i < nComp; i++)
        {
            compositions[i][L] = 0;
            phi[i][L - 1] = 0;
        }
    }
}


int myLimitIndex(const char* max_or_min, DoubleVec& composition, DoubleVec& property, int k)
{
    int i;
    double lastProperty;
    double tol = 1e-6;
    int index = -1; /* 初始化返回值，未找到有效元素时返回 -1 */

    if (strcmp(max_or_min, "max") == 0)
    {
        lastProperty = -DBL_MAX; /* 非常小的初始值 */
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
        lastProperty = DBL_MAX; /* 非常大的初始值 */
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

void InitialBeta(int eqVapor, int eqAqueous, double TK, double PBar, int NonZeroNumGases,
    DoubleVec& z, DoubleVec& TCr, DoubleVec& PCr, DoubleVec& MWgas,
    DoubleVec& beta, int NumGases, int MaxBeta, BoolVec& betaexists)
{
    int m, j;
    double sumBeta = 0;

    /* 初始猜测 */
    if (MaxBeta >= 2)
    {
        beta[0] = 0.5;
        beta[1] = 0.5;
    }

    /* 检查水是否存在 (最后一个组分) */
    if (z[NumGases - 1] > 0.0)
    {
        if (TK < TCr[NumGases - 1] && (PBar > PsatH2O(TK)))
        { /* 调用外部函数 PsatH2O */
            if (MaxBeta > 2)
            {
                beta[0] = (1.0 - z[NumGases - 1]) / 2.0;
                beta[1] = (1.0 - z[NumGases - 1]) / 2.0;
                beta[MaxBeta - 1] = z[NumGases - 1];
            }
            else
            {
                beta[0] = 1.0 - z[NumGases - 1];
                beta[MaxBeta - 1] = z[NumGases - 1];
            }
        }
        else
        {
            if (MaxBeta > 2)
            {
                beta[0] = (1.0 - z[NumGases - 1]) / 2.0 + z[NumGases - 1];
                beta[1] = (1.0 - z[NumGases - 1]) / 2.0;
                beta[MaxBeta - 1] = 0.0;
            }
            else
            {
                beta[0] = 1.0;
                beta[MaxBeta - 1] = 0.0;
            }
        }

        if (NonZeroNumGases == 2)
        {
            /* VB 中 ReDim beta(MaxBeta) 等价于重置，此处仅重新赋值 */
            if (MaxBeta >= 2)
            {
                beta[0] = 0.5;
                beta[MaxBeta - 1] = 0.5;
            }
        }
    }
    for (j = 0; j < MaxBeta; j++)
        sumBeta += beta[j];
    for (j = 0; j < MaxBeta; j++)
        beta[j] /= sumBeta;
    for (j = 0; j < MaxBeta; j++)
        if (beta[j] > 0)
            betaexists[j] = true;
}


void Initialization(double zGlobalWater, int NonZeroNumGases, int MaxBeta, const char* EOS,
    DoubleVec& MWgas, DoubleVec& TCr, DoubleVec& PCr, DoubleVec& Omega,
    double TK, double PBar, DoubleMatrix& logphi3phase, int NumGases)
{
    int i, j;
    double correction;
    DoubleVec lnKi(NumGases, 0);

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

        /*
         * SSP special corrections (Henry's law) are commented out
         * Example:
         * logphi3phase[0][MaxBeta-1] = logKHenry_CH4(TK) - log(PBar);
         */
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
}



void EQCalculation(DoubleVec& beta, DoubleMatrix& lnPHI,
    DoubleVec& zGlobal, DoubleVec& E, double* Q, int NumGases, int MaxBeta)
{
    int i, k;
    double Qphase = 0.0;
    double Qcomp = 0.0;

    E.resize(NumGases, 0);


    /* 计算辅助向量 E[i] */
    for (i = 0; i < NumGases; ++i)
    {
        E[i] = 0.0;
        for (k = 0; k < MaxBeta; ++k)
        {
            if (beta[k] > 0.0 && zGlobal[i] > 0.0)
                E[i] += beta[k] / exp(lnPHI[i][k]); /* lnPHI[i][k] -> e^(-lnPHI) */
        }
    }

    /* 计算 Qphase */
    for (k = 0; k < MaxBeta; ++k)
        Qphase += beta[k];

    /* 计算 Qcomp */
    for (i = 0; i < NumGases; ++i)
        if (zGlobal[i] > 0.0)
            Qcomp += zGlobal[i] * log(E[i]);

    *Q = Qphase - Qcomp;
}


void Gauss(DoubleMatrix& a, DoubleVec& b, int n, DoubleVec& xs)
{
    int i, j, k, kMax;
    DoubleMatrix AB(n, DoubleVec(n + 1, 0));
    DoubleVec rowmax(n + 1, 0);
    double max, sum;

    // 动态分配增广矩阵 AB[n][n+1]
    //AB = (double**)malloc(n * sizeof(double*));
    //for (i = 0; i < n; i++)
    //    AB[i] = (double*)malloc((n + 1) * sizeof(double));

    //rowmax = (double*)malloc((n + 1) * sizeof(double));

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
        kMax = 1;

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
    for (j = n - 2; j >= 0; j--)
    {
        sum = 0.0;
        for (k = j + 1; k < n; k++)
        {
            sum += AB[j][k] * xs[k];
        }
        xs[j] = (AB[j][n] - sum) / AB[j][j];
    }

}


double myMin(DoubleVec& property, int n)
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


void Equilibrium(DoubleVec& zGlobal, DoubleMatrix& logphi3phase, DoubleVec& E, double* Q,
    DoubleVec& beta, BoolVec& betaexists, DoubleMatrix& compositions,
    bool* continue_flash, int* counterEquilibrium_final, int* iter_final, int* counter_final, int NumGases, int MaxBeta)
{
    bool converged = false;
    long iter = 0, itermax = 1000, counter = 0, counterEquilibrium = 0;

    DoubleVec g(MaxBeta, 0);
    DoubleVec gcopy(MaxBeta, 0);
    DoubleVec gsolve(MaxBeta, 0);
    DoubleVec alpha0(MaxBeta, 0);
    DoubleVec betanew(MaxBeta, 0);
    DoubleVec w(MaxBeta, 0);
    DoubleVec Enew(NumGases, 0);
    DoubleMatrix Hessian(MaxBeta, DoubleVec(MaxBeta, 0));


    // 错误处理标签
    int error_code = 0;

    double alpha0used, Er2, tol, Difference, epsilon, Qnew;
    //*continue_flash = true;

    while (!converged && counterEquilibrium < 100)
    {
        Er2 = 1.0;
        tol = 1e-12;
        iter = 0;
        itermax = 1000;
        counterEquilibrium++;

        if (counterEquilibrium > 2000) {
            error_code = 3;  // 迭代超时
            goto ErrHandler;
        }

        while (Er2 > tol && iter < itermax)
        {
            iter++;
            if (iter > itermax) {
                error_code = 4;  // 内部迭代超时
                goto ErrHandler;
            }
            if (counter > 10000) {
                error_code = 5;  // 计数器超限
                goto ErrHandler;
            }

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

            while (Difference > epsilon && counter < 10000)
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
            }

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

        }

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

    }

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

    return;
ErrHandler:
    *continue_flash = false;
    *counterEquilibrium_final = counterEquilibrium;
    *iter_final = iter;
    *counter_final = counter;

    if (error_code != 0) {
        *continue_flash = false;

        // 如果还没有计算结果，尽量计算compositions
        if (error_code != 1) {  // 如果不是内存分配失败
            for (int i = 0; i < NumGases; ++i) {
                compositions[i][0] = zGlobal[i];
                for (int k = 0; k < MaxBeta; ++k) {
                    if (zGlobal[i] > 0.0) {
                        compositions[i][k + 1] = zGlobal[i] / E[i] / exp(logphi3phase[i][k]);
                    }
                }
            }
        }

        *counterEquilibrium_final = counterEquilibrium;
        *iter_final = iter;
        *counter_final = counter;
    }
    return;
}



int isVapor(int eqVapor, DoubleVec& composition, DoubleVec& globalComposition, double TK, double PBar,
    DoubleVec& TCr, DoubleVec& PCr, int lightest, int NumGases)
{
    double composition_ratio = 1.0; // 组分比率阈值
    double tol = 0.00001;           // 容差
    double total_dry = 0.0;
    double total_comp_dry = 0.0;
    double dry_composition_Lightest = 0.0;
    double dry_global_comp_Lightest = 0.0;

    // 计算除水或最轻组分外的干物质总量
    for (int i = 0; i < NumGases - 1; i++)
    {
        total_dry += globalComposition[i];
        total_comp_dry += composition[i];
    }

    dry_global_comp_Lightest = globalComposition[lightest] / total_dry;
    dry_composition_Lightest = composition[lightest] / total_comp_dry;

    if ((dry_composition_Lightest - composition_ratio * dry_global_comp_Lightest) > 0.0)
    {
        return 1; // 气相
    }
    else if (eqVapor && fabs(dry_composition_Lightest - composition_ratio * dry_global_comp_Lightest) < tol)
    {
        return 1; // 气相
    }

    return 0; // 非气相
}


void normalize2Darray(DoubleMatrix& x, int nRows, int nCols)
{
    int i, j;
    DoubleVec sum(nCols, 0);

    // 计算每列的和
    for (j = 0; j < nCols; j++)
        for (i = 0; i < nRows; i++)
            sum[j] += x[i][j];

    // 对每列进行归一化
    for (j = 0; j < nCols; j++)
        if (sum[j] > 0.0)
            for (i = 0; i < nRows; i++)
                x[i][j] /= sum[j];
}


bool equalArrays(double tolerance, DoubleVec& array1, DoubleVec& array2, int nElements)
{
    int i;
    double D1, D2;
    double sumError = 0.0;
    int sumResultArray = 0;

    for (i = 0; i < nElements; i++)
    {
        D1 = array1[i];
        D2 = array2[i];

        if (D1 > D2)
        {
            if (fabs(D1 - D2) / fabs(D1) < tolerance)
                sumResultArray += 0;
            else
                sumResultArray += 1;
            sumError += fabs(D1 - D2) / fabs(D1);
        }
        else if (D1 < D2)
        {
            if (fabs(D1 - D2) / fabs(D2) < tolerance)
                sumResultArray += 0;
            else
                sumResultArray += 1;
            sumError += fabs(D1 - D2) / fabs(D2);
        }
        else
            sumResultArray += 0;
    }

    if (sumResultArray == 0 || (sumError / nElements) < tolerance)
        return true;
    else
        return false;
}


double myMax(DoubleVec& property, int nElements)
{
    int i;
    double maxVal = property[0]; // C语言数组从0开始
    for (i = 1; i < nElements; i++)
    {
        if (property[i] > maxVal)
            maxVal = property[i];
    }
    return maxVal;
}


void normalizeWithoutWater(DoubleMatrix& x, int nRows, int nCols)
{
    int i, j;
    DoubleVec sum(nCols, 0);

    // 计算每列（不包括最后一行）的和
    for (j = 0; j < nCols; j++)
    {
        for (i = 0; i < nRows - 1; i++)
            sum[j] += x[i][j];
    }

    // 对每列（不包括最后一行）进行归一化
    for (j = 0; j < nCols; j++)
    {
        if (sum[j] > 0.0)
        {
            for (i = 0; i < nRows - 1; i++)
                x[i][j] /= sum[j];
        }
    }
}


void indexBubbleSrt(bool ascending, DoubleVec& ArrayIn, IntVec& index, int nElements)
{
    int i, j;
    double tempVal;
    int tempIdx;

    if (ascending) {
        for (i = 0; i < nElements - 1; i++) {
            for (j = i + 1; j < nElements; j++) {
                if (ArrayIn[i] > ArrayIn[j]) {
                    // 交换 ArrayIn 元素
                    tempVal = ArrayIn[j];
                    ArrayIn[j] = ArrayIn[i];
                    ArrayIn[i] = tempVal;

                    // 交换对应索引
                    tempIdx = index[j];
                    index[j] = index[i];
                    index[i] = tempIdx;
                }
            }
        }
    }
    else {
        for (i = 0; i < nElements - 1; i++) {
            for (j = i + 1; j < nElements; j++) {
                if (ArrayIn[i] < ArrayIn[j]) {
                    // 交换 ArrayIn 元素
                    tempVal = ArrayIn[j];
                    ArrayIn[j] = ArrayIn[i];
                    ArrayIn[i] = tempVal;

                    // 交换对应索引
                    tempIdx = index[j];
                    index[j] = index[i];
                    index[i] = tempIdx;
                }
            }
        }
    }
}


void Fix_Position_Phases(bool AuroraCalculation, int lightest, int heaviest, bool eqAqueous, int NonZeroComponents,
    DoubleVec& beta, DoubleVec& density, DoubleVec& Compr, DoubleMatrix& compositions,
    DoubleMatrix& logPHI, StrVec& phase, int* nPhases, StrVec& phaseName,
    int maxPhases, int nComponents)
{
    int i, j, k, m, n;
    int aqueousPhase = -1;
    double composition_ratio = 1.001;
    double tol = composition_ratio - 1.0;

    std::vector<BoolVec> equalPhases(maxPhases, BoolVec(maxPhases, false));

    DoubleVec relativeComposition_Water(maxPhases, 0);
    DoubleVec compoPhaseA(nComponents, 0);
    DoubleVec compoPhaseB(nComponents, 0);

    BoolVec uniquePhase(maxPhases, false);

    IntVec phaseIndex(maxPhases, 0);
    IntVec phaseOrder(maxPhases, 0);
    DoubleVec variable_sorted(maxPhases, 0);

    double sumBeta = 0.0;
    bool vaporExists = false;
    bool aqueousExists = false;
    int nonAqueousPhases;

    DoubleVec beta_prime(maxPhases, 0);
    DoubleVec density_prime(maxPhases, 0);
    DoubleVec Compr_prime(maxPhases, 0);

    StrVec phase_prime(maxPhases, "");

    //compo_prime(nComponents, maxPhases + 1), logPHI_prime(nComponents, maxPhases)
    DoubleMatrix compo_prime(nComponents, DoubleVec(maxPhases + 1, 0));
    DoubleMatrix logPHI_prime(nComponents, DoubleVec(maxPhases, 0));

    for (int k = 0; k < maxPhases; k++)
    {
        // C 语言中，VB的 beta(k) 翻译为 beta[k]
        if (beta[k] <= 0)
        {
            beta[k] = 0;
            density[k] = 0;
            Compr[k] = 0;

            // VB的 For i = 1 To nComponents 翻译为 i=0 To nComponents-1
            for (int i = 0; i < nComponents; i++)
            {
                // C 语言中，VB的 compositions(i, k + 1) 翻译为 compositions[i][k + 1]
                compositions[i][k + 1] = 0;
                // C 语言中，VB的 logPHI(i, k) 翻译为 logPHI[i][k]
                logPHI[i][k] = 0;
            }
        }
    }

    // 复制 compositions 到 relativeCompositions 并归一化
    DoubleMatrix relativeCompositions(nComponents, DoubleVec(maxPhases + 1, 0));
    for (i = 0; i < nComponents; i++)
        for (j = 0; j <= maxPhases; j++)
            relativeCompositions[i][j] = compositions[i][j];

    normalize2Darray(relativeCompositions, nComponents, maxPhases + 1);

    for (i = 0; i < maxPhases; i++)
        for (j = 0; j < maxPhases; j++)
            equalPhases[i][j] = false;
    // 判断相是否相等
    double tolerance = 0.03; // 3%
    int totalEqualPhases = 0;
    for (m = 0; m < maxPhases; m++)
    {
        for (n = 0; n < m; n++)
        {
            for (i = 0; i < nComponents; i++)
            {
                compoPhaseA[i] = relativeCompositions[i][m + 1];
                compoPhaseB[i] = relativeCompositions[i][n + 1];
            }
            if (equalArrays(tolerance, compoPhaseA, compoPhaseB, nComponents))
            {
                equalPhases[m][n] = true;
                equalPhases[n][m] = true;
                if (m != n)
                    totalEqualPhases++;
            }
        }
        equalPhases[m][m] = true;
    }

    // 找水相
    if (compositions[nComponents - 1][0] > 0)
    {
        for (j = 0; j < maxPhases; j++)
            relativeComposition_Water[j] = relativeCompositions[nComponents - 1][j + 1];

        int maxWater_Phase = myLimitIndex("max", relativeComposition_Water, relativeComposition_Water, maxPhases);
        if (maxWater_Phase >= 0)
        {
            if (relativeComposition_Water[maxWater_Phase] > 0.5 && density[maxWater_Phase] > 0.5)
            {
                if (relativeComposition_Water[maxWater_Phase] - composition_ratio * relativeCompositions[nComponents - 1][0] > 0)
                    aqueousPhase = maxWater_Phase;
                else if (relativeComposition_Water[maxWater_Phase] > 0.9)
                    aqueousPhase = maxWater_Phase;
                else if (eqAqueous && fabs(myMax(relativeComposition_Water, maxPhases) - relativeCompositions[nComponents - 1][0]) <= tol * relativeCompositions[nComponents - 1][0])
                    aqueousPhase = maxWater_Phase;
            }
        }
        // 除水外归一化
        normalizeWithoutWater(relativeCompositions, nComponents, maxPhases + 1);
    }

    // 判断非零组分数
    NonZeroComponents = 0;
    for (i = 0; i < nComponents; i++)
        if (compositions[i][0] > 0)
            NonZeroComponents++;

    for (j = 0; j < maxPhases; j++)
    {
        if (aqueousPhase > 0 && !aqueousExists)
            if (beta[j] > 0 && equalPhases[j][aqueousPhase])
            {
                aqueousExists = true;
                aqueousPhase = j;
            }
    }

    // 唯一相检测
    for (n = 0; n < maxPhases; n++)
    {
        bool found = false;
        int m_min = -1;
        for (m = 0; m <= n; m++)
        {
            if (equalPhases[n][m] && beta[m] > 0)
            {
                if (!found)
                {
                    uniquePhase[m] = true;
                    found = true;
                    m_min = m;
                }
                else
                {
                    beta[m_min] += beta[m];
                }
            }
        }
    }

    *nPhases = 0;
    //Finds the existing unique phases, creates a copy the current properties values
    for (j = 0; j < maxPhases; j++)
    {
        /* This loop nullifies fugacity coefficients for non-existing components */
        for (i = 0; i < nComponents; i++)
        {
            if (compositions[i][0] == 0.0)
            {
                logPHI[i][j] = 0.0;
                compositions[i][j + 1] = 0.0;
            }
        }

        if (uniquePhase[j] == true)
        {
            (*nPhases)++;
            beta_prime[j] = beta[j];
            density_prime[j] = density[j];
            Compr_prime[j] = Compr[j];
            phase_prime[j] = phase[j];
            //strcpy(phase_prime[j], phase[j]);

            for (i = 0; i < nComponents; i++)
            {
                compo_prime[i][j + 1] = compositions[i][j + 1];
                compo_prime[i][0] = compositions[i][0];
                logPHI_prime[i][j] = logPHI[i][j];
            }
        }
    }
    //relativeCompositionsPrime = compo_prime

    DoubleMatrix relativeCompositionsPrime(nComponents, DoubleVec(maxPhases + 1, 0));
    for (i = 0; i < nComponents; i++)
        for (j = 0; j < maxPhases + 1; j++)
            relativeCompositionsPrime[i][j] = compo_prime[i][j];
    normalize2Darray(relativeCompositionsPrime, nComponents, maxPhases + 1);
    //rel_Compo_Lightest(maxPhases), rel_Compo_Heaviest(maxPhases)
    DoubleVec rel_Compo_Lightest(maxPhases, 0);
    DoubleVec rel_Compo_Heaviest(maxPhases, 0);
    for (j = 0; j < maxPhases; j++)
    {
        rel_Compo_Lightest[j] = relativeCompositionsPrime[lightest][j + 1];
        rel_Compo_Heaviest[j] = relativeCompositionsPrime[heaviest][j + 1];
    }

    // 不明意义的Redim行为，暂忽略
    //ReDim beta(maxPhases), density(maxPhases), Compr(maxPhases)
    //,compositions(nComponents, maxPhases + 1), logPHI(nComponents, maxPhases)

    //Finds the proper order of the existing phases in a sorted desired property value
    for (j = 0; j < maxPhases; j++)
        phaseIndex[j] = j;
    bool ascending = false;
    for (i = 0; i < maxPhases; i++)
        variable_sorted[i] = rel_Compo_Lightest[i];

    indexBubbleSrt(ascending, variable_sorted, phaseIndex, maxPhases);
    for (j = 0; j < maxPhases; j++)
        phaseOrder[j] = phaseIndex[j];

    //nonAqueousPhases = *nPhases - abs(CLng(aqueousExists))
    if (aqueousExists)
        nonAqueousPhases = *nPhases - 1;
    else
        nonAqueousPhases = *nPhases;

    if (nonAqueousPhases == 1 && density_prime[phaseOrder[0]] < 0.2)
        vaporExists = true;
    else if (*nPhases == 1 && aqueousExists)
    {
        /* Do nothing (same as VB's empty ElseIf) */
    }
    else
    {
        if (phase_prime[phaseOrder[0]] == "vapor" &&
            density_prime[phaseOrder[0]] < 0.2)
        {
            vaporExists = true;
        }
        else if ((rel_Compo_Lightest[phaseOrder[0]] *
            rel_Compo_Heaviest[phaseOrder[1]]) >
            1000.0 *
            (rel_Compo_Lightest[phaseOrder[1]] *
                rel_Compo_Heaviest[phaseOrder[0]]) &&
            phaseOrder[0] != aqueousPhase)
        {
            vaporExists = true;
        }
        else if (2.0 * density_prime[phaseOrder[0]] < density_prime[phaseOrder[1]] &&
            phaseOrder[1] != aqueousPhase)
        {
            vaporExists = true;
        }
        else
        {
            vaporExists = false;
        }
    }
    //z(nComponents)
    DoubleVec z(nComponents);
    for (i = 0; i < nComponents; i++)
        z[i] = compo_prime[i][0];

    k = 1;
    m = 0;

    DoubleVec x(nComponents);
    for (j = 0; j < maxPhases; j++)
    {
        if (beta_prime[phaseOrder[j]] > 0.0)
        {
            /* If k == 1 And Not vaporExists And nPhases < maxPhases */
            if (k == 1 && !vaporExists && (*nPhases < maxPhases))
            {
                k = 2;
            }

            m = k;

            for (i = 0; i < nComponents; i++)
                x[i] = compo_prime[i][phaseOrder[j] + 1];

            if ((density_prime[phaseOrder[j]] > 0.5 &&
                aqueousExists &&
                isAqueous_2016(true, x, z, nComponents))
                ||
                (NonZeroComponents == 2 &&
                    x[nComponents - 1] > 0.5 &&
                    phase[j] == "liquid"))
            {
                k = maxPhases;
                m = m - 1;
            }

            beta[k - 1] = beta_prime[phaseOrder[j]];
            density[k - 1] = density_prime[phaseOrder[j]];
            Compr[k - 1] = Compr_prime[phaseOrder[j]];
            phase[k - 1] = phase_prime[phaseOrder[j]];
            //strcpy(phase[k - 1], phase_prime[phaseOrder[j]]);

            for (i = 0; i < nComponents; i++)
            {
                compositions[i][k] = compo_prime[i][phaseOrder[j] + 1];
                logPHI[i][k - 1] = logPHI_prime[i][phaseOrder[j]];
            }

            k = m + 1;
        }
    }

    /* Normalizing phase fractions */
    sumBeta = 0.0;
    for (j = 0; j < maxPhases; j++)
        sumBeta += beta[j];
    for (j = 0; j < maxPhases; j++)
        beta[j] = beta[j] / sumBeta;


    /* Naming phases (Vapor, Liquid, Aqueous) */

    m = 1;

    if (*nPhases == 1)
    {
        if (vaporExists)
            phaseName[0] = "Steam";
        else if (!aqueousExists)
            phaseName[1] = "Liquid";
        else
            phaseName[maxPhases - 1] = "Aqueous";
    }
    else
    {
        for (j = 0; j < maxPhases; j++)
        {
            if (beta[j] > 0.0)
            {
                for (i = 0; i < nComponents; i++)
                    x[i] = compositions[i][j + 1];

                if (aqueousExists && isAqueous_2016(eqAqueous, x, z, nComponents))
                    phaseName[j] = "Aqueous";
                else
                {
                    /* VB: "Phase " & CStr(m)  (C 要自己格式化字符串) */
                    phaseName[j] = "Phase " + std::to_string(m);
                    m = m + 1;
                }
            }
        }
    }
}


void Flash(bool* eqVapor, bool* eqAqueous, const char* EOS, double TK, double PBar, int NonZeroNumGases, DoubleVec& zGlobal,
    DoubleVec& gNeut, double aH2O, DoubleVec& TCr, DoubleVec& PCr, DoubleVec& Omega,
    DoubleVec& MWgas, DoubleMatrix& kPr, DoubleVec& c0, DoubleVec& c1,
    double* Q, int* Numbeta, StrVec& phaseName, DoubleVec& beta, DoubleMatrix& compositions,
    DoubleMatrix& phi, DoubleVec& Compr, DoubleVec& density,
    int* counterEquilibrium_final, int* iter_final, int* counter_final, DoubleVec& zOutput,
    int NumGases, int MaxBeta)
{

    int i, j, k;
    int iter = 0, itermax = 30;
    int lightest, heaviest;
    bool continue_flash = true;
    bool vaporFound;
    bool highErr = true;
    double pseudoTCr_mix, pseudoPCr_mix;
    double sumorganic, Er, tol, MWsum;
    double Sz, errNet, Er_old, errRate;
    double sumOld, sumNew;
    double maxErrorLogPhi;
    double Gibbs;
    double temp;
    bool testTemp = false;

    const double RBar = 83.14;

    // 获取组分和相数
    //NumGases = NonZeroNumGases; // UBound(zGlobal) → 非零组分数
    //MaxBeta = *Numbeta;         // UBound(beta)

    // 分配主要数组
    DoubleVec Ki(NumGases, 0);

    DoubleMatrix logphi3phase(NumGases, DoubleVec(MaxBeta, 0));

    DoubleMatrix logphi3phase_new(NumGases, DoubleVec(MaxBeta, 0));

    DoubleVec logphi3phase_newphase(NumGases, 0);

    DoubleVec x(NumGases, 0);
    DoubleVec z(NumGases, 0);

    DoubleVec everyLogPHI(NumGases, 0);
    DoubleMatrix errorLogPhi(NumGases, DoubleVec(MaxBeta, 0));

    BoolVec betaexists(MaxBeta, false);

    StrVec phase(MaxBeta, "");

    // 归一化全局摩尔组成
    Sz = 0.0;
    for (i = 0; i < NumGases; i++) Sz += zGlobal[i];
    for (i = 0; i < NumGases; i++) {
        z[i] = zGlobal[i] / Sz;
        zOutput[i] = z[i];
    }

    // 找出最轻/最重组分 (排除水)
    lightest = myLimitIndex("min", z, TCr, NumGases - 1);
    heaviest = myLimitIndex("max", z, TCr, NumGases - 1);

    // 初始相分数估计
    InitialBeta(*eqVapor, *eqAqueous, TK, PBar, NonZeroNumGases, z, TCr, PCr, MWgas, beta, NumGases, MaxBeta, betaexists);

    // 初始 logphi 计算
    Initialization(zGlobal[NumGases - 1], NonZeroNumGases, MaxBeta, EOS, MWgas, TCr, PCr, Omega, TK, PBar, logphi3phase, NumGases);

    // 避免数值不稳定：把零组分的 logphi 设置为 0
    for (j = 0; j < MaxBeta; j++) {
        for (i = 0; i < NumGases; i++) {
            if (zGlobal[i] == 0.0) {
                logphi3phase[i][j] = 0.0;
            }
        }
    }

    // 迭代初值
    highErr = true;
    Er = 100.0;
    errNet = 100.0;
    Er_old = 1.0;
    errRate = 100.0;
    tol = 0.0002;

    iter = 0;
    continue_flash = true;

    if (PBar < 1.02) {
        testTemp = true;
    }

    // 主循环
    while ((errNet > tol && iter < 30) &&
        (errRate > 5 && highErr || Er > 200 * tol)) {
        iter++;

        EQCalculation(beta, logphi3phase, z, E, Q, NumGases, MaxBeta);

        Equilibrium(z, logphi3phase, E, Q, beta, betaexists, compositions,
            &continue_flash, counterEquilibrium_final, iter_final, counter_final, NumGases, MaxBeta);

        vaporFound = false;

        for (k = 0; k < MaxBeta; k++) {
            //明明只需要MaxBeta个，为什么要Redim成NumGases长度？？;
            everyLogPHI[k] = 0.0;
            for (i = 0; i < NumGases; i++) {
                x[i] = compositions[i][k + 1];
                everyLogPHI[k] += fabs(logphi3phase[i][k]);
            }

            if (everyLogPHI[k] < 0.1) {
                phase[k] = "vapor";
                vaporFound = true;
            }
            else if (isAqueous_2016(*eqAqueous, x, z, NumGases)) {
                phase[k] = "liquid";
            }
            else if (isVapor(*eqVapor, x, z, TK, PBar, TCr, PCr, lightest, NumGases)) {
                phase[k] = "vapor";
                vaporFound = true;
            }
            else if (!vaporFound && (NonZeroNumGases == 2) && (z[NumGases - 1] > 0) &&
                ((z[0] > 0 || z[1] > 0 || z[2] > 0 || z[3] > 0) ||
                    (z[NumGases - 2] > 0 && MWgas[NumGases - 2] < 29))) {
                phase[k] = "vapor";
                vaporFound = true;
            }
            else {
                phase[k] = "liquid";
            }

            phi_calc(*eqVapor, *eqAqueous, EOS, phase[k], TK, PBar, x, z, gNeut, aH2O,
                TCr, PCr, Omega, c0, c1, kPr, logphi3phase_newphase, &Compr[k], NumGases);

            for (i = 0; i < NumGases; i++) {
                logphi3phase_new[i][k] = logphi3phase_newphase[i];
            }
        }

        if (!continue_flash) {
            goto exit_flash;
        }

        errRate = 100.0 * fabs(Er - Er_old) / Er_old;
        if (Er > 0.2) highErr = true;
        Er_old = Er;
        Er = 0.0;

        sumOld = 0.0;
        sumNew = 0.0;
        errNet = 0.0;

        for (i = 0; i < NumGases; i++) {

            for (k = 0; k < MaxBeta; k++) {
                if (z[i] > 0.00001) {
                    Er += fabs(logphi3phase_new[i][k] - logphi3phase[i][k]);
                    sumNew += fabs(logphi3phase_new[i][k]);
                    sumOld += fabs(logphi3phase[i][k]);

                    if (fabs(logphi3phase_new[i][k]) > 0) {
                        errorLogPhi[i][k] = fabs((logphi3phase_new[i][k] - logphi3phase[i][k]) /
                            logphi3phase_new[i][k]);
                    }
                }

                logphi3phase[i][k] = logphi3phase_new[i][k];
                phi[i][k] = exp(logphi3phase[i][k]);
            }
        }

        errNet = fabs(sumOld - sumNew);
    }

exit_flash:
    // 计算各相分子量与密度
    for (j = 0; j < MaxBeta; j++) {
        if (beta[j] > 0) {
            MWsum = 0.0;
            for (i = 0; i < NumGases; i++) {
                MWsum += MWgas[i] * compositions[i][j + 1];
                //printf("\n %.16f \n", compositions[i][j + 1]);
            }
            if (Compr[j] > 0) {
                density[j] = (MWsum * PBar) / RBar / Compr[j] / TK;
            }
            else {
                density[j] = 0.0;
                beta[j] = 0.0;
            }
        }
        else {
            density[j] = 0.0;
            Compr[j] = 0.0;
        }
    }

    // 修正相的排列与赋值
    Fix_Position_Phases(true, lightest, heaviest, *eqAqueous,
        NonZeroNumGases, beta, density, Compr,
        compositions, phi, phase, Numbeta, phaseName, MaxBeta, NumGases);
}


//Compr(MaxBeta), beta(MaxBeta), density(MaxBeta), phi(NumGases, MaxBeta), compositions(NumGases, MaxBeta + 1), zOutput(NumGases)
void TrueFlash(const char* EOS, double TK, double PBar, int NonZeroNumGases, DoubleVec& zGlobal, DoubleVec& gNeut,
    double aH2O, DoubleVec& TCr, DoubleVec& PCr, DoubleVec& Omega, DoubleVec& MWgas,
    DoubleMatrix& kPr, DoubleVec& c0, DoubleVec& c1, int* Numbeta,
    StrVec& phaseName, DoubleVec& beta, DoubleMatrix& compositions, DoubleMatrix& phi,
    DoubleVec& Compr, DoubleVec& density, int* counterEquilibrium_final, int* iter_final, int* counter_final,
    DoubleVec& zOutput, int NumGases, int* MaxBeta)
{
    int i, j;
    //int NumGases, MaxBeta;
    double Q = 0.0, Q_alt = 0.0;

    // Alternative solution containers
    int temp_MaxBeta = *MaxBeta;
    int numBeta_alt = 0;
    //char** phaseName_alt = (char**)malloc(temp_MaxBeta * sizeof(char*));
    //for (i = 0; i < temp_MaxBeta; i++)
    //    phaseName_alt[i] = (char*)malloc(15 * sizeof(char));
    StrVec phaseName_alt(temp_MaxBeta, "");  //全部初始化为空
    //
    //double* beta_alt = NULL;
    //double** compositions_alt = NULL;
    //double** phi_alt = NULL;
    //double* zOutput_alt = NULL;
    //double* Compr_alt = NULL;
    //double* density_alt = NULL;
    int counterEquilibrium_final_alt = 0, iter_final_alt = 0, counter_final_alt = 0;

    // sizes
    //NumGases;       // In VB: UBound(zGlobal)
    DoubleVec Compr_alt;
    DoubleVec beta_alt;
    DoubleVec density_alt;
    DoubleMatrix phi_alt;
    DoubleMatrix compositions_alt;
    DoubleVec zOutput_alt;

    // 重新调整大小并初始化为0
    Compr_alt.resize(temp_MaxBeta, 0.0);
    beta_alt.resize(temp_MaxBeta, 0.0);
    density_alt.resize(temp_MaxBeta, 0.0);

    // 调整phi_alt大小：NumGases行，每行temp_MaxBeta列
    // 调整compositions_alt大小：NumGases行，每行temp_MaxBeta+1列
    phi_alt.resize(NumGases);
    compositions_alt.resize(NumGases);
    for (i = 0; i < NumGases; i++) {
        phi_alt[i].resize(temp_MaxBeta, 0.0);
        compositions_alt[i].resize(temp_MaxBeta + 1, 0.0);
    }

    zOutput_alt.resize(NumGases, 0.0);

    //Compr_alt = (double*)malloc(temp_MaxBeta * sizeof(double));
    //beta_alt = (double*)malloc(temp_MaxBeta * sizeof(double));
    //density_alt = (double*)malloc(temp_MaxBeta * sizeof(double));
    //phi_alt = (double**)malloc(NumGases * sizeof(double*));
    //compositions_alt = (double**)malloc(NumGases * sizeof(double*));
    //zOutput_alt = (double*)malloc(NumGases * sizeof(double));
    //for (i = 0; i < NumGases; i++)
    //{
    //    phi_alt[i] = (double*)malloc(temp_MaxBeta * sizeof(double));
    //    compositions_alt[i] = (double*)malloc((temp_MaxBeta + 1) * sizeof(double));
    //}
    ////初始化
    //memset(Compr_alt, 0, temp_MaxBeta * sizeof(double));
    //memset(beta_alt, 0, temp_MaxBeta * sizeof(double));
    //memset(density_alt, 0, temp_MaxBeta * sizeof(double));
    //memset(Compr_alt, 0, temp_MaxBeta * sizeof(double));
    //memset(beta_alt, 0, temp_MaxBeta * sizeof(double));
    //memset(Compr_alt, 0, temp_MaxBeta * sizeof(double));
    //for (i = 0; i < NumGases; i++)
    //{
    //    zOutput_alt[i] = 0;
    //    memset(phi_alt[i], 0, temp_MaxBeta * sizeof(double));
    //    memset(compositions_alt[i], 0, (temp_MaxBeta + 1) * sizeof(double));
    //    memset(phi_alt[i], 0, temp_MaxBeta * sizeof(double));
    //    memset(compositions_alt[i], 0, (temp_MaxBeta + 1) * sizeof(double));
    //}

    // conditions
    bool eqVapor = false, eqAqueous = false;
    if ((zGlobal[NumGases - 1] > 0.99) || ((zGlobal[NumGases - 1] > 0) && (NonZeroNumGases < 4) && (NumGases < 7))) {
        eqVapor = false;
        eqAqueous = true;
    }


    // First Flash
    Flash(&eqVapor, &eqAqueous, EOS, TK, PBar, NonZeroNumGases, zGlobal, gNeut, aH2O,
        TCr, PCr, Omega, MWgas, kPr, c0, c1,
        &Q, Numbeta, phaseName, beta, compositions, phi, Compr, density,
        counterEquilibrium_final, iter_final, counter_final, zOutput,
        NumGases, temp_MaxBeta);
    if (*Numbeta < 3) {
        eqVapor = true;
        eqAqueous = false;
        if ((zGlobal[NumGases - 1] > 0.99) || ((zGlobal[NumGases - 1] > 0) && (NonZeroNumGases < 4) && (NumGases < 7))) {
            eqVapor = false;
            eqAqueous = false;
        }

        Flash(&eqVapor, &eqAqueous, EOS, TK, PBar, NonZeroNumGases, zGlobal, gNeut, aH2O,
            TCr, PCr, Omega, MWgas, kPr, c0, c1,
            &Q_alt, &numBeta_alt, phaseName_alt, beta_alt, compositions_alt, phi_alt, Compr_alt, density_alt,
            &counterEquilibrium_final_alt, &iter_final_alt, &counter_final_alt, zOutput_alt,
            NumGases, temp_MaxBeta);

        if (Q > Q_alt) {
            *MaxBeta = 3;

            //重新分配空间
            //ReDim Compr(MaxBeta), beta(MaxBeta), density(MaxBeta), phi(NumGases, MaxBeta), compositions(NumGases, MaxBeta + 1), zOutput(NumGases)
            Compr.assign(*MaxBeta, 0);
            beta.assign(*MaxBeta, 0);
            density.assign(*MaxBeta, 0);
            zOutput.assign(*MaxBeta, 0);
            // 当需要重新分配时：
            phi.assign(NumGases, DoubleVec(*MaxBeta, 0.0));
            compositions.assign(NumGases, DoubleVec(*MaxBeta + 1, 0.0));
            *Numbeta = numBeta_alt;
            for (j = 0; j < *MaxBeta; j++) {
                //phaseName[j] = phaseName_alt[j];
                if (j != *MaxBeta - 1)
                {
                    phaseName[j] = phaseName_alt[j];
                    beta[j] = beta_alt[j];
                    Compr[j] = Compr_alt[j];
                    density[j] = density_alt[j];
                }
                else
                {
                    beta[j] = 0;
                    Compr[j] = 0;
                    density[j] = 0;
                }

                for (i = 0; i < NumGases; i++) {
                    if (j != *MaxBeta - 1)
                    {
                        compositions[i][j + 1] = compositions_alt[i][j + 1];
                        phi[i][j] = phi_alt[i][j];
                        zOutput[i] = zOutput_alt[i];
                    }
                    else
                    {
                        compositions[i][j + 1] = 0;
                        phi[i][j] = 0;
                        zOutput[i] = 0;
                    }
                }
            }
        }
    }
}


// MultiPhaseFlash。  变量mf_Vg不知道有什么用，暂时当作局部变量。
void MultiPhaseFlash(bool& mf_ParametersWereRead, DoubleVec& TCr, DoubleVec& PCr, DoubleVec& Omega,
    DoubleVec& MWgas, DoubleMatrix& kPr, DoubleVec& c0, DoubleVec& c1,
    double TK, double PBar, double total_moles, DoubleVec& z, DoubleVec& gNeut, double aH2O, DoubleVec& density, DoubleMatrix& compositions,
    DoubleMatrix& phi, DoubleVec& Compr, DoubleVec& beta, double* zOutput,
    DoubleVec& mass_phase, DoubleVec& MW_Phase, int* No_Phases)
{
    /*
    Dim max_NumGases As Long, NonZeroNumGases As Long, mf_NumGases As Long, MaxBeta As Long, i As Long, j As Long, k As Long, L As Long, m As Long
Dim mf_TK As Double, mf_PBar As Double, mf_Vg As Double, mf_mass_w_GO As Double
Dim zInput() As Double, zOut() As Double, mf_compositions() As Double, mf_gNeut() As Double, mf_beta() As Double, mf_Compr() As Double
Dim mf_density() As Double, mf_phi() As Double
Dim counterEquilibrium_final As Long, iter_final As Long, counter_final As Long
Dim EOS As String
Dim eqAqueous As Boolean
Dim eqVapor As Boolean
Dim x() As Double, mf_aH2O As Double, mf_nTCO2_GO As Double, mf_nTH2S_GO As Double
Dim iPure As Long
Dim zGlobal() As Double
Dim Numbeta As Long
Dim iFlash() As Long, pure As Long, Sz As Double
Dim logphipure() As Double, PsatPure As Double
Dim phase As String
Dim phaseName() As String
Dim logphipureL As Double, ComprL As Double, logphipureV As Double, ComprV As Double
    */
    char EOS[] = "PR"; //建立了状态方程为 Peng-Robinson: PR，或 Soave-RK: SRK
    int max_NumGases = 15; //这确定了组件的数量和可能的最大阶段数

    const double my_RBar = 83.14;//Const RBar As Double = 83.14    ' cm3 bar / (mol K)

    //变量定义
    int NonZeroNumGases, mf_NumGases, MaxBeta, i, j, k, L, m = 0;
    double mf_Vg, mf_mass_w_GO;
    int counterEquilibrium_final, counter_final;
    bool eqAqueous, eqVapor;
    double mf_nTCO2_GO, mf_nTH2S_GO;
    int iPure;
    int Numbeta;
    int pure;
    double Sz;

    DoubleVec logphipure(15, 0);
    double PsatPure;

    std::string phase = "vapor";
    StrVec phaseName;
    double logphipureL = 0, ComprL = 0, logphipureV = 0, ComprV = 0;

    //这表明，如果没有水，就不应该存在水相
    z[max_NumGases - 1] > 0 ? MaxBeta = 3 : MaxBeta = 2;
    phaseName = StrVec(MaxBeta, "");

    DoubleVec mf_gNeut(2, 0);
    DoubleVec zInput(max_NumGases, 0);
    //将值从外部变体类型变量输出传递到双精度类型变量
    double mf_TK = TK, mf_PBar = PBar, mf_aH2O = aH2O;
    for (i = 0; i < 2; i++)
        mf_gNeut[i] = gNeut[i];
    for (i = 0; i < max_NumGases; i++)
        zInput[i] = z[i];

    mf_NumGases = 6;
    NonZeroNumGases = 0;

    for (i = 0; i < max_NumGases - 2; i++)
    {
        if (zInput[i] > 0)
        {
            NonZeroNumGases++;
            //pure = i;   知道pure是用于判断而非数组索引，因此和vb保持一致，从1开始，需要加1
            pure = i + 1;
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
        pure = max_NumGases; //vb原语句，不需要偏移进行偏移操作
    }
    if (zInput[max_NumGases - 1] > 0 && NonZeroNumGases == 2 && (zInput[0] > 0 || zInput[1] > 0 || zInput[2] > 0 ||
        zInput[3] > 0 || zInput[max_NumGases - 2] > 0))
    {
        MaxBeta = 2;
    }
    if (NonZeroNumGases == 1)MaxBeta = 3;
    //赋值后，就可以开辟一些需要的数组了。

    DoubleVec zGlobal(mf_NumGases, 0);

    IntVec iFlash(max_NumGases, 0);

    DoubleVec zOut(mf_NumGases, 0);

    DoubleVec mf_beta(MaxBeta, 0);
    DoubleVec mf_Compr(MaxBeta, 0);
    DoubleVec mf_density(MaxBeta, 0);

    DoubleMatrix mf_phi(mf_NumGases, DoubleVec(MaxBeta, 0));
    DoubleMatrix mf_compositions(mf_NumGases, DoubleVec(MaxBeta + 1, 0));

    //ReDim mf_Compr(MaxBeta), mf_density(MaxBeta), mass_phase(MaxBeta)，注，参数列表参数mass_phase的ReDim行为
    mass_phase = DoubleVec(MaxBeta, 0);
    //double* x = (double*)malloc(mf_NumGases * sizeof(double));
    DoubleVec x(mf_NumGases, 0);

    //获取参数值、缩小组成以适应现有化合物并使其标准化的子程序
    mf_ParametersWereRead = false;
    //Call InitialPreparationSSP(mf_ParametersWereRead, EOS, zInput, iFlash, zGlobal, TCr, PCr, Omega, MWgas, kPr, c0, c1)
    InitialPreparationSSP(mf_ParametersWereRead, EOS, zInput, iFlash, zGlobal, TCr, PCr, Omega,
        MWgas, kPr, c0, c1, max_NumGases, mf_NumGases);

    counterEquilibrium_final = 0;
    int iter_final = 0;
    counter_final = 0;//初始化防止传入空值

    //计算多次闪蒸平衡的主要子程序
    if (NonZeroNumGases > 1)
    {
        Numbeta = 0;//原代码无这个，这里赋值防止传入未初始化的值。

        TrueFlash(EOS, mf_TK, mf_PBar, NonZeroNumGases, zGlobal, mf_gNeut, mf_aH2O, TCr, PCr, Omega, MWgas, kPr, c0, c1,
            &Numbeta, phaseName, mf_beta, mf_compositions, mf_phi, mf_Compr, mf_density,
            &counterEquilibrium_final, &iter_final, &counter_final, zOut, mf_NumGases, &MaxBeta);
    }
    else if (NonZeroNumGases == 1)//纯组分情况下的性质计算
    {
        Numbeta = 1;
        for (i = 0; i < mf_NumGases; i++)
            if (iFlash[i] == pure)iPure = i;
        if (pure < 3 || (pure == max_NumGases))
        {
            phase = "vapor";
            k = 1;
            if (TK < TCr[iPure])
            {
                if (pure == 1)PsatPure = 1;//CH4的位置
                if (pure == 2)PsatPure = PsatCO2(mf_TK);//CO2的位置
                if (pure == max_NumGases)PsatPure = PsatH2O(mf_TK);//水的位置，因为默认情况下水是最后一个

                if (PBar > PsatPure && PsatPure != 0)
                {
                    phase = "liquid";
                    if (pure = max_NumGases)k = 3;
                }
            }
        }
        //Call phi_calc(False, False, EOS, phase, mf_TK, mf_PBar, zGlobal, zGlobal, mf_gNeut, mf_aH2O, TCr, PCr, Omega, c0, c1, kPr, logphipure, mf_Compr(k))
        //mf_NumGases
        phi_calc(false, false, EOS, phase, mf_TK, mf_PBar, zGlobal, zGlobal, mf_gNeut, mf_aH2O, TCr, PCr,
            Omega, c0, c1, kPr, logphipure, &mf_Compr[k - 1], mf_NumGases);//NumGases
        if (phase == "vapor")
            Compr[k - 1] = ComprV; // vb源码是Compr(k) = ComprV，k并非从0开始，因此做偏移
        // 上面注意到，iPure是和循环相关的变量，从0开始的，因此不做偏移
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
                mf_Compr[k - 1] = ComprL;
                logphipure[iPure] = logphipureL;
            }
            else
            {
                k = 1;
                mf_Compr[k - 1] = ComprV;
                logphipure[iPure] = logphipureV;
            }
        }
        mf_phi[iPure][k - 1] = exp(logphipure[iPure]);//Exp
        mf_compositions[iPure][k] = 1.0;

        //相密度计算
        if (mf_Compr[k - 1] > 0)
        {
            mf_density[k - 1] = MWgas[iPure] * PBar / my_RBar / mf_Compr[k - 1] / TK;
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
            DoubleVec mf_beta_temp = mf_beta;
            DoubleVec mf_density_temp = mf_density;
            DoubleVec mf_Compr_temp = mf_Compr;

            //给空间重新定义
            mf_beta.assign(MaxBeta, 0);
            mf_density.assign(MaxBeta, 0);
            mf_Compr.assign(MaxBeta, 0);
            //将原来的数据拿回来
            for (i = 0; i < temMaxBeta; i++)
            {
                mf_beta[i] = mf_beta_temp[i];
                mf_density[i] = mf_density_temp[i];
                mf_Compr[i] = mf_Compr_temp[i];
            }
            if (zInput[max_NumGases - 1] > 0)
            {
                mf_beta[2] = mf_beta_temp[1];
                mf_density[2] = mf_density_temp[1];
                mf_Compr[2] = mf_Compr_temp[1];
            }

            //二维数组同理，先保存
            DoubleMatrix mf_compositions_temp = mf_compositions;
            DoubleMatrix mf_phi_temp = mf_phi;

            //然后让它们的长度进行变化
            for (auto& row : mf_compositions)
                row.resize(MaxBeta);
            for (auto& row : mf_phi)
                row.resize(MaxBeta);
            // 从临时数组复制回主数组
            for (i = 0; i < mf_NumGases; i++)
                for (j = 0; j < temMaxBeta; j++)
                {
                    mf_compositions[i][j] = mf_compositions_temp[i][j];
                    mf_phi[i][j] = mf_phi_temp[i][j];
                }
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
        }
    }
    //ReDim MW_Phase(MaxBeta), mass_phase(MaxBeta)
    MW_Phase.resize(MaxBeta, 0);
    mass_phase.resize(MaxBeta, 0);

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
        MW_Phase[j] = weighted_mean(x, MWgas, mf_NumGases);
        mass_phase[j] = total_moles * mf_beta[j] * MW_Phase[j] / 1000;
        density[j] = mf_density[j];
        Compr[j] = mf_Compr[j];
    }
    // 重新排序输出（CALEB）
    OrderPhases(compositions, phi, beta, MW_Phase, mass_phase, density, Compr, max_NumGases, MaxBeta);
    *No_Phases = Numbeta;
    mf_Vg = Compr[0] * total_moles * beta[0] * 0.00008314461 * mf_TK / PBar; // 气体体积 m^3

}


void Get_EOS_Parameters(int NumGases, std::string& EOS, DoubleVec& MWgas, DoubleVec& TCr,
    DoubleVec& PCr, DoubleVec& Omega, DoubleVec& mf_c0, DoubleVec& mf_c1,
    DoubleMatrix& kPr)
{
    int i, j;

    //for (i = 0; i < NumGases; i++) {
    //    //MWgas(i) = Worksheets("Input").Cells(4 + i, 24)
    //    //    TCr(i) = Worksheets("Input").Cells(4 + i, 26)
    //    //    PCr(i) = Worksheets("Input").Cells(4 + i, 27)
    //    //    Omega(i) = Worksheets("Input").Cells(4 + i, 28)
    //    if (EOS == "PR") {
    //        //mf_c0(i) = Worksheets("Input").Cells(4 + i, 29)
    //        //    mf_c1(i) = Worksheets("Input").Cells(4 + i, 30)
    //    }
    //    if (EOS == "SRK") {
    //        //mf_c0(i) = Worksheets("Input").Cells(4 + i, 31)
    //        //    mf_c1(i) = Worksheets("Input").Cells(4 + i, 32)
    //    }

    //    for (j = 0; j < i; j++) {
    //        //(*kPr)[i][j] = Worksheets("Input").Cells(4 + i, 34 + j)
    //    }
    //}
    MWgas[0] = 16.043;
    TCr[0] = 190.6;
    PCr[0] = 46;
    Omega[0] = 0.008;

    MWgas[1] = 44.01;
    TCr[1] = 304.2;
    PCr[1] = 73.76;
    Omega[1] = 0.225;

    MWgas[2] = 34.08;
    TCr[2] = 373.5;
    PCr[2] = 89.63;
    Omega[2] = 0.094;

    MWgas[3] = 30.07;
    TCr[3] = 305.5;
    PCr[3] = 48.84;
    Omega[3] = 0.098;

    MWgas[4] = 44.094;
    TCr[4] = 369.9;
    PCr[4] = 42.46;
    Omega[4] = 0.152;

    MWgas[5] = 58.124;
    TCr[5] = 408.2;
    PCr[5] = 36.48;
    Omega[5] = 0.176;

    MWgas[6] = 58.124;
    TCr[6] = 425.3;
    PCr[6] = 38;
    Omega[6] = 0.193;

    MWgas[7] = 72.151;
    TCr[7] = 460.45;
    PCr[7] = 33.84;
    Omega[7] = 0.227;

    MWgas[8] = 72.151;
    TCr[8] = 469.6;
    PCr[8] = 33.74;
    Omega[8] = 0.251;

    MWgas[9] = 86.178;
    TCr[9] = 507.5;
    PCr[9] = 29.69;
    Omega[9] = 0.296;

    MWgas[10] = 108.46355910141;
    TCr[10] = 575.5;
    PCr[10] = 26.26;
    Omega[10] = 0.545;

    MWgas[11] = 203.768522560204;
    TCr[11] = 708.6;
    PCr[11] = 16.76;
    Omega[11] = 0.857;

    MWgas[12] = 321.359072235909;
    TCr[12] = 945.6;
    PCr[12] = 13.36;
    Omega[12] = 1.248;

    MWgas[13] = 28.013;
    TCr[13] = 126.2;
    PCr[13] = 33.94;
    Omega[13] = 0.040;

    MWgas[14] = 18;
    TCr[14] = 647.1;
    PCr[14] = 220.55;
    Omega[14] = 0.345;

    if (EOS == "PR") {
        mf_c0[0] = -5.19998589998324;
        mf_c1[0] = 0.0;

        mf_c0[1] = -1.90887591231191;
        mf_c1[1] = 0.0;

        mf_c0[2] = -3.92147265302369;
        mf_c0[3] = -5.7950356085208;
        mf_c0[4] = -6.35369532428396;
        mf_c0[5] = -7.18062223841032;
        mf_c0[6] = -6.48762839907065;
        mf_c0[7] = -6.19849948277172;
        mf_c0[8] = -5.12105162468565;
        mf_c0[9] = -3.4814313848538;
        mf_c0[10] = -14.8622158330974;
        mf_c0[11] = -29.8260575190994;
        mf_c0[12] = 5.77517376259963;
        mf_c0[13] = -4.23208315053643;
        mf_c0[14] = 3.00773525559068;

        mf_c1[2] = 0.0;
        mf_c1[3] = 0.0;
        mf_c1[4] = 0.0;
        mf_c1[5] = 0.0;
        mf_c1[6] = 0.0;
        mf_c1[7] = 0.0;
        mf_c1[8] = 0.0;
        mf_c1[9] = 0.0;
        mf_c1[10] = -0.00560031997876096;
        mf_c1[11] = 0.101896336756916;
        mf_c1[12] = 0.279177061769215;
        mf_c1[13] = 0.0;
        mf_c1[14] = 0.00873186552687633;
    }
    if (EOS == "SRK") {
        //??
    }

    kPr[0][0] = 0;

    kPr[1][0] = 0.12;

    kPr[2][0] = 0.08;
    kPr[2][1] = 0.12;

    kPr[3][1] = 0.15;
    kPr[3][2] = 0.07;

    kPr[4][1] = 0.15;
    kPr[4][2] = 0.07;

    kPr[5][1] = 0.15;
    kPr[5][2] = 0.06;

    kPr[6][1] = 0.15;
    kPr[6][2] = 0.06;

    kPr[7][1] = 0.15;
    kPr[7][2] = 0.06;

    kPr[8][1] = 0.15;
    kPr[8][2] = 0.06;

    kPr[9][1] = 0.15;
    kPr[9][2] = 0.05;

    kPr[10][1] = 0.15;
    kPr[10][2] = 0.05;

    kPr[11][1] = 0.15;
    kPr[11][2] = 0.05;

    kPr[12][1] = 0.15;
    kPr[12][2] = 0.05;

    kPr[13][0] = 0.02;
    kPr[13][3] = 0.06;
    kPr[13][4] = 0.08;
    kPr[13][5] = 0.08;
    kPr[13][6] = 0.08;
    kPr[13][7] = 0.08;
    kPr[13][8] = 0.08;
    kPr[13][9] = 0.08;
    kPr[13][10] = 0.08;
    kPr[13][11] = 0.08;
    kPr[13][12] = 0.08;

    kPr[14][0] = 0.45;
    kPr[14][3] = 0.45;
    kPr[14][4] = 0.53;
    kPr[14][5] = 0.52;
    kPr[14][6] = 0.52;
    kPr[14][7] = 0.5;
    kPr[14][8] = 0.5;
    kPr[14][9] = 0.5;
    kPr[14][10] = 0.5;
    kPr[14][11] = 0.5;
    kPr[14][12] = 0.5;
}


void true_composition(double TK, double PBar, double mol_HC, double mol_W, double temp_aH2O, DoubleVec& gNeut, double nTCO2, double nTH2S,
    int useEOS, DoubleVec& reservoir_Composition, DoubleVec& feed_Composition, double* total_moles)
{
    int i, j, NumGases = 15;
    /*
Dim composition_G() As Double
Dim MWgas() As Double, TCr() As Double, PCr() As Double, Omega() As Double, mf_c0() As Double, mf_c1() As Double, kPr() As Double 'lnphi_Gas() As Double,
'Dim xH2O_aq As Double
Dim lnphiH2O_Aq As Double, composition_Aq() As Double, lnphi_Aq() As Double, Compr_composition_Aq As Double
Dim EOS As String, mf_gNeut(2) As Double
Const RBar As Double = 0.000083144621 ' m3*Bar/mol*K
    */
    double lnphiH2O_Aq = 0;
    double Compr_composition_Aq = 0;

    const double RBar = 0.000083144621; // 注意：这里源码的RBar和全局的那个重名了，但数值不同
    std::string EOS = "PR";
    DoubleVec mf_gNeut(2, 0);

    mf_gNeut[0] = gNeut[0];
    mf_gNeut[1] = gNeut[1];

    /*
        , mf_feed_Composition(NumGases)  这两个数组没有被任何函数调用过，在这里凭空产生ReDim，不明确什么意思，此处已经省略

ReDim TCr(NumGases), PCr(NumGases), MWgas(NumGases), Omega(NumGases), mf_c0(NumGases), mf_c1(NumGases), kPr(NumGases, NumGases), lnphi_Gas(5)
ReDim composition_Aq(NumGases), lnphi_Aq(NumGases)
*/

// 分配内存

    DoubleVec MWgas(NumGases, 0);
    DoubleVec TCr(NumGases, 0);
    DoubleVec PCr(NumGases, 0);
    DoubleVec Omega(NumGases, 0);
    DoubleVec mf_c0(NumGases, 0);
    DoubleVec mf_c1(NumGases, 0);
    DoubleMatrix kPr(NumGases, DoubleVec(NumGases, 0));


    DoubleVec lnphi_Gas(5, 0);
    DoubleVec composition_Aq(NumGases, 0);
    DoubleVec lnphi_Aq(NumGases, 0);

    // 获取EOS参数
    Get_EOS_Parameters(NumGases, EOS, MWgas, TCr, PCr, Omega, mf_c0, mf_c1, kPr);

    // 初始化feed_Composition数组
    for (i = 0; i < NumGases; i++) {
        feed_Composition[i] = 0.0;
    }

    // 计算各组分的摩尔数（不包括水）
    // 注意：VB中reservoir_Composition索引从1到NumGases-1，对应C中的0到NumGases-2
    for (i = 0; i < (NumGases - 1); i++) {
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
    for (i = 0; i < NumGases; i++) {
        // if UseEOS=2 then total_moles are different from that calc in partD
        *total_moles += feed_Composition[i];
    }

    // 计算摩尔分数
    for (i = 0; i < NumGases; i++) {
        feed_Composition[i] = feed_Composition[i] / (*total_moles);
    }
}


void pseudo_composition(double API, double SGG, double VgTP, double mol_opd, double mol_wpd, double TK,
    double PBar, double aH2O, DoubleVec& gNeut, double nTCO2, double nTH2S, double temp_yCO2, double temp_yH2S, double YH2O,
    double* total_moles, DoubleVec& temp_feed_Composition, double* mol_HC) {
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
    int i, j;
    int NumGases = 15;

    //-----参考true_composition

    DoubleVec MWgas(NumGases, 0);
    DoubleVec TCr(NumGases, 0);
    DoubleVec PCr(NumGases, 0);
    DoubleVec Omega(NumGases, 0);
    DoubleVec mf_c0(NumGases, 0);
    DoubleVec mf_c1(NumGases, 0);


    // 分配二维数组kPr
    DoubleMatrix kPr(NumGases, DoubleVec(NumGases, 0));

    DoubleVec lnphi_Aq(NumGases, 0);
    DoubleVec lnphi_Gas(NumGases, 0);
    DoubleVec composition_G(NumGases, 0);
    DoubleVec composition_Aq(NumGases, 0);
    DoubleVec lnphi_Water(NumGases, 0);



    double SGO = 0, sgLightO = 0, sgHeavyO = 0, mwLightO = 0, mwHeavyO = 0, mwLightG = 0, mwHeavyG = 0, mwAir = 0;
    double xLightO = 0, mwO = 0, mwG = 0, yLightG = 0, yHeavyG = 0;
    double mol_gpd = 0, mol_LightOpd = 0, mol_HeavyOpd = 0, mol_LightGpd = 0;
    double mol_HeavyGpd = 0;
    double lnphiH2O_Aq = 0, Compr_composition_Aq = 0;
    double Compr_G = 0;

    std::string EOS = "PR";
    DoubleVec mf_gNeut(2, 0);

    const double RBar = 0.000083144621;  // m3*Bar/mol*K



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


    // calculation of mole fractions of hexane and toluene equivalent in oil phase (from API gravity)

    xLightO = mwHeavyO * (sgHeavyO - SGO) / (mwLightO * (sgHeavyO / sgLightO) * (SGO - sgLightO) + mwHeavyO * (sgHeavyO - SGO));  // mole fraction of hexane in the oil

    mwO = xLightO * mwLightO + (1 - xLightO) * mwHeavyO;  // calculated molecular weight of oil

    // calculation of mole fractions of methane and butane equivalent in the gas phase (from SG of gas)

    mwG = SGG * mwAir;  // calculated molar mass of gas

    // Subroutine that gets critical properties and parameters for the components
    Get_EOS_Parameters(NumGases, EOS, MWgas, TCr, PCr, Omega, mf_c0, mf_c1, kPr);

    // gNeut[0] 有问题 vb 1.3378450925412 CPP 1.337801774221661
    mf_gNeut[0] = gNeut[0];
    mf_gNeut[1] = gNeut[1];
    /*     VB 中注释了           */
    //  Sub phi_calc(eqVapor,eqAqueous, EOS As String, phase As String, 
    // TK As Double, PBar As Double, x() As Double, xGlobal() As Double, gNeut() As Double,
    //  aH2O As Double, TCr() As Double, PCr() As Double, Omega() As Double, c0() As Double, c1() As Double, kPr() As Double, lnphi() As Double, z As Double)
    // Call phi_calc(False, False, EOS, "liquid", CDbl(TK), CDbl(PBar),
    //  composition_Aq, composition_Aq, mf_gNeut, CDbl(aH2O), TCr, PCr, Omega, mf_c0, mf_c1, kPr, lnphi_Water, Compr_composition_Aq)  'commented out by Amy ???
    // lnphiH2O_Aq = lnphi_Water(NumGases)
    /*     VB 中注释了           */

    yLightG = (mwHeavyG - mwG - mwHeavyG * (simContext.yCO2 + simContext.yH2S + YH2O) + simContext.yCO2 * 44.01 + simContext.yH2S * 34.08 + YH2O * 18.01528) / (mwHeavyG - mwLightG);  // mole fraction of methane in the gas phase
    if (yLightG < 0) yLightG = 0;

    yHeavyG = 1 - (yLightG + simContext.yCO2 + simContext.yH2S + YH2O);
    if (yHeavyG < 0) yHeavyG = 0;
    // This calculates the compressibility factor of the vapor phase
    for (i = 3; i < NumGases - 1; i++) {  // 0-based: i=3 to 12 for components 4 to 13
        composition_G[i] = 0.0;
    }
    double denom = yLightG + simContext.yCO2 + simContext.yH2S + yHeavyG + YH2O;
    composition_G[0] = yLightG / denom;  // methane
    composition_G[1] = simContext.yCO2 / denom;  // CO2
    composition_G[2] = simContext.yH2S / denom;  // H2S
    composition_G[9] = yHeavyG / denom;  // n-butane (index 10 in 1-based)
    composition_G[NumGases - 1] = YH2O;  // water (index 15 in 1-based)

    phi_calc(0, 0, EOS, "vapor", TK, PBar, composition_G, composition_G, mf_gNeut, aH2O, TCr, PCr, Omega, mf_c0, mf_c1, kPr, lnphi_Gas, &Compr_G, NumGases);
    mol_gpd = VgTP * PBar / (Compr_G * RBar * TK);  // moles of gas per day -- R is in m3/bar/K/g-mol

    mol_LightOpd = xLightO * mol_opd;  // moles of hexane equivalent per day

    mol_HeavyOpd = (1 - xLightO) * mol_opd;  // moles of toluene equivalent per day

    mol_LightGpd = yLightG * mol_gpd;  // moles of methane equivalent per day

    mol_HeavyGpd = (1 - yLightG) * mol_gpd;  // moles of butane equivalent per day
    // total_moles = mol_gpd + mol_wpd + mol_opd  ' TOTAL MOLES PER DAY (GAS + OIL + WATER)

    for (i = 0; i < NumGases; i++) {
        temp_feed_Composition[i] = 0.0;
    }

    temp_feed_Composition[0] = mol_LightGpd;  // mole fraction of methane in the feed
    temp_feed_Composition[1] = nTCO2;  // mole fraction of CO2 in the feed
    temp_feed_Composition[2] = nTH2S;  // mole fraction of H2S in the feed
    temp_feed_Composition[6] = mol_HeavyGpd;  // mole fraction of n-butane in the feed (index 7 in 1-based)
    temp_feed_Composition[9] = mol_LightOpd;  // mole fraction of n-hexane in the feed (index 10)
    temp_feed_Composition[12] = mol_HeavyOpd;  // mole fraction of toluene in the feed (index 13)
    temp_feed_Composition[14] = mol_wpd;  // mole fraction of water in the feed (index 15)

    *total_moles = 0;
    for (i = 0; i < NumGases; i++) {
        *total_moles += temp_feed_Composition[i];
    }
    *mol_HC = *total_moles - mol_wpd;
    for (i = 0; i < NumGases; i++) {
        temp_feed_Composition[i] /= *total_moles;
    }
}


void Flash_Input_Processing(char* ioSheet, char* ioCol, double eosProps[][6], double kij[][15], DoubleVec& zInput, int nComp,
    DoubleVec& zi,        // 输出数组地址（函数内malloc）
    IntVec& idx_CompA,    // 输出索引数组地址
    int& nCompA, IntVec& idx_Henry,
    int& idx_Water,
    IntVec& phaseID,      // 输出相标识数组地址
    int& nPhase
)
{
    int i, k;

    // 初始化
    for (i = 0; i < 5; ++i) idx_Henry[i] = 0;
    idx_Water = 0;
    nCompA = 0;

    // -------------------------
    // 1. 计算活性（非零）组分数量
    // -------------------------
    for (i = 0; i < nComp; ++i) {
        if (zInput[i] > 0.0)
            nCompA++;

        // 识别轻气体 (H2, CO2, H2S, CH4)
        if (i < 4 && zInput[i] > 0.0)
            idx_Henry[i] = i;    // 原语句是idx_Henry(i) = i，，后面用到idx_Henry(i)的时候，请务必确认此处的idx_Henry(i)和vb的值不一样了。
    }

    // 判断 N2 是否存在（nComp-2位置）
    if (zInput[nComp - 2] > 0.0)
        idx_Henry[4] = nCompA - 2;//原语句：idx_Henry(5) = nCompA - 1，但上方可知，我们的idx_Henry是有偏移的，所以改成了-2

    // 判断 H2O 是否存在（最后一个）
    if (zInput[nComp - 1] > 0.0)
        idx_Water = nCompA;

    // -------------------------
    // 2. 填充活性组分向量
    // -------------------------
    zi = DoubleVec(nCompA);
    idx_CompA = IntVec(nCompA);

    k = 0;
    for (i = 0; i < nComp; ++i) {
        if (zInput[i] > 0.0) {
            zi[k] = zInput[i];
            idx_CompA[k] = i; // C 用0-based索引
            k++;
        }
    }

    // -------------------------
    // 3. 归一化组分
    // -------------------------
    //NormalizeVector(*zi, *nCompA);

    // -------------------------
    // 4. 若 EOS 参数尚未加载
    // -------------------------
    if (eosProps[0][0] == 0.0) {
        // 从 sheet 读取参数
        //sheet2vba(ioSheet, ioCol, eosProps, kij, nComp);
        // 重新排序、修剪非零组分
        //compReorder_c2a(eosProps, kij, nComp, *idx_CompA, *nCompA);
    }

    // -------------------------
    // 5. 判断存在几相
    // -------------------------
    if (idx_Water != 0) {
        nPhase = 3;
        phaseID = IntVec(nPhase, 0);
        phaseID[0] = 1;
        phaseID[1] = 2;
        phaseID[2] = 3;
    }
    else {
        nPhase = 2;
        phaseID = IntVec(nPhase, 0);
        phaseID[0] = 1;
        phaseID[1] = 2;
    }
}


void MultiPhaseFlash_CS(char* iosheet, char* ioCol, double eosProps[][6], double kij[][15], double TK, double PBar,
    double total_moles, DoubleVec& zInput, DoubleVec& gNeut, double aH2O, DoubleVec& density, DoubleMatrix& compositions,
    DoubleMatrix& phi, DoubleVec& Compr, DoubleVec& beta, double* zOutput,
    DoubleVec& mass_phase, DoubleVec& MW_Phase, int* No_Phases)
{
    int i, j, k, L;
    int FlashType = 1;  /* modified Wilson PT-flash */
    int EOS = 2;        /* Peng-Robinson EOS */
    int SSPactive = 1;

    int nCompMax = 15;  /* 最大组分数 */

    /* ========= 输入预处理：调用外部函数 ========== */
    DoubleVec zi;
    IntVec idx_CompA; int nCompA = 0;
    IntVec idx_Henry(5, 0); int idx_Water = 0;
    IntVec phaseID; int nPhaseMax = 0;
    Flash_Input_Processing(iosheet, ioCol, eosProps, kij, zInput, nCompMax,
        zi, idx_CompA, nCompA, idx_Henry, idx_Water,
        phaseID, nPhaseMax);

    /* ========= 将输入转为C数组 ========= */
    double T = TK;
    double P = PBar;
    DoubleVec MW(nCompA);
    double gamma[2];
    gamma[0] = gNeut[0];
    gamma[1] = gNeut[1];
    double aW = aH2O;

    for (i = 0; i < nCompA; i++) {
        MW[i] = eosProps[i][0]; /* VB的 eosProps(i,1) → C的 eosProps[i][0] */
    }

    /* ========= 调用 ShellFlash ========= */
    DoubleVec beta_mol(nPhaseMax);
    DoubleVec Zk(nPhaseMax);

    DoubleMatrix yik(nCompA, DoubleVec(nPhaseMax, 0));
    DoubleMatrix logPHI(nCompA, DoubleVec(nPhaseMax, 0));

    //ShellFlash(FlashType, EOS, T, P, zi, eosProps, kij, nCompA, nPhaseMax,
    //    phaseID, beta_mol, Zk, yik, logPHI,
    //    SSPactive, idx_Water, idx_Henry, aW, gamma);

    /* ========= 排序相并计算性质 ========= */
    int nPhase = 0;
    double MWk[3], vol_cm3mol[3], dens_gcm3[3];

    //为什么这里的参数和OrderPhases定义的参数不同？
    //C语言：OrderPhases(double** compositions, double** phi, double* beta, double* MW_Phase,
    //    double* mass_phase, double* density, double* Compr, int* nComp, int* nPhase)

    //vb原定义：Sub OrderPhases(compositions, phi, beta, MW_Phase, mass_phase, density, Compr, nComp, nPhase)
    // 参数列表里是数组类型的，但vb源码传的是标量
    //vb调用：Call OrderPhases(T, P, MW, nCompA, nPhaseMax, idx_Water, beta_mol, Zk, MWk, vol_cm3mol, dens_gcm3, yik, logPHI)
    //OrderPhases(T, P, MW, nCompA, nPhaseMax, idx_Water, beta_mol, Zk,
        //MWk, vol_cm3mol, dens_gcm3, yik, logPHI);


    //下面不知道是在做什么
    /*
         ' output properties to variant type variables
     ReDim density(nPhaseMax), Compr(nPhaseMax), beta(nPhaseMax), mass_phase(nPhaseMax), MW_Phase(nPhaseMax)
     ReDim zOutput(nCompMax), compositions(nCompMax, nPhaseMax + 1), phi(nCompMax, nPhaseMax)
    */

    /* ========= 输出到调用方变量 ========= */
    for (i = 0; i < nCompMax; i++) {
        compositions[i][0] = zInput[i];
        zOutput[i] = zInput[i];
    }
    nPhase = 0;
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
}


/**
 * @brief fMeSSpeciation 的迭代子过程
 * @param im   金属索引
 * @param igas 气体标志
 */
void MeSWhileloop(int im, int igas, double* HSOld, double ZP1, double ZP2, double ZP3, double ZP4, double ZP5, double ZP6, double ZP7, double ZP9)
{
    *HSOld = simContext.HS;
    double coef1 = 1.0;
    double coef2 = (mc[iZn] * ZP3 + mc[iPb] * ZP6) / (mc[iZn] * ZP4 + mc[iPb] * ZP7);
    double coef3 = (hydHS + mc[iFe] * ZP1) / (mc[iZn] * ZP4 + mc[iPb] * ZP7);
    double coef4 = -TH2Saq / (mc[iZn] * ZP4 + mc[iPb] * ZP7);

    double a = coef1 * pow(simContext.HS, 3) / coef4;
    double b = coef2 * pow(simContext.HS, 2) / coef4;

    if (fabs(a) < 1e-4) {
        if (fabs(b) < 1e-4) {
            simContext.HS = -coef4 / coef3;
        }
        else {
            QuadraticRoots(coef2, coef3, coef4, root1, root2);
            simContext.HS = (root2 > root1) ? root2 : root1;
        }
        if (simContext.HS <= 0) return;
    }
    else {
        CubicRoots(coef1, coef2, coef3, coef4, root1, root2, root3);
        simContext.HS = root1;
        if (root2 > simContext.HS) simContext.HS = root2;
        if (root3 > simContext.HS) simContext.HS = root3;
        if (simContext.HS <= 0) {
            simContext.HS = 0;
            return;
        }
    }

    // 更新 mc
    if (im == iFe && igas == 3)
        mc[iFe] = (TFe - ppt) / (1.0 + ZP1 * simContext.HS);
    else
        mc[iFe] = TFe / (1.0 + ZP1 * simContext.HS);

    if (im == iZn && igas == 3)
        mc[iZn] = (TZn - ppt) / (1.0 + ZP2 + ZP3 * pow(simContext.HS, 2) + ZP4 * pow(simContext.HS, 3) + ZP9 * pow(OH, 2));
    else
        mc[iZn] = TZn / (1.0 + ZP2 + ZP3 * pow(simContext.HS, 2) + ZP4 * pow(simContext.HS, 3) + ZP9 * pow(OH, 2));

    if (im == iPb && igas == 3)
        mc[iPb] = (TPb - ppt) / (1.0 + ZP5 + ZP6 * pow(simContext.HS, 2) + ZP7 * pow(simContext.HS, 3));
    else
        mc[iPb] = TPb / (1.0 + ZP5 + ZP6 * pow(simContext.HS, 2) + ZP7 * pow(simContext.HS, 3));
}

/**
 * @brief 金属硫化物水溶液物种化计算
 *        （只考虑水相，不包括与气相、油相的平衡；Cl- 物种化未考虑）
 *
 * @param im   金属索引 (如 iFe, iZn, iPb)
 * @param igas 气体标志 (3 表示硫化物沉淀)
 * @return double (未定义时返回 NAN)
 */
double fMeSSpeciation(int im, int igas) {
    double HSOld = simContext.HS;
    double FeOld = TFe;
    double ZnOld = TZn;
    double PbOld = TPb;

    // 预计算常用项
    double cl_term = gDot[iClDot] * ma[iCl];
    double hs_term_sq = pow(gDot[iHSDot], 2); // HS的平方项
    double hs_term_cb = pow(gDot[iHSDot], 3); // HS的立方项

    // 计算ZP1 - FeSaq配合物参数
    double ZP1 = (KstFeSaq * gAn[iHS] * gNAn[iHS] * gCat[iFe] * gNCat[iFe]) /
        (gNeut[iFeSaq] * gNNeut[iFeSaq] * aH);

    // 计算ZP2 - Zn-Cl配合物参数
    double ZP2 = BetaDot[iZnCl] * cl_term * gDot[iZnDot] / gDot[iZnCl] +
        BetaDot[iZnCl2] * pow(cl_term, 2) * gDot[iZnDot] / gDot[iZnCl2] +
        BetaDot[iZnCl3] * pow(cl_term, 3) * gDot[iZnDot] / gDot[iZnCl3] +
        BetaDot[iZnCl4] * pow(cl_term, 4) * gDot[iZnDot] / gDot[iZnCl4];

    // 计算ZP3 - Zn-HS配合物参数
    double ZP3 = BetaDot[iZnHS2] * gDot[iZnDot] * hs_term_sq / gDot[iZnHS2];

    // 计算ZP4 - Zn-HS配合物参数
    double ZP4 = BetaDot[iZnHS3] * gDot[iZnDot] * hs_term_cb / gDot[iZnHS3];

    // 计算ZP5 - Pb-Cl配合物参数
    double ZP5 = BetaDot[iPbCl] * cl_term * gDot[iPbDot] / gDot[iPbCl] +
        BetaDot[iPbCl2] * pow(cl_term, 2) * gDot[iPbDot] / gDot[iPbCl2] +
        BetaDot[iPbCl3] * pow(cl_term, 3) * gDot[iPbDot] / gDot[iPbCl3] +
        BetaDot[iPbCl4] * pow(cl_term, 4) * gDot[iPbDot] / gDot[iPbCl4];

    // 计算ZP6 - Pb-HS配合物参数
    double ZP6 = BetaDot[iPbHS2] * gDot[iPbDot] * hs_term_sq / gDot[iPbHS2];

    // 计算ZP7 - Pb-HS配合物参数
    double ZP7 = BetaDot[iPbHS3] * gDot[iPbDot] * hs_term_cb / gDot[iPbHS3];

    double ZP9 = BetaDot[iZnOH2] * gDot[iZnDot] * gDot[iOHDot] * gDot[iOHDot] / gDot[iZnOH2];

    /* ------- 根据是否有沉淀 (ppt) 来计算初始溶解相浓度 mc[] ------- */
       /* Fe */
    if (im == iFe && igas == 3) { /* FeS 正在沉淀 */
        FeOld = TFe - ppt;
        mc[iFe] = (TFe - ppt) / (1.0 + ZP1 * simContext.HS);
    }
    else {
        mc[iFe] = TFe / (1.0 + ZP1 * simContext.HS);
    }

    /* Zn  -  修改： 与excel保持一致 */
    if (im == iZn && igas == 3) { /* ZnS 正在沉淀 */
        ZnOld = TZn - ppt;
        mc[iZn] = (TZn - ppt) / (1.0 + ZP2 + ZP3 * pow(simContext.HS, 2.0) + ZP4 * pow(simContext.HS, 3.0) + ZP9 * pow(OH, 2.0));
    }
    else {
        mc[iZn] = TZn / (1.0 + ZP2 + ZP3 * pow(simContext.HS, 2.0) + ZP4 * pow(simContext.HS, 3.0) + ZP9 * pow(OH, 2.0));
    }

    /* Pb */
    if (im == iPb && igas == 3) { /* PbS 正在沉淀 */
        PbOld = TPb - ppt;
        mc[iPb] = (TPb - ppt) / (1.0 + ZP5 + ZP6 * pow(simContext.HS, 2.0) + ZP7 * pow(simContext.HS, 3.0));
    }
    else {
        mc[iPb] = TPb / (1.0 + ZP5 + ZP6 * pow(simContext.HS, 2.0) + ZP7 * pow(simContext.HS, 3.0));
    }

    /* ------- 有 Zn 或 Pb 时，求解多项式（可能为三次或二次）以求 simContext.HS ------- */
    if (mc[iZn] > 0.0 || mc[iPb] > 0.0) {
        double coef1 = 1.0;
        double denom = (mc[iZn] * ZP4 + mc[iPb] * ZP7);
        /* 若 denom = 0 可能造成除零，按 VB 的原逻辑直接除；如需防护可添加 epsilon 检查 */
        double coef2 = (mc[iZn] * ZP3 + mc[iPb] * ZP6) / denom;
        double coef3 = (hydHS + mc[iFe] * ZP1) / denom;
        double coef4 = -TH2Saq / denom;
        //if (igas == 3) {
        //    coef4 = -(TH2Saq - ppt) / denom;
        //}
        //else {
        //    coef4 = -TH2Saq / denom;
        //}

        double a = coef1 * pow(simContext.HS, 3.0) / coef4;
        double b = coef2 * pow(simContext.HS, 2.0) / coef4;

        if (fabs(a) < 1e-4) {
            if (fabs(b) < 1e-4) {
                simContext.HS = -coef4 / coef3;
            }
            else {
                QuadraticRoots(coef2, coef3, coef4, root1, root2);
                simContext.HS = (root2 > root1) ? root2 : root1;
            }
            if (simContext.HS <= 0.0) return NAN; /* VB: GoTo 1000 */
        }
        else {
            CubicRoots(coef1, coef2, coef3, coef4, root1, root2, root3);
            simContext.HS = root1;
            if (root2 > simContext.HS) simContext.HS = root2;
            if (root3 > simContext.HS) simContext.HS = root3;
            if (simContext.HS <= 0.0) {
                simContext.HS = 0.0;
                return NAN; /* VB: GoTo 1000 */
            }
        }
    }
    /* ------- 若只有 Fe 存在，则直接解析 simContext.HS ------- */
    else if (mc[iFe] > 0.0) {
        //if (igas == 3) {
        //    simContext.HS = (TH2Saq - ppt) / (hydHS + ZP1 * mc[iFe]);
        //}
        //else {
        //    simContext.HS = TH2Saq / (hydHS + ZP1 * mc[iFe]);
        //}
        simContext.HS = TH2Saq / (hydHS + ZP1 * mc[iFe]);
    }
    /* 若都不存在可解的金属，则返回（VB: GoTo 1000） */
    else {
        return NAN;
    }

    /* ========== 下面展开所有 VB 中的不同组合情形的迭代收敛循环 ========== */
    /* 1) If TFe > 0 And simContext.HS > 0 And TZn > 0 And TPb > 0 Then 'Fe,Zn,Pb sulfide */
    if (TFe > 0 && simContext.HS > 0 && TZn > 0 && TPb > 0) {
        while (fabs((FeOld / (1.0 + ZP1 * simContext.HS) - mc[iFe]) / mc[iFe]) > 0.001
            || fabs((PbOld / (1.0 + ZP5 + ZP6 * pow(simContext.HS, 2) + ZP7 * pow(simContext.HS, 3)) - mc[iPb]) / mc[iPb]) > 0.001
            || fabs((ZnOld / (1.0 + ZP2 + ZP3 * pow(simContext.HS, 2) + ZP4 * pow(simContext.HS, 3) + ZP9 * pow(OH, 2)) - mc[iZn]) / mc[iZn]) > 0.001
            || fabs((HSOld - simContext.HS) / simContext.HS) > 0.001)
        {
            MeSWhileloop(im, igas, &HSOld, ZP1, ZP2, ZP3, ZP4, ZP5, ZP6, ZP7, ZP9);
        }
    }

    /* 2) If TFe > 0 And simContext.HS > 0 And TZn = 0 And TPb > 0 Then 'Fe,Pb sulfide */
    if (TFe > 0 && simContext.HS > 0 && TZn == 0 && TPb > 0) {
        while (fabs((FeOld / (1.0 + ZP1 * simContext.HS) - mc[iFe]) / mc[iFe]) > 0.001
            || fabs((PbOld / (1.0 + ZP5 + ZP6 * pow(simContext.HS, 2) + ZP7 * pow(simContext.HS, 3)) - mc[iPb]) / mc[iPb]) > 0.001
            || fabs((HSOld - simContext.HS) / simContext.HS) > 0.001)
        {
            MeSWhileloop(im, igas, &HSOld, ZP1, ZP2, ZP3, ZP4, ZP5, ZP6, ZP7, ZP9);
        }
    }

    /* 3) If TFe > 0 And simContext.HS > 0 And TZn > 0 And TPb = 0 Then 'Fe,Zn sulfide */
    if (TFe > 0 && simContext.HS > 0 && TZn > 0 && TPb == 0) {
        while (fabs((FeOld / (1.0 + ZP1 * simContext.HS) - mc[iFe]) / mc[iFe]) > 0.001
            || fabs((ZnOld / (1.0 + ZP2 + ZP3 * pow(simContext.HS, 2) + ZP4 * pow(simContext.HS, 3) + ZP9 * pow(OH, 2)) - mc[iZn]) / mc[iZn]) > 0.001
            || fabs((HSOld - simContext.HS) / simContext.HS) > 0.001)
        {
            MeSWhileloop(im, igas, &HSOld, ZP1, ZP2, ZP3, ZP4, ZP5, ZP6, ZP7, ZP9);
        }
    }

    /* 4) If TFe = 0 And simContext.HS > 0 And TZn > 0 And TPb > 0 Then 'Zn,Pb sulfide */
    if (TFe == 0 && simContext.HS > 0 && TZn > 0 && TPb > 0) {
        /* Xin_10/14/2021, correction for Pb speciation */
        while (fabs((PbOld / (1.0 + ZP5 + ZP6 * pow(simContext.HS, 2) + ZP7 * pow(simContext.HS, 3)) - mc[iPb]) / mc[iPb]) > 0.001
            || fabs((ZnOld / (1.0 + ZP2 + ZP3 * pow(simContext.HS, 2) + ZP4 * pow(simContext.HS, 3) + ZP9 * pow(OH, 2)) - mc[iZn]) / mc[iZn]) > 0.001
            || fabs((HSOld - simContext.HS) / simContext.HS) > 0.001)
        {
            MeSWhileloop(im, igas, &HSOld, ZP1, ZP2, ZP3, ZP4, ZP5, ZP6, ZP7, ZP9);
        }
    }

    /* 5) If TFe = 0 And simContext.HS > 0 And TZn > 0 And TPb = 0 Then 'Zn sulfide */
    if (TFe == 0 && simContext.HS > 0 && TZn > 0 && TPb == 0) {
        /* Xin_ZnS 5/05/2021 */
        while (fabs((ZnOld / (1.0 + ZP2 + ZP3 * pow(simContext.HS, 2) + ZP4 * pow(simContext.HS, 3) + ZP9 * pow(OH, 2)) - mc[iZn]) / mc[iZn]) > 0.001
            || fabs((HSOld - simContext.HS) / simContext.HS) > 0.001)
        {
            MeSWhileloop(im, igas, &HSOld, ZP1, ZP2, ZP3, ZP4, ZP5, ZP6, ZP7, ZP9);
        }
    }

    /* Case G: TFe > 0 And simContext.HS > 0 And TZn = 0 And TPb = 0
       这是只有 Fe 存在的情况，需要交替更新 mc(iFe) 与 simContext.HS 直到收敛 */
    if (TFe > 0.0 && simContext.HS > 0.0 && TZn == 0.0 && TPb == 0.0) {
        while (fabs((FeOld - mc[iFe]) / mc[iFe]) > 0.001
            || fabs((HSOld - simContext.HS) / simContext.HS) > 0.001) {
            /* VB: FeOld = mc(iFe): HSOld = simContext.HS: */
            FeOld = mc[iFe];
            HSOld = simContext.HS;

            if (igas == 3) {
                mc[iFe] = (TFe - ppt) / (1.0 + ZP1 * simContext.HS);
                simContext.HS = (TH2Saq - ppt) / (hydHS + ZP1 * mc[iFe]);
            }
            else {
                mc[iFe] = TFe / (1.0 + ZP1 * simContext.HS);
                simContext.HS = TH2Saq / (hydHS + ZP1 * mc[iFe]);
            }
        }
    }

    /* VB 末尾: 1000 fMeSSpeciation = Null */
    return NAN;
}

void C4_SSPEquilCalcs(int ppt_or_not, int im, int igas, double Ksp)
{
    //ppt_or_not：“0”表示不发生沉淀；“1”表示只发生沉淀；“2”表示允许同时发生沉淀和溶解

    int RunPpt = (ppt_or_not == 1) ? 1 : 0;

    ppt = 0.0;
    pHHigh = 14.0;
    pHLow = 0.0;
    double pptOld, Vg_ZRT_Old;
    double s1, s2, s3;
    double coef1, coef2, coef3, coef4;
    double pptNew;
    int ISilicate;
    double fppt, dfppt;

    DoubleVec gNeutAq = { gNeut[iCO2aq], gNeut[iH2Saq] }; // Array for MultiPhaseFlash

    // i=5 pH 不同
    // i=6 mc [7,8,9]不同
    for (int i = 0; i < 30; i++) {
        pH = (pHHigh + pHLow) / 2.0;
        aH = pow(10.0, -pH);
        H = aH / (gCat[iH] * gNCat[iH]);
        OH = KH2O / (aH * gAn[iOH] * gNAn[iOH]);

        ppt = 0.0;
        pptOld = 1e42;
        Vg_ZRT = 0.0;
        Vg_ZRT_Old = 1e42;
        Iteration = 0;

        while (fabs((ppt - pptOld + 1e-21) / (ppt + pptOld + 1e-16)) +
            fabs((Vg_ZRT - Vg_ZRT_Old + 1e-21) / (Vg_ZRT + Vg_ZRT_Old + 1e-16)) > 0.0001 && Iteration < 100) {
            Iteration++;
            pptOld = ppt;
            Vg_ZRT_Old = Vg_ZRT;

            if (useEOS == 0) {
                t1 = gGas[iCH4g] * (RatioOilBPoints * KgoCH4 * Mass_o / gL[iCH4o] +
                    KgwCH4 * mass_w / (gNeut[iCH4aq] * gNNeut[iCH4aq]));
                t2 = gGas[iCO2g] * (RatioOilBPoints * KgoCO2 * Mass_o / gL[iCO2o] +
                    KgwCO2 * mass_w * (1.0 / (gNeut[iCO2aq] * gNNeut[iCO2aq]) +
                        (K1H2CO3 * aH2O) / (aH * gAn[iHCO3] * gNAn[iHCO3]) +
                        (K1H2CO3 * aH2O) * K2HCO3 / (aH * aH * gAn[iCO3] * gNAn[iCO3])));
                t3 = gGas[iH2Sg] * (RatioOilBPoints * KgoH2S * Mass_o / gL[iH2So] +
                    KgwH2S * mass_w * (1.0 / (gNeut[iH2Saq] * gNNeut[iH2Saq]) +
                        K1H2S / (aH * gAn[iHS] * gNAn[iHS]) +
                        K1H2S * K2HS / (aH * aH * gAn[iSion] * gNAn[iSion])));
                s1 = nTCH4 / Ppsia;
                s2 = nTCO2 / Ppsia;
                if (igas == 2) s2 = (nTCO2 - ppt * mass_w) / Ppsia;
                s3 = nTH2S / Ppsia;
                if (igas == 3) s3 = (nTH2S - ppt * mass_w) / Ppsia;
                coef1 = 1.0;
                coef2 = t1 + t2 + t3 - (s1 + s2 + s3);
                coef3 = t1 * t2 + t1 * t3 + t2 * t3 - (s1 * t2 + s1 * t3 + s3 * t1 + s3 * t2 + s2 * t1 + s2 * t3);
                coef4 = t1 * t2 * t3 - (t1 * t2 * s3 + t1 * t3 * s2 + t2 * t3 * s1);
                CubicRoots(coef1, coef2, coef3, coef4, root1, root2, root3);
                Vg_ZRT = root1;
                if (root2 > Vg_ZRT) Vg_ZRT = root2;
                if (root3 > Vg_ZRT) Vg_ZRT = root3;
                if (Vg_ZRT < 0.0) Vg_ZRT = 0.0;
            }

            if (ppt_or_not == 0) goto label223;

            if (igas == 2) {
                if (useEOS == 0) {
                    a = 1.0;
                    b = -(mc[im] * mass_w + nTCO2) / mass_w;
                    cc = (mc[im] * nTCO2 - Ksp * (Vg_ZRT + t2) * (aH * aH) / ((K1H2CO3 * aH2O) * K2HCO3 * KgwCO2 * gGas[iCO2g] * gCat[im] * gNCat[im])) / mass_w;
                }
                else {
                    total_moles_Temp = Total_moles_before_precipitation - ppt;
                    zTemp[0] = z_before_precipitation[0] * Total_moles_before_precipitation / total_moles_Temp;
                    zTemp[1] = (z_before_precipitation[1] * Total_moles_before_precipitation - ppt * mass_w) / total_moles_Temp;
                    if (zTemp[1] < 0.0) {
                        ppt = 0.0;
                        goto label223;
                    }
                    for (int iz = 2; iz < 15; iz++) {
                        zTemp[iz] = z_before_precipitation[iz] * Total_moles_before_precipitation / total_moles_Temp;
                    }
                    if (simContext.RunH2SGUI != 1) {

                        MultiPhaseFlash(mf_ParametersWereRead, mf_TCr, mf_PCr, mf_Omega, mf_MWgas, mf_kPr, mf_c0, mf_c1, TK,
                            PBar, total_moles_Temp, zTemp, gNeutAq, aH2O, density, compositions, phi, Compr, beta, zOutput, mass_phase, MW_Phase, &No_Phases);
                    }
                    else {
                        MultiPhaseFlash_CS(myiosheet, myiocol, eosProps, kij, TK, PBar, total_moles_Temp, zTemp, gNeutAq,
                            aH2O, density, compositions, phi, Compr, beta, zOutput, mass_phase, MW_Phase, &No_Phases);
                    }
                    PCO2 = compositions[1][3] * phi[1][2] * Ppsia / gGas[iCO2g];
                    a = 1.0;
                    b = -(mc[im] * mass_w + nTCO2) / mass_w;
                    cc = (mc[im] * nTCO2 - Ksp * (nTCO2 / PCO2) * (aH * aH) / ((K1H2CO3 * aH2O) * K2HCO3 * KgwCO2 * gGas[iCO2g] * gCat[im] * gNCat[im])) / mass_w;
                }
                QuadraticRoots(a, b, cc, ppt, root2);
            }

            if (igas == 3) {
                if (useEOS == 0) {
                    a = 1.0;
                    b = -(mc[im] * mass_w + nTH2S) / mass_w;
                    cc = (mc[im] * nTH2S - Ksp * (Vg_ZRT + t3) * (aH * aH) / (K1H2S * KgwH2S * gGas[iH2Sg] * gCat[im] * gNCat[im])) / mass_w;
                }
                else {
                    total_moles_Temp = Total_moles_before_precipitation - ppt;
                    for (int iz = 0; iz < 2; iz++) {
                        zTemp[iz] = z_before_precipitation[iz] * Total_moles_before_precipitation / total_moles_Temp;
                    }
                    zTemp[2] = (z_before_precipitation[2] * Total_moles_before_precipitation - ppt * mass_w) / total_moles_Temp;
                    if (zTemp[2] < 0.0) {
                        ppt = 0.0;
                        goto label223;//When ppt is too large, HCS value is too high, i.e. pH is too high
                    }
                    for (int iz = 3; iz < 15; iz++) {
                        zTemp[iz] = z_before_precipitation[iz] * Total_moles_before_precipitation / total_moles_Temp;
                    }
                    if (simContext.RunH2SGUI != 1) {

                        MultiPhaseFlash(mf_ParametersWereRead, mf_TCr, mf_PCr, mf_Omega, mf_MWgas, mf_kPr, mf_c0, mf_c1, TK,
                            PBar, total_moles_Temp, zTemp, gNeutAq, aH2O, density, compositions, phi, Compr, beta, zOutput, mass_phase, MW_Phase, &No_Phases);
                    }
                    else {

                        MultiPhaseFlash_CS(myiosheet, myiocol, eosProps, kij, TK, PBar, total_moles_Temp, zTemp, gNeutAq,
                            aH2O, density, compositions, phi, Compr, beta, zOutput, mass_phase, MW_Phase, &No_Phases);
                    }
                    PH2S = compositions[2][3] * phi[2][2] * Ppsia / gGas[iH2Sg];
                    a = 1.0;
                    b = -(mc[im] * mass_w + nTH2S) / mass_w;
                    cc = (mc[im] * nTH2S - Ksp * (nTH2S / PH2S) * (aH * aH) / (K1H2S * KgwH2S * gGas[iH2Sg] * gCat[im] * gNCat[im])) / mass_w;
                }
                QuadraticRoots(a, b, cc, ppt, root2);
            }

            if (igas == 8) {
                ppt = mc[iCa] - KspCaOH2 / (OH * OH * gCat[iCa] * gAn[iOH] * gAn[iOH]);
            }
            if (igas == 9) {
                ppt = mc[iMg] - KspMgOH2 / (OH * OH * gCat[iMg] * gAn[iOH] * gAn[iOH]);
            }

            if (igas == 12) {
                double kc = KspGreenalite * pow(aH, 6) / (pow(gCat[iFe], 3) * pow(gNeut[iH4SiO4aq], 2));
                if (mc[iFe] / 3.0 >= TH4SiO4 / 2.0) pptNew = (TH4SiO4 - pow(kc / pow(mc[iFe], 3), 0.5)) / 2.0;
                else pptNew = (mc[iFe] - pow(kc / pow(TH4SiO4, 2), 0.3333)) / 3.0;
                ISilicate = 0;
                while (ISilicate < 1000 && fabs((pptNew - ppt) / pptNew) > 0.00001) {
                    ISilicate++;
                    ppt = pptNew;
                    fppt = pow(mc[iFe] - 3.0 * ppt, 3) * pow(TH4SiO4 - 2.0 * ppt, 2) - kc;
                    dfppt = -1.0 * (9.0 * pow(mc[iFe] - 3.0 * ppt, 2) * pow(TH4SiO4 - 2.0 * ppt, 2) + 4.0 * pow(mc[iFe] - 3.0 * ppt, 3) * (TH4SiO4 - 2.0 * ppt));
                    if ((mc[iFe] - 3.0 * ppt) == 0.0 || (TH4SiO4 - 2.0 * ppt) == 0.0) goto label312;
                    pptNew = ppt - fppt / dfppt;
                    if (ISilicate == 999) simContext.errmsg[5] = 6;
                }
            label312:
                TH4SiO4_Greenalite = TH4SiO4 - 2.0 * ppt;
                ppt = ppt * 6.0 / 2.0;
            }

            if (igas == 13) {
                double kc = KspDiopside * pow(aH, 4) / (gCat[iCa] * gCat[iMg] * pow(gNeut[iH4SiO4aq], 2));
                if (mc[iCa] < mc[iMg] && mc[iCa] < TH4SiO4 / 2.0) pptNew = mc[iCa] - kc / (mc[iMg] * pow(TH4SiO4, 2));
                else if (mc[iMg] < mc[iCa] && mc[iMg] < TH4SiO4 / 2.0) pptNew = mc[iMg] - kc / (mc[iCa] * pow(TH4SiO4, 2));
                else pptNew = (TH4SiO4 - pow(kc / (mc[iCa] * mc[iMg]), 0.5)) / 2.0;
                ISilicate = 0;
                while (ISilicate < 1000 && fabs((pptNew - ppt) / pptNew) > 0.00001) {
                    ISilicate++;
                    ppt = pptNew;
                    fppt = (mc[iCa] - ppt) * (mc[iMg] - ppt) * pow(TH4SiO4 - 2.0 * ppt, 2) - kc;
                    dfppt = -1.0 * ((mc[iMg] - ppt) * pow(TH4SiO4 - 2.0 * ppt, 2) + (mc[iCa] - ppt) * pow(TH4SiO4 - 2.0 * ppt, 2) + 4.0 * (mc[iCa] - ppt) * (mc[iMg] - ppt) * (TH4SiO4 - 2.0 * ppt));
                    if ((mc[iMg] - ppt) == 0.0 || (TH4SiO4 - 2.0 * ppt) == 0.0 || (mc[iCa] - ppt) == 0.0) goto label313;
                    pptNew = ppt - fppt / dfppt;
                    if (ISilicate == 999) simContext.errmsg[5] = 6;
                }
            label313:
                ppt = ppt * 4.0 / 2.0;
                TH4SiO4_Dipside = TH4SiO4 - 2.0 * ppt;
            }

            if (igas == 14) {
                double kc = KspChrysotile * pow(aH, 6) / (pow(gCat[iMg], 3) * pow(gNeut[iH4SiO4aq], 2));
                if (mc[iMg] / 3.0 >= TH4SiO4 / 2.0) pptNew = (TH4SiO4 - pow(kc / pow(mc[iMg], 3), 0.5)) / 2.0;
                else pptNew = (mc[iMg] - pow(kc / pow(TH4SiO4, 2), 0.3333)) / 3.0;
                ISilicate = 0;
                while (ISilicate < 1000 && fabs((pptNew - ppt) / pptNew) > 0.00001) {
                    ISilicate++;
                    ppt = pptNew;
                    fppt = pow(mc[iMg] - 3.0 * ppt, 3) * pow(TH4SiO4 - 2.0 * ppt, 2) - kc;
                    if ((mc[iMg] - 3.0 * ppt) == 0.0 || (TH4SiO4 - 2.0 * ppt) == 0.0) goto label314;
                    dfppt = -1.0 * (9.0 * pow(mc[iMg] - 3.0 * ppt, 2) * pow(TH4SiO4 - 2.0 * ppt, 2) + 4.0 * pow(mc[iMg] - 3.0 * ppt, 3) * (TH4SiO4 - 2.0 * ppt));
                    pptNew = ppt - fppt / dfppt;
                    if (ISilicate == 999) simContext.errmsg[5] = 6;
                }
            label314:
                ppt = ppt * 6.0 / 2.0;
                TH4SiO4_Chrysotile = TH4SiO4 - 2.0 * ppt;
            }

            if (ppt_or_not == 1 && ppt < 0.0) ppt = 0.0;
            if (Iteration == 1 && igas > 3) goto label223;
        }
    label223:
        if (useEOS == 0) {
            PCH4 = nTCH4 / (Vg_ZRT + t1);
            PCO2 = nTCO2 / (Vg_ZRT + t2);
            if (igas == 2) PCO2 = (nTCO2 - ppt * mass_w) / (Vg_ZRT + t2);
            PH2S = nTH2S / (Vg_ZRT + t3);
            if (igas == 3) PH2S = (nTH2S - ppt * mass_w) / (Vg_ZRT + t3);
        }
        else {
            PCO2 = compositions[1][3] * phi[1][2] * Ppsia / gGas[iCO2g];
            PH2S = compositions[2][3] * phi[2][2] * Ppsia / gGas[iH2Sg];
        }
        // CO2aq 计算结果不同  PCO2 不同 
        simContext.CO2aq = KgwCO2 * PCO2 * gGas[iCO2g] / (gNeut[iCO2aq] * gNNeut[iCO2aq]);
        simContext.H2Saq = KgwH2S * PH2S * gGas[iH2Sg] / (gNeut[iH2Saq] * gNNeut[iH2Saq]);
        // HCO3计算结果不同 
        simContext.HCO3 = (K1H2CO3 * aH2O) * simContext.CO2aq * gNeut[iCO2aq] * gNNeut[iCO2aq] / (aH * gAn[iHCO3] * gNAn[iHCO3]);
        simContext.CO3 = K2HCO3 * simContext.HCO3 * gAn[iHCO3] * gNAn[iHCO3] / (aH * gAn[iCO3] * gNAn[iCO3]);
        simContext.HS = K1H2S * simContext.H2Saq * gNeut[iH2Saq] * gNNeut[iH2Saq] / (aH * gAn[iHS] * gNAn[iHS]);
        TH2Saq = simContext.H2Saq + simContext.HS;
        hydHS = aH * gAn[iHS] * gNAn[iHS] / (K1H2S * gNeut[iH2Saq] * gNNeut[iH2Saq]) + 1.0;
        if (nTH2S > 0.0) {

            //mt = fMeSSpeciation(im, igas) 'speciate mc(iFe) and simContext.HS
            fMeSSpeciation(im, igas);// 'speciate mc(iFe) and simContext.HS
        }

        S = K2HS * simContext.HS * gAn[iHS] * gNAn[iHS] / (aH * gAn[iSion] * gNAn[iSion]);
        hydAc = aH * gAn[iAc] * gNAn[iAc] / (KHAc * gNeut[iHAcaq] * gNNeut[iHAcaq]) + 1.0;
        AC = TAc / hydAc;
        HAcaq = TAc - AC;
        hydH2BO3 = aH * gAn[iH2BO3] * gNAn[iH2BO3] / (KH3BO3 * gNeut[iH3BO3] * gNNeut[iH3BO3]) + 1.0;
        H2BO3 = TH3BO3 / hydH2BO3;
        hydNH3 = aH * gNeut[iNH3] * gNNeut[iNH3] / (KNH4 * gCat[iNH4] * gNCat[iNH4]) + 1.0;
        NH3 = TNH4 / hydNH3;

        hydH2SiO4 = (aH * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq])) +
            (aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4])) + 1.0;
        H2SiO4 = TH4SiO4 / hydH2SiO4;
        H3SiO4 = H2SiO4 * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]);
        H4SiO4 = H2SiO4 * aH * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]);

        if (igas == 10) {
            H4SiO4 = KspAmSilica / gNeut[iH4SiO4aq];
            H3SiO4 = KH4SiO4 * KspAmSilica / (aH * gAn[iH3SiO4]);
            H2SiO4 = KH4SiO4 * KH3SiO3 * KspAmSilica / (aH * aH * gAn[iH2SiO4]);
            pptAmSilica = TH4SiO4 - (H2SiO4 + H3SiO4 + H4SiO4);
        }
        if (igas == 11) {
            H4SiO4 = KspQuartz;
            H3SiO4 = KH4SiO4 * KspQuartz / aH;
            H2SiO4 = KH4SiO4 * KH3SiO3 * KspQuartz / (aH * aH); // Note: VB has KspAmSilica, probably typo, changed to KspQuartz
        }

        if (igas == 12) {
            H2SiO4 = TH4SiO4_Greenalite / hydH2SiO4;
            H3SiO4 = H2SiO4 * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]);
            H4SiO4 = H2SiO4 * aH * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]);
        }
        if (igas == 13) {
            H2SiO4 = TH4SiO4_Dipside / hydH2SiO4;
            H3SiO4 = H2SiO4 * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]);
            H4SiO4 = H2SiO4 * aH * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]);
        }
        if (igas == 14) {
            H2SiO4 = TH4SiO4_Chrysotile / hydH2SiO4;
            H3SiO4 = H2SiO4 * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]);
            H4SiO4 = H2SiO4 * aH * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]);
        }
        // faH计算不同 HCO3不同 CO3 不同
        faH = Alk - 2.0 * ppt - (simContext.HCO3 + 2.0 * simContext.CO3 + simContext.HS + 2.0 * S + AC + NH3 + H2BO3 + H3SiO4 + 2.0 * H2SiO4 + OH - H);
        ypH[i] = faH;
        xpH[i] = pH;
        xppt[i] = ppt;
        if (faH > 0.0) pHLow = pH;
        if (faH < 0.0) pHHigh = pH;
    }

    // Newton-Raphson
    if (fabs(ypH[29] / ((ypH[29] - ypH[28]) / (xpH[29] - xpH[28]))) < 0.001) {
        pH = pH - ypH[29] / ((ypH[29] - ypH[28]) / (xpH[29] - xpH[28]));
    }
    else {
        pH = fabs((xpH[29] + xpH[28])) / 2.0;
        if (RunPpt == 0) simContext.errmsg[4] = 5; // VB simContext.errmsg(5)=5, 0-based [4]
    }

    aH = pow(10.0, -pH);
    H = aH / (gCat[iH] * gNCat[iH]);
    OH = KH2O / (aH * gAn[iOH] * gNAn[iOH]);

    ppt = 0.0;
    pptOld = 1e42;
    Vg_ZRT = 0.0;
    Vg_ZRT_Old = 1e42;
    Iteration = 0;

    while (fabs((ppt - pptOld + 1e-21) / (ppt + pptOld + 1e-16)) +
        fabs((Vg_ZRT - Vg_ZRT_Old + 1e-21) / (Vg_ZRT + Vg_ZRT_Old + 1e-16)) > 0.0001 && Iteration < 100) {
        Iteration++;
        pptOld = ppt;
        Vg_ZRT_Old = Vg_ZRT;

        if (useEOS == 0) {
            t1 = gGas[iCH4g] * (RatioOilBPoints * KgoCH4 * Mass_o / gL[iCH4o] +
                KgwCH4 * mass_w / (gNeut[iCH4aq] * gNNeut[iCH4aq]));
            t2 = gGas[iCO2g] * (RatioOilBPoints * KgoCO2 * Mass_o / gL[iCO2o] +
                KgwCO2 * mass_w * (1.0 / (gNeut[iCO2aq] * gNNeut[iCO2aq]) +
                    (K1H2CO3 * aH2O) / (aH * gAn[iHCO3] * gNAn[iHCO3]) +
                    (K1H2CO3 * aH2O) * K2HCO3 / (aH * aH * gAn[iCO3] * gNAn[iCO3])));
            t3 = gGas[iH2Sg] * (RatioOilBPoints * KgoH2S * Mass_o / gL[iH2So] +
                KgwH2S * mass_w * (1.0 / (gNeut[iH2Saq] * gNNeut[iH2Saq]) +
                    K1H2S / (aH * gAn[iHS] * gNAn[iHS]) +
                    K1H2S * K2HS / (aH * aH * gAn[iSion] * gNAn[iSion])));
            s1 = nTCH4 / Ppsia;
            s2 = nTCO2 / Ppsia;
            if (igas == 2) s2 = (nTCO2 - ppt * mass_w) / Ppsia;
            s3 = nTH2S / Ppsia;
            if (igas == 3) s3 = (nTH2S - ppt * mass_w) / Ppsia;
            coef1 = 1.0;
            coef2 = t1 + t2 + t3 - (s1 + s2 + s3);
            coef3 = t1 * t2 + t1 * t3 + t2 * t3 - (s1 * t2 + s1 * t3 + s3 * t1 + s3 * t2 + s2 * t1 + s2 * t3);
            coef4 = t1 * t2 * t3 - (t1 * t2 * s3 + t1 * t3 * s2 + t2 * t3 * s1);
            CubicRoots(coef1, coef2, coef3, coef4, root1, root2, root3);
            Vg_ZRT = root1;
            if (root2 > Vg_ZRT) Vg_ZRT = root2;
            if (root3 > Vg_ZRT) Vg_ZRT = root3;
            if (Vg_ZRT < 0.0) Vg_ZRT = 0.0;
        }

        if (ppt_or_not == 0) goto label323;

        if (igas == 2) {
            if (useEOS == 0) {
                a = 1.0;
                b = -(mc[im] * mass_w + nTCO2) / mass_w;
                cc = (mc[im] * nTCO2 - Ksp * (Vg_ZRT + t2) * (aH * aH) / ((K1H2CO3 * aH2O) * K2HCO3 * KgwCO2 * gGas[iCO2g] * gCat[im] * gNCat[im])) / mass_w;
            }
            else {
                total_moles_Temp = Total_moles_before_precipitation - ppt;
                zTemp[0] = z_before_precipitation[0] * Total_moles_before_precipitation / total_moles_Temp;
                zTemp[1] = (z_before_precipitation[1] * Total_moles_before_precipitation - ppt * mass_w) / total_moles_Temp;
                for (int iz = 2; iz < 15; iz++) {
                    zTemp[iz] = z_before_precipitation[iz] * Total_moles_before_precipitation / total_moles_Temp;
                }
                if (simContext.RunH2SGUI != 1) {
                    MultiPhaseFlash(mf_ParametersWereRead, mf_TCr, mf_PCr, mf_Omega, mf_MWgas, mf_kPr, mf_c0, mf_c1, TK,
                        PBar, total_moles_Temp, zTemp, gNeutAq, aH2O, density, compositions, phi, Compr, beta, zOutput, mass_phase, MW_Phase, &No_Phases);
                }
                else {

                    MultiPhaseFlash_CS(myiosheet, myiocol, eosProps, kij, TK, PBar, total_moles_Temp, zTemp, gNeutAq,
                        aH2O, density, compositions, phi, Compr, beta, zOutput, mass_phase, MW_Phase, &No_Phases);
                }
                simContext.CO2aq = compositions[1][3] * beta[2] * total_moles_Temp / mass_w / (gNeut[iCO2aq] * gNNeut[iCO2aq]); // Assume beta[2] water?
                a = 1.0;
                b = -(mc[im] + simContext.CO2aq * gNeut[iCO2aq] * gNNeut[iCO2aq] * (K1H2CO3 * aH2O) / (aH * gAn[iHCO3] * gNAn[iHCO3]));
                cc = (mc[im] * simContext.CO2aq * gNeut[iCO2aq] * gNNeut[iCO2aq] * (K1H2CO3 * aH2O) / (aH * gAn[iHCO3] * gNAn[iHCO3]) - Ksp * aH / K2HCO3);
            }
            QuadraticRoots(a, b, cc, ppt, root2);
        }

        if (igas == 3) {
            if (useEOS == 0) {
                a = 1.0;
                b = -(mc[im] * mass_w + nTH2S) / mass_w;
                cc = (mc[im] * nTH2S - Ksp * (Vg_ZRT + t3) * (aH * aH) / (K1H2S * KgwH2S * gGas[iH2Sg] * gCat[im] * gNCat[im])) / mass_w;
            }
            else {
                total_moles_Temp = Total_moles_before_precipitation - ppt;
                for (int iz = 0; iz < 2; iz++) {
                    zTemp[iz] = z_before_precipitation[iz] * Total_moles_before_precipitation / total_moles_Temp;
                }
                zTemp[2] = (z_before_precipitation[2] * Total_moles_before_precipitation - ppt * mass_w) / total_moles_Temp;
                for (int iz = 3; iz < 15; iz++) {
                    zTemp[iz] = z_before_precipitation[iz] * Total_moles_before_precipitation / total_moles_Temp;
                }
                if (simContext.RunH2SGUI != 1) {

                    MultiPhaseFlash(mf_ParametersWereRead, mf_TCr, mf_PCr, mf_Omega, mf_MWgas, mf_kPr, mf_c0, mf_c1, TK,
                        PBar, total_moles_Temp, zTemp, gNeutAq, aH2O, density, compositions, phi, Compr, beta, zOutput, mass_phase, MW_Phase, &No_Phases);
                }
                else {

                    MultiPhaseFlash_CS(myiosheet, myiocol, eosProps, kij, TK, PBar, total_moles_Temp, zTemp, gNeutAq,
                        aH2O, density, compositions, phi, Compr, beta, zOutput, mass_phase, MW_Phase, &No_Phases);
                }
                simContext.H2Saq = compositions[2][3] * beta[2] * total_moles_Temp / mass_w / gNeut[iH2Saq] / gNNeut[iH2Saq];
                a = 1.0;
                b = -(mc[im] + simContext.H2Saq * gNeut[iH2Saq] * gNNeut[iH2Saq] * K1H2S / (aH * gAn[iHS] * gNAn[iHS]));
                cc = (mc[im] * simContext.H2Saq * gNeut[iH2Saq] * gNNeut[iH2Saq] * K1H2S / (aH * gAn[iHS] * gNAn[iHS]) - Ksp * aH);
            }
            QuadraticRoots(a, b, cc, ppt, root2);
        }

        if (igas == 8) {
            ppt = mc[iCa] - KspCaOH2 / (OH * OH * gCat[iCa] * gAn[iOH] * gAn[iOH]);
        }
        if (igas == 9) {
            ppt = mc[iMg] - KspMgOH2 / (OH * OH * gCat[iMg] * gAn[iOH] * gAn[iOH]);
        }

        if (igas == 12) {
            double kc = KspGreenalite * pow(aH, 6) / (pow(gCat[iFe], 3) * pow(gNeut[iH4SiO4aq], 2));
            if (mc[iFe] / 3.0 >= TH4SiO4 / 2.0) pptNew = (TH4SiO4 - pow(kc / pow(mc[iFe], 3), 0.5)) / 2.0;
            else pptNew = (mc[iFe] - pow(kc / pow(TH4SiO4, 2), 0.3333)) / 3.0;
            ISilicate = 0;
            while (ISilicate < 1000 && fabs((pptNew - ppt) / pptNew) > 0.00001) {
                ISilicate++;
                ppt = pptNew;
                fppt = pow(mc[iFe] - 3.0 * ppt, 3) * pow(TH4SiO4 - 2.0 * ppt, 2) - kc;
                dfppt = -1.0 * (9.0 * pow(mc[iFe] - 3.0 * ppt, 2) * pow(TH4SiO4 - 2.0 * ppt, 2) + 4.0 * pow(mc[iFe] - 3.0 * ppt, 3) * (TH4SiO4 - 2.0 * ppt));
                if ((mc[iFe] - 3.0 * ppt) == 0.0 || (TH4SiO4 - 2.0 * ppt) == 0.0) goto label412;
                pptNew = ppt - fppt / dfppt;
                if (ISilicate == 999) simContext.errmsg[5] = 6;
            }
        label412:
            TH4SiO4_Greenalite = TH4SiO4 - 2.0 * ppt;
            ppt = ppt * 6.0 / 2.0;
        }

        if (igas == 13) {
            double kc = KspDiopside * pow(aH, 4) / (gCat[iCa] * gCat[iMg] * pow(gNeut[iH4SiO4aq], 2));
            if (mc[iCa] < mc[iMg] && mc[iCa] < TH4SiO4 / 2.0) pptNew = mc[iCa] - kc / (mc[iMg] * pow(TH4SiO4, 2));
            else if (mc[iMg] < mc[iCa] && mc[iMg] < TH4SiO4 / 2.0) pptNew = mc[iMg] - kc / (mc[iCa] * pow(TH4SiO4, 2));
            else pptNew = (TH4SiO4 - pow(kc / (mc[iCa] * mc[iMg]), 0.5)) / 2.0;
            ISilicate = 0;
            while (ISilicate < 1000 && fabs((pptNew - ppt) / pptNew) > 0.00001) {
                ISilicate++;
                ppt = pptNew;
                fppt = (mc[iCa] - ppt) * (mc[iMg] - ppt) * pow(TH4SiO4 - 2.0 * ppt, 2) - kc;
                dfppt = -1.0 * ((mc[iMg] - ppt) * pow(TH4SiO4 - 2.0 * ppt, 2) + (mc[iCa] - ppt) * pow(TH4SiO4 - 2.0 * ppt, 2) + 4.0 * (mc[iCa] - ppt) * (mc[iMg] - ppt) * (TH4SiO4 - 2.0 * ppt));
                if ((mc[iMg] - ppt) == 0.0 || (TH4SiO4 - 2.0 * ppt) == 0.0 || (mc[iCa] - ppt) == 0.0) goto label413;
                pptNew = ppt - fppt / dfppt;
                if (ISilicate == 999) simContext.errmsg[5] = 6;
            }
        label413:
            ppt = ppt * 4.0 / 2.0;
            TH4SiO4_Dipside = TH4SiO4 - 2.0 * ppt;
        }

        if (igas == 14) {
            double kc = KspChrysotile * pow(aH, 6) / (pow(gCat[iMg], 3) * pow(gNeut[iH4SiO4aq], 2));
            if (mc[iMg] / 3.0 >= TH4SiO4 / 2.0) pptNew = (TH4SiO4 - pow(kc / pow(mc[iMg], 3), 0.5)) / 2.0;
            else pptNew = (mc[iMg] - pow(kc / pow(TH4SiO4, 2), 0.3333)) / 3.0;
            ISilicate = 0;
            while (ISilicate < 1000 && fabs((pptNew - ppt) / pptNew) > 0.00001) {
                ISilicate++;
                ppt = pptNew;
                fppt = pow(mc[iMg] - 3.0 * ppt, 3) * pow(TH4SiO4 - 2.0 * ppt, 2) - kc;
                if ((mc[iMg] - 3.0 * ppt) == 0.0 || (TH4SiO4 - 2.0 * ppt) == 0.0) goto label414;
                dfppt = -1.0 * (9.0 * pow(mc[iMg] - 3.0 * ppt, 2) * pow(TH4SiO4 - 2.0 * ppt, 2) + 4.0 * pow(mc[iMg] - 3.0 * ppt, 3) * (TH4SiO4 - 2.0 * ppt));
                pptNew = ppt - fppt / dfppt;
                if (ISilicate == 999) simContext.errmsg[5] = 6;
            }
        label414:
            ppt = ppt * 6.0 / 2.0;
            TH4SiO4_Chrysotile = TH4SiO4 - 2.0 * ppt;
        }

        if (ppt_or_not == 1 && ppt < 0.0) ppt = 0.0;
        if (Iteration == 1 && igas > 3) goto label323;
    }
label323:
    if (useEOS == 0) {
        PCH4 = nTCH4 / (Vg_ZRT + t1);
        PCO2 = nTCO2 / (Vg_ZRT + t2);
        if (igas == 2) PCO2 = (nTCO2 - ppt * mass_w) / (Vg_ZRT + t2);
        PH2S = nTH2S / (Vg_ZRT + t3);
        if (igas == 3) PH2S = (nTH2S - ppt * mass_w) / (Vg_ZRT + t3);
    }
    else {
        PCO2 = compositions[1][3] * phi[1][2] * Ppsia / gGas[iCO2g];
        PH2S = compositions[2][3] * phi[2][2] * Ppsia / gGas[iH2Sg];
    }
    simContext.CO2aq = KgwCO2 * PCO2 * gGas[iCO2g] / (gNeut[iCO2aq] * gNNeut[iCO2aq]);
    simContext.H2Saq = KgwH2S * PH2S * gGas[iH2Sg] / (gNeut[iH2Saq] * gNNeut[iH2Saq]);
    simContext.HCO3 = (K1H2CO3 * aH2O) * simContext.CO2aq * gNeut[iCO2aq] * gNNeut[iCO2aq] / (aH * gAn[iHCO3] * gNAn[iHCO3]);
    simContext.CO3 = K2HCO3 * simContext.HCO3 * gAn[iHCO3] * gNAn[iHCO3] / (aH * gAn[iCO3] * gNAn[iCO3]);
    simContext.HS = K1H2S * simContext.H2Saq * gNeut[iH2Saq] * gNNeut[iH2Saq] / (aH * gAn[iHS] * gNAn[iHS]);

    //if (TH2Saq > 0.0) mt = fMeSSpeciation(im, igas, TZn, TPb, hydHS, *ppt, &root1, &root2, &root3, BetaDot, gDot, KstFeSaq, aH);
    if (TH2Saq > 0.0) fMeSSpeciation(im, igas);//output mc(iFe),mc(iZn), mc(iPb), simContext.HS

    mn[iFeSaq] = KstFeSaq * mc[iFe] * simContext.HS * gAn[iHS] * gNAn[iHS] * gCat[iFe] * gNCat[iFe] / (gNeut[iFeSaq] * gNNeut[iFeSaq] * aH);
    S = K2HS * simContext.HS * gAn[iHS] * gNAn[iHS] / (aH * gAn[iSion] * gNAn[iSion]);
    hydAc = aH * gAn[iAc] * gNAn[iAc] / (KHAc * gNeut[iHAcaq] * gNNeut[iHAcaq]) + 1.0;
    AC = TAc / hydAc;
    HAcaq = TAc - AC;
    hydH2BO3 = aH * gAn[iH2BO3] * gNAn[iH2BO3] / (KH3BO3 * gNeut[iH3BO3] * gNNeut[iH3BO3]) + 1.0;
    H2BO3 = TH3BO3 / hydH2BO3;
    hydNH3 = aH * gNeut[iNH3] * gNNeut[iNH3] / (KNH4 * gCat[iNH4] * gNCat[iNH4]) + 1.0;
    NH3 = TNH4 / hydNH3;

    hydH2SiO4 = (aH * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq])) +
        (aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4])) + 1.0;
    H2SiO4 = TH4SiO4 / hydH2SiO4;
    H3SiO4 = H2SiO4 * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]);
    H4SiO4 = H2SiO4 * aH * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]);

    if (igas == 10) {
        H4SiO4 = KspAmSilica / gNeut[iH4SiO4aq];
        H3SiO4 = KH4SiO4 * KspAmSilica / (aH * gAn[iH3SiO4]);
        H2SiO4 = KH4SiO4 * KH3SiO3 * KspAmSilica / (aH * aH * gAn[iH2SiO4]);
        pptAmSilica = TH4SiO4 - (H2SiO4 + H3SiO4 + H4SiO4);
    }
    if (igas == 11) {
        H4SiO4 = KspQuartz;
        H3SiO4 = KH4SiO4 * KspQuartz / aH;
        H2SiO4 = KH4SiO4 * KH3SiO3 * KspAmSilica / (aH * aH); // VB has KspAmSilica, perhaps typo
    }

    if (igas == 12) {
        H2SiO4 = TH4SiO4_Greenalite / hydH2SiO4;
        H3SiO4 = H2SiO4 * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]);
        H4SiO4 = H2SiO4 * aH * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]);
    }
    if (igas == 13) {
        H2SiO4 = TH4SiO4_Dipside / hydH2SiO4;
        H3SiO4 = H2SiO4 * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]);
        H4SiO4 = H2SiO4 * aH * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]);
    }
    if (igas == 14) {
        H2SiO4 = TH4SiO4_Chrysotile / hydH2SiO4;
        H3SiO4 = H2SiO4 * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]);
        H4SiO4 = H2SiO4 * aH * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] / (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]);
    }

    faH = Alk - 2.0 * ppt - (simContext.HCO3 + 2.0 * simContext.CO3 + simContext.HS + 2.0 * S + AC + NH3 + H2BO3 + H3SiO4 + 2.0 * H2SiO4 + OH - H);

    if (igas == 10 || igas == 11) {
        ppt = TH4SiO4 - (H2SiO4 + H3SiO4 + H4SiO4);
        if (ppt < 0.0) ppt = 0.0;
    }
}



void QualityControlCalculations(int kk, int j)
{
    //double root1 = 0, root2 = 0, root3 = 0;

    double pHMeterReading_from_QC = 0;
    double Alk_from_QC, NaQC;

    if (simContext.RunStat == 1) {
        /*
                AlkMix(kk) = Worksheets(mySheet).Cells(24, j + 2).Value / (61019 * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#))
            AlkMix(kk) = AlkMix(kk) + 2 * Worksheets(mySheet).Cells(25, j + 2).Value / (60019 * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#))
            TAcMix(kk) = Worksheets(mySheet).Cells(30, j + 2).Value / 59.044 'convert to sum of carboxylic acid in meq/L
            TAcMix(kk) = TAcMix(kk) + Worksheets(mySheet).Cells(31, j + 2).Value / 73.07
            TAcMix(kk) = TAcMix(kk) + Worksheets(mySheet).Cells(32, j + 2).Value / 87.098
            TAcMix(kk) = TAcMix(kk) + Worksheets(mySheet).Cells(33, j + 2).Value / 87.11
            TAcMix(kk) = TAcMix(kk) + Worksheets(mySheet).Cells(34, j + 2).Value / 101.13
            TAcMix(kk) = TAcMix(kk) + Worksheets(mySheet).Cells(35, j + 2).Value / 101.13
            TAcMix(kk) = TAcMix(kk) + Worksheets(mySheet).Cells(36, j + 2).Value / 115.16
            TAcMix(kk) = TAcMix(kk) + Worksheets(mySheet).Cells(37, j + 2).Value / 115.16
            TAcMix(kk) = TAcMix(kk) + Worksheets(mySheet).Cells(38, j + 2).Value / 129.178
            TAcMix(kk) = TAcMix(kk) * 59.044 'convert to mg/L as acetate
            TAcMix(kk) = TAcMix(kk) / (59044 * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#))
        */
    }
    else {
        if (simContext.UseMolal == 0) {
            /*
            AlkMix(kk) = Worksheets(mySheet).Cells(24, j + 2).Value / (61019 * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#))
            AlkMix(kk) = AlkMix(kk) + 2 * Worksheets(mySheet).Cells(25, j + 2).Value / (60019 * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#))
            AlkMix(kk) = AlkMix(kk) + Worksheets(mySheet).Cells(48, j + 2).Value / (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#)
            AlkMix(kk) = AlkMix(kk) - Worksheets(mySheet).Cells(47, j + 2).Value / (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#)
            TAcMix(kk) = Worksheets(mySheet).Cells(26, j + 2).Value / (59044 * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#))
            */
        }
        else {
            /*
            AlkMix(kk) = Worksheets(mySheet).Cells(24, j + 2).Value
            AlkMix(kk) = AlkMix(kk) + 2 * Worksheets(mySheet).Cells(25, j + 2).Value
            AlkMix(kk) = AlkMix(kk) + Worksheets(mySheet).Cells(48, j + 2).Value
            AlkMix(kk) = AlkMix(kk) - Worksheets(mySheet).Cells(47, j + 2).Value
            TAcMix(kk) = Worksheets(mySheet).Cells(26, j + 2).Value
            */
        }

    }
    //double UseH2Sgas;
    UseH2Sgas = simContext.UseH2SgasMix[kk];
    SumOfCations = 0.0000001;
    for (int c = 0; c < NumCat; c++)
        SumOfCations += ChCat[c] * mc[c];

    SumOfAnions = -0.0000001;
    for (int a = 0; a < NumAn; a++)
        SumOfAnions += ChAn[a] * ma[a];

    //1、使用测量的 pH 值和 Alk 值来计算输入表 QC 部分中的 P-CO2。
    pH = simContext.pHMeterStpMix[kk] + DpHj;
    aH = pow(10, -pH);


    if (simContext.RunStat == 1) {
        //simContext.yCO2 = Worksheets(mySheet).Cells(26, j + 2).Value / 100
        //TH2Saq = Worksheets(mySheet).Cells(27, j + 2).Value / (34080 * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#))
        simContext.yH2S = 0;
    }
    else {
        // simContext.yCO2 = Worksheets(mySheet).Cells(31, j + 2).Value / 100

        if (simContext.UseH2SgasMix[kk] == 1) {
            //simContext.yH2S = Worksheets(mySheet).Cells(33, j + 2).Value / 100
            TH2Saq = 0;
        }
        if (simContext.UseH2SgasMix[kk] == 0) {
            if (simContext.UseMolal == 0) {
                //TH2Saq = Worksheets(mySheet).Cells(33, j + 2).Value / (34080 * (rho25CMix[kk] - CalculatedTDSMix[kk] / 1000000.0));
                simContext.yH2S = 0;
            }
            else {
                //TH2Saq = Worksheets(mySheet).Cells(33, j + 2).Value
                simContext.yH2S = 0;
            }
        }

    }


    if (useEOS == 0) {
        if (UseH2Sgas == 0) {
            // Calculate TH2Saq from yH2Sstp and pH. If TH2Saq is given, use it, otherwise use PH2S.
            simContext.HS = TH2Saq / (aH * gAn[iHS] * gNAn[iHS] / (K1H2S * gNeut[iH2Saq] * gNNeut[iH2Saq]) + 1);

            if (TH2Saq > 0) {
                // speciation for simContext.HS
                // fMeSSpeciation(int im, int igas,double TZn, double TPb, double hydHS, double ppt, double* root1, double* root2, double* root3, double* BetaDot, double* gDot, double KstFeSaq, double aH)
                //mt = fMeSSpeciation(2, 2) 'speciation for simContext.HS
                fMeSSpeciation(2, 2);
                simContext.H2Saq = aH * simContext.HS * gAn[iHS] * gNAn[iHS] / (K1H2S * gNeut[iH2Saq] * gNNeut[iH2Saq]);
                simContext.yH2S = simContext.H2Saq * gNeut[iH2Saq] * gNNeut[iH2Saq] / (KgwH2S * Ppsia * gGas[iH2Sg]);

                if (simContext.yH2S > 1) {
                    simContext.errmsg[2] = 3; // 数组下标从0开始，3对应索引2
                }
            }
            else {
                simContext.yH2S = 0;
            }
        }

        // Calculate the P-H2S, or simContext.yH2S, from pH and TH2Saq.
        if (UseH2Sgas == 1) {
            if (simContext.yH2S > 0) {
                simContext.H2Saq = KgwH2S * Ppsia * simContext.yH2S * gGas[iH2Sg] / gNeut[iH2Saq] / gNNeut[iH2Saq];
                simContext.HS = K1H2S * simContext.H2Saq * gNeut[iH2Saq] * gNNeut[iH2Saq] / (aH * gAn[iHS] * gNAn[iHS]);

                // 计算各种配合物形成的ZP参数
                double ZP1 = (KstFeSaq * gAn[iHS] * gNAn[iHS] * gCat[iFe] * gNCat[iFe]) /
                    (gNeut[iFeSaq] * gNNeut[iFeSaq] * aH);

                double ZP2 = BetaDot[iZnCl] * (gDot[iClDot] * ma[iCl]) * gDot[iZnDot] / gDot[iZnCl] +
                    BetaDot[iZnCl2] * pow((gDot[iClDot] * ma[iCl]), 2) * gDot[iZnDot] / gDot[iZnCl2] +
                    BetaDot[iZnCl3] * pow((gDot[iClDot] * ma[iCl]), 3) * gDot[iZnDot] / gDot[iZnCl3] +
                    BetaDot[iZnCl4] * pow((gDot[iClDot] * ma[iCl]), 4) * gDot[iZnDot] / gDot[iZnCl4];

                double ZP3 = BetaDot[iZnHS2] * gDot[iZnDot] * pow(gDot[iHSDot], 2) / gDot[iZnHS2];
                double ZP4 = BetaDot[iZnHS3] * gDot[iZnDot] * pow(gDot[iHSDot], 3) / gDot[iZnHS3];

                double ZP5 = BetaDot[iPbCl] * (gDot[iClDot] * ma[iCl]) * gDot[iPbDot] / gDot[iPbCl] +
                    BetaDot[iPbCl2] * pow((gDot[iClDot] * ma[iCl]), 2) * gDot[iPbDot] / gDot[iPbCl2] +
                    BetaDot[iPbCl3] * pow((gDot[iClDot] * ma[iCl]), 3) * gDot[iPbDot] / gDot[iPbCl3] +
                    BetaDot[iPbCl4] * pow((gDot[iClDot] * ma[iCl]), 4) * gDot[iPbDot] / gDot[iPbCl4];

                double ZP6 = BetaDot[iPbHS2] * gDot[iPbDot] * pow(gDot[iHSDot], 2) / gDot[iPbHS2];
                double ZP7 = BetaDot[iPbHS3] * gDot[iPbDot] * pow(gDot[iHSDot], 3) / gDot[iPbHS3];

                // 计算金属离子浓度，考虑配合物形成
                mc[iFe] = TFe / (1 + ZP1 * simContext.HS);
                mc[iZn] = TZn / (1 + ZP2 + ZP3 * pow(simContext.HS, 2) + ZP4 * pow(simContext.HS, 3));
                mc[iPb] = TPb / (1 + ZP5 + ZP6 * pow(simContext.HS, 2) + ZP7 * pow(simContext.HS, 3));

                // 计算总H2S浓度
                TH2Saq = simContext.H2Saq + (simContext.HS + mc[iFe] * ZP1 * simContext.HS +
                    mc[iZn] * (ZP3 * pow(simContext.HS, 2) + ZP4 * pow(simContext.HS, 3)) +
                    mc[iPb] * (ZP6 * pow(simContext.HS, 2) + ZP7 * pow(simContext.HS, 3)));
            }
            else {
                TH2Saq = 0;
            }
        }

        if (TH2Saq == 0 && simContext.yH2S == 0) {
            simContext.HS = 0;
        }

        // 计算其他化学物种的浓度
        H = aH / gCat[iH] / gNCat[iH];
        OH = KH2O / (aH * gAn[iOH] * gNAn[iOH]);

        AC = TAc / (aH * gAn[iAc] * gNAn[iAc] / (KHAc * gNeut[iHAcaq] * gNNeut[iHAcaq]) + 1); // Note gNNeut(iHACaq)=1

        double hydH2BO3 = aH * gAn[iH2BO3] * gNAn[iH2BO3] / (KH3BO3 * gNeut[iH3BO3] * gNNeut[iH3BO3]) + 1;
        H2BO3 = TH3BO3 / hydH2BO3;

        double hydNH3 = aH * gNeut[iNH3] * gNNeut[iNH3] / (KNH4 * gCat[iNH4] * gNCat[iNH4]) + 1;
        NH3 = TNH4 / hydNH3;

        double hydH2SiO4 = pow(aH, 2) * gAn[iH2SiO4] * gNAn[iH2SiO4] /
            (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]) +
            aH * gAn[iH2SiO4] * gNAn[iH2SiO4] /
            (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]) + 1;
        H2SiO4 = TH4SiO4 / hydH2SiO4;
        H3SiO4 = H2SiO4 * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] /
            (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]);
        H4SiO4 = H2SiO4 * pow(aH, 2) * gAn[iH2SiO4] * gNAn[iH2SiO4] /
            (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]);

        // 计算CO2相关的参数
        double tHCO3 = (K1H2CO3 * aH2O) * KgwCO2 * Ppsia * gGas[iCO2g] /
            (aH * gAn[iHCO3] * gNAn[iHCO3]);
        double tCO3 = (K1H2CO3 * aH2O) * K2HCO3 * KgwCO2 * Ppsia * gGas[iCO2g] /
            (pow(aH, 2) * gAn[iCO3] * gNAn[iCO3]);

        // 计算yCO2
        simContext.yCO2 = (Alk + H - AC - NH3 - H2BO3 - simContext.HS - H3SiO4 - 2 * H2SiO4 - OH) /
            (tHCO3 + 2 * tCO3);

        simContext.HCO3 = (K1H2CO3 * aH2O) * KgwCO2 * Ppsia * simContext.yCO2 * gGas[iCO2g] /
            (aH * gAn[iHCO3] * gNAn[iHCO3]);
        simContext.CO3 = (K1H2CO3 * aH2O) * K2HCO3 * KgwCO2 * Ppsia * simContext.yCO2 * gGas[iCO2g] /
            (pow(aH, 2) * gAn[iCO3] * gNAn[iCO3]);

        if (simContext.yCO2 > 1) {
            simContext.errmsg[0] = 1;
        }

        if (simContext.yCO2 < 0) {
            simContext.errmsg[1] = 2; // 数组下标从0开始，2对应索引1
            simContext.yCO2 = 0;  // Caps %CO2 at 0
        }

        // 输出结果到相应的工作表
        //源码：If simContext.RunStat = 1 And CaseCount(1) = 1
        if (simContext.RunStat == 1 && simContext.CaseCount[0] == 1) {
            // for StatQC produced water
            //Worksheets("Output data sheet").Cells(6, 4) = simContext.yCO2 * 100  '假设 pH 和 Alk 正确，QC 检查 PCO2 是多少？

        }
        else if (simContext.RunStat == 1 && simContext.CaseCount[0] == 9) {
            // for StatQC fresh water
            //Worksheets("Output data sheet").Cells(18, 4) = simContext.yCO2 * 100  'QC check for assuming that pH and Alk are correct, what would PCO2 be?
        }
        else if (simContext.RunH2SGUI == 1) {
            //Worksheets(mySheet).Cells(83, j + 2) = simContext.yH2S * 100

            if (simContext.UseMolal == 0) {
                //Worksheets(mySheet).Cells(84, j + 2) = TH2Saq * (34080 * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#))
            }
            else {
                //Worksheets(mySheet).Cells(84, j + 2) = TH2Saq
            }

            //Worksheets(mySheet).Cells(86, j + 2) = simContext.yCO2 * 100 'QC check for assuming that pH and Alk are correct, what would PCO2 be?
        }
        else {
            //源码的kk从1开始，[1, nob_Input]，对应c语言的kk < nob_Input [0,nob_Input-1] = [0,nob_Input)
            if (kk < simContext.nob_Input) {
                //Worksheets("Input").Cells(14, j + 15) = simContext.yH2S * 100

                if (simContext.UseMolal == 0) {
                    //Worksheets("Input").Cells(15, j + 15) = TH2Saq * (34080 * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#))
                }
                else {
                    //Worksheets("Input").Cells(15, j + 15) = TH2Saq
                }

                //Worksheets("Input").Cells(17, j + 15) = simContext.yCO2 * 100  'QC check for assuming that pH and Alk are correct, what would PCO2 be?
            }
            //注意下方源码： If kk > nob_Input And kk <= nob_Input + nob_InputII Then
            if (kk >= simContext.nob_Input && kk < simContext.nob_Input + simContext.nob_InputII) {
                // Worksheets("Input II").Cells(83, j + 2) = simContext.yH2S * 100

                if (simContext.UseMolal == 0) {
                    //Worksheets("Input II").Cells(84, j + 2) = TH2Saq * (34080 * (rho25CMix(i + nob_Input) - CalculatedTDSMix(i + nob_Input) / 1000000#))
                }
                else {
                    //Worksheets("Input II").Cells(84, j + 2) = TH2Saq
                }

                //Worksheets("Input II").Cells(86, j + 2) = simContext.yCO2 * 100 'QC check for assuming that pH and Alk are correct, what would PCO2 be?
            }
        }


        //2、使用 STP P-CO2 和碱度来计算 pH 值。
        if (UseH2Sgas == 0) {
            // 然后，必须指定yCO2stp。需要计算pH值和yH2Sstp。

            // 读取yCO2值（从相应的单元格）
            if (simContext.RunStat == 1) {
                // simContext.yCO2 = Worksheets(mySheet).Cells(26, j + 2).Value / 100
            }
            else {
                //simContext.yCO2 = Worksheets(mySheet).Cells(31, j + 2).Value / 100
            }

            // 半区间搜索根查找器对于 pH 方程非常有效。
            pHHigh = 14.0;
            pHLow = 0.0;
            //double hydHS, hydAc, hydH2BO3, hydNH3, hydH2SiO4;
            int k;

            // 二分法求解pH值（30次迭代）
            for (k = 0; k < 30; k++) {
                pH = (pHHigh + pHLow) / 2.0;
                aH = pow(10.0, -pH);

                // 计算H+和OH-浓度
                H = aH / (gCat[iH] * gNCat[iH]);
                OH = KH2O / (aH * gAn[iOH] * gNAn[iOH]);

                // 计算CO2系统物种浓度
                simContext.CO2aq = KgwCO2 * Ppsia * simContext.yCO2 * gGas[iCO2aq] / (gNeut[iCO2aq] * gNNeut[iCO2aq]);
                simContext.HCO3 = (K1H2CO3 * aH2O) * simContext.CO2aq * gNeut[iCO2aq] * gNNeut[iCO2aq] /
                    (aH * gAn[iHCO3] * gNAn[iHCO3]);
                simContext.CO3 = K2HCO3 * simContext.HCO3 * gAn[iHCO3] * gNAn[iHCO3] /
                    (aH * gAn[iCO3] * gNAn[iCO3]);

                // 计算H2S系统物种浓度
                hydHS = aH * gAn[iHS] * gNAn[iHS] / (K1H2S * gNeut[iH2Saq] * gNNeut[iH2Saq]) + 1.0;
                if (TH2Saq > 0) {
                    simContext.HS = TH2Saq / hydHS;
                    //fMeSSpeciation(im, igas, simContext.HS, TFe, TZn, TPb, TH2Saq, hydHS, *ppt, & root1, & root2, & root3);
                    fMeSSpeciation(2, 2); // speciation for simContext.HS
                }
                else {
                    simContext.HS = 0.0;
                }

                // 计算乙酸系统物种浓度
                hydAc = aH * gAn[iAc] * gNAn[iAc] / (KHAc * gNeut[iHAcaq] * gNNeut[iHAcaq]) + 1.0; // Note gNNeut(iHACaq)=1
                AC = TAc / hydAc;

                // 计算硼酸系统物种浓度
                hydH2BO3 = aH * gAn[iH2BO3] * gNAn[iH2BO3] / (KH3BO3 * gNeut[iH3BO3] * gNNeut[iH3BO3]) + 1.0;
                H2BO3 = TH3BO3 / hydH2BO3;

                // 计算氨系统物种浓度
                hydNH3 = aH * gNeut[iNH3] * gNNeut[iNH3] / (KNH4 * gCat[iNH4] * gNCat[iNH4]) + 1.0;
                NH3 = TNH4 / hydNH3;

                // 计算硅酸系统物种浓度
                hydH2SiO4 = pow(aH, 2) * gAn[iH2SiO4] * gNAn[iH2SiO4] /
                    (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]) +
                    aH * gAn[iH2SiO4] * gNAn[iH2SiO4] /
                    (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]) + 1.0;
                H2SiO4 = TH4SiO4 / hydH2SiO4;
                H3SiO4 = H2SiO4 * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] /
                    (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]);
                H4SiO4 = H2SiO4 * pow(aH, 2) * gAn[iH2SiO4] * gNAn[iH2SiO4] /
                    (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]);

                // 计算碱度平衡误差
                faH = Alk - (simContext.HCO3 + 2 * simContext.CO3 + simContext.HS + AC + NH3 + H2BO3 + H3SiO4 + 2 * H2SiO4 + OH - H);

                // 调整pH搜索区间
                if (faH > 0) {
                    pHLow = pH;
                }
                if (faH < 0) {
                    pHHigh = pH;
                }
            }

            // 计算pH计读数
            pHMeterReading_from_QC = pH - DpHj;

            if (simContext.yH2S > 1) simContext.errmsg[2] = 3; // 注意：yH2S在此代码段中未计算

            // 根据条件设置最终pH计读数
            if (simContext.Run_Seawater_Mixing == 0 || simContext.Run_MixingTwoWells == 0) {
                pHMeterReading_from_QC = pH - DpHj;
            }
        }


        if (UseH2Sgas == 1) {
            // Use PCO2 and Alk to calculate the pH value at STP
            // 如果使用 PH2S，简化输入将跳过此部分

            // 读取yCO2和yH2S值
            //待定：yCO2 = Worksheets(mySheet).Cells(31, j + 2).Value / 100
            //待定：yH2S = Worksheets(mySheet).Cells(33, j + 2).Value / 100

            // 二分法求解pH值
            pHHigh = 14.0;
            pHLow = 0.0;
            //double faH;
            //double hydAc, hydH2BO3, hydNH3, hydH2SiO4;
            int k;

            // 这使得 pH 值能够收敛到 8 位有效数字；
            // 如果沉淀物仅占碱度的一小部分，有时就需要这样做。
            for (k = 0; k < 30; k++) {
                pH = (pHHigh + pHLow) / 2.0;
                aH = pow(10.0, -pH);

                // 计算H+和OH-浓度
                H = aH / (gCat[iH] * gNCat[iH]);
                OH = KH2O / (aH * gAn[iOH] * gNAn[iOH]);

                // 计算CO2系统物种浓度
                simContext.CO2aq = KgwCO2 * Ppsia * simContext.yCO2 * gGas[iCO2g] / (gNeut[iCO2aq] * gNNeut[iCO2aq]);
                simContext.HCO3 = (K1H2CO3 * aH2O) * simContext.CO2aq * gNeut[iCO2aq] * gNNeut[iCO2aq] /
                    (aH * gAn[iHCO3] * gNAn[iHCO3]);
                simContext.CO3 = K2HCO3 * simContext.HCO3 * gAn[iHCO3] * gNAn[iHCO3] /
                    (aH * gAn[iCO3] * gNAn[iCO3]);

                // 计算H2S系统物种浓度
                simContext.H2Saq = KgwH2S * Ppsia * simContext.yH2S * gGas[iH2Sg] / gNeut[iH2Saq] / gNNeut[iH2Saq];
                simContext.HS = K1H2S * simContext.H2Saq * gNeut[iH2Saq] * gNNeut[iH2Saq] /
                    (aH * gAn[iHS] * gNAn[iHS]);

                // 计算乙酸系统物种浓度
                hydAc = aH * gAn[iAc] * gNAn[iAc] / (KHAc * gNeut[iHAcaq] * gNNeut[iHAcaq]) + 1.0;
                AC = TAc / hydAc;

                // 计算硼酸系统物种浓度
                hydH2BO3 = aH * gAn[iH2BO3] * gNAn[iH2BO3] / (KH3BO3 * gNeut[iH3BO3] * gNNeut[iH3BO3]) + 1.0;
                H2BO3 = TH3BO3 / hydH2BO3;

                // 计算氨系统物种浓度
                hydNH3 = aH * gNeut[iNH3] * gNNeut[iNH3] / (KNH4 * gCat[iNH4] * gNCat[iNH4]) + 1.0;
                NH3 = TNH4 / hydNH3;

                // 计算硅酸系统物种浓度
                hydH2SiO4 = pow(aH, 2) * gAn[iH2SiO4] * gNAn[iH2SiO4] /
                    (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]) +
                    aH * gAn[iH2SiO4] * gNAn[iH2SiO4] /
                    (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]) + 1.0;
                H2SiO4 = TH4SiO4 / hydH2SiO4;
                H3SiO4 = H2SiO4 * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] /
                    (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]);
                H4SiO4 = H2SiO4 * pow(aH, 2) * gAn[iH2SiO4] * gNAn[iH2SiO4] /
                    (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]);

                // 计算碱度平衡误差
                faH = Alk - (simContext.HCO3 + 2 * simContext.CO3 + simContext.HS + AC + NH3 + H2BO3 + H3SiO4 + 2 * H2SiO4 + OH - H);

                // 调整pH搜索区间
                if (faH > 0) {
                    pHLow = pH;
                }
                if (faH < 0) {
                    pHHigh = pH;
                }
            }

            // 计算pH计读数
            pHMeterReading_from_QC = pH - DpHj;
        }

        // 输出结果到相应的工作表
        if (simContext.RunStat == 1 && simContext.CaseCount[0] == 1) {
            // for StatQC produced water
            //Worksheets("Output data sheet").Cells(5, 4) = pHMeterReading_from_QC  'QC check given Alk and PCO2, the calculated pH for meter reading.
        }
        else if (simContext.RunStat == 1 && simContext.CaseCount[0] == 9) {
            // for StatQC fresh water
            //Worksheets("Output data sheet").Cells(17, 4) = pHMeterReading_from_QC  'QC check given Alk and PCO2, the calculated pH for meter reading.
        }
        else if (simContext.RunH2SGUI == 1) {
            //Worksheets(mySheet).Cells(85, j + 2) = pHMeterReading_from_QC  'QC check given Alk and PCO2, the calculated pH for meter reading.
        }
        else {
            if (kk < simContext.nob_Input) {
                //Worksheets("Input").Cells(16, j + 15) = pHMeterReading_from_QC  'QC check given Alk and PCO2, the calculated pH for meter reading.
            }
            if (kk >= simContext.nob_Input && kk < simContext.nob_Input + simContext.nob_InputII) {
                //Worksheets("Input II").Cells(85, j + 2) = pHMeterReading_from_QC  'QC check given Alk and PCO2, the calculated pH for meter reading.
            }
        }
    }
    else {
        //对应于 useEOS<>0 QC 中的大多数参数已在 ReadInDataD sub 的末尾计算
        if (simContext.RunStat == 1 && simContext.CaseCount[0] == 1) {
            // for StatQC produced water
            //Worksheets("Output data sheet").Cells(6, 4) = compositions(2, 2) * 100
            //    Worksheets("Output data sheet").Cells(5, 4) = pHMeterReading
        }
        else if (simContext.RunStat == 1 && simContext.CaseCount[0] == 9) {
            // for StatQC fresh water
            //Worksheets("Output data sheet").Cells(18, 4) = compositions(2, 2) * 100
                //Worksheets("Output data sheet").Cells(17, 4) = pHMeterReading
        }
        else if (simContext.RunH2SGUI == 1) {
            //Worksheets(mySheet).Cells(83, j + 2) = compositions(3, 2) * 100

            if (simContext.UseMolal == 0) {
                //Worksheets(mySheet).Cells(84, j + 2) = (simContext.H2Saq + simContext.HS + S) * (34080 * (simContext.rho25c - TDS / 1000000#))
            }
            else {
                //Worksheets(mySheet).Cells(84, j + 2) = (simContext.H2Saq + simContext.HS + S)
            }

            //Worksheets(mySheet).Cells(85, j + 2) = pHMeterReading
            //    Worksheets(mySheet).Cells(86, j + 2) = compositions(2, 2) * 100 'QC check for assuming that pH and Alk are correct, what would PCO2 be?
        }
        else {
            if (kk < simContext.nob_Input) {
                //Worksheets("Input").Cells(14, j + 15) = compositions(3, 2) * 100

                if (simContext.UseMolal == 0) {
                    //Worksheets("Input").Cells(15, j + 15) = (simContext.H2Saq + simContext.HS + S) * (34080 * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#))
                }
                else {
                    //Worksheets("Input").Cells(15, j + 15) = TH2Saq
                }
                //Worksheets("Input").Cells(16, j + 15) = pHMeterReading
                //    Worksheets("Input").Cells(17, j + 15) = compositions(2, 2) * 100  'QC check for assuming that pH and Alk are correct, what would PCO2 be?
            }

            if (kk >= simContext.nob_Input && kk < simContext.nob_Input + simContext.nob_InputII) {
                //    Worksheets("Input II").Cells(83, j + 2) = compositions(3, 2) * 100

                if (simContext.UseMolal == 0) {
                    //Worksheets("Input II").Cells(84, j + 2) = (simContext.H2Saq + simContext.HS + S) * (34080 * (rho25CMix(i + nob_Input) - CalculatedTDSMix(i + nob_Input) / 1000000#))
                }
                else {
                    //Worksheets("Input II").Cells(84, j + 2) = (simContext.H2Saq + simContext.HS + S)
                }

                //Worksheets("Input II").Cells(85, j + 2) = pHMeterReading
                //    Worksheets("Input II").Cells(86, j + 2) = compositions(2, 2) * 100 'QC check for assuming that pH and Alk are correct, what would PCO2 be?
            }
        }
    }

    //3、根据报告的pH值和STP P-CO2计算碱度。注意，H2Saq值已在上文计算过。
    pHMeterReading = simContext.pHMeterStpMix[kk];
    simContext.yCO2 = simContext.yCO2Mix[kk];

    // 计算pH和相关参数
    pH = simContext.pHMeterStpMix[kk] + DpHj;
    aH = pow(10.0, -pH);

    H = aH / (gCat[iH] * gNCat[iH]);
    OH = KH2O / (aH * gAn[iOH] * gNAn[iOH]);

    // 根据useEOS选择计算CO2aq和H2Saq的方法
    if (useEOS == 0) {
        simContext.CO2aq = KgwCO2 * Ppsia * simContext.yCO2 * gGas[iCO2g] / (gNeut[iCO2aq] * gNNeut[iCO2aq]);
        // 注意：在useEOS=0时，H2Saq可能需要其他计算方法
    }
    else {
        // Unit Kg/ per mole of phase 3
        // 注意：VB中compositions(2,4)对应C中compositions[1][3]
        // compositions(15,4)对应compositions[14][3]
        if (compositions[14][3] != 0) { // 避免除以零
            simContext.CO2aq = compositions[1][3] / (0.01801528 * compositions[14][3]);
            simContext.H2Saq = compositions[2][3] / (0.01801528 * compositions[14][3]);
        }
        else {
            simContext.CO2aq = 0;
            simContext.H2Saq = 0;
        }
    }

    // 计算碳酸系统物种
    simContext.HCO3 = (K1H2CO3 * aH2O) * simContext.CO2aq * gNeut[iCO2aq] * gNNeut[iCO2aq] /
        (aH * gAn[iHCO3] * gNAn[iHCO3]);
    simContext.CO3 = K2HCO3 * simContext.HCO3 * gAn[iHCO3] * gNAn[iHCO3] /
        (aH * gAn[iCO3] * gNAn[iCO3]);

    // 计算硫化物系统物种
    simContext.HS = K1H2S * simContext.H2Saq * gNeut[iH2Saq] * gNNeut[iH2Saq] /
        (aH * gAn[iHS] * gNAn[iHS]);

    // 计算乙酸系统物种
    hydAc = aH * gAn[iAc] * gNAn[iAc] / (KHAc * gNeut[iHAcaq] * gNNeut[iHAcaq]) + 1;
    AC = TAc / hydAc;

    // 计算硼酸系统物种
    hydH2BO3 = aH * gAn[iH2BO3] * gNAn[iH2BO3] / (KH3BO3 * gNeut[iH3BO3] * gNNeut[iH3BO3]) + 1;
    H2BO3 = TH3BO3 / hydH2BO3;

    // 计算氨系统物种
    hydNH3 = aH * gNeut[iNH3] * gNNeut[iNH3] / (KNH4 * gCat[iNH4] * gNCat[iNH4]) + 1;
    NH3 = TNH4 / hydNH3;

    // 计算硅酸系统物种
    hydH2SiO4 = pow(aH, 2) * gAn[iH2SiO4] * gNAn[iH2SiO4] /
        (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]) +
        (aH)*gAn[iH2SiO4] * gNAn[iH2SiO4] /
        (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]) + 1;
    H2SiO4 = TH4SiO4 / hydH2SiO4;
    H3SiO4 = H2SiO4 * aH * gAn[iH2SiO4] * gNAn[iH2SiO4] /
        (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]);
    H4SiO4 = H2SiO4 * pow(aH, 2) * gAn[iH2SiO4] * gNAn[iH2SiO4] /
        (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]);

    // 计算碱度和钠平衡
    Alk_from_QC = (simContext.HCO3 + 2 * simContext.CO3 + simContext.HS + AC + NH3 + H2BO3 + H3SiO4 + 2 * H2SiO4 + OH - H);
    NaQC = (-SumOfAnions - (SumOfCations - simContext.NaMix[kk]));

    // 输出结果到相应的工作表
    if (simContext.RunStat == 1 && simContext.CaseCount[0] == 1) {
        /*
        *     Worksheets("Output Data Sheet").Cells(7, 4) = Alk_from_QC * (61019 * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#))
            Worksheets("Output Data Sheet").Cells(8, 4).Value = SumOfCations * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#) 'Convert output from equiv/kg-water to equiv/l-solution
            Worksheets("Output Data Sheet").Cells(9, 4).Value = SumOfAnions * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#) 'Convert output from equiv/kg-water to equiv/l-solution
            Worksheets("Output Data Sheet").Cells(10, 4).Value = rho25CMix(kk)
            Worksheets("Output Data Sheet").Cells(11, 4).Value = CalculatedTDSMix(kk)
            Worksheets("Output Data Sheet").Cells(12, 4).Value = NaQC * (22990 * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#))
        */
    }
    else if (simContext.RunStat == 1 && simContext.CaseCount[0] == 9) {
        /*
        Worksheets("Output Data Sheet").Cells(19, 4) = Alk_from_QC * (61019 * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#))
        Worksheets("Output Data Sheet").Cells(20, 4).Value = SumOfCations * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#) 'Convert output from equiv/kg-water to equiv/l-solution
        Worksheets("Output Data Sheet").Cells(21, 4).Value = SumOfAnions * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#) 'Convert output from equiv/kg-water to equiv/l-solution
        Worksheets("Output Data Sheet").Cells(22, 4).Value = rho25CMix(kk)
        Worksheets("Output Data Sheet").Cells(23, 4).Value = CalculatedTDSMix(kk)
        Worksheets("Output Data Sheet").Cells(24, 4).Value = NaQC * (22990 * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#))
        */
    }
    else if (simContext.RunH2SGUI == 1) {
        if (simContext.UseMolal == 0) {
            /*
            Worksheets(mySheet).Cells(87, j + 2) = Alk_from_QC * (61019 * (simContext.rho25c - TDS / 1000000#))
            Worksheets(mySheet).Cells(88, j + 2).Value = SumOfCations * (simContext.rho25c - TDS / 1000000#) 'Convert output from equiv/kg-water to equiv/l-solution
            Worksheets(mySheet).Cells(89, j + 2).Value = SumOfAnions * (simContext.rho25c - TDS / 1000000#) 'Convert output from equiv/kg-water to equiv/l-solution
            Worksheets(mySheet).Cells(91, j + 2).Value = NaQC * (22990 * (simContext.rho25c) - TDS / 1000000#)
            */
        }
        else {
            /*
            Worksheets(mySheet).Cells(87, j + 2) = Alk_from_QC
            Worksheets(mySheet).Cells(88, j + 2).Value = SumOfCations
            Worksheets(mySheet).Cells(89, j + 2).Value = SumOfAnions
            Worksheets(mySheet).Cells(91, j + 2).Value = NaQC
            */
        }
        //Worksheets(mySheet).Cells(90, j + 2).Value = TDS
    }
    else {
        if (kk < simContext.nob_Input) {
            if (simContext.UseMolal == 0) {
                /*
                Worksheets("Input").Cells(18, j + 15) = Alk_from_QC * (61019 * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#))
                Worksheets("Input").Cells(19, j + 15).Value = SumOfCations * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#) 'Convert output from equiv/kg-water to equiv/l-solution
                Worksheets("Input").Cells(20, j + 15).Value = SumOfAnions * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#) 'Convert output from equiv/kg-water to equiv/l-solution
                Worksheets("Input").Cells(22, j + 15).Value = NaQC * (22990 * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#))

                */
            }
            else {
                /*
                Worksheets("Input").Cells(18, j + 15) = Alk_from_QC
                Worksheets("Input").Cells(19, j + 15).Value = SumOfCations
                Worksheets("Input").Cells(20, j + 15).Value = SumOfAnions
                Worksheets("Input").Cells(22, j + 15).Value = NaQC
                */
            }
            //Worksheets("Input").Cells(21, j + 15).Value = CalculatedTDSMix(kk)
        }

        if (kk >= simContext.nob_Input && kk < simContext.nob_Input + simContext.nob_InputII) {
            if (simContext.UseMolal == 0) {
                /*
                Worksheets("Input II").Cells(87, j + 2) = Alk_from_QC * (61019 * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#))
                Worksheets("Input II").Cells(88, j + 2).Value = SumOfCations * (rho25CMix(i + nob_Input) - CalculatedTDSMix(i + nob_Input) / 1000000#) 'Convert output from equiv/kg-water to equiv/l-solution
                Worksheets("Input II").Cells(89, j + 2).Value = SumOfAnions * (rho25CMix(i + nob_Input) - CalculatedTDSMix(i + nob_Input) / 1000000#) 'Convert output from equiv/kg-water to equiv/l-solution
                Worksheets("Input II").Cells(91, j + 2).Value = NaQC * (22990 * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#))
                */
            }
            else {
                /*
                Worksheets("Input II").Cells(87, j + 2) = Alk_from_QC
                Worksheets("Input II").Cells(88, j + 2).Value = SumOfCations
                Worksheets("Input II").Cells(89, j + 2).Value = SumOfAnions
                Worksheets("Input II").Cells(91, j + 2).Value = NaQC

                */
            }
            //Worksheets("Input II").Cells(90, j + 2).Value = CalculatedTDSMix(i + nob_Input)
        }
    }

}

// 新增函数 by 黄志缘
/********************************************************************************************************/


/**
 * @brief 计算密度和pH（D2_CalcDensitypH）
 *
 * 此函数计算溶液的密度和pH值，通过迭代方法调整浓度以实现密度收敛。
 * 包括离子强度计算、pH调整、密度迭代、浓度修正，以及热力学平衡常数和活度系数的计算。
 * 适用于电解质溶液的密度和pH计算，处理CO2、H2S、FeSaq等物种的物种分布。
 *
 */
void D2_CalcDensitypH(int i, int Iteration, double* mt, int use_pH) {
    // Call CalcIonicStrength 'before CO2, H2S, FeSaq speciation
    CalcIonicStrength();

    pH = simContext.pHMeterStpMix[i] + DpHj;

    if (ISt >= 25) {
        printf("The calculated ionic strength is %.2f. This is greater than 20 m (moles of salt/kg of water), the upper limit. The calculation will be terminated. It is suggested that you check the input, conc unit, and retry the calculation.\n", ISt);
        return;  // End
    }
    Iteration = 0, rhoSSE = 1;
    while (rhoSSE > 0.00000001 && Iteration < 30) {
        *mt = fTPFunc(0);  // iTP=0 T=77F, P=14.696 psi: iTP=1 T=TVol, P=Pvol;iTP=2 T=TpH, P=PpH

        simContext.rho25c = CalcRhoTP(TK, TC, PBar, Patm);
        // *simContext.rho25c = CalcRhoTP(TK, TC, PBar, Patm);

        int iden;
        for (iden = 0; iden < NumCat; iden++) {
            mc[iden] *= (rhoOld - TDS / 1000000.0) / (simContext.rho25c - TDS / 1000000.0);
        }
        for (iden = 0; iden < NumAn; iden++) {
            ma[iden] *= (rhoOld - TDS / 1000000.0) / (simContext.rho25c - TDS / 1000000.0);
        }
        for (iden = 0; iden < NumNeut; iden++) {
            mn[iden] *= (rhoOld - TDS / 1000000.0) / (simContext.rho25c - TDS / 1000000.0);
        }
        Alk *= (rhoOld - TDS / 1000000.0) / (simContext.rho25c - TDS / 1000000.0);
        TAc *= (rhoOld - TDS / 1000000.0) / (simContext.rho25c - TDS / 1000000.0);
        TCO2 *= (rhoOld - TDS / 1000000.0) / (simContext.rho25c - TDS / 1000000.0);
        TNH4 *= (rhoOld - TDS / 1000000.0) / (simContext.rho25c - TDS / 1000000.0);
        TH3BO3 *= (rhoOld - TDS / 1000000.0) / (simContext.rho25c - TDS / 1000000.0);
        TH2Saq *= (rhoOld - TDS / 1000000.0) / (simContext.rho25c - TDS / 1000000.0);
        TH4SiO4 *= (rhoOld - TDS / 1000000.0) / (simContext.rho25c - TDS / 1000000.0);
        TFe *= (rhoOld - TDS / 1000000.0) / (simContext.rho25c - TDS / 1000000.0);

        rhoSSE = pow(simContext.rho25c - rhoOld, 2);
        rhoOld = simContext.rho25c;
        (Iteration)++;
    }

    xMeOH = 0;
    xMEG = 0;
    IStCosolvent = ISt;

    // Call CalcIonicStrength 'before CO2, H2S, FeSaq speciation
    CalcIonicStrength();

    pH = simContext.pHMeterStpMix[i] + DpHj;

    // If use_pH = 0 Then mt = fTPFunc(0) 'Option0 77F, 14.696 psi: Option1 T=TVol, P=Pvol; Option2 T=TpH, P=PpH
    // If use_pH = 1 Then mt = fTPFunc(2) 'Option0 77F, 14.696 psi: Option1 T=TVol, P=Pvol; Option2 T=TpH, P=PpH
    // If use_pH = 2 Or use_pH = 3 Then mt = fTPFunc(0) 'Option0 77F, 14.696 psi: Option1 T=TVol, P=Pvol; Option2 T=TpH, P=PpH

    // Call C1_ThermodynamicEquilConsts
    C1_ThermodynamicEquilConsts();

    // Call C2_PitzerActCoefs_T_P_ISt(gNeut, aH2O, TK, TC, PBar, Patm)
    C2_PitzerActCoefs_T_P_ISt(gNeut, aH2O, TK, TC, PBar, Patm);

    // Call PengRobinson3
    PengRobinson3();

    // Call C5_CalcpHPCO2PH2SSTP 'CO2, H2S, FeSaq speciation
    C5_CalcpHPCO2PH2SSTP(use_pH, UseH2Sgas, useEOS);

    mc[iH] = H;
    ma[iOH] = OH;
    ma[iAc] = AC;
    mn[iNH3] = NH3;
    ma[iH2BO3] = H2BO3;
    ma[iHCO3] = simContext.HCO3;
    ma[iCO3] = simContext.CO3;
    ma[iHS] = simContext.HS;
    ma[iH3SiO4] = H3SiO4;
    ma[iH2SiO4] = H2SiO4;
    mn[iH4SiO4aq] = H4SiO4;
    mn[iNH3] = NH3;
    mn[iH3BO3] = TH3BO3 - H2BO3;
    mn[iCO2aq] = simContext.CO2aq;
    mn[iH2Saq] = simContext.H2Saq;
    mn[iHAcaq] = HAcaq;

    if (simContext.useEOSmix[simContext.kk] == 1) {
        mc[iH] = 0.0000001;
        ma[iOH] = 0.0000001;
        ma[iAc] = TAc;
        ma[iHCO3] = Alk;
        ma[iHS] = TH2Saq;
        mn[iH4SiO4aq] = TH4SiO4;
        mn[iH3BO3] = TH3BO3;
        mc[iNH4] = TNH4;
        mn[iNH3] = 0;
        ma[iH2BO3] = 0;
        ma[iCO3] = 0;
        ma[iH3SiO4] = 0;
        ma[iH2SiO4] = 0;
        mn[iH3BO3] = 0;
        mn[iCO2aq] = 0;
        mn[iH2Saq] = 0;
        mn[iHAcaq] = 0;
    }
}


/**
 * @brief 计算密度（D1_CalcDensity）
 *
 * 此函数计算混合物i的密度，通过迭代调整TDS和浓度以实现收敛。
 * 支持摩尔浓度（UseMolal=0）和直接计算（UseMolal=1）两种模式。
 * 在摩尔浓度模式下，进行TDS迭代；在直接模式下，计算离子强度、活度系数和物种分布。
 *
 */
void D1_CalcDensity(int i, int* Iteration2, double* mt) {

    // HCO3stpMix(i) = AlkMix(i): ACstpMix(i) = TAcMix(i): HstpMix(i) = 0.000001: OHstpMix(i) = 0.0000001: CO3stpMix(i) = 0: _
    // HSstpMix(i) = 0: NH4STPMix(i) = TNH4Mix(i): H2BO3stpMix(i) = 0
    simContext.HCO3stpMix[i] = simContext.AlkMix[i];
    simContext.ACstpMix[i] = simContext.TAcMix[i];
    simContext.HstpMix[i] = 0.000001;
    simContext.OHstpMix[i] = 0.0000001;
    simContext.CO3stpMix[i] = 0;
    simContext.HSstpMix[i] = 0;
    simContext.NH4STPMix[i] = simContext.TNH4Mix[i];
    simContext.H2BO3stpMix[i] = 0;

    *Iteration2 = 0;

    mc[iH] = simContext.HstpMix[i];
    mc[iNa] = simContext.NaMix[i];  // 注意：NaMix等需从全局或参数中获取，假设已定义
    mc[iK] = simContext.KMix[i];
    mc[iMg] = simContext.MgMix[i];
    mc[iCa] = simContext.CaMix[i];
    TCa = mc[iCa];
    mc[iSr] = simContext.SrMix[i];
    mc[iBa] = simContext.BaMix[i];
    mc[iFe] = simContext.FeMix[i];
    mc[iZn] = simContext.ZnMix[i];
    mc[iPb] = simContext.PbMix[i];

    ma[iOH] = simContext.OHstpMix[i];
    ma[iCl] = simContext.ClMix[i];
    ma[iAc] = simContext.ACstpMix[i];
    mc[iNH4] = simContext.NH4STPMix[i];
    ma[iH2BO3] = simContext.H2BO3stpMix[i];
    ma[iHCO3] = simContext.HCO3stpMix[i];
    ma[iCO3] = simContext.CO3stpMix[i];

    ma[iH3SiO4] = 0;
    ma[iH2SiO4] = 0;

    ma[iSO4] = simContext.SO4Mix[i];
    ma[iHS] = simContext.HSstpMix[i];
    ma[intF] = simContext.FMix[i];
    ma[iBr] = simContext.BrMix[i];

    Alk = simContext.AlkMix[i];
    TAc = simContext.TAcMix[i];
    TH2Saq = simContext.TH2SaqMix[i];
    TH4SiO4 = simContext.TH4SiO4Mix[i];
    TH3BO3 = simContext.TH3BO3Mix[i];
    TNH4 = simContext.TNH4Mix[i];

    mn[iNH3] = 0;
    mn[iH3BO3] = TH3BO3;
    mn[iH4SiO4aq] = TH4SiO4;
    TFe = mc[iFe];

    // If use_pH = 3 Then TCO2 = TCO2Mix(i)
    // TDSOld = TDSMix(i):  rhoOld = rho_Mix(i): TDS = TDSMix(i): TDSSSE = 1:
    TDSOld = simContext.TDSMix[i];
    rhoOld = simContext.rho_Mix[i];  // 假设rhoOld在D2_CalcDensitypH中处理
    TDS = simContext.TDSMix[i];
    TDSSSE = 1;

    // simContext.yCO2 = yCO2Mix(i): simContext.yH2S = yH2SMix(i): yCH4 = 1 - simContext.yCO2 - simContext.yH2S
    //simContext.yCO2 = yCO2Mix[i];
    //simContext.yH2S = yH2SMix[i];
    //yCH4 = 1 - simContext.yCO2 - simContext.yH2S;

    if (simContext.UseMolal == 0) {
        while (TDSSSE > 0.00000001 && *Iteration2 < 20) {
            // Call D2_CalcDensitypH(i) 'Calculate ISt, density, and simContext.HCO3, AC, simContext.HS speciation from TDS
            D2_CalcDensitypH(i, Iteration, mt, use_pH);
            TDS = 0;
            CalculateTDSDen = 0;  // Calculate TDS from density
            int iden;
            for (iden = 1; iden < NumCat; iden++) {
                TDS += 1000 * (simContext.rho25c) * mc[iden] * MWCat[iden];  // =Sum of mg salt/L*(Kg soln/Kg H2O)
                CalculateTDSDen += 0.001 * mc[iden] * MWCat[iden];  // =Sum of Kg salt/Kg H2O
            }
            for (iden = 1; iden < NumAn; iden++) {
                TDS += 1000 * (simContext.rho25c) * ma[iden] * MWAn[iden];
                CalculateTDSDen += 0.001 * ma[iden] * MWAn[iden];
            }
            for (iden = 1; iden < NumNeut; iden++) {
                TDS += 1000 * (simContext.rho25c) * mn[iden] * MWNeut[iden];
                CalculateTDSDen += 0.001 * mn[iden] * MWNeut[iden];
            }
            TDS /= (1 + CalculateTDSDen);  // denometer=(1+Kgsalt/KgH2O)=(Kgsoln/KgH2O)

            for (iden = 1; iden < NumCat; iden++) {  // Calculate molality from new TDS
                mc[iden] *= ((simContext.rho25c) - (TDSOld) / 1000000.0) / ((simContext.rho25c) - (TDS) / 1000000.0);
            }
            for (iden = 1; iden < NumAn; iden++) {
                ma[iden] *= ((simContext.rho25c) - (TDSOld) / 1000000.0) / ((simContext.rho25c) - (TDS) / 1000000.0);
            }
            for (iden = 1; iden < NumNeut; iden++) {
                mn[iden] *= ((simContext.rho25c) - (TDSOld) / 1000000.0) / ((simContext.rho25c) - (TDS) / 1000000.0);
            }
            Alk *= ((simContext.rho25c) - (TDSOld) / 1000000.0) / ((simContext.rho25c) - (TDS) / 1000000.0);
            TAc *= ((simContext.rho25c) - (TDSOld) / 1000000.0) / ((simContext.rho25c) - (TDS) / 1000000.0);
            TCO2 *= (simContext.rho25c - TDSOld / 1000000.0) / (simContext.rho25c - TDS / 1000000.0);
            TNH4 *= ((simContext.rho25c) - (TDSOld) / 1000000.0) / ((simContext.rho25c) - (TDS) / 1000000.0);
            TH3BO3 *= ((simContext.rho25c) - (TDSOld) / 1000000.0) / ((simContext.rho25c) - (TDS) / 1000000.0);
            TH2Saq *= ((simContext.rho25c) - (TDSOld) / 1000000.0) / ((simContext.rho25c) - (TDS) / 1000000.0);
            TH4SiO4 *= ((simContext.rho25c) - (TDSOld) / 1000000.0) / ((simContext.rho25c) - (TDS) / 1000000.0);
            TFe *= ((simContext.rho25c) - (TDSOld) / 1000000.0) / ((simContext.rho25c) - (TDS) / 1000000.0);

            // Call D2_CalcDensitypH(i) 'Calculate ISt, density, and simContext.HCO3, AC, simContext.HS speciation from TDS
            D2_CalcDensitypH(i, Iteration, mt, use_pH);

            if (TDSOld == 0) goto label10;
            TDSSSE = pow((TDS / TDSOld) - 1, 2);
            TDSOld = TDS;
            (*Iteration2)++;
        }
    }
    else {
        // 未定义函数：计算离子强度
        CalcIonicStrength();

        xMeOH = 0;
        xMEG = 0;
        IStCosolvent = ISt;

        *mt = fTPFunc(0);  // iTP=0 T=77F, P=14.696 psi: iTP=1 T=TVol, P=Pvol;iTP=2 T=TpH, P=PpH

        // simContext.rho25c = CalcRhoTP(TK, TC, PBar, Patm) 'Function subroutine  // 未定义函数：计算密度
        simContext.rho25c = CalcRhoTP(TK, TC, PBar, Patm);  // 参数：TK, TC, PBar, Patm -> double

        pH = simContext.pHMeterStpMix[i] + DpHj;

        // amy check ????????????????????????????????????????????????????
        // If use_pH = 0 Then mt = fTPFunc(0) 'Option0 77F, 14.696 psi: Option1 T=TVol, P=Pvol; Option2 T=TpH, P=PpH
        // If use_pH = 1 Then mt = fTPFunc(2) 'Option0 77F, 14.696 psi: Option1 T=TVol, P=Pvol; Option2 T=TpH, P=PpH
        // If use_pH = 2 Or use_pH = 3 Then mt = fTPFunc(0) 'Option0 77F, 14.696 psi: Option1 T=TVol, P=Pvol; Option2 T=TpH, P=PpH

        // Call C1_ThermodynamicEquilConsts  // 未定义函数：计算热力学平衡常数
        C1_ThermodynamicEquilConsts();  // 参数：无（使用全局TK, PBar等）

        // Call C2_PitzerActCoefs_T_P_ISt(gNeut, aH2O, TK, TC, PBar, Patm)  // 已部分定义
        C2_PitzerActCoefs_T_P_ISt(gNeut, aH2O, TK, TC, PBar, Patm);  // 参数：gNeut[], aH2O*, TK, TC, PBar, Patm

        PengRobinson3();

        C5_CalcpHPCO2PH2SSTP(use_pH, UseH2Sgas, useEOS);
        mc[iH] = H;
        ma[iOH] = OH;
        ma[iAc] = AC;
        mn[iNH3] = NH3;
        ma[iH2BO3] = H2BO3;
        ma[iHCO3] = simContext.HCO3;
        ma[iCO3] = simContext.CO3;
        ma[iHS] = simContext.HS;
        ma[iH3SiO4] = H3SiO4;
        ma[iH2SiO4] = H2SiO4;
        mn[iH4SiO4aq] = H4SiO4;
        mn[iNH3] = NH3;
        mn[iH3BO3] = TH3BO3 - H2BO3;
        mn[iCO2aq] = simContext.CO2aq;
        mn[iH2Saq] = simContext.H2Saq;
        mn[iHAcaq] = HAcaq;

        if (simContext.useEOSmix[simContext.kk] == 1) {
            mc[iH] = 0.0000001;
            ma[iOH] = 0.0000001;
            ma[iAc] = TAc;
            ma[iHCO3] = Alk;
            ma[iHS] = 0;
            mn[iH4SiO4aq] = TH4SiO4;
            mn[iH3BO3] = TH3BO3;
            mc[iNH4] = TNH4;
        }

        *mt = fTPFunc(0);  // iTP=0 T=77F, P=14.696 psi: iTP=1 T=TVol, P=Pvol;iTP=2 T=TpH, P=PpH
        CalcIonicStrength();

        //C5_CalcpHPCO2PH2SSTP(use_pH, UseH2Sgas, useEOS);

        simContext.rho25c = CalcRhoTP(TK, TC, PBar, Patm);

        // If simContext.yCO2 + simContext.yH2S <= 1 Then  ' UseTPpH is chosen, the gas composition is calculated at T,P of pH
        // yCH4 = 1 - (simContext.yCO2 + simContext.yH2S)
        // Else
        // yCH4 = 0
        // End If
    }

label10:;
}


/**
 * @brief 读取输入部分C（ReadInputPartC）
 *
 * 此函数读取并设置混合物kk的输入参数，包括标准条件下的浓度、pH选项、气体组成等。
 * 计算密度和TDS，更新各种浓度数组，并处理摩尔浓度到TDS的转换。
 * 适用于电解质溶液混合物的输入处理，支持多种运行模式（如测试案例、海水混合等）。
 *
 */
void ReadInputPartC(int kk, double* mt, int* Iteration2) {
    *mt = fTPFunc(0);  // Densitym TDS, and m calculated at STP condition
    if (simContext.UseTPpHMix[kk] == 1) *mt = fTPFunc(2);

    if (simContext.Run10TestCases == 1 && simContext.Loop10 > 1) goto label100;
    if (simContext.Run_Seawater_Mixing == 1 && simContext.LoopMixing > 1) goto label100;
    if (simContext.Run_MixingTwoWells == 1 && simContext.LoopMixing > 1) goto label100;
    if (simContext.RunMultiMix == 1 && simContext.LoopResChem > 1) goto label100;
    if (simContext.RunStatMix == 1 && simContext.LoopMixing > 1) goto label100;

    simContext.HCO3stpMix[kk] = simContext.AlkMix[kk];
    simContext.ACstpMix[kk] = simContext.TAcMix[kk];
    simContext.HstpMix[kk] = 0.000001;
    simContext.OHstpMix[kk] = 0.0000001;
    simContext.CO3stpMix[kk] = 0;
    simContext.HSstpMix[kk] = 0;
    simContext.NH4STPMix[kk] = simContext.TNH4Mix[kk];
    simContext.H2BO3stpMix[kk] = 0;
    TDS = 0;

    *Iteration2 = 0;

    mc[iH] = simContext.HstpMix[kk];
    mc[iNa] = simContext.NaMix[kk];
    mc[iK] = simContext.KMix[kk];
    mc[iMg] = simContext.MgMix[kk];
    mc[iCa] = simContext.CaMix[kk];
    TCa = mc[iCa];
    mc[iSr] = simContext.SrMix[kk];
    mc[iBa] = simContext.BaMix[kk];
    mc[iFe] = simContext.FeMix[kk];
    mc[iZn] = simContext.ZnMix[kk];
    mc[iPb] = simContext.PbMix[kk];
    mc[iRa] = simContext.RaMix[kk];

    ma[iOH] = simContext.OHstpMix[kk];
    ma[iCl] = simContext.ClMix[kk];
    ma[iAc] = simContext.ACstpMix[kk];
    mc[iNH4] = simContext.NH4STPMix[kk];
    ma[iH2BO3] = simContext.H2BO3stpMix[kk];
    ma[iHCO3] = simContext.HCO3stpMix[kk];
    ma[iCO3] = simContext.CO3stpMix[kk];

    ma[iH3SiO4] = 0;
    ma[iH2SiO4] = 0;

    ma[iSO4] = simContext.SO4Mix[kk];
    ma[iHS] = simContext.HSstpMix[kk];
    ma[intF] = simContext.FMix[kk];
    ma[iBr] = simContext.BrMix[kk];

    Alk = simContext.AlkMix[kk];
    TAc = simContext.TAcMix[kk];
    TH2Saq = simContext.TH2SaqMix[kk];
    TH4SiO4 = simContext.TH4SiO4Mix[kk];
    TH3BO3 = simContext.TH3BO3Mix[kk];
    TNH4 = simContext.TNH4Mix[kk];

    TFe = simContext.FeMix[kk];
    TPb = simContext.PbMix[kk];
    TZn = simContext.ZnMix[kk];

    mn[iNH3] = 0;
    mn[iH3BO3] = TH3BO3;
    mn[iH4SiO4aq] = TH4SiO4;

    use_pH = simContext.usepHmix[kk];
    UseH2Sgas = simContext.UseH2SgasMix[kk];

    if (use_pH == 3) {
        TCO2 = simContext.TCO2Mix[kk];
    }
    else {
        TCO2 = 0;
    }

    simContext.yCO2 = simContext.yCO2Mix[kk];
    simContext.yH2S = simContext.yH2SMix[kk];  // yCH4 = 1 - simContext.yCO2 - simContext.yH2S 'assume dry gas

    // If yCH4 < 0 Then yCH4 = 0  '?????????????????????????????????????????????????????????????

    if (simContext.useEOSmix[kk] == 1 && simContext.yCO2 == 0 && simContext.SumofZMix[kk] > 0) simContext.yCO2 = simContext.zMix[kk][1];  // if UseEOS=1 then set YCO2 and YH2S to reservoir condition to calculate density and TDS only if the reservoir fluid comp is given
    if (simContext.useEOSmix[kk] == 1 && simContext.yH2S == 0 && simContext.SumofZMix[kk] > 0) simContext.yH2S = simContext.zMix[kk][2];

    // Call C1_ThermodynamicEquilConsts  'Only function of T,P, does not recalculate in D1_CalcDensity
    // Call PengRobinson3

    // Call D1_CalcDensity(kk) 'Calculate TDS, density and fix molality based on the predicted density and TDS
    D1_CalcDensity(kk, Iteration2, mt);

    simContext.CalculatedTDSMix[kk] = TDS;

    simContext.HstpMix[kk] = mc[iH];
    simContext.NaMix[kk] = mc[iNa];
    simContext.KMix[kk] = mc[iK];
    simContext.MgMix[kk] = mc[iMg];
    simContext.CaMix[kk] = mc[iCa];
    TCa = mc[iCa];

    simContext.SrMix[kk] = mc[iSr];
    simContext.BaMix[kk] = mc[iBa];
    simContext.FeMix[kk] = TFe;
    simContext.ZnMix[kk] = mc[iZn];
    simContext.PbMix[kk] = mc[iPb];
    simContext.RaMix[kk] = mc[iRa];

    simContext.OHstpMix[kk] = ma[iOH];
    simContext.ClMix[kk] = ma[iCl];
    simContext.ACstpMix[kk] = ma[iAc];
    simContext.NH4STPMix[kk] = mc[iNH4];
    simContext.H2BO3stpMix[kk] = ma[iH2BO3];
    simContext.HCO3stpMix[kk] = ma[iHCO3];
    simContext.CO3stpMix[kk] = ma[iCO3];

    simContext.SO4Mix[kk] = ma[iSO4];
    simContext.HSstpMix[kk] = ma[iHS];
    simContext.FMix[kk] = ma[intF];
    simContext.BrMix[kk] = ma[iBr];

    simContext.rho25CMix[kk] = simContext.rho25c;

    simContext.H3SiO4Mix[kk] = ma[iH3SiO4];
    simContext.H2SiO4Mix[kk] = ma[iH2SiO4];

    simContext.NH3Mix[kk] = mn[iNH3];
    simContext.H4SiO4Mix[kk] = mn[iH4SiO4aq];
    simContext.H3BO3Mix[kk] = mn[iH3BO3];

    simContext.CO2aqMix[kk] = mn[iCO2aq];
    simContext.H2SaqMix[kk] = mn[iH2Saq];
    simContext.HACaqMix[kk] = mn[iHAcaq];

    simContext.AlkMix[kk] = Alk;
    simContext.TAcMix[kk] = TAc;
    simContext.TH2SaqMix[kk] = TH2Saq;
    simContext.TH4SiO4Mix[kk] = TH4SiO4;
    simContext.TNH4Mix[kk] = TNH4;
    simContext.TH3BO3Mix[kk] = TH3BO3;

    simContext.yCO2Mix[kk] = simContext.yCO2;
    simContext.yH2SMix[kk] = simContext.yH2S;
    simContext.yCH4Mix[kk] = 1 - simContext.yCO2 - simContext.yH2S;  // Set YCO2 and YH2S to the calculated value if pH option is used.

    if (simContext.UseTPVolMix[kk] == 0) simContext.WaterDensityMix[kk] = simContext.rho25c;

    if (simContext.UseMolal == 1) {
        TDS = 0;
        CalculateTDSDen = 0;  // Calculate TDS from density

        int iden;
        for (iden = 1; iden < NumCat; iden++) {
            CalculateTDSDen += 0.001 * mc[iden] * MWCat[iden];  // =Sum of g salt/Kg H2O
        }
        for (iden = 1; iden < NumAn; iden++) {
            CalculateTDSDen += 0.001 * ma[iden] * MWAn[iden];
        }
        for (iden = 1; iden < NumNeut; iden++) {
            CalculateTDSDen += 0.001 * mn[iden] * MWNeut[iden];
        }

        TDS = CalculateTDSDen / (1 + CalculateTDSDen) * simContext.rho25c * 1000000.0;  // TDS in unit of mg/L,  numerator=(Kgsalt/KgH2O), denometer=(1+Kgsalt/KgH2O)=(Kgsoln/KgH2O);density Kgsoln/Lsoln
        simContext.CalculatedTDSMix[kk] = TDS;
    }
label100:;
}
/********************************************************************************************************/

void ReadInputPartD(int kk, int j, SampleData* data)
{
    int mol_w3;
    double mol_W, mol_o, mol_g, mol_HC;

    DoubleVec feed_Composition(15, 0.0);

    double mw_oil;

    DoubleVec tempgNeut = { gNeut[iCO2aq], gNeut[iH2Saq] };

    //-----------------------------------------------------
    int iNG = 0;
    for (iNG = 0; iNG < 15; iNG++) {
        zOutput[iNG] = 0;
        z[iNG] = 0;
        simContext.zMix[kk][iNG] = 0;
    }
    for (iNG = 0; iNG < 2; iNG++)
        density[iNG] = 0;

    int usedryHC = 0;

    if (simContext.RunStat == 0) {
        // vb 中 RunMultiMix, RunWhatIf, Run_Seawater_Mixing 未赋初值
        if (simContext.Run10TestCases == 1 || Run_MassTransfer == 1 || simContext.Run_Seawater_Mixing == 1 || simContext.Run_MixingTwoWells == 1 || simContext.RunWhatIf == 1 || simContext.RunMultiMix == 1) {
            if (simContext.nob_Input + simContext.nob_InputII == 1) {
                //待定：Worksheets(mySheet).Cells(55, j + 2).Value = 0
                data->Option_EoS = 0;
            }
            else {
                //待定：Worksheets(MySheetMix).Cells(55, 8).Value = Empty
                data->Option_EoS = 0;
            }
        }
        if (simContext.nob_Input + simContext.nob_InputII == 1) {
            //待定：useEOSmix(kk) = Worksheets(mySheet).Cells(55, j + 2).Value
            // 待定值
            simContext.useEOSmix[kk] = data->Option_EoS;
        }
        else {
            //待定：useEOSmix(kk) = Worksheets(MySheetMix).Cells(55, 8).Value 'UseEOSMix  values for mixed solution IS set on Column H
            // 待定值
            simContext.useEOSmix[kk] = data->Option_EoS;
        }


        //usedryHC = Worksheets(mySheet).Cells(56, j + 2).Value
        usedryHC = data->Option_Water_HC;

        if (RunShellMultiflash == 1) simContext.useEOSmix[kk] = 0;

        /*
        * 待定，注：useEOSmix是整形的，为什么会用""来判断？？？？？
        if (useEOSmix[kk] == "") useEOSmix(kk) = 0
        */
        if (simContext.RunMultiMix == 1) simContext.useEOSmix[kk] = 0;

        //注意，由于 ZMix (kk,i) 在下面重新分配，因此值将从此 Sub 开头的输入表中重新读取。
        //for (iNG = 0; iNG < 14; iNG++) {
        //    //待定：zMix(kk, iNG) = Worksheets(mySheet).Cells(65 + iNG, j + 2) / 100:
        //}
        simContext.zMix[kk][0] = data->C1_o / 100.0;           simContext.zMix[kk][1] = data->CO2_o / 100.0;
        simContext.zMix[kk][2] = data->H2S_o / 100.0;          simContext.zMix[kk][3] = data->C2_o / 100.0;
        simContext.zMix[kk][4] = data->C3_o / 100.0;           simContext.zMix[kk][5] = data->iC4_o / 100.0;
        simContext.zMix[kk][6] = data->nC4_o / 100.0;          simContext.zMix[kk][7] = data->iC5_o / 100.0;
        simContext.zMix[kk][8] = data->nC5_o / 100.0;          simContext.zMix[kk][9] = data->C6_o / 100.0;
        simContext.zMix[kk][10] = data->C7_C12_o / 100.0;      simContext.zMix[kk][11] = data->C13_C25_o / 100.0;
        simContext.zMix[kk][12] = data->C26_C80_o / 100.0;     simContext.zMix[kk][13] = data->N2_o / 100.0;
        simContext.zMix[kk][14] = 0.0;
    }

    simContext.SumofZMix[kk] = 0.0;
    for (iNG = 0; iNG < 14; iNG++) {
        if (simContext.zMix[kk][iNG] < 1e-7) simContext.zMix[kk][iNG] = 0;
        simContext.SumofZMix[kk] += simContext.zMix[kk][iNG];
    }
    if (simContext.SumofZMix[kk] > 0) {
        for (iNG = 0; iNG < 14; iNG++)//Normalized z[0] to z[13]
            simContext.zMix[kk][iNG] /= simContext.SumofZMix[kk];
    }

    //为该子程序分配局部变量
    mc[iH] = simContext.HstpMix[kk]; mc[iNa] = simContext.NaMix[kk]; mc[iK] = simContext.KMix[kk]; mc[iMg] = simContext.MgMix[kk]; mc[iCa] = simContext.CaMix[kk];
    double TCa = mc[iCa];
    mc[iSr] = simContext.SrMix[kk]; mc[iBa] = simContext.BaMix[kk]; mc[iFe] = simContext.FeMix[kk]; mc[iZn] = simContext.ZnMix[kk]; mc[iPb] = simContext.PbMix[kk]; mc[iRa] = simContext.RaMix[kk];

    ma[iOH] = simContext.OHstpMix[kk]; ma[iCl] = simContext.ClMix[kk]; ma[iAc] = simContext.ACstpMix[kk]; mc[iNH4] = simContext.NH4STPMix[kk]; ma[iH2BO3] = simContext.H2BO3stpMix[kk];
    ma[iHCO3] = simContext.HCO3stpMix[kk]; ma[iCO3] = simContext.CO3stpMix[kk];
    ma[iSO4] = simContext.SO4Mix[kk]; ma[iHS] = simContext.HSstpMix[kk]; ma[intF] = simContext.FMix[kk]; ma[iBr] = simContext.BrMix[kk];

    Alk = simContext.AlkMix[kk]; TAc = simContext.TAcMix[kk]; TNH4 = simContext.TNH4Mix[kk]; TH3BO3 = simContext.TH3BO3Mix[kk]; TH2Saq = simContext.TH2SaqMix[kk]; TCO2 = simContext.TCO2Mix[kk];

    VW = simContext.VwMix[kk]; VgTP = simContext.VgTPMix[kk]; VO = simContext.VoMix[kk]; VMeOH = simContext.VMeOHMix[kk]; VMEG = simContext.VMEGMix[kk]; mass_MeOH = simContext.mass_MeOH_mix[kk]; mass_MEG = simContext.mass_MEG_mix[kk];

    simContext.yCO2 = simContext.yCO2Mix[kk], simContext.yH2S = simContext.yH2SMix[kk], yCH4 = simContext.yCH4Mix[kk];   // Local variable values; in this loop only.

    useEOS = simContext.useEOSmix[kk]; use_pH = simContext.usepHmix[kk]; UseH2Sgas = simContext.UseH2SgasMix[kk];

    TFe = mc[iFe];

    for (iNG = 0; iNG < 14; iNG++) {
        z[iNG] = simContext.zMix[kk][iNG];
        if (useEOS == 3 || useEOS == 0)
            z[iNG] = 0;
    }

    SGG = simContext.gasSpGravMix[kk]; data->API = simContext.oilAPIgravMix[kk];

    //根据Cragoe, C.S.1929《石油产品的热力学性质》（美国商务部标准局杂项出版物第97号）计算石油的平均分子量。
    if (data->API > 20 && data->API < 80) mw_oil = 6084.0 / (data->API - 5.9);
    if (data->API <= 20) mw_oil = 6084 / (20 - 5.9);
    if (data->API >= 80) mw_oil = 6084 / (80 - 5.9);

    //根据流体温度和压力计算水、油的质量以及水、油和气体的摩尔数

    if (simContext.UseTPVolMix[kk] == 1) fTPFunc(1);//检查 T、P 以进行 mass_W、mass_O、mol_g 计算
    if (simContext.UseTPVolMix[kk] == 0) fTPFunc(0);

    CalcIonicStrength();
    pH = simContext.pHMeterStpMix[kk] + DpHj;
    RatioOilBPoints = fRatioOilBPoints(data->API);
    C1_ThermodynamicEquilConsts();
    C2_PitzerActCoefs_T_P_ISt(gNeut, aH2O, TK, TC, PBar, Patm); //计算在标准温度和压力下水合物抑制剂存在下的作用系数
    PengRobinson3();
    //注意，如果 useEOSmix(kk)<>0，则省略此 pH 计算步骤。换句话说，如果此处 useEOS<>0，则不会运行 pH 和形态分析。形态分析已在 ReadInputPartC 中完成
    //C5_CalcpHPCO2PH2SSTP();
    C5_CalcpHPCO2PH2SSTP(use_pH, UseH2Sgas, useEOS);
    //重新分配 CO2aq、HCO3、CO3、H2Saq、HS 并重新计算离子强度
    //mt = fmn();
    fmn();

    CalcIonicStrength();
    C2_PitzerActCoefs_T_P_ISt(gNeut, aH2O, TK, TC, PBar, Patm);
    C5_CalcpHPCO2PH2SSTP(use_pH, UseH2Sgas, useEOS);

    //使用温度 TVol 和 PVol 下计算的水密度
    if (simContext.UseTPVolMix[kk] == 1) simContext.WaterDensityMix[kk] = CalcRhoTP(TK, TC, PBar, Patm);

    simContext.mass_o_Mix[kk] = 159 * simContext.VoMix[kk] * simContext.OilDensityMix[kk]; // 这些用于在 useEOS=0 和 useEOS=3 中计算 nTCO2 和 nTH2S。这些值在 useEOS=1 或 2 中重新计算。
    simContext.mass_w_Mix[kk] = 159 * simContext.VwMix[kk] * simContext.WaterDensityMix[kk]; // 换算成公斤盐水
    simContext.mass_w_Mix[kk] = simContext.mass_w_Mix[kk] * (1 - simContext.CalculatedTDSMix[kk] / simContext.rho25CMix[kk] * 0.000001);

    Mass_o = simContext.mass_o_Mix[kk];
    mass_w = simContext.mass_w_Mix[kk];

    //mt = fTotalCO2H2Smoles(); //计算每种气体成分的总摩尔数（无论是在气体中还是在油或水中）
    fTotalCO2H2Smoles();
    //请注意，VgTP 的单位是 m^3，当 Vg 的单位是 MMCF 时，829 是从 He 的系数 78740 转换而来的。

    simContext.yCO2Mix[kk] = simContext.yCO2;
    simContext.yCH4Mix[kk] = yCH4;
    simContext.yH2SMix[kk] = simContext.yH2S;

    simContext.nTCO2Mix[kk] = nTCO2;
    simContext.nTCH4Mix[kk] = nTCH4;
    simContext.nTH2SMix[kk] = nTH2S;// vb 536.894123326537 cpp 537.8953812802223

    double YH2O = PsatH2O(TK) / PBar;

    double mol_w_Orig = 1000 * mass_w / 18.01528; // moles of water per day

    if (usedryHC == 1)
        mol_W = mol_w_Orig;
    else
        mol_W = 1000 * mass_w / 18.01528 + PsatH2O(TK) * VgTP * 1000 / Znew / RBar / TK; //包括气相中的 H2O，仅用作使用 EOS > 0 时的初始猜测



    //***************步骤 1  碳氢化合物调节**************

    //If useEOS <> 0 Then
    if (useEOS != 0)
    {
        if (useEOS == 3)
        {
            mol_o = 1000 * Mass_o / mw_oil; // moles of oil per day
            // vb 9860366.98058253 计算有误 VgTP不同
            mol_g = VgTP * PBar * 1000 / (Znew * RBar * TK); // moles of gas per day
            mol_HC = mol_o + mol_g;

            if ((simContext.VoMix[kk] == (1.0 / 159.0 / 1000.0)) && (VgTP == 0.000001)) // 当气体和油都为零并混合水时
            {
                if (simContext.nob == 1)
                {
                    simContext.errmsg[13] = 14;
                    useEOS = 0;
                    goto label_3003;
                }
                else
                {
                    z_before_precipitation[1] = simContext.nTCO2Mix[kk] / mol_W;
                    z_before_precipitation[2] = simContext.nTH2SMix[kk] / mol_W;
                    z_before_precipitation[14] = 1 - z_before_precipitation[1] - z_before_precipitation[2];
                    simContext.Total_molesMix[kk] = (mol_W);
                    simContext.nTCO2MixEOS[kk] = nTCO2;
                    simContext.nTH2SMixEOS[kk] = nTH2S;
                    //mt = fTPFunc(0);
                    fTPFunc(0);

                    if (simContext.usepHmix[kk] == 1) pH = simContext.pHMeterStpMix[kk] + DpHj;

                    CalcIonicStrength();
                    RatioOilBPoints = fRatioOilBPoints(data->API);
                    C1_ThermodynamicEquilConsts();
                    C2_PitzerActCoefs_T_P_ISt(gNeut, aH2O, TK, TC, PBar, Patm); // 计算在标准温度和压力下水合物抑制剂存在下的作用系数
                    PengRobinson3();

                    // ***注意，如果 useEOSmix[kk]!=0，则省略此 pH 计算步骤。换句话说，如果此处 useEOS!=0，则不会运行 pH 和形态分析。形态分析已在 ReadInputPartC 中完成。
                    C5_CalcpHPCO2PH2SSTP(use_pH, UseH2Sgas, useEOS);

                    // ******重新分配 CO2aq、HCO3、CO3、H2Saq、HS 并重新计算离子强度
                    fmn();
                    CalcIonicStrength();
                    C2_PitzerActCoefs_T_P_ISt(gNeut, aH2O, TK, TC, PBar, Patm);

                    C5_CalcpHPCO2PH2SSTP(use_pH, UseH2Sgas, useEOS);

                    pHMeterReading = pH - DpHj;
                    goto label_3003;
                }
            }
        }
        if (useEOS == 1 || useEOS == 2)
        {
            if (simContext.SumofZMix[kk] == 0)
            {       // 仅当给出碳氢化合物成分时，才运行 HC 调节
                if (simContext.Run_Seawater_Mixing == 1 || simContext.Run_MixingTwoWells == 1)
                {
                    goto label_3002;
                }
                else
                {
                    if (simContext.nob == 1)
                    {
                        simContext.errmsg[13] = 14;  // 数组下标从0开始，14对应索引13
                        useEOS = 0;
                        goto label_3003;
                    }
                    else
                    {
                        if (mol_W > 0)
                        {
                            z_before_precipitation[1] = simContext.nTCO2Mix[kk] / mol_W;  // 原下标2对应索引1
                            z_before_precipitation[2] = simContext.nTH2SMix[kk] / mol_W;  // 原下标3对应索引2
                            z_before_precipitation[14] = 1 - z_before_precipitation[1] - z_before_precipitation[2];  // 原下标15对应索引14
                            simContext.Total_molesMix[kk] = mol_W;
                            simContext.nTCO2MixEOS[kk] = nTCO2;
                            simContext.nTH2SMixEOS[kk] = nTH2S;
                            goto label_3003;
                        }
                        else
                        {
                            simContext.errmsg[13] = 14;  // 数组下标从0开始，14对应索引13
                            useEOS = 0;
                            goto label_3003;
                        }
                    }
                }
            }
            else
            {
                total_moles = 1;
                //Compr(3)  beta(3)
                MultiPhaseFlash(mf_ParametersWereRead, mf_TCr, mf_PCr, mf_Omega, mf_MWgas, mf_kPr, mf_c0, mf_c1, TK, PBar, total_moles, z, tempgNeut,
                    aH2O, density, compositions, phi, Compr, beta, zOutput, mass_phase, MW_Phase, &No_Phases);

                if (beta[0] > 0 && beta[1] > 0)
                {
                    if (simContext.VoMix[kk] > 1.0 / 159.0 / 1000.0)
                    { // key to oil if vo is greater than 0
                        simContext.mass_o_Mix[kk] = 159 * simContext.VoMix[kk] * density[1]; // Kg
                        mol_o = simContext.mass_o_Mix[kk] * 1000 / MW_Phase[1];
                        mol_g = mol_o / beta[1] * beta[0];
                    }
                    else if (VgTP > 0.000001)
                    {    // 仅当 Vg 大于 0 且 Vo=0 时才允许
                        mol_g = VgTP * PBar * 1000 / (Compr[0] * RBar * TK);
                        mol_o = mol_g / beta[0] * beta[1];
                    }
                    else
                    { // 当没有给出油气体积时，只允许计算混合
                        if (simContext.Run_Seawater_Mixing == 1 || simContext.Run_MixingTwoWells == 1)
                        {
                            goto label_3002;
                        }
                        else
                        {
                            if (simContext.nob == 1)
                            {
                                mol_g = 0;
                                mol_o = 0;
                                simContext.errmsg[15] = 16;  // 数组下标从0开始，16对应索引15
                                useEOS = 0;
                                goto label_3003;
                            }
                            else
                            {
                                if (mol_W > 0)
                                {
                                    z_before_precipitation[1] = simContext.nTCO2Mix[kk] / mol_W;
                                    z_before_precipitation[2] = simContext.nTH2SMix[kk] / mol_W;
                                    z_before_precipitation[14] = 1 - z_before_precipitation[1] - z_before_precipitation[2];
                                    simContext.Total_molesMix[kk] = mol_W;
                                    simContext.nTCO2MixEOS[kk] = nTCO2;
                                    simContext.nTH2SMixEOS[kk] = nTH2S;
                                    goto label_3003;
                                }
                                else
                                {
                                    simContext.errmsg[15] = 16;  // 数组下标从0开始，16对应索引15
                                    useEOS = 0;
                                    goto label_3003;
                                }
                            }
                        }
                    }
                }

                if (beta[0] > 0 && beta[1] == 0)
                { // 仅存在气相
                    if (simContext.VoMix[kk] > 1.0 / 159.0 / 1000.0)
                    { // 如果存在石油量，则石油的关键
                        if ((z[7] + z[8] + z[9] + z[10] + z[11] + z[12]) > 0) { // 只有当 HC 大于 C4 时，才是石油的关键
                            simContext.mass_o_Mix[kk] = 159 * simContext.VoMix[kk] * density[0]; // Kg
                            mol_o = simContext.mass_o_Mix[kk] * 1000 / (MW_Phase)[0];
                            mol_g = 0;
                        }
                        else
                        {
                            mol_g = VgTP * PBar * 1000 / (Compr[0] * RBar * TK); // 如果不存在大于 C4 的 HC，则为气体键
                            mol_o = 0;
                        }
                    }
                    else if (VgTP > 0.000001) { // key to gas if vol of oil=0
                        mol_g = VgTP * PBar * 1000 / (Compr[0] * RBar * TK);
                        mol_o = 0;
                    }
                    else
                    { // 如果没有给出石油或天然气产量
                        if (simContext.Run_Seawater_Mixing == 1 || simContext.Run_MixingTwoWells == 1)
                        {
                            goto label_3002;
                        }
                        else
                        {
                            if (simContext.nob == 1)
                            {
                                mol_g = 0;
                                mol_o = 0;
                                simContext.errmsg[15] = 16;  // 数组下标从0开始，16对应索引15
                                useEOS = 0;
                                goto label_3003;
                            }
                            else
                            {
                                if (mol_W > 0)
                                {
                                    z_before_precipitation[1] = simContext.nTCO2Mix[kk] / mol_W;
                                    z_before_precipitation[2] = simContext.nTH2SMix[kk] / mol_W;
                                    z_before_precipitation[14] = 1 - z_before_precipitation[1] - z_before_precipitation[2];
                                    simContext.Total_molesMix[kk] = mol_W;
                                    simContext.nTCO2MixEOS[kk] = nTCO2;
                                    simContext.nTH2SMixEOS[kk] = nTH2S;
                                    goto label_3003;
                                }
                                else
                                {
                                    simContext.errmsg[15] = 16;  // 数组下标从0开始，16对应索引15
                                    useEOS = 0;
                                    goto label_3003;
                                }
                            }
                        }
                    }
                }

                if (beta[0] == 0 && beta[1] > 0)
                {
                    if (simContext.VoMix[kk] > 1.0 / 159.0 / 1000.0)
                    { // key to oil if vol of oil is present
                        simContext.mass_o_Mix[kk] = 159 * simContext.VoMix[kk] * density[1]; // Kg
                        mol_o = simContext.mass_o_Mix[kk] * 1000 / MW_Phase[1];
                        mol_g = 0;
                    }
                    else if (VgTP > 0.000001) { // key to gas if vol of oil=0
                        mol_g = VgTP * PBar * 1000 / (Compr[1] * RBar * TK);
                        mol_o = 0;
                    }
                    else
                    { // if neither oil or gas vol is given
                        if (simContext.Run_Seawater_Mixing == 1 || simContext.Run_MixingTwoWells == 1)
                        {
                            goto label_3002;
                        }
                        else
                        {
                            if (simContext.nob == 1)
                            {
                                simContext.errmsg[15] = 16;  // 数组下标从0开始，16对应索引15
                                useEOS = 0;
                                goto label_3003;
                            }
                            else
                            {
                                if (mol_W > 0)
                                {
                                    z_before_precipitation[1] = simContext.nTCO2Mix[kk] / mol_W;
                                    z_before_precipitation[2] = simContext.nTH2SMix[kk] / mol_W;
                                    z_before_precipitation[14] = 1 - z_before_precipitation[1] - z_before_precipitation[2];
                                    simContext.Total_molesMix[kk] = mol_W;
                                    simContext.nTCO2MixEOS[kk] = nTCO2;
                                    simContext.nTH2SMixEOS[kk] = nTH2S;
                                    goto label_3003;
                                }
                                else
                                {
                                    simContext.errmsg[15] = 16;  // 数组下标从0开始，16对应索引15
                                    useEOS = 0;
                                    goto label_3003;
                                }
                            }
                        }
                    }
                }

                mol_HC = mol_o + mol_g;

                if (useEOS == 1) {
                    nTCO2EOS = z[1] * mol_HC;  // 原z(2)对应z[1]
                    nTH2sEOS = z[2] * mol_HC;  // 原z(3)对应z[2]
                    nTCO2 = nTCO2EOS;
                    nTH2S = nTH2sEOS;
                }

                if (useEOS == 2) {
                    nTCO2EOS = nTCO2;
                    nTH2sEOS = nTH2S;
                    double zHC = z[0]; // 通过在 STP 条件下分解出（HCO3、CO3 和 HS）的总摩尔数来重新调整 nTCO2EOS 和 nTH2SEOS

                    for (iNG = 3; iNG < 14; iNG++) zHC = zHC + z[iNG]; // 总量中 z(碳氢化合物-CO2-H2s) 的摩尔分数


                    mol_HC = mol_HC * zHC + nTCO2EOS + nTH2sEOS;  // 如果 useEOS=2 nTCO2EOS，则修改 HC 的总摩尔数，计算气体、油和 simContext.CO2aq 中的 CO2 摩尔数（不包括 simContext.HCO3 和 CO3）

                    z[0] = z[0] * mol_HC;
                    z[1] = nTCO2EOS;
                    z[2] = nTH2sEOS;

                    for (iNG = 3; iNG < 14; iNG++) z[iNG] = z[iNG] * mol_HC;

                    for (iNG = 0; iNG < 14; iNG++) z[iNG] = z[iNG] / mol_HC;
                }
            } // end SumofZ>0
        }

        if (VW == 1.0 / 159.0 / 1000.0)
        {               // if water volume=0, skip Flash
            if (simContext.Run_Seawater_Mixing == 1 || simContext.Run_MixingTwoWells == 1)
            {
                goto label_3002;
            }
            else
            {
                simContext.errmsg[13] = 14;  // 数组下标从0开始，14对应索引13
                useEOS = 0;
                goto label_3003;
            }
        }
        //'*****使用EOS选项1或2
        //仅当给定实际的G/O/W体积时才允许进行Flash计算，例外情况是混合两口井和海水
        if (useEOS == 1 || useEOS == 2)
        {
            // 首先利用气相中的纯水蒸气计算储层成分
            if (usedryHC == 0)
            { // 湿烃情况，总水量=输入水量+HC中的水量
                mol_w3 = 0;
                while (mol_w3 <= 0)
                { // 确保添加足够的水以开始迭代
                    mol_W = mol_W * 1.1;

                    true_composition(TK, PBar, mol_HC, mol_W, aH2O, tempgNeut, nTCO2EOS, nTH2sEOS, useEOS, z, feed_Composition, &total_moles);

                    if (feed_Composition[0] < 0.0000001)
                    {  // 原feed_Composition(1)对应feed_Composition[0]
                        feed_Composition[0] = 0;
                    }

                    MultiPhaseFlash(mf_ParametersWereRead, mf_TCr, mf_PCr, mf_Omega, mf_MWgas, mf_kPr, mf_c0, mf_c1, TK, PBar, total_moles, feed_Composition,
                        tempgNeut, aH2O, density, compositions, phi, Compr, beta, zOutput, mass_phase, MW_Phase, &No_Phases);

                    if (compositions[14][3] > 0.5)
                    {  // 原compositions(15,4)对应compositions[14][3]
                        mol_w3 = total_moles * beta[2] * compositions[14][3]; // H2O moles in aqueous phase
                    }
                    else
                    {
                        mol_w3 = 0;
                    }
                }

                if (mol_w3 > 0)
                {
                    if (compositions[14][3] > 0.5) {
                        while (((mol_w3 - mol_w_Orig) / mol_w_Orig) * ((mol_w3 - mol_w_Orig) / mol_w_Orig) > 0.0001) {
                            mol_W = mol_W - mol_w3 + mol_w_Orig;
                            total_moles = mol_HC + mol_W;
                            true_composition(TK, PBar, mol_HC, mol_W, aH2O, tempgNeut, nTCO2EOS, nTH2sEOS, useEOS, z, feed_Composition, &total_moles);

                            MultiPhaseFlash(mf_ParametersWereRead, mf_TCr, mf_PCr, mf_Omega, mf_MWgas, mf_kPr, mf_c0, mf_c1, TK, PBar, total_moles, feed_Composition,
                                tempgNeut, aH2O, density, compositions, phi, Compr, beta, zOutput, mass_phase, MW_Phase, &No_Phases);
                            mol_w3 = total_moles * beta[2] * compositions[14][3];
                        }
                    }
                    else
                    {
                        simContext.errmsg[13] = 14;  // 数组下标从0开始，14对应索引13
                        useEOS = 0;
                        goto label_3003;
                    }
                }
                else
                {
                    simContext.errmsg[13] = 14;  // 数组下标从0开始，14对应索引13
                    useEOS = 0;
                    goto label_3003;
                }
            }
            else if (usedryHC == 1)
            { // dry hydrocarbon, only water from Input
                true_composition(TK, PBar, mol_HC, mol_W, aH2O, tempgNeut, nTCO2EOS, nTH2sEOS, useEOS, z, feed_Composition, &total_moles);

                if (feed_Composition[0] < 0.0000001)
                {
                    feed_Composition[0] = 0;
                }

                MultiPhaseFlash(mf_ParametersWereRead, mf_TCr, mf_PCr, mf_Omega, mf_MWgas, mf_kPr, mf_c0, mf_c1, TK, PBar, total_moles, feed_Composition,
                    tempgNeut, aH2O, density, compositions, phi, Compr, beta, zOutput, mass_phase, MW_Phase, &No_Phases);
            }

            simContext.nTCO2MixEOS[kk] = total_moles * zOutput[1];  // 原zOutput(2)对应zOutput[1]
            simContext.nTH2SMixEOS[kk] = total_moles * zOutput[2];  // 原zOutput(3)对应zOutput[2]
            mass_w = total_moles * beta[2] * compositions[14][3] * 0.01801528; // mass_w is only the aqueous phase H2O
            simContext.Total_molesMix[kk] = total_moles;
            simContext.mass_w_Mix[kk] = mass_w;

            if (No_Phases == 3)
            { // Output oil and gas density from flash calculation if useEOS=1 or useEOS=2
                simContext.GasDensityMix[kk] = density[0] * 1000;  // 原density(1)对应density[0]
                simContext.OilDensityMix[kk] = density[1];         // 原density(2)对应density[1]
            }

            if (No_Phases == 2)
            { // Only g/o, G/W or O/W phases exist
                if (beta[2] == 0)
                {
                    //本例中不存在水相。增加水量并重新运行。程序中止。
                    printf("Aqueous phase does not existed in this case. Increase water volume and run again. Program abort.\n"); // Aqueous phase does not exist
                    exit(1);  // 替代VB的End
                }
                else
                {
                    if (density[0] > 0 && density[0] < 0.3) simContext.GasDensityMix[kk] = density[0] * 1000;// phase 1 is gas

                    else if (density[0] > 0 && density[0] > 0.3) simContext.OilDensityMix[kk] = density[0]; // phase 1 is oil

                    else if (density[1] > 0 && density[1] < 0.3) simContext.GasDensityMix[kk] = density[1] * 1000;// Phase 2 is gas

                    else if (density[1] > 0 && density[1] > 0.3) simContext.OilDensityMix[kk] = density[1];

                }
            }

            if (No_Phases == 1)
            {
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

        if (useEOS == 3)
        {
            // 仅在给定天然气和石油产量时运行 useEOS；第一次迭代涵盖 usedryHC=0 或 1
            pseudo_composition(data->API, SGG, VgTP, mol_o, mol_W, TK, PBar, aH2O, tempgNeut, nTCO2, nTH2S, simContext.yCO2,
                simContext.yH2S, YH2O, &total_moles, feed_Composition, &mol_HC);

            if (feed_Composition[0] < 0.0000001) {
                feed_Composition[0] = 0;
            }

            SumofZ = feed_Composition[0];
            for (iNG = 1; iNG < 14; iNG++) {
                if (feed_Composition[iNG] < 0.0000001) feed_Composition[iNG] = 0;
                SumofZ = SumofZ + feed_Composition[iNG];
            }

            if (SumofZ < 0.000001) {
                simContext.errmsg[13] = 14;
                useEOS = 0;
            }

            MultiPhaseFlash(mf_ParametersWereRead, mf_TCr, mf_PCr, mf_Omega, mf_MWgas, mf_kPr, mf_c0, mf_c1, TK, PBar, total_moles, feed_Composition, tempgNeut,
                aH2O, density, compositions, phi, Compr, beta, zOutput, mass_phase, MW_Phase, &No_Phases);

            if (usedryHC == 0)
            {
                mol_w3 = total_moles * beta[2] * compositions[14][3];

                while (mol_w3 <= 0)
                { // 确保添加足够的水以开始迭代
                    mol_W = mol_W * 1.1;
                    pseudo_composition(data->API, SGG, VgTP, mol_o, mol_W, TK, PBar, aH2O, tempgNeut, nTCO2, nTH2S, simContext.yCO2,
                        simContext.yH2S, YH2O, &total_moles, feed_Composition, &mol_HC);

                    MultiPhaseFlash(mf_ParametersWereRead, mf_TCr, mf_PCr, mf_Omega, mf_MWgas, mf_kPr, mf_c0, mf_c1, TK, PBar, total_moles, feed_Composition,
                        tempgNeut, aH2O, density, compositions, phi, Compr, beta, zOutput, mass_phase, MW_Phase, &No_Phases);

                    if (compositions[14][3] > 0.5)  mol_w3 = total_moles * beta[2] * compositions[14][3]; // 水相中的 H2O 摩尔数
                    else mol_w3 = 0;

                }

                if (mol_w3 > 0)
                {
                    while (((mol_w3 - mol_w_Orig) / mol_w_Orig) * ((mol_w3 - mol_w_Orig) / mol_w_Orig) > 0.0001)
                    {
                        mol_W = mol_W - mol_w3 + mol_w_Orig;
                        total_moles = total_moles - mol_w3 + mol_w_Orig;

                        // recalculate z(i) into reservoir composition to use true_composition sub
                        for (iNG = 0; iNG < 14; iNG++)
                        {
                            z[iNG] = zOutput[iNG] / (1 - zOutput[14]);
                        }

                        pseudo_composition(data->API, SGG, VgTP, mol_o, mol_W, TK, PBar, aH2O, tempgNeut,
                            nTCO2, nTH2S, simContext.yCO2, simContext.yH2S, YH2O, &total_moles, feed_Composition, &mol_HC);

                        MultiPhaseFlash(mf_ParametersWereRead, mf_TCr, mf_PCr, mf_Omega, mf_MWgas, mf_kPr, mf_c0, mf_c1, TK, PBar, total_moles, feed_Composition, tempgNeut,
                            aH2O, density, compositions, phi, Compr, beta, zOutput, mass_phase, MW_Phase, &No_Phases);
                        mol_w3 = total_moles * beta[2] * compositions[14][3];
                    }
                }

                simContext.nTCO2MixEOS[kk] = total_moles * zOutput[1];
                simContext.nTH2SMixEOS[kk] = total_moles * zOutput[2];
                mass_w = total_moles * beta[2] * compositions[14][3] * 0.01801528; // mass_w 仅为水相 H2O
                simContext.Total_molesMix[kk] = total_moles;
                simContext.mass_w_Mix[kk] = mass_w;
            }
            else {
                simContext.Total_molesMix[kk] = total_moles;
            }
        } // Correspond to useEOS=3



        // ***** Run speciation and pH calculation for useEOS>0
        nTCO2_before_precipitation = nTCO2EOS;
        nTH2S_before_precipitation = nTH2sEOS;
        Total_moles_before_precipitation = total_moles;

        for (int iz = 0; iz < 15; iz++) {
            z_before_precipitation[iz] = zOutput[iz];
            simContext.zMix[kk][iz] = zOutput[iz];
        }

        if (usedryHC == 1)
        {
            // 如果假设碳氢化合物为干碳氢化合物，则重新计算 T、P 处的真实水浓度
            for (int c = 0; c < NumCat; c++) mc[c] = mc[c] * mol_w_Orig * 0.01801528 / mass_w;

            for (int a = 0; a < NumAn; a++)  ma[a] = ma[a] * mol_w_Orig * 0.01801528 / mass_w;

            for (int n = 0; n < NumNeut; n++) mn[n] = mn[n] * mol_w_Orig * 0.01801528 / mass_w;

            Alk = Alk * mol_w_Orig * 0.01801528 / mass_w;
            TAc = TAc * mol_w_Orig * 0.01801528 / mass_w;
            TNH4 = TNH4 * mol_w_Orig * 0.01801528 / mass_w;
            TH3BO3 = TH3BO3 * mol_w_Orig * 0.01801528 / mass_w;
            TH2Saq = TH2Saq * mol_w_Orig * 0.01801528 / mass_w;
            TH4SiO4 = TH4SiO4 * mol_w_Orig * 0.01801528 / mass_w;

            simContext.HstpMix[kk] = mc[iH];
            simContext.NaMix[kk] = mc[iNa];
            simContext.KMix[kk] = mc[iK];
            simContext.MgMix[kk] = mc[iMg];
            simContext.CaMix[kk] = mc[iCa];
            TCa = mc[iCa];
            simContext.SrMix[kk] = mc[iSr];
            simContext.BaMix[kk] = mc[iBa];
            simContext.ZnMix[kk] = mc[iZn];
            simContext.PbMix[kk] = mc[iPb];
            simContext.RaMix[kk] = mc[iRa];

            simContext.OHstpMix[kk] = ma[iOH];
            simContext.ClMix[kk] = ma[iCl];
            simContext.ACstpMix[kk] = ma[iAc];
            simContext.NH4STPMix[kk] = mc[iNH4];
            simContext.H2BO3stpMix[kk] = ma[iH2BO3];
            simContext.HCO3stpMix[kk] = ma[iHCO3];
            simContext.CO3stpMix[kk] = ma[iCO3];
            simContext.SO4Mix[kk] = ma[iSO4];
            simContext.HSstpMix[kk] = ma[iHS];
            simContext.FMix[kk] = ma[intF];
            simContext.BrMix[kk] = ma[iBr];

            simContext.rho25CMix[kk] = simContext.rho25c;

            simContext.H3SiO4Mix[kk] = ma[iH3SiO4];
            simContext.H2SiO4Mix[kk] = ma[iH2SiO4];

            simContext.NH3Mix[kk] = mn[iNH3];
            simContext.H4SiO4Mix[kk] = mn[iH4SiO4aq];
            simContext.H3BO3Mix[kk] = mn[iH3BO3];

            simContext.CO2aqMix[kk] = mn[iCO2aq];
            simContext.H2SaqMix[kk] = mn[iH2Saq];
            simContext.HACaqMix[kk] = mn[iHAcaq];

            simContext.AlkMix[kk] = Alk;
            simContext.TAcMix[kk] = TAc;
            simContext.TH2SaqMix[kk] = TH2Saq;
            simContext.TH4SiO4Mix[kk] = TH4SiO4;
            simContext.TNH4Mix[kk] = TNH4;
            simContext.TH3BO3Mix[kk] = TH3BO3;
        }

        // **** flash and convert to STP condition
        // Call C4_EOS_TCO2_SSPEquilCalcs(0, 2, 2, KspCalcite) 'calculate Separator or STP pH if useEOS>0
        // *********
        // mt = fmn  'Assign neutral species and simContext.HCO3, simContext.CO3, simContext.HS and recalculate ionic strength, Pitzer Coeff and pH
        // Call CalcIonicStrength
        // Call C2_PitzerActCoefs_T_P_ISt(gNeut, aH2O, TK, TC, PBar, Patm)
        // Call C4_EOS_TCO2_SSPEquilCalcs(0, 5, 2, KspCalcite) '(ppt_or_not, iMetals, iGas, Ksp)

        //mt = fTPFunc(0);
        fTPFunc(0);
        C4_SSPEquilCalcs(0, 5, 2, KspCalcite);  // (ppt_or_not, iMetals, iGas, Ksp)
        fmn();        // simContext.CO2aq, simContext.HCO3, simContext.CO3, simContext.H2Saq, simContext.HS
        C2_PitzerActCoefs_T_P_ISt(gNeut, aH2O, TK, TC, PBar, Patm);
        C4_SSPEquilCalcs(0, 5, 2, KspCalcite);  // (ppt_or_not, iMetals, iGas, Ksp)

        pHMeterReading = pH - DpHj; // at separator T&P or STP
        rhoTP = CalcRhoTP(TK, TC, PBar, Patm); // at separator T&P or STP

        //下面三个数组在B2中初始化：
        //molc(NumCat, nob_Input + nob_InputII),mola(NumAn, nob_Input + nob_InputII), moln(NumNeut, nob_Input + nob_InputII) ' Convert all ions to moles
        for (int c = 0; c < NumCat; c++)
            molc[c][kk] = mc[c] * mass_w;

        for (int a = 0; a < NumAn; a++)
            mola[a][kk] = ma[a] * mass_w;

        for (int n = 0; n < NumNeut; n++)
            moln[n][kk] = mn[n] * mass_w;

        simContext.molAlk[kk] = Alk * mass_w;
        simContext.molTAC[kk] = TAc * mass_w;
        simContext.molTNH4[kk] = TNH4 * mass_w;
        simContext.molTH3BO3[kk] = TH3BO3 * mass_w;
        simContext.molTH2Saq[kk] = TH2Saq * mass_w;
        simContext.molTH4SiO4[kk] = TH4SiO4 * mass_w;

    }

label_3002:
    if (simContext.Run_Seawater_Mixing == 1 || simContext.Run_MixingTwoWells == 1) {
        if (simContext.Run_Seawater_Mixing == 1) {
            // 源代码是kk == 1，但我们的Kk是从0开始的
            if (simContext.LoopMixing == 1 && kk == 0) { // when Fr of brine#1=0
                for (iNG = 0; iNG < 15; iNG++)
                    z_before_precipitation[iNG] = 0;

                simContext.Total_molesMix[kk] = 0;
                simContext.nTCO2MixEOS[kk] = 0;
                simContext.nTH2SMixEOS[kk] = 0;
                simContext.mass_w_Mix[kk] = 0;
                useEOS = 0;
                goto label_3003;
            }

            if (simContext.LoopMixing == 11 && kk == 1) { // when Fr of seawater=0
                for (iNG = 0; iNG < 15; iNG++)
                    z_before_precipitation[iNG] = 0;

                simContext.Total_molesMix[kk] = 0;
                simContext.nTCO2MixEOS[kk] = 0;
                simContext.nTH2SMixEOS[kk] = 0;
                simContext.mass_w_Mix[kk] = 0;
                goto label_3003;
            }

            if (kk == 1) { // Set EOS parameter for seawater
                simContext.Total_molesMix[kk] = mol_W;
                z_before_precipitation[1] = simContext.nTCO2Mix[kk] / simContext.Total_molesMix[kk];
                z_before_precipitation[2] = simContext.nTH2SMix[kk] / simContext.Total_molesMix[kk];
                z_before_precipitation[14] = 1 - z_before_precipitation[1] - z_before_precipitation[2];
                goto label_3003;
            }
        }

        if (simContext.Run_MixingTwoWells == 1) {
            if (simContext.LoopMixing == 1 && kk == 0) { // when Fr Brine#1 =0
                for (iNG = 0; iNG < 15; iNG++)
                    z_before_precipitation[iNG] = 0;

                simContext.Total_molesMix[kk] = 0;
                simContext.nTCO2MixEOS[kk] = 0;
                simContext.nTH2SMixEOS[kk] = 0;
                simContext.mass_w_Mix[kk] = 0;
                goto label_3003;
            }

            if (simContext.LoopMixing == 11 && kk == 1) { // When Fr brine#2 =0
                for (iNG = 0; iNG < 15; iNG++)
                    z_before_precipitation[iNG] = 0;

                simContext.Total_molesMix[kk] = 0;
                simContext.nTCO2MixEOS[kk] = 0;
                simContext.nTH2SMixEOS[kk] = 0;
                simContext.mass_w_Mix[kk] = 0;
                goto label_3003;
            }

            if (VgTP == 0.000001 && VO == 1.0 / 159.0 / 1000.0 && VW > 1.0 / 159.0 / 1000.0) { // when only water present for either brine 1 or 2
                simContext.Total_molesMix[kk] = mol_W;
                z_before_precipitation[1] = simContext.nTCO2Mix[kk] / simContext.Total_molesMix[kk];
                z_before_precipitation[2] = simContext.nTH2SMix[kk] / simContext.Total_molesMix[kk];
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
        simContext.zMix[kk][iNG] = z_before_precipitation[iNG];


    // ***** QC calculation
    if (simContext.RunStat == 1)
    {
        if (simContext.RunQualityControlChecks == 1) {
            // Worksheets(mySheet).Activate 
            QualityControlCalculations(kk, j);
        }
    }
    else
    {
        // kk在vb中是1开始的，这里我们的kk是从0开始，导致了 nob_input=1时，kk为0和1都符合条件，因此手动进行偏移
        if (simContext.RunQualityControlChecks == 1 && kk + 1 <= simContext.nob_Input) { // only run QC if requested from Input Sheet.
            // Worksheets(mySheet).Activate
            QualityControlCalculations(kk, j);
            if (kk == simContext.nob) {
                goto label_123;   // Last Brine QC has been printed to Input Sheet; exit calculations.
            }
        }
        if (RunQualityControlChecks_II == 1 && kk + 1 >= simContext.nob_Input) { // only run QC if requested from InputII Sheet.
            // Worksheets(mySheet).Activate
            QualityControlCalculations(kk, j);
            if (kk == simContext.nob_Input + simContext.nob_InputII) {
                exit(0);//vb源码 ： End 
            }
        }
    }

    // *****End of QC calculation
    // ***** 将 pH、PCO2、PH2S、TH2Saq 打印到输入表
    // If RunShellMultiflash <> 1 Then
    if (use_pH == 0 && useEOS == 0) {        // 然后，必须指定yCO2stp。需要计算pH值和yH2Sstp。
        if (simContext.Run_Seawater_Mixing == 0 && simContext.Run_MixingTwoWells == 0 && simContext.RunMultiMix == 0 && simContext.Run10TestCases == 0 && simContext.RunWhatIf == 0 && simContext.RunStat == 0) {
            if (simContext.Run_CalcConcFactor == 1) {
                // 在C中，可能需要调用特定的输出函数
                //待定：Worksheets(MySheetMix).Cells(34, 8) = pHMeterReading
                int bbba = pHMeterReading;
            }
            else {
                //待定： Worksheets(mySheet).Cells(34, j + 2) = pHMeterReading
                int bbba = pHMeterReading;
            }
        }
    }

    if (use_pH == 1 && useEOS == 0) {                            // Use the pH measured at STP.
        if (simContext.RunWhatIf != 1 && simContext.RunStat == 0) {
            // 待定： Worksheets(mySheet).Cells(31, j + 2) = simContext.yCO2 * 100

        }
    }

    if (use_pH == 2 && useEOS == 0) {                            // Use the pH measured at STP.
        if (simContext.RunWhatIf != 1 && simContext.RunStat == 0) {
            // 待定： Worksheets(mySheet).Cells(24, j + 2) = Alk * (61019 * (simContext.rho25c - TDS / 1000000.0))
        }
        simContext.AlkMix[kk] = Alk;
    }

    if (use_pH == 3 && useEOS == 0) {         // use TCO2 to calculate pH
        if (simContext.Run_Seawater_Mixing == 0 && simContext.Run_MixingTwoWells == 0 && simContext.RunMultiMix == 0 && simContext.Run10TestCases == 0 && simContext.RunWhatIf == 0 && simContext.RunStat == 0) {
            pHMeterReading = pH - DpHj;
            if (simContext.Run_CalcConcFactor == 1) {
                // 待定： Worksheets(MySheetMix).Cells(34, 8) = pHMeterReading
            }
            else {
                //待定： Worksheets(mySheet).Cells(34, j + 2) = pHMeterReading
            }
        }
        if (simContext.RunWhatIf != 1) {
            //待定：  Worksheets(mySheet).Cells(31, j + 2) = simContext.yCO2 * 100
        }
    }

    if (useEOS != 0) {
        if (simContext.Run_Seawater_Mixing == 0 && simContext.Run_MixingTwoWells == 0 && simContext.RunMultiMix == 0 && simContext.Run10TestCases == 0 && simContext.RunWhatIf == 0 && simContext.RunStat == 0) {
            if (simContext.Run_CalcConcFactor == 1) {
                //待定：  Worksheets(MySheetMix).Cells(34, 8) = pHMeterReading
            }
            else {
                //待定：  Worksheets(mySheet).Cells(34, j + 2) = pHMeterReading
            }
        }
    }

    if (simContext.UseMolal == 1) {
        //待定：  Worksheets(mySheet).Cells(29, j + 2) = CalculatedTDSMix(kk)
    }

    if (mass_MeOH > 0) {
        xMeOH = (mass_MeOH / 32.0) / (mass_MeOH / 32.0 + mass_w / 18.0); // ?????????????????????????????
    }

    if (mass_MEG > 0) {
        xMEG = (mass_MEG / 62.07) / (mass_MEG / 62.07 + mass_w / 18.0);
    }

    if (mass_MeOH > 0) {
        IStCosolvent = ISt * mass_w / (mass_w + mass_MeOH);
    }

    if (mass_MEG > 0) {
        IStCosolvent = ISt * mass_w / (mass_w + mass_MEG);
    }

label_123:
    // End Sub - 在C中，这可能是函数的结束
    // 假设这是一个函数的结尾
    return;


}




void C2_PitzerActCoefsConstants()
{
    b0[iH][iBr] = 0.2085; b0[iNa][iAc] = 0.1426; b0[iNa][iHS] = -0.103;
    b0[iNa][iBr] = 0.0973; b0[iNa][iH2BO3] = -0.0427;

    b0[iK][iAc] = 0.15298; b0[iK][iBr] = 0.0569; b0[iK][iH2BO3] = 0.035;

    b0[iMg][iHS] = 0.466; b0[iMg][iBr] = 0.5769 * 3.0 / 4.0;

    b0[iCa][iOH] = -0.1747; b0[iCa][iCO3] = 0.16; b0[iCa][iHS] = 0.069;
    b0[iCa][iBr] = 0.5088 * 3.0 / 4.0;

    b0[iSr][iHCO3] = 0.12; b0[iSr][iCO3] = 0;
    b0[iSr][iBr] = 0.4415 * 3.0 / 4.0;

    b0[iBa][iOH] = 0.17175; b0[iBa][iBr] = 0.4194 * 3.0 / 4.0;

    b0[iFe][iOH] = 0.17175; b0[iFe][iCl] = 0.35011;
    b0[iFe][iAc] = 0.28725; b0[iFe][iCO3] = 1.919;
    b0[iFe][iSO4] = 0.2568;

    b0[iZn][iOH] = 0.17175; b0[iZn][iCl] = 0.0887;
    b0[iZn][iSO4] = 0.18404;

    b0[iZn][iBr] = 0.6213 * 3.0 / 4.0;
    b0[iZn][iHS] = -0.5; b0[iZn][iCl] = -0.5;

    b1[iH][iBr] = 0.3477; b1[iNa][iAc] = 0.3237; b1[iNa][iHS] = 0.884;
    b1[iNa][iBr] = 0.2791; b1[iNa][iH2BO3] = 0.089;

    b1[iK][iAc] = 0.34195; b1[iK][iBr] = 0.2212; b1[iK][iH2BO3] = 0.14;

    b1[iMg][iHS] = 2.264; b1[iMg][iBr] = 2.337 * 3.0 / 4.0;

    b1[iCa][iOH] = -0.2303; b1[iCa][iCO3] = 2.1;
    b1[iCa][iHS] = 2.264; b1[iCa][iBr] = 2.151 * 3.0 / 4.0;

    b1[iSr][iBr] = 2.282 * 3.0 / 4.0;

    b1[iBa][iOH] = 1.2; b1[iBa][iBr] = 2.093 * 3.0 / 4.0;

    b1[iFe][iCl] = 1.40092; b1[iFe][iHCO3] = 14.76;
    b1[iFe][iCO3] = -5.134; b1[iFe][iSO4] = 3.063;

    b1[iZn][iSO4] = 3.031; b1[iZn][iBr] = 2.179 * 3.0 / 4.0;
    b1[iZn][iOH] = 4.72441418675659e-02;
    b1[iZn][iCl] = -6.03067094634178; b1[iZn][iHS] = -7.44795278611779;

    b2[iCa][iOH] = -5.72; b2[iCa][iCO3] = -69;
    b2[iFe][iCO3] = -274; b2[iFe][iSO4] = -42;
    b2[iZn][iSO4] = -27.709;

    CPhi[iH][iBr] = 0.00152;

    CPhi[iNa][iAc] = -0.00629; CPhi[iNa][iBr] = 0.00116;
    CPhi[iNa][iH2BO3] = 0.0114;

    CPhi[iK][iAc] = -0.00474; CPhi[iK][iHCO3] = -0.008;
    CPhi[iK][iCO3] = -0.0015; CPhi[iK][iCO3] = 0.0005;

    CPhi[iK][iBr] = -0.0018;

    CPhi[iMg][iBr] = 0.00589 * 0.5303;

    CPhi[iFe][iCl] = -0.01412; CPhi[iFe][iSO4] = 0.0209;

    CPhi[iZn][iSO4] = 0.03286;
    CPhi[iZn][iHS] = 0.5; CPhi[iZn][iCl] = -0.5;

    CPhi[iCa][iBr] = -0.00485 * 0.5303;
    CPhi[iSr][iBr] = 0.00231 * 0.5303;
    CPhi[iBa][iBr] = -0.03009 * 0.5303;
    CPhi[iZn][iBr] = -0.2035 * 0.5303;

    Tccp[iH][iNa] = 0.036; Tccp[iNa][iH] = 0.036;
    Tccp[iH][iK] = 0.005; Tccp[iK][iH] = 0.005;
    Tccp[iH][iMg] = 0.005; Tccp[iMg][iH] = 0.005;
    Tccp[iH][iCa] = 0.092; Tccp[iCa][iH] = 0.092;
    Tccp[iH][iSr] = 0.07; Tccp[iSr][iH] = 0.07;

    Tccp[iNa][iK] = -0.012; Tccp[iK][iNa] = -0.012;
    Tccp[iNa][iMg] = 0.07; Tccp[iMg][iNa] = 0.07;

    Tccp[iK][iCa] = 0.032; Tccp[iCa][iK] = 0.032;
    Tccp[iMg][iCa] = 0.1244; Tccp[iCa][iMg] = 0.1244;

    Yccpa[iH][iNa][iCl] = -0.004; Yccpa[iNa][iH][iCl] = -0.004;
    Yccpa[iH][iK][iCl] = -0.011; Yccpa[iK][iH][iCl] = -0.011;
    Yccpa[iH][iK][iSO4] = 0.197; Yccpa[iK][iH][iSO4] = 0.197;
    Yccpa[iH][iMg][iCl] = -0.011; Yccpa[iMg][iH][iCl] = -0.011;
    Yccpa[iH][iMg][iSO4] = -0.0178; Yccpa[iMg][iH][iSO4] = -0.0178;
    Yccpa[iH][iCa][iCl] = -0.015; Yccpa[iCa][iH][iCl] = -0.015;
    Yccpa[iH][iSr][iCl] = 0.01; Yccpa[iSr][iH][iCl] = 0.01;
    Yccpa[iH][iSr][iSO4] = 0.03; Yccpa[iSr][iH][iSO4] = 0.03;

    Yccpa[iNa][iK][iCl] = -0.0018; Yccpa[iK][iNa][iCl] = -0.0018;
    Yccpa[iNa][iK][iHCO3] = -0.003; Yccpa[iK][iNa][iHCO3] = -0.003;
    Yccpa[iNa][iK][iCO3] = 0.003; Yccpa[iK][iNa][iCO3] = 0.003;
    Yccpa[iNa][iK][iSO4] = -0.01; Yccpa[iK][iNa][iSO4] = -0.01;

    Yccpa[iNa][iMg][iCl] = -0.008; Yccpa[iMg][iNa][iCl] = -0.008;
    Yccpa[iNa][iMg][iSO4] = -0.015; Yccpa[iMg][iNa][iSO4] = -0.015;

    Yccpa[iNa][iCa][iSO4] = -0.055; Yccpa[iCa][iNa][iSO4] = -0.055;

    Yccpa[iK][iMg][iCl] = -0.022; Yccpa[iMg][iK][iCl] = -0.022;
    Yccpa[iK][iMg][iSO4] = -0.048; Yccpa[iMg][iK][iSO4] = -0.048;

    Yccpa[iK][iCa][iCl] = -0.025; Yccpa[iCa][iK][iCl] = -0.025;

    Yccpa[iMg][iCa][iCl] = -0.0238; Yccpa[iCa][iMg][iCl] = -0.0238;
    Yccpa[iMg][iCa][iSO4] = 0.024; Yccpa[iCa][iMg][iSO4] = 0.024;

    Taap[iOH][iCl] = -0.05; Taap[iCl][iOH] = -0.05;
    Taap[iOH][iCO3] = 0.1; Taap[iCO3][iOH] = 0.1;
    Taap[iOH][iSO4] = -0.013; Taap[iSO4][iOH] = -0.013;

    Taap[iCl][iHCO3] = 0.0359; Taap[iHCO3][iCl] = 0.0359;
    Taap[iCl][iCO3] = -0.053; Taap[iCO3][iCl] = -0.053;
    Taap[iHCO3][iCO3] = -0.04; Taap[iCO3][iHCO3] = -0.04;

    Taap[iHCO3][iSO4] = 0.01; Taap[iSO4][iHCO3] = 0.01;
    Taap[iCO3][iSO4] = 0.02; Taap[iSO4][iCO3] = 0.02;

    Taap[iH2BO3][iCl] = -0.065; Taap[iCl][iH2BO3] = -0.065;
    Taap[iH2BO3][iSO4] = -0.012; Taap[iSO4][iH2BO3] = -0.012;

    Yaapc[iOH][iCl][iNa] = -0.006; Yaapc[iCl][iOH][iNa] = -0.006;
    Yaapc[iOH][iCl][iK] = -0.006; Yaapc[iCl][iOH][iK] = -0.006;

    Yaapc[iOH][iCl][iCa] = -0.025; Yaapc[iCl][iOH][iCa] = -0.025;

    Yaapc[iOH][iCO3][iNa] = -0.017; Yaapc[iCO3][iOH][iNa] = -0.017;
    Yaapc[iOH][iCO3][iK] = -0.01; Yaapc[iCO3][iOH][iK] = -0.01;

    Yaapc[iOH][iSO4][iNa] = -0.009; Yaapc[iSO4][iOH][iNa] = -0.009;
    Yaapc[iOH][iSO4][iK] = -0.05; Yaapc[iSO4][iOH][iK] = -0.05;

    Yaapc[iCl][iHCO3][iK] = -0.006; Yaapc[iHCO3][iCl][iK] = -0.006;
    Yaapc[iCl][iHCO3][iMg] = -0.025; Yaapc[iHCO3][iCl][iMg] = -0.025;

    Yaapc[iCl][iCO3][iNa] = 0.016; Yaapc[iCO3][iCl][iNa] = 0.016;

    Yaapc[iCl][iCO3][iK] = 0.004; Yaapc[iCO3][iCl][iK] = 0.004;

    Yaapc[iCl][iSO4][iMg] = -0.004; Yaapc[iSO4][iCl][iMg] = -0.004;
    Yaapc[iCl][iSO4][iCa] = -0.018; Yaapc[iSO4][iCl][iCa] = -0.018;

    Yaapc[iHCO3][iCO3][iNa] = 0.002; Yaapc[iCO3][iHCO3][iNa] = 0.002;
    Yaapc[iHCO3][iCO3][iK] = 0.012; Yaapc[iCO3][iHCO3][iK] = 0.012;

    Yaapc[iHCO3][iSO4][iNa] = -0.005; Yaapc[iSO4][iHCO3][iNa] = -0.005;
    Yaapc[iHCO3][iSO4][iMg] = -0.161; Yaapc[iSO4][iHCO3][iMg] = -0.161;

    Yaapc[iCO3][iSO4][iNa] = -0.005; Yaapc[iSO4][iCO3][iNa] = -0.005;
    Yaapc[iCO3][iSO4][iK] = -0.009; Yaapc[iSO4][iCO3][iK] = -0.009;

    Yaapc[iH2BO3][iCl][iNa] = -0.0073; Yaapc[iCl][iH2BO3][iNa] = -0.0073;

    Lnc[iH3BO3][iNa] = -0.097;

    Lna[iHAcaq][iCl] = 0.076;
    Lna[iH3BO3][iCl] = 0.091;
    Lna[iH3BO3][iSO4] = 0.018;
}



void B1_InitializeIndices()
{
    mf_ParametersWereRead = false;
    //nComponents = 15;
    int i;

    NumCat = 12;
    NumAn = 13;
    NumNeut = 8;
    int NumMean = 6;

    /* ----------- Cation charge array ----------- */
    ChCat[iH] = 1;
    ChCat[iNa] = 1;
    ChCat[iK] = 1;
    ChCat[iMg] = 2;
    ChCat[iCa] = 2;
    ChCat[iSr] = 2;
    ChCat[iBa] = 2;
    ChCat[iFe] = 2;
    ChCat[iZn] = 2;
    ChCat[iPb] = 2;
    ChCat[iNH4] = 1;
    ChCat[iRa] = 2;

    /* ----------- Anion charge array ----------- */
    ChAn[iOH] = -1;
    ChAn[iCl] = -1;
    ChAn[iAc] = -1;
    ChAn[iHCO3] = -1;
    ChAn[iCO3] = -2;
    ChAn[iSO4] = -2;
    ChAn[iHS] = -1;
    ChAn[intF] = -1;
    ChAn[iBr] = -1;
    ChAn[iH2BO3] = -1;
    ChAn[iH3SiO4] = -1;
    ChAn[iH2SiO4] = -2;
    ChAn[iSion] = -2;

    /* ----------- Neutral aquatic indexes ----------- */
    /* (You already offset these so direct assignment is OK) */

    /* ----------- Set multiplicities gNCat, gNAn etc. ----------- */
    for (i = 0; i < NumCat; i++)  gNCat[i] = 1;
    for (i = 0; i < NumAn; i++)  gNAn[i] = 1;
    for (i = 0; i < NumNeut; i++)  gNNeut[i] = 1;
    for (i = 0; i < NumMean; i++)  gNMean[i] = 1;

    dSIMeOHBar = 0;
    dSIMEGBar = 0;
    dSIMeOHcal = 0;
    dSIMEGcal = 0;
    dSIMeOHHal = 0;
    dSIMEGHal = 0;

    aNH2O = 1;

    C2_PitzerActCoefsConstants();

    gL[iCO2o] = 1;
    gL[iCH4o] = 1;
    gL[iH2So] = 1;

    /* ----------- Water EOS constants (IAPWS IF-97) ----------- */
    rgas_water = 461.526;
    tc_water = 647.096;
    pc_water = 220.64;
    dc_water = 322.0;


    /* ---------- MW of ions ---------- */
    MWCat[iH] = 1.008;
    MWCat[iNa] = 22.99;
    MWCat[iK] = 39.102;
    MWCat[iMg] = 24.305;
    MWCat[iCa] = 40.08;
    MWCat[iSr] = 87.62;
    MWCat[iBa] = 137.33;
    MWCat[iFe] = 55.847;
    MWCat[iZn] = 65.38;
    MWCat[iPb] = 207.2;
    MWCat[iNH4] = 18.039;
    MWCat[iRa] = 226.0254;

    MWAn[iOH] = 17.007;
    MWAn[iCl] = 35.45;
    MWAn[iAc] = 59.054;
    MWAn[iHCO3] = 61.019;
    MWAn[iCO3] = 60.019;
    MWAn[iSO4] = 96.064;
    MWAn[iHS] = 33.073;
    MWAn[intF] = 18.998;
    MWAn[iBr] = 79.904;
    MWAn[iH2BO3] = 60.825;
    MWAn[iH3SiO4] = 59.07;
    MWAn[iH2SiO4] = 58.06;
    MWAn[iSion] = 32.065;

    MWNeut[iCH4aq] = 16;
    MWNeut[iCO2aq] = 62.03;
    MWNeut[iH2Saq] = 34.08;
    MWNeut[iHAcaq] = 60.054;
    MWNeut[iH4SiO4aq] = 60.08;
    MWNeut[iFeSaq] = 87.912;
    MWNeut[iH3BO3] = 61.833;
    MWNeut[iNH3] = 17.031;

    /* ---------- Partial molar volumes ---------- */
    V0_c[iH] = 0;
    V0_c[iNa] = -1.21;
    V0_c[iK] = 9.02;
    V0_c[iMg] = -21.17;
    V0_c[iCa] = -17.85;
    V0_c[iNH4] = 17.87;
    V0_c[iSr] = -18.16;
    V0_c[iBa] = -12.47;
    V0_c[iFe] = -24.7;
    V0_c[iZn] = -21.6;
    V0_c[iPb] = -15.5;
    V0_c[iRa] = -12.47;

    V0_a[iOH] = -4.04;
    V0_a[iCl] = 17.83;
    V0_a[iAc] = 40.46;
    V0_a[iHCO3] = 23.4;
    V0_a[iCO3] = -4.3;
    V0_a[iSO4] = 35.67;
    V0_a[iHS] = 20.71;
    V0_a[intF] = -1.16;
    V0_a[iBr] = 24.71;
    V0_a[iH2BO3] = 21.84;
    V0_a[iSion] = -8.2;

    V0_n[iHAcaq] = 51.94;
    V0_n[iCO2aq] = 50.78;
    V0_n[iH2Saq] = 35.71;
    V0_n[iCH4aq] = 0;
    V0_n[iH4SiO4aq] = 0;
    V0_n[iNH3] = 19.74;
    V0_n[iH3BO3] = 31.6;
    V0_n[iFeSaq] = 35.71;

    /* ---------- bi(...) assignments ---------- */

    for (i = 0; i < 11; i++) {
        bi[i][1] = -1.09852;
        bi[i][2] = 947.26074;
        bi[i][3] = 0.11604;
        bi[i][4] = 0.10637;
    }

    for (i = 11; i < 14; i++) {
        bi[i][1] = -1.09052;
        bi[i][2] = 1128.54873;
        bi[i][3] = 0.16117;
        bi[i][4] = 0.1664;
    }

    bi[14][1] = -0.88889;
    bi[14][2] = 1060.77063;
    bi[14][3] = 0.0058;
    bi[14][4] = 0.14347;

    for (i = 15; i < 18; i++) {
        bi[i][1] = -1.91121;
        bi[i][2] = 2112.54073;
        bi[i][3] = -0.03065;
        bi[i][4] = 0.1278;
    }

    bi[0][0] = -1.54021;
    bi[1][0] = -1.34794;
    bi[2][0] = -1.3758;
    bi[3][0] = -1.34399;
    bi[4][0] = -1.12207;
    bi[5][0] = -1.34495;
    bi[6][0] = -2.1112;
    bi[7][0] = -1.85606;
    bi[8][0] = -1.63406;
    bi[9][0] = -1.64298;
    bi[10][0] = -1.31222;
    bi[11][0] = -2.40229;
    bi[12][0] = -2.60203;
    bi[13][0] = -2.51562;
    bi[14][0] = -1.69223;
    bi[15][0] = -2.46339;
    bi[16][0] = -2.68192;
    bi[17][0] = -2.46894;
    bi[18][0] = -1.512;


    ci[0][0] = -0.8592;   ci[0][1] = -1.871;   ci[0][2] = 702.5264;  ci[0][3] = 0.32;  ci[0][4] = 0.13;
    ci[1][0] = -0.7779;   ci[1][1] = -1.255;   ci[1][2] = 645.7227;  ci[1][3] = 0.16;  ci[1][4] = 0.27;
    ci[2][0] = -0.9406;   ci[2][1] = -0.943;   ci[2][2] = 319.4692;  ci[2][3] = 0.26;  ci[2][4] = 0.29;
    ci[3][0] = -5.2447;   ci[3][1] = -0.8143;  ci[3][2] = 1775.0323; ci[3][3] = 0.22;  ci[3][4] = 0.27;
    ci[5][0] = -1.1879;   ci[5][1] = -1.9919;  ci[5][2] = 956.1444;  ci[5][3] = 0.32;  ci[5][4] = 0.14;
    ci[11][0] = -1.9439;  ci[11][1] = -1.4624; ci[11][2] = 1299.3775; ci[11][3] = 0.07; ci[11][4] = 0.27;
    ci[12][0] = -2.4594;  ci[12][1] = -1.5937; ci[12][2] = 1651.1601; ci[12][3] = 0.06; ci[12][4] = 0.33;
    ci[14][0] = -0.2386;  ci[14][1] = -1.87;   ci[14][2] = 980.1616;  ci[14][3] = 0.02; ci[14][4] = 0.2;
    ci[15][0] = -0.3512;  ci[15][1] = -1.8329; ci[15][2] = 1027.0341; ci[15][3] = 0.02; ci[15][4] = 0.2;
    ci[16][0] = -0.0514;  ci[16][1] = -1.1846; ci[16][2] = 479.2197;  ci[16][3] = 0.02; ci[16][4] = 0.2;
    ci[10][0] = -0.9734;  ci[10][1] = -1.871;  ci[10][2] = 702.5264;  ci[10][3] = 0.32; ci[10][4] = 0.13;

    ai[3][0] = -1.981;
    ai[3][1] = -1.552;
    ai[3][2] = 968.6;
    ai[3][3] = 0.102;
    ai[3][4] = 0.137;

    gi[3][0] = -1.4815;
    gi[3][1] = 0;
    gi[3][2] = 241.8075;
    gi[3][3] = 0.0908;
    gi[3][4] = 0.1843;

    celi[2][0] = -1.478;
    celi[2][1] = -1.845;
    celi[2][2] = 1218.91;

    celi[11][0] = 0.1994;
    celi[11][1] = -1.846;
    celi[11][2] = 602.42;

    celi[16][0] = -1.1569;
    celi[16][1] = -2.212;
    celi[16][2] = 1160.87;

    EaRInh[0] = 8064.9;   LnAInh[0] = 6.761;
    EaRInh[1] = 5524.5;   LnAInh[1] = 0.01131;
    EaRInh[2] = 8064.9;   LnAInh[2] = 6.761;
    EaRInh[3] = 5524.5;   LnAInh[3] = 0.01131;
    EaRInh[5] = 4872.3;   LnAInh[5] = -1.469;
    EaRInh[11] = 9936.7;  LnAInh[11] = 8.327;
    EaRInh[12] = 12120;   LnAInh[12] = 14.51;
    EaRInh[14] = 8064.9;  LnAInh[14] = 6.761;
    EaRInh[15] = 4872.3;  LnAInh[15] = -1.469;
    EaRInh[16] = 2822.5;  LnAInh[16] = -8.586;
    EaRInh[18] = 8064.9;  LnAInh[18] = 6.761;

}


void B2_ReadinAllData(SampleData* data)
{
    //

    // 声明并初始化所有需要的变量
    double mt = 0.0;
    simContext.Run10TestCases = 0, simContext.Loop10 = 0, simContext.Run_Seawater_Mixing = 0, simContext.LoopMixing = 0;
    simContext.Run_MixingTwoWells = 0, simContext.RunMultiMix = 0, simContext.LoopResChem = 0;
    UseH2Sgas = 0;
    //---------------------

    // D2_CalcDensitypH
    rhoOld = 0;
    rhoSSE = 0.00000002;
    Iteration = 0;

    double TDS = 0.0;
    simContext.yH2S = 0.0, simContext.yCO2 = 0.0;
    int Iteration2 = 0;
    //

    if (simContext.Read_InputII == 1) { simContext.nob = simContext.nob_Input; }
    if (simContext.Run1000Cases == 1) { simContext.nob = simContext.nob_Input; }

    for (int iRead = 0; iRead < simContext.nob; iRead++) {
        printf("%d", iRead);
        if (simContext.Run1000Cases == 1) {
            simContext.CaseCount[iRead] = simContext.LoopTP1000Cases;
        }
        simContext.j = simContext.CaseCount[iRead];
        simContext.kk = iRead;
        ReadInputPartA(simContext.kk, data);
    }

    for (int iRead = 0; iRead < simContext.nob; iRead++) {
        simContext.j = simContext.CaseCount[iRead];
        simContext.kk = iRead;
        ReadInputPartB(simContext.kk, data);
    }

    for (int iRead = 0; iRead < simContext.nob; iRead++) {
        simContext.j = simContext.CaseCount[iRead];
        simContext.kk = iRead;
        ReadInputPartC(simContext.kk, &mt, &Iteration2);
        if (simContext.RunStat == 0)
        {
            if (simContext.RunH2SGUI != 1)
            {
                if (simContext.Run_CalcConcFactor == 1)
                {
                    // Worksheets(MySheetMix).Cells(30, 8).Value = simContext.rho25c
                }
                else
                {
                    // Worksheets(mySheet).Cells(30, j + 2).Value = simContext.rho25c
                }
            }
        }
    }
    for (int iRead = 0; iRead < simContext.nob; iRead++)
    {
        simContext.j = simContext.CaseCount[iRead];
        simContext.kk = iRead;
        ReadInputPartD(simContext.kk, simContext.j, data);
    }

    if (simContext.nob > 0)
    {
        if (simContext.nob + simContext.nob_InputII > 1 && simContext.RunQualityControlChecks == 0)
        {
            simContext.j = 6;
        }
        simContext.TBH = data->T_initial;
        if (simContext.RunMultiMix == 1)
        {
            // simContext.TBH = Worksheets("MultiMix").Cells(2 + simContext.LoopResChem, 3).Value
        }
        if (simContext.RunWhatIf == 1)
        {
            simContext.TBH = simContext.TWIInit;
        }

        if (simContext.UseSI == 1)
        {
            simContext.TBH = simContext.TBH * 9 / 5 + 32;
        }

        if (simContext.TBH < 24.8 && simContext.Run1000Cases != 1)
        {
            // MsgBox("Initial temperature is below 24.8F or -4C. Initial temperature will be set to 77 F or 25 C.")
        }

        if (simContext.TBH < 24.8)
        {
            if (simContext.UseSI == 1)
            {
                // Worksheets(mySheet).Cells(39, j + 2) = 250;
            }
            else
            {
                // Worksheets(mySheet).Cells(39, j + 2) = 482;
            }
        }

        if (simContext.TBH < 24.8)
        {
            simContext.TBH = 77;
        }

        if (simContext.TBH > 482 && simContext.Run1000Cases != 1)
        {
            // MsgBox ("Initial temperature is above 482 F. Initail temperature will be set to 482 F or 250 C.")
        }
        if (simContext.TBH > 482)
        {
            if (simContext.UseSI == 1)
            {
                // Worksheets(mySheet).Cells(39, j + 2) = 250
            }
            else
            {
                // Worksheets(mySheet).Cells(39, j + 2) = 482
            }
        }
        if (simContext.TBH > 482)
        {
            simContext.TBH = 482;
        }

        simContext.TWH = data->T_final;
        if (simContext.RunMultiMix == 1)
        {
            // simContext.TWH = Worksheets("MultiMix").Cells(2 + simContext.LoopResChem, 3).Value
        }
        if (simContext.RunWhatIf == 1)
        {
            simContext.TWH = simContext.TWIInit;
        }
        if (simContext.UseSI == 1)
        {
            simContext.TWH = simContext.TWH * 9 / 5 + 32;
        }

        if (simContext.TWH < 24.8 && simContext.Run1000Cases != 1)
        {
            // MsgBox("Final temperature is below 24.8F or -4C. Final temperature will be set to 77 F or 25 C.")
        }

        if (simContext.TWH < 24.8)
        {
            if (simContext.UseSI == 1)
            {
                // Worksheets(mySheet).Cells(40, j + 2) = 25
            }
            else
            {
                // Worksheets(mySheet).Cells(40, j + 2) = 77
            }
        }
        if (simContext.TWH < 24.8)
        {
            simContext.TWH = 77;
        }
        if (simContext.TWH > 482 && simContext.Run1000Cases != 1)
        {
            // MsgBox ("Final temperature is above 482 F. Final temperature will be set to 482 F or 250 C.")
        }
        if (simContext.TWH > 482)
        {
            if (simContext.UseSI == 1)
            {
                // Worksheets(mySheet).Cells(40, j + 2) = 250
            }
            else
            {
                // Worksheets(mySheet).Cells(40, j + 2) = 482
            }
        }
        if (simContext.TWH > 482)
        {
            simContext.TWH = 482;
        }

        // simContext.PBH = Worksheets(mySheet).Cells(41, j + 2);
        simContext.PBH = data->P_initial;

        if (simContext.RunMultiMix == 1)
        {
            // simContext.PBH = Worksheets("MultiMix").Cells(2 + simContext.LoopResChem, 4).Value;
        }

        if (simContext.RunWhatIf == 1)
        {
            simContext.PBH = simContext.PWIInit;
        }

        if (simContext.UseSI == 1)
        {
            simContext.PBH = simContext.PBH * 14.503774; // note that simContext.TBH, simContext.TWH, simContext.PBH, simContext.PWH is defaulted to F, and psia
        }
        if (simContext.PBH < 12 && simContext.Run1000Cases != 1)
        {
            // MsgBox("Initial pressure is below 12 psia or 0.827 bar. Initial pressure will be set to 14.7 psi or 1.01325 bar.");
        }
        if (simContext.PBH < 12)
        {
            if (simContext.UseSI == 1)
            {
                // Worksheets(mySheet).Cells(41, j + 2) = 1.01325;
                simContext.PBH = 1.01325;
            }
            else
            {
                // Worksheets(mySheet).Cells(41, j + 2) = 14.7;
                simContext.PBH = 14.7;
            }
        }

        if (simContext.PBH < 12)
        {
            simContext.PBH = 14.7; // If P is zero, probably it is gauge and therefore add 1 atm
        }
        if (simContext.PBH > 30000 && simContext.Run1000Cases != 1)
        {
            // MsgBox("Initial pressure is above 30,000 psi. Initial pressure will be set to 30,000 psia or 2,068.4 bar.");
        }
        if (simContext.PBH > 30000)
        {
            if (simContext.UseSI == 1)
            {
                // Worksheets(mySheet).Cells(41, j + 2) = 30000 / 14.503774;
                simContext.PBH = 30000 / 14.503774;
            }
            else
            {
                // Worksheets(mySheet).Cells(41, j + 2) = 30000;
                simContext.PBH = 30000;
            }
        }

        if (simContext.PBH > 30000)
        {
            simContext.PBH = 30000;
        }

        // simContext.PWH = Worksheets(mySheet).Cells(42, j + 2);
        simContext.PWH = data->P_final;

        if (simContext.RunMultiMix == 1)
        {
            // simContext.PWH = Worksheets("MultiMix").Cells(2 + simContext.LoopResChem, 4).Value;
        }

        if (simContext.RunWhatIf == 1)
        {
            simContext.PWH = simContext.PWIInit;
        }

        if (simContext.UseSI == 1)
        {
            simContext.PWH = simContext.PWH * 14.503774;
        }

        if (simContext.PWH < 12 && simContext.Run1000Cases != 1)
        {
            // MsgBox("Final pressure is below 12 psia or 0.827 bar. Final pressure will be set to 14.7 psia or 1.01325 bar.");
        }

        if (simContext.PWH < 12)
        {
            if (simContext.UseSI == 1)
            {
                // Worksheets(mySheet).Cells(42, j + 2) = 1.01325;
                simContext.PWH = 1.01325;
            }
            else
            {
                // Worksheets(mySheet).Cells(42, j + 2) = 14.7;
                simContext.PWH = 14.7;
            }
        }

        if (simContext.PWH < 12)
        {
            simContext.PWH = 14.7; // If P is zero, probably it is gauge and therefore add 1 atm
        }

        if (simContext.PWH > 30000 && simContext.Run1000Cases != 1)
        {
            // MsgBox("Final pressure is above 30,000 psi. Final pressure will be set to 30000 psia or 2068.4 bar.");
        }

        if (simContext.PWH > 30000)
        {
            if (simContext.UseSI == 1)
            {
                // Worksheets(mySheet).Cells(42, j + 2) = 30000 / 14.503774;
                simContext.PWH = 30000 / 14.503774;
            }
            else
            {
                // Worksheets(mySheet).Cells(42, j + 2) = 30000;
                simContext.PWH = 30000;
            }
        }

        if (simContext.PWH > 30000)
        {
            simContext.PWH = 30000;
        }
    }

    if (simContext.tInh == 0)
    {
        simContext.tInh = 1;
    }
    // SelectInh = Worksheets("Input").Range("O4").Value
    // InhNo = Worksheets("Input").Range("O5").Value
    // InhNo1 = Worksheets("input").Range("o7").Value
    // FracInhNo1 = Worksheets("Input").Range("o8").Value / 100
    // InhNo2 = Worksheets("Input").Range("o9").Value
    simContext.InhNo2 = 18;

    for (int i = 0; i < 10; i++)
    {
        // InhName[i] = Worksheets("Input").Cells(1 + i, 18);
    }
    for (int i = 10; i < 20; i++)
    {
        // InhName(i) = Worksheets("Input").Cells(1 + i - 10, 20)
    }

    simContext.NoRiskcalc = 3;
    for (int i = 0; i < (simContext.NoRiskcalc - 1); i++)
    {
        simContext.ConcInhBarRisk[i] = simContext.MaxInh / (simContext.NoRiskcalc - 1) * i;
        simContext.ConcInhCalRisk[i] = simContext.MaxInh / (simContext.NoRiskcalc - 1) * i;
        simContext.ConcInhAnRisk[i] = simContext.MaxInh / (simContext.NoRiskcalc - 1) * i;
        simContext.ConcInhGypRisk[i] = simContext.MaxInh / (simContext.NoRiskcalc - 1) * i;

        simContext.ConcInhCelRisk[i] = simContext.MaxInh / (simContext.NoRiskcalc - 1) * i;
    };
}


void InitializeOptionClearCellContent()
{
    // 基本变量初始化（按 VBA 逻辑）
    UseSR = 0; simContext.UseTPCalciteSheet = 0;
    simContext.Read_InputII = 0;
    //int NCount_II = 0;  //NCount_II:InitializeO、CountNOB;   - 但在前者中仅赋值，在后者中赋值再使用，因此注释掉
    useEOS = 0;
    simContext.LoopMixing = 0; simContext.UseMolal = 0; //iTP = 0
    RunShellMultiflash = 0; H2Oevap = 0;
    simContext.Run_CalcConcFactor = 0;

    // 这些变量应为全局变量
    //int RunGoalSeek, RunStatGoalSeek; Run_MassTransfer;
    //RunH2SGUI, RunSMT, RunStat, RunH2SPartition;

    RunGoalSeek = (RunGoalSeek == 1 ? 1 : 0);
    RunStatGoalSeek = (RunStatGoalSeek == 1 ? 1 : 0);
    Run_MassTransfer = (Run_MassTransfer == 1 ? 1 : 0);

    // If Worksheets("Input").Range("I11") = "Saturation Ratio values" Then UseSR = 1;
    //没看懂，发现总是不会进入UseSR=1;先注释掉了

    usePTB = 1;   //目前来看一定为 1
    /*
    if (simContext.RunH2SGUI != 1 && RunSMT != 1 && simContext.RunStat != 1) {
         if Worksheets("Calcite").Range("F3") == "mg/L"
             usePTB = 0;
         else
             usePTB = 1;
    }
    */

    // ============== 以下所有 Excel 清空操作都保留为注释 ==============
    if (RunSMT != 1) {
        if (simContext.RunH2SGUI != 1) {

            // Worksheets("Input").Range("P14:T22").Value = Null;
            // Worksheets("input II").Range("C83:Cx91").Value = Null;
            // Worksheets("input II").Range("C80:Cx80").Value = Null;
            // Worksheets("input").Range("H10:H37").Value = Null;
            // Worksheets("input").Range("H43:H49").Value = Null;
            // Worksheets("input").Range("H62:H79").Value = Null;

            // Worksheets("Well#1+Water").Range("B5:AC15") = Null;
            // Worksheets("Well#1+Water").Range("B40:AG50") = Null;

            // Worksheets("Mixing Two Wells").Range("B5:AC15") = Null;
            // Worksheets("Mixing Two Wells").Range("B40:AG50") = Null;

            // Worksheets("Calcite").Range("D16:G25") = Null;
            // Worksheets("Calcite").Range("P20:P27") = Null;

            // Worksheets("Barite").Range("C16:E25") = Null;
            // Worksheets("Halite").Range("C16:E25") = Null;

            // Worksheets("Other SO4s").Range("C17:J26") = Null;
            // Worksheets("Sulfides,Fluorite,Carbonates").Range("C17:J26") = Null;

            // Worksheets("Use Mass Transfer").Range("A12:H21") = Null;

            // Deposition Prediction
            // Worksheets("Deposition Prediction").Range("C25:R34") = Null;
            // Worksheets("Deposition Prediction").Range("C41:R50") = Null;

            // Corrosion
            // Worksheets("Corrosion").Range("C25:N34") = Null;
            // Worksheets("Corrosion").Range("C41:N50") = Null;

            // Activity coefficients
            // Worksheets("Calcite").Range("Q4:AG13") = Null;
            // Worksheets("Barite").Range("Q4:Z13") = Null;
            // Worksheets("Other SO4s").Range("Q4:AA13") = Null;
            // Worksheets("Silicates").Range("X2:Z13") = Null;
        }
    }

    //double feed_Composition[15];这个仅在partD中使用，所以设为局部的了，此处忽略

    DpHj = 0;


    // ============== 第二段 Excel 清空（条件判断）==============
    if (RunSMT != 1) {
        if (RunH2SPartition != 1) {

            // For i = 1 to 17:
            for (int i = 0; i < 17; i++) {
                // Worksheets("Input").Cells(11+2*i, 9).Value = Null;
                // Worksheets("Input").Cells(11+2*i,10).Value = Null;
            }

            // For i = 1 to 7:
            for (int i = 0; i < 7; i++) {
                // Worksheets("Input").Cells(46+2*i,9) = Null;
                // Worksheets("Input").Cells(46+2*i,10) = Null;
            }

            // Worksheets("Multiple Ppt").Unprotect("eesi");
            // Worksheets("Multiple Ppt").Range("A4:P100") = Null;
            // Worksheets("Multiple Ppt").Protect("eesi");
        }
    }
}


void CountNOB()
{

    //nob, Ncount, NCount_II, nob_Input, nob_InputII;
    //extern int RunStatReservoirCalc, RunStatGoalSeek, RunStatSICalcSSP;
    //extern int RunStatMix, Run_Seawater_Mixing, RunMultiMix, RunMultiMixSlb;
    //extern int Read_InputII, Run1000Cases;

    // 初始化
    //simContext.nob = 0;
    Ncount = 0;
    NCount_II = 0;
    //simContext.nob_Input = 0;
    //simContext.nob_InputII = 0;
    simContext.Run1000Cases = 0;
    RunStatReservoirCalc = 0;//没见过给这个值赋值？？？

    if (RunStatReservoirCalc == 1)
    {
        simContext.nob = 2;
        simContext.CaseCount[0] = 1;   // CaseCount(1)
        simContext.CaseCount[1] = 9;   // CaseCount(2)
    }
    else if (RunStatGoalSeek == 1)
    {
        simContext.nob = 2;
        /*
        If Worksheets("Halite module").Range("H18").Value = "Use produced water composition" Then
            CaseCount(1) = 1: CaseCount(2) = 9:
            Worksheets("Output data sheet").Range("Q8") = "Produced water"
        End If

        If Worksheets("Halite module").Range("H18").Value = "Use calculated reservoir water composition" Then
            CaseCount(1) = 5: CaseCount(2) = 9:
            Worksheets("Output data sheet").Range("Q8") = "Calculated reservoir water"
        End If
        */
    }

    // ====================== RunStatSICalcSSP = 1 (Halite 1 water) ======================
    else if (RunStatSICalcSSP == 1)
    {
        /*
        If Worksheets("Halite module").Range("G33").Value = "Use produced water composition" Then
            nob = 1: CaseCount(1) = 1
            Worksheets("Halite analysis").Range("C3") = "Produced water"
            Worksheets("Halite analysis").Range("b18") = "Figure 1. Plots of halite SI and pressure versus temperature based on calculated reservoir water composition."
        End If
            If Worksheets("Halite module").Range("G33").Value = "Use calculated reservoir water composition" Then
            nob = 1: CaseCount(1) = 5
            Worksheets("Halite analysis").Range("c3") = "Calculated reservoir water"
            Worksheets("Halite analysis").Range("b18") = "Figure 1. Plots of halite SI and pressure versus temperature based on calculated reservoir water composition."
        End If
        If Worksheets("Halite module").Range("G33").Value = "Use fresh water composition" Then
            nob = 1: CaseCount(1) = 9
            Worksheets("Halite analysis").Range("c3") = "Fresh water"
            Worksheets("Halite analysis").Range("b18") = "Figure 1. Plots of halite SI and pressure versus temperature based on calculated reservoir water composition."
        End If
        */
    }
    else if (RunStatSICalcSSP == 2)
    {
        simContext.nob = 2;
        simContext.CaseCount[1] = 9;  // CaseCount(2)

        /*
        If Worksheets("Halite module").Range("I33").Value = "Use produced water composition" Then
            CaseCount(1) = 1
            Worksheets("Halite analysis").Range("K3") = "Produced water and fresh water"
            Worksheets("Halite analysis").Range("J18") = "Figure 2. Plots of halite SI/SR and pressure versus temperature based on mixing produced water with fresh water.
        End If
        If Worksheets("Halite module").Range("I33").Value = "Use calculated reservoir water composition" Then
            CaseCount(1) = 5
            Worksheets("Halite analysis").Range("K3") = "Calculated reservoir water and fresh water"
            Worksheets("Halite analysis").Range("J18") = "Plots of halite SI/SR and pressure versus temperature based on mixing produced water with fresh water."
        End If
        */

    }

    // ====================== RunStatSICalcSSP = 3 ======================
    else if (RunStatSICalcSSP == 3)
    {
        /*
        If Worksheets("Halite module").Range("H18").Value = "Use produced water composition" Then
            nob = 1: CaseCount(1) = 1
            Worksheets("Output data sheet").Range("J3") = "Produced water"
        End If
        If Worksheets("Halite module").Range("H18").Value = "Use calculated reservoir water composition" Then
            nob = 1: CaseCount(1) = 5
            Worksheets("Output data sheet").Range("J3") = "Calculated reservoir water"
        End If
        */
    }
    else if (simContext.RunStatMix == 1)
    {
        simContext.nob = 2;
        simContext.CaseCount[1] = 9;
        /*
        If Worksheets("Halite module").Range("I33").Value = "Use produced water composition" Then
            CaseCount(1) = 1
            Worksheets("Halite analysis").Range("T3") = "Produced water and fresh water"
            Worksheets("Halite analysis").Range("R18") = "Figure 3. Plot of halite SI/SR versus mixing ratio at bottom and surface T and P based on mixing produced water with fresh water."
        End If
        If Worksheets("Halite module").Range("I33").Value = "Use calculated reservoir water composition" Then
            CaseCount(1) = 5
            Worksheets("Halite analysis").Range("T3") = "Calculated reservoir water and fresh water"
            Worksheets("Halite analysis").Range("R18") = "Figure 3. Plot of halite SI/SR versus mixing ratio at bottom and surface T and P based on mixing calculated reservoir water with fresh water."
        End If
        */
    }
    // ====================== 最后 else: 读取 Input Sheet 的勾选框 ======================
    else
    {
        //阅读输入表顶部的复选框，以确定要使用哪些井。
        for (int i = 0; i < 5; i++)
        {

            // if Worksheets("Input").Cells(3, 2+i+1) = True
                //nob++;
                //CaseCount[i - Ncount] = i;
            //else   
                //Ncount++;
        }

        //nob_Input = nob;

        //阅读输入表顶部的复选框，以确定要使用哪些井。
        for (int i = 0; i < 100; i++)
        {
            /*
            If Worksheets("Input II").Cells(3, 2 + i).Value = True Then
                nob = nob + 1
                CaseCount_II(i - NCount_II) = i
            Else
                NCount_II = NCount_II + 1
            End If
            */

            //nob_InputII = nob - nob_Input;
        }

        //if (nob_InputII > 0)
        //    Read_InputII = 1;
        //if (nob == 0) {
        //    printf("No box on Row 3 of the Input sheet is selected.\n");
        //    exit(0);
        //}
    }

    // ====================== Seawater mixing ======================
    if (simContext.Run_Seawater_Mixing == 1)
    {
        simContext.CaseCount[0] = 1;
        simContext.CaseCount[1] = 2;
        simContext.nob = 2;
    }

    // ====================== Multi-mix ======================
    if (simContext.RunMultiMix == 1 || RunMultiMixSlb == 1)
    {

        Ncount = 0;
        simContext.nob = 0;

        for (int i = 0; i < 5; i++) {
            /*
            If Worksheets("Input").Cells(3, 2 + i).Value = True Then
                nob = nob + 1
                CaseCount(i - Ncount) = i
            Else
                Ncount = Ncount + 1
            End If
            */
        }

        //nob_InputII = 0;
        //nob_Input = nob;

        //if (nob == 0) {
        //    printf("You must select Brines from the Input Sheet to run.\n");
        //    exit(0);
        //}
    }
}


//B3中会修改data的API，所以这里直接引用传递
void B3_CalcConcs(double& API)
{
    /*  局部变量声明 */
    int i, c, a, n, iNG, iden;
    double molAlkF = 0, molTACF = 0, molTNH4F = 0;
    double molTH3BO3F = 0, molTH2SaqF = 0, molTH4SiO4F = 0;
    double molTFeF;

    xMeOH = 0;
    xMEG = 0;
    Alk = 0;
    TAc = 0;
    TH2Saq = 0;
    TDS = 0;
    TH4SiO4 = 0;
    TNH4 = 0;
    TH3BO3 = 0;
    TFe = 0;

    nTCO2 = 0;
    nTCH4 = 0;
    nTH2S = 0;
    mass_w = 0;
    Mass_o = 0;
    mass_MeOH = 0;
    mass_MEG = 0;

    API = 0;
    SGG = 0;
    QTotal = 0;
    TC = 0;
    VgTP = 0;
    VO = 0;

    total_moles = 0;
    VW = 0;
    VO = 0;
    VgTP = 0;
    SGG = 0;

    nTCO2EOS = 0;
    nTH2sEOS = 0;

    /* ---- Reset arrays ---- */
    //下面这三个数组仅在B3中使用
    double molcF[15]; double molaF[15]; double molnF[10];
    for (c = 0; c < NumCat; c++)
    {
        mc[c] = 0;
        molcF[c] = 0;
    }
    for (a = 0; a < NumAn; a++)
    {
        ma[a] = 0;
        molaF[a] = 0;
    }
    for (n = 0; n < NumNeut; n++)
    {
        mn[n] = 0;
        molnF[n] = 0;
    }
    for (i = 0; i < 15; i++)
        z[i] = 0;

    /* ---- Sum basic volumes and masses ---- */
    for (i = 0; i < simContext.nob; i++)
    {
        VW += simContext.VwMix[i];
        VgTP += simContext.VgTPMix[i];
        VO += simContext.VoMix[i];
        mass_w += simContext.mass_w_Mix[i];
        Mass_o += simContext.mass_o_Mix[i];
    }

    mass_w_0 = mass_w;

    /* ---- Determine mixing fractions ---- */
    if (simContext.Run_Seawater_Mixing == 0 && simContext.Run_MixingTwoWells == 0 &&
        simContext.RunStatMix == 0 && RunShellMultiflash != 1)
    {

        for (i = 0; i < simContext.nob; i++)
        {
            simContext.MixFrac[i] = simContext.mass_w_Mix[i] / mass_w;
            simContext.MixFracOil[i] = simContext.VoMix[i] / VO;
            simContext.MixFracGas[i] = simContext.VgTPMix[i] / VgTP;
        }
    }

    if (simContext.Run_Seawater_Mixing == 1 || simContext.RunStatMix == 1)
    {
        simContext.MixFracOil[0] = 1; simContext.MixFracOil[1] = 0;
        simContext.MixFracGas[0] = 1; simContext.MixFracGas[1] = 0;

        simContext.MixFrac[0] = simContext.mass_w_Mix[0] / mass_w;
        simContext.MixFrac[1] = 1 - simContext.MixFrac[0];
    }

    double MixFracTwoWells[11] = { 0 };
    if (simContext.Run_MixingTwoWells == 1)
    {
        for (i = 0; i < simContext.nob; i++)
        {
            MixFracTwoWells[0] = simContext.mass_w_Mix[0] / mass_w;
            MixFracTwoWells[1] = 1 - MixFracTwoWells[0];
            simContext.MixFracOil[i] = simContext.MixFrac[i];
            simContext.MixFracGas[i] = simContext.MixFrac[i];
        }
    }

    /* ---- Component mixing ---- */
    for (i = 0; i < simContext.nob; i++)
    {
        /*
        注意vb是判断Run_MixingTwoWells是否为1
        若为1，这些数组都会乘一个MixFracTwoWells[i]

        否则乘的是MixFrac[i],所以这里先用一个MF来保存要乘的值，避免重复索引耗时
        */
        double MF = (simContext.Run_MixingTwoWells == 1) ? MixFracTwoWells[i] : simContext.MixFrac[i];

        /* cations */
        mc[iH] += MF * simContext.HstpMix[i];
        mc[iNa] += MF * simContext.NaMix[i];
        mc[iK] += MF * simContext.KMix[i];
        mc[iMg] += MF * simContext.MgMix[i];
        mc[iCa] += MF * simContext.CaMix[i];
        mc[iSr] += MF * simContext.SrMix[i];
        mc[iBa] += MF * simContext.BaMix[i];
        mc[iFe] += MF * simContext.FeMix[i];
        mc[iZn] += MF * simContext.ZnMix[i];
        mc[iPb] += MF * simContext.PbMix[i];
        mc[iNH4] += MF * simContext.NH4STPMix[i];

        //不为1时才会赋值这个iRa
        if (simContext.Run_MixingTwoWells != 1)
            mc[iRa] += MF * simContext.RaMix[i];

        /* anions */
        ma[iOH] += MF * simContext.OHstpMix[i];
        ma[iCl] += MF * simContext.ClMix[i];
        ma[iAc] += MF * simContext.ACstpMix[i];
        ma[iH2BO3] += MF * simContext.H2BO3stpMix[i];
        ma[iHCO3] += MF * simContext.HCO3stpMix[i];
        ma[iCO3] += MF * simContext.CO3stpMix[i];
        ma[iSO4] += MF * simContext.SO4Mix[i];
        ma[iBr] += MF * simContext.BrMix[i];
        ma[intF] += MF * simContext.FMix[i];

        /* neutral species */
        TH4SiO4 += MF * simContext.TH4SiO4Mix[i];
        Alk += MF * simContext.AlkMix[i];
        TAc += MF * simContext.TAcMix[i];
        TNH4 += MF * simContext.TNH4Mix[i];
        TH3BO3 += MF * simContext.TH3BO3Mix[i];

        /* cosolvents */
        mass_MeOH += simContext.mass_MeOH_mix[i];
        mass_MEG += simContext.mass_MEG_mix[i];

        /* densities */
        API += simContext.OilDensityMix[i] * simContext.MixFracOil[i];
        SGG += simContext.GasDensityMix[i] * simContext.MixFracGas[i];

        /* gas */
        for (iNG = 0; iNG < 15; iNG++)
            z[iNG] += simContext.zMix[i][iNG] * simContext.Total_molesMix[i];

        total_moles += simContext.Total_molesMix[i];
        nTCO2 += simContext.nTCO2Mix[i];
        nTCH4 += simContext.nTCH4Mix[i];
        nTH2S += simContext.nTH2SMix[i];

        nTCO2EOS += simContext.nTCO2MixEOS[i];
        nTH2sEOS += simContext.nTH2SMixEOS[i];
    }

    TFe = mc[iFe];

    /* ---- EOS logic ---- */
    if (useEOS > 0)
    {
        nTCO2 = nTCO2EOS;
        nTH2S = nTH2sEOS;

        if (total_moles > 0)
        {
            for (iNG = 0; iNG < 15; iNG++)
            {
                z[iNG] /= total_moles;
                if (z[iNG] < 1e-7) z[iNG] = 0;
                z_before_precipitation[iNG] = z[iNG];
            }

            double zHC = z[0];
            for (iNG = 3; iNG < 14; iNG++)
                zHC += z[iNG];

            if (zHC == 0)
            {
                useEOS = 0;
                if (simContext.Run_Seawater_Mixing != 1 &&
                    simContext.Run_MixingTwoWells != 1 && simContext.Run10TestCases != 1)
                    simContext.errmsg[13] = 14;
            }

            SumofZ = 0;
            for (iNG = 0; iNG < 15; iNG++)
                SumofZ += z[iNG];
        }

        if ((simContext.Run_Seawater_Mixing == 1 ||
            simContext.Run_MixingTwoWells == 1) &&
            simContext.LoopMixing == 1)
            useEOS = 0;

        /* mixing moles */
        if (simContext.Run_MixingTwoWells == 1)
        {
            for (i = 0; i < simContext.nob; i++)
            {
                double mixFracTwoWells_i = simContext.MixFrac[i];

                for (c = 0; c < NumCat; c++)
                    molcF[c] += molc[c][i] * mixFracTwoWells_i;

                for (a = 0; a < NumAn; a++)
                    molaF[a] += mola[a][i] * mixFracTwoWells_i;

                for (n = 0; n < NumNeut; n++)
                    molnF[n] += moln[n][i] * mixFracTwoWells_i;

                molAlkF += simContext.molAlk[i] * mixFracTwoWells_i;
                molTACF += simContext.molTAC[i] * mixFracTwoWells_i;
                molTNH4F += simContext.molTNH4[i] * mixFracTwoWells_i;
                molTH3BO3F += simContext.molTH3BO3[i] * mixFracTwoWells_i;
                molTH2SaqF += simContext.molTH2Saq[i] * mixFracTwoWells_i;
                molTH4SiO4F += simContext.molTH4SiO4[i] * mixFracTwoWells_i;
            }

        }
        else {

            for (i = 0; i < simContext.nob; i++)
            {
                double mixFrac_i = simContext.MixFrac[i];

                for (c = 0; c < NumCat; c++)
                    molcF[c] += molc[c][i] * mixFrac_i;

                for (a = 0; a < NumAn; a++)
                    molaF[a] += mola[a][i] * mixFrac_i;

                for (n = 0; n < NumNeut; n++)
                    molnF[n] += moln[n][i] * mixFrac_i;



                molAlkF += simContext.molAlk[i] * mixFrac_i;
                molTACF += simContext.molTAC[i] * mixFrac_i;
                molTNH4F += simContext.molTNH4[i] * mixFrac_i;
                molTH3BO3F += simContext.molTH3BO3[i] * mixFrac_i;
                molTH2SaqF += simContext.molTH2Saq[i] * mixFrac_i;
                molTH4SiO4F += simContext.molTH4SiO4[i] * mixFrac_i;
            }
        }

        molTFeF = molcF[iFe];
    }

    /* ---- cosolvent mole fractions ---- */
    if (mass_MeOH > 0)
        xMeOH = (mass_MeOH / 32.0) / ((mass_MeOH / 32.0) + (mass_w / 18.0));

    if (mass_MEG > 0)
        xMEG = (mass_MEG / 62.07) / ((mass_MEG / 62.07) + (mass_w / 18.0));

    CalcIonicStrength();

    if (mass_MeOH > 0)
        IStCosolvent = ISt * mass_w / (mass_w + mass_MeOH);
    if (mass_MEG > 0)
        IStCosolvent = ISt * mass_w / (mass_w + mass_MEG);

    /* --------------------------------------- */
    /* =========== Non-EOS branch ============= */
    /* --------------------------------------- */

    if (useEOS == 0)
    {
        if (simContext.useTPVol == 1)
            fTPFunc(1);
        else
            fTPFunc(0);

        SGG = SGG / (Patm * 28.97 / (0.08206 * TK));
        API = 141.5 / (API / (fH2ODensity(TK, PBar) / 1000.0)) - 131.5;

        C1_ThermodynamicEquilConsts();
        C2_PitzerActCoefs_T_P_ISt(gNeut, aH2O, TK, TC, PBar, Patm);

        nTCO2_before_precipitation = nTCO2;
        nTH2S_before_precipitation = nTH2S;

        PengRobinson3();
        RatioOilBPoints = fRatioOilBPoints(API);

        //Call C4_SSPEquilCalcs(0, 5, 2, KspCalcite)
        C4_SSPEquilCalcs(0, 4, 2, KspCalcite);

        fmn();

        CalcIonicStrength();
        C2_PitzerActCoefs_T_P_ISt(gNeut, aH2O, TK, TC, PBar, Patm);
        ////Call C4_SSPEquilCalcs(0, 5, 2, KspCalcite)
        C4_SSPEquilCalcs(0, 4, 2, KspCalcite);

        rhoTP = CalcRhoTP(TK, TC, PBar, Patm);

        yCH4 = PCH4 / Ppsia;
        simContext.yCO2 = PCO2 / Ppsia;
        simContext.yH2S = PH2S / Ppsia;

        /* ---- calculate TDS at STP ---- */
        if (TC == 25)
            simContext.rho25c = rhoTP;
        else
        {

            fTPFunc(0);
            C1_ThermodynamicEquilConsts();
            C2_PitzerActCoefs_T_P_ISt(gNeut, aH2O, TK, TC, PBar, Patm);

            nTCO2_before_precipitation = nTCO2;
            nTH2S_before_precipitation = nTH2S;

            PengRobinson3();
            RatioOilBPoints = fRatioOilBPoints(API);
            //Call C4_SSPEquilCalcs(0, 5, 2, KspCalcite)
            C4_SSPEquilCalcs(0, 4, 2, KspCalcite);

            fmn();
            CalcIonicStrength();
            C2_PitzerActCoefs_T_P_ISt(gNeut, aH2O, TK, TC, PBar, Patm);
            //Call C4_SSPEquilCalcs(0, 5, 2, KspCalcite)
            C4_SSPEquilCalcs(0, 4, 2, KspCalcite);

            simContext.rho25c = CalcRhoTP(TK, TC, PBar, Patm);
        }
    }

    /* --------------------------------------- */
    /* ============== EOS branch ============== */
    /* --------------------------------------- */

    if (useEOS != 0) {
        if (simContext.useTPVol == 1)
            fTPFunc(1);
        else
            fTPFunc(0);

        nTCO2_before_precipitation = nTCO2EOS;
        nTH2S_before_precipitation = nTH2sEOS;
        Total_moles_before_precipitation = total_moles;

        bool temp_ParametersWereRead = false;
        DoubleVec temp_gNeut(2);
        temp_gNeut[0] = gNeut[iCO2aq];
        temp_gNeut[1] = gNeut[iH2Saq];
        MultiPhaseFlash(temp_ParametersWereRead,
            mf_TCr, mf_PCr, mf_Omega, mf_MWgas, mf_kPr, mf_c0, mf_c1,
            TK, PBar, total_moles, z,
            temp_gNeut,
            aH2O, density, compositions, phi, Compr, beta,
            zOutput, mass_phase, MW_Phase, &No_Phases);

        /* compute water mass */
        mass_w = total_moles * beta[2] * compositions[14][3] * 0.01801528;

        if (mass_w == 0) {
            simContext.errmsg[13] = 14;
            useEOS = 0;
            goto L4000;
        }

        ////Call C4_SSPEquilCalcs(0, 5, 2, KspCalcite)
        C4_SSPEquilCalcs(0, 4, 2, KspCalcite);
        fmn();

        C2_PitzerActCoefs_T_P_ISt(gNeut, aH2O, TK, TC, PBar, Patm);
        //Call C4_SSPEquilCalcs(0, 5, 2, KspCalcite)
        C4_SSPEquilCalcs(0, 4, 2, KspCalcite);

        /* Recompute API + SGG */
        fTPFunc(3);

        //bool temp_ParametersWereRead = false;
        //DoubleVec temp_gNeut(2);
        temp_gNeut[0] = gNeut[iCO2aq];
        temp_gNeut[1] = gNeut[iH2Saq];
        MultiPhaseFlash(temp_ParametersWereRead,
            mf_TCr, mf_PCr, mf_Omega, mf_MWgas, mf_kPr, mf_c0, mf_c1,
            TK, PBar, total_moles, z,
            temp_gNeut,
            aH2O, density, compositions, phi, Compr, beta,
            zOutput, mass_phase, MW_Phase, &No_Phases);

        if (beta[0] > 0 && density[0] < 0.3)
            SGG = density[0] * 1000.0 / (Patm * 28.97 / (0.08206 * TK));
        else
            SGG = 0.6;

        if (beta[1] > 0 && density[1] > 0.3)
            API = 141.5 / (density[1] / (fH2ODensity(TK, PBar) / 1000.0)) - 131.5;
        else
            API = 30;
    }

L4000:
    TDS = 0;
    CalculateTDSDen = 0;

    for (iden = 1; iden < NumCat; iden++)
        CalculateTDSDen += 0.001 * mc[iden] * MWCat[iden];

    for (iden = 1; iden < NumAn; iden++)
        CalculateTDSDen += 0.001 * ma[iden] * MWAn[iden];

    for (iden = 1; iden < NumNeut; iden++)
        CalculateTDSDen += 0.001 * mn[iden] * MWNeut[iden];

    TDS = CalculateTDSDen / (1 + CalculateTDSDen) * simContext.rho25c * 1000000.0;

    for (c = 0; c < NumCat; c++)
        mcInit[c] = mc[c];

    for (a = 0; a < NumAn; a++)
        maInit[a] = ma[a];

    AlkInit = Alk;
    TH4SiO4Init = TH4SiO4;
    TNH4Init = TNH4;
    TH3BO3Init = TH3BO3;
    TAcInit = TAc;
    TFeInit = TFe;

    fmn();
}


void B4_CalcFinalBrine()
{

    if (simContext.nob == 1 && simContext.Run_CalcConcFactor == 0) {

        for (int i = 8; i < 36; i++) {
            // Worksheets("Input").Cells(i + 1, 8).Value = "";
        }
        for (int i = 41; i < 48; i++) {
            // Worksheets("Input").Cells(i + 1, 8).Value = "";
        }
        for (int i = 51; i < 53; i++) {
            // Worksheets("Input").Cells(i + 1, 8).Value = "";
        }
        for (int i = 56; i < 78; i++) {
            // Worksheets("Input").Cells(i + 1, 8).Value = "";
        }

    }
    else
    {

        if (simContext.Run_Seawater_Mixing == 1 || simContext.Run_MixingTwoWells == 1)
            goto label100;

        if (simContext.UseMolal == 0)
        {

            //Worksheets("Input").Range("h10").Value = mc[iNa] * (22990 * (simContext.rho25c - TDS / 1000000.0));

            //Worksheets("Input").Range("h11").Value = mc[iK] * (39102 * (simContext.rho25c - TDS / 1000000.0));

            //Worksheets("Input").Range("h12").Value = mc[iMg] * (24305 * (simContext.rho25c - TDS / 1000000.0));
            //Worksheets("Input").Range("h13").Value = mc[iCa] * (40080 * (simContext.rho25c - TDS / 1000000.0));
            //Worksheets("Input").Range("h14").Value = mc[iSr] * (87620 * (simContext.rho25c - TDS / 1000000.0));
            //Worksheets("Input").Range("h15").Value = mc[iBa] * (137330 * (simContext.rho25c - TDS / 1000000.0));
            //Worksheets("Input").Range("h16").Value = TFe* (55847 * (simContext.rho25c - TDS / 1000000.0));
            //Worksheets("Input").Range("h17").Value = mc[iZn] * (65380 * (simContext.rho25c - TDS / 1000000.0));
            //Worksheets("Input").Range("h18").Value = mc[iPb] * (207200 * (simContext.rho25c - TDS / 1000000.0));

            //Worksheets("Input").Range("h19").Value = ma[iCl] * (35450 * (simContext.rho25c - TDS / 1000000.0));
            //Worksheets("Input").Range("h20").Value = ma[iSO4] * (96064 * (simContext.rho25c - TDS / 1000000.0));
            //Worksheets("Input").Range("h21").Value = ma[intF] * (18998 * (simContext.rho25c - TDS / 1000000.0));
            //Worksheets("Input").Range("h22").Value = ma[iBr] * (79904 * (simContext.rho25c - TDS / 1000000.0));
            //
            //Worksheets("Input").Range("h23").Value = TH4SiO4* (28085 * (simContext.rho25c - TDS / 1000000.0));
            //Worksheets("Input").Range("h24").Value = Alk* (61019 * (simContext.rho25c - TDS / 1000000.0));
            //Worksheets("Input").Range("h26").Value = TAc* (60054 * (simContext.rho25c - TDS / 1000000.0));
            //Worksheets("Input").Range("h27").Value = TNH4* (17031 * (simContext.rho25c - TDS / 1000000.0));
            //Worksheets("Input").Range("h28").Value = TH3BO3* (10811 * (simContext.rho25c - TDS / 1000000.0));
            //
            //Worksheets("Input").Range("h62").Value = mc[iRa] * 1e12 * (226 * (simContext.rho25c - TDS / 1000000.0));

        }
        else
        {
            //Worksheets("Input").Range("h10").Value = mc[iNa];

            //Worksheets("Input").Range("h11").Value = mc[iK];
            //Worksheets("Input").Range("h12").Value = mc[iMg];
            //Worksheets("Input").Range("h13").Value = mc[iCa];
            //Worksheets("Input").Range("h14").Value = mc[iSr];
            //Worksheets("Input").Range("h15").Value = mc[iBa];
            //Worksheets("Input").Range("h16").Value = TFe;
            //Worksheets("Input").Range("h17").Value = mc[iZn];
            //Worksheets("Input").Range("h18").Value = mc[iPb];

            //Worksheets("Input").Range("h19").Value = ma[iCl];
            //Worksheets("Input").Range("h20").Value = ma[iSO4];
            //Worksheets("Input").Range("h21").Value = ma[intF];
            //Worksheets("Input").Range("h22").Value = ma[iBr];
            //Worksheets("Input").Range("h23").Value = TH4SiO4;
            //Worksheets("Input").Range("h24").Value = Alk;
            //Worksheets("Input").Range("h26").Value = TAc;
            //Worksheets("Input").Range("h27").Value = TNH4;
            //Worksheets("Input").Range("h28").Value = TH3BO3;
        }

        //Worksheets("Input").Range("h29").Value = TDS;

        //Worksheets("Input").Range("h30").Value = rho25c;

        //Worksheets("Input").Range("h31").Value = yCO2 * 100;

        //Worksheets("Input").Range("h33").Value = (H2Saq + HS)* (34080 * (rho25c - TDS / 1000000.0));

        //Worksheets("Input").Range("h34").Value = pH - DpHj;

        //Worksheets("Input").Range("h35").Value = VgTP / 28.31685;
        if (simContext.UseSI == 1)
            //Worksheets("Input").Range("h35").Value =  = VgTP / 1000.0;

        //Worksheets("Input").Range("h36").Value = VO;
            if (simContext.UseSI == 1) //Worksheets("Input").Range("h36").Value  = VO * 0.159;

                //Worksheets("Input").Range("h37").Value = VW;
                if (simContext.UseSI == 1) //Worksheets("Input").Range("h37").Value  = VW * 0.159;

                    //Worksheets("Input").Range("h45").Value =  mass_MeOH / 0.7914 / 159.0;
                    if (simContext.UseSI == 1) //Worksheets("Input").Range("h45").Value = tmp45 * 0.159;

                        //Worksheets("Input").Range("h46").Value = mass_MEG / 1.1135 / 159.0;
                        if (simContext.UseSI == 1)
                            //Worksheets("Input").Range("h45").Value = mass_MeOH / 1.1135 / 159 * 0.159
                            double aaa123 = 0;// 占位语句，确认后请删除
    }

label100:
    return;
    // 跳转目标
}


void OutputActivityCoefficients()
{
    /*
    Worksheets("Calcite").Range("x1:ax13").Value = ""
    Worksheets("Barite").Range("w1:z13").Value = ""
    Worksheets("Halite").Range("O1:R13").Value = ""
    Worksheets("Other SO4s").Range("x1:aa13").Value = ""
    Worksheets("Sulfides,Fluorite,Carbonates").Range("z1:ap13").Value = ""
    */

    if (simContext.OutPutActCoefs == 1)
    {
        /*
            Worksheets("Calcite").Range("x2") = "Activity coefficients"
            Worksheets("Calcite").Range("z3") = "aH2O"
            Worksheets("Calcite").Range("aa3") = "gH"
            Worksheets("Calcite").Range("ab3") = "gOH"
            Worksheets("Calcite").Range("ac3") = "gCa"
            Worksheets("Calcite").Range("ad3") = "gMg"
            Worksheets("Calcite").Range("ae3") = "gHCO3"
            Worksheets("Calcite").Range("af3") = "gCO3"
            Worksheets("Calcite").Range("ag3") = "gAc"


            Worksheets("Barite").Range("w2") = "Activity coefficients"

            Worksheets("Barite").Range("y3") = "gBa"
            Worksheets("Barite").Range("z3") = "gSO4"

            Worksheets("Halite").Range("O2") = "Activity coefficients"

            Worksheets("Halite").Range("Q3") = "gNa"
            Worksheets("Halite").Range("R3") = "gCl"

            Worksheets("Other SO4s").Range("x2") = "Activity coefficients"

            Worksheets("Other SO4s").Range("z3") = "gSr"
            Worksheets("Other SO4s").Range("aa3") = "gSO4"

            Worksheets("Sulfides,Fluorite,Carbonates").Range("z2") = "Activity coefficients"

            Worksheets("Sulfides,Fluorite,Carbonates").Range("ab3") = "gFe"
            Worksheets("Sulfides,Fluorite,Carbonates").Range("ac3") = "gZn"
            Worksheets("Sulfides,Fluorite,Carbonates").Range("ad3") = "gF"
            Worksheets("Sulfides,Fluorite,Carbonates").Range("ae3") = "gHS"
            Worksheets("Sulfides,Fluorite,Carbonates").Range("af3") = "gPb"

            Worksheets("Silicates").Range("x2") = "Activity coefficients" 'Dai added 2020
            Worksheets("Silicates").Range("z3") = "gSilicates" 'Dai added 2020
        */
    }
    if (simContext.UseSI == 0)
    {
        /*
            Worksheets("Calcite").Range("x3") = "T(F)"
            Worksheets("Calcite").Range("y3") = "P(psia)"
            Worksheets("Barite").Range("w3") = "T(F)"
            Worksheets("Barite").Range("x3") = "P(psia)"
            Worksheets("Halite").Range("O3") = "T(F)"
            Worksheets("Halite").Range("P3") = "P(psia)"
            Worksheets("Other SO4s").Range("x3") = "T(F)"
            Worksheets("Other SO4s").Range("y3") = "P(psia)"
            Worksheets("Sulfides,Fluorite,Carbonates").Range("z3") = "T(F)"
            Worksheets("Sulfides,Fluorite,Carbonates").Range("aa3") = "P(psia)"
            Worksheets("Silicates").Range("x3") = "T(F)" 'Dai added 2020
            Worksheets("Silicates").Range("y3") = "P(psia)" 'Dai added 2020
        */
    }
    else
    {
        /*
            Worksheets("Calcite").Range("Q3") = "T(C)"
            Worksheets("Calcite").Range("R3") = "P(Bar)"
            Worksheets("Barite").Range("O3") = "T(C)"
            Worksheets("Barite").Range("P3") = "P(Bar)"
            Worksheets("Halite").Range("O3") = "T(C)"
            Worksheets("Halite").Range("P3") = "P(Bar)"
            Worksheets("Other SO4s").Range("x3") = "T(C)"
            Worksheets("Other SO4s").Range("y3") = "P(Bar)"
            Worksheets("Sulfides,Fluorite,Carbonates").Range("z3") = "T(C)"
            Worksheets("Sulfides,Fluorite,Carbonates").Range("aa3") = "P(Bar)"
            Worksheets("Silicates").Range("x3") = "T(C)" 'Dai added 2020
            Worksheets("Silicates").Range("y3") = "P(Bar)" 'Dai added 2020
        */

    }

}


void B5_CalculateSIvalues(double& API)
{
    if (LoopTP == 10) {
        pH = pH;
    }

    // --- Initialize SI values (use 0 to represent Null)
    SICal = 0; SIDol = 0; SISid = 0; SIBar = 0; SIGyp = 0;
    SIHemi = 0; SIAn = 0; SICel = 0; SIHal = 0; SICaF2 = 0;
    SIFeS = 0; SIFeSAm = 0; SITrot = 0; SIZnS = 0; SIBaCO3 = 0;
    SISrCO3 = 0; SIPbS = 0; SIZnCO3 = 0; SIAmSilica = 0;
    SIQuartz = 0; SIDiopside = 0; SIChrysotile = 0;
    SIGreenalite = 0;

    // --- Copy mcInit → mc, maInit → ma
    for (int c = 0; c < NumCat; c++) {
        mc[c] = mcInit[c];
    }
    for (int a = 0; a < NumAn; a++) {
        ma[a] = maInit[a];
    }

    Alk = AlkInit;
    TH4SiO4 = TH4SiO4Init;
    TNH4 = TNH4Init;
    TH3BO3 = TH3BO3Init;
    TAc = TAcInit;
    TFe = TFeInit;

    if (Run_MassTransfer == 1 && LoopTP > 1) {
        mc[iCa] = mcMT[iCa];
        mc[iBa] = mcMT[iBa];
        ma[iSO4] = maMT[iSO4];
        Alk = AlkMT;

        mc[iFe] = mcMT[iFe];
        TFe = TFe_MT;
    }

    // --- Core thermodynamic calculations
    CalcIonicStrength();
    C1_ThermodynamicEquilConsts();
    C2_PitzerActCoefs_T_P_ISt(gNeut, aH2O, TK, TC, PBar, Patm);

    if (useEOS == 0)
    {
        PengRobinson3();
        RatioOilBPoints = fRatioOilBPoints(API);

        //Call C4_SSPEquilCalcs(0, 5, 2, KspCalcite)
        C4_SSPEquilCalcs(0, 4, 2, KspCalcite);
        fmn();

        C2_PitzerActCoefs_T_P_ISt(gNeut, aH2O, TK, TC, PBar, Patm);
        //Call C4_SSPEquilCalcs(0, 5, 2, KspCalcite)
        C4_SSPEquilCalcs(0, 4, 2, KspCalcite);

        PBubblePt = nTCH4 / t1 + nTCO2 / t2 + nTH2S / t3;

        if (Run_MassTransfer == 1) {
            if (Ppsia > PBubblePt)
                VgTP_MT = 0;
            else
                VgTP_MT = Vg_ZRT * Znew * RAtm * 14.696 * TK / 1000.0;

            VO_MT = VO * 159 / 1000.0;
            VW_MT = VW * 159 / 1000.0;
            QBrineFlow = VW_MT * 1e6 / 86400.0;
        }
    }

    // ================================
    //        useEOS != 0 分支
    // ================================

    if (useEOS != 0)
    {
        if (a1 < 20.0 / ISt) a1 = 2;
        else a1 = 20.0 / ISt;

        for (int c = 0; c < NumCat; c++) mc[c] = mcInit[c] * a1;
        for (int a = 0; a < NumAn; a++) ma[a] = maInit[a] * a1;

        Alk = AlkInit * a1;
        TH4SiO4 = TH4SiO4Init * a1;
        TNH4 = TNH4Init * a1;
        TH3BO3 = TH3BO3Init * a1;
        TAc = TAcInit * a1;
        TFe = TFeInit * a1;

        CalcIonicStrength();
        C1_ThermodynamicEquilConsts();
        C2_PitzerActCoefs_T_P_ISt(gNeut, aH2O, TK, TC, PBar, Patm);

        DoubleVec temp_gNeut(2);
        temp_gNeut[0] = gNeut[iCO2aq]; temp_gNeut[1] = gNeut[iH2Saq];
        MultiPhaseFlash(mf_ParametersWereRead, mf_TCr, mf_PCr, mf_Omega, mf_MWgas, mf_kPr, mf_c0, mf_c1, TK, PBar, total_moles, z, temp_gNeut,
            aH2O, density, compositions, phi, Compr, beta, zOutput, mass_phase, MW_Phase, &No_Phases);

        mass_w = total_moles * beta[2] * compositions[14][3] * 0.01801528;

        if (compositions[14][3] < 0.5 || mass_w < 1e-7)
        {
            //当此标记被触发时，水会蒸发并退出计算。
            simContext.errmsg[7] = 8;
            ISt = 0; pH = 0; rhoTP = 0; H2Oevap = 1;
            goto exit_label_500;
        }

        if (mass_w / mass_w_0 < 0.05) {
            simContext.errmsg[7] = 8;
            ISt = 0; pH = 0; rhoTP = 0; H2Oevap = 1;
            goto exit_label_500;
        }

        // ---- rescale concentrations
        for (int c = 0; c < NumCat; c++)
            mc[c] = mcInit[c] * mass_w_0 / mass_w;

        for (int a = 0; a < NumAn; a++)
            ma[a] = maInit[a] * mass_w_0 / mass_w;

        Alk = AlkInit * mass_w_0 / mass_w;
        TH4SiO4 = TH4SiO4Init * mass_w_0 / mass_w;
        TNH4 = TNH4Init * mass_w_0 / mass_w;
        TH3BO3 = TH3BO3Init * mass_w_0 / mass_w;
        TAc = TAcInit * mass_w_0 / mass_w;
        double mass_wOld = mass_w;

        // ---- mass transfer override
        if (Run_MassTransfer == 1 && LoopTP > 1)
        {
            mc[iCa] = mcMT[iCa] * mass_w_0 / mass_w;
            mc[iBa] = mcMT[iBa] * mass_w_0 / mass_w;
            ma[iSO4] = maMT[iSO4] * mass_w_0 / mass_w;
            Alk = AlkMT * mass_w_0 / mass_w;
        }

        CalcIonicStrength();
        if (ISt > 25) {
            simContext.errmsg[7] = 8;
            ISt = 0; pH = 0; rhoTP = 0; H2Oevap = 1;
            goto exit_label_500;
        }

        C2_PitzerActCoefs_T_P_ISt(gNeut, aH2O, TK, TC, PBar, Patm);
        temp_gNeut[0] = gNeut[iCO2aq]; temp_gNeut[1] = gNeut[iH2Saq];
        MultiPhaseFlash(mf_ParametersWereRead, mf_TCr, mf_PCr, mf_Omega, mf_MWgas, mf_kPr, mf_c0, mf_c1, TK, PBar, total_moles, z, temp_gNeut,
            aH2O, density, compositions, phi, Compr, beta, zOutput, mass_phase, MW_Phase, &No_Phases);

        mass_w = total_moles * beta[2] * compositions[14][3] * 0.01801528;
        Iteration = 0;

        while ((mass_wOld - mass_w) * (mass_wOld - mass_w) > 0.0001 * mass_w * mass_w)
        {
            if (compositions[14][3] < 0.5 || mass_w < 1e-7)
            {
                simContext.errmsg[7] = 8; ISt = 0; pH = 0; rhoTP = 0; H2Oevap = 1;
                goto exit_label_500;
            }

            for (int c = 0; c < NumCat; c++)
                mc[c] = mc[c] * mass_wOld / mass_w;
            for (int a = 0; a < NumAn; a++)
                ma[a] = ma[a] * mass_wOld / mass_w;

            Alk *= mass_wOld / mass_w;
            TH4SiO4 *= mass_wOld / mass_w;
            TNH4 *= TNH4Init * mass_wOld / mass_w;
            TH3BO3 *= TH3BO3Init * mass_wOld / mass_w;
            TAc *= TAcInit * mass_wOld / mass_w;
            mass_wOld = mass_w;

            CalcIonicStrength();
            C2_PitzerActCoefs_T_P_ISt(gNeut, aH2O, TK, TC, PBar, Patm);
            temp_gNeut[0] = gNeut[iCO2aq]; temp_gNeut[1] = gNeut[iH2Saq];
            MultiPhaseFlash(mf_ParametersWereRead, mf_TCr, mf_PCr, mf_Omega, mf_MWgas, mf_kPr, mf_c0, mf_c1, TK, PBar, total_moles, z, temp_gNeut,
                aH2O, density, compositions, phi, Compr, beta, zOutput, mass_phase, MW_Phase, &No_Phases);

            mass_w = total_moles * beta[2] * compositions[14][3] * 0.01801528;
            mass_w = 0.5 * (mass_wOld + mass_w);

            Iteration++;
            if (Iteration > 10)
            {
                simContext.errmsg[7] = 8; ISt = 0; pH = 0; rhoTP = 0; H2Oevap = 1;
                goto exit_label_500;
            }
        }

        // more calculations
        //Call C4_SSPEquilCalcs(0, 5, 2, KspCalcite)
        C4_SSPEquilCalcs(0, 4, 2, KspCalcite);
        fmn();
        C2_PitzerActCoefs_T_P_ISt(gNeut, aH2O, TK, TC, PBar, Patm);
        //Call C4_SSPEquilCalcs(0, 5, 2, KspCalcite)
        C4_SSPEquilCalcs(0, 4, 2, KspCalcite);

        // Bubble point evaluation
        QPBubblePt = 1;
        if (No_Phases == 3) QPBubblePt = 0;
        else if (No_Phases == 2) {
            if ((beta[0] > 0 && density[0] < 0.3) ||
                (beta[1] > 0 && density[1] < 0.3))
                QPBubblePt = 0;
        }

        if (Run_MassTransfer == 1)
        {
            if (QPBubblePt == 1) VgTP_MT = 0;
            else VgTP_MT = beta[0] * total_moles * Compr[0] * RBar * TK / PBar / 1000.0;

            VO_MT = VO * 159 / 1000.0;
            VW_MT = VW * 159 / 1000.0;
        }
    }

    // ================================
    //   ShellMultiFlash not used
    // ================================
    if (RunShellMultiflash != 1)
    {
        // 计算饱和度指数
        SICal = log10(mc[iCa] * simContext.HCO3 * gCat[iCa] * gNCat[iCa] *
            gAn[iHCO3] * gNAn[iHCO3] * K2HCO3 / (aH * KspCalcite));
        // vb 0.565472 cpp 0.564576 gCat vb 0.062964 cpp 0.0628701
        // CO3 vb 4.183682E-06 cpp 4.182635E-06
        SIDol = log10(mc[iCa] * mc[iMg] * pow(simContext.CO3, 2) *
            gCat[iCa] * gCat[iMg] * pow(gAn[iCO3], 2) / KspDol);

        SISid = log10(mc[iFe] * simContext.HCO3 * gCat[iFe] * gAn[iHCO3] * K2HCO3 / (aH * KspSiderite));

        SIBar = log10(mc[iBa] * ma[iSO4] * gCat[iBa] * gAn[iSO4] / KspBarite) + dSIMeOHBar + dSIMEGBar;

        SIGyp = log10(mc[iCa] * ma[iSO4] * gCat[iCa] * gAn[iSO4] * pow(aH2O, 2) *
            pow(gNMean[iCaSO42H2O], 2) * pow(aNH2O, 2) / KspGypsum);

        SIHemi = log10(mc[iCa] * ma[iSO4] * gCat[iCa] * gAn[iSO4] * pow(aH2O, 0.5) *
            pow(gNMean[ihemiCaSO4], 2) * pow(aNH2O, 0.5) / KspHemihydrate);

        SIAn = log10(mc[iCa] * ma[iSO4] * gCat[iCa] * gAn[iSO4] * pow(gNMean[iCaSO4], 2) / KspAnhydrite);

        SICel = log10(mc[iSr] * ma[iSO4] * gCat[iSr] * gAn[iSO4] * pow(gNMean[iSrSO4], 2) / KspCelestite);

        SIHal = log10(mc[iNa] * ma[iCl] * gCat[iNa] * gAn[iCl] / KspHalite) + dSIMeOHHal + dSIMEGHal;

        SICaF2 = log10((mc[iCa] * pow(ma[intF], 2) * gCat[iCa] * pow(gAn[intF], 2) / KspCaF2));

        SIBaCO3 = log10(mc[iBa] * simContext.HCO3 * gCat[iBa] * gNCat[iBa] * gAn[iHCO3] * gNAn[iHCO3] * K2HCO3 / (aH * KspBaCO3));

        SISrCO3 = log10(mc[iSr] * simContext.HCO3 * gCat[iSr] * gNCat[iSr] * gAn[iHCO3] * gNAn[iHCO3] * K2HCO3 / (aH * KspSrCO3));

        // If Use_ZnCl2Const = 1 Then 注释掉的条件语句
        SIZnS = log10(mc[iZn] * simContext.HS * gCat[iZn] * gAn[iHS] * gNAn[iHS] / aH / KspZnS);

        SIPbS = log10(mc[iPb] * simContext.HS * gCat[iPb] * gAn[iHS] * gNAn[iHS] / aH / KspPbS);

        SIZnCO3 = log10(mc[iZn] * simContext.CO3 * gCat[iZn] * gAn[iCO3] / KspZnCO3);

        // added by Dai 2016
        SIFeSAm = log10(mc[iFe] * simContext.HS * gCat[iFe] * gAn[iHS] * gNAn[iHS] / aH / KspFeSAm);
        SIFeS = log10(mc[iFe] * simContext.HS * gCat[iFe] * gAn[iHS] * gNAn[iHS] / aH / KspFeS);
        SITrot = log10(mc[iFe] * simContext.HS * gCat[iFe] * gAn[iHS] * gNAn[iHS] / aH / KspTrot);

        if (xMeOH > 0 && mc[iCa] * simContext.HCO3 > 0)
        {
            double temp = mc[iCa] * simContext.HCO3;
            double logTerm = log10(temp);
            double denominator = 1.0 + pow(logTerm, 2);
            SICal = SICal + (10.0432 * logTerm / denominator) * xMeOH;
            SIHal = SIHal - (0.0696 * exp(xMeOH)) * mc[iCa] * ma[iCl] * xMeOH;
        }

        if (xMEG > 0 && mc[iCa] * simContext.HCO3 > 0)
        {
            double temp = mc[iCa] * simContext.HCO3;
            double logTerm = log10(temp);
            double denominator = 1.0 + pow(logTerm, 2);
            SICal = SICal + (8.4131 * logTerm / denominator) * xMEG;
        }

        // Silicate SI values. No MeOH/MEG is considered at this time.
        // Also, all SiO2 is assumed to remain as H4SiO4 and not ionized.
        SIAmSilica = log10(H4SiO4 * gNeut[iH4SiO4aq] / KspAmSilica);
        SIQuartz = log10(H4SiO4 * gNeut[iH4SiO4aq] / KspQuartz);
        // vb -0.630987 cpp -0.0.631902
        SIDiopside = log10((mc[iCa] * mc[iMg] * pow(H4SiO4, 2) * gCat[iCa] *
            gCat[iMg] * pow(gNeut[iH4SiO4aq], 2) / pow(aH, 4)) / KspDiopside);

        SIChrysotile = log10((pow(mc[iMg], 3) * pow(H4SiO4, 2) *
            pow(gCat[iMg], 3) * pow(gNeut[iH4SiO4aq], 2) / pow(aH, 6)) / KspChrysotile);
        // vb -0.01688 cpp -0.017644
        SIGreenalite = log10((pow(mc[iFe], 3) * pow(H4SiO4, 2) *
            pow(gCat[iFe], 3) * pow(gNeut[iH4SiO4aq], 2) / pow(aH, 6)) / KspGreenalite);

        OH = KH2O / (pow(10, -pH) * gAn[iOH]);

        SIMgOH2 = log10((mc[iMg] * gCat[iMg] * pow(OH * gAn[iOH], 2)) / KspMgOH2);
        SICaOH2 = log10((mc[iCa] * gCat[iCa] * pow(OH * gAn[iOH], 2)) / KspCaOH2);

        rhoTP = CalcRhoTP(TK, TC, PBar, Patm);
    }

exit_label_500:
    return;
}



double fbInhBar(int InhNo, double SI)
{
    return bi[InhNo][0]
        + bi[InhNo][1] * SI
        + bi[InhNo][2] / TK
        + bi[InhNo][3] * log10(1.0 / aH)
        + bi[InhNo][4] * fabs(log10(mc[iBa] / ma[iSO4]));
}


double flogT0Bar(double SI)
{
    double value;

    value = -3.153194285
        + (-0.92635504 / SI)
        + (716.694987 / TK)
        + (1879.905802 / (SI * TK))
        + 0.189075542 * fabs(log10(mc[iBa] / ma[iSO4]));   // Zhaoyi Dai 2020

    // Amy 2013 correction for MeOH
    if (xMeOH > 0.0)
        value += 1.1136 * SI * (1.0 - 1.2976 * xMeOH) * xMeOH;

    // Amy 2013 correction for MEG
    if (xMEG > 0.0)
        value += 4.8464 * xMEG;

    return value;
}


double fbInhCal(int InhNo, double SI)
{
    return
        ci[InhNo][0]                         // ci(InhNo,1)
        + ci[InhNo][1] * SI                  // ci(InhNo,2) * SI
        + ci[InhNo][2] / TK                  // ci(InhNo,3) / TK
        + ci[InhNo][3] * log10(1.0 / aH)     // ci(InhNo,4) * Log10(1/aH)
        + ci[InhNo][4]                       // ci(InhNo,5)
        * fabs(log10(mc[iCa] / simContext.HCO3));     // Abs(Log10(mc(iCa)/HCO3))
}


double fCinhCal(double SI, double tInh)
{
    double bInhCalMixed;

    if (InhNoCal == 20)
    {
        // 注意：所有 VB 下标均需 -1
        bInhCal[simContext.InhNo1] = pow(10.0, fbInhCal(simContext.InhNo1, SI));
        bInhCal[simContext.InhNo2] = pow(10.0, fbInhCal(simContext.InhNo2, SI));

        bInhCalMixed = simContext.FracInhNo1 * bInhCal[simContext.InhNo1]
            + (1.0 - simContext.FracInhNo1) * bInhCal[simContext.InhNo2];

        return (1.0 / bInhCalMixed)
            * log10(fSafetyCal * tInh / t0Cal);
    }
    else
    {
        bInhCal[InhNoCal] = pow(10.0, fbInhCal(InhNoCal, SI));

        return (1.0 / bInhCal[InhNoCal])
            * log10(fSafetyCal * tInh / t0Cal);
    }
}

double fCinhBar(double SI, double tInh)
{
    double bInhBarMixed;

    if (InhNoBar == 20)
    {
        bInhBar[simContext.InhNo1] = pow(10.0, fbInhBar(simContext.InhNo1, SI));

        bInhBar[simContext.InhNo2 - 1] = pow(10.0, fbInhBar(simContext.InhNo2, SI));

        bInhBarMixed =
            simContext.FracInhNo1 * bInhBar[simContext.InhNo1 - 1] +
            (1.0 - simContext.FracInhNo1) * bInhBar[simContext.InhNo2 - 1];

        return (1.0 / bInhBarMixed) * log10(fSafetyBar * tInh / t0Bar);
    }
    else
    {
        bInhBar[InhNoBar] = pow(10.0, fbInhBar(InhNoBar, SI));

        return (1.0 / bInhBar[InhNoBar]) * log10(fSafetyBar * tInh / t0Bar);
    }
}

double flogT0Cal(double SI)
{
    return
        -5.36
        + 1.5 / SI
        + 1779.17 / TK
        + 0.95 * fabs(log10(mc[iCa] / ma[iHCO3]));
}


double flogT0Gyp(double SI)
{
    return (-6.2971 - 0.2212 / SI + 2171.2067 / (TK * pow(SI, 0.2852)) + 1.715 / (1.0 + sqrt(ISt)));//ISt ^ 0.5
}


double fbInhGyp(int InhNo, double SI)
{
    return gi[InhNo - 1][0] + gi[InhNo - 1][2] / (TK * SI) + gi[InhNo - 1][3] *
        log10(1.0 / aH) + gi[InhNo - 1][4] * abs(log10(mc[iCa] / ma[iSO4]));
}

double fCinhGyp(double SI, double tInh)
{

    bInhGyp[InhNoGyp - 1] = pow(10.0, fbInhGyp(InhNoGyp, SI));

    return (1.0 / bInhGyp[InhNoGyp - 1]) *
        log10(fSafetyGyp * tInh / t0Gyp);
}


double flogT0An(double SI)
{
    return (2.15 - 2.83 / SI - 885.8 / TK + 1766.3 / (SI * TK));
}


double fbInhAn(int InhNo, double SI)
{
    return ai[InhNo - 1][0] + ai[InhNo - 1][1] * SI + ai[InhNo - 1][2] / TK + ai[InhNo - 1][3]
        * log10(1.0 / aH) + ai[InhNo - 1][4] * abs(log10(mc[iCa] / ma[iSO4]));
}


double fCinhAn(double SI, double tInh)
{
    bInhAn[InhNoAn - 1] = pow(10, fbInhAn(InhNoAn, SI));
    return (1 / bInhAn[InhNoAn - 1]) * log10(fsafetyAn * tInh / t0An);
}


double flogT0Cel(double SI)
{
    return  -1.713 - 3.411 / SI + 2646.1 / (SI * TK);
}


double fbInhCel(double InhNo, double SI)
{
    return celi[InhNo - 1][0] + celi[InhNo - 1][1] * SI + celi[InhNo - 1][2] / TK;
}


double fCinhCel(double SI, double tInh)
{
    bInhCel[InhNoCel - 1] = pow(10, fbInhCel(InhNoCel, SI));
    return (1.0 / bInhCel[InhNoCel - 1]) * log10(fSafetyCel * tInh / t0Cel);
}



void B6_InhibitorNeeded()
{
    /* ---------- Barite ---------- */
    fSafetyBar = 1.0;
    ConcInhBar = 0.0;

    if (mc[iBa] * ma[iSO4] > 0.0)
    {
        if (LoopTP == 1 || simContext.RunWhatIf == 1 || simContext.Run1000Cases == 1 || simContext.LoopResChem == 1)
        {
            InhNoBar = simContext.InhNo;

            if (simContext.SelectInh == 1)
            {
                //此处C语言未作下标偏移，因为涉及到InhNoBar等
                //因此，对于数组[i]，要进行数组[i - 1], 传入i时不能-1，因为在fbInhBar里会做偏移
                int iMaxBar = 1;
                double bInhBarMax = fbInhBar(iMaxBar, SIBar);

                for (int i = 1; i <= 19; i++)
                {
                    bInhBar[i - 1] = fbInhBar(i, SIBar);
                    if (bInhBar[i - 1] > bInhBarMax)
                    {
                        iMaxBar = i;
                        bInhBarMax = bInhBar[i - 1];
                    }
                }
                InhNoBar = iMaxBar;
            }

            if (simContext.SelectInh == 1)
                InhNoBar = 2;   /* VB: always choose BHPMP */

            /* Excel writes → remove → optional store */
            /*
            Worksheets("Barite").Cells(2, 6) = InhName(InhNoBar)
            Worksheets("Input").Cells(49, 10) = InhName(InhNoBar)
            */
        }

        if (SIBar > 0.001)
        {
            double BarExpon10 = flogT0Bar(SIBar);

            if (log10(tInh) > BarExpon10 && BarExpon10 < 8.0)
            {
                t0Bar = pow(10.0, BarExpon10);
                ConcInhBar = fCinhBar(SIBar, tInh);
            }
        }
    }


    /* ---------- Calcite ---------- */
    fSafetyCal = 1.0;
    ConcInhCal = 0.0;
    ConcInhGyp;

    if (mc[iCa] * simContext.HCO3 > 0.0)
    {
        if (LoopTP == 1 || simContext.RunWhatIf == 1 || simContext.Run1000Cases == 1 || simContext.LoopResChem == 1)
        {
            InhNoCal = simContext.InhNo;

            if (simContext.SelectInh == 1)
            {
                int iMaxCal = 1;
                double bInhCalMax = fbInhCal(iMaxCal, SICal);

                for (int i = 1; i <= 19; i++)
                {
                    bInhCal[i - 1] = fbInhCal(i, SICal);
                    if (bInhCal[i - 1] > bInhCalMax)
                    {
                        iMaxCal = i;
                        bInhCalMax = bInhCal[i - 1];
                    }
                }
                InhNoCal = iMaxCal;
            }

            if (simContext.SelectInh == 1)
                InhNoCal = 1;
            /*
            Worksheets("Calcite").Cells(2, 7) = InhName(InhNoCal)
            Worksheets("Input").Cells(47, 10) = InhName(InhNoCal)
            */

        }

        if (SICal > 0.001)
        {
            double CalExpon10 = flogT0Cal(SICal);

            if (log10(tInh) > CalExpon10 && CalExpon10 < 8.0)
            {
                double t0Cal = pow(10.0, CalExpon10);
                ConcInhCal = fCinhCal(SICal, tInh);
            }
        }
    }


    /* ---------- Gypsum ---------- */
    fSafetyGyp = 1.0;
    ConcInhGyp = 0.0;

    if (mc[iCa] * ma[iSO4] > 0.0)
    {
        InhNoGyp = 4;

        //Worksheets("Input").Cells(51, 10) = InhName(InhNoGyp) 

        if (SIGyp > 0.1)
        {
            double GypExpon10 = flogT0Gyp(SIGyp);

            if (log10(tInh) > GypExpon10 && GypExpon10 < 8.0)
            {
                if (TK < 373.0)
                {
                    t0Gyp = pow(10.0, GypExpon10);
                    ConcInhGyp = fCinhGyp(SIGyp, tInh);
                }
                else
                {
                    ConcInhGyp = 0;  /* VB Null */
                }
            }
        }
    }

    /* ---------- Anhydrite ---------- */
    fsafetyAn = 1.0;
    ConcInhAn = 0.0;

    if (mc[iCa] * ma[iSO4] > 0.0)
    {
        InhNoAn = 4;
        // Worksheets("Input").Cells(53, 10) = InhName(InhNoAn)

        if (SIAn > 0.1)
        {
            double AnExpon10 = flogT0An(SIAn);

            if (log10(tInh) > AnExpon10 && AnExpon10 < 8.0)
            {
                if (TK > 373.0)
                {
                    t0An = pow(10.0, AnExpon10);
                    ConcInhAn = fCinhAn(SIAn, tInh);
                }
                else
                    ConcInhAn = 0;
            }
        }
    }

    /* ---------- Celestite (SrSO4) ---------- */
    fSafetyCel = 1.0;
    ConcInhCel = 0.0;

    if (mc[iSr] * ma[iSO4] > 0.0)
    {
        if (simContext.InhNo <= 11)
            InhNoCel = 3;      /* phosphonate → DTPMP */
        else if (simContext.InhNo <= 14)
            InhNoCel = 12;     /* carboxylates → PPCA */
        else
            InhNoCel = 17;     /* others → PVS */

        if (simContext.SelectInh == 1)
            InhNoCel = 3;      /* always choose DTPMP */

        //Worksheets("Input").Cells(55, 10) = InhName(InhNoCel)

        if (SICel > 0.001)
        {
            double CelExpon10 = flogT0Cel(SICel);

            if (log10(tInh) > CelExpon10 && CelExpon10 < 8.0)
            {
                t0Cel = pow(10.0, CelExpon10);
                ConcInhCel = fCinhCel(SICel, tInh);
            }
        }
    }
}


void B7_ScaleRisk()
{
    // Barite MIC=0: Barite index=1, calcite index=2
    if (mc[iBa] > 1e-7 && ma[iSO4] > 1e-7)
    {
        double BarExpon10 = log10(tInh);

        double a = -1.1136 * (xMeOH - 1.2976 * pow(xMeOH, 2)) * TK;
        double b = 716.694987 + (BarExpon10 - -3.153194285
            - 0 * mc[iCa]
            - 0.189075542 * fabs(log10(mc[iBa] / ma[iSO4]))
            - 4.8404 * xMEG) * TK;
        double cc = (-0.92635504 * TK - 1879.905802);

        if (xMeOH == 0)
        {
            SIRisk[InhNoBar][LoopTP][1][1] = -cc / b;
        }
        else
        {
            double qroot = -0.5 * (b - sqrt(b * b - 4 * a * cc));
            double root1 = cc / qroot;
            qroot = -0.5 * (b + sqrt(b * b - 4 * a * cc));
            double root2 = cc / qroot;

            SIRisk[InhNoBar][LoopTP][1][1] = root2;
        }

        if (SIRisk[InhNoBar][LoopTP][1][1] < 0)
            SIRisk[InhNoBar][LoopTP][1][1] = 0;

        double tInhRisk = tInh;

        for (int iRisk = 1; iRisk <= simContext.NoRiskcalc - 1; iRisk++)
        {
            double SIRiskLow = SIRisk[InhNoBar][LoopTP][iRisk][1];
            double SIRiskHigh = 5.0;
            double SIBarRisk = 0;

            for (int k = 1; k <= 10; k++)
            {
                SIBarRisk = (SIRiskLow + SIRiskHigh) / 2.0;
                double t0Bar = pow(10, flogT0Bar(SIBarRisk));

                if (t0Bar == 0) {
                    SIBarRisk = 0.0;
                    break;
                }

                double Dep1 = fCinhBar(SIBarRisk, tInhRisk) - simContext.ConcInhBarRisk[iRisk];

                if (Dep1 > 0)
                    SIRiskHigh = SIBarRisk;
                else
                    SIRiskLow = SIBarRisk;
            }

            SIRisk[InhNoBar][LoopTP][iRisk + 1][1] = SIBarRisk;
            if (SIRisk[InhNoBar][LoopTP][iRisk + 1][1] < 0)
                SIRisk[InhNoBar][LoopTP][iRisk + 1][1] = 0;
        }
    }

    // Calcite MIC = 0
    if (mc[iCa] > 1e-7 && simContext.HCO3 > 1e-7)
    {
        double CalExpon10 = log10(tInh);
        double b = 1876.4 + (CalExpon10 - 4.22) * TK;
        double cc = (13.8 * TK - 6259.6);

        double root2 = -cc / b;
        SIRisk[InhNoCal][LoopTP][1][2] = root2;

        if (SIRisk[InhNoCal][LoopTP][1][2] < 0)
            SIRisk[InhNoCal][LoopTP][1][2] = 0;

        double tInhRisk = tInh;

        for (int iRisk = 1; iRisk <= simContext.NoRiskcalc - 1; iRisk++)
        {
            double SIRiskLow = SIRisk[InhNoCal][LoopTP][iRisk][2];
            double SIRiskHigh = 5.0;
            double SICalRisk = 0;

            for (int k = 1; k <= 10; k++)
            {
                SICalRisk = (SIRiskLow + SIRiskHigh) / 2.0;
                double t0Cal = pow(10, flogT0Cal(SICalRisk));

                if (t0Cal == 0) {
                    SICalRisk = 0.0;
                    break;
                }

                double Dep1 = fCinhCal(SICalRisk, tInhRisk) - simContext.ConcInhCalRisk[iRisk];
                if (Dep1 > 0)
                    SIRiskHigh = SICalRisk;
                else
                    SIRiskLow = SICalRisk;
            }

            if (InhNoCal < 6 || InhNoCal > 9)
            {
                SIRisk[InhNoCal][LoopTP][iRisk + 1][2] = SICalRisk;
                if (SIRisk[InhNoCal][LoopTP][iRisk + 1][2] < 0)
                    SIRisk[InhNoCal][LoopTP][iRisk + 1][2] = 0;
            }

            if (InhNoCal == 14)
            {
                if (simContext.InhNo1 > 5 && simContext.InhNo1 < 10)
                    SIRisk[InhNoCal][LoopTP][iRisk + 1][2] = NAN;

                if (simContext.InhNo2 > 5 && simContext.InhNo2 < 10)
                    SIRisk[InhNoCal][LoopTP][iRisk + 1][2] = NAN;
            }
        }
    }

    // Anhydrite MIC = 0
    if (mc[iCa] > 1e-5 && ma[iSO4] > 1e-4)
    {
        if (TK > 373)
        {
            double AnExpon10 = log10(tInh);
            double b = 885.8 + (AnExpon10 - 2.15) * TK;
            double cc = 2.83 * TK - 1766.3;

            double root2 = -cc / b;

            SIRisk[InhNoAn][LoopTP][1][3] = root2;
            if (root2 < 0)
                SIRisk[InhNoAn][LoopTP][1][3] = 0;

            double tInhRisk = tInh;

            for (int iRisk = 1; iRisk <= simContext.NoRiskcalc - 1; iRisk++)
            {
                double SIRiskLow = SIRisk[InhNoAn][LoopTP][iRisk][3];
                double SIRiskHigh = 5.0;
                double SIAnRisk = 0;

                for (int k = 1; k <= 10; k++)
                {
                    SIAnRisk = (SIRiskLow + SIRiskHigh) / 2;
                    double t0An = pow(10, flogT0An(SIAnRisk));

                    if (t0An == 0)
                    {
                        SIAnRisk = 0.0;
                        break;
                    }

                    double Dep1 = fCinhAn(SIAnRisk, tInhRisk) - simContext.ConcInhAnRisk[iRisk];
                    if (Dep1 > 0)
                        SIRiskHigh = SIAnRisk;
                    else
                        SIRiskLow = SIAnRisk;
                }

                SIRisk[InhNoAn][LoopTP][iRisk + 1][3] = SIAnRisk;
                if (SIRisk[InhNoAn][LoopTP][iRisk + 1][3] < 0)
                    SIRisk[InhNoAn][LoopTP][iRisk + 1][3] = 0;
            }
        }
    }
}


void LoopTPSI(double& API)
{
    B5_CalculateSIvalues(API);
    pH_before_precipitation = pH;  // 保存沉淀前 pH，用于 SqSoft

    if (H2Oevap != 1)
        B6_InhibitorNeeded();   // 计算抑制剂量

    // Calcite: 保存 LoopTP = 1 时的基准 SI
    if (LoopTP == 1)
    {
        SICalBh = SICal;
        SIBarBH = SIBar;
        SIGypBH = SIGyp;
        SIHemiBH = SIHemi;
        SIAnBH = SIAn;
        SICelBH = SICel;
        SIHalBH = SIHal;
        SICaF2BH = SICaF2;
        SISidBH = SISid;
        SIFeSBH = SIFeS;
        SIZnSBH = SIZnS;
        SIBaCO3BH = SIBaCO3;
        SISrCO3BH = SISrCO3;
    }

    // ---- SqueezeSoftPitzer 需要的输出 ----
    double TCO2BH, pHBH;
    double BHConcInhCal, BHConcInhBar, WHConcInhCal, WHConcInhBar, SICalWH, SIBarWH;
    if (RunGoalSeek != 1)
    {
        if (LoopTP == 1)
        {
            TCO2BH = simContext.HCO3 + simContext.CO3 + simContext.CO2aq;
            pHBH = pH_before_precipitation;
            BHConcInhCal = ConcInhCal;
            BHConcInhBar = ConcInhBar;
        }
        if (LoopTP == 10)
        {
            WHConcInhCal = ConcInhCal;
            WHConcInhBar = ConcInhBar;
            SICalWH = SICal;
            SIBarWH = SIBar;
        }
    }
    else { // RunGoalSeek == 1
        if (LoopTP == 1) {
            TCO2BH = simContext.HCO3 + simContext.CO3 + simContext.CO2aq;
            pHBH = pH_before_precipitation;
            BHConcInhCal = ConcInhCal;
            BHConcInhBar = ConcInhBar;
        }
        if (LoopTP == 2) {
            WHConcInhCal = ConcInhCal;
            WHConcInhBar = ConcInhBar;
            SICalWH = SICal;
            SIBarWH = SIBar;
        }
    }

    // ----- 计算 ΔSI -----
    double dSIGyp, dSIHemi, dSIAn, dSICel, dSICaF2;

    if (RunGoalSeek != 1)
    {
        dSICal = SICal - SICalBh;
        dSIBar = SIBar - SIBarBH;
        dSIGyp = SIGyp - SIGypBH;
        dSIHemi = SIHemi - SIHemiBH;
        dSIAn = SIAn - SIAnBH;
        dSICel = SICel - SICelBH;
        dSIHal = SIHal - SIHalBH;
        dSICaF2 = SICaF2 - SICaF2BH;
        dSISid = SISid - SISidBH;
        dSIFeS = SIFeS - SIFeSBH;

        if (RunSSP == 1 && H2Oevap != 1)
        {
            B7_ScaleRisk();   // Compute risk
        }
    }

    // ----- 全部沉淀量清零 -----

    pptCal = pptBar = pptSid = pptGyp = pptHemi = pptAn = pptCel = pptHal = 0;
    pptFeS = pptZnS = pptZnCO3 = pptPbS = pptCaF2 = 0;
    pptMgOH2 = pptCaOH2 = pptFeSAm = pptTrot = 0;
    pptAmSilica = pptQuartz = pptGreenalite = pptDiopside = pptChrysotile = 0;

    pptFeCO3_NoMassTransfer = pptFeCO3_MassTransfer = 0;
    pptCalcite_MassTransfer = pptCalcite_NoMassTransfer = 0;
    pptBarite_MassTransfer = pptBarite_NoMassTransfer = 0;
    pptCelestite_MassTransfer = pptCelestite_NoMassTransfer = 0;
    pptFeS_MassTransfer = pptFeS_NoMassTransfer = 0;
}

double ViscWatIst, CpPerMl;

void LoopTPVisHeatCap() {
    // Density calculation
    double dens0 = 999.83952;
    a1 = 16.945176; double a2 = -0.0079870401, a3 = -0.000046170461;
    double a4 = 0.00000010556302, a5 = -2.8054253e-10, a6 = 0.01687985;
    double a7 = 0.043922, a8 = -0.000076, a9 = 0.000000126, a10 = 0.00000000319;

    double dens1atm = (dens0 + a1 * TC + a2 * pow(TC, 2) + a3 * pow(TC, 3) +
        a4 * pow(TC, 4) + a5 * pow(TC, 5)) / (1 + a6 * TC);
    dens = dens1atm + (a7 + a8 * TC + a9 * pow(TC, 2) + a10 * pow(TC, 3)) * (Patm - 1);

    // Salt fraction
    double WtFrSalt = ISt / (ISt + 20.0);
    double SumM = 31.83 * WtFrSalt / (1 - WtFrSalt);

    a1 = 0.0046581; a2 = 0.00014712; a3 = -0.0000017413; a4 = 0.0000000044668;
    a5 = 0.0026539; a6 = -0.000073598; a7 = 0.0000009171; a8 = -0.0000000027177;

    double Den0overDenIst = (1 - WtFrSalt) * (1 + 0.001 * dens1atm * SumM * (
        (a1 + a2 * TC + a3 * pow(TC, 2) + a4 * pow(TC, 3)) +
        (a5 + a6 * TC + a7 * pow(TC, 2) + a8 * pow(TC, 3)) * sqrt(ISt)
        ));
    double DenIst = 0.001 * dens / Den0overDenIst;

    // Kinematic viscosity of pure water
    double ViscWatKinam0 = pow(10, -1.464 + 205.4 / (TK - 153.0));
    if (TC > 150) {
        ViscWatKinam0 = pow(10, -1.416 + 175.2 / (TK - 177.2));
    }

    // Relative viscosity calculation
    a1 = 0.18788; a2 = -0.0037532; a3 = 0.0019142; a4 = -0.007951;
    a5 = -0.69447; a6 = 0.013646; a7 = -0.001413; a8 = -0.0054073;

    double RelativeKinViscWater = pow(10, (a1 * ISt + a2 * pow(ISt, 2) + a3 * pow(ISt, 3) + a4 * ISt * TC) +
        (a5 * ISt + a6 * pow(ISt, 2) + a7 * pow(ISt, 3) + a8 * ISt * TC) * log10(ViscWatKinam0));
    double ViscWatKinamIst = RelativeKinViscWater * ViscWatKinam0;
    ViscWatIst = ViscWatKinamIst * DenIst * exp(0.0000567 * Ppsia); // dynamic viscosity in cP

    // Heat capacity calculation
    double RCp = 8.314;

    // Polynomial coefficients for Cp
    double Cp01 = (1.42910599e-12 * pow(TK, 6) - 3.0405349e-9 * pow(TK, 5) + 2.68432001e-6 * pow(TK, 4) -
        0.00125862425 * pow(TK, 3) + 0.330571467 * pow(TK, 2) - 46.1155258 * TK + 2678.82979) * RCp;
    double Cp02 = (-4.75805453e-11 * pow(TK, 6) + 1.0269399e-7 * pow(TK, 5) - 9.21125253e-5 * pow(TK, 4) +
        0.043953835 * pow(TK, 3) - 11.7702434 * pow(TK, 2) + 1677.55304 * TK - 99455.8841) * RCp;
    double aj = (-7.62506036e-13 * pow(TK, 6) + 1.70111051e-9 * pow(TK, 5) - 1.56757109e-6 * pow(TK, 4) +
        0.000765652877 * pow(TK, 3) - 0.209217557 * pow(TK, 2) + 30.360998 * TK - 1827.83946) * RCp;
    double b0Cp = (-6.43096664e-11 * pow(TK, 6) + 1.3960701e-7 * pow(TK, 5) - 0.000125978879 * pow(TK, 4) +
        0.0604944374 * pow(TK, 3) - 16.3067026 * pow(TK, 2) + 2340.14685 * TK - 139737.502) * 1e-6;
    double b1Cp = (-1.723664970418e-10 * pow(TK, 6) + 3.741076393426e-7 * pow(TK, 5) - 0.0003375079792185 * pow(TK, 4) +
        0.1620205208996 * pow(TK, 3) - 43.65554038758 * pow(TK, 2) + 6261.081592614 * TK - 373498.3550791) * 1e-6;
    double cj = (6.65778657e-12 * pow(TK, 6) - 1.44552177e-8 * pow(TK, 5) + 1.30444127e-5 * pow(TK, 4) -
        0.00626293561 * pow(TK, 3) + 1.68757513 * pow(TK, 2) - 242.00911 * TK + 14433.8114) * 1e-6;

    // Pressure adjustment
    double Cp01400 = (1.1151399759e-12 * pow(TK, 6) - 2.3524142054e-9 * pow(TK, 5) + 2.0567964964e-6 * pow(TK, 4) -
        0.00095381633879 * pow(TK, 3) + 0.24738105724 * pow(TK, 2) - 34.015476356 * TK + 1945.6172051) * RCp;
    double Cp02400 = (-3.2956183298e-11 * pow(TK, 6) + 7.0838438096e-8 * pow(TK, 5) - 6.3240160287e-5 * pow(TK, 4) +
        0.030014769748 * pow(TK, 3) - 7.9888685643 * pow(TK, 2) + 1130.8904696 * TK - 66540.444671) * RCp;
    double AJ400 = (-3.368047257e-12 * pow(TK, 6) + 7.1124201481e-9 * pow(TK, 5) - 6.2311623578e-6 * pow(TK, 4) +
        0.0028997187171 * pow(TK, 3) - 0.75599605111 * pow(TK, 2) + 104.72624175 * TK - 6022.1138681) * RCp;
    double b0400 = (-3.1301714299e-11 * pow(TK, 6) + 6.792050067e-8 * pow(TK, 5) - 6.1270715171e-5 * pow(TK, 4) +
        0.029417824308 * pow(TK, 3) - 7.9306534931 * pow(TK, 2) + 1138.6367071 * TK - 68057.583888) * 1e-6;
    double b1400 = (-1.7188829774e-10 * pow(TK, 6) + 3.7305662383e-7 * pow(TK, 5) - 0.00033655325035 * pow(TK, 4) +
        0.16156170942 * pow(TK, 3) - 43.532509995 * pow(TK, 2) + 6243.6263496 * TK - 372474.54806) * 1e-6;
    double CJ400 = (4.3505452423e-12 * pow(TK, 6) - 9.431670242e-9 * pow(TK, 5) + 8.5005436339e-6 * pow(TK, 4) -
        0.004077397537 * pow(TK, 3) + 1.097972865 * pow(TK, 2) - 157.41127409 * TK + 9388.9944053) * 1e-6;

    // Linear pressure interpolation
    double Cp01PT = Cp01 - (Cp01 - Cp01400) * (Patm - 1) / (400 - 1);
    double Cp02PT = Cp02 - (Cp02 - Cp02400) * (Patm - 1) / (400 - 1);
    double AJPT = aj - (aj - AJ400) * (Patm - 1) / (400 - 1);
    double b0PT = b0Cp - (b0Cp - b0400) * (Patm - 1) / (400 - 1);
    double b1PT = b1Cp - (b1Cp - b1400) * (Patm - 1) / (400 - 1);
    double CJPT = cj - (cj - CJ400) * (Patm - 1) / (400 - 1);

    double bCp = 1.2;
    double neu = 2;
    double hISt = log(1 + bCp * sqrt(ISt)) / (2 * bCp);
    double alphaCp = 2;
    double zM = 1;
    double zX = -1;
    double BJNaCl = b0PT + 2 * b1PT * (1 - (1 + alphaCp * sqrt(ISt)) * exp(-alphaCp * sqrt(ISt))) / (alphaCp * alphaCp * ISt);

    double n1Water = 1000.0 / 18.0153;
    double CpPhi = Cp02PT + neu * fabs(zM * zX) * AJPT * hISt - 2 * 1 * 1 * RCp * TK * TK * (ISt * BJNaCl + ISt * ISt * 1 * 1 * CJPT);
    double Cp0 = Cp01PT * 55.508;
    double CpCp = n1Water * Cp01PT + ISt * CpPhi;
    double massCp = 1000 + ISt * 58.45;
    double SaltMass = ISt * 58.45;
    double VtotalCp = (1000 + SaltMass) / simContext.rho25c;
    double CpPerGram = CpCp / massCp;
    CpPerMl = (CpCp / VtotalCp) * 0.239; // cal/(ml*C)
}


void MassTransferCoefficients(double& API)
{
    /*
    Note:this function has  a major problem: the variable `n` is used in the variables <thick_lsl> and
    shear_lsl>. However, n is used in loops in many places.
    */
    int n = 9;

    // Convert cm to m
    double PipeID_m = PipeID / 100.0;
    double PipeL_m = PipeL / 100.0;

    // Superficial velocities
    VgTP_MT = fmax(0.000001, VgTP_MT);
    double V_s_g = VgTP_MT / 86400.0 / (pi * pow(PipeID_m / 2.0, 2));
    double V_s_l = (VO_MT + VW_MT) / 86400.0 / (pi * pow(PipeID_m / 2.0, 2));
    double V_s_m = V_s_g + V_s_l;

    // Reduced gas properties
    double P_reduced_g = Ppsia / (708.75 - 57.7 * SGG);
    double T_reduced_g = (TF + 460.0) / (169.0 + 314.0 * SGG);

    // Ratio_gas_visc coefficients
    a1 = -2.462; double  a2 = 2.97, a3 = -0.2862, a4 = 0.008054, a5 = 2.808, a6 = -3.498;
    double a7 = 0.3603, a8 = -0.01044, a9 = -0.7933, a10 = 1.396, a11 = -0.1491, a12 = 0.00441;
    double a13 = 0.08393, a14 = -0.1864, a15 = 0.02033, a16 = -0.0006095;

    double Ratio_gas_visc = a1 + a2 * P_reduced_g + a3 * pow(P_reduced_g, 2) + a4 * pow(P_reduced_g, 3)
        + T_reduced_g * (a5 + a6 * P_reduced_g + a7 * pow(P_reduced_g, 2) + a8 * pow(P_reduced_g, 3))
        + pow(T_reduced_g, 2) * (a9 + a10 * P_reduced_g + a11 * pow(P_reduced_g, 2) + a12 * pow(P_reduced_g, 3))
        + pow(T_reduced_g, 3) * (a13 + a14 * P_reduced_g + a15 * pow(P_reduced_g, 2) + a16 * pow(P_reduced_g, 3));

    double Visc_g = ((1.709 / 100000.0 - 2.062 / 1000000.0 * SGG) * TF + 8.188 / 1000.0 - 6.15 / 1000.0 * log10(SGG))
        / T_reduced_g * exp(Ratio_gas_visc);

    // Liquid properties
    double SGO = 141.5 / (131.5 + API);
    double SGL = (VO_MT * SGO + VW_MT * rhoTP * 1000.0 / dens) / (VO_MT + VW_MT);

    double rho_gas_MT = SGG * (Patm * 28.97 / (0.08206 * TK));
    double rho_liquid_MT = SGL * dens;

    // Surface tension
    double surf_tension_og = (39.0 - 0.2571 * API - (TF - 68.0) * 1.5 / 32.0) / 1000.0;
    double surf_tension_wg = -0.11538 * pow(1 - TK / 647.096, 2.471) + 0.2151 * pow(1 - TK / 647.096, 1.233);
    double surf_tension_lg = (VO_MT * surf_tension_og + VW_MT * surf_tension_wg) / (VO_MT + VW_MT);

    // Dimensionless numbers
    double N_v_l = V_s_l * pow(rho_liquid_MT / surf_tension_lg / 9.8, 0.25);
    double N_v_g = V_s_g * pow(rho_liquid_MT / surf_tension_lg / 9.8, 0.25);
    double N_dia = PipeID_m * sqrt(rho_liquid_MT * 9.8 / surf_tension_lg);
    double N_visc_l = ViscWatIst / 1000.0 * pow(9.8 / pow(surf_tension_lg, 3) / rho_liquid_MT, 0.25);

    // Liquid hold up
    double X_1 = log10(N_visc_l + 3.0);
    double CN_l = pow(10, -2.69851 + 0.15841 * X_1 - 0.551 * pow(X_1, 2) + 0.54785 * pow(X_1, 3) - 0.12195 * pow(X_1, 4));
    double X_2 = N_v_l * CN_l / (pow(N_v_g, 0.575) * N_dia) * pow(PBar / 1.01325, 0.1);
    double X_3 = N_v_g * pow(N_visc_l, 0.38) / pow(N_dia, 2.14);
    double a_sign = (X_3 - 0.012) / fabs(X_3 - 0.012);
    X_3 = (1 - a_sign) / 2 * 0.012 + (1 + a_sign) / 2 * X_3;
    double fai_MT = 0.91163 - 4.82176 * X_3 + 1232.25 * pow(X_3, 2) - 22253.6 * pow(X_3, 3) + 116174.3 * pow(X_3, 4);
    if (X_3 <= 0.01) fai_MT = 1.0;

    Hold_l = fai_MT * (-0.10307 + 0.61777 * (log10(X_2) + 6) - 0.63295 * pow(log10(X_2) + 6, 2)
        + 0.29598 * pow(log10(X_2) + 6, 3) - 0.0401 * pow(log10(X_2) + 6, 4));

    if (VgTP_MT < (VO_MT + VW_MT) * 0.018)
        Hold_l = 1.0;
    else if (Hold_l<0.01 || VgTP_MT>(VO_MT + VW_MT) * 1000 || Hold_l >= 1.0)
        Hold_l = 1.0 - 0.5 * (1.0 + V_s_m / 0.8 - sqrt(pow(1 + V_s_m / 0.8, 2) - 4.0 * V_s_g / 0.8));
    else if (Hold_l < 0)
        Hold_l = V_s_l / V_s_m;

    // Flow pattern
    if (V_s_g > 3.1 * pow(9.8 * surf_tension_lg * (rho_liquid_MT - rho_gas_MT) / pow(rho_gas_MT, 2), 0.25))
    {
        Flow_Pattern = "Annular";
        Eff_D = (1 - sqrt(1 - Hold_l)) * PipeID_m;
        Eff_V = V_s_l / Hold_l;
    }
    else if (V_s_g > 0.429 * V_s_l + 0.546 * pow(9.8 * surf_tension_lg * (rho_liquid_MT - rho_gas_MT) / pow(rho_liquid_MT, 2), 0.25))
    {
        Flow_Pattern = "SlugChurn";
        Eff_D = PipeID_m;
        Eff_V = V_s_l / Hold_l;
    }
    else
    {
        Flow_Pattern = "Bubbly";
        Eff_D = PipeID_m;
        Eff_V = (VO_MT + VW_MT + VgTP_MT) / 86400.0 / (pi * pow(Eff_D / 2.0, 2));
    }

    // Diffusion coefficients
    const double com_factor_for_DiffC = 10000.0 / pi / (ViscWatIst / 1000.0) * TK / 4.71 * 1.3806488e-23;
    //原句：DiffC[iCa] = 10000.0 / (radiusC[iCa] * 1e-10) / pi / (ViscWatIst / 1000.0) * TK / 4.71 * 1.3806488e-23;
    // 将重复的部分抽取出来
    DiffC[iCa] = com_factor_for_DiffC / (radiusC[iCa] * 1e-10);
    DiffC[iFe] = com_factor_for_DiffC / (radiusC[iFe] * 1e-10);
    DiffC[iBa] = com_factor_for_DiffC / (radiusC[iBa] * 1e-10);
    DiffC[iSr] = com_factor_for_DiffC / (radiusC[iSr] * 1e-10);

    const double com_factor_for_DiffA = 10000.0 / pi / (ViscWatIst / 1000.0) * TK / 4.0 * 1.3806488e-23;
    //原句：DiffA[iCO3] = 10000.0 / (radiusA[iCO3] * 1e-10) / pi / (ViscWatIst / 1000.0) * TK / 4.0 * 1.3806488e-23;
    DiffA[iCO3] = com_factor_for_DiffA / (radiusA[iCO3] * 1e-10);
    DiffA[iSO4] = com_factor_for_DiffA / (radiusA[iSO4] * 1e-10);
    DiffA[iHCO3] = com_factor_for_DiffA / (radiusA[iHCO3] * 1e-10);
    DiffA[iHS] = com_factor_for_DiffA / (radiusA[iHS] * 1e-10);

    Diff[iCaHCO32] = 3.0 / (1.0 / DiffC[iCa] + 2.0 / DiffA[iHCO3]);
    Diff[iFeHCO32] = 3.0 / (1.0 / DiffC[iFe] + 2.0 / DiffA[iHCO3]);
    Diff[iFeHS2] = 3.0 / (1.0 / DiffC[iFe] + 2.0 / DiffA[iHS]);
    Diff[iBaSO4] = 2.0 / (1.0 / DiffC[iBa] + 1.0 / DiffA[iSO4]);
    Diff[iSrSO4] = 2.0 / (1.0 / DiffC[iSr] + 1.0 / DiffA[iSO4]);
    Diff[iCaSO4] = 2.0 / (1.0 / DiffC[iCa] + 1.0 / DiffA[iSO4]);

    // ----------------- Mass transfer coefficients (continued) -----------------

    // Reynolds number and friction factor
    double ReNO = rhoTP * 1000.0 * Eff_V * Eff_D / (ViscWatIst / 1000.0);
    double f_Churchill = pow(pow(8.0 / ReNO, 12.0) + pow(pow(2.457 * log(1.0 / (pow(7.0 / ReNO, 0.9) + 0.27 * PipeRough / 1000.0 / Eff_D)), 16.0) + pow(37530.0 / ReNO, 16.0), -3.0 / 2.0), 1.0 / 12.0) * 8.0;
    if (ReNO < 7.0) f_Churchill = 64.0 / ReNO;

    // ----------------- Barite (BaSO4) -----------------
    double ScNO = ViscWatIst / 1000.0 / (rhoTP * 1000.0) / (Diff[iBaSO4] / 10000.0);
    double LengthEntry = LoopTP * PipeL_m;
    double GzNO = ReNO * ScNO * (Eff_D / LengthEntry);
    double ShNO_La = 3.66 + 0.0668 * GzNO / (1.0 + 0.04 * pow(GzNO, 2.0 / 3.0));
    double ShNO;

    if (ReNO < 2300.0) {
        Flow_Regime = "Laminar";
        ShNO = ShNO_La;
    }
    else if (ReNO < 3100.0) {
        Flow_Regime = "Transitional";
        double ShNO_lc = 3.66 + 0.0668 * (GzNO / ReNO * 2100.0) / (1.0 + 0.04 * pow(GzNO / ReNO * 2100.0, 2.0 / 3.0));
        ShNO = pow(pow(ShNO_La, 10.0) + pow(1.0 / (pow(ShNO_lc, 2.0)) + pow(4.8 + 0.079 * ReNO * pow(f_Churchill / 8.0, 0.5) * ScNO / pow(1 + pow(ScNO, 0.8), 5.0 / 6.0), -2.0), -5.0), 1.0 / 10.0);
    }
    else {
        Flow_Regime = "Turbulent";
        ShNO = (f_Churchill / 8.0) * (ReNO - 1000.0) * ScNO / (1.0 + 12.7 * pow(f_Churchill / 8.0, 0.5) * (pow(ScNO, 2.0 / 3.0) - 1.0));
    }
    km[iBaSO4] = ShNO * (Diff[iBaSO4] / 10000.0) / Eff_D * 100.0;

    // ----------------- Calcite (CaHCO3) -----------------
    ScNO = ViscWatIst / 1000.0 / (rhoTP * 1000.0) / (Diff[iCaHCO32] / 10000.0);
    LengthEntry = LoopTP * PipeL_m;
    GzNO = ReNO * ScNO * (Eff_D / LengthEntry);
    ShNO_La = 3.66 + 0.0668 * GzNO / (1.0 + 0.04 * pow(GzNO, 2.0 / 3.0));

    if (ReNO < 2300.0) {
        ShNO = ShNO_La;
    }
    else if (ReNO < 3100.0) {
        double ShNO_lc = 3.66 + 0.0668 * (GzNO / ReNO * 2100.0) / (1.0 + 0.04 * pow(GzNO / ReNO * 2100.0, 2.0 / 3.0));
        ShNO = pow(pow(ShNO_La, 10.0) + pow(1.0 / (pow(ShNO_lc, 2.0)) + pow(4.8 + 0.079 * ReNO * pow(f_Churchill / 8.0, 0.5) * ScNO / pow(1 + pow(ScNO, 0.8), 5.0 / 6.0), -2.0), -5.0), 1.0 / 10.0);
    }
    else {
        ShNO = (f_Churchill / 8.0) * (ReNO - 1000.0) * ScNO / (1.0 + 12.7 * pow(f_Churchill / 8.0, 0.5) * (pow(ScNO, 2.0 / 3.0) - 1.0));
    }
    km[iCaHCO32] = ShNO * (Diff[iCaHCO32] / 10000.0) / Eff_D * 100.0;

    // ----------------- Siderite (FeHCO3) -----------------
    ScNO = ViscWatIst / 1000.0 / (rhoTP * 1000.0) / (Diff[iFeHCO32] / 10000.0);
    LengthEntry = LoopTP * PipeL_m;
    GzNO = ReNO * ScNO * (Eff_D / LengthEntry);
    ShNO_La = 3.66 + 0.0668 * GzNO / (1.0 + 0.04 * pow(GzNO, 2.0 / 3.0));

    if (ReNO < 2300.0)
    {
        ShNO = ShNO_La;
    }
    else if (ReNO < 3100.0)
    {
        double ShNO_lc = 3.66 + 0.0668 * (GzNO / ReNO * 2100.0) / (1.0 + 0.04 * pow(GzNO / ReNO * 2100.0, 2.0 / 3.0));
        ShNO = pow(pow(ShNO_La, 10.0) + pow(1.0 / (pow(ShNO_lc, 2.0)) + pow(4.8 + 0.079 * ReNO * pow(f_Churchill / 8.0, 0.5) * ScNO / pow(1 + pow(ScNO, 0.8), 5.0 / 6.0), -2.0), -5.0), 1.0 / 10.0);
    }
    else
    {
        ShNO = (f_Churchill / 8.0) * (ReNO - 1000.0) * ScNO / (1.0 + 12.7 * pow(f_Churchill / 8.0, 0.5) * (pow(ScNO, 2.0 / 3.0) - 1.0));
    }
    km[iFeHCO32] = ShNO * (Diff[iFeHCO32] / 10000.0) / Eff_D * 100.0;
    // Surface-limited rate adjustment
    km[iFeHCO32] = 1.0 / (1.0 / km[iFeHCO32] + 1.0 / (exp(21.6) * exp(-70400.0 / TK / 8.31) / 10.0));

    // ----------------- Laminar sublayer stability -----------------
    double lgt_barite_MT = 1000.0;
    if (SIBar > 0.001)
    {
        lgt_barite_MT = fmin(flogT0Bar(SIBar) + pow(10.0, fbInhBar(i_Inhibitor_MT, SIBar)) * c_Inhibitor_MT, 1000.0);
    }

    double lgt_calcite_MT = 1000.0;
    if (SICal > 0.001)
    {

        lgt_calcite_MT = fmin(flogT0Cal(SICal) + pow(10, fbInhCal(i_Inhibitor_MT, SICal)) * c_Inhibitor_MT, 1000.0);

        if (i_Inhibitor_MT == 5 || i_Inhibitor_MT == 7 || i_Inhibitor_MT == 8 ||
            i_Inhibitor_MT == 9 || i_Inhibitor_MT == 10 || i_Inhibitor_MT == 14 ||
            i_Inhibitor_MT == 18)
        {
            lgt_calcite_MT = fmin(flogT0Cal(SICal), 1000.0);

            if (LoopTP == 1 && Iter_MT_WI == 1) {
                /*
                MsgBox ("SSP does not include inhibitor impacts on calcite batch crystallization induction time
                for 5-EDTMP, 7-HEBMP, 8-PhnAA, 9-PAPEMP, 10-PBTC, 14-PASP, and 18-AA/AMPS.
                Click OK to continue deposition prediction.")
                */
            }
        }
    }

    Time_lsl_burst = 0;   // initialize

    if (ReNO > 2300) {

        double n_power = 0.8391 * log(ReNO) - 2.5488;

        double thick_lsl = pow(
            pow(11.6, 2) * ViscWatIst / 1000.0 /
            (rhoTP * 1000.0) * pow(Eff_D / 2.0, 1.0 / n) / Eff_V * n,
            n / (n + 1.0)
        );

        double shear_lsl = Eff_V / pow(Eff_D / 2.0, 1.0 / n) / n
            * pow(thick_lsl, 1.0 / n - 1.0) * ViscWatIst / 1000.0;

        Time_lsl_burst = 300.0 * ViscWatIst / 1000.0 / shear_lsl;

        if (Use_lsl_stability == 1) {

            if (SIBar > 0) {
                if (lgt_barite_MT > log10(Time_lsl_burst))
                    km[iBaSO4] = 0;
            }

            if (SICal > 0) {
                if (lgt_calcite_MT > log10(Time_lsl_burst))
                    km[iCaHCO32] = 0;
            }
        }
    }

    /* ******** inhibitor molar weights ******** */
    double MW_Inh[20];
    MW_Inh[0] = 299; MW_Inh[1] = 685; MW_Inh[2] = 573; MW_Inh[3] = 492; MW_Inh[4] = 436;
    MW_Inh[5] = 206; MW_Inh[6] = 249; MW_Inh[7] = 270; MW_Inh[8] = 600; MW_Inh[9] = 270;
    MW_Inh[10] = 389; MW_Inh[11] = 3800; MW_Inh[12] = 116; MW_Inh[13] = 114;
    MW_Inh[14] = 205.5; MW_Inh[15] = 2897; MW_Inh[16] = 107; MW_Inh[17] = 279;

    /* backup original */
    int i_Inhibitor_MT_original = i_Inhibitor_MT;
    double c_Inhibitor_MT_original = c_Inhibitor_MT;


    DoubleVec Keq_inhibitor_dep_barite(20, 0);
    DoubleVec Theta0_inhibitor_dep_barite(20, 0);
    DoubleVec Keq_inhibitor_dep_calcite(20, 0);
    DoubleVec Theta0_inhibitor_dep_calcite(20, 0);
    /* --- barite inhibition constants --- */
    Keq_inhibitor_dep_barite[2] =
        exp(-902.2 + 47533.3 / TK + 131.02 * log(TK));
    Theta0_inhibitor_dep_barite[2] = 0;

    Keq_inhibitor_dep_barite[11] =
        exp(-637.7 + 32441.7 / TK + 93.01 * log(TK));
    Theta0_inhibitor_dep_barite[11] =
        0.00204 * exp(1784.0 / TK);

    Keq_inhibitor_dep_barite[15] =
        exp(-347.1 + 15692.3 / TK + 51.66 * log(TK));
    Theta0_inhibitor_dep_barite[15] =
        0.001021 * exp(1792.0 / TK);

    /* inhibitor group selection */
    if (i_Inhibitor_MT == 1 || i_Inhibitor_MT == 2 || i_Inhibitor_MT == 3 ||
        i_Inhibitor_MT == 4 || i_Inhibitor_MT == 5 || i_Inhibitor_MT == 6 ||
        i_Inhibitor_MT == 7 || i_Inhibitor_MT == 8 || i_Inhibitor_MT == 9 ||
        i_Inhibitor_MT == 10 || i_Inhibitor_MT == 11)
    {
        c_Inhibitor_MT = c_Inhibitor_MT / MW_Inh[i_Inhibitor_MT - 1] * MW_Inh[2];
        i_Inhibitor_MT = 3;
    }
    else if (i_Inhibitor_MT == 12 || i_Inhibitor_MT == 13 ||
        i_Inhibitor_MT == 14 || i_Inhibitor_MT == 15)
    {
        c_Inhibitor_MT = c_Inhibitor_MT / MW_Inh[i_Inhibitor_MT - 1] * MW_Inh[11];
        i_Inhibitor_MT = 12;
    }
    else {
        c_Inhibitor_MT = c_Inhibitor_MT / MW_Inh[i_Inhibitor_MT - 1] * MW_Inh[15];
        i_Inhibitor_MT = 16;
    }

    km[iBaSO4] = km[iBaSO4] *
        (1 - Keq_inhibitor_dep_barite[i_Inhibitor_MT - 1] * c_Inhibitor_MT /
            (1 + Keq_inhibitor_dep_barite[i_Inhibitor_MT - 1] * c_Inhibitor_MT) *
            (1 - Theta0_inhibitor_dep_barite[i_Inhibitor_MT - 1]));

    /* restore for calcite */
    i_Inhibitor_MT = i_Inhibitor_MT_original;
    c_Inhibitor_MT = c_Inhibitor_MT_original;

    /* --- calcite inhibition --- */
    Keq_inhibitor_dep_calcite[2] = exp(-4127.0 / TK + 14.923);
    Theta0_inhibitor_dep_calcite[2] = 0;

    if (i_Inhibitor_MT != 3) {
        c_Inhibitor_MT = c_Inhibitor_MT / MW_Inh[i_Inhibitor_MT - 1] * MW_Inh[2];
        i_Inhibitor_MT = 3;
    }

    km[iCaHCO32] = km[iCaHCO32] *
        (1 - Keq_inhibitor_dep_calcite[i_Inhibitor_MT - 1] * c_Inhibitor_MT /
            (1 + Keq_inhibitor_dep_calcite[i_Inhibitor_MT - 1] * c_Inhibitor_MT) *
            (1 - Theta0_inhibitor_dep_calcite[i_Inhibitor_MT - 1]));

    /* restore final */
    i_Inhibitor_MT = i_Inhibitor_MT_original;
    c_Inhibitor_MT = c_Inhibitor_MT_original;

}

// ==================== 主转换函数 ====================
void LoopTPWrite(ExcelWriter* excelWriter = nullptr) {
    // 创建数据存储对象
    OutputData outputData;

    // 第一部分：条件写入 location 和 depth
    if (simContext.RunStat == 1) {
        if (RunStatSICalcSSP == 1) {
            // VB: Worksheets("Halite analysis").Cells(6 + StatTPnob - LoopTP, 2) = location
            // VB: Worksheets("Halite analysis").Cells(6 + StatTPnob - LoopTP, 3) = depth
            int row = calculateRowWithStat(6, StatTPnob, LoopTP);
            outputData.addData("Halite analysis", row, 2, location);
            outputData.addData("Halite analysis", row, 3, depth);

        }
        else if (RunStatSICalcSSP == 2) {
            if (RunStatGoalSeek == 1) {
                // VB: Worksheets("Output data sheet").Cells(11 + StatTPnob - LoopTP, 16) = location
                // VB: Worksheets("Output data sheet").Cells(11 + StatTPnob - LoopTP, 17) = depth
                int row = calculateRowWithStat(11, StatTPnob, LoopTP);
                outputData.addData("Output data sheet", row, 16, location);
                outputData.addData("Output data sheet", row, 17, depth);
            }
            else {
                // VB: Worksheets("Halite analysis").Cells(6 + StatTPnob - LoopTP, 10) = location
                // VB: Worksheets("Halite analysis").Cells(6 + StatTPnob - LoopTP, 11) = depth
                int row = calculateRowWithStat(6, StatTPnob, LoopTP);
                outputData.addData("Halite analysis", row, 10, location);
                outputData.addData("Halite analysis", row, 11, depth);
            }
        }
        else if (RunStatSICalcSSP == 3) {
            // VB: Worksheets("Output data sheet").Cells(6 + StatTPnob - LoopTP, 8) = location
            // VB: Worksheets("Output data sheet").Cells(6 + StatTPnob - LoopTP, 9) = depth
            int row = calculateRowWithStat(6, StatTPnob, LoopTP);
            outputData.addData("Output data sheet", row, 8, location);
            outputData.addData("Output data sheet", row, 9, depth);
        }
    }

    // 第二部分：当 RunGoalSeek 不为1时的输出
    if (RunGoalSeek != 1) {
        if (simContext.UseSI == 0) {
            // 使用英制单位输出温度和压力
            // Calcite 工作表
            outputData.addData("Calcite", 14 - LoopTP, 1, TF);
            outputData.addData("Calcite", 14 - LoopTP, 2, Ppsia);

            // Barite 工作表
            outputData.addData("Barite", 14 - LoopTP, 1, TF);
            outputData.addData("Barite", 14 - LoopTP, 2, Ppsia);

            // Other SO4s 工作表
            outputData.addData("Other SO4s", 14 - LoopTP, 1, TF);
            outputData.addData("Other SO4s", 14 - LoopTP, 2, Ppsia);

            // Halite 工作表
            outputData.addData("Halite", 14 - LoopTP, 1, TF);
            outputData.addData("Halite", 14 - LoopTP, 2, Ppsia);

            // Sulfides,Fluorite,Carbonates 工作表
            outputData.addData("Sulfides,Fluorite,Carbonates", 14 - LoopTP, 1, TF);
            outputData.addData("Sulfides,Fluorite,Carbonates", 14 - LoopTP, 2, Ppsia);

            // Use Mass Transfer 工作表
            outputData.addData("Use Mass Transfer", 22 - LoopTP, 1, TF);
            outputData.addData("Use Mass Transfer", 22 - LoopTP, 2, Ppsia);

            // Dai 2020 deposition - Run_MassTransfer = 1
            if (Run_MassTransfer == 1) {
                double distance_ft = (LoopTP - 1) * PipeL * 0.0328084; // 转换为英尺

                outputData.addData("Deposition Prediction", 35 - LoopTP, 3, distance_ft);
                outputData.addData("Deposition Prediction", 35 - LoopTP, 4, TF);
                outputData.addData("Deposition Prediction", 35 - LoopTP, 5, Ppsia);

                outputData.addData("Deposition Prediction", 51 - LoopTP, 3, distance_ft);
                outputData.addData("Deposition Prediction", 51 - LoopTP, 4, TF);
                outputData.addData("Deposition Prediction", 51 - LoopTP, 5, Ppsia);

                // 输出流动参数
                outputData.addData("Deposition Prediction", 35 - LoopTP, 6, Flow_Pattern);
                outputData.addData("Deposition Prediction", 35 - LoopTP, 7, Flow_Regime);
                outputData.addData("Deposition Prediction", 35 - LoopTP, 8, ReNO);
                outputData.addData("Deposition Prediction", 35 - LoopTP, 9, Hold_l);
                outputData.addData("Deposition Prediction", 35 - LoopTP, 10, Time_lsl_burst);
                outputData.addData("Deposition Prediction", 35 - LoopTP, 11, lgt_calcite_MT);
                outputData.addData("Deposition Prediction", 35 - LoopTP, 12, lgt_barite_MT);
            }

            // Dai 2020 deposition - Run_MassTransfer_WhatIf = 1
            if (simContext.Run_MassTransfer_WhatIf == 1) {
                double distance_ft = (LoopTP - 1) * PipeL * 0.0328084; // 转换为英尺
                int row_offset = (Iter_MT_WI - 1) * 32;

                outputData.addData("Deposition Prediction_WhatIf", 35 - LoopTP + row_offset, 3, distance_ft);
                outputData.addData("Deposition Prediction_WhatIf", 35 - LoopTP + row_offset, 4, TF);
                outputData.addData("Deposition Prediction_WhatIf", 35 - LoopTP + row_offset, 5, Ppsia);

                outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 3, distance_ft);
                outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 4, TF);
                outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 5, Ppsia);

                // 输出流动参数
                outputData.addData("Deposition Prediction_WhatIf", 35 - LoopTP + row_offset, 6, Flow_Pattern);
                outputData.addData("Deposition Prediction_WhatIf", 35 - LoopTP + row_offset, 7, Flow_Regime);
                outputData.addData("Deposition Prediction_WhatIf", 35 - LoopTP + row_offset, 8, ReNO);
                outputData.addData("Deposition Prediction_WhatIf", 35 - LoopTP + row_offset, 9, Hold_l);
                outputData.addData("Deposition Prediction_WhatIf", 35 - LoopTP + row_offset, 10, Time_lsl_burst);
                outputData.addData("Deposition Prediction_WhatIf", 35 - LoopTP + row_offset, 11, lgt_calcite_MT);
                outputData.addData("Deposition Prediction_WhatIf", 35 - LoopTP + row_offset, 12, lgt_barite_MT);
            }

            // Mg(OH)2,Ca(OH)2 工作表
            outputData.addData("Mg(OH)2,Ca(OH)2", 14 - LoopTP, 1, TF);
            outputData.addData("Mg(OH)2,Ca(OH)2", 14 - LoopTP, 2, Ppsia);

            // Silicates 工作表
            outputData.addData("Silicates", 14 - LoopTP, 1, TF);
            outputData.addData("Silicates", 14 - LoopTP, 2, Ppsia);

            // 根据 RunStat 条件输出
            if (simContext.RunStat == 1) {
                if (RunStatSICalcSSP == 1) {
                    int row = calculateRowWithStat(6, StatTPnob, LoopTP);
                    outputData.addData("Halite analysis", row, 4, TF);
                    outputData.addData("Halite analysis", row, 5, Ppsia);
                }
                else if (RunStatSICalcSSP == 2) {
                    if (RunStatGoalSeek == 1) {
                        int row = calculateRowWithStat(11, StatTPnob, LoopTP);
                        outputData.addData(myname, row, 18, TF);
                        outputData.addData(myname, row, 19, Ppsia);
                    }
                    else {
                        int row = calculateRowWithStat(6, StatTPnob, LoopTP);
                        outputData.addData("Halite analysis", row, 12, TF);
                        outputData.addData("Halite analysis", row, 13, Ppsia);
                    }
                }
                else if (RunStatSICalcSSP == 3) {
                    int row = calculateRowWithStat(6, StatTPnob, LoopTP);
                    outputData.addData("Output data sheet", row, 10, TF);
                    outputData.addData("Output data sheet", row, 11, Ppsia);
                }
            }

        }
        else if (simContext.UseSI == 1) {
            // 使用SI单位制（摄氏度，Bar）
            double TC = (TF - 32) * 5.0 / 9.0;    // 华氏度转摄氏度
            double Pbar = Ppsia / 14.503774;      // psi转bar

            // Calcite 工作表
            outputData.addData("Calcite", 14 - LoopTP, 1, TC);
            outputData.addData("Calcite", 14 - LoopTP, 2, Pbar);

            // Barite 工作表
            outputData.addData("Barite", 14 - LoopTP, 1, TC);
            outputData.addData("Barite", 14 - LoopTP, 2, Pbar);

            // Other SO4s 工作表
            outputData.addData("Other SO4s", 14 - LoopTP, 1, TC);
            outputData.addData("Other SO4s", 14 - LoopTP, 2, Pbar);

            // Halite 工作表
            outputData.addData("Halite", 14 - LoopTP, 1, TC);
            outputData.addData("Halite", 14 - LoopTP, 2, Pbar);

            // Sulfides,Fluorite,Carbonates 工作表
            outputData.addData("Sulfides,Fluorite,Carbonates", 14 - LoopTP, 1, TC);
            outputData.addData("Sulfides,Fluorite,Carbonates", 14 - LoopTP, 2, Pbar);

            // Use Mass Transfer 工作表
            outputData.addData("Use Mass Transfer", 22 - LoopTP, 1, TC);
            outputData.addData("Use Mass Transfer", 22 - LoopTP, 2, Pbar);

            // Dai 2020 deposition - Run_MassTransfer = 1 (SI单位)
            if (Run_MassTransfer == 1) {
                double distance_m = (LoopTP - 1) * PipeL / 100.0; // 转换为米

                outputData.addData("Deposition Prediction", 35 - LoopTP, 3, distance_m);
                outputData.addData("Deposition Prediction", 35 - LoopTP, 4, TC);
                outputData.addData("Deposition Prediction", 35 - LoopTP, 5, Pbar);

                outputData.addData("Deposition Prediction", 51 - LoopTP, 3, distance_m);
                outputData.addData("Deposition Prediction", 51 - LoopTP, 4, TC);
                outputData.addData("Deposition Prediction", 51 - LoopTP, 5, Pbar);

                // 输出流动参数
                outputData.addData("Deposition Prediction", 35 - LoopTP, 6, Flow_Pattern);
                outputData.addData("Deposition Prediction", 35 - LoopTP, 7, Flow_Regime);
                outputData.addData("Deposition Prediction", 35 - LoopTP, 8, ReNO);
                outputData.addData("Deposition Prediction", 35 - LoopTP, 9, Hold_l);
                outputData.addData("Deposition Prediction", 35 - LoopTP, 10, Time_lsl_burst);
            }

            // Dai 2020 deposition - Run_MassTransfer_WhatIf = 1 (SI单位)
            if (simContext.Run_MassTransfer_WhatIf == 1) {
                double distance_m = (LoopTP - 1) * PipeL / 100.0; // 转换为米
                int row_offset = (Iter_MT_WI - 1) * 32;

                outputData.addData("Deposition Prediction_WhatIf", 35 - LoopTP + row_offset, 3, distance_m);
                outputData.addData("Deposition Prediction_WhatIf", 35 - LoopTP + row_offset, 4, TC);
                outputData.addData("Deposition Prediction_WhatIf", 35 - LoopTP + row_offset, 5, Pbar);

                outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 3, distance_m);
                outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 4, TC);
                outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 5, Pbar);

                // 输出流动参数
                outputData.addData("Deposition Prediction_WhatIf", 35 - LoopTP + row_offset, 6, Flow_Pattern);
                outputData.addData("Deposition Prediction_WhatIf", 35 - LoopTP + row_offset, 7, Flow_Regime);
                outputData.addData("Deposition Prediction_WhatIf", 35 - LoopTP + row_offset, 8, ReNO);
                outputData.addData("Deposition Prediction_WhatIf", 35 - LoopTP + row_offset, 9, Hold_l);
                outputData.addData("Deposition Prediction_WhatIf", 35 - LoopTP + row_offset, 10, Time_lsl_burst);
            }

            // Mg(OH)2,Ca(OH)2 工作表
            outputData.addData("Mg(OH)2,Ca(OH)2", 14 - LoopTP, 1, TC);
            outputData.addData("Mg(OH)2,Ca(OH)2", 14 - LoopTP, 2, Pbar);

            // Silicates 工作表
            outputData.addData("Silicates", 14 - LoopTP, 1, TC);
            outputData.addData("Silicates", 14 - LoopTP, 2, Pbar);

            // 根据 RunStat 条件输出
            if (simContext.RunStat == 1) {
                if (RunStatSICalcSSP == 1) {
                    int row = calculateRowWithStat(6, StatTPnob, LoopTP);
                    outputData.addData("Halite analysis", row, 4, TC);
                    outputData.addData("Halite analysis", row, 5, Pbar);
                }
                else if (RunStatSICalcSSP == 2) {
                    if (RunStatGoalSeek == 1) {
                        int row = calculateRowWithStat(11, StatTPnob, LoopTP);
                        outputData.addData(myname, row, 18, TC);
                        outputData.addData(myname, row, 19, Pbar);
                    }
                    else {
                        int row = calculateRowWithStat(6, StatTPnob, LoopTP);
                        outputData.addData("Halite analysis", row, 12, TC);
                        outputData.addData("Halite analysis", row, 13, Pbar);
                    }
                }
                else if (RunStatSICalcSSP == 3) {
                    int row = calculateRowWithStat(6, StatTPnob, LoopTP);
                    outputData.addData("Output data sheet", row, 10, TC);
                    outputData.addData("Output data sheet", row, 11, Pbar);
                }
            }
        }

        // ==================== 检查Amy部分 ====================
        // pH输出
        outputData.addData("Calcite", 14 - LoopTP, 3, pH_before_precipitation - DpHj);

        // 抑制剂浓度清空逻辑
        // 注意：在C++中，我们通常不直接"清空"变量，而是设置特定值或使用可选类型
        // 这里模拟VB中的Empty值设置
        if (InhNoCal == 5 || InhNoCal == 7 || InhNoCal == 8 || InhNoCal == 9 ||
            InhNoCal == 10 || InhNoCal == 14 || InhNoCal == 18) {
            ConcInhCal = std::numeric_limits<double>::quiet_NaN(); // 表示"空"值
        }

        if (InhNoCal == 20) {
            if (simContext.InhNo1 == 5 || simContext.InhNo1 == 7 || simContext.InhNo1 == 9 || simContext.InhNo1 == 8 ||
                simContext.InhNo1 == 10 || simContext.InhNo1 == 14 || simContext.InhNo1 == 18) {
                ConcInhCal = std::numeric_limits<double>::quiet_NaN();
            }
            if (simContext.InhNo2 == 5 || simContext.InhNo2 == 7 || simContext.InhNo2 == 8 || simContext.InhNo2 == 9 ||
                simContext.InhNo2 == 10 || simContext.InhNo2 == 11 || simContext.InhNo2 == 14 || simContext.InhNo2 == 18) {
                ConcInhCal = std::numeric_limits<double>::quiet_NaN();
            }
        }

        // 输出抑制剂浓度
        outputData.addData("Calcite", 14 - LoopTP, 7, ConcInhCal);
        outputData.addData("Calcite", 14 - LoopTP, 8, pHaftercalciteppt - DpHj);
        outputData.addData("Calcite", 14 - LoopTP, 11, ISt * (simContext.rho25c - TDS / 1000000.0));

        outputData.addData("Barite", 14 - LoopTP, 6, ConcInhBar);
        outputData.addData("Halite", 14 - LoopTP, 6, rhoTP);

        outputData.addData("Silicates", 14 - LoopTP, 5, pHafterAmsilicappt - DpHj);
        outputData.addData("Silicates", 14 - LoopTP, 8, pHafterQuartzppt - DpHj);
        outputData.addData("Silicates", 14 - LoopTP, 17, pHafterGreenaliteppt - DpHj);
        outputData.addData("Silicates", 14 - LoopTP, 14, pHafterDiopsideppt - DpHj);
        outputData.addData("Silicates", 14 - LoopTP, 11, pHafterChrysotileppt - DpHj);

        outputData.addData("Mg(OH)2,Ca(OH)2", 14 - LoopTP, 3, pH_before_precipitation - DpHj);
        outputData.addData("Mg(OH)2,Ca(OH)2", 14 - LoopTP, 6, pHafterMgOH2ppt - DpHj);
        outputData.addData("Mg(OH)2,Ca(OH)2", 14 - LoopTP, 9, pHafterCaOH2ppt - DpHj);

        // ==================== 沉淀量计算和输出 ====================
        double density_factor = simContext.rho25c - TDS * 0.000001;
        double halite_density_factor = Rho25cHalite - TDSHalite / 1000000.0;

        if (usePTB == 0) {
            // 使用 mg/L 单位
            // Calcite 沉淀量
            double pptCalcite_mgL = pptCalcite_NoMassTransfer * 100091.0 * density_factor;
            outputData.addData("Calcite", 14 - LoopTP, 6, pptCalcite_mgL);

            // Barite 沉淀量
            double pptBarite_mgL = pptBarite_NoMassTransfer * 233390.0 * density_factor;
            outputData.addData("Barite", 14 - LoopTP, 5, pptBarite_mgL);

            // 其他硫酸盐沉淀量
            outputData.addData("Other SO4s", 14 - LoopTP, 4, pptGyp * 172172.0 * density_factor);
            outputData.addData("Other SO4s", 14 - LoopTP, 6, pptHemi * 145148.0 * density_factor);
            outputData.addData("Other SO4s", 14 - LoopTP, 8, pptAn * 136140.0 * density_factor);
            outputData.addData("Other SO4s", 14 - LoopTP, 10, pptCel * 183680.0 * density_factor);

            // Halite 沉淀量
            double pptHal_mgL = pptHal * 58443.0 * halite_density_factor;
            outputData.addData("Halite", 14 - LoopTP, 5, pptHal_mgL);

            // 硫化物、氟化物、碳酸盐沉淀量
            outputData.addData("Sulfides,Fluorite,Carbonates", 14 - LoopTP, 4, pptFeSAm * 87910.0 * density_factor);
            outputData.addData("Sulfides,Fluorite,Carbonates", 14 - LoopTP, 6, pptFeS_NoMassTransfer * 87910.0 * density_factor);
            outputData.addData("Sulfides,Fluorite,Carbonates", 14 - LoopTP, 8, pptTrot * 87910.0 * density_factor);
            outputData.addData("Sulfides,Fluorite,Carbonates", 14 - LoopTP, 10, pptZnS * 97440.0 * density_factor);
            outputData.addData("Sulfides,Fluorite,Carbonates", 14 - LoopTP, 12, pptCaF2 * 78080.0 * density_factor);
            outputData.addData("Sulfides,Fluorite,Carbonates", 14 - LoopTP, 15, pptFeCO3_NoMassTransfer * 115861.0 * density_factor);
            outputData.addData("Sulfides,Fluorite,Carbonates", 14 - LoopTP, 17, pptZnCO3 * 125417.0 * density_factor);
            outputData.addData("Sulfides,Fluorite,Carbonates", 14 - LoopTP, 19, pptPbS * 239265.0 * density_factor);
            outputData.addData("Sulfides,Fluorite,Carbonates", 14 - LoopTP, 21, pptSrCO3 * 147639.0 * density_factor);
            outputData.addData("Sulfides,Fluorite,Carbonates", 14 - LoopTP, 23, pptBaCO3 * 197349.0 * density_factor);

            // 硅酸盐沉淀量
            outputData.addData("Silicates", 14 - LoopTP, 4, pptAmSilica * 60084.0 * density_factor);
            outputData.addData("Silicates", 14 - LoopTP, 7, pptQuartz * 60084.0 * density_factor);
            outputData.addData("Silicates", 14 - LoopTP, 10, pptChrysotile * 277110.0 * density_factor / 2.0);
            outputData.addData("Silicates", 14 - LoopTP, 13, pptDiopside * 216550.0 * density_factor / 2.0);
            outputData.addData("Silicates", 14 - LoopTP, 16, pptGreenalite * 371770.0 * density_factor / 3.0);

            // 氢氧化物沉淀量
            outputData.addData("Mg(OH)2,Ca(OH)2", 14 - LoopTP, 5, pptMgOH2 * 58321.0 * density_factor);
            outputData.addData("Mg(OH)2,Ca(OH)2", 14 - LoopTP, 8, pptCaOH2 * 74094.0 * density_factor);

        }
        else if (usePTB == 1) {
            // 使用 PTB 单位（磅/千桶）
            const double PTB_FACTOR = 0.35051;

            // Calcite 沉淀量
            double pptCalcite_ptb = pptCalcite_NoMassTransfer * 100091.0 * density_factor * PTB_FACTOR;
            outputData.addData("Calcite", 14 - LoopTP, 6, pptCalcite_ptb);

            // Barite 沉淀量
            double pptBarite_ptb = pptBarite_NoMassTransfer * 233390.0 * density_factor * PTB_FACTOR;
            outputData.addData("Barite", 14 - LoopTP, 5, pptBarite_ptb);

            // 其他硫酸盐沉淀量
            outputData.addData("Other SO4s", 14 - LoopTP, 4, pptGyp * 172172.0 * density_factor * PTB_FACTOR);
            outputData.addData("Other SO4s", 14 - LoopTP, 6, pptHemi * 145148.0 * density_factor * PTB_FACTOR);
            outputData.addData("Other SO4s", 14 - LoopTP, 8, pptAn * 136140.0 * density_factor * PTB_FACTOR);
            outputData.addData("Other SO4s", 14 - LoopTP, 10, pptCel * 183680.0 * density_factor * PTB_FACTOR);

            // Halite 沉淀量
            double pptHal_ptb = pptHal * 58443.0 * halite_density_factor * PTB_FACTOR;
            outputData.addData("Halite", 14 - LoopTP, 5, pptHal_ptb);

            // 硫化物、氟化物、碳酸盐沉淀量
            outputData.addData("Sulfides,Fluorite,Carbonates", 14 - LoopTP, 4, pptFeSAm * 87910.0 * density_factor * PTB_FACTOR);
            outputData.addData("Sulfides,Fluorite,Carbonates", 14 - LoopTP, 6, pptFeS_NoMassTransfer * 87910.0 * density_factor * PTB_FACTOR);
            outputData.addData("Sulfides,Fluorite,Carbonates", 14 - LoopTP, 8, pptTrot * 87910.0 * density_factor * PTB_FACTOR);
            outputData.addData("Sulfides,Fluorite,Carbonates", 14 - LoopTP, 10, pptZnS * 97440.0 * density_factor * PTB_FACTOR);
            outputData.addData("Sulfides,Fluorite,Carbonates", 14 - LoopTP, 12, pptCaF2 * 78080.0 * density_factor * PTB_FACTOR);
            outputData.addData("Sulfides,Fluorite,Carbonates", 14 - LoopTP, 15, pptFeCO3_NoMassTransfer * 115861.0 * density_factor * PTB_FACTOR);
            outputData.addData("Sulfides,Fluorite,Carbonates", 14 - LoopTP, 17, pptZnCO3 * 125417.0 * density_factor * PTB_FACTOR);
            outputData.addData("Sulfides,Fluorite,Carbonates", 14 - LoopTP, 19, pptPbS * 239265.0 * density_factor * PTB_FACTOR);
            outputData.addData("Sulfides,Fluorite,Carbonates", 14 - LoopTP, 21, pptSrCO3 * 147639.0 * density_factor * PTB_FACTOR);
            outputData.addData("Sulfides,Fluorite,Carbonates", 14 - LoopTP, 23, pptBaCO3 * 197349.0 * density_factor * PTB_FACTOR);

            // 硅酸盐沉淀量
            outputData.addData("Silicates", 14 - LoopTP, 4, pptAmSilica * 60084.0 * density_factor * PTB_FACTOR);
            outputData.addData("Silicates", 14 - LoopTP, 7, pptQuartz * 60084.0 * density_factor * PTB_FACTOR);
            outputData.addData("Silicates", 14 - LoopTP, 10, pptChrysotile * 277110.0 * density_factor / 2.0 * PTB_FACTOR);
            outputData.addData("Silicates", 14 - LoopTP, 13, pptDiopside * 216550.0 * density_factor / 2.0 * PTB_FACTOR);
            outputData.addData("Silicates", 14 - LoopTP, 16, pptGreenalite * 371770.0 * density_factor / 3.0 * PTB_FACTOR);

            // 氢氧化物沉淀量
            outputData.addData("Mg(OH)2,Ca(OH)2", 14 - LoopTP, 5, pptMgOH2 * 58321.0 * density_factor * PTB_FACTOR);
            outputData.addData("Mg(OH)2,Ca(OH)2", 14 - LoopTP, 8, pptCaOH2 * 74094.0 * density_factor * PTB_FACTOR);
        }

        // 根据RunStat条件输出Halite数据
        if (simContext.RunStat == 1) {
            double pptHal_output = 0.0;
            if (usePTB == 0) {
                pptHal_output = pptHal * 58443.0 * halite_density_factor;
            }
            else {
                pptHal_output = pptHal * 58443.0 * halite_density_factor * 0.35051;
            }

            if (RunStatSICalcSSP == 1) {
                int row = calculateRowWithStat(6, StatTPnob, LoopTP);
                outputData.addData("Halite analysis", row, 7, pptHal_output);
                outputData.addData("Halite analysis", row, 8, rhoTP);
            }
            else if (RunStatSICalcSSP == 2) {
                if (RunStatGoalSeek == 1) {
                    int row = calculateRowWithStat(11, StatTPnob, LoopTP);
                    outputData.addData(myname, row, 22, rhoTP);
                    outputData.addData(myname, row, 21, pptHal_output);
                }
                else {
                    int row = calculateRowWithStat(6, StatTPnob, LoopTP);
                    outputData.addData("Halite analysis", row, 15, pptHal_output);
                    outputData.addData("Halite analysis", row, 16, rhoTP);
                }
            }
            else if (RunStatSICalcSSP == 3) {
                int row = calculateRowWithStat(6, StatTPnob, LoopTP);
                outputData.addData("Output data sheet", row, 13, pptHal_output);
                outputData.addData("Output data sheet", row, 14, rhoTP);
            }
        }

        // ==================== 气泡点检测输出 ====================
        if (useEOS == 0) {
            std::string bubbleStatus = (Ppsia > PBubblePt) ? "Yes" : "No";
            outputData.addData("Calcite", 14 - LoopTP, 9, bubbleStatus);
        }
        else {
            std::string bubbleStatus = (QPBubblePt == 1) ? "Yes" : "No";
            outputData.addData("Calcite", 14 - LoopTP, 9, bubbleStatus);
        }

        // ==================== 传质相关输出 ====================
        // Use Mass Transfer 工作表
        if (usePTB == 0) {
            outputData.addData("Use Mass Transfer", 22 - LoopTP, 4,
                pptCalcite_NoMassTransfer * 100091.0 * density_factor);
            outputData.addData("Use Mass Transfer", 22 - LoopTP, 7,
                pptBarite_NoMassTransfer * 233390.0 * density_factor);
        }
        else {
            outputData.addData("Use Mass Transfer", 22 - LoopTP, 4,
                pptCalcite_NoMassTransfer * 100091.0 * density_factor * 0.35051);
            outputData.addData("Use Mass Transfer", 22 - LoopTP, 7,
                pptBarite_NoMassTransfer * 233390.0 * density_factor * 0.35051);
        }

        // Dai 2020 deposition - 传质输出
        if (Run_MassTransfer == 1) {
            // Calcite
            if (usePTB == 0) {
                outputData.addData("Deposition Prediction", 51 - LoopTP, 7,
                    pptCalcite_NoMassTransfer * 100091.0 * density_factor);
                outputData.addData("Deposition Prediction", 51 - LoopTP, 11,
                    pptBarite_NoMassTransfer * 233390.0 * density_factor);
                outputData.addData("Deposition Prediction", 51 - LoopTP, 15,
                    pptFeCO3_NoMassTransfer * 233390.0 * density_factor);
            }
            else {
                outputData.addData("Deposition Prediction", 51 - LoopTP, 7,
                    pptCalcite_NoMassTransfer * 100091.0 * density_factor * 0.35051);
                outputData.addData("Deposition Prediction", 51 - LoopTP, 11,
                    pptBarite_NoMassTransfer * 233390.0 * density_factor * 0.35051);
                outputData.addData("Deposition Prediction", 51 - LoopTP, 15,
                    pptFeCO3_NoMassTransfer * 233390.0 * density_factor * 0.35051);
            }
        }

        // WhatIf 情景
        if (simContext.Run_MassTransfer_WhatIf == 1) {
            if (usePTB == 0) {
                outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + (Iter_MT_WI - 1) * 32, 7,
                    pptCalcite_NoMassTransfer * 100091.0 * density_factor);
                outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + (Iter_MT_WI - 1) * 32, 11,
                    pptBarite_NoMassTransfer * 233390.0 * density_factor);
            }
            else {
                outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + (Iter_MT_WI - 1) * 32, 7,
                    pptCalcite_NoMassTransfer * 100091.0 * density_factor * 0.35051);
                outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + (Iter_MT_WI - 1) * 32, 11,
                    pptBarite_NoMassTransfer * 233390.0 * density_factor * 0.35051);
            }
        }

        // 如果有传质，输出传质沉淀量
        if (Run_MassTransfer == 1) {
            if (usePTB == 0) {
                outputData.addData("Use Mass Transfer", 22 - LoopTP, 5,
                    pptCalcite_MassTransfer * 100091.0 * density_factor);
                outputData.addData("Use Mass Transfer", 22 - LoopTP, 8,
                    pptBarite_MassTransfer * 233390.0 * density_factor);
            }
            else {
                outputData.addData("Use Mass Transfer", 22 - LoopTP, 5,
                    pptCalcite_MassTransfer * 100091.0 * density_factor * 0.35051);
                outputData.addData("Use Mass Transfer", 22 - LoopTP, 8,
                    pptBarite_MassTransfer * 233390.0 * density_factor * 0.35051);
            }

            // Deposition Prediction 工作表
            if (usePTB == 0) {
                outputData.addData("Deposition Prediction", 51 - LoopTP, 8,
                    pptCalcite_MassTransfer * 100091.0 * density_factor);
                outputData.addData("Deposition Prediction", 51 - LoopTP, 12,
                    pptBarite_MassTransfer * 233390.0 * density_factor);
                outputData.addData("Deposition Prediction", 51 - LoopTP, 16,
                    pptFeCO3_MassTransfer * 233390.0 * density_factor);
            }
            else {
                outputData.addData("Deposition Prediction", 51 - LoopTP, 8,
                    pptCalcite_MassTransfer * 100091.0 * density_factor * 0.35051);
                outputData.addData("Deposition Prediction", 51 - LoopTP, 12,
                    pptBarite_MassTransfer * 233390.0 * density_factor * 0.35051);
                outputData.addData("Deposition Prediction", 51 - LoopTP, 16,
                    pptFeCO3_MassTransfer * 233390.0 * density_factor * 0.35051);
            }

            // 沉积速率
            if (simContext.UseSI == 1)
            {
                outputData.addData("Deposition Prediction", 51 - LoopTP, 9, pptCalcite_MassTransfer_V); // cm/yr
                outputData.addData("Deposition Prediction", 51 - LoopTP, 13, pptBarite_MassTransfer_V); // cm/yr
                outputData.addData("Deposition Prediction", 51 - LoopTP, 17, pptFeCO3_MassTransfer_V); // cm/yr
                outputData.addData("Deposition Prediction", 51 - LoopTP, 18, m_CR_selected); // mm/yr
            }
            else {
                outputData.addData("Deposition Prediction", 51 - LoopTP, 9, pptCalcite_MassTransfer_V * 0.393701); // inch/yr
                outputData.addData("Deposition Prediction", 51 - LoopTP, 13, pptBarite_MassTransfer_V * 0.393701); // inch/yr
                outputData.addData("Deposition Prediction", 51 - LoopTP, 17, pptFeCO3_MassTransfer_V * 0.393701); // inch/yr
                outputData.addData("Deposition Prediction", 51 - LoopTP, 18, m_CR_selected * 39.37); // mil/yr
            }

            if (Use_Corr_in_Deposition == 0) {
                outputData.addData("Deposition Prediction", 51 - LoopTP, 18, std::string(""));
            }

            // WhatIf 情景
            if (simContext.Run_MassTransfer_WhatIf == 1) {
                int row_offset = (Iter_MT_WI - 1) * 32;

                if (usePTB == 0) {
                    outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 8,
                        pptCalcite_MassTransfer * 100091.0 * density_factor);
                    outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 12,
                        pptBarite_MassTransfer * 233390.0 * density_factor);
                }
                else {
                    outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 8,
                        pptCalcite_MassTransfer * 100091.0 * density_factor * 0.35051);
                    outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 12,
                        pptBarite_MassTransfer * 233390.0 * density_factor * 0.35051);
                }

                if (simContext.UseSI == 1) {
                    outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 9,
                        pptCalcite_MassTransfer_V); // cm/yr
                    outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 13,
                        pptBarite_MassTransfer_V); // cm/yr
                }
                else {
                    outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 9,
                        pptCalcite_MassTransfer_V * 0.393701); // inch/yr
                    outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 13,
                        pptBarite_MassTransfer_V * 0.393701); // inch/yr
                }
            }
        }

        // ==================== 饱和指数(SI)输出 ====================
        if (simContext.useSR == 0) {
            // 输出 SI 值（对数形式）
            outputData.addData("Calcite", 14 - LoopTP, 4, SICal);
            outputData.addData("Calcite", 14 - LoopTP, 5, dSICal);
            outputData.addData("Calcite", 14 - LoopTP, 10, SIDol);

            outputData.addData("Barite", 14 - LoopTP, 3, SIBar);
            outputData.addData("Barite", 14 - LoopTP, 4, dSIBar);

            outputData.addData("Other SO4s", 14 - LoopTP, 3, SIGyp);
            outputData.addData("Other SO4s", 14 - LoopTP, 5, SIHemi);
            outputData.addData("Other SO4s", 14 - LoopTP, 7, SIAn);
            outputData.addData("Other SO4s", 14 - LoopTP, 9, SICel);

            outputData.addData("Halite", 14 - LoopTP, 3, SIHal);
            outputData.addData("Halite", 14 - LoopTP, 4, dSIHal);

            outputData.addData("Sulfides,Fluorite,Carbonates", 14 - LoopTP, 3, SIFeSAm);
            outputData.addData("Sulfides,Fluorite,Carbonates", 14 - LoopTP, 5, SIFeS);
            outputData.addData("Sulfides,Fluorite,Carbonates", 14 - LoopTP, 7, SITrot);
            outputData.addData("Sulfides,Fluorite,Carbonates", 14 - LoopTP, 9, SIZnS);
            outputData.addData("Sulfides,Fluorite,Carbonates", 14 - LoopTP, 11, SICaF2);
            outputData.addData("Sulfides,Fluorite,Carbonates", 14 - LoopTP, 13, SISid);
            outputData.addData("Sulfides,Fluorite,Carbonates", 14 - LoopTP, 14, dSISid);
            outputData.addData("Sulfides,Fluorite,Carbonates", 14 - LoopTP, 16, SIZnCO3);
            outputData.addData("Sulfides,Fluorite,Carbonates", 14 - LoopTP, 18, SIPbS);
            outputData.addData("Sulfides,Fluorite,Carbonates", 14 - LoopTP, 20, SISrCO3);
            outputData.addData("Sulfides,Fluorite,Carbonates", 14 - LoopTP, 22, SIBaCO3);

            outputData.addData("Mg(OH)2,Ca(OH)2", 14 - LoopTP, 7, SICaOH2);
            outputData.addData("Mg(OH)2,Ca(OH)2", 14 - LoopTP, 4, SIMgOH2);

            outputData.addData("Silicates", 14 - LoopTP, 3, SIAmSilica);
            outputData.addData("Silicates", 14 - LoopTP, 6, SIQuartz);
            outputData.addData("Silicates", 14 - LoopTP, 9, SIChrysotile);
            outputData.addData("Silicates", 14 - LoopTP, 12, SIDiopside);
            outputData.addData("Silicates", 14 - LoopTP, 15, SIGreenalite);

            outputData.addData("Use Mass Transfer", 22 - LoopTP, 3, SICal);
            outputData.addData("Use Mass Transfer", 22 - LoopTP, 6, SIBar);

            // Deposition Prediction
            if (Run_MassTransfer == 1) {
                outputData.addData("Deposition Prediction", 51 - LoopTP, 6, SICal);
                outputData.addData("Deposition Prediction", 51 - LoopTP, 10, SIBar);
                outputData.addData("Deposition Prediction", 51 - LoopTP, 14, SISid);
            }

            // WhatIf 情景
            if (simContext.Run_MassTransfer_WhatIf == 1) {
                int row_offset = (Iter_MT_WI - 1) * 32;
                outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 6, SICal);
                outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 10, SIBar);
            }

            // RunStat 条件下的输出
            if (simContext.RunStat == 1) {
                if (RunStatSICalcSSP == 1) {
                    int row = calculateRowWithStat(6, StatTPnob, LoopTP);
                    outputData.addData("Halite analysis", row, 6, SIHal);
                }
                else if (RunStatSICalcSSP == 2) {
                    if (RunStatGoalSeek == 1) {
                        int row = calculateRowWithStat(11, StatTPnob, LoopTP);
                        outputData.addData(myname, row, 20, SIHal);
                    }
                    else {
                        int row = calculateRowWithStat(6, StatTPnob, LoopTP);
                        outputData.addData("Halite analysis", row, 14, SIHal);
                    }
                }
                else if (RunStatSICalcSSP == 3) {
                    int row = calculateRowWithStat(6, StatTPnob, LoopTP);
                    outputData.addData("Output data sheet", row, 12, SIHal);
                }
            }

        }
        else {
            // 输出饱和度比（10^SI）
            outputData.addData("Calcite", 14 - LoopTP, 4, pow(10.0, SICal));
            outputData.addData("Calcite", 14 - LoopTP, 5, pow(10.0, dSICal));
            outputData.addData("Calcite", 14 - LoopTP, 10, pow(10.0, SIDol));

            outputData.addData("Barite", 14 - LoopTP, 3, pow(10.0, SIBar));
            outputData.addData("Barite", 14 - LoopTP, 4, pow(10.0, SIBar - SIBarBH));

            // 需要调用 outputData.addData
            // Worksheets("Other SO4s").Cells(14 - LoopTP, 3) = 10 ^ SIGyp
            // Worksheets("Other SO4s").Cells(14 - LoopTP, 5) = 10 ^ SIHemi
            // Worksheets("Other SO4s").Cells(14 - LoopTP, 7) = 10 ^ SIAn
            // Worksheets("Other SO4s").Cells(14 - LoopTP, 9) = 10 ^ SICel
            // Worksheets("Halite").Cells(14 - LoopTP, 3) = 10 ^ SIHal
            // Worksheets("Halite").Cells(14 - LoopTP, 4) = 10 ^ dSIHal
            // Worksheets("Sulfides,Fluorite,Carbonates").Cells(14 - LoopTP, 3) = 10 ^ SIFeSAm
            // Worksheets("Sulfides,Fluorite,Carbonates").Cells(14 - LoopTP, 5) = 10 ^ SIFeS
            // Worksheets("Sulfides,Fluorite,Carbonates").Cells(14 - LoopTP, 7) = 10 ^ SITrot
            // Worksheets("Sulfides,Fluorite,Carbonates").Cells(14 - LoopTP, 9) = 10 ^ SIZnS
            // Worksheets("Sulfides,Fluorite,Carbonates").Cells(14 - LoopTP, 11) = 10 ^ SICaF2
            // Worksheets("Sulfides,Fluorite,Carbonates").Cells(14 - LoopTP, 13) = 10 ^ SISid
            // Worksheets("Sulfides,Fluorite,Carbonates").Cells(14 - LoopTP, 14) = 10 ^ dSISid
            // Worksheets("Sulfides,Fluorite,Carbonates").Cells(14 - LoopTP, 16) = 10 ^ SIZnCO3
            // Worksheets("Sulfides,Fluorite,Carbonates").Cells(14 - LoopTP, 18) = 10 ^ SIPbS
            // Worksheets("Sulfides,Fluorite,Carbonates").Cells(14 - LoopTP, 20) = 10 ^ SISrCO3
            // Worksheets("Sulfides,Fluorite,Carbonates").Cells(14 - LoopTP, 22) = 10 ^ SIBaCO3
            // Worksheets("Mg(OH)2,Ca(OH)2").Cells(14 - LoopTP, 7) = 10 ^ SICaOH2
            // Worksheets("Mg(OH)2,Ca(OH)2").Cells(14 - LoopTP, 4) = 10 ^ SIMgOH2
            // Worksheets("Silicates").Cells(14 - LoopTP, 3) = 10 ^ SIAmSilica
            // Worksheets("Silicates").Cells(14 - LoopTP, 6) = 10 ^ SIQuartz
            // Worksheets("Silicates").Cells(14 - LoopTP, 9) = 10 ^ SIChrysotile
            // Worksheets("Silicates").Cells(14 - LoopTP, 12) = 10 ^ SIDiopside
            // Worksheets("Silicates").Cells(14 - LoopTP, 15) = 10 ^ SIGreenalite
            // Worksheets("Use Mass Transfer").Cells(22 - LoopTP, 3) = 10 ^ SICal
            // Worksheets("Use Mass Transfer").Cells(22 - LoopTP, 6) = 10 ^ SIBar

            // Deposition Prediction
            if (Run_MassTransfer == 1) {
                outputData.addData("Deposition Prediction", 51 - LoopTP, 6, pow(10.0, SICal));
                outputData.addData("Deposition Prediction", 51 - LoopTP, 10, pow(10.0, SIBar));
                outputData.addData("Deposition Prediction", 51 - LoopTP, 14, pow(10.0, SISid));
            }

            // WhatIf 情景
            if (simContext.Run_MassTransfer_WhatIf == 1) {
                int row_offset = (Iter_MT_WI - 1) * 32;
                outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 6, pow(10.0, SICal));
                outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 10, pow(10.0, SIBar));
            }
        }

        // ==================== 抑制剂风险评估输出 ====================
        if (H2Oevap != 1) {
            // Calcite 风险评估
            if (mc.size() > 0 && mc[0] > 0.0000001 && simContext.HCO3 > 0.0000001) {
                // 假设 mc[0] 对应 iCa，需要根据实际索引调整
                int row = 14 - LoopTP;

                outputData.addData("Calcite", row, 16, outputData.worksheets["Calcite"].numericData[{row, 1}]);
                outputData.addData("Calcite", row, 17, outputData.worksheets["Calcite"].numericData[{row, 2}]);
                outputData.addData("Calcite", row, 18, outputData.worksheets["Calcite"].numericData[{row, 4}]);

                // 清空 SIRisk 数组中的某些值
                if (InhNoCal == 5 || InhNoCal == 7 || InhNoCal == 8 || InhNoCal == 9 ||
                    InhNoCal == 10 || InhNoCal == 14 || InhNoCal == 18) {
                    // 模拟VB中的 Empty 赋值
                    if (SIRisk.size() > InhNoCal && SIRisk[InhNoCal].size() > LoopTP) {
                        SIRisk[InhNoCal][LoopTP][2][2] = std::numeric_limits<double>::quiet_NaN();
                        SIRisk[InhNoCal][LoopTP][3][2] = std::numeric_limits<double>::quiet_NaN();
                    }
                }

                if (InhNoCal == 20) {
                    if (simContext.InhNo1 == 5 || simContext.InhNo1 == 7 || simContext.InhNo1 == 8 || simContext.InhNo1 == 9 ||
                        simContext.InhNo1 == 10 || simContext.InhNo1 == 14 || simContext.InhNo1 == 18 ||
                        simContext.InhNo2 == 5 || simContext.InhNo2 == 7 || simContext.InhNo2 == 8 || simContext.InhNo2 == 9 ||
                        simContext.InhNo2 == 10 || simContext.InhNo2 == 14 || simContext.InhNo2 == 18) {
                        if (SIRisk.size() > InhNoCal && SIRisk[InhNoCal].size() > LoopTP) {
                            SIRisk[InhNoCal][LoopTP][2][2] = std::numeric_limits<double>::quiet_NaN();
                            SIRisk[InhNoCal][LoopTP][3][2] = std::numeric_limits<double>::quiet_NaN();
                        }
                    }
                }

                // 输出风险评估值
                if (simContext.useSR == 0) {
                    if (SIRisk.size() > InhNoCal && SIRisk[InhNoCal].size() > LoopTP) {
                        outputData.addData("Calcite", row, 19, SIRisk[InhNoCal][LoopTP][1][2]);
                        outputData.addData("Calcite", row, 20, SIRisk[InhNoCal][LoopTP][2][2]);
                        outputData.addData("Calcite", row, 21, SIRisk[InhNoCal][LoopTP][3][2]);
                    }
                }
                else {
                    if (SIRisk.size() > InhNoCal && SIRisk[InhNoCal].size() > LoopTP) {
                        outputData.addData("Calcite", row, 19, pow(10.0, SIRisk[InhNoCal][LoopTP][1][2]));
                        outputData.addData("Calcite", row, 20, pow(10.0, SIRisk[InhNoCal][LoopTP][2][2]));
                        outputData.addData("Calcite", row, 21, pow(10.0, SIRisk[InhNoCal][LoopTP][3][2]));
                    }
                }
            }

            // Barite 风险评估（类似处理，为简洁起见省略）
            // If mc(iBa) > 0.0000001 And ma(iSO4) > 0.0000001 Then
            //     Worksheets("Barite").Cells(14 - LoopTP, 16) = Worksheets("Barite").Cells(14 - LoopTP, 1)
            //     Worksheets("Barite").Cells(14 - LoopTP, 17) = Worksheets("Barite").Cells(14 - LoopTP, 2)
            //     Worksheets("Barite").Cells(14 - LoopTP, 18) = Worksheets("Barite").Cells(14 - LoopTP, 3)
            //     If UseSR = 0 Then
            //     Worksheets("Barite").Cells(14 - LoopTP, 19) = SIRisk(InhNoBar, LoopTP, 1, 1)
            //     Worksheets("Barite").Cells(14 - LoopTP, 20) = SIRisk(InhNoBar, LoopTP, 2, 1)
            //     Worksheets("Barite").Cells(14 - LoopTP, 21) = SIRisk(InhNoBar, LoopTP, 3, 1)
            //     Else
            //     Worksheets("Barite").Cells(14 - LoopTP, 19) = 10 ^ SIRisk(InhNoBar, LoopTP, 1, 1)
            //     Worksheets("Barite").Cells(14 - LoopTP, 20) = 10 ^ SIRisk(InhNoBar, LoopTP, 2, 1)
            //     Worksheets("Barite").Cells(14 - LoopTP, 21) = 10 ^ SIRisk(InhNoBar, LoopTP, 3, 1)
            //     End If
            //     End If

            // // Other SO4s 风险评估
            //     If mc(iCa) > 0.00001 And ma(iSO4) > 0.0001 Then
            //     If TK > 373.15 Then
            //         If UseSI = 0 Then
            //         Worksheets("Other SO4s").Cells(15 - LoopTP, 16) = TF
            //         Else
            //         Worksheets("Other SO4s").Cells(15 - LoopTP, 16) = (TF - 32) * 5 / 9
            //         End If
            //         Worksheets("Other SO4s").Cells(15 - LoopTP, 17) = Worksheets("Other SO4s").Cells(14 - LoopTP, 2)
            //         Worksheets("Other SO4s").Cells(15 - LoopTP, 18) = Worksheets("Other SO4s").Cells(14 - LoopTP, 7)
            //         If UseSR = 0 Then
            //         Worksheets("Other SO4s").Cells(15 - LoopTP, 19) = SIRisk(InhNoAn, LoopTP, 1, 3)
            //         Worksheets("Other SO4s").Cells(15 - LoopTP, 20) = SIRisk(InhNoAn, LoopTP, 2, 3)
            //         Worksheets("Other SO4s").Cells(15 - LoopTP, 21) = SIRisk(InhNoAn, LoopTP, 3, 3)
            //         Else
            //         Worksheets("Other SO4s").Cells(15 - LoopTP, 19) = 10 ^ SIRisk(InhNoAn, LoopTP, 1, 3)
            //         Worksheets("Other SO4s").Cells(15 - LoopTP, 20) = 10 ^ SIRisk(InhNoAn, LoopTP, 2, 3)
            //         Worksheets("Other SO4s").Cells(15 - LoopTP, 21) = 10 ^ SIRisk(InhNoAn, LoopTP, 3, 3)
            //         End If
            //     Else
            //         If UseSI = 0 Then
            //         Worksheets("Other SO4s").Cells(15 - LoopTP - 1, 16) = TF
            //         Else
            //         Worksheets("Other SO4s").Cells(15 - LoopTP - 1, 16) = (TF - 32) * 5 / 9
            //         End If
            //     Worksheets("Other SO4s").Cells(15 - LoopTP - 1, 17) = Worksheets("Other SO4s").Cells(14 - LoopTP, 2)
            //     Worksheets("Other SO4s").Cells(15 - LoopTP - 1, 18) = Worksheets("Other SO4s").Cells(14 - LoopTP, 3)
            //     Worksheets("Other SO4s").Cells(15 - LoopTP - 1, 19) = SIRisk(InhNoAn, LoopTP, 1, 3)
            //     Worksheets("Other SO4s").Cells(15 - LoopTP - 1, 20) = SIRisk(InhNoAn, LoopTP, 2, 3)
            //     Worksheets("Other SO4s").Cells(15 - LoopTP - 1, 21) = SIRisk(InhNoAn, LoopTP, 3, 3)
            //     End If
            //     End If

            //     // 'Dai 2020 update inhibition model for celestite
            //     If mc(iSr) > 0.0000001 And ma(iSO4) > 0.0000001 Then
            //         Worksheets("Other SO4s").Cells(26 - LoopTP, 16) = Worksheets("Other SO4s").Cells(14 - LoopTP, 1)
            //         Worksheets("Other SO4s").Cells(26 - LoopTP, 17) = Worksheets("Other SO4s").Cells(14 - LoopTP, 2)
            //         Worksheets("Other SO4s").Cells(26 - LoopTP, 18) = Worksheets("Other SO4s").Cells(14 - LoopTP, 9)
            //         If UseSR = 0 Then
            //             Worksheets("Other SO4s").Cells(26 - LoopTP, 19) = SIRisk(InhNoCel, LoopTP, 1, 4)
            //             Worksheets("Other SO4s").Cells(26 - LoopTP, 20) = SIRisk(InhNoCel, LoopTP, 2, 4)
            //             Worksheets("Other SO4s").Cells(26 - LoopTP, 21) = SIRisk(InhNoCel, LoopTP, 3, 4)
            //         Else
            //             Worksheets("Other SO4s").Cells(26 - LoopTP, 19) = 10 ^ SIRisk(InhNoCel, LoopTP, 1, 4)
            //             Worksheets("Other SO4s").Cells(26 - LoopTP, 20) = 10 ^ SIRisk(InhNoCel, LoopTP, 2, 4)
            //             Worksheets("Other SO4s").Cells(26 - LoopTP, 21) = 10 ^ SIRisk(InhNoCel, LoopTP, 3, 4)
            //         End If
            //     End If
            // ... 其他物质的风险评估
        }
    }

    // ==================== LoopTP特殊条件输出 ====================
    // 第一部分：LoopTP = 1 时的输出
    if (LoopTP == 1) {
        if (simContext.useSR == 0) {
            // 输出各种矿物的饱和指数到Input工作表
            outputData.addData("Input", 13, 9, SICal);
            outputData.addData("Input", 15, 9, SIBar);
            outputData.addData("Input", 17, 9, SIHal);
            outputData.addData("Input", 19, 9, SIGyp);
            outputData.addData("Input", 21, 9, SIHemi);
            outputData.addData("Input", 23, 9, SIAn);
            outputData.addData("Input", 25, 9, SICel);
            outputData.addData("Input", 27, 9, SIFeS);
            outputData.addData("Input", 29, 9, SIZnS);
            outputData.addData("Input", 31, 9, SICaF2);
            outputData.addData("Input", 33, 9, SISid);
            outputData.addData("Input", 35, 9, SIAmSilica);
            outputData.addData("Input", 37, 9, SIQuartz);
            outputData.addData("Input", 39, 9, SIChrysotile);
            outputData.addData("Input", 41, 9, SIDiopside);
            outputData.addData("Input", 43, 9, SIGreenalite);
        }
        else {
            // 输出饱和度比
            outputData.addData("Input", 13, 9, pow(10.0, SICal));
            outputData.addData("Input", 15, 9, pow(10.0, SIBar));
            outputData.addData("Input", 17, 9, pow(10.0, SIHal));
            outputData.addData("Input", 19, 9, pow(10.0, SIGyp));
            outputData.addData("Input", 21, 9, pow(10.0, SIHemi));
            outputData.addData("Input", 23, 9, pow(10.0, SIAn));
            outputData.addData("Input", 25, 9, pow(10.0, SICel));
            outputData.addData("Input", 27, 9, pow(10.0, SIFeS));
            outputData.addData("Input", 29, 9, pow(10.0, SIZnS));
            outputData.addData("Input", 31, 9, pow(10.0, SICaF2));
            outputData.addData("Input", 33, 9, pow(10.0, SISid));
            outputData.addData("Input", 35, 9, pow(10.0, SIAmSilica));
            outputData.addData("Input", 37, 9, pow(10.0, SIQuartz));
            outputData.addData("Input", 39, 9, pow(10.0, SIChrysotile));
            outputData.addData("Input", 41, 9, pow(10.0, SIDiopside));
            outputData.addData("Input", 43, 9, pow(10.0, SIGreenalite));
        }

        // 输出其他参数到Input工作表
        outputData.addData("Input", 45, 9, pH_before_precipitation - DpHj);
        outputData.addData("Input", 48, 9, ConcInhCal);
        outputData.addData("Input", 50, 9, ConcInhBar);
        outputData.addData("Input", 52, 9, ConcInhGyp);
        outputData.addData("Input", 54, 9, static_cast<double>(ConcInhAn));
        outputData.addData("Input", 56, 9, ConcInhCel);
        outputData.addData("Input", 58, 9, ViscWatIst);
        outputData.addData("Input", 60, 9, CpPerMl);
        outputData.addData("Input", 62, 9, rhoTP);
    }

    // 第二部分：LoopTP = 10 或 LoopTP = 2（根据RunGoalSeek）时的输出
    if (RunGoalSeek != 1)
    {
        if (LoopTP == 10)
        {
            // 类似LoopTP=1的处理，但输出到第10列
            if (simContext.useSR == 0)
            {
                outputData.addData("Input", 13, 10, SICal);
                outputData.addData("Input", 15, 10, SIBar);
                outputData.addData("Input", 19, 10, SIGyp);
                outputData.addData("Input", 21, 10, SIHemi);
                outputData.addData("Input", 23, 10, SIAn);
                outputData.addData("Input", 25, 10, SICel);
                outputData.addData("Input", 17, 10, SIHal);
                outputData.addData("Input", 31, 10, SICaF2);
                outputData.addData("Input", 33, 10, SISid);
                outputData.addData("Input", 27, 10, SIFeS);
                outputData.addData("Input", 29, 10, SIZnS);
                outputData.addData("Input", 35, 10, SIAmSilica);
                outputData.addData("Input", 37, 10, SIQuartz);
                outputData.addData("Input", 39, 10, SIChrysotile);
                outputData.addData("Input", 41, 10, SIDiopside);
                outputData.addData("Input", 43, 10, SIGreenalite);
            }
            else {
                // 输出饱和度比到第10列
                outputData.addData("Input", 13, 10, pow(10.0, SICal));
                outputData.addData("Input", 15, 10, pow(10.0, SIBar));
                outputData.addData("Input", 19, 10, pow(10.0, SIGyp));
                outputData.addData("Input", 21, 10, pow(10.0, SIHemi));
                outputData.addData("Input", 23, 10, pow(10.0, SIAn));
                outputData.addData("Input", 25, 10, pow(10.0, SICel));
                outputData.addData("Input", 17, 10, pow(10.0, SIHal));
                outputData.addData("Input", 31, 10, pow(10.0, SICaF2));
                outputData.addData("Input", 33, 10, pow(10.0, SISid));
                outputData.addData("Input", 27, 10, pow(10.0, SIFeS));
                outputData.addData("Input", 29, 10, pow(10.0, SIZnS));
                outputData.addData("Input", 35, 10, pow(10.0, SIAmSilica));
                outputData.addData("Input", 37, 10, pow(10.0, SIQuartz));
                outputData.addData("Input", 39, 10, pow(10.0, SIChrysotile));
                outputData.addData("Input", 41, 10, pow(10.0, SIDiopside));
                outputData.addData("Input", 43, 10, pow(10.0, SIGreenalite));
            }

            outputData.addData("Input", 45, 10, pH_before_precipitation - DpHj);
            outputData.addData("Input", 48, 10, ConcInhCal);
            outputData.addData("Input", 50, 10, ConcInhBar);
            outputData.addData("Input", 52, 10, ConcInhGyp);
            outputData.addData("Input", 54, 10, static_cast<double>(ConcInhAn));
            outputData.addData("Input", 56, 10, ConcInhCel);
            outputData.addData("Input", 58, 10, ViscWatIst);
            outputData.addData("Input", 60, 10, CpPerMl);
            outputData.addData("Input", 62, 10, rhoTP);
        }
    }
    else {
        // RunGoalSeek = 1 的情况
        if (LoopTP == 2) {
            // 类似处理，输出到第10列
            if (simContext.useSR == 0) {
                outputData.addData("Input", 13, 10, SICal);
                outputData.addData("Input", 15, 10, SIBar);
                outputData.addData("Input", 19, 10, SIGyp);
                outputData.addData("Input", 21, 10, SIHemi);
                outputData.addData("Input", 23, 10, SIAn);
                outputData.addData("Input", 25, 10, SICel);
                outputData.addData("Input", 17, 10, SIHal);
                outputData.addData("Input", 31, 10, SICaF2);
                outputData.addData("Input", 33, 10, SISid);
                outputData.addData("Input", 27, 10, SIFeS);
                outputData.addData("Input", 29, 10, SIZnS);
                outputData.addData("Input", 35, 10, SIAmSilica);
                outputData.addData("Input", 37, 10, SIQuartz);
                outputData.addData("Input", 39, 10, SIChrysotile);
                outputData.addData("Input", 41, 10, SIDiopside);
                outputData.addData("Input", 43, 10, SIGreenalite);
            }
            else {
                // 输出饱和度比
                outputData.addData("Input", 13, 10, pow(10.0, SICal));
                outputData.addData("Input", 15, 10, pow(10.0, SIBar));
                outputData.addData("Input", 19, 10, pow(10.0, SIGyp));
                outputData.addData("Input", 21, 10, pow(10.0, SIHemi));
                outputData.addData("Input", 23, 10, pow(10.0, SIAn));
                outputData.addData("Input", 25, 10, pow(10.0, SICel));
                outputData.addData("Input", 17, 10, pow(10.0, SIHal));
                outputData.addData("Input", 31, 10, pow(10.0, SICaF2));
                outputData.addData("Input", 33, 10, pow(10.0, SISid));
                outputData.addData("Input", 27, 10, pow(10.0, SIFeS));
                outputData.addData("Input", 29, 10, pow(10.0, SIZnS));
                outputData.addData("Input", 35, 10, pow(10.0, SIAmSilica));
                outputData.addData("Input", 37, 10, pow(10.0, SIQuartz));
                outputData.addData("Input", 39, 10, pow(10.0, SIChrysotile));
                outputData.addData("Input", 41, 10, pow(10.0, SIDiopside));
                outputData.addData("Input", 43, 10, pow(10.0, SIGreenalite));
            }

            outputData.addData("Input", 45, 10, pH_before_precipitation - DpHj);
            outputData.addData("Input", 48, 10, ConcInhCal);
            outputData.addData("Input", 50, 10, ConcInhBar);
            outputData.addData("Input", 52, 10, ConcInhGyp);
            outputData.addData("Input", 54, 10, static_cast<double>(ConcInhAn));
            outputData.addData("Input", 56, 10, ConcInhCel);
            outputData.addData("Input", 58, 10, ViscWatIst);
            outputData.addData("Input", 60, 10, CpPerMl);
            outputData.addData("Input", 62, 10, rhoTP);
        }
    }

    // ==================== 表头和其他信息输出 ====================
    if (RunGoalSeek != 1) {
        // 设置表头
        outputData.addData("calcite", 3, 20, "SI(Inh=" + std::to_string(simContext.MaxInh / 2) + "mg/L)");
        outputData.addData("Calcite", 3, 21, "SI (Inh=" + std::to_string(simContext.MaxInh) + "mg/L)");
        outputData.addData("Barite", 3, 20, "SI (Inh=" + std::to_string(simContext.MaxInh / 2) + "mg/L)");
        outputData.addData("Barite", 3, 21, "SI (Inh=" + std::to_string(simContext.MaxInh) + "mg/L)");
        outputData.addData("Other SO4s", 3, 20, "SI (Inh=" + std::to_string(simContext.MaxInh / 2) + "mg/L)");
        outputData.addData("Other SO4s", 3, 21, "SI (Inh=" + std::to_string(simContext.MaxInh) + "mg/L)");
        outputData.addData("Other SO4s", 15, 20, "SI (Inh=" + std::to_string(simContext.MaxInh / 2) + "mg/L)");
        outputData.addData("Other SO4s", 15, 21, "SI (Inh=" + std::to_string(simContext.MaxInh) + "mg/L)");

        // 输出井名
        int nob1 = simContext.nob;
        if (simContext.nob > 2) nob1 = 2; // 只输出2个井名

        for (int i = 0; i < nob1 && i < WellNameMix.size(); i++) {
            outputData.addData("Calcite", 1, i + 2, WellNameMix[i]);
            outputData.addData("Barite", 1, i + 2, WellNameMix[i]);
            outputData.addData("Halite", 1, i + 2, WellNameMix[i]);
            outputData.addData("Other SO4s", 1, i + 2, WellNameMix[i]);
            outputData.addData("Sulfides,Fluorite,Carbonates", 1, i + 2, WellNameMix[i]);
            outputData.addData("Silicates", 1, i + 2, WellNameMix[i]);
            outputData.addData("Mg(OH)2,Ca(OH)2", 1, i + 2, WellNameMix[i]);
        }

        // 输出MeOH和MEG数据
        if (simContext.UseSI == 0) {
            outputData.addData("Calcite", 1, 15, mass_MeOH / 0.7914 / 159.0);
            outputData.addData("Calcite", 1, 12, mass_MEG / 1.1135 / 159.0);
            outputData.addData("Barite", 2, 13, mass_MeOH / 0.7914 / 159.0);
            outputData.addData("Barite", 3, 13, mass_MEG / 1.1135 / 159.0);
            outputData.addData("Halite", 2, 13, mass_MeOH / 0.7914 / 159.0);
            outputData.addData("Halite", 3, 13, mass_MEG / 1.1135 / 159.0);
            outputData.addData("Other SO4s", 2, 14, mass_MeOH / 0.7914 / 159.0);
            outputData.addData("Other SO4s", 3, 14, mass_MEG / 1.1135 / 159.0);
        }
        else {
            outputData.addData("Calcite", 1, 15, mass_MeOH / 0.7914 / 159.0 * 0.159);
            outputData.addData("Calcite", 1, 10, mass_MEG / 1.1135 / 159.0 * 0.159);
            outputData.addData("Barite", 2, 13, mass_MeOH / 0.7914 / 159.0 * 0.159);
            outputData.addData("Barite", 3, 13, mass_MEG / 1.1135 / 159.0 * 0.159);
            outputData.addData("Halite", 2, 13, mass_MeOH / 0.7914 / 159.0 * 0.159);
            outputData.addData("Halite", 3, 13, mass_MEG / 1.1135 / 159.0 * 0.159);
            outputData.addData("Other SO4s", 2, 14, mass_MeOH / 0.7914 / 159.0 * 0.159);
            outputData.addData("Other SO4s", 3, 14, mass_MEG / 1.1135 / 159.0 * 0.159);
        }

        // ==================== 活度系数输出 ====================
        if (simContext.OutPutActCoefs == 1) {
            // Calcite 工作表中的活度系数
            if (simContext.UseSI == 0) {
                outputData.addData("Calcite", 14 - LoopTP, 24, TF);
                outputData.addData("Calcite", 14 - LoopTP, 25, Ppsia);
            }
            else {
                outputData.addData("Calcite", 14 - LoopTP, 24, (TF - 32) * 5.0 / 9.0);
                outputData.addData("Calcite", 14 - LoopTP, 25, Ppsia / 14.503774);
            }
            // 输出各种活度系数
            if (gCat.size() > 0 && gAn.size() > 0 && gNeut.size() > 0) {
                // 需要调整索引值以匹配C++的0基索引
                // H+ 活度系数 iH
                if (gCat.size() > 0) outputData.addData("Calcite", 14 - LoopTP, 27, gCat[0]);
                // OH- 活度系数 iOH
                if (gAn.size() > 0) outputData.addData("Calcite", 14 - LoopTP, 28, gAn[0]);
                // Ca2+ 活度系数 iCa
                if (gCat.size() > 1) outputData.addData("Calcite", 14 - LoopTP, 29, gCat[1]);
                // Mg2+ 活度系数 iMg
                if (gCat.size() > 2) outputData.addData("Calcite", 14 - LoopTP, 30, gCat[2]);
                // HCO3- 活度系数 iHCO3
                if (gAn.size() > 1) outputData.addData("Calcite", 14 - LoopTP, 31, gAn[1]);
                // CO3^2- 活度系数 iCO3
                if (gAn.size() > 2) outputData.addData("Calcite", 14 - LoopTP, 32, gAn[2]);
                // Ac- 活度系数 iAc
                if (gAn.size() > 3) outputData.addData("Calcite", 14 - LoopTP, 33, gAn[3]);

                // 水的活度
                outputData.addData("Calcite", 14 - LoopTP, 26, aH2O);

                // 其他工作表的活度系数输出，类似处理
                // Worksheets("Barite").Cells(14 - LoopTP, 23) = TF
                // Worksheets("Barite").Cells(14 - LoopTP, 24) = Ppsia
                // If UseSI = 1 Then Worksheets("Barite").Cells(14 - LoopTP, 23) = (TF - 32) * 5 / 9
                // If UseSI = 1 Then Worksheets("Barite").Cells(14 - LoopTP, 24) = Ppsia / 14.503774
                // Worksheets("Barite").Cells(14 - LoopTP, 25) = gCat(iBa)
                // Worksheets("Barite").Cells(14 - LoopTP, 26) = gAn(iSO4)
                // Worksheets("Other SO4s").Cells(14 - LoopTP, 24) = TF
                // Worksheets("Other SO4s").Cells(14 - LoopTP, 25) = Ppsia
                // If UseSI = 1 Then Worksheets("Other SO4s").Cells(14 - LoopTP, 24) = (TF - 32) * 5 / 9
                // If UseSI = 1 Then Worksheets("Other SO4s").Cells(14 - LoopTP, 25) = Ppsia / 14.503774
                // Worksheets("Other SO4s").Cells(14 - LoopTP, 26) = gCat(iSr)
                // Worksheets("Other SO4s").Cells(14 - LoopTP, 27) = gAn(iSO4)
                // Worksheets("Halite").Cells(14 - LoopTP, 15) = TF
                // Worksheets("Halite").Cells(14 - LoopTP, 16) = Ppsia
                // If UseSI = 1 Then Worksheets("Halite").Cells(14 - LoopTP, 15) = (TF - 32) * 5 / 9
                // If UseSI = 1 Then Worksheets("Halite").Cells(14 - LoopTP, 16) = Ppsia / 14.503774
                // Worksheets("Halite").Cells(14 - LoopTP, 17) = gCat(iNa)
                // Worksheets("Halite").Cells(14 - LoopTP, 18) = gAn(iCl)
                // Worksheets("Sulfides,Fluorite,Carbonates").Cells(14 - LoopTP, 26) = TF
                // Worksheets("Sulfides,Fluorite,Carbonates").Cells(14 - LoopTP, 27) = Ppsia
                // If UseSI = 1 Then Worksheets("Sulfides,Fluorite,Carbonates").Cells(14 - LoopTP, 26) = (TF - 32) * 5 / 9
                // If UseSI = 1 Then Worksheets("Sulfides,Fluorite,Carbonates").Cells(14 - LoopTP, 27) = Ppsia / 14.503774
                // Worksheets("Sulfides,Fluorite,Carbonates").Cells(14 - LoopTP, 28) = gCat(iFe)
                // Worksheets("Sulfides,Fluorite,Carbonates").Cells(14 - LoopTP, 29) = gCat(iZn)
                // Worksheets("Sulfides,Fluorite,Carbonates").Cells(14 - LoopTP, 30) = gAn(intF)
                // Worksheets("Sulfides,Fluorite,Carbonates").Cells(14 - LoopTP, 31) = gAn(iHS)
                // Worksheets("Sulfides,Fluorite,Carbonates").Cells(14 - LoopTP, 32) = gCat(iPb)

                // Worksheets("Silicates").Cells(14 - LoopTP, 24) = TF
                // Worksheets("Silicates").Cells(14 - LoopTP, 25) = Ppsia
                // If UseSI = 1 Then Worksheets("Silicates").Cells(14 - LoopTP, 24) = (TF - 32) * 5 / 9
                // If UseSI = 1 Then Worksheets("Silicates").Cells(14 - LoopTP, 25) = Ppsia / 14.503774
                // Worksheets("Silicates").Cells(14 - LoopTP, 26) = gNeut(iH4SiO4aq)
            }
        }
    }

    // ==================== 将数据写入Excel（如果有excelWriter）====================
    if (excelWriter != nullptr) {
        // 遍历所有工作表和数据，调用excelWriter的方法写入
        for (auto& sheetPair : outputData.worksheets) {
            const std::string& sheetName = sheetPair.first;
            auto& worksheet = sheetPair.second;

            // 写入数值数据
            for (auto& dataPair : worksheet.numericData) {
                int row = dataPair.first.first;
                int col = dataPair.first.second;
                double value = dataPair.second;
                excelWriter->writeCell(sheetName, row, col, value);
            }

            // 写入文本数据
            for (auto& dataPair : worksheet.textData) {
                int row = dataPair.first.first;
                int col = dataPair.first.second;
                const std::string& value = dataPair.second;
                excelWriter->writeCell(sheetName, row, col, value);
            }

            // 写入布尔数据
            for (auto& dataPair : worksheet.booleanData) {
                int row = dataPair.first.first;
                int col = dataPair.first.second;
                bool value = dataPair.second;
                excelWriter->writeCell(sheetName, row, col, value);
            }
        }
    }
    else {
        // 如果不使用 excelWriter，就将数据保存到文件或数据库
        // 例如：保存为CSV文件或JSON格式
        std::cout << "数据已收集，但未写入Excel（未提供ExcelWriter实例）" << std::endl;
        std::cout << "共收集了 " << outputData.worksheets.size() << " 个工作表的数据" << std::endl;
    }
}

void LoopTPppt() {
    // 初始化所有沉淀量为0
    pptCal = 0.0; pptBar = 0.0; pptSid = 0.0; pptGyp = 0.0; pptHemi = 0.0;
    pptAn = 0.0; pptCel = 0.0; pptHal = 0.0; pptFeS = 0.0; pptZnS = 0.0;
    pptZnCO3 = 0.0; pptPbS = 0.0; pptCaF2 = 0.0; pptMgOH2 = 0.0; pptCaOH2 = 0.0;
    pptFeSAm = 0.0; pptTrot = 0.0; pptAmSilica = 0.0; pptQuartz = 0.0;
    pptGreenalite = 0.0; pptDiopside = 0.0; pptChrysotile = 0.0;
    pptFeCO3_NoMassTransfer = 0.0; pptFeCO3_MassTransfer = 0.0;
    pptCalcite_MassTransfer = 0.0; pptCalcite_NoMassTransfer = 0.0;
    pptBarite_MassTransfer = 0.0; pptBarite_NoMassTransfer = 0.0;
    pptCelestite_MassTransfer = 0.0; pptCelestite_NoMassTransfer = 0.0;
    pptBarite_MassTransfer_V = 0.0; pptCalcite_MassTransfer_V = 0.0;
    pptFeCO3_MassTransfer_V = 0.0;
    pptFeS_MassTransfer = 0.0; pptFeS_NoMassTransfer = 0.0;

    // 初始化pH值
    pHaftercalciteppt = pH_before_precipitation;
    pHafterMgOH2ppt = pH_before_precipitation;
    pHafterCaOH2ppt = pH_before_precipitation;
    pHafterAmsilicappt = pH_before_precipitation;
    pHafterQuartzppt = pH_before_precipitation;
    pHafterGreenaliteppt = pH_before_precipitation;
    pHafterDiopsideppt = pH_before_precipitation;
    pHafterChrysotileppt = pH_before_precipitation;

    // 如果是第一次循环且有传质，初始化传质浓度
    if (LoopTP == 1 && Run_MassTransfer == 1) {
        // 确保mcMT和maMT有足够大小
        if (mcMT.size() < mc.size()) mcMT.resize(mc.size());
        if (maMT.size() < ma.size()) maMT.resize(ma.size());
        
        mcMT[iCa] = mc[iCa];
        mcMT[iBa] = mc[iBa];
        maMT[iSO4] = ma[iSO4];
        AlkMT = 0.0; // 碱度需要从其他地方获取
        
        mcMT[iFe] = mc[iFe];
        TFe_MT = 0.0; // 总铁需要从其他地方获取
    }

    // ==================== Barite (BaSO4) 沉淀计算 ====================
    double ppt = 0.0;
    if (SIBar > 0) {
        double a = 1.0;
        double b = -(mc[iBa] + ma[iSO4]);
        double c = mc[iBa] * ma[iSO4] - KspBarite / (gCat[iBa] * gAn[iSO4] * std::pow(gNMean[iBaSO4], 2));
        double root2 = 0.0;
        QuadraticRoots(a, b, c, ppt, root2);
        
        // 计算平衡浓度
        if (mcEq.size() <= iBa) mcEq.resize(iBa + 1);
        mcEq[iBa] = mc[iBa] - ppt;
        
        // 镭共沉淀计算
        RaBaSO4ppt = mc[iRa] / (1.0 + ((mc[iBa] - ppt) * gCat[iBa] * GammaSolidRaBarite / 
                                      (ppt * DRaBarite * gCat[iRa])));
    }
    pptBarite_NoMassTransfer = ppt;
    
    // Barite 传质计算
    if (SIBar > 0 && Run_MassTransfer == 1) {
        // 确保km数组有足够大小
        if (km.size() > 0) {
            double km_BaSO4 = km[0]; // 假设索引0对应BaSO4的传质系数
            
            // 计算动力学浓度
            mcKin[iBa] = (std::exp(-km_BaSO4 * SArea / QBrineFlow) * (mc[iBa] - mcEq[iBa]) + mcEq[iBa]);
            pptBarite_MassTransfer = mc[iBa] - mcKin[iBa];
            
            // 更新传质浓度
            mcMT[iBa] = mcKin[iBa];
            maMT[iSO4] = ma[iSO4] - pptBarite_MassTransfer;
            
            // Dai 2020 deposition - 计算沉积速率
            double density_factor = simContext.rho25c - TDS * 0.000001;
            pptBarite_MassTransfer_V = pptBarite_MassTransfer * 233390.0 * density_factor * 
                                     QBrineFlow / 1000.0 * 3600.0 * 24.0 * 365.0 / 
                                     (4480.0 * pi * PipeID * PipeL); // cm/yr
        }
    }

    // ==================== Gypsum (CaSO4·2H2O) 沉淀计算 ====================
    ppt = 0.0;
    if (SIGyp > 0) {
        double a = 1.0;
        double b = -(mc[iCa] + ma[iSO4]);
        double c = mc[iCa] * ma[iSO4] - KspGypsum / 
                   (std::pow(aH2O, 2) * gCat[iCa] * gAn[iSO4] * std::pow(aNH2O, 2) * std::pow(gNMean[iCa], 2));
        
        double root2 = 0.0;
        QuadraticRoots(a, b, c, ppt, root2);
    }
    pptGyp = ppt;

    // ==================== Hemihydrate (CaSO4·0.5H2O) 沉淀计算 ====================
    ppt = 0.0;
    if (SIHemi > 0) {
        double a = 1.0;
        double b = -(mc[iCa] + ma[iSO4]);
        double c = mc[iCa] * ma[iSO4] - KspHemihydrate / 
                   (std::pow(aH2O, 0.5) * gCat[iCa] * gAn[iSO4] * std::pow(aNH2O, 0.5) * std::pow(gNMean[iCa], 2));
        
        double root2 = 0.0;
        QuadraticRoots(a, b, c, ppt, root2);
    }
    pptHemi = ppt;

    // ==================== Anhydrite (CaSO4) 沉淀计算 ====================
    ppt = 0.0;
    if (SIAn > 0) {
        double a = 1.0;
        double b = -(mc[iCa] + ma[iSO4]);
        double c = mc[iCa] * ma[iSO4] - KspAnhydrite / 
                   (gCat[iCa] * gAn[iSO4] * std::pow(gNMean[iCa], 2));
        
        double root2 = 0.0;
        QuadraticRoots(a, b, c, ppt, root2);
        
        // 镭共沉淀计算
        RaCaSO4ppt = mc[iRa] / (1.0 + ((mc[iCa] - ppt) * gCat[iCa] * GammaSolidRaAnhydrite / 
                                      (ppt * DRaAnhydrite * gCat[iRa])));
    }
    pptAn = ppt;

    // ==================== Celestite (SrSO4) 沉淀计算 ====================
    ppt = 0.0;
    if (SICel > 0) {
        double a = 1.0;
        double b = -(mc[iSr] + ma[iSO4]);
        double c = mc[iSr] * ma[iSO4] - KspCelestite / 
                   (gCat[iSr] * gAn[iSO4] * std::pow(gNMean[iSr], 2));
        
        double root2 = 0.0;
        QuadraticRoots(a, b, c, ppt, root2);
        
        // 计算平衡浓度
        if (mcEq.size() <= iSr) mcEq.resize(iSr + 1);
        mcEq[iSr] = mc[iSr] - ppt;
        
        // 镭共沉淀计算
        RaSrSO4ppt = mc[iRa] / (1.0 + ((mc[iSr] - ppt) * gCat[iSr] * GammaSolidRaCelestite / 
                                      (ppt * DRaCelestite * gCat[iRa])));
    }
    pptCel = ppt;
    pptCelestite_NoMassTransfer = ppt;

    // ==================== Halite (NaCl) 沉淀计算 ====================
    ppt = 0.0;
    if (SIHal > 0) {
        // 保存原始值
        double TNa = mc[iNa];
        double TCl = ma[iCl];
        double TKTemp = TK;
        double PBarTemp = PBar;
        double TCTemp = TC;
        double PatmTemp = Patm;
        
        // 第一次计算
        double a = 1.0;
        double b = -(mc[iNa] + ma[iCl]);
        double c = mc[iNa] * ma[iCl] - KspHalite / 
                   (gCat[iNa] * gAn[iCl] * std::pow(gNMean[iNa], 2));
        
        double root2 = 0.0;
        QuadraticRoots(a, b, c, ppt, root2);
        
        // 更新浓度
        mc[iNa] = mc[iNa] - ppt;
        ma[iCl] = ma[iCl] - ppt;
        
        // 迭代计算直到收敛
        int iteration = 0;
        const int maxIterations = 100;
        
        while (std::abs(ppt) >= 0.001 && iteration < maxIterations) {
            // 重新计算离子强度和活度系数
            CalcIonicStrength();
            C1_ThermodynamicEquilConsts();
            C2_PitzerActCoefs_T_P_ISt(gNeut, aH2O, TK, TC, PBar, Patm);
            
            // 重新计算沉淀量
            a = 1.0;
            b = -(mc[iNa] + ma[iCl]);
            c = mc[iNa] * ma[iCl] - KspHalite / 
                (gCat[iNa] * gAn[iCl] * std::pow(gNMean[iNa], 2));
            
            QuadraticRoots(a, b, c, ppt, root2);
            
            // 更新浓度
            mc[iNa] = mc[iNa] - ppt;
            ma[iCl] = ma[iCl] - ppt;
            
            iteration++;
        }
        
        // 计算实际沉淀量
        pptHal = TNa - mc[iNa];
        
        // 在标准条件下计算TDS和密度
        TK = 298.15;
        TC = 25.0;
        PBar = 1.01325;
        Patm = 1.0;
        
        CalcIonicStrength();
        C2_PitzerActCoefs_T_P_ISt(gNeut, aH2O, TK, TC, PBar, Patm);
        
        // 计算密度
        Rho25cHalite = CalcRhoTP(TK, TC, PBar, Patm);
        
        // 计算TDS（简化版本）
        TDSHalite = 0.0;
        double CalculateTDSDen = 0.0;
        
        // 阳离子贡献
        for (size_t i = 1; i < mc.size(); i++) {
            CalculateTDSDen += 0.001 * mc[i] * MWCat[i];
        }
        
        // 阴离子贡献
        for (size_t i = 1; i < ma.size(); i++) {
            CalculateTDSDen += 0.001 * ma[i] * MWAn[i];
        }
        
        // 中性物质贡献
        for (size_t i = 1; i < mn.size(); i++) {
            CalculateTDSDen += 0.001 * mn[i] * MWNeut[i];
        }
        
        // 恢复原始值
        mc[iNa] = TNa;
        ma[iCl] = TCl;
        TK = TKTemp;
        TC = TCTemp;
        PBar = PBarTemp;
        Patm = PatmTemp;
        Ppsia = PBar * 14.503774;
        
        // 重新计算活度系数
        CalcIonicStrength();
        C1_ThermodynamicEquilConsts();
        C2_PitzerActCoefs_T_P_ISt(gNeut, aH2O, TK, TC, PBar, Patm);
    } else {
        TDSHalite = TDS;
        Rho25cHalite = simContext.rho25c;
    }

    // ==================== Fluorite (CaF2) 沉淀计算 ====================
    ppt = 0.0;
    if (SICaF2 > 0) {
        // 立方方程系数
        double coef1 = 1.0;
        double coef2 = -1.0 * (mc[iCa] + ma[intF]);
        double coef3 = (4.0 * mc[iCa] * ma[intF] + std::pow(ma[intF], 2)) / 4.0;
        double coef4 = -1.0 * (mc[iCa] * std::pow(ma[intF], 2) - 
                             KspCaF2 / (gCat[iCa] * std::pow(gAn[intF], 2))) / 4.0;
        
        double root1 = 0.0, root2 = 0.0, root3 = 0.0;
        CubicRoots(coef1, coef2, coef3, coef4, root1, root2, root3);
        
        // 选择最小的正根
        ppt = root1;
        if (root2 != 0.0 && root2 < ppt) ppt = root2;
        if (root3 != 0.0 && root3 < ppt) ppt = root3;
    }
    pptCaF2 = ppt;

    // ==================== Calcite (CaCO3) 沉淀计算 ====================
    ppt = 0.0;
    pH = pH_before_precipitation;
    
    if (SICal > 0) {
        // 调用SSP平衡计算 vb代码中注释了
        // if (useEOS == 0) {
            C4_SSPEquilCalcs(1, 5, 2, KspCalcite);  // ppt_or_not=1, iMetals=5(Ca), iGas=2(CO2)
        // } else {
        //     C4_EOS_TCO2_SSPEquilCalcs(1, 5, 2, KspCalcite);
        // }
        
        pptCalcite_NoMassTransfer = ppt;
        
        // 计算平衡浓度
        mcEq[iCa] = mc[iCa] - ppt;
    }
    pHaftercalciteppt = pH;

    // Calcite 传质计算
    if (SICal > 0 && Run_MassTransfer == 1) {
        // 假设km[1]对应CaHCO3+的传质系数
        if (km.size() > 1) {
            double km_CaHCO3 = km[1];
            
            // 计算动力学浓度
            mcKin[iCa] = std::exp(-km_CaHCO3 * SArea / QBrineFlow) * (mc[iCa] - mcEq[iCa]) + mcEq[iCa];
            pptCalcite_MassTransfer = mc[iCa] - mcKin[iCa];
            
            // 更新传质浓度
            mcMT[iCa] = mcKin[iCa];
            AlkMT = AlkMT - 2.0 * pptCalcite_MassTransfer;  // 假设AlkMT已初始化
            
            // Dai 2020 deposition - 计算沉积速率
            double density_factor = simContext.rho25c - TDS * 0.000001;
            pptCalcite_MassTransfer_V = pptCalcite_MassTransfer * 100091.0 * density_factor * 
                                      QBrineFlow / 1000.0 * 3600.0 * 24.0 * 365.0 / 
                                      (2710.0 * pi * PipeID * PipeL); // cm/yr
        }
    }

    // ==================== Siderite (FeCO3) 沉淀计算 ====================
    ppt = 0.0;
    pH = pH_before_precipitation;
    
    if (SISid > 0) {
        // 如果有方解石沉淀，使用方解石沉淀后的pH
        // if (SICal > 0) pH = pHaftercalciteppt;
        
        // 调用SSP平衡计算
        // if (useEOS == 0) {
            C4_SSPEquilCalcs(1, 8, 2, KspSiderite);  // ppt_or_not=1, iMetals=8(Fe), iGas=2(CO2)
        // } else {
        //     C4_EOS_TCO2_SSPEquilCalcs(1, 8, 2, KspSiderite);
        // }
        
        pptFeCO3_NoMassTransfer = ppt;
        
        // 计算平衡浓度
        mcEq[iFe] = mc[iFe] - ppt;
    }
    pHafterFeCO3ppt = pH;

    // Siderite 传质计算
    if (SISid > 0 && Run_MassTransfer == 1) {
        // 假设km[2]对应FeHCO3+的传质系数
        if (km.size() > 2) {
            double km_FeHCO3 = km[2];
            
            // 计算动力学浓度
            mcKin[iFe] = std::exp(-km_FeHCO3 * SArea / QBrineFlow) * (mc[iFe] - mcEq[iFe]) + mcEq[iFe];
            pptFeCO3_MassTransfer = mc[iFe] - mcKin[iFe];
            
            // 更新传质浓度
            mcMT[iFe] = mcKin[iFe];
            AlkMT = AlkMT - 2.0 * pptFeCO3_MassTransfer;
            
            // Dai 2021 corrosion deposition - 计算沉积速率
            double density_factor = simContext.rho25c - TDS * 0.000001;
            pptFeCO3_MassTransfer_V = pptFeCO3_MassTransfer * 100091.0 * density_factor * 
                                     QBrineFlow / 1000.0 * 3600.0 * 24.0 * 365.0 / 
                                     (3960.0 * pi * PipeID * PipeL); // cm/yr
        }
    }

    // ==================== 计算总沉积速率 ====================
    // Dai 2020 deposition
    pptTotal_MassTransfer_V = pptBarite_MassTransfer_V + pptCalcite_MassTransfer_V;
    if (pptTotal_MassTransfer_V > pptMax_MassTransfer_V) {
        pptMax_MassTransfer_V = pptTotal_MassTransfer_V;
    }
    
    // Dai 2021 corrosion deposition
    pptTotal_MassTransfer_V = pptTotal_MassTransfer_V + pptFeCO3_MassTransfer_V;
    if (pptTotal_MassTransfer_V > pptMax_MassTransfer_V) {
        pptMax_MassTransfer_V = pptTotal_MassTransfer_V;
    }

    // ==================== Zinc Carbonate (ZnCO3) 沉淀计算 ====================
    ppt = 0.0;
    pH = pH_before_precipitation;
    
    if (SIZnCO3 > 0) {
        // if (useEOS == 0) {
            C4_SSPEquilCalcs(1, 9, 2, KspZnCO3);  // ppt_or_not=1, iMetals=9(Zn), iGas=2(CO2)
        // } else {
        //     C4_EOS_TCO2_SSPEquilCalcs(1, 9, 2, KspZnCO3);
        // }
        pptZnCO3 = ppt;
    }

    // ==================== Amorphous Iron Sulfide (FeS) 沉淀计算 ====================
    ppt = 0.0;
    pH = pH_before_precipitation;
    
    if (SIFeSAm > 0) {
        // if (useEOS == 0) {
            C4_SSPEquilCalcs(1, 8, 3, KspFeSAm);  // ppt_or_not=1, iMetals=8(Fe), iGas=3(H2S)
        // } else {
        //     C4_EOS_TCO2_SSPEquilCalcs(1, 8, 3, KspFeSAm);
        // }
        pptFeSAm = ppt;
    }

    // ==================== Mackinawite (FeS) 沉淀计算 ====================
    ppt = 0.0;
    pH = pH_before_precipitation;
    
    if (SIFeS > 0) {
        // if (useEOS == 0) {
            C4_SSPEquilCalcs(1, 8, 3, KspFeS);  // ppt_or_not=1, iMetals=8(Fe), iGas=3(H2S)
        // } else {
        //     C4_EOS_TCO2_SSPEquilCalcs(1, 8, 3, KspFeS);
        // }
        pptFeS_NoMassTransfer = ppt;
    }

    // ==================== Troilite (FeS) 沉淀计算 ====================
    ppt = 0.0;
    pH = pH_before_precipitation;
    
    if (SITrot > 0) {
        // if (useEOS == 0) {
            C4_SSPEquilCalcs(1, 8, 3, KspTrot);  // ppt_or_not=1, iMetals=8(Fe), iGas=3(H2S)
        // } else {
        //     C4_EOS_TCO2_SSPEquilCalcs(1, 8, 3, KspTrot);
        // }
        pptTrot = ppt;
    }

    // ==================== Pyrrotite (FeS) 沉淀计算 ====================
    ppt = 0.0;
    pH = pH_before_precipitation;
    
    if (simContext.SIPyrr > 0) {
        // if (useEOS == 0) {
            C4_SSPEquilCalcs(1, 8, 3, KspPyrr);  // ppt_or_not=1, iMetals=8(Fe), iGas=3(H2S)
        // } else {
        //     C4_EOS_TCO2_SSPEquilCalcs(1, 8, 3, KspPyrr);
        // }
        // pptPyrr = ppt;  // 变量未在顶部声明
    }

    // ==================== Zinc Sulfide (ZnS) 沉淀计算 ====================
    ppt = 0.0;
    pH = pH_before_precipitation;
    
    if (SIZnS > 0) {
        // if (useEOS == 0) {
            C4_SSPEquilCalcs(1, 9, 3, KspZnS);  // ppt_or_not=1, iMetals=9(Zn), iGas=3(H2S)
        // } else {
        //     C4_EOS_TCO2_SSPEquilCalcs(1, 9, 3, KspZnS);
        // }
        pptZnS = ppt;
    }

    // ==================== Lead Sulfide (PbS) 沉淀计算 ====================
    ppt = 0.0;
    pH = pH_before_precipitation;
    
    if (SIPbS > 0) {
        // if (useEOS == 0) {
            C4_SSPEquilCalcs(1, 10, 3, KspPbS);  // ppt_or_not=1, iMetals=10(Pb), iGas=3(H2S)
        // } else {
        //     C4_EOS_TCO2_SSPEquilCalcs(1, 10, 3, KspPbS);
        // }
        pptPbS = ppt;
    }

    // ==================== Calcium Hydroxide (Ca(OH)2) 沉淀计算 ====================
    ppt = 0.0;
    pH = pH_before_precipitation;
    
    if (SICaOH2 > 0) {
        // if (useEOS == 0) {
            C4_SSPEquilCalcs(1, 5, 8, KspCaOH2);  // ppt_or_not=1, iMetals=5(Ca), iGas=8(信号不改变气体)
        // } else {
        //     C4_EOS_TCO2_SSPEquilCalcs(1, 5, 8, KspCaOH2);
        // }
        pptCaOH2 = ppt;
    }
    pHafterCaOH2ppt = pH;

    // ==================== Magnesium Hydroxide (Mg(OH)2) 沉淀计算 ====================
    ppt = 0.0;
    pH = pH_before_precipitation;
    
    if (SIMgOH2 > 0) {
        // if (useEOS == 0) {
            C4_SSPEquilCalcs(1, 4, 9, KspMgOH2);  // ppt_or_not=1, iMetals=4(Mg), iGas=9(信号不改变气体)
        // } else {
        //     C4_EOS_TCO2_SSPEquilCalcs(1, 4, 9, KspMgOH2);
        // }
        pptMgOH2 = ppt;
    }
    pHafterMgOH2ppt = pH;

    // ==================== Amorphous Silica (SiO2) 沉淀计算 ====================
    ppt = 0.0;
    pH = pH_before_precipitation;
    
    if (SIAmSilica > 0) {
        // if (useEOS == 0) {
            C4_SSPEquilCalcs(1, 8, 10, KspAmSilica);  // iMetals=8(Fe), iGas=10(信号)
        // } else {
        //     C4_EOS_TCO2_SSPEquilCalcs(1, 8, 10, KspAmSilica);
        // }
        pptAmSilica = ppt;
    }
    pHafterAmsilicappt = pH;

    // ==================== Quartz (SiO2) 沉淀计算 ====================
    ppt = 0.0;
    pH = pH_before_precipitation;
    
    if (SIQuartz > 0) {
        // if (useEOS == 0) {
            C4_SSPEquilCalcs(1, 8, 11, KspQuartz);  // iMetals=8(Fe), iGas=11(信号)
        // } else {
        //     C4_EOS_TCO2_SSPEquilCalcs(1, 8, 11, KspQuartz);
        // }
        pptQuartz = ppt;
    }
    pHafterQuartzppt = pH;

    // ==================== Greenalite (Fe3Si2O5(OH)4) 沉淀计算 ====================
    ppt = 0.0;
    pH = pH_before_precipitation;
    
    if (SIGreenalite > 0) {
        // if (useEOS == 0) {
            C4_SSPEquilCalcs(1, 8, 12, KspGreenalite);  // iMetals=8(Fe), iGas=12(信号)
        // } else {
        //     C4_EOS_TCO2_SSPEquilCalcs(1, 8, 12, KspGreenalite);
        // }
        pptGreenalite = ppt;
    }
    pHafterGreenaliteppt = pH;

    // ==================== Diopside (CaMgSi2O6) 沉淀计算 ====================
    ppt = 0.0;
    pH = pH_before_precipitation;
    
    if (SIDiopside > 0) {
        // if (useEOS == 0) {
            C4_SSPEquilCalcs(1, 8, 13, KspDiopside);  // iMetals=8(Fe), iGas=13(信号)
        // } else {
        //     C4_EOS_TCO2_SSPEquilCalcs(1, 8, 13, KspDiopside);
        // }
        pptDiopside = ppt;
    }
    pHafterDiopsideppt = pH;

    // ==================== Chrysotile (Mg3Si2O5(OH)4) 沉淀计算 ====================
    ppt = 0.0;
    pH = pH_before_precipitation;
    
    if (SIChrysotile > 0) {
        // if (useEOS == 0) {
            C4_SSPEquilCalcs(1, 8, 14, KspChrysotile);  // iMetals=8(Fe), iGas=14(信号)
        // } else {
        //     C4_EOS_TCO2_SSPEquilCalcs(1, 8, 14, KspChrysotile);
        // }
        pptChrysotile = ppt;
    }
    pHafterChrysotileppt = pH;

    // ==================== Strontium Carbonate (SrCO3) 沉淀计算 ====================
    ppt = 0.0;
    pH = pH_before_precipitation;
    
    if (SISrCO3 > 0) {
        // if (useEOS == 0) {
            C4_SSPEquilCalcs(1, 6, 2, KspSrCO3);  // ppt_or_not=1, iMetals=6(Sr), iGas=2(CO2)
        // } else {
        //     C4_EOS_TCO2_SSPEquilCalcs(1, 6, 2, KspSrCO3);
        // }
        pptSrCO3 = ppt;
    }

    // ==================== Barium Carbonate (BaCO3) 沉淀计算 ====================
    ppt = 0.0;
    pH = pH_before_precipitation;
    
    if (SIBaCO3 > 0) {
        // if (useEOS == 0) {
            C4_SSPEquilCalcs(1, 6, 2, KspBaCO3);  // ppt_or_not=1, iMetals=6(Ba), iGas=2(CO2)
        // } else {
        //     C4_EOS_TCO2_SSPEquilCalcs(1, 6, 2, KspBaCO3);
        // }
        pptBaCO3 = ppt;
    }

    // ==================== 多重沉淀计算 ====================
    if (MultiplePpt == 1) {
        // Call LoopTPMultiple();  // 需要实现此函数
    }
}

int A1_Start_ScaleSoftPitzer(SampleData& data)
{
    RunSSP = 1;
    simContext.RunStat = 0;
    simContext.RunQualityControlChecks = 0;
    simContext.RunH2SGUI = 0; simContext.RunNORM = 0;

    /*
    Worksheets("Calcite").Range("B1:J1").Value = Null
    Worksheets("Calcite").Range("C4:K13").Value = Null
    Worksheets("Calcite").Range("p4:u13").Value = Null
    Worksheets("Barite").Range("C4:F13").Value = Null
    Worksheets("Barite").Range("B1:J1").Value = Null
    Worksheets("Barite").Range("p4:u13").Value = Null
    Worksheets("Other SO4s").Range("c4:J13").Value = Null
    Worksheets("Other SO4s").Range("B1:J1").Value = Null
    Worksheets("Other SO4s").Range("P4:U14").Value = Null
    Worksheets("Other SO4s").Range("P16:U25").Value = Null
    Worksheets("Halite").Range("C4:F13").Value = Null
    Worksheets("Halite").Range("B1:J1").Value = Null
    Worksheets("Sulfides,Fluorite,Carbonates").Range("C4:W13").Value = Null
    Worksheets("Sulfides,Fluorite,Carbonates").Range("B1:J1").Value = Null
    Worksheets("Silicates").Range("C4:Q13").Value = Null
    Worksheets("Silicates").Range("B1:J1").Value = Null
    Worksheets("Mg(OH)2,Ca(OH)2").Range("C4:I13").Value = Null
    Worksheets("Mg(OH)2,Ca(OH)2").Range("B1:J1").Value = Null
    */

    InitializeOptionClearCellContent();

    // ReDim CaseCount(5), CaseCount_II(110) 
    CaseCount_II.resize(110, 0);

    CountNOB();
    B1_InitializeIndices();
    B2_ReadinAllData(&data);
    B3_CalcConcs(data.API);
    B4_CalcFinalBrine();
    OutputActivityCoefficients();


    // deltaT = (TBH - TWH) / 9
    double deltaT = (simContext.TBH - simContext.TWH) / 9.0;
    double deltaP = (simContext.PBH - simContext.PWH) / 9.0;
    /*
    Worksheets("Calcite").Cells(13, 1) = TBH
    Worksheets("Calcite").Cells(4, 1) = TWH
    Worksheets("Calcite").Cells(13, 2) = PBH
    Worksheets("Calcite").Cells(4, 2) = PWH
   */
   // SI 单位转换逻辑
    if (simContext.UseSI == 1) {
        simContext.TBH = (simContext.TBH - 32.0) * 5.0 / 9.0;
        simContext.TWH = (simContext.TWH - 32.0) * 5.0 / 9.0;
        simContext.PBH /= 14.503774;
        simContext.PWH /= 14.503774;
    }
    /*
      Worksheets("Deposition Prediction").Range("F7") = TBH
      Worksheets("Deposition Prediction").Range("H7") = TWH
      Worksheets("Deposition Prediction").Range("G7") = PBH
      Worksheets("Deposition Prediction").Range("I7") = PWH
      Worksheets("Deposition Prediction").Range("J7") = VgTP / 28.31685 'Mscf/d
      Worksheets("Deposition Prediction").Range("K7") = VO 'bpd
      Worksheets("Deposition Prediction").Range("L7") = VW 'bpd

      Worksheets("Deposition Prediction_WhatIf").Range("F7") = TBH
      Worksheets("Deposition Prediction_WhatIf").Range("H7") = TWH
      Worksheets("Deposition Prediction_WhatIf").Range("G7") = PBH
      Worksheets("Deposition Prediction_WhatIf").Range("I7") = PWH
      Worksheets("Deposition Prediction_WhatIf").Range("J7") = VgTP / 28.31685 'Mscf/d
      Worksheets("Deposition Prediction_WhatIf").Range("K7") = VO 'bpd
      Worksheets("Deposition Prediction_WhatIf").Range("L7") = VW 'bpd

      'dai 2021 corrosion
      Worksheets("Corrosion").Range("F7") = TBH
      Worksheets("Corrosion").Range("H7") = TWH
      Worksheets("Corrosion").Range("G7") = PBH
      Worksheets("Corrosion").Range("I7") = PWH
      Worksheets("Corrosion").Range("J7") = VgTP / 28.31685 'Mscf/d
      Worksheets("Corrosion").Range("K7") = VO 'bpd
      Worksheets("Corrosion").Range("L7") = VW 'bpd

      Worksheets("Corrosion").Range("F11") = TBH
      Worksheets("Corrosion").Range("H11") = TWH
      Worksheets("Corrosion").Range("G11") = PBH
      Worksheets("Corrosion").Range("I11") = PWH
      Worksheets("Corrosion").Range("J11") = VgTP / 28.31685 'Mscf/d
      Worksheets("Corrosion").Range("K11") = VO 'bpd
      Worksheets("Corrosion").Range("L11") = VW 'bpd

    */

    if (Iter_MT_WI == 0) Iter_MT_WI = 1;

    if (simContext.UseSI == 1)
    {
        /*
        Worksheets("Deposition Prediction").Range("F7") = (TBH - 32) * 5 / 9    'default T as TF, P as psia, Output TBH, TWH etc to calcite sheet
        Worksheets("Deposition Prediction").Range("H7") = (TWH - 32) * 5 / 9
        Worksheets("Deposition Prediction").Range("G7") = PBH / 14.503774
        Worksheets("Deposition Prediction").Range("I7") = PWH / 14.503774
        Worksheets("Deposition Prediction").Range("J7") = VgTP / 1000 'Mm3/d
        Worksheets("Deposition Prediction").Range("K7") = VO * 0.159 'm3/d
        Worksheets("Deposition Prediction").Range("L7") = VW * 0.159 'm3/d

        Worksheets("Deposition Prediction_WhatIf").Range("F7") = (TBH - 32) * 5 / 9    'default T as TF, P as psia, Output TBH, TWH etc to calcite sheet
        Worksheets("Deposition Prediction_WhatIf").Range("H7") = (TWH - 32) * 5 / 9
        Worksheets("Deposition Prediction_WhatIf").Range("G7") = PBH / 14.503774
        Worksheets("Deposition Prediction_WhatIf").Range("I7") = PWH / 14.503774
        Worksheets("Deposition Prediction_WhatIf").Range("J7") = VgTP / 1000 'Mm3/d
        Worksheets("Deposition Prediction_WhatIf").Range("K7") = VO * 0.159 'm3/d
        Worksheets("Deposition Prediction_WhatIf").Range("L7") = VW * 0.159 'm3/d

        'dai 2021 corrosion
        Worksheets("Corrosion").Range("F7") = (TBH - 32) * 5 / 9
        Worksheets("Corrosion").Range("H7") = (TWH - 32) * 5 / 9
        Worksheets("Corrosion").Range("G7") = PBH / 14.503774
        Worksheets("Corrosion").Range("I7") = PWH / 14.503774
        Worksheets("Corrosion").Range("J7") = VgTP / 1000 'Mm3/d
        Worksheets("Corrosion").Range("K7") = VO * 0.159 'm3/d
        Worksheets("Corrosion").Range("L7") = VW * 0.159 'm3/d

        Worksheets("Corrosion").Range("F11") = (TBH - 32) * 5 / 9
        Worksheets("Corrosion").Range("H11") = (TWH - 32) * 5 / 9
        Worksheets("Corrosion").Range("G11") = PBH / 14.503774
        Worksheets("Corrosion").Range("I11") = PWH / 14.503774
        Worksheets("Corrosion").Range("J11") = VgTP / 1000 'Mm3/d
        Worksheets("Corrosion").Range("K11") = VO * 0.159 'm3/d
        Worksheets("Corrosion").Range("L11") = VW * 0.159 'm3/d
        */
    }

    if (RunGoalSeek != 1)
    {
        for (LoopTP = 1; LoopTP <= 10; LoopTP++)
        {
            if (simContext.UseTPCalciteSheet == 0)
            {
                TF = simContext.TBH - deltaT * (LoopTP - 1);
                TK = (TF - 32) * 5.0 / 9.0 + 273.15;

                Ppsia = simContext.PBH - deltaP * (LoopTP - 1);
                Patm = Ppsia / 14.696;
            }
            else
            {
                // TF = Worksheets("Calcite").Cells(14 - LoopTP, 1)
                // Ppsia = Worksheets("Calcite").Cells(14 - LoopTP, 2)

                if (simContext.UseSI == 1)
                {
                    //TF = Worksheets("Calcite").Cells(14 - LoopTP, 1) * 9 / 5 + 32
                    //Ppsia = Worksheets("Calcite").Cells(14 - LoopTP, 2) * 14.503774
                }
                TK = (TF - 32) * 5.0 / 9.0 + 273.15;
                Patm = Ppsia / 14.696;
            }

            PBar = Ppsia / 14.503774;
            TC = (TF - 32) * 5.0 / 9.0;

            LoopTPSI(data.API);

            if (H2Oevap != 1)
            {
                if (ISt != 0)
                {
                    LoopTPVisHeatCap();

                    if (Run_MassTransfer == 1)
                        MassTransferCoefficients(data.API);
                }
                LoopTPppt();
            }
            LoopTPWrite();
            if (H2Oevap == 1)
                goto label201;
        }
    }
    else // RunGoalSeek == 1
    {
        for (LoopTP = 0; LoopTP < 2; LoopTP++) {
            if (LoopTP == 0) {
                TF = simContext.TBH;
                TK = (TF - 32) * 5.0 / 9.0 + 273.15;
                TC = (TF - 32) * 5.0 / 9.0;

                Ppsia = simContext.PBH;
                Patm = Ppsia / 14.696;
                PBar = Ppsia / 14.503774;
            }
            if (LoopTP == 1) {
                TF = simContext.TWH;
                TK = (TF - 32) * 5.0 / 9.0 + 273.15;
                TC = (TF - 32) * 5.0 / 9.0;

                Ppsia = simContext.PWH;
                Patm = Ppsia / 14.696;
                PBar = Ppsia / 14.503774;
            }
            LoopTPSI(data.API);
            //LoopTPWrite();
            if (H2Oevap == 1)
                goto label201;
        }
    }
    // mytime1 = Time
    //nettime = (mytime1 - mytime)

    if (RunSSP == 1) {
        // Worksheets("Input").Range("R28") = nettime * 24 * 3600
    }

label201:
    // Application.ScreenUpdating = True

    if (RunGoalSeek == 0 && Run_MassTransfer == 0)
    {
        //ErrorMsgBox();
        if (simContext.UseSI == 1) {
            simContext.TBH = (simContext.TBH - 32) * 5.0 / 9.0;
            simContext.TWH = (simContext.TWH - 32) * 5.0 / 9.0;
            simContext.PBH = simContext.PBH / 14.503774;
            simContext.PWH = simContext.PWH / 14.503774;

            /*
            MsgBox "Calculation is finished. No. of brine mixed = " & nob & Chr(13) & Chr(13) & "The initial and final temperatures are " & TBH & " and " & TWH & " C;" & Chr(13) _
            & "Initial and final pressures are " & PBH & " and " & PWH & " Bar." & Chr(13) & Chr(13) & "Flash calculation option =" & useEOS & "." & Chr(13) & Chr(13) _
            & "Calculation Option(s) on Row 51 used for Input 1 is " & usepHmix(1)  '

            */

        }
        else
        {
            /*
            MsgBox "Calculation is finished. No. of brine mixed = " & nob & Chr(13) & Chr(13) & "The initial and final temperatures are " & TBH & " and " & TWH & " F;" & Chr(13) _
            & "Initial and final pressures are " & PBH & " and " & PWH & " psia." & Chr(13) & Chr(13) & "Flash calculation option =" & useEOS & "." & Chr(13) & Chr(13) _
            & "Calculation Option(s) on Row 51 used for Input 1 is " & Chr(13) & usepHmix(1)
            */
        }
    }

    if (Run_MassTransfer != 1 || simContext.Run10TestCases != 1 ||
        simContext.RunWhatIf != 1 || simContext.Run_Seawater_Mixing != 1 ||
        simContext.Run_MixingTwoWells != 1)
    {
        // Worksheets(mySheet).Activate
        // Worksheets(mySheet).Range("A1").Select
    }

    if (RunGoalSeek != 1 && Run_MassTransfer != 1)
    {
        return 0;  //End
    }
    return 0;
}


int main()
{
    SampleData data;

    initData();
    pointerInit_pf();
    mockData_sheetInput(&data);
    A1_Start_ScaleSoftPitzer(data);
    //B1_InitializeIndices();  //  初始化
    //B2_ReadinAllData(&data);

    printf("Calc Result: TDS = %f, simContext.yCO2 = %f, simContext.yH2S = %f\n", TDS, simContext.yCO2, simContext.yH2S);
    printf("Alkalinity Alk = %f, Tatal Ac TAc = %f\n", Alk, TAc);
    // cleanMemory_pf();
    printf("hello world\n");
    return 0;
}


