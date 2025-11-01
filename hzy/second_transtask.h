#ifndef SECOND_TARNSTASK_H
#define SECOND_TRANSTASK_H

/**
 * @brief 计算离子强度及Pitzer理论常用系统项和函数
 * 
 * 计算离子强度(Ist)、总摩尔浓度(mtotal)、摩尔电荷(MoleCharge)、
 * 阳离子总和(SumOfCations)、阴离子总和(SumOfAnions)和DpHj值。
 * 用于Pitzer电解质溶液理论的相关计算。
 * 
 * @param Ist [输出] 离子强度
 * @param mtotal [输出] 总摩尔浓度
 * @param MoleCharge [输出] 摩尔电荷
 * @param SumOfCations [输出] 阳离子总和
 * @param SumOfAnions [输出] 阴离子总和
 * @param DpHj [输出] pH相关参数
 */
// void CalcIonicStrength(publicpara_m  *glob_var) ;
void CalcIonicStrength(double *Ist, double *mtotal, double *MoleCharge,
                      double *SumOfCations, double *SumOfAnions, double *DpHj);
                      
// 没有找到VB代码
// double CalcRhoTP(publicpara_m  *glob_var);

// 闫师兄的翻译
void C1_ThermodynamicEquilConsts(void);

/*
 * @brief 计算标准部分摩尔体积（V0TP）
 *
 * 此函数计算各种离子的标准部分摩尔体积（单位：cm³/mol），
 * 基于温度（TK）和压力（PBar）的依赖性，使用fV0函数进行计算。
 * 适用于阳离子（Na、K、Mg、Ca、Ba、Sr）和阴离子（Cl、SO4）。
 *
 * @param TK 绝对温度（K）
 * @param PBar 压力（bar）
 * @param V0_c 阳离子标准部分摩尔体积数组 [NumCat]
 * @param V0_a 阴离子标准部分摩尔体积数组 [NumAnion]
 * @param iNa Na+ 索引
 * @param iK K+ 索引
 * @param iMg Mg²+ 索引
 * @param iCa Ca²+ 索引
 * @param iBa Ba²+ 索引
 * @param iSr Sr²+ 索引
 * @param iCl Cl⁻ 索引
 * @param iSO4 SO₄²⁻ 索引
 */
// void V0TP(publicpara_m  *glob_var)
void V0TP(double TK, double PBar, double* V0_c, double* V0_a);

/** 
 * @brief 计算Pitzer模型参数（C2_Pitzer2019）
 * 
 * 此函数计算Pitzer模型中的各种二元和三元相互作用参数，包括b0、b1、b2、CPhi等， 
 * 基于温度（TK）和压力（PBar）依赖性，适用于电解质溶液的热力学建模。 
 * 包括Na、K、Ca、Mg、Ba、Sr、H、OH、Cl、SO4、HCO3、CO3、HS等离子的参数计算。 
 * 
 * @param TK 绝对温度（K） 
 * @param PBar 压力（bar） 
 * @param TC 摄氏温度（°C） 
 * @param Patm 大气压力（bar） 
 * @param b0 二元参数数组 [NumCat+NumAnion][NumCat+NumAnion] 
 * @param b1 二元参数数组 [NumCat+NumAnion][NumCat+NumAnion] 
 * @param b2 二元参数数组 [NumCat+NumAnion][NumCat+NumAnion] 
 * @param CPhi 三元参数数组 [NumCat+NumAnion][NumCat+NumAnion] 
 * @param Lnc 中性物种与阳离子参数 [NumSpecies][NumCat] 
 * @param Lna 中性物种与阴离子参数 [NumSpecies][NumAnion] 
 * @param zeta 三元参数 [NumSpecies][NumCat][NumAnion] 
 * @param Taap 三元参数 [NumCat+NumAnion][NumCat+NumAnion] 
 * @param Yaapc 三元参数 [NumCat+NumAnion][NumCat+NumAnion][NumCat] 
 * @param Yccpa 三元参数 [NumCat][NumCat][NumAnion] 
 * @param Tccp 三元参数 [NumCat][NumCat] 
 * @param NumCat 阳离子数量 
 * @param NumAnion 阴离子数量 
 * @param NumSpecies 物种数量 
 * @param iH H+ 索引 
 * @param iCl Cl- 索引 
 * @param iSO4 SO4^2- 索引 
 * @param iNa Na+ 索引 
 * @param iOH OH- 索引 
 * @param iHCO3 HCO3- 索引 
 * @param iCO3 CO3^2- 索引 
 * @param iK K+ 索引 
 * @param iHS HS- 索引 
 * @param iCa Ca^2+ 索引 
 * @param iMg Mg^2+ 索引 
 * @param iBa Ba^2+ 索引 
 * @param iSr Sr^2+ 索引 
 * @param iFe Fe^2+ 索引 
 * @param iAc Ac- 索引 
 * @param iSion S离子索引 
 * @param iCO2aq CO2(aq) 索引 
 * @param iH2Saq H2S(aq) 索引 
 * @param Tr 参考温度（298.15 K） 
 */
// void C2_Pitzer2019(publicpara_m  *glob_var);
void C2_Pitzer2019(double TK, double PBar, double TC, double Patm, double** b0, double** b1, double** b2, double** CPhi, 
                   double** Lnc, double** Lna, double*** zeta, double** Taap, 
                   double*** Yaapc, double*** Yccpa, double** Tccp, 
                   int NumCat, int NumAnion, int NumSpecies);

/**
 * @brief 计算Pitzer活性系数（温度、压力、离子强度版本）
 *
 * 此函数使用温度(TK)、温度(摄氏度 TC)、压力(PBar)和大气压(Patm)计算阴离子、阳离子和中性化合物的活性系数。
 * 基于Pitzer离子相互作用模型，包含二元、三元和中性项的计算。还计算水的渗透系数和活性，以及特定络合物的Bdot参数。
 * 注意：函数依赖外部数组和全局变量（如mc, ma, b0, b1等），这些需在调用前定义。
 * 特殊处理：对于某些离子对（如Na-SO4, K-SO4, Ca-SO4）使用自定义alpha值。
 *
 * @param gNeut [输出] 中性化合物的活性系数数组
 * @param aH2O [输出] 水的活性
 * @param TK 绝对温度 (K)
 * @param TC 温度 (摄氏度)
 * @param PBar 压力 (bar)
 * @param Patm 大气压 (atm)
 */
void C2_PitzerActCoefs_T_P_ISt(double *gNeut, double *aH2O, double TK, double TC, double PBar, double Patm);

/**
 * @brief Peng-Robinson状态方程计算气体逸度系数
 * 
 * 此函数使用Peng-Robinson状态方程计算CH4、CO2和H2S气体的逸度系数。
 * 基于温度、压力和气体组成，计算混合气体的压缩因子和各组分的逸度系数。
 * 
 * @param TK 温度 (K)
 * @param PBar 压力 (bar)
 * @param yCH4 甲烷摩尔分数
 * @param yCO2 二氧化碳摩尔分数  
 * @param yH2S 硫化氢摩尔分数
 * @param phiCH4 [输出] 甲烷逸度系数
 * @param phiCO2 [输出] 二氧化碳逸度系数
 * @param phiH2S [输出] 硫化氢逸度系数
 * @param Znew [输出] 气体压缩因子
 */
void PengRobinson3(double TK, double PBar, double yCH4, double yCO2, double yH2S,
                   double *phiCH4, double *phiCO2, double *phiH2S, double *Znew);

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
 * @param yCO2 CO2摩尔分数
 * @param yH2S H2S摩尔分数
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
 * @param errmsg [输出] 错误信息数组
 */
void C5_CalcpHPCO2PH2SSTP(int use_pH, int UseH2Sgas, int useEOS,
                         double TK, double Ppsia, double *yCO2, double *yH2S,
                         double Alk, double TAc, double TH2Saq, double TFe, double TCO2,
                         double TNH4, double TH3BO3, double TH4SiO4,
                         double *pH, double *pHMeterReading, int *errmsg);

/**
 * @brief 计算密度和pH（D2_CalcDensitypH）
 *
 * 此函数计算溶液的密度和pH值，通过迭代方法调整浓度以实现密度收敛。
 * 包括离子强度计算、pH调整、密度迭代、浓度修正，以及热力学平衡常数和活度系数的计算。
 * 适用于电解质溶液的密度和pH计算，处理CO2、H2S、FeSaq等物种的物种分布。
 *
 * @param i 当前混合物索引
 * @param pH [输出] pH值
 * @param Ist [输入/输出] 离子强度
 * @param rhoOld [输入/输出] 旧密度值
 * @param Iteration [输入/输出] 迭代次数
 * @param rhoSSE [输入/输出] 密度平方误差
 * @param TK 绝对温度（K）
 * @param TC 摄氏温度（°C）
 * @param PBar 压力（bar）
 * @param Patm 大气压力（bar）
 * @param mc 阳离子摩尔浓度数组 [MaxSpecies]
 * @param ma 阴离子摩尔浓度数组 [MaxSpecies]
 * @param mn 中性物种摩尔浓度数组 [MaxSpecies]
 * @param Alk [输入/输出] 碱度
 * @param TAc [输入/输出] 总醋酸浓度
 * @param TCO2 [输入/输出] 总CO2浓度
 * @param TNH4 [输入/输出] 总NH4浓度
 * @param TH3BO3 [输入/输出] 总H3BO3浓度
 * @param TH2Saq [输入/输出] 总H2S(aq)浓度
 * @param TH4SiO4 [输入/输出] 总H4SiO4浓度
 * @param TFe [输入/输出] 总Fe浓度
 * @param TDS 总溶解固体（mg/L）
 * @param NumCat 阳离子数量
 * @param NumAn 阴离子数量
 * @param NumNeut 中性物种数量
 * @param iH H+ 索引
 * @param iOH OH- 索引
 * @param iAc Ac- 索引
 * @param iNH3 NH3 索引
 * @param iH2BO3 H2BO3- 索引
 * @param iHCO3 HCO3- 索引
 * @param iCO3 CO3^2- 索引
 * @param iHS HS- 索引
 * @param iH3SiO4 H3SiO4- 索引
 * @param iH2SiO4 H2SiO4^2- 索引
 * @param iH4SiO4aq H4SiO4(aq) 索引
 * @param iNH4 NH4+ 索引
 * @param iCO2aq CO2(aq) 索引
 * @param iH2Saq H2S(aq) 索引
 * @param iHAcaq HAc(aq) 索引
 * @param useEOSmix [数组] EOS混合标志 [MaxMix]
 * @param kk 当前EOS索引
 * @param gNeut 中性活度系数数组 [MaxNeut]
 * @param aH2O [输出] 水活度
 * @param DpHj pH校正
 * @param H [输入] H+ 浓度（从C5_CalcpHPCO2PH2SSTP）
 * @param OH [输入] OH- 浓度
 * @param AC [输入] Ac- 浓度
 * @param NH3 [输入] NH3 浓度
 * @param H2BO3 [输入] H2BO3- 浓度
 * @param HCO3 [输入] HCO3- 浓度
 * @param CO3 [输入] CO3^2- 浓度
 * @param HS [输入] HS- 浓度
 * @param H3SiO4 [输入] H3SiO4- 浓度
 * @param H2SiO4 [输入] H2SiO4^2- 浓度
 * @param H4SiO4 [输入] H4SiO4 浓度
 * @param CO2aq [输入] CO2(aq) 浓度
 * @param H2Saq [输入] H2S(aq) 浓度
 * @param HAcaq [输入] HAc(aq) 浓度
 * @param xMeOH [输出] 甲醇摩尔分数
 * @param xMEG [输出] 乙二醇摩尔分数
 * @param IStCosolvent [输出] 共溶剂离子强度
 * @param mt [输出] 温度压力函数值
 * @param use_pH pH使用选项
 * @param pHMeterStpMix pH 测量值
 * @param rho25c [输入/输出] 25°C密度
 */
void D2_CalcDensitypH(int i, double *pH, double *Ist, double *rhoOld, int *Iteration, double *rhoSSE, double TK, double TC, double PBar, double Patm,
                      double* mc, double* ma, double* mn, double *Alk, double *TAc, double *TCO2, double *TNH4, double *TH3BO3,
                      double *TH2Saq, double *TH4SiO4, double *TFe, double TDS, int NumCat, int NumAn, int NumNeut,
                      int* useEOSmix, int kk, double* gNeut, double *aH2O, double DpHj,
                      double H, double OH, double AC, double NH3, double H2BO3, double HCO3, double CO3, double HS, double H3SiO4,
                      double H2SiO4, double H4SiO4, double CO2aq, double H2Saq, double HAcaq, double *xMeOH, double *xMEG, double *IStCosolvent,
                      double *mt, int use_pH, double *pHMeterStpMix, double *rho25c);

/**
 * @brief 计算密度（D1_CalcDensity）
 *
 * 此函数计算混合物i的密度，通过迭代调整TDS和浓度以实现收敛。
 * 支持摩尔浓度（UseMolal=0）和直接计算（UseMolal=1）两种模式。
 * 在摩尔浓度模式下，进行TDS迭代；在直接模式下，计算离子强度、活度系数和物种分布。
 *
 * @param i 当前混合物索引
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
 * @param Iteration2 [输入/输出] 迭代2计数
 * @param mc 阳离子摩尔浓度数组 [MaxSpecies]
 * @param ma 阴离子摩尔浓度数组 [MaxSpecies]
 * @param mn 中性物种摩尔浓度数组 [MaxSpecies]
 * @param iH H+ 索引
 * @param iNa Na+ 索引
 * @param iK K+ 索引
 * @param iMg Mg^2+ 索引
 * @param iCa Ca^2+ 索引
 * @param TCa [输出] 总Ca浓度
 * @param iSr Sr^2+ 索引
 * @param iBa Ba^2+ 索引
 * @param iFe Fe^2+ 索引
 * @param iZn Zn^2+ 索引
 * @param iPb Pb^2+ 索引
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
 * @param iNH3 NH3 索引
 * @param iH3BO3 H3BO3 索引
 * @param iH4SiO4aq H4SiO4(aq) 索引
 * @param TFe [输入/输出] 总Fe浓度
 * @param TDSMix [数组] TDS混合 [MaxMix]
 * @param TDSOld [输入/输出] 旧TDS值
 * @param TDS [输入/输出] 总溶解固体（mg/L）
 * @param TDSSSE [输入/输出] TDS平方误差
 * @param UseMolal 使用摩尔浓度标志
 * @param rho25c [输入/输出] 25°C密度
 * @param CalculateTDSDen [输入/输出] 计算TDS密度
 * @param NumCat 阳离子数量
 * @param NumAn 阴离子数量
 * @param NumNeut 中性物种数量
 * @param MWCat 阳离子分子量数组 [MaxCat]
 * @param MWAn 阴离子分子量数组 [MaxAnion]
 * @param MWNeut 中性物种分子量数组 [MaxNeut]
 * @param xMeOH 甲醇摩尔分数
 * @param xMEG 乙二醇摩尔分数
 * @param IStCosolvent 共溶剂离子强度
 * @param Ist [输入/输出] 离子强度
 * @param mt [输出] 温度压力函数值
 * @param pH [输出] pH值
 * @param pHMeterStpMix pH 测量值
 * @param DpHj pH校正
 * @param gNeut 中性活度系数数组 [MaxNeut]
 * @param aH2O [输出] 水活度
 * @param TK 绝对温度（K）
 * @param TC 摄氏温度（°C）
 * @param PBar 压力（bar）
 * @param Patm 大气压力（bar）
 * @param use_pH pH使用选项
 * @param useEOSmix [数组] 使用EOS混合 [MaxMix]
 * @param kk 当前EOS索引
 * @param H [输入] H+ 浓度（从C5_CalcpHPCO2PH2SSTP）
 * @param OH [输入] OH- 浓度
 * @param AC [输入] Ac- 浓度
 * @param NH3 [输入] NH3 浓度
 * @param H2BO3 [输入] H2BO3- 浓度
 * @param HCO3 [输入] HCO3- 浓度
 * @param CO3 [输入] CO3^2- 浓度
 * @param HS [输入] HS- 浓度
 * @param H3SiO4 [输入] H3SiO4- 浓度
 * @param H2SiO4 [输入] H2SiO4^2- 浓度
 * @param H4SiO4 [输入] H4SiO4 浓度
 * @param CO2aq [输入] CO2(aq) 浓度
 * @param H2Saq [输入] H2S(aq) 浓度
 * @param HAcaq [输入] HAc(aq) 浓度
 * @param iCO2aq CO2(aq) 索引
 * @param iH2Saq H2S(aq) 索引
 * @param iHAcaq HAc(aq) 索引
 * @param yCO2 [输入/输出] CO2气体摩尔分数
 * @param yH2S [输入/输出] H2S气体摩尔分数
 * @param yCH4 [输入/输出] CH4气体摩尔分数
 * @param rho_Mix [数组] 密度混合 [MaxMix]
 */
void D1_CalcDensity(int i, double* HCO3stpMix, double* AlkMix, double* ACstpMix, double* TAcMix, double* HstpMix, double* OHstpMix,
                    double* CO3stpMix, double* HSstpMix, double* NH4STPMix, double* TNH4Mix, double* H2BO3stpMix, int *Iteration2,
                    double* mc, double* ma, double* mn, 
                    double *Alk, double *TAc, double *TH2Saq, double *TH4SiO4, double* TCa,
                    double *TH3BO3, double *TNH4,
                    double *TFe, double* TDSMix, double *TDSOld,
                    double *TDS, double *TDSSSE, int UseMolal, double *rho25c, double *CalculateTDSDen, int NumCat, int NumAn,
                    int NumNeut, double* MWCat, double* MWAn, double* MWNeut, double *xMeOH, double *xMEG, double *IStCosolvent,
                    double *Ist, double *mt, double *pH, double *pHMeterStpMix,  double DpHj, double gNeut[], double *aH2O,
                    double TK, double TC, double PBar, double Patm, int use_pH, int* useEOSmix, int kk, double H, double OH, double AC,
                    double NH3, double H2BO3, double HCO3, double CO3, double HS, double H3SiO4, double H2SiO4, double H4SiO4,
                    double CO2aq, double H2Saq, double HAcaq,
                    double *yCO2, double *yH2S,
                    double *yCH4, double* rho_Mix);

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
void ReadInputPartC(int kk, double *mt, int* UseTPpHMix, int Run10TestCases, int Loop10, int Run_Seawater_Mixing, int LoopMixing,
                    int Run_MixingTwoWells, int RunMultiMix, int LoopResChem, int RunStatMix,
                    double* HCO3stpMix, double* AlkMix, double* ACstpMix, double* TAcMix, double* HstpMix, double* OHstpMix,
                    double* CO3stpMix, double* HSstpMix, double* NH4STPMix, double* TNH4Mix, double* H2BO3stpMix, double *TDS,
                    double *yH2S, double *yCO2, int *Iteration2,
                    double* mc, double* ma, double* mn,
                    double *Alk, double *TAc, double *TH2Saq, double *TH4SiO4, double *TH3BO3, double *TNH4, double *TFe, double *TPb,
                    double *TZn, int *use_pH, int* usepHmix, int *UseH2Sgas, int* UseH2SgasMix,
                    double *TCO2, double* TCO2Mix, double* yCO2Mix, double* yH2SMix, int* useEOSmix, double* SumofZMix,
                    double** zMix, double* CalculatedTDSMix, double* NaMix, double* KMix, double* MgMix, double* CaMix,
                    double *TCa, double* SrMix, double* BaMix, double* FeMix, double* ZnMix, double* PbMix, double* RaMix,
                    double* ClMix, double* SO4Mix, double* FMix, double* BrMix, double* rho25CMix, double* H3SiO4Mix, double* H2SiO4Mix,
                    double* NH3Mix, double* H4SiO4Mix, double* H3BO3Mix, double* CO2aqMix, double* H2SaqMix, double* HACaqMix,
                    double* AlkMix2, double* TAcMix2, double* TH2SaqMix, double* TH4SiO4Mix2, double* TNH4Mix2, double* TH3BO3Mix2,
                    double* yCH4Mix, int* UseTPVolMix, double* WaterDensityMix, double rho25c, int UseMolal, double *CalculateTDSDen,
                    int NumCat, int NumAn, int NumNeut, double* MWCat, double* MWAn, double* MWNeut
                    /*int iCO2aq, int iH2Saq, int iHAcaq*/
                );

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
void Get_EOS_Parameters(int NumGases, char* EOS, double* MWgas, double* TCr, double* PCr, double* Omega, double* mf_c0, double* mf_c1, double** kPr);

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
void pseudo_composition(double API, double SGG, double VgTP, double mol_opd,
    double mol_wpd, double TK, double PBar, double aH2O,
    double gNeut[], double nTCO2, double nTH2S, double yCO2,
    double yH2S, double YH2O, double *total_moles,
    double feed_Composition[], double *mol_HC);


#endif // SECOND_TRANSTASK_H