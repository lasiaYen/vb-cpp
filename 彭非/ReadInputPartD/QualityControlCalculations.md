# QualityControlCalculations

> code

```cpp
void QualityControlCalculations(int kk, int j,
    int* UseH2Sgas, double* pH, double* aH, double* yH2S, int useEOS, double K1H2S, double KgwH2S,
    double Ppsia, double KstFeSaq, double TZn, double TPb, double hydHS, double ppt, double* BetaDot, double* gDot,
    double KH2O, double KHAc, double KH3BO3, double KNH4, double KH4SiO4, double KH3SiO3, double K1H2CO3, double KgwCO2,
    double K2HCO3,
    double* yCO2,
    double* pHMeterReading
)
{
    double root1 = 0, root2 = 0, root3 = 0 ;

    double pHMeterReading_from_QC = 0;

    if (RunStat == 1) {
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
        if (UseMolal == 0){
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
    *UseH2Sgas = globalvar.UseH2SgasMix[kk];
    SumOfCations = 0.0000001;
    for (int c = 0; c < NumCat; c++)
        SumOfCations += ChCat[c] * mc[c];

    SumOfAnions = -0.0000001;
    for (int a = 0; a < NumCat; a++)
        SumOfAnions += ChAn[a] * ma[a];

    //1、使用测量的 pH 值和 Alk 值来计算输入表 QC 部分中的 P-CO2。
    *pH = globalvar.pHMeterStpMix[kk] + DpHj;
    *aH = pow(10, -(*pH));


    if (RunStat == 1) {
        //yCO2 = Worksheets(mySheet).Cells(26, j + 2).Value / 100
        //TH2Saq = Worksheets(mySheet).Cells(27, j + 2).Value / (34080 * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#))
        *yH2S = 0;
    }
    else {
        // yCO2 = Worksheets(mySheet).Cells(31, j + 2).Value / 100

        if (globalvar.UseH2SgasMix[kk] == 1) {
            //yH2S = Worksheets(mySheet).Cells(33, j + 2).Value / 100
            TH2Saq = 0;
        }
        else {
            // TH2Saq = Worksheets(mySheet).Cells(33, j + 2).Value
            *yH2S = 0;
        }

    }


    if (useEOS == 0) {
        if (*UseH2Sgas == 0) {
            // Calculate TH2Saq from yH2Sstp and pH. If TH2Saq is given, use it, otherwise use PH2S.
            HS = TH2Saq / (*aH * gAn[iHS] * gNAn[iHS] / (K1H2S * gNeut[iH2Saq] * gNNeut[iH2Saq]) + 1);

            if (TH2Saq > 0) {
                // speciation for HS
                // fMeSSpeciation(int im, int igas,double TZn, double TPb, double hydHS, double ppt, double* root1, double* root2, double* root3, double* BetaDot, double* gDot, double KstFeSaq, double aH)
                fMeSSpeciation(2, 2, TZn, TPb, hydHS, ppt, &root1, &root2, &root3, BetaDot, gDot, KstFeSaq, *aH);
                H2Saq = *aH * HS * gAn[iHS] * gNAn[iHS] / (K1H2S * gNeut[iH2Saq] * gNNeut[iH2Saq]);
                *yH2S = H2Saq * gNeut[iH2Saq] * gNNeut[iH2Saq] / (KgwH2S * Ppsia * gGas[iH2Sg]);

                if (*yH2S > 1) {
                    errmsg[2] = 3; // 数组下标从0开始，3对应索引2
                }
            }
            else {
                *yH2S = 0;
            }
        }

        // Calculate the P-H2S, or yH2S, from pH and TH2Saq.
        if (*UseH2Sgas == 1) {
            if (*yH2S > 0) {
                H2Saq = KgwH2S * Ppsia * (*yH2S) * gGas[iH2Sg] / gNeut[iH2Saq] / gNNeut[iH2Saq];
                HS = K1H2S * H2Saq * gNeut[iH2Saq] * gNNeut[iH2Saq] / (*aH * gAn[iHS] * gNAn[iHS]);

                // 计算各种配合物形成的ZP参数
                double ZP1 = (KstFeSaq * gAn[iHS] * gNAn[iHS] * gCat[iFe] * gNCat[iFe]) /
                    (gNeut[iFeSaq] * gNNeut[iFeSaq] * *aH);

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
                mc[iFe] = TFe / (1 + ZP1 * HS);
                mc[iZn] = TZn / (1 + ZP2 + ZP3 * pow(HS, 2) + ZP4 * pow(HS, 3));
                mc[iPb] = TPb / (1 + ZP5 + ZP6 * pow(HS, 2) + ZP7 * pow(HS, 3));

                // 计算总H2S浓度
                TH2Saq = H2Saq + (HS + mc[iFe] * ZP1 * HS +
                    mc[iZn] * (ZP3 * pow(HS, 2) + ZP4 * pow(HS, 3)) +
                    mc[iPb] * (ZP6 * pow(HS, 2) + ZP7 * pow(HS, 3)));
            }
            else {
                TH2Saq = 0;
            }
        }

        if (TH2Saq == 0 && *yH2S == 0) {
            HS = 0;
        }

        // 计算其他化学物种的浓度
        H = *aH / gCat[iH] / gNCat[iH];
        OH = KH2O / (*aH * gAn[iOH] * gNAn[iOH]);

        AC = TAc / (*aH * gAn[iAc] * gNAn[iAc] / (KHAc * gNeut[iHAcaq] * gNNeut[iHAcaq]) + 1); // Note gNNeut(iHACaq)=1

        double hydH2BO3 = *aH * gAn[iH2BO3] * gNAn[iH2BO3] / (KH3BO3 * gNeut[iH3BO3] * gNNeut[iH3BO3]) + 1;
        H2BO3 = TH3BO3 / hydH2BO3;

        double hydNH3 = *aH * gNeut[iNH3] * gNNeut[iNH3] / (KNH4 * gCat[iNH4] * gNCat[iNH4]) + 1;
        NH3 = TNH4 / hydNH3;

        double hydH2SiO4 = pow(*aH, 2) * gAn[iH2SiO4] * gNAn[iH2SiO4] /
            (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]) +
            (*aH) * gAn[iH2SiO4] * gNAn[iH2SiO4] /
            (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]) + 1;
        H2SiO4 = TH4SiO4 / hydH2SiO4;
        H3SiO4 = H2SiO4 * (*aH) * gAn[iH2SiO4] * gNAn[iH2SiO4] /
            (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]);
        H4SiO4 = H2SiO4 * pow(*aH, 2) * gAn[iH2SiO4] * gNAn[iH2SiO4] /
            (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]);

        // 计算CO2相关的参数
        double tHCO3 = (K1H2CO3 * aH2O) * KgwCO2 * Ppsia * gGas[iCO2g] /
            (*aH * gAn[iHCO3] * gNAn[iHCO3]);
        double tCO3 = (K1H2CO3 * aH2O) * K2HCO3 * KgwCO2 * Ppsia * gGas[iCO2g] /
            (pow(*aH, 2) * gAn[iCO3] * gNAn[iCO3]);

        // 计算yCO2
        *yCO2 = (Alk + H - AC - NH3 - H2BO3 - HS - H3SiO4 - 2 * H2SiO4 - OH) /
            (tHCO3 + 2 * tCO3);

        HCO3 = (K1H2CO3 * aH2O) * KgwCO2 * Ppsia * (*yCO2) * gGas[iCO2g] /
            (*aH * gAn[iHCO3] * gNAn[iHCO3]);
        CO3 = (K1H2CO3 * aH2O) * K2HCO3 * KgwCO2 * Ppsia * (*yCO2) * gGas[iCO2g] /
            (pow(*aH, 2) * gAn[iCO3] * gNAn[iCO3]);

        if ((*yCO2) > 1) {
            errmsg[0] = 1;
        }

        if (*yCO2 < 0) {
            errmsg[1] = 2; // 数组下标从0开始，2对应索引1
            *yCO2 = 0;  // Caps %CO2 at 0
        }

        // 输出结果到相应的工作表
        if (RunStat == 1 && CaseCount[1] == 1) {
            // for StatQC produced water
            //Worksheets("Output data sheet").Cells(6, 4) = yCO2 * 100  '假设 pH 和 Alk 正确，QC 检查 PCO2 是多少？

        }
        else if (RunStat == 1 && CaseCount[1] == 9) {
            // for StatQC fresh water
            //Worksheets("Output data sheet").Cells(18, 4) = yCO2 * 100  'QC check for assuming that pH and Alk are correct, what would PCO2 be?
        }
        else if (RunH2SGUI == 1) {
            //Worksheets(mySheet).Cells(83, j + 2) = yH2S * 100

            if (UseMolal == 0) {
                //Worksheets(mySheet).Cells(84, j + 2) = TH2Saq * (34080 * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#))
            }
            else {
                //Worksheets(mySheet).Cells(84, j + 2) = TH2Saq
            }

            //Worksheets(mySheet).Cells(86, j + 2) = yCO2 * 100 'QC check for assuming that pH and Alk are correct, what would PCO2 be?
        }
        else {
            if (kk <= nob_Input) {
                //Worksheets("Input").Cells(14, j + 15) = yH2S * 100

                if (UseMolal == 0) {
                    //Worksheets("Input").Cells(15, j + 15) = TH2Saq * (34080 * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#))
                }
                else {
                    //Worksheets("Input").Cells(15, j + 15) = TH2Saq
                }

                //Worksheets("Input").Cells(17, j + 15) = yCO2 * 100  'QC check for assuming that pH and Alk are correct, what would PCO2 be?
            }

            if (kk > nob_Input && kk <= nob_Input + nob_InputII) {
                // Worksheets("Input II").Cells(83, j + 2) = yH2S * 100

                if (UseMolal == 0) {
                    //Worksheets("Input II").Cells(84, j + 2) = TH2Saq * (34080 * (rho25CMix(i + nob_Input) - CalculatedTDSMix(i + nob_Input) / 1000000#))
                }
                else {
                    //Worksheets("Input II").Cells(84, j + 2) = TH2Saq
                }

                //Worksheets("Input II").Cells(86, j + 2) = yCO2 * 100 'QC check for assuming that pH and Alk are correct, what would PCO2 be?
            }
        }


        //2、使用 STP P-CO2 和碱度来计算 pH 值。
        if (UseH2Sgas == 0) {
            // 然后，必须指定yCO2stp。需要计算pH值和yH2Sstp。

            // 读取yCO2值（从相应的单元格）
            if (RunStat == 1) {
                // yCO2 = Worksheets(mySheet).Cells(26, j + 2).Value / 100
            }
            else {
                //yCO2 = Worksheets(mySheet).Cells(31, j + 2).Value / 100
            }

            // 半区间搜索根查找器对于 pH 方程非常有效。
            double pHHigh = 14.0;
            double pHLow = 0.0;
            double faH;
            double hydHS, hydAc, hydH2BO3, hydNH3, hydH2SiO4;
            int k;

            // 二分法求解pH值（30次迭代）
            for (k = 0; k < 30; k++) {
                *pH = (pHHigh + pHLow) / 2.0;
                *aH = pow(10.0, -(*pH));

                // 计算H+和OH-浓度
                H = *aH / (gCat[iH] * gNCat[iH]);
                OH = KH2O / (*aH * gAn[iOH] * gNAn[iOH]);

                // 计算CO2系统物种浓度
                CO2aq = KgwCO2 * Ppsia * (*yCO2) * gGas[iCO2aq] / (gNeut[iCO2aq] * gNNeut[iCO2aq]);
                HCO3 = (K1H2CO3 * aH2O) * CO2aq * gNeut[iCO2aq] * gNNeut[iCO2aq] /
                    (*aH * gAn[iHCO3] * gNAn[iHCO3]);
                CO3 = K2HCO3 * HCO3 * gAn[iHCO3] * gNAn[iHCO3] /
                    (*aH * gAn[iCO3] * gNAn[iCO3]);

                // 计算H2S系统物种浓度
                hydHS = *aH * gAn[iHS] * gNAn[iHS] / (K1H2S * gNeut[iH2Saq] * gNNeut[iH2Saq]) + 1.0;
                if (TH2Saq > 0) {
                    HS = TH2Saq / hydHS;
                    //fMeSSpeciation(im, igas, HS, TFe, TZn, TPb, TH2Saq, hydHS, *ppt, & root1, & root2, & root3);
                    fMeSSpeciation(2, 2, TZn, TPb, hydHS, ppt, &root1, &root2, &root3, BetaDot, gDot, KstFeSaq, *aH); // speciation for HS
                }
                else {
                    HS = 0.0;
                }

                // 计算乙酸系统物种浓度
                hydAc = *aH * gAn[iAc] * gNAn[iAc] / (KHAc * gNeut[iHAcaq] * gNNeut[iHAcaq]) + 1.0; // Note gNNeut(iHACaq)=1
                AC = TAc / hydAc;

                // 计算硼酸系统物种浓度
                hydH2BO3 = *aH * gAn[iH2BO3] * gNAn[iH2BO3] / (KH3BO3 * gNeut[iH3BO3] * gNNeut[iH3BO3]) + 1.0;
                H2BO3 = TH3BO3 / hydH2BO3;

                // 计算氨系统物种浓度
                hydNH3 = *aH * gNeut[iNH3] * gNNeut[iNH3] / (KNH4 * gCat[iNH4] * gNCat[iNH4]) + 1.0;
                NH3 = TNH4 / hydNH3;

                // 计算硅酸系统物种浓度
                hydH2SiO4 = pow(*aH, 2) * gAn[iH2SiO4] * gNAn[iH2SiO4] /
                    (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]) +
                    *aH * gAn[iH2SiO4] * gNAn[iH2SiO4] /
                    (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]) + 1.0;
                H2SiO4 = TH4SiO4 / hydH2SiO4;
                H3SiO4 = H2SiO4 * (*aH) * gAn[iH2SiO4] * gNAn[iH2SiO4] /
                    (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]);
                H4SiO4 = H2SiO4 * pow(*aH, 2) * gAn[iH2SiO4] * gNAn[iH2SiO4] /
                    (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]);

                // 计算碱度平衡误差
                faH = Alk - (HCO3 + 2 * CO3 + HS + AC + NH3 + H2BO3 + H3SiO4 + 2 * H2SiO4 + OH - H);

                // 调整pH搜索区间
                if (faH > 0) {
                    pHLow = *pH;
                }
                if (faH < 0) {
                    pHHigh = *pH;
                }
            }

            // 计算pH计读数
            pHMeterReading_from_QC = *pH - DpHj;

            if (*yH2S > 1) errmsg[2] = 3; // 注意：yH2S在此代码段中未计算

            // 根据条件设置最终pH计读数
            if (Run_Seawater_Mixing == 0 || Run_MixingTwoWells == 0) {
                pHMeterReading_from_QC = *pH - DpHj;
            }
        }


        if (*UseH2Sgas == 1) {
            // Use PCO2 and Alk to calculate the pH value at STP
            // 如果使用 PH2S，简化输入将跳过此部分

            // 读取yCO2和yH2S值
            //待定：yCO2 = Worksheets(mySheet).Cells(31, j + 2).Value / 100
            //待定：yH2S = Worksheets(mySheet).Cells(33, j + 2).Value / 100

            // 二分法求解pH值
            double pHHigh = 14.0;
            double pHLow = 0.0;
            double faH;
            double hydAc, hydH2BO3, hydNH3, hydH2SiO4;
            int k;

            // 这使得 pH 值能够收敛到 8 位有效数字；
            // 如果沉淀物仅占碱度的一小部分，有时就需要这样做。
            for (k = 0; k < 30; k++) {
                *pH = (pHHigh + pHLow) / 2.0;
                *aH = pow(10.0, -(*pH));

                // 计算H+和OH-浓度
                H = *aH / (gCat[iH] * gNCat[iH]);
                OH = KH2O / (*aH * gAn[iOH] * gNAn[iOH]);

                // 计算CO2系统物种浓度
                CO2aq = KgwCO2 * Ppsia * (*yCO2) * gGas[iCO2g] / (gNeut[iCO2aq] * gNNeut[iCO2aq]);
                HCO3 = (K1H2CO3 * aH2O) * CO2aq * gNeut[iCO2aq] * gNNeut[iCO2aq] /
                    (*aH * gAn[iHCO3] * gNAn[iHCO3]);
                CO3 = K2HCO3 * HCO3 * gAn[iHCO3] * gNAn[iHCO3] /
                    (*aH * gAn[iCO3] * gNAn[iCO3]);

                // 计算H2S系统物种浓度
                H2Saq = KgwH2S * Ppsia * (*yH2S) * gGas[iH2Sg] / gNeut[iH2Saq] / gNNeut[iH2Saq];
                HS = K1H2S * H2Saq * gNeut[iH2Saq] * gNNeut[iH2Saq] /
                    (*aH * gAn[iHS] * gNAn[iHS]);

                // 计算乙酸系统物种浓度
                hydAc = *aH * gAn[iAc] * gNAn[iAc] / (KHAc * gNeut[iHAcaq] * gNNeut[iHAcaq]) + 1.0;
                AC = TAc / hydAc;

                // 计算硼酸系统物种浓度
                hydH2BO3 = *aH * gAn[iH2BO3] * gNAn[iH2BO3] / (KH3BO3 * gNeut[iH3BO3] * gNNeut[iH3BO3]) + 1.0;
                H2BO3 = TH3BO3 / hydH2BO3;

                // 计算氨系统物种浓度
                hydNH3 = *aH * gNeut[iNH3] * gNNeut[iNH3] / (KNH4 * gCat[iNH4] * gNCat[iNH4]) + 1.0;
                NH3 = TNH4 / hydNH3;

                // 计算硅酸系统物种浓度
                hydH2SiO4 = pow(*aH, 2) * gAn[iH2SiO4] * gNAn[iH2SiO4] /
                    (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]) +
                    *aH * gAn[iH2SiO4] * gNAn[iH2SiO4] /
                    (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]) + 1.0;
                H2SiO4 = TH4SiO4 / hydH2SiO4;
                H3SiO4 = H2SiO4 * (*aH) * gAn[iH2SiO4] * gNAn[iH2SiO4] /
                    (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]);
                H4SiO4 = H2SiO4 * pow(*aH, 2) * gAn[iH2SiO4] * gNAn[iH2SiO4] /
                    (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]);

                // 计算碱度平衡误差
                faH = Alk - (HCO3 + 2 * CO3 + HS + AC + NH3 + H2BO3 + H3SiO4 + 2 * H2SiO4 + OH - H);

                // 调整pH搜索区间
                if (faH > 0) {
                    pHLow = *pH;
                }
                if (faH < 0) {
                    pHHigh = *pH;
                }
            }

            // 计算pH计读数
            pHMeterReading_from_QC = *pH - DpHj;
        }

        // 输出结果到相应的工作表
        if (RunStat == 1 && CaseCount[1] == 1) {
            // for StatQC produced water
            //Worksheets("Output data sheet").Cells(5, 4) = pHMeterReading_from_QC  'QC check given Alk and PCO2, the calculated pH for meter reading.
        }
        else if (RunStat == 1 && CaseCount[1] == 9) {
            // for StatQC fresh water
            //Worksheets("Output data sheet").Cells(17, 4) = pHMeterReading_from_QC  'QC check given Alk and PCO2, the calculated pH for meter reading.
        }
        else if (RunH2SGUI == 1) {
            //Worksheets(mySheet).Cells(85, j + 2) = pHMeterReading_from_QC  'QC check given Alk and PCO2, the calculated pH for meter reading.
        }
        else {
            if (kk <= nob_Input) {
                //Worksheets("Input").Cells(16, j + 15) = pHMeterReading_from_QC  'QC check given Alk and PCO2, the calculated pH for meter reading.
            }
            if (kk > nob_Input && kk <= nob_Input + nob_InputII) {
                //Worksheets("Input II").Cells(85, j + 2) = pHMeterReading_from_QC  'QC check given Alk and PCO2, the calculated pH for meter reading.
            }
        }
    }
    else {
        //对应于 useEOS<>0 QC 中的大多数参数已在 ReadInDataD sub 的末尾计算
        if (RunStat == 1 && CaseCount[1] == 1) {
            // for StatQC produced water
            //Worksheets("Output data sheet").Cells(6, 4) = compositions(2, 2) * 100
            //    Worksheets("Output data sheet").Cells(5, 4) = pHMeterReading
        }
        else if (RunStat == 1 && CaseCount[1] == 9) {
            // for StatQC fresh water
            //Worksheets("Output data sheet").Cells(18, 4) = compositions(2, 2) * 100
                //Worksheets("Output data sheet").Cells(17, 4) = pHMeterReading
        }
        else if (RunH2SGUI == 1) {
            //Worksheets(mySheet).Cells(83, j + 2) = compositions(3, 2) * 100

            if (UseMolal == 0) {
                //Worksheets(mySheet).Cells(84, j + 2) = (H2Saq + HS + S) * (34080 * (rho25c - TDS / 1000000#))
            }
            else {
                //Worksheets(mySheet).Cells(84, j + 2) = (H2Saq + HS + S)
            }

            //Worksheets(mySheet).Cells(85, j + 2) = pHMeterReading
            //    Worksheets(mySheet).Cells(86, j + 2) = compositions(2, 2) * 100 'QC check for assuming that pH and Alk are correct, what would PCO2 be?
        }
        else {
            if (kk <= nob_Input) {
                //Worksheets("Input").Cells(14, j + 15) = compositions(3, 2) * 100

                if (UseMolal == 0) {
                    //Worksheets("Input").Cells(15, j + 15) = (H2Saq + HS + S) * (34080 * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#))
                }
                else {
                    //Worksheets("Input").Cells(15, j + 15) = TH2Saq
                }
                //Worksheets("Input").Cells(16, j + 15) = pHMeterReading
                //    Worksheets("Input").Cells(17, j + 15) = compositions(2, 2) * 100  'QC check for assuming that pH and Alk are correct, what would PCO2 be?
            }

            if (kk > nob_Input && kk <= nob_Input + nob_InputII) {
                //    Worksheets("Input II").Cells(83, j + 2) = compositions(3, 2) * 100

                if (UseMolal == 0) {
                    //Worksheets("Input II").Cells(84, j + 2) = (H2Saq + HS + S) * (34080 * (rho25CMix(i + nob_Input) - CalculatedTDSMix(i + nob_Input) / 1000000#))
                }
                else {
                    //Worksheets("Input II").Cells(84, j + 2) = (H2Saq + HS + S)
                }

                //Worksheets("Input II").Cells(85, j + 2) = pHMeterReading
                //    Worksheets("Input II").Cells(86, j + 2) = compositions(2, 2) * 100 'QC check for assuming that pH and Alk are correct, what would PCO2 be?
            }
        }
    }

    //3、根据报告的pH值和STP P-CO2计算碱度。注意，H2Saq值已在上文计算过。
    *pHMeterReading = globalvar.pHMeterStpMix[kk];
    (*yCO2) = yCO2Mix[kk];

    // 计算pH和相关参数
    *pH = globalvar.pHMeterStpMix[kk] + DpHj;
    *aH = pow(10.0, - (*pH));

    H = *aH / (gCat[iH] * gNCat[iH]);
    OH = KH2O / (*aH * gAn[iOH] * gNAn[iOH]);

    // 根据useEOS选择计算CO2aq和H2Saq的方法
    if (useEOS == 0) {
        CO2aq = KgwCO2 * Ppsia * (*yCO2) * gGas[iCO2g] / (gNeut[iCO2aq] * gNNeut[iCO2aq]);
        // 注意：在useEOS=0时，H2Saq可能需要其他计算方法
    }
    else {
        // Unit Kg/ per mole of phase 3
        // 注意：VB中compositions(2,4)对应C中compositions[1][3]
        // compositions(15,4)对应compositions[14][3]
        if (compositions[14][3] != 0) { // 避免除以零
            CO2aq = compositions[1][3] / (0.01801528 * compositions[14][3]);
            H2Saq = compositions[2][3] / (0.01801528 * compositions[14][3]);
        }
        else {
            CO2aq = 0;
            H2Saq = 0;
        }
    }

    // 计算碳酸系统物种
    HCO3 = (K1H2CO3 * aH2O) * CO2aq * gNeut[iCO2aq] * gNNeut[iCO2aq] /
        (*aH * gAn[iHCO3] * gNAn[iHCO3]);
    CO3 = K2HCO3 * HCO3 * gAn[iHCO3] * gNAn[iHCO3] /
        (*aH * gAn[iCO3] * gNAn[iCO3]);

    // 计算硫化物系统物种
    HS = K1H2S * H2Saq * gNeut[iH2Saq] * gNNeut[iH2Saq] /
        (*aH * gAn[iHS] * gNAn[iHS]);

    // 计算乙酸系统物种
    double hydAc = *aH * gAn[iAc] * gNAn[iAc] / (KHAc * gNeut[iHAcaq] * gNNeut[iHAcaq]) + 1;
    AC = TAc / hydAc;

    // 计算硼酸系统物种
    double hydH2BO3 = *aH * gAn[iH2BO3] * gNAn[iH2BO3] / (KH3BO3 * gNeut[iH3BO3] * gNNeut[iH3BO3]) + 1;
    H2BO3 = TH3BO3 / hydH2BO3;

    // 计算氨系统物种
    double hydNH3 = *aH * gNeut[iNH3] * gNNeut[iNH3] / (KNH4 * gCat[iNH4] * gNCat[iNH4]) + 1;
    NH3 = TNH4 / hydNH3;

    // 计算硅酸系统物种
    double hydH2SiO4 = pow(*aH, 2) * gAn[iH2SiO4] * gNAn[iH2SiO4] /
        (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]) +
        (*aH) * gAn[iH2SiO4] * gNAn[iH2SiO4] /
        (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]) + 1;
    H2SiO4 = TH4SiO4 / hydH2SiO4;
    H3SiO4 = H2SiO4 * (*aH) * gAn[iH2SiO4] * gNAn[iH2SiO4] /
        (KH3SiO3 * gAn[iH3SiO4] * gNAn[iH3SiO4]);
    H4SiO4 = H2SiO4 * pow(*aH, 2) * gAn[iH2SiO4] * gNAn[iH2SiO4] /
        (KH4SiO4 * KH3SiO3 * gNeut[iH4SiO4aq] * gNNeut[iH4SiO4aq]);

    // 计算碱度和钠平衡
    double Alk_from_QC = (HCO3 + 2 * CO3 + HS + AC + NH3 + H2BO3 + H3SiO4 + 2 * H2SiO4 + OH - H);
    double NaQC = (-SumOfAnions - (SumOfCations - globalvar.NaMix[kk]));

    // 输出结果到相应的工作表
    if (RunStat == 1 && CaseCount[1] == 1) {
    /*
    *     Worksheets("Output Data Sheet").Cells(7, 4) = Alk_from_QC * (61019 * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#))
        Worksheets("Output Data Sheet").Cells(8, 4).Value = SumOfCations * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#) 'Convert output from equiv/kg-water to equiv/l-solution
        Worksheets("Output Data Sheet").Cells(9, 4).Value = SumOfAnions * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#) 'Convert output from equiv/kg-water to equiv/l-solution
        Worksheets("Output Data Sheet").Cells(10, 4).Value = rho25CMix(kk)
        Worksheets("Output Data Sheet").Cells(11, 4).Value = CalculatedTDSMix(kk)
        Worksheets("Output Data Sheet").Cells(12, 4).Value = NaQC * (22990 * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#))
    */
    }
    else if (RunStat == 1 && CaseCount[1] == 9) {
        /*
        Worksheets("Output Data Sheet").Cells(19, 4) = Alk_from_QC * (61019 * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#))
        Worksheets("Output Data Sheet").Cells(20, 4).Value = SumOfCations * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#) 'Convert output from equiv/kg-water to equiv/l-solution
        Worksheets("Output Data Sheet").Cells(21, 4).Value = SumOfAnions * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#) 'Convert output from equiv/kg-water to equiv/l-solution
        Worksheets("Output Data Sheet").Cells(22, 4).Value = rho25CMix(kk)
        Worksheets("Output Data Sheet").Cells(23, 4).Value = CalculatedTDSMix(kk)
        Worksheets("Output Data Sheet").Cells(24, 4).Value = NaQC * (22990 * (rho25CMix(kk) - CalculatedTDSMix(kk) / 1000000#))
        */
    }
    else if (RunH2SGUI == 1) {
        if (UseMolal == 0) {
            /*
            Worksheets(mySheet).Cells(87, j + 2) = Alk_from_QC * (61019 * (rho25c - TDS / 1000000#))
            Worksheets(mySheet).Cells(88, j + 2).Value = SumOfCations * (rho25c - TDS / 1000000#) 'Convert output from equiv/kg-water to equiv/l-solution
            Worksheets(mySheet).Cells(89, j + 2).Value = SumOfAnions * (rho25c - TDS / 1000000#) 'Convert output from equiv/kg-water to equiv/l-solution
            Worksheets(mySheet).Cells(91, j + 2).Value = NaQC * (22990 * (rho25c) - TDS / 1000000#)    
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
        if (kk <= nob_Input) {
            if (UseMolal == 0) {
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

        if (kk > nob_Input && kk <= nob_Input + nob_InputII) {
            if (UseMolal == 0) {
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
```