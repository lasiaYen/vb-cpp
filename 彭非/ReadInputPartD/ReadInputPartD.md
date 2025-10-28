# ReadInputPartD

> code

/*
-----------------------------------！！！！！可能存在的重大错误！！！！！---------------------------------

观察代码，其中kk常被用于数组访问下标。  例如：数组arr调用 arr(kk)，
此时，vb的arr(kk) 在c语言中实际为 arr[kk-1]，因为c语言的数组下标是从0开始的，
但我并没有把kk换为kk-1。 

解决方案：
法一： 我们设置传入的kk从0开始，  同时   我们要去代码中把诸如 kk==1 的代码换成 kk==0

法二：传入kk时，依旧按照vb源码的结果传入（也就是kk肯定不为0）   我们把 数组访问的 kk 换成 kk - 1（偏移offset）


暂时没有确定选哪一个，等待后续的确定方案

2025-10月15号
*/

```cpp

void ReadInputPartD(int kk, int j,
    int Run_MassTransfer, int RunShellMultiflash, double* VMeOH, double* VMEG, double* mw_oil, double* pH, double* RatioOilBPoints,
    double* YH2O, double* mol_w_Orig, int* usedryHC, double* mol_W, double Znew, double* mol_o, double* wm_oil, double* mol_g, double* mol_HC,
    bool* mf_ParametersWereRead, double* mf_TCr, double* mf_PCr, double* mf_Omega, double* mf_MWgas, double** mf_kPr, 
    double* mf_c0, double* mf_c1, double** mass_phase, double** MW_Phase, int* No_Phases,
    double** MWgas, double** TCr, double** PCr, double** Omega, double*** kPr,
    double** composition_G, double** lnphi_Gas, double** composition_Aq, double** lnphi_Aq, double* Compr_composition_Aq, double** mf_reservoir_Composition,
    double** mf_feed_Composition, int* NumGases,
    double* rhoTP, double* yCO2, double* yH2S, double* yCH4, int RunQualityControlChecks_II, double* aH, 
    double K1H2S, double KgwH2S, double Ppsia, double KstFeSaq, double TZn, double TPb, double hydHS, double ppt, double* BetaDot, double* gDot,
    double KH2O, double KHAc, double KH3BO3, double KNH4, double KH4SiO4, double KH3SiO3, double K1H2CO3, double KgwCO2, double K2HCO3,
    double* pHMeterReading
) {
    
    int mol_w3;

    double feed_Composition[15];    

    double tempgNeut[2] = { gNeut[iCO2aq], gNeut[iH2Saq] };

    //-----------------------------------------------------
	int iNG = 0;
	for (; iNG < 15; iNG++) {
		zOutput[iNG] = 0;
		z[iNG] = 0;
	}
	for (; iNG < 2; iNG++)
		globalvar.density[iNG] = 0;

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

		if (RunShellMultiflash == 1) globalvar.useEOSmix[kk] = 0;

		/*
		* 待定，注：useEOSmix是整形的，为什么会用""来判断？？？？？
		if (useEOSmix[kk] == "") useEOSmix(kk) = 0
		*/
		if (RunMultiMix == 1) globalvar.useEOSmix[kk] = 0;
		
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
	mc[iH] = HstpMix[kk]; mc[iNa] = globalvar.NaMix[kk]; mc[iK] = globalvar.KMix[kk]; mc[iMg] = globalvar.MgMix[kk]; mc[iCa] = globalvar.CaMix[kk];
	int TCa = mc[iCa];
	mc[iSr] = globalvar.SrMix[kk]; mc[iBa] = globalvar.BaMix[kk]; mc[iFe] = globalvar.FeMix[kk]; mc[iZn] = globalvar.ZnMix[kk]; mc[iPb] = globalvar.PbMix[kk]; mc[iRa] = globalvar.RaMix[kk];

	ma[iOH] = OHstpMix[kk]; ma[iCl] = globalvar.ClMix[kk]; ma[iAc] = ACstpMix[kk]; mc[iNH4] = NH4STPMix[kk]; ma[iH2BO3] = H2BO3stpMix[kk];
	ma[iHCO3] = HCO3stpMix[kk]; ma[iCO3] = CO3stpMix[kk];
	ma[iSO4] = globalvar.SO4Mix[kk]; ma[iHS] = HSstpMix[kk]; ma[intF] = globalvar.FMix[kk]; ma[iBr] = globalvar.BrMix[kk];

	Alk = globalvar.AlkMix[kk]; TAc = globalvar.TAcMix[kk]; TNH4 = TNH4Mix[kk]; TH3BO3 = TH3BO3Mix[kk]; TH2Saq = TH2SaqMix[kk]; TCO2 = TCO2Mix[kk];

	VW = globalvar.VwMix[kk]; VgTP = globalvar.VgTPMix[kk]; VO = globalvar.VoMix[kk]; *VMeOH = globalvar.VMeOHMix[kk]; *VMEG = globalvar.VMEGMix[kk]; mass_MeOH = mass_MeOH_mix[kk]; mass_MEG = mass_MEG_mix[kk];

    *yCO2 = yCO2Mix[kk], *yH2S = yH2SMix[kk], *yCH4 = yCH4Mix[kk];   // Local variable values; in this loop only.

	int useEOS = globalvar.useEOSmix[kk]; int use_pH = usepHmix[kk]; int UseH2Sgas = globalvar.UseH2SgasMix[kk];

    void* mt;


	TFe = mc[iFe];

	for (iNG = 0; iNG < 14; iNG++) {
		z[iNG] = zMix[kk][iNG];
		if (useEOS == 3 || useEOS == 0)
			z[iNG] = 0;
	}

	SGG = globalvar.gasSpGravMix[kk]; API = globalvar.oilAPIgravMix[kk];

	//根据Cragoe, C.S.1929《石油产品的热力学性质》（美国商务部标准局杂项出版物第97号）计算石油的平均分子量。
	if (API > 20 && API < 80) *mw_oil = 6084.0 / (API - 5.9);
	if (API <= 20) *mw_oil = 6084 / (20 - 5.9);
	if (API >= 80) *mw_oil = 6084 / (80 - 5.9);

	//根据流体温度和压力计算水、油的质量以及水、油和气体的摩尔数

	if (globalvar.UseTPVolMix[kk] == 1) mt = fTPFunc(1);//检查 T、P 以进行 mass_W、mass_O、mol_g 计算
	if (globalvar.UseTPVolMix[kk] == 0) mt = fTPFunc(0);

	CalcIonicStrength();
	*pH = globalvar.pHMeterStpMix[kk] + DpHj;
	*RatioOilBPoints = fRatioOilBPoints();
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
	if (globalvar.UseTPVolMix[kk] == 1) WaterDensityMix[kk] = CalcRhoTP(TK, TC, PBar, Patm);

	mass_o_Mix[kk] = 159 * globalvar.VoMix[kk] * OilDensityMix[kk]; // 这些用于在 useEOS=0 和 useEOS=3 中计算 nTCO2 和 nTH2S。这些值在 useEOS=1 或 2 中重新计算。
	mass_w_Mix[kk] = 159 * globalvar.VwMix[kk] * WaterDensityMix[kk]; // 换算成公斤盐水
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

			if ((globalvar.VoMix[kk] == (1.0 / 159.0 / 1000.0)) && (VgTP == 0.000001)) { // 当气体和油都为零并混合水时
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

					if (usepHmix[kk] == 1) *pH = globalvar.pHMeterStpMix[kk] + DpHj;

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
                    aH2O, globalvar.density, compositions, phi, Compr, beta, zOutput, mass_phase, MW_Phase, No_Phases);

                if (beta[0] > 0 && beta[1] > 0) { 
                    if (globalvar.VoMix[kk] > 1.0 / 159.0 / 1000.0) { // key to oil if vo is greater than 0
                        mass_o_Mix[kk] = 159 * globalvar.VoMix[kk] * globalvar.density[1]; // Kg
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
                    if (globalvar.VoMix[kk] > 1.0 / 159.0 / 1000.0) { // 如果存在石油量，则石油的关键
                        if ((z[7] + z[8] + z[9] + z[10] + z[11] + z[12]) > 0) { // 只有当 HC 大于 C4 时，才是石油的关键
                            mass_o_Mix[kk] = 159 * globalvar.VoMix[kk] * globalvar.density[0]; // Kg
                            (*mol_o) = mass_o_Mix[kk] * 1000 / (*MW_Phase)[0];
                            (*mol_g) = 0;
                        }
                        else {
                            (*mol_g )= VgTP * PBar * 1000 / (Compr[0] * RBar * TK); // 如果不存在大于 C4 的 HC，则为气体键
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
                    if (globalvar.VoMix[kk] > 1.0 / 159.0 / 1000.0) { // key to oil if vol of oil is present
                        mass_o_Mix[kk] = 159 * globalvar.VoMix[kk] * globalvar.density[1]; // Kg
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
                        tempgNeut, aH2O, globalvar.density, compositions, phi, Compr, beta, zOutput, mass_phase, MW_Phase, No_Phases);

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
                                tempgNeut, aH2O, globalvar.density, compositions, phi, Compr, beta, zOutput, mass_phase, MW_Phase, No_Phases);
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
                    tempgNeut, aH2O, globalvar.density, compositions, phi, Compr, beta, zOutput, mass_phase, MW_Phase, No_Phases);
            }

            nTCO2MixEOS[kk] = total_moles * zOutput[1];  // 原zOutput(2)对应zOutput[1]
            nTH2SMixEOS[kk] = total_moles * zOutput[2];  // 原zOutput(3)对应zOutput[2]
            mass_w = total_moles * beta[2] * compositions[14][3] * 0.01801528; // mass_w is only the aqueous phase H2O
            Total_molesMix[kk] = total_moles;
            mass_w_Mix[kk] = mass_w;

            if (*No_Phases == 3) { // Output oil and gas density from flash calculation if useEOS=1 or useEOS=2
                GasDensityMix[kk] = globalvar.density[0] * 1000;  // 原density(1)对应density[0]
                OilDensityMix[kk] = globalvar.density[1];         // 原density(2)对应density[1]
            }

            if (*No_Phases == 2) { // Only g/o, G/W or O/W phases exist
                if (beta[2] == 0) {
                    //本例中不存在水相。增加水量并重新运行。程序中止。
                    printf("Aqueous phase does not existed in this case. Increase water volume and run again. Program abort.\n"); // Aqueous phase does not exist
                    exit(1);  // 替代VB的End
                }
                else {
                    
                    if (globalvar.density[0] > 0 && globalvar.density[0] < 0.3) GasDensityMix[kk] = globalvar.density[0] * 1000;// phase 1 is gas

                    else if (globalvar.density[0] > 0 && globalvar.density[0] > 0.3) OilDensityMix[kk] = globalvar.density[0]; // phase 1 is oil

                    else if (globalvar.density[1] > 0 && globalvar.density[1] < 0.3)  GasDensityMix[kk] = globalvar.density[1] * 1000;// Phase 2 is gas

                    else if (globalvar.density[1] > 0 && globalvar.density[1] > 0.3) OilDensityMix[kk] = globalvar.density[1];

                }
            }

            if (*No_Phases == 1) {
                if (globalvar.density[0] > 0 && compositions[14][1] < 0.8) { // hydrocarbon phase
                    printf("Aqueous phase does not existed in this case. Increase water volume and run again. Program abort.\n");
                    exit(1);
                }
                else if (globalvar.density[1] > 0 && compositions[14][2] < 0.8) {
                    printf("Aqueous phase does not existed in this case. Increase water volume and run again. Program abort.\n");
                    exit(1);
                }
                else if (globalvar.density[2] > 0 && compositions[14][3] < 0.8) {
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
                aH2O, globalvar.density, compositions, phi, Compr, beta, zOutput, mass_phase, MW_Phase, No_Phases);

            if (*usedryHC == 0) {
                mol_w3 = total_moles * beta[2] * compositions[14][3];

                while (mol_w3 <= 0) { // 确保添加足够的水以开始迭代
                    (*mol_W) = (*mol_W) * 1.1;
                    pseudo_composition(API, SGG, VgTP, *mol_o, *mol_W, TK, PBar, aH2O, tempgNeut, nTCO2, nTH2S, *yCO2, *yH2S, *YH2O, &total_moles, feed_Composition, mol_HC,
                        MWgas, TCr, PCr, Omega, &mf_c0, &mf_c1, kPr,
                        composition_G, lnphi_Gas, composition_Aq, lnphi_Aq, Compr_composition_Aq, mf_reservoir_Composition, mf_feed_Composition, NumGases);
                    MultiPhaseFlash(mf_ParametersWereRead, mf_TCr, mf_PCr, mf_Omega, mf_MWgas, mf_kPr, mf_c0, mf_c1, TK, PBar, total_moles, feed_Composition, 
                        tempgNeut, aH2O, globalvar.density, compositions, phi, Compr, beta, zOutput, mass_phase, MW_Phase, No_Phases);

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
                            aH2O, globalvar.density, compositions, phi, Compr, beta, zOutput, mass_phase, MW_Phase, No_Phases);
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
            globalvar.NaMix[kk] = mc[iNa];
            globalvar.KMix[kk] = mc[iK];
            globalvar.MgMix[kk] = mc[iMg];
            globalvar.CaMix[kk] = mc[iCa];
            TCa = mc[iCa];
            globalvar.SrMix[kk] = mc[iSr];
            globalvar.BaMix[kk] = mc[iBa];
            globalvar.ZnMix[kk] = mc[iZn];
            globalvar.PbMix[kk] = mc[iPb];
            globalvar.RaMix[kk] = mc[iRa];

            OHstpMix[kk] = ma[iOH];
            globalvar.ClMix[kk] = ma[iCl];
            ACstpMix[kk] = ma[iAc];
            NH4STPMix[kk] = mc[iNH4];
            H2BO3stpMix[kk] = ma[iH2BO3];
            HCO3stpMix[kk] = ma[iHCO3];
            CO3stpMix[kk] = ma[iCO3];
            globalvar.SO4Mix[kk] = ma[iSO4];
            HSstpMix[kk] = ma[iHS];
            globalvar.FMix[kk] = ma[intF];
            globalvar.BrMix[kk] = ma[iBr];

            rho25CMix[kk] = rho25c;

            H3SiO4Mix[kk] = ma[iH3SiO4];
            H2SiO4Mix[kk] = ma[iH2SiO4];

            NH3Mix[kk] = mn[iNH3];
            H4SiO4Mix[kk] = mn[iH4SiO4aq];
            H3BO3Mix[kk] = mn[iH3BO3];

            globalvar.CO2aqMix[kk] = mn[iCO2aq];
            globalvar.H2SaqMix[kk] = mn[iH2Saq];
            HACaqMix[kk] = mn[iHAcaq];

            globalvar.AlkMix[kk] = Alk;
            globalvar.TAcMix[kk] = TAc;
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

```