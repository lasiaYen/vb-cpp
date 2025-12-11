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
    simContext.nob = 0;
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
    if (useEOS > 0) {

        nTCO2 = nTCO2EOS;
        nTH2S = nTH2sEOS;

        if (total_moles > 0) {
            for (iNG = 0; iNG < 15; iNG++) {
                z[iNG] /= total_moles;
                //可能会导致c和vb行为不同的代码
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

        fmn;

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

            fmn;
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
        fmn;

        C2_PitzerActCoefs_T_P_ISt(gNeut, aH2O, TK, TC, PBar, Patm);
        //Call C4_SSPEquilCalcs(0, 5, 2, KspCalcite)
        C4_SSPEquilCalcs(0, 4, 2, KspCalcite);

        /* Recompute API + SGG */
        fTPFunc(3);

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

    fmn;
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
        fmn;

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
        fmn;
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
        // 需要先定义 SICal 变量
        SICal = 0.0; // 假设初始值

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

        SIDiopside = log10((mc[iCa] * mc[iMg] * pow(H4SiO4, 2) * gCat[iCa] *
            gCat[iMg] * pow(gNeut[iH4SiO4aq], 2) / pow(aH, 4)) / KspDiopside);

        SIChrysotile = log10((pow(mc[iMg], 3) * pow(H4SiO4, 2) *
            pow(gCat[iMg], 3) * pow(gNeut[iH4SiO4aq], 2) / pow(aH, 6)) / KspChrysotile);

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


void B6_InhibitorNeeded()
{
    /* ---------- Barite ---------- */
    double fSafetyBar = 1.0;
    ConcInhBar = 0.0;

    if (mc[iBa] * ma[iSO4] > 0.0)
    {
        if (LoopTP == 1 || RunWhatIf == 1 || Run1000Cases == 1 || LoopResChem == 1)
        {
            InhNoBar = InhNo;

            if (SelectInh == 1)
            {
                int iMaxBar = 1;
                double bInhBarMax = fbInhBar(iMaxBar, SIBar);

                for (int i = 1; i <= 19; i++)
                {
                    bInhBar[i] = fbInhBar(i, SIBar);
                    if (bInhBar[i] > bInhBarMax)
                    {
                        iMaxBar = i;
                        bInhBarMax = bInhBar[i];
                    }
                }
                InhNoBar = iMaxBar;
            }

            if (SelectInh == 1)
                InhNoBar = 2;   /* VB: always choose BHPMP */

            /* Excel writes → remove → optional store */
            InhNameSelected_Barite = InhName[InhNoBar];
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
    double fSafetyCal = 1.0;
    ConcInhCal = 0.0;

    if (mc[iCa] * HCO3 > 0.0)
    {
        if (LoopTP == 1 || RunWhatIf == 1 || Run1000Cases == 1 || LoopResChem == 1)
        {
            InhNoCal = InhNo;

            if (SelectInh == 1)
            {
                int iMaxCal = 1;
                double bInhCalMax = fbInhCal(iMaxCal, SICal);

                for (int i = 1; i <= 19; i++)
                {
                    bInhCal[i] = fbInhCal(i, SICal);
                    if (bInhCal[i] > bInhCalMax)
                    {
                        iMaxCal = i;
                        bInhCalMax = bInhCal[i];
                    }
                }
                InhNoCal = iMaxCal;
            }

            if (SelectInh == 1)
                InhNoCal = 1;   /* always choose NTMP */

            InhNameSelected_Calcite = InhName[InhNoCal];
        }

        if (SICal > 0.001)
        {
            double CalExpon10 = flogT0Cal(SICal);

            if (log10(tInh) > CalExpon10 && CalExpon10 < 8.0)
            {
                t0Cal = pow(10.0, CalExpon10);
                ConcInhCal = fCinhCal(SICal, tInh);
            }
        }
    }


    /* ---------- Gypsum ---------- */
    double fSafetyGyp = 1.0;
    ConcInhGyp = 0.0;

    if (mc[iCa] * ma[iSO4] > 0.0)
    {
        InhNoGyp = 4;

        InhNameSelected_Gypsum = InhName[InhNoGyp];

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
                    ConcInhGyp = NAN;  /* VB Null */
                }
            }
        }
    }


    /* ---------- Anhydrite ---------- */
    double fSafetyAn = 1.0;
    ConcInhAn = 0.0;

    if (mc[iCa] * ma[iSO4] > 0.0)
    {
        InhNoAn = 4;
        InhNameSelected_Anhydrite = InhName[InhNoAn];

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
                {
                    ConcInhAn = NAN;
                }
            }
        }
    }


    /* ---------- Celestite (SrSO4) ---------- */
    double fSafetyCel = 1.0;
    ConcInhCel = 0.0;

    if (mc[iSr] * ma[iSO4] > 0.0)
    {
        if (InhNo <= 11)
            InhNoCel = 3;      /* phosphonate → DTPMP */
        else if (InhNo <= 14)
            InhNoCel = 12;     /* carboxylates → PPCA */
        else
            InhNoCel = 17;     /* others → PVS */

        if (SelectInh == 1)
            InhNoCel = 3;      /* always choose DTPMP */

        InhNameSelected_Celestite = InhName[InhNoCel];

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


int A1_Start_ScaleSoftPitzer()
{
    SampleData data;
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

        for (LoopTP = 0; LoopTP < 10; LoopTP++)
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

            LoopTPSI();

            if (H2Oevap != 1) {

                if (ISt != 0) {
                    LoopTPVisHeatCap();

                    if (Run_MassTransfer == 1)
                        MassTransferCoefficients();
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

            LoopTPSI();
            LoopTPWrite();

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

        ErrorMsgBox();

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
        return;  //End
    }
}