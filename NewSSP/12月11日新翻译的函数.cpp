void InitializeOptionClearCellContent()
{

    UseSR = 0; simContext.UseTPCalciteSheet = 0;
    simContext.Read_InputII = 0;
    //int NCount_II = 0; 
    useEOS = 0;
    simContext.LoopMixing = 0; simContext.UseMolal = 0; //iTP = 0
    RunShellMultiflash = 0; H2Oevap = 0;
    simContext.Run_CalcConcFactor = 0;

    //int RunGoalSeek, RunStatGoalSeek; Run_MassTransfer;
    //RunH2SGUI, RunSMT, RunStat, RunH2SPartition;

    RunGoalSeek = (RunGoalSeek == 1 ? 1 : 0);
    RunStatGoalSeek = (RunStatGoalSeek == 1 ? 1 : 0);
    Run_MassTransfer = (Run_MassTransfer == 1 ? 1 : 0);

    // If Worksheets("Input").Range("I11") = "Saturation Ratio values" Then UseSR = 1;


    usePTB = 1; 


    /*
    if (simContext.RunH2SGUI != 1 && RunSMT != 1 && simContext.RunStat != 1) {
         if Worksheets("Calcite").Range("F3") == "mg/L"
             usePTB = 0;
         else
             usePTB = 1;
    }
    */
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

    //double feed_Composition[15];
    DpHj = 0;

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
    simContext.nob = 0;
    Ncount = 0;
    NCount_II = 0;
    //simContext.nob_Input = 0;
    //simContext.nob_InputII = 0;
    simContext.Run1000Cases = 0;
    RunStatReservoirCalc = 0;


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
    else
    {
        for (int i = 0; i < 5; i++)
        {

            // if Worksheets("Input").Cells(3, 2+i+1) = True
                //nob++;
                //CaseCount[i - Ncount] = i;
            //else   
                //Ncount++;
        }

        //nob_Input = nob;
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



void B3_CalcConcs(double& API)
{
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
                            double aaa123 = 0;
    }
label100:
    return;
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
        SIZnS = log10(mc[iZn] * simContext.HS * gCat[iZn] * gAn[iHS] * gNAn[iHS] / aH / KspZnS);

        SIPbS = log10(mc[iPb] * simContext.HS * gCat[iPb] * gAn[iHS] * gNAn[iHS] / aH / KspPbS);

        SIZnCO3 = log10(mc[iZn] * simContext.CO3 * gCat[iZn] * gAn[iCO3] / KspZnCO3);
        SICal = 0.0; 

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



double fbInhBar(int InhNo, double SI)
{
    return bi[InhNo - 1][0]
        + bi[InhNo - 1][1] * SI
        + bi[InhNo - 1][2] / TK
        + bi[InhNo - 1][3] * log10(1.0 / aH)
        + bi[InhNo - 1][4] * fabs(log10(mc[iBa] / ma[iSO4]));
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
        ci[InhNo - 1][0]                         // ci(InhNo,1)
        + ci[InhNo - 1][1] * SI                  // ci(InhNo,2) * SI
        + ci[InhNo - 1][2] / TK                  // ci(InhNo,3) / TK
        + ci[InhNo - 1][3] * log10(1.0 / aH)     // ci(InhNo,4) * Log10(1/aH)
        + ci[InhNo - 1][4]                       // ci(InhNo,5)
        * fabs(log10(mc[iCa] / simContext.HCO3));     // Abs(Log10(mc(iCa)/HCO3))
}


double fCinhCal(double SI, double tInh, double& fSafetyCal)
{
    double bInhCalMixed;

    if (InhNoCal == 20)
    {
        // 注意：所有 VB 下标均需 -1
        bInhCal[simContext.InhNo1 - 1] = pow(10.0, fbInhCal(simContext.InhNo1, SI));
        bInhCal[simContext.InhNo2 - 1] = pow(10.0, fbInhCal(simContext.InhNo2, SI));

        bInhCalMixed = simContext.FracInhNo1 * bInhCal[simContext.InhNo1 - 1]
            + (1.0 - simContext.FracInhNo1) * bInhCal[simContext.InhNo2 - 1];

        return (1.0 / bInhCalMixed)
            * log10(fSafetyCal * tInh / t0Cal);
    }
    else
    {
        bInhCal[InhNoCal - 1] = pow(10.0, fbInhCal(InhNoCal, SI));

        return (1.0 / bInhCal[InhNoCal - 1])
            * log10(fSafetyCal * tInh / t0Cal);
    }
}

double fCinhBar(double SI, double tInh, double& fSafetyBar)
{
    double bInhBarMixed;

    if (InhNoBar == 20)
    {
        bInhBar[simContext.InhNo1 - 1] = pow(10.0, fbInhBar(simContext.InhNo1, SI));

        bInhBar[simContext.InhNo2 - 1] = pow(10.0, fbInhBar(simContext.InhNo2, SI));

        bInhBarMixed =
            simContext.FracInhNo1 * bInhBar[simContext.InhNo1 - 1] +
            (1.0 - simContext.FracInhNo1) * bInhBar[simContext.InhNo2 - 1];

        return (1.0 / bInhBarMixed) * log10(fSafetyBar * tInh / t0Bar);
    }
    else
    {
        bInhBar[InhNoBar - 1] = pow(10.0, fbInhBar(InhNoBar, SI));

        return (1.0 / bInhBar[InhNoBar - 1]) * log10(fSafetyBar * tInh / t0Bar);
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

double fCinhGyp(double SI, double tInh, double& fSafetyGyp, double& t0Gyp)
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


double fCinhAn(double SI, double tInh, double& fsafetyAn, double& t0An)
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


double fCinhCel(double SI, double tInh, double& fSafetyCel, double& t0Cel)
{
    bInhCel[InhNoCel - 1] = pow(10, fbInhCel(InhNoCel, SI));
    return (1.0 / bInhCel[InhNoCel - 1]) * log10(fSafetyCel * tInh / t0Cel);
}



void B6_InhibitorNeeded()
{
    /* ---------- Barite ---------- */
    double fSafetyBar = 1.0;
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

            InhNameSelected_Barite = InhName[InhNoBar];
        }

        if (SIBar > 0.001)
        {
            double BarExpon10 = flogT0Bar(SIBar);

            if (log10(tInh) > BarExpon10 && BarExpon10 < 8.0)
            {
                t0Bar = pow(10.0, BarExpon10);
                ConcInhBar = fCinhBar(SIBar, tInh, fSafetyBar);
            }
        }
    }


    /* ---------- Calcite ---------- */
    double fSafetyCal = 1.0;
    ConcInhCal = 0.0;
    double ConcInhGyp;

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
                ConcInhCal = fCinhCal(SICal, tInh, fSafetyCal);
            }
        }
    }


    /* ---------- Gypsum ---------- */
    double fSafetyGyp = 1.0;
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
                    double t0Gyp = pow(10.0, GypExpon10);
                    ConcInhGyp = fCinhGyp(SIGyp, tInh, fSafetyGyp, t0Gyp);
                }
                else
                {
                    ConcInhGyp = 0;  /* VB Null */
                }
            }
        }
    }


    /* ---------- Anhydrite ---------- */
    double fsafetyAn = 1.0;
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
                    double t0An = pow(10.0, AnExpon10);
                    ConcInhAn = fCinhAn(SIAn, tInh, fsafetyAn, t0An);
                }
                else
                    ConcInhAn = 0;
            }
        }
    }

    /* ---------- Celestite (SrSO4) ---------- */
    double fSafetyCel = 1.0;
    ConcInhCel = 0.0;

    if (mc[iSr] * ma[iSO4] > 0.0)
    {
        if (simContext.InhNo <= 11)
            InhNoCel = 3;      /* phosphonate → DTPMP */
        else if (simContext.InhNo <= 14)
            InhNoCel = 12;     /* carboxylates → PPCA */
        else
            InhNoCel = 17;     

        if (simContext.SelectInh == 1)
            InhNoCel = 3;      /* always choose DTPMP */

        //Worksheets("Input").Cells(55, 10) = InhName(InhNoCel)

        if (SICel > 0.001)
        {
            double CelExpon10 = flogT0Cel(SICel);

            if (log10(tInh) > CelExpon10 && CelExpon10 < 8.0)
            {
                double t0Cel = pow(10.0, CelExpon10);
                ConcInhCel = fCinhCel(SICel, tInh, fSafetyCel, t0Cel);
            }
        }
    }
}


void LoopTPSI(double& API)
{
    B5_CalculateSIvalues(API);
    pH_before_precipitation = pH;  

    if (H2Oevap != 1)
        B6_InhibitorNeeded();   

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

    // ----- ���� ��SI -----
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

    // ----- ȫ������������ -----
=======

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
   // SI ��λת���߼�
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

            LoopTPSI(data.API);

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

            LoopTPSI(data.API);
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

void B7_ScaleRisk()
{
    // 缁撴瀯浣撳叏灞�鍙橀噺
    // double mc[], ma[], TK, tInh, xMeOH, xMEG, HCO3;
    // double SIRisk[ ... ];
    // double ConcInhBarRisk[], ConcInhCalRisk[], ConcInhAnRisk[];
    // int InhNoBar, InhNoCal, InhNoAn;
    // int InhNo1, InhNo2;
    // int LoopTP, NoRiskcalc;

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

        for (int iRisk = 1; iRisk <= NoRiskcalc - 1; iRisk++)
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

                double Dep1 = fCinhBar(SIBarRisk, tInhRisk) - ConcInhBarRisk[iRisk];

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
    if (mc[iCa] > 1e-7 && HCO3 > 1e-7)
    {
        double CalExpon10 = log10(tInh);
        double b = 1876.4 + (CalExpon10 - 4.22) * TK;
        double cc = (13.8 * TK - 6259.6);

        double root2 = -cc / b;
        SIRisk[InhNoCal][LoopTP][1][2] = root2;

        if (SIRisk[InhNoCal][LoopTP][1][2] < 0)
            SIRisk[InhNoCal][LoopTP][1][2] = 0;

        double tInhRisk = tInh;

        for (int iRisk = 1; iRisk <= NoRiskcalc - 1; iRisk++)
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

                double Dep1 = fCinhCal(SICalRisk, tInhRisk) - ConcInhCalRisk[iRisk];
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
                if (InhNo1 > 5 && InhNo1 < 10)
                    SIRisk[InhNoCal][LoopTP][iRisk + 1][2] = NAN;

                if (InhNo2 > 5 && InhNo2 < 10)
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

            for (int iRisk = 1; iRisk <= NoRiskcalc - 1; iRisk++)
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

                    double Dep1 = fCinhAn(SIAnRisk, tInhRisk) - ConcInhAnRisk[iRisk];
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


// VB原始代码参考：
// Sub LoopTPWrite()
// ... (完整VB代码如上)
// End Sub

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <memory>

// ==================== 全局变量声明（根据VB代码中的变量）====================
// 注意：这些变量需要在实际项目中正确定义和初始化
int RunStat = 0;
int RunStatSICalcSSP = 0;
int RunStatGoalSeek = 0;
int RunGoalSeek = 0;
int LoopTP = 0;
int StatTPnob = 0;
int UseSI = 0;
int Run_MassTransfer = 0;
int Run_MassTransfer_WhatIf = 0;
int Iter_MT_WI = 0;
int usePTB = 0;
int useEOS = 0;
int OutPutActCoefs = 0;
int QPBubblePt = 0;
int useSR = 0;
int H2Oevap = 0;
int RunSimpSSP = 0;
int InhNoCal = 0;
int InhNoBar = 0;
int InhNoAn = 0;
int InhNoCel = 0;
int InhNo1 = 0;
int InhNo2 = 0;
int MaxInh = 0;
int nob = 0;
int Use_Corr_in_Deposition = 0;
int Flow_Pattern = 0;
int Flow_Regime = 0;
int FlowType = 0;

double location = 0.0;
double depth = 0.0;
double TF = 0.0;
double Ppsia = 0.0;
double PipeL = 0.0;
double pH_before_precipitation = 0.0;
double DpHj = 0.0;
double ConcInhCal = 0.0;
double ConcInhBar = 0.0;
double rhoTP = 0.0;
double ISt = 0.0;
double rho25c = 0.0;
double TDS = 0.0;
double TDSHalite = 0.0;
double Rho25cHalite = 0.0;
double pHaftercalciteppt = 0.0;
double pHafterAmsilicappt = 0.0;
double pHafterQuartzppt = 0.0;
double pHafterGreenaliteppt = 0.0;
double pHafterDiopsideppt = 0.0;
double pHafterChrysotileppt = 0.0;
double pHafterMgOH2ppt = 0.0;
double pHafterCaOH2ppt = 0.0;
double pptCalcite_NoMassTransfer = 0.0;
double pptBarite_NoMassTransfer = 0.0;
double pptGyp = 0.0;
double pptHemi = 0.0;
double pptAn = 0.0;
double pptCel = 0.0;
double pptHal = 0.0;
double pptFeSAm = 0.0;
double pptFeS_NoMassTransfer = 0.0;
double pptTrot = 0.0;
double pptZnS = 0.0;
double pptCaF2 = 0.0;
double pptFeCO3_NoMassTransfer = 0.0;
double pptZnCO3 = 0.0;
double pptPbS = 0.0;
double pptSrCO3 = 0.0;
double pptBaCO3 = 0.0;
double pptAmSilica = 0.0;
double pptQuartz = 0.0;
double pptChrysotile = 0.0;
double pptDiopside = 0.0;
double pptGreenalite = 0.0;
double pptMgOH2 = 0.0;
double pptCaOH2 = 0.0;
double pptCalcite_MassTransfer = 0.0;
double pptBarite_MassTransfer = 0.0;
double pptFeCO3_MassTransfer = 0.0;
double pptCalcite_MassTransfer_V = 0.0;
double pptBarite_MassTransfer_V = 0.0;
double pptFeCO3_MassTransfer_V = 0.0;
double m_CR_selected = 0.0;
double PBubblePt = 0.0;
double SICal = 0.0;
double dSICal = 0.0;
double SIDol = 0.0;
double SIBar = 0.0;
double dSIBar = 0.0;
double SIBarBH = 0.0;
double SIHal = 0.0;
double dSIHal = 0.0;
double SIGyp = 0.0;
double SIHemi = 0.0;
double SIAn = 0.0;
double SICel = 0.0;
double SIFeSAm = 0.0;
double SIFeS = 0.0;
double SITrot = 0.0;
double SIZnS = 0.0;
double SICaF2 = 0.0;
double SISid = 0.0;
double dSISid = 0.0;
double SIZnCO3 = 0.0;
double SIPbS = 0.0;
double SISrCO3 = 0.0;
double SIBaCO3 = 0.0;
double SICaOH2 = 0.0;
double SIMgOH2 = 0.0;
double SIAmSilica = 0.0;
double SIQuartz = 0.0;
double SIChrysotile = 0.0;
double SIDiopside = 0.0;
double SIGreenalite = 0.0;
double ReNO = 0.0;
double Hold_l = 0.0;
double Time_lsl_burst = 0.0;
double lgt_calcite_MT = 0.0;
double lgt_barite_MT = 0.0;
double ViscWatIst = 0.0;
double CpPerMl = 0.0;
double aH2O = 0.0;
double TK = 0.0;
double HCO3 = 0.0;
double mass_MeOH = 0.0;
double mass_MEG = 0.0;

// 化学物质浓度数组
std::vector<double> mc;  // 阳离子浓度
std::vector<double> ma;  // 阴离子浓度

// 活度系数数组
std::vector<double> gCat;   // 阳离子活度系数
std::vector<double> gAn;    // 阴离子活度系数
std::vector<double> gNeut;  // 中性分子活度系数
std::vector<double> gGas;   // 气体活度系数

// 三维数组（假设的结构）
// SIRisk[抑制剂编号][LoopTP][参数索引][物质类型]
std::vector<std::vector<std::vector<std::vector<double>>>> SIRisk;

// 井名数组
std::vector<std::string> WellNameMix;

// 工作表名称
std::string myname;

// 传质系数数组
std::vector<double> km;

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

// ==================== 主转换函数 ====================
void LoopTPWrite(ExcelWriter* excelWriter = nullptr) {
    // 创建数据存储对象
    OutputData outputData;
    
    // 第一部分：条件写入 location 和 depth
    if (RunStat == 1) {
        if (RunStatSICalcSSP == 1) {
            // VB: Worksheets("Halite analysis").Cells(6 + StatTPnob - LoopTP, 2) = location
            // VB: Worksheets("Halite analysis").Cells(6 + StatTPnob - LoopTP, 3) = depth
            int row = calculateRowWithStat(6, StatTPnob, LoopTP);
            outputData.addData("Halite analysis", row, 2, location);
            outputData.addData("Halite analysis", row, 3, depth);
            
        } else if (RunStatSICalcSSP == 2) {
            if (RunStatGoalSeek == 1) {
                // VB: Worksheets("Output data sheet").Cells(11 + StatTPnob - LoopTP, 16) = location
                // VB: Worksheets("Output data sheet").Cells(11 + StatTPnob - LoopTP, 17) = depth
                int row = calculateRowWithStat(11, StatTPnob, LoopTP);
                outputData.addData("Output data sheet", row, 16, location);
                outputData.addData("Output data sheet", row, 17, depth);
            } else {
                // VB: Worksheets("Halite analysis").Cells(6 + StatTPnob - LoopTP, 10) = location
                // VB: Worksheets("Halite analysis").Cells(6 + StatTPnob - LoopTP, 11) = depth
                int row = calculateRowWithStat(6, StatTPnob, LoopTP);
                outputData.addData("Halite analysis", row, 10, location);
                outputData.addData("Halite analysis", row, 11, depth);
            }
        } else if (RunStatSICalcSSP == 3) {
            // VB: Worksheets("Output data sheet").Cells(6 + StatTPnob - LoopTP, 8) = location
            // VB: Worksheets("Output data sheet").Cells(6 + StatTPnob - LoopTP, 9) = depth
            int row = calculateRowWithStat(6, StatTPnob, LoopTP);
            outputData.addData("Output data sheet", row, 8, location);
            outputData.addData("Output data sheet", row, 9, depth);
        }
    }
    
    // 第二部分：当 RunGoalSeek 不为1时的输出
    if (RunGoalSeek != 1) {
        if (UseSI == 0) {
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
                outputData.addData("Deposition Prediction", 35 - LoopTP, 6, static_cast<double>(Flow_Pattern));
                outputData.addData("Deposition Prediction", 35 - LoopTP, 7, static_cast<double>(Flow_Regime));
                outputData.addData("Deposition Prediction", 35 - LoopTP, 8, ReNO);
                outputData.addData("Deposition Prediction", 35 - LoopTP, 9, Hold_l);
                outputData.addData("Deposition Prediction", 35 - LoopTP, 10, Time_lsl_burst);
                outputData.addData("Deposition Prediction", 35 - LoopTP, 11, lgt_calcite_MT);
                outputData.addData("Deposition Prediction", 35 - LoopTP, 12, lgt_barite_MT);
            }
            
            // Dai 2020 deposition - Run_MassTransfer_WhatIf = 1
            if (Run_MassTransfer_WhatIf == 1) {
                double distance_ft = (LoopTP - 1) * PipeL * 0.0328084; // 转换为英尺
                int row_offset = (Iter_MT_WI - 1) * 32;
                
                outputData.addData("Deposition Prediction_WhatIf", 35 - LoopTP + row_offset, 3, distance_ft);
                outputData.addData("Deposition Prediction_WhatIf", 35 - LoopTP + row_offset, 4, TF);
                outputData.addData("Deposition Prediction_WhatIf", 35 - LoopTP + row_offset, 5, Ppsia);
                
                outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 3, distance_ft);
                outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 4, TF);
                outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 5, Ppsia);
                
                // 输出流动参数
                outputData.addData("Deposition Prediction_WhatIf", 35 - LoopTP + row_offset, 6, static_cast<double>(Flow_Pattern));
                outputData.addData("Deposition Prediction_WhatIf", 35 - LoopTP + row_offset, 7, static_cast<double>(Flow_Regime));
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
            if (RunStat == 1) {
                if (RunStatSICalcSSP == 1) {
                    int row = calculateRowWithStat(6, StatTPnob, LoopTP);
                    outputData.addData("Halite analysis", row, 4, TF);
                    outputData.addData("Halite analysis", row, 5, Ppsia);
                } else if (RunStatSICalcSSP == 2) {
                    if (RunStatGoalSeek == 1) {
                        int row = calculateRowWithStat(11, StatTPnob, LoopTP);
                        outputData.addData(myname, row, 18, TF);
                        outputData.addData(myname, row, 19, Ppsia);
                    } else {
                        int row = calculateRowWithStat(6, StatTPnob, LoopTP);
                        outputData.addData("Halite analysis", row, 12, TF);
                        outputData.addData("Halite analysis", row, 13, Ppsia);
                    }
                } else if (RunStatSICalcSSP == 3) {
                    int row = calculateRowWithStat(6, StatTPnob, LoopTP);
                    outputData.addData("Output data sheet", row, 10, TF);
                    outputData.addData("Output data sheet", row, 11, Ppsia);
                }
            }
            
        } else if (UseSI == 1) {
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
                outputData.addData("Deposition Prediction", 35 - LoopTP, 6, static_cast<double>(Flow_Pattern));
                outputData.addData("Deposition Prediction", 35 - LoopTP, 7, static_cast<double>(Flow_Regime));
                outputData.addData("Deposition Prediction", 35 - LoopTP, 8, ReNO);
                outputData.addData("Deposition Prediction", 35 - LoopTP, 9, Hold_l);
                outputData.addData("Deposition Prediction", 35 - LoopTP, 10, Time_lsl_burst);
            }
            
            // Dai 2020 deposition - Run_MassTransfer_WhatIf = 1 (SI单位)
            if (Run_MassTransfer_WhatIf == 1) {
                double distance_m = (LoopTP - 1) * PipeL / 100.0; // 转换为米
                int row_offset = (Iter_MT_WI - 1) * 32;
                
                outputData.addData("Deposition Prediction_WhatIf", 35 - LoopTP + row_offset, 3, distance_m);
                outputData.addData("Deposition Prediction_WhatIf", 35 - LoopTP + row_offset, 4, TC);
                outputData.addData("Deposition Prediction_WhatIf", 35 - LoopTP + row_offset, 5, Pbar);
                
                outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 3, distance_m);
                outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 4, TC);
                outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 5, Pbar);
                
                // 输出流动参数
                outputData.addData("Deposition Prediction_WhatIf", 35 - LoopTP + row_offset, 6, static_cast<double>(Flow_Pattern));
                outputData.addData("Deposition Prediction_WhatIf", 35 - LoopTP + row_offset, 7, static_cast<double>(Flow_Regime));
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
            if (RunStat == 1) {
                if (RunStatSICalcSSP == 1) {
                    int row = calculateRowWithStat(6, StatTPnob, LoopTP);
                    outputData.addData("Halite analysis", row, 4, TC);
                    outputData.addData("Halite analysis", row, 5, Pbar);
                } else if (RunStatSICalcSSP == 2) {
                    if (RunStatGoalSeek == 1) {
                        int row = calculateRowWithStat(11, StatTPnob, LoopTP);
                        outputData.addData(myname, row, 18, TC);
                        outputData.addData(myname, row, 19, Pbar);
                    } else {
                        int row = calculateRowWithStat(6, StatTPnob, LoopTP);
                        outputData.addData("Halite analysis", row, 12, TC);
                        outputData.addData("Halite analysis", row, 13, Pbar);
                    }
                } else if (RunStatSICalcSSP == 3) {
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
            if (InhNo1 == 5 || InhNo1 == 7 || InhNo1 == 9 || InhNo1 == 8 || 
                InhNo1 == 10 || InhNo1 == 14 || InhNo1 == 18) {
                ConcInhCal = std::numeric_limits<double>::quiet_NaN();
            }
            if (InhNo2 == 5 || InhNo2 == 7 || InhNo2 == 8 || InhNo2 == 9 || 
                InhNo2 == 10 || InhNo2 == 11 || InhNo2 == 14 || InhNo2 == 18) {
                ConcInhCal = std::numeric_limits<double>::quiet_NaN();
            }
        }
        
        // 输出抑制剂浓度
        outputData.addData("Calcite", 14 - LoopTP, 7, ConcInhCal);
        outputData.addData("Calcite", 14 - LoopTP, 8, pHaftercalciteppt - DpHj);
        outputData.addData("Calcite", 14 - LoopTP, 11, ISt * (rho25c - TDS / 1000000.0));
        
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
        double density_factor = rho25c - TDS * 0.000001;
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
            
        } else if (usePTB == 1) {
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
        if (RunStat == 1) {
            double pptHal_output = 0.0;
            if (usePTB == 0) {
                pptHal_output = pptHal * 58443.0 * halite_density_factor;
            } else {
                pptHal_output = pptHal * 58443.0 * halite_density_factor * 0.35051;
            }
            
            if (RunStatSICalcSSP == 1) {
                int row = calculateRowWithStat(6, StatTPnob, LoopTP);
                outputData.addData("Halite analysis", row, 7, pptHal_output);
                outputData.addData("Halite analysis", row, 8, rhoTP);
            } else if (RunStatSICalcSSP == 2) {
                if (RunStatGoalSeek == 1) {
                    int row = calculateRowWithStat(11, StatTPnob, LoopTP);
                    outputData.addData(myname, row, 22, rhoTP);
                    outputData.addData(myname, row, 21, pptHal_output);
                } else {
                    int row = calculateRowWithStat(6, StatTPnob, LoopTP);
                    outputData.addData("Halite analysis", row, 15, pptHal_output);
                    outputData.addData("Halite analysis", row, 16, rhoTP);
                }
            } else if (RunStatSICalcSSP == 3) {
                int row = calculateRowWithStat(6, StatTPnob, LoopTP);
                outputData.addData("Output data sheet", row, 13, pptHal_output);
                outputData.addData("Output data sheet", row, 14, rhoTP);
            }
        }
        
        // ==================== 气泡点检测输出 ====================
        if (useEOS == 0) {
            std::string bubbleStatus = (Ppsia > PBubblePt) ? "Yes" : "No";
            outputData.addData("Calcite", 14 - LoopTP, 9, bubbleStatus);
        } else {
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
        } else {
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
            } else {
                outputData.addData("Deposition Prediction", 51 - LoopTP, 7,
                    pptCalcite_NoMassTransfer * 100091.0 * density_factor * 0.35051);
                outputData.addData("Deposition Prediction", 51 - LoopTP, 11,
                    pptBarite_NoMassTransfer * 233390.0 * density_factor * 0.35051);
                outputData.addData("Deposition Prediction", 51 - LoopTP, 15,
                    pptFeCO3_NoMassTransfer * 233390.0 * density_factor * 0.35051);
            }
        }
        
        // WhatIf 情景
        if (Run_MassTransfer_WhatIf == 1) {
            if (usePTB == 0) {
                outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + (Iter_MT_WI - 1) * 32, 7,
                    pptCalcite_NoMassTransfer * 100091.0 * density_factor);
                outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + (Iter_MT_WI - 1) * 32, 11,
                    pptBarite_NoMassTransfer * 233390.0 * density_factor);
            } else {
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
            } else {
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
            } else {
                outputData.addData("Deposition Prediction", 51 - LoopTP, 8,
                    pptCalcite_MassTransfer * 100091.0 * density_factor * 0.35051);
                outputData.addData("Deposition Prediction", 51 - LoopTP, 12,
                    pptBarite_MassTransfer * 233390.0 * density_factor * 0.35051);
                outputData.addData("Deposition Prediction", 51 - LoopTP, 16,
                    pptFeCO3_MassTransfer * 233390.0 * density_factor * 0.35051);
            }
            
            // 沉积速率
            if (UseSI == 1) {
                outputData.addData("Deposition Prediction", 51 - LoopTP, 9, pptCalcite_MassTransfer_V); // cm/yr
                outputData.addData("Deposition Prediction", 51 - LoopTP, 13, pptBarite_MassTransfer_V); // cm/yr
                outputData.addData("Deposition Prediction", 51 - LoopTP, 17, pptFeCO3_MassTransfer_V); // cm/yr
                outputData.addData("Deposition Prediction", 51 - LoopTP, 18, m_CR_selected); // mm/yr
            } else {
                outputData.addData("Deposition Prediction", 51 - LoopTP, 9, pptCalcite_MassTransfer_V * 0.393701); // inch/yr
                outputData.addData("Deposition Prediction", 51 - LoopTP, 13, pptBarite_MassTransfer_V * 0.393701); // inch/yr
                outputData.addData("Deposition Prediction", 51 - LoopTP, 17, pptFeCO3_MassTransfer_V * 0.393701); // inch/yr
                outputData.addData("Deposition Prediction", 51 - LoopTP, 18, m_CR_selected * 39.37); // mil/yr
            }
            
            if (Use_Corr_in_Deposition == 0) {
                outputData.addData("Deposition Prediction", 51 - LoopTP, 18, std::string(""));
            }
            
            // WhatIf 情景
            if (Run_MassTransfer_WhatIf == 1) {
                int row_offset = (Iter_MT_WI - 1) * 32;
                
                if (usePTB == 0) {
                    outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 8,
                        pptCalcite_MassTransfer * 100091.0 * density_factor);
                    outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 12,
                        pptBarite_MassTransfer * 233390.0 * density_factor);
                } else {
                    outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 8,
                        pptCalcite_MassTransfer * 100091.0 * density_factor * 0.35051);
                    outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 12,
                        pptBarite_MassTransfer * 233390.0 * density_factor * 0.35051);
                }
                
                if (UseSI == 1) {
                    outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 9, 
                        pptCalcite_MassTransfer_V); // cm/yr
                    outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 13, 
                        pptBarite_MassTransfer_V); // cm/yr
                } else {
                    outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 9, 
                        pptCalcite_MassTransfer_V * 0.393701); // inch/yr
                    outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 13, 
                        pptBarite_MassTransfer_V * 0.393701); // inch/yr
                }
            }
        }
        
        // ==================== 饱和指数(SI)输出 ====================
        if (useSR == 0) {
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
            if (Run_MassTransfer_WhatIf == 1) {
                int row_offset = (Iter_MT_WI - 1) * 32;
                outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 6, SICal);
                outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 10, SIBar);
            }
            
            // RunStat 条件下的输出
            if (RunStat == 1) {
                if (RunStatSICalcSSP == 1) {
                    int row = calculateRowWithStat(6, StatTPnob, LoopTP);
                    outputData.addData("Halite analysis", row, 6, SIHal);
                } else if (RunStatSICalcSSP == 2) {
                    if (RunStatGoalSeek == 1) {
                        int row = calculateRowWithStat(11, StatTPnob, LoopTP);
                        outputData.addData(myname, row, 20, SIHal);
                    } else {
                        int row = calculateRowWithStat(6, StatTPnob, LoopTP);
                        outputData.addData("Halite analysis", row, 14, SIHal);
                    }
                } else if (RunStatSICalcSSP == 3) {
                    int row = calculateRowWithStat(6, StatTPnob, LoopTP);
                    outputData.addData("Output data sheet", row, 12, SIHal);
                }
            }
            
        } else {
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
            if (Run_MassTransfer_WhatIf == 1) {
                int row_offset = (Iter_MT_WI - 1) * 32;
                outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 6, pow(10.0, SICal));
                outputData.addData("Deposition Prediction_WhatIf", 51 - LoopTP + row_offset, 10, pow(10.0, SIBar));
            }
        }
        
        // ==================== 抑制剂风险评估输出 ====================
        if (H2Oevap != 1) {
            // Calcite 风险评估
            if (mc.size() > 0 && mc[0] > 0.0000001 && HCO3 > 0.0000001) {
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
                    if (InhNo1 == 5 || InhNo1 == 7 || InhNo1 == 8 || InhNo1 == 9 || 
                        InhNo1 == 10 || InhNo1 == 14 || InhNo1 == 18 || 
                        InhNo2 == 5 || InhNo2 == 7 || InhNo2 == 8 || InhNo2 == 9 || 
                        InhNo2 == 10 || InhNo2 == 14 || InhNo2 == 18) {
                        if (SIRisk.size() > InhNoCal && SIRisk[InhNoCal].size() > LoopTP) {
                            SIRisk[InhNoCal][LoopTP][2][2] = std::numeric_limits<double>::quiet_NaN();
                            SIRisk[InhNoCal][LoopTP][3][2] = std::numeric_limits<double>::quiet_NaN();
                        }
                    }
                }
                
                // 输出风险评估值
                if (useSR == 0) {
                    if (SIRisk.size() > InhNoCal && SIRisk[InhNoCal].size() > LoopTP) {
                        outputData.addData("Calcite", row, 19, SIRisk[InhNoCal][LoopTP][1][2]);
                        outputData.addData("Calcite", row, 20, SIRisk[InhNoCal][LoopTP][2][2]);
                        outputData.addData("Calcite", row, 21, SIRisk[InhNoCal][LoopTP][3][2]);
                    }
                } else {
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
        if (useSR == 0) {
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
        } else {
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
        outputData.addData("Input", 54, 9, ConcInhAn);
        outputData.addData("Input", 56, 9, ConcInhCel);
        outputData.addData("Input", 58, 9, ViscWatIst);
        outputData.addData("Input", 60, 9, CpPerMl);
        outputData.addData("Input", 62, 9, rhoTP);
    }
    
    // 第二部分：LoopTP = 10 或 LoopTP = 2（根据RunGoalSeek）时的输出
    if (RunGoalSeek != 1) {
        if (LoopTP == 10) {
            // 类似LoopTP=1的处理，但输出到第10列
            if (useSR == 0) {
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
            } else {
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
            outputData.addData("Input", 54, 10, ConcInhAn);
            outputData.addData("Input", 56, 10, ConcInhCel);
            outputData.addData("Input", 58, 10, ViscWatIst);
            outputData.addData("Input", 60, 10, CpPerMl);
            outputData.addData("Input", 62, 10, rhoTP);
        }
    } else {
        // RunGoalSeek = 1 的情况
        if (LoopTP == 2) {
            // 类似处理，输出到第10列
            if (useSR == 0) {
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
            } else {
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
            outputData.addData("Input", 54, 10, ConcInhAn);
            outputData.addData("Input", 56, 10, ConcInhCel);
            outputData.addData("Input", 58, 10, ViscWatIst);
            outputData.addData("Input", 60, 10, CpPerMl);
            outputData.addData("Input", 62, 10, rhoTP);
        }
    }
    
    // ==================== 表头和其他信息输出 ====================
    if (RunGoalSeek != 1) {
        // 设置表头
        outputData.addData("calcite", 3, 20, "SI(Inh=" + std::to_string(MaxInh / 2) + "mg/L)");
        outputData.addData("Calcite", 3, 21, "SI (Inh=" + std::to_string(MaxInh) + "mg/L)");
        outputData.addData("Barite", 3, 20, "SI (Inh=" + std::to_string(MaxInh / 2) + "mg/L)");
        outputData.addData("Barite", 3, 21, "SI (Inh=" + std::to_string(MaxInh) + "mg/L)");
        outputData.addData("Other SO4s", 3, 20, "SI (Inh=" + std::to_string(MaxInh / 2) + "mg/L)");
        outputData.addData("Other SO4s", 3, 21, "SI (Inh=" + std::to_string(MaxInh) + "mg/L)");
        outputData.addData("Other SO4s", 15, 20, "SI (Inh=" + std::to_string(MaxInh / 2) + "mg/L)");
        outputData.addData("Other SO4s", 15, 21, "SI (Inh=" + std::to_string(MaxInh) + "mg/L)");
        
        // 输出井名
        int nob1 = nob;
        if (nob > 2) nob1 = 2; // 只输出2个井名
        
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
        if (UseSI == 0) {
            outputData.addData("Calcite", 1, 15, mass_MeOH / 0.7914 / 159.0);
            outputData.addData("Calcite", 1, 12, mass_MEG / 1.1135 / 159.0);
            outputData.addData("Barite", 2, 13, mass_MeOH / 0.7914 / 159.0);
            outputData.addData("Barite", 3, 13, mass_MEG / 1.1135 / 159.0);
            outputData.addData("Halite", 2, 13, mass_MeOH / 0.7914 / 159.0);
            outputData.addData("Halite", 3, 13, mass_MEG / 1.1135 / 159.0);
            outputData.addData("Other SO4s", 2, 14, mass_MeOH / 0.7914 / 159.0);
            outputData.addData("Other SO4s", 3, 14, mass_MEG / 1.1135 / 159.0);
        } else {
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
        if (OutPutActCoefs == 1) {
            // Calcite 工作表中的活度系数
            if (UseSI == 0) {
                outputData.addData("Calcite", 14 - LoopTP, 24, TF);
                outputData.addData("Calcite", 14 - LoopTP, 25, Ppsia);
            } else {
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
    } else {
        // 如果不使用 excelWriter，就将数据保存到文件或数据库
        // 例如：保存为CSV文件或JSON格式
        std::cout << "数据已收集，但未写入Excel（未提供ExcelWriter实例）" << std::endl;
        std::cout << "共收集了 " << outputData.worksheets.size() << " 个工作表的数据" << std::endl;
    }
}

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
        
        // CSV格式：行,列,值
        fileStreams[filename] << row << "," << col << "," << value << std::endl;
    }
    
    void writeCell(const std::string& sheetName, int row, int col, const std::string& value) override {
        std::string filename = basePath + "/" + sheetName + ".csv";
        ensureFileOpen(filename);
        
        // 处理可能包含逗号的字符串
        std::string escapedValue = value;
        std::replace(escapedValue.begin(), escapedValue.end(), ',', ';');
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

// 方案2：使用第三方库（如OpenXLSX）的具体实现
#ifdef USE_OPENXLSX
#include "OpenXLSX.hpp"
class OpenXLSXWriter : public ExcelWriter {
private:
    OpenXLSX::XLDocument doc;
    
public:
    OpenXLSXWriter(const std::string& filename) {
        doc.create(filename);
    }
    
    void writeCell(const std::string& sheetName, int row, int col, double value) override {
        try {
            auto wks = doc.workbook().worksheet(sheetName);
            wks.cell(row, col).value() = value;
        } catch (...) {
            // 如果工作表不存在，创建它
            doc.workbook().addWorksheet(sheetName);
            auto wks = doc.workbook().worksheet(sheetName);
            wks.cell(row, col).value() = value;
        }
    }
    
    // ... 其他方法类似实现
};
#endif

// 方案3：使用Windows COM（仅Windows）
#ifdef _WIN32
#include <windows.h>
#include <oleauto.h>
#include <comdef.h>
class ExcelCOMWriter : public ExcelWriter {
private:
    IDispatch* excelApp = nullptr;
    IDispatch* workbooks = nullptr;
    IDispatch* workbook = nullptr;
    std::map<std::string, IDispatch*> worksheets;
    
public:
    ExcelCOMWriter() {
        // 初始化COM和Excel应用程序
        CoInitialize(nullptr);
        // ... 创建Excel实例的代码
    }
    
    ~ExcelCOMWriter() {
        // 清理COM对象
        CoUninitialize();
    }
    
    void writeCell(const std::string& sheetName, int row, int col, double value) override {
        // 使用COM接口写入单元格
        // ... 实现细节
    }
    
    // ... 其他方法
};
#endif

// ==================== 主程序示例 ====================
int main() {
    // 示例用法
    CSVExcelWriter csvWriter("./output");
    
    // 初始化全局变量（实际应用中从文件或数据库读取）
    RunStat = 1;
    LoopTP = 1;
    TF = 150.0;  // 华氏度
    Ppsia = 5000.0;
    
    // 调用转换后的函数
    LoopTPWrite(&csvWriter);
    
    std::cout << "数据输出完成" << std::endl;
    
    return 0;
}