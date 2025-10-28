# fTotalCO2H2Smoles

> code

```cpp
void* fTotalCO2H2Smoles(double yCO2, double* RatioOilBPoints, double Znew) {
    // nTCO2 等于气体中的 CO2 加上油中的 CO2 加上水中的 CO2 以及水中的 HCO3 和 CO3。829= (1/(R*14.696) 和 R=atm m^3/(K*mol)
    // 829 = (1/(R*14.696) and R = atm m^3/(K*mol)
    nTCO2 = Patm * 14.696 * yCO2 * (829.0 * VgTP / (Znew * TK) +
        gGas[iCO2g] * (KgwCO2 * mass_w / (gNeut[iCO2aq] * gNNeut[iCO2aq]) +
            RatioOilBPoints * KgoCO2 * Mass_o / gL[iCO2o])) +
        (HCO3 + CO3) * mass_w;

    nTCO2EOS = nTCO2;

    nTCH4 = Patm * 14.696 * yCH4 * (829.0 * VgTP / (Znew * TK) +
        gGas[iCH4g] * (KgwCH4 * mass_w / gNeut[iCH4aq] +
            RatioOilBPoints * KgoCH4 * Mass_o / gL[iCH4g]));

    nTH2S = Patm * 14.696 * yH2S * (829.0 * VgTP / (Znew * TK) +
        gGas[iH2Sg] * (KgwH2S * mass_w / (gNeut[iH2Saq] * gNNeut[iH2Saq]) +
            RatioOilBPoints * KgoH2S * Mass_o / gL[iH2So])) +
        HS * (1 + KstFeSaq * mc[iFe] * HS * gAn[iHS] * gNAn[iHS] * gCat[iFe] * gNCat[iFe] /
            (gNeut[iFeSaq] * gNNeut[iFeSaq] * aH)) * mass_w;

    nTH2sEOS = nTH2S;

    return NULL;  // 对应VB的Null
}

```