# fmn
> code

```cpp
void fmn() {
    // 初始化中性物种浓度
    mn[0] = 0;           // 原mn(1)对应mn[0]
    mn[1] = CO2aq;       // 原mn(2)对应mn[1]
    mn[2] = H2Saq;       // 原mn(3)对应mn[2]
    mn[3] = HAcaq;       // 原mn(4)对应mn[3]
    mn[4] = H4SiO4;      // 原mn(5)对应mn[4]
    mn[5] = NH3;         // 原mn(6)对应mn[5]
    mn[6] = 0;           // 原mn(7)对应mn[6]

    // 计算FeSaq浓度（基于平衡常数）
    mn[7] = KstFeSaq * mc[iFe] * HS * gAn[iHS] * gNAn[iHS] * gCat[iFe] * gNCat[iFe] /
        (gNeut[iFeSaq] * gNNeut[iFeSaq] * aH);  // 原mn(8)对应mn[7]

    // 设置阴离子浓度
    ma[0] = OH;          // 原ma(1)对应ma[0]
    ma[2] = AC;          // 原ma(3)对应ma[2]
    ma[3] = HCO3;        // 原ma(4)对应ma[3]
    ma[4] = CO3;         // 原ma(5)对应ma[4]
    ma[6] = HS;          // 原ma(7)对应ma[6]
    ma[9] = H2BO3;       // 原ma(10)对应ma[9]
    ma[10] = H3SiO4;     // 原ma(11)对应ma[10]
    ma[11] = H2SiO4;     // 原ma(12)对应ma[11]

    // 设置阳离子浓度
    mc[0] = H;           // 原mc(1)对应mc[0]
}

```