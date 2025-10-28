# CalcIonicStrength

> code
此子程序设置单井 nTCO2、nTCH4、nTH2O 和 QC 计算的参数

此子程序在两个计算例程中通过 B2_ReadinAllData 例程为每个被检查的井调用。

```cpp
void CalcIonicStrength() {
    mtotal = 0.0;
    MoleCharge = 0.0;
    Ist = 0.0;
    SumOfCations = 0.0;

    // 处理阳离子
    for (int c = 0; c < NumCat; c++) {
        double charge = ChCat[c];
        double concentration = mc[c];

        MoleCharge += charge * concentration;
        SumOfCations += charge * concentration;
        Ist += charge * charge * concentration;  // 比pow()更高效
        mtotal += concentration;
    }

    SumOfAnions = 0.0;

    // 处理阴离子
    for (int a = 0; a < NumAn; a++) {
        double charge = ChAn[a];
        double abs_charge = fabs(charge);
        double concentration = ma[a];

        MoleCharge += abs_charge * concentration;
        SumOfAnions += charge * concentration;
        Ist += charge * charge * concentration;  // 比pow()更高效
        mtotal += concentration;
    }

    if (Ist <= 0.0) {
        Ist = 2.0 * 1e-7;  // 0.0000001
    }

    Ist *= 0.5;  // Ist = Ist / 2
    DpHj = 0.129 * sqrt(Ist);
}


```