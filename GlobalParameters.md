# global parameters

```c
// constants.h
#ifndef CONSTANTS_H
#define CONSTANTS_H

// 数学常数
#define pi          3.14159265358979     // π (圆周率)

// 物理常数
#define NAv         6.0221367E+23        // 阿伏伽德罗常数 (mol⁻¹)
#define eElec       1.60217733E-19       // 电子电荷 (C)
#define eps0        8.854187818E-12      // 真空介电常数 (F/m)
#define kBoltz      1.380658E-23         // 玻尔兹曼常数 (J/K)

// 气体常数（不同单位）
#define RBar        0.083144             // Gas constant (L·bar/(K·mol))
#define R           83.144               // Gas constant (cm³·bar/(K·mol))
#define RAtm        0.082057             // Gas constant (L·atm/(K·mol))

#endif  // CONSTANTS_H

// Cation indexes: 阳离子索引（宏定义，编译时替换为对应整数）
#define iH      1   // H+ 的索引为 1
#define iNa     2   // Na+ 的索引为 2
#define iK      3   // K+ 的索引为 3
#define iMg     4   // Mg2+ 的索引为 4
#define iCa     5   // Ca2+ 的索引为 5
#define iSr     6   // Sr2+ 的索引为 6
#define iBa     7   // Ba2+ 的索引为 7
#define iFe     8   // Fe2+/3+ 的索引为 8（根据实际价态调整注释）
#define iZn     9   // Zn2+ 的索引为 9
#define iPb     10  // Pb2+ 的索引为 10
#define iNH4    11  // NH4+（铵根离子）的索引为 11
#define iRa     12  // Ra2+ 的索引为 12

// Anion indexes: 阴离子索引（宏定义，对应数组下标）
#define iOH         1   // OH- (氢氧根离子)
#define iCl         2   // Cl- (氯离子)
#define iAc         3   // Ac- (醋酸根离子，CH3COO-)
#define iHCO3       4   // HCO3- (碳酸氢根离子)
#define iCO3        5   // CO3^2- (碳酸根离子)
#define iSO4        6   // SO4^2- (硫酸根离子)
#define iHS         7   // HS- (硫氢根离子)
#define intF        8   // F- (氟离子，名称可能为"iF"的笔误，暂保留原拼写)
#define iBr         9   // Br- (溴离子)
#define iH2BO3      10  // H2BO3- (硼酸一氢根离子)
#define iH3SiO4     11  // H3SiO4- (原硅酸一氢根离子)
#define iH2SiO4     12  // H2SiO4^2- (原硅酸二氢根离子，或偏硅酸根，视具体体系而定)

// 特殊索引：离子总数/离子种类数（推测为 "ion species count" 或 "special ion index"）
#define iSion       13  // 可能表示阴离子总数，或特殊离子的索引（如"固体离子"、"表面离子"等）

// Neutral aquatic indexes
#define iCH4aq       1
#define iCO2aq       2
#define iH2Saq       3
#define iHAcaq       4
#define iH4SiO4aq    5
#define iNH3         6   // 注意：原数据中该名称未带aq后缀，保持原样
#define iH3BO3       7   // 同上
#define iFeSaq       8

// Oil phase indexes
#define iCH4o        1
#define iCO2o        2
#define iH2So        3

// Gas indexes
#define iCH4g        1
#define iCO2g        2
#define iH2Sg        3
#define iC2g         4
#define iC3g         5
#define iC4ig        6
#define iC4ng        7
#define iC5ig        8
#define iC5ng        9   // 原数据中此处序号为9，保持连续
#define iC6g         10
#define iC7_12g      11
#define iC13_25g     12
#define iC26_80g     13
#define iN2g         14  // 原数据中iN2g在iH2Og之后，调整为按序号排序
#define iH2Og        15

// 最后一组未命名类别（推测为矿物/盐类索引）
#define iBaSO4       1
#define iCaSO42H2O   2
#define iSrSO4       3
#define ihemiCaSO4   4
#define iCaSO4       5
#define iNaCl        6
#define iCaHCO32     7
#define iFeHCO32     8
#define iFeHS2       9

/* bDot物质种类索引定义（用于标识不同络合物种或离子的编号） */
#define iZnDot     1   // Zn²⁺离子（锌离子）
#define iPbDot     2   // Pb²⁺离子（铅离子）
#define iHSDot     3   // HS⁻离子（硫氢根离子）
#define iClDot     4   // Cl⁻离子（氯离子）
#define iZnCl      5   // ZnCl⁺络离子（一氯合锌离子）
#define iZnCl2     6   // ZnCl₂中性分子（二氯合锌）
#define iZnCl3     7   // ZnCl₃⁻络离子（三氯合锌离子）
#define iZnCl4     8   // ZnCl₄²⁻络离子（四氯合锌离子）
#define iZnHS2     9   // Zn(HS)₂中性分子（二硫氢合锌）
#define iZnHS3     10  // Zn(HS)₃⁻络离子（三硫氢合锌离子）
#define iPbCl      11  // PbCl⁺络离子（一氯合铅离子）
#define iPbCl2     12  // PbCl₂中性分子（二氯合铅）
#define iPbCl3     13  // PbCl₃⁻络离子（三氯合铅离子）
#define iPbCl4     14  // PbCl₄²⁻络离子（四氯合铅离子）
#define iPbHS2     15  // Pb(HS)₂中性分子（二硫氢合铅）
#define iPbHS3     16  // Pb(HS)₃⁻络离子（三硫氢合铅离子）
```
