#include <stdio.h>
#include <math.h>
// 只是象征性的计算一下，很多函数都不需要测试，因为公式固定，计算值也固定
double fPP(double TK, double p0, double p1, double p2,
           double p3, double p4, double p5)
{
    const double Tr = 298.15;
    double value = p0
                 + p1 * (TK - Tr)
                 + p2 * (TK * TK - Tr * Tr)
                 + p3 * (1.0 / TK - 1.0 / Tr)
                 + p4 * log(TK / Tr)
                 + p5 * (1.0 / (TK * TK) - 1.0 / (Tr * Tr));
    return value;
}

void test_fPP() {
    double TK = 300.0;
    double p0 = 1.0, p1 = 0.01, p2 = 0.001;
    double p3 = -10.0, p4 = 2.0, p5 = 0.5;

    double result = fPP(TK, p0, p1, p2, p3, p4, p5);

    printf("Test fPP with TK=%.2f\n", TK);
    printf("Parameters: p0=%.2f, p1=%.2f, p2=%.2f, p3=%.2f, p4=%.2f, p5=%.2f\n",
           p0, p1, p2, p3, p4, p5);
    printf("Result = %.6f\n", result);

    double Tr = 298.15;
    double expected = p0
                    + p1 * (TK - Tr)
                    + p2 * (TK * TK - Tr * Tr)
                    + p3 * (1.0 / TK - 1.0 / Tr)
                    + p4 * log(TK / Tr)
                    + p5 * (1.0 / (TK * TK) - 1.0 / (Tr * Tr));

    printf("Expected = %.6f\n", expected);
    printf("Difference = %.6e\n\n", fabs(result - expected));
}

double fHolmes(double TK, double PBar,
               double q1, double q2, double q3, double q4, double q5,
               double q6, double q7, double q8, double q9, double q10,
               double q11, double q12, double q13, double q14, double q15,
               double q16, double q17) 
{
    double term1 = q1 
                 + 0.5 * q2 * TK
                 + (1.0/6.0) * q3 * (TK * TK)
                 + (1.0/12.0) * q4 * (TK * TK * TK)
                 + (1.0/6.0) * q5 * (TK * TK) * (log(TK) - 5.0/6.0);

    double term2 = q6 * (TK / 2.0
                 + (3.0 * 227.0 * 227.0) / (2.0 * TK)
                 + 227.0 * (TK - 227.0) / TK * log(TK - 227.0));

    double term3 = q7 * ((2.0 * (647.0 - TK) / TK + 1.0) * log(647.0 - TK));

    double term4 = PBar * (q8 
                 + q9 / TK 
                 + q10 * TK 
                 + q11 * (TK * TK) 
                 + q12 / (TK - 227.0) 
                 + q13 / (647.0 - TK));

    double term5 = (PBar * PBar) * (q14 
                 + q15 / TK 
                 + q16 * TK 
                 + q17 * (TK * TK));

    return term1 + term2 + term3 + term4 + term5;
}

void test_fHolmes() {
    double TK = 300.0;
    double PBar = 1.0;

    double q1=1, q2=0.1, q3=0.01, q4=0.001, q5=0.0005;
    double q6=0.05, q7=0.02, q8=0.3, q9=-0.2, q10=0.001;
    double q11=1e-5, q12=0.5, q13=-0.1, q14=0.01, q15=0.02;
    double q16=-0.0001, q17=5e-6;

    double result = fHolmes(TK, PBar,
                            q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,
                            q11,q12,q13,q14,q15,q16,q17);

    printf("Test fHolmes with TK=%.2f, PBar=%.2f\n", TK, PBar);
    printf("Result = %.6f\n", result);

    double expected = fHolmes(TK, PBar,
                              q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,
                              q11,q12,q13,q14,q15,q16,q17);
    printf("Expected = %.6f\n", expected);
    printf("Difference = %.6e\n\n", fabs(result - expected));
}

double fDuan1(double TK, double PBar,
              double q1, double q2, double q3, double q4,
              double q5, double q6, double q7, double q8)
{
    double term1 = q1 
                 + q2 * TK
                 + q3 / TK
                 + q4 / (TK - 210.0)
                 + q5 / (647.0 - TK);

    double term2 = q6 * pow(TK - 443.0, 3.0) / 3.0;

    double term3 = q7 * (PBar - 1.0)
                 + q8 * pow(PBar - 1.0, 2.0) / 2.0;

    return term1 + term2 + term3;
}

void test_fDuan1() {
    double TK   = 350.0;
    double PBar = 5.0;

    double q1=1.0, q2=0.01, q3=50.0, q4=20.0;
    double q5=10.0, q6=1e-4, q7=0.5, q8=0.05;

    double result = fDuan1(TK, PBar, q1,q2,q3,q4,q5,q6,q7,q8);

    printf("Test fDuan1 with TK=%.2f, PBar=%.2f\n", TK, PBar);
    printf("Parameters: q1=%.2f q2=%.2f q3=%.2f q4=%.2f q5=%.2f q6=%.6f q7=%.2f q8=%.2f\n",
           q1,q2,q3,q4,q5,q6,q7,q8);
    printf("Result   = %.6f\n", result);

    double term1 = q1 
                 + q2 * TK
                 + q3 / TK
                 + q4 / (TK - 210.0)
                 + q5 / (647.0 - TK);

    double term2 = q6 * pow(TK - 443.0, 3.0) / 3.0;

    double term3 = q7 * (PBar - 1.0)
                 + q8 * pow(PBar - 1.0, 2.0) / 2.0;

    double expected = term1 + term2 + term3;

    printf("Expected = %.6f\n", expected);
    printf("Difference = %.6e\n\n", fabs(result - expected));
}

double fDuan2(double TK, double PBar,
              double q1, double q2, double q3, double q4, double q5,
              double q6, double q7, double q8, double q9,
              double q10, double q11)
{
    double denom = 630.0 - TK;

    double term1 = q1 
                 + q2 * TK 
                 + q3 / TK 
                 + q4 * (TK * TK) 
                 + q5 / denom;

    double term2 = q6 * PBar
                 + q7 * PBar * log(TK)
                 + q8 * PBar / TK
                 + q9 / denom;

    double term3 = q10 * (PBar * PBar) / (denom * denom)
                 + q11 * TK * log(PBar);

    return term1 + term2 + term3;
}

void test_fDuan2() {
    double TK   = 400.0;
    double PBar = 10.0;

    double q1=1.0, q2=0.01, q3=5.0, q4=1e-5, q5=20.0;
    double q6=0.05, q7=0.001, q8=-0.02, q9=10.0;
    double q10=0.005, q11=0.0001;

    double result = fDuan2(TK, PBar,
                           q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11);

    printf("Test fDuan2 with TK=%.2f, PBar=%.2f\n", TK, PBar);
    printf("Parameters: q1=%.2f q2=%.2f q3=%.2f q4=%.5f q5=%.2f ... q11=%.5f\n",
           q1,q2,q3,q4,q5,q11);
    printf("Result   = %.6f\n", result);

    double expected = fDuan2(TK, PBar,
                             q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11);
    printf("Expected = %.6f\n", expected);
    printf("Difference = %.6e\n\n", fabs(result - expected));
}

double fPZ_DL_6(double TK, double Patm,
                double z1, double z2, double z3,
                double z4, double z5, double z6)
{
    double Pbar = Patm * 1.013254;
    double value = z1
                 + z2 * TK
                 + z3 / (TK - 210.0)
                 + z4 / (Pbar + 100.0)
                 + z5 * TK * log(Pbar)
                 + z6 * TK * Pbar * log(Pbar);
    return value;
}

void test_fPZ_DL_6() {
    double TK   = 350.0;
    double Patm = 2.0;

    double z1=1.0, z2=0.01, z3=50.0, z4=20.0, z5=1e-4, z6=1e-6;

    double result = fPZ_DL_6(TK, Patm, z1,z2,z3,z4,z5,z6);

    printf("Test fPZ_DL_6 with TK=%.2f, Patm=%.2f\n", TK, Patm);
    printf("Parameters: z1=%.2f z2=%.2f z3=%.2f z4=%.2f z5=%.6f z6=%.6f\n",
           z1,z2,z3,z4,z5,z6);
    printf("Result   = %.6f\n", result);

    double Pbar = Patm * 1.013254;
    double expected = z1
                    + z2 * TK
                    + z3 / (TK - 210.0)
                    + z4 / (Pbar + 100.0)
                    + z5 * TK * log(Pbar)
                    + z6 * TK * Pbar * log(Pbar);

    printf("Expected = %.6f\n", expected);
    printf("Difference = %.6e\n\n", fabs(result - expected));
}

double fV0(double TK, double PBar,
           double q1, double q2, double q3, double q4, double q5, double q6,
           double q7, double q8, double q9, double q10, double q11, double q12,
           double q13, double q14, double q15, double q16, double q17, double q18)
{
    double denom1 = TK - 227.0;
    double denom2 = 647.0 - TK;

    double term1 = q1
                 + q2 / TK
                 + q3 * TK
                 + q4 * (TK * TK)
                 + q5 / denom1
                 + q6 / denom2;

    double term2 = PBar * (q7
                 + q8 / TK
                 + q9 * TK
                 + q10 * (TK * TK)
                 + q11 / denom1
                 + q12 / denom2);

    double term3 = (PBar * PBar) * (q13
                 + q14 / TK
                 + q15 * TK
                 + q16 * (TK * TK)
                 + q17 / denom1
                 + q18 / denom2);

    return term1 + term2 + term3;
}

void test_fV0() {
    double TK   = 350.0;
    double PBar = 5.0;

    double q1=1.0, q2=10.0, q3=0.01, q4=1e-5, q5=20.0, q6=30.0;
    double q7=0.5, q8=5.0, q9=0.001, q10=1e-6, q11=2.0, q12=3.0;
    double q13=0.05, q14=0.5, q15=1e-4, q16=1e-7, q17=0.2, q18=0.3;

    double result = fV0(TK, PBar,
                        q1,q2,q3,q4,q5,q6,
                        q7,q8,q9,q10,q11,q12,
                        q13,q14,q15,q16,q17,q18);

    printf("Test fV0 with TK=%.2f, PBar=%.2f\n", TK, PBar);
    printf("Result   = %.6f\n", result);

    double denom1 = TK - 227.0;
    double denom2 = 647.0 - TK;

    double term1 = q1
                 + q2 / TK
                 + q3 * TK
                 + q4 * (TK * TK)
                 + q5 / denom1
                 + q6 / denom2;

    double term2 = PBar * (q7
                 + q8 / TK
                 + q9 * TK
                 + q10 * (TK * TK)
                 + q11 / denom1
                 + q12 / denom2);

    double term3 = (PBar * PBar) * (q13
                 + q14 / TK
                 + q15 * TK
                 + q16 * (TK * TK)
                 + q17 / denom1
                 + q18 / denom2);

    double expected = term1 + term2 + term3;
    printf("Expected = %.6f\n", expected);
    printf("Difference = %.6e\n\n", fabs(result - expected));
}

// 定义外部数组
double nreg1[34];
int ireg1[34];
int jreg1[34];

double fH2ODensity(double TK, double PBar) {
    const double rgas_water = 0.461526;
    double tau, piDen, gammapireg1;
    double volreg1;
    int igammapireg;

    tau = 1386.0 / TK;
    piDen = 0.1 * PBar / 16.53;

    gammapireg1 = 0.0;
    for (igammapireg = 0; igammapireg < 34; igammapireg++) {
        gammapireg1 -= nreg1[igammapireg] * ireg1[igammapireg] *
                       pow(7.1 - piDen, ireg1[igammapireg] - 1) *
                       pow(tau - 1.222, jreg1[igammapireg]);
    }

    volreg1 = rgas_water * TK * piDen * gammapireg1 / (PBar * 100000.0);

    if (volreg1 != 0.0) {
        return 1.0 / volreg1;
    } else {
        return -1.0;
    }
}

void test_fH2ODensity() {
    for (int i = 0; i < 34; i++) {
        nreg1[i] = 0.1 * (i + 1);
        ireg1[i] = (i % 3) + 1;
        jreg1[i] = (i % 4) + 1;
    }

    double TK = 300.0;
    double PBar = 1.0;

    double density = fH2ODensity(TK, PBar);

    printf("Test fH2ODensity with TK=%.2f K, PBar=%.2f bar\n", TK, PBar);
    printf("Density = %.6f kg/m^3\n", density);
}

double fPZ6(double TK, double z1, double z2, double z3, double z4, double z5, double z6) {
    double result;

    result = z1
           + z2 * (100.0 / TK)
           + z3 * (0.01 * TK)
           + z4 * (0.0001 * TK * TK)
           + z5 * (10.0 / (TK - 227.0))
           + z6 * log(TK);

    return result;
}

void test_fPZ6() {
    double TK = 300.0;

    double z1=1.0, z2=2.0, z3=0.5, z4=0.01, z5=20.0, z6=0.1;

    double result = fPZ6(TK, z1, z2, z3, z4, z5, z6);

    printf("Test fPZ6 with TK=%.2f\n", TK);
    printf("Parameters: z1=%.2f z2=%.2f z3=%.2f z4=%.2f z5=%.2f z6=%.2f\n",
           z1,z2,z3,z4,z5,z6);
    printf("Result   = %.6f\n", result);

    double expected = z1
                    + z2 * (100.0 / TK)
                    + z3 * (0.01 * TK)
                    + z4 * (0.0001 * TK * TK)
                    + z5 * (10.0 / (TK - 227.0))
                    + z6 * log(TK);

    printf("Expected = %.6f\n", expected);
    printf("Difference = %.6e\n\n", fabs(result - expected));
}

void CubicRoots(double coef2, double coef3, double coef4,
                double* root1, double* root2, double* root3) {
    const double pi = 3.14159265358979323846;
    *root1 = 0;
    *root2 = 0;
    *root3 = 0;

    double QCubic = (coef2 * coef2 - 3.0 * coef3) / 9.0;
    double RCubic = (2.0 * coef2 * coef2 * coef2 - 9.0 * coef2 * coef3 + 27.0 * coef4) / 54.0;

    if (RCubic * RCubic - QCubic * QCubic * QCubic > 0) {
        double A = pow(fabs(RCubic) + sqrt(RCubic * RCubic - QCubic * QCubic * QCubic), 1.0 / 3.0);
        if (RCubic < 0) A = -A;
        *root1 = - (A + QCubic / A) - coef2 / 3.0;
    } else {
        double Xcubic = RCubic / sqrt(QCubic * QCubic * QCubic);
        double Theta = atan(-Xcubic / sqrt(1 - Xcubic * Xcubic)) + pi / 2.0;
        double sqrtQ = 2.0 * sqrt(QCubic);
        *root1 = - sqrtQ * cos(Theta / 3.0) - coef2 / 3.0;
        *root2 = - sqrtQ * cos((Theta + 2.0 * pi) / 3.0) - coef2 / 3.0;
        *root3 = - sqrtQ * cos((Theta + 4.0 * pi) / 3.0) - coef2 / 3.0;
    }
}

void test_CubicRoots() {
    // 测试一组三次方程: x^3 - 6x^2 + 11x - 6 = 0, 根为 1, 2, 3
    double coef2 = -6.0;
    double coef3 = 11.0;
    double coef4 = -6.0;

    double r1, r2, r3;
    CubicRoots(coef2, coef3, coef4, &r1, &r2, &r3);

    printf("Test CubicRoots for x^3 - 6x^2 + 11x - 6 = 0\n");
    printf("Roots: r1=%.6f, r2=%.6f, r3=%.6f\n", r1, r2, r3);

    double f1 = r1*r1*r1 + coef2*r1*r1 + coef3*r1 + coef4;
    double f2 = r2*r2*r2 + coef2*r2*r2 + coef3*r2 + coef4;
    double f3 = r3*r3*r3 + coef2*r3*r3 + coef3*r3 + coef4;

    printf("Check: f(r1)=%.6e, f(r2)=%.6e, f(r3)=%.6e\n\n", f1, f2, f3);
}

void CalcIonicStrength(int NumCat, int NumAn,
                       const double mc[], const double ma[],
                       const int ChCat[], const int ChAn[],
                       double *mtotal, double *MoleCharge,
                       double *Ist, double *SumOfCations,
                       double *SumOfAnions, double *DpHj) 
{
    *mtotal = 0.0;
    *MoleCharge = 0.0;
    *Ist = 0.0;
    *SumOfCations = 0.0;

    for (int c = 0; c < NumCat; ++c) {
        *MoleCharge += ChCat[c] * mc[c];
        *SumOfCations += ChCat[c] * mc[c];
        *Ist += (double)(ChCat[c] * ChCat[c]) * mc[c];
        *mtotal += mc[c];
    }

    *SumOfAnions = 0.0;

    for (int a = 0; a < NumAn; ++a) {
        *MoleCharge += fabs((double)ChAn[a]) * ma[a];
        *SumOfAnions += ChAn[a] * ma[a];
        *Ist += (double)(ChAn[a] * ChAn[a]) * ma[a];
        *mtotal += ma[a];
    }

    if (*Ist <= 0.0) {
        *Ist = 2.0 * 1.0e-7;
    }

    *Ist /= 2.0;

    *DpHj = 0.129 * sqrt(*Ist);
}

void test_CalcIonicStrength() {
    double mc[2] = {0.1, 0.05};
    int ChCat[2] = {1, 2};

    double ma[2] = {0.1, 0.025};
    int ChAn[2] = {-1, -2};

    double mtotal, MoleCharge, Ist, SumC, SumA, DpHj;

    CalcIonicStrength(2, 2, mc, ma, ChCat, ChAn,
                      &mtotal, &MoleCharge, &Ist, &SumC, &SumA, &DpHj);

    printf("Test CalcIonicStrength:\n");
    printf("Total concentration (mtotal) = %.6f\n", mtotal);
    printf("MoleCharge                  = %.6f\n", MoleCharge);
    printf("Ionic strength (Ist)        = %.6f\n", Ist);
    printf("SumOfCations                = %.6f\n", SumC);
    printf("SumOfAnions                 = %.6f\n", SumA);
    printf("DpHj                        = %.6f\n", DpHj);
}

int main() {
    test_fPP();
    test_fHolmes();
    test_fDuan1();
    test_fDuan2();
    test_fPZ_DL_6();
    test_fV0();
    test_fH2ODensity();
    test_fPZ6();
    test_CubicRoots();
    test_CalcIonicStrength();
    return 0;
}