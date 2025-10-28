# fRatioOilBPoints

> code

```cpp
double fRatioOilBPoints() {
    double base = SGG / 0.6;
    double rankine = TF + 460.0; // 转换为兰金温度
    double api_diff = API - 30.0;

    if (API >= 30.0) {
        return pow(base * exp(23.93 * api_diff / rankine), 1.0 / 1.187);
    }
    else {
        return pow(base * exp(25.724 * api_diff / rankine), 1.0 / 1.0937);
    }
}

```