# BubbleSort

> 代码部分

```cpp
#include <stdio.h>

/**
 * @brief 冒泡排序（返回排序后的索引数组）
 *
 * @param Ascend 排序方式：1 升序，其它值降序
 * @param ArrayIn 输入数组，长度为 n
 * @param index   索引数组，初始值应为 0,1,...,n-1
 * @param n       数组长度
 * @return        排序后的索引数组
 */
int* BubbleSort(int Ascend, double ArrayIn[], int index[], int n)
{
    double SrtTemp;
    int indSrtTemp;
    int i, j;

    if (Ascend == 1)
    {
        for (i = 0; i < n; i++)
        {
            for (j = i + 1; j < n; j++)
            {
                if (ArrayIn[i] > ArrayIn[j])
                {
                    SrtTemp = ArrayIn[j];
                    ArrayIn[j] = ArrayIn[i];
                    ArrayIn[i] = SrtTemp;

                    indSrtTemp = index[j];
                    index[j] = index[i];
                    index[i] = indSrtTemp;
                }
            }
        }
    }
    else
    {
        for (i = 0; i < n; i++)
        {
            for (j = i + 1; j < n; j++)
            {
                if (ArrayIn[i] < ArrayIn[j])
                {
                    SrtTemp = ArrayIn[j];
                    ArrayIn[j] = ArrayIn[i];
                    ArrayIn[i] = SrtTemp;

                    indSrtTemp = index[j];
                    index[j] = index[i];
                    index[i] = indSrtTemp;
                }
            }
        }
    }
    return index;
}
```

> 测试函数

```cpp
void testBubbleSort() {
    double arr[] = { 3.2, 1.5, 4.8, 2.1, 9.0 };
    int n = sizeof(arr) / sizeof(arr[0]);
    int index[5];
    int i;

    // 初始化索引数组
    for (i = 0; i < n; i++)
        index[i] = i;

    printf("原始数组: ");
    for (i = 0; i < n; i++)
        printf("%.2f ", arr[i]);
    printf("\n");

    // 升序排序
    BubbleSort(1, arr, index, n);

    printf("升序排序结果: ");
    for (i = 0; i < n; i++)
        printf("%.2f ", arr[i]);
    printf("\n");

    printf("升序索引顺序: ");
    for (i = 0; i < n; i++)
        printf("%d ", index[i]);
    printf("\n");

    // 降序排序
    double arr2[] = { 3.2, 1.5, 4.8, 2.1, 9.0 };
    int index2[5];
    for (i = 0; i < n; i++)
        index2[i] = i;

    BubbleSort(0, arr2, index2, n);

    printf("降序排序结果: ");
    for (i = 0; i < n; i++)
        printf("%.2f ", arr2[i]);
    printf("\n");

    printf("降序索引顺序: ");
    for (i = 0; i < n; i++)
        printf("%d ", index2[i]);
    printf("\n");

}
```

> 预期结果

```text
原始数组:       3.20 1.50 4.80 2.10 9.00
升序排序结果:   1.50 2.10 3.20 4.80 9.00
升序索引顺序:   1 3 0 2 4
降序排序结果:   9.00 4.80 3.20 2.10 1.50
降序索引顺序:   4 2 0 3 1
```
