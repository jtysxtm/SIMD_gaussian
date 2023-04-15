#include<immintrin.h>
#include<iostream>
#include<windows.h>
#include<stdlib.h>
#include<malloc.h>
using namespace std;
#define NUM 1000

float(*A)[NUM] = (float(*)[NUM])_aligned_malloc(NUM * NUM * sizeof(float),16);

void gaussian_elimination(float A[][NUM], int n) {
    // 使用 aligned_alloc 函数进行内存对齐操作，以 16 字节对齐

     for (int k = 0; k < NUM; k++) {
        __m128 v1 = _mm_set1_ps(A[k][k]); // 把A[k][k]装入一个SSE变量v1中
        int j;
        for (j = k + 1; j + 4 <= NUM; j += 4) {
            __m128 va = _mm_load_ps(A[k] + j); // 把A[k][j]~A[k][j+7]装入一个SSE变量va中
            va = _mm_div_ps(va, v1); // 执行va/v1操作
            _mm_store_ps(A[k] + j, va); // 把va中的值写回A[k][j]~A[k][j+3]中
        }
        for (; j < NUM; j++)
            A[k][j] /= A[k][k];
        A[k][k] = 1.0; // A[k][k]=1
        for (int i = k + 1; i < NUM; i++) {
            __m128 vaik = _mm_set1_ps(A[i][k]); // 把A[i][k]装入一个SSE变量vaik中
            for (j = k + 1; j + 8 <= NUM; j += 8) {
                __m128 vaij = _mm_load_ps(A[i] + j); // 把A[i][j]~A[i][j+7]装入一个SSE变量vaij中
                __m128 vakj = _mm_load_ps(A[k] + j); // 把A[k][j]~A[k][j+7]装入一个SSE变量vakj中
                __m128 vx = _mm_mul_ps(vaik, vakj); // 把vaik和vakj中的值取出来相乘
                vaij = _mm_sub_ps(vaij, vx); // 执行vaij-vx操作
                _mm_store_ps(A[i] + j, vaij); // 把vaij中的值写回A[i][j]~A[i][j+7]中
            }
            for (; j < NUM; j++)
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            A[i][k] = 0; // A[i][k]=0
        }
    }

}







void m_reset()//测试样例生成
{
    for (int i = 0; i < NUM; i++)
    {
        A[i][i] = 1.0;
        for (int j = i + 1; j < NUM; j++)
            A[i][j] = rand();
    }
    for (int k = 0; k < NUM; k++)
        for (int i = k + 1; i < NUM; i++)
            for (int j = 0; j < NUM; j++)
                A[i][j] += A[k][j];
}

int main()
{
    long long head, tail, freq; // timers
    int n = 10;
    m_reset();//生成测试用例
    // similar to CLOCKS_PER_SEC
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    // start time
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    while (n--)
        gaussian_elimination(A, NUM);
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    cout << "Col: " << (tail - head) * 1000.0 / freq << endl;
    _aligned_free(A);
    return 0;
}