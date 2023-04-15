#include<immintrin.h>
#include<iostream>
#include<windows.h>
#include<cstdlib>
using namespace std;
#define NUM 10

float A[NUM][NUM];

void gaussian_elimination(float A[][NUM], int n)
{

    for (int k = 0; k < NUM; k++) {
        __m128 v1 = _mm_set1_ps(A[k][k]); // ��A[k][k]װ��һ��SSE����v1��
        int j;
        for (j = k + 1; j + 4 <= NUM; j += 4) {
            __m128 va = _mm_loadu_ps(A[k] + j); // ��A[k][j]~A[k][j+7]װ��һ��SSE����va��
            va = _mm_div_ps(va, v1); // ִ��va/v1����
            _mm_storeu_ps(A[k] + j, va); // ��va�е�ֵд��A[k][j]~A[k][j+3]��
        }
        for (; j < NUM; j++)
            A[k][j] /= A[k][k];
        A[k][k] = 1.0; // A[k][k]=1
        for (int i = k + 1; i < NUM; i++) {
            __m128 vaik = _mm_set1_ps(A[i][k]); // ��A[i][k]װ��һ��SSE����vaik��
            for (j = k + 1; j + 8 <= NUM; j += 8) {
                __m128 vaij = _mm_loadu_ps(A[i] + j); // ��A[i][j]~A[i][j+7]װ��һ��SSE����vaij��
                __m128 vakj = _mm_loadu_ps(A[k] + j); // ��A[k][j]~A[k][j+7]װ��һ��SSE����vakj��
                __m128 vx = _mm_mul_ps(vaik, vakj); // ��vaik��vakj�е�ֵȡ�������
                vaij = _mm_sub_ps(vaij, vx); // ִ��vaij-vx����
                _mm_storeu_ps(A[i] + j, vaij); // ��vaij�е�ֵд��A[i][j]~A[i][j+7]��
            }
            for (; j < NUM; j++)
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            A[i][k] = 0; // A[i][k]=0
        }
    }

}







void m_reset()//������������
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
    int n = 1000;
    m_reset();//���ɲ�������
    // similar to CLOCKS_PER_SEC
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    // start time
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    while(n--)
       gaussian_elimination(A, NUM);
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    cout << "Col: " << (tail - head) * 1000.0 / freq << endl;
    return 0;
}