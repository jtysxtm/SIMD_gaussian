#include<iostream>
#include<windows.h>
#include<cstdlib>
using namespace std;

#define N 1000

float A[N][N];


void gaussian_elimination(float* A, int n) {
	for (int k = 0; k < n; k++) {
		for (int j = k + 1; j < n; j++) {
			A[k * n + j] /= A[k * n + k];
		}
		A[k * n + k] = 1.0f;
		for (int i = k + 1; i < n; i++) {
			for (int j = k + 1; j < n; j++) {
				A[i * n + j] -= A[i * n + k] * A[k * n + j];
			}
			A[i * n + k] = 0.0f;
		}
	}
}

void m_reset()//测试样例生成
{
	for (int i = 0; i < N; i++)
	{
		A[i][i] = 1.0;
		for (int j = i + 1; j < N; j++)
			A[i][j] = rand();
	}
	for (int k = 0; k < N; k++)
		for (int i = k + 1; i < N; i++)
			for (int j = 0; j < N; j++)
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
		gaussian_elimination(*A, N);
	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << "Col: " << (tail - head) * 1000.0 / freq << endl;
	return 0;
}