#include <iostream>
#include <cstdlib>
#include <sys/time.h>
#include <arm_neon.h>
#include<time.h>

using namespace std;

#define NUM 50

float(*A)[NUM] = (float(*)[NUM])malloc(NUM * NUM * sizeof(float));

void gaussian_elimination(float A[][NUM], int n)
{
    for (int k = 0; k < NUM; k++)
    {
        float v1_0 = A[k][k];
        float32x4_t v1 = vdupq_n_f32(v1_0);
        int j;
        for (j = k + 1; j + 4 <= NUM; j += 4)
        {
            float32x4_t va = vld1q_f32(A[k] + j);
  float32x4_t vt=vrecpeq_f32(v1);
            va = vmulq_f32(va, v1);
            vst1q_f32(A[k] + j, va);
        }
        for (; j < NUM; j++)
            A[k][j] /= A[k][k];
        A[k][k] = 1.0;
        for (int i = k + 1; i < NUM; i++)
        {
            float vaik_0 = A[i][k];
            float32x4_t vaik = vdupq_n_f32(vaik_0);
            for (j = k + 1; j + 8 <= NUM; j += 8)
            {
                float32x4_t vaij_0 = vld1q_f32(A[i] + j);
                float32x4_t vakj_0 = vld1q_f32(A[k] + j);
                float32x4_t vx_0 = vmulq_f32(vaik, vakj_0);
                float32x4_t vaij = vsubq_f32(vaij_0, vx_0);

                float32x4_t vaij_1 = vld1q_f32(A[i] + j + 4);
                float32x4_t vakj_1 = vld1q_f32(A[k] + j + 4);
                float32x4_t vx_1 = vmulq_f32(vaik, vakj_1);
                float32x4_t vaij_2 = vsubq_f32(vaij_1, vx_1);

                vst1q_f32(A[i] + j, vaij);
                vst1q_f32(A[i] + j + 4, vaij_2);
            }
            for (; j < NUM; j++)
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            A[i][k] = 0;
        }
    }
}

void m_reset()
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
    m_reset();
    struct timespec sts,ets;
    timespec_get(&sts, TIME_UTC);

     for (int i = 0; i < 50; i++)
        gaussian_elimination(A, NUM);
    
    timespec_get(&ets, TIME_UTC);
    time_t dsec=ets.tv_sec-sts.tv_sec;
    long dnsec=ets.tv_nsec-sts.tv_nsec;
    if(dnsec<0){
        dsec--;
        dnsec+=1000000000ll;
     }

    cout << "Col: " << dsec<<"."<<dnsec << endl;
    free(A);
    return 0;
}

