#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <immintrin.h>
#include <windows.h>
const int col = 3000, row = 60000;
int n, m1, m2;

using namespace std;
struct bitset {
	unsigned int bits[col];
	int lp = -1;
	void addNum(int num) {
		int n1 = num / 32, n2 = num % 32;
		bits[n1] += (0x80000000 >> n2);
		lp = max(lp, num);
		n = max(n, n1 + 1);
	}
};
bitset R[row], E[row];//R为消元子，E为被消元行
void print(bitset& b, ofstream& of) {
	if (b.lp == -1) {
		of << endl;
		return;
	}
	for (int i = n; i >= 0; --i) {
		for (unsigned int temp = 1, j = 31; temp != 0; temp <<= 1, --j)
			if (temp & b.bits[i])
				of << i * 32 + j << " ";
	}
	of << endl;
}
void Xor(bitset& b1, bitset& b2) //逐位进行异或操作
{
	int i = 0;
	for (; i + 4 <= n; i += 4)
	{
		__m128i v1 = _mm_loadu_si128((__m128i*) & b1.bits[i]);
		__m128i v2 = _mm_loadu_si128((__m128i*) & b2.bits[i]);
		__m128i result = _mm_xor_si128(v1, v2);
		_mm_storeu_si128((__m128i*) & b1.bits[i], result);
	}
	for (; i <= n; ++i)b1.bits[i] ^= b2.bits[i];
	int k = n, j = 31;
	while (b1.bits[k] == 0 && k >= 0)--k;
	if (k < 0) {
		b1.lp = -1;
		return;
	}
	unsigned temp = 1;
	while ((b1.bits[k] & temp) == 0 && j >= 0)temp <<= 1, --j;
	b1.lp = 32 * k + j;

}
void grobner_gaussian_elimination() {
	for (int i = 0; i < m2; ++i) {
		while (E[i].lp != -1) {
			if (R[E[i].lp].lp != -1) {
				Xor(E[i], R[E[i].lp]);
			}

			else {
				R[E[i].lp] = E[i];
				break;
			}
		}
	}
}
void readData() {
	fstream fs1("消元子3.txt");//鲲鹏
	string line;
	while (getline(fs1, line)) {
		stringstream ss(line);
		int index, num;
		ss >> index;
		R[index].addNum(index);
		while (ss >> num)
			R[index].addNum(num);
	}
	fs1.close();
	fs1.open("被消元行1.txt");//鲲鹏
	while (getline(fs1, line)) {
		stringstream ss(line);
		int num;
		while (ss >> num)
			E[m2].addNum(num);
		++m2;
	}
	fs1.close();
}
int main() {
	readData();

	long long head, tail, freq;
	QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
	QueryPerformanceCounter((LARGE_INTEGER*)&head);
	grobner_gaussian_elimination();
	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << (tail - head) * 1000.0 / freq  << endl;



	return 0;
}