#include<iostream>
#include<fstream>
#include<sstream>
#include<windows.h>
#include<stdlib.h>
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
bitset R[row], E[row];//RΪ��Ԫ�ӣ�EΪ����Ԫ��
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

void Xor(bitset& b1, bitset& b2) //��λ����������
{
	for (int i = n; i >= 0; --i) {
		b1.bits[i] ^= b2.bits[i];
	}
	int i = n, j = 31;
	while (b1.bits[i] == 0 && i >= 0)--i;
	if (i < 0) {
		b1.lp = -1;
		return;
	}
	unsigned temp = 1;
	while ((b1.bits[i] & temp) == 0)temp <<= 1, --j;
	b1.lp = 32 * i + j;
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
	fstream fs1("��Ԫ��9.txt");//��������1
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
	fs1.open("����Ԫ��9.txt");//��������2
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
	cout << (tail - head) * 1000.0 / freq << endl;

	//����
	//ofstream of("output.txt");
	//for (int i = 0; i < m2; ++i)
	//	print(E[i], of);
	//of.close();
	//ifstream if1("output.txt");
	//ifstream if2("3.txt");//����
	//int a = 0, b = 0, sum = 0;
	//while (if1 >> a || if2 >> b) {
	//	sum += (a - b) * (a - b);
	//}
	//cout << sum;
	//if1.close();
	//if2.close();

	//ofstream of_result("../result.txt", ios::app);

	//of_result << "special gauss given m2 =" << m2 << "time:" << (tail - head) * 1000000.0 / freq << endl << "s";
	//of_result.close();
	return 0;
}