#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>

using namespace std;

void showMatrix(vector<vector<double>>matrix) {
	for (int i = 0; i < matrix.size(); i++) {
		for (int j = 0; j < matrix.size(); j++)
			cout << '\t' << matrix[i][j];
		cout << endl;
	}
	cout << endl;
}
vector<vector<double>> enterMatrix(int n) {
	vector<vector<double>> A(n, vector<double>(n));
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			cin >> A[i][j];
	return A;
}
double normMatrix(vector<vector<double>>matrix)
{
	int n = matrix.size();
	vector<double>sums(n);
	for (int i = 0; i < n; i++) {
		double sum = 0;
		for (int j = 0; j < n; j++)
			sum += fabs(matrix[i][j]);
		sums[i] = sum;
	}
	double max = *max_element(sums.begin(), sums.end());
	return max;
}
vector<double> gauss(vector<vector<double>> matrix_A, vector<double> vector_b)
{
	int p, n = matrix_A.size();
	double r, c, s;
	vector<double>x(n), b(n);
	vector<vector<double>>a(n, vector<double>(n));
	a = matrix_A;
	b = vector_b;
	for (int k = 0; k < n; k++) {
		p = k;
		for (int m = k + 1; m < n; m++)
			if (abs(a[p][k]) < abs(a[m][k])) //поиск максимального ведущего элемента
				p = m;
		for (int j = k; j < n; j++) {
			r = a[k][j];
			a[k][j] = a[p][j];   //перестановка строк
			a[p][j] = r;
		}
		r = b[k];
		b[k] = b[p];   //перестановка свободных членов
		b[p] = r;
		for (int m = k + 1; m < n; m++) {
			c = a[m][k] / a[k][k];
			b[m] = b[m] - c * b[k]; //приведение матрицы к верхнетреугольному виду
			for (int i = k; i < n; i++)
				a[m][i] = a[m][i] - c * a[k][i];
		}
	}
	x[n - 1] = b[n - 1] / a[n - 1][n - 1];
	for (int k = n - 1; k >= 0; k--) {
		s = 0;
		for (int i = k + 1; i < n; i++)	 //обратный ход метода Гаусса
			s = s + a[k][i] * x[i];
		x[k] = (b[k] - s) / a[k][k];
	}
	return x;
}
vector<vector<double>> reverseMatrix(vector<vector<double>> matrix)
{
	int n = matrix.size();
	vector<vector<double>> reverse(n, vector<double>(n));
	vector<double>b(n);
	vector<double>solve(n);
	for (int i = 0; i < n; i++)
	{
		b[i] = 1;
		solve = gauss(matrix, b);
		for (int j = 0; j < n; j++)
			reverse[j][i] = solve[j];
		b[i] = 0;
	}
	return reverse;
}
double сonditionNumberOfMatrix(vector<vector<double>> matrix) {
	return normMatrix(reverseMatrix(matrix)) * normMatrix(matrix);
}


int main()
{
	setlocale(LC_ALL, "Russian");
	int n = 0;
	double mu = 0;
	cout << "Введите n = ";
	cin >> n;
	vector<vector<double>> A(n,vector<double>(n));
	
	return 0;
}


