#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

// gauss func
vector<double> gauss(vector<vector<double>>& matrix)
{
    int n = matrix.size();

    for (int i = 0; i < n; i++)
    {
        //maxElement
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (abs(matrix[k][i]) > abs(matrix[maxRow][i])) {
                maxRow = k;
            }
        }
        //pomenyat mestami
        swap(matrix[i], matrix[maxRow]);

        //treug vid
        for (int k = i + 1; k < n; k++) {
            double factor = matrix[k][i] / matrix[i][i];
            for (int j = i; j < n + 1; j++) {
                matrix[k][j] -= factor * matrix[i][j];
            }
        }
    }

    //gauss obratnii
    vector<double> solution(n, 0);
    for (int i = n - 1; i >= 0; i--) {
        solution[i] = matrix[i][n];
        for (int j = i + 1; j < n; j++) {
            solution[i] -= matrix[i][j] * solution[j];
        }
        solution[i] /= matrix[i][i];
    }

    return solution;
}

//funkciya dlya aproksimacii metodom naimensh kvadratov
void slau(const vector<double>& F, const vector<double>& v, double& c, double& e)
{
    int n = F.size();

    // postroenie SLU dlya MNK
    vector<vector<double>> matrix(2, vector<double>(3, 0));

    for (int i = 0; i < n; i++) {
        matrix[0][0] += log(F[i]);
        matrix[0][1] += log(F[i]) * log(F[i]);
        matrix[0][2] += log(F[i]) * v[i];
        matrix[1][0] += log(F[i]);
        matrix[1][1] += log(F[i]);
        matrix[1][2] += v[i];
    }

    //gaussmetod
    vector<double> coefficients = gauss(matrix);

    //izvlechenie coefficientov
    e = -coefficients[0];
    c = exp(coefficients[1]);

    cout << "Коэффициент c: " << c << endl;
    cout << "Коэффициент e: " << e << endl;
}

int main()
{
    setlocale(LC_ALL, "rus");

    vector<double> F = { 1.1, 1.4, 1.7, 2.1, 2.6, 4.7, 6.1, 7.0, 10.0, 12.8, 16.5, 20.8, 40.6 };
    vector<double> v = { 25.0, 22.7, 22.1, 19.8, 17.0, 12.3, 10.7, 10.0, 8.2, 6.7, 5.6, 5.0, 3.5 };

    double c, e;
    slau(F, v, c, e);//approksimaciya MNK

    return 0;
}