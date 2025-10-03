#include <vector>
#include <string>
#include <iostream>
#include "linearalgebra.h"


#define EPS 1e-10
#define INF 10

using namespace std;

namespace la{

void printMatrix(const vector<vector<double>>& matrix) //Auxiliary print function
{
    for (int i = 0; i < matrix.size(); i++)
    {
        for (int j = 0; j < matrix[i].size(); j++)
        {
            cout << matrix[i][j] << ' ';
        }

        cout << "\n";
        
    }
}

int gauss (vector < vector<double> > a, vector<double> & ans) { //Gauss elimination... kinda slow
    int n = (int) a.size();
    int m = (int) a[0].size() - 1;

    vector<int> where (m, -1);
    for (int col=0, row=0; col<m && row<n; ++col) {
        int sel = row;
        for (int i=row; i<n; ++i)
            if (abs (a[i][col]) > abs (a[sel][col]))
                sel = i;
        if (abs (a[sel][col]) < EPS)
            continue;
        for (int i=col; i<=m; ++i)
            swap (a[sel][i], a[row][i]);
        where[col] = row;

        for (int i=0; i<n; ++i)
            if (i != row) {
                double c = a[i][col] / a[row][col];
                for (int j=col; j<=m; ++j)
                    a[i][j] -= a[row][j] * c;
            }
        ++row;
    }

    ans.assign (m, 0);
    for (int i=0; i<m; ++i)
        if (where[i] != -1)
            ans[i] = a[where[i]][m] / a[where[i]][i];
    for (int i=0; i<n; ++i) {
        double sum = 0;
        for (int j=0; j<m; ++j)
            sum += ans[j] * a[i][j];
        if (abs (sum - a[i][m]) > EPS)
            return 0;
    }

    for (int i=0; i<m; ++i)
        if (where[i] == -1)
            return INF;
    return 1;
}

vector<vector<double>> transpose(const vector<vector<double>>& matrix)
{
    int rows = matrix.size();
    int cols = matrix[0].size();

    vector<vector<double>> tMat(cols, vector<double>(rows));

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            tMat[j][i] = matrix[i][j];
        }
    }

    return tMat;
}
}