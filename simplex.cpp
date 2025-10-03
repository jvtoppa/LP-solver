#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <utility>
#include <chrono>
#include <fstream>
#include <string>
#include "linearalgebra.h"

#define EPS 1e-10
#define INF 1e9

using namespace std;

void phase1(vector<double>& costVector, 
    vector<vector<double>>& restrictionsM, 
    vector<int>& basicColumns,
    std::vector<int>& nonBasicColumns)
{
    int n = restrictionsM.size();    // Rows
    int m = restrictionsM[0].size(); // Columns

    
    vector<bool> rowsUsed(n, false); 
    vector<bool> zeroRows(n, false);
    basicColumns.clear();
    nonBasicColumns.clear();
    cout << "\n";


    for (int row = 0; row < n; row++)
    {
        bool isZeroRow = true;
        for (int col = 0; col < m; col++)
        {
            if(abs(restrictionsM[row][col]) > EPS)
            {
                isZeroRow = false;
            }
        }
        zeroRows[row] = isZeroRow;
    }
    



    for (int col = 0; col < m; col++)
    {
    bool basic = true;
    int oneRow = -1;
    for (int row = 0; row < n; row++)
    {
        if (abs(restrictionsM[row][col] - 1.0) < EPS)
        {
            if (oneRow != -1)
            {
                basic = false;
                break;
            }
            oneRow = row;
            
        }
        else if (abs(restrictionsM[row][col]) > EPS)
        {
            basic = false;
            break;
        }
    }

    if (basic && oneRow != -1 && !rowsUsed[oneRow])
        {
            basicColumns.push_back(col);
            rowsUsed[oneRow] = true;
        }
    else
    {
        nonBasicColumns.push_back(col);
    }
    }

    double M = INF;

    for (int i = 0; i < n; i++) {
        if (!rowsUsed[i]) {
            int newColIndex = restrictionsM[0].size();
    
            for (int row = 0; row < n; row++) {
                if (row == i) {
                    restrictionsM[row].push_back(1.0);
                } else {
                    restrictionsM[row].push_back(0.0);
                }
            }

            costVector.push_back(M);
    
            basicColumns.push_back(newColIndex);
            rowsUsed[i] = true;
    
        }
    }
    
    return;
}

struct relativeCost
{
    int lowestNonBColumn{};
    double lowestRC{};
};

// -x- Auxiliary phase 2 functions

relativeCost relativeCosts(const vector<vector<double>>& basicPartition, 
    const vector<int>& basicColumns, 
    const vector<double>& costVector,
    std::vector<std::vector<double>>& nonBasicPartition,
    std::vector<int>& nonBasicColumns)
{
 //Calculate SMV
    vector<vector<double>> basicPartitionT = la::transpose(basicPartition);

    vector<double> basicCosts;
    for (int i = 0; i < basicColumns.size(); i++)
    {
        basicCosts.push_back(costVector[basicColumns[i]]);
    }

    for (int row = 0; row < basicPartitionT.size(); row++)
    {
        basicPartitionT[row].push_back(basicCosts[row]);
    }
    
    vector<double> lambda;

    int result = la::gauss(basicPartitionT, lambda);
    
    if(result == INF)
    {
        relativeCost r;
        r.lowestNonBColumn = -1;
        r.lowestRC = INF;
        cout << "Infeasible Solution in Phase 2. \n";
        return r;
    }

 //Calculate relative costs
    vector<double> relativeC;

    for (int col = 0; col < nonBasicColumns.size(); ++col) //Non basic columns are the indices of the non basic costVector
    {
        double cost_j = costVector[nonBasicColumns[col]];
    
        double dotProduct = 0.0;
        for (int row = 0; row < lambda.size(); ++row)
        {
            dotProduct += lambda[row] * nonBasicPartition[row][col];
        }
    
        relativeC.push_back(cost_j - dotProduct);
    }
    /*
    cout << "Relative costs: \n";
    for (int i = 0; i < relativeC.size(); i++)
    {
        cout << relativeC[i] << " ";
    }
    cout << "\n";
    */
    
 
 //Choose variable to enter basis
     int mv = 0;
     double lowestRC = relativeC[0];

     for(int i = 0; i < relativeC.size(); ++i)
     {
        
        if(relativeC[i] < relativeC[mv])
        {
            mv = i;
            lowestRC = relativeC[i];
            
        }
    }

    relativeCost rc;
    rc.lowestNonBColumn = nonBasicColumns[mv];
    rc.lowestRC = lowestRC;
    return rc;
    
}

vector<double> simplexDirection(const vector<vector<double>>& basicPartition,
    const vector<int>& nonBasicColumns,
    const vector<vector<double>>& nonBasicPartition,
    int k)
{

vector<double> y;


int rowSize = basicPartition.size();

vector<vector<double>> augmentedBasicPartition = basicPartition;

for (int i = 0; i < rowSize; i++) // i == rows :p
{
    
augmentedBasicPartition[i].push_back(nonBasicPartition[i][k]);

}

la::gauss(augmentedBasicPartition, y);

return y;
}

void rebuildPartitions(const vector<vector<double>>& A,
                       const vector<int>& basicCols,
                       const vector<int>& nonBasicCols,
                       vector<vector<double>>& basicPart,
                       vector<vector<double>>& nonBasicPart)
{
    basicPart.clear();
    nonBasicPart.clear();
    int n = A.size();
    for (int i = 0; i < n; i++) {
        vector<double> bRow, nbRow;
        for (int j : basicCols) {
            bRow.push_back(A[i][j]);
        }
        for (int j : nonBasicCols) {
            nbRow.push_back(A[i][j]);
        }
        basicPart.push_back(bRow);
        nonBasicPart.push_back(nbRow);
    }
}


vector<double> basicSolution(vector<vector<double>> basicPartition, const vector<int>& basicColumns, const vector<double>& costVector, const vector<double>& resources)
{
    vector<double> xB;
    int rowSize = basicPartition.size();
    for (int row = 0; row < rowSize; row++)
    {
        basicPartition[row].push_back(resources[row]);
    }
    la::gauss(basicPartition, xB);
    vector<double> basicCosts;
    for (int i = 0; i < basicColumns.size(); i++)
    {
        basicCosts.push_back(costVector[basicColumns[i]]);
    }
    double functionAnswer = 0;
    for (int i = 0; i < xB.size(); i++)
    {
        functionAnswer += xB[i]*basicCosts[i];
    }
    /*
    for (int i = 0; i < basicCosts.size(); i++)
    {
        cout << basicCosts[i] << " ";
    }
    */
    
        cout << "f(x) = " << functionAnswer << '\n';
    return xB;
}

int outBasis(const vector<double>& y, const vector<double>& xB, const vector<int>& basicColumns)
{
    vector<pair<double, int>> ratios;

    for(int i = 0; i < y.size(); i++)
    {
        if(y[i] > 0)
        {
            double ratio = xB[i] / y[i];
            ratios.push_back({ratio, i});
        }
    }

    if (ratios.empty()) {
        cout << "No valid leaving variable found (unbounded).\n";
    }
    /*
    cout << "Ratio values = ";
    for (int i = 0; i < ratios.size(); i++)
    {
        cout << ratios[i].first << ' ';
    }
    cout << "\n";
    */

    auto minPair = *min_element(ratios.begin(), ratios.end(),
                                [](const pair<double,int>& a, const pair<double,int>& b) {
                                    return a.first < b.first;
                                });

    int leavingIndex = minPair.second;
    /*
    cout << "y = ";
    for (int i = 0; i < y.size(); i++)
    cout << y[i] << " ";
    cout << "\n";
    
    cout << "xB = ";
    for (int i = 0; i < xB.size(); i++)
    cout << xB[i] << " ";
    cout << "\n";
    
    cout << "Minimum ratio: " << minPair.first << " at basicColumns index: " << leavingIndex << "\n";
    */

    return basicColumns[leavingIndex];
}


int updateBasis(vector<int>& basicColumns, vector<int>& nonBasicColumns, int enteringCol, int leavingColumn)
{
    for (int i = 0; i < basicColumns.size(); i++)
    {
        if(basicColumns[i] == leavingColumn)
        {
            basicColumns[i] = enteringCol;
        }

    }
    for (int i = 0; i < nonBasicColumns.size(); i++)
    {
        if(nonBasicColumns[i] == enteringCol)
        {
            nonBasicColumns[i] = leavingColumn;
        }
    }
    
    return 1;
}

int phase2(const vector<vector<double>>& basicPartition,
     vector<int>& basicColumns,
     const vector<double>& costVector,
     const vector<vector<double>>& restrictionsM,
     vector<vector<double>>& nonBasicPartition,
     vector<int>& nonBasicColumns,
     const vector<double>& resources)
{
   relativeCost rc = relativeCosts(basicPartition, basicColumns, costVector, nonBasicPartition, nonBasicColumns);
    
    
    if(rc.lowestRC >= -EPS)
    {
        cout << "Exited. ";
        return -1;
    }
    vector<double> xB = basicSolution(basicPartition, basicColumns, costVector, resources);
    cout << "\nEnters Non basic column: " <<rc.lowestNonBColumn << "\nLowest Relative Cost: " << rc.lowestRC << "\n";
    
    int enteringCol = rc.lowestNonBColumn;

    int k = -1;
    for (int i = 0; i < nonBasicColumns.size(); ++i) {
        
        if (nonBasicColumns[i] == enteringCol) {
            k = i;
            break;
        }
    }
    /*
    cout << "k = " <<  k << " \n";
    */
    vector<double> y = simplexDirection(basicPartition, nonBasicColumns, nonBasicPartition, k);
    //Test for finite sol
    bool foundPositive = false;
    for (int i = 0; i < y.size(); i++)
    {
        if (y[i] > EPS)
        {
            foundPositive = true;
            break;
        }
    }

    if (!foundPositive)
    {
        cout << "Problem is unbounded (no positive y).\n";
        return -2;
    }
    int leavingColumn = outBasis(y, xB, basicColumns);
    cout << "Basic column that leaves basis: " << leavingColumn << "\n";
    updateBasis(basicColumns, nonBasicColumns, enteringCol, leavingColumn);
    //cout << "Returned 0\n";
    return 0;
}


vector<double> simplex(const int& restQtt, 
                       const int& varQtt, 
                       vector<double> obFunction, 
                       vector<double> resources, 
                       vector<vector<double>> restrictionsM)
{
    //Phase 1
    vector<vector<double>> basicPartition;
    vector<int> basicColumns;
    vector<int> nonBasicColumns;
    vector<vector<double>> nonBasicPartition;

    phase1(obFunction, restrictionsM, basicColumns, nonBasicColumns);

    //Phase 2
    int maxIterations = 10000; // safeguard
    int status = 0;
    int it = 0;
    while (status == 0 && maxIterations--)
    {

        
        cout << "---\nIteration " << it << "\n---\n";
        it++;
        rebuildPartitions(restrictionsM, basicColumns, nonBasicColumns, basicPartition, nonBasicPartition);

        status = phase2(basicPartition, basicColumns, obFunction, restrictionsM, nonBasicPartition, nonBasicColumns, resources);
    }

    if (status == -1)
    {
        cout << "Optimal solution found.\n";
    }
    else if (status == -2)
    {
        cout << "The problem is unbounded.\n";
        return {};
    }
    else if (maxIterations < 2)
    {
        cout << "Max iterations reached. Possibly cycling.\n";
        
    }
    rebuildPartitions(restrictionsM, basicColumns, nonBasicColumns, basicPartition, nonBasicPartition);
    vector<double> xB = basicSolution(basicPartition, basicColumns, obFunction, resources);

    vector<double> fullSolution(varQtt, 0.0);

    for (int i = 0; i < basicColumns.size(); ++i)
    {
        if (basicColumns[i] < varQtt)
        {
            fullSolution[basicColumns[i]] = xB[i];
        }
    }
 
    return fullSolution;
}

void receiveInput(vector<double>& costVector,
    vector<vector<double>>& A,
    vector<double>& resources,
    int& constraints,
    int& variables)
{
    cout << "Insert the number of variables (How many Xs/Size of c): ";
    cin >> variables;
    cout << "Insert the number of constraints (How many rows in A/ Size of b): ";
    cin >> constraints;
    cout << "Insert the matrix of coefficients:\n";
    A.resize(constraints, vector<double>(variables));
    costVector.resize(variables);
    resources.resize(constraints);
    for (int i = 0; i < constraints; i++)
    {
        for (int j = 0; j < variables; j++)
        {
            cin >> A[i][j];
        }
        
    }
    cout << "\nInsert c: \n";

    for (int i = 0; i < variables; i++)
    {
        cin >> costVector[i];
    }
    cout << "\nInsert b: \n";
    for (int i = 0; i < constraints; i++)
    {
        cin >> resources[i];
    }
    return;

}

int receiveInputDIMACS(const string& filename,
    vector<double>& costVector,
    vector<vector<double>>& A,
    vector<double>& resources,
    int& constraints,
    int& variables)
{
    ifstream file(filename);
    if (!file.is_open()) {
        cout << "Failed to open file." << endl;
        return;
    }

    string line;

    while (getline(file, line))
    {
        if(line[0] == 'c')
        {
            continue;
        }
    }
    

}

int main()
{
    
    vector<double> costVector{};
    vector<vector<double>> A{};
    vector<double> resources{};
    int constraints = -1;
    int variables = -1;
    receiveInput(costVector, A, resources, constraints, variables);

    cout << "-------\nMin of \n";
    for (int i = 0; i < costVector.size(); i++)
    {
        
        cout << costVector[i] << "x" << i;
        if(i != costVector.size() - 1)
        {
            cout << " + ";
        }
    }
    cout << "\n\n";
    cout << "Subjected to:\n\n";
    for (int i = 0; i < A.size(); i++)
    {
        for (int j = 0; j < A[i].size(); j++)
        {
            cout << A[i][j] << ' ';
        }

        cout << " = " << resources[i] << "\n";
        
    }
    cout << "\n---";

   
    auto start = chrono::high_resolution_clock::now(); 
    vector<double> solution = simplex(constraints, variables, costVector, resources, A);
    auto end = chrono::high_resolution_clock::now(); 

    if(solution.size() > 0) cout << "\nFinal solution:\n";
    for (int i = 0; i < solution.size(); ++i)
    {
        cout << "x[" << i << "] = " << solution[i] << "\n";
    }


    cout << "Simplex Method - Time elapsed: [" << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "ms]\n";
    return 0;
}


//TODO:
// - Testes e Experimentacao
// - Otimizar gauss (sera que isso é necessario mesmo? tipo, o codigo ja ta super rapido)