#include <bits/stdc++.h>

using namespace std;

void print(int n, int** mat)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << mat[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

int** assignMatrix(int n)
{
    int* data = (int*)malloc(n * n * sizeof(int));
    int** array = (int**)malloc(n * sizeof(int*));
    for (int i = 0; i < n; i++)
    {
        array[i] = &(data[n * i]);
    }
    return array;
}

void populateMat(int n, int**& mat)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            mat[i][j] = rand() % 5;
        }
    }
}

void MatrixClear(int n, int** mat)
{
    free(mat[0]);
    free(mat);
}

int** naiveMatMult(int n, int** mat1, int** mat2)
{
    int** prod = assignMatrix(n);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            prod[i][j] = 0;
            for (int k = 0; k < n; k++)
            {
                prod[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }

    return prod;
}

int** openSL(int n, int** mat, int offseti, int offsetj)
{
    int m = n / 2;
    int** slice = assignMatrix(m);
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < m; j++)
        {
            slice[i][j] = mat[offseti + i][offsetj + j];
        }
    }
    return slice;
}

int** MatSum(int n, int** mat1, int** mat2, bool add)
{
    int** result = assignMatrix(n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (add)
                result[i][j] = mat1[i][j] + mat2[i][j];
            else
                result[i][j] = mat1[i][j] - mat2[i][j];
        }
    }

    return result;
}

int** MatJoin(int m, int** c11, int** c12, int** c21, int** c22)
{
    int n = 2 * m;
    int** result = assignMatrix(n);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i < m && j < m)
                result[i][j] = c11[i][j];
            else if (i < m)
                result[i][j] = c12[i][j - m];
            else if (j < m)
                result[i][j] = c21[i - m][j];
            else
                result[i][j] = c22[i - m][j - m];
        }
    }

    return result;
}

int** StrassenMatMult(int n, int** mat1, int** mat2)
{

    if (n <= 32)
    {
        return naiveMatMult(n, mat1, mat2);
    }

    int m = n / 2;

    int** a = openSL(n, mat1, 0, 0);
    int** b = openSL(n, mat1, 0, m);
    int** c = openSL(n, mat1, m, 0);
    int** d = openSL(n, mat1, m, m);
    int** e = openSL(n, mat2, 0, 0);
    int** f = openSL(n, mat2, 0, m);
    int** g = openSL(n, mat2, m, 0);
    int** h = openSL(n, mat2, m, m);

    int** bds = MatSum(m, b, d, false);
    int** gha = MatSum(m, g, h, true);
    int** s1 = StrassenMatMult(m, bds, gha);
    MatrixClear(m, bds);
    MatrixClear(m, gha);

    int** ada = MatSum(m, a, d, true);
    int** eha = MatSum(m, e, h, true);
    int** s2 = StrassenMatMult(m, ada, eha);
    MatrixClear(m, ada);
    MatrixClear(m, eha);

    int** acs = MatSum(m, a, c, false);
    int** efa = MatSum(m, e, f, true);
    int** s3 = StrassenMatMult(m, acs, efa);
    MatrixClear(m, acs);
    MatrixClear(m, efa);

    int** aba = MatSum(m, a, b, true);
    int** s4 = StrassenMatMult(m, aba, h);
    MatrixClear(m, aba);
    MatrixClear(m, b);

    int** fhs = MatSum(m, f, h, false);
    int** s5 = StrassenMatMult(m, a, fhs);
    MatrixClear(m, fhs);
    MatrixClear(m, a);
    MatrixClear(m, f);
    MatrixClear(m, h);

    int** ges = MatSum(m, g, e, false);
    int** s6 = StrassenMatMult(m, d, ges);
    MatrixClear(m, ges);
    MatrixClear(m, g);

    int** cda = MatSum(m, c, d, true);
    int** s7 = StrassenMatMult(m, cda, e);
    MatrixClear(m, cda);
    MatrixClear(m, c);
    MatrixClear(m, d);
    MatrixClear(m, e);

    int** s1s2a = MatSum(m, s1, s2, true);
    int** s6s4s = MatSum(m, s6, s4, false);
    int** c11 = MatSum(m, s1s2a, s6s4s, true);
    MatrixClear(m, s1s2a);
    MatrixClear(m, s6s4s);
    MatrixClear(m, s1);

    int** c12 = MatSum(m, s4, s5, true);
    MatrixClear(m, s4);

    int** c21 = MatSum(m, s6, s7, true);
    MatrixClear(m, s6);

    int** s2s3s = MatSum(m, s2, s3, false);
    int** s5s7s = MatSum(m, s5, s7, false);
    int** c22 = MatSum(m, s2s3s, s5s7s, true);
    MatrixClear(m, s2s3s);
    MatrixClear(m, s5s7s);
    MatrixClear(m, s2);
    MatrixClear(m, s3);
    MatrixClear(m, s5);
    MatrixClear(m, s7);

    int** prod = MatJoin(m, c11, c12, c21, c22);

    MatrixClear(m, c11);
    MatrixClear(m, c12);
    MatrixClear(m, c21);
    MatrixClear(m, c22);

    return prod;
}

int main()
{
    int n;
    cout << "\nEnter Matrix Size in power of 2 (2,4,8,16,32,64,128,256,....16384): ";
    cin >> n;

    int** mat1 = assignMatrix(n);
    populateMat(n, mat1);

    int** mat2 = assignMatrix(n);
    populateMat(n, mat2);

    clock_t start, end;
    start = clock();
    int** prod1 = naiveMatMult(n, mat1, mat2);
    end = clock();
    double time = double(end - start) / double(CLOCKS_PER_SEC);
    cout << "\nSequential runtime of " << n << " X " << n << " matrix multiplication is: " << time << " seconds\n";

    start = clock();
    int** prod2 = StrassenMatMult(n, mat1, mat2);
    end = clock();
    time = double(end - start) / double(CLOCKS_PER_SEC);
    cout << "\nStrassen Sequential runtime of " << n << " X " << n << " matrix multiplication is: " << time << " seconds\n";

    cout << endl;
    return 0;
}
