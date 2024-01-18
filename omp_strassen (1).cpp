#include <omp.h>
#include <bits/stdc++.h>

using namespace std;

void print(int n, int** mat)
{
	int i;
	int j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
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

void populateMatrix(int n, int**& mat)
{
	int i;
	int j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
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
    int i;
	int j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
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
	int i;
	int j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
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

int** strassen(int n, int** mat1, int** mat2)
{

    int m = n / 2;

    int** a = openSL(n, mat1, 0, 0);
    int** b = openSL(n, mat1, 0, m);
    int** c = openSL(n, mat1, m, 0);
    int** d = openSL(n, mat1, m, m);
    int** e = openSL(n, mat2, 0, 0);
    int** f = openSL(n, mat2, 0, m);
    int** g = openSL(n, mat2, m, 0);
    int** h = openSL(n, mat2, m, m);

    int** s1;
    #pragma omp task shared(s1)
    {
        int** bds = MatSum(m, b, d, false);
        int** gha = MatSum(m, g, h, true);
        s1 = strassen(m, bds, gha);
        MatrixClear(m, bds);
        MatrixClear(m, gha);
    }

    int** s2;
    #pragma omp task shared(s2)
    {
        int** ada = MatSum(m, a, d, true);
        int** eha = MatSum(m, e, h, true);
        s2 = strassen(m, ada, eha);
        MatrixClear(m, ada);
        MatrixClear(m, eha);
    }

    int** s3;
    #pragma omp task shared(s3)
    {
        int** acs = MatSum(m, a, c, false);
        int** efa = MatSum(m, e, f, true);
        s3 = strassen(m, acs, efa);
        MatrixClear(m, acs);
        MatrixClear(m, efa);
    }

    int** s4;
    #pragma omp task shared(s4)
    {
        int** aba = MatSum(m, a, b, true);
        s4 = strassen(m, aba, h);
        MatrixClear(m, aba);
    }

    int** s5;
    #pragma omp task shared(s5)
    {
        int** fhs = MatSum(m, f, h, false);
        s5 = strassen(m, a, fhs);
        MatrixClear(m, fhs);
    }

    int** s6;
    #pragma omp task shared(s6)
    {
        int** ges = MatSum(m, g, e, false);
        s6 = strassen(m, d, ges);
        MatrixClear(m, ges);
    }

    int** s7;
    #pragma omp task shared(s7)
    {
        int** cda = MatSum(m, c, d, true);
        s7 = strassen(m, cda, e);
        MatrixClear(m, cda);
    }

    #pragma omp taskwait

    MatrixClear(m, a);
    MatrixClear(m, b);
    MatrixClear(m, c);
    MatrixClear(m, d);
    MatrixClear(m, e);
    MatrixClear(m, f);
    MatrixClear(m, g);
    MatrixClear(m, h);

    int** c11;
    #pragma omp task shared(c11)
    {
        int** s1s2a = MatSum(m, s1, s2, true);
        int** s6s4s = MatSum(m, s6, s4, false);
        c11 = MatSum(m, s1s2a, s6s4s, true);
        MatrixClear(m, s1s2a);
        MatrixClear(m, s6s4s);
    }

    int** c12;
    #pragma omp task shared(c12)
    {
        c12 = MatSum(m, s4, s5, true);
    }

    int** c21;
    #pragma omp task shared(c21)
    {
        c21 = MatSum(m, s6, s7, true);
    }

    int** c22;
    #pragma omp task shared(c22)
    {
        int** s2s3s = MatSum(m, s2, s3, false);
        int** s5s7s = MatSum(m, s5, s7, false);
        c22 = MatSum(m, s2s3s, s5s7s, true);
        MatrixClear(m, s2s3s);
        MatrixClear(m, s5s7s);
    }

    #pragma omp taskwait

    MatrixClear(m, s1);
    MatrixClear(m, s2);
    MatrixClear(m, s3);
    MatrixClear(m, s4);
    MatrixClear(m, s5);
    MatrixClear(m, s6);
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
    populateMatrix(n, mat1);

    int** mat2 = assignMatrix(n);
    populateMatrix(n, mat2);

    double startParStrassen = omp_get_wtime();
    int** prod;

    omp_set_num_threads(8);

    #pragma omp parallel
    {
    #pragma omp single
        {
            prod = strassen(n, mat1, mat2);
        }
    }
    double endParStrassen = omp_get_wtime();
    cout << "\nParallel runtime with matrix " << n << " X " << n << " is: ";
    cout << setprecision(5) << endParStrassen - startParStrassen << endl;

    cout << endl;

    return 0;
}
