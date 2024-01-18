#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <bits/stdc++.h>

using namespace std;

void print(int n, int** mat)
{
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << mat[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

int* assignMatrix(int n)
{
    int* data = (int*)malloc(n * n * sizeof(int));
    return data;
}

int** assignMatrix2D(int n)
{
    int* data = (int*)malloc(n * n * sizeof(int));
    int** array = (int**)malloc(n * sizeof(int*));
    for (int i = 0; i < n; i++)
    {
        array[i] = &(data[n * i]);
    }
    return array;
}

void populateMat(int n, int*& mat)
{
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            mat[i * n + j] = rand() % 5;
        }
    }
}

void populateMat2D(int n, int** &mat)
{
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            mat[i][j] = rand() % 5;
        }
    }
}

int** openSL(int n, int** mat, int offseti, int offsetj)
{
    int m = n / 2;
    int** slice = assignMatrix2D(m);
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
    int** result = assignMatrix2D(n);
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
    int** result = assignMatrix2D(n);

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

void MatrixClear(int n, int* mat)
{
    free(mat);
}

void MatrixClear2D(int n, int** mat)
{
    free(mat[0]);
    free(mat);
}

__global__ void multiply(int* mat1, int* mat2, int* product, int n)
{
    int prod = blockIdx.x * blockDim.x + threadIdx.x;
    int i = prod / n;
    int j = prod % n;
    for (int k = 0; k < n; k++) {
        product[i * n + j] += mat1[i * n + k] * mat2[k * n + j];
    }
}

int** cudaNaive(int n, int** mat1, int** mat2)
{
    int* h_mat1 = assignMatrix(n);
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            h_mat1[i*n + j] = mat1[i][j];
        }
    }

    int* h_mat2 = assignMatrix(n);
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            h_mat2[i*n + j] = mat2[i][j];
        }
    }

    int* h_product = assignMatrix(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            h_product[i * n + j] = 0;
        }
    }

    size_t bytes = n * n * sizeof(int);

    int *d_mat1, *d_mat2, *d_product;

    cudaMalloc(&d_mat1, bytes);
    cudaMalloc(&d_mat2, bytes);
    cudaMalloc(&d_product, bytes);

    cudaMemcpy(d_mat1, h_mat1, bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_mat2, h_mat2, bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_product, h_product, bytes, cudaMemcpyHostToDevice);

    int threads = min(1024, n);
    int blocks = (n * n) / threads;
    dim3 gridSize(blocks, 1, 1);
    dim3 blockSize(threads, 1, 1);

    multiply<<<gridSize, blockSize>>>(d_mat1, d_mat2, d_product, n);
    cudaDeviceSynchronize();

    cudaMemcpy(h_product, d_product, bytes, cudaMemcpyDeviceToHost);

    cudaFree(d_mat1);
    cudaFree(d_mat2);
    cudaFree(d_product);

    MatrixClear(n, h_mat1);
    MatrixClear(n, h_mat2);

    int** product = assignMatrix2D(n);
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            product[i][j] = h_product[i*n + j];
        }
    }
    return product;
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

    int** bds = MatSum(m, b, d, false);
    int** gha = MatSum(m, g, h, true);
    int** s1 = cudaNaive(m, bds, gha);
    MatrixClear2D(m, bds);
    MatrixClear2D(m, gha);

    int** ada = MatSum(m, a, d, true);
    int** eha = MatSum(m, e, h, true);
    int** s2 = cudaNaive(m, ada, eha);
    MatrixClear2D(m, ada);
    MatrixClear2D(m, eha);

    int** acs = MatSum(m, a, c, false);
    int** efa = MatSum(m, e, f, true);
    int** s3 = cudaNaive(m, acs, efa);
    MatrixClear2D(m, acs);
    MatrixClear2D(m, efa);

    int** aba = MatSum(m, a, b, true);
    int** s4 = cudaNaive(m, aba, h);
    MatrixClear2D(m, aba);
    MatrixClear2D(m, b);

    int** fhs = MatSum(m, f, h, false);
    int** s5 = cudaNaive(m, a, fhs);
    MatrixClear2D(m, fhs);
    MatrixClear2D(m, a);
    MatrixClear2D(m, f);
    MatrixClear2D(m, h);

    int** ges = MatSum(m, g, e, false);
    int** s6 = cudaNaive(m, d, ges);
    MatrixClear2D(m, ges);
    MatrixClear2D(m, g);

    int** cda = MatSum(m, c, d, true);
    int** s7 = cudaNaive(m, cda, e);
    MatrixClear2D(m, cda);
    MatrixClear2D(m, c);
    MatrixClear2D(m, d);
    MatrixClear2D(m, e);

    int** s1s2a = MatSum(m, s1, s2, true);
    int** s6s4s = MatSum(m, s6, s4, false);
    int** c11 = MatSum(m, s1s2a, s6s4s, true);
    MatrixClear2D(m, s1s2a);
    MatrixClear2D(m, s6s4s);
    MatrixClear2D(m, s1);

    int** c12 = MatSum(m, s4, s5, true);
    MatrixClear2D(m, s4);

    int** c21 = MatSum(m, s6, s7, true);
    MatrixClear2D(m, s6);

    int** s2s3s = MatSum(m, s2, s3, false);
    int** s5s7s = MatSum(m, s5, s7, false);
    int** c22 = MatSum(m, s2s3s, s5s7s, true);
    MatrixClear2D(m, s2s3s);
    MatrixClear2D(m, s5s7s);
    MatrixClear2D(m, s2);
    MatrixClear2D(m, s3);
    MatrixClear2D(m, s5);
    MatrixClear2D(m, s7);

    int** prod = MatJoin(m, c11, c12, c21, c22);

    MatrixClear2D(m, c11);
    MatrixClear2D(m, c12);
    MatrixClear2D(m, c21);
    MatrixClear2D(m, c22);

    return prod;
}

int main()
{
    int n;
    cout << "\nEnter Matrix Size in power of 2 (2,4,8,16,32,64,128,256,....16384): ";
    cin >> n;

    //n = 1024

    int** mat1 = assignMatrix2D(n);
    populateMat2D(n, mat1);

    int** mat2 = assignMatrix2D(n);
    populateMat2D(n, mat2);

    clock_t start, end;
    start = clock();

    int** prod = strassen(n, mat1, mat2);

    end = clock();
    double time = double(end - start) / double(CLOCKS_PER_SEC);
    cout << endl;
    cout << "\nParallel CUDA runtime with matrix " << n << " X " << n << " is: " << time;
    cout << endl;
    cout << endl;
    return 0;
}
