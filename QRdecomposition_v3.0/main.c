#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <getopt.h>

#include <cblas.h>

#include <lapacke.h>

typedef double* matrix;

matrix create(size_t m, size_t n)
{
	matrix mat;

	mat = malloc(m * n * sizeof(double));
	for(size_t i = 0; i < m*n; ++i)
	{
        mat[i] = 0;
	}

	return mat;
}

matrix copy_matrix(const matrix A, size_t m, size_t n)
{
    matrix mat;
	mat = malloc(m * n * sizeof(double));

	for(size_t i = 0; i < m*n; ++i)
	{
		mat[i] = A[i];
	}

	return mat;

}

matrix create_rand_matrix(size_t m, size_t n)
{
	matrix mat;

    srand(time(0));

	mat = malloc(m * n * sizeof(double));
	for(size_t i = 0; i < m*n; ++i)
	{
		mat[i] = (double) rand() / RAND_MAX;
	}

	return mat;
}

void print_matrix(const matrix mat, size_t m, size_t n)
{
	for(size_t i = 0; i < m; ++i)
	{
		for(size_t j = 0; j < n; ++j)
			printf("%lf ", mat[i*n + j]);
		printf("\n");
	}
}

matrix create_I(size_t m, size_t n)
{
    matrix matrix2;

 	matrix2 = malloc(m * n * sizeof(double));
	for(size_t i = 0; i < m; ++i)
	{
		for(size_t j = 0; j < n; ++j)
			if (i == j)
                matrix2[i*n + j] = 1;
            else
                matrix2[i*n + j] = 0;
	}

    return matrix2;
}

matrix create_O(size_t m, size_t n)
{
    matrix mat;

 	mat = malloc(m * n * sizeof(double));
	for(size_t i = 0; i < m*n; ++i)
	{
        mat[i] = 0;
	}

    return mat;
}

double error_of_qr(const matrix A, const matrix newA, size_t n)
{
    double er = 0;
    double norm = 0;
    for(size_t i = 0; i < n*n; ++i)
    {
        er = er + (A[i] - newA[i]) * (A[i] - newA[i]);
        norm = norm + A[i] * A[i];
    }
    norm = sqrt(norm);
    er = sqrt(er) / norm;

    return er;
}

matrix* QRdecomposition(const matrix A, size_t N)
{
	matrix* QR = malloc( 2 * sizeof(matrix));

    double *vT, *vrT, *vq;
	vT = malloc(N * sizeof(double));
    vrT = malloc(N * sizeof(double));
    vq = malloc(N * sizeof(double));

    matrix Q, R, P;
    Q = create_I(N, N);
    R = copy_matrix(A, N, N);
    P = create_O(N, N);

    // A = Q*R
	// P = I - (2/vT*v) * v*vT - householder matrix, v - householder
    // Pn*Pn-1*...*P1*A = R
    // P1*P2*...*Pn = Q
    for(size_t k = 0; k < N; ++k)
	{
        size_t size;
        size = N - k;
        for(size_t i = 0; i < size; ++i)
            vT[i] = R[(k + i) * N + k];

        double w1;
        w1 = cblas_dnrm2(size, vT, 1);

        if (R[k * N + k] > 0)
            vT[0] = vT[0] + w1;
        else
            vT[0] = vT[0] - w1;

        cblas_dgemv(CblasRowMajor, CblasTrans, size, size, 1, R + k*N + k, N, vT, 1, 0, vrT, 1);
        double w2;
        w2 = cblas_ddot(size, vT, 1, vT, 1);

        cblas_dgemv(CblasRowMajor, CblasNoTrans, N, size, 1, Q + k, N, vT, 1, 0, vq, 1);

        w2 = 2.0 / w2;

        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, size, size, 1, -w2, vT, N, vrT, N, 1, R + k*N + k, N);

        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, N, size, 1, -w2, vq, N, vT, N, 1, Q + k, N);
    }

    QR[0] = Q;
    QR[1] = R;

    free(vT);
    free(vrT);
    free(vq);
    free(P);

    return QR;
}

int main(int argc, char *argv[])
{
    // Set size of matrix
    size_t size;
    size = 2048;

    matrix A;
	A = create_rand_matrix(size, size);

    // QR - decomposition
    matrix* QR;
    clock_t start = clock();
    QR = QRdecomposition(A, size);
    clock_t end = clock();
    printf("Time: %lf\n",(double) (end - start) / CLOCKS_PER_SEC );

    // Compute Error of QR decomposition
    matrix newA;
    newA = create(size, size);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, size, size, size, 1, QR[0], size, QR[1], size, 0, newA, size);
    printf("Error: %.16lf\n", error_of_qr(A, newA, size));

    free(A);
    free(newA);
    free(QR[0]);
    free(QR[1]);
    free(QR);

	return 0;
}




