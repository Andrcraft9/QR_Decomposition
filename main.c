#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

typedef double** matrix;


matrix create(size_t m, size_t n)
{
	matrix mat;
		
	mat = malloc(m * sizeof(double*));
	for(size_t i = 0; i < m; ++i)
	{
		mat[i] = malloc(n * sizeof(double));
		for(size_t j = 0; j < n; ++j)
			mat[i][j] = 0;
	}
	
	return mat;
}


matrix scan_matrix(size_t m, size_t n)
{
	matrix mat;
	
	mat = malloc(m * sizeof(double*));
	for(size_t i = 0; i < m; ++i)
	{
		mat[i] = malloc(n * sizeof(double));
		for(size_t j = 0; j < n; ++j)
		{
			double d;
			scanf("%lf", &d);		
			mat[i][j] = d;
		}	
	}
	
	return mat;
}

matrix create_seq_matrix(size_t n)
{
	matrix mat;
	int value = 0;
	mat = malloc(n * sizeof(double*));
	for(size_t i = 0; i < n; ++i)
	{
		mat[i] = malloc(n * sizeof(double));
		for(size_t j = 0; j < n; ++j)
			mat[i][j] = ++value;
	}
	
    return mat;
}

void print_matrix(matrix mat, size_t m, size_t n)
{
	for(size_t i = 0; i < m; ++i)
	{
		for(size_t j = 0; j < n; ++j)
			printf("%lf ", mat[i][j]);
		printf("\n");
	}
}

void delete_matrix(matrix mat, size_t m)
{
	for(size_t i = 0; i < m; ++i)
	{
		free(mat[i]);
	}
	free(mat);
}

matrix create_from(matrix matrix1, size_t m, size_t n)
{
    matrix matrix2;

 	matrix2 = malloc(m * sizeof(double*));
	for(size_t i = 0; i < m; ++i)
	{
		matrix2[i] = malloc(n * sizeof(double));
		for(size_t j = 0; j < n; ++j)
			matrix2[i][j] = matrix1[i][j];
	}   

    return matrix2;
}

matrix create_I(size_t m, size_t n)
{
    matrix matrix2;

 	matrix2 = malloc(m * sizeof(double*));
	for(size_t i = 0; i < m; ++i)
	{
		matrix2[i] = malloc(n * sizeof(double));
		for(size_t j = 0; j < n; ++j)
			if ( i == j)
                matrix2[i][j] = 1;
            else
                matrix2[i][j] = 0;
	}   

    return matrix2;
}

matrix* QRdecomposition(const matrix A, size_t N)
{
	matrix* QR = malloc( 2 * sizeof(matrix));
    
    matrix Q, R;
	Q = create_I(N, N);
	R = create_from(A, N, N);
    QR[0] = Q;
    QR[1] = R;

    double *vT, *vrT, *vq;
	vT = malloc(N * sizeof(double));
    vrT = malloc(N * sizeof(double));
    vq = malloc(N * sizeof(double));

    // A = Q*R
	// P = I - (2/vT*v) * v*vT - householder matrix, v - householder 
    // Pn*Pn-1*...*P1*A = R
    // P1*P2*...*Pn = Q
    for(size_t n = 0; n < N; ++n)
	{
        // Compute Rn = Pn * Rn-1 and Qn = Qn-1 * Pn
		
        // Set size our vectors
		size_t size = N - n;

		// Compute Householder vector vT
		double w1 = 0;
		for(size_t i = 0; i < size; ++i)
		{
			vT[i] = R[n+i][n];
			w1 += vT[i] * vT[i];
		}	
		w1 = sqrt(w1);
		if (R[n][n] > 0)
			vT[0] = vT[0] + w1;
		else
			vT[0] = vT[0] - w1;

		// Multiply vrT = vT * R and compute w2 = (v, v) and compute vq = Q * v
		double w2 = 0;
		for(size_t i = 0; i < size; ++i)
		{
			// compute w2
            w2 += vT[i] * vT[i];

            // compute vrT and vq
			vrT[i] = 0;
            vq[i] = 0;
			for(size_t s = 0; s < size; ++s)
            {
				vrT[i] += vT[s] * R[n + s][n + i]; 
		        vq[i] += Q[n + i][n + s] * vT[s];
            }
        }
		
		// Compute w2 = 2/(v, v)
		w2 = 2 / w2;

		// Compute new R = R - w2(v * vrT) and new Q = Q - w2(vq * vT)
		for(size_t i = 0; i < size; ++i)
			for(size_t j = 0; j < size; ++j)
            {
				R[n + i][n + j] = R[n + i][n + j] - w2 * vT[i] * vrT[j];
	            Q[n + i][n + j] = Q[n + i][n + j] - w2 * vq[i] * vT[j];
            }
    
    }
    
    return QR;
}

int main(void)
{
	size_t size;
    printf("Write dimension\n");
    scanf("%zu", &size);
    printf("\n");

    matrix mat1;
	matrix* QR;
	// Generating matrix
    mat1 = scan_matrix(size, size);
    // mat1 = create_seq_matrix(size);

    clock_t start = clock();
    QR = QRdecomposition(mat1, size);
    clock_t end = clock();
    printf("%lf\n",(double) end - start);

    printf("Q: \n");
	print_matrix(QR[0], size, size);
    printf("R: \n");
    print_matrix(QR[1], size, size);
    
    printf("new A:\n ");
    double a;
	for(size_t i = 0; i < size; ++i)
	{
		for(size_t j = 0; j < size; ++j)
        {
			a = 0;
            for(size_t s = 0; s < size; ++s)
                a += QR[0][i][s] * QR[1][s][j];
            printf("%lf ", a);
        }
        printf("\n");
	}

	//delete_matrix(mat1, size);
	
	return 0;
}




