#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <getopt.h>

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
            if (scanf("%lf", &d) != 1)
                exit(1);
            mat[i][j] = d;
        }
    }

    return mat;
}

matrix create_rand_matrix(size_t m, size_t n)
{
    matrix mat;

    srand(time(0));

    mat = malloc(m * sizeof(double*));
    for(size_t i = 0; i < m; ++i)
    {
        mat[i] = malloc(n * sizeof(double));
        for(size_t j = 0; j < n; ++j)
            mat[i][j] = (double) rand() / RAND_MAX;
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

void print_matrix(const matrix mat, size_t m, size_t n)
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

matrix create_from(const matrix matrix1, size_t m, size_t n)
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

matrix matrix_multiply_square(const matrix A, const matrix B, size_t n)
{
    matrix mat;
    mat = create(n, n);

    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
        {
            mat[i][j] = 0;
            for (size_t s = 0; s < n; ++s)
                mat[i][j] += A[i][s] * B[s][j];
        }

    return mat;

}

// Generate bad matrix Ak = H1 * diag(yk) * H2,
// where y1 / yn = 10**k and y1 > y2 ... yn-1 > yn,
// H1, H2 - Householder matrix
matrix create_rand_bad_matrix(size_t n, size_t k)
{
    srand(time(0));

    // Generate random Householder vector v1 and v2 and compute w = (v,v)
    double *v1, *v2;
    double w1 = 0, w2 = 0;
    v1 = malloc(n * sizeof(double));
    v2 = malloc(n * sizeof(double));
    for (size_t i = 0; i < n; ++i)
    {
        v1[i] = (double) rand() / RAND_MAX;
        w1 = w1 + v1[i] * v1[i];
        v2[i] = (double) rand() / RAND_MAX;
        w2 = w2 + v2[i] * v2[i];
    }

    w1 = 2.0 / w1;
    w2 = 2.0 / w2;

    // Create Householder matrix H1 and H2
    matrix H1, H2;
    H1 = create_I(n, n);
    H2 = create_I(n, n);
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
        {
            H1[i][j] = H1[i][j] - w1 * v1[i] * v1[j];
            H2[i][j] = H2[i][j] - w2 * v2[i] * v2[j];
        }

    // Create diag(yk)
    matrix Y;
    Y = create_I(n, n);
    double h = (double) k / (n - 1);
    double deg = 0;
    for (size_t i = 0; i < n; ++i)
    {
        Y[n - i - 1][n - i - 1] = pow(10, deg);
        deg = deg + h;
    }

    //print_matrix(Y, n, n);

    // Create Ak = H1 * Y * H2
    matrix Ak, P;
    P = matrix_multiply_square(H1, Y, n);
    Ak = matrix_multiply_square(P, H2, n);

    free(v1);
    free(v2);
    delete_matrix(Y, n);
    delete_matrix(P, n);
    delete_matrix(H1, n);
    delete_matrix(H2, n);

    //print_matrix(Ak, n, n);

    return Ak;

}

double matrix_Frobenius_norm(const matrix A, size_t m, size_t n)
{
    long double norm = 0;
    for(size_t i = 0; i < m; ++i)
        for(size_t j = 0; j < n; ++j)
            norm += A[i][j] * A[i][j];

    norm = sqrt(norm);

    return norm;
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
        // Compute Rn = Pn * Rn-1 and Qn = Qn-1 * PnT

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

        // Compute vrT = vT * R and compute w2 = (v, v)
        double w2 = 0;
        for(size_t i = 0; i < size; ++i)
        {
            // compute w2
            w2 += vT[i] * vT[i];

            // compute vrT
            vrT[i] = 0;
            for(size_t s = 0; s < size; ++s)
                vrT[i] += vT[s] * R[n + s][n + i];
        }

        //Compute vq = Q * v
        for(size_t i = 0; i < N; ++i)
        {
            vq[i] = 0;
            for(size_t s = 0; s < size; ++s)
                vq[i] += Q[i][n + s] * vT[s];
        }

        // Compute w2 = 2/(v, v)
        w2 = 2.0 / w2;

        // Compute new R = R - w2(v * vrT)
        for(size_t i = 0; i < size; ++i)
            for(size_t j = 0; j < size; ++j)
                R[n + i][n + j] = R[n + i][n + j] - w2 * vT[i] * vrT[j];

        // Compute new Q = Q - w2(vq * vT)
        for(size_t i = 0; i < N; ++i)
            for(size_t j = 0; j < size; ++j)
                Q[i][n + j] = Q[i][n + j] - w2 * vq[i] * vT[j];
    }

    free(vT);
    free(vrT);
    free(vq);

    return QR;
}

int main(int argc, char *argv[])
{
    struct globalArgs_t
    {
        int print_matrix;   // -a print matrix newA, Q, R, QQT
        int print_time;     // -t print time of QR decomposition computing
        int print_err;      // -e print error of computing
        int print_A;        // -A print matrix A
        int rand;           // -r random generating of matrix A
        char* bad_rand;       // -R k bad random generating of matrix A. k is degree
        char* size;   // -s size of matrix A

    } globalArgs;

    const char *optString = "ateArR:s:";

    // Initialization structure of arguments
    globalArgs.print_matrix = 0;
    globalArgs.print_time = 0;
    globalArgs.print_err = 0;
    globalArgs.print_A = 0;
    globalArgs.rand = 0;
    globalArgs.bad_rand = NULL;
    globalArgs.size = NULL;

    // Analyze command string
    int opt;
    opt = getopt(argc, argv, optString);
    while( opt != -1 )
    {
        switch( opt )
        {
            case 'a':
                globalArgs.print_matrix = 1; /* true */
                break;
            case 't':
                globalArgs.print_time = 1;
                break;
            case 'e':
                globalArgs.print_err = 1;
                break;
            case 'r':
                globalArgs.rand = 1;
                break;
            case 'R':
                globalArgs.bad_rand = optarg;
                break;
            case 'A':
                globalArgs.print_A = 1;
                break;
            case 's':
                globalArgs.size = optarg;
                break;
        }
        opt = getopt(argc, argv, optString);
    }

    // Set size of matrix
    size_t size;
    if (globalArgs.size != NULL)
        size = strtoul(globalArgs.size, NULL, 10);
    else
    {
        printf("Write dimension\n");
        if (scanf("%zu", &size) != 1)
            exit(1);
        printf("\n");
    }

    matrix A;
    matrix* QR;

    // Generating matrix A
    if (globalArgs.rand)
    {
        A = create_rand_matrix(size, size);
    }
    else if (globalArgs.bad_rand)
    {
        A = create_rand_bad_matrix(size, strtoul(globalArgs.bad_rand, NULL, 10));
    }
    else
    {
        A = scan_matrix(size, size);
    }

    // QR - decomposition
    clock_t start = clock();
    QR = QRdecomposition(A, size);
    clock_t end = clock();

    if (globalArgs.print_time)
        // Print time of computing
        printf("Time: %lf\n",(double) (end - start) / CLOCKS_PER_SEC );


    matrix newA = NULL;
    matrix QQT = NULL;

    if (globalArgs.print_matrix || globalArgs.print_err)
    {
        // Building new matrix A = Q * R
        newA = matrix_multiply_square(QR[0], QR[1], size);

        // Building QQT = Q * QT
        QQT = create(size, size);
        for (size_t i = 0; i < size; ++i)
            for (size_t j = 0; j < size; ++j)
            {
                QQT[i][j] = 0;
                for(size_t s = 0; s < size; ++s)
                    QQT[i][j] += QR[0][i][s] * QR[0][j][s];
            }
    }

    if (globalArgs.print_err)
    {
        // Compute errA = ||QR - A|| / ||A||; by Frobenius norm. Print error
        long double errA = 0;
        for(size_t i = 0; i < size; ++i)
            for(size_t j = 0; j < size; ++j)
                errA += (newA[i][j] - A[i][j]) * (newA[i][j] - A[i][j]);
        errA = sqrt(errA);
        errA = errA / matrix_Frobenius_norm(A, size, size);
        printf("Error of matrix A = Q * R: %.20Lf\n", errA);

        // Compute errQQT = ||QQT - I|| / ||I||
        //double errQQT = 0;
        //for(size_t i = 0; i < size; ++i)
        //    errQQT += (1 - QQT[i][i])
        // Compute errA = ||A|| - ||newA||; by Frobenius norm. Print error
        //double errA = fabs(matrix_Frobenius_norm(A, size, size) - matrix_Frobenius_norm(newA, size, size));
        //printf("Error of matrix A = Q * R: %.20lf\n", errA);

        // Compute errQQT = ||I|| - ||QQT||; by Frobenius norm. Print error
        //double errQQT = fabs(sqrt(size) - matrix_Frobenius_norm(QQT, size, size));
        //printf("Error of matrix QQT = Q * QT: %.20lf\n", errQQT);
    }

    if (globalArgs.print_matrix)
    {
        // Print all information
        printf("Q:\n");
        print_matrix(QR[0], size, size);
        printf("R:\n");
        print_matrix(QR[1], size, size);
        printf("new A:\n");
        print_matrix(newA, size, size);
        printf("QQT:\n");
        print_matrix(QQT, size, size);
    }

    if (globalArgs.print_A)
    {
        printf("A:\n");
        print_matrix(A, size, size);
    }

    // Clear memory
    delete_matrix(A, size);
    delete_matrix(QR[0], size);
    delete_matrix(QR[1], size);
    free(QR);
    if (newA)
        delete_matrix(newA, size);
    if (QQT)
        delete_matrix(QQT, size);

    return 0;
}




