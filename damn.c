#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#define MAX 1000000
#define eps 1e-10


void WhereIsX(int n, double x, double* arr_x, int* m) {
    double h_min = MAX;
    for (int i = 0; i < n; i++) {
        if (((x - arr_x[i]) < h_min) && (x - arr_x[i]) >= 0) {
            h_min = x - arr_x[i];
            *m = i;
        }
    }
}

void steps(double* arr_x, double* arr_y, double* h, int n, double* alpha, double* c, double* b, double* d) {
    double* l = (double*)calloc(n + 1, sizeof(*l));
    double* u = (double*)calloc(n + 1, sizeof(*u));
    double* z = (double*)calloc(n + 1, sizeof(*z));

    for (int i = 0; i <= n - 1; ++i)
        h[i] = arr_x[i + 1] - arr_x[i];

    for (int i = 1; i <= n - 1; i++)
        alpha[i] = (3 * (arr_y[i + 1] - arr_y[i]) / h[i]) - 3 * (arr_y[i] - arr_y[i - 1]) / h[i - 1];

    l[0] = 1;
    u[0] = 0;
    z[0] = 0;

    for (int i = 1; i <= n - 1; ++i) {
        l[i] = 2 * (arr_x[i + 1] - arr_x[i - 1]) - h[i - 1] * u[i - 1];
        u[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    l[n - 1] = 1;
    u[n - 1] = 0;
    z[n - 1] = 0;

    for (int j = n - 1; j >= 0; --j) {
        c[j] = z[j] - u[j] * c[j + 1];
        b[j] = (arr_y[j + 1] - arr_y[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3;
        d[j] = (c[j + 1] - c[j]) / (3 * h[j]);
    }
}

void spline(double* arr_y, double* b, double x, double* arr_x, double* c, double* d, int m, double* znach) {
    *znach = arr_y[m] + b[m] * (x - arr_x[m]) + c[m] * (x - arr_x[m]) * (x - arr_x[m]) + d[m] * (x - arr_x[m]) * (x - arr_x[m]) * (x - arr_x[m]);
}

void equal(double* arr_x1, double* arr_x2, double* arr_y1, double* arr_y2, int n, int k, int* g) {
    int count = 0;
    for (int i = 0; i < n + 1; i++) {
        for (int j = 0; j < k + 1; j++) {
            if ((arr_x1[i] == arr_x2[j]) && (arr_y1[i] == arr_y2[j]))
                count += 1;
        }
    }
    if (count == k + 1 || count == n + 1)
        *g = 1;
    else
        *g = 0;
}


struct {
    double first;
    double second;
}   root[3] = { eps, eps, eps, eps, eps, eps };


double croot(double x) {
    if (x < 0)
        return -pow(-x, 1.0 / 3.0);
    return pow(x, 1.0 / 3.0);
}

void square(double a, double b, double c) {
    double discr = b * b - 4 * a * c;

    if (discr < 0) {
        for (int i = 0; i < 3; i++)
            root[i].second = eps;
    }

    else if (discr == 0) {
        root[0].first = ((-b + sqrt(discr)) / (2 * a));
        root[0].second = 0;
        for (int i = 1; i < 3; i++) { root[i].second = eps; }
    }

    else {
        printf("x = %f and x = %f \n", ((-b + sqrt(discr)) / (2 * a)), ((-b - sqrt(discr)) / (2 * a)));
        for (int i = 0; i < 2; i++) {
            root[i].second = 0;
            root[i].first = ((-b + pow((-1), i) * sqrt(discr)) / (2 * a));
        }
        root[2].second = eps;
    }
}


int main() {
    double x1, x2, left, right, y1, y2;
    int m1, m2, n, k;
    printf("Hi!\nEnter the number of points of first spline ");
    scanf("%d", &n);
    n = n - 1;
    printf("Enter the number of points of second spline ");
    scanf("%d", &k);
    k = k - 1;

    double* arr_x1 = (double*)calloc(n + 1, sizeof(*arr_x1));
    double* arr_y1 = (double*)calloc(n + 1, sizeof(*arr_y1));
    double* h1 = (double*)calloc(n, sizeof(*h1));
    double* alpha1 = (double*)calloc(n, sizeof(*alpha1));
    double* c1 = (double*)calloc(n + 1, sizeof(*c1));
    double* b1 = (double*)calloc(n, sizeof(*b1));
    double* d1 = (double*)calloc(n, sizeof(*d1));

    double* arr_x2 = (double*)calloc(k + 1, sizeof(*arr_x2));
    double* arr_y2 = (double*)calloc(k + 1, sizeof(*arr_y2));
    double* h2 = (double*)calloc(k, sizeof(*h2));
    double* alpha2 = (double*)calloc(k, sizeof(*alpha2));
    double* c2 = (double*)calloc(k + 1, sizeof(*c2));
    double* b2 = (double*)calloc(k, sizeof(*b2));
    double* d2 = (double*)calloc(k, sizeof(*d2));

    printf("Enter x of first spline: ");
    for (int i = 0; i < n + 1; i++)
        scanf("%lf", &arr_x1[i]);

    printf("Enter x of second spline: ");
    for (int i = 0; i < k + 1; i++)
        scanf("%lf", &arr_x2[i]);

    printf("Enter y of first spline: ");
    for (int i = 0; i < n + 1; ++i)
        scanf("%lf", &arr_y1[i]);

    printf("Enter y of second spline: ");
    for (int i = 0; i < k + 1; ++i)
        scanf("%lf", &arr_y2[i]);

    int flag1 = 0;
    printf("\nEnter the desired point of first spline: ");
    scanf("%lf", &x1);

    // check if the desired point is already set
    for (int i = 0; i < n + 1; i++) {
        if (arr_x1[i] == x1) {
            printf("\nS1(%f) = %f\n", x1, arr_y1[i]);
            flag1 = 1;
            break;
        }
    }
    //next to what is the given point ?
    WhereIsX(n + 1, x1, arr_x1, &m1);

    int flag2 = 0;
    printf("Enter the desired point of second spline: ");
    scanf("%lf", &x2);

    // check if the desired point is already set
    for (int i = 0; i < k + 1; i++) {
        if (arr_x2[i] == x2) {
            printf("\nS2(%f) = %f\n", x2, arr_y2[i]);
            flag2 = 1;
            break;
        }
    }
    //next to what is the given point ?
    WhereIsX(k + 1, x2, arr_x2, &m2);

    //coefficient search
    steps(arr_x1, arr_y1, h1, n, alpha1, c1, b1, d1);
    steps(arr_x2, arr_y2, h2, k, alpha2, c2, b2, d2);

    printf("\n--------------------coefficients--------------------\n");

    printf("%2s %8s %8s %8s %8s of 1st spline\n", "i", "ai", "bi", "ci", "di");
    for (int i = 0; i <= n - 1; i++)
        printf("%2d %8.2f %8.2f %8.2f %8.2f\n", i, arr_y1[i], b1[i], c1[i], d1[i]);
    printf("\n");
    printf("%2s %8s %8s %8s %8s of 2nd spline\n", "i", "ai", "bi", "ci", "di");
    for (int i = 0; i <= k - 1; i++)
        printf("%2d %8.2f %8.2f %8.2f %8.2f\n", i, arr_y2[i], b2[i], c2[i], d2[i]);


    if (flag1 == 0) {
        spline(arr_y1, b1, x1, arr_x1, c1, d1, m1, &y1);
        printf("\nS1(%f)= %f \n", x1, y1);
    }
    if (flag2 == 0) {
        spline(arr_y2, b2, x2, arr_x2, c2, d2, m2, &y2);
        printf("\nS2(%f)= %f \n", x2, y2);
    }

    //do the splines match?
    int g;
    equal(arr_x1, arr_x2, arr_y1, arr_y2, n, k, &g);
    if (g == 1) {
        printf("\nSplines are equivalent\n");
        return 0;
    }

    //search for an intersection
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            left = fmax(arr_x1[i], arr_x2[j]);
            right = fmin(arr_x1[i + 1], arr_x2[j + 1]);
            if (left < right) {
                double A = d1[i] - d2[j];
                double B = c1[i] - 3 * d1[i] * arr_x1[i] - c2[j] + 3 * d2[j] * arr_x2[j];
                double C = b1[i] - 2 * c1[i] * arr_x1[i] + 3 * d1[i] * arr_x1[i] * arr_x1[i] - b2[j] + 2 * c2[j] * arr_x2[j] - 3 * d2[j] * arr_x2[j] * arr_x2[j];
                double D = arr_y1[i] - b1[i] * arr_x1[i] + c1[i] * arr_x1[i] * arr_x1[i] - d1[i] * arr_x1[i] * arr_x1[i] * arr_x1[i] - arr_y2[j] + b2[j] * arr_x2[j] - c2[j] * arr_x2[j] * arr_x2[j] + d2[j] * arr_x2[j] * arr_x2[j] * arr_x2[j];

                if (A == 0 && B != 0)
                    square(B, C, D);

                else if (A == 0 && B == 0) {// Cx + D = 0
                    root[0].first = (-D / C);
                    root[0].second = 0;
                }
                else {//cubic (Cardano)
                    double p = (3.0 * A * C - B * B) / (3.0 * A * A);
                    double q = (2.0 * B * B * B - 9.0 * A * B * C + 27.0 * A * A * D) / (27.0 * A * A * A);
                    double S = (q * q / 4.0) + (p * p * p / 27.0);

                    double F;
                    if (q == 0)
                        F = M_PI / 2.0;
                    if (q < 0)
                        F = atan(-2.0 * sqrt(-S) / q);
                    if (q > 0)
                        F = atan(-2.0 * sqrt(-S) / q) + M_PI;

                    for (int i = 0; i < 3; i++)
                        root[i].first = root[i].second = 0;

                    if (S < 0) {
                        root[0].first = 2.0 * sqrt(-p / 3.0) * cos(F / 3.0) - B / (3.0 * A);
                        root[1].first = 2.0 * sqrt(-p / 3.0) * cos((F / 3.0) + 2.0 * M_PI / 3.0) - B / (3.0 * A);
                        root[2].first = 2.0 * sqrt(-p / 3.0) * cos((F / 3.0) + 4.0 * M_PI / 3.0) - B / (3.0 * A);
                    }

                    if (S == 0) {
                        root[0].first = 2.0 * croot(-q / 2.0) - B / (3.0 * A);
                        root[1].first = -croot(-q / 2.0) - B / (3.0 * A);
                        root[2].first = -croot(-q / 2.0) - B / (3.0 * A);
                    }

                    if (S > 0) {
                        double temp1 = croot((-q / 2.0) + sqrt(S)) + croot((-q / 2.0) - sqrt(S));
                        double temp2 = croot((-q / 2.0) + sqrt(S)) - croot((-q / 2.0) - sqrt(S));
                        root[0].first = temp1 - B / (3.0 * A);
                        root[1].first = -temp1 / 2.0 - B / (3.0 * A);
                        root[1].second = sqrt(3) * temp2 / 2.0;
                        root[2].first = -temp1 / 2.0 - B / (3.0 * A);
                        root[2].second = -sqrt(3) * temp2 / 2.0;
                    }
                }

            }
        }
    }

    printf("\nIntersection in: \n");

    //intersection point in the spline definition area
    int flag = 0;
    int num = 1;
    for (int i = 0; i < 3; i++) {
        if (!(root[i].first >= fmin(arr_x1[0], arr_x2[0]) && root[i].first <= fmax(arr_x1[n], arr_x2[k]))) {
            root[i].second = eps;
        }
    }

    //checking that the roots don't match
    if (root[0].first == root[1].first)
        root[0].first = eps;
    if (root[0].first == root[2].first)
        root[0].first = eps;
    if (root[1].first == root[2].first)
        root[1].first = eps;

    //checking the intersection point that it exists

    for (int i = 0; i < 3; i++) {
        if (root[i].second == 0 && root[i].first != eps) {
            double tmp = root[i].first;
            WhereIsX(n + 1, tmp, arr_x1, &m1);
            spline(arr_y1, b1, tmp, arr_x1, c1, d1, m1, &y1);
            printf("%d.  x = %f, y = %f\n", num, tmp, y1);
            flag = 1;
            num++;
        }
    }

    //minimum distance
    if (flag == 0) {
        printf("nowhere :c\n");

        double distnc = MAX;
        for (double i = arr_x1[0]; i <= arr_x1[n]; i += 0.01) {
            for (double j = arr_x2[0]; j <= arr_x2[k]; j += 0.01) {
                WhereIsX(n + 1, i, arr_x1, &m1);
                spline(arr_y1, b1, i, arr_x1, c1, d1, m1, &y1);
                WhereIsX(k + 1, j, arr_x2, &m2);
                spline(arr_y2, b2, j, arr_x2, c2, d2, m2, &y2);
                if (sqrt(pow((i - j), 2) + pow((y1 - y2), 2)) < distnc)
                    distnc = sqrt(pow((i - j), 2) + pow((y1 - y2), 2));
            }
        }
        printf("But we know the distance between the splines: %f\n", distnc);
    }
    return 0;
}
