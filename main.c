#include <stdio.h>

int chebyshev(int i, int j) {
    if(i > j && j >= 2) {
        return 2 * chebyshev(i - 1, j - 1) - chebyshev(i - 2, j);
    }
    else if(i >= 3 && j == i) {
        return 2 * chebyshev(i - 1, j - 1);
    }
    else if(i >= 3 && j == 1) {
        return -1 * chebyshev(i - 2, 1);
    }
    else {
        return 1;
    }
}

double p_coefficients(int k, double g) {
    double p = 0;
    for(int a = 0; a <= k; a++){
        double val = chebyshev(2 * k + 1, 2 * a + 1) * sqrt(2) / M_PI *
            gamma_identity(a) * pow(a + g + 0.5, -a - 0.5) * exp(a + g + 0.5);
        p += val;
    }
    return p;
}


int main() {
    printf("%d\n", chebyshev(5, 3));
    return 0;
}
