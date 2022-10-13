#include "./gaunt-table.h"
#include "stdlib.h"
#include "stdio.h"


int
main(int argc, const char *argv[])
{
    int l1 = atoi(argv[1]);
    int l2 = atoi(argv[2]);
    int l3 = atoi(argv[3]);
    int m1 = atoi(argv[4])+ L_MAX;
    int m2 = atoi(argv[5]) + L_MAX;
    int m3 = atoi(argv[6]) + L_MAX;

    fprintf(stdout,"<%d %d | %d %d | %d %d > = %f", l1, m1-L_MAX, l2, m2-L_MAX, l3, m3-L_MAX, slater::gaunt_table[l1][l2][l3][m1][m2][m3]);

    return 0;

}
