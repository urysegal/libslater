#pragma once
#define L_MAX 6



#define M_MAX L_MAX

#define M_MIN (-L_MAX)

// M goes from -L to L. So 2L+1 possibilities
#define M_COUNT ((L_MAX*2) + 1)



namespace slater {

extern double gaunt_table[L_MAX+1][L_MAX+1][L_MAX+1][M_COUNT][M_COUNT][M_COUNT];

}