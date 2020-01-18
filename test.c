/* test.c - test `cd_vector.h` */

#include "./cd_vector.h"
#include <stdio.h>

int main(int argc, char** argv) {

    v2 x = V2_X, y = V2_Y;
    v2 xy = v2_add(x, y), xyu = v2_unit(v2_add(x, y));


    printf(
        "X: " V2_FMT "\n" 
        "Y: " V2_FMT "\n" 
        "X+Y: " V2_FMT "\n" 
        "unit(X+Y): " V2_FMT "\n"
        , 
        V2__(x), V2__(y),
        V2__(xy), V2__(xyu)
    );
    

    m2x2 A = m2x2_mul(M2X2_I, M2X2_1);

    printf("A: " M2X2_FMT "\n", M2X2__(A));



    return 0;
}

