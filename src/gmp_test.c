#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>

mpfr_ptr gmp_ret_f(void) {
    mpfr_ptr x = (mpfr_ptr) malloc(sizeof(__mpfr_struct));
    mpfr_init2(x, 64);
    mpfr_set_d(x, 0.0, MPFR_RNDN);
    mpfr_add_d(x, x, 5.0, MPFR_RNDN);
    mpfr_add(x, x, x, MPFR_RNDN);
    return x;
}
 
void gmp_test(void)
{
    mpfr_ptr k = gmp_ret_f();
    mpfr_printf ("%.6Rf", k);
    
    
    
    mpz_t x;
    mpz_t y;
    mpz_t result;

    mpz_init(x);
    mpz_init(y);
    mpz_init(result);

    mpz_set_str(x, "7612058254738945", 10);
    mpz_set_str(y, "9263591128439081", 10);

    mpz_mul(result, x, y);
    gmp_printf("\n    %Zd\n*\n    %Zd\n--------------------\n%Zd\n\n", x, y, result);

    /* free used memory */
    mpz_clear(x);
    mpz_clear(y);
    mpz_clear(result);
}

