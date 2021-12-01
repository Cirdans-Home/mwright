/* BSD 3-Clause License

 Copyright (c) 2021, Fabio Durastante, Lidia Aceto
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this
    list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

 3. Neither the name of the copyright holder nor the names of its
    contributors may be used to endorse or promote products derived from
    this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "include/wright_config.h"
#include <stdio.h>
#include <math.h>

void printversion(){
  printf("Wright Library Version %d.%d\n",wright_VERSION_MAJOR,wright_VERSION_MINOR);
  flint_printf("Working with arb-%s\n", arb_version);
}

/**
 * WRIGHT computes the Wright function using the series formulation and the ARB
 * library with interval arithmetic in high-precision. This version is used to
 * compare the results obtained with the inversion of the Laplace transform for
 * some of the cases for which we do not have a closed form expression.
 *
 * :param: x The `x` in the argument :math:`-|x|t^\lambda`
 * :param: t: The `t` in the argument :math:`-|x|t^\lambda` ( :math:`t > 0` )
 * :param: lambda First parameter of the Wright function
 * :param: prec: Precision parameter for the ARB library
 * :param: n Number of terms of the series to be used
 * :param: w Variable containing the output of the computation
 * :returns: info 0 for correct run, 1 otherwise
 */
int wrightseries(arb_t x, arb_t lambda, arb_t mu, slong prec, slong n, arb_t w){
  int info;
  arb_t c,arbi,arbifact;
  char stringi[2*n];

  info = 0;

  arb_poly_t truncatedseries;
  arb_poly_init(truncatedseries);

  arb_init(c);
  arb_init(arbi);
  arb_init(arbifact);


  for(int i = 0; i <= n; i++){
    sprintf(stringi,"%d",i);
    arb_set_str(arbi,stringi,prec);
    arb_fac_ui(arbifact,i,prec);     // i!
#ifdef DEBUG
    fprintf(stdout, "i = %d : i! = %s\n",i,arb_get_str(arbifact,5,0));
#endif
    arb_mul(arbi,lambda,arbi,prec);  // i <- λ i
#ifdef DEBUG
    fprintf(stdout, "i = %d : λ i = %s\n",i,arb_get_str(arbi,5,0));
#endif
    arb_add(c,arbi,mu,prec);         // c <- i + \mu = λ i + μ
#ifdef DEBUG
    fprintf(stdout, "i = %d : λ i + μ = %s\n",i,arb_get_str(c,5,0));
#endif
    arb_gamma(c,c,2*prec);             // c <- Γ(c) = Γ(λ i + μ)
#ifdef DEBUG
    fprintf(stdout, "i = %d : Γ(λ i + μ) = %s\n",i,arb_get_str(c,5,0));
#endif
    if( arb_is_finite(c) ){
      arb_mul(c,arbifact,c,2*prec);      // c <- i! c = i! Γ(λ i + μ)
#ifdef DEBUG
      fprintf(stdout, "i = %d : i! Γ(λ i + μ) = %s\n",i,arb_get_str(c,5,0));
#endif
      arb_inv(c,c,2*prec);               // c <- 1/c = 1 / (i! Γ(λ i + μ))
#ifdef DEBUG
      fprintf(stdout, "i = %d : 1 / (i! Γ(λ i + μ)) = %s\n",i,arb_get_str(c,5,0));
#endif
    }else{
      arb_zero(c);
    }
#ifdef DEBUG
    fprintf(stdout, "i = %d : 1 / (i! Γ(λ i + μ)) = %s\n",i,arb_get_str(c,5,0));
#endif
    arb_poly_set_coeff_arb(truncatedseries,i,c);
  }

#ifdef DEBUG
  arb_poly_printd(truncatedseries,6);
#endif

  arb_poly_evaluate(w,truncatedseries,x,prec);

  return(info);
}
