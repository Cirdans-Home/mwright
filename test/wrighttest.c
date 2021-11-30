/* Test for the computation of the Wrigt Function with the Series */
// BSD 3-Clause License
//
// Copyright (c) 2021, Fabio Durastante, Lidia Aceto
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#include "include/wright_config.h"
#include <stdio.h>

/*
This test program contains the test for the computation of the Wright Function
using the series representation and the intervall arithmetic as implementend in
the ARB libray.
*/
int main(){

  double x,lambda,mu;
  arb_t arx, arlambda, armu, w;
  char input[100];
  slong prec = 128;
  int n;

  arb_init(arx);
  arb_init(arlambda);
  arb_init(armu);
  arb_init(w);

  x = -1.0;
  lambda = -1.0/3.0;
  mu = 2.0/3.0;
  n = 4;

  sprintf(input,"%f",x);
  arb_set_str(arx,input,prec);

  sprintf(input,"%f",lambda);
  arb_set_str(arlambda,input,prec);

  sprintf(input,"%f",mu);
  arb_set_str(armu,input,prec);

  wright(arx, arlambda, armu, prec, n, w);

  fprintf(stdout, "\nW_{%f,%f}(%f) = %s\n",lambda,mu,x,arb_get_str(w,prec,0));

  return 0;
}
