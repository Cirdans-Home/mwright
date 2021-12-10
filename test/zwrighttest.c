/* Test for the computation of the Wright function with high precision series
representation. */
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

#define BUFFER_LENGTH 50

/*
This test program contains the test for the computation of the Wright Function
using the series representation and the intervall arithmetic as implementend in
the ARB libray.
*/

/**
 * This test program contains the test for the computation of the Wright
 * Function using the series representation and the intervall arithmetic as
 * implementend in the ARB libray. The input has to be given from stdin as
 * :code:`./zwrighttest λ Re(μ) Im(μ) n prec inputfile`.
 *
 * :param: :math:`\lambda` with :math:`\lambda \in (-1,0)`,
 * :param: :math:`\Re(\mu)` with :math:`\mu \in \mathbb{C}`
 * :param: :math:`\Im(\mu)` with :math:`\mu \in \mathbb{C}`
 * :param: :math:`n` number of terms in the series
 * :param: prec: ARB precision parameter
 * :param: inputfile path of the file containing the input
 *
 * :returns: File wrighttest.out with the computed values
 *
 */
int main(int argc, char *argv[]){

  double x,lambda,remu,immu;
  arb_t arlambda, arremu, arimmu, xread;
  acb_t armu;
  acb_ptr xvec,wvec;
  slong prec;
  int n,numberofx;
  int info = -1;
  char buffer[BUFFER_LENGTH];
  FILE *inputfile, *outputfile;

  if( argc < 7 || argc > 7){
    fprintf(stderr, "./wrighttest λ Re(μ) Im(μ) n prec inputfile\n");
    fprintf(stderr, "\t-1 < λ ≦ 0\n");
    fprintf(stderr, "\tRe(μ) ∈ R\n");
    fprintf(stderr, "\tIm(μ) ∈ R\n");
    fprintf(stderr, "\tprec precision for ARB library");
    fprintf(stderr, "\tPath to the input file formatted as:\n");
    fprintf(stderr, "\t\t 3 ! Number of values followed by x values\n");
    fprintf(stderr, "\t\t 0\n\t\t-1.\n\t\t-2\n");
    return(info);
  }

  // Reading n and precision
  n = atoi(argv[4]);
  prec = atoi(argv[5]);
  // Reading λ
  lambda = atof(argv[1]);
  arb_init(arlambda);
  arb_set_str(arlambda,argv[1],prec);
  // Reading μ
  remu = atof(argv[2]);
  immu = atof(argv[3]);
  arb_init(arremu);
  arb_set_str(arremu,argv[2],prec);
  arb_init(arimmu);
  arb_set_str(arimmu,argv[3],prec);
  acb_init(armu);
  acb_set_arb_arb(armu,arremu,arimmu);

  fprintf(stdout, "Computing Wright Function W_{%1.2f,%1.2f + %1.2fi}(x)\n",lambda,remu,immu);

  fprintf(stdout, "\t Reading data from file %s...\n", argv[5]);
  inputfile = fopen(argv[6],"r");
  if (inputfile == NULL){
    fprintf(stderr, "\t ERROR: File not found\n");
    return(info);
  }

  fscanf(inputfile,"%d",&numberofx);
  fprintf(stdout, "\t Number of x values is %d\n", numberofx);

  xvec = _acb_vec_init(numberofx);
  arb_init(xread);
  for (int i = 0; i < numberofx; i++){
    fgets(buffer, BUFFER_LENGTH, inputfile);
    arb_set_str(xread, buffer, prec);
    acb_set_arb(xvec+i,xread);
  }

  fprintf(stdout, "\t Reading data complete, closing %s.\n", argv[5]);
  fclose(inputfile);

  fprintf(stdout, "\t Computing Wright function values");
  wvec = _acb_vec_init(numberofx);
  for (int i = 0; i < numberofx; i++){
    info = zwrightseries(xvec+i, arlambda, armu, prec, n, wvec+i);
    if( info != 0){
      fprintf(stderr, "ERROR: Error in the %dth Wright computation\n",i);
      return(info);
    }
  }
  fprintf(stdout, "\t Computation complete.\n");

  fprintf(stdout, "\t Writing output to file\n");
  outputfile = fopen("wrighttest.out", "w");
  for (int i = 0; i < numberofx; i++)
  {
      acb_fprintd(outputfile, wvec + i, 16);
      fprintf(outputfile, "\n");    // or any whitespace character
  }
  fclose(outputfile);
  fprintf(stdout, "\t Writing complete\n");

  return(info);
}
