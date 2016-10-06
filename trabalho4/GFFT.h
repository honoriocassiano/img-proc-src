/*
 * GFFT.h
 *
 *  Created on: Oct 1, 2016
 *      Author: Ronaldo
 */

#ifndef SRC_GFFT_H_
#define SRC_GFFT_H_

#include "DanielsonLanczos.h"

template<unsigned P,
         typename T=double>
class GFFT {
   enum { N = 1<<P };
   DanielsonLanczos<N,T> recursion;
public:
   void fft(T* data) {
      scramble(data,N);
      recursion.apply(data);
   }

   void scramble(T* data, unsigned n) {
	   unsigned j=1;
	   unsigned m;

	   for (unsigned i = 1; i < n; i += 2) {
		   if (j > i) {
			   std::swap(data[j-1], data[i-1]);
			   std::swap(data[j], data[i]);
		   }
		   m = n;
		   while (m>=2 && j>m) {
			   j -= m;
			   m >>= 1;
		   }
		   j += m;
	   };
   }
};

#endif /* SRC_GFFT_H_ */
