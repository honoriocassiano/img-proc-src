/*
 * DanielsonLanczos.h
 *
 *  Created on: Oct 1, 2016
 *      Author: Ronaldo
 */

#ifndef SRC_DANIELSONLANCZOS_H_
#define SRC_DANIELSONLANCZOS_H_

#include <cmath>

using std::sin;

template<unsigned N, typename T=double>
class DanielsonLanczos {
   DanielsonLanczos<N/2,T> next;
public:
   void apply(T* data) {
      next.apply(data);
      next.apply(data+N);

      T wtemp,tempr,tempi,wr,wi,wpr,wpi;
      wtemp = sin(M_PI/N);
      wpr = -2.0*wtemp*wtemp;
      wpi = -sin(2*M_PI/N);
      wr = 1.0;
      wi = 0.0;
      for (unsigned i=0; i<N; i+=2) {
        tempr = data[i+N]*wr - data[i+N+1]*wi;
        tempi = data[i+N]*wi + data[i+N+1]*wr;
        data[i+N] = data[i]-tempr;
        data[i+N+1] = data[i+1]-tempi;
        data[i] += tempr;
        data[i+1] += tempi;

        wtemp = wr;
        wr += wr*wpr - wi*wpi;
        wi += wi*wpr + wtemp*wpi;
      }
   }
};

template<typename T>
class DanielsonLanczos<4,T> {
public:
   void apply(T* data) {
      T tr = data[2];
      T ti = data[3];
      data[2] = data[0]-tr;
      data[3] = data[1]-ti;
      data[0] += tr;
      data[1] += ti;
      tr = data[6];
      ti = data[7];
      data[6] = data[5]-ti;
      data[7] = tr-data[4];
      data[4] += tr;
      data[5] += ti;

      tr = data[4];
      ti = data[5];
      data[4] = data[0]-tr;
      data[5] = data[1]-ti;
      data[0] += tr;
      data[1] += ti;
      tr = data[6];
      ti = data[7];
      data[6] = data[2]-tr;
      data[7] = data[3]-ti;
      data[2] += tr;
      data[3] += ti;
   }
};

template<typename T>
class DanielsonLanczos<2,T> {
public:
   void apply(T* data) {
      T tr = data[2];
      T ti = data[3];
      data[2] = data[0]-tr;
      data[3] = data[1]-ti;
      data[0] += tr;
      data[1] += ti;
   }
};

template<typename T>
class DanielsonLanczos<1,T> {
public:
   void apply(T* data) { }
};

#endif /* SRC_DANIELSONLANCZOS_H_ */
