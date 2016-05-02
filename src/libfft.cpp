#include <math.h>

#ifndef PI
#define PI 3.14159265358979323844
#endif

void four1(double data[], int nn, int sign)
/*
   FOUR1 replaces data by its discrete Fourier transform,
   if isigen is input as 1; or replaces data by nn times its
   inverse discrete Fourier transform, if sign is input as -1.
   Data is a complex array of length nn or equivalently, a real
   array of length 2*nn. nn must be an integer power of 2.
*/
{
    register double wr, wi, wpr, wpi;
    double wtemp, theta, tempr, tempi;
    register int i, j;
    int n, m, max, step;

    n = 2*nn;
    j = 0;
    for (i = 0; i < n; i += 2) {
        if (j>i) {
            tempr = data[j];
            tempi = data[j+1];
            data[j] = data[i];
            data[j+1] = data[i+1];
            data[i] = tempr;
            data[i+1] = tempi;
        }
        m = n/2;
        while (m >= 2 && j >= m) {
            j -= m;
            m /= 2;
        }
        j += m;
    }

    max = 2;
    while (n>max) {
        step = 2*max;
        theta = 2.0*PI/(sign*max);
        wpr = -2.0*sin(0.5*theta)*sin(0.5*theta);
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        for (m = 0; m < max; m += 2) {
            for (i = m; i < n; i += step) {
                j = i+max;
                tempr = wr*data[j]-wi*data[j+1];
                tempi = wr*data[j+1]+wi*data[j];
                data[j] = data[i]-tempr;
                data[j+1] = data[i+1]-tempi;
                data[i] = data[i]+tempr;
                data[i+1] = data[i+1]+tempi;
            }
            wtemp = wr;
            wr += wr*wpr-wi*wpi;
            wi += wi*wpr+wtemp*wpi;
        }
        max = step;
    }
}

void twofft(double data1[], double data2[], double fft1[], double fft2[], int n)
/*
   Given two real input arrays data1 and data2, each of length n, this routine
   calls four1 and returns two complex output arrays, fft1 and fft2, each of
   complex length n (i.e. real length 2*n), which contain the discrete Fourier
   transforms of the respective data's. n must be an integer power of 2.
*/
{
    int nn2, j, jj;
    double rep, rem, aip, aim;

    for (j = 0; j < n; j++) {
        jj = 2*j;
        fft1[jj] = data1[j];
        fft1[jj+1] = data2[j];
    }
    four1(fft1, n, 1);
    fft2[0] = fft1[1];
    fft1[1] = fft2[1] = 0.0;
    nn2 = 2*n;
    for (j = 1; j <= n/2; j++) {
        jj = 2*j;
        rep = 0.5*(fft1[jj]+fft1[nn2-jj]);
        rem = 0.5*(fft1[jj]-fft1[nn2-jj]);
        aip = 0.5*(fft1[jj+1]+fft1[nn2-jj+1]);
        aim = 0.5*(fft1[jj+1]-fft1[nn2-jj+1]);
        fft1[jj] = rep;
        fft1[jj+1] = aim;
        fft1[nn2-jj] = rep;
        fft1[nn2-jj+1] = -aim;
        fft2[jj] = aip;
        fft2[jj+1] = -rem;
        fft2[nn2-jj] = aip;
        fft2[nn2-jj+1] = rem;
    }
}

void realft(double data[], int n, int sign)
/*
   Calculates the Fourier transform of a set of 2n real-valued
   data points. Replaces this data (array data) by the positive
   frequency half of its complex Fourier transform. The real-valued
   first and last components of the complex transform are returned
   as elements data[0] and data[1] respectively. n must be a power
   of 2. This routine also calculates the inverse transform of a
   complex data array if it is the transform of real data. The
   result in this case must be multiplied by 1/n.
*/
{
    double theta, c1, c2, wr, wi, wpr, wpi, h1r, h1i, h2r, h2i, wtemp;
    int i, i1, i2, i3, i4;

    theta = PI/n;
    c1 = 0.5;
    if (sign == 1) {
        c2 = -0.5;
        four1 (data,n,sign);
    } else {
        c2 = 0.5;
        theta = -theta;
    }

    wpr = -2.0*sin(0.5*theta)*sin(0.5*theta);
    wpi = sin(theta);
    wr = 1.0+wpr;
    wi = wpi;
    for (i = 1; i < n/2+1; i++) {
        i1 = i+i; i2 = i1+1; i3 = n+n-i1; i4 = i3+1;
        h1r =  c1*(data[i1]+data[i3]);
        h1i =  c1*(data[i2]-data[i4]);
        h2r = -c2*(data[i2]+data[i4]);
        h2i =  c2*(data[i1]-data[i3]);
        data[i1] =  h1r+wr*h2r-wi*h2i;
        data[i2] =  h1i+wr*h2i+wi*h2r;
        data[i3] =  h1r-wr*h2r+wi*h2i;
        data[i4] = -h1i+wr*h2i+wi*h2r;
        wtemp = wr;
        wr += wr*wpr-wi*wpi;
        wi += wi*wpr+wtemp*wpi;
    }

    if (sign == -1) {
        h1r = data[0];
        data[0] = c1*(h1r+data[1]);
        data[1] = c1*(h1r-data[1]);
        four1 (data,n,sign);
    } else {
        h1r = data[0];
        data[0] = h1r+data[1];
        data[1] = h1r-data[1];
    }
}

void cosft(double y[], int n, int sign)
/*
   Calculates the cosine transform of a set y of n real-valued data points.
   The transformed data replace the original data in array y. n must be a
   power of 2. Set sign to +1 for a transform, and to -1 for an inverse
   transform. For an inverse transform, the output array should be multiplied
   by 2/n.
*/
{
    double wr, wi, wpr, wpi, wtemp, theta;
    double y1, y2, sum, enfo, sumo, sume, even, odd;
    int m, j;

    theta = PI/n;
    wr = 1.0; wi = 0.0;
    wpr = -2.0*pow(sin(0.5*theta),2.0);  wpi = sin(theta);
    sum = y[0];
    m = n/2;
    for (j = 1; j < m; j++) {
        wtemp = wr;
        wr += wr*wpr-wi*wpi;  wi += wi*wpr+wtemp*wpi;
        y1 = 0.5*(y[j]+y[n-j]);  y2 = (y[j]-y[n-j]);
        y[j] = y1-wi*y2;  y[n-j] = y1+wi*y2;
        sum += wr*y2;
    }
    realft (y,m,1);
    y[1] = sum;
    for (j = 3; j < n; j += 2) {
        sum += y[j];
        y[j] = sum;
    }
    if (sign == -1) {
        even = y[0];  odd = y[1];
        for (j = 2; j < n; j+=2) {
            even += y[j];  odd += y[j+1];
        }
        enfo = 2.0*(even-odd);
        sumo = y[0]-enfo;  sume = 2.0*odd/n-sumo;
        y[0] = 0.5*enfo;  y[1] -= sume;
        for (j = 2; j < n; j += 2) {
            y[j] -= sumo;  y[j+1] -= sume;
        }
    }
}

void sinft(double y[], int n)
/*
   Calculates the sine transform of a set y of n real-valued data points.
   The transformed data replace the original data in array y. n must be a
   power of 2. For an inverse transform, the output array should be multiplied
   by 2/n.
*/
{
    double wr, wi, wpr, wpi, wtemp, theta;
    double y1, y2, sum;
    int m, j;

    theta = PI/n;
    wr = 1.0; wi = 0.0;
    wpr = -2.0*pow(sin(0.5*theta),2.0);  wpi = sin(theta);
    y[0] = 0.0;
    m = n/2;
    for (j = 1; j < m; j++) {
        wtemp = wr;
        wr += wr*wpr-wi*wpi;  wi += wi*wpr+wtemp*wpi;
        y1 = wi*(y[j]+y[n-j]);  y2 = 0.5*(y[j]-y[n-j]);
        y[j] = y1+y2;  y[n-j] = y1-y2;
    }
    realft (y,m,1);
    sum = 0.0;
    y[0] = 0.5*y[1];
    y[1] = 0.0;
    for (j = 0; j < n; j += 2) {
        sum += y[j];
        y[j] = y[j+1];
        y[j+1] = sum;
    }
}
