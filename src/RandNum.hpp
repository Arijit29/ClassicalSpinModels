/** Implementation of Numerical Recipes random number generators **/
#pragma once
#include <cstdlib>

/*--------------- ran0 ---------------*/
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

inline float ran0(long *idum) {
/* “Minimal” random number generator of Park and Miller. Returns a uniform random deviate
between 0.0 and 1.0. Set or reset idum to any integer value (except the unlikely value MASK)
to initialize the sequence; idum must not be altered between calls for successive deviates in
a sequence. */

    long k;
    float ans;

    *idum ^= MASK;              // XORing with MASK allows use of zero and other simple bit patterns for idum.
    k=(*idum)/IQ;
    *idum=IA*(*idum-k*IQ)-IR*k; // Compute idum=(IA*idum) % IM without overflows by Schrage's method.
    if (*idum < 0) *idum += IM;
    ans=AM*(*idum);             // Convert idum to a floating result.
    *idum ^= MASK;              // Unmask before return.
    return ans;
}
/*======================================*/


/*------------------ ran1 --------------*/
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

inline float ran1(long *idum) {
/* “Minimal” random number generator of Park and Miller with Bays-Durham shuffle and added
safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint
values). Call with idum a negative integer to initialize; thereafter, do not alter idum between
successive deviates in a sequence. RNMX should approximate the largest floating value that is
less than 1. */

    int j;
    long k;
    static long iy=0;
    static long iv[NTAB];
    float temp;

    if (*idum <= 0 || !iy) {                // Initialize.
        if (-(*idum) < 1) *idum=1;          // Be sure to prevent idum = 0.
        else *idum = -(*idum);
        for (j=NTAB+7;j>=0;j--) {           // Load the shuffle table (after 8 warm-ups).
            k=(*idum)/IQ;
            *idum=IA*(*idum-k*IQ)-IR*k;
            if (*idum < 0) *idum += IM;
            if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ;                           // Start here when not initializing.
    *idum=IA*(*idum-k*IQ)-IR*k;             // Compute idum=(IA*idum) % IM without overflows by Schrage’s method.
    if (*idum < 0) *idum += IM;
    j=iy/NDIV;                              // Will be in the range 0..NTAB-1.
    iy=iv[j];                               // Output previously stored value and refill the shuffle table.
    iv[j] = *idum;
    if ((temp=AM*iy) > RNMX) return RNMX;   // Because users don’t expect endpoint values.
    else return temp;
}
/*======================================*/

/*------------------ ran3 --------------*/
#define IM1 2147483563
#define IM2 2147483399
#define AM1 (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV1 (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

inline float ran2(long *idum) {
/* Long period (> 2 × 1018 ) random number generator of L’Ecuyer with Bays-Durham shuffle
and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
the endpoint values). Call with idum a negative integer to initialize; thereafter, do not alter
idum between successive deviates in a sequence. RNMX should approximate the largest floating
value that is less than 1. */
    int j;
    long k;
    static long idum2=123456789;
    static long iy=0;
    static long iv[NTAB];
    float temp;
    if (*idum <= 0) {                           // Initialize.
        if (-(*idum) < 1) *idum=1;              // Be sure to prevent idum = 0.
        else *idum = -(*idum);
        idum2=(*idum);
        for (j=NTAB+7;j>=0;j--) {               // Load the shuffle table (after 8 warm-ups).
            k=(*idum)/IQ1;
            *idum=IA1*(*idum-k*IQ1)-k*IR1;
            if (*idum < 0) *idum += IM1;
            if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ1;                              // Start here when not initializing.
    *idum=IA1*(*idum-k*IQ1)-k*IR1;              // Compute idum=(IA1*idum) % IM1 without overflows by Schrage’s method.
    if (*idum < 0) *idum += IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;              // Compute idum2=(IA2*idum) % IM2 likewise.
    if (idum2 < 0) idum2 += IM2;
    j=iy/NDIV1;                                  // Will be in the range 0..NTAB-1.
    iy=iv[j]-idum2;                             // Here idum is shuffled, idum and idum2 are combined to generate output.
    iv[j] = *idum;
    if (iy < 1) iy += IMM1;
    if ((temp=AM1*iy) > RNMX) return RNMX;       // Because users don’t expect endpoint values.
    else return temp;
}
/* ====================================== */

/*----------------- ran3 -----------------*/
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

// According to Knuth, any large MBIG, and any smaller (but still large) MSEED can be substituted for the above values.

inline float ran3(long *idum) {
// Returns a uniform random deviate between 0.0 and 1.0. Set idum to any negative value to initialize or reinitialize the sequence.
    static int inext,inextp;
    static long ma[56];                         // The value 56 (range ma[1..55]) is special and
    static int iff=0;                           // should not be modified; see Knuth.
    long mj,mk;
    int i,ii,k;
    if (*idum < 0 || iff == 0) {                // Initialization.
        iff=1;
        mj=labs(MSEED-labs(*idum));             // Initialize ma[55] using the seed idum and the
        mj %= MBIG;                             // large number MSEED.
        ma[55]=mj;
        mk=1;
        for (i=1;i<=54;i++) {                   // Now initialize the rest of the table, in a slightly random order, with numbers that are not especially random.
            ii=(21*i) % 55;
            ma[ii]=mk;
            mk=mj-mk;
            if (mk < MZ) mk += MBIG;
            mj=ma[ii];
        }
        for (k=1;k<=4;k++) {                    // We randomize them by “warming up the generator.”
            for (i=1;i<=55;i++) {
                ma[i] -= ma[1+(i+30) % 55];
                if (ma[i] < MZ) ma[i] += MBIG;
            }
        }
        inext=0;                            // Prepare indices for our first generated number.
        inextp=31;                          // The constant 31 is special; see Knuth.
        *idum=1;
    }
    // Here is where we start, except on initialization.
    if (++inext == 56) inext=1;                     // Increment inext and inextp, wrapping around 56 to 1.
    if (++inextp == 56) inextp=1;
    mj=ma[inext]-ma[inextp];                        // Generate a new random number subtractively.
    if (mj < MZ) mj += MBIG;                        // Be sure that it is in range.
    ma[inext]=mj;                                   // Store it,  and output the derived uniform deviate.
    return mj*FAC;
}
/* ========================================= */
