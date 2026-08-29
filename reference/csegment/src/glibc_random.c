/*
 * glibc_random.c
 *
 * A byte-exact port of glibc's TYPE_3 (default) random()/srandom().
 *
 * WHY THIS FILE EXISTS
 * --------------------
 * region.c uses  #define flip() (random() & 01)  and never seeds.  The golden
 * fixtures were produced on Linux/glibc, whose unseeded random() runs from
 * state seeded with 1.  Apple's libc random() is a *different* generator, so a
 * stock macOS build cannot reproduce the golden output.  Defining random() and
 * srandom() here makes the program's own definitions win at link time.
 *
 * Algorithm (glibc stdlib/random_r.c, TYPE_3: DEG=31, SEP=3):
 *   r[0] = seed
 *   r[i] = 16807 * r[i-1] mod 2147483647    (Schrage, i = 1..30)
 *   fptr = &r[3], rptr = &r[0]
 *   discard the first 10*31 = 310 outputs
 *   *fptr += *rptr  (uint32 wraparound); result = (uint32)*fptr >> 1
 */

#include <stdint.h>

#define DEG 31
#define SEP 3

static int32_t  gr_state[DEG];
static int      gr_f;           /* index of fptr */
static int      gr_r;           /* index of rptr */
static int      gr_inited = 0;

static long int glibc_random_internal(void);

void
srandom(unsigned int seed)
{
    int32_t  word;
    int      i;
    long int hi, lo;

    /* glibc: "We must make sure the seed is not 0." */
    if (seed == 0)
        seed = 1;

    gr_state[0] = (int32_t) seed;

    word = (int32_t) seed;
    for (i = 1; i < DEG; ++i) {
        /*
         * Compute  word = (16807 * word) % 2147483647  using Schrage's
         * method to avoid overflowing a signed 32-bit intermediate.
         *   hi = word / 127773
         *   lo = word % 127773
         *   word = 16807 * lo - 2836 * hi
         *   if (word < 0) word += 2147483647
         */
        hi = word / 127773;
        lo = word % 127773;
        word = (int32_t) (16807 * lo - 2836 * hi);
        if (word < 0)
            word += 2147483647;
        gr_state[i] = word;
    }

    gr_f = SEP;
    gr_r = 0;
    gr_inited = 1;

    for (i = 0; i < DEG * 10; ++i)
        (void) glibc_random_internal();
}

static long int
glibc_random_internal(void)
{
    uint32_t val;

    /* *fptr += *rptr, wrapping in uint32 */
    val = (uint32_t) gr_state[gr_f] + (uint32_t) gr_state[gr_r];
    gr_state[gr_f] = (int32_t) val;

    if (++gr_f >= DEG)
        gr_f = 0;
    if (++gr_r >= DEG)
        gr_r = 0;

    return (long int) (val >> 1);
}

long int
random(void)
{
    if (!gr_inited)
        srandom(1);             /* glibc's unseeded default state */
    return glibc_random_internal();
}
