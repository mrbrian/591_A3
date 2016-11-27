
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <ctime>
#include <vector>
#include <string>
#include <unordered_map>
#include <random>
#include <cstdint>
#include <algorithm>



// -----------------------------------------
//    RANDOM NUMBER GENERATION
// -----------------------------------------

// Mersenne Twister 19937 generator
// A Mersenne Twister pseudo-random generator of 32-bit numbers with a state size of 19937 bits.
// http://www.cplusplus.com/reference/random/mt19937/

std::mt19937 MersenneTwisterPRNG;

// Uniform real distribution
// Random number distribution that produces floating-point values according to a uniform distribution, 
// which is described by the following probability density function (...) - more details here:
// http://www.cplusplus.com/reference/random/uniform_real_distribution/

std::uniform_real_distribution<double> m_URD;

// Hemisphere sampling (i.e., for diffuse reflection) 
// function "m_Vector HemisphereSampling(m_Vector m_normal)" below
// calls both m_RND_1 and m_RND_2
#define m_RND_1 (2.0*m_URD(MersenneTwisterPRNG)-1.0)
#define m_RND_2 (m_URD(MersenneTwisterPRNG))

int main()
{
for (int i =0 ; i < 100; i++)
	printf("%f\n", m_RND_2);
	return 0;
}
