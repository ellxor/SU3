#pragma once
#define _GNU_SOURCE
#include <complex.h>
#include <math.h>
#include <stdint.h>

/*
 * Minimal pseudo-random number generator (PRNG) library for quantum trajectories.
 */

const uint64_t HighBit = (uint64_t)1 << 63;
const float TwoPi = 2.0f * (float)M_PI;

thread_local uint64_t _random_seed = HighBit;

void set_random_seed(uint64_t seed) {
	_random_seed = seed | HighBit;  // seed must be non zero!
}

// Generate a random 64 bit integer using the xorshift* algorithm (G. Marsaglia)
// Periodicity of (2^64 - 1), never outputs zero.
uint64_t random_u64() {
	_random_seed ^= _random_seed >> 12;
	_random_seed ^= _random_seed << 25;
	_random_seed ^= _random_seed >> 27;
	return _random_seed * 0x2545f4914f6cdd1d;
}

// Generate a uniform random variable from (0,1)
float random_uniform() {
	return (float)random_u64() / (float)UINT64_MAX;
}

