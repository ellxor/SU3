#define _GNU_SOURCE
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "parameters.h"
#include "utility.h"

#define triangle(n) ((n)*((n)+1)/2)


enum { Dimension = 2048 > (N*N*N/8) ? 2048 : (N*N*N/8) };
enum { Levels = 3 };

typedef const int Tau[Levels];

typedef uint32_t index_t;
typedef uint32_t sector_t;

enum EnergyLevel {
	Level1,
	Level2,
	Level3,
};

enum BoxChange {
	RemoveBox = -1,
	InsertBox = +1,
};

struct GT {
	int M[triangle(Levels)];
};

struct HamiltonianTableEntry {
	index_t target;
	float   factor;
};

struct LindbladTableEntry {
	sector_t sector;
	struct HamiltonianTableEntry entries[];
};

struct LindbladTableEntry2 {
	sector_t sector;
	struct HamiltonianTableEntry entries[2];
};

struct LindbladTableEntry4 {
	sector_t sector;
	struct HamiltonianTableEntry entries[4];
};

struct Combination {
	uint8_t insert_index;
	uint8_t remove_index;
};


sector_t sectorof(struct GT W_nu) {
	int nu1 = W_nu.M[3];
	int nu2 = W_nu.M[4];

	static_assert(N < UINT16_MAX, "number of molecules does not fit into 16 bits");
	return nu1 | nu2 << 16; // nu3 = N - nu1 - nu2
}

void read_sector(sector_t sector, int *nu1, int *nu2, int *nu3) {
	*nu1 =  sector        & UINT16_MAX;
	*nu2 = (sector >> 16) & UINT16_MAX;
	*nu3 = N - *nu1 - *nu2;
}

index_t indexof6(int n1, int n2, int n3, int nu1, int nu2, int nu3) {
	index_t index = (long)(n1 - nu2)
	              + (long)(n2 - nu3) * (long)(nu1 - nu2 + 1)
	              + (long)(n3 - n2)  * (long)(nu1 - nu2 + 1) * (long)(nu2 - nu3 + 1)
		      + 1;

	assert(index < Dimension);
	return index;
}

index_t indexof(struct GT W_nu) {
	return indexof6(W_nu.M[1], W_nu.M[2], W_nu.M[0], W_nu.M[3], W_nu.M[4], W_nu.M[5]);
}


thread_local struct HamiltonianTableEntry Transition12[Dimension][4];
thread_local struct HamiltonianTableEntry Transition21[Dimension][4];
thread_local struct HamiltonianTableEntry Transition23[Dimension][2];
thread_local struct HamiltonianTableEntry Transition32[Dimension][2];
thread_local struct HamiltonianTableEntry Transition31[Dimension][2];

thread_local struct LindbladTableEntry4 LindbladJump21[Dimension][Levels*Levels];
thread_local struct LindbladTableEntry2 LindbladJump32[Dimension][Levels*Levels];

const size_t TotalTableBytes
	= sizeof(Transition12) + sizeof(Transition21)
	+ sizeof(Transition23) + sizeof(Transition32)
	+ sizeof(LindbladJump21) + sizeof(LindbladJump32);

Tau SymmetryChanges[Levels][6] = {
	[Level1] = {{1,1,1}, {1,2,1}, {1,1,2}, {1,2,2}, {1,1,3}, {1,2,3}},
	[Level2] = {{0,1,1}, {0,2,1}, {0,1,2}, {0,2,2}, {0,1,3}, {0,2,3}},
	[Level3] = {{0,0,1}, {0,0,2}, {0,0,3}},
};

int p(struct GT W_mu, int j, int k) {
	int index = triangle(k - 1) + j - 1;
	return W_mu.M[index] + k - j;
}

double A(struct GT W_mu, int l, Tau tau) {
	int sign = (tau[l-2] >= tau[l-1]) ? 1 : -1;
	double numerator, denominator, factor = 1;

	for (int k = 1; k <= l; ++k) {
		if (k != tau[l-1]) {
			numerator = p(W_mu, tau[l-2], l-1) - p(W_mu, k, l) + 1;
			denominator = p(W_mu, tau[l-1], l) - p(W_mu, k, l);
			assert(denominator);
			factor *= numerator / denominator;
		}
	}

	for (int k = 1; k <= l - 1; ++k) {
		if (k != tau[l-2]) {
			numerator = p(W_mu, tau[l-1], l) - p(W_mu, k, l-1);
			denominator = p(W_mu, tau[l-2],l-1) - p(W_mu, k, l-1) + 1;
			assert(denominator);
			factor *= numerator / denominator;
		}
	}

	assert(factor >= 0);
	return sign * sqrt(factor);
}

double zeta(struct GT W_mu, enum EnergyLevel i, Tau tau) {
	int numerator = 1;
	int denominator = 1;

	for (int k = 1; k <= i; ++k) {
		numerator *= p(W_mu, tau[i], i+1) - p(W_mu, k, i);
	}

	for (int k = 1; k <= i + 1; ++k) {
		if (k != tau[i]) {
			denominator *= p(W_mu, tau[i], i+1) - p(W_mu, k, i+1);
		}
	}

	assert(numerator * denominator > 0);
	double factor = sqrt((double)numerator / (double)denominator);

	for (int l = i+2; l <= Levels; ++l) {
		factor *= A(W_mu, l, tau);
	}

	return factor;
}

double log_f_factor(struct GT W_nu) {
	long nu1 = W_nu.M[3];
	long nu2 = W_nu.M[4];
	long nu3 = W_nu.M[5];

	double log_numerator = log((nu1 - nu2 + 1)*(nu2 - nu3 + 1)*(nu1 - nu3 + 2));
	double log_denominator = lgamma(nu1 + 3) + lgamma(nu2 + 2) + lgamma(nu3 + 1);

	return log_numerator - log_denominator;
}

double ratio(struct GT W_nu, struct GT W_mu) {
	return exp(log_f_factor(W_mu) - log_f_factor(W_nu));
}

struct GT compute_delta(struct GT W_nu, Tau tau, enum BoxChange change) {
	struct GT W_mu = W_nu;

	for (int k = 0; k < Levels; ++k) {
		if (!tau[k]) continue;

		int index = triangle(k) + tau[k] - 1;
		W_mu.M[index] += change;
	}

	return W_mu;
}

bool check_valid_swt(struct GT W_nu) {
	int n1 = W_nu.M[1];
	int n2 = W_nu.M[2];
	int n3 = W_nu.M[0];
	int nu1 = W_nu.M[3];
	int nu2 = W_nu.M[4];
	int nu3 = W_nu.M[5];

	return (nu1 >= n1) && (n1 >= nu2) && (nu2 >= n2)
	    && (n1 >= n3) && (n3 >= n2) && (n2 >= nu3) && (nu3 >= 0);
}


void hamiltonian_clebsh_gordan(struct GT W_nu, enum EnergyLevel a, enum EnergyLevel b, double R[Levels],
                               const struct Combination combs[], struct HamiltonianTableEntry table[], size_t CombinationCount)
{
	size_t count = CombinationCount / Levels;
	size_t mod = count / Levels;

	for (int k = 0; k < count; ++k) {
		const int *tau_a = SymmetryChanges[a][combs[k].insert_index];
		const int *tau_b = SymmetryChanges[b][combs[k].remove_index];

		struct GT W_mu = compute_delta(W_nu, tau_b, RemoveBox);
		struct GT W_la = compute_delta(W_mu, tau_a, InsertBox);

		if (!check_valid_swt(W_mu)) continue;
		if (!check_valid_swt(W_la)) continue;

		table[k % mod].target  = indexof(W_la);
		table[k % mod].factor += R[tau_b[Levels - 1] - 1] * zeta(W_mu, a, tau_a) * zeta(W_mu, b, tau_b);
	}
}


void lindblad_clebsh_gordan(struct GT W_nu, enum EnergyLevel a, enum EnergyLevel b, double R[Levels],
	                    const struct Combination combs[], struct LindbladTableEntry *table, size_t CombinationCount)
{
	size_t mod = CombinationCount / (Levels * Levels);

	for (int k = 0; k < CombinationCount; ++k) {
		const int *tau_a = SymmetryChanges[a][combs[k].insert_index];
		const int *tau_b = SymmetryChanges[b][combs[k].remove_index];

		struct GT W_mu = compute_delta(W_nu, tau_b, RemoveBox);
		struct GT W_la = compute_delta(W_mu, tau_a, InsertBox);

		if (!check_valid_swt(W_mu)) continue;
		if (!check_valid_swt(W_la)) continue;

		size_t i = k / mod, j = k % mod;

		size_t entry_size = (mod == 4) ? sizeof(struct LindbladTableEntry4) : sizeof(struct LindbladTableEntry2);
		struct LindbladTableEntry *entry = (struct LindbladTableEntry *)((char *)table + i * entry_size);

		entry->sector |= sectorof(W_la);
		entry->entries[j].target = indexof(W_la);
		entry->entries[j].factor = sqrt(R[tau_b[Levels - 1] - 1]) * zeta(W_mu, a, tau_a) * zeta(W_mu, b, tau_b);
	}
}


void update_tables(int nu1, int nu2, int nu3) {
	static const struct Combination GroundToGroundCombinations[] = {
		{0,0}, {0,1}, {1,0}, {1,1}, {2,2}, {2,3}, {3,2}, {3,3}, {4,4}, {4,5}, {5,4}, {5,5},
		{0,2}, {0,3}, {1,2}, {1,3}, {0,4}, {0,5}, {1,4}, {1,5}, 
		{2,0}, {2,1}, {3,0}, {3,1}, {2,4}, {2,5}, {3,4}, {3,5}, 
		{4,0}, {4,1}, {5,0}, {5,1}, {4,2}, {4,3}, {5,2}, {5,3}, 
	};

	static const struct Combination GroundToExcitedCombinations[] = {
		{0,0}, {0,1}, {1,2}, {1,3}, {2,4}, {2,5},
		{0,2}, {0,3}, {0,4}, {0,5}, {1,0}, {1,1},
		{1,4}, {1,5}, {2,0}, {2,1}, {2,2}, {2,3},
	};

	static const struct Combination ExcitedToGroundCombinations[] = {
		{0,0}, {1,0}, {2,1}, {3,1}, {4,2}, {5,2},
		{2,0}, {3,0}, {4,0}, {5,0}, {0,1}, {1,1},
		{4,1}, {5,1}, {0,2}, {1,2}, {2,2}, {3,2},
	};

	static_assert(countof(GroundToGroundCombinations) == 6 * 6, "size mismatch");
	static_assert(countof(GroundToExcitedCombinations) == 6 * 3, "size mismatch");
	static_assert(countof(ExcitedToGroundCombinations) == 3 * 6, "size mismatch");

	memset(LindbladJump21, 0, sizeof(LindbladJump21));
	memset(LindbladJump32, 0, sizeof(LindbladJump32));
	memset(Transition12, 0, sizeof(Transition12));
	memset(Transition21, 0, sizeof(Transition21));
	memset(Transition23, 0, sizeof(Transition23));
	memset(Transition32, 0, sizeof(Transition32));

	struct GT W_nu = {{ [3] = nu1, [4] = nu2, [5] = nu3 }};
	double ratios[Levels] = {0};

	for (int k = 0; k < Levels; ++k) {
		struct GT W_mu = W_nu;
		W_mu.M[triangle(Levels - 1) + k] -= 1;
		ratios[k] = ratio(W_nu, W_mu);
	}

	for (int n1 = nu2; n1 <= nu1; ++n1) {
		for (int n2 = nu3; n2 <= nu2; ++n2) {
			for (int n3 = n2; n3 <= n1; ++n3) {
				W_nu.M[0] = n3;
				W_nu.M[1] = n1;
				W_nu.M[2] = n2;

				index_t index = indexof6(n1, n2, n3, nu1, nu2, nu3);

				//                              |to>    <from|
				hamiltonian_clebsh_gordan(W_nu, Level2, Level1, ratios, GroundToGroundCombinations,  Transition21[index], countof(GroundToGroundCombinations) );
				hamiltonian_clebsh_gordan(W_nu, Level1, Level2, ratios, GroundToGroundCombinations,  Transition12[index], countof(GroundToGroundCombinations) );
				hamiltonian_clebsh_gordan(W_nu, Level3, Level2, ratios, GroundToExcitedCombinations, Transition32[index], countof(GroundToExcitedCombinations));
				hamiltonian_clebsh_gordan(W_nu, Level2, Level3, ratios, ExcitedToGroundCombinations, Transition23[index], countof(ExcitedToGroundCombinations));
				hamiltonian_clebsh_gordan(W_nu, Level3, Level1, ratios, GroundToExcitedCombinations, Transition32[index], countof(GroundToExcitedCombinations));

				//                           |to>    <from|
				lindblad_clebsh_gordan(W_nu, Level2, Level1, ratios, GroundToGroundCombinations, (struct LindbladTableEntry *)LindbladJump21[index], countof(GroundToGroundCombinations));
				lindblad_clebsh_gordan(W_nu, Level3, Level2, ratios, GroundToExcitedCombinations, (struct LindbladTableEntry *)LindbladJump32[index], countof(GroundToExcitedCombinations));
			}
		}
	}
}
