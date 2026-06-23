#define _GNU_SOURCE
#define _POSIX_C_SOURCE 200809L
#define NDEBUG

#ifdef __clang__
#pragma clang diagnostic ignored "-Wnon-literal-null-conversion"
#endif

#include <complex.h>
#include <math.h>
#include <pthread.h>
#include <stdalign.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

enum { nullptr = 0 };

#include "clebsh_gordan.c"
#include "parameters.h"
#include "random.h"
#include "utility.h"


typedef complex float WaveVector[Dimension];


void effective_hamiltonian_step(WaveVector dst, WaveVector src, int nu1, int nu2, int nu3, float time_step)
{
	memset(dst, 0, sizeof(WaveVector));

	for (int n1 = nu2; n1 <= nu1; ++n1) {
		for (int n2 = nu3; n2 <= nu2; ++n2) {
			for (int n3 = n2; n3 <= n1; ++n3) {
				int level1_count = n3;
				int level2_count = n1 + n2 - n3;
				int level3_count = N - n1 - n2;

				index_t index = indexof6(n1, n2, n3, nu1, nu2, nu3);
				dst[index] -= src[index] * time_step/2 * (cnormf(alpha) * level1_count + cnormf(beta) * level2_count);

				complex float coeff = src[index] * I * time_step;
				dst[index] -= coeff * h * (level1_count - level3_count);

				for (int i = 0; i < countof(Transition21[index]); ++i) {
					index_t target = Transition21[index][i].target;
					float factor = Transition21[index][i].factor;

					for (int j = 0; j < countof(Transition12[target]); ++j) {
						index_t target2 = Transition12[target][j].target;
						float factor2 = Transition12[target][j].factor;
						dst[target2] += factor2 * factor * coeff * g_plus * g_plus / N;
					}

					for (int j = 0; j < countof(Transition23[target]); ++j) {
						index_t target2 = Transition23[target][j].target;
						float factor2 = Transition23[target][j].factor;
						dst[target2] += factor2 * factor * coeff * g_minus * g_plus / N;
					}
				}

				for (int i = 0; i < countof(Transition32[index]); ++i) {
					index_t target = Transition32[index][i].target;
					float factor = Transition32[index][i].factor;

					for (int j = 0; j < countof(Transition12[target]); ++j) {
						index_t target2 = Transition12[target][j].target;
						float factor2 = Transition12[target][j].factor;
						dst[target2] += factor2 * factor * coeff * g_plus * g_minus / N;
					}

					for (int j = 0; j < countof(Transition23[target]); ++j) {
						index_t target2 = Transition23[target][j].target;
						float factor2 = Transition23[target][j].factor;
						dst[target2] += factor2 * factor * coeff * g_minus * g_minus / N;
					}
				}
			}
		}
	}
}


// Higher order exponetial evolution with Wiener fluctuation for QSD
//
void evolve_under_H_eff(WaveVector wave, int nu1, int nu2, int nu3, float time_step) {
	// In this case of an exponential and linear Hamiltonian, the Runge-Kutta method
	// is identical to a Taylor series expansion, so this is performed directly for efficiency.

	static thread_local WaveVector _a, _b; // Create two temporary vectors using double-buffering technique.
	complex float *a = wave; // Controlled by pointers which are cheap to swap.
	complex float *b = _b;

	int factorial = 1;

	for (int i = 1; i <= RungeKuttaPoly; ++i) {
		effective_hamiltonian_step(b, a, nu1, nu2, nu3, time_step); // |b> = -i Heff dt |a>
		factorial *= i;

		// accumulate Taylor series expansion
		float factor = 1.0f / factorial;

		for (int n1 = nu2; n1 <= nu1; ++n1) {
			for (int n2 = nu3; n2 <= nu2; ++n2) {
				for (int n3 = n2; n3 <= n1; ++n3) {
					index_t index = indexof6(n1, n2, n3, nu1, nu2, nu3);
					wave[index] += factor * b[index];
				}
			}
		}

		if (i == 1) a = _a; // a is temporarily set to wave for first iteration to avoid a copy

		// perform double-buffering: swap pointers
		complex float *tmp = a;
		a = b;
		b = tmp;
	}
}


void collective_lindblad_jump(WaveVector wave, int nu1, int nu2, int nu3)
{
	static thread_local WaveVector copy;
	memcpy(copy, wave, sizeof(WaveVector));
	memset(wave, 0, sizeof(WaveVector));

	for (int n1 = nu2; n1 <= nu1; ++n1) {
		for (int n2 = nu3; n2 <= nu2; ++n2) {
			for (int n3 = n2; n3 <= n1; ++n3) {
				index_t index = indexof6(n1, n2, n3, nu1, nu2, nu3);
				complex float coeff = copy[index];

				for (int i = 0; i < countof(Transition21[index]); ++i) {
					index_t target = Transition21[index][i].target;
					float factor = Transition21[index][i].factor;
					wave[target] += alpha_collective * factor * coeff;
				}

				for (int i = 0; i < countof(Transition32[index]); ++i) {
					index_t target = Transition32[index][i].target;
					float factor = Transition32[index][i].factor;
					wave[target] += beta_collective * factor * coeff;
				}
			}
		}
	}
}


void lindblad_jump(WaveVector wave, size_t choice, int *nu1, int *nu2, int *nu3)
{
	alignas(64) static thread_local WaveVector copy;
	memcpy(copy, wave, sizeof(WaveVector));
	memset(wave, 0, sizeof(WaveVector));

	sector_t next_sector = 0;

	for (int n1 = *nu2; n1 <= *nu1; ++n1) {
		for (int n2 = *nu3; n2 <= *nu2; ++n2) {
			for (int n3 = n2; n3 <= n1; ++n3) {
				index_t index = indexof6(n1, n2, n3, *nu1, *nu2, *nu3);

				struct LindbladTableEntry4 *alpha_jump = &LindbladJump21[index][choice];
				struct LindbladTableEntry2 *beta_jump = &LindbladJump32[index][choice];

				for (int i = 0; i < countof(alpha_jump->entries); ++i)
					wave[alpha_jump->entries[i].target] += copy[index] * alpha * alpha_jump->entries[i].factor;

				for (int i = 0; i < countof(beta_jump->entries); ++i)
					wave[beta_jump->entries[i].target] += copy[index] * beta * beta_jump->entries[i].factor;

				if (next_sector) {
					assert(!alpha_jump->sector || alpha_jump->sector == next_sector);
					assert(!beta_jump->sector || beta_jump->sector == next_sector);
				}

				next_sector |= alpha_jump->sector;
				next_sector |= beta_jump->sector;
			}
		}
	}

	read_sector(next_sector, nu1, nu2, nu3);
	assert(*nu1 + *nu2 + *nu3 == N);

	update_tables(*nu1, *nu2, *nu3);
}


float normalize(WaveVector wave, int nu1, int nu2, int nu3) {
	float norm = 0;

	for (int n1 = nu2; n1 <= nu1; ++n1) {
		for (int n2 = nu3; n2 <= nu2; ++n2) {
			for (int n3 = n2; n3 <= n1; ++n3) {
				index_t index = indexof6(n1, n2, n3, nu1, nu2, nu3);
				norm += cnormf(wave[index]);
			}
		}
	}

	assert(norm);
	float scale = 1.0f / sqrtf(norm);

	for (int n1 = nu2; n1 <= nu1; ++n1) {
		for (int n2 = nu3; n2 <= nu2; ++n2) {
			for (int n3 = n2; n3 <= n1; ++n3) {
				index_t index = indexof6(n1, n2, n3, nu1, nu2, nu3);
				wave[index] *= scale;
			}
		}
	}

	return norm;
}


void run_simulation(int thread_index)
{
	set_random_seed(thread_index);

	struct GT initial = {{ (InitialStateLevel == 1) ? N : 0, (InitialStateLevel == 3) ? 0 : N,0, N,0,0 }}; // fully symmetric sector
	assert(check_valid_swt(initial));

	int nu1 = initial.M[3];
	int nu2 = initial.M[4];
	int nu3 = initial.M[5];

	alignas(64) static thread_local WaveVector wave;
	memset(wave, 0, sizeof(WaveVector));
	wave[indexof(initial)] = 1;
	update_tables(nu1, nu2, nu3);

	float time = 0;
	int jump_hist[Levels*Levels+1] = {0};

	char filepath[100];
	snprintf(filepath, sizeof(filepath), "%s/log-%d.txt", OutputDirectory, thread_index);

	FILE *log = fopen(filepath, "w");

	while (time < TotalTime) {
		float expected_level1 = 0;
		float expected_level2 = 0;
		float expected_level3 = 0;
		float expected_alpha = 0;
		float expected_beta  = 0;

		for (int n1 = nu2; n1 <= nu1; ++n1) {
			for (int n2 = nu3; n2 <= nu2; ++n2) {
				for (int n3 = n2; n3 <= n1; ++n3) {
					index_t index = indexof6(n1, n2, n3, nu1, nu2, nu3);
					assert(!isnan(wave[index]));

					float norm = cnormf(wave[index]);
					expected_level1 += n3 * norm;
					expected_level2 += (n1 + n2 - n3) * norm;

					for (int i = 0; i < countof(Transition21[index]); ++i) {
						float factor = Transition21[index][i].factor;
						expected_alpha += factor*factor * norm;
					}

					for (int i = 0; i < countof(Transition32[index]); ++i) {
						float factor = Transition32[index][i].factor;
						expected_beta += factor*factor * norm;
					}
				}
			}
		}

		expected_level3 = N - expected_level2 - expected_level1;
		if (expected_level3 < 0) expected_level3 = 0; // floating point error correction

		fprintf(log, "%g\t%g\t%g\t%g\n", time, expected_level1, expected_level2, expected_level3);

		float jump_probability = cnormf(alpha) * expected_level1 + cnormf(beta) * expected_level2;
		float jump_probability_collective = cnormf(alpha_collective) * expected_alpha + cnormf(beta_collective) * expected_beta;

		float time_step = 0.05f / fmaxf(fmaxf(jump_probability, jump_probability_collective), 1.0f);

		jump_probability *= time_step;
		jump_probability_collective *= time_step;
		time += time_step;

		float r = random_uniform();

		if (r > jump_probability + jump_probability_collective) {
			evolve_under_H_eff(wave, nu1, nu2, nu3, time_step);
			++jump_hist[Levels*Levels];
		}

		else if (r > jump_probability) {
			collective_lindblad_jump(wave, nu1, nu2, nu3);
		}

		else {
			float sub_jump_table[Levels*Levels] = {0};

			for (int n1 = nu2; n1 <= nu1; ++n1) {
				for (int n2 = nu3; n2 <= nu2; ++n2) {
					for (int n3 = n2; n3 <= n1; ++n3) {
						index_t index = indexof6(n1, n2, n3, nu1, nu2, nu3);

						for (int choice = 0; choice < Levels*Levels; ++choice) {
							struct LindbladTableEntry4 *alpha_jump = &LindbladJump21[index][choice];
							struct LindbladTableEntry2 *beta_jump = &LindbladJump32[index][choice];

							for (int i = 0; i < countof(alpha_jump->entries); ++i)
								sub_jump_table[choice] += cnormf(wave[index] * alpha) * powf(alpha_jump->entries[i].factor, 2);

							for (int i = 0; i < countof(beta_jump->entries); ++i)
								sub_jump_table[choice] += cnormf(wave[index] * beta) * powf(beta_jump->entries[i].factor, 2);
						}
					}
				}
			}

			int choice = 0;
			float accumulator = 0;

			for (; choice < countof(sub_jump_table); ++choice) {
				accumulator += sub_jump_table[choice] * time_step;
				if (accumulator >= r) break;
			}

			assert(choice < countof(sub_jump_table));
			++jump_hist[choice];
			lindblad_jump(wave, choice, &nu1, &nu2, &nu3);
		}

		normalize(wave, nu1, nu2, nu3);
	}

	fclose(log);

	// for (int i = 0; i <= Levels*Levels; ++i) {
	// 	printf("%d: %d\n", i, jump_hist[i]);
	// }
}


atomic(int) thread_pool;
atomic(int) threads_done = 0;
atomic(size_t) millis = 0;

void *run_simulation_thread_wrapper() {
	int next;

	while ((next = atomic_fetch_add(&thread_pool, -1)) > 0) {
		double begin = get_time_from_os();
		run_simulation(next);
		double end = get_time_from_os();

		millis += (size_t)((end - begin) * 1000);
		int complete = atomic_fetch_add(&threads_done, 1) + 1;

		fprintf(stderr, "Trajectory [%d/%d] completed in %.3f seconds.\n", complete, TrajectoryCount, end - begin);
	}

	return nullptr;
}


int main() {
	fprintf(stderr, "Allocating tables: %zu KB.\n", (TotalTableBytes * ThreadCount) >> 10);

	thread_pool = TrajectoryCount;
	threads_done = 0;

	pthread_t threads[ThreadCount];

	for (int i = 0; i < ThreadCount; ++i) pthread_create(&threads[i], nullptr, run_simulation_thread_wrapper, nullptr);
	for (int i = 0; i < ThreadCount; ++i) pthread_join(threads[i], nullptr);

	fprintf(stderr, "Time per trajectory: %zu ms\n", millis / TrajectoryCount);
}
