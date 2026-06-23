// CONFIG FILE
#pragma once

enum {
	N = 32,			// number of molecules (algorithm scales as N^4, i.e. doubling N takes 16 times as long)
	ThreadCount = 8,	// number of threads to use (set equal to number of CPU cores available)
	TrajectoryCount = 200,  // number of Monte Carlo simulations to average over
	RungeKuttaPoly = 4,	// Runge-Kutta integration order (can try reducing to see if results change, lower is faster)
	InitialStateLevel = 1,  // what level to put all initial states in (i.e. in fully symmetric sector)
				// > for superposition, this can be changed in main.c using Weyl Tableaux irreps
};

char OutputDirectory[] = "local"; // remember to avoid overwriting previous data!!!
                                  // (and remember to create it before running)

// Model Parameters
// Hamiltonian
const float h = 1;
const float g_plus = 1;
const float g_minus = 1;

// Lindblad channels
const complex float alpha = 1;
const complex float beta = 1;
const complex float alpha_collective = 0;
const complex float beta_collective = 0;

// Integration parameters
float TotalTime = 10; // total time in integration t-space
float JumpTolerance = 0.05f; // decrease this if results don't quite look right (lower is more accurate, but slower)
