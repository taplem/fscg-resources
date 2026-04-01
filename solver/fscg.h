// definitions shared by C and Lua code.

typedef struct FSCGFile FSCGFile;

enum {
	FSCG_MODE_SERIAL_PRIMAL   = 'p',
	FSCG_MODE_SERIAL_GREEDY   = 'g',
	FSCG_MODE_PARALLEL_PRIMAL = 'P',
	FSCG_MODE_PARALLEL_DUAL   = 'D'
};

typedef struct FSCGModel {
	int I;          // number of units
	int J;          // number of schedules
	int J0;         // number of schedules, including pruned
	int K;          // number of items
	int T;          // number of transports
	int C;          // number of general constraints
	int CX;         // number of nonzeros in objective+constraint matrix
	int Ca;         // start of artificial constraints (artificial z-variable definitions)
	int Z;          // number of z-variables (including artificial)
	int Za;         // start of artificial z-variables
	int *Kt;        // number of transports for each item    (K entries)
	uint32_t *Ij;   // unit -> start of schedules            (I+1 entries)
	double *ITu;    // (unit, transport) -> utility          (I*T entries)
	uint32_t *Jx;   // schedule -> start of nonzeros         (J+1 entries)
	double *Jo;     // schedule -> objective                 (J entries)
	uint8_t *Xk;    // nonzero -> item                       (X entries)
	double *Xv;     // nonzero -> value                      (X entries)
	int *Cnnz;      // constraint -> number of nonzeros      (C+1 entries)
	char *Cs;       // constraint -> sense                   (C+1 entries)
	double *Cr;     // constraint -> rhs                     (C+1 entries)
	int *CXc;       // constraint nonzero -> column          (CX entries)
	double *CXv;    // constraint nonzero -> value           (CX entries)
	char *err;      // error message
	// config:
	double scale;   // model scaling
	double epsilon; // error tolerance
	double gamma;   // new schedule acceptance threshold relative to best schedule
	int threads;    // how many threads can fscg use?
	int s_W;        // solver work heap size
	int s_N0;       // solver additional schedule pool size (master problem)
	int s_N1;       // solver additional schedule pool size (subproblem)
	char mode;      // see FSCG_MODE_* constants
} FSCGModel;
