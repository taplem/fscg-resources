#define _GNU_SOURCE // for madvise, qsort_q, vasprintf

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <pthread.h>
#include <x86intrin.h>

#ifdef FSCG_GUROBI
#include <gurobi_c.h>
#endif
#ifdef FSCG_HIGHS
#include <highs_c_api.h>
#endif

#ifdef FSCG_TRACE
#define trace(fmt, ...) fprintf(stderr, fmt "\n", ##__VA_ARGS__)
#else
#define trace(...)
#endif

#define DEFAULT_SUBEPSILON 1e-2
#define DEFAULT_EPSILON    1e-6
#define DEFAULT_GAMMA      0.20

#define LIKELY(x)    __builtin_expect(!!(x), 1)
#define UNLIKELY(x)  __builtin_expect(!!(x), 0)

// LP row numbering:
//   0..T-1     transport row t (implicit item k):   Tt = sum_ij w_ij*r_ijk
//   T..T+I-1   area constraint for unit i:          sum_ij w_ij = 1
//
// LP column numbering:
//   0..T-1     t-variables
//   T..T+Z-1   z-variables
//    (...)     slacks phase 1
//   last S     schedule pool

#define trow(t)     (t)

#include "fscg.h"

typedef struct LP {
#ifdef FSCG_GUROBI
	GRBmodel *model;
	int *j;
	int *bc, *br;
	int sizej;
	int sizebc, sizebr;
	int nrow, ncol;
	char cmd;
#endif
#ifdef FSCG_HIGHS
	void *model;
	int *i, *j;
	int *bc, *br;
	double *sol;
	double *x, *y, *z, *a, *t;
	int nrow, ncol;
	int nxyzi, sizexyzi;
	int naj, sizeaj;
	int sizet;
	int sizebc, sizebr;
	int errcount;
	char cmd;
#endif
} LP;

// an unit can be in exactly one of the following two states:
// (1) it has one uninstantiated schedule and no area constraint;
// (2) it has at least two instantiated schedules and an area constraint.
typedef struct Unit {
	int i;          // global unit
	int ns;         // number of schedules
	int sx;         // xor of schedules
	int r;          // area constraint or -1 if uninstantiated
	uint32_t it;    // last iteration this unit changed
} Unit;

typedef union Fp32Bits { float f32; uint32_t u32; } Fp32Bits;
#define b2fp(x) (((Fp32Bits){.u32=(x)}).f32)
#define fp2b(x) (((Fp32Bits){.f32=(x)}).u32)
#define sem2fp(s,e,m) b2fp(((s)<<31)|((e)<<23)|(m))

typedef uint64_t WSchd;

#define WSCHD_K          (0xffffffff00000000ull)
#define wschd_k(w)       ((w)>>32)
#define wschd_fpk(fp)    (~fp2b((float)(fp)))
#define wschd_knew(k)    ((uint64_t)(k)<<32)
#define wschd_ss(w)      ((uint16_t)((w)>>16))
#define wschd_ssnew(ss)  ((uint64_t)(ss)<<16)
#define wschd_s(w)       ((uint16_t)(w))
#define wschd_snew(s)    ((uint64_t)(s))
#define wschd_new(k,ss,s) (wschd_knew(k) | wschd_ssnew(ss) | wschd_snew(s))

// active schedule
typedef struct ASchd {
	int s;          // pool schedule number (active) or -1 (pending deletion)
	uint32_t it;    // last iteration this schedule was basic
} ASchd;

enum {
	SCLOCK_INIT,
	SCLOCK_SCAN,
	SCLOCK_MASTER,
	SCLOCK__MAX
};

#define SOLVER_CLOCK CLOCK_MONOTONIC
// #define SOLVER_CLOCK CLOCK_PROCESS_CPUTIME_ID

typedef struct SchPool {
	int N;          // size of this schedule pool
	int n;          // number of schedules in use
	int nextfree;   // freelist head
	int *l;         // schedule -> local unit                    (N entries)
	size_t *j;      // schedule -> global schedule number        (N entries)
	int *Kt;        // (schedule, item) -> transport             (N entries)
	int *free;      // freelist, if in use                       (N entries or NULL)
} SchPool;

typedef struct ScanState {
	struct Solver *solver;
	SchPool pool;   // schedule pool
	uint16_t no;    // scan state number
	int l;          // first local unit managed by this scan state
	int L;          // number of local units managed by this scan state
	double wd;      // score of the best schedule in the work queue
	double lag;     // lagrangian bound
	WSchd *Ws;      // work buffer                               (solver->W entries)
	int *t;         // best transport work array                 (model->K entries)
	double *d;      // best item duals work array                (model->K entries)
} ScanState;

enum {
	SS_PHASE1,
	SS_PHASE2,
	SS_EXIT
};

typedef struct Solver {
	LP lp;
	SchPool pool;
	FSCGModel *model;
	FILE *log;
	double epsilon; // optimality tolerance
	double gamma;   // scan accept tolerance
	double scale;   // area constraint right hand size
	double obj;     // current objective value
	double lag;     // current lagrangian bound
	double ub;      // best known objective upper bound
	double maxgap;  // max optimality gap
	bool serial;    // true if LP should not use parallelism
	uint8_t state;  // solver state
	uint8_t no;     // subproblem number if this is a subproblem solver
	uint8_t nsub;   // total number of subsolvers
	uint32_t iter;  // slave of calls to lp solver
	int nss;        // number of scanning threads
	int L;          // number of local units
	int W;          // max number of new schedules
	int Tmaskbytes; // ((model->T+63)>>6)<<3
	int R0;         // start of area constraints (model->T + model->C)
	int Rmax;       // maximum number of area constraints = S-L
	int A0;         // start of instantiated schedules (lp column)
	int r;          // number of instantiated units (area constraints)
	int a;          // number of instantiated schedules
	int r1;         // number of new units to instantiate
	int a1;         // number of schedules queued for instantiation
	ASchd *As;      // active list                               (S entries)
	int *as1;       // schedule instantiation buffer             (S entries)
	Unit *unit;     // units                                     (L entries)
	int *rl;        // area constraint -> unit                   (Rmax entries)
	int *rl1;       // unit instantiation buffer                 (Rmax entries)
	double *Tr;     // transport constraints rhs                 (model->T entries)
	uint64_t *Trx;  // transport constraint rhs changed?         (model->T bits)
	double *dual;   // current dual                              (model->T+model->C+Rmax entries)
	int *tmpi;      // integer work array                        (sizetmp entries)
	double *tmpx;   // double work array                         (sizetmp entries)
	uint64_t *tmpq; // u64 work array                            (sizetmp entries)
	ScanState *ss;  // scan thread states                        (nss entries)
	pthread_t *st;  // scan threads, if using threaded scanning  (nss entries)
	pthread_barrier_t bar; // scan thread barrier
	uint64_t *mask; // subsolver finish mask
	double clock[SCLOCK__MAX]; // clocks
	struct timespec tp;
	struct timespec tp0;
} Solver;

/* ---- System -------------------------------------------------------------- */

static int fscg_cpus(void)
{
	cpu_set_t set;
	CPU_ZERO(&set);
	if (sched_getaffinity(0, sizeof(set), &set))
		return 0;
	return CPU_COUNT(&set);
}

static double fscg_time(void)
{
	struct timespec tp;
	clock_gettime(SOLVER_CLOCK, &tp);
	return tp.tv_sec + tp.tv_nsec*1e-9;
}

/* ---- Model management ---------------------------------------------------- */

static void fscg_model_destroy(FSCGModel *model)
{
	free(model->Kt);
	free(model->Ij);
	free(model->ITu);
	free(model->Jx);
	free(model->Jo);
	// free(model->Xk);
	// free(model->Xv);
	munmap(model->Xk, (size_t)model->J0*model->K*sizeof(*model->Xk));
	munmap(model->Xv, (size_t)model->J0*model->K*sizeof(*model->Xv));
	free(model->Cnnz);
	free(model->Cs);
	free(model->Cr);
	free(model->CXc);
	free(model->CXv);
	free(model->err);
}

static bool model_err(FSCGModel *model, const char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	free(model->err);
	vasprintf(&model->err, fmt, ap);
	va_end(ap);
	return false;
}

/* ---- LP (Gurobi) --------------------------------------------------------- */

#define LP_BARRIER 0x0
#define LP_SIMPLEX 0x1
#define LP_SERIAL  0x2
#define LP_WARM    0x4

#ifdef FSCG_GUROBI

#define LP_BASIC  0
#define LP_LOWER (-1)

__attribute__((cold, noinline))
static void lp_assertfail(LP *lp, int err)
{
	if (lp)
		fprintf(stderr, "Gurobi call failed: %s (%d)\n", GRBgetmerrormsg(lp->model), err);
	else
		fprintf(stderr, "Gurobi call failed (%d)\n", err);
	exit(1);
}

__attribute__((always_inline))
static inline void lp_assertstatus(LP *lp, int err)
{
	if (UNLIKELY(err))
		lp_assertfail(lp, err);
}

static void lp_create(LP *lp, int flags)
{
	memset(lp, 0, sizeof(*lp));
	GRBenv *env;
	lp_assertstatus(NULL, GRBloadenv(&env, NULL));
	lp_assertstatus(NULL, GRBsetintparam(env, "Method", (flags & LP_SIMPLEX) ? 1 : 2));
	lp_assertstatus(NULL, GRBsetintparam(env, "OutputFlag", 0));
	if (flags & LP_SERIAL)
		lp_assertstatus(NULL, GRBsetintparam(env, "Threads", 1));
	if (flags & LP_WARM)
		lp_assertstatus(NULL, GRBsetintparam(env, "LpWarmStart", 2));
	lp_assertstatus(NULL, GRBnewmodel(env, &lp->model, NULL, 0, NULL, NULL, NULL, NULL, NULL));
	lp_assertstatus(lp, GRBsetintattr(lp->model, "ModelSense", -1)); // max
}

static void lp_destroy(LP *lp)
{
	lp_assertstatus(NULL, GRBfreemodel(lp->model));
	free(lp->j);
	free(lp->bc);
	free(lp->br);
}

static void lp_cmd(LP *lp, char s)
{
	if (s != lp->cmd) {
		lp_assertstatus(lp, GRBupdatemodel(lp->model));
		lp->cmd = s;
	}
}

static bool lp_growbufsz(int *sz, int n)
{
	int s = *sz;
	if (n <= s)
		return false;
	if (!s)
		 s = 32;
	while (s < n)
		s <<= 1;
	*sz = s;
	return true;
}

static void lp_growj(LP *lp, int n)
{
	if (lp_growbufsz(&lp->sizej, n))
		lp->j = realloc(lp->j, lp->sizej*sizeof(*lp->j));
}

static void lp_growbasis(LP *lp, int ncol, int nrow)
{
	if (lp_growbufsz(&lp->sizebc, ncol)) {
		lp->bc = realloc(lp->bc, lp->sizebc*sizeof(*lp->bc));
		// lp->sol = realloc(lp->sol, lp->sizebc*sizeof(*lp->sol));
		memset(lp->bc+ncol, 0, (lp->sizebc-ncol)*sizeof(*lp->bc));
	}
	if (lp_growbufsz(&lp->sizebr, nrow)) {
		lp->br = realloc(lp->br, lp->sizebr*sizeof(*lp->br));
		memset(lp->br+nrow, 0, (lp->sizebr-nrow)*sizeof(*lp->br));
	}
}

static void lp_add_col(LP *lp, int nnz, int *nzi, double *nzx, double lb, double ub, double obj)
{
	// fprintf(stderr, "lp_add_col: %d [%g %g %g]\n", nnz, lb, ub, obj);
	lp_cmd(lp, 'c');
	lp_assertstatus(lp, GRBaddvar(lp->model, nnz, nzi, nzx, obj, lb, ub, 'C', NULL));
	if (lp->bc) {
		lp_growbasis(lp, lp->ncol+1, lp->nrow);
		lp->bc[lp->ncol] = LP_LOWER;
	}
	lp->ncol++;
}

static void qsorti(int *xs, size_t n)
{
	for (;;) {
		if (n < 32) {
			// insertion sort for small arrays
			for (size_t i=1; i<n; i++) {
				int x = xs[i];
				size_t j = i;
				for (; j; j--) {
					if (xs[j-1] <= x)
						break;
					xs[j] = xs[j-1];
				}
				xs[j] = x;
			}
			return;
		}
		// pivot = pseudo-median of xs[0], xs[n/2], xs[n-1]
		int a = xs[0];
		int b = xs[n>>1];
		int c = xs[n-1];
		int p = a < b ? (b < c ? b : (a < c ? c : a)) : (a < c ? a : (b < c ? c : b));
		size_t i = 0;
		size_t j = n-1;
		// partition left <= p <= right
		for (;;) {
			while (xs[i] <= p && i < j) i++;
			while (xs[j] >= p && i < j) j--;
			if (i == j) break;
			int x = xs[i];
			xs[i] = xs[j];
			xs[j] = x;
		}
		int *ys;
		if (i < (n>>1)) {
			// left interval is smaller
			ys = xs;
			xs += i;
			n -= i;
		} else {
			// right interval is smaller
			ys = xs+i;
			j = n-i;
			n = i;
		}
		qsorti(ys, j);
	}
}

static void lp_delidx(int *xs, int *del, int n, int ndel)
{
	if (!ndel)
		return;
	int p = del[0];
	for (int i=p+1, j=1; i<n; i++) {
		if (j < ndel && i == del[j])
			j++;
		else
			xs[p++] = xs[i];
	}
	assert(p == n-ndel);
}

static void lp_del_cols(LP *lp, int ndel, int *is)
{
	lp_cmd(lp, 0);
	qsorti(is, ndel);
	lp_assertstatus(lp, GRBdelvars(lp->model, ndel, is));
	lp_assertstatus(lp, GRBupdatemodel(lp->model));
	if (lp->bc)
		lp_delidx(lp->bc, is, lp->ncol, ndel);
	lp->ncol -= ndel;
}

static void lp_add_row(LP *lp, int nnz, int *nzi, double *nzx, char sense, double rhs)
{
	// fprintf(stderr, "lp_add_row: %d [%c %g]\n", nnz, sense, rhs);
	lp_cmd(lp, 'r');
	lp_assertstatus(lp, GRBaddconstr(lp->model, nnz, nzi, nzx, sense, rhs, NULL));
	if (lp->br) {
		lp_growbasis(lp, lp->ncol, lp->nrow+1);
		lp->br[lp->nrow] = LP_LOWER;
	}
	lp->nrow++;
}

static void lp_del_rows(LP *lp, int ndel, int *is)
{
	lp_cmd(lp, 0);
	qsorti(is, ndel);
	lp_assertstatus(lp, GRBdelconstrs(lp->model, ndel, is));
	lp_assertstatus(lp, GRBupdatemodel(lp->model));
	if (lp->br)
		lp_delidx(lp->br, is, lp->nrow, ndel);
	lp->nrow -= ndel;
}

static void lp_change_rhs(LP *lp, int nnz, int *nzi, double *nzx)
{
	lp_assertstatus(lp, GRBsetdblattrlist(lp->model, "RHS", nnz, nzi, nzx));
}

static void lp_set_obj(LP *lp, int nnz, int *nzi, double *nzx)
{
	lp_assertstatus(lp, GRBsetdblattrlist(lp->model, "Obj", nnz, nzi, nzx));
}

static double lp_get_obj(LP *lp)
{
	double obj;
	lp_assertstatus(lp, GRBgetdblattr(lp->model, "ObjVal", &obj));
	return obj;
}

static void lp_get_primal(LP *lp, int start, int num, double *x)
{
	lp_assertstatus(lp, GRBgetdblattrarray(lp->model, "X", start, num, x));
}

static void lp_get_dual(LP *lp, int start, int num, double *y)
{
	lp_assertstatus(lp, GRBgetdblattrarray(lp->model, "Pi", start, num, y));
}

static int lp_get_basis(LP *lp, int *is)
{
	lp_growj(lp, lp->ncol);
	lp_assertstatus(lp, GRBgetintattrarray(lp->model, "VBasis", 0, lp->ncol, lp->j));
	size_t n = 0;
	for (size_t i=0, c=lp->ncol; i<c; i++) {
		if (lp->j[i] == 0)
			is[n++] = i;
	}
	return n;
}

static bool lp_solve(LP *lp)
{
	lp_assertstatus(lp, GRBupdatemodel(lp->model));
	lp_assertstatus(lp, GRBoptimize(lp->model));
	int status;
	lp_assertstatus(lp, GRBgetintattr(lp->model, "Status", &status));
	return status == 2;
}

#endif

/* ---- LP (HiGHS) ---------------------------------------------------------- */

#if FSCG_HIGHS

#define LP_BASIC kHighsBasisStatusBasic
#define LP_LOWER kHighsBasisStatusLower

__attribute__((cold, noinline))
static void lp_assertfail(void)
{
	fputs("HiGHS API call failed", stderr);
	exit(1);
}

__attribute__((always_inline))
static inline void lp_assertstatus(int err)
{
	if (UNLIKELY(err < 0))
		lp_assertfail();
}

static void lp_create(LP *lp, int flags)
{
	memset(lp, 0, sizeof(*lp));
	lp->model = Highs_create();
	lp_assertstatus(Highs_setOptionValue(lp->model, "solver", flags & LP_SIMPLEX ? "simplex" : "ipm"));
	lp_assertstatus(Highs_setOptionValue(lp->model, "output_flag", "false"));
	lp_assertstatus(Highs_changeObjectiveSense(lp->model, kHighsObjSenseMaximize));
	lp->errcount = 3;
}

static void lp_destroy(LP *lp)
{
	Highs_destroy(lp->model);
	free(lp->i);
	free(lp->j);
	free(lp->x);
	free(lp->y);
	free(lp->z);
	free(lp->a);
	free(lp->t);
}

static bool lp_growbufsz(int *sz, int n)
{
	int s = *sz;
	if (n <= s)
		return false;
	if (!s)
		 s = 32;
	while (s < n)
		s <<= 1;
	*sz = s;
	return true;
}

static void lp_growbasis(LP *lp, int ncol, int nrow)
{
	if (lp_growbufsz(&lp->sizebc, ncol)) {
		lp->bc = realloc(lp->bc, lp->sizebc*sizeof(*lp->bc));
		lp->sol = realloc(lp->sol, lp->sizebc*sizeof(*lp->sol));
		memset(lp->bc+ncol, 0, (lp->sizebc-ncol)*sizeof(*lp->bc));
	}
	if (lp_growbufsz(&lp->sizebr, nrow)) {
		lp->br = realloc(lp->br, lp->sizebr*sizeof(*lp->br));
		memset(lp->br+nrow, 0, (lp->sizebr-nrow)*sizeof(*lp->br));
	}
}

static void lp_cmd(LP *lp, char s)
{
	if (s == lp->cmd)
		return;
	switch (lp->cmd) {
		case 'c':
			lp_assertstatus(Highs_addCols(
					lp->model,
					lp->nxyzi,    // nnew
					lp->x,        // obj [nnew]
					lp->y,        // lb  [nnew]
					lp->z,        // ub  [nnew]
					lp->naj,      // nnz
					lp->i,        // nzi [nnew]
					lp->j,        // nzr [nnz]
					lp->a         // nzx [nnz]
			));
			lp->ncol += lp->nxyzi;
			if (lp->bc)
				lp_growbasis(lp, lp->ncol, lp->nrow);
			break;
		case 'r':
			lp_assertstatus(Highs_addRows(
				lp->model,
				lp->nxyzi,      // nnew
				lp->x,          // lb  [nnew]
				lp->y,          // ub  [nnew]
				lp->naj,        // nnz
				lp->i,          // nzi [nnew]
				lp->j,          // nzc [nnz]
				lp->a           // nzx [nnz]
			));
			lp->nrow += lp->nxyzi;
			if (lp->br) {
				lp_growbasis(lp, lp->ncol, lp->nrow);
				for (int i=lp->nrow-lp->nxyzi; i<lp->nrow; i++)
					lp->br[i] = kHighsBasisStatusLower;
			}
			break;
	}
	lp->nxyzi = lp->naj = 0;
	lp->cmd = s;
}

static void lp_growxyzi(LP *lp, int n)
{
	if (lp_growbufsz(&lp->sizexyzi, lp->nxyzi+=n)) {
		lp->x = realloc(lp->x, lp->sizexyzi*sizeof(*lp->x));
		lp->y = realloc(lp->y, lp->sizexyzi*sizeof(*lp->y));
		lp->z = realloc(lp->z, lp->sizexyzi*sizeof(*lp->z));
		lp->i = realloc(lp->i, lp->sizexyzi*sizeof(*lp->i));
	}
}

static void lp_growaj(LP *lp, int n)
{
	if (lp_growbufsz(&lp->sizeaj, lp->naj+=n)) {
		lp->a = realloc(lp->a, lp->sizeaj*sizeof(*lp->a));
		lp->j = realloc(lp->j, lp->sizeaj*sizeof(*lp->j));
	}
}

static void lp_add(LP *lp, char cmd, int nnz, int *nzi, double *nzx, double x, double y, double z)
{
	lp_cmd(lp, cmd);
	int p = lp->nxyzi, nzp = lp->naj;
	lp_growxyzi(lp, 1);
	lp_growaj(lp, nnz);
	lp->x[p] = x;
	lp->y[p] = y;
	lp->z[p] = z;
	lp->i[p] = nzp;
	memcpy(lp->j+nzp, nzi, nnz*sizeof(*nzi));
	memcpy(lp->a+nzp, nzx, nnz*sizeof(*nzx));
}

static void lp_add_col(LP *lp, int nnz, int *nzi, double *nzx, double lb, double ub, double obj)
{
	lp_add(lp, 'c', nnz, nzi, nzx, obj, lb, ub);
}

static void qsorti(int *xs, size_t n)
{
	for (;;) {
		if (n < 32) {
			// insertion sort for small arrays
			for (size_t i=1; i<n; i++) {
				int x = xs[i];
				size_t j = i;
				for (; j; j--) {
					if (xs[j-1] <= x)
						break;
					xs[j] = xs[j-1];
				}
				xs[j] = x;
			}
			return;
		}
		// pivot = pseudo-median of xs[0], xs[n/2], xs[n-1]
		int a = xs[0];
		int b = xs[n>>1];
		int c = xs[n-1];
		int p = a < b ? (b < c ? b : (a < c ? c : a)) : (a < c ? a : (b < c ? c : b));
		size_t i = 0;
		size_t j = n-1;
		// partition left <= p <= right
		for (;;) {
			while (xs[i] <= p && i < j) i++;
			while (xs[j] >= p && i < j) j--;
			if (i == j) break;
			int x = xs[i];
			xs[i] = xs[j];
			xs[j] = x;
		}
		int *ys;
		if (i < (n>>1)) {
			// left interval is smaller
			ys = xs;
			xs += i;
			n -= i;
		} else {
			// right interval is smaller
			ys = xs+i;
			j = n-i;
			n = i;
		}
		qsorti(ys, j);
	}
}

static void lp_delidx(int *xs, int *del, int n, int ndel)
{
	if (!ndel)
		return;
	int p = del[0];
	for (int i=p+1, j=1; i<n; i++) {
		if (j < ndel && i == del[j])
			j++;
		else
			xs[p++] = xs[i];
	}
	assert(p == n-ndel);
}

static void lp_delidx2(double *xs, int *del, int n, int ndel)
{
	if (!ndel)
		return;
	int p = del[0];
	for (int i=p+1, j=1; i<n; i++) {
		if (j < ndel && i == del[j])
			j++;
		else
			xs[p++] = xs[i];
	}
	assert(p == n-ndel);
}

static void lp_del_cols(LP *lp, int ndel, int *is)
{
	lp_cmd(lp, 0);
	qsorti(is, ndel);
	lp_assertstatus(Highs_deleteColsBySet(lp->model, ndel, is));
	if (lp->bc) {
		lp_delidx(lp->bc, is, lp->ncol, ndel);
		lp_delidx2(lp->sol, is, lp->ncol, ndel);
	}
	lp->ncol -= ndel;
}

static void lp_del_rows(LP *lp, int ndel, int *is)
{
	lp_cmd(lp, 0);
	qsorti(is, ndel);
	lp_assertstatus(Highs_deleteRowsBySet(lp->model, ndel, is));
	if (lp->br)
		lp_delidx(lp->br, is, lp->nrow, ndel);
	lp->nrow -= ndel;
}

static void lp_change_rhs(LP *lp, int nnz, int *nzi, double *nzx)
{
	// hack, we only ever change the rhs of equality constraints
	lp_cmd(lp, 0);
	int nrow;
	int nrownnz;
	lp_assertstatus(Highs_getRowsBySet(
		lp->model,
		nnz,        // number of rows to get
		nzi,        // row indices
		&nrow,      // number of rows returned
		lp->x,      // lb [nrow]
		lp->y,      // ub [nrow]
		&nrownnz,   // number of nonzero elements
		NULL,       // nzi
		NULL,       // nzc
		NULL        // nzx
	));
	for (int r=0; r<nnz; r++) {
		if (lp->x[r] == -1/0.0) {
			// '<' row
			lp->y[r] = nzx[r];
		} else if (lp->y[r] == 1/0.0) {
			// '>' row
			lp->x[r] = nzx[r];
		} else {
			// '=' row
			lp->x[r] = nzx[r];
			lp->y[r] = nzx[r];
		}
	}
	lp_assertstatus(Highs_changeRowsBoundsBySet(lp->model, nnz, nzi, lp->x, lp->y));
}

// static int lp_get_col(LP *lp, int col, int nnz, int *nzi, double *nzx, double *obj)
// {
// 	int ncol, start;
// 	lp_assertstatus(
// 		Highs_getColsByRange(lp->model, col, col, &ncol, obj, NULL, NULL, &nnz, &start, nzi, nzx));
// 	return nnz;
// }

static void lp_add_row(LP *lp, int nnz, int *nzi, double *nzx, char sense, double rhs)
{
	double lb, ub;
	switch (sense) {
		case '<':
			lb = -1/0.0;
			ub = rhs;
			break;
		case '>':
			lb = rhs;
			ub = 1/0.0;
			break;
		case '=':
			lb = ub = rhs;
			break;
		default:
			lp_assertfail();
			return;
	}
	lp_add(lp, 'r', nnz, nzi, nzx, lb, ub, /* unused: */ lb);
}

static void lp_set_obj(LP *lp, int nnz, int *nzi, double *nzx)
{
	lp_cmd(lp, 0);
	lp_assertstatus(Highs_changeColsCostBySet(lp->model, nnz, nzi, nzx));
}

static double lp_get_obj(LP *lp)
{
	double obj;
	lp_assertstatus(Highs_getDoubleInfoValue(lp->model, "objective_function_value", &obj));
	return obj;
}

static void lp_getsolution(LP *lp, bool primal, int start, int num, double *x)
{
	if (lp_growbufsz(&lp->sizet, primal ? lp->ncol : lp->nrow))
		lp->t = realloc(lp->t, lp->sizet*sizeof(*lp->t));
	if (primal) {
		lp_assertstatus(Highs_getSolution(lp->model, lp->t, NULL, NULL, NULL));
	} else {
		lp_assertstatus(Highs_getSolution(lp->model, NULL, NULL, NULL, lp->t));
	}
	memcpy(x, lp->t+start, num*sizeof(*x));
}

static void lp_get_primal(LP *lp, int start, int num, double *x)
{
	lp_getsolution(lp, true, start, num, x);
}

static void lp_get_dual(LP *lp, int start, int num, double *y)
{
	lp_getsolution(lp, false, start, num, y);
}

static int lp_get_basis(LP *lp, int *is)
{
	lp_growaj(lp, lp->ncol+lp->nrow);
	lp_assertstatus(Highs_getBasis(lp->model, lp->j, lp->j+lp->ncol));
	size_t n = 0;
	for (size_t i=0, c=lp->ncol; i<c; i++) {
		if (lp->j[i] == 1)
			is[n++] = i;
	}
	return n;
}

static bool lp_solve(LP *lp)
{
	lp_cmd(lp, 0);
	int ok = Highs_run(lp->model);
	int status = Highs_getModelStatus(lp->model);
	if (ok < 0 && status == kHighsModelStatusSolveError) {
		trace("[warn] HiGHS solve error, switching to simplex");
		lp_assertstatus(Highs_setHighsOptionValue(lp->model, "solver", "simplex"));
		status = kHighsModelStatusNotset;
	}
	if (ok < 0 && (status == kHighsModelStatusNotset || status == kHighsModelStatusUnknown)) {
		trace("[warn] HiGHS solve failed, re-solving model from scratch");
		Highs_clearSolver(lp->model);
		ok = Highs_run(lp->model);
		status = Highs_getModelStatus(lp->model);
	}
	if (status == kHighsModelStatusUnboundedOrInfeasible || status == kHighsModelStatusUnbounded)
		return false;
	lp_assertstatus(ok);
	return true;
}

#endif

/* ---- Clocks -------------------------------------------------------------- */

static void solver_clock_start(Solver *solver)
{
	clock_gettime(SOLVER_CLOCK, &solver->tp);
}

static void solver_clock_stop(Solver *solver, int clock)
{
	struct timespec tp;
	clock_gettime(SOLVER_CLOCK, &tp);
	double t = (tp.tv_sec-solver->tp.tv_sec) + (tp.tv_nsec-solver->tp.tv_nsec)*1e-9;
	solver->clock[clock] += t;
}

/* ---- Schedule pool management -------------------------------------------- */

static void pool_reset(SchPool *pool)
{
	pool->n = 0;
	if (pool->free) {
		for (int n=0; n<pool->N-1; n++)
			pool->free[n] = n+1;
		pool->free[pool->N-1] = -1;
	}
}

static void pool_create(SchPool *pool, FSCGModel *model, int N, bool fl)
{
	pool->N = N;
	pool->l = malloc(N*sizeof(*pool->l));
	pool->j = malloc(N*sizeof(*pool->j));
	pool->Kt = malloc(model->K*N*sizeof(*pool->Kt));
	pool->free = fl ? malloc(N*sizeof(*pool->free)) : NULL;
	pool->nextfree = 0;
	pool_reset(pool);
}

static void pool_destroy(SchPool *pool)
{
	free(pool->l);
	free(pool->j);
	free(pool->Kt);
	free(pool->free);
}

static int pool_alloc(SchPool *pool)
{
	int n = pool->n;
	if (n == pool->N)
		return -1;
	pool->n++;
	if (pool->free) {
		int s = pool->nextfree;
		assert(s >= 0);
		pool->nextfree = pool->free[s];
		pool->free[s] = -1;
		return s;
	} else {
		return n;
	}
}

static void pool_free(SchPool *pool, int s)
{
	if (pool->free) {
		assert(pool->free[s] == -1);
		pool->free[s] = pool->nextfree;
		pool->nextfree = s;
		pool->n--;
	}
}

static size_t pool_transfer(FSCGModel *model, SchPool *dst, SchPool *src, size_t ss, size_t lofs)
{
	size_t ds = pool_alloc(dst);
	assert(ds < (size_t)dst->N);
	assert(ss < (size_t)src->N);
	dst->l[ds] = src->l[ss] + lofs;
	dst->j[ds] = src->j[ss];
	memcpy(dst->Kt+model->K*ds, src->Kt+model->K*ss, model->K*sizeof(*dst->Kt));
	return ds;
}

static bool pool_cmp(FSCGModel *model, SchPool *p1, size_t s1, SchPool *p2, size_t s2)
{
	size_t j = p1->j[s1];
	if (j != p2->j[s2])
		return false;
	size_t K = model->K;
	int *t1 = p1->Kt + s1*K;
	int *t2 = p2->Kt + s2*K;
	for (size_t x=model->Jx[j], X=model->Jx[j+1]; x<X; x++) {
		size_t k = model->Xk[x];
		if (t1[k] != t2[k])
			return false;
	}
	return true;
}

/* ---- Sync ---------------------------------------------------------------- */

// queue an active schedule for deinstantiation and return the schedule number to the schedule pool
static void solver_deinst(Solver *solver, size_t a)
{
	assert(a < (size_t)solver->a);
	int s = solver->As[a].s;
	assert(s >= 0);
	pool_free(&solver->pool, s);
	solver->As[a].s = -1;
	Unit *uni = &solver->unit[solver->pool.l[s]];
	uni->ns--;
	uni->sx ^= s;
	// it's not possible for uni->ns to ever fall to zero, because only nonbasic schedules can be
	// deinstantiated, and at least one schedule will always be basic.
	assert(uni->ns >= 1);
}

static void sync_modify_rhs(Solver *solver, size_t l, double sign)
{
	assert(solver->unit[l].ns == 1);
	FSCGModel *model = solver->model;
	double scale = sign * solver->scale;
	size_t s = solver->unit[l].sx;
	size_t j = solver->pool.j[s];
	for (size_t x=model->Jx[j], X=model->Jx[j+1]; x<X; x++) {
		size_t t = solver->pool.Kt[s*model->K+model->Xk[x]];
		solver->Trx[t>>6] |= 1ULL << (t&0x3f);
		solver->Tr[t] += scale*model->Xv[x];
	}
}

// queue a schedule number for instantiation
static void solver_inst(Solver *solver, size_t s)
{
	size_t l = solver->pool.l[s];
	Unit *uni = &solver->unit[l];
	uni->it = solver->iter;
	solver->as1[solver->a1++] = s;
	if (uni->r < 0) {
		// unit switches from uninstantiated to instantiated. queue it for instantiation.
		solver->as1[solver->a1++] = ~uni->sx; // basic
		solver->rl1[solver->r1++] = l;
		sync_modify_rhs(solver, l, 1);
	}
	uni->sx ^= s;
	uni->ns++;
}

static void sync_deinst_schd(Solver *solver)
{
	int *del = solver->tmpi;
	size_t ndel = 0;
	size_t ap = 0;
	for (size_t a=0, A=solver->a; a<A; a++) {
		ASchd as = solver->As[a];
		if (as.s < 0 || solver->unit[solver->pool.l[as.s]].ns == 1) {
			// pending deletion or unit will be deleted
			del[ndel++] = solver->A0 + a;
		} else {
			solver->As[ap++] = as;
		}
	}
	solver->a = ap;
	lp_del_cols(&solver->lp, ndel, del);
}

static void sync_deinst_units(Solver *solver)
{
	int *del = solver->tmpi;
	size_t ndel = 0;
	size_t rp = 0;
	for (size_t r=0, R=solver->r; r<R; r++) {
		size_t l = solver->rl[r];
		assert(solver->unit[l].ns >= 1);
		if (solver->unit[l].ns == 1) {
			sync_modify_rhs(solver, l, -1);
			solver->unit[l].r = -1;
			del[ndel++] = solver->R0 + r;
		} else {
			solver->unit[l].r = rp;
			solver->rl[rp++] = l;
		}
	}
	solver->r = rp;
	lp_del_rows(&solver->lp, ndel, del);
}

static void sync_inst_units(Solver *solver)
{
	size_t rp = solver->r;
	for (size_t r=0, R=solver->r1; r<R; r++) {
		size_t l = solver->rl1[r];
		solver->unit[l].r = rp;
		solver->rl[rp++] = l;
		lp_add_row(&solver->lp, 0, NULL, NULL, '=', solver->scale);
	}
	solver->r = rp;
	solver->r1 = 0;
	assert(solver->r <= solver->Rmax);
}

static void sync_update_rhs(Solver *solver)
{
	FSCGModel *model = solver->model;
	size_t nnz = 0;
	int *nzi = solver->tmpi;
	double *nzx = solver->tmpx;
	// for (size_t t=0, tx=0, T=model->T; t<T; t+=64, tx++) {
	// 	for (uint64_t m=solver->Trx[tx]; m; m&=m-1) {
	// 		size_t idx = t + __builtin_ctzll(m);
	// 		nzi[nnz] = trow(idx);
	// 		nzx[nnz] = solver->Tr[idx];
	// 		nnz++;
	// 	}
	// }
	for (size_t t=0; t<(size_t)model->T; t++) {
		nzi[nnz] = trow(t);
		nzx[nnz] = solver->Tr[t];
		nnz++;
	}
	lp_change_rhs(&solver->lp, nnz, nzi, nzx);
	memset(solver->Trx, 0, solver->Tmaskbytes);
}

// rhs_t = -scale * ∑_{i ∈ uninst} r_{i,j_i,k->t}
static void sync_recompute_rhs(Solver *solver)
{
	FSCGModel *model = solver->model;
	memset(solver->Tr, 0, model->T*sizeof(*solver->Tr));
	for (size_t l=0, L=solver->L; l<L; l++) {
		if (solver->unit[l].r < 0)
			sync_modify_rhs(solver, l, -1);
	}
	sync_update_rhs(solver);
}

static void sync_inst_schd(Solver *solver)
{
	FSCGModel *model = solver->model;
	size_t T = model->T;
	int *nzi = solver->tmpi;
	double *nzx = solver->tmpx;
	size_t ap = solver->a;
	for (size_t a=0, A=solver->a1; a<A; a++) {
		int as = solver->as1[a];
		size_t s = as >= 0 ? as : ~as;
		assert(s < (size_t)solver->pool.N);
		solver->As[ap].s = s;
		solver->As[ap].it = as < 0 ? solver->iter : (solver->iter-1);
		ap++;
		size_t l = solver->pool.l[s];
		size_t j = solver->pool.j[s];
		size_t i = solver->unit[l].i;
		assert(l < (size_t)solver->L);
		assert(i < (size_t)model->I);
		assert(j >= (size_t)model->Ij[i] && j < (size_t)model->Ij[i+1]);
		nzi[0] = solver->R0 + solver->unit[l].r;
		nzx[0] = 1;
		size_t nnz = 1;
		double obj = model->Jo[j];
		for (size_t x=model->Jx[j], X=model->Jx[j+1]; x<X; x++) {
			size_t t = solver->pool.Kt[s*model->K+model->Xk[x]];
			double v = model->Xv[x];
			nzi[nnz] = trow(t);
			nzx[nnz] = v;
			nnz++;
			obj += model->ITu[i*T+t]*v;
			assert(t < T);
			assert((size_t)model->Xk[x] < (size_t)model->K);
			assert(v > 0 || model->Kt[model->Xk[x]] == 1);
		}
		assert(nnz < (size_t)model->K+1);
		lp_add_col(&solver->lp, nnz, nzi, nzx, 0, 1/0.0, solver->state == SS_PHASE1 ? 0 : obj);
	}
	solver->a = ap;
	solver->a1 = 0;
}

static void sync_dual(Solver *solver)
{
	lp_get_dual(&solver->lp, 0, solver->R0+solver->r, solver->dual);
}

static void sync_primal(Solver *solver)
{
	int *basis = solver->tmpi;
	size_t n = lp_get_basis(&solver->lp, basis);
	for (size_t i=0; i<n; i++) {
		int v = basis[i];
		if (v >= solver->A0)
			solver->As[v-solver->A0].it = solver->iter;
	}
}

static double solver_offsetobj(Solver *solver)
{
	FSCGModel *model = solver->model;
	size_t T = model->T;
	double obj = 0;
	for (size_t l=0, L=solver->L; l<L; l++) {
		if (solver->unit[l].ns == 1) {
			size_t s = solver->unit[l].sx;
			size_t j = solver->pool.j[s];
			size_t i = solver->unit[solver->pool.l[s]].i;
			obj += model->Jo[j];
			for (size_t x=model->Jx[j], X=model->Jx[j+1]; x<X; x++) {
				size_t t = solver->pool.Kt[s*model->K+model->Xk[x]];
				obj += model->ITu[i*T+t]*model->Xv[x];
			}
		}
	}
	return solver->scale*obj;
}

static bool solver_sync(Solver *solver)
{
	sync_deinst_schd(solver);
	sync_deinst_units(solver);
	sync_inst_units(solver);
	sync_inst_schd(solver);
	// sync_update_rhs(solver);
	sync_recompute_rhs(solver);
	solver_clock_start(solver);
	if (!lp_solve(&solver->lp)) return false;
	double ofs = solver->state == SS_PHASE1 ? 0 : solver_offsetobj(solver);
	solver->obj = lp_get_obj(&solver->lp) + ofs;
	solver_clock_stop(solver, SCLOCK_MASTER);
	solver->iter++;
	sync_dual(solver);
	sync_primal(solver);
	return true;
}

static void solver_syncinit(Solver *solver, size_t warm)
{
	for (size_t n=0, N=solver->pool.n; n<N; n++) {
		Unit *uni = &solver->unit[solver->pool.l[n]];
		uni->sx ^= n;
		uni->ns++;
	}
	for (size_t n=0, N=solver->pool.n; n<N; n++) {
		Unit *uni = &solver->unit[solver->pool.l[n]];
		if (uni->ns > 1) {
			// treat warm schedules as basic here so cold schedules get yeeted out first
			solver->as1[solver->a1++] = n < warm ? ~n : n;
		}
	}
	for (size_t l=0, L=solver->L; l<L; l++) {
		if (solver->unit[l].ns > 1)
			solver->rl1[solver->r1++] = l;
		else
			sync_modify_rhs(solver, l, -1);
	}
}

static void solver_syncobj(Solver *solver)
{
	FSCGModel *model = solver->model;
	size_t T = model->T;
	int *nzi = solver->tmpi;
	double *nzx = solver->tmpx;
	memcpy(nzi, model->CXc, model->Cnnz[0]*sizeof(*nzi));
	memcpy(nzx, model->CXv, model->Cnnz[0]*sizeof(*nzx));
	size_t nnz = model->Cnnz[0];
	for (size_t a=0, A=solver->a; a<A; a++) {
		size_t s = solver->As[a].s;
		size_t j = solver->pool.j[s];
		size_t i = solver->unit[solver->pool.l[s]].i;
		double obj = model->Jo[j];
		for (size_t x=model->Jx[j], X=model->Jx[j+1]; x<X; x++) {
			size_t t = solver->pool.Kt[s*model->K+model->Xk[x]];
			obj += model->ITu[i*T+t]*model->Xv[x];
		}
		nzi[nnz] = solver->A0 + a;
		nzx[nnz] = obj;
		nnz++;
	}
	lp_set_obj(&solver->lp, nnz, nzi, nzx);
}


/* ---- Scanning ------------------------------------------------------------ */

static void scan_create(ScanState *ss, Solver *solver, int no, int l, int L)
{
	ss->solver = solver;
	pool_create(&ss->pool, solver->model, solver->W, false);
	ss->no = no;
	ss->l = l;
	ss->L = L;
	ss->Ws = malloc(solver->W*sizeof(*ss->Ws));
	ss->t = malloc(solver->model->K*sizeof(*ss->t));
	ss->d = malloc(solver->model->K*sizeof(*ss->d));
}

static void scan_destroy(ScanState *ss)
{
	pool_destroy(&ss->pool);
	free(ss->Ws);
	free(ss->t);
	free(ss->d);
}

static void scan_transport_phase1(Solver *solver, double *dx, int *tx)
{
	FSCGModel *model = solver->model;
	int *Kt = model->Kt;
	// phase 1: utility=0, so the dual for the transport rows is the same for all units.
	for (size_t k=0, t=0, K=model->K; k<K; k++) {
		double dumin = 1/0.0;
		size_t dumint = 0;
		for (size_t ti=0, ktn=Kt[k]; ti<ktn; ti++, t++) {
			double du = solver->dual[trow(t)];
			dumint = du < dumin ? t : dumint;
			dumin = du < dumin ? du : dumin;
		}
		dx[k] = dumin;
		tx[k] = dumint;
	}
}

static void scan_transport_phase2(Solver *solver, double *dx, int *tx, size_t l)
{
	FSCGModel *model = solver->model;
	size_t T = model->T;
	int *Kt = model->Kt;
	size_t i = solver->unit[l].i;
	// phase 2: dual depends on unit because of the utilities.
	for (size_t k=0, t=0, K=model->K; k<K; k++) {
		// this loop is too small to be worth vectorizing so we're just going to brute force it.
		size_t ktn = Kt[k];
		double dumin;
		size_t dumint;
		if (ktn < 4) {
			dumin = solver->dual[trow(t)] - model->ITu[i*T+t];
			dumint = t;
			t++;
			for (size_t ti=1; ti<ktn; ti++, t++) {
				double du = solver->dual[trow(t)] - model->ITu[i*T+t];
				dumint = du < dumin ? t : dumint;
				dumin = du < dumin ? du : dumin;
			}
		} else {
			double dum1 = solver->dual[trow(t+0)] - model->ITu[i*T+t+0];
			double dum2 = solver->dual[trow(t+1)] - model->ITu[i*T+t+1];
			double dum3 = solver->dual[trow(t+2)] - model->ITu[i*T+t+2];
			double dum4 = dum3;
			size_t dut1 = t+0;
			size_t dut2 = t+1;
			size_t dut3 = t+2;
			size_t dut4 = dut3;
			size_t ti = ktn&3;
			t += ti;
			for (; ti<ktn; ti+=4, t+=4) {
				double du1 = solver->dual[trow(t+0)] - model->ITu[i*T+t+0];
				double du2 = solver->dual[trow(t+1)] - model->ITu[i*T+t+1];
				double du3 = solver->dual[trow(t+2)] - model->ITu[i*T+t+2];
				double du4 = solver->dual[trow(t+3)] - model->ITu[i*T+t+3];
				dut1 = du1 < dum1 ? t+0 : dut1;
				dum1 = du1 < dum1 ? du1 : dum1;
				dut2 = du2 < dum2 ? t+1 : dut2;
				dum2 = du2 < dum2 ? du2 : dum2;
				dut3 = du3 < dum3 ? t+2 : dut3;
				dum3 = du3 < dum3 ? du3 : dum3;
				dut4 = du4 < dum4 ? t+3 : dut4;
				dum4 = du4 < dum4 ? du4 : dum4;
			}
			dut1 = dum1 < dum2 ? dut1 : dut2;
			dum1 = dum1 < dum2 ? dum1 : dum2;
			dut3 = dum3 < dum4 ? dut3 : dut4;
			dum3 = dum3 < dum4 ? dum3 : dum4;
			dumint = dum1 < dum3 ? dut1 : dut3;
			dumin = dum1 < dum3 ? dum1 : dum3;
		}
		dx[k] = dumin;
		tx[k] = dumint;
	}
}

static size_t scan_minschd(ScanState *ss, size_t l, double *pdumin)
{
	Solver *solver = ss->solver;
	FSCGModel *model = solver->model;
	double dumin = 1/0.0;
	size_t duminj = 0;
	size_t i = solver->unit[l].i;
	for (size_t j=model->Ij[i], J=model->Ij[i+1]; j<J; j++) {
		// double du = solver->state == SS_PHASE1 ? 0 : -model->Jo[j];
		// for (size_t x=model->Jx[j], X=model->Jx[j+1]; x<X; x++)
		// 	du += ss->d[model->Xk[x]]*model->Xv[x];
		double du1 = (solver->state == SS_PHASE1) ? 0 : -model->Jo[j];
		double du2 = 0;
		double du3 = 0;
		double du4 = 0;
		size_t x = model->Jx[j], X = model->Jx[j+1];
		for (; x+4<=X; x+=4) {
			du1 += ss->d[model->Xk[x+0]]*model->Xv[x+0];
			du2 += ss->d[model->Xk[x+1]]*model->Xv[x+1];
			du3 += ss->d[model->Xk[x+2]]*model->Xv[x+2];
			du4 += ss->d[model->Xk[x+3]]*model->Xv[x+3];
		}
		switch (X-x) {
			case 3: du3 += ss->d[model->Xk[x+2]]*model->Xv[x+2]; // fallthrough
			case 2: du2 += ss->d[model->Xk[x+1]]*model->Xv[x+1]; // fallthrough
			case 1: du1 += ss->d[model->Xk[x+0]]*model->Xv[x+0];
		}
		double du = (du1+du2)+(du3+du4);
		duminj = du < dumin ? j : duminj;
		dumin = du < dumin ? du : dumin;
	}
	*pdumin = dumin;
	return duminj;
}

// static size_t scan_minschd(ScanState *ss, size_t l, double *pdumin)
// {
// 	Solver *solver = ss->solver;
// 	FSCGModel *model = solver->model;
// 	double dumin = 1/0.0;
// 	size_t duminj = 0;
// 	size_t i = solver->unit[l].i;
// 	size_t j=model->Ij[i], J=model->Ij[i+1];
// 	size_t x=model->Jx[j];
// 	// handle length-1 schedules
// 	for (; j<J && x+1 == model->Jx[j+1]; j++, x++) {
// 		double du = (solver->state == SS_PHASE1) ? 0 : -model->Jo[j];
// 		du += ss->d[model->Xk[x]]*model->Xv[x];
// 		duminj = du < dumin ? j : duminj;
// 		dumin = du < dumin ? du : dumin;
// 	}
// 	// handle length>1 schedules
// 	for (; j<J; j++) {
// 		// double du = solver->phase1 ? 0 : -model->Jo[j];
// 		// for (size_t x=model->Jx[j], X=model->Jx[j+1]; x<X; x++)
// 		// 	du += solver->d[model->Xk[x]]*model->Xv[x];
// 		double du1 = (solver->state == SS_PHASE1) ? 0 : -model->Jo[j];
// 		double du2 = 0;
// 		double du3 = 0;
// 		double du4 = 0;
// 		size_t X = model->Jx[j+1];
// 		for (; x+4<=X; x+=4) {
// 			du1 += ss->d[model->Xk[x+0]]*model->Xv[x+0];
// 			du2 += ss->d[model->Xk[x+1]]*model->Xv[x+1];
// 			du3 += ss->d[model->Xk[x+2]]*model->Xv[x+2];
// 			du4 += ss->d[model->Xk[x+3]]*model->Xv[x+3];
// 		}
// 		switch (X-x) {
// 			case 3: du3 += ss->d[model->Xk[x+2]]*model->Xv[x+2]; // fallthrough
// 			case 2: du2 += ss->d[model->Xk[x+1]]*model->Xv[x+1]; // fallthrough
// 			case 1: du1 += ss->d[model->Xk[x+0]]*model->Xv[x+0];
// 		}
// 		x = X;
// 		double du = (du1+du2)+(du3+du4);
// 		duminj = du < dumin ? j : duminj;
// 		dumin = du < dumin ? du : dumin;
// 	}
// 	*pdumin = dumin;
// 	return duminj;
// }


// static size_t scan_minschd(ScanState *ss, size_t l, double *pdumin)
// {
// 	Solver *solver = ss->solver;
// 	FSCGModel *model = solver->model;
// 	double dumin1 = 1/0.0, dumin2 = 1/0.0;
// 	size_t duminj1 = 0, duminj2 = 0;
// 	size_t i = solver->unit[l].i;
// 	for (size_t j=model->Ij[i], J=model->Ij[i+1]; j+1<J; j+=2) {
// 		double du = solver->state == SS_PHASE1 ? 0 : -model->Jo[j];
// 		for (size_t x=model->Jx[j], X=model->Jx[j+1]; x<X; x++)
// 			du += ss->d[model->Xk[x]]*model->Xv[x];
// 		duminj1 = du < dumin1 ? j : duminj1;
// 		dumin1 = du < dumin1 ? du : dumin1;
// 		du = solver->state == SS_PHASE1 ? 0 : -model->Jo[j+1];
// 		for (size_t x=model->Jx[j+1], X=model->Jx[j+2]; x<X; x++)
// 			du += ss->d[model->Xk[x]]*model->Xv[x];
// 		duminj2 = du < dumin2 ? j+1 : duminj2;
// 		dumin2 = du < dumin2 ? du : dumin2;
// 	}
// 	if ((model->Ij[i+1]-model->Ij[i]) & 1) {
// 		size_t j = model->Ij[i+1]-1;
// 		double du = solver->state == SS_PHASE1 ? 0 : -model->Jo[j];
// 		for (size_t x=model->Jx[j], X=model->Jx[j+1]; x<X; x++)
// 			du += ss->d[model->Xk[x]]*model->Xv[x];
// 		duminj1 = du < dumin1 ? j : duminj1;
// 		dumin1 = du < dumin1 ? du : dumin1;
// 	}
// 	*pdumin = dumin1 < dumin2 ? dumin1 : dumin2;
// 	return dumin1 < dumin2 ? duminj1 : duminj2;
// }

// static size_t scan_minschd(ScanState *ss, size_t l, double *pdumin)
// {
// 	Solver *solver = ss->solver;
// 	FSCGModel *model = solver->model;
// 	double dumin = 1/0.0;
// 	size_t duminj = 0;
// 	size_t i = solver->unit[l].i;
// 	size_t j = model->Ij[i];
// 	size_t J = model->Ij[i+1];
// 	size_t j1=j, j2=j+1<J?j+1:j, j3=j2+1<J?j2+1:j2, j4=j3+1<J?j3+1:j3;
// 	size_t x1=model->Jx[j1], x2=model->Jx[j2], x3=model->Jx[j3], x4=model->Jx[j4];
// 	size_t X1=model->Jx[j1+1], X2=model->Jx[j2+1], X3=model->Jx[j3+1], X4=model->Jx[j4+1];
// 	double du1 = solver->state == SS_PHASE1 ? 0 : -model->Jo[j1];
// 	double du2 = solver->state == SS_PHASE1 ? 0 : -model->Jo[j2];
// 	double du3 = solver->state == SS_PHASE1 ? 0 : -model->Jo[j3];
// 	double du4 = solver->state == SS_PHASE1 ? 0 : -model->Jo[j4];
// 	for (;;) {
// 		assert(X1-x1 <= X2-x2); // invariant
// 		if (x1 == X1) {
// 			duminj = du1 < dumin ? j1 : duminj;
// 			dumin = du1 < dumin ? du1 : dumin;
// 			if (j1 == j2) break;
// 			x1 = x2; X1 = X2; j1 = j2; du1 = du2;
// 			x2 = x3; X2 = X3; j2 = j3; du2 = du3;
// 			x3 = x4; X3 = X4; j3 = j4; du3 = du4;
// 			j4 = j4+1 < J ? j4+1 : j4;
// 			x4 = model->Jx[j4];
// 			X4 = model->Jx[j4+1];
// 			du4 = solver->state == SS_PHASE1 ? 0 : -model->Jo[j4];
// 			continue;
// 		}
// 		du1 += ss->d[model->Xk[x1]]*model->Xv[x1]; x1++;
// 		du2 += ss->d[model->Xk[x2]]*model->Xv[x2]; x2++;
// 		du3 += ss->d[model->Xk[x3]]*model->Xv[x3]; x3++;
// 		du4 += ss->d[model->Xk[x4]]*model->Xv[x4]; x4++;
// 	}
// 	*pdumin = dumin;
// 	return duminj;
// }

static WSchd scan_put_heap(WSchd *Ws, size_t W, size_t h, uint32_t k)
{
	WSchd wh = Ws[h] & ~WSCHD_K;
	for (;;) {
		size_t l = 2*h + 1;
		size_t r = 2*h + 2;
		size_t m = h;
		uint32_t mk = k;
		if (l < W && wschd_k(Ws[l]) > mk) { m = l; mk = wschd_k(Ws[l]); }
		if (r < W && wschd_k(Ws[r]) > mk) m = r;
		if (m == h) {
			return Ws[h] = wh | wschd_knew(k);
		} else {
			Ws[h] = Ws[m];
			h = m;
		}
	}
}

// transports are in ss->t
static void scan_put(ScanState *ss, size_t l, size_t j, double du)
{
	Solver *solver = ss->solver;
	if (du > solver->gamma*ss->wd)
		return; // quick reject because it won't be selected anyway
	uint32_t k = wschd_fpk(du);
	WSchd *Ws = ss->Ws;
	size_t s;
	if (ss->pool.n < ss->pool.N) {
		// worklist is linear, work indices match pool indices
		s = pool_alloc(&ss->pool);
		assert(s == (size_t)ss->pool.n-1);
		Ws[s] = wschd_new(k, ss->no, s);
		// pool just became full?
		if (ss->pool.n == ss->pool.N) {
			// heapify work list
			size_t w = ss->pool.N >> 1;
			for (;;) {
				scan_put_heap(Ws-1, ss->pool.N, w+1, wschd_k(Ws[w]));
				if (!w--) break;
			}
		}
	} else {
		// work space is heapified
		if (k >= wschd_k(Ws[0]))
			return; // offered schedule is worse than the worst one we already have
		s = wschd_s(scan_put_heap(Ws-1, ss->pool.N, 1, k));
	}
	if (du < ss->wd)
		ss->wd = du;
	ss->pool.l[s] = l;
	ss->pool.j[s] = j;
	size_t K = solver->model->K;
	memcpy(ss->pool.Kt+s*K, ss->t, K*sizeof(*ss->t));
}

static void scan_schd(ScanState *ss, size_t l)
{
	Solver *solver = ss->solver;
	FSCGModel *model = solver->model;
	size_t i = solver->unit[l].i;
	// determine best schedule
	double dumin;
	size_t duminj = scan_minschd(ss, l, &dumin);
	// compute area constraint dual
	double wdu;
	if (solver->unit[l].r >= 0) {
		wdu = solver->dual[solver->R0+solver->unit[l].r];
	} else {
		assert(solver->unit[l].ns == 1);
		size_t s = solver->unit[l].sx;
		size_t j = solver->pool.j[s];
		size_t T = model->T;
		if (solver->state == SS_PHASE1) {
			wdu = 0;
			for (size_t x=model->Jx[j], X=model->Jx[j+1]; x<X; x++) {
				size_t t = solver->pool.Kt[s*model->K+model->Xk[x]];
				wdu -= model->Xv[x]*solver->dual[trow(t)];
			}
		} else {
			wdu = model->Jo[j];
			for (size_t x=model->Jx[j], X=model->Jx[j+1]; x<X; x++) {
				size_t t = solver->pool.Kt[s*model->K+model->Xk[x]];
				wdu += model->Xv[x]*(model->ITu[i*T+t] - solver->dual[trow(t)]);
			}
		}
	}
	dumin += wdu;
	// add schedule if it's improving
	if (dumin < -solver->epsilon) {
		ss->lag += dumin;
		scan_put(ss, l, duminj, dumin);
	}
}

static void scan_thread_scan(ScanState *ss)
{
	Solver *solver = ss->solver;
	ss->lag = 0;
	ss->wd = 0;
	pool_reset(&ss->pool);
	if (solver->state == SS_PHASE1)
		scan_transport_phase1(solver, ss->d, ss->t);
	for (size_t l=ss->l, L=l+ss->L; l<L; l++) {
		if (solver->state == SS_PHASE2)
			scan_transport_phase2(solver, ss->d, ss->t, l);
		scan_schd(ss, l);
	}
}

static void *scan_thread_run(ScanState *ss)
{
	Solver *solver = ss->solver;
	for (;;) {
		pthread_barrier_wait(&solver->bar);  // wait for command
		if (solver->state == SS_EXIT)
			return NULL;
		scan_thread_scan(ss);
		pthread_barrier_wait(&solver->bar);  // wait until all work threads are done
	}
}

static void qselect64(uint64_t *xs, size_t u, size_t n)
{
	for (;;) {
		if (u <= n)
			return; // nothing to do
		if (u < 16 || n < 5) {
			// selection sort for small arrays
			for (size_t i=0; i<n; i++) {
				size_t m = xs[i];
				size_t mi = i;
				for (size_t j=i+1; j<u; j++) {
					mi = m < xs[j] ? mi : j;
					m = m < xs[j] ? m : xs[j];
				}
				xs[mi] = xs[i];
				xs[i] = m;
			}
			return;
		}
		// pivot = pseudo-median of xs[0], xs[u/2], xs[u-1]
		uint64_t a = xs[0];
		uint64_t b = xs[u>>1];
		uint64_t c = xs[u-1];
		uint64_t p = a < b ? (b < c ? b : (a < c ? c : a)) : (a < c ? a : (b < c ? c : b));
		size_t i = 0;
		size_t j = u-1;
		// partition left <= p <= right
		for (;;) {
			while (xs[i] <= p && i < j) i++;
			while (xs[j] >= p && i < j) j--;
			if (i == j) break;
			uint64_t xi = xs[i];
			xs[i] = xs[j];
			xs[j] = xi;
		}
		if (i < n) {
			// all elements on the left are selected
			xs += i;
			u -= i;
			n -= i;
		} else {
			// no elements on the right are selected
			u = i;
		}
	}
}

static bool solver_exists(ScanState *ss, SchPool *pool, size_t s)
{
	Solver *solver = ss->solver;
	FSCGModel *model = solver->model;
	size_t l = pool->l[s];
	bool exists;
	if (solver->unit[l].ns == 1) {
		size_t ss = solver->unit[l].sx;
		exists = pool_cmp(model, &solver->pool, ss, pool, s);
	} else {
		exists = false;
		for (size_t a=0, A=solver->a; a<A; a++) {
			size_t ss = solver->As[a].s;
			if (pool_cmp(model, &solver->pool, ss, pool, s)) {
				exists = true;
				break;
			}
		}
	}
	// fprintf(stderr, "CHECK %d %d (%d) [", (int)l, (int)j, solver->unit[l].ns);
	// for (size_t k=0; k<K; k++)
	// 	fprintf(stderr, " %d", pool->Kt[s*K+k]);
	// fprintf(stderr, " ] -> %d\n", (int)exists);
	return exists;
}

static size_t solver_commit(Solver *solver)
{
	// collect all schedules with du<=gamma*min(wd) into solver->tmpq
	double wd = 0;
	for (size_t s=0, S=solver->nss; s<S; s++)
		wd = solver->ss[s].wd < wd ? solver->ss[s].wd : wd;
	uint32_t k = wschd_fpk(solver->gamma*wd);
	bool check = wd > -0.01;
	size_t ns = 0;
	for (size_t s=0, S=solver->nss; s<S; s++) {
		// either the pool isn't filled, and the schedules are in linear order,
		// or the pool is filled, and every entry contains a schedule.
		// either way, all schedules in 0..n are valid.
		ScanState *ss = &solver->ss[s];
		for (size_t w=0, W=ss->pool.n; w<W; w++) {
			if (wschd_k(ss->Ws[w]) <= k && !(check && solver_exists(ss, &ss->pool, wschd_s(ss->Ws[w]))))
				solver->tmpq[ns++] = ss->Ws[w];
		}
	}
	// we currently have up to `space` slots available for instantiation,
	// and we would like to instantiate up to `inst` slots.
	size_t W = solver->W;
	size_t space = solver->pool.N - solver->pool.n;
	size_t inst = ns < W ? ns : W;
	// need to make space?
	if (inst > space) {
		// deinstantiate up to `space-inst` slots
		uint64_t *de = solver->tmpq + ns;
		size_t nde = 0; // up to `nde` slots are eligible for deinstantiation
		for (size_t a=0, A=solver->a; a<A; a++) {
			if (solver->As[a].it < solver->iter)
				de[nde++] = a | ((uint64_t)solver->As[a].it << 32);
		}
		if (nde > inst-space) {
			// more eligible slots than we need to deinstantiate, pick the worst ones
			qselect64(de, nde, inst-space);
			nde = inst-space;
		}
		// deinstantiate `nde` slots
		for (size_t d=0; d<nde; d++)
			solver_deinst(solver, (uint32_t) de[d]);
		space += nde;
	}
	// trace("eligible=%d  want=%d  space=%d  wd=%g", (int)ns, (int)inst, (int)space, wd);
	// instantiate up to `space` slots
	if (inst > space)
		inst = space;
	if (inst < ns) {
		// we can't instantiate everything, choose the best ones
		qselect64(solver->tmpq, ns, inst);
		ns = inst;
	}
	// instantiate `ns` slots
	for (size_t w=0; w<ns; w++) {
		size_t ws = solver->tmpq[w];
		int s = pool_transfer(solver->model, &solver->pool, &solver->ss[wschd_ss(ws)].pool,
			wschd_s(ws), 0);
		solver_inst(solver, s);
	}
	return ns;
}

static size_t solver_scan(Solver *solver)
{
	solver_clock_start(solver);
	if (solver->nss > 1) {
		pthread_barrier_wait(&solver->bar);  // start scanning
		pthread_barrier_wait(&solver->bar);  // wait until scanning finishes
	} else {
		scan_thread_scan(solver->ss);
	}
	solver_clock_stop(solver, SCLOCK_SCAN);
	// collect lagrangian bounds
	solver->lag = 0;
	for (size_t s=0, S=solver->nss; s<S; s++)
		solver->lag += solver->ss[s].lag;
	double ub = solver->obj - solver->scale*solver->lag;
	if (ub < solver->ub) solver->ub = ub;
	return solver_commit(solver);
}

/* ---- RNG ----------------------------------------------------------------- */

typedef uint32_t Rng;

static uint32_t rng_next(Rng *rng)
{
	uint32_t x = *rng;
	*rng = x * 0x93d765dd;
	return x;
}

static float rng_rand12(Rng *rng)
{
	return sem2fp(0, 127, rng_next(rng) >> (32-23));
}

float rng_unif(Rng *rng, float a, float b)
{
	return (2*a-b) + (b-a)*rng_rand12(rng);
}

/* ---- Initialization ------------------------------------------------------ */

static double solver_gettime(Solver *solver)
{
	struct timespec tp;
	clock_gettime(SOLVER_CLOCK, &tp);
	return (tp.tv_sec-solver->tp0.tv_sec) + (tp.tv_nsec-solver->tp0.tv_nsec)*1e-9;
}

static void solver_init(Solver *solver)
{
	solver_clock_start(solver);
	solver->state = SS_PHASE2;
	FSCGModel *model = solver->model;
	ScanState *ss = solver->ss;
	size_t T = model->T, C = model->C;
	pool_reset(&solver->pool);
	memset(solver->pool.Kt, 0, solver->L*model->K*sizeof(*solver->pool.Kt));
	for (size_t l=0, L=solver->L; l<L; l++) {
		size_t s = pool_alloc(&solver->pool);
		assert(s == l); (void)s;
		solver->pool.l[l] = l;
		solver->pool.j[l] = 0;
	}
	int *SKtbuf = malloc(solver->L*model->K*sizeof(*SKtbuf));
	int *lbuf = malloc(solver->L*sizeof(*lbuf));
	size_t *jbuf = malloc(solver->L*sizeof(*jbuf));
	int *tabu = calloc(solver->L, sizeof(*tabu));
	int *order = malloc(solver->L*sizeof(*order));
	for (size_t l=0, L=solver->L; l<L; l++)
		order[l] = l;
	// solver->dual is partitioned as:
	//   0..T-1     transport (T) duals
	//   T..T+C-1   general constraint (C) duals
	memset(solver->dual+T, 0, C*sizeof(*solver->dual));
	// work arrays:
	double *tx = calloc(T, sizeof(*tx));
	double *tx2 = malloc(T*sizeof(*tx2));
	double *dual = calloc(C, sizeof(*dual));
	double *dual2 = malloc(C*sizeof(*dual2));
	double *rho = malloc(2*C*sizeof(*rho));
	LP lp;
	lp_create(&lp, LP_SIMPLEX);
	// first Z columns are z-variables
	for (int z=0; z<model->Z; z++)
		lp_add_col(&lp, 0, NULL, NULL, z >= model->Za ? -1/0.0 : 0, 1/0.0, 0);
	// first C rows are general constraint rows (z-variables only)
	for (size_t c=1, cx=model->Cnnz[0]; c<=C; cx+=model->Cnnz[c++]) {
		size_t nnz = 0;
		for (size_t x=0, X=model->Cnnz[c]; x<X; x++) {
			size_t col = model->CXc[cx+x];
			if (col >= T) {
				// only z-variables
				solver->tmpi[nnz] = col-T;
				solver->tmpx[nnz] = model->CXv[cx+x];
				nnz++;
			}
		}
		lp_add_row(&lp, nnz, solver->tmpi, solver->tmpx, model->Cs[c], 0);
	}
	// slacks
	// double rho = 1e2;
	double rho0 = 1e3;
	// double rhomin = 1e2;
	size_t nslack = 0;
	double objk = 1;
	for (int c=1, cx=model->Cnnz[0]; c<=model->C; cx+=model->Cnnz[c++]) {
			// slack only rows that either mention a t-variable or are artificial z-variable definitions
			if (c < model->Ca) {
				for (int x=0; x<model->Cnnz[c]; x++) {
					if (model->CXc[cx+x] < model->T)
						goto slack;
				}
				continue; // no t-variable
			}
slack:
		char s = model->Cs[c];
		if (s == '>' || s == '=') {
			rho[nslack] = rho0;
			lp_add_col(&lp, 1, (int[]){c-1}, (double[]){1}, 0, 1/0.0, -rho[nslack]);
			nslack++;
		}
		if (s == '<' || s == '=') {
			rho[nslack] = rho0;
			lp_add_col(&lp, 1, (int[]){c-1}, (double[]){-1}, 0, 1/0.0, -rho[nslack]);
			nslack++;
		}
	}
	// objective (z-variables)
	size_t onnz = 0;
	for (int x=0; x<model->Cnnz[0]; x++) {
		size_t col = model->CXc[x];
		if (col >= T) {
			solver->tmpi[onnz] = col-T;
			solver->tmpx[onnz] = objk*model->CXv[x];
			onnz++;
		}
	}
	lp_set_obj(&lp, onnz, solver->tmpi, solver->tmpx);
	int noimprove = 0;
	// int noimprovelimit = 100;
	int noimprovelimit = solver->L / model->C;
	if (noimprovelimit < 100) noimprovelimit = 100;
	int addstep = model->C;
	if (addstep > solver->L/100) addstep = solver->L/100;
	if (addstep < 10) addstep = 10;
	double scale = solver->scale;
	size_t cursor = 0;
	double sobj = 0, obj = -1/0.0;
	int txupd = 0;
	double dulim = -0.01;
	Rng rng = 314159;
	for (int it=0; /*it<200*/; it++) {
		// compute T duals
		memset(solver->dual, 0, T*sizeof(*solver->dual));
		for (size_t c=1, cx=model->Cnnz[0]; c<=C; cx+=model->Cnnz[c++]) {
			for (size_t x=0, X=model->Cnnz[c]; x<X; x++) {
				size_t col = model->CXc[cx+x];
				if (col < T)
					solver->dual[col] += solver->dual[T+c-1]*model->CXv[cx+x];
			}
		}
		// solve pricing problem
		size_t maxadd = it == 0 ? solver->L : addstep;
		size_t nadd = 0;
		double sobj2 = sobj;
		double dusum = 0;
		memcpy(tx2, tx, T*sizeof(*tx));
		for (size_t ll=0, L=solver->L; ll<L; ll++) {
			size_t l = order[cursor++];
			if (cursor == L) {
				// reshuffle
				for (size_t li=1; li<L; li++) {
					size_t lj = (size_t)(li+1)*rng_unif(&rng, 0, 1); // lj ∈ [0, li]
					int tmp = order[li];
					order[li] = order[lj];
					order[lj] = tmp;
				}
				cursor = 0;
			}
			if (tabu[l] > it)
				continue;
			tabu[l] = it+1; // save some time by not redoing this unit on this iteration
			size_t i = solver->unit[l].i;
			// scan transports (TODO: call scan_*)
			for (size_t k=0, t=0, K=model->K; k<K; k++) {
				double dumin = 1/0.0;
				size_t dumint = 0;
				for (size_t ti=0, ktn=model->Kt[k]; ti<ktn; ti++, t++) {
					double du = solver->dual[trow(t)] - objk*model->ITu[i*T+t];
					dumint = du < dumin ? t : dumint;
					dumin = du < dumin ? du : dumin;
				}
				ss->d[k] = dumin;
				ss->t[k] = dumint;
			}
			// scan schedules
			double dumin;
			size_t duminj = scan_minschd(ss, l, &dumin);
			// compute area constraint dual
			size_t j = solver->pool.j[l];
			double wdu = objk*model->Jo[j];
			for (size_t x=model->Jx[j], X=model->Jx[j+1]; x<X; x++) {
				size_t t = solver->pool.Kt[l*model->K+model->Xk[x]];
				wdu += model->Xv[x]*(objk*model->ITu[i*T+t] - solver->dual[trow(t)]);
			}
			// trace("cursor = %d; dumin+wdu = %g", (int)cursor, dumin+wdu);
			if (dumin+wdu < dulim || it == 0) {
				dusum += dumin+wdu;
				// it's improving
				if (it > 0) {
					// remove old schedule
					sobj2 -= scale*model->Jo[j];
					for (size_t x=model->Jx[j], X=model->Jx[j+1]; x<X; x++) {
						size_t t = solver->pool.Kt[l*model->K+model->Xk[x]];
						tx2[t] -= scale*model->Xv[x];
						sobj2 -= scale*model->Xv[x]*model->ITu[i*T+t];
					}
				}
				// insert new schedule
				lbuf[nadd] = l;
				jbuf[nadd] = duminj;
				memcpy(SKtbuf+nadd*model->K, ss->t, model->K*sizeof(*SKtbuf));
				sobj2 += scale*model->Jo[duminj];
				for (size_t x=model->Jx[duminj], X=model->Jx[duminj+1]; x<X; x++) {
					size_t t = ss->t[model->Xk[x]];
					tx2[t] += scale*model->Xv[x];
					sobj2 += scale*model->Xv[x]*model->ITu[i*T+t];
				}
				tabu[l] = it+10;
				if (++nadd >= maxadd) break;
			}
		}
		// recompute rhs & solve z-variables
		memcpy(solver->tmpx, model->Cr+1, C*sizeof(*solver->tmpx));
		for (size_t c=1, cx=model->Cnnz[0]; c<=C; cx+=model->Cnnz[c++]) {
			solver->tmpi[c-1] = c-1;
			for (size_t x=0, X=model->Cnnz[c]; x<X; x++) {
				size_t col = model->CXc[cx+x];
				if (col < T)
					solver->tmpx[c-1] -= tx2[col]*model->CXv[cx+x];
			}
		}
		lp_change_rhs(&lp, model->C, solver->tmpi, solver->tmpx);
		lp_solve(&lp);
		lp_get_dual(&lp, 0, model->C, dual2);
		// did objective improve?
		double lpobj = lp_get_obj(&lp);
		double alpha = 0.001;
		if (sobj2+lpobj > obj) {
			// improved, replace center and schedules
			sobj = sobj2;
			obj = sobj2+lpobj;
			memcpy(dual, dual2, C*sizeof(*dual));
			for (size_t s=0; s<nadd; s++) {
				size_t l = lbuf[s];
				tabu[l] = 0;
				solver->pool.j[l] = jbuf[s];
				memcpy(solver->pool.Kt+l*model->K, SKtbuf+s*model->K, model->K*sizeof(*solver->pool.Kt));
			}
			noimprove = 0;
			if (++txupd == 100) {
				txupd = 0;
				// recompute tx for numerical stability
				memset(tx, 0, T*sizeof(*tx));
				for (size_t l=0, L=solver->L; l<L; l++) {
					size_t j = solver->pool.j[l];
					for (size_t x=model->Jx[j], X=model->Jx[j+1]; x<X; x++)
						tx[solver->pool.Kt[l*model->K+model->Xk[x]]] += scale*model->Xv[x];
				}
			} else {
				memcpy(tx, tx2, T*sizeof(*tx));
			}
			double residual = 0; // used for logging only
			lp_get_primal(&lp, model->Z, nslack, solver->tmpx);
			for (size_t s=0; s<nslack; s++)
				residual += solver->tmpx[s];
			double time = solver_gettime(solver);
			trace("%-4d %4ds  obj: %.10e  resid: %.10e  step: %d",
				it, (int)time, obj, residual, addstep);
			alpha *= 0.1;
			dulim = nadd > 0 ? 0.1*dusum/nadd : -0.01;
		} else {
			memcpy(dual, dual2, C*sizeof(*dual));
			dulim *= 0.95;
			if (++noimprove == noimprovelimit) {
				noimprove = 0;
				addstep /= 2;
				// addstep = 2*addstep/3;
				if (addstep <= 2)
					break;
			}
		}
		if (dulim > -0.01) dulim = -0.01;
		// update dual
		for (size_t c=0; c<(size_t)model->C; c++)
			solver->dual[T+c] = alpha*dual[c] + (1-alpha)*solver->dual[T+c];
		// update penalty
		// lp_get_primal(&lp, model->Z, nslack, solver->tmpx);
		// for (size_t s=0; s<nslack; s++) {
		// 	rho[s] = alpha*solver->tmpx[s] + (1-alpha)*rho[s];
		// 	if (rho[s] < rhomin) rho[s] = rhomin;
		// 	solver->tmpi[s] = model->Z+s;
		// 	solver->tmpx[s] = -rho[s];
		// }
		// lp_set_obj(&lp, nslack, solver->tmpi, solver->tmpx);
		// trace("it %d [obj %e %e %e] [best %e] [step %d] (%d)",
		// 	it, sobj2+lpobj, sobj2, lpobj, obj, addstep, noimprove);
	}
	solver_clock_stop(solver, SCLOCK_INIT);
	free(SKtbuf);
	free(lbuf);
	free(jbuf);
	free(tabu);
	free(order);
	free(tx);
	free(tx2);
	free(dual);
	free(dual2);
	free(rho);
}

/* ---- Solver -------------------------------------------------------------- */

static void solver_log(Solver *solver, int nnew, double time)
{
	if (solver->log) {
		fprintf(solver->log, "%d\t%d\t%.10e\t%.10e\t%.10e\t%g", solver->iter-1, nnew, solver->obj,
			solver->ub-solver->obj, -solver->scale*solver->lag, time);
		for (int c=0; c<SCLOCK__MAX; c++)
			fprintf(solver->log, "\t%g", solver->clock[c]);
		fputc('\n', solver->log);
	}
}

static bool solver_checkstop(Solver *solver)
{
	if (solver->maxgap <= 0 || solver->ub-solver->obj > solver->obj*solver->maxgap)
		return false;
	// we are close enough to optimal to stop, but are the sibling subsolvers still working?
	if (solver->mask) {
		uint64_t m = __atomic_or_fetch(solver->mask, 1ULL << solver->no, __ATOMIC_RELAXED);
		if (m != ((1ULL << solver->nsub)-1))
			return false; // siblings are not done yet
	}
	// we can stop now
	return true;
}

static bool solver_run(Solver *solver, int iterlimit)
{
	solver->ub = 1/0.0;
	for (int it=0; it!=iterlimit; it++) {
		bool ok = solver_sync(solver);
		assert(ok); (void) ok;
		assert(solver->lp.ncol == solver->A0+solver->a);
		assert(solver->lp.nrow == solver->R0+solver->r);
		if (iterlimit == -2 && solver_checkstop(solver))
			return true;
		int nnew = solver_scan(solver);
		double time = solver_gettime(solver);
		trace("%-4d %4ds %.10e %.10e %.10e [+%d]", solver->iter-1, (int)time,
			solver->obj, solver->ub-solver->obj, -solver->scale*solver->lag, nnew);
		solver_log(solver, nnew, time);
		if (!nnew) {
			solver->ub = solver->obj;
			return true; // optimal
		}
	}
	return false;
}

#define SOLVER_LOG    0x1
#define SOLVER_APPEND 0x2
#define SOLVER_SERIAL 0x4

static void solver_create(Solver *solver, FSCGModel *model, int L, int NN, int flags, int nss)
{
	memset(solver, 0, sizeof(*solver));
	solver->iter = 1;
	clock_gettime(SOLVER_CLOCK, &solver->tp0);
	lp_create(
		&solver->lp,
		LP_BARRIER
		| LP_WARM
		| ((flags & SOLVER_SERIAL) ? LP_SERIAL : 0)
	);
	solver->model = model;
	solver->epsilon = model->epsilon;
	solver->gamma = model->gamma;
	solver->scale = model->scale * ((double)model->I / (double)L);
	solver->L = L;
	solver->nss = nss > 1 ? nss : 1;
	int N = L + (NN ? NN : 30*sqrt(L));
	pool_create(&solver->pool, model, N, true);
	solver->W = model->s_W ? model->s_W : 2*(solver->pool.N-solver->L)/3;
	if (solver->W > L) solver->W = L;
	solver->Tmaskbytes = ((model->T+63)>>6)<<3;
	solver->R0 = model->T + model->C;
	solver->Rmax = solver->pool.N - solver->L;
	solver->As = malloc(solver->pool.N*sizeof(*solver->As));
	solver->as1 = malloc(solver->pool.N*sizeof(*solver->as1));
	solver->unit = calloc(solver->L, sizeof(*solver->unit));
	solver->rl = malloc(solver->Rmax*sizeof(*solver->rl));
	solver->rl1 = malloc(solver->Rmax*sizeof(*solver->rl1));
	solver->Tr = calloc(model->T, sizeof(*solver->Tr));
	solver->Trx = calloc(solver->Tmaskbytes, 1);
	solver->dual = malloc((model->T+model->C+solver->Rmax)*sizeof(*solver->dual));
	// TODO figure this out
	size_t sizetmp = solver->pool.N + solver->nss*solver->W + model->C + model->T + model->Cnnz[0];
	solver->tmpi = malloc(sizetmp*sizeof(*solver->tmpi));
	solver->tmpx = malloc(sizetmp*sizeof(*solver->tmpx));
	solver->tmpq = malloc(sizetmp*sizeof(*solver->tmpq));
	solver->ss = malloc(solver->nss*sizeof(*solver->ss));
	solver->st = solver->nss > 1 ? malloc(solver->nss*sizeof(*solver->st)) : NULL;
	if (solver->nss > 1)
		pthread_barrier_init(&solver->bar, NULL, solver->nss+1);
	for (int s=0, S=solver->nss, l=0, L=solver->L/S; s<S; s++, l+=L) {
		scan_create(&solver->ss[s], solver, s, l, s == solver->nss-1 ? solver->L-l : L);
		if (solver->nss > 1)
			pthread_create(&solver->st[s], NULL, (void *)scan_thread_run, &solver->ss[s]);
	}
	solver->obj = -1/0.0;
	for (int l=0; l<L; l++)
		solver->unit[l].r = -1;
	// first T columns are t-variables
	for (int t=0; t<model->T; t++)
		lp_add_col(&solver->lp, 0, NULL, NULL, -1/0.0, 1/0.0, 0);
	// next Z columns are z-variables
	for (int z=0; z<model->Z; z++)
		lp_add_col(&solver->lp, 0, NULL, NULL, z >= model->Za ? -1/0.0 : 0, 1/0.0, 0);
	// first T rows are t-variable definitions
	for (int t=0; t<model->T; t++)
		lp_add_row(&solver->lp, 1, (int[]){t}, (double[]){-1}, '=', 0);
	// nect C rows are general constraints
	for (int c=1, cx=model->Cnnz[0]; c<=model->C; cx+=model->Cnnz[c++])
		lp_add_row(&solver->lp, model->Cnnz[c], model->CXc+cx, model->CXv+cx, model->Cs[c],
			model->Cr[c]);
	if (flags & SOLVER_LOG) {
		const char *logname = getenv("FSCG_LOG");
		if (logname)
			solver->log = fopen(logname, (flags & SOLVER_APPEND) ? "a" : "w");
	}
}

static void solver_destroy(Solver *solver)
{
	if (solver->nss > 1) {
		// send exit command to threads
		solver->state = SS_EXIT;
		pthread_barrier_wait(&solver->bar);
		for (size_t s=0, S=solver->nss; s<S; s++)
			pthread_join(solver->st[s], NULL);
		pthread_barrier_destroy(&solver->bar);
	}
	for (size_t s=0, S=solver->nss; s<S; s++)
		scan_destroy(&solver->ss[s]);
	lp_destroy(&solver->lp);
	pool_destroy(&solver->pool);
	if (solver->log)
		fclose(solver->log);
	free(solver->As);
	free(solver->as1);
	free(solver->unit);
	free(solver->rl);
	free(solver->rl1);
	free(solver->Tr);
	free(solver->Trx);
	free(solver->dual);
	free(solver->tmpi);
	free(solver->tmpx);
	free(solver->tmpq);
	free(solver->ss);
	free(solver->st);
}

static void solver_initschd(Solver *solver, size_t l)
{
	size_t s = pool_alloc(&solver->pool);
	FSCGModel *model = solver->model;
	size_t r = s;
	size_t K = model->K;
	for (size_t k=0, t=0; k<K; t+=model->Kt[k++])
		solver->pool.Kt[s*K+k] = t + (r++)%model->Kt[k];
	solver->pool.l[s] = l;
	Unit *uni = &solver->unit[l];
	size_t i = uni->i;
	solver->pool.j[s] = model->Ij[i] + r%(model->Ij[i+1]-model->Ij[i]);
}

static int solver_initgreedy(Solver *solver, size_t l)
{
	size_t s = pool_alloc(&solver->pool);
	FSCGModel *model = solver->model;
	size_t i = solver->unit[l].i;
	size_t T = model->T;
	int *Kt = model->Kt;
	for (size_t k=0, t=0, K=model->K; k<K; k++) {
		double maxobj = -1/0.0;
		size_t maxt = 0;
		for (size_t ti=0, ktn=Kt[k]; ti<ktn; ti++, t++) {
			double obj = model->ITu[i*T+t];
			maxt = obj > maxobj ? t : maxt;
			maxobj = obj > maxobj ? obj : maxobj;
		}
		solver->tmpx[k] = maxobj;
		solver->pool.Kt[s*K+k] = maxt;
	}
	double maxobj = -1/0.0;
	size_t maxj = 0;
	for (size_t j=model->Ij[i], J=model->Ij[i+1]; j<J; j++) {
		double obj = model->Jo[j];
		for (size_t x=model->Jx[j], X=model->Jx[j+1]; x<X; x++)
			obj += solver->tmpx[model->Xk[x]]*model->Xv[x];
		maxj = obj > maxobj ? j : maxj;
		maxobj = obj > maxobj ? obj : maxobj;
	}
	solver->pool.l[s] = l;
	solver->pool.j[s] = maxj;
	return s;
}

static int solver_initwarm(Solver *solver, size_t l)
{
	scan_transport_phase2(solver, solver->tmpx, solver->tmpi, l);
	size_t s = pool_alloc(&solver->pool);
	FSCGModel *model = solver->model;
	for (int k=0; k<model->K; k++)
		solver->pool.Kt[s*model->K+k] = solver->tmpi[k];
	double dumin = 1/0.0;
	size_t duminj = 0;
	size_t i = solver->unit[l].i;
	for (size_t j=model->Ij[i], J=model->Ij[i+1]; j<J; j++) {
		double du = -model->Jo[j];
		for (size_t x = model->Jx[j], X = model->Jx[j+1];x<X; x++)
			du += solver->tmpx[model->Xk[x]]*model->Xv[x];
		duminj = du < dumin ? j : duminj;
		dumin = du < dumin ? du : dumin;
	}
	solver->pool.l[s] = l;
	solver->pool.j[s] = duminj;
	return s;
}

enum {
	SOLVER_WARMPRIMAL = 0x1,  // primal start point is provided
	SOLVER_WARMDUAL   = 0x2,  // dual start point is provided
	SOLVER_WARMDUALEQ = 0x4,  // use dual equality instead of inequality
	SOLVER_PHASE0     = 0x8,  // run phase 0
	SOLVER_PHASE1     = 0x10, // run phase 1
	SOLVER_PHASE2     = 0x20, // run phase 2
};

static bool solver_solve(Solver *solver, int mode)
{
	FSCGModel *model = solver->model;
	// add slacks if running phase 1
	int slack = model->T + model->Z;
	int nslack = 0;
	if (mode & (SOLVER_PHASE0|SOLVER_PHASE1)) {
		for (int c=1, cx=model->Cnnz[0]; c<=model->C; cx+=model->Cnnz[c++]) {
			// slack only rows that either mention a t-variable or are artificial z-variable definitions
			if (c < model->Ca) {
				for (int x=0; x<model->Cnnz[c]; x++) {
					if (model->CXc[cx+x] < model->T)
						goto slack;
				}
				continue; // no t-variable
			}
slack:
			char s = model->Cs[c];
			double lo, hi;
			if (mode & SOLVER_PHASE0) {
				double du = solver->dual[model->T+c-1];
				double radius = (mode & SOLVER_WARMDUALEQ) ? 0 : 0.05;
				double scale = radius*fabs(du);
				if (scale < radius) scale = radius;
				lo = du - scale;
				hi = -(du + scale);
			} else {
				lo = hi = -1;
			}
			if (s == '>' || s == '=') {
				lp_add_col(&solver->lp, 1, (int[]){model->T+c-1}, (double[]){1}, 0, 1/0.0, lo);
				nslack++;
			}
			if (s == '<' || s == '=') {
				lp_add_col(&solver->lp, 1, (int[]){model->T+c-1}, (double[]){-1}, 0, 1/0.0, hi);
				nslack++;
			}
		}
	}
	solver->A0 = slack + nslack;
	solver->state = SS_PHASE2;
	// if a primal start is provided, skip the first L inits
	if (mode & SOLVER_WARMDUAL) {
		// init up to L warm schedules if a dual start is provided
		// (NOTE: assuming warm primal is schedules 0..n)
		int *warmpri;
		if (mode & SOLVER_WARMPRIMAL) {
			warmpri = calloc(solver->L, sizeof(*warmpri));
			for (int p=0; p<solver->pool.n; p++)
				warmpri[solver->pool.l[p]] = p;
		} else {
			warmpri = NULL;
		}
		for (int l=0; l<solver->L && solver->pool.n<solver->pool.N; l++) {
			int s = solver_initwarm(solver, l);
			// if we also have a warm primal, make sure we didn't just regenerate the warm schedule
			// and waste a slot.
			if (mode & SOLVER_WARMPRIMAL) {
				int p = warmpri[l];
				if (solver->pool.j[p] == solver->pool.j[s] && !memcmp(solver->pool.Kt+p*model->K,
						solver->pool.Kt+s*model->K, model->K*sizeof(*solver->pool.Kt))) {
					// trace("skip wasted slot");
					pool_free(&solver->pool, s);
				}
			}
		}
		free(warmpri);
	}
	int nwarm = solver->pool.n;
	// mix in at most L greedy schedules
	for (int l=0; l<solver->L && solver->pool.n<solver->pool.N; l++)
		solver_initgreedy(solver, l);
	// fill the remaining slots with random schedules
	while (solver->pool.n<solver->pool.N)
		solver_initschd(solver, solver->pool.n%solver->L);
	solver_syncinit(solver, nwarm);
	if (mode & SOLVER_PHASE0) {
		// run phase 0
		solver_syncobj(solver);
		solver_run(solver, 20);
	}
	if (mode & SOLVER_PHASE1) {
		// run phase 1
		if (mode & SOLVER_PHASE0) {
			// reset objective if we ran phase 0
			for (int x=0; x<model->Cnnz[0]; x++) {
				solver->tmpi[x] = model->CXc[x];
				solver->tmpx[x] = 0;
			}
			lp_set_obj(&solver->lp, model->Cnnz[0], solver->tmpi, solver->tmpx);
			for (int a=0, A=solver->a; a<A; a++) {
				solver->tmpi[a] = solver->A0+a;
				solver->tmpx[a] = 0;
			}
			lp_set_obj(&solver->lp, solver->a, solver->tmpi, solver->tmpx);
			for (int s=0; s<nslack; s++) {
				solver->tmpi[s] = slack+s;
				solver->tmpx[s] = -1;
			}
			lp_set_obj(&solver->lp, nslack, solver->tmpi, solver->tmpx);
		}
		solver->state = SS_PHASE1;
		solver_run(solver, -1);
		if (lp_get_obj(&solver->lp) != 0) {
			model_err(model, "infeasible problem");
			return false;
		}
	}
	if (mode & (SOLVER_PHASE0|SOLVER_PHASE1)) {
		// remove slacks
		for (int t=0; t<nslack; t++)
			solver->tmpi[t] = slack+t;
		lp_del_cols(&solver->lp, nslack, solver->tmpi);
	}
	if (mode & SOLVER_PHASE2) {
		// recompute objective
		solver->A0 = slack;
		solver_syncobj(solver);
		// run phase2
		solver->state = SS_PHASE2;
		solver_run(solver, solver->maxgap > 0 ? -2 : -1);
	}
	return true;
}

static void solver_add_solution(Solver *dst, Solver *src, int Iofs)
{
	FSCGModel *model = dst->model;
	for (size_t a=0, A=src->a; a<A; a++) {
		if (src->As[a].it == src->iter)
			pool_transfer(model, &dst->pool, &src->pool, src->As[a].s, Iofs);
	}
	for (size_t l=0, L=src->L; l<L; l++) {
		Unit *uni = &src->unit[l];
		if (uni->ns == 1)
			pool_transfer(model, &dst->pool, &src->pool, uni->sx, Iofs);
	}
}

static bool model_solve_serial(FSCGModel *model)
{
	Solver solver;
	solver_create(&solver, model, model->I, model->s_N0, SOLVER_LOG, model->threads);
	for (int i=0; i<model->I; i++)
		solver.unit[i].i = i;
	int flags;
	switch (model->mode) {
	case 'p': // primal heuristic initialization
		solver_init(&solver);
		flags = SOLVER_WARMPRIMAL | SOLVER_PHASE1 | SOLVER_PHASE2;
		break;
	case 'g': // primal greedy initialization
		flags = SOLVER_PHASE1 | SOLVER_PHASE2;
		break;
	default:
		flags = 0; // silence warning
		assert(0);
	}
	bool ok = solver_solve(&solver, flags);
	solver_destroy(&solver);
	return ok;
}

static void *solver_thread_main(Solver *solver)
{
	solver_init(solver);
	if (solver_solve(solver, SOLVER_WARMPRIMAL | SOLVER_PHASE1 | SOLVER_PHASE2))
		trace("subproblem objective: %.10e solved in %d iterations", solver->obj, solver->iter-1);
	return NULL;
}

static bool model_solve_parallel(FSCGModel *model)
{
	// create solver for master problem
	// for very small problems with primal sampling, we need to ensure that we can fit all schedules
	// from the subproblem.
	int master_N = model->s_N0 ? model->s_N0 : 30*sqrt(model->I);;
	if (model->mode == 'P' && master_N < model->threads*model->C/2) {
		// TODO: this is still a quite optimistic estimate.
		// it works in practice but may fail for pathological cases.
		// the actual proper way to do this would be to compute N after solving the subproblems.
		master_N = model->threads*model->C/2;
		// trace("adjusted master pool size %d -> %d", model->s_N0, master_N);
	}
	Solver solver;
	solver_create(&solver, model, model->I, master_N, SOLVER_LOG, model->threads);
	size_t P = model->threads > 0 ? model->threads : fscg_cpus();
	solver_clock_start(&solver);
	struct Thread {
		Solver solver;
		pthread_t phtread;
	} *threads = malloc(P*sizeof(*threads));
	uint64_t mask = 0;
	// start subproblems
	for (size_t p=0; p<P; p++) {
		struct Thread *thr = &threads[p];
		int I = (model->I+P-p-1)/P;
		solver_create(&thr->solver, model, I, model->s_N1, SOLVER_SERIAL, 0);
		thr->solver.epsilon = DEFAULT_SUBEPSILON;
		thr->solver.maxgap = 1e-7;
		thr->solver.mask = &mask;
		thr->solver.nsub = P;
		thr->solver.no = p;
		for (int i=0; i<I; i++)
			thr->solver.unit[i].i = i*P+p;
		pthread_create(&thr->phtread, NULL, (void *) solver_thread_main, &thr->solver);
	}
	int flags;
	switch (model->mode) {
	case 'D': // dual parallel initialization
		memset(solver.dual, 0, (model->T+model->C)*sizeof(*solver.dual));
		for (int i=0; i<model->I; i++)
			solver.unit[i].i = i;
		for (size_t p=0; p<P; p++) {
			struct Thread *thr = &threads[p];
			// solver_thread_main(&thr->solver);
			pthread_join(thr->phtread, NULL);
			for (size_t d=0, D=model->T+model->C; d<D; d++)
				solver.dual[d] += (1/(double)P) * thr->solver.dual[d];
		}
		flags = SOLVER_WARMDUAL | SOLVER_PHASE0 | SOLVER_PHASE1 | SOLVER_PHASE2;
		break;
	case 'P': // primal parallel initialization
		for (size_t p=0, i=0; p<P; p++) {
			struct Thread *thr = &threads[p];
			pthread_join(thr->phtread, NULL);
			solver_add_solution(&solver, &thr->solver, i);
			for (int ii=0; ii<thr->solver.L; ii++)
				solver.unit[i+ii].i = thr->solver.unit[ii].i;
			i += thr->solver.L;
		}
		flags = SOLVER_WARMPRIMAL | SOLVER_PHASE2;
		break;
	default:
		flags = 0; // silence warning
		assert(0);
	}
	for (size_t p=0; p<P; p++)
		solver_destroy(&threads[p].solver);
	free(threads);
	solver_clock_stop(&solver, SCLOCK_INIT);
	trace("solving master problem");
	bool ok = solver_solve(&solver, flags);
	return ok;
}

static bool fscg_model_solve(FSCGModel *model)
{
	if (!model->epsilon) model->epsilon = DEFAULT_EPSILON;
	if (!model->gamma) model->gamma = DEFAULT_GAMMA;
	switch (model->mode) {
	case 'p': case 'g':
		return model_solve_serial(model);
	case 'P': case 'D':
		return model_solve_parallel(model);
	default:
		model_err(model, "bad mode: %c", model->mode);
		return false;
	}
}

/* ---- Problem reading ----------------------------------------------------- */

struct FSCGFile {
	char *map;
	char *p;
	char *e;
	int fd;
};

static FSCGFile *fscg_file_open(const char *fname)
{
	int fd = open(fname, O_RDONLY);
	if (fd < 0)
		return NULL;
	struct stat sb;
	if (fstat(fd, &sb) < 0)
		goto fail;
	void *map = mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
	if (map == MAP_FAILED)
		goto fail;
	madvise(map, sb.st_size, MADV_SEQUENTIAL);
	FSCGFile *fp = malloc(sizeof(*fp));
	fp->map = map;
	fp->p = map;
	fp->e = map + sb.st_size;
	fp->fd = fd;
	return fp;
fail:
	close(fd);
	return NULL;
}

static const char *fscg_file_line(FSCGFile *fp, size_t *n)
{
	char *p = fp->p;
	char *pp = memchr(p, '\n', fp->e-p);
	if (pp) {
		*n = pp - p;
		fp->p = pp+1;
	} else {
		*n = 0;
	}
	return p;
}

static int fscg_file_read(FSCGFile *fp, int n, double *buf)
{
	char b[33];
	int i = 0;
	char *p = fp->p;
	char *e = fp->e;
	while (i < n) {
		// skip leading space
		while (p < e && isspace(*p)) p++;
		if (p >= e) break;
		// scan number
		char *ee = e-p < 32 ? e : p+32;
		size_t s = 0;
		while (p < ee && !isspace(*p)) b[s++] = *p++;
		p++; // skip space
		// convert number
		b[s] = 0;
		buf[i++] = atof(b);
	}
	fp->p = p;
	return i;
}

static void fscg_file_close(FSCGFile *fp)
{
	munmap(fp->map, fp->e-fp->map);
	close(fp->fd);
}

static int cmpschedule(double *a, double *b, int ncol, int *col2k)
{
	for (int col=0; col<ncol; col++) {
		if (col2k[col] >= 0 && a[col] != b[col])
			return a[col] < b[col] ? -1 : 1;
	}
	return 0;
}

typedef struct SortState {
	double *buf;
	double *obj;
	int *nnz;
	int *col2k;
	int ncol;
} SortState;

static int ss_cmpschedule(const void *ip, const void *jp, void *ssp)
{
	int i = *(int *) ip;
	int j = *(int *) jp;
	SortState *ss = (SortState *) ssp;
	// sort by nnz first, this improves branch prediction when computing duals,
	// and also saves some comparisons here.
	if (ss->nnz[i] != ss->nnz[j]) {
		return (ss->nnz[i] & 3) == (ss->nnz[j] & 3)
			? (ss->nnz[i] - ss->nnz[j])
			: ((ss->nnz[i] & 3) - (ss->nnz[j] & 3));
		// if (ss->nnz[i] == 1) return -1;
		// if (ss->nnz[j] == 1) return 1;
		// // return ss->nnz[i] - ss->nnz[j];
		// // sort by tail. this improves branch prediction for the remainder switch in `scan_minschd`.
		// // `c&-c` isolates the lowest differing bit, the number with that bit set goes first.
		// int c = ss->nnz[i] ^ ss->nnz[j];
		// return (ss->nnz[i] & c & -c) ? -1 : 1;
	} else {
		int c = cmpschedule(ss->buf+ss->ncol*i, ss->buf+ss->ncol*j, ss->ncol, ss->col2k);
		return c ? c : (ss->obj[i] > ss->obj[j] ? -1 : (ss->obj[i] < ss->obj[j] ? 1 : 0));
	}
}

static bool fscg_model_read_xda(FSCGModel *model, const char *fname, int ncol, int *col2k,
	double *col2o)
{
	if (!model->scale) model->scale = 1;
	double scale1 = 1/model->scale;
	FSCGFile *fp = fscg_file_open(fname);
	if (!fp)
		return model_err(model, "failed to read `%s'", fname);
	int sizebuf = 0; // size of schedule buffer
	double *buf = NULL; // schedule rows for current unit
	double *obj = NULL; // schedule objective for current unit
	int *snnz = NULL;   // schedule nonzero count buffer
	int *order = NULL; // schedule order for current
	int J = 0; // number of schedules actually included in the model
	size_t nnz = 0; // number of schedule nonzeros in the model
	// address space is cheap. allocate more than we could possibly need.
	// tcmalloc doesn't let us malloc() this much at once, and if we use realloc() we get killed,
	// so we have to use mmap() here.
	uint8_t *Xk = model->Xk = mmap(NULL, (size_t)model->J*model->K*sizeof(*model->Xk),
		PROT_READ|PROT_WRITE, MAP_ANONYMOUS|MAP_PRIVATE|MAP_NORESERVE, -1, 0);
	double *Xv = model->Xv = mmap(NULL, (size_t)model->J*model->K*sizeof(*model->Xv),
		PROT_READ|PROT_WRITE, MAP_ANONYMOUS|MAP_PRIVATE|MAP_NORESERVE, -1, 0);
	uint32_t *Jx = model->Jx = malloc((model->J+1)*sizeof(*model->Jx));
	double *Jo = model->Jo = malloc(model->J*sizeof(*model->Jo));
	double xmax = 0;
	bool ok;
	for (int i=0; i<model->I; i++) {
		int nj = model->Ij[i];
		model->Ij[i] = J;
		if (nj > sizebuf) {
			sizebuf = nj;
			buf = realloc(buf, sizebuf*ncol*sizeof(*buf));
			obj = realloc(obj, sizebuf*sizeof(*obj));
			snnz = realloc(snnz, sizebuf*sizeof(*snnz));
			order = realloc(order, sizebuf*sizeof(*order));
		}
		int nread = fscg_file_read(fp, nj*ncol, buf);
		if (nread != nj*ncol) {
			model_err(model, "file contains too few values: got %d, expected %d", nread, nj*ncol);
			goto fail;
		}
		// compute objective & nonzeros
		for (int j=0; j<nj; j++) {
			double o = 0;
			int nz = 0;
			for (int col=0; col<ncol; col++) {
				o += col2o[col]*buf[j*ncol+col];
				if (col2k[col] >= 0)
					nz += buf[j*ncol+col] != 0;
			}
			obj[j] = o;
			snnz[j] = nz;
		}
		// sort schedules
		for (int j=0; j<nj; j++)
			order[j] = j;
		qsort_r(order, nj, sizeof(*order), ss_cmpschedule,
			&(SortState){.buf=buf, .obj=obj, .nnz=snnz, .col2k=col2k, .ncol=ncol});
		// add non-duplicates to model
		for (int jj=0; jj<nj; jj++) {
			size_t j = order[jj];
			if (jj == 0 || cmpschedule(buf+j*ncol, buf+order[jj-1]*ncol, ncol, col2k)) {
				Jx[J] = nnz;
				Jo[J] = obj[j]*scale1;
				J++;
				double *row = buf + j*ncol;
				for (int col=0; col<ncol; col++) {
					int k = col2k[col];
					double v = row[col];
					if (k >= 0 && v != 0) {
						Xk[nnz] = k;
						Xv[nnz] = v*scale1;
						xmax = v > xmax ? v : xmax;
						nnz++;
					}
				}
			}
		}
		model->Ij[model->I] = J;
		model->Jx[J] = nnz;
		model->J0 = model->J;
		model->J = J;
	}
	ok = true;
out:
	free(buf);
	free(obj);
	free(snnz);
	free(order);
	fscg_file_close(fp);
	return ok;
fail:
	ok = false;
	goto out;
}

/* ---- Lua API ------------------------------------------------------------- */

#include <lua.h>
#include <lauxlib.h>

static const char FSCG_LUASRC[] = {
#embed "fscg.lua"
	,0
};

static const char FSCG_LUAH[] = {
#embed "fscg.h"
};

static const char FSCG_SOLVER[] =
#ifdef FSCG_GUROBI
"gurobi"
#endif
#ifdef FSCG_HIGHS
"highs"
#endif
;

// order must match fscg.lua!
static const void *FSCG_FTAB[] = {
	malloc,
	fscg_file_open,
	fscg_file_line,
	fscg_file_read,
	fscg_file_close,
	fscg_model_read_xda,
	fscg_model_solve,
	fscg_model_destroy,
	fscg_cpus,
	fscg_time
};

int luaopen_fscg(lua_State *L)
{
	luaL_loadstring(L, FSCG_LUASRC);
	lua_pushlightuserdata(L, FSCG_FTAB);
	lua_pushlstring(L, FSCG_LUAH, sizeof(FSCG_LUAH));
	lua_pushlstring(L, FSCG_SOLVER, sizeof(FSCG_SOLVER)-1);
	lua_call(L, 3, 1);
	return 1;
}
