/* 
   Data-parallel sum or factorial implementation; for performance measurements
   This version uses arbitrary precision arithmetic from the GMP library and MPI for parallelism.
   Same as psum_gmp.c (chunking rather than striding) but now uses the 2nd command line argument 
   to specify the chunk size.
   Master performs work as well, at the end of sending out one set of data items.

   Basic GMP API usage is taken 
   from: http://www.cs.colorado.edu/~srirams/classes/doku.php/gmp_usage_tutorial
   Retrieved: 25/1/2013 
   GMP API: http://gmplib.org/manual/

   CPP flags:
    DO_FACT ... perform factorial computation, rather than sum computation (default: sum computation)
    DEBUG ... enable debugging (lots of messages)
    CHECK ... check the result of the parallel computation against seq (or tabular) results
    PRINT ... print the result value

   Compile: mpicc -UDEBUG -DCHECK -UPRINT -DDO_FACT -lm -lgmp -O2 -o psum3_gmp psum3_gmp.c
   Run:     mpirun -n 6 psum3_gmp 1M 100k

   Batch job of measurements + visualisation:
    > for ((i=1;i<8;i++)) ; do echo "PEs: $i">>L0; mpirun -n ${i} psum3_gmp 1M 100k 1>>L0 2>&1 ; done
    > cat L0 | sed -ne '/^PEs/H;/^Elapsed time:/H;${x;p}' | sed -e 's/^PEs: \([0-9]*\).*$/\1/' | sed -e 's/Elapsed time: \([.0-9]*\) secs /\1/' |  sed -e '/[.]/a\X' | sed ':a;N;$!ba;s/\n/ /g' | sed -e 's/X/\n/g'  > rt0.dat
    > echo "set term x11; plot 'rt0.dat' with lines ; pause 10" | gnuplot
*/

#include "gmp.h"
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "string.h"

// size of a string, containing the result of type mpq_t produced by one processor
// Q: Is this in general large enough? If not, how can you fix it?
#define GMP_STR_SIZE 9000000

/* types */
typedef unsigned long int  ui;

/* prototypes */
static inline void sum(mpq_t res, ui m, ui n);
static inline void prod(mpq_t res, ui m, ui n);
static inline void fact(mpq_t res, ui n);
static inline void pow_int(mpq_t res, ui z, ui n) ;
static inline void power_e(mpq_t res, ui z, ui d);
ui readIntK (char *str);

/* variables */
double startTime, stopTime;

/* aux functions */
// res = m+(m+1)+...+n
static inline
void sum(mpq_t res, ui m, ui n){
  ui i;
  mpz_t p, t;
 
  mpz_init_set_ui(p,0); /* p = 0 */
  mpz_init_set_ui(t,0);
  for (i=m; i <= n ; ++i){
    mpz_set(t,p);      /* t = p */ 
    mpz_add_ui(p,t,i); /* p = p + i */
  }
  mpq_set_z(res, p);
  mpz_clear(p);
  mpz_clear(t);
  // return res;
}

/* aux functions */
// res = m*(m+1)*...*n
static inline
void prod(mpq_t res, ui m, ui n){
  ui i;
  mpz_t p, t;
 
  mpz_init_set_ui(p,1); /* p = 0 */
  mpz_init_set_ui(t,1);
  for (i=m; i <= n ; ++i){
    mpz_set(t,p);      /* t = p */ 
    mpz_mul_ui(p,p,i); /* p = p * i */
  }
  mpq_set_z(res, p);
  mpz_clear(p);
  // return res;
}

static inline 
void fact(mpq_t res, ui n) {
  prod(res, 1ul, n);
}  

static inline
ui fact_ui(ui n) { return (n==1) ? 1ul : n*fact_ui(n-1); }

// res = z^n
static inline
void pow_int(mpq_t res, ui z, ui n) {
  mpz_t tmp, zq;
  mpz_init(tmp);
  mpz_init(zq);
  mpz_set_ui(zq,z);
  mpz_pow_ui(tmp, zq, n);
  mpq_set_z(res, tmp);
  mpz_clear(tmp);
  mpz_clear(zq);
  // return res;
}

ui getSize (mpq_t op){ return (ui)mpz_sizeinbase(mpq_numref(op), 10) + (ui)mpz_sizeinbase(mpq_denref(op), 10) + 3ul ; }

// reading large integers, eg. 80k, and store them into a long, if possible
ui readIntK (char *str) {
  ui multiplyBy, maxVal, val;
  int len;
  
  len = (int)strlen(str);
#if defined(DEBUG)
  printf("readIntK: last char of argument 1 '%s' is '%c'\n", str, str[len-1]);
#endif
  switch (str[len-1]) {
  case 'k':
  case 'K': multiplyBy = 1000ul; str[len-1] = '\0'; break;
  case 'm':
  case 'M': multiplyBy = 1000000ul; str[len-1] = '\0'; break;
  case 'g':
  case 'G': multiplyBy = 1000000000ul; str[len-1] = '\0'; break;
  default: multiplyBy = 1ul; break;
  }
  maxVal = (unsigned long int)(-1);
  val = atol(str);
#if defined(DEBUG)
  printf("readIntK: base value is %ul\n", val);
#endif
  if (val > (maxVal / multiplyBy)) {
    // printf(stderr, "Input too large; maximum value: %ul\n", maxVal);
    exit(3);
  } 
  return (val * multiplyBy);
}

int main (int argc, char **argv) {
  int p, p0, id;
  ui m, n, z;
  mpq_t res, res1, tmp1;
  double zz;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);

#if defined(DEBUG)
  fprintf(stderr, "[%d] Init: %d processors\n", id, p);
#endif

  /*
  if (p<=1) {
    fprintf(stderr,"[%d] Need at least 2 processors: 1 master and 1 worker", id); 
    MPI_Abort(MPI_COMM_WORLD, 2);
  }
  */

  if (id == 0) { /*  master   --------------------------------- */

    if (argc<1) {
      fprintf(stderr, "Usage: psum3_gmp <n> <chunksize> ... computes sum over 1 .. n ");
      exit(1);
    }
    if (argc==1) { // just for testing
      n = readIntK(argv[1]); // atol(argv[1]);
      fprintf(stderr, "Usage: psum3_gmp <n> <chunksize> ... computes sum over 1 .. n ");
      exit(1);
      // d = 6; // default precision
    } else {
      n = readIntK(argv[1]);  // input for n!  //  atol(argv[1]);
      z = readIntK(argv[2]);  // chunksize     //  atol(argv[2]);
      //fprintf(stderr, "WARNING: chunk size argument currently unused; it is automatically calculated from number of workers (here: %d) ...\n", n, p-1);
    }
    /* MAIN ------------------- */
#if defined(DO_FACT)
    fprintf(stderr, "Computing factorial of %ul with chunksize %d...\n", n, z);
#else
    fprintf(stderr, "Computing sum [1..%ul] with chunksize %d...\n", n, z);
#endif
    fprintf(stderr, "Using 1 master (also acting as worker) and %d workers ...\n", p-1);

    int i, chunks;
    char* res_str; // was: char res_str[GMP_STR_SIZE]; 
    // double res, res1, res2, res_check;
    mpq_t result, res, tmp_res, tres, resq_check;
    mpz_t resz_check;
    mpf_t resf;
    ui res_check, res_exp;
    double elapsed_time;
    double t1;
    double time;

    mpq_init(result);
    mpq_init(res);
    mpq_init(tmp_res);
    mpq_init(result);
    
    // start the timer
    MPI_Barrier(MPI_COMM_WORLD);
#if defined(DEBUG)
    fprintf(stderr, ".. Master has passed barrier\n");
#endif
    elapsed_time = - MPI_Wtime();

    /* General case: dole out intervals of size t */

    long from, to, last_n, max_n; 
    long send_buf[2];
    int len, l, r;
    double *times;
    MPI_Status status;

    /* compute chunk size z and number of large chunks r */
    if (z == 0l) { // chunksize 0 means compute chunksize based on input an processors (as in psum2_gmp.c)
      z = n/(p-1); r = n%(p-1);
    } else {
      r = p; // r==p means chunksize z has been chosen from commandline => we have more than #workers chunks
    } 

#if defined(DEBUG)
    fprintf(stderr, "[%d] block size %d with %d large chunks, gives in total %ul (expected %ul)\n", id, z, r, r*(z+1)+(p-1-r)*z, n);
#endif

    // initialise variable for overall result
#if defined(DO_FACT)
    mpq_set_ui(result, 1l, 1l);
#else
    mpq_set_ui(result, 0l, 1l);
#endif

    /* init array for timing info per worker */
    times = (double *)malloc(p*sizeof(double));
    memset(times, 0, p*sizeof(double));

    /* Distribute work: send intervals to all workers */
    from=1l;
    to=1l;
    // chunks=0;
    i = 1; 
    while (to<n) { // ((r==p && from<=n) || (r<p && i<p)) { 

      /* send one set of intervals to all workers */
      for (i=1; i<p; i++) { 
	to = (r==p) ? from+z-1 : ((i<=r) ? from+z : from+z-1);
	to = (to>n) ? n : to;
	//// MPI_Send(&z, 1, MPI_LONG, i, 0, MPI_COMM_WORLD);
	//// MPI_Send(&d, 1, MPI_LONG, i, 0, MPI_COMM_WORLD);
	MPI_Send(&from, 1, MPI_LONG, i, 0, MPI_COMM_WORLD);
	MPI_Send(&to, 1, MPI_LONG, i, 0, MPI_COMM_WORLD);
	// chunks++;
	from = to+1;
	if (to==n) { p0 = i+1; break; } // done sending
      }

      /* master computes last interval of the set itself */
      if (to!=n) { 
	to = (r==p) ? from+z-1 : ((i<=r) ? from+z : from+z-1);
	to = (to>n) ? n : to;
#if defined(DEBUG)
	fprintf(stderr, "[%d] MASTER: from=%ul, to=%ul\t", id, from, to);
	// fprintf(stderr, "[%d] computing sum ...\n", id);
#endif
	mpq_init(tres);
	// startTime = clock();
	t1 = - MPI_Wtime();
	prod(tres, from, to);
	t1 += MPI_Wtime();

	times[0] += t1;
#if defined(DO_FACT)
	mpq_set(tmp_res, result);
	mpq_mul(result, tmp_res, tres); // factorial instead of sum
#else
	mpq_set(tmp_res, result);
	mpq_add(result, tmp_res, tres); // do sum; default
#endif
	mpq_clear(tres);
	from = to+1;
      }

      /* Collect result: receive string-encoded GMP values from all workers; any order; we expect #chunks results */

      /* receive one set of results from all workers */
      for (i=1; i<((to==n)?p0:p); i++) { 
	MPI_Recv(&len, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
	if (len==0) { break ; } // done
	// fprintf(stderr, "[%d] len=%d\n", id, len);
	res_str = (char*)malloc((len+1) * sizeof(char));
	if (res_str==NULL) {
	  fprintf(stderr,"[%d] Failed to malloc a buffer of %d char for receiving the worker's results\n", id, GMP_STR_SIZE); 
	  MPI_Abort(MPI_COMM_WORLD, 7);
	}
	MPI_Recv(res_str, len, MPI_CHAR, i, 0, MPI_COMM_WORLD, &status);
	// fprintf(stderr, "[%d] res_str=%s\n", id, res_str);
	MPI_Recv(&time, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
	// fprintf(stderr, "[%d] time=%f\n", id, time);
	/* MPI_Recv(&l1, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &status); */
	/* MPI_Recv(&eps1_str, l1, MPI_CHAR, 1, 0, MPI_COMM_WORLD, &status); */
	res_str[len] = '\0';
#if defined(DEBUG)
	// fprintf(stderr, "[%d] FROM PE %d (as string) res=%s\n", id, i, res_str);
	if (len > GMP_STR_SIZE) {
	  fprintf(stderr,"[%d] Buffer overflow when receiving result from PE %d ", id, i); 
	  MPI_Abort(MPI_COMM_WORLD, 6);
	}
#endif
	/* unmarshall the GMP data */
	if (gmp_sscanf(res_str,"%Qd",&res)==0) {
	  fprintf(stderr,"[%d] Error in gmp_sscanf", id); 
	  MPI_Abort(MPI_COMM_WORLD, 2);
	}
#if defined(DEBUG)
	//gmp_fprintf(stderr, "[%d] FROM PE %d res=%Qd\n", id, i, res);
#endif
	// ress[i] = res;
	times[i] += time;
	// max_n = (last_n>max_n) ? last_n : max_n;
#if defined(DO_FACT)
	mpq_set(tmp_res, result);
	mpq_mul(result, tmp_res, res); // factorial instead of sum
#else
	mpq_set(tmp_res, result);
	mpq_add(result, tmp_res, res); // do sum; default
#endif
	free(res_str);
	// chunks--; // that many chunks left to collect
      } // for

      /*
	if (len==0 && chunks>0) { 
	fprintf(stderr,"[%d] PANIC: left-over chunks after having received finish (len=%d,chunks=%d)\n", len, chunks);
	MPI_Abort(MPI_COMM_WORLD, 8);
	}
      */

    }  // while

    from = 0; to = 0; // end
    for (i=1; i<p; i++) {
      MPI_Send(&from, 1, MPI_LONG, i, 0, MPI_COMM_WORLD);
      MPI_Send(&to, 1, MPI_LONG, i, 0, MPI_COMM_WORLD);
    }

    // stop the timer
    elapsed_time += MPI_Wtime();

    // mpf_set_q(resf,eps2);
#if defined(DEBUG)
    gmp_printf("\n[%d] Finished computation ------------------------------------------\n", id);
    // gmp_printf("\n[%d] Last n was %d (last %d done on master), last epsilon was %.18Ff", id, to2, m, resf);
#endif

#if defined(PRINT)
    mpf_init(resf);
    gmp_fprintf(stderr, "Result = %Qd\t", result);
    mpf_set_q(resf,result);
    // gmp_printf(" = %Ff", resf);
    // mpq_out_str(stdout, 10, result);
    mpf_clear(resf);
#else
    printf("\ndone");
#endif
#if defined(CHECK)
#endif
    printf("\nElapsed time: %f secs \n by PEs: ", 
	   /*   n, res, res_check,*/ elapsed_time);
    for (i=0; i<p; i++) { 
      printf(", PE %d: %f secs", i, times[i]);
    }
    putchar('\n');

#if defined(CHECK)
# if defined(DO_FACT)
    if (n<21ul) { // compute factorial sequentially to check result
      /* check result: do a sequential computation, using integer arithmetic, if the input is small enough */
      mpq_init(resq_check);
      res_exp = fact_ui(n);
      mpq_set_ui(resq_check, res_exp, 1);
#  if defined(PRINT)
      gmp_fprintf(stderr, "expected result of factorial %ul = %Qd\n", n, resq_check);
#  endif
      if (mpq_cmp(result, resq_check)==0) {
        printf("++ Result OK\n");
      } else {
        printf("** Result WRONG\n");
      }
      mpq_clear(resq_check);
    } else if (n<50ul) { // compute factorial sequentially to check result
      /* check result: do a sequential computation, using GMP */
      mpq_init(resq_check);
      // startTime = clock();
      elapsed_time = - MPI_Wtime();
      prod(resq_check, 1, n);
      elapsed_time += MPI_Wtime();
      // stopTime = clock();
#  if defined(PRINT)
      gmp_fprintf(stderr, "expected result of factorial %ul = %Qd\n", n, resq_check);
#  endif
      if (mpq_cmp(result, resq_check)==0) {
	printf("++ Result OK\n");
      } else {
	printf("** Result WRONG\n");
      }
      printf("\nSequential time: %f secs \n", elapsed_time);
      mpq_clear(resq_check);
    } else if (n==50ul) { // table lookup to check result
      mpq_init(resq_check);
      mpq_set_str(resq_check, "50", 10);
#  if defined(PRINT)
      gmp_fprintf(stderr, "expected result of factorial %ul = %Qd\n", n, resq_check);
#  endif
      if (mpq_cmp(result, resq_check)==0) {
         printf("++ Result OK\n");
      } else {
         printf("** Result WRONG\n");
      }
       mpq_clear(resq_check);
    } else if (n==1000ul) { // table lookup to check result
      mpq_init(resq_check);
      mpq_set_str(resq_check, "1000", 10);
#  if defined(PRINT)
      gmp_fprintf(stderr, "expected result of factorial %ul = %Qd\n", n, resq_check);
#  endif
      if (mpq_cmp(result, resq_check)==0) {
         printf("++ Result OK\n");
      } else {
         printf("** Result WRONG\n");
      }
       mpq_clear(resq_check);
    } else {
      printf("\nFYI: Using mpz_fac_ui to check result for input larger than 50\n");
      elapsed_time = - MPI_Wtime();
      mpz_init(resz_check);
      mpq_init(resq_check);
      mpz_fac_ui(resz_check, n);
      mpq_set_z(resq_check, resz_check);
      elapsed_time += MPI_Wtime();
      if (mpq_cmp(result, resq_check)==0) {
	printf("++ Result OK\n");
      } else {
	printf("** Result WRONG\n");
      }
      printf("\nmpz_fac_ui sequential time: %f secs \n ", elapsed_time);
      mpq_clear(resq_check);
      mpz_clear(resz_check);
    }
# else
    /* check result: do a sequential computation */
    mpq_init(resq_check);
    // startTime = clock();
    elapsed_time = - MPI_Wtime();
    sum(resq_check, 1, n);
    elapsed_time += MPI_Wtime();
    // stopTime = clock();
    if (mpq_cmp(result, resq_check)==0) {
      printf("++ Result OK\n");
    } else {
      printf("** Result WRONG\n");
    }
    printf("\nSequential time: %f secs \n ", elapsed_time);

    //gmp_printf("\ne^%d = %Qd", n, res);
    //mpq_out_str(stdout, 10, res);
    res_exp = (n % 2) ? ((n+1ul)/2)*n : (n/2)*(n+1); // beware of overflows!
    mpq_init(resq_check);
    mpq_set_ui(resq_check, res_exp, 1);
    // if (n<=1000000000ul) { gmp_fprintf(stdout, "expected result of sum [1..%ul] = %Qd\n", n, resq_check); } 
    if (mpq_cmp(result, resq_check)==0) {
      printf("++ Result OK\n");
    } else {
      printf("** Result WRONG\n");
    }
    mpq_clear(resq_check);
# endif

    // mpf_set_q(resf,res);
    // gmp_printf("\ne^%d = %Ff", n, resf);
    // printf("\n%d! = ", n);
    // mpf_out_str(stdout, 10, d, resf);
    //printf("\nElapsed time: %f secs\n", 
    //   (stopTime-startTime)/CLOCKS_PER_SEC);
#endif

  } else { /*  worker   ----------------------------------------------------------------------------------- */
    MPI_Status status;
    int len, len0, max_n;
    ui from, to, sz;
    double elapsed_time;
    mpq_t res, eps;
    // was: char buf[GMP_STR_SIZE], buf0[GMP_STR_SIZE]; // FIXME: EVIL CONSTANT
    char *buf, *buf0; // FIXME: EVIL CONSTANT

    /*
    buf = (char*)malloc(GMP_STR_SIZE * sizeof(char));
    buf0 = (char*)malloc(GMP_STR_SIZE * sizeof(char));
    if (buf==NULL || buf0==NULL) {
      fprintf(stderr,"[%d] Failed to malloc a buffer of %d char for sending the worker's results\n", id, GMP_STR_SIZE); 
      MPI_Abort(MPI_COMM_WORLD, 7);
    }
    */

    mpq_init(res);
    mpq_init(eps);

    // start the timer
    MPI_Barrier(MPI_COMM_WORLD);
    // elapsed_time = - MPI_Wtime();

#if defined(DEBUG)
    fprintf(stderr, ".. Worker %d has passed barrier\n", id);
#endif

    while (1) {
      /* Receive interval to calculate */
      //MPI_Recv(&z, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD, &status);
      //MPI_Recv(&d, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(&from, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(&to, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD, &status);

      if (from==0ul && to==0ul) { break; }
      
      // start the timer
      elapsed_time = - MPI_Wtime();

#if defined(DEBUG)
      fprintf(stderr, "[%d] from=%ul, to=%ul\t", id, from, to);
      // fprintf(stderr, "[%d] computing sum ...\n", id);
#endif
      /* Calculate own result */
#if defined(DO_FACT)
      prod(res, from, to);           // result of our interval
#else
      sum(res, from, to);            // result of our interval
#endif
      sz = getSize(res);
      buf = (char*)malloc((sz+1) * sizeof(char));
#if defined(DEBUG)
      // get size:  mpz_sizeinbase (mpq_numref(op), base) + mpz_sizeinbase (mpq_denref(op), base) + 3
      if (sz > GMP_STR_SIZE) {
	fprintf(stderr,"[%d] PANIC: size of result %ul is higher than GMP_STR_SIZE %ul\n", id, sz, GMP_STR_SIZE); 
	MPI_Abort(MPI_COMM_WORLD, 7);
      }
      // gmp_fprintf(stderr, "[%d] res=%Qd", id, res);
      //fprintf(stderr, "[%d] computing sum ...\n", id);
#endif
      if ((len = gmp_sprintf(buf,"%Qd", res))==0 || len>=GMP_STR_SIZE) {      // marshall into a string
	fprintf(stderr,"[%d] Error in gmp_sprintf", id); 
	MPI_Abort(MPI_COMM_WORLD, 2);
      }
#if defined(CHECK)
      fprintf(stderr, "[%d] sending result of size %d (GMP_STR_SIZE %d)\n", id, len, GMP_STR_SIZE);
      // fprintf(stderr, "[%d] last_n=%d, res=%s\n", id, max_n, buf);
      // fprintf(stderr, "[%d] last_n=%d\n", id, max_n);
#endif

      // stop the timer
      elapsed_time += MPI_Wtime();

      /* Send own result (and timing info) */
      MPI_Send(&len, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Send(buf, len, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
      // MPI_Send(&max_n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Send(&elapsed_time, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
      /* MPI_Send(&len0, 1, MPI_INT, 0, 0, MPI_COMM_WORLD); */
      /* MPI_Send(buf0, len0, MPI_CHAR, 0, 0, MPI_COMM_WORLD); */

    } // while

    /*
    if (!(from==0ul && len==0ul)) {
	fprintf(stderr,"[%d] PANIC: worker finished before receiving FINISH token", id); 
	MPI_Abort(MPI_COMM_WORLD, 9);
    }
    */

    //MPI_Send(&len, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

    free(buf);
  } /* end worker */

  MPI_Finalize();
  //return 0;
  exit(0);
}
 
