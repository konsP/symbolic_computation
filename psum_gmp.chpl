/*
  Parallel factorial computation, matching the C+MPI code in psum_gmp.c

  Setup needs to be run in top-level chapel dir. 
  CHPL_COMM=gasnet for multi-node exec

  Setup:   source util/setchplenv.bash 
  Compile: chpl -O -o psum_gmp psum_gmp.chpl
  Run:     ./psum_gmp 
*/

use GMP;

/* factorial, builtin */
var a = new BigInt(); // initialize a GMP value, set it to zero
a.fac_ui(100); // set a to 100!
writeln(a); // output 100!
delete a; // free memory used by the GMP value

/* var b = new BigInt("48473822929893829847"); // initialize from a decimal string */
/* b.add_ui(b, 1); // add one to b */
/* delete b; // free memory used by b */

/* ----------------------------------------------------------------------------- */
/* factorial, as a loop over multiplication */
var i: c_ulong;
var n: c_ulong = 200;
var x = new BigInt("1");
for i in 1..n do
	x.mul_ui(x, i);

writeln("Factorial of "); writeln(n); writeln(x);

var x0 = new BigInt(); // initialize a GMP value, set it to zero
x0.fac_ui(n);          // set a to 100!
writeln("Chapel says: Factorial of "); writeln(n); writeln(x0);
if (x==x0)
then 
  writeln("++ Result OK");
else 
  writeln("** Result WRONG");
delete x0; // free memory used by the GMP value
delete x; // free memory used by the GMP value


/* ----------------------------------------------------------------------------- */
/* factorial, in parallel over intervals; same structure as in the C+MPI version */

/*
var ProblemSpace = all;
var RgLocales = {0..#numLocales};
var AllLocales: [RgLocales] locale = reshape(Locales, RgLocales);
*/
// var D : domain(1) dmapped Cyclic(1) = all;

const LocaleSpace = {0..numLocales-1};
//const res : [LocaleSpace] BigInt;           // local result

var result = new BigInt("1"); // global variable for result; needs to be synchronised
//var p : int;
const nn : int = 20;
var p,z,r : int;
/* compute chunk size z and number of large chunks r */
p = 5; // numLocales; // number of processors/locales to use
z = nn/(p-1); 
r = nn%(p-1);
// writeln(stderr, "[%d] block size %d with %d large chunks, gives in total %ul (expected %ul)\n", id, z, r, r*(z+1)+(p-1-r)*z, n);

/* Distribute work: send intervals to all workers */
forall i in 1..p-1 { // or: LocaleSpace { // WAS: for (i=1; i<p; i++) {   // should be done in parallel
  var res = new BigInt("1");                                           // local variable for local result
  var from: int = if (i<=r) then (i-1)*(z+1) else r*(z+1)+(i-r-1)*z;  // lower bound for our interval
  var to: int = if (i<=r) then from+z else from+z-1;                  // upper bound for our interval
  for val in from..to do  // should be done sequentially
	res.mul_si(res, val);
  result.mul(result,res); // multiply global with local result; needs to be synchronised
}
// or: result = BigInt.mul reduce res;
// writeln("[%d] Master now waiting for results from %d workers ...\n", id, p-1);
writeln("Result of parallel factorial of "); writeln(nn); writeln(result);
delete result;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* pure GMP code: */

/* aux functions */
// res = m*(m+1)*...*n
/* static inline */
/* void prod(mpq_t res, ui m, ui n){ */
/*   ui i; */
/*   mpz_t p, t; */
 
/*   mpz_init_set_ui(p,1); /\* p = 0 *\/ */
/*   mpz_init_set_ui(t,1); */
/*   for (i=m; i <= n ; ++i){ */
/*     mpz_set(t,p);      /\* t = p *\/  */
/*     mpz_mul_ui(p,p,i); /\* p = p * i *\/ */
/*   } */
/*   mpq_set_z(res, p); */
/*   mpz_clear(p); */
/*   // return res; */
/* } */

