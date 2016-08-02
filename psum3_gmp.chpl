/*
  Parallel factorial computation, matching the C+MPI code in psum3_gmp.c:
  naive factorial, as a fold over mult, done by n workers in strides of n
  Same as psum_gmp.chpl (chunking rather than striding) but now uses the z 
  config var to specify the chunk size.
  Master performs work as well, same as in the C+MPI version.

  Setup needs to be run in top-level chapel dir. 
  CHPL_COMM=gasnet for multi-node exec

  Setup:   source util/setchplenv.bash 
  Compile: chpl -O -o psum3_gmp psum3_gmp.chpl
  Run:     ./psum3_gmp --numLocales=2 --n=100000 --n4=100000 --z=50000
*/

use GMP;

use Time;
//
var t: Timer;

/* ----------------------------------------------------------------------------- */
/* factorial, builtin */
/*
var a = new BigInt(); // initialize a GMP value, set it to zero
a.fac_ui(100); // set a to 100!
writeln(a); // output 100!
delete a; // free memory used by the GMP value
*/
/* var b = new BigInt("48473822929893829847"); // initialize from a decimal string */
/* b.add_ui(b, 1); // add one to b */
/* delete b; // free memory used by b */

/* ----------------------------------------------------------------------------- */
/* factorial, as a loop over multiplication */
/*
var i: c_ulong;
var n0: c_ulong = 200;
var x = new BigInt("1");
for i in 1..n0 do
	x.mul_ui(x, i);

writeln("Factorial of "); writeln(n0); writeln(x);

var x0 = new BigInt(); // initialize a GMP value, set it to zero
x0.fac_ui(n0);          // set a to 100!
writeln("Chapel says: Factorial of "); writeln(n0); writeln(x0);
if (x.cmp(x0)==0)
then 
  writeln("++ Result OK");
else 
  writeln("** Result WRONG");
delete x0; // free memory used by the GMP value
delete x; // free memory used by the GMP value
*/

/* ----------------------------------------------------------------------------- */
/* factorial, in parallel over intervals; same structure as in the C+MPI version */

/*
var ProblemSpace = all;
var RgLocales = {0..#numLocales};
var AllLocales: [RgLocales] locale = reshape(Locales, RgLocales);
*/
// var D : domain(1) dmapped Cyclic(1) = all;

proc prod (from: int, to: int, stride: int): BigInt {
  var res = new BigInt("1");                                           // local variable for local result
  var i: int;
  for i in from..to by stride {
    res.mul_si(res, i);
  }
  return res;
}

proc prod0 (from: int, to: int): BigInt {
  var res = new BigInt("1");                                           // local variable for local result
  var i: int;
  for i in from..to  {
    res.mul_si(res, i);
  }
  return res;
}

const LocaleSpace = {0..numLocales-1};
const res : [LocaleSpace] BigInt;           // local result

var result = new BigInt("1"); // global variable for result; needs to be synchronised
var tres = new BigInt("1"); // tmp result
var done : bool = false;
var from : int = 1;
config const n : int = 200;
config const n4 : c_ulong = 200; // TODO: use variable 'n' and cast it to c_ulong
config const z : int = 20;

writeln("Computing factorial of ", n, " in chunks of ", z, " on ", numLocales, " locales");

for loc in Locales do
  res[loc.id] =  new BigInt("1"); 

t.start();
done=false;
from=1;
while (!done) {
  // TODO: handle the case of z*numLocales doesn't divide n
  for loc in Locales do
    on loc {
      tres =  prod(from+loc.id*z,from+loc.id*z+z-1,numLocales);
      res[loc.id].mul(res[loc.id],tres);
    }
  
  from = from+numLocales*z;
  if (from > n) {
    done = true;
  }
 }

for loc in Locales do
    result.mul(result, res[loc.id]);

t.stop();
// writeln("Result of a 2 Locale, parallel factorial of "); writeln(n); writeln(result);
writeln("Elapsed time: ", t.elapsed());

writeln("Computing factorial (sequential in strides of ", numLocales, ")  ...");
// var t0: Timer;
var tmp = new BigInt("1"); // global variable for result; needs to be synchronised
var result4 = new BigInt("1"); // global variable for result; needs to be synchronised
t.start();
for i in 1..numLocales {
  tmp = prod(i,n,numLocales);
  result4.mul(result4, tmp);
}
t.stop();
writeln("Sequential time: ", t.elapsed());

writeln("Using mpz_fac_ui to compute factorial ...");
var t0: Timer;
var test = new BigInt("1"); // global variable for result; needs to be synchronised
t0.start();
test.fac_ui(n4);
t0.stop();
if (result.cmp(test)==0) { 
  writeln("++ OK");
} else {
  writeln("** WRONG");
}
writeln("fac_ui elapsed time: ", t0.elapsed());
delete test;
delete result;

/*
const n4: c_ulong = 20;
var result4 = new BigInt("1"); // global variable for result; needs to be synchronised
var test = new BigInt("1"); // global variable for result; needs to be synchronised
test.fac_ui(n4);
result4 = prod(1,n3,1);
writeln("Testing prod over 1..", n4, " stride 1: ", result4);
writeln("same as fact ", n4, "? :", test);
if (test.cmp(result4)==0) {
  writeln("++ OK");
} else {
  writeln("** WRONG");
}
delete result4;
*/
/*
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
*/

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

