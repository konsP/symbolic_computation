use GMP;
use CyclicDist, BlockCycDist;
use Time;
var t: Timer;

config debug =false;
config const n : int = 10;
config const n4 : c_ulong = 10; 
config const z : int = 2;


var result = new BigInt("1"); 
var tres = new BigInt("1"); 

const Space={1..n};
const BCSpace= Space dmapped BlockCyclic(startIdx=1, blocksize=z);

const LocaleSpace = {0..numLocales-1};
const res : [LocaleSpace] BigInt; 



for loc in Locales do
  res[loc.id] =  new BigInt("1"); 

t.start();
forall i in BCSpace {
	res[here.id].mul_si(res[here.id],i);
}
for loc in Locales do
    result.mul(result, res[loc.id]);
t.stop();
writeln("Elapsed time: ", t.elapsed());



writeln("===================================================");
if debug then writeln("Using mpz_fac_ui to compute factorial ...");
var t0: Timer;
var test = new BigInt("1"); // global variable, stores result; requires sync
t0.start();
test.fac_ui(n4);
t0.stop();
if (result.cmp(test)==0) { 
  writeln("++ OK");
} else {
  writeln("** WRONG");
}
if debug then writeln("fac_ui result = ", test );
if debug then writeln("fac_ui elapsed time: ", t0.elapsed());
delete test;
delete result;


/*
 use GMP;
 var a = new BigInt(); 	// initialize a GMP value, set it to zero
 a.fac_ui(100); 	// set a to 100!
 writeln(a); 		// output 100!
 delete a; 		// free memory used by the GMP value
 var b = new BigInt("48473822929893829847"); // initialize from a decimal string
 b.add_ui(b, 1); 	// add one to b
 delete b; 		// free memory used by b
*/
