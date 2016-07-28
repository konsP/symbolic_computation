/*
Using a Cyclic distribution of elements.
Calculating from..to range doesnt require a chunk, 
it is done automatically.
Cyclic gives a strided domain by stride = numLocales so no index out of bound exceptions.

under CHPL_HOME
export CHPL_COMM=gasnet;gmake;export GASNET_SPAWNFN=L;export GASNET_ROUTE_OUTPUT=0;source util/setchplenv.bash

Compile: chpl --fast -o factorial3 factorial3.chpl
Run: ./factorial3 -nl 5 --n=100 --n4=100

*/
use GMP;
use Time;
var t: Timer;
use CyclicDist;

var result = new BigInt("1"); 	// global variable for result; needs to be synchronised
var tres = new BigInt("1"); 	// tmp result
var done : bool = false;
var from : int = 1;
config const n : int = 200;
config const n4 : c_ulong = 200; 	// TODO: use variable 'n' and cast it to c_ulong
config const z : int = 20;
config debug = false;


proc prod (from: int, to: int): BigInt {
  var res = new BigInt("1");       	// local variable for local result
  var i: int;
  for i in from..to {
    res.mul_si(res, i);
  }
  return res;
}



const LocaleSpace = {0..numLocales-1};
const res : [LocaleSpace] BigInt;           // local result

const Space = {1..n by numLocales};
const CyclicSpace= Space dmapped Cyclic(startIdx=1);


for loc in Locales do
  res[loc.id] =  new BigInt("1"); 

t.start();
done=false;
from=1;

while (!done) {
  for loc in Locales do
    on loc {
	
      var myFrom=from+loc.id*numLocales;
      var myTo= from+loc.id*numLocales+numLocales-1;
      tres =  prod(myFrom,myTo);
  
      if (myTo <= n) then
      	res[loc.id].mul(res[loc.id],tres);
  }
  from = from+(numLocales-1)*numLocales+numLocales;
  if (from >= n) {
    done = true;
    if debug then writeln("done");
  }
 }

for loc in Locales do
    result.mul(result, res[loc.id]);

t.stop();
writeln("Elapsed time: ", t.elapsed());
if debug then writeln("result ", result);


if debug then writeln("Using mpz_fac_ui to compute factorial ...");

var test = new BigInt("1"); // global variable for result; needs to be synchronised
test.fac_ui(n4);
if (result.cmp(test)==0) { 
	if debug then  writeln("++ OK");
} else {
  if debug then  writeln("** WRONG");
}
if debug then writeln("fac ui :", test);
delete test;
delete result;



