/* factorial5.chpl
*  Low level SPMD-style implementation using BLockCyclic distribution of elements	
*  Chunksize z is passed as argument
*  Calculates range for each worker
*
*  USE:
*  under directory: chapel-1.8.0:
*  on one node:
*  export CHPL_COMM=gasnet;gmake;export GASNET_SPAWNFN=L;export GASNET_ROUTE_OUTPUT=0;source util/setchplenv.bash
*   or:
*  on multiple nodes
*  export CHPL_COMM=gasnet;gmake;export GASNET_SPAWNFN=S;export SSH_SERVERS="bwlf01 bwlf02 .. bwlf32";
*  export SSH_CMD=ssh;export SSH_OPTIONS=-x;export GASNET_ROUTE_OUTPUT=0;source util/setchplenv.bash;
*
*  Compile: chpl --fast -o factorial5 factorial5.chpl
*  Run: ./factorial5 -nl 32 --n=17000 --n4=17000 --z=3000
*/
use GMP;
use CyclicDist, BlockCycDist;
use Time;
var t: Timer;


config const n : int = 10; 		//upper bound of input range
config const n4 : c_ulong = 10; 	// -//- for result validation
config const z : int = 2;		//chunk size
config debug = false;
var result = new BigInt("1"); 		//final result

const LocaleSpace = {0..numLocales-1};
const res : [LocaleSpace] BigInt; 

const Space={1..n};
const BCSpace= Space dmapped BlockCyclic(startIdx=1, blocksize=z); 	//BlockCyclic distribution of elements


proc prod0 (from: int, to: int): BigInt { 				//SPMD part of the calculation -- calculates local result
  var res = new BigInt("1");                                        
  var i: int;
  for i in from..to  {
    res.mul_si(res, i);
  }
  return res;
}

proc prodS (from: int, to: int): BigInt { 				//SPMD part of the calculation -- calculates local result
  var syncr$ :sync bool = true;				
  var res = new BigInt("1");                                        
  var i: int;
  forall i in from..to  {
    syncr$;
    res.mul_si(res, i);
    syncr$=true; 
  }
  return res;
}


for loc in Locales do
	res[loc.id] =  new BigInt("1"); 			

t.start(); 							//start timer



coforall locid in BCSpace.dist._value.targetLocDom do		//iteration over locales
  on BCSpace.dist._value.targetLocales[locid] {
 
	
	var myStarts =BCSpace._value.locDoms[locid].myStarts;	//returns a strided range with the starting index of each block 
								//eg. n = 10, z = 2, numLocales = 2
								//locale 0: myStarts= {1..10 by 4}
								//locale 1: myStarts= {3..10 by 4}
	var syncronize$ :sync bool = true;

	forall i in myStarts {
		var j=i+(z-1);					//calculate last index: start + (chunksize -1)	
		if(j>n) then j=n;				//last index in range else, go till n
			var tres = new BigInt("1"); 
			tres =  prod0(i,j);
			//tres = prodS(i, j);
			syncronize$;				//acquire lock to write in the array
     			res[locid].mul(res[locid],tres);
			syncronize$=true;			//release lock
		
		
	}
    }

for loc in Locales do						//reduce to final result
    result.mul(result, res[loc.id]);


t.stop();							//stop timer
writeln( t.elapsed());


var test = new BigInt("1"); 					//check result
test.fac_ui(n4);
if (result.cmp(test)==0) { 
  if debug then writeln("++ OK");
} else {
  if debug then writeln("** WRONG");
}


delete test;	
						//free memory
delete result;


