// Worker script that calculates LCA
// NOTE: When changing this script make sure the browser is not loading old cached version!
importScripts('lca.js'); 

onmessage = function(e) {

  let P = e.data[0];
  let dag = e.data[1];
  let total = P.length;

  let results = P.map( (p,i) => [p[0], p[1], LCA(dag,p[0],p[1])] );
  
  // If third argument exists, pass it along with the results
  // We use this to pass the name of the experiment in the multiple experiment case
  if (typeof e.data[2] == 'undefined') {
      postMessage(results);
  } else {
      postMessage([results,e.data[2]]);
  }


}
