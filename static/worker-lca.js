// Worker script that calculates LCA
// NOTE: When changing this script make sure the browser is not loading old cached version!
importScripts('lca.js'); 

onmessage = function(e) {

  let P = e.data[0];
  let dag = e.data[1];
  let total = P.length;

  let results = P.map( (p,i) => [p[0], p[1], LCA(dag,p[0],p[1])] );
  
  postMessage(results);
}
