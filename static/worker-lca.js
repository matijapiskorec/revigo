// Worker script that calculates LCA
// NOTE: When changing this script make sure the browser is not loading old cached version!
importScripts('lca.js'); 

onmessage = function(e) {
  console.log('Data received from main script');
  // console.log(e);
  let P = e.data[0];
  let dag = e.data[1];
  let total = P.length;

  let results = P.map( (p,i) => [p[0], p[1], LCA(dag,p[0],p[1])] );
 
  // TODO: For now just reporting progress!
  // let results = [];
  // for (const [i,p] of P.entries()) {
  	  // results.push([p[0], p[1], LCA(dag,p[0],p[1])]);
          // if (i%1000==0) { postMessage(100*(i/total).toPrecision(3)); }
	  // // if (i%1000==0) {console.log(i);}
  // }
  
  console.log('Posting message back to main script');
  postMessage(results);
}
