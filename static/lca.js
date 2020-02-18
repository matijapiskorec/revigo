
// Potentially faster version of LCA, without calculating full ancestor paths
function LCA(dag,u1,u2) {
  
  if (u1==u2) {return u1;}
  
  let A1 = new Map(); A1.set(u1,0);
  let A2 = new Map(); A2.set(u2,0);
  let depth = 1;
  let lcaCand = new Map();
  let curNodes1 = new Array([u1]);
  let curNodes2 = new Array([u2]);
  
  while ( (lcaCand.size==0) && (curNodes1.length!=0 || curNodes2.length!=0) ) {
  
    let nextNodes1 = new Array();
    let nextNodes2 = new Array();

    for (const v1 of curNodes1) {
      // for (const p1 of dag[v1].parents) {
      for (const p1 of dag[v1]) {
        if (A2.has(p1)) {
          lcaCand.set(p1,depth);
        } else {
          nextNodes1.push(p1);
          A1.set(p1,depth);
        }
      }
    }
    curNodes1 = nextNodes1;
    nextNodes1 = new Array();
    
    for (const v2 of curNodes2) {
      // for (const p2 of dag[v2].parents) {
      for (const p2 of dag[v2]) {
        if (A1.has(p2)) {
          lcaCand.set(p2,depth);
        } else {
          nextNodes2.push(p2);
          A2.set(p2,depth);
        }
      }
    }
    curNodes2 = nextNodes2;
    nextNodes2 = new Array();
    
    depth = depth + 1;

  }
  // return Array.from(lcaCand).reduce( (res,x) => x[1] < res ? x[1] : res )[0];
  return Array.from(lcaCand).reduce( (res,x) => x[1] < res[1] ? x : res, ["",Infinity])[0];
}
