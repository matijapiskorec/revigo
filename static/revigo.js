
// var app = new Vue({
  // el: '#app',
  // data: {
    // message: 'Hello Vue!'
  // }
// });


// Generate all pairs of elements of given array, without repetition and without pairs with identical elements
function generate_pairs(N) {
  let result = [];
  let i = -1;
  for (const a of N) {
    i = i + 1;
    for (const b of N.slice(i+1,)) {
       result.push([a,b]);
    }
  }
  return result;
}

function getUrl(url) {
return fetch(url)
  .then((response) => {
    return response.json();
  })
};

// Prepare canvas
let canvas;

// Worker which will calculate LCA in a separate thread
// NOTE: When updating this file make sure your browser does not load old cached version!
var lcaWorker = new Worker('static/worker-lca.js');

// Wait for all data to load
Promise.all(["data/go-dag-molecular-function.json",
             "data/go-dag-cellular-component.json",
	     "data/go-dag-biological-process.json",
	     "data/go-terms-logcount-goa.json",
	     "data/revigo-enrichments1.json"].map(url=>getUrl(url)))
  .then(([dagMol,dagCell,dagBio,terms,enrichments]) => {
	  
  let dag = dagBio; // dagBio, dagCell, dagMol
  let ontologyName = "biological processes"; // biological_processes, cellular component, molecular function

  // Selecting only enrichements from a specific ontology namespace
  console.log("Selecting only enrichements from "+ontologyName+" namespace.");
  let enrichmentsSelected = Object.entries(enrichments).filter( x => dag.hasOwnProperty(x[0]) );
  console.log("Using "+enrichmentsSelected.length+" out of "+Object.entries(enrichments).length+" enrichments");

  // We are only interested in small p-values
  enrichmentsSelected = enrichmentsSelected.filter( x => x[1] < 0.01);
  console.log("Using "+enrichmentsSelected.length+" enrichments with p-value < 0.01.");
  let P = generate_pairs(enrichmentsSelected.map(x => x[0]));
  let termsToId = new Map(enrichmentsSelected.map( (x,i) => [x[0],i] ));

  // Calculation of LCA in the worker
  let t0 = performance.now();
  lcaWorker.postMessage([P,dag]);
  console.log('Sent pairs and dag to worker!');
  $("#progressBox").html('Calculating semantic similarities between '+enrichmentsSelected.length+' selected '+ontologyName+' GO terms with p<0.01 <img src="static/ajax-loader.gif">');

  // Wait for the worker to return LCA results
  lcaWorker.onmessage = function(e) {
    let t1 = performance.now();
    $("#progressBox").html("Calculated semantic similaritites between "+enrichmentsSelected.length+" selected "+ontologyName+" GO terms with p<0.01 GO terms in "+(t1-t0)+" miliseconds!");
    console.log("Calculated "+P.length+" LCA's for all pairs of "+enrichmentsSelected.length+" GO terms in " + (t1-t0) + " miliseconds (on the worker!).");
    let results = e.data; 
    let similarity = results.map( (x) => [x[0], x[1], terms[x[2]]]);

    // Select only pairs which have actual values
    similarity = similarity.filter( x => x[2] );
    console.log("Selecting "+similarity.length+" GO term pairs with defined LCA's (out of "+P.length+" pairs)");
    s = similarity;

    // List of GO terms for visualization
    let goTerms = similarity.map(x=>[x[0],x[1]]).flat().unique()
    console.log(goTerms);

  // Default max value which we use for t-sne distance matrix - 10x the max distance in data
  console.log("Creating a "+goTerms.length+"x"+goTerms.length+" distance (dissimilarity) matrix.");
  t0 = performance.now();
  let maxValue = Math.max.apply(Math, s.map(x=>x[2]));
  let distMat = [...Array(goTerms.length)].map(e => Array(goTerms.length).fill(maxValue));
  for (const x of similarity) {
	  let x0 = termsToId.get(x[0]);
	  let x1 = termsToId.get(x[1]);
	  distMat[x0][x1] = x[2];
	  distMat[x1][x0] = x[2];
  }
  t1 = performance.now();
  console.log(distMat);
  console.log("Distance matrix created in " + (t1-t0) + " miliseconds.");

  // Prepare canvas
  // canvas = document.querySelector("canvas"),
  canvas = $("#tsne-plot")[0],
  context = canvas.getContext("2d"),
  width = canvas.width,
  height = canvas.height;

  console.log("Starting t-SNE embedding...");

  // Calculating t-SNE
  var opt = {}
  opt.epsilon = 10; // epsilon is learning rate (10 = default)
  opt.perplexity = 20; // roughly how many neighbors each point influences (30 = default)
  opt.dim = 2; // dimensionality of the embedding (2 = default)
  
  var tsne = new tsnejs.tSNE(opt); // create a tSNE instance
  tsne.initDataDist(distMat);

  // Initial iterations before we start dynamic visualization
  for(var k = 0; k < 10; k++) {
     tsne.step(); // every time you call this, solution gets better
  }

  // Returns control to the browser so that canvas can be redrawn
  // Time intrval is set to 0, with no delay between redrawing
  setInterval(function() {
      for(var k = 0; k < 1; k++) {
         tsne.step(); // every time you call this, solution gets better
      }
      Y = tsne.getSolution();
      Y0min = Math.min(...Y.map(x=>x[0]));
      Y0max = Math.max(...Y.map(x=>x[0]));
      Y1min = Math.min(...Y.map(x=>x[1]));
      Y1max = Math.max(...Y.map(x=>x[1]));
      context.clearRect(0, 0, canvas.width, canvas.height);
      Y.forEach((y,i)=>drawNode(y,Y0min,Y0max,Y1min,Y1max,
	                        canvas.width,canvas.height,
	                        enrichmentsSelected[i][0].substring(3)));
  },0);

  }

function drawNode(d,Y0min,Y0max,Y1min,Y1max,w,h,label) {
  let ds = scaleToCanvas(d,Y0min,Y0max,Y1min,Y1max,w,h);
  context.beginPath();
  context.fillStyle = "#9ecae1";
  context.strokeStyle = "#000";
  context.arc(ds[0],ds[1],3,0,2*Math.PI);
  context.fill();
  context.stroke();
  context.font = "9px Helvetica";
  context.fillStyle = "#000";
  context.fillText(label,ds[0]+5,ds[1]+3);
}

// [min,max] -> [a,b]
// [min,max] -> [0,canvas.width]
// f(x) = (b-a)*(x-min) / (max-min) + a
// f(x) = canvas.width*x / (max-min)
function scaleToCanvas(x,Y0min,Y0max,Y1min,Y1max,w,h) {
  return [(w*(x[0]-Y0min)) / (Y0max-Y0min),(h*(x[1]-Y1min)) / (Y1max-Y1min)];
}

});


// Custom unique function for arrays
Array.prototype.unique = function() {
  return this.filter(function (value, index, self) { 
    return self.indexOf(value) === index;
  });
}

