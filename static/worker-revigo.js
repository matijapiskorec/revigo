
// Web worker script that calculates non-redundant GO terms with Revigo algorithm

// Script for calculation of Lowest Common Ancestor
// NOTE: When changing this script make sure the browser is not loading old cached version!
importScripts('lca.js'); 

onmessage = function(e) {

    console.log("Revigo web worker: I received data!");
    console.log(e);

    let algorithmType = e.data[0]; // LCA or NMF

    let enrichmentsFiltered = e.data[1];
    let terms = e.data[2];
    let dag = e.data[3];

    // Experiment is passed as either fourth (LCA case) or fifth (NMF case) argument
    let experiment;
    if (algorithmType=='NMF') {
        experiment = e.data[5];
    } else {
        experiment = e.data[4];
    }

    let go_terms = enrichmentsFiltered.map(x=>x[0]);
    let totalAnnotations = Object.entries(terms).reduce((sum,x)=>sum+x[1],0);

    let semanticSimilarityMatrix = [];

    // TODO: Map where remaining (non-redundant) terms are keys and their redundant counterparts are values.
    // let revigoTerms = {};
    let revigoTerms = new Map(go_terms.map(x=>[x,[]]))

    if (algorithmType=='NMF') {

        let nmfEmbedding = e.data[4];

        // Filter the NMF embeddings to contain only terms in user's dataset
        let nmfEmbeddingFiltered = Object.entries(nmfEmbedding).filter(x=>go_terms.includes(x[0]));

        let W = nmfEmbeddingFiltered.map(x=>x[1]["W"]);
        let H = nmfEmbeddingFiltered.map(x=>x[1]["H"]);

        semanticSimilarityMatrix = multiply(W,H);

    } else {

        let P = generate_pairs(go_terms);
        let total = P.length;
        let results = P.map( (p,i) => [p[0], p[1], LCA(dag,p[0],p[1])] );
        let termsToId = new Map(go_terms.map( (x,i) => [x,i] ));

        // SimRel method: sim(t1,t2) = 2*IC(LCA)*(1-p(LCA)) / (IC(t1)+IC(t2))
        let similarity =  results.map( (x) => [x[0], 
                                    x[1], 
                                    2*(-Math.log((terms[x[2]]||1)/totalAnnotations))*
                                    (1 - ((terms[x[2]]||1)/totalAnnotations))/ 
                                    (-Math.log((terms[x[0]]||1)/totalAnnotations)
                                     -Math.log((terms[x[1]]||1)/totalAnnotations)) ]);

        // Default max value which we use for t-sne distance matrix 
        // Custom max function because Math.max is recursive and fails for large arrays
        let maxValue = similarity.map(x=>x[2]).reduce((max, v) => max >= v ? max : v, -Infinity);
        semanticSimilarityMatrix = [...Array(go_terms.length)].map(e=>Array(go_terms.length)
                                                                 .fill(maxValue));
        for (const x of similarity) {
            let x0 = termsToId.get(x[0]);
            let x1 = termsToId.get(x[1]);

            semanticSimilarityMatrix[x0][x1] = x[2];
            semanticSimilarityMatrix[x1][x0] = x[2];
        }

        // Set diagonal to zero
        semanticSimilarityMatrix = semanticSimilarityMatrix
            .map( (x,i) => x.map( (y,j) => j == i ? 0 : y ) );

    }

    // Find the max element in 2D array and its indices in both arrays
    let max_element = findMax2D(semanticSimilarityMatrix);

    let remainingMatrix = semanticSimilarityMatrix;
    let removed_terms = []; // Terms removed as redundant

    // A set of non-redundant (remaining) terms which will be returned as a result
    // Also used as a reference to track terms in similarity matrix!
    let remaining_terms = go_terms; 

    // Index of the removed term in the remaining_terms array, not the original go_terms!
    let removed_index; 

    // While the most similar pair is above a predefined similarity threshold
    // TODO: This threshold should be set as a tunable parameter!
    // TODO: This is not neccessarily the same for LCA and NMF-based approaches!
    while (max_element[0] > 1.0) {

        // Remove one of the terms depending on several criteria:
        // 1. If one term has very broad interpretation (freq. > 5%) remove it!
        // 2. If one term has lower (up to a threshold) enrichment (p-value) remove it!
        // 3. If one term is a child of the other remove it!
        // 4. If none of the criteria is matched then remove one of the terms at random.

        // The two terms, one of which should be removed
        let term1_id = max_element[1][0];
        let term2_id = max_element[1][1];

        let term1 = remaining_terms[term1_id];
        let term2 = remaining_terms[term2_id];

        // Calculate term frequencies - used in step 1.
        let freq1 = terms[term1] / totalAnnotations;
        let freq2 = terms[term2] / totalAnnotations;

        // Fetch term enrichments - used in step 2.
        // TODO: enrichmentsFiltered is currently an Array, consider making it a 
        //       Map so that you can access it more easily here!
        let pvalue1 = enrichmentsFiltered.filter( x => x[0] == term1 );
        let pvalue2 = enrichmentsFiltered.filter( x => x[0] == term2 );

        // If (only) one of the terms has frequency of more than 5%
        // TODO: This frequency should be a tunable parameter!
        if ((freq1 > 0.05)^(freq2 > 0.05)) {

            // If the first term has frequency > 5% remove it, otherwise remove the second!
            removed_index = (freq1>0.05) ? term1_id : term2_id;

        // 2. One term has lower (up to a threshold) enrichment
        // TODO: This threshold should be a tunable parameter!
        } else if ( Math.abs(pvalue1-pvalue2) > 0.005 ) {
            
            // If the first term has lower enrichment remove it!
            removed_index = (pvalue1<pvalue2) ? term1_id : term2_id; 
            
        // 3. If either one term is a parent of another then we remove the child term! 
        } else if (dag[term1].includes(term2) || dag[term2].includes(term1)) {

            // If term2 is parent of term1 then remove term1, otherwise remove term2
            removed_index = dag[term1].includes(term2) ? term1_id : term2_id;

        // 4. If none of the criteria is matched then remove one of the terms at random.
        // TODO: Original Revigo algorithm uses deterministic choice for reproducibility! 
        //       1. Math.random() does not allow setting a seed unfortunatelly:-(
        //       2. You can choose a term by using a deterministic function of their ID's!
        } else {

            removed_index = Math.random() > 0.5 ? term1_id : term2_id;

        }

        removed_terms.push(remaining_terms[removed_index]); 

        // Record that remaining term is a representative of a (redundant) removed term
        // NOTE: Both remaining_index and removed_index are indices in the remaining_term array!
        // TODO: I guess we can use revigoTerms as the main output variable?
        remaining_index = (removed_index == term1_id) ? term2_id : term1_id;
        term_remaining = remaining_terms[remaining_index]
        term_removed = remaining_terms[removed_index]
        revigoTerms.set(term_remaining,
                        revigoTerms.get(term_remaining).concat(term_removed));
        revigoTerms.delete(term_removed);
        // if (remaining_terms[remaining_index] in revigoTerms) {
        //     revigoTerms[remaining_terms[remaining_index]].push(remaining_terms[removed_index]);
        // } else {
        //     revigoTerms[remaining_terms[remaining_index]] = [remaining_terms[removed_index]]; 
        // }

        // Send the progress as an intermediary output from web worker! 
        if (typeof experiment == 'undefined') {
            postMessage(["progress",removed_terms.length,go_terms.length]);
        } else {
            // Send experiment name as well in case of multiple experiments!
            postMessage(["progress",experiment,removed_terms.length,go_terms.length]);
        }

        // Remove the term from the reference go term array
        remaining_terms = remaining_terms.filter( (x,i) => i!=removed_index );

        // Remove corresponding row and column from similarity matrix
        remainingMatrix = remainingMatrix.map( arr => arr.filter( (x,i) => i!=removed_index ) )
                                         .filter( (x,i) => i!=removed_index );

        // Find the next most similar pair
        max_element = findMax2D(remainingMatrix);
    }

    console.log(revigoTerms);
    console.log(remaining_terms);

    // TODO: Union of all remaining and removed terms is not neccessarily equal to all GO terms!
    //       All other GO terms are considered non-redundant by default!
    // let temp = new Set(Object.values(revigoTerms)
    //                          .reduce((arr1,arr2)=>arr1.concat(arr2),[])
    //                          .concat(Object.keys(revigoTerms)));
    // console.log("All remaining and removed terms...");
    // console.log(temp);

    // If experiment argument exists, pass it along with the results
    // We use this to pass the name of the experiment in the multiple experiment case
    if (typeof experiment == 'undefined') {
      // postMessage(remaining_terms);
      postMessage(revigoTerms);
    } else {
      // postMessage([remaining_terms,experiment]);
      postMessage([revigoTerms,experiment]);
    }

}

// Custom function for finding maximum value and its index in 2D array
const findMax2D = (arr) => arr.map( arr => arr.reduce( (max,x,i) => ( x > max[0] ? [x,i] : max ), [0,0] ) )
                              .reduce( (max,x,i) => ( x[0] > max[0] ? [x[0],[i,x[1]]] : max ), [0,[0,0]] );

// Generate all pairs of elements of given array, without repetition 
const generate_pairs = 
    function(N) {
        let result = [];
        let i = -1;
        for (const a of N) {
            i = i + 1;
            for (const b of N.slice(i+1,)) {
                result.push([a,b]);
            }
        }
        return result;
    };

// Multiply two matrices stored as 2D arrays
// https://stackoverflow.com/questions/27205018/multiply-2-matrices-in-javascript
const multiply =
    function(A, B) {
        var result = new Array(A.length).fill(0).map(row => new Array(B[0].length).fill(0));

        return result.map((row, i) => {
            return row.map((val, j) => {
                return A[i].reduce((sum, elm, k) => sum + (elm*B[k][j]) ,0)
            })
        })
    };

