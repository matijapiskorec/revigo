var temp;

var app = new Vue({
    
    el: '#revigo',

    data: {
        message: 'Hello Vue!',
        dag: [],
        dagBio: [],
        dagMol: [],
        dagCell: [],
        terms: [],
        enrichments: [],
        dataLoaded: false,
        // biological_processes, cellular component, molecular function
        ontologyName: "biological processes", 
        lcaWorker: [],
        canvas: [],
        context: [],
        termsToId: [],
        results: [],
        distmat: [],
        enrichmentsSelected: [],
        loadingStatus: [true, 0, 0, 0, ""], // [statusLoading, length, t0, t1, label]
        similarityMeasure: 'resnik' // resnik, lin, rel, jiang
    },

    created: function() {

        // Fetch data as soon as component is created
        this.fetchData();

        // Worker which will calculate LCA in a separate thread
        // NOTE: make sure your browser does not load old cached version!
        this.lcaWorker = new Worker('static/worker-lca.js');
    },

    watch: {

        dataLoaded: function(newDataLoaded,oldDataLoaded) {

            // Sometimes this reference gets confused, especially when used in a callback
            // So it is safe to define it explicitly
            let vm = this;
            if (newDataLoaded) {

                console.log("revigo/watch/dataLoaded: All data successfully loaded!");

                // Selecting only enrichements from a specific ontology namespace
                // Selecting only small p-values
                // NOTE: It is passed to child component so we assign it only once!
                this.enrichmentsSelected = Object.entries(this.enrichments)
                                                 .filter( x => this.dag.hasOwnProperty(x[0]) )
                                                 .filter( x => x[1] < 0.01);
                this.calculateLCA();
            }
        },

        ontologyName: function(newOntologyName, oldOntologyName) {
            console.log("revigo/watch/ontologyName: Changed value to "+newOntologyName);
            if (newOntologyName=="biological processes") {
                this.dag = this.dagBio;
            } else if (newOntologyName=="molecular function") {
                this.dag = this.dagMol;
            } else if (newOntologyName=="cellular component") {
                this.dag = this.dagCell;
            }
            this.enrichmentsSelected = Object.entries(this.enrichments)
                                             .filter( x => this.dag.hasOwnProperty(x[0]) )
                                             .filter( x => x[1] < 0.01);
            this.calculateLCA();
        },

        similarityMeasure: function(newSimilarityMeasure, oldSimilarityMeasure) {
            console.log("revigo/watch/similarityMeasure: Changed value to "+newSimilarityMeasure);
            this.calculateDistanceMatrix();
        },

        results: function(newResults, oldResults) {
            this.calculateDistanceMatrix();
        }
    },

    methods: {

        // Fetch all needed data - DAG, term counts, enrichments
        fetchData: function() {

            // Needed because => functions have no defined this property
            var vm = this; 

            // Wait for all data to load
            Promise.all(["data/go-dag-molecular-function.json",
               "data/go-dag-cellular-component.json",
               "data/go-dag-biological-process.json",
               "data/go-terms-count-goa.json",
               "data/revigo-enrichments1.json"].map(url=>vm.getUrl(url)))
               .then(([dagMol,dagCell,dagBio,terms,enrichments]) => {
                    vm.dag = dagBio;
                    vm.dagBio = dagBio;
                    vm.dagMol = dagMol;
                    vm.dagCell = dagCell;
                    vm.terms = terms;
                    vm.enrichments = enrichments;
                    vm.dataLoaded = true;
           });
        },

        // Fetch generic url
        getUrl: function(url) {
            return fetch(url)
                .then((response) => {
                return response.json();
            })
        },

        calculateLCA: function() {

            let vm = this;

            let P = this.generate_pairs(this.enrichmentsSelected.map(x => x[0]));
            this.termsToId = new Map(this.enrichmentsSelected.map( (x,i) => [x[0],i] ));

            // Calculation of LCA in the worker
            let t0 = performance.now();
            this.lcaWorker.postMessage([P,this.dag]);

            // [statusLoading (true/false), length, t0, t1, label]
            vm.loadingStatus = [true, vm.enrichmentsSelected.length, 0, 0, vm.ontologyName];

            console.log('revigo/calculateLCA: Sent pairs and dag to worker!');

            // Wait for the worker to return LCA results
            this.lcaWorker.onmessage = function(e) {
                
                let t1 = performance.now();
                vm.loadingStatus = [false, // loading status (true/false)
                                    vm.enrichmentsSelected.length, // number of GO terms
                                    t0, 
                                    t1, 
                                    vm.ontologyName]; // ontology namespace of GO terms
                vm.results = e.data; 
                console.log("revigo/calculateLCA/onmessage: Calculated LCA in the worker!");

            };
        },

        // Generate all pairs of elements of given array, without repetition 
        generate_pairs: function(N) {
            let result = [];
            let i = -1;
            for (const a of N) {
                i = i + 1;
                for (const b of N.slice(i+1,)) {
                    result.push([a,b]);
                }
            }
            return result;
        },

        csvURI2: function() {
            return this.distmat.map(e=>e.toString()).join("\n");
        },

        calculateDistanceMatrix: function() {

            let similarity = this.semanticSimilarity(this.results,this.terms,this.similarityMeasure);

            // List of GO terms for visualization
            let goTerms = similarity.map(x=>[x[0],x[1]]).flat().unique()

            // Default max value which we use for t-sne distance matrix 
            let maxValue = Math.max.apply(Math, similarity.map(x=>x[2]));
            let distMat = [...Array(goTerms.length)].map(e=>Array(goTerms.length).fill(maxValue));
            for (const x of similarity) {
                let x0 = this.termsToId.get(x[0]);
                let x1 = this.termsToId.get(x[1]);
                distMat[x0][x1] = x[2];
                distMat[x1][x0] = x[2];
            }
            // Final distance matrix passed to child components is assigned only once!
            this.distmat = distMat;
            console.log("revigo/methods/calculateDistanceMatrix: Created distance matrix!");

        },

        // Calculate semantic similarity between terms
        // GO terms that do not have any annotations in GOA are considered to have at least one!
        semanticSimilarity: function(results,terms,type) {

            // results: an array of GO term triplets - [term1,term2,lca]
            // terms: number of annotations for each GO term
            // Type: resnik, lin, rel, jiang

            let totalAnnotations = Object.entries(terms).reduce((sum,x)=>sum+x[1],0);

            // If terms contains raw counts we have to calculate p(t) and IC(t) from scratch
            // IC(t) = -log(p(t))
            // p(t) = n(t)/N

            // TODO: Some GO terms in this.terms are missing, for example GO:0009208!

            switch (type) {

                case 'resnik':

                    // Resnik method: sim(t1,t2) = IC(LCA)
                    return results.map( (x) => [x[0], x[1], -Math.log((terms[x[2]]||1)/totalAnnotations)] );
                    break;

                case 'lin':

                    // Lin method: sim(t1,t2) = 2*IC(LCA) / (IC(t1)+IC(t2))
                    return results.map( (x) => [x[0], x[1], 2*(-Math.log((terms[x[2]]||1)/totalAnnotations))/(-Math.log((terms[x[0]]||1)/totalAnnotations)-Math.log((terms[x[1]]||1)/totalAnnotations)) ]);
                    break;

                case 'rel':

                    // Rel method: sim(t1,t2) = 2*IC(LCA)*(1-p(LCA)) / (IC(t1)+IC(t2))
                    return results.map( (x) => [x[0], x[1], 2*(-Math.log((terms[x[2]]||1)/totalAnnotations)) * (1 - ((terms[x[2]]||1)/totalAnnotations)) / (-Math.log((terms[x[0]]||1)/totalAnnotations)-Math.log((terms[x[1]]||1)/totalAnnotations)) ]);
                    break;

                case 'jiang':

                    // Jiang method: sim(t1,t2) = 1 - min(1,IC(t1)+IC(t2)-2*IC(LCA))
                    return results.map( (x) => [x[0], x[1], 1 - Math.min( 1, (-Math.log((terms[x[0]]||1)/totalAnnotations)-Math.log((terms[x[1]]||1)/totalAnnotations)) - 2*(-Math.log((terms[x[2]]||1)/totalAnnotations)) )] );
                    break;

                default:
                    return [];
            }

        }
    }

});

Vue.component('progress-box', {
    props: ["data"], // [statusLoading, length, t0, t1, label]
    data: function() {
        return;
    },
    computed: {
        message: function() {
            if (this.data[0]) {
                return 'Calculating semantic similarities between '+this.data[1]+' selected '+this.data[4]+' GO terms with p<0.01 <img src="static/ajax-loader.gif">';
             } else if (!this.data[0]) {
                return "Calculated semantic similaritites between "+this.data[1]+" selected "+this.data[4]+" GO terms with p<0.01 GO terms in "+(this.data[3]-this.data[2])+" miliseconds!";
             }
        }
    },
    template: '<p><span v-html="message"></span></p>'
});

Vue.component('scatter-plot', {
    props: ["distmat","enrichments"],
        data: function () {
            return {
                canvas: null,
                intervalID: null,
                distmatLocal: null, // local version of parent prop distance matrix
                enrichmentsLocal: null // local version of parent prop enrichments
            }
    },
    computed: {
        // Will run whenever either of the props changes in the parent element
        // Have to check that both props are updated to the same set of GO terms
        dataLoaded: function() {
            console.log("scatter-plot/computed/dataLoaded: Someone changed dataLoaded!");
            // TODO: Not entirely correct - if both old and new data have same number
            // of GO terms this will pass but the GO terms will not match!
            if (this.distmat.length!=0 && this.enrichments.length!=0 &&
                this.distmat.length==this.enrichments.length) {
                // Drawing will be done in watch expression 
                // TODO: Whenever distmat or enrichments change we want to send an unique value to watch! 
                return Date.now();
            }
        }
    },
    watch: {
        dataLoaded(newDataLoaded, oldDataLoaded) {
            console.log("scatter-plot/watch/dataLoaded: All data successfully loaded!");
            let vm = this;
            // Only when both distmat and enrichments are loaded we update the local version simultaneously
            this.distmatLocal = this.distmat;
            this.enrichmentsLocal = this.enrichments;
            console.log("scatter-plot/watch/dataLoaded: distmatLocal.lenght="+this.distmatLocal.length+", enrichmentsLocal.length="+this.enrichmentsLocal.length);
            this.drawPlot();
        }
    },
    methods: {
        drawPlot: function() {
            let vm = this;

            // Stop any previous setInterval() function
            clearInterval(this.intervalID);

            // Prepare canvas
            this.canvas = vm.$el; 
            this.context = this.canvas.getContext("2d");
            width = this.canvas.width,
            height = this.canvas.height;
            vm.context.clearRect(0, 0, width, height);

            // Calculating t-SNE
            var opt = {}
            opt.epsilon = 10; // learning rate 
            opt.perplexity = 10; // roughly how many neighbors each point influences 
            opt.dim = 2; // dimensionality of the embedding 

            var tsne = new tsnejs.tSNE(opt); // create a tSNE instance
            tsne.initDataDist(this.distmatLocal);

            // Initial iterations before we start dynamic visualization
            for(var k = 0; k < 10; k++) {
                tsne.step(); 
            }

            // Returns control to the browser so that canvas can be redrawn
            // Time interval is set to 0, with no delay between redrawing
            this.intervalID = setInterval(function() {
                for(var k = 0; k < 1; k++) {
                    tsne.step(); // every time you call this, solution gets better
                }
                Y = tsne.getSolution();
                Y0min = Math.min(...Y.map(x=>x[0]));
                Y0max = Math.max(...Y.map(x=>x[0]));
                Y1min = Math.min(...Y.map(x=>x[1]));
                Y1max = Math.max(...Y.map(x=>x[1]));
                vm.context.clearRect(0, 0, width, height);

                // TODO: In some cases there is an error that vm.enrichmentsLocal[i] is not defined!

                Y.forEach(function(y,i) {
                    vm.drawNode(y,Y0min,Y0max,Y1min,Y1max,
                            width,height,vm.context,
                            vm.enrichmentsLocal[i][0].substring(3));
                });
            },0);
        },

        drawNode: function(d,Y0min,Y0max,Y1min,Y1max,w,h,context,label) {
             let ds = this.scaleToCanvas(d,Y0min,Y0max,Y1min,Y1max,w,h);
             context.beginPath();
             context.fillStyle = "#9ecae1";
             context.strokeStyle = "#000";
             context.arc(ds[0],ds[1],3,0,2*Math.PI);
             context.fill();
             context.stroke();
             context.font = "9px Helvetica";
             context.fillStyle = "#000";
             context.fillText(label,ds[0]+5,ds[1]+3);
        },

        // [min,max] -> [a,b]
        // [min,max] -> [0,canvas.width]
        // f(x) = (b-a)*(x-min) / (max-min) + a
        // f(x) = canvas.width*x / (max-min)
        scaleToCanvas: function(x,Y0min,Y0max,Y1min,Y1max,w,h) {
            return [(w*(x[0]-Y0min)) / (Y0max-Y0min),(h*(x[1]-Y1min)) / (Y1max-Y1min)];
        }
    },
    template: '<canvas></canvas>'
});

// Custom unique function for arrays
Array.prototype.unique = function() {
    return this.filter(function (value, index, self) { 
        return self.indexOf(value) === index;
    });
}

Vue.component('download-csv', {
    props: ["distmat","ontology","similarity"], 
    data: function() {
        return;
    },
    computed: {
        message: function() {
            return this.distmat.length > 0 ? "data:text/csv," + encodeURIComponent(this.distmat.map(e=>e.toString()).join("\n")) : 'javascript:void(0);';
        },
        matrixShape: function() {
            return this.distmat.length+'x'+this.distmat.length;
        },
        similarityFullName: function() {
            let dict = {'resnik':'Resnik','lin':'Lin','rel':'SimRel','jiang':'Jiang'}
            return dict[this.similarity];
        }
    },
    template: '<a v-bind:href="message" download="distmat.csv">Download {{similarityFullName}} distance matrix for {{ontology}} ({{matrixShape}})</a>'
});

// A custom function for applying map to all object properties
const objectMap = (obj, fn) =>
    Object.fromEntries(
        Object.entries(obj).map(
            ([k, v], i) => [k, fn(v, k, i)]
        )
    )

