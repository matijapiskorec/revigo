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
        termsData: [], // data on selected GO terms
        loadingStatus: [true, 0, 0, 0, ""], // [statusLoading, length, t0, t1, label]
        similarityMeasure: 'resnik' // resnik, lin, rel, jiang
    },

    created: function() {

        // Fetch data as soon as component is created
        this.fetchData();

        // Worker which will calculate LCA in a separate thread
        // NOTE: make sure your browser does not load old cached version!
        this.lcaWorker = new Worker('static/worker-lca.js');

        console.log('revigo/created: New worker created!');
    },

    watch: {

        dataLoaded: function(newDataLoaded,oldDataLoaded) {

            if (newDataLoaded) {
                console.log("revigo/watch/dataLoaded: All data successfully loaded!");
                this.filterEnrichments(0.01);
                this.calculateLCA();

                // TODO: Reset the dataLoaded so that we can set it to true when needed
                // Reset will again trigger this dataLoaded watch but will have no effect as else clause is empty
                this.dataLoaded = false;
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
            this.filterEnrichments(0.01);

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

        receiveDataFromChild: function(value) {
            console.log("revigo/methods/receiveDataFromChild: Received data from input box!");

            // TODO: This is not enough to start new LCA calculation!
            this.enrichments = value;

            // TODO: This is enough to start the new LCA calulation because dataLoaded
            // is properly reset to false after loading new data
            this.dataLoaded = true;
        },

        // Select only enrichements from a current ontology space and with small p-value
        filterEnrichments: function(pvalue) {
            var vm = this; 
            // NOTE: It is passed to child component so we assign it only once!
            this.enrichmentsSelected = Object.entries(this.enrichments)
                                             .filter( x => this.dag.hasOwnProperty(x[0]) )
                                             .filter( x => x[1] < 0.01);
            
            // TODO: Create one dataset with all information on GO terms which will be passed to visuals!
            this.termsData = new Map(Object.entries(this.enrichmentsSelected).map( function(d,i) {
                var key = d[1][0];
                var value = d[1][1];
                return [key,
                        {'pvalue': value, 'frequency': vm.terms[d[1][0]] || 1}
                ];
            }));
            // console.log(this.termsData);

        },

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
                return 'Calculating semantic similarities between '+this.data[1]+
                       ' selected '+this.data[4]+
                       ' GO terms with p<0.01 <img src="static/ajax-loader.gif">';
             } else if (!this.data[0]) {
                return "Calculated semantic similaritites between "+this.data[1]+
                       " selected "+this.data[4]+
                       " GO terms with p<0.01 GO terms in "+(this.data[3]-this.data[2])+" miliseconds!";
             }
        }
    },
    template: '<p><span v-html="message"></span></p>'
});

Vue.component('input-box', {
    props: [], 
    data: function() {
        return {
            inputData: null,
            examplesLinks: null
        }
    },
    created: function() {

        // Example enrichments dat is in csv format 
        this.examplesLinks = ["data/revigo-enrichments1.csv",
                              "data/revigo-enrichments2.csv",
                              "data/revigo-enrichments3.csv"];

    },
    methods: {
        sendDataToParent: function() {
            console.log("input-box/methods/sendDataToParent");

            // Data is in text/csv format so we have to convert it to Object
            // TODO: Consider using Map instead of Object for enrichments data
            var data = Object.fromEntries(
                this.inputData
                    .split('\n')
                    .filter(x=>(x.substring(0,1)!='%')&&(x.substring(0,1)!=''))
                    .map(x=>[x.split(/ |\t/)[0],Number(x.split(/ |\t/)[1])])
            );
            this.$emit('clicked',data);
        },

        // Fetch all needed data - DAG, term counts, enrichments
        fetchData: function(link) {

            // Needed because => functions have no defined this property
            var vm = this; 

            // TODO: Make so that example data is originally in text not json format!
            fetch(link)
                .then(response => response.text()) // Use .json() if enrichment data is in json format!
                .then(data => {
                    vm.inputData = data;
                });

        }

    },
    template: '<div>'+
              'Examples: '+
              '<span v-for="(link,index) of examplesLinks">'+
              '<a class="exampleLink" @click="fetchData(link)">#{{index+1}}</a>&nbsp;'+
              '</span>'+
              '</br>'+
              '<textarea v-model="inputData" cols="70" rows="10"></textarea>'+
              '</br>'+
              '<button @click="sendDataToParent">Submit</button>'+
              '</div>'
});

Vue.component('scatter-plot', {
    props: ["distmat","enrichments","termsdata"],
        data: function () {
            return {
                canvas: null,
                width: 500, 
                height: 500, 
                intervalID: null,
                distmatLocal: null, // local version of parent prop distance matrix
                enrichmentsLocal: null, // local version of parent prop enrichments
                termsDataLocal: null,
                custom: null,
                qtree: null,
                scaleX: null,
                scaleY: null,
                scaleR: null,
                coordToData: null
            }
    },
    mounted: function() {

        // TODO: Consider putting this in separate css file!
        // Needed so that tooltips (which have absolute position) are positioned correctly
        d3.select(this.$el).style('position','relative');

        // Prepare canvas
        this.canvas = d3.select(this.$el)
                        .append('canvas')
                        .classed('mainCanvas', true)
                        .attr('width',this.width)
                        .attr('height',this.height);

		var customBase = document.createElement('custom');
		this.custom = d3.select(customBase); 

        this.scaleX = d3.scaleLinear().range([0, this.width]);
        this.scaleY = d3.scaleLinear().range([0, this.height]);
        this.scaleR = d3.scaleLog().range([4, 10]);
        this.scaleP = d3.scaleLog().range(['red','steelblue']);
        
        d3.select(this.$el).append("div").attr("id","tooltip");

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
            this.termsDataLocal = this.termsdata;

            console.log("scatter-plot/watch/dataLoaded: distmatLocal.lenght="+this.distmatLocal.length+
                        ", enrichmentsLocal.length="+this.enrichmentsLocal.length+
                        ", termsDataLocal.length="+this.termsDataLocal.size);

            // Set scale for radius which depends on the pvalues of newly loaded GO terms
            fmax = Math.max(...Array.from(this.termsDataLocal.values()).map(x=>x.frequency));
            fmin = Math.min(...Array.from(this.termsDataLocal.values()).map(x=>x.frequency));
            vm.scaleR.domain([fmin, fmax]);

            // Set scale for radius which depends on the pvalues of newly loaded GO terms
            pmax = Math.max(...Array.from(this.termsDataLocal.values()).map(x=>x.pvalue));
            pmin = Math.min(...Array.from(this.termsDataLocal.values()).map(x=>x.pvalue));
            vm.scaleP.domain([pmin, pmax]);

            this.drawCanvas();
        }
    },
    methods: {

        drawCanvas: function() {
            let vm = this;

            // Stop any previous setInterval() function
            console.log("scatter-plot/methods/drawCanvas: Clearing setInterval with ID: "+this.intervalID);
            clearInterval(this.intervalID);

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

            var goterms = Array.from(vm.termsDataLocal.keys());
            
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

                vm.scaleX.domain([Y0min, Y0max]);
                vm.scaleY.domain([Y1min, Y1max]);

                vm.qtree = d3.quadtree()
                  .addAll(Y.map(d=>[vm.scaleX(d[0]),vm.scaleY(d[1])]));

                // Connect point coordinates with the information on the corresponding GO term
                vm.coordToData = new Map(Y.map( function(d,i) {
                    
                    // We have to stringify coordinates to use them as key, arrays won't work!
                    var key = String([vm.scaleX(d[0]),vm.scaleY(d[1])]);

                    // TODO: Possible error "value is not defined"
                    // This happens because termsLocalData updates while goterms still has old value
                    var value = vm.termsDataLocal.get(goterms[i]);
                    value.name = goterms[i];
                    value.x = Y[i][0];
                    value.y = Y[i][1];
                    return [key,value];

                }));


                // Bind data to visual elements (does not draw anything yet!)
                // TODO: Maybe move variable declaration outside setInterval?
                var join = vm.custom.selectAll('custom.circle')
                    .data(Array.from(vm.coordToData.values()))
                    .enter()
                    .append('custom')
                    .attr('class', 'circle')
                    .attr('x', function(d, i) { return vm.scaleX(d.x); })
                    .attr('y', function(d, i) { return vm.scaleY(d.y); })
                    .attr('radius', function(d){ return vm.scaleR(d.frequency); })
                    .attr('label',function(d){ return d.name; }) 
                    .attr('fillStyle', function(d) { return vm.scaleP(d.pvalue); }); 

                // Draw data that was bound to visual elements
                // TODO: Not sure this is really needed?!
                var exitSel = join.exit()
                    .transition()
                    .attr('radius', 0)
                    .remove();

                if (Y.length==vm.enrichmentsLocal.length) {
                    vm.drawNodes(vm.canvas,false);
                }

            },0);


            d3.select('.mainCanvas').on('mousemove', function() {

                var xy = d3.mouse(this);

                var xyTooltip = vm.qtree.find(xy[0],xy[1]);
                var closestNode = vm.coordToData.get(String(xyTooltip));

                // If mouse cursor is close enough to the closest point show tooltip
                if (Math.abs(xy[0]-xyTooltip[0])<=vm.scaleR(closestNode.frequency) && 
                    Math.abs(xy[1]-xyTooltip[1])<=vm.scaleR(closestNode.frequency)) {

                    // Show the tooltip only when there is nodeData found by the mouse
                    d3.select('#tooltip')
                      .style('opacity', 0.8)
                      .style('top', xy[1]+1+'px') 
                      .style('left', xy[0]+1+'px') 
                      .html(function() { 
                          return closestNode.name+
                                 "</br>GOA annotations: "+closestNode.frequency+
                                 "</br>p-value: "+closestNode.pvalue; 
                      });
                } else {
                    // Hide the tooltip when there our mouse doesn't find nodeData
                    d3.select('#tooltip')
                      .style('opacity', 0);
                }
            }); 

        },

        drawNodes: function(canvas) {
            let vm = this;
            let context = canvas.node().getContext("2d");
            context.clearRect(0, 0, vm.width, vm.height);

            var elements = vm.custom.selectAll('custom.circle');

            elements.remove();

            // Draw all nodes
            elements.each(function(d,i) { // for each virtual/custom element...
                var node = d3.select(this);
                context.beginPath();
                context.fillStyle = node.attr('fillStyle'); 
                context.strokeStyle = "#000";
                context.arc(node.attr('x'), node.attr('y'), node.attr('radius'),0,2*Math.PI) 
                context.fill();
                context.stroke();
            });

            // Draw all node labels
            // TODO: Labels are disabled because all relevant information is in the tooltip!
            // elements.each(function(d,i) { // for each virtual/custom element...
                // var node = d3.select(this);
                // context.font = "9px Helvetica";
                // context.fillStyle = "#000";
                // context.fillText(node.attr('label'),
                                    // parseFloat(node.attr('x'))+parseFloat(node.attr('radius')-1),
                                    // parseFloat(node.attr('y'))-parseFloat(node.attr('radius'))+1);
            // });
        },

    },
    template: '<div></div>'
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

