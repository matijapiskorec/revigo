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
        dataLoaded: false, // for user data (enrichments)
        basicDataLoaded: false, // for server data (gene ontology and annotations)
        // biological_processes, cellular component, molecular function
        ontologyName: "biological processes", 
        lcaWorker: [],
        canvas: [],
        context: [],
        termsToId: [],
        results: [],
        distmat: [], 
        termsData: [], // data on selected GO terms
        loadingStatus: {'start': true, 
                        'loading': true, 
                        'length': 0, 
                        't0': 0, 
                        't1': 0, 
                        'label': '', 
                        'warning':''}, 
        similarityMeasure: 'resnik', // resnik, lin, rel, jiang
        showExperimentSelector: false,
        experimentFeatures: [],
        selectedFeatures: [],
        experiments: [],
        selectedExperiment: [],
        experimentEnrichments: [],
        multipleExperimentSetting: false,
        thresholdPvalue: 0.01
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

        // For server data - gene ontology and annotations
        basicDataLoaded: function(newBasicDataLoaded,oldBasicDataLoaded) {
            if (newBasicDataLoaded) {
                console.log("revigo/watch/basicDataLoaded: Loaded basic data!");

                // A hacky way to force mounting of scatter plot.
                this.termsData = [];
            }
        },

        // For user data - enrichments
        dataLoaded: function(newDataLoaded,oldDataLoaded) {

            if (newDataLoaded) {
                console.log("revigo/watch/dataLoaded: All data successfully loaded!");
                this.filterEnrichments(this.thresholdPvalue);
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
            this.filterEnrichments(this.thresholdPvalue);

            this.calculateLCA();
        },

        similarityMeasure: function(newSimilarityMeasure, oldSimilarityMeasure) {
            console.log("revigo/watch/similarityMeasure: Changed value to "+newSimilarityMeasure);
            this.calculateDistanceMatrix();
        },

        results: function(newResults, oldResults) {
            this.calculateDistanceMatrix();
        },

        selectedFeatures: function(newSelectedFeatures, oldSelectedFeatures) {
            console.log('Somebody selected new subset of features: ' + newSelectedFeatures);

            // TODO: I'm still not sure how to implement selecting GO terms by experiment features.
            //       The problem is that some feature combinations will not correspond to any experiment,
            //       and some will correspond to multiple experiments and so it will not be clear which
            //       p-value from which experiment to use for node color.
            //       This is why for now I only implemented selection by the whole experiment.

        },

        selectedExperiment: function(newSelectedExperiment, oldSelectedExperiment) {
            var vm = this; 
            console.log('Somebody selected a new experiment: ' + newSelectedExperiment);

            // Update termsdata with the information from the newly selected experiment.
            // An assumption is that termsdata already contains all GO terms which exist in our DAG
            // and are lower than threshold p-value. GO terms visibility will reflect 'selected' attribute.
            this.termsData = new Map([...this.termsData].map( function(d) {
                var key = d[0];
                var value = d[1];
                var temp = vm.experimentEnrichmentsFiltered.get(d[0]);
                value['pvalue'] = temp.hasOwnProperty(newSelectedExperiment) ? 
                                  temp[newSelectedExperiment] :
                                  vm.thresholdPvalue, 
                value['selected'] = vm.experimentEnrichmentsFiltered.get(d[0])
                                      .hasOwnProperty(newSelectedExperiment);
                return [key,value];
            }));

        }
    },

    methods: {

        receiveDataFromChild: function(value) {

            var vm = this; 

            // TODO: Enrichments are passed to this function as a Map, although for a single experiment
            //       we are working with the Object later. This is because it is easier to check whether
            //       experiment is single or multiple if we now at least the type of variable.
            //       In the single case the Map is actually converted to Object for further processing.

            if (typeof(value.entries().next().value[1])=='number') {
                // If data has a simple format then it is a single experiment
                //{"GO:00001": 0.0001}
                this.multipleExperimentSetting = false;
            } else if (typeof(value.entries().next().value[1])=='object') {
                // If data has complex format it is a multiple experiment
                // {"GO:00001": [{"enrichment": 0.0001}]}
                this.multipleExperimentSetting = true;
            }

            if (!this.multipleExperimentSetting) {

                console.log("revigo/methods/receiveDataFromChild: Received data for a single experiment");

                // Hide experiment selector
                this.showExperimentSelector = false;

                // TODO: Hack to convert Map back to Object, ideally enrichments should be in one or the other!
                // this.enrichments = value;
                this.enrichments = Object.fromEntries([...value]);

                this.dataLoaded = true;
                this.loadingStatus['warning'] = '';

            } else {

                console.log("revigo/methods/receiveDataFromChild: Received data for multiple experiments");

                // Show experiment selector
                this.showExperimentSelector = true;

                // Experiment features are extracted from the first experiment of the first GO term
                // TODO: Consider whether column labels can be extracted and passed earlier!
                var temp = Object.keys([...value][0][1][0]);
                temp.splice(temp.indexOf('enrichment'),1); // remove 'enrichment' (modifies array in place)
                this.experimentFeatures = temp;

                // Extract all feature combinations for all experiments in enrichments dataset
                // These experiment codes will be used for experiment selectors
                this.experiments = [...value].map( x => // Map to Array
                    x[1].map( y => vm.experimentFeatures // first element is a list of experiments
                                     .filter( a => y[a] != '' ) // filter only non-empty feature names
                                     .map( a => y[a] ).join('_') )).flat().unique();

                // Select first experiment by default by setting the selectedExperiment variable
                this.selectedExperiment = this.experiments[0];

                // Only GO terms that contain enrichment for a particular experiment
                // TODO: We still have all experiments for those particular GO terms!
                //       We have to filter experiments as well!
                //       Check whether this is really needed!
                var valueFiltered = new Map([...value].filter( x =>  
                    x[1].map( y => vm.experimentFeatures
                                     .filter( a => y[a] != '' )
                                     .map( a => y[a] ).join('_') ).some( y => y == this.experiments[2] )
                ));

                // Stores all GO terms, codes of experiments where they appear and corresponding enrichments
                // TODO: Should this replace termsData? A master format for all enrichments?
                //       There are properties which are not tied to experiments - e.g. frequency.
                this.experimentEnrichments = new Map([...value].map( x =>
                    [x[0], Object.fromEntries( x[1].map( y => 
                        [ vm.experimentFeatures
                            .filter( z => y[z] != '' )
                            .map( z => y[z] ).join('_'),
                          y['enrichment']] ))]
                ));

                // TODO: Hacky solution to make multiple experiment format identical to the single experiment.
                //       I just picked first experiment's enrichment from the filtered Map.
                //       Ideally the specific experiment will be selected with the menu, with the first
                //       one selected by default. I can even have a separate function for selection.
                var enrichments = Object.fromEntries([...value].map( x=> [x[0], x[1][0]['enrichment']] ));
                console.log(enrichments);

                // TODO: Check whether enrichment variable is used anywhere anymore!

                // If there are zero p-values, replace them all with the smalles p-value which is not zero!
                // Ideally this should not happen, but it messes our color log scale (zero is infinity)
                // and so this is a pragmatic way to deal with it, along with reporting a warning message.

                var minPvalue = [...this.experimentEnrichments]
                                        .map(x=>Object.values(x[1]))
                                        .flat()
                                        .reduce((min,x) => min <= x ? min : x, Infinity);

                // If there are zero p-values replace them all with the smallest non-zero p-value
                if (minPvalue==0.0) {

                    // Minimum non-zero p-value
                    var minNonZeroPvalue = [...this.experimentEnrichments]
                                                    .map(x=>Object.values(x[1]))
                                                    .flat()
                                                    .reduce((min,x) => x<min && x!==0 ? x : min, 1.0);

                    // TODO: We now use two variables to store enrichments in different ways.
                    //       The experimentEnrichments is more general as it stores all pvalues from all exp.

                    // Replace all zero p-values with the smallest non-zero p-value
                    enrichments = Object.fromEntries(
                                    Object.entries(enrichments)
                                          .map( x => [x[0],x[1] == 0.0 ? minNonZeroPvalue : x[1]] )
                                  );

                    // Replace all zero p-values with the smallest non-zero p-value
                    this.experimentEnrichments = new Map( 
                        [...this.experimentEnrichments]
                                .map( x => [ x[0],
                                             Object.fromEntries(
                                                 Object.entries(x[1])
                                                       .map(y=>[y[0],y[1]==0.0?minNonZeroPvalue:y[1]]))
                                           ] )
                    );

                    this.loadingStatus['warning'] = 'All zero p-values replaced with the smallest '+
                                                    'non-zero p-value!';
                }

                // Pass the enrichments from mutliple experiments to be visualized
                this.enrichments = enrichments;
                this.dataLoaded = true;

            }
        },

        // Select only enrichements from a current ontology space and with small p-value
        // TODO: Later we can use the same function for Fran's redundancy reduction algorithm
        //       where GO terms which are too general are filtered out.
        filterEnrichments: function(pvalue) {
            var vm = this; 

            // Single experiment setting
            var enrichmentsSelected = Object.entries(this.enrichments)
                                            .filter( x => this.dag.hasOwnProperty(x[0]) )
                                            .filter( x => x[1] < this.thresholdPvalue);

            if (!this.multipleExperimentSetting) {

                // Create one dataset with all information on GO terms which will be passed to visuals!
                this.termsData = new Map(enrichmentsSelected.map( function(d,i) {
                    var key = d[0];
                    var value = d[1];
                    return [key,
                            {'pvalue': value, 
                             'frequency': vm.terms[d[0]] || 1,
                             'selected': true}
                    ];
                }));

            } else {

                // In the multiple experiment setup the first experiment is selected by default,
                // but we draw all GO terms which satisfy p-value threshold, regardless of whether
                // they appear in the first experiment, although the missing GO terms will be drawn
                // as transparent.

                // Store all p-values of all experiments but only those that are above p-value threshold. 
                // All of these GO terms will be embedded, but only some will be actually visible!
                this.experimentEnrichmentsFiltered = 
                    new Map( [...vm.experimentEnrichments]
                        .filter( x => vm.dag.hasOwnProperty(x[0]) ) // make sure GO term exists in our DAG
                        .map( x => [ x[0],
                                     Object.fromEntries(
                                        Object.entries(x[1]).filter(y => y[1] < vm.thresholdPvalue ) 
                                     ) 
                                   ] )
                        .filter( x => Object.keys(x[1]).length != 0 )
                    );
                
                // We use a new variable experimentEnrichmentsFiltered which already has only those
                // GO terms that exist in our DAG and only enrichments lower than threshold p-value.
                // First experiment is selected by default!
                this.termsData = 
                    new Map( [...vm.experimentEnrichmentsFiltered] 
                        .map( x => [ x[0],
                                     {'pvalue': x[1].hasOwnProperty(vm.experiments[0]) ? 
                                                x[1][vm.experiments[0]] : vm.thresholdPvalue, 
                                      'frequency': vm.terms[x[0]] || 1,
                                      'selected': x[1].hasOwnProperty(vm.experiments[0])} ] )
                    );

            }
            
        },

        // Fetch all needed data - DAG, term counts, enrichments
        fetchData: function() {

            // Needed because => functions have no defined this property
            var vm = this; 

            // Wait for all data to load
            // NOTE: We are not loading example enrichments data anymore, the visualization starts empty!
            Promise.all(["data/go-dag-molecular-function.json",
               "data/go-dag-cellular-component.json",
               "data/go-dag-biological-process.json",
               "data/go-terms-count-goa.json"].map(url=>vm.getUrl(url))) 
               .then(([dagMol,dagCell,dagBio,terms]) => { 
                    vm.dag = dagBio;
                    vm.dagBio = dagBio;
                    vm.dagMol = dagMol;
                    vm.dagCell = dagCell;
                    vm.terms = terms; 
                    vm.basicDataLoaded = true;
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

            let P = this.generate_pairs([...this.termsData.keys()]);
            this.termsToId = new Map([...this.termsData.keys()].map( (x,i) => [x,i] ));

            // Calculation of LCA in the worker
            let t0 = performance.now();
            this.lcaWorker.postMessage([P,this.dag]);

            vm.loadingStatus = {'start':false,
                                'loading':true, 
                                'length':vm.termsData.size, 
                                't0':0, 
                                't1':0, 
                                'label':vm.ontologName,
                                'warning':vm.loadingStatus['warning']};

            console.log('revigo/calculateLCA: Sent pairs and dag to worker!');

            // Wait for the worker to return LCA results
            this.lcaWorker.onmessage = function(e) {
                
                let t1 = performance.now();

                vm.loadingStatus = {'start':false,
                                    'loading':false, 
                                    'length':vm.termsData.size, 
                                    't0':t0, 
                                    't1':t1, 
                                    'label':vm.ontologName,
                                    'warning':vm.loadingStatus['warning']};

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
            // Custom max function because Math.max is recursive and fails for large arrays
            let maxValue = similarity.map(x=>x[2]).reduce((max, v) => max >= v ? max : v, -Infinity);
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
                    return results.map( (x) => [x[0], 
                                                x[1], 
                                                -Math.log((terms[x[2]]||1)/totalAnnotations)] );
                    break;

                case 'lin':

                    // Lin method: sim(t1,t2) = 2*IC(LCA) / (IC(t1)+IC(t2))
                    return results.map( (x) => [x[0], 
                                                x[1], 
                                                2*(-Math.log((terms[x[2]]||1)/totalAnnotations))/
                                                (-Math.log((terms[x[0]]||1)/totalAnnotations)
                                                 -Math.log((terms[x[1]]||1)/totalAnnotations)) ]);
                    break;

                case 'rel':

                    // Rel method: sim(t1,t2) = 2*IC(LCA)*(1-p(LCA)) / (IC(t1)+IC(t2))
                    return results.map( (x) => [x[0], 
                                                x[1], 
                                                2*(-Math.log((terms[x[2]]||1)/totalAnnotations))*
                                                (1 - ((terms[x[2]]||1)/totalAnnotations))/ 
                                                (-Math.log((terms[x[0]]||1)/totalAnnotations)
                                                 -Math.log((terms[x[1]]||1)/totalAnnotations)) ]);
                    break;

                case 'jiang':

                    // Jiang method: sim(t1,t2) = 1 - min(1,IC(t1)+IC(t2)-2*IC(LCA))
                    return results.map( (x) => [x[0], 
                                                x[1], 
                                                1 - Math.min(1, (-Math.log((terms[x[0]]||1)/totalAnnotations)
                                                -Math.log((terms[x[1]]||1)/totalAnnotations))
                                                -2*(-Math.log((terms[x[2]]||1)/totalAnnotations)) )] );
                    break;

                default:
                    return [];
            }

        }
    }

});

Vue.component('progress-box', {
    // Fields in loadingStatus are "loading", "length", "t0", "t1", "label", "warning"
    props: ["loadingstatus","thresholdpvalue"], 
    data: function() {
        return;
    },
    computed: {
        message: function() {
            if (this.loadingstatus['start']) {
                return '';
            } else if (this.loadingstatus['loading']) {
                return 'Calculating semantic similarities between '+this.loadingstatus['length']+
                       ' selected '+this.loadingstatus['label']+
                       ' GO terms with p<'+this.thresholdpvalue+' <img src="static/ajax-loader.gif">';
             } else if (!this.loadingstatus['loading']) {
                return "Calculated semantic similaritites between "+this.loadingstatus['length']+
                       " selected "+this.loadingstatus['label']+
                       " GO terms with p<"+this.thresholdpvalue+" GO terms in "+
                       (this.loadingstatus['t1']-this.loadingstatus['t0']).toFixed(0)+" miliseconds!"+
                       '<span style="color:red">'+
                       ((this.loadingstatus['warning']!='')?(' '+this.loadingstatus['warning']):'')+
                       '</span>';
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

        // Example enrichments data in CSV format which we will present to the user
        this.examplesLinks = ["data/revigo-enrichments1.csv",
                              "data/revigo-enrichments2.csv",
                              "data/revigo-enrichments3.csv"]; 

    },
    methods: {

        // Parse CSV data to JSON and send it to the main Revigo app
        // Handles both single and multiple experiment cases!
        // GO term is the key, value is an array of experiment features along with the enrichment
        sendDataToParent: function() {
            console.log("input-box/methods/sendDataToParent");

            // Allowed delimiters for the CSV
            // NOTE: Only commas allow usage of empty column values (for experiment features)
            // TODO: Consider allowing semicolon as well!
            var allowedDelimiters = / |\t|,|;/;

            // TODO: Check whether dataset is for single or multiple experiments!
            var numberOfColumns;
            for (const d of this.inputData.split('\n')) {
                // Extract first data line in CSV to count the number of columns
                if ( (d.substring(0,1)!='%') && (d.substring(0,1)!='') ) {
                    numberOfColumns = d.split(allowedDelimiters).length;
                    console.log('Number of columns in CSV file is '+numberOfColumns);
                    break;
                }
            }

            // If there is more than two columns in CSV file this means we have a multiple experiment data
            if (numberOfColumns==2) {

                // Parsing CSV to JSON - single experiment setup
                // Data is in text/csv format so we have to convert it to Object
                // GO term is key, enrichment is value
                // TODO: This will be converted to Object later for further processing!
                var data = new Map(
                    this.inputData
                        .split('\n')
                        .filter(x=>(x.substring(0,1)!='%')&&(x.substring(0,1)!=''))
                        .map(x=>[x.split(/ |\t/)[0],Number(x.split(/ |\t/)[1])])
                );

            } else if (numberOfColumns>2) {

                // Column labels should be in the last header line (which begin with '%')
                // If column names are defined, it is assumed that the last two are GO terms and enrichments
                var columnLabels = '';
                for (const d of this.inputData.split('\n')) {
                    if (d.substring(0,1)=='%') {
                        columnLabels = d;
                    } else {
                        columnLabels = columnLabels.substring(1).trim();
                        columnLabels = columnLabels.split(allowedDelimiters);
                        console.log('Column names -> '+columnLabels); 
                        break;
                    }
                }

                var data = new Map();
                for (const d of this.inputData.split('\n')
                                    .filter(x=>(x.substring(0,1)!='%')&&(x.substring(0,1)!=''))) {
                    var row = d.split(allowedDelimiters);

                    // Assumption is that the last two columns are GO term and enrichment
                    var goTerm = row[row.length-2];
                    var enrichment = Number(row[row.length-1]);

                    var experimentInfo = row.slice(0,-2);
                    var featureLabels = columnLabels.slice(0,-2);

                    // Assumption is that the number of columns corresponds to the column label header
                    // TODO: We do not cover case when header is not defined at all!
                    if (experimentInfo.length==featureLabels.length) {

                        // Prepare the experiment info with the enrichment value of GO term 
                        var goTermData = Object.fromEntries(zip(featureLabels,experimentInfo));
                        goTermData['enrichment'] = enrichment;

                        if (data.has(goTerm)) {
                            // Append a new enrichment (in a new experiment) for existing GO term
                            data.get(goTerm).push(goTermData);
                        }
                        else {
                            // GO term is first encountered - create an entry with a single-element array
                            data.set(goTerm,[goTermData]);
                        }
                    } else {
                        console.log('ERROR in parsing CSV - Number of columns does not match header labels!');
                        break;
                    }

                }

            } else {
                // TODO: If there are less than two columns this should probably be reported as error!
            }

            // TODO: Emitting data should be done here, just make sure you unify parents receiving function!
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

        },

        // Load local CSV file to input box
        // TODO: Consider using drag and drop (to input box, or maybe to canvas) to load local files!
        fileChange: function(e) {

            var vm = this;

            // TODO: We expect only one file, but in theory we can load multiple files as well!
            var file = e.target.files[0];

            const reader = new FileReader();
            reader.onload = function(event) {
                var contents = event.target.result; // event.target points to FileReader
                vm.inputData = contents;
            };
            reader.readAsText(file);

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
              '<span style="display:inline-block; width:5px;"></span>'+
              '<input type="file" @change="fileChange"/>'+ 
              '</div>'
});

Vue.component('scatter-plot', {
    props: ["distmat","termsdata"],
        data: function () {
            return {
                canvas: null,
                width: 500, 
                height: 500, 
                intervalID: null,
                intervalCheckID: null,
                custom: null,
                qtree: null,
                scaleX: null,
                scaleY: null,
                scaleR: null,
                coordToData: null,
                coordinatesReady: false
            }
    },
    mounted: function() {

        // TODO: Consider putting this in separate css file!
        // Needed so that tooltips (which have absolute position) are positioned correctly
        d3.select(this.$el)
          .style('position','relative')
          .style('width','fit-content');

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

        // TODO: Legend for radius 
        d3.select(this.$el).append("div")
                           .attr("id","legendR")
                           .style("float","right")

        // TODO: Legend for p-values
        d3.select(this.$el).append("div")
                           .attr("id","legendP")
                           .style("float","right")

        // Put a placeholder text for empty canvas!

        let context = this.canvas.node().getContext("2d");
        context.clearRect(0, 0, this.width, this.height);

        context.font = "24px Helvetica";
        context.fillStyle = "#777";
        context.fillText('Please load enrichment data! :-)',0.2*this.width,0.5*this.height);

    },
    watch: {
        // Trigger when termsdata change.
        // For now we assume that node positions did not change, only node properties (enrichment-color,
        // selected-transparency) so we don't recalculate the whole embedding.
        // This watch is not triggered if we modify just the pvalue in termsdata - we have to reassign
        // the whole variable! 
        termsdata(newTermsdata,oldTermsdata) {
            // console.log("scatter-plot/watch/termsdata: Someone changed termsdata and now I react!");

            let vm = this;

            // Set scale for radius which depends on the frequency of newly loaded GO terms
            fmax = Math.max(...Array.from(this.termsdata.values()).map(x=>x.frequency));
            fmin = Math.min(...Array.from(this.termsdata.values()).map(x=>x.frequency));
            vm.scaleR.domain([fmin, fmax]);

            // Set scale for radius which depends on the pvalues of newly loaded GO terms
            pmax = Math.max(...Array.from(this.termsdata.values()).map(x=>x.pvalue));
            pmin = Math.min(...Array.from(this.termsdata.values()).map(x=>x.pvalue));
            vm.scaleP.domain([pmin, pmax]);

            // TODO: If there are undefined values in termsdata set coordinatesReady to false!
            //       This should disable drawing of nodes with undefined coordinates, but it also
            //       freezes the last frame of the previous embedding, which is worse.
            // var undefinedCoordinates = [...vm.termsdata].some( x => 
                // typeof(x['x'])=='undefined' || typeof(x['y'])=='undefined'
            // );
            // this.coordinatesReady = !undefinedCoordinates;

            // Bind data to visual elements and redraw nodes - embedding did not change so node
            // positions do not change as well, only their properties (color, visibility).

            // Bind data to visual elements
            vm.bindDataToVisuals();

            // Draw data that was bound to visual elements
            vm.drawNodes(vm.canvas);
        },
        distmat(newDataLoaded, oldDataLoaded) {
            console.log("scatter-plot/watch/distmat: Distmat changed, will check for termsdata as well!!");
            let vm = this;

            // Only when both distmat and termsdata are loaded we proceed with the watch expression
            // This avoids having a separate computed property which checks whether both termsdata and
            // distmat changed - this expression checks that both are properly loaded.
            // TODO: Not entirely correct - if both old and new data have same number
            //       of GO terms this will pass but the GO terms will not match!
            if (this.distmat.length!=0 && this.termsdata.size!=0 &&
                this.distmat.length==this.termsdata.size) {

                console.log('scatter-plot/watch/distmat: termsdata also changed, trigger redrawing!');

                // Redraw the whole scatter plot - recalculate embedding and redraw all nodes
                this.drawCanvas();

            }
        }
    },
    methods: {

        // TODO: This function is technically not drawing canvas but just calculating all elements
        //       and their properties which have to be drawn. Actual drawing is done in drawNodes()!
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
            tsne.initDataDist(this.distmat);

            // Initial iterations before we start dynamic visualization
            for(var k = 0; k < 10; k++) {
                tsne.step(); 
            }

            var goterms = Array.from(vm.termsdata.keys());
            
            // Stop the tsne calculation after 7 seconds
            // TODO: Check for the convergence of coordinates, rather than some specified time!
            //       You can get current coordinates with tsne.getSolution()
            clearInterval(this.intervalCheckID);
            this.intervalCheckID = setInterval(function() {
                clearInterval(vm.intervalID);
            },10000);

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
                // This is only used to identify node closest to the mouse pointer with quad tree
                vm.coordToData = new Map(Y.map( function(d,i) {
                    
                    // We have to stringify coordinates to use them as key, arrays won't work!
                    var key = String([vm.scaleX(d[0]),vm.scaleY(d[1])]);

                    // value inherits all fields from termsdata, including "frequency" and "pvalue"

                    // TODO: Sometimes we get an error "Cannot set property 'name' of undefined"
                    //       This happens because termsdata updates while goterms still has old value
                    //       Check if this is still an issue!
                    var value = vm.termsdata.get(goterms[i]);
                    value.name = goterms[i];
                    value.x = Y[i][0];
                    value.y = Y[i][1];
                    return [key,value];

                }));

                // Update termsdata with the node coordinates. We will use it in binding data to visuals.
                // Note that that this will trigger termsdata watcher!
                vm.termsdata = new Map([...vm.termsdata].map( function(d) {
                    var key = d[0];
                    var value = d[1];
                    var indexOfGOterm = goterms.indexOf(d[0]);
                    value['x'] = Y[indexOfGOterm][0];
                    value['y'] = Y[indexOfGOterm][1];
                    value['name'] = d[0];
                    return [key,value];
                }));

                // Now the data coordinates are ready and we will check this variable before actual drawing.
                vm.coordinatesReady = true;

                // Bind data to visual elements (does not draw anything yet!)
                vm.bindDataToVisuals();

                // Draw data that was bound to visual elements
                vm.drawNodes(vm.canvas);

            },0);


            d3.select('.mainCanvas').on('mousemove', function() {

                var xy = d3.mouse(this);

                // Finding the closest node based on the coordinates of the mouse is the only 
                // purpose of the coordToData, as it uses node coordinates as keys.

                // TODO: Sometimes we get an error here "Cannot read property 'find' of null"
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

        // Bind data to visual elements using D3. Actual drawing is done in drawNodes() function!
        bindDataToVisuals: function() {

            let vm = this;

            // Multiple experiments are implemented by changing the fillStyle of each node based
            // on its enrichment in the currently selected experiment - this is "pvalue" field.
            // Color of the node reflects the newly changed pvalue.
            // Transparency of the node reflects which features are selected.

            // TODO: Maybe define globally within the component?
            var selectedAlpha = 0.2;

            // Bind data to visual elements (does not draw anything yet!)
            // The color of non-selected nodes is not rendered, as their p-value is set to some default value.
            var join = vm.custom.selectAll('custom.circle')
                .data(Array.from(vm.termsdata.values())) 
                .enter()
                .append('custom')
                .attr('class', 'circle')
                .attr('x', function(d, i) { return vm.scaleX(d.x); })
                .attr('y', function(d, i) { return vm.scaleY(d.y); })
                .attr('radius', function(d){ return vm.scaleR(d.frequency); })
                .attr('label',function(d){ return d.name; }) 
                .attr('strokeStyle', function(d) { 
                    return d.selected ? 
                           'rgba(0,0,0,1)' : 
                           'rgba(0,0,0,'+String(selectedAlpha)+')'})
                .attr('fillStyle', function(d) { 
                    return d.selected ? 
                           vm.scaleP(d.pvalue) : 
                           'rgba(255,255,255,'+String(selectedAlpha)+')'});
            
        },

        // Draw data that was bound to visual elements with D3
        drawNodes: function(canvas) {

            let vm = this;

            // Draw legend for p-values (node color)
            vm.drawPlegend("#legendP", vm.scaleP);

            // Draw legend for GOA annotations (circle radius)
            vm.drawRlegend("#legendR", vm.scaleR);


            // Check whether coordinates of the nodes are ready
            if (vm.coordinatesReady) {

                let context = canvas.node().getContext("2d");
                context.clearRect(0, 0, vm.width, vm.height);

                var elements = vm.custom.selectAll('custom.circle');

                elements.remove();

                // Draw all nodes
                elements.each(function(d,i) { // for each virtual/custom element...
                    var node = d3.select(this);
                    context.beginPath();
                    context.fillStyle = node.attr('fillStyle'); 
                    context.strokeStyle = node.attr('strokeStyle'); 
                    context.arc(node.attr('x'), node.attr('y'), node.attr('radius'),0,2*Math.PI) 
                    context.fill();
                    context.stroke();
                });
            }

        },

        // Generate legend with continuous colors from a prespecified scale
        // Strangely, this is not built in D3 by default!?
        // http://bl.ocks.org/syntagmatic/e8ccca52559796be775553b467593a9f
        drawPlegend: function(selector_id, colorscale) {

          var legendheight = 200,
              legendwidth = 80,
              margin = {top: 20, right: 60, bottom: 10, left: 2};
          
          d3.select(selector_id).selectAll("*").remove();

          var canvas = d3.select(selector_id)
            .style("height", legendheight + "px")
            .style("width", legendwidth + "px")
            .style("position", "relative")
            .append("canvas")
            .attr("height", legendheight - margin.top - margin.bottom)
            .attr("width", 1)
            .style("height", (legendheight - margin.top - margin.bottom) + "px")
            .style("width", (legendwidth - margin.left - margin.right) + "px")
            .style("border", "1px solid #000")
            .style("position", "absolute")
            .style("top", (margin.top) + "px")
            .style("left", (margin.left) + "px")
            .node();

          var ctx = canvas.getContext("2d");

          var legendscale = d3.scaleLog()
            .range([1, legendheight - margin.top - margin.bottom])
            .domain(colorscale.domain());

          // Generate image with continuous scale colors. If too slow see faster solution bellow!
          // http://stackoverflow.com/questions/4899799
          //       /whats-the-best-way-to-set-a-single-pixel-in-an-html5-canvas   
          d3.range(legendheight).forEach(function(i) {
            ctx.fillStyle = colorscale(legendscale.invert(i));
            ctx.fillRect(0,i,1,1);
          });

          // TODO: Possibly a faster way to generate scale image.
          //       http://bl.ocks.org/mbostock/048d21cf747371b11884f75ad896e5a5
          // var image = ctx.createImageData(1, legendheight);
          // d3.range(legendheight).forEach(function(i) {
            // var c = d3.rgb(colorscale(legendscale.invert(i)));
            // image.data[4*i] = c.r;
            // image.data[4*i + 1] = c.g;
            // image.data[4*i + 2] = c.b;
            // image.data[4*i + 3] = 255;
          // });
          // ctx.putImageData(image, 0, 0);

          var legendaxis = d3.axisRight()
            .scale(legendscale)
            .tickSize(6)
            .ticks(8);

          var svg = d3.select(selector_id)
            .append("svg")
            .attr("height", (legendheight) + "px")
            .attr("width", (legendwidth) + "px")
            .style("position", "absolute")
            .style("left", "0px")
            .style("top", "0px");

          svg
            .append("g")
            .attr("class", "axis")
            .attr("transform", "translate(" + (legendwidth - margin.left - margin.right + 3) + 
                                        "," + (margin.top) + ")")
            .call(legendaxis);

          svg.append("text")
            .attr("x", 0)
            .attr("y", 12)
            .attr("fill", "currentColor")
            .attr("text-anchor", "start")
            .text('p-value'); 

        },

        // https://www.youtube.com/watch?v=XmVPHq4NhMA
        drawRlegend: function(selector_id, sizescale) {

          var legendheight = 150,
              legendwidth = 80,
              margin = {top: 20, right: 60, bottom: 10, left: 2};

          d3.select(selector_id).selectAll("*").remove();

          d3.select(selector_id)
            .style("height", legendheight + "px")
            .style("width", legendwidth + "px")
            .style("position", "relative")

          var svg = d3.select(selector_id)
            .append("svg")
            .attr("height", (legendheight) + "px")
            .attr("width", (legendwidth) + "px")
            .style("position", "absolute")
            .style("left", "0px")
            .style("top", "0px");

          svg
            .append("g")
            .attr("class", "axis")
            .attr("transform", "translate(" + (legendwidth - margin.left - margin.right + 3) + 
                                        "," + (margin.top) + ")")
            .call( function(selection) {

                var groups = selection.selectAll('g').data([1,10,100,1000,10000,100000]);
                var groupsEnter = groups.enter().append('g');

                groupsEnter
                    .merge(groups)
                    .attr('transform',(d,i) => 'translate(0,'+i*20+')');

                groups.exit().remove();

                groupsEnter
                      .append('circle')
                      .merge(groups.select('circle'))
                      .attr('r',sizescale)
                      .attr('stroke','black')
                      .attr('fill','white')

                groupsEnter
                      .append('text')
                      .merge(groups.select('text'))
                      .text(d => d)
                      .attr('dy','0.32em')
                      .attr('x',10)

              svg.append("text")
                .attr("x", 0)
                .attr("y", 12)
                .attr("fill", "currentColor")
                .attr("text-anchor", "start")
                .text('Annotations'); 

            });
        }
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
            return this.distmat.length > 0 ? 
                      'data:text/csv,' + encodeURIComponent(this.distmat.map(e=>e.toString()).join("\n")) :
                      'javascript:void(0);';
        },
        matrixShape: function() {
            return this.distmat.length+'x'+this.distmat.length;
        },
        similarityFullName: function() {
            let dict = {'resnik':'Resnik','lin':'Lin','rel':'SimRel','jiang':'Jiang'}
            return dict[this.similarity];
        }
    },
    template: '<a v-bind:href="message" download="distmat.csv">'+
              'Download {{similarityFullName}} distance matrix for {{ontology}} ({{matrixShape}})'+
              '</a>'
});

// Custom zip function
const zip = (arr1, arr2) => arr1.map((k, i) => [k, arr2[i]]); // custom zip function

