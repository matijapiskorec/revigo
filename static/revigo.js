var app = new Vue({
    
    el: '#revigo',

    data: {

        // Input data
        dag: [],
        dagBio: [],
        dagMol: [],
        dagCell: [],
        terms: [],
        names: [],
        enrichments: [],

        // Status of the application
        dataLoaded: false, // for user data (enrichments)
        basicDataLoaded: false, // for server data (gene ontology and annotations)
        lcaStatus: {'algorithm': 'LCA',
                    'start': true, 
                    'loading': true, 
                    'termCount': 0, 
                    'duration': 0,
                    'namespace': '', 
                    'result': '', // If we need to report on the results of the calculation in some way
                    'warning': ''}, 
        revigoStatus: {'algorithm': 'Revigo',
                       'start': true, 
                       'loading': true, 
                       'progress': 0, // In percentages (LCA calculation doesn't have this!)
                       'experiment': '', // Only in multiple experiment case!
                       'termCount': 0, 
                       'duration': 0,
                       'namespace': '', 
                       'result': '', // If we need to report on the results of the calculation in some way
                       'warning': ''}, 
        umapdisabled: false, // for disabling umap radio button if needed
        filterRevigoDone: false, // for Revigo algorithm
        resetEmbeddingType: true,
        revigoStartTime: '',

        // Intermediary data
        lcaWorker: [],
        revigoWorker: [],
        termsToId: [],
        results: [],
        experimentFeatures: [],
        selectedFeatures: [],
        experiments: [],
        selectedExperiment: [],
        enrichmentsFiltered: [], // TODO: Currently an Array, consider making it a Map!
        enrichmentsFilteredDAG: [],
        experimentEnrichments: [],
        experimentEnrichmentsFiltered: [],
        experimentEnrichmentsFilteredDAG: [], // TODO: Still testing this!
        experimentDataTable: [],
        // experimentDataTableUnfiltered: [], 

        // Calculated data
        distmat: [], 
        termsData: [], // data on selected GO terms
        experimentData: [], // data on experiments
        nodeData: [], // data on nodes which will be visualized
        simGICmatrix: [],
        experimentHierarchy: [],
        semanticSimilarityMatrix: [],
        remaining_terms_all: new Map(),
        final_remaining_terms: [],
        redundancyMap: new Map(), // Which term is represented by some other term in the Revigo hierarchy!

        // Visualization type
        multipleExperimentSetting: false, // Defaults to GO-GO visualization context
        visualizationContext: 'initial', // go-go, go-ex, ex-ex, initial
        embeddingType: '', // tsne, mds, umap, biplot, hierarchical

        // Visualization parameters
        similarityMeasure: 'resnik', // resnik, lin, rel, jiang
        selectedPvalue: 0.01, // for SimGIC
        pvalueThresholds: [0.01,0.05,0.1,0.5], // p-value thresholds for SimGIC
        ontologyName: "biological processes", // biological processes, cellular component, molecular function

        // TODO: Testing and development parameters!
        revigoAlgorithmType: 'LCA'
        // revigoAlgorithmType: 'NMF'

    },

    created: function() {

        // Fetch data as soon as component is created
        this.fetchData();

        // Worker which will calculate LCA in a separate thread
        // NOTE: make sure your browser does not load old cached version!
        this.lcaWorker = new Worker('static/worker-lca.js');

        // Use different web worker depending whether we are using NMF embedding for Revigo or not!
        this.revigoWorker = new Worker('static/worker-revigo.js');

        console.log('revigo/created: New worker created!');
    },

    watch: {

        // For server data - gene ontology and annotations
        basicDataLoaded: function(newBasicDataLoaded,oldBasicDataLoaded) {
            if (newBasicDataLoaded) {
                console.log("revigo/watch/basicDataLoaded: Loaded basic data!");

                // A hacky way to force mounting of scatter plot.
                this.nodeData = [];
            }
        },

        // For user data - enrichments
        dataLoaded: function(newDataLoaded,oldDataLoaded) {

            if (newDataLoaded) {
                var vm = this; 
                console.log("revigo/watch/dataLoaded: All data successfully loaded!");

                this.filterEnrichmentsDAG();

                // Call a web worker for calculation so results are asynchronous!
                // Results are processed in watch expressions and callback functions!
                this.filterRevigo(); 

            }
        },

        // Revigo's redundancy elimination algorithm is executed asynchronously in a worker,
        // so this watcher is triggered when results are done!
        filterRevigoDone: function(newFilterRevigoDone,oldFilterRevigoDone) {
            var vm = this; 

            if (newFilterRevigoDone) {
                console.log("revigo/watch/filterRevigoDone: Revigo algorithm finished!");

                this.updateTermsData();

                // Reset the dataLoaded so that we can set it to true when needed
                // Autotriggers this dataLoaded watch but will have no effect as else clause is empty
                this.dataLoaded = false;

                let revigoEndTime = performance.now();
                vm.revigoStatus['algorithm'] = 'Revigo';
                vm.revigoStatus['start'] = false;
                vm.revigoStatus['loading'] = false; 
                vm.revigoStatus['termCount'] = vm.multipleExperimentSetting ? 
                                                    vm.experimentEnrichmentsFilteredDAG.size : 
                                                    vm.enrichmentsFilteredDAG.length;
                vm.revigoStatus['duration'] = revigoEndTime-vm.revigoStartTime;
                vm.revigoStatus['namespace'] = vm.ontologyName;
                vm.revigoStatus['result'] = '';

                if (this.visualizationContext=='initial') { 
                    
                    // Change the visualization context from initial to GO-GO
                    this.visualizationContext = 'go-go';
                }

                if (vm.visualizationContext=='go-go') {

                    // Calculation of distance matrix is automatically called in the watcher 
                    // for results from LCA!
                    vm.calculateLCA(); 

                    // TODO: Check if this would be better put in filterEnrichments, where we already have it
                    //       for a single experiment setting!
                    //       In multiple experiment setting nodeData stores information on experiments, 
                    //       not terms!
                    // TODO: Consider moving out of if (doesn't apply to EX-EX case!)
                    vm.nodeData = vm.termsData;  
                }

                if (vm.visualizationContext=='go-ex') {
                    // TODO: Consider moving out of if (doesn't apply to EX-EX case!)
                    vm.nodeData = vm.termsData; 
                }

                // TODO: For this case to work properly I first have to decouple filtering by ontology space from
                //       filtering by p-value. The SimGIC algorithm is doing selection of enriched/non-enriched
                //       terms and so it makes sense to pass all non-redundat GO terms. The redundancy will be
                //       calculated with Revigo's algorithm, on a specific ontology space.
                if (vm.visualizationContext=='ex-ex') {

                    vm.calculateSimGIC();

                    // TODO: Sets the nodeData with experiment information instead of term information
                    vm.passInformationOnExperiments();

                    if (vm.embeddingType=='hierarchical') {
                        vm.calculateHierarchicalClustering();
                    } 
                }

            }

            // Reset status variable (won't trigger this watcher because there is if clause at beggining!)
            vm.filterRevigoDone = false;

        },

        ontologyName: function(newOntologyName, oldOntologyName) {

            var vm = this; 
            console.log("revigo/watch/ontologyName: Changed value to "+newOntologyName);

            if (newOntologyName=="biological processes") {
                vm.dag = vm.dagBio;
                vm.nmfEmbedding = vm.nmfEmbeddingBio;
            } else if (newOntologyName=="molecular function") {
                vm.dag = vm.dagMol;
                vm.nmfEmbedding = vm.nmfEmbeddingMol;
            } else if (newOntologyName=="cellular component") {
                vm.dag = vm.dagCell;
                vm.nmfEmbedding = vm.nmfEmbeddingCell;
            }

            // Reset the embedding type!
            vm.resetEmbeddingType = false;

            this.filterEnrichmentsDAG();

            // Runs asynchronously so the results are processed in a separate watcher!
            this.filterRevigo(); 

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

            // NOTE: If embedding type is not defined this continues as normal, but in the nodedata watcher 
            //       in the scatter plot component the visualization will be skipped!
            
            // Update termsdata with the information from the newly selected experiment.
            // An assumption is that termsdata already contains all GO terms which exist in our DAG
            // and are lower than threshold p-value. GO terms visibility will reflect 'selected' attribute.
            this.termsData = new Map([...this.termsData].map( function(d) {
                var key = d[0];
                var value = d[1];
                var temp = vm.experimentEnrichmentsFilteredDAG.get(d[0]); 
                value['pvalue'] = temp.hasOwnProperty(newSelectedExperiment) ? 
                                  temp[newSelectedExperiment] :
                                  undefined; 

                value['selected'] = vm.experimentEnrichmentsFilteredDAG.get(d[0]) 
                                      .hasOwnProperty(newSelectedExperiment);
                return [key,value];
            }));

            vm.nodeData = vm.termsData;

        },

        // p-value threshold for SimGIC calculation
        // SimGIC is used as a similarity measure between experiments
        // TODO: How to trigger this watcher for the default p-value which is preselected!
        selectedPvalue: function(newSelectedPvalue, oldSelectedPvalue) {
            var vm = this;
            console.log('Somebody selected a new p-value:' + vm.selectedPvalue);

            if (vm.visualizationContext=='ex-ex') {

                // Calculate SimGIC with the chosen p-value
                // This will also update t-SNE, MDS and UMAP embeddings automatically if they are selected
                vm.calculateSimGIC();
                vm.passInformationOnExperiments();

                // In case hierarchical clustering is needed, calculate it here
                if (vm.embeddingType=='hierarchical') {
                    vm.calculateHierarchicalClustering();
                } 
            }

        },

        embeddingType: function(newEmbeddingType, oldEmbeddingType) {
            var vm = this;
            console.log("revigo/watch/embeddingType: Someone changed embedding type to " + vm.embeddingType);

            // Do nothing if embedding type is empty!
            if (vm.embeddingType=='') {
                return;
            }

            if (vm.visualizationContext=='ex-ex') {
                console.log('revigo/watch/embeddingType: Visualize experiment embedding in EX-EX context! ');

                vm.calculateSimGIC();
                vm.passInformationOnExperiments();

                // Hierarchical clustering is only available in multiple experiments visualization context,
                // so we don't have to explicitly check for the visualization context!
                if (vm.embeddingType=='hierarchical') {
                    vm.calculateHierarchicalClustering();
                } 

            }

            // TODO: If the embedding type changes back to GO-GO, check whether semantic similarity is
            //       already calculated and if yes, assign semanticSimilarityMatrix to the distmat which
            //       will then be passed to the scatter-plot component.

        },

        visualizationContext: function(newVisualizationContext, oldVisualizationContext) {
            let vm = this;
            console.log("revigo/watch/visualizationContext: Someone changed visualization context to "+
                vm.visualizationContext);

            // Reset the embedding type when you switch the visualization context.
            vm.embeddingType = '';

            if (vm.visualizationContext=='go-go') {

                // I need this to set results properly, because they are used in calculateDistanceMatrix!
                vm.calculateLCA();

                vm.nodeData = vm.termsData;

                // Typeset the MathJax equations in semantic similarity menu after it is revealed!
                this.$nextTick(function() {
                    MathJax.typeset();
                 });

            }

            if (vm.visualizationContext=='go-ex') {
                vm.nodeData = vm.termsData;
            }

            if (vm.visualizationContext=='ex-ex') {

                // TODO: We calculate SimGIC as soon as we switch to EX-EX visualization context!
                //       We do almost the same thing in the watcher for embedding type, check that!

                // Calculate SimGIC as soon as the context is switched!
                vm.calculateSimGIC();
                vm.passInformationOnExperiments();

            }
        },

        remaining_terms_all(newRemainingTerms,oldRemainingTerms) {
            let vm = this;
            console.log("revigo/watch/remaining_terms_all: Someone changed remaining_terms_all!");

            // NOTE: In single experiment setting the remaining_terms_all is an Array, while in
            //       multiple experiment setting it is a Map!

            // TODO: If the Map is empty, this means that we just initialized it, so we can skip the rest!
            if (vm.remaining_terms_all.size==0) {
                return;
            }

            // Single experiment setting
            if (!vm.multipleExperimentSetting) {

                // In single experiment setting there is no aggregation of non-redundant terms!
                // This is passed to the data table for visualization!
                vm.final_remaining_terms = vm.remaining_terms_all;

                // Filter enrichments based on the Revigo's algorithm!
                vm.enrichmentsFiltered = vm.enrichmentsFilteredDAG 
                                           .filter( x => vm.final_remaining_terms.includes(x[0]) );

                // Revigo runs asynchronously so we need to set status variable so that the main app
                // knows that algorithm finished.
                vm.filterRevigoDone = true;

            // Multiple experiment setting
            } else {

                if (vm.remaining_terms_all.size == vm.experiments.length) {
                    console.log('revigo/watch/remaining_terms_all: Terms from all experiments collected!');

                    // Aggregate results from all experiments by choosing all terms which are non-redundant
                    // in at least one experiment, basically taking union of all non-redundant terms accross
                    // all experiments.
                    // This is passed to the data table for visualization!
                    vm.final_remaining_terms = 
                        [...vm.remaining_terms_all].reduce((sum,x)=>sum.concat(x[1]),[]).unique();

                    // Filter the enrichments based on the final set of non-redundant terms!
                    vm.experimentEnrichmentsFiltered = 
                        new Map( [...vm.experimentEnrichmentsFilteredDAG] 
                               .filter( x => vm.final_remaining_terms.includes(x[0]) )
                        );

                    // Set the status variable so that the app knows that Revigo algorithm finished!
                    vm.filterRevigoDone = true;

                }
            }

        }

    },

    methods: {

        receiveStatusFromScatterPlot: function(statusMessage) {
            let vm = this;
            console.log("revigo/methods/receiveStatusFromScatterPlot: "+statusMessage.message);

            if (statusMessage.umapdisabled) {
                vm.umapdisabled = true;
            } else {
                vm.umapdisabled = false;
            }

        },

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

                // TODO: Hack to convert Map back to Object, ideally enrichments should be in one or the other!
                this.enrichments = Object.fromEntries([...value]);

                this.dataLoaded = true;
                this.lcaStatus['warning'] = '';

            } else {

                console.log("revigo/methods/receiveDataFromChild: Received data for multiple experiments");

                // Show experiment selector
                this.multipleExperimentSetting = true;

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
                    [x[0], 
                     Object.fromEntries( x[1].map( y => 
                        [ vm.experimentFeatures.filter( z => y[z] != '' )
                                               .map( z => y[z] )
                                               .join('_'),
                          y['enrichment']
                        ] 
                     ))
                    ]
                ));

                // TODO: Hacky solution to make multiple experiment format identical to the single experiment.
                //       I just picked first experiment's enrichment from the filtered Map.
                //       Ideally the specific experiment will be selected with the menu, with the first
                //       one selected by default. I can even have a separate function for selection.
                var enrichments = Object.fromEntries([...value].map( x=> [x[0], x[1][0]['enrichment']] ));

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

                    this.lcaStatus['warning'] = 'All zero p-values replaced with the smallest '+
                                                    'non-zero p-value!';
                }

                // // TODO: ------- MOVED FROM UPDATE TERMS DATA!
                // // BUG: PCA BIPLOT DOES NOT WORK ANYMORE!

                // // // TODO: MOVED TO filterEnrichmentsDAG!
                // // // If GO term is not present in some experiment give it p-value of 1.0
                // // // TODO: Consider removing!
                // // this.experimentDataTable = [...this.experimentEnrichmentsFiltered]
                // //     .map( x => vm.experiments.map( e => e in x[1] ? x[1][e] : 1.0 ) );  

                // // Experiment data table where all GO terms are present!
                // // TODO: This is the one we should keep!
                // this.experimentDataTableUnfiltered = [...this.experimentEnrichments]
                //     .map( x => vm.experiments.map( e => e in x[1] ? x[1][e] : 1.0 ) );  

                // // TODO: -------

                // Pass the enrichments from multiple experiments to be visualized
                this.enrichments = enrichments;
                this.dataLoaded = true;

            }
        },

        // Filter all GO terms belonging to the current ontology namespace (stored as DAG)
        filterEnrichmentsDAG: function() {
            var vm = this; 
            console.log("revigo/methods/filterEnrichmentsDAG: Selecting terms from specific ontology space!");

            if (!this.multipleExperimentSetting) {

                // TODO: NOT NEEDED ANYMORE?!
                // // Single experiment setting
                // this.enrichmentsFiltered = Object.entries(this.enrichments)
                //                                  .filter( x => this.dag.hasOwnProperty(x[0]) );

                // For use in data table where we show all GO terms from current ontology namespace,
                // regardless of further filtering with Revigo algorithm!
                this.enrichmentsFilteredDAG = Object.entries(this.enrichments)
                                                 .filter( x => this.dag.hasOwnProperty(x[0]) );

            } else {

                // TODO: NOT NEEDED ANYMORE?!
                // // Store p-values of all GO terms of all experiments - all of these will be embedded but
                // // only those belonging to the currently selected experiment will actually be visible!
                // this.experimentEnrichmentsFiltered = 
                //     new Map( [...vm.experimentEnrichments]
                //         .filter( x => vm.dag.hasOwnProperty(x[0]) ) // make sure GO term exists in our DAG
                //     );
                
                // Store p-values of all GO terms of all experiments - all of these will be embedded but
                // only those belonging to the currently selected experiment will actually be visible!
                this.experimentEnrichmentsFilteredDAG = 
                    new Map( [...vm.experimentEnrichments]
                        .filter( x => vm.dag.hasOwnProperty(x[0]) ) // make sure GO term exists in our DAG
                    );

                // // TODO: MOVED FROM receiveDataFromChild (ORIGINALY FROM updateTermsData)
                // // If GO term is not present in some experiment give it p-value of 1.0
                // // TODO: Consider removing!
                // this.experimentDataTable = [...this.experimentEnrichmentsFiltered]
                //     .map( x => vm.experiments.map( e => e in x[1] ? x[1][e] : 1.0 ) );  

            }

        },

        filterRevigo: function() {
            var vm = this; 
            console.log("revigo/methods/filterRevigo: Revigo's redundancy elimination algorithm!");

            // Reset the embedding type, except if the resetEmbeddingType is false, then reset it to true!
            if (vm.resetEmbeddingType) {
                vm.embeddingType = '';
            } else {
                vm.resetEmbeddingType = true;
            }

            // TODO: Reinitialize the Map before you run the Revigo algorithm!
            //       In single experiment case it is actually an Array, but then it will be simply overwritten!
            vm.remaining_terms_all = new Map(); 

            // End time for Revigo is taken in filterRevigoDone watcher!
            vm.revigoStartTime = performance.now();

            vm.revigoStatus['algorithm'] = 'Revigo';
            vm.revigoStatus['start'] = false;
            vm.revigoStatus['loading'] = true;
            vm.revigoStatus['termCount'] = vm.multipleExperimentSetting ? 
                                                vm.experimentEnrichmentsFilteredDAG.size : 
                                                vm.enrichmentsFilteredDAG.length;
            vm.revigoStatus['duration'] = 0;
            vm.revigoStatus['namespace'] = vm.ontologyName;
            vm.revigoStatus['result'] = '';

            // Reset the LCA progress bar so it doesn't display old information when Revigo is calculated
            vm.lcaStatus['start'] = true;

            // Single experiment setting
            if (!this.multipleExperimentSetting) {

                if (vm.revigoAlgorithmType=='LCA') {
                    this.revigoWorker.postMessage([vm.revigoAlgorithmType,
                                                   vm.enrichmentsFilteredDAG, 
                                                   vm.terms,
                                                   vm.dag]);
                } else {
                    this.revigoWorker.postMessage([vm.revigoAlgorithmType,
                                                   vm.enrichmentsFiltered, 
                                                   vm.terms,
                                                   vm.dag,
                                                   vm.nmfEmbedding]);
                }

                this.revigoWorker.onmessage = function(e) {

                    // If the output from web worker is progress then set the progress status variable.
                    if (e.data[0]=="progress") {
                        vm.revigoStatus['progress'] = ((e.data[1]/e.data[2])*100).toFixed(0);
                        return;
                    }

                    // TODO: Non-redundant terms are stored as keys!
                    // vm.remaining_terms_all = e.data; 
                    vm.remaining_terms_all = Array.from(e.data.keys());

                    // TODO: Redundancy Map, which term is represented by some other term!
                    //       Generalize this to the multiple experiment case as well!
                    vm.redundancyMap = new Map([...e.data].map(x=>x[1].map(y=>[y,x[0]]))
                                                          .reduce((arr,t)=>arr.concat(t),[]));
                };

            // Multiple experiment setting
            } else {

                // Send pairs of terms for LCA calculation to a worker, in batches for each experiment
                for (const experiment of vm.experiments) {

                    // let enrichmentsFromExperiment = [...vm.experimentEnrichmentsFiltered]
                    let enrichmentsFromExperiment = [...vm.experimentEnrichmentsFilteredDAG] // TODO: TESTING!
                        .map( x => x[1].hasOwnProperty(experiment) ? [x[0],x[1][experiment]] : null )
                        .filter( x => x !== null );

                    if (vm.revigoAlgorithmType=='LCA') {

                        // Send data to web worker that runs Revigo algorithm
                        // In multiple experiment setting we also pass the name of the experiment
                        this.revigoWorker.postMessage([vm.revigoAlgorithmType,
                                                       enrichmentsFromExperiment,
                                                       vm.terms,
                                                       vm.dag,
                                                       experiment]); // Pass the name of the experiment as well!
                    } else {

                        // Send data to web worker that runs Revigo algorithm
                        this.revigoWorker.postMessage([vm.revigoAlgorithmType,
                                                       enrichmentsFromExperiment,
                                                       vm.terms,
                                                       vm.dag,
                                                       vm.nmfEmbedding,
                                                       experiment]); // Pass the name of the experiment as well!
                    }
                }

                // Wait for the worker to return non-reduntant GO terms from Revigo algorithm
                this.revigoWorker.onmessage = function(e) {
                    
                    // If the output from web worker is progress then set the progress status variable.
                    if (e.data[0]=="progress") {
                        vm.revigoStatus['progress'] = ((e.data[2]/e.data[3])*100).toFixed(0);
                        vm.revigoStatus['experiment'] = e.data[1];
                        return;
                    }

                    // Im multiple experiment setting we get both the results and the name of the experiment
                    // TODO: Non-redundant terms are stored as keys!
                    // let remaining_terms = e.data[0];
                    let remaining_terms = Array.from(e.data[0].keys());
                    let experiment = e.data[1];

                    console.log("revigo/methods/filterRevigo: Remaining terms for the "+experiment+"...");
                    console.log(remaining_terms);

                    // A hack to force triggering a watch expression on remaining_terms_all, as Vue 2.0
                    // is not reactive on changes to Map and Set (although it is on their assignment)
                    vm.remaining_terms_all = new Map(vm.remaining_terms_all.set(experiment,remaining_terms));

                    // Results are aggregated asynchronously in remaining_terms_all watcher!

                    // Redundancy Map, which term is represented by some other term!
                    // This is a generalization to the multiple experiment case, a Map of Maps,
                    // every experiment is stored as a entry in a Map!
                    vm.redundancyMap = new Map(vm.redundancyMap.set(
                        experiment,
                        new Map([...e.data[0]].map(x=>x[1].map(y=>[y,x[0]]))
                                                          .reduce((arr,t)=>arr.concat(t),[]))
                    ));
                    // console.log(vm.redundancyMap);

                };

            }

        },

        updateTermsData: function() {
            var vm = this; 
            console.log("revigo/methods/updateTermsData: Updating terms data!");

            if (!this.multipleExperimentSetting) {

                // Create one dataset with all information on GO terms which will be passed to visuals!
                // GO terms with missing annotations are assigned one annotation by default
                this.termsData = new Map(this.enrichmentsFiltered.map( function(d,i) {
                    var key = d[0];
                    var value = d[1];
                    return [key,
                            {'pvalue': value, 
                             'frequency': vm.terms[d[0]] || 1, // assign default value to missing GO terms
                             'description': vm.names[d[0]],
                             'selected': true}
                    ];
                }));

                // TODO: We are passing terms data here but not in the multiple experiment setting?
                vm.nodeData = vm.termsData;

            } else {

                // We use a new variable experimentEnrichmentsFiltered which already has only those
                // GO terms that exist in our DAG and which are non-redundant by Revigo algorithm.
                // First experiment is selected by default!
                this.termsData = 
                    new Map( [...vm.experimentEnrichmentsFiltered] 
                        .map( x => [ x[0],
                                     {'pvalue': x[1].hasOwnProperty(vm.experiments[0]) ? 
                                                x[1][vm.experiments[0]] : undefined, 
                                      'frequency': vm.terms[x[0]] || 1,
                                      'description': vm.names[x[0]],
                                      'selected': x[1].hasOwnProperty(vm.experiments[0])} ] )
                    );

                // TODO: THESE CAN BE DEFINED EVEN EARLIER, NOT ONLY AFTER REVIGO FINISHES?!
                // BUG: MOVED TO receiveDataFromChild AND NOW THE PCA BIPLOT IS NOT WORKING ANYMORE!

                // If GO term is not present in some experiment give it p-value of 1.0
                // TODO: Consider removing!
                this.experimentDataTable = [...this.experimentEnrichmentsFiltered]
                // this.experimentDataTable = [...this.experimentEnrichmentsFilteredDAG] // TODO: TESTING!
                    .map( x => vm.experiments.map( e => e in x[1] ? x[1][e] : 1.0 ) );  

                // // Experiment data table where all GO terms are present!
                // // TODO: This is the one we should keep!
                // // this.experimentDataTableUnfiltered = [...this.experimentEnrichments]
                // this.experimentDataTableUnfiltered = [...this.experimentEnrichmentsFilteredDAG] // TODO: TESTING!
                //     .map( x => vm.experiments.map( e => e in x[1] ? x[1][e] : 1.0 ) );  

            }
            
        },

        // Fetch all needed data - DAG, term counts, enrichments
        fetchData: function() {

            // Needed because => functions have no defined this property
            var vm = this; 

            // Wait for all data to load - the visualization starts empty!
            Promise.all(["data/go-dag-molecular-function.json",
               "data/go-dag-cellular-component.json",
               "data/go-dag-biological-process.json",
               "data/go-terms-count-goa.json",
               "data/go-names.json",
               "data/nmf_embedding_biological_process.json",
               "data/nmf_embedding_molecular_function.json",
               "data/nmf_embedding_cellular_component.json"].map(url=>vm.getUrl(url))) 
               .then(([dagMol,dagCell,dagBio,terms,names,nmfEmbeddingBio,nmfEmbeddingMol,nmfEmbeddingCell]) => { 
                    vm.dag = dagBio;
                    vm.dagBio = dagBio;
                    vm.dagMol = dagMol;
                    vm.dagCell = dagCell;
                    vm.terms = terms; 
                    vm.names = names;
                    vm.nmfEmbedding = nmfEmbeddingBio;
                    vm.nmfEmbeddingBio = nmfEmbeddingBio;
                    vm.nmfEmbeddingMol = nmfEmbeddingMol;
                    vm.nmfEmbeddingCell = nmfEmbeddingCell;
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

            vm.lcaStatus = {'algorithm': 'LCA',
                            'start':false,
                            'loading':true, 
                            'termCount':vm.termsData.size, 
                            'duration': 0,
                            'namespace':vm.ontologyName,
                            'result': '',
                            'warning':vm.lcaStatus['warning']};

            console.log('revigo/calculateLCA: Sent pairs and dag to worker!');

            // Wait for the worker to return LCA results
            this.lcaWorker.onmessage = function(e) {
                
                let t1 = performance.now();

                vm.lcaStatus = {'algorithm': 'LCA',
                                'start':false,
                                'loading':false, 
                                'termCount':vm.termsData.size, 
                                'duration': t1-t0,
                                'namespace':vm.ontologyName,
                                'result': '',
                                'warning':vm.lcaStatus['warning']};

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
            let vm = this;

            let similarity = this.semanticSimilarity(this.results,this.terms,this.similarityMeasure);

            // List of GO terms for visualization
            let goTerms = similarity.map(x=>[x[0],x[1]]).flat().unique()

            // Default max value which we use for t-sne distance matrix 
            // Custom max function because Math.max is recursive and fails for large arrays
            let maxValue = similarity.map(x=>x[2]).reduce((max, v) => max >= v ? max : v, -Infinity);
            vm.semanticSimilarityMatrix = [...Array(goTerms.length)].map(e=>Array(goTerms.length)
                                                                    .fill(maxValue));
            for (const x of similarity) {
                let x0 = this.termsToId.get(x[0]);
                let x1 = this.termsToId.get(x[1]);
                vm.semanticSimilarityMatrix[x0][x1] = x[2];
                vm.semanticSimilarityMatrix[x1][x0] = x[2];
            }
            // Final distance matrix passed to child components is assigned only once!
            this.distmat = vm.semanticSimilarityMatrix;
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

        },

        calculateSimGIC: function() {
            var vm = this;
            console.log("Calculating SimGIC with threshold p-value: " + vm.selectedPvalue);

            // Formula for SimGIC:
            //
            // Information content: IC(t) = -log(p(t))
            // Probability of a term: p(t) = n(t)/N
            //
            // SimGIC(p,q) is then simply the ratio of:
            // - the sum of IC's for terms in the intersection of the GO sets p and q
            // - and the sum of the IC's for terms in the union of the GO sets p and q.

            // experimentEnrichments is a Map which contains all GO terms and enrichments for their experiments
            // TODO: Zero p-values are replaced with the smallest non-zero p-value
            //       Some GO terms have missing experiments! Give them value of 1.0 there!
            let goterms = [...vm.experimentEnrichments.keys()];

            let annotations = new Map(Object.entries(vm.terms));

            // We are only interested in annotations of GO terms that actually appear in our dataset
            // Those that do not appear get the default value of one annotation!
            let annotationsFiltered = goterms.map( x => annotations.get(x) || 1);

            let totalAnnotations = [...annotationsFiltered].reduce((sum,x)=>sum+x,0);

            // Information content of each GO term in multiple experiment dataset
            let infoContentGO = [...annotationsFiltered].map( (x) => -Math.log((x)/totalAnnotations) );

            // Experiment x experiment intersection, union and similarity matrices
            let intersectionMatrix = [...Array(vm.experiments.length)]
                                            .map(e=>Array(vm.experiments.length)
                                            .fill(0));
            let unionMatrix = [...Array(vm.experiments.length)]
                                            .map(e=>Array(vm.experiments.length)
                                            .fill(0));
            vm.simGICmatrix = [...Array(vm.experiments.length)]
                                          .map(e=>Array(vm.experiments.length)
                                          .fill(0));

            // Marks the significant GO terms in each experiment
            // We use experimentDataTableUnfiltered to determine which GO terms are significant or
            // not significant in certain experiments, depending on their p-value.
            // let temp = zip(vm.experimentDataTableUnfiltered,infoContentGO).map( function(x) {
            let temp = zip([...this.experimentEnrichmentsFilteredDAG],infoContentGO).map( function(x) {
                return [x[0].map( function(y) {
                    return y <= vm.selectedPvalue;
                }),x[1]];
            });

            // Calculate intersection and union matrices
            temp.forEach( function(x) {
                for (i = 0; i < vm.experiments.length; i++) {
                    for (j = 0; j < vm.experiments.length; j++) {
                       if (x[0][i]||x[0][j]) {
                           unionMatrix[i][j] += x[1];
                       }
                       if (x[0][i]&&x[0][j]) {
                           intersectionMatrix[i][j] += x[1];
                       }
                    }
                }
            });

            // Divide intersection and union matrix elementwise to get SimGIC
            for (i = 0; i < vm.experiments.length; i++) {
                for (j = 0; j < vm.experiments.length; j++) {
                   vm.simGICmatrix[i][j] = intersectionMatrix[i][j] / unionMatrix[i][j];
                }
            }

            // TODO: I decided to pass information to components in a separate function!

        },
        
        passInformationOnExperiments: function() {
            var vm = this;
            console.log("revigo/methods/passInformationOnExperiments: Passing information on experiments!");

            // We will embbed experiments in the same way that we embbed GO terms (t-SNE, MDS, UMAP)
            // so we can use the same variables which we used when embbeding GO terms to pass this
            // information to the scatterplot component.

            // Pass the SimGIC matrix through the distance matrix to the scatterplot component
            vm.distmat = vm.simGICmatrix;

            // Store experiment data in a separate variable
            vm.experimentData = new Map(vm.experiments.map( function(d,i) {
                return [
                    d,
                    {'pvalue': vm.selectedPvalue,  // TODO: Makes no sense in experiment setting!
                     'frequency': 1, // TODO: Makes no sense in experiment setting!
                     'count': [...vm.experimentEnrichments].reduce(
                                function(x0,x1) {
                                    if (x1[1].hasOwnProperty(d)) {return x0+1;} 
                                    else{ return x0;}
                                },0),  // How many GO terms are in each experiment
                     'selected': true} 
                    ];
            }));

            // We are using nodedata to pass node annotation data to scatter plot component!
            vm.nodeData = vm.experimentData;

        },

        calculateHierarchicalClustering: function() {
            var vm = this;
            console.log("Calculating hierarchical clustering with threshold p-value: " + vm.selectedPvalue);

            // Hierarchical clustering with Agnes algorithm (agglomerative nesting)
            // Part of the ml-hclust library https://github.com/mljs/hclust
            let tree = agnes(vm.simGICmatrix,{method:'ward',isDistanceMatrix:true});

            // Add names to leaf nodes which correspond to experiments
            // TODO: You can also name non-leaf nodes (clusters) but for now I keep it empty!
            tree.traverse(function(x){ 
                if (x.isLeaf) { 
                    x.name = vm.experiments[x.index]; 
                } else { 
                    x.name = ""} 
            });

            vm.experimentHierarchy = tree;

        }

    }

});

Vue.component('progress-box', {
    // Fields in loadingStatus are 'algorithm', 'start', 'loading', 'progress' (only for Revigo), 
    //                             'experiment', 'termCount', 'duration', 'namespace', 'result', 'warning'.
    props: ["loadingstatus"], 
    data: function() {
        return;
    },
    computed: {
        message: function() {
            if (this.loadingstatus['start']) {
                return '';
            } else if (this.loadingstatus['loading']) {
                return 'Calculating '+this.loadingstatus['algorithm']+
                       ' between '+this.loadingstatus['termCount']+
                       ' '+this.loadingstatus['namespace']+
                       ((this.loadingstatus['algorithm']=='Revigo') ?
                            ' GO terms.</br>' + this.loadingstatus['progress'] + '% redundant terms removed' +
                            ((this.loadingstatus['experiment']!=='') ? 
                                    ' for experiment ' + this.loadingstatus['experiment'] + '.' :
                                    '.') :
                            ' GO terms <img src="static/ajax-loader.gif">');
             } else if (!this.loadingstatus['loading']) {
                return "Calculated "+this.loadingstatus['algorithm']+
                       " between "+this.loadingstatus['termCount']+
                       " "+this.loadingstatus['namespace']+
                       " GO terms in "+
                       this.loadingstatus['duration'].toFixed(0)+" miliseconds!"+
                       '</br><span style="color:red">'+
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
            console.log('CSV file loaded!');

            var vm = this;

            // TODO: We expect only one file, but in theory we can load multiple files as well!
            var file = e.target.files[0];

            // Hack to reset the value of the file input element so that @change event will be
            // triggered even if we decide to load the same file twice in a row (e.g. after file change).
            e.target.value = "";

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

    props: ["lcastatus",
            "distmat",
            "nodedata",
            "embeddingtype",
            "visualizationcontext",
            "experimentdatatable",
            "experiments",
            "experimenthierarchy"], 

    data: function () {

        return {
            canvas: null,
            context: null,
            width: 500, 
            height: 500, 
            margin: {top: 15, right: 15, bottom: 15, left: 15}, // TODO: We also have legend margin!
            intervalID: null,
            intervalCheckID: null,
            custom: null,
            qtree: null,
            scaleX: null,
            scaleY: null,
            scaleR: null,
            coordToData: null,
            coordinatesReady: false,
            loadings: null, // PCA loadings for the PCA biplot
            nNeighbors: 3, // for UMAP
            nEpochs: 400 // for UMAP
        }

    },

    mounted: function() {

        // TODO: Consider putting this in separate css file!
        // Needed so that tooltips (which have absolute position) are positioned correctly
        d3.select(this.$el)
          .style('position','relative')
          .style('width','fit-content')
          .style('width','-moz-fit-content'); // For Firefox!

        // NOTE: lower() is used to position element as the first child of its parent.
        //       As we have a form for choosing embedding type in our scatter plot component
        //       we want all elements to be positioned above it.

        // Prepare canvas
        this.canvas = d3.select(this.$el)
                        .append('canvas')
                        .lower()
                        .classed('mainCanvas', true)
                        .attr('width',this.width)
                        .attr('height',this.height);

		var customBase = document.createElement('custom');
		this.custom = d3.select(customBase); 

        // Set margins for scatter plot
        this.scaleX = d3.scaleLinear().range([this.margin.left, this.width-this.margin.right]);
        this.scaleY = d3.scaleLinear().range([this.margin.top, this.height-this.margin.bottom]);

        this.scaleR = d3.scaleLog().range([4, 10]);
        this.scaleP = d3.scaleLog().range(['red','steelblue']);
        
        d3.select(this.$el).append("div")
                           .lower()
                           .attr("id","tooltip");

        // Legend for p-values
        d3.select(this.$el).append("div")
                           .lower()
                           .attr("id","legendP")
                           .style("float","right")

        // Legend for radius 
        d3.select(this.$el).append("div")
                           .lower()
                           .attr("id","legendR")
                           .style("float","right")

        // Put a placeholder text for empty canvas!

        this.context = this.canvas.node().getContext("2d");
        this.context.clearRect(0, 0, this.width, this.height);

        this.context.font = "24px Helvetica";
        this.context.fillStyle = "#777";
        this.context.fillText('Please load enrichment data! :-)',0.2*this.width,0.5*this.height);

    },
    watch: {
        lcastatus() {
            console.log('scatter-plot/watch/lcastatus: LCA is being calculated, better clean the plot!');
            let vm = this;
            
            vm.cleanPlot();

        },
        // Trigger when nodedata change.
        // For now we assume that node positions did not change, only node properties (enrichment-color,
        // selected-transparency) so we don't recalculate the whole embedding.
        // This watch is not triggered if we modify just the pvalue in nodedata - we have to reassign
        // the whole variable! 
        nodedata(newNodedata,oldNodedata) {
            // console.log('scatter-plot/watch/nodedata: nodedata changed!'); // Fires too much times!
            let vm = this;

            // If the embedding type is not defined do not continue with the visualization!
            if (vm.embeddingtype=='') {
                return;
            }

            // If we are visualizing GO terms or experiments in scatter plot layout

            // Initialize domains for the p-value and radius legends
            // This is needed only when we visualize scatter plot of GO terms!
            // TODO: Do we need this when we visualize scatter plot of experiments as well?
            if (['go-go','go-ex'].includes(vm.visualizationcontext)) {

                // Set scale for radius which depends on the frequency of newly loaded GO terms
                // NOTE: We include even the terms that are not part of the currently shown experiment!
                fmax = Math.max(...Array.from(this.nodedata.values()).map(x=>x.frequency));
                fmin = Math.min(...Array.from(this.nodedata.values()).map(x=>x.frequency));
                vm.scaleR.domain([fmin, fmax]);

                // Set scale for color which depends on the pvalues of newly loaded GO terms
                // NOTE: We do not include terms which are not part of the currently shown experiment!
                pmax = Math.max(...Array.from(this.nodedata.values()).filter(x=>x.selected).map(x=>x.pvalue));
                pmin = Math.min(...Array.from(this.nodedata.values()).filter(x=>x.selected).map(x=>x.pvalue));
                vm.scaleP.domain([pmin, pmax]);

                // TODO: MOVED HERE FROM DRAW NODES!
                // Draw legend for p-values (node color)
                vm.drawPlegend("#legendP", vm.scaleP, "p-value");
                // Draw legend for GOA annotations (circle radius)
                vm.drawRlegend("#legendR", vm.scaleR, "Annotations");

            }

            // TODO: Initialize domains for p-value and radius scales but using something meaningful
            //       for the experiments! Also, pass the legend label!
            // if (['ex-ex'].includes(vm.visualizationcontext)) {
                // // Set scale for radius which depends on the frequency of newly loaded GO terms
                // fmax = Math.max(...Array.from(this.nodedata.values()).map(x=>x.count));
                // fmin = Math.min(...Array.from(this.nodedata.values()).map(x=>x.count));
                // vm.scaleR.domain([fmin, fmax]);
            // }

            // Bind data to visual elements and redraw nodes - embedding did not change so node
            // positions do not change as well, only their properties (color, visibility).

            // Bind data to visual elements
            vm.bindDataToVisuals();

            // Draw data that was bound to visual elements
            vm.drawNodes(vm.canvas);

            // If PCA biplot is selected as embedding you should make sure that arrows are drawn at the end!
            // Without this the last call to drawNodes would erase the arrows from the plot!
            if (vm.embeddingtype=='biplot') {
                console.log('scatter-plot/watch/nodedata: Drawing PCA biplot arrows!');

                // TODO: Calculation of loadings is in the drawCanvas method which is not invoked here!
                
                this.loadings.forEach(function(x){
                    vm.drawLineWithArrows(vm.context,
                                          0.5*vm.width,0.5*vm.height,
                                          vm.scaleX(x[0]),vm.scaleY(x[1]),
                                          5,8,false,true);
                });

                // Add labels for arrows loadings!
                zip(this.loadings,this.experiments).forEach(function(x){
                    vm.context.font = "10px Arial";
                    vm.context.fillStyle = "black";
                    vm.context.fillText(x[1], vm.scaleX(x[0][0])+5, vm.scaleY(x[0][1])+5);
                });
            }
        },

        distmat(newDataLoaded, oldDataLoaded) {
            console.log("scatter-plot/watch/distmat: Distmat changed, will check for nodedata as well!");
            let vm = this;

            // Check if there is a potential error with UMAP calculation due to small number of data points.
            // This is emitted back to the main component through the @status attribute in the scatter-plot
            // element, whose change triggers receiveStatusFromScatterPlot method in the main component.
            if (vm.distmat.length<=vm.nNeighbors) {
                this.$emit('status',{message: "Potential error if you choose UMAP as embedding!", 
                                     umapdisabled:true});
            } else {
                this.$emit('status',{message: "There is no potential error if you choose UMAP as embedding!",
                                     umapdisabled:false});
            }

            // Only when both distmat and nodedata are loaded we proceed with the watch expression
            // This avoids having a separate computed property which checks whether both nodedata and
            // distmat changed - this expression checks that both are properly loaded.
            // TODO: Not entirely correct - if both old and new data have same number
            //       of GO terms this will pass but the GO terms will not match!
            if (this.distmat.length!=0 && this.nodedata.size!=0 &&
                this.distmat.length==this.nodedata.size) {

                console.log('scatter-plot/watch/distmat: nodedata also changed, trigger redrawing!');

                // Redraw the whole scatter plot - recalculate embedding and redraw all nodes
                this.drawCanvas();

            }
        },

        embeddingtype(newEmbeddingType,oldEmbeddingType) {
            console.log("scatter-plot/watch/embeddingType: Embedding type changed, redrawing canvas!");
            let vm = this;

            // If we are drawing PCA biplot we need experimentdatatable and nodedata to correspond!
            // if (vm.embeddingtype=='biplot') {
            if (['biplot','network'].includes(vm.embeddingtype)) { // TODO: FOR NETWORK AND BIPLOT LAYOUT!

                if (this.experimentdatatable!=0 && this.nodedata.size!=0 &&
                    this.experimentdatatable.length==this.nodedata.size) {

                    console.log("scatter-plot/watch/embeddingType: All conditions for PCA biplot satisfied!");
                    vm.drawCanvas();

                }
            }

            // Check that distance matrix and data on terms is defined and that they are of equal size!
            if (this.distmat.length!=0 && this.nodedata.size!=0 &&
                this.distmat.length==this.nodedata.size) {

                vm.drawCanvas();
            }


        },

        visualizationcontext(newVisualizationContext,oldVisualizationContext) {
            console.log("scatter-plot/watch/visualizationcontext: Visualization context changed!");
            let vm = this;

            vm.cleanPlot();

        },

        experimentdatatable(newExperimentDataTable,oldExperimentDataTable) {
            let vm = this;
            console.log("scatter-plot/watch/experimentdatatable: Experiment data table changed!");

            vm.drawCanvas();

        },

        experimenthierarchy(newExperimentHierarchy,oldExperimentHierarchy) {
            console.log("scatter-plot/watch/experimentHierarchy: Someone changed experiment hierarchy!");
            let vm = this;

            vm.cleanPlot();

            vm.drawExperimentHierarchy();

        }

    },
    methods: {

        // TODO: This function is technically not drawing canvas but just calculating all elements
        //       and their properties which have to be drawn. Actual drawing is done in drawNodes()!
        // TODO: Also, t-SNE embedding is calculated here, although it would make sense to separate it
        //       to a separate function in order to be able to call different embeddings.
        drawCanvas: function() {
            let vm = this;

            // Stop any previous setInterval() function
            console.log("scatter-plot/methods/drawCanvas: Clearing setInterval with ID: "+this.intervalID);
            clearInterval(this.intervalID);

            // TODO: Not sure whether cleaning of all plots should be here or somewhere else?
            vm.cleanPlot();

            // Stop the tsne calculation after 100 seconds
            // TODO: Check for the convergence of coordinates, rather than some specified time!
            //       You can get current coordinates with tsne.getSolution()
            clearInterval(this.intervalCheckID);
            this.intervalCheckID = setInterval(function() {
                clearInterval(vm.intervalID);
            },10000);

            if (vm.embeddingtype=='mds') {
                console.log("scatter-plot/methods/drawCanvas: Calculating MDS embedding!");

                var Y = mds.classic(this.distmat);
                vm.updateCoordinates(Y);

            } else if (vm.embeddingtype=='tsne') {
                console.log("scatter-plot/methods/drawCanvas: Calculating t-SNE embedding!");

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

                // Returns control to the browser so that canvas can be redrawn
                // Time interval is set to 0, with no delay between redrawing
                this.intervalID = setInterval(function() {
                    for(var k = 0; k < 1; k++) { 
                        tsne.step(); // every time you call this, solution gets better
                    }
                    var Y = tsne.getSolution();
                    vm.updateCoordinates(Y);

                },0);

            } else if (vm.embeddingtype=='umap') {
                console.log("scatter-plot/methods/drawCanvas: Calculating UMAP embedding!");

                // TODO: If the number of data points is very small - equal or less then nNeighbors,
                //       we will not be able to calculate UMAP. This is not likely to happen with GO terms
                //       but it may happen with the experiments whose numbers are typically much smaller!
                //       We should check this explicitly and report error to the user!

                // For each GO term find k most similar terms and save their similarity and index
                let neighborsDistances = this.distmat.map(function(row,i){
                             return row.map(function(x,i){
                                 return {'val':x,'ind':i}})
                               .sort(function(x, y){return x.val > y.val ? 1 : x.val == y.val ? 0 : -1})
                               .slice(0,vm.nNeighbors)
                });

                let umap = new UMAP({
                  nComponents: 2,
                  nEpochs: vm.nEpochs,
                  nNeighbors: vm.nNeighbors,
                });

                // Indices of neirest neighbors
                let knnIndices = neighborsDistances.map(function(x,i){
                    return x.map(function(x,i){
                        return x['ind']});
                });

                // Distances to the neirest neighbors
                let knnDistances = neighborsDistances.map(function(x,i){
                    return x.map(function(x,i){
                        return x['val']});
                });

                // Dummy data which we need to initialize UMAP, but will never actually be used
                // as we use data on neirest neighbors instead!
                let data = new Array(neighborsDistances.length).fill([1,2,3,4,2,4]);

                // Set precomputed knn indices and distances and initialize UMAP with dummy data
                umap.setPrecomputedKNN(knnIndices,knnDistances);

                // Actually, the button for UMAP is disabled if there are fewer data points
                // (either GO terms or experiments) then nNeighbors, so user will not be able to select
                // calculation of UMAP at all!
                try {
                    umap.initializeFit(data);
                } 
                catch {
                    console.log("Error - There is too few data points for calculation of UMAP! "+
                                "Need more than "+vm.nNeighbors+"!");
                    this.$emit('status',"Error - There is too few data points for calculation of UMAP! "+
                                "Need more than "+vm.nNeighbors+"!");
                    return;
                }

                // Initial iterations before we start dynamic visualization
                for (let i = 0; i < 10; i++) {
                  umap.step();
                }

                // Returns control to the browser so that canvas can be redrawn
                // Time interval is set to 0, with no delay between redrawing
                // NOTE: The same intervalID is used for t-sne embedding as well, but they
                //       are never rendered at the same time so I guess this is ok!
                this.intervalID = setInterval(function() {
                    for(var k = 0; k < 1; k++) { 
                        umap.step(); // every time you call this, solution gets better
                    }
                    let Y = umap.getEmbedding();
                    vm.updateCoordinates(Y);

                },0);

            } else if (vm.embeddingtype=='biplot') {
                console.log("scatter-plot/methods/drawCanvas: Calculating PCA biplot embedding!");

                // TODO: Experiment data table is calculated in the main app!

                // TODO: revigo.js:1100 Uncaught TypeError: Cannot read property '0' of undefined
                //       This happens even if currently selected embedding type is not biplot?!
                // PCA analysis - https://github.com/bitanath/pca
                let vectors = PCA.getEigenVectors(this.experimentdatatable);

                this.loadings = zip(vectors[0].vector,vectors[1].vector);

                // Adjusted data - GO terms projected to PCA components
                let adData = PCA.computeAdjustedData(this.experimentdatatable,vectors[0],vectors[1]);

                // Coordinates need to be in [[x,y],...] and not in [[x,...],[y,...]]
                var Y = zip(adData.adjustedData[0],adData.adjustedData[1]); 
                vm.updateCoordinates(Y);

                // TODO: Ideally, arrows for the biplot would be drawn here, but then they are overdrawn
                //       by the drawNodes in the nodedata watcher! So arrows are drawn there!
                //       The PCA loadings need to be accessible within the whole component for this to work!

            } else if (vm.embeddingtype=='network') {
                console.log("scatter-plot/methods/drawCanvas: Calculating network embedding!");

                // TODO: Experiment data table is calculated in the main app!
                console.log(this.experimentdatatable);

                // TODO: FINISH THIS! CALCULATE CORRELATION MATRIX AND DRAW NETWORK!

                // // Coordinates need to be in [[x,y],...] and not in [[x,...],[y,...]]
                // var Y = zip(adData.adjustedData[0],adData.adjustedData[1]); 
                // vm.updateCoordinates(Y);

            }

            // TODO: Maybe separate this is another function? This makes sense if we will reuse code for
            //       tooltip in different contexts, for example when drawing scatter plot of experiments.
            d3.select('.mainCanvas').on('mousemove', function() {

                var xy = d3.mouse(this);

                // Finding the closest node based on the coordinates of the mouse is the only 
                // purpose of the coordToData, as it uses node coordinates as keys.

                // TODO: Sometimes we get an error here "Cannot read property 'find' of null"
                var xyTooltip = vm.qtree.find(xy[0],xy[1]);
                var closestNode = vm.coordToData.get(String(xyTooltip));

                let closestNodeRadius = ['go-go','go-ex'].includes(vm.visualizationcontext) ?
                                            vm.scaleR(closestNode.frequency) :
                                            10;

                // If mouse cursor is close enough to the closest point show tooltip
                if (Math.abs(xy[0]-xyTooltip[0])<=closestNodeRadius && 
                    Math.abs(xy[1]-xyTooltip[1])<=closestNodeRadius) {

                    // Show the tooltip only when there is nodeData found by the mouse
                    d3.select('#tooltip')
                      .style('opacity', 0.8)
                      .style('top', xy[1]+1+'px') 
                      .style('left', xy[0]+1+'px') 
                      .html(function() { 
                          if (vm.visualizationcontext=='ex-ex') {
                              return closestNode.name+
                                     "</br>Count: "+closestNode.count;
                          } else {
                              return closestNode.name+
                                     "</br>GOA annotations: "+closestNode.frequency+
                                     "</br>description: "+closestNode.description+
                                     ((typeof closestNode.pvalue !== 'undefined') ? 
                                         "</br>p-value: "+closestNode.pvalue : 
                                         ""); 
                          }
                      });
                } else {
                    // Hide the tooltip when there our mouse doesn't find nodeData
                    d3.select('#tooltip')
                      .style('opacity', 0);
                }
            }); 

        },

        updateCoordinates: function(Y) {

            let vm = this;

            let goterms = Array.from(vm.nodedata.keys());

            let Y0min = Math.min(...Y.map(x=>x[0]));
            let Y0max = Math.max(...Y.map(x=>x[0]));
            let Y1min = Math.min(...Y.map(x=>x[1]));
            let Y1max = Math.max(...Y.map(x=>x[1]));

            vm.scaleX.domain([Y0min, Y0max]);
            vm.scaleY.domain([Y1min, Y1max]);

            vm.qtree = d3.quadtree()
              .addAll(Y.map(d=>[vm.scaleX(d[0]),vm.scaleY(d[1])]));

            // Connect point coordinates with the information on the corresponding GO term
            // This is only used to identify node closest to the mouse pointer with quad tree
            vm.coordToData = new Map(Y.map( function(d,i) {
                
                // We have to stringify coordinates to use them as key, arrays won't work!
                let key = String([vm.scaleX(d[0]),vm.scaleY(d[1])]);

                // value inherits all fields from nodedata, including "frequency" and "pvalue"

                // TODO: Sometimes we get an error "Cannot set property 'name' of undefined"
                //       This happens because nodedata updates while goterms still has old value
                //       Check if this is still an issue!

                let value = vm.nodedata.get(goterms[i]);
                value.name = goterms[i];
                value.x = Y[i][0];
                value.y = Y[i][1];
                return [key,value];

            }));

            // Update nodedata with the node coordinates. We will use it in binding data to visuals.
            // Note that that this will trigger nodedata watcher!
            vm.nodedata = new Map([...vm.nodedata].map( function(d) {
                var key = d[0];
                var value = d[1];
                var indexOfGOterm = goterms.indexOf(d[0]);
                value['x'] = Y[indexOfGOterm][0]; // TODO: Warning - Cannot read property '0' of undefined!
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
                .data(Array.from(vm.nodedata.values())) 
                .enter()
                .append('custom')
                .attr('class', 'circle')
                .attr('x', function(d, i) { return vm.scaleX(d.x); })
                .attr('y', function(d, i) { return vm.scaleY(d.y); })
                .attr('radius', function(d){ return ['go-go','ex-ex'].includes(vm.visualizationcontext) ? vm.scaleR(d.frequency) : 10; }) // .attr('radius', function(d){ return vm.scaleR(d.frequency); })
                .attr('label',function(d){ return d.name; }) 
                .attr('strokeStyle', function(d) { 
                    return d.selected ? 
                           'rgba(0,0,0,1)' : 
                           'rgba(0,0,0,'+String(selectedAlpha)+')'})
                .attr('fillStyle', function(d) { 
                    return ['go-go','go-ex'].includes(vm.visualizationcontext) ?
                                d.selected ? 
                                    vm.scaleP(d.pvalue) : 
                                    'rgba(255,255,255,'+String(selectedAlpha)+')' : 
                                    'rgba(255,255,255,'+String(selectedAlpha)+')'});
                    // return d.selected ? 
                           // vm.scaleP(d.pvalue) : 
                           // 'rgba(255,255,255,'+String(selectedAlpha)+')'});
            
        },

        // Draw data that was bound to visual elements with D3
        drawNodes: function(canvas) {

            let vm = this;

            // // TODO: MOVED WHERE WE DEFINE LEGEND SCALES IN NODEDATA WATCHER!
            // //       CONSIDER WHETHER WE NEED IT IN UPDATECOORDINATES AS WELL!
            // // p-value and annotations legend only makes sense when nodes are GO terms, not experiments!
            // if (['go-go','go-ex'].includes(vm.visualizationcontext)) {

            //     // Draw legend for p-values (node color)
            //     vm.drawPlegend("#legendP", vm.scaleP, "p-value");


            //     // Draw legend for GOA annotations (circle radius)
            //     vm.drawRlegend("#legendR", vm.scaleR, "Annotations");
            // }


            // Check whether coordinates of the nodes are ready
            if (vm.coordinatesReady) {

                this.context.clearRect(0, 0, vm.width, vm.height);

                var elements = vm.custom.selectAll('custom.circle');

                elements.remove();

                // Draw all nodes
                elements.each(function(d,i) { // for each virtual/custom element...
                    var node = d3.select(this);
                    vm.context.beginPath();
                    vm.context.fillStyle = node.attr('fillStyle'); 
                    vm.context.strokeStyle = node.attr('strokeStyle'); 
                    vm.context.arc(node.attr('x'), node.attr('y'), node.attr('radius'),0,2*Math.PI) 
                    vm.context.fill();
                    vm.context.stroke();
                });

            }

        },

        // Generate legend with continuous colors from a prespecified scale
        // Strangely, this is not built in D3 by default!?
        // http://bl.ocks.org/syntagmatic/e8ccca52559796be775553b467593a9f
        drawPlegend: function(selector_id, colorscale, legendLabel) {

          var legendheight = 200,
              legendwidth = 80,
              margin = {top: 25, right: 60, bottom: 10, left: 2};
          
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
          // http://bl.ocks.org/mbostock/048d21cf747371b11884f75ad896e5a5
          // http://stackoverflow.com/questions/4899799
          //       /whats-the-best-way-to-set-a-single-pixel-in-an-html5-canvas   
          d3.range(legendheight).forEach(function(i) {
            ctx.fillStyle = colorscale(legendscale.invert(i));
            ctx.fillRect(0,i,1,1);
          });

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
            .text(legendLabel); 

        },

        // https://www.youtube.com/watch?v=XmVPHq4NhMA
        drawRlegend: function(selector_id, sizescale, legendLabel) {

          var legendheight = 150,
              legendwidth = 110,
              margin = {top: 30, right: 100, bottom: 10, left: 2};

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
                .text(legendLabel); 

            });
        },

        drawExperimentHierarchy: function() {

            console.log("scatter-plot/methods/drawExperimentHierarchy: Drawing experiment hierarchy!");
            let vm = this;

            var width = 500,
                height = 500,
                nodeRadius = 4.5;

            let outerRadius = width / 2;
            let innerRadius = outerRadius - 170;

            // TODO: Use this to delete svg visualization when data changes!
            // d3.select('#experimentHierarchy').selectAll("*").remove();

            var svg = d3.select(this.$el)
                .append("svg")
                .lower()
                .attr('id','experimentHierarchy')
                .attr("viewBox", [-outerRadius, -outerRadius, width, width])
                .attr("font-family", "sans-serif")
                .attr("font-size", 10)
                .attr('width',width)
                .attr('height',height)
                .style('position','absolute');

            let root = d3.hierarchy(vm.experimenthierarchy, d => d.children)
                  .sum(d => d.children ? 0 : 1)
                  .sort((a, b) => (a.value - b.value) || d3.ascending(a.data.length, b.data.length));

            let cluster = d3.cluster()
                .size([360, innerRadius])
                .separation((a, b) => 1);

            cluster(root);

            // Compute the maximum cumulative length of any node in the tree.
            function maxLength(d) {
              return d.data.length + (d.children ? d3.max(d.children, maxLength) : 0);
            };

            // Set the radius of each node by recursively summing and scaling the distance from the root.
            function setRadius(d, y0, k) {
              d.radius = (y0 += d.data.length) * k;
              if (d.children) d.children.forEach(d => setRadius(d, y0, k));
            };

            let color = d3.scaleOrdinal() // .domain(["Bacteria", "Eukaryota", "Archaea"])
                .range(d3.schemeCategory10)

            // Set the color of each node by recursively inheriting.
            function setColor(d) {
              var name = d.data.name;
              d.color = color.domain().indexOf(name) >= 0 ? color(name) : d.parent ? d.parent.color : null;
              if (d.children) d.children.forEach(setColor);
            };

            setRadius(root, root.data.length = 0, innerRadius / maxLength(root));
            setColor(root);

            function linkVariable(d) {
              return linkStep(d.source.x, d.source.radius, d.target.x, d.target.radius);
            }

            function linkConstant(d) {
              return linkStep(d.source.x, d.source.y, d.target.x, d.target.y);
            }

            function linkExtensionVariable(d) {
              return linkStep(d.target.x, d.target.radius, d.target.x, innerRadius);
            }

            function linkExtensionConstant(d) {
              return linkStep(d.target.x, d.target.y, d.target.x, innerRadius);
            }

            function linkStep(startAngle, startRadius, endAngle, endRadius) {
              const c0 = Math.cos(startAngle = (startAngle - 90) / 180 * Math.PI);
              const s0 = Math.sin(startAngle);
              const c1 = Math.cos(endAngle = (endAngle - 90) / 180 * Math.PI);
              const s1 = Math.sin(endAngle);
              return "M" + startRadius * c0 + "," + startRadius * s0
                  + (endAngle === startAngle ? "" : "A" + startRadius + "," + startRadius + " 0 0 " + (endAngle > startAngle ? 1 : 0) + " " + startRadius * c1 + "," + startRadius * s1)
                  + "L" + endRadius * c1 + "," + endRadius * s1;
            }

              svg.append("style").text(`

            .link--active {
              stroke: #000 !important;
              stroke-width: 1.5px;
            }

            .link-extension--active {
              stroke-opacity: .6;
            }

            .label--active {
              font-weight: bold;
            }

            `);

              const linkExtension = svg.append("g")
                  .attr("fill", "none")
                  .attr("stroke", "#000")
                  .attr("stroke-opacity", 0.25)
                .selectAll("path")
                .data(root.links().filter(d => !d.target.children))
                .join("path")
                  .each(function(d) { d.target.linkExtensionNode = this; })
                  .attr("d", linkExtensionConstant);

              const link = svg.append("g")
                  .attr("fill", "none")
                  .attr("stroke", "#000") // .attr("stroke-opacity", 0.25) // TODO: I added this! 
                .selectAll("path")
                .data(root.links())
                .join("path")
                  .each(function(d) { d.target.linkNode = this; })
                  .attr("d", linkConstant)
                  .attr("stroke", d => d.target.color);

              svg.append("g")
                .selectAll("text")
                .data(root.leaves())
                .join("text")
                  .attr("dy", ".31em")
                  .attr("transform", d => `rotate(${d.x - 90}) translate(${innerRadius + 4},0)${d.x < 180 ? "" : " rotate(180)"}`)
                  .attr("text-anchor", d => d.x < 180 ? "start" : "end")
                  .text(d => d.data.name.replace(/_/g, " "))
                  .on("mouseover", mouseovered(true))
                  .on("mouseout", mouseovered(false));

              function mouseovered(active) {
                return function(event, d) {
                  d3.select(this).classed("label--active", active);
                  d3.select(d.linkExtensionNode).classed("link-extension--active", active).raise();
                  do d3.select(d.linkNode).classed("link--active", active).raise();
                  while (d = d.parent);
                };
              }



        },

        // TODO: A single function for removing all plot elements from the DOM!
        //       Should be called before drawing new plot.
        //       Still not sure whether this is the best way to implement this?
        cleanPlot: function() {

            console.log("scatter-plot/methods/cleanPlot: Cleaning all plots!");
            let vm = this;

            // Reset the interval for dynamic calculation of embeddings.
            clearInterval(vm.intervalID);

            // Clean all the legends. 
            // SVG elements will remain. We create them at mount time so we cannot remove them completelly!
            d3.select('#legendP').selectAll('*').remove();
            d3.select('#legendR').selectAll('*').remove();

            // Clear the canvas! The canvas element remains as it is created only once at mount time!
            vm.context.clearRect(0, 0, vm.width, vm.height);

            // Remove SVG for the experiment hierarchy.
            // Otherwise remains on top of canvas element and messes with the mouseover listener for tooltip.
            // We can remove the whole element because it is recreated in the drawExperimentHierarchy()
            d3.select('#experimentHierarchy').remove();

        },

        // x0,y0: the line's starting point
        // x1,y1: the line's ending point
        // width: the distance the arrowhead perpendicularly extends away from the line
        // height: the distance the arrowhead extends backward from the endpoint
        // arrowStart: true/false directing to draw arrowhead at the line's starting point
        // arrowEnd: true/false directing to draw arrowhead at the line's ending point
        drawLineWithArrows: function(ctx,x0,y0,x1,y1,aWidth,aLength,arrowStart,arrowEnd){
            var dx=x1-x0;
            var dy=y1-y0;
            var angle=Math.atan2(dy,dx);
            var length=Math.sqrt(dx*dx+dy*dy);

            ctx.translate(x0,y0);
            ctx.rotate(angle);
            ctx.beginPath();
            ctx.moveTo(0,0);
            ctx.lineTo(length,0);
            if(arrowStart){
                ctx.moveTo(aLength,-aWidth);
                ctx.lineTo(0,0);
                ctx.lineTo(aLength,aWidth);
            }
            if(arrowEnd){
                ctx.moveTo(length-aLength,-aWidth);
                ctx.lineTo(length,0);
                ctx.lineTo(length-aLength,aWidth);
            }

            ctx.strokeStyle = "black"; 
            ctx.stroke();
            ctx.setTransform(1,0,0,1,0,0);
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
    props: ["distmat","ontology","similarity","visualizationcontext"], 
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
    template: '<a v-if="visualizationcontext==\'go-go\'" v-bind:href="message" download="distmat.csv">'+
              'Download {{similarityFullName}} distance matrix for {{ontology}} ({{matrixShape}})'+
              '</a>'
});

Vue.component('data-table', {

    props: ["experimentenrichments",
            "experiments",
            "multipleexperimentsetting",
            "enrichmentsfiltereddag",
            "names",
            "remainingterms",
            "representatives",
            "experiment"],

    data: function() {
        return {
            table: null
        }
    },

    mounted: function() {

        let vm = this;

        vm.scaleP = d3.scaleLog().range(['red','steelblue']);

    },

    computed: {

        // Whether both experiment enrichments and non-redundant terms are ready!
        // Watches for the changes in both remaining terms and experiment enrichments!
        dataWatcher: function() {

            // return [this.remainingterms,this.experimentenrichments,this.enrichmentsfiltereddag];

            // console.log([this.remainingterms,this.experimentenrichments,this.enrichmentsfiltereddag,this.experiment]);

            // TODO: IT IS BETTER TO SIMPLY REDRAW THE WHOLE TABLE WHEN EXPERIMENT CHANGES!  
            return [this.remainingterms,this.experimentenrichments,this.enrichmentsfiltereddag,this.experiment];

        }
    },

    watch: {

        // // TODO: IT IS BETTER TO SIMPLY REDRAW THE WHOLE TABLE WHEN EXPERIMENT CHANGES!  
        // // experiment: function(newData,oldData) {
        //     let vm = this;
        //     console.log("data-table/watch/experiment: Experiment changed, redrawing hierarchy in the table!");

        //     // TODO: When selected experiment change redraw Revigo redundancy hierarchy in the table!

        // },

        dataWatcher: function(newData,oldData) {
            let vm = this;
            console.log("data-table/watch/dataWatcher: Checking data for data table!");

            // Single experiment setting
            if (!vm.multipleexperimentsetting) {

                if (vm.remainingterms.length!=0 && vm.enrichmentsfiltereddag.length!=0) {

                    // TODO: MODIFY TABLEDATA SO THAT REDUNDANT TERMS ARE PLACED RIGHT BELLOW THEIR
                    //       NON-REDUNDAT REPRESENTATIVES!
                    //       Hmmm, unfortunatelly we don't have this connection in the remaining terms? :-(

                    // console.log(vm.representatives);

                    // let tabledata = [...vm.nodedata].map(function(d){return [d[0],d[1].pvalue]});
                    // let tabledata = vm.enrichmentsfiltereddag;

                    // TODO: Add column with representative terms!
                    // let tabledata = vm.enrichmentsfiltereddag.map(x=>[x[0],vm.names[x[0]],x[1]]);
                    // let tabledata = vm.enrichmentsfiltereddag.map(x=>[vm.representatives.get(x[0]),x[0],vm.names[x[0]],x[1]]);
                    let tabledata = vm.enrichmentsfiltereddag.map(
                        x => [vm.representatives.has(x[0]) ? vm.representatives.get(x[0]) : x[0],
                              x[0],
                              vm.names[x[0]],x[1]]
                    );

                    // Existing datatable cannot be reinitialized, so if it already exist destroy it before
                    if ($.fn.dataTable.isDataTable( this.$refs.datatable )) {
                        console.log("data-table/watch/dataWatcher: Destroying existing table!");
                        vm.table.destroy(); // Not enough, we need to empty table as well!
                        $(vm.$refs.datatable).empty();
                    }

                    // Set scale for color which depends on the pvalues of newly loaded GO terms

                    // TODO: The p-value color scale should be defined only for the non-redundant terms, 
                    //       while here we are calculating it for all terms! (Also in multiple experiment case!)
                    
                    // let pmax = Math.max(...vm.enrichmentsfiltereddag.map(x=>x[1]));
                    // let pmin = Math.min(...vm.enrichmentsfiltereddag.map(x=>x[1]));

                    let pvalues = vm.enrichmentsfiltereddag
                                    .filter(x=>vm.remainingterms.includes(x[0]))
                                    .map(x=>x[1]);
                    let pmax = Math.max(...pvalues);
                    let pmin = Math.min(...pvalues);

                    vm.scaleP.domain([pmin, pmax]);

                    vm.table = $(this.$refs.datatable).DataTable( {
                            scrollY:        "200px",
                            scrollCollapse: true,
                            paging:         false,
                            deferRender:    true, // Faster on large tables!
                            // scroller:       true,  // Requires pagination to be disabled (paging: false)! 
                            // autoWidth:      false,
                            // fixedColumns: true,
                            data: tabledata,
                            columns: [

                                { title: "Representative" , 
                                  width: "70px", 
                                  defaultContent: "",
                                  render: function(data,type,row,meta) {
                                    return "<div style='color:"+
                                            // (vm.representatives.has(row[0]) ?  "black" : "gray")+"'>"+
                                            ( (row[0] != row[1]) ? "black" : "LightGray" )+"'>"+
                                            row[0] + 
                                            ( (row[0] == row[1]) ? "*" : "" ) + 
                                            "<div>";
                                  }
                                }, 

                                { title: "GO term" , width: "70px"},
                                // { title: "description", width: "100px" },

                                { title: "description"},

                                { title: "p-value",
                                  width: "70px",
                                  className: 'dt-body-center',
                                  render: function(data,type,row,meta) {
                                    // let term = row[0]; // TODO: row[1] if there is representative column!
                                    let term = row[1]; // TODO: row[0] if there is no representative column!
                                    let pvalue = data;
                                    // return "<div style='background-color:"+vm.scaleP(data)+"'>"+
                                    return "<div style='background-color:"+
                                           ((vm.remainingterms.includes(term)) ?
                                               vm.scaleP(pvalue) :
                                               "white")+"'>"+
                                           Number.parseFloat(data).toExponential(1)+
                                           "<div>";
                                  }
                                }
                            ]
                        } );
                }

            } else {
                // Multiple experiment setting

                // TODO: Check whether remaining terms is filled with all experiments!

                if (vm.remainingterms.length!=0 && vm.experimentenrichments.length!=0) {

                    // console.log('READY TO RENDER TABLE!');
                    // console.log(vm.remainingterms);
                    // console.log(vm.experimentenrichments);
                    // console.log(vm.representatives)

                    // Extract representatives for a current experiment and redraw the whole table
                    let representativesExperiment = vm.representatives.get(vm.experiment);
                    let tabledata = [...vm.experimentenrichments].map(
                        x => [representativesExperiment.has(x[0]) ? 
                                representativesExperiment.get(x[0]) : 
                                x[0], // representatives
                              x[0], // GO term id
                              vm.names[x[0]]].concat(vm.experiments.map(y=>x[1][y])) // GO term description
                    );

                    // Existing datatable cannot be reinitialized, so if it already exist destroy it before
                    if ($.fn.dataTable.isDataTable( this.$refs.datatable )) {
                        console.log("data-table/watch/experimentenrichments: Destroying existing table!");
                        vm.table.destroy(); // Not enough, we need to empty table as well!
                        $(vm.$refs.datatable).empty();
                    }

                    // Construct an array of p-value scales - one for each column/experiments!
                    // The p-value color scale should be defined only for the non-redundant terms.

                    // Applies a function to all p-values from each experiment, and outputs the result as 
                    // a Map along with the experiment name
                    let maxMin = function(f) {
                        return new Map(vm.experiments.map(e => [e,
                        f(...[...vm.experimentenrichments] // f will later be either Math.min or Math.max
                            .filter(x=>vm.remainingterms.includes(x[0])) // Only non-redundant terms!
                            .map(y => y[1].hasOwnProperty(e) ? 
                                                   y[1][e] : 
                                                   null) // These will be filtered out later
                            .filter(x=>x) // Filter out null values from previous map
                        )] 
                    ))};

                    // Apply Math.min and Math.max functions to find min and max p-values for each experiment
                    let minValues = maxMin(Math.min);
                    let maxValues = maxMin(Math.max);

                    // An array of scales for color which depends on the pvalues of newly loaded GO terms
                    let experimentScales = new Map(vm.experiments.map( e => [ e, 
                        vm.scaleP.domain([minValues.get(e),maxValues.get(e)]) ] ));

                    vm.table = $(this.$refs.datatable).DataTable( {
                            scrollY:        "400px",
                            scrollCollapse: true,
                            deferRender:    true, // Faster on large tables!
                            // scroller:       true,  // Requires pagination to be disabled (paging: false)! 
                            paging: false,
                            // autoWidth:      false,
                            data: tabledata,
                            columns: [ 

                                { title: "Representative" , 
                                  width: "70px", 
                                  defaultContent: "",
                                  render: function(data,type,row,meta) {
                                    return "<div style='color:"+
                                            // (vm.representatives.has(row[0]) ?  "black" : "gray")+"'>"+
                                            ( (row[0] != row[1]) ? "black" : "LightGray" )+"'>"+
                                            row[0] + 
                                            ( (row[0] == row[1]) ? "*" : "" ) + 
                                            "<div>";
                                  }
                                }, 

                               { title: "Term ID" , defaultContent: "" },
                               // { title: "description", defaultContent: "", width: "100px"} ].concat( 

                               { title: "description", defaultContent: ""} ].concat( 
                                vm.experiments.map( function(d,i) {
                                    return {title: String(i), 
                                        defaultContent: "", // Not needed because we are handling it in render!
                                        width: "50px",
                                        className: 'dt-body-center',
                                        // TODO: Past versions of handling p-value color:
                                        //       experimentScales.get(d)(data)+"'>"+
                                        //       ((data<0.01)?experimentScales.get(d)(data):"white")+"'>"+
                                        //       ((data<0.01)?vm.scaleP(data):"white")+"'>"+
                                        render: function(data,type,row,meta) {
                                            let term = row[0];
                                            let experiment = d;
                                            let pvalue = data;
                                            return "<div style='background-color:"+
                                                   ((vm.remainingterms.includes(term)) ?
                                                       experimentScales.get(experiment)(pvalue) :
                                                       "white")+"'>"+
                                                   ((typeof data !== 'undefined') ?
                                                        Number.parseFloat(pvalue).toExponential(1) :
                                                        "")+
                                                   "<div>";

                                        }
                                    };
                                })
                            )
                        } );

                    // Putting tooltips on headers
                    vm.table
                      .columns()
                      .header()
                      .to$()
                      .attr('data-toggle','tooltip')
                      .attr('title',function(){
                          // Assumption is that header titles are simple integers!
                          if (Number.isInteger(Number($(this).text()))) {
                              return vm.experiments[Number($(this).text())];
                          }
                      });

                }

            }

        },

    },
    // With the enclosing div in the template we have to access the table through a reference!
    // With enclosing div: this.$refs.datatable
    // Without enclosing div: this.$el
    template: '<div style="width:600px"><table ref="datatable" class="display compact"></table></div>'
});

// Custom zip function
const zip = (arr1, arr2) => arr1.map((k, i) => [k, arr2[i]]); // custom zip function

