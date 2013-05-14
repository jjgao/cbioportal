/**
 * Constructor for the network (sbgn) visualization class.
 *
 * @param divId     target div id for this visualization.
 * @constructor
 */
function NetworkSbgnVis(divId)
{
	// call the parent constructor
	NetworkVis.call(this, divId);
	this.networkTabsSelector = "#" + this.divId + " #network_tabs_sbgn";
    this.filteringTabSelector = "#" + this.divId + " #filtering_tab_sbgn";
    this.genesTabSelector = "#" + this.divId + " #genes_tab_sbgn";
    this.detailsTabSelector = "#" + this.divId + " #element_details_tab_sbgn";

	// node glyph class constants
	this.MACROMOLECULE = "macromolecule";
	this.PROCESS = "process";
	this.COMPARTMENT = "compartment";
	this.COMPLEX = "complex";
	this.HUGOGENES = new Array();
}

//this simulates NetworkSbgnVis extends NetworkVis (inheritance)
NetworkSbgnVis.prototype = new NetworkVis("");

//update constructor
NetworkSbgnVis.prototype.constructor = NetworkSbgnVis;

//TODO override necessary methods (filters, inspectors, initializers, etc.) to have a proper UI.

//Genomic data parser method
NetworkSbgnVis.prototype.parseGenomicData = function(genomicData)
{
	var hugoToGene 		= "hugo_to_gene_index";
	var geneData   		= "gene_data";
	var cna 	   	= "cna";
	var hugo 	   	= "hugo";
	var mrna	   	= "mrna";
	var mutations  		= "mutations";
	var rppa	   	= "rppa";
	var percent_altered 	= "percent_altered";

	//first extend node fields to support genomic data
	this.addGenomicFields();

	// iterate for every hugo gene symbol in incoming data
	for(var hugoSymbol in genomicData[hugoToGene])
	{
		var geneDataIndex 	= genomicData[hugoToGene][hugoSymbol];		// gene data index for hugo gene symbol
		var _geneData 	= genomicData[geneData][geneDataIndex];		// corresponding gene data

		// Arrays and percent altered data 
		var cnaArray   		= _geneData[cna];
		var mrnaArray  	= _geneData[mrna];
		var mutationsArray 	= _geneData[mutations];
		var rppaArray	  	= _geneData[rppa];
		var percentAltered 	= _geneData[percent_altered];
		
		// corresponding cytoscape web node
		var vis = this._vis;
		var targetNodes = findNode(hugoSymbol, vis );

		this.calcCNAPercents(cnaArray, targetNodes);
		this.calcMutationPercent(mutationsArray, targetNodes);
		this.calcRPPAorMRNAPercent(mrnaArray, mrna, targetNodes);
		this.calcRPPAorMRNAPercent(rppaArray, rppa, targetNodes);
		var alterationPercent = parseInt(percentAltered.split('%'),10)/100;		
		var alterationData =  {PERCENT_ALTERED: alterationPercent };
		this._vis.updateData("nodes",targetNodes, alterationData);
	}
};


//Searches an sbgn node whose label fits with parameter hugoSymbol
function findNode(hugoSymbol, vis)
{
	var nodeArray = vis.nodes();
	var nodes = new Array();
	for ( var i = 0; i < nodeArray.length; i++) 
	{
		if( nodeArray[i].data.glyph_label_text == hugoSymbol)
		{
			nodes.push(nodeArray[i].data.id);
		}
	}
	return nodes;
}


//calculates cna percents ands adds them to target node
NetworkSbgnVis.prototype.calcCNAPercents = function(cnaArray, targetNodes)
{  
	var amplified	= "AMPLIFIED";
	var gained    	= "GAINED";
	var hemiDeleted = "HEMIZYGOUSLYDELETED";
	var homoDeleted	= "HOMODELETED";

	var percents = {};
	percents[amplified] = 0;
	percents[gained] = 0;
	percents[hemiDeleted] = 0;
	percents[homoDeleted] = 0;

	var increment = 1/cnaArray.length;

	for(var i = 0; i < cnaArray.length; i++)
	{
		if(cnaArray[i] != null)
			percents[cnaArray[i]] += increment; 

	}

	var ampl = { PERCENT_CNA_AMPLIFIED:percents[amplified] };
	var gain = { PERCENT_CNA_GAINED: percents[gained]};
	var hem =  { PERCENT_CNA_HEMIZYGOUSLY_DELETED: percents[hemiDeleted] };
	var hom =  { PERCENT_CNA_HOMOZYGOUSLY_DELETED: percents[homoDeleted]};
	
	this._vis.updateData("nodes",targetNodes, ampl);
	this._vis.updateData("nodes",targetNodes, gain);
	this._vis.updateData("nodes",targetNodes, hem);
	this._vis.updateData("nodes",targetNodes, hom);

};

//calculates rppa or mrna percents ands adds them to target node, data indicator determines which data will be set
NetworkSbgnVis.prototype.calcRPPAorMRNAPercent = function(dataArray, dataIndicator, targetNodes)
{  
	var up		= "UPREGULATED";
	var down   	= "DOWNREGULATED";
	
	var upData = null;
	var DownData = null;

	var percents = {};
	percents[up] = 0;
	percents[down] = 0;

	var increment = 1/dataArray.length;

	for(var i = 0; i < dataArray.length; i++)
	{
		if(dataArray[i] != null)
			percents[dataArray[i]] += increment; 
	}

	if (dataIndicator == "mrna") 
	{
		upData =  {PERCENT_MRNA_UP: percents[up]};
		downData = {PERCENT_MRNA_DOWN: percents[down]};
	} 
	else if(dataIndicator == "rppa") 
	{
		upData =   {PERCENT_RPPA_UP: percents[up]};
		downData = {PERCENT_RPPA_DOWN: percents[down]};
	}
	
	this._vis.updateData("nodes",targetNodes, upData);
	this._vis.updateData("nodes",targetNodes, downData);
};

//calculates mutation percents ands adds them to target node
NetworkSbgnVis.prototype.calcMutationPercent = function(mutationArray, targetNodes)
{  
	var percent = 0;
	var increment = 1/mutationArray.length
	for(var i = 0; i < mutationArray.length; i++)
	{
		if(mutationArray[i] != null)
			percent += increment;  
	}
	var mutData = {PERCENT_MUTATED: percent};
	this._vis.updateData("nodes",targetNodes, mutData);
};

//extends node fields by adding new fields according to genomic data
NetworkSbgnVis.prototype.addGenomicFields = function()
{
	var cna_amplified 	= {name:"PERCENT_CNA_AMPLIFIED", type:"number", defValue: 0};
	var cna_gained		= {name:"PERCENT_CNA_GAINED", type:"number"};
	var cna_homodel 	= {name:"PERCENT_CNA_HOMOZYGOUSLY_DELETED", type:"number", defValue: 0};
	var cna_hemydel		= {name:"PERCENT_CNA_HEMIZYGOUSLY_DELETED", type:"number", defValue: 0};

	var mrna_up 		= {name:"PERCENT_MRNA_UP", type:"number", defValue: 0};
	var mrna_down 		= {name:"PERCENT_MRNA_DOWN", type:"number", defValue: 0};

	var rppa_up 		= {name:"PERCENT_RPPA_UP", type:"number", defValue: 0};
	var rppa_down 		= {name:"PERCENT_RPPA_DOWN", type:"number", defValue: 0};

	var mutated			= {name:"PERCENT_MUTATED", type:"number", defValue: 0};
	var altered			= {name:"PERCENT_ALTERED", type:"number", defValue: 0};

	var label			= {name:"label", type:"text", defValue: ""};


	this._vis.addDataField(cna_amplified);
	this._vis.addDataField(cna_gained);
	this._vis.addDataField(cna_homodel);
	this._vis.addDataField(cna_hemydel);

	this._vis.addDataField(mrna_down);
	this._vis.addDataField(mrna_up);

	this._vis.addDataField(rppa_down);
	this._vis.addDataField(rppa_up);

	this._vis.addDataField(mutated);
	this._vis.addDataField(altered);
	//this._vis.addDataField(label);
};

/**
 * Initializes all necessary components. This function should be invoked, before
 * calling any other function in this script.
 *
 * @param vis	CytoscapeWeb.Visualization instance associated with this UI
 */
NetworkSbgnVis.prototype.initNetworkUI = function(vis, genomicData)
{
    this._vis = vis;
    this._linkMap = this._xrefArray();

    // init filter arrays
    this._alreadyFiltered = new Array();
    this._filteredBySlider = new Array();
    this._filteredByDropDown = new Array();
    this._filteredByIsolation = new Array();
    this._edgeTypeVisibility = this._edgeTypeArray();
    this._edgeSourceVisibility = this._edgeSourceArray();
    // parse and add genomic data to cytoscape nodes
    this.parseGenomicData(genomicData); 
    // this.setInitialData();

    this._geneWeightMap = this._geneWeightArray(0, null);
    this._geneWeightThreshold = this.ALTERATION_PERCENT;
    this._maxAlterationPercent = this._maxAlterValNonSeed(this._geneWeightMap);

    this._resetFlags();

    this._initControlFunctions();
    this._initLayoutOptions();

    this._initMainMenu();

    this._initDialogs();
    this._initPropsUI();
    this._initSliders();
    this._initDropDown();
    this._initTooltipStyle();

    // add listener for the main tabs to hide dialogs when user selects
    // a tab other than the Network tab

    var self = this;

    var hideDialogs = function(evt, ui){
        self.hideDialogs(evt, ui);
    };

    $("#tabs").bind("tabsshow", hideDialogs);

    // this is required to prevent hideDialogs function to be invoked
    // when clicked on a network tab
    $(this.networkTabsSelector).bind("tabsshow", false);

    // init tabs
    $(this.networkTabsSelector).tabs();
    $(this.networkTabsSelector + " .network-tab-ref").tipTip(
        {defaultPosition: "top", delay:"100", edgeOffset: 10, maxWidth: 200});

    this._initGenesTab();
    this._refreshGenesTab();
    this._refreshRelationsTab();

    // adjust things for IE
    this._adjustIE();

    // make UI visible
    this._setVisibility(true);

};


/*
/**
 * Calculates weight values for each gene by using the formula:
 *
 * weight = Max[(Total Alteration of a node),
 *    Max(Total Alteration of its neighbors) * coeff] * 100
 *
 * @param coeff	coefficient value used in the weight function
 * @returns		a map (array) containing weight values for each gene
 */


NetworkSbgnVis.prototype._geneWeightArray = function(flag, selected)
{
	var parents = new Array();
	var pId = new Array();
	var weights = new Array();
	var processes = new Array();
	var leaves = new Array();
	
	var nodes = this._vis.nodes();
	// just the visible nodes
	var i = 0;
	while (i < nodes.length)
	{
		if (!this.currentVisibility(nodes[i]))
		{
			weights[nodes[i].data.id] = -1;
			nodes.splice(i, 1);
		}
		else
		{
			i++;
		}
	}

	if (flag == 0) // first time calling this method
	{
		// A0: initialization
		for (var i = 0; i < nodes.length; i++)
		{
			var glyph = nodes[i].data.glyph_class;
			if (glyph == this.MACROMOLECULE)
			{
				// first update hugogenes to hold one node of every macromolecule hugotype
				// these are either proteins or genes
				var label = this.geneLabel(nodes[i].data);
				var check = 0;

				for (var j = 0; j < this.HUGOGENES.length; j++)
				{
					if (label == this.geneLabel(this.HUGOGENES[j].data))
					{
						check = 1;
						break;
					}
				}
				// if its a new hugolabel add it to the hugogenes
				if (check == 0)
				{
						this.HUGOGENES.push(nodes[i]);
				}

			}
			else
			{ 
				// set initial weight of other nodes as 0
				weights[nodes[i].data.id] = 0;
				
				// make a list of processes for latter update of weights
				if (glyph == this.PROCESS)
				{
					processes.push(nodes[i]);
				}
				
			}
			// initialize the parent ID
			pId[nodes[i].data.id] = -1;
			
			// update leaves array
			if ( this._vis.childNodes(nodes[i].data.id).length == 0)
			{				
				leaves.push(nodes[i]);
			}
			
			if (nodes[i].data["PERCENT_ALTERED"] != 0)
			{
				weights[nodes[i].data.id] = nodes[i].data["PERCENT_ALTERED"]  *  100;
				
			}
			else
			{
				weights[nodes[i].data.id] = 0;
			}
		}
	}
	else // flag == 1 for show and hide methods
	{
		for (var i = 0; i < nodes.length; i++)
		{
			if (nodes[i].data.glyph_class == this.PROCESS)
			{
				processes.push(nodes[i]);
			}
			// update leaves array
			if ( this._vis.childNodes(nodes[i].data.id).length == 0)
			{				
				leaves.push(nodes[i]);
			}
			
			weights[nodes[i].data.id] = 0;
		}
		for (var i = 0; i < selected.length; i++)
		{
			weights[selected[i].data.id] = 1;
		//	if (selected[i].data.glyph_class == this.MACROMOLECULE)
		//	{
		//		weights[selected[i].data.id] =  selected[i].data["PERCENT_ALTERED"]  *  100;
		//	}
		}
		
	}
	
	// update parents array
	var k = 0;
	for (var i = 0; i < nodes.length; i++)
	{
		if ( this._vis.childNodes(nodes[i]).length > 0)
		{
			var children = this._vis.childNodes(nodes[i]);
			for (var j = 0; j < children.length; j++)
			{
				pId[children[j].data.id] = k;
				
			}
			parents[k] = nodes[i].data.id;
			k++;
		}
	}
	
	 // A1: update process weights based on neighbors
	// for each process, set the initial weight the maximum of its neighbors
	for (var i = 0; i < processes.length; i++)
	{
		var max = 0;
		var neighbors = this._vis.firstNeighbors([processes[i]]).neighbors;
		for(var j = 0; j < neighbors.length; j++)
		{
			var nID = neighbors[j].data.id;
			if (weights[nID] > max)
			{
				max = weights[nID];
			}
		}
		if (weights[processes[i].data.id] < max)
		{
			weights[processes[i].data.id] = max;
		}
	}
	
	// update all neighbors of processes to have the weight of the process
	for (var i = 0; i < processes.length; i++)
	{
		var w = weights[processes[i].data.id] ;
		var neighbors = this._vis.firstNeighbors([processes[i]]).neighbors;
		var complexNeighbors = new Array();
		for(var j = 0; j < neighbors.length; j++)
		{
			if (weights[neighbors[j].data.id]  < w)
			{
				weights[neighbors[j].data.id] = w;
			}
		}
	}
	
	// make the parent nodes hold the maximum weight of its children
	for (var i = 0; i < leaves.length; i++)
	{
		var node = leaves[i];
		
		while (true)
		{
			// see if we can go higher one level
			if (pId[node.data.id] == -1)
				break;
			var parentID = parents[pId[node.data.id]];
			// if the weight of the parent is less than the child, update its weight by the child
			if (weights[node.data.id] > weights[parentID] )
			{
				 weights[parentID]  = weights[node.data.id];
			}
			else
			{
				// if the parent is higher than the child what should be done?
				// by default it breaks
				// in our case else will never happen
				// kept for debugging or future changes
				break;
			}
			// go up one level
			node = this._vis.node(parentID);
		}
	}
	

	// A3: propogate max values to parents from leaves to root
	for (var i = 0; i < leaves.length; i++)
	{
		var node = leaves[i];
		var nodeID = node.data.id;
		var pCheck= pId[nodeID];
		while (pCheck > -1)
		{
			var parent = this._vis.node(parents[pCheck]);
			var parentID = parent.data.id;
			if (weights[parentID] < weights[nodeID])
			{
				weights[parentID]  = weights[nodeID];
			}
			pCheck = pId[parentID];
			node = parent;
		}
	}
	
	
	// make sure all complex nodes 
	// A4: propogate max values of complex hierarchies down to leaves
	for (var i = 0; i < nodes.length; i++)
	{
		var n = nodes[i];
		if (n.data.glyph_class == this.COMPLEX && weights[n.data.id] > 0)
		{
			var children = this._vis.childNodes(n);
			while (children.length > 0)
			{
				var nextGeneration = new Array();
				for(var j = 0; j < children.length; j++)
				{
					weights[children[j].data.id] = weights[n.data.id];
					if (children[j].data.glyph_class == this.COMPLEX)
					{
						nextGeneration.push(children[j]);
					}
				}
				children = nextGeneration;
			}	
		}
	}
	
	return weights;
};


/**
 * Creates an array of visible (i.e. non-filtered) genes.
 *
 * @return		array of visible genes
 */
NetworkSbgnVis.prototype._visibleGenes = function()
{
	// set the genes to be shown in the gene list
	var genes = this.HUGOGENES;

	// sort genes by glyph class (alphabetically)
	genes.sort(_geneSort);
	return genes;
};

/**
 * Updates selected genes when clicked on a gene on the Genes Tab. This function
 * helps the synchronization between the genes tab and the visualization.
 *
 * @param evt	target event that triggered the action
 */
NetworkSbgnVis.prototype.updateSelectedGenes = function(evt)
{
	// this flag is set to prevent updateGenesTab function to update the tab
	// when _vis.select function is called.
	
	this._selectFromTab = true;
	var hugoIds = new Array();

	// deselect all nodes
	this._vis.deselect("nodes");


	// collect id's of selected node's on the tab
	$(this.geneListAreaSelector + " select option").each(
        function(index)
        {
            if ($(this).is(":selected"))
            {
                var nodeId = $(this).val();
                hugoIds.push(nodeId);
            }
        });

	var nodes = this._vis.nodes();
	
	// array of nodes to select
    	var nodeIds = new Array();
	var check = 0;
	for (var i = 0; i < nodes.length; i++)
	{
		for ( var j = 0; j < hugoIds.length; j++) 
		{
			if (this.geneLabel(nodes[i].data) == this.geneLabel(this._vis.node(hugoIds[j]).data))
			{
				nodeIds.push(nodes[i].data.id);
				break;
			}
		}
	}

	// select all checked nodes
	this._vis.select("nodes", nodeIds);

	// reset flag
	this._selectFromTab = false;
};
// used for selecting genes
/**
 * returns all nodes from HUGOGENES array that have the same label as
 * the nodes in the elements list.
 */
NetworkSbgnVis.prototype.hugoGenes = function(elements)
{
	var hugoElements = new Array();
	
	for (var i=0; i < elements.length; i++)
	{
		if (elements[i].data.glyph_class != this.MACROMOLECULE)
		{
			continue;
		}
		for(var j=0; j<this.HUGOGENES.length; j++)
		{
			if (this.geneLabel(this.HUGOGENES[j].data) == this.geneLabel(elements[i].data))
			{
				hugoElements.push(this.HUGOGENES[j]);
				break;
			}
		}
	}
	return hugoElements;
}
/**
 * returns all nodes in the graph that have the same label as
 * the nodes in the elements list.
 */
NetworkSbgnVis.prototype.sameHugoGenes = function(elements)
{
	var sameElements = new Array();
	var nodes = this._vis.nodes();
	for (var i=0; i < elements.length; i++)
	{
		if (elements[i].data.glyph_class == this.MACROMOLECULE)
		{
			for(var j=0; j<nodes.length; j++)
			{
				if (this.geneLabel(nodes[j].data) == this.geneLabel(elements[i].data))
				{
					sameElements.push(nodes[j]);
				}
			}
		}
		else
		{
			sameElements.push(elements[i]);
		}
		
	}
	return sameElements;
}
/**
 * Updates the gene tab if at least one node is selected or deselected on the
 * network. This function helps the synchronization between the genes tab and
 * visualization.
 * for now whenever a gene is selected the row associated with the glyph label
 * is highlighted and all genes with same property are selected.
 *
 * @param evt	event that triggered the action
 */
NetworkSbgnVis.prototype.updateGenesTab = function(evt)
{
    var selected = this._vis.selected("nodes");

    // do not perform any action on the gene list,
    // if the selection is due to the genes tab
    if(!this._selectFromTab)
    {
	
        if (_isIE())
        {
            this._setComponentVis($(this.geneListAreaSelector + " select"), false);
        }

        // deselect all options
        $(this.geneListAreaSelector + " select option").each(
            function(index)
            {
                $(this).removeAttr("selected");
            });
	var nodes = this._vis.nodes();
        
	// select all nodes with same label
	var sameNodes =  this.sameHugoGenes(selected);
	this._vis.select("nodes", sameNodes);
	// select all nodes with same glyph label text in the gene tab list
	var hugos = this.hugoGenes(selected);
        for (var i=0; i < hugos.length; i++)
        {
		$(this.geneListAreaSelector + " #" +  _safeProperty(hugos[i].data.id)).attr(
				     "selected", "selected");
        }

        if (_isIE())
        {
            this._setComponentVis($(this.geneListAreaSelector + " select"), true);
        }
    }
    // also update Re-submit button
    if (selected.length > 0)
    {
        // enable the button
        $(this.genesTabSelector + " #re-submit_query").button("enable");
    }
    else
    {
        // disable the button
        $(this.genesTabSelector + " #re-submit_query").button("disable");
    }
};

/**
 * Comparison function to sort genes alphabetically.
 *
 * @param node1	node to compare to node2
 * @param node2 node to compare to node1
 * @return 		positive integer if node1 is alphabetically greater than node2
 * 				negative integer if node2 is alphabetically greater than node1
 * 				zero if node1 and node2 are alphabetically equal
 */
function _geneSort (node1, node2)
{
    if (node1.data.glyph_label_text > node2.data.glyph_label_text)
    {
        return 1;
    }
    else if (node1.data.glyph_label_text < node2.data.glyph_label_text)
    {
        return -1;
    }
    else
    {
        return 0;
    }
}

// make a function to give the label so that we can neasily over ride in the sbgn
NetworkSbgnVis.prototype.geneLabel = function(data)
{
	return data.glyph_label_text;
};

NetworkSbgnVis.prototype.geneDetailsCheck = function (selected)
{
	// if all the selected nodes have the same glyph label text show only one
	// to be decided upon by Ugur Hoca
	var check = 1;
	for(var i=0; i < selected.length; i++)
	{
		if(selected[i].data.glyph_class == this.COMPARTMENT || selected[i].data.glyph_class == this.COMPLEX)
		{
			check = -1;
			break;
		}
		if(this.geneLabel(selected[i].data) != this.geneLabel(selected[0].data))
		{
			check = 0;
			break;
		}
	}
	return check;
};


/**
 * Filters out all non-selected nodes.
 */
NetworkSbgnVis.prototype.filterNonSelected = function()
{
	var self = this;

	// this is required to pass "this" instance to the listener
	var geneVisibility = function(element) {
	return self.geneVisibility(element);
	};

	// update selected elements map
	var selected = this._vis.selected("nodes");
	var map = new Array();
	var weights = this._geneWeightArray(1, selected);
	var nodes = this._vis.nodes();

	for (var i=0; i < nodes.length; i++)
	{
		if (weights[nodes[i].data.id] == 1)
		{
			var key = nodes[i].data.id;
			map[key] = nodes[i];
		}
	}
	this._selectedElements = map;

	// filter out non-selected elements
	this._vis.filter('nodes', geneVisibility);

	// also, filter disconnected nodes if necessary
	this._filterDisconnected();

	// refresh Genes tab
	this._refreshGenesTab();
	this.updateGenesTab();

	// visualization changed, perform layout if necessary
	this._visChanged();
};
/**
 * Filters out all selected genes.
 */
NetworkSbgnVis.prototype.filterSelectedGenes = function()
{
	var self = this;

	// this is required to pass "this" instance to the listener
	var selectionVisibility = function(element) {
	return self.selectionVisibility(element);
	};

	// update selected elements map
	var selected = this._vis.selected("nodes");
	var map = new Array();
	var weightsSelected = this._geneWeightArray(1, selected);
	var nodes = this._vis.nodes();
	var complement = new Array();
	for (var i=0; i < nodes.length; i++)
	{
		if (weightsSelected[nodes[i].data.id] == 0)
		{
			complement.push(nodes[i]);
		}
	}
	weightsSelected = this._geneWeightArray(1, complement);
	for (var i=0; i < nodes.length; i++)
	{
		if (weightsSelected[nodes[i].data.id] == 0)
		{
			var key = nodes[i].data.id;
			map[key] = nodes[i];
		}
	}
	this._selectedElements = map;

	// filter out selected elements
	this._vis.filter("nodes", selectionVisibility);

	// also, filter disconnected nodes if necessary
	this._filterDisconnected();

	// refresh genes tab
	this._refreshGenesTab();

	// visualization changed, perform layout if necessary
	this._visChanged();
};

/**
 * Initializes the gene filter sliders.
 */
NetworkSbgnVis.prototype._initSliders = function()
{
    var self = this;

    var keyPressListener = function(evt) {
      self._keyPressListener(evt);
    };

    var weightSliderStop = function(evt, ui) {
        self._weightSliderStop(evt, ui);
    };

    var weightSliderMove = function(evt, ui) {
        self._weightSliderMove(evt, ui);
    };


    // add key listeners for input fields

    $(this.filteringTabSelector + " #weight_slider_field").keypress(keyPressListener);
    $(this.filteringTabSelector + " #affinity_slider_field").keypress(keyPressListener);

    // show gene filtering slider
    $(this.filteringTabSelector + " #weight_slider_bar").slider(
        {value: this.ALTERATION_PERCENT,
            stop: weightSliderStop,
            slide: weightSliderMove});

};

/**
 * Updates the contents of the details tab according to
 * the currently selected elements.
 *
 * @param evt
 */
NetworkSbgnVis.prototype.updateDetailsTab = function(evt)
{
    // TODO also consider selected edges?
	var selected = this._vis.selected("nodes");
    var data;

    var self = this;

    // clean previous content
    $(self.detailsTabSelector + " div").empty();
	$(self.detailsTabSelector + " .error").hide();

    if (selected.length == 1)
    {
        data = selected[0].data;
    }
    else if (selected.length > 1)
    {
        //$(self.detailsTabSelector + " div").empty();
        $(self.detailsTabSelector + " .error").append(
            "Currently more than one node is selected. Please, select only one node to see details.");
	    $(self.detailsTabSelector + " .error").show();
        return;
    }
    else
    {
        //$(self.detailsTabSelector + " div").empty();
        $(self.detailsTabSelector + " .error").append(
            "Currently there is no selected node. Please, select a node to see details.");
	    $(self.detailsTabSelector + " .error").show();
        return;
    }

    var handler = function(queryResult) {
        // check the initial conditions, if they don't match do nothing
	    // (they may not match because of a delay in the ajax request)

	    selected = self._vis.selected("nodes");

	    if (selected.length != 1 ||
	        data.id != selected[0].data.id)
	    {
		    return;
	    }

        // update tab content
        $(self.detailsTabSelector + " .biogene-content").empty();

        if (queryResult.returnCode != "SUCCESS")
        {
            $(self.detailsTabSelector + " .error").append(
	            "Error retrieving data: " + queryResult.returnCode);
	        $(self.detailsTabSelector + " .error").show();
        }
        else
        {
	        if (queryResult.count > 0)
	        {
		        // generate the view by using backbone
		        var biogeneView = new BioGeneView(
			        {el: self.detailsTabSelector + " .biogene-content",
					data: queryResult.geneInfo[0]});
	        }
	        else
	        {
		        $(self.detailsTabSelector + " .error").append(
			        "No additional information available for the selected node.");
	        }
        }

	    // generate view for genomic profile data
	    var genomicProfileView = new GenomicProfileView(
		    {el: self.detailsTabSelector + " .genomic-profile-content",
			data: data});
    };


	if (data.type == this.DRUG)
	{
		// update tab content
		$(self.detailsTabSelector + " div").empty();

		var drugView = new DrugInfoView({el: this.detailsTabSelector + " .drug-info-content",
			data: data,
			linkMap: this._linkMap,
			idPlaceHolder: this.ID_PLACE_HOLDER,
			edges: this._vis.edges()});
	}
	// send biogene query for only genes
	else
	{
		var queryParams = {"query": data.label,
			"org": "human",
			"format": "json"};

		$(self.detailsTabSelector + " .biogene-content").append(
			'<img src="images/ajax-loader.gif">');

		$.post("bioGeneQuery.do",
			queryParams,
			handler);
	}
};
