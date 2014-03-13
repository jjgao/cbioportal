/*
 * Basic DC PieChart Component.
 * 
 * @param _param -- Input object
 *                  chartDivClass: currently only accept class name for DIV.chartDiv,
 *                                  (TODO: Add more specific parameters later) 
 *                  chartID: the current pie chart ID which is treated as
 *                           identifier using in global,
 *                  attrId: the attribute name, 
 *                  displayName: the display content of this attribute, 
 *                  transitionDuration: this will be used for initializing
 *                                      DC Pie Chart,
 *                  ndx: crossfilter dimension, used by initializing DC Pie Chart
 *                  chartColors: color schema
 *                  
 * @interface: getChart -- return DC Pie Chart Object.
 * @interface: getCluster -- return the cluster of DC Pie Chart.
 * @interface: pieLabelClickCallbackFunction -- pass a function to be called when
 *                                              the pie label been clicked.
 * @interface: scatterPlotCallbackFunction -- pass a function to connect with
 *                                            Scatter Plot after filtering DC
 *                                            Pie Chart.
 * @interface: postFilterCallbackFunc -- pass a function to be called after DC Pie
 *                                       Chart filtered.
 *                                       
 * @authur: Hongxin Zhang
 * @date: Mar. 2014
 * 
 */


var PieChart = function(){
    var pieChart, cluster;
    
    //All DIV ID names are organized based on the structure rule, initialized 
    //in initParam function
    var DIV = {
        parentID : "",
        mainDiv : "",
        chartDiv : "",
        labelTableID : "",
        labelTableTdID : ""
    };
    
    var chartID, 
        className, 
        selectedAttr, 
        selectedAttrDisplay, 
        transitionDuration, 
        ndx, 
        chartColors;
    
    var label =[],
        labelSize = 10,
        fontSize = labelSize +1;
            
    var postFilterCallback,
        pieLabelClickCallback,
        scatterPlotCallback;
    
    
    
    //This function is designed to draw Pie Labels based on current color the
    //Pie Chart has. Pagging function will be added when the number of labels
    //bigger than 5.
    function addPieLabels() {
       
        $('#' + DIV.mainDiv + ' .study-view-pie-label').html("");
          
        initLabelInfo();
        
        if(label.length > 6){
            bigLabelFunction();
        }else{
            smallLabelFunction();
        }
        
        addPieLabelEvents();
    }
    
    //Called when the number of label biggen than 6, used by addPieLabels()
    function bigLabelFunction() {
        var _totalTableNum = 1;
        
        if( label.length % 5 === 0 )
            _totalTableNum = label.length / 5;
        else
            _totalTableNum = parseInt( label.length / 5 ) + 1;
        
        
        for(var j = 0 ; j < label.length ; j+=5){
            var _innerID = 0,
                _tableId = parseInt(j/5),
                _showTableStyle = '',
                _leftArrowColor = 'blue',
                _rightArrowColor = 'blue',
                _currentTableDivId = DIV.labelTableID + "-" + _tableId;

            if( _tableId !== 0 )  
                _showTableStyle = 'style="display:none"';

            $('#' + DIV.mainDiv)
                    .find('.study-view-pie-label')
                    .append("<table id='"+_currentTableDivId+"' "+
                        _showTableStyle + " ></table>");

            for(var i=j; i< j+5; i++){
                if(i<label.length){
                    var _tmpName = label[i].name;
                    if(_tmpName.length > 9)
                       _tmpName = _tmpName.substring(0,5) + " ...";
                    if(( i - _tableId ) % 2 === 0 ){
                        $('#' + DIV.mainDiv)
                                .find('#' + _currentTableDivId)
                                .append("<tr id="+ _innerID +" width='150px'></tr>");
                        _innerID++;
                    } 
                    $('#' + DIV.mainDiv)
                            .find('#'+_currentTableDivId+
                                ' tr:nth-child(' + _innerID +')')
                            .append('<td class="pieLabel" id="' +
                                DIV.labelTableTdID+label[i].id + "-" + i +
                                '" style="font-size:' + fontSize +
                                'px"><svg width="75" height="13"><rect width="' +
                                labelSize+'" height="'+ labelSize +
                                '" style="fill:' + label[i].color + 
                                ';" /><text x="15" y="10">' + _tmpName +
                                '</text></svg></td>');

                    //Only add qtip when the length of pie label bigger than 9
                    if(label[i].name.length > 9){
                        var _qtip = jQuery.extend(true, {}, StudyViewBoilerplate.pieLabelQtip);
                        
                        _qtip.content.text = label[i].name;
                        $('#'+DIV.labelTableTdID+label[i].id+'-'+i).qtip(_qtip);
                    }
                }else{
                    if(( i - _tableId ) % 2 === 0){
                        $('#' + DIV.mainDiv)
                            .find('#' + _currentTableDivId)
                            .append("<tr id="+ _innerID +" width='150px'></tr>");
                        _innerID++;
                    } 
                    $('#' + DIV.mainDiv )
                        .find('#'+ _currentTableDivId +
                            ' tr:nth-child(' + _innerID +')')
                        .append('<td style="width="75px" height="15px"'+
                                ' font-size:'+fontSize+'px"></td>');
                }
            }

            if(_tableId === 0)
                _leftArrowColor = 'grey';
            if(_tableId+1 === _totalTableNum)
                _rightArrowColor = 'grey';

            $('#' + DIV.mainDiv)
                .find('#'+_currentTableDivId+
                    ' tr:nth-child(' + _innerID +')')
                .append('<td id="pieLabel-pagging-' + DIV.chartDiv + "-" +
                    _tableId + '" style=" width="75px" height="16px"' +
                    ' border="0px" font-size:' + fontSize+'px">' +
                    '<svg  width="75" height="13">' +
                    '<path class="pie-label-left-pagging" ' +
                    'd="M5 1 L0 11 L10 11 Z" fill="' + _leftArrowColor +
                    '"/><text x=15 y=10 fill="black">' +
                    ( _tableId + 1 ) + '/' + _totalTableNum + '</text>' +
                    '<path class="pie-label-right-pagging"' +
                    'd="M45 11 L40 1 L50 1 Z" fill="' + _rightArrowColor +
                    '"/></svg></td>');
        }
    }
    
    function addPieLabelEvents() {
        $('#' + DIV.mainDiv + ' .pieLabel').mouseenter(function(){
            var idArray = $(this).attr('id').split('-'),
                childID = Number(idArray[idArray.length-2])+1,
                fatherID = Number(idArray[idArray.length-3]);

            $('#' + DIV.chartDiv + ' svg>g>g:nth-child(' + childID+')').css({
                'fill-opacity': '.5',
                'stroke-width': '3'
            });
            
            drawMarker(childID,fatherID);
        });

        $('#' + DIV.mainDiv + ' .pieLabel').mouseleave(function(){
            var idArray = $(this).attr('id').split('-');
            var childID = Number(idArray[idArray.length-2])+1;
            var fatherID = Number(idArray[idArray.length-3]);
            var arcID = fatherID+"-"+(Number(childID)-1);

            $("#" + DIV.chartDiv + " svg g #arc-" + arcID).remove();

            $('#' + DIV.chartDiv + ' svg>g>g:nth-child(' + childID+')').css({
                'fill-opacity': '1',
                'stroke-width': '1px'
            });
        });

        $('#' + DIV.mainDiv + ' .pieLabel').click(function(){
            var idArray = $(this).attr('id').split('-');

            var childaLabelID = Number(idArray[idArray.length-1]),
                childID = Number(idArray[idArray.length-2])+1,
                chartID = Number(idArray[idArray.length-3]);

            var arcID = chartID+"-"+(Number(childID)-1);

            pieChart.onClick({
                key: label[childaLabelID].name, 
                value: label[childaLabelID].value
            });       

            pieChart.redraw;           

            $("#" + DIV.chartDiv + " svg g #" + arcID).remove();

            $('#' + DIV.chartDiv + ' svg>g>g:nth-child(' + childID+')').css({
                'fill-opacity': '1',
                'stroke-width': '1px'
            });
        });

        $('#' + DIV.mainDiv + ' .pie-label-left-pagging').click(function(){
            var tmpValue = $(this).parent().parent().attr('id').split('-');
            var currentTableID = Number(tmpValue[tmpValue.length-1]);
            if(currentTableID !== 0){
                var nextTableID = currentTableID-1;
                $('#'+ DIV.labelTableID + '-'+currentTableID).css('display','none');            
                $('#'+DIV.labelTableID+'-'+nextTableID).css('display','block');
            }
        });
        $('#' + DIV.mainDiv + ' .pie-label-right-pagging').click(function(){
            var fill = $(this).attr('fill');
            if(fill === 'blue'){
                var tmpValue = $(this).parent().parent().attr('id').split('-');
                var currentTableID = Number(tmpValue[tmpValue.length-1]);
                var nextTableID = currentTableID+1;

                $('#'+DIV.labelTableID+'-'+currentTableID).css('display','none');            
                $('#'+DIV.labelTableID+'-'+nextTableID).css('display','block');
            }
        });
    }
    
    //This function is designed to add functions like click, on, or other
    //other functions added after initializing this Pie Chart.
    function addFunctions() {
        if(selectedAttr !== 'CASE_ID'){
            pieChart.on("filtered", function(chart,filter){
                var _currentFilters = pieChart.filters();
                
                if(_currentFilters.length === 0){
                    $("#" + DIV.mainDiv + " .study-view-dc-chart-change")
                                .css('display','none');
                    $("#" + DIV.mainDiv)
                            .css({'border-width':'1px', 'border-style':'solid'});
                }else{
                    $("#" + DIV.mainDiv + " .study-view-dc-chart-change")
                                .css('display','block');
                    $("#" + DIV.mainDiv)
                            .css({'border-width':'2px', 'border-style':'inset'});
                }

                updateScatterPlot(_currentFilters);
                removeMarker();
                postFilterCallback();
            });
            pieChart.on("postRedraw",function(chart){
                addPieLabels();
            });
            pieChart.on("postRender",function(chart){
                addPieLabels();
            });
        }
    }
    
    //Initialize HTML tags which will be used for current Pie Chart.
    function createDiv() {
        if(selectedAttr === 'CASE_ID'){
            $(DIV.parentID)
                    .append("<div id=\"" + DIV.mainDiv + 
                    "\" class='study-view-dc-chart study-view-pie-main' "+
                    "style='display:none'><div id=\"" +
                    DIV.chartDiv + "\"></div></div>");
        }else
            $("#"+DIV.parentID).append("<div id=\"" + DIV.mainDiv +
                "\" class='study-view-dc-chart study-view-pie-main'>"+
                "<div id=\"" + DIV.chartDiv + "\" class='" + 
                className + "'  value='"+ selectedAttr + "," + 
                selectedAttrDisplay + ",pie'>"+
                "<div style='width:180px; float:right; text-align:center;'>"+
                "<span class='study-view-dc-chart-delete'>x</span>"+
                "<a href='javascript:varChart["+chartID+"].getChart().filterAll();dc.redrawAll();'>" +
                "<span title='Reset Chart' class='study-view-dc-chart-change' "+
                "style='font-size:10px;'>RESET</span></a></div>"+
                "<div style='width:180px;float:left;text-align:center'><chartTitleH4>" +
                selectedAttrDisplay + "</chartTitleH4></div></div>"+
                "<div class='study-view-pie-label'></div>"+
                "<div style='width:180px; text-align:center;float:left;'></div></div>");
        
    }
    
    //This function is designed to draw Pie Slice Marker(Arc) based on the
    //selected pie slice color.
    function drawMarker(_childID,_fatherID) {
        var _pointsInfo = 
                $('#' + DIV.chartDiv + ' svg>g>g:nth-child(' + _childID+')')
                    .find('path')
                    .attr('d')
                    .split(/[\s,MLHVCSQTAZ]/);
        
        var _pointsInfo1 = 
                $('#' + DIV.chartDiv + ' svg>g>g:nth-child(' + _childID+')')
                    .find('path')
                    .attr('d')
                    .split(/[A]/);
        
        var _fill = 
                $('#' + DIV.chartDiv + ' svg>g>g:nth-child(' + _childID+')')
                    .find('path')
                    .attr('fill');    

        var _x1 = Number(_pointsInfo[1]),
            _y1 = Number(_pointsInfo[2]),
            _largeArc = Number(_pointsInfo[6]),
            _x2 = Number(_pointsInfo[8]),
            _y2 = Number(_pointsInfo[9]),
            _r = Number(_pointsInfo[3]);

        if((_x1 - _x2!==0 || _y1 - _y2!==0) && _pointsInfo1.length === 2){
            var _xm = (_x1 + _x2) /2,
                _ym = (_y1 + _y2) /2;

            var m = Math.sqrt((Math.pow(_xm,2)+Math.pow(_ym,2)));

            var _tmpX = (_r + 3) / m * _xm,
                _tmpY = (_r + 3) / m * _ym;

            if(_xm === 0 && _ym === 0){
                _tmpY = 0;
                _tmpX = _r + 6;
            }

            if(_largeArc === 1 && Math.abs(_x1 - _x2) >0.1) {
                _tmpX = -_tmpX;
                _tmpY = -_tmpY;
            }

            var _pointOne = Math.atan2(_y1,_x1);
            var _pointTwo = Math.atan2(_y2,_x2);

            if(_pointOne < -Math.PI/2){
                _pointOne = Math.PI/2 + Math.PI *2 +_pointOne;
            }else{
                _pointOne = Math.PI/2 +_pointOne;
            }

            if(_pointTwo < -Math.PI/2){
                _pointTwo = Math.PI/2 + Math.PI*2 +_pointTwo;
            }else{
                _pointTwo = Math.PI/2 +_pointTwo;
            }

            var _arcID = "arc-" +_fatherID+"-"+(Number(_childID)-1);

            var _arc = d3.svg.arc()
                            .innerRadius(_r + 3)
                            .outerRadius(_r + 5)
                            .startAngle(_pointOne)
                            .endAngle(_pointTwo);
            
            d3.select("#" + DIV.chartDiv + " svg g").append("path")
                .attr("d", _arc)
                .attr('fill',_fill)
                .attr('id',_arcID)
                .attr('class','mark');
        }
    }
    
    //Initialize PieChart in DC.js
    function initDCPieChart() {
        var _pieWidth = 130,
            _pieRadius = (_pieWidth - 20) /2;

        
        pieChart = dc.pieChart("#" + DIV.chartDiv);
        
        cluster = ndx.dimension(function (d) {
            if(!d[selectedAttr] || d[selectedAttr].toLowerCase()==="unknown" 
                    || d[selectedAttr].toLowerCase()==="none")
                return "NA";
            return d[selectedAttr];
        });
        
        pieChart
            .width(_pieWidth)
            .height(_pieWidth)
            .radius(_pieRadius)
            .dimension(cluster)
            .group(cluster.group())
            .transitionDuration(600)
            .ordinalColors(chartColors)
            .label(function (d) {
                return d.value;
            })
            .ordering(function(d){ return d.key;});
    }
    
    //Initial Label Information stored in `label` array
    function initLabelInfo() {
        var _labelID = 0;

        label = [];
        
        $('#' + DIV.chartDiv + '>svg>g>g').each(function(){
            var _labelDatum = {},
                _labelName = $(this).find('title').text().split(':'),
                _color = $(this).find('path').attr('fill'),            
                _pointsInfo = $(this).find('path').attr('d').split(/[\s,MLHVCSQTAZ]/);            
            
            var _x1 = Number( _pointsInfo[1] ),
                _y1 = Number( _pointsInfo[2] ),
                _x2 = Number( _pointsInfo[8] ),
                _y2 = Number( _pointsInfo[9] );
                
            if(_x1 !== _x2 || _y1 !== _y2){
                _labelDatum.id = _labelID;
                _labelDatum.name = _labelName[0];
                _labelDatum.color = _color;
                _labelDatum.parentID = DIV.chartDiv;
                _labelDatum.value = _labelName[1];
                label.push(_labelDatum);
            }
            _labelID++;            
        });
    }
    
    //Initial global parameters by using passed object .
    function initParam(_param) {
        var _baseID = _param.baseID;
        
        className = _param.chartDivClass,
        chartID = _param.chartID;
        selectedAttr = _param.attrID;
        selectedAttrDisplay = _param.displayName;
        transitionDuration = _param.transitionDuration;
        ndx = _param.ndx;
        chartColors = _param.chartColors;

        DIV.mainDiv = _baseID + "-dc-chart-main-" + chartID;
        DIV.chartDiv = _baseID + "-dc-chart-" + chartID;
        DIV.labelTableID = "table-" + _baseID + "-dc-chart-" + chartID;
        DIV.labelTableTdID = "pieLabel-" + _baseID + "-dc-chart-" + chartID + "-";
        DIV.parentID = _baseID + "-charts";
    }

    //Remove drawed Pie Markder.
    function removeMarker() {
        $("#" + DIV.chartDiv).find('svg g .mark').remove();
    }
    
    //Called when the number of label biggen than 6, used by addPieLabels()
    function smallLabelFunction() {
        var _innerID = 0;
        $('#' + DIV.mainDiv)
                .find('.study-view-pie-label')
                .append("<table id="+DIV.labelTableID+"-0></table>");

        for(var i=0; i< label.length; i++){
            var _tmpName = label[i].name;
            
            if(_tmpName.length > 9){
                _tmpName = _tmpName.substring(0,5) + " ...";
            }
            
            if(i % 2 === 0){
                $('#' + DIV.mainDiv)
                        .find('table')
                        .append("<tr id="+ _innerID +" width='150px'></tr>");
                _innerID++;
            } 
            
            $('#' + DIV.mainDiv)
                    .find('table tr:nth-child(' + _innerID +')')
                    .append('<td class="pieLabel" id="'+
                        DIV.labelTableTdID +label[i].id+'-'+i+
                        '"  style="font-size:'+fontSize+'px">'+
                        '<svg width="75" height="13"><rect width="'+
                        labelSize+'" height="'+labelSize+'" style="fill:'+
                        label[i].color + ';" /><text x="15" y="10">'+
                        _tmpName+'</text></svg></td>');

            //Only add qtip when the length of pie label bigger than 9
            if(label[i].name.length > 9){
                var _qtip = jQuery.extend(true,{},StudyViewBoilerplate.pieLabelQtip);
                
                _qtip.content.text = label[i].name;
                $('#'+DIV.labelTableTdID +label[i].id+'-'+i).qtip(_qtip);
            }
        }
    }
    
    //Pie Chart will have communications with ScatterPlot, this function is used
    //to call the callback function.
    function updateScatterPlot(_currentFilters) {
        scatterPlotCallback(_currentFilters);
    }
    
    return {
        init: function(_param){
            initParam(_param);
            createDiv();
            initDCPieChart();
            addFunctions();
        },

        getChart : function(){
            return pieChart;
        },
        
        getCluster: function(){
            return cluster;
        },
        
        drawMarker: drawMarker,
        
        pieLabelClickCallbackFunction: function(_callback){
            pieLabelClickCallback = _callback;
        },
        
        removeMarker: removeMarker,
        
        scatterPlotCallbackFunction: function (_callback){
            scatterPlotCallback = _callback;
        },
        
        postFilterCallbackFunc: function(_callback){
            postFilterCallback = _callback;
        }
    };
};