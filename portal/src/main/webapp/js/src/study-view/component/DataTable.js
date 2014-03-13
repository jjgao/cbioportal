/*
 * DataTable component.
 * 
 * @param _data -- TWO PROPERTIES
 *                 attr: data attributes, each attribute is an object which
 *                       include attr_id, datatype, description and dispaly_name
 *                 arr: data contents, each object of this array stand for one
 *                      case. Each case should include all attributes as key and
 *                      their relative value. If attribute of this case doesn't
 *                      exist, the NA value will still give to this attribute.
 *                  
 * @interface: getDataTable -- return DataTable Object.
 * @interface: updateTable -- giving filteredResultList, the dateTable will be
 *                            refreshed.
 * @interface: rowClickCallback -- pass a function to dataTable. It will be
 *                                 called when one row is clicked.
 * @interface: rowShiftClickCallback -- pass a function to dataTable. It will be
 *                                      called when one row is clicked and
 *                                      ShiftKeys is pressed at same time.
 * @interface: resizeTable: will be used when width of dataTable changed.
 *                                       
 * @authur: Hongxin Zhang
 * @date: Mar. 2014
 * 
 */
        
var DataTable = function() {
    
    var attr,
        arr,
        attrLength,
        arrLength,
        aoColumnsLength,
        aaDataLength,
        dataTable,
        dataType = [],
        dataTableNumericFilter = [],
        disableFiltId = [0],
        aoColumns = [], //DataTable Title Data
        aaData = [], //DataTable Content Data
        columnIndexMappingColumnId = [];
    
    var rowClickCallback,
        rowShiftClickCallback;

    var parObject = {
        studyId: "",
        caseIds: "",
        cnaProfileId: "",
        mutationProfileId: "",
        caseSetId: ""
    };
    
    function initParam(_param, _data) {
        var i;
        
        attr = _data.attr;
        arr = _data.arr;
        
        attrLength = attr.length;
        arrLength = arr.length;
        
        parObject.studyId = _param.studyId;
        parObject.caseIds = _param.caseIds;
        parObject.cnaProfileId = _param.cnaProfileId;
        parObject.mutationProfileId = _param.mutationProfileId;
        parObject.caseSetId = _param.caseSetId;
        
        for( i = 0; i < attrLength; i++ ){
            if(attr[i]["datatype"] === "NUMBER"){
                dataType[attr[i]["attr_id"]] = 'allnumeric';
            }else{
                dataType[attr[i]["attr_id"]] = 'string';
            }
            dataTableNumericFilter[i] = '';
        }
        
        initColumnsTitleData();
        initContentData();
    }
    
    //Initialize aoColumns Data
    function initColumnsTitleData() {
        var i;
        
        aoColumns.length = 0;
        
        aoColumns.push({sTitle:"CASE ID",sType:'string'});

        for( i = 0; i < attr.length; i++ ){
            if( attr[i].attr_id !== 'CASE_ID' ){
                var _tmp = {};
                
                if(attr[i].attr_id === 'COPY_NUMBER_ALTERATIONS'){
                    _tmp.sTitle = 'CNA';
                }else{
                    _tmp.sTitle = attr[i].attr_id.replace(/[_]/g,' ');
                }
                
                _tmp.sType = dataType[attr[i].attr_id];
                aoColumns.push(_tmp);
            }
        }
        
        aoColumnsLength = aoColumns.length;
    }
    
    //Initialize aaData Data
    function initContentData() {
        aaData.length = 0;
        
        $.each(arr, function(key,value){ 
            
            aaData[key] = [];
            
            $.each(aoColumns, function(key1,value1){
                var _selectedString,
                    _specialCharLength,
                    _tmpValue ='',
                    _specialChar = ['(',')','/','?','+'];

                if(value1.sTitle === 'CNA'){
                    _tmpValue = value['COPY_NUMBER_ALTERATIONS'];                
                }else if ( value1.sTitle === 'COMPLETE (ACGH, MRNA, SEQUENCING)'){
                    _tmpValue = value[value1.sTitle];
                }else if ( value1.sTitle === 'CASE ID'){
                    _tmpValue = "<a href='tumormap.do?case_id=" + 
                    value['CASE_ID'] + "&cancer_study_id=" +
                    parObject.studyId + "' target='_blank'>" + 
                    value['CASE_ID'] + "</a></strong>";
                }else{
                    _tmpValue = value[value1.sTitle.replace(/[ ]/g,'_')];
                }
                if(!isNaN(_tmpValue) && (_tmpValue % 1 !== 0)){
                    _tmpValue = cbio.util.toPrecision(Number(_tmpValue),3,0.01);
                }
                
               
                _selectedString = _tmpValue.toString();
                _specialCharLength = _specialChar.length;
                
                if ( value1.sTitle !== 'CASE ID'){
                    var j;
                    
                    for( j = 0; j < _specialCharLength; j++){
                        if(_selectedString.indexOf(_specialChar[j]) !== -1){
                            var _re = new RegExp("\\" + _specialChar[j], "g");
                            
                            _selectedString = _selectedString.replace(_re, _specialChar[j] + " ");
                        } 
                    }
                }
                aaData[key].push(_selectedString);
            });
        });
        
        aaDataLength = aaData.length;
    }
    
    //Initialize the basic dataTable component by using jquery.dataTables.min.js
    function initDataTable() {
        dataTable = $('#dataTable').dataTable({
            "sScrollX": "1200px",
            "sScrollY": "300px",
            "bPaginate": false,
            "bScrollCollapse": true,
            "aoColumns": aoColumns,
            "aaData": aaData,
            "sDom": '<"dataTableTop"Ci<"dataTableReset">f>rt',
            "fnDrawCallback": function(oSettings){
                if($(".ColVis_collection.TableTools_collection").css('display') === 'block'){
                    var _currentIndex= 0 ;
                    
                    columnIndexMappingColumnId = {};
                    $.each(dataTable.fnSettings().aoColumns, function(c){
                        if(dataTable.fnSettings().aoColumns[c].bVisible === true){
                            columnIndexMappingColumnId[_currentIndex] = c;
                            _currentIndex++;
                        }
                    });
                    $("#dataTable_filter label input").attr("value","");
                    $.fn.dataTableExt.afnFiltering = [];
                    disableFiltId = [0];
                    resizeLeftColumn();         
                    refreshSelectionInDataTable();
                    $(".dataTableReset span").css('display','none');
                }
            }
        });
    }
    
    //Add th tags based on number of attributes
    function initDataTableTfoot() {
        for( var i = 0; i < aoColumnsLength; i++ ){
            $("#dataTable tfoot tr").append("<th></th>");
            columnIndexMappingColumnId[i] = i;
        }
    }
    
    //Add all HTML events by using JQUERY
    function addEvents() {
        
        $(".ColVis_MasterButton").css({height:'auto', padding:'0 4px'});
        
        $(".ColVis_MasterButton").click(function() {
            $('.ColVis_collection.TableTools_collection')
                .find('button')
                .first()
                .prop('disabled',true);
        
            $('.ColVis_collection.TableTools_collection')
                .find('button')
                .first()
                .find('input')
                .prop('disabled',true);
         
        });
        
        $(".dataTableReset")
            .append("<a><span class='hidden' title='Reset Chart'>RESET</span></a>");
    
        $(".dataTableReset span").click(function(){
            $("#dataTable_filter label input").attr("value","");
            $.fn.dataTableExt.afnFiltering = [];
            updateTable([]);
            resizeLeftColumn();            
            refreshSelectionInDataTable();
            dataTable.fnAdjustColumnSizing();
            showDataTableReset();
        });
        
        $("#dataTable_filter label input").keyup(function() {
            showDataTableReset(dataTable);
            resizeLeftColumn();
        });
        
        $("#dataTable tbody").mousedown(function(event){
            var _selectedRowCaseId = [],
                _deSelect = false;
                
            if(event.shiftKey){
                event.preventDefault();

                if($(event.target.parentNode).hasClass('row_selected')){
                    $(event.target.parentNode).removeClass('row_selected');
                    if($(event.target.parentNode).hasClass('odd')){
                       $(event.target.parentNode).css('background-color','#E2E4FF'); 
                    }else{
                        $(event.target.parentNode).css('background-color','white');
                    }
                }else{
                    $(event.target.parentNode).addClass('row_selected');
                    $(event.target.parentNode).css('background-color','lightgray');
                }
                
                _selectedRowCaseId = getRowSelectedCases();
                rowShiftClickCallback(_selectedRowCaseId);
                
            }else{
                if($(event.target.parentNode).hasClass('row_selected')){
                    $(event.target.parentNode).removeClass('row_selected');
                    if($(event.target.parentNode).hasClass('odd')){
                       $(event.target.parentNode).css('background-color','#E2E4FF'); 
                    }else{
                        $(event.target.parentNode).css('background-color','white');
                    }
                    _deSelect = true;
                }else{
                    $(dataTable.fnSettings().aoData).each(function (){
                        if($(this.nTr).hasClass('row_selected')){
                            $(this.nTr).removeClass('row_selected');
                            if($(this.nTr).hasClass('odd')){
                               $(this.nTr).css('background-color','#E2E4FF'); 
                            }else{
                                $(this.nTr).css('background-color','white');
                            }
                        }
                    });

                    $(event.target.parentNode).addClass('row_selected');
                    $(event.target.parentNode).css('background-color','lightgray');
                }
                
                _selectedRowCaseId = getRowSelectedCases();
                rowClickCallback(_deSelect, _selectedRowCaseId);
            }
        });
    }
    
    function getRowSelectedCases() {
        var _selectedRowCaseIds = [],
            _returnValue = fnGetSelected(),
            _returnValueLength = _returnValue.length;
        
        for( var i = 0; i < _returnValueLength; i++ ){
            _selectedRowCaseIds.push($(_returnValue[i]).find('td').first().text());
        }
        
        return _selectedRowCaseIds;
    }
    
    function updateTable(_exceptionColumns){
        var _oSettings = dataTable.fnSettings();
        
        for(var iCol = 0; iCol < _oSettings.aoPreSearchCols.length; iCol++) {
            if(_exceptionColumns.indexOf(iCol) === -1){
                _oSettings.aoPreSearchCols[ iCol ].sSearch = '';
            }
        }
        _oSettings.oPreviousSearch.sSearch = '';
        dataTable.fnDraw();
    }
    
    //Return the selected nodes
    function fnGetSelected(){
        var i,
            _aReturn = [],
            _aTrs = dataTable.fnGetNodes();

        for ( i = 0; i < _aTrs.length; i++ ){
                if ( $(_aTrs[i]).hasClass('row_selected') ){
                        _aReturn.push( _aTrs[i] );
                }
        }
        return _aReturn;
    }
    
    //Create Regular Selector or Numeric Selector based on data type.
    function fnCreateSelect( aData, index ){
        var _isNumericArray = true,
            _hasNullValue = false;
    
        for(var i=0;i<aData.length;i++){
            if(isNaN(aData[i])){
                if(aData[i] !== 'NA'){
                    _isNumericArray = false;
                    break;
                }else{
                    _hasNullValue = true;
                }
            }
        }
        
        if(_isNumericArray && _hasNullValue){
            var _index = aData.indexOf("NA");
            if (_index > -1) {
                aData.splice(_index, 1);
            }
        }
        
        if(_isNumericArray){            
            aData.sort(function(a,b) {
                return Number(a) - Number(b);
            });
        }else{
            aData.sort();
        }
        if(!_isNumericArray || aData.length === 0){
            var r='<select><option value=""></option>', i, iLen=aData.length;
            for ( i=0 ; i<iLen ; i++ )
            {
                r += '<option value="'+aData[i]+'">'+aData[i]+'</option>';
            }
            return r+'</select>';
        }else{
            var _min = aData[0],
                _max = aData[aData.length-1],
                _x1 = 5,
                _x2 = 65;
            
            var _leftTriangelCoordinates = (_x1-5) + ",2 "+ (_x1+5)+",2 "+_x1+",10",
                _rightTriangelCoordinates = (_x2-5) + ",2 "+ (_x2+5)+",2 "+_x2+",10",
                _leftText = "x='"+(_x1-3)+"' y='20'",
                _rightText = "x='"+(_x2-10)+"' y='20'",
                _resetText = "x='"+(_x2+15)+"' y='20'";
           
            var _svgLine = "<svg width='110' height='30' start='"+ _min +"' end='"+ _max +"'>" + 
                    "<g><line x1='"+ _x1 +"' y1='10' x2='"+ _x2 +"' y2='10' style='stroke:black;stroke-width:2' /></g>"+
                    "<g id='dataTable-"+ index +"-left' class='clickable left'><polygon points='"+_leftTriangelCoordinates+"' style='fill:grey'></polygon>"+
                    "<text "+_leftText+" fill='black' style='font-size:8'>"+ _min +"</text></g>" + 
                    "<g id='dataTable-"+ index +"-right' class='clickable right'><polygon points='"+_rightTriangelCoordinates+"' style='fill:grey'></polygon>"+
                    "<text "+_rightText+" fill='black' style='font-size:8'>"+ _max +"</text></g>" +
                    "<text "+ _resetText +" id='dataTable-"+ index +"-reset' class='clickable hidden'  fill='black' style='font-size:8'>RESET</text>" + 
                    "</svg>";

            return _svgLine;
        }
    }
    
    //This function will be called when the dataTable has been resized
    function resizeLeftColumn(){
        var _heightBody = $(".dataTables_scrollBody").css('height'),
            _heightTable = $('.dataTables_scroll').css('height'),
            _widthBody = $("#dataTable tbody>tr>td:nth-child(1)").css('width');
        
        _widthBody = _widthBody.slice(0,_widthBody.length-2);
        _widthBody = Number(_widthBody) + 20;
        _widthBody = _widthBody.toString() + 'px';
        
        $(".DTFC_LeftBodyLiner").css('height',_heightBody);
        $(".DTFC_LeftBodyWrapper").css('height',_heightBody); 
        $(".DTFC_LeftWrapper").css('width',_widthBody);
        $(".DTFC_LeftBodyLiner").css('width',_widthBody);
        $('.DTFC_ScrollWrapper').css('height',_heightTable);            
    }
    
    function showDataTableReset( ){
        var _showedColumnNumber = dataTable.fnSettings().fnRecordsDisplay(),
            _totalColumnNumber = dataTable.fnSettings().fnRecordsTotal();
        
        if(_showedColumnNumber !== _totalColumnNumber){
            $(".dataTableReset span").css('display','block');
        }else{
            $(".dataTableReset span").css('display','none');
            disableFiltId = [0];
            refreshSelectionInDataTable();
        }
    }
    
    function updateDataTableNumericFilter(){
        var i,
            _dataTableNumericFilterLength = dataTableNumericFilter.length;
        
        $.fn.dataTableExt.afnFiltering = [];
        for( i = 0; i < _dataTableNumericFilterLength; i++ ){
            if(dataTableNumericFilter[i] !== ''){
                $.fn.dataTableExt.afnFiltering.push(dataTableNumericFilter[i]);
            }
        }
        dataTable.fnDraw();
    }
    
    function refreshSelectionInDataTable(){
        $(".dataTables_scrollFoot tfoot th").each( function ( i ) {
            if(disableFiltId.indexOf(i) === -1){                
                $(this).css('z-index','1500');
                this.innerHTML = fnCreateSelect( dataTable.fnGetColumnData(columnIndexMappingColumnId[i]), i);
                
                var _drag = d3.behavior.drag()
                        .on("drag", dragMove)
                        .on("dragend",dragEnd);
                
                d3.select("#dataTable-" + i + "-left")
                                .call(_drag);
                d3.select("#dataTable-" + i + "-right")
                                .call(_drag);
                        
                $("#dataTable-" + i + "-reset").unbind('click');
                $("#dataTable-" + i + "-reset").click(function(){
                    dataTableNumericFilter[columnIndexMappingColumnId[i]] = '';
                    updateDataTableNumericFilter();
                    disableFiltId.splice(disableFiltId.indexOf(i),1);
                    resizeLeftColumn();
                    showDataTableReset();
                    $("#dataTable-" + i + "-reset").css('display','none');
                    refreshSelectionInDataTable();
                    dataTable.fnAdjustColumnSizing();
                });
                
                $('select', this).change( function () {
                    if($(this).val() === ''){
                        dataTable.fnFilter($(this).val(), columnIndexMappingColumnId[i]);
                        disableFiltId.splice(disableFiltId.indexOf(i),1);
                    }else{
                        var j,
                            _selectedString = $(this).val().toString(),
                            _specialChar = ['(',')','/','?','+'],
                            _specialCharLength = _specialChar.length;
                        for( j = 0; j < _specialCharLength; j++){
                            if(_selectedString.indexOf(_specialChar[j]) !== -1){
                                var re = new RegExp("\\" + _specialChar[j],"g");
                                _selectedString = _selectedString.replace(re ,"\\"+ _specialChar[j]);
                            } 
                        }
                        dataTable.fnFilter("^"+_selectedString+"$", columnIndexMappingColumnId[i], true);
                        disableFiltId.push(i);
                    }
                    
                    resizeLeftColumn();
                    showDataTableReset();
                    refreshSelectionInDataTable();
                    dataTable.fnAdjustColumnSizing();
                });
                
                //Used by numeric selector
                function dragMove(d) {
                    var _start = Number($(this).parent().attr('start')),
                        _end = Number($(this).parent().attr('end')),
                        _xMoved = d3.event.x - 5,
                        _lineLength = Number($(this).parent().find('line').attr('x2')) - Number($(this).parent().find('line').attr('x1'));
                    
                    if(_start > _end){
                        var _tmp = _end;
                        
                        _end = _start;
                        _start = _tmp;
                    }
                    
                    if(_xMoved >= 0 && _xMoved <= _lineLength){
                        var _text = (_end-_start) * _xMoved / _lineLength + _start;
                        
                        _text = cbio.util.toPrecision(Number(_text),3,0.1);
                        
                        if($(this).attr('id').toString().indexOf('left') !== -1){
                            d3.select(this)
                                .attr("transform", "translate(" +_xMoved +",0)");
                        }else{
                            _xMoved -= _lineLength;
                            d3.select(this)
                                .attr("transform", "translate(" +_xMoved +",0)");
                        }
                        $(this).find('text').text(_text);
                    }
                }
                
                //Used by numeric selector
                function dragEnd(d) {
                    var _min = Number($(this).parent().find('g.left').find('text').text()),
                        _max = Number($(this).parent().find('g.right').find('text').text());
                    
                    if(_min > _max){
                        var _tmp = _max;
                        
                        _max = _min;
                        _min = _tmp;
                    }
                    
                    dataTableNumericFilter[columnIndexMappingColumnId[i]] = function( oSettings, aData, iDataIndex ) {
                        var _iMin = _min,
                            _iMax = _max,
                            _iCurrent = aData[columnIndexMappingColumnId[i]];
                    
                        if ( _iMin === "" && _iMax === "" ){
                                return true;
                        }else if ( _iMin === "" && _iCurrent < _iMax ){
                                return true;
                        }else if ( _iMin < _iCurrent && "" === _iMax ){
                                return true;
                        }else if ( _iMin < _iCurrent && _iCurrent < _iMax ){
                                return true;
                        }
                        
                        return false;
                    };
                        
                    updateDataTableNumericFilter();
                    dataTable.fnSort([ [columnIndexMappingColumnId[i],'asc']]);
                    disableFiltId.push(i);
                    resizeLeftColumn();
                    showDataTableReset();
                    $("#dataTable-" + i + "-reset").css('display','block');
                    refreshSelectionInDataTable();
                    dataTable.fnAdjustColumnSizing();
                }
            }
        });
    }
    
    function deleteChartResetDataTable(_filteredResult) {
        var _filterArray = [];

        for(var i=0 ; i<_filteredResult.length ; i++){
            _filterArray.push(_filteredResult[i].CASE_ID);
        }

        $.fn.dataTableExt.afnFiltering = [function( oSettings, aData, iDataIndex ) {
            var _data = aData[0],
                _dataContent = $(_data).text();
                
            if ( _filterArray.indexOf(_dataContent) !== -1){
                return true;
            }
            return false;
        }];
        dataTable.fnDraw();
    }
    
    function resizeTable(){
        //Before resize data table, the window should be showed first
        $('#dc-plots-loading-wait').hide();
        $('#study-view-main').show();
        
        refreshSelectionInDataTable();
        dataTable.fnAdjustColumnSizing();
        new FixedColumns(dataTable);
        $(".DTFC_LeftBodyLiner").css("overflow-y","hidden");
        $(".dataTables_scroll").css("overflow-x","scroll");
        $(".DTFC_LeftHeadWrapper").css("background-color","white");
    }
    
    return {
        init: function(_param, _data) {
            initParam(_param, _data);
            initDataTableTfoot();
            initDataTable();
            addEvents();
        },
        
        getDataTable: function() {
            return dataTable;
        },
        
        updateTable: function(_filteredResult) {
            deleteChartResetDataTable(_filteredResult);
            resizeLeftColumn();            
            refreshSelectionInDataTable();
            dataTable.fnAdjustColumnSizing();
            $("#dataTable_filter label input").attr("value","");
            showDataTableReset();
        },
        
        rowClickCallback: function(_callback) {
            rowClickCallback = _callback;
        },
        
        rowShiftClickCallback: function(_callback) {
            rowShiftClickCallback = _callback;
        },
        
        resizeTable: resizeTable
    };
};