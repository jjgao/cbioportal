<html>
<head>
<link href="css/global_portal.css" type="text/css" rel="stylesheet" />

<style type="text/css">
textarea {
  border:1px solid #999999;
  font-family:Consolas,Monaco,Lucida Console,Liberation Mono,DejaVu Sans Mono,Bitstream Vera Sans Mono,Courier New, monospace;
}
</style>
</head>
<body>

<div align="left">
<h2>OncoViz</h2>
<form action="onco_viz.json" METHOD="POST">

    <p>Format is as follows:</p>
    <pre>
CASE_ID CHR START END REFERENCE_ALLELE OBSERVED_ALLELE COLOR
    </pre>

    <p>
    <b>Warning! </b>  This is alpha code.  If you forget a column
        or format your data incorrectly, you are likely to get bad a messy NullPointerException!
    </p>
    
    Enter a gene:
    <input type="text" name="hugoGeneSymbol" value="PIK3CA">

    <P>Paste Tab-Delim Mutations (restricted to single gene!):

    <BR>
    <textarea name="mutations" rows="20" cols="60">
TCGA-AP-A051	3	178916537	178916537	G	T	red
TCGA-BS-A0UM	3	178916614	178916614	A	G	blue
TCGA-B5-A0JY	3	178916648	178916648	G	A	red
TCGA-A5-A0GI	3	178916725	178916725	C	A	blue
TCGA-AX-A0J0	3	178916725	178916725	C	T	red
TCGA-BS-A0TC	3	178916725	178916725	C	T	red
TCGA-A5-A0GA	3	178916726	178916726	G	A	blue
TCGA-AP-A0LD	3	178916726	178916726	G	A	blue
TCGA-D1-A163	3	178916726	178916726	G	A	blue
TCGA-B5-A11R	3	178916728	178916728	G	A	blue
TCGA-B5-A11X	3	178916728	178916728	G	A	green
TCGA-D1-A103	3	178916728	178916728	G	A	red
TCGA-BG-A0M4	3	178916729	178916729	A	G	blue
TCGA-AX-A0J0	3	178916823	178916823	C	T	red
TCGA-AP-A052	3	178916854	178916854	G	A	green
TCGA-A5-A0GP	3	178916876	178916876	G	A	red
TCGA-A5-A0VP	3	178916876	178916876	G	A	red
TCGA-AP-A056	3	178916876	178916876	G	A	red
TCGA-AX-A06B	3	178916876	178916876	G	A	green
TCGA-AX-A06F	3	178916876	178916876	G	A	red
TCGA-BG-A0VX	3	178916876	178916876	G	A	blue
TCGA-BS-A0TJ	3	178916876	178916876	G	A	blue
TCGA-BS-A0UV	3	178916876	178916876	G	A	red
TCGA-D1-A16F	3	178916876	178916876	G	A	blue
TCGA-D1-A16R	3	178916876	178916876	G	A	green
TCGA-D1-A174	3	178916876	178916876	G	A	blue
TCGA-D1-A17T	3	178916876	178916876	G	A	green
TCGA-AX-A0IZ	3	178916890	178916890	C	T	red
TCGA-B5-A0K2	3	178916890	178916890	C	T	blue
TCGA-AP-A0LP	3	178916891	178916891	G	A	blue
TCGA-B5-A11I	3	178916891	178916891	G	A	green
TCGA-B5-A11Y	3	178916891	178916891	G	A	red
TCGA-BG-A0MQ	3	178916891	178916891	G	A	blue
TCGA-BG-A187	3	178916891	178916891	G	A	green
TCGA-D1-A16E	3	178916917	178916919	ATT	-	green
TCGA-BG-A18C	3	178916928	178916933	AGGCAA	-	green
TCGA-D1-A175	3	178916930	178916930	G	A	red
TCGA-BG-A0M9	3	178916933	178916933	A	G	green
TCGA-A5-A0GB	3	178916936	178916936	G	A	blue
TCGA-AP-A0LD	3	178916936	178916936	G	A	blue
TCGA-B5-A0K3	3	178916936	178916936	G	A	green
TCGA-BS-A0TC	3	178916936	178916936	G	A	red
TCGA-AP-A0LQ	3	178916938	178916940	GAA	-	green
TCGA-B5-A0K6	3	178916938	178916940	GAA	-	blue
TCGA-B5-A11G	3	178916938	178916940	GAA	-	blue
TCGA-AP-A053	3	178916944	178916944	A	G	green
TCGA-AX-A064	3	178916944	178916944	A	G	blue
TCGA-BS-A0V6	3	178916944	178916944	A	G	green
TCGA-BG-A0VZ	3	178916946	178916946	G	T	blue
TCGA-AX-A06D	3	178916957	178916957	G	C	blue
TCGA-A5-A0RA	3	178916957	178916957	G	T	green
TCGA-A5-A0GP	3	178917478	178917478	G	A	red
TCGA-AX-A06J	3	178917478	178917478	G	A	green
TCGA-B5-A11X	3	178917478	178917478	G	A	green
TCGA-BS-A0WQ	3	178917478	178917478	G	A	green
TCGA-D1-A17F	3	178917478	178917478	G	A	blue
TCGA-D1-A17Q	3	178917478	178917478	G	A	red
TCGA-AX-A0IZ	3	178917496	178917496	C	A	red
TCGA-AP-A0LM	3	178917634	178917634	A	G	red
TCGA-AP-A0LM	3	178921353	178921353	C	A	red
TCGA-B5-A11R	3	178921520	178921520	C	T	blue
TCGA-AP-A059	3	178921544	178921544	C	A	red
TCGA-AP-A0LE	3	178921548	178921548	G	A	blue
TCGA-B5-A11Y	3	178921548	178921548	G	A	red
TCGA-AX-A060	3	178921549	178921549	T	C	blue
TCGA-AP-A05D	3	178921552	178921552	A	T	green
TCGA-B5-A0K0	3	178921553	178921553	T	A	green
TCGA-BG-A186	3	178921553	178921553	T	A	green
TCGA-BS-A0V7	3	178921553	178921553	T	A	green
TCGA-AP-A0LM	3	178921566	178921566	G	A	red
TCGA-B5-A11H	3	178921567	178921567	A	G	red
TCGA-A5-A0GB	3	178922306	178922306	G	C	blue
TCGA-AP-A0LH	3	178922321	178922321	G	A	green
TCGA-BG-A0M0	3	178922324	178922324	G	A	blue
TCGA-B5-A0K2	3	178922363	178922363	T	C	blue
TCGA-A5-A0R6	3	178922364	178922364	G	A	green
TCGA-AX-A06J	3	178922364	178922364	G	T	green
TCGA-BS-A0UF	3	178927439	178927439	G	A	red
TCGA-BG-A187	3	178927980	178927980	T	C	green
TCGA-D1-A15X	3	178927980	178927980	T	C	red
TCGA-D1-A16N	3	178927980	178927980	T	C	green
TCGA-D1-A17H	3	178927980	178927980	T	C	blue
TCGA-B5-A0JY	3	178928068	178928068	C	T	red
TCGA-BK-A0C9	3	178928079	178928079	G	C	blue
TCGA-AP-A0LM	3	178928082	178928082	G	T	red
TCGA-BG-A0MT	3	178928110	178928115	GATCAA	-	green
TCGA-BG-A0M0	3	178928226	178928226	C	T	blue
TCGA-BG-A18B	3	178928297	178928297	C	T	blue
TCGA-B5-A0JY	3	178936023	178936023	A	C	red
TCGA-A5-A0GG	3	178936082	178936082	G	A	blue
TCGA-A5-A0GM	3	178936082	178936082	G	A	green
TCGA-A5-A0GX	3	178936082	178936082	G	A	green
TCGA-A5-A0VQ	3	178936082	178936082	G	A	blue
TCGA-AX-A062	3	178936082	178936082	G	A	green
TCGA-B5-A0JV	3	178936082	178936082	G	A	blue
TCGA-BK-A0CB	3	178936082	178936082	G	A	green
TCGA-D1-A0ZP	3	178936082	178936082	G	A	green
TCGA-D1-A0ZU	3	178936082	178936082	G	A	green
TCGA-D1-A17K	3	178936082	178936082	G	A	green
TCGA-BK-A0C9	3	178936082	178936082	G	C	blue
TCGA-BG-A0MC	3	178936083	178936083	A	C	green
TCGA-D1-A16S	3	178936083	178936083	A	C	green
TCGA-A5-A0VP	3	178936083	178936083	A	T	red
TCGA-AP-A0LN	3	178936091	178936091	G	A	green
TCGA-AP-A0LV	3	178936091	178936091	G	A	green
TCGA-B5-A0K7	3	178936091	178936091	G	A	green
TCGA-BG-A0RY	3	178936091	178936091	G	A	green
TCGA-BG-A0W2	3	178936091	178936091	G	A	green
TCGA-BS-A0TG	3	178936091	178936091	G	A	green
TCGA-D1-A168	3	178936091	178936091	G	A	green
TCGA-D1-A16O	3	178936091	178936091	G	A	green
TCGA-BS-A0UT	3	178936092	178936092	A	C	green
TCGA-D1-A17L	3	178936092	178936092	A	G	green
TCGA-BG-A0M6	3	178936093	178936093	G	T	green
TCGA-D1-A102	3	178936093	178936093	G	T	green
TCGA-AX-A05W	3	178936094	178936094	C	A	green
TCGA-B5-A11O	3	178936094	178936094	C	A	green
TCGA-D1-A16G	3	178936094	178936094	C	A	green
TCGA-A5-A0GN	3	178936095	178936095	A	C	green
TCGA-AP-A0LJ	3	178936095	178936095	A	C	green
TCGA-BG-A0VX	3	178936095	178936095	A	C	blue
TCGA-BS-A0UJ	3	178936095	178936095	A	C	red
TCGA-D1-A0ZO	3	178936095	178936095	A	G	blue
TCGA-D1-A15V	3	178936095	178936095	A	G	green
TCGA-D1-A17A	3	178936095	178936095	A	G	blue
TCGA-A5-A0G2	3	178937024	178937024	C	A	red
TCGA-AX-A05Z	3	178937046	178937046	C	A	red
TCGA-AX-A060	3	178937422	178937422	T	C	blue
TCGA-D1-A103	3	178937461	178937461	C	T	red
TCGA-A5-A0G1	3	178937826	178937826	C	A	red
TCGA-D1-A17F	3	178938931	178938931	G	A	blue
TCGA-AP-A0LM	3	178943785	178943785	C	T	red
TCGA-A5-A0G2	3	178947119	178947119	G	A	red
TCGA-AP-A054	3	178947161	178947161	T	G	red
TCGA-D1-A175	3	178947200	178947200	A	G	red
TCGA-D1-A16R	3	178947827	178947827	G	T	green
TCGA-D1-A17T	3	178947827	178947827	G	T	green
TCGA-A5-A0G2	3	178947833	178947833	G	A	red
TCGA-BS-A0V6	3	178948044	178948044	A	G	green
TCGA-D1-A161	3	178948044	178948044	A	G	green
TCGA-B5-A11I	3	178951957	178951957	G	T	green
TCGA-B5-A11S	3	178951964	178951964	G	C	green
TCGA-BS-A0UF	3	178951992	178951992	T	G	red
TCGA-AP-A056	3	178952007	178952007	A	G	red
TCGA-BG-A0MQ	3	178952007	178952007	A	G	blue
TCGA-AX-A05Z	3	178952018	178952018	A	G	red
TCGA-AX-A0J0	3	178952018	178952018	A	T	red
TCGA-A5-A0R6	3	178952072	178952072	A	G	green
TCGA-B5-A11S	3	178952072	178952072	A	G	green
TCGA-BG-A0W1	3	178952072	178952072	A	G	green
TCGA-AP-A054	3	178952074	178952074	G	A	red
TCGA-B5-A11F	3	178952074	178952074	G	A	green
TCGA-D1-A16J	3	178952074	178952074	G	A	blue
TCGA-BG-A0VZ	3	178952077	178952077	T	G	blue
TCGA-BS-A0UM	3	178952084	178952084	C	T	blue
TCGA-AP-A0LT	3	178952085	178952085	A	G	blue
TCGA-AX-A05T	3	178952085	178952085	A	G	green
TCGA-B5-A0JS	3	178952085	178952085	A	G	green
TCGA-B5-A0JT	3	178952085	178952085	A	G	green
TCGA-BG-A18A	3	178952085	178952085	A	G	green
TCGA-BS-A0UL	3	178952085	178952085	A	G	blue
TCGA-D1-A0ZV	3	178952085	178952085	A	G	green
TCGA-D1-A162	3	178952085	178952085	A	G	green
TCGA-D1-A16B	3	178952085	178952085	A	G	green
TCGA-D1-A179	3	178952085	178952085	A	G	green
TCGA-D1-A17M	3	178952085	178952085	A	G	blue
TCGA-D1-A17S	3	178952085	178952085	A	G	green
TCGA-AP-A0L9	3	178952085	178952085	A	T	green
TCGA-AP-A0LG	3	178952085	178952085	A	T	blue
TCGA-BK-A0CA	3	178952085	178952085	A	T	green
TCGA-BS-A0V7	3	178952085	178952085	A	T	green
TCGA-D1-A16D	3	178952085	178952085	A	T	green
TCGA-B5-A11Y	3	178952086	178952086	T	A	red
TCGA-BG-A0M8	3	178952088	178952088	A	G	green
TCGA-B5-A11M	3	178952090	178952090	G	C	green
TCGA-A5-A0G2	3	178952119	178952119	C	T	red
    </textarea>
    <P>
    <input type="submit">
</form>
</div>
</body>
</html>