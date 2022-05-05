$title Check Consistency of the WiNDC Accounting Models (MGE and MCP)

$if not set ds $set ds WiNDC_bluenote_cps_census_2017.gdx

* $if not set gdxdir $set gdxdir ..%system.dirsep%..%system.dirsep%..%system.dirsep%bmk_data%system.dirsep%
$if not set gdxdir $set gdxdir ..%system.dirsep%ltran%system.dirsep%

*	Solve sequence is either mcpmge or mgemcp

$if not set sequence $set sequence mgemcp

*	Read the dataset:

*$include windc_coredata
$include readdata_alt.gms

*	Read the MGE model:
*$exit

$include mgemodel

* $exit

*	Read the MCP model:

* $include mcpmodel

*	Replicate the benchmark equilibrium in both formats:

mgemodel.workspace = 100;
mgemodel.iterlim = 0;

*mgemodel.iterlim = 100000;
$INCLUDE MGEMODEL.GEN
SOLVE mgemodel using mcp;
ABORT$(mgemodel.objval > 1e-4) "Error in benchmark calibration of the MGE model.";

te(r,"col","ele") = 0.3;
* te(r,"gas",s) = 0.1;
* te(r,"oil",s) = 0.1;

* aeeir(r,"gas") = 0.1;

* tfp(r,"ele") = -0.1;

parameter rep;

rep(r,o,s,"LOS","BMK") = LOS.l(r,o,s);

rep(r,"all",s,"LDS","BMK") = sum(o,LOS.l(r,o,s));
rep(r,o,"all","LDO","BMK") = sum(s,LOS.l(r,o,s));

rep(r,o,s,"LSHR","BMK")$[(sum(s.local,ld0_o(r,o,s)))] = LOS.l(r,o,s)/sum(s.local,ld0_o(r,o,s));

rep(r,o,s,"LSHRbmk","BMK") = lshr(r,o,s);

mgemodel.iterlim = 100000;
$INCLUDE MGEMODEL.GEN
SOLVE mgemodel using mcp;
ABORT$(mgemodel.objval > 1e-4) "Error in benchmark calibration of the MGE model.";

rep(r,o,s,"LOS","CF") = LOS.l(r,o,s);

rep(r,"all",s,"LDS","CF") = sum(o,LOS.l(r,o,s));
rep(r,o,"all","LDO","CF") = sum(s,LOS.l(r,o,s));


rep(r,o,s,"LSHR","CF")$[(sum(s.local,ld0_o(r,o,s)))] = LOS.l(r,o,s)/sum(s.local,ld0_o(r,o,s));


rep(r,o,s,"LOSdiffraw","CF") =
	rep(r,o,s,"LOS","CF") - rep(r,o,s,"LOS","BMK");

rep(r,o,s,"LOSdiffpct","CF")$[rep(r,o,s,"LOS","BMK")] =
	100*(rep(r,o,s,"LOS","CF") - rep(r,o,s,"LOS","BMK"))/ rep(r,o,s,"LOS","BMK");

rep(r,"all",s,"LDSdiffraw","CF") =
	rep(r,"all",s,"LDS","CF") - rep(r,"all",s,"LDS","BMK");

rep(r,"all",s,"LDSdiffpct","CF")$[rep(r,"all",s,"LDS","BMK")] =
	(rep(r,"all",s,"LDS","CF") - rep(r,"all",s,"LDS","BMK"))/ rep(r,"all",s,"LDS","BMK");

rep(r,o,"all","LDOdiffraw","CF") =
	rep(r,o,"all","LDO","CF") - rep(r,o,"all","LDO","BMK");

rep(r,o,"all","LDOdiffpct","CF")$[rep(r,o,"all","LDO","BMK")] =
	(rep(r,o,"all","LDO","CF") - rep(r,o,"all","LDO","BMK"))/ rep(r,o,"all","LDO","BMK");

rep(r,o,s,"LSHRdiffraw","CF") =
	rep(r,o,s,"LSHR","CF") - rep(r,o,s,"LSHR","BMK");

rep(r,o,s,"LSHRdiffpct","CF")$[rep(r,o,s,"LSHR","BMK")] =
	(rep(r,o,s,"LSHR","CF") - rep(r,o,s,"LSHR","BMK"))/ rep(r,o,s,"LSHR","BMK");



rep("all",o,s,"LOSdiffraw","CF") =
	sum(r,rep(r,o,s,"LOS","CF") - rep(r,o,s,"LOS","BMK"));

rep("all",o,s,"LOSdiffpct","CF")$[sum(r,rep(r,o,s,"LOS","BMK"))] =
	100*(sum(r,(rep(r,o,s,"LOS","CF"))) - sum(r,rep(r,o,s,"LOS","BMK")))/ sum(r,rep(r,o,s,"LOS","BMK"));

execute_unload "ltran_out.gdx";

$exit
