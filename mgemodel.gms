$title MPSGE Model with Pooled National Markets
 
* set o used for origin sector
alias(s,o,oo);

parameter lshr(r,o,s) 		labor transition rate from origin (o) to destination sector (s);
parameter numberl(r,o,s)	randomization of transition rates;
parameter popr(r,o)			household disagg share by o;

* Random number for toy transition rates
numberl(r,o,s)$ld0(r,s) = uniform(0.9,1.1);
* Rescale according to the size of sectoral labor demand
numberl(r,o,s)$(ld0(r,s)) = numberl(r,o,s)*ld0(r,s)/sum((s.local),ld0(r,s));

* Assume 75% same sector and distribute non same sector across remaining 25%
loop((r,o),
lshr(r,o,s)$(not sameas(o,s)) =
	0.25*numberl(r,o,s)$(not sameas(o,s))/sum(s.local,numberl(r,o,s)$(not sameas(o,s)));
);

lshr(r,o,s)$(sameas(o,s)) = 0.75;

* Solve for benchmark labor supply in origin sector o using transition rate and labor demand
variables
	VFO(r,o)	labor supply in origin sector o
	MN			dummy
;

equations
	main(r,s)	
	dummy
;

main(r,s).. 	ld0(r,s) =e= sum(o, VFO(r,o)*lshr(r,o,s));
dummy..			MN =e= 0;

model getvfo	/main, dummy/;

MN.l = 0;
VFO.lo(r,o) = 0;

solve getvfo minimizing MN using lp;

* used to disaggregate househols
parameter	oshr(r,o)	share of origin sector o in region r;
oshr(r,o) = VFO.l(r,o)/sum((o.local),VFO.l(r,o));
popr(r,o) = oshr(r,o);

* Origin specific labor demand enters model to inform supply and demand
parameter ld0_o(r,o,s)	origin specific labor demand;
ld0_o(r,o,s) = VFO.l(r,o)*lshr(r,o,s);

* Preprocessing for household disaggregation by o
parameter bmk_trev(r)	benchmark tax revenues;

bmk_trev(r) = sum((s,g)$y_(r,s), ty0(r,s) * ys0(r,s,g))
	+ sum((g)$a_(r,g),   ta0(r,g)*a0(r,g) + tm0(r,g)*m0(r,g));
*	+ taxrevL(r) + sum(s,tk0(r)*kd0(r,s));

parameter tpo(r,o)	benchmark transfers;

tpo(r,o) = (c0(r) - (
	sum(g,yh0(r,g))
	+ bopdef0(r)
	+ hhadj(r)
	- sum(g,i0(r,g))
	+ sum(s,ld0(r,s))
	+ sum(s,kd0(r,s))
	))*popr(r,o)
;

parameter gov0	government deficit;
gov0 = sum(r,sum(g,g0(r,g)) + sum(o,tpo(r,o)) - bmk_trev(r));


* policy counterfactuals
parameter te(r,g,s)	energy goods tax;

te(r,g,s) = 0;

parameter aeeir(r,g)	aeei rate;

aeeir(r,g) = 0;


parameter tfp(r,s)	tfp shock;

tfp(r,s) = 0;


* $exit

																						
* * New/updated model variables
* $sectors:
* 	C(r,o)$popr(r,o)		! Aggregate final demand
* 	LS(r,o)$popr(r,o)		! Labor Supply

* $commodities:
* 	PC(r,o)$popr(r,o)		! Consumer price index
* 	PL(r,o,s)$ld0_o(r,o,s)	! Wage
* 	PELL(r,o)$popr(r,o)		! Opportunity cost of labor
	
* $consumer:
* 	RA(r,o)$popr(r,o)		! Representative agent
* 	GOVT					! Government

* $auxiliary:
* 	TRANS					! Transfers

* $prod:Y(r,s)$y_(r,s)  s:0 va:1
* 	o:PY(r,g)		q:(ys0(r,s,g)*(1-tfp(r,s)))		a:GOVT 	t:ty(r,s)	p:(1-ty0(r,s))
* 	i:PA(r,g)		q:(id0(r,g,s)*(1-aeeir(r,g)))	a:GOVT	t:te(r,g,s)
* 	i:PK(r,s)		q:kd0(r,s)		va:	
* 	i:PL(r,o,s)		q:ld0_o(r,o,s)	va:

* $prod:LS(r,o)$popr(r,o)	t:1
* 	o:PL(r,o,s)		q:ld0_o(r,o,s)
* 	i:PELL(r,o)		q:(sum(s,ld0_o(r,o,s)))

* $prod:C(r,o)$popr(r,o)  s:1
*     o:PC(r,o)		q:(c0(r)*popr(r,o))
* 	i:PA(r,g)		q:(cd0(r,g)*popr(r,o))

* $demand:RA(r,o)$popr(r,o)
* 	d:PC(r,o)		q:(c0(r)*popr(r,o))
* 	e:PY(r,g)		q:(yh0(r,g)*popr(r,o))
* 	e:PFX			q:((bopdef0(r) + hhadj(r))*popr(r,o))
* 	e:PA(r,g)		q:(-(i0(r,g))*popr(r,o))
* 	e:PELL(r,o)		q:((sum((s,o.local),ld0_o(r,o,s)))*popr(r,o))
* 	e:PK(r,s)		q:(kd0(r,s)*popr(r,o))
* 	e:PFX			q:tpo(r,o)	r:TRANS
	
* $demand:GOVT
* 	d:PA(r,g)		q:g0(r,g)
* 	e:PFX			q:gov0
* 	e:PFX			q:(-sum((r,o),tpo(r,o)))	r:TRANS

* $constraint:TRANS
* 	GOVT =e= sum((r,g),PA(r,g)*g0(r,g));

* $report:
* 	v:LOS(r,o,s)	o:PL(r,o,s)	prod:LS(r,o)



$ONTEXT
$model:mgemodel

$sectors:
	Y(r,s)$y_(r,s)		!	Production
	X(r,g)$x_(r,g)		!	Disposition
	A(r,g)$a_(r,g)		!	Absorption
	C(r,o)$popr(r,o)			!	Aggregate final demand
	MS(r,m)			!	Margin supply
	LS(r,o)$popr(r,o)

$commodities:
	PA(r,g)$a0(r,g)		!	Regional market (input)
	PY(r,g)$s0(r,g)		!	Regional market (output)
	PD(r,g)$xd0(r,g)	!	Local market price
	PN(g)			!	National market
	PL(r,o,s)$ld0_o(r,o,s)			!	Wage rate
	PELL(r,o)$popr(r,o)
	PK(r,s)$kd0(r,s)	!	Rental rate of capital
 	PM(r,m)			!	Margin price
	PC(r,o)$popr(r,o)			!	Consumer price index
	PFX			!	Foreign exchange

$consumer:
	RA(r,o)$popr(r,o)			!	Representative agent
	GOVT

$auxiliary:
	TRANS

$prod:Y(r,s)$y_(r,s)  s:0 va:1
	o:PY(r,g)	q:(ys0(r,s,g)*(1-tfp(r,s)))            a:GOVT t:ty(r,s)    p:(1-ty0(r,s))
	i:PA(r,g)	q:(id0(r,g,s)*(1-aeeir(r,g)))	a:GOVT	t:te(r,g,s)
	i:PL(r,o,s)		q:ld0_o(r,o,s)	va:
	i:PK(r,s)	q:kd0(r,s)	va:	

$prod:X(r,g)$x_(r,g)  t:4
	o:PFX		q:(x0(r,g)-rx0(r,g))
	o:PN(g)		q:xn0(r,g)
	o:PD(r,g)	q:xd0(r,g)
	i:PY(r,g)	q:s0(r,g)

$prod:A(r,g)$a_(r,g)  s:0 dm:2  d(dm):4
 	o:PA(r,g)	q:a0(r,g)		a:GOVT	t:ta(r,g)	p:(1-ta0(r,g))
	o:PFX		q:rx0(r,g)
	i:PN(g)		q:nd0(r,g)	d:
	i:PD(r,g)	q:dd0(r,g)	d:
	i:PFX		q:m0(r,g)	dm: 	a:GOVT	t:tm(r,g) 	p:(1+tm0(r,g))
	i:PM(r,m)	q:md0(r,m,g)

$prod:MS(r,m)
	o:PM(r,m)	q:(sum(gm, md0(r,m,gm)))
	i:PN(gm)	q:nm0(r,gm,m)
	i:PD(r,gm)	q:dm0(r,gm,m)

$prod:C(r,o)$popr(r,o)  s:1
    o:PC(r,o)		q:(c0(r)*popr(r,o))
	i:PA(r,g)	q:(cd0(r,g)*popr(r,o))

$prod:LS(r,o)$popr(r,o)	t:1
	o:PL(r,o,s)		q:ld0_o(r,o,s)
	i:PELL(r,o)	q:(sum(s,ld0_o(r,o,s)))

$demand:RA(r,o)$popr(r,o)
	d:PC(r,o)		q:(c0(r)*popr(r,o))
	e:PY(r,g)	q:(yh0(r,g)*popr(r,o))
	e:PFX		q:((bopdef0(r) + hhadj(r))*popr(r,o))
	e:PA(r,g)	q:(-(i0(r,g))*popr(r,o))
	e:PELL(r,o)		q:((sum((s,o.local),ld0_o(r,o,s)))*popr(r,o))
	e:PK(r,s)	q:(kd0(r,s)*popr(r,o))
	e:PFX		q:tpo(r,o)	r:TRANS
	
$demand:GOVT
	d:PA(r,g)	q:g0(r,g)
	e:PFX		q:gov0
	e:PFX		q:(-sum((r,o),tpo(r,o)))	r:TRANS

$constraint:TRANS
	GOVT =e= sum((r,g),PA(r,g)*g0(r,g));

$report:
	v:LOS(r,o,s)	o:PL(r,o,s)	prod:LS(r,o)

$OFFTEXT
$SYSINCLUDE mpsgeset mgemodel -mt=0

TRANS.l = 1;
PFX.fx = 1;



execute_unload "ltran_out.gdx";