/* HSQMBC - Experiment for efficient measurement fo long-range H-X couplings
From MRC 2000, 38: p 265 - 273, Williamson, Marquez, Gerwich, and Kover */

/* CPMG added per Williamson's suggestion */
/* J Magn Reson 181 p 89-97 (2006)  Kover, Batta & Feher */
/* R Crouch 18 Sept 2006 */

#include <standard.h>
#include <chempack.h> 

extern int dps_flag;

static int	ph1[2] = {0,2},
		ph2[4] = {0,0,2,2},
		ph3[8] = {0,0,0,0,2,2,2,2},
		ph5[8] = {1,1,1,1,3,3,3,3},
		ph4[8] = {0,2,0,2,2,0,2,0};
		
pulsesequence()

{

double	pwxlvl    = getval("pwxlvl"),
	pwxlvlcp  = getval("pwxlvlcp"),
	tpwrcp    = getval("tpwrcp"),
	pwx       = getval("pwx"),
	pwcp      = getval("pwcp"),
	pwxcp     = getval("pwxcp"),
        tpwr_cf   = getval("tpwr_cf"),
	gt1       = getval("gt1"),
	pwxlvl_cf = getval("pwxlvl_cf"),
	gzlvlE    = getval("gzlvlE"),
	pwxlvl180 = getval("pwxlvl180"),
	pwx180    = getval("pwx180"),
	gtE       = getval("gtE"),
	EDratio   = getval("EDratio"),
	hsgt      = getval("hsgt"),
	hsglvl    = getval("hsglvl"),
	gstab     = getval("gstab"),
	jnxh      = getval("jnxh"),
	j1xh      = getval("j1xh"),
	phase     = getval("phase"),
	corfacth  = getval("corfacth"),
	corfactx  = getval("corfactx"),
	tau,taug,tauB;
	
int	iphase,
	icosel;
	
char	sspul[MAXSTR],pwx180ad[MAXSTR],cmd[MAXSTR],cpmg[MAXSTR],hshape[MAXSTR],xshape[MAXSTR];

tauB = 1/(2*j1xh);   /*not yet used */
tau  = 1/(2*jnxh);
taug = gtE + gstab + 2*GRADIENT_DELAY;

getstr("sspul",sspul);
getstr("pwx180ad",pwx180ad);
getstr("cpmg",cpmg);
getstr("hshape",hshape);
getstr("xshape",xshape);

iphase = (int) (phase + 0.5);
   icosel = -1;

/* This phase cycle exactly as for HSQMBC */
   settable(t1,2,ph1);
   settable(t2,4,ph2);
   settable(t3,8,ph3);
   settable(t5,8,ph5);
   settable(t4,8,ph4);
   
   /* Setup FAD and 2nd dimension */
   
   getelem(t1,ct,v2);
   getelem(t4,ct,oph);
   
   initval(2.0*(double)(((int)(d2*getval("sw1")+0.5)%2)),v14);
   
   if (iphase == 2)
     icosel = 1;

   add(v2,v14,v2);
   add(oph,v14,oph);

   status(A);
   decpower(pwxlvl);
   obspower(tpwr);
    
   rcvroff();
   if (sspul[0] == 'y')
   {
         zgradpulse(1.2*hsglvl,hsgt);
         rgpulse(pw,zero,rof1,rof1);
         zgradpulse(1.2*hsglvl,hsgt);
   }

   delay(d1);
   rgpulse(pw,zero,rof1,rof1);
   if (cpmg[0] == 'y')
   {
    obspower(tpwrcp);
    decpower(pwxlvlcp);
    obspwrf(4095.0*corfacth);
    decpwrf(4095.0*corfactx);
    obsprgon(hshape,pwcp,1.0);
    decprgon(xshape,pwxcp,1.0);
    decunblank();
    obsunblank();
    decon();
    xmtron();
    delay(tau);
    decoff();
    xmtroff();
    obsblank();
    decblank();
    obsprgoff();
    decprgoff();
    obspwrf(4095.0);
    decpwrf(4095.0);
   }
  else
   {
    delay(tau);
   }
   status(B);
 
     obspower(tpwr);
     obspwrf(4095.0);
     decpower(pwxlvl);
     decpwrf(4095.0);
     
     rgpulse(pw,one,rof1,rof1);

     zgradpulse(hsglvl,gt1);
     delay(gstab);
     decrgpulse(pwx,v2,0.0,0.0);
     
     delay(d2/2);
     rgpulse(2*pw,t2,rof1,rof1);
     delay(d2/2);
 
     zgradpulse(gzlvlE,gtE);
     delay(gstab);
     if (pwx180 != 0.0)
     {
       decpower(pwxlvl180);
       decshaped_pulse(pwx180ad,pwx180,t5,rof1,rof1);
       decpower(pwxlvl);
     }
     else
     {
       decrgpulse(pwx,t5,3e-6,3e-6);
       decrgpulse(2.0*pwx,t3,0.0,0.0);
       decrgpulse(pwx,t5,3e-6,3e-6);
     }
     delay(taug);
     decrgpulse(pwx,t3,3e-6,3e-6);
     zgradpulse(0.8*hsglvl,gt1);
     delay(gstab);
     rgpulse(pw,zero,3e-6,3e-6);
     delay(taug);
     rgpulse(2.0*pw,zero,3e-6,3e-6);
     zgradpulse(gzlvlE*icosel/EDratio,gtE);
     delay(gstab);
     decpower(dpwr);
     rcvron();
   status(C);
}  
     
