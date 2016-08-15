/*PIP_HSQMBC Version 1*/
/*Williamson,  and Parella */
/* Version 1.3 R Crouch 29 July 2014 */
/* simple presat and wet supported but not purge presat */


#include <standard.h>
#include <chempack.h>

static int	ph1[1] = {0},
		ph2[1] = {1},
		ph3[2] = {0,2},
		ph4[1] = {0},
		ph5[1] = {0},
		ph6[1] = {0},
		ph10[2] = {0,2};
		

pulsesequence()

{

double	pwxlvl = getval("pwxlvl"),
	pwx = getval("pwx"),
	trim = getval("trim"),
	gzlvl1 = getval("gzlvl1"),
	gt1 = getval("gt1"),
	gzlvlE = getval("gzlvlE"),
	gtE = getval("gtE"),
	EDratio = getval("EDratio"),
	hsgt = getval("hsgt"),
	hsglvl = getval("hsglvl"),
        gstab = getval("gstab"),
        jnxh = getval("jnxh"),
        phase = getval("phase"),
        pwx180 = getval("pwx180"),
	pwxlvl180 = getval("pwxlvl180"),
        zqfpw1 = getval("zqfpw1"),
        zqfpwr1 = getval("zqfpwr1"),
        gzlvlzq1 = getval("gzlvlzq1"),
	evolcorr,
        grad1, grad2,
	tau, taug,
	delta, delta2, delta4;
	
int	iphase, ipap,
	icosel;
	
char	sspul[MAXSTR],
        bipflg[MAXSTR],
        zqfpat1[MAXSTR],
        pwx180ad[MAXSTR];       
  
   ipap = getval("ipap");
   getstr("sspul",sspul);
   getstr("bipflg",bipflg);
   getstr("zqfpat1",zqfpat1);
   getstr("pwx180ad", pwx180ad);

   pwxlvl180=getval("dnbippwr");
   pwx180=getval("dnbippw");
   getstr("dnbipshp",pwx180ad);

   tau  = 1/(4*jnxh);
   delta = gtE + gstab + (2*pw);
   delta2=tau - pwx180/2;
   delta4=gtE + gstab;


   iphase = (int) (phase + 0.5);
   icosel = 1; 
   
   settable(t1,1,ph1);
   settable(t2,1,ph2);
   settable(t3,2,ph3);
   settable(t7,2,ph10);
   settable(t4,1,ph4);
   settable(t5,1,ph5);
   settable(t6,1,ph6);

   /* Setup FAD and 2nd dimension */
   
   getelem(t3,ct,v3);
   getelem(t6,ct,v6);
   getelem(t7,ct,oph);
   
   initval(2.0*(double)(((int)(d2*getval("sw1")+0.5)%2)),v14);
   
    if (iphase == 2)
     icosel = -1;

   add(v3,v14,v3);
   add(v6,v14,v6);
   add(oph,v14,oph);
  
   grad1 = gzlvlE;
   grad2 =  gzlvlE/EDratio;


status(A);
     decpower(pwxlvl);
     obspower(tpwr);
     if (sspul[0] == 'y')
   {
         zgradpulse(hsglvl,hsgt);
         rgpulse(pw,zero,rof1,rof1);
         zgradpulse(hsglvl,hsgt);
   }

      if (satmode[0] == 'y')
     {
        if ((d1-satdly) > 0.02)
                delay(d1-satdly);
        else
                delay(0.02);
        if (getflag("slpsat"))
           {
                shaped_satpulse("relaxD",satdly,zero);
           }
        else
           {
                satpulse(satdly,zero,rof1,rof1);
           }
     }
   else
        delay(d1);

   if (getflag("wet"))
     wet4(zero,one);

status(B);
      rgpulse(pw,t1,rof1,rof1);
 
       delay(delta2);
       decpower(pwxlvl180);
       simshaped_pulse("",pwx180ad,pw*2,pwx180,t1,v6,2e-6,2e-6); 
       decpower(pwxlvl);    
       delay(delta2);
 
      rgpulse(pw,t2,2e-6,2e-6);
      zgradpulse(gzlvl1,gt1);
      delay(gstab);
      decpulse(pwx,v3);

      delay(d2/2.0);
      rgpulse(2.0*pw,t5,rof1,rof1);
      delay(d2/2.0);
       
         decpower(pwxlvl180);
         decshaped_pulse(pwx180ad,pwx180,t4,rof1,rof1);
// Correct here for every detail more or less
//	 evolcorr=POWER_DELAY + pwx180/PI + 2*rof1;
         delay(POWER_DELAY + 3*rof1 + 2*(pwx180/3.1416));
         zgradpulse(grad1*icosel,gtE);
         delay(gstab);
         decshaped_pulse(pwx180ad,pwx180,t4,rof1,rof1);
         decpower(pwxlvl);
      
      delay(delta);
      pulse(pw,t1);
      decpulse(pwx,t4);

      if (ipap == 0)
       {
        
        delay(delta2);
        
          decpower(pwxlvl180);
          simshaped_pulse("",pwx180ad,pw*2,pwx180,t1,t1,2e-6,2e-6); 
          decpower(pwxlvl);
       
        delay(delta2);
       }
      else
       {
        delay(tau);
        pulse(2.0*pw,t2);
        delay(tau);
       }

      if (ipap == 0)
       {
        pulse(pw,t2);
       }
      else
       {
        pulse(pw,t1);
       }
      if (getflag("Gzqfilt"))
        {
         delay(gstab);
         obspower(zqfpwr1);
         rgradient('z',gzlvlzq1);
         delay(100.0e-6);
         shaped_pulse(zqfpat1,zqfpw1,t1,rof1,rof1);
         delay(100.0e-6);
         rgradient('z',0.0);
         delay(gstab);
	 obspower(tpwr);
        }
      pulse(pw,t1);
      delay(delta4);
      pulse(2.0*pw,t1);
      delay(delta4);
      pulse(pw,t2);
      delay(delta4);
      pulse(2.0*pw,t1);
      zgradpulse(grad2,gtE);
      delay(gstab);
      decpower(dpwr);
status(C);
}
