/*LR_HSQMBC Version 1*/
/*Williamson, Buevich, Martin, and Parella JOC 79 p3887-3894 (2014)*/
/* Version 1 R Crouch 21 May 2014 */


#include <standard.h>
#include <chempack.h>

static int	ph1[1] = {0},
		ph2[1] = {1},
		ph3[2] = {0,2},
		ph4[8] = {0,0,0,0,2,2,2,2},
		ph5[4] = {0,0,2,2},
		ph6[1] = {0},
		ph7[1] = {0},
		ph8[1] = {0},
		ph9[1] = {2},
		ph10[8] = {0,2,0,2,2,0,2,0};
		

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
        j1xh = getval("j1xh"),
        phase = getval("phase"),
        pwx180 = getval("pwx180"),
	pwxlvl180 = getval("pwxlvl180"),
	pwx180r = getval("pwx180r"),
	pwxlvl180r = getval("pwxlvl180r"),
        tpwr180 = getval("tpwr180"),
        pw180 = getval("pw180"),
        grad1, grad2,
	tau, taug, tauB,
	delta, delta1, delta2, delta3;
	
	
	
int	iphase,
	icosel;
	
char	sspul[MAXSTR],
        shaped[MAXSTR],
        pwx180adR[MAXSTR],
        pwx180ad[MAXSTR];
        
  
   
   getstr("sspul",sspul);
   getstr("shaped",shaped);
   getstr("pwx180ad", pwx180ad);
   getstr("pwx180adR", pwx180adR);


tauB = 1/(2*j1xh);
tau  = 1/(4*jnxh);

delta = gtE + gstab + (2*pw) + 4*rof1 + 4e-6;
delta1=tau;
delta2=tau;
delta3=tau;


  if (pw>pwx)
  {
   delta1 = tau - pw - rof1;;
   delta2 = tau - pw - (2*pw/3.141);
   delta3 = tau - gtE - pw - 4*rof1;
  }

  if (pwx>pw)
  {
   delta1 = tau - pwx - rof1;
   delta2 = tau - pwx - (2*pw/3.141);
   delta3 = tau - gtE - pwx - 4*rof1;
  }

   iphase = (int) (phase + 0.5);
   icosel = -1;

   settable(t3,2,ph3);
   settable(t7,8,ph10);
   settable(t4,8,ph4);
   settable(t5,4,ph5);
   
   /* Setup FAD and 2nd dimension */
   
   getelem(t3,ct,v3);
   getelem(t7,ct,oph);
   
   initval(2.0*(double)(((int)(d2*getval("sw1")+0.5)%2)),v14);
   
    if (iphase == 2)
     icosel = 1;

   add(v3,v14,v3);
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

      delay(d1);
status(B);
      rgpulse(pw,zero,rof1,rof1);
      
      delay(delta1);
      zgradpulse(gzlvl1,gt1);
      delay(gstab);
      pulse(pw,zero);
      
      if (shaped[0] == 'y')
       { 
        delay(tauB - pwx180/2);
        decpower(pwxlvl180);
        simshaped_pulse("",pwx180ad,pw*2,pwx180,zero,zero,2e-6,2e-6); 
        delay(tauB - pwx180/2);
        decpower(pwxlvl);
       }
      else
       { 
        delay(tauB);
        simpulse(2.0*pw,2.0*pwx,zero,zero,2e-6,2e-6);
        delay(tauB);
       }
      pulse(pw,two);
      delay(gstab);
      zgradpulse(gzlvl1,gt1);
      delay(delta1);
     
      rgpulse(pw,one,rof1,2e-6);
      zgradpulse(gzlvl1*1.3,gt1);
      delay(gstab);
      decrgpulse(pwx,v3,rof1,rof1);
      
      delay(d2/2.0);
      rgpulse(2.0*pw,t5,rof1,rof1);
      delay(d2/2.0);
      
      zgradpulse(grad1*icosel,gtE);
      delay(gstab);
       if (shaped[0] == 'y')
         /* { 
         decpower(pwxlvl180);
         decshaped_pulse(pwx180ad,pwx180,t4,rof1,rof1);
         delay(delta - pwx180*2 - 2*POWER_DELAY);
         decshaped_pulse(pwx180ad,pwx180,t4,rof1,rof1);
         decpower(pwxlvl);
        }  */
        { 
         decrgpulse(2.0*pwx,t4,rof1,rof1); 
         delay(delta - 2.5*pwx);
        }  
       else
        { 
         decrgpulse(2.0*pwx,t4,rof1,rof1); 
         delay(delta - 2.5*pwx);
        }

       simpulse(pw,pwx,zero,t4,2e-6,2e-6); 
      
       if (shaped[0] == 'y')
        { 
          decpower(pwxlvl180);
          delay(delta2 - pwx180/2 + pw - POWER_DELAY);
          simshaped_pulse("",pwx180ad,pw*2,pwx180,zero,zero,rof1,rof1); 
          zgradpulse(grad2,gtE);
          decpower(pwxlvl);
          delay(delta3 - pwx180/2 + pw - POWER_DELAY);
        }
       else
        { 
         delay(delta2);
         simpulse(2.0*pw,2.0*pwx,zero,zero,rof1,rof1);
         zgradpulse(grad2,gtE);
         delay(delta3);
        }
      decrgpulse(pwx,zero,rof1,rof2); 
      decpower(dpwr);
status(C);
}
