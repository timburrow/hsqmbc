//
//  MEselHSQMBC.c
//  Multiplicity-edited selHSQMBC with options for TOCSY transfer
//
//  Ref: Gary Martin Journal of Natural Products 2015 78 (9), 2236-2241
//
//  VnmrJ version Created by Timothy Burrow on 2016-08-12.
//

#include <standard.h>
#include <chempack.h>

/* Bruker phases
 ph1=0
 ph2=1
 ph3=0 2
 ph4=0
 ph5=0
 ph6=0
 ph7=0
 ph8=0
 ph23=3 : ph9
 ph25=1 : ph10
 ph30=0 : ph11
 ph31=0 2 : ph12 (oph)
*/

/* Bruker pulse to Agilent 
 p1 : pw
 p2 : pw*2
 p3 : pwx
 p4 : pwx*2
 p6 : (TOCSY)
 p12 : selpwA
 p14 : 13C 180deg inversion (500 µs)
 p16 : hsgt homospoil gradient
 p31 : 13C 180deg adiabatic matched sweep (1730 µs @ 600 MHz)
 p32 : 1H 180deg adiabatic (30 ms)
 
 */

/* Bruker delay to Agilent
 d0 : 2d delay (NA)
 d1 : d1 relaxation delay
 d4 : taumb 1/(4*JnCH)
 d9 : mixT (TOCSY mix time)
 d11 : NA (diskIO)
 d16 : gstab
 d21 : tau 1/(2*J1CH)
 
 */
 

static int	ph1[1] = {0},
            ph2[1] = {1},
            ph3[2] = {0,2},
            ph4[1] = {0},
            ph5[1] = {0},
            ph6[1] = {0},
            ph7[1] = {0},
            ph8[1] = {0},
            ph9[1] = {3},       // ph23
            ph10[1] = {1},      // ph25
            ph11[1] = {0},      // ph30
            ph12[2] = {0,2};    // ph31 (oph)

pulsesequence()
{
    
    double  jnxh = getval("jnxh"),
            j1xh = getval("j1xh"),
            selpwA = getval("selpwA"), // length of selective H pulse
            selpwrA = getval("selpwrA"), // power of selective H pulse
            pwx180 = getval("pwx180"),
            pwxlvl180 = getval("pwxlvl180"),
            mult = getval("mult"),
            gzlvl0 = getval("gzlvl0"),
            gt0 = getval("gt0"),
            gzlvlE = getval("gzlvlE"),
            gtE = getval("gtE"),
            EDratio = getval("EDratio"),
            gstab = getval("gstab"),
            hsgt = getval("hsgt"),
            hsglvl = getval("hsglvl"),
            phase = getval("phase"),
            mixT = getval("mixT"),
            gzlvl1 = getval("gzlvl1"), //?
            gt1 = getval("gt1"),
            gzlvl3 = getval("gzlvl3"),
            gt3 = getval("gt3"),
            tau, taumb,
            grad1,
            grad2,
            grad3,
            grad4
            grad5
            delta, delta1, delta2, delta4, delta5,
            factor1;
    
    char    tocsyflg[MAXSTR],
            slpatT[MAXSTR],
            pwx180ad[MAXSTR],
            selshapeA[MAXSTR];

    int     iphase,
            icosel; // gradient switches sign
    
    getstr("tocsyflg", tocsyflg);
    getstr("pwx180ad", pwx180ad);
    getstr("selshapeA", selshapeA);

    taumb  = 1/(4*jnxh);
    tau = 1/(2*j1xh);

    // DELTA=p16+d16+50u+p2+d0*2
    // DELTA1=50u+p16+d16
    // DELTA2=d4-larger(p12,p14)/2-50u-p16-d16
    // DELTA4=d21-p2+p3*2/PI
    // DELTA5=d21-p16-d16-p2-d0*2-p3*2/PI
    // FACTOR1=(d9/(p6*115.112))/2+0.5
    
    delta = gstab + gtE + pwx;
    delta1 = pw + gtE;
    delta2 = taumb - pwx180/2;
    delta4 = taumb - fmax(pwx180,)
    iphase = (int) (phase + 0.5);
    icosel = 1;
    
    settable(t1,1,ph1);
    settable(t2,1,ph2);
    settable(t3,2,ph3);
    settable(t4,1,ph4);
    settable(t5,1,ph5);
    settable(t6,1,ph6);
    settable(t7,1,ph7);
    settable(t8,1,ph8);
    settable(t9,1,ph9);
    settable(t10,1,ph10);
    settable(t11,1,ph11);
    settable(t12,2,ph12);

    /* Setup axial displacement and 2nd dimension */
    
    getelem(t3,ct,v3);
    getelem(t6,ct,v6);
    getelem(t12,ct,oph);
    
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

    delay(5.0e-5);
    if (getflag("sspul"))
        steadystate();
    
    if (satmode[0] == 'y')
    {
        if ((d1-satdly) > 0.02)
            delay(d1-satdly);
        else
            delay(0.02);
        if (getflag("slpsat"))
        {
            shaped_satpulse("relaxD",satdly,zero);
            if (getflag("prgflg"))
                shaped_purge(v6,zero,v18,v19);
        }
        else
        {
            satpulse(satdly,zero,rof1,rof1);
            if (getflag("prgflg"))
                purge(v6,zero,v18,v19);
        }
    }
    else
        delay(d1);
    
    if (getflag("wet"))
        wet4(zero,one);

status(B);
    rgpulse(pw,t1,rof1,rof1);
    zgradpulse(gzlvl1,gt1);
    
    obspower(selpwrA);
    decpower(pwxlvl180);
    
    delay(delta2);
    simshaped_pulse("selshapeA","pwx180ad",selpwA,pwx180,t1,v6,2e-6,2e-6);
    decpower(pwxlvl);
    obspower(tpwr);

    delay(delta2);
    
    zgradpulse(gzlvl1,gt1);
    delay(gstab);
    
    pulse(pw,t2);
    decpulse(pwx,t3);
    
    delay(d2/2.0);
    pulse(2.0*pw,t5);
    delay(d2/2.0);
    
    zgradpulse(grad1*icosel,gtE);
    decpower(pwxlvl180);
    delay(delta5);
    simshaped_pulse("","pwx180ad",pw*2,pwx180,t1,v6,2e-6,2e-6);
    delay(delta5);
    
    decshaped_pulse("pwx180ad",pwx180,t1,rof1,rof1);
    
    decpower(pwxlvl);
    obspower(tpwr);
    
    simpulse(pw,pwx,t2,t4,rof1,0.0);
  
    obspower(selpwrA);
    decpower(pwxlvl180);
    zgradpulse(gzlvl1,gt1); // fix
    delay(gstab);
    
    delay(delta2);
    simshaped_pulse("selshapeA","pwx180ad",selpwA,pwx180,t1,v6,2e-6,2e-6);
    delay(delta2);
    
    obspower(tpwr);
    
    zgradpulse(gzlvl1,gt1); // fix
    delay(gstab);
    zgradpulse(gzlvl2,gt2); // fix
    delay(gstab);

    pulse(pw,t1);
    
    decpower(dpwr);
status(C);
    
}


