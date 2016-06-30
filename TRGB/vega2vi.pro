PRO vega2vi,f606,f814,vmag,imag
;takes VEGA F606 and F814 magnitudes and converts to V&I

;RDFLOAT,'ngc0055-disk/ngc0055-disk_psfphot.dat',x,dx,y,dy,f606,f606e,f814,f814e
nstar=N_ELEMENTS(f814)

f606=DOUBLE(f606)
f814=DOUBLE(f814)

;have to remove zeropoints to use Sirianni's transformation
f606zeropoint=26.398
f814zeropoint=25.501
f606=f606-f606zeropoint
f814=f814-f814zeropoint

;transformations from Sirianni using BPGC stellar atlas
vblue_c0=26.394d0 & vblue_c0e=0.005
vblue_c1=0.153d0 & vblue_c1e=0.018
vblue_c2=0.096d0 & vblue_c2e=0.085
vred_c0=26.331d0 & vred_c0e=0.008
vred_c1=0.340d0 & vred_c1e=0.008
vred_c2=-0.038d0 & vred_c2e=0.002

iblue_c0=25.489d0 & iblue_c0e=0.013
iblue_c1=0.041d0 & iblue_c1e=0.211
iblue_c2=-0.093d0 & iblue_c2e=0.803
ired_c0=25.496d0 & ired_c0e=0.010
ired_c1=-0.014d0 & ired_c1e=0.013
ired_c2=0.015d0 & ired_c2e=0.003


vminusi=f606-f814
vold=f606
iold=f814
test=0
diff=1.
miniter=10 ;minimum number of iterations before quiting
vmag=DBLARR(nstar)
imag=DBLARR(nstar)
WHILE (diff GT 0.0001) DO BEGIN
    FOR i=0l,nstar-1 DO BEGIN
        print, diff
        ;first do V transformation
        CASE 1 OF
            (vminusi[i] GE 0.5): vmag[i]=f606[i] + vred_c0 + (vminusi[i]*vred_c1) + (vminusi[i]^2*vred_c2)
            (vminusi[i] LE 0.3): vmag[i]=f606[i] + vblue_c0 + (vminusi[i]*vblue_c1) + (vminusi[i]^2*vblue_c2)
            ELSE: BEGIN
                ;smooth ~0.01 mag jump in V transformation eq.
                redfrac=(vminusi[i]-0.3)*5.
                bluefrac=1.-redfrac
                vmag[i]=bluefrac*(f606[i] + vblue_c0 + (vminusi[i]*vblue_c1) + (vminusi[i]^2*vblue_c2))+redfrac*(f606[i] + vred_c0 + (vminusi[i]*vred_c1) + (vminusi[i]^2*vred_c2))
            END
        ENDCASE
        ;Now do I transformation
        CASE 1 OF
            (vminusi[i] GE 0.2): imag[i]=f814[i] + ired_c0 + (vminusi[i]*ired_c1) + (vminusi[i]^2*ired_c2)
            (vminusi[i] LE 0.0): imag[i]=f814[i] + iblue_c0 + (vminusi[i]*iblue_c1) + (vminusi[i]^2*iblue_c2)
            ELSE: BEGIN
                ;smooth ~0.003 mag jump in I transformation eq.
                redfrac=vminusi[i]*5.
                bluefrac=1.-redfrac
                imag[i]=bluefrac*(f814[i] + iblue_c0 + (vminusi[i]*iblue_c1) + (vminusi[i]^2*iblue_c2))+redfrac*(f814[i] + ired_c0 + (vminusi[i]*ired_c1) + (vminusi[i]^2*ired_c2))
            END
        ENDCASE
    ENDFOR
    
    ;now compare to previous and set up for next iteration
    diff=MAX([ABS(vmag-vold),ABS(imag-iold)],maxpos)
    ind=WHERE(ABS(vmag-vold) GT 0.0001)
;    plothist,vminusi[ind],bin=0.001
;    plothist,ABS(vmag-vold),bin=0.0001,/ylog,yrange=[0.1,1e6],xrange=[-0.002,0.010]
;    print, diff
    vminusi=vmag-imag & vold=vmag & iold=imag
    test=test+1
    IF (test LT miniter) THEN diff=1.
ENDWHILE

END
