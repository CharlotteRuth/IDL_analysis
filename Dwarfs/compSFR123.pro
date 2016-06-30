pro compSFR123
;This graphs the (star fromation rate)/(gass mass) over time for galaxies with
;different resolution.  One graph per galactic mass is produced.

r=['5E1R', '1E2R', '1E3R','1E4R','1E5R','1E6R_Rok']
resN = [50,100,1000,10000,100000]
nRes = 5
m=['9M', '10M', '11M','12M','13M']
mN = [1.e9,1.e10,1.e11,1.e12,1.e13]
nMass = 4
trials = ['1','2','3','4','5','6','7','8','9']
integrate=fltarr(nRes)
aveSFR = fltarr(nRes,nMass)
binsize=1.e8
timeunit=1e9
massunit=2.325e5
;set_plot,'ps'
set_plot,'x'

for mct=0,nMass-1 do begin ;Finds ymax
    first = 0
    start = 0
    ymax=0   
;    device,filename=m[mct]+'sfrMulti.eps'
    for rct=0,nRes-1 do begin
        if (rct lt 3) then begin trial = 9. ; Allows for different number of trials at different resolutions
        endif else if (rct lt 5) then begin  trial = 3.
        endif else if (rct lt 6) then begin trial = 1.
        endif else begin trial = 0
        endelse
        cd,r[rct]+'/'+m[mct]; Moves to the correct directory
        s={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,METALS:0,TFORM:0,EPS:0,PHI:0};Initializes the star structure
        for t = 0, trial-1 do begin
            ssub={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,METALS:0,TFORM:0,EPS:0,PHI:0}
            file='o'+m[mct]+'__'+trials[t]+'.00300'
            rtipsy,file,h,g,d,ssub
            if (t ne 0) then concat_structs,s,ssub,ssum else ssum = ssub ;Makes an array of all stars formed in all trials
            s = ssum
        ENDFOR
        if (s[0].mass ne 0)then begin ;If some stars actually formed . . .
            ind=WHERE(s.tform GT 0.0);For all stars where tform is greater than 0,
            tform=s[ind].tform*timeunit ;tform=tform*timeunit and mass = mass*massunit
            mass=s[ind].mass*massunit
            mintime=MIN(tform) ;Find min and max star formation times and create bins
            maxtime=MAX(tform)
            nbins=FIX((maxtime-mintime)/binsize)+1
            timearray=FINDGEN(nbins)*binsize+mintime
            sfr=FLTARR(nbins)
            FOR i=0,nbins-1 DO BEGIN ;For each bin, figure out the mass of For each bin, figure out the mass of stars were formed in that time bin and make into a SFR 
                inbin=WHERE(tform GE timearray[i] AND tform LT (timearray[i]+binsize), ninbin)
                IF (ninbin EQ 0) THEN sfr[i]=0 ELSE sfr[i]=TOTAL(mass[inbin])/binsize/10.^(mct+9)*10./trial ;Check to see that m is really measuring the baryonic matter mass, not dark matter
            ENDFOR
            ymaxtemp=MAX(sfr)
 ;           print,ymaxtemp
            if (ymaxtemp GT ymax)then ymax=ymaxtemp ; Finds the maximum SFR to be used for the y-axes
            integrate[rct]=TOTAL(mass) ;Total mass of stars formed
        endif else integrate[rct]=0.0000
        cd,'../..'
    ENDFOR

    ymax=1.01*ymax ;Scale the y-axes

    for rct=0,nRes-1 do begin
        if (rct lt 3) then begin trial = 9 ; Allows for different number of trials at different resolutions
        endif else if (rct lt 5) then begin  trial = 3
        endif else if (rct lt 6) then begin trial = 1
        endif else begin trial = 0
        endelse
        cd,r[rct]+'/'+m[mct]; Moves to the correct directory
        s={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,METALS:0,TFORM:0,EPS:0,PHI:0};Initializes the star structure
        for t = 0, trial-1 do begin
            ssub={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,METALS:0,TFORM:0,EPS:0,PHI:0}
            file='o'+m[mct]+'__'+trials[t]+'.00300'
            rtipsy,file,h,g,d,ssub
              if (t ne 0) then concat_structs,s,ssub,ssum else ssum = ssub ;makes an array of all  stars formed
            s = ssum
        ENDFOR
        if (s[0].mass ne 0)then begin
            lines = nRES-first-start
            gm = mct+9
            IF first EQ 0 then begin
                first = 1 ; If this is the first line plotted, don't overplot
;call the function sfr.pro which will graph the star formation rate
                sfr,s,a,massu=2.325e5,time=1e9,LINESTYLE=lines,OVERPLOT=0,xrange=[0,3],binsize=binsize,gmass=gm,ntrial = trial,YMAX = ymax
            endif else begin
                first = first+1 ; If it is not the first line plotted, overplot
                sfr,s,a,massu=2.325e5,time=1e9,LINESTYLE=lines,OVERPLOT=1,binsize=binsize,gmass= gm, ntrial = trial        
            endelse
            aveSFR[rct,mct]=a/mn[mct]
        endif else start = start+1
        cd,'../..'
    endfor
    legend,['1E5R '+trim(integrate[4]),'1E4R '+trim(integrate[3]),'1E3R '+trim(integrate[2]),'1E2E '+trim(integrate[1]),'2E5R '+trim(integrate[0])],linestyle=[0,1,2,3,4],/right
;    device,/close
endfor

set_plot,'ps'
device,filename='MeanSFR_vs_Res.eps'

FOR mct=0,nMass-1 do begin
    print,aveSFR[*,mct]
    IF mct EQ 0 THEN plot,resN,aveSFR[*,mct],xtitle = 'Number of Particles',YTITLE = 'Average SFR over 3 Gigayears',TITLE = 'SFR/Mass of Galaxy as a Function of Resolution',psym= mct*(-1) - 1,yrange=[0,max(aveSFR)],xrange=[50,100000],/xlog ELSE oplot,resN,aveSFR[*,mct],psym= mct*(-1) - 1
ENDFOR
legend,['M: 1e9 '+trim(integrate[4]),'M: 1e10 '+trim(integrate[3]),'M: 1e11 '+trim(integrate[2]),'M: 1e12 '+trim(integrate[1]),'M: 1e13 '+trim(integrate[0])],psym=[-1,-2,-3,-4,-5],/right
device,/close
stop
END
