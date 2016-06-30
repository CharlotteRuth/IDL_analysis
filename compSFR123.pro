pro compSFR123
;This graphs the (star fromation rate)/(gass mass) over time for galaxies with
;different resolution.  One graph per galactic mass is produced.

r=['5E1R', '1E2R', '1E3R','1E4R','1E5R']
m=['9M', '10M', '11M','12M','13M']
trials = ['1','2','3','4','5','6','7','8','9']
integrate=fltarr(5)
binsize=1.e8
timeunit=1e9
massunit=2.325e5
;set_plot,'ps'
set_plot,'x'

for mct=0,4 do begin
    first = 0
    start = 0
    ymax=0   
;    device,filename=m[mct]+'sfr.eps'
    for rct=0,4 do begin
        ; Allows for different number of trials at different resolutions
        if (rct lt 3) then begin trial = 9 
        endif else if (rct lt 5) then begin  trial = 3
        endif else if (rct lt 6) then begin trial = 1
        endif else begin trial = 0
        endelse
        ; Moves to the correct directory
        cd,r[rct]    
        cd,m[mct]
        ;Initializes the star structure
        s={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,METALS:0,TFORM:0,EPS:0,PHI:0}
        for t = 0, trial-1 do begin
            ssub={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,METALS:0,TFORM:0,EPS:0,PHI:0}
            file='o'+m[mct]+'__'+trials[t]+'.00300'
            rtipsy,file,h,g,d,ssub
            ;Makes an array of all stars formed in all trials
            if (t ne 0) then concat_structs,s,ssub,ssum else ssum = ssub 
            s = ssum
        ENDFOR
        if (s[0].mass ne 0)then begin
            ;If some stars actually formed . . .
            ind=WHERE(s.tform GT 0.0)
            ;For all stars where tform is greater than 0, 
            ;tform=tform*timeunit and mass = mass*massunit
            tform=s[ind].tform*timeunit
            mass=s[ind].mass*massunit
            ;Find min and max star formation times and create bins
            mintime=MIN(tform)
            maxtime=MAX(tform)
            nbins=FIX((maxtime-mintime)/binsize)+1
            timearray=FINDGEN(nbins)*binsize+mintime
            sfr=FLTARR(nbins)
            FOR i=0,nbins-1 DO BEGIN
                                ;For each bin, figure out the mass of
                                ;stars were formed in that time bin
                                ;and make into a SFR 
                inbin=WHERE(tform GE timearray[i] AND tform LT (timearray[i]+binsize), ninbin)
                IF (ninbin EQ 0) THEN sfr[i]=0 ELSE sfr[i]=TOTAL(mass[inbin])/binsize/10.^(m+9)*10.
                                ;Check to see that m is really
                                ;measuring the baryonic matter mass,
                                ;not dark matter
            ENDFOR
  ;          print,sfr
            ymaxtemp=MAX(sfr)
 ;           print,ymaxtemp
            if (ymaxtemp GT ymax)then ymax=ymaxtemp
            ; Finds the maximum SFR to be used for the y-axes
            print,ymax
            integrate[rct]=TOTAL(mass)
            ;Total mass of stars formed
            print,10.^m[mct]
        endif else integrate[rct]=0.0000
        cd,'../..'
    ENDFOR

    ;Scale the y-axes
    ymax=1.01*ymax

    for rct=0,4 do begin
        ; Allows for different number of trials at different resolutions
        if (rct lt 3) then begin trial = 9 
        endif else if (rct lt 5) then begin  trial = 3
        endif else if (rct lt 6) then begin trial = 1
        endif else begin trial = 0
        endelse
        ; Moves to the correct directory
        cd,r[rct]    
        cd,m[mct]
        ;Initializes the star structure
        s={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,METALS:0,TFORM:0,EPS:0,PHI:0}
        for t = 0, trial-1 do begin
            ssub={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,METALS:0,TFORM:0,EPS:0,PHI:0}
            file='o'+m[mct]+'__'+trials[t]+'.00300'
            rtipsy,file,h,g,d,ssub
            ;makes an array of all a stars formed
            if (t ne 0) then concat_structs,s,ssub,ssum else ssum = ssub 
            s = ssum
        ENDFOR
        if (s[0].mass ne 0)then begin
  ;          print,sfr
   ;         print,s
            IF first EQ 0 then begin
                ; If this is the first line plotted, don't overplot
                first = 1
                                ;call the function sfr.pro which will
                                ;graph the star formation rate
                sfr,s,massu=2.325e5,time=1e9,LINESTYLE=(5-first-start),OVERPLOT=0,xrange=[0,3],binsize=binsize,ymax=ymax,gmass=m[mct]
            endif else begin
                ; If it is not the first line plotted, overplot
                first = first+1
                sfr,s,massu=2.325e5,time=1e9,LINESTYLE=(5-first-start),OVERPLOT=1,binsize=binsize,ymax=ymax,gmass=m[mct]            
            endelse
        endif else start = start+1
        cd,'../..'
    endfor
    legend,['1E5R '+trim(integrate[4]),'1E4R '+trim(integrate[3]),'1E3R '+trim(integrate[2]),'1E2E '+trim(integrate[1]),'2E5R '+trim(integrate[0])],linestyle=[0,1,2,3,4],/right
    ;device,/close
    stop
endfor
END
