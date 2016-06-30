pro compSFR
;This graphs the (star fromation rate)/(gass mass) over time for galaxies with
;different resolution.  One graph per galactic mass is produced.

r=['5E1R', '1E2R', '1E3R', '1E4R', '1E5R']
m=['7M', '8M', '9M', '10M', '11M','12M', '13M']
integrate=fltarr(5)
binsize=1.e8
timeunit=1e9
massunit=2.325e5
set_plot,'ps'
;set_plot,'x'

for mct=0,6 do begin
    first = 0
    start = 0
    ymax=0   
    device,filename=m[mct]+'sfr.eps'
    for rct=0,4 do begin
        s={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,METALS:0,TFORM:0,EPS:0,PHI:0}
        cd,r[rct]    
        cd,m[mct]
        file='o'+m[mct]+'.00300'
        rtipsy,file,h,g,d,s
        if (s[0].mass ne 0)then begin
            ind=WHERE(s.tform GT 0.0)
            tform=s[ind].tform*timeunit
            mass=s[ind].mass*massunit
            mintime=MIN(tform)
            maxtime=MAX(tform)
            nbins=FIX((maxtime-mintime)/binsize)+1
            timearray=FINDGEN(nbins)*binsize+mintime
            sfr=FLTARR(nbins)
            FOR i=0,nbins-1 DO BEGIN
                inbin=WHERE(tform GE timearray[i] AND tform LT (timearray[i]+binsize),ninbin)
                IF (ninbin EQ 0) THEN sfr[i]=0 ELSE sfr[i]=TOTAL(mass[inbin])/binsize/10.^m[mct]*10
            ENDFOR
            ymaxtemp=MAX(sfr)
            if (ymaxtemp GT ymax)then ymax=ymaxtemp
            integrate[rct]=TOTAL(mass)
            print,10.^m[mct]
        endif else integrate[rct]=0.0000
        cd,'../..'
    ENDFOR

    ymax=1.01*ymax

    for rct=0,4 do begin
        s={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,METALS:0,TFORM:0,EPS:0,PHI:0}
        cd,r[rct]    
        cd,m[mct]
        file='o'+m[mct]+'.00300'
        rtipsy,file,h,g,d,s
       if (s[0].mass ne 0)then begin
            IF first EQ 0 then begin
                first = 1
                sfr,s,massu=2.325e5,time=1e9,LINESTYLE=(5-first-start),OVERPLOT=0,xrange=[0,3],binsize=binsize,ymax=ymax,gmass=m[mct]
            endif else begin
                first = first+1
                sfr,s,massu=2.325e5,time=1e9,LINESTYLE=(5-first-start),OVERPLOT=1,binsize=binsize,ymax=ymax,gmass=m[mct]            
            endelse
        endif else start = start+1
        cd,'../..'
    endfor
legend,['1E5R '+trim(integrate[4]),'1E4R '+trim(integrate[3]),'1E3R '+trim(integrate[2]),'1E2E '+trim(integrate[1]),'2E5R '+trim(integrate[0])],linestyle=[0,1,2,3,4],/right
    device,/close
;    stop
endfor
END
