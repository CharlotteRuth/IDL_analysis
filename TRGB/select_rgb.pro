FUNCTION SELECT_RGB,ORIGINAL=original,INTERPOL=INTERPOL
;get rid of red clump and AGB so we can interpolate between isochrones
;ORIGINAL -- set keyword to use the old set of Padova isochrones
;INTERPOL -- set this to filename of new proc_interp.pl'ed file.  

IF KEYWORD_SET(original) THEN BEGIN
    oiso=READ_PADOVA()
    ind=WHERE(oiso.logage EQ 10.0)
    iso=oiso[ind]
ENDIF ELSE BEGIN
    ;read in interpolated 10 Gyr Padova Isochrones
    iso=READ_PADOVA(/interpol)
ENDELSE

uniqz=iso[UNIQ(iso.z)].z
nz=N_ELEMENTS(uniqz)

goodind=[0l] ;use to store indices at each z
FOR i=0,nz-1 DO BEGIN
    atz=WHERE(iso.z EQ uniqz[i],natz)
    ind=INDGEN(natz-1)
    diff=iso[atz[ind]].f814mag-iso[atz[ind+1]].f814mag
    min=MIN(diff,minpos)
;    plot,(iso.f606mag-iso.f814mag)[atz],iso[atz].f814mag,xrange=[0,2],yrange=[0,-6],title=STRTRIM(uniqz[i]),color=100
;    oplot,(iso.f606mag-iso.f814mag)[atz[0:minpos]],iso[atz[0:minpos]].f814mag
    goodind=[goodind,atz[0:minpos]]
ENDFOR
goodind=goodind[1:*]
iso=iso[goodind]

;plot,iso.f606mag-iso.f814mag,iso.f814mag,/nodata,xrange=[0,4],yrange=[0,-6],xtitle='F606W-F814W',ytitle='F814W'
loadct,12
FOR i=0,nz-1 DO BEGIN
    atz=WHERE(iso.z EQ uniqz[i],natz)
;    oplot,(iso.f606mag-iso.f814mag)[atz],iso[atz].f814mag,color=200
ENDFOR


RETURN,iso

END
