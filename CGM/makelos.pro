
;impact_param = [2.0,4.0,8.0,14.0,20.0,50.0,100.0,150.0,200.0,250.0,300.0]
;kpcunit = 50000.
;nangle = 32
;outfile = '1_11.83_10.62_90' halonumber_logDM_logstar_angle
;outfile = '1_11.48_9.89_90'
;;;;LOSfile = 'LOS.zoom.10kpc.90.36angles'
;LOSfile = 'LOS.h603.g14HMbwK.00512.1'

PRO makeLOS,LOSfile,outfile,impact_param,kpcunit,nangle

ni = n_elements(impact_param)
x = fltarr(ni*nangle + 1)
y = fltarr(ni*nangle + 1)
z = fltarr(ni*nangle + 1)
axis = fltarr(ni*nangle + 1) + 2

FOR i = 0, ni - 1 DO BEGIN
    angles = findgen(nangle)*2*!PI/nangle
    x[i*nangle + 1 : nangle + i*nangle] = cos(angles)*impact_param[i]/kpcunit
    y[i*nangle + 1 : nangle + i*nangle] = sin(angles)*impact_param[i]/kpcunit
ENDFOR
xstring = strtrim(round(x*kpcunit*100)/100.0,2)
FOR i = 0, n_elements(xstring) - 1 DO xstring[i] = strmid(xstring[i],0,strlen(xstring[i]) - 3)
pos = where(x GT 0)
xstring[pos] = '+' + xstring[pos]
xstring = 'x'+xstring+'kpc'

ystring = strtrim(round(y*kpcunit*100)/100.0,2)
FOR i = 0, n_elements(ystring) - 1 DO ystring[i] = strmid(ystring[i],0,strlen(ystring[i]) - 3)
pos = where(y GT 0)
ystring[pos] = '+' + ystring[pos]
ystring = 'y'+ystring+'kpc'

filename = outfile+xstring+'_'+ystring
writecol,LOSfile,x,y,z,axis,filename,format = '(F,F,F,I,A40)'
END
