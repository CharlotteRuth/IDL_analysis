
;impact_param = 10.
;kpcunit = 50000.
;nangle = 32
;filebase = '1_11.83_10.62_90'
;outfile = 'LOS.zoom.10kpc.90.36angles'

PRO makeLOS,outfile,filebase,impact_param,kpcunit,nangle

angles = findgen(nangle)*2*!PI/nangle
x = cos(angles)*impact_param/kpcunit
y = sin(angles)*impact_param/kpcunit
z = angles*0
axis = fix(angles*0) + 2 ;2 means the line of sight goes down the z axis

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

filename = filebase+xstring+'_'+ystring
writecol,outfile,x,y,z,axis,filename,format = '(F,F,F,I,A40)'
END
