;use with cooltable_read.c

fileext = '_lowt_0.txt'
fileext = '_0.txt'
coolname = 'cooltable' + fileext
heatname = 'heattable' + fileext

readcol,coolname,temp
nnH= fix(temp[0])
nHminlog = temp[1]
dnH = temp[2]
nt = fix(temp[3])
tminlog = temp[4]
dt = temp[5]
nHarray = findgen(nnH)*dnH + nHminlog
tarray = findgen(nt)*dt + tminlog

coolarray = temp[6:n_elements(temp) - 1]
cooltable = transpose(reform(coolarray,nt,nnH))

loadct,39
formatplot
min = -80 ;70
max = -45
window,0
contour,cooltable,nHarray,tarray,/fill,nlevels = 254,xrange = [-9,4],yrange = [1,9],xstyle = 1,ystyle = 1,xtitle = 'log Density [amu/cc]',ytitle = 'log Temperature [K]',min_value = min, max_value = max,title = 'log Cool';,xrange = [min(nHarray),max(nHarray)],yrange = [min(tarray),max(tarray)

endcell = 115
xnHlog = nnH - endcell
wnHlog1 = xnHlog
wnHlog0 = 1-xnhlog
estimate_coolarray = wnHlog0*cooltable[endcell,*] + wnHlog1*cooltable[endcell + 1,*]
print,transpose([[tarray],$
                 [reform(cooltable[nnH - 1,*])],$
                 [reform(estimate_coolarray)],$
                 [reform(estimate_coolarray - cooltable[nnH - 1,*])],$
                 [reform(exp(estimate_coolarray) - exp(cooltable[nnH - 1,*]))/exp(cooltable[nnH - 1,*])]])
print,minmax(estimate_coolarray - cooltable[nnH - 1,*])
print,minmax((exp(estimate_coolarray) - exp(cooltable[nnH - 1,*]))/exp(cooltable[nnH - 1,*]))

;-------------------------------------------------------------------
readcol,heatname,temp
nnH= fix(temp[0])
nHminlog = temp[1]
dnH = temp[2]
nt = fix(temp[3])
tminlog = temp[4]
dt = temp[5]
nHarray = findgen(nnH)*dnH + nHminlog
tarray = findgen(nt)*dt + tminlog

heatarray = temp[6:n_elements(temp) - 1]
heattable = transpose(reform(heatarray,nt,nnH))

min = -80
max = -45
window,1
contour,heattable,nHarray,tarray,/fill,nlevels = 254,xrange = [-9,4],yrange = [1,9],xstyle = 1,ystyle = 1,xtitle = 'log Density [amu/cc]',ytitle = 'log Temperature [K]',min_value = min, max_value = max, title = 'log Heat';,xrange = [min(nHarray),max(nHarray)],yrange = [min(tarray),max(tarray)

endcell = 115
xnHlog = nnH - endcell
wnHlog1 = xnHlog
wnHlog0 = 1-xnhlog
estimate_heatarray = wnHlog0*heattable[endcell,*] + wnHlog1*heattable[endcell + 1,*]
print,transpose([[tarray],$
                 [reform(heattable[nnH - 1,*])],$
                 [reform(estimate_heatarray)],$
                 [reform(estimate_heatarray - heattable[nnH - 1,*])],$
                 [reform(exp(estimate_heatarray) - exp(heattable[nnH - 1,*]))/exp(heattable[nnH - 1,*])]])
print,minmax(estimate_heatarray - heattable[nnH - 1,*])
print,minmax((exp(estimate_heatarray) - exp(heattable[nnH - 1,*]))/exp(heattable[nnH - 1,*]))

;--------------------------------------------------------
print,transpose([[tarray],[reform(((exp(estimate_heatarray) - exp(estimate_coolarray)) - (exp(heattable[nnH - 1,*]) - exp(cooltable[nnH - 1,*])))/(exp(heattable[nnH - 1,*]) - exp(cooltable[nnH - 1,*])))]])
print,minmax(((exp(estimate_heatarray) - exp(estimate_coolarray)) - (exp(heattable[nnH - 1,*]) - exp(cooltable[nnH - 1,*])))/(exp(heattable[nnH - 1,*]) - exp(cooltable[nnH - 1,*])))
