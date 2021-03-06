pro cooltable_read,outplot= outplot

formatplot,outplot = outplot

M_H = 1.6719999999999997e-24
Y_H = 0.68649999884516
cm_per_kpc = 3.0857d21
kpc_per_syslength = 1e5
h = 2.60000e-07*cm_per_kpc*kpc_per_syslength
loadct,0
;filename_c = '/astro/net/scratch2/christensen/MolecH/cooltable.txt'
;filename_h = '/astro/net/scratch2/christensen/MolecH/heattable.txt'
filename_c = '/home/christensen/Code/changa_uw_bk_NR/cool_z0.txt'
filename_h = '/home/christensen/Code/changa_uw_bk_NR/heat_z0.txt'
openr, lun_c, filename_c, /GET_LUN
readf,lun_c, nnH, nT
openr, lun_h, filename_h, /GET_LUN
readf,lun_h, dummy

;nnH = 121
nHminlog = -9
nHmaxlog = 3
dnH = 0.1

;nt = 142
tminlog = 2
tmaxlog = 9.05
dt = 0.05

h_array = findgen(nnH)*dnH + nHminlog
temp_array = findgen(nt)*dt + tminlog
h_table = fltarr(nnH,nt)
for i = 0, N_ELEMENTS(temp_array)-1 do h_table[*,i] = 10^h_array
;cooltable = fltarr(nt,nnH)
;heattable = fltarr(nt,nnH)

cooltable = ''
WHILE NOT EOF(lun_c) DO BEGIN
    readf,lun_c,line
    cooltable = [cooltable,line]
ENDWHILE
;readf,1,cooltable
close,lun_c
cooltable = cooltable[1:n_elements(cooltable)-1]
cooltable = reform(cooltable,nnh,nt)
;cooltable2[*,0] = same temp, changing density
;cooltable2[0,*] = same density, changing temp

heattable = ''
WHILE NOT EOF(lun_h) DO BEGIN
    readf,lun_h,line
    heattable = [heattable,line]
ENDWHILE
;readf,1,heattable
close,lun_h
heattable = heattable[1:n_elements(heattable)-1]
heattable = reform(heattable,nnh,nt)

;cooltable = rotate(cooltable,4)
;heattable = rotate(heattable,4)
cooltablefinal = exp(cooltable);*Y_H/M_H*h_table
heattablefinal = exp(heattable);*Y_H/M_H*h_table

;***************************************** Cooling ********************************

IF keyword_set(outplot) THEN device,filename=outplot + 'metalcooling.eps',/color,bits_per_pixel= 8,/times ELSE window,0
contour,alog(cooltablefinal),H_array,temp_array,ytitle="Log(Temperature)",xtitle = 'Log(nH)',/FILL,nlevels = 60,title = 'Log(Metal Cooling Coeff)'
oplot,-1.0*alog(cooltablefinal[120,*]) - MIN(alog(cooltablefinal[120,*])),temp_array
oplot,H_array,alog(cooltablefinal[*,141])/6.0 - MIN(alog(cooltablefinal[*,141]))/6.0+2.0
IF keyword_set(outplot) THEN device,/close

extrap_H_min = -2
extrap_H_max = 4
extrap_temp_min = 0
extrap_temp_max = 3.0 ;5
N_el = 100
poly = 2
sub_cooltablefinal = cooltablefinal[71:120,0:30]
sub_H_array = H_array[71:120]
sub_temp_array = temp_array[0:30]

z = reform(alog(sub_cooltablefinal),50*31)
input = fltarr(3,50*31)
iH = 0
iT = 0
for i = 0, 50*31 -1 DO BEGIN
    input[0,i] = sub_H_array[iH]
    input[1,i] = sub_temp_array[iT]
    input[2,i] = z[i]
    iH = iH + 1
    if (iH eq 50) then begin
        iH = 0
        iT = iT + 1
    endif
ENDFOR
results = reform(SFIT(input,/IRREGULAR,poly,kx = kx),50,31)

extrap_H_array = findgen(n_el)*(extrap_H_max - extrap_H_min)/n_el + extrap_H_min
extrap_temp_array = findgen(n_el)*(extrap_temp_max - extrap_temp_min)/n_el + extrap_temp_min
extrap = fltarr(n_el,n_el)
sub_2 = fltarr(N_ELEMENTS(sub_H_array),N_ELEMENTS(sub_temp_array))
print,'Metal Cooling Fit'
print,kx
FOR iH = 0, n_el - 1 DO BEGIN
    FOR iT = 0, n_el - 1 DO BEGIN
        FOR cH = 0, poly  DO BEGIN
            FOR cT = 0, poly  DO BEGIN
                extrap[iH,iT] = extrap[iH,iT] + kx[cT,cH]*extrap_H_array[iH]^cH*extrap_temp_array[iT]^cT
            ENDFOR
        ENDFOR
    ENDFOR
ENDFOR

FOR cH = 0, poly  DO BEGIN
    FOR cT = 0, poly  DO BEGIN
        FOR iH = 0, N_ELEMENTS(sub_H_array) - 1 DO BEGIN
            FOR iT = 0, N_ELEMENTS(sub_temp_array) - 1 DO BEGIN
                sub_2[iH,iT] = sub_2[iH,iT] + kx[cT,cH]*sub_H_array[iH]^cH*sub_temp_array[iT]^cT
            ENDFOR
        ENDFOR
    ENDFOR
ENDFOR

IF keyword_set(outplot) THEN device,filename=outplot+'metalcooling_fit.eps',/color,bits_per_pixel= 8,/times ELSE window,1
contour,alog(sub_cooltablefinal),sub_H_array,sub_temp_array,ytitle="Log(Temperature)",xtitle = 'Log(nH)',/FILL,nlevels = 60,title = 'Log(Metal Cooling Coeff)',xrange = [extrap_H_min, extrap_H_max],yrange = [extrap_temp_min,extrap_temp_max]
contour,extrap,extrap_H_array,extrap_temp_array,/overplot,nlevels =12,c_labels = [1,1,1,1,1,1,1,1,1,1,1,1]
IF keyword_set(outplot) THEN device,/close

loadct,39
IF keyword_set(outplot) THEN device,filename=outplot+'metalcooling_temp.eps',/color,bits_per_pixel= 8,/times ELSE window,2
plot,temp_array,alog(cooltablefinal[120,*]),xrange = [extrap_temp_min,extrap_temp_max],xtitle = 'Log(Temperature)',yrange = [-75,-55],ytitle = 'Log(Metal Cooling Coeff)'
y = MIN(ABS(extrap_H_array - 3),x)
oplot,extrap_temp_array,extrap[x,*],linestyle = 1

oplot,temp_array,alog(cooltablefinal[120,*]),color = 254
oplot,extrap_temp_array,(extrap_temp_array - tminlog)*(alog(cooltablefinal[120,1]) - alog(cooltablefinal[120,0]))/dt + alog(cooltablefinal[120,0]),color = 254,linestyle = 1

oplot,temp_array,alog(cooltablefinal[110,*]),color = 220
oplot,extrap_temp_array,(extrap_temp_array - tminlog)*(alog(cooltablefinal[110,1]) - alog(cooltablefinal[110,0]))/dt + alog(cooltablefinal[110,0]),color = 220,linestyle = 1

oplot,temp_array,alog(cooltablefinal[100,*]),color = 200
oplot,extrap_temp_array,(extrap_temp_array - tminlog)*(alog(cooltablefinal[100,1]) - alog(cooltablefinal[100,0]))/dt + alog(cooltablefinal[100,0]),color = 200,linestyle = 1

oplot,temp_array,alog(cooltablefinal[90,*]),color = 180
oplot,extrap_temp_array,(extrap_temp_array - tminlog)*(alog(cooltablefinal[90,1]) - alog(cooltablefinal[90,0]))/dt + alog(cooltablefinal[90,0]),color = 180,linestyle = 1

oplot,temp_array,alog(cooltablefinal[80,*]),color = 160
oplot,extrap_temp_array,(extrap_temp_array - tminlog)*(alog(cooltablefinal[80,1]) - alog(cooltablefinal[80,0]))/dt + alog(cooltablefinal[80,0]),color = 160,linestyle = 1

oplot,[1.0,1.0],[-100,100]
legend,['Original','Fit'],linestyle = [0,1]
print,h_array[80],h_array[90],h_array[100],h_array[110],h_array[120]
IF keyword_set(outplot) THEN device,/close
stop

IF keyword_set(outplot) THEN device,filename=outplot+'metalcooling_rho.eps',/color,bits_per_pixel= 8,/times ELSE window,3
plot,H_array,alog(cooltablefinal[*,0]),xrange = [extrap_H_min,extrap_H_max],xtitle = 'Log(Density)',yrange = [-6,2],ytitle = 'Log(Metal Cooling Coeff)'
y = MIN(ABS(extrap_temp_array - 2),x)
oplot,extrap_H_array,extrap[*,x],linestyle = 1
legend,['Original','Fit'],linestyle = [0,1]
IF keyword_set(outplot) THEN device,/close else stop

;*********************** Heating *******************************

IF keyword_set(outplot) THEN device,filename=outplot + 'metalheating.eps',/color,bits_per_pixel= 8,/times ELSE window,0
contour,alog(heattablefinal),H_array,temp_array,ytitle="Log(Temperature)",xtitle = 'Log(nH)',/FILL,nlevels = 60,title = 'Log(Metal Heating Coeff)'
oplot,-1.0*alog(heattablefinal[120,*])/2 - MIN(alog(heattablefinal[120,*]))/2 - 7.3866,temp_array
oplot,H_array,alog(heattablefinal[*,141])/6.0 - MIN(alog(heattablefinal[*,141]))/6.0+2.0
IF keyword_set(outplot) THEN device,/close

extrap_H_min = -2
extrap_H_max = 4
extrap_temp_min = 0
extrap_temp_max = 3.5
N_el = 100
poly = 2
sub_heattablefinal = heattablefinal[71:120,0:30]
sub_H_array = H_array[71:120]
sub_temp_array = temp_array[0:30]

z = reform(alog(sub_heattablefinal),50*31)
input = fltarr(3,50*31)
iH = 0
iT = 0
for i = 0, 50*31 -1 DO BEGIN
    input[0,i] = sub_H_array[iH]
    input[1,i] = sub_temp_array[iT]
    input[2,i] = z[i]
    iH = iH + 1
    if (iH eq 50) then begin
        iH = 0
        iT = iT + 1
    endif
ENDFOR
results = reform(SFIT(input,/IRREGULAR,poly,kx = kx),50,31)

extrap_H_array = findgen(n_el)*(extrap_H_max - extrap_H_min)/n_el + extrap_H_min
extrap_temp_array = findgen(n_el)*(extrap_temp_max - extrap_temp_min)/n_el + extrap_temp_min
extrap = fltarr(n_el,n_el)
sub_2 = fltarr(N_ELEMENTS(sub_H_array),N_ELEMENTS(sub_temp_array))
print,'Metal Heating Fit'
print,kx
FOR iH = 0, n_el - 1 DO BEGIN
    FOR iT = 0, n_el - 1 DO BEGIN
        FOR cH = 0, poly  DO BEGIN
            FOR cT = 0, poly  DO BEGIN
                extrap[iH,iT] = extrap[iH,iT] + kx[cT,cH]*extrap_H_array[iH]^cH*extrap_temp_array[iT]^cT
            ENDFOR
        ENDFOR
    ENDFOR
ENDFOR

FOR cH = 0, poly  DO BEGIN
    FOR cT = 0, poly  DO BEGIN
        FOR iH = 0, N_ELEMENTS(sub_H_array) - 1 DO BEGIN
            FOR iT = 0, N_ELEMENTS(sub_temp_array) - 1 DO BEGIN
                sub_2[iH,iT] = sub_2[iH,iT] + kx[cT,cH]*sub_H_array[iH]^cH*sub_temp_array[iT]^cT
            ENDFOR
        ENDFOR
    ENDFOR
ENDFOR

IF keyword_set(outplot) THEN device,filename=outplot+'metalheating_fit.eps',/color,bits_per_pixel= 8,/times ELSE window,1
contour,alog(sub_heattablefinal),sub_H_array,sub_temp_array,ytitle="Log(Temperature)",xtitle = 'Log(nH)',/FILL,nlevels = 60,title = 'Log(Metal Heating Coeff)',xrange = [extrap_H_min, extrap_H_max],yrange = [extrap_temp_min,extrap_temp_max]
contour,extrap,extrap_H_array,extrap_temp_array,/overplot,nlevels =12,c_labels = [1,1,1,1,1,1,1,1,1,1,1,1]
IF keyword_set(outplot) THEN device,/close

IF keyword_set(outplot) THEN device,filename=outplot+'metalheating_T.eps',/color,bits_per_pixel= 8,/times ELSE window,2
plot,temp_array,alog(heattablefinal[120,*]),xrange = [extrap_temp_min,extrap_temp_max],xtitle = 'Log(Temperature)',ytitle = "Log(Metal Heating Coeff)",yrange = [-11,-10]
y = MIN(ABS(extrap_H_array - 3),x)
oplot,extrap_temp_array,extrap[x,*],linestyle = 1
legend,['Original','Fit'],linestyle = [0,1]
IF keyword_set(outplot) THEN device,/close

IF keyword_set(outplot) THEN device,filename=outplot+'metalheating_rho.eps',/color,bits_per_pixel= 8,/times ELSE window,3
plot,H_array,alog(heattablefinal[*,0]),xrange = [extrap_H_min,extrap_H_max],xtitle = 'Log(Density)',ytitle = "Log(Metal Heating Coeff)",yrange = [-15,9]
y = MIN(ABS(extrap_temp_array - 2),x)
oplot,extrap_H_array,extrap[*,x],linestyle = 1
legend,['Original','Fit'],linestyle = [0,1]
IF keyword_set(outplot) THEN device,/close ELSE stop
end
