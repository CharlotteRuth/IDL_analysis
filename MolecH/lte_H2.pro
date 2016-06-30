
function S_H2,yH2,h
omega_H2 = 0.2
x = yH2*h/5d14
return, (1 - omega_H2)/(1 + x)^2 + omega_H2/SQRT(1 + x)*exp(-0.00085*SQRT(1 + x))
end

function S_H22,yH2,h
return, (yH2*h/1e14)^(-0.75)
end

function S_d,yHI,yH2,h;,z
ZSOLAR = 0.0130215 ;0.0177 ;0.0130215
sigma_d = 2d-21
;return, exp(-1.0*sigma_d*z/ZSOLAR*(yHI*h + 2.0*yH2*h))
return, exp(-1.0*sigma_d*(yHI*h + 2.0*yH2*h))
end

function lte_H2,rho,zmetal,h,lw,yH
yHI = 1d0
yH2 = 1d-6
CL_B_gm  = 6.022d23*(938.7830/931.494)
ZSOLAR = 0.0130215 ;0.0177 ;0.0130215
eps = 1e-5
cp = 10.0
;print,rho,zmetal,h,lw,yH

for  i=0, 20 do begin
    yHI_old   = yHI             
    yH2_old   = yH2             
   
    S_H2 = S_H2(yH2*rho*yH,h) ;self shielding
    S_H22 = S_H22(yH2*rho*yH,3e21)   ;self shielding
    S_d  = S_d(yHI*rho*yH,yH2*rho*yH,h) ;shielding from dust
    DustForm = 3.5e-17*zmetal/zsolar*cp
    PhotDissoc = lw*1.175577d-16;1.1d8*12.87*1.6021746d-12*4.0*!PI*4.13e-5 

    IF(s_d eq 0 OR s_H2 eq 0) THEN BEGIN
        yHI = 0 
 ;       fHI = -1
    ENDIF ELSE  BEGIN
        fHI = 2.0*DustForm*rho*yH/(s_d*S_H22*PhotDissoc) ;2.0*(DustForm*yHI_old)/(PhotDissoc*yH*yH2_old*S_d*S_H2)
;    IF(s_d eq 0 OR s_H2 eq 0) THEN fHI = 1
        yHI = yH/(1 + fHI)
    ENDELSE
    yH2 = (yH - yHI)/2.

   ; if ( abs(yHI_old-yHI) lt EPS * yHI ) then  break
    if (YH2 lt 1e-20) then YH2 = 1e-20                 ;
    if (YHI lt 1e-20) then yHI = 1e-20 ;
    if (YH2 gt yH/2.0) then YH2 = yH/2.0 
    if (YHI gt yH) then YHI = yH 
;    print,fHI,s_H2,s_d,DustForm,PhotDissoc,yHI_old,yH2_old
;    print,yHI,yH2
;    stop
endfor
Y_HI = yHI                      
Y_H2 = yH2                      
IF(~FINITE(Y_HI) OR ~FINITE(Y_H2)) THEN stop
return,[Y_H2,yHI,s_H2,s_d,2.0*DustForm*rho*yH,s_d*S_H2*PhotDissoc,S_H22]
end

pro check_lte_H2
loadct,39


file = '/astro/net/scratch2/christensen/MolecH/11M/Disk_Iso_1e5_zsol/lte.txt'
rtipsy,'/astro/net/scratch2/christensen/MolecH/11M/Disk_Iso_1e5_zsol/MW_disk_zsol.std',h,g,d,s
readcol,file,en_B,ye,H2,HI,HII,HeI,HeII,HeIII,zmetal,CorreLength,s_dust,s_self,Coll_HI,Coll_HeI,Coll_HeII,Coll_e_H2,Coll_H_H2,Coll_H2_H2,Coll_HII_H2,Radr_HII,Totr_HeII,Radr_HeII,Radr_HeIII,DustForm_H2,Phot_HI,Phot_HeI,Phot_HeII,Phot_H2

fHII = (Radr_HII*en_B*ye)/(Phot_HI*s_dust + Coll_HI*en_B*ye)          ;  
fHI = 2.0*(DustForm_H2*en_B*(HI + 2.0*H2))/((Coll_e_H2*en_B*ye + Coll_H_H2*en_B*HI + Coll_H2_H2*en_B*H2 + Phot_H2)*s_dust*s_self)    ;
rfH  =  1 / ( 1 + fHII * (1 + fHI) ) ; 
yH = (HI + 2.0*H2 + HII)
yHII =  yH * rfH                ;
yHI  =  yH * rfH * fHII         ;
yH2  = (yH - yHI - yHII)/2.0    ;

ind = where(H2 gt 1e-7)
stop
lteH2 = fltarr(N_ELEMENTS(ind))
FOR i = 0, N_ELEMENTS(ind)-1 DO lteH2[i] = lte_H2(en_B[ind[i]],zmetal[ind[i]],CorreLength[ind[i]],Phot_H2[ind[i]]/1.175577d-16,(HI[ind[i]] + 2.0*H2[ind[i]] + HII[ind[i]]))
plot,alog10(CorreLength[ind]*en_B[ind]),H2[ind]/(H2[ind] + HI[ind])*2, xtitle = textoidl("N_{HI} + 2N_{H_2} (cm^{-2})"),ytitle = textoidl('f_{H_2}'),yrange = [1e-6,1.0],/ylog,xstyle=1,ystyle=1,xrange = [19,23],psym = 3
oplot,alog10(CorreLength[ind]*en_B[ind]),lteH2/(H2[ind] + HI[ind])*2,psym = 3,color = 240


end

pro plot_lte_H2,outfile = outfile
loadct,39
if (KEYWORD_SET(outfile)) then begin
    !Y.STYLE = 3
    !X.STYLE = 3
    !P.THICK = 3.5
    !P.CHARTHICK=4
    !X.THICK=4
    !Y.THICK=4
    !P.charsize=1.0
    !x.charsize=2.25
    !y.charsize=2.25
    set_plot,'ps' 
;    if (KEYWORD_SET(outfile)) then device,filename = outfile+'_fH2_col.eps',/color,bits_per_pixel= 8,/times,xsize = 19,ysize = 12,xoffset =  2,yoffset =  2
endif else begin
    set_plot,'x'
    !P.THICK = 3.5
    !P.CHARTHICK=1.5
    !X.THICK=1.5
    !Y.THICK=1.5
    window,0
endelse

column_den= 10^(findgen(5000)/5000*(23-19) + 19.0)
z = 0.0177 ;10^(-2.8)                 ;MEAN(g[where(frac_ind gt 1e-6)].zmetal)
ZSOLAR = 0.0177 ;0.0130215
total_H = 0.70
sigma_d = 2d-21
shield = exp(-1.0*sigma_d*z/ZSOLAR*(column_den))


plot,alog10(column_Den),1.0-shield, xtitle = textoidl("N_{HI} + 2N_{H_2} (cm^{-2})"),ytitle = textoidl('f_{H_2}'),yrange = [1e-6,1.0],xstyle=1,ystyle=1,xrange = [19,23]
shield = exp(-1.0*sigma_d*(column_den))
oplot,alog10(column_Den),1.0-shield
lte_H2 = column_den
s_H2 = column_den
s_H22 = column_den
s_d = column_den
form = column_den
dis = column_den
   
if (KEYWORD_SET(outfile)) then device,filename = outfile+'_fH2_col_516.eps',/color,bits_per_pixel= 8,/times,xsize = 19,ysize = 12,xoffset =  2,yoffset =  2
lw = 1e5
length = 1e20
z = 0.016*0.0177
FOR i = 0, N_ELEMENTS(column_den) - 1 DO BEGIN
  result = lte_H2(column_den[i]/length,z,length,lw,total_H)
  lte_H2[i] = result[0]
  s_H2[i] = result[2]
  s_d[i] = result[3]
  form[i] = result[4]
  dis[i] = result[5]
  s_H22[i] = result[6]
ENDFOR
oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 2,color =240

lw = 1e6
z = 0.1*0.0177
length = 1e20
FOR i = 0, N_ELEMENTS(column_den) - 1 DO BEGIN
  result = lte_H2(column_den[i]/length,z,length,lw,total_H)
  lte_H2[i] = result[0]
  s_H2[i] = result[2]
  s_d[i] = result[3]
  form[i] = result[4]
  dis[i] = result[5]
  s_H22[i] = result[6]
ENDFOR
oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 0,color =240
if (KEYWORD_SET(outfile)) then device,/close ELSE stop
   if (KEYWORD_SET(outfile)) then device,filename = outfile+'_fH2_col.eps',/color,bits_per_pixel= 8,/times,xsize = 19,ysize = 12,xoffset =  2,yoffset =  2
;-----------------------------------------------------------
plot,alog10(column_Den),1.0-shield, xtitle = textoidl("N_{HI} + 2N_{H_2} (cm^{-2})"),ytitle = textoidl('f_{H_2}'),yrange = [1e-6,1.0],/ylog,xstyle=1,ystyle=1,xrange = [19,23]
shield = exp(-1.0*sigma_d*(column_den))
oplot,alog10(column_Den),1.0-shield
lte_H2 = column_den
s_H2 = column_den
s_H22 = column_den
s_d = column_den
form = column_den
dis = column_den

lw = 1e5
;length = 1e19
length = 1e20
z = 0.0177*0.1
FOR i = 0, N_ELEMENTS(column_den) - 1 DO BEGIN
  result = lte_H2(column_den[i]/length,z,length,lw,total_H)
  lte_H2[i] = result[0]
  s_H2[i] = result[2]
  s_d[i] = result[3]
  form[i] = result[4]
  dis[i] = result[5]
  s_H22[i] = result[6]
ENDFOR
oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 4,color =240
;length = 1e20
z = 0.0177*0.016
FOR i = 0, N_ELEMENTS(column_den) - 1 DO BEGIN
    result = lte_H2(column_den[i]/length,z,length,lw,total_H)
    lte_H2[i] = result[0]
    s_H2[i] = result[2]
    s_d[i] = result[3]
    form[i] = result[4]
    dis[i] = result[5]
  s_H22[i] = result[6]
ENDFOR
oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 2,color =240
;length = 1e21
z = 0.0177*0.1
z = 0.0177*0.01

FOR i = 0, N_ELEMENTS(column_den) - 1 DO BEGIN
    result = lte_H2(column_den[i]/length,z,length,lw,total_H)
    lte_H2[i] = result[0]
    s_H2[i] = result[2]
    s_d[i] = result[3]
    form[i] = result[4]
    dis[i] = result[5]
  s_H22[i] = result[6]
ENDFOR
oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 1,color =240
oplot,alog10(column_Den),1 - s_d*s_H2,linestyle = 0,color = 240

lw = 1e7
length = 1e19
length = 1e20
z = 0.0177
FOR i = 0, N_ELEMENTS(column_den) - 1 DO BEGIN
    result = lte_H2(column_den[i]/length,z,length,lw,total_H)
    lte_H2[i] = result[0]
    s_H2[i] = result[2]
    s_d[i] = result[3]
    form[i] = result[4]
    dis[i] = result[5]
  s_H22[i] = result[6]
ENDFOR
oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 4,color =190
;length = 1e20
z = 0.0177*0.3
FOR i = 0, N_ELEMENTS(column_den) - 1 DO BEGIN
    result = lte_H2(column_den[i]/length,z,length,lw,total_H)
    lte_H2[i] = result[0]
    s_H2[i] = result[2]
    s_d[i] = result[3]
    form[i] = result[4]
    dis[i] = result[5]
  s_H22[i] = result[6]
ENDFOR
oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 2,color =190
;length = 1e21
z = 0.0177*0.1
z = 0.0177*0.01
FOR i = 0, N_ELEMENTS(column_den) - 1 DO BEGIN
    result = lte_H2(column_den[i]/length,z,length,lw,total_H)
    lte_H2[i] = result[0]
    s_H2[i] = result[2]
    s_d[i] = result[3]
    form[i] = result[4]
    dis[i] = result[5]
    s_H22[i] = result[6]
ENDFOR
oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 1,color =190
oplot,alog10(column_Den),1 - s_d*s_H2,linestyle = 0,color = 190
;length = 1e22
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 3,color =240

lw = 1e9
length = 1e19
length = 1e20
z = 0.0177
FOR i = 0, N_ELEMENTS(column_den) - 1 DO BEGIN
    result = lte_H2(column_den[i]/length,z,length,lw,total_H)
    lte_H2[i] = result[0]
    s_H2[i] = result[2]
    s_d[i] = result[3]
    form[i] = result[4]
    dis[i] = result[5]
    s_H22[i] = result[6]
ENDFOR
oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 4,color =150
;length = 1e20
z = 0.0177*0.3
FOR i = 0, N_ELEMENTS(column_den) - 1 DO BEGIN
    result = lte_H2(column_den[i]/length,z,length,lw,total_H)
    lte_H2[i] = result[0]
    s_H2[i] = result[2]
    s_d[i] = result[3]
    form[i] = result[4]
    dis[i] = result[5]
    s_H22[i] = result[6]
ENDFOR
oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 2,color =150
;length = 1e21
z = 0.0177*0.1
z = 0.0177*0.01
FOR i = 0, N_ELEMENTS(column_den) - 1 DO BEGIN
    result = lte_H2(column_den[i]/length,z,length,lw,total_H)
    lte_H2[i] = result[0]
    s_H2[i] = result[2]
    s_d[i] = result[3]
    form[i] = result[4]
    dis[i] = result[5]
    s_H22[i] = result[6]
ENDFOR
oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 1,color =150
oplot,alog10(column_Den),1 - s_d*s_H2,linestyle = 0,color = 150
if (KEYWORD_SET(outfile)) then begin 
    device,/close 
    device,filename = outfile+'_fH2_den.eps',/color,bits_per_pixel= 8,/times,xsize = 19,ysize = 12,xoffset =  2,yoffset =  2    
endif else window,1

;length = 1e22
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 3,color =220

;length = 1e20
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = lte_H2(column_den[i]/length,z,length,1e11,total_H)
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 2,color =200
;length = 1e21
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = lte_H2(column_den[i]/length,z,length,1e11,total_H)
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 1,color =200
;length = 1e22
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = lte_H2(column_den[i]/length,z,length,1e11,total_H)
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 3,color =200

;lw = 1e13
;length = 1e19
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 4,color =180
;length = 1e20
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 2,color =180
;length = 1e21
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 1,color =180
;length = 1e22
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 3,color =180

;length = 1e20
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;lte_H2[i] = lte_H2(column_den[i]/length,z,length,1e15,total_H)
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 2,color =160
;length = 1e21
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = lte_H2(column_den[i]/length,z,length,1e15,total_H)
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 1,color =160
;length = 1e22
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = lte_H2(column_den[i]/length,z,length,1e15,total_H)
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 3,color =160

;lw = 1e17
;length = 1e19
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 4,color =140
;length = 1e20
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 2,color =140
;length = 1e21
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 1,color =140
;length = 1e22
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 3,color =140

;lw = 1e21
;length = 1e19
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 4,color =100
;length = 1e20
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 2,color =100
;length = 1e21
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 1,color =100
;length = 1e22
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 3,color =100

;window,1
plot,alog10(column_Den/length),1.0-shield, xtitle = textoidl("n_H (cm^-3)"),ytitle = textoidl('f_{H_2}'),yrange = [0,1.0],xstyle=1,ystyle=1,xrange = [-1,3]

lw = 1e5
length = 1e19
length = 1e20
z = 0.0177
FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
oplot,alog10(column_Den/length),lte_H2/total_H*2.0,linestyle = 4,color =240
;length = 1e20
z = 0.0177*0.3
FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
oplot,alog10(column_Den/length),lte_H2/total_H*2.0,linestyle = 2,color =240
;length = 1e21
z = 0.0177*0.1
z = 0.0177*0.01
FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
oplot,alog10(column_Den/length),lte_H2/total_H*2.0,linestyle = 1,color =240

lw = 1e7
;length = 1e19
length = 1e20
z = 0.0177
FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
oplot,alog10(column_Den/length),lte_H2/total_H*2.0,linestyle = 4,color =190
;length = 1e20
z = 0.0177*0.3
FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
oplot,alog10(column_Den/length),lte_H2/total_H*2.0,linestyle = 2,color =190
;length = 1e21
z = 0.0177*0.1
z = 0.0177*0.01
FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
oplot,alog10(column_Den/length),lte_H2/total_H*2.0,linestyle = 1,color =190
;length = 1e22
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den/length),lte_H2/total_H*2.0,linestyle = 3,color =240

lw = 1e9
;length = 1e19
length = 1e20
z = 0.0177
FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
oplot,alog10(column_Den/length),lte_H2/total_H*2.0,linestyle = 4,color =150
;length = 1e20
z = 0.0177*0.3
FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
oplot,alog10(column_Den/length),lte_H2/total_H*2.0,linestyle = 2,color =150
;length = 1e21
z = 0.0177*0.1
z = 0.0177*0.01
FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
oplot,alog10(column_Den/length),lte_H2/total_H*2.0,linestyle = 1,color =150
;length = 1e22
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = lte_H2(column_den[i]/length,z,length,1e9,total_H)
;oplot,alog10(column_Den/length),lte_H2/total_H*2.0,linestyle = 3,color =220

;length = 1e20
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = lte_H2(column_den[i]/length,z,length,1e11,total_H)
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 2,color =200
;length = 1e21
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = lte_H2(column_den[i]/length,z,length,1e11,total_H)
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 1,color =200
;length = 1e22
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = lte_H2(column_den[i]/length,z,length,1e11,total_H)
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 3,color =200

;length = 1e20
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = lte_H2(column_den[i]/length,z,length,1e13,total_H)
;oplot,alog10(column_Den/length),lte_H2/total_H*2.0,linestyle = 2,color =180
;length = 1e21
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = lte_H2(column_den[i]/length,z,length,1e13,total_H)
;oplot,alog10(column_Den/length),lte_H2/total_H*2.0,linestyle = 1,color =180
;length = 1e22
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = lte_H2(column_den[i]/length,z,length,1e13,total_H)
;oplot,alog10(column_Den/length),lte_H2/total_H*2.0,linestyle = 3,color =180

;length = 1e20
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;lte_H2[i] = lte_H2(column_den[i]/length,z,length,1e15,total_H)
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 2,color =160
;length = 1e21
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = lte_H2(column_den[i]/length,z,length,1e15,total_H)
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 1,color =160
;length = 1e22
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = lte_H2(column_den[i]/length,z,length,1e15,total_H)
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 3,color =160

;lw = 1e17
;length = 1e19
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den/length),lte_H2/total_H*2.0,linestyle = 4,color =140
;length = 1e20
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den/length),lte_H2/total_H*2.0,linestyle = 2,color =140
;length = 1e21
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den/length),lte_H2/total_H*2.0,linestyle = 1,color =140
;length = 1e22
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den/length),lte_H2/total_H*2.0,linestyle = 3,color =140

;length = 1e19
;lw = 1e21
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den/length),lte_H2/total_H*2.0,linestyle = 4,color =100
;length = 1e20
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den/length),lte_H2/total_H*2.0,linestyle = 2,color =100
;length = 21
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den/length),lte_H2/total_H*2.0,linestyle = 1,color =100
;length = 1e22
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den/length),lte_H2/total_H*2.0,linestyle =
;3,color =100
if (KEYWORD_SET(outfile)) then begin 
    device,/close 
    device,filename = outfile+'_fH2_den_lin.eps',/color,bits_per_pixel= 8,/times,xsize = 19,ysize = 12,xoffset =  2,yoffset =  2    
endif else window,2

plot,alog10(column_Den),1.0-shield, xtitle = textoidl("N_{HI} + 2N_{H_2} (cm^{-2})"),ytitle = textoidl('f_{H_2}'),yrange = [1e-6,1.0],xstyle=1,ystyle=1,xrange = [19,23]

lw = 1e5
length = 1e19
length = 1e20
z = 0.0177
FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 4,color =240
;length = 1e20
z = 0.0177*0.3
FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 2,color =240
;length = 1e21
z = 0.0177*0.1
z = 0.0177*0.01
FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 1,color =240
;length = 1e22
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 3,color =240

lw = 1e7
length = 1e19
length = 1e20
z = 0.0177
FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 4,color =190
;length = 1e20
z = 0.0177*0.3
FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 2,color =190
;length = 1e21
z = 0.0177*0.1
z = 0.0177*0.01
FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 1,color =190
;length = 1e22
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 3,color =240

lw = 1e9
length = 1e19
length = 1e20
z = 0.0177
FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 4,color =150
;length = 1e20
z = 0.0177*0.3
FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 2,color =150
;length = 1e21
z = 0.0177*0.1
z = 0.0177*0.01
FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 1,color =150
;length = 1e22
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 3,color =220

;length = 1e20
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = lte_H2(column_den[i]/length,z,length,1e11,total_H)
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 2,color =200
;length = 1e21
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = lte_H2(column_den[i]/length,z,length,1e11,total_H)
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 1,color =200
;length = 1e22
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = lte_H2(column_den[i]/length,z,length,1e11,total_H)
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 3,color =200

;lw = 1e13
;length = 1e19
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 4,color =180
;length = 1e20
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 2,color =180
;length = 1e21
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 1,color =180
;length = 1e22
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 3,color =180

;length = 1e20
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;lte_H2[i] = lte_H2(column_den[i]/length,z,length,1e15,total_H)
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 2,color =160
;length = 1e21
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = lte_H2(column_den[i]/length,z,length,1e15,total_H)
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 1,color =160
;length = 1e22
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = lte_H2(column_den[i]/length,z,length,1e15,total_H)
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 3,color =160

;lw = 1e17
;length = 1e19
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 4,color =140
;length = 1e20
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 2,color =140
;length = 1e21
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 1,color =140
;length = 1e22
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 3,color =140

;lw = 1e21
;length = 1e19
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 4,color =100
;length = 1e20
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 2,color =100
;length = 1e21
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 1,color =100
;length = 1e22
;FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
;  lte_H2[i] = (lte_H2(column_den[i]/length,z,length,lw,total_H))[0]
;oplot,alog10(column_Den),lte_H2/total_H*2.0,linestyle = 3,color =100
if (KEYWORD_SET(outfile)) then  device,/close  else stop

END
