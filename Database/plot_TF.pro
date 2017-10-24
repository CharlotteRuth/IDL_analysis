;This is the massive plotting program to make all sort so Tully Fisher graphs
;plot_TF,'leoTF'

pro plot_TF, outputfile
;msol = 2.362e5
msol = 2.0e17 ; mass of Sum in grams /1e16
dMsolUnit = 1.69875e16 ; Solar mass in system units
kpc = 3.085 ; km per kpc /1e16
grav = 6.67e-23
dKpcUnit = 50000.0 ;Kpc in system units
SYS_AMUCC = 0.0096302124  ;  conversion between system units and amu/c
hubble = 0.70
loadct,39

filename1 = 'Data_files/X2X2randv.dat'
filename2 = 'Data_files/X3X3randv.dat'
filename3 = 'Data_files/X5X5randv.dat'
filename4 = 'Data_files/X2X2nfrandv.dat'

readcol,filename1,indicies1,r2001, v2001, s2001, bm2001, r11, v11, s11, bm11, r21, v21, s21, bm21,format = 'A,F,F,F,F,F,F,F,F,F,F,F,F',/silent
readcol,filename2,indicies2,r2002, v2002, s2002, bm2002, r12, v12, s12, bm12, r22, v22, s22, bm22,format = 'A,F,F,F,F,F,F,F,F,F,F,F,F',/silent
readcol,filename3,indicies3,r2003, v2003, s2003, bm2003, r13, v13, s13, bm13, r23, v23, s23, bm23,format = 'A,F,F,F,F,F,F,F,F,F,F,F,F',/silent
readcol,filename4,indicies4,r2004, v2004, s2004, bm2004, r14, v14, s14, bm14, r24, v24, s24, bm24,format = 'A,F,F,F,F,F,F,F,F,F,F,F,F',/silent

list1 = 'Data_files/halo_listX2X2.dat'
readcol, list1, format = 'A60', files1
mag_file1 = '../Fabio_Runs/h1201-h1755-X2X2g1bwK/h1201-h1755-X2X2g1bwK.00512.amiga.halos.Mv.leo.fits'
print,"Reading ",mag_file1
halo_mags1 = MRDFITS(mag_file1,1)
num_halos1 = N_ELEMENTS(halo_mags1)
indices1 = STRTRIM(sindgen(num_halos1),2)
imag1 = fltarr(n_elements(files1))
kmag1 = fltarr(n_elements(files1))
id1 = fltarr(n_elements(files1))

for ct = 0, n_elements(files1)-1 do begin
;    print,files1[ct]
    halo_num = (strsplit(files1[ct],'.',/extract))[4]
    IF (strpos(halo_num,'0') eq 0) THEN halo_num =  STRMID(halo_num,1,STRLEN(halo_num)-1)
;    print,halo_num
;    IF (strpos(halo_num,'0') eq 0) THEN halo_num =  STRMID(halo_num,1,STRLEN(halo_num)-1)
;    print,halo_num
;    IF (strpos(halo_num,'0') eq 0) THEN halo_num =  STRMID(halo_num,1,STRLEN(halo_num)-1)
    list = strcmp(halo_num,indices1)
    halo_index = where(list eq 1)
    if(halo_index[0] ne -1)then  begin
        imag1[ct] = halo_mags1[halo_index].i 
        kmag1[ct] = halo_mags1[halo_index].k
        id1[ct] = halo_mags1[halo_index].id 
;        print,halo_num,halo_index,halo_mags1[halo_index].id        
    endif else begin
        imag1[ct] = -1
        kmag1[ct] = -1
        id1[ct] = -1
        print,halo_num,halo_index
        stop
    endelse
endfor

list2 = 'Data_files/halo_listX3X3.dat'
readcol, list2, format = 'A60', files2
mag_file2 = '../Fabio_Runs/h1201-h1755-X3X3g1bwK/h1201-h1755-X3X3g1bwK.00512.amiga.halos.Mv.leo.fits'
print,"Reading ",mag_file2
halo_mags2 = MRDFITS(mag_file2,1)
num_halos2 = N_ELEMENTS(halo_mags2)
indices2 = STRTRIM(sindgen(num_halos2),2)
imag2 = fltarr(n_elements(files2))
kmag2 = fltarr(n_elements(files2))
id2 = fltarr(n_elements(files2))

for ct = 0, n_elements(files2)-1 do begin
;    print,files2[ct]
    halo_num = (strsplit(files2[ct],'.',/extract))[4]
    IF (strpos(halo_num,'0') eq 0) THEN halo_num = STRMID(halo_num,1,STRLEN(halo_num)-1)
    list = strcmp(halo_num,indices2)
    halo_index = where(list eq 1)
    if(halo_index[0] ne -1)then begin
        imag2[ct] = halo_mags2[halo_index].i 
        kmag2[ct] = halo_mags2[halo_index].k
        id2[ct] = halo_mags2[halo_index].id 
;        print,halo_num,halo_index,halo_mags2[halo_index].id          
    endif else begin
        imag2[ct] = -1
        kmag2[ct] = -1
        id2[ct] = -1     
        print,halo_num,halo_index   
        stop
    endelse
endfor

list3 = 'Data_files/halo_listX5X5.dat'
readcol, list3, format = 'A60', files3
mag_file3 = '../Fabio_Runs/h1201-h1755-X5X5g1bwK/h1201-h1755-X5X5g1bwK.00512.amiga.halos.Mv.leo.fits'
print,"Reading ",mag_file3
halo_mags3 = MRDFITS(mag_file3,1)
num_halos3 = N_ELEMENTS(halo_mags3)
indices3 = STRTRIM(sindgen(num_halos3),2)
imag3 = fltarr(n_elements(files3))
kmag3 = fltarr(n_elements(files3))
id3 = fltarr(n_elements(files3))

for ct = 0, n_elements(files3)-1 do begin
;   print,files3[ct]
    halo_num = (strsplit(files3[ct],'.',/extract))[4]
    IF (strpos(halo_num,'0') eq 0) THEN halo_num =   STRMID(halo_num,1,STRLEN(halo_num)-1)
    list = strcmp(halo_num,indices3)
    halo_index = where(list eq 1)
    if(halo_index[0] ne -1)then begin    
        imag3[ct] = halo_mags3[halo_index].i
        kmag3[ct] = halo_mags3[halo_index].k
        id3[ct] = halo_mags3[halo_index].id 
;        print,halo_num,halo_index,halo_mags3[halo_index].id   
    endif else begin
        imag3[ct] = -1
        kmag3[ct] = -1
        id3[ct] = -1 
        print,halo_num,halo_index   
        stop
    endelse         
endfor

list4 = 'Data_files/halo_listX2X2nf.dat'
readcol, list4, format = 'A60', files4
mag_file4 = '../Fabio_Runs/h1201-h1755-X2X2g3bwK/h1201-h1755-X2X2g3bwK.00512.amiga.halos.Mv.leo.fits'
print,"Reading ",mag_file4
halo_mags4 = MRDFITS(mag_file4,1)
num_halos4 = N_ELEMENTS(halo_mags4)
indices4 = STRTRIM(sindgen(num_halos4),2)
imag4 = fltarr(n_elements(files4))
kmag4 = fltarr(n_elements(files4))
id4 = fltarr(n_elements(files4))

for ct = 0, n_elements(files4)-1 do begin
;   print,files4[ct]
    halo_num = (strsplit(files4[ct],'.',/extract))[4]
    IF (strpos(halo_num,'0') eq 0) THEN halo_num = STRMID(halo_num,1,STRLEN(halo_num)-1)
    list = strcmp(halo_num,indices4)
    halo_index = where(list eq 1)
    if(halo_index[0] ne -1)then begin  
        imag4[ct] = halo_mags4[halo_index].i
        kmag4[ct] = halo_mags4[halo_index].k
        id4[ct] = halo_mags4[halo_index].id
;        print,halo_num,halo_index,halo_mags4[halo_index].id 
    endif else begin
        imag4[ct] = -1
        kmag4[ct] = -1
        id4[ct] = -1
        print,halo_num,halo_index 
        stop
    endelse     
endfor

close,1
openw,1,'X2X2_data.dat'
printf,1,'ID ','v200 ','v20 ','v10 ','Ileo ','Kleo '
printf,1,TRANSPOSE([[id1],[v2001],[v11],[v21],[imag1],[kmag1]])
close,1

close,1
openw,1,'X3X3_data.dat'
printf,1,'ID ','v200 ','v20 ','v10 ','Ileo ','Kleo '
printf,1,TRANSPOSE([[id2],[v2002],[v12],[v22],[imag2],[kmag2]])
close,1

close,1
openw,1,'X5X5_data.dat'
printf,1,'ID ','v200 ','v20 ','v10 ','Ileo ','Kleo '
printf,1,TRANSPOSE([[id3],[v2003],[v13],[v23],[imag3],[kmag3]])
close,1

openw,1,'X2X2_data_nf.dat'
printf,1,'ID ','v200 ','v20 ','v10 ','Ileo ','Kleo '
printf,1,TRANSPOSE([[id4],[v2004],[v14],[v24],[imag4],[kmag4]])
close,1
;stop


!Y.style = 1
!X.style = 1
mags = findgen((24.0-8.0)/0.1)*0.1 -24.0
slope = -0.128412
velocities = 10.0^(slope*(mags + 15.5) + 1.5)
;Fit given in Eke, Navarro and Steimetz 2001


;*************************Comparing Magnitudes*************************************************************
set_plot,'x'
plot,[1e8,1e12],[0,0],xrange = [1e8,1e12],yrange=[-5,0.5],xtitle="M200 -- Baryonic", ytitle = "M/L = 2.5 - Leo I Mags",title = "Magnitude Difference -- Different Radii",/xlog
oplot,bm2001,imag1 - 4.08+2.5*ALOG10(s11/1.5),psym = 1, color = 250
oplot,bm2001,imag1 - 4.08+2.5*ALOG10(s21/1.5),psym = 1, color = 150
oplot,bm2001,imag1 - 4.08+2.5*ALOG10(s2001/1.5),psym = 1, color = 50
legend,['R200','20/20*v200/220','10/220*v200/h','Goal'],color=[250,150,50,0],psym = [1,1,1,-3],linestyle=[0,0,0,0],/bottom,/left
set_plot,'ps'
device,filename='comp_mags_leo_ML1.eps',/color,bits_per_pixel=8
plot,[1e8,1e12],[0,0],xrange = [1e8,1e12],ytitle = "M/L = 1.5 - Leo I Mags", xtitle="M200 -- Baryonic",title = "Magnitude Difference -- Different Res",/xlog,yrange=[-5,0.5]
oplot,bm2001,imag1 - 4.08+2.5*ALOG10(s11/1.5),psym = 1, color = 250
oplot,bm2001,imag1 - 4.08+2.5*ALOG10(s21/1.5),psym = 1, color = 150
oplot,bm2001,imag1 - 4.08+2.5*ALOG10(s2001/1.5),psym = 1, color = 50
legend,['R200','20/20*v200/220','10/220*v200/h','Goal'],color=[250,150,50,0],psym = [1,1,1,-3],linestyle=[0,0,0,0],/bottom,/left
device,/close
stop

set_plot,'x'
plot,[1e8,1e12],[0,0], xrange = [1e8,1e12],ytitle = "M/L = 1.5 - Leo I Mags", xtitle = "M200 -- Baryonic",title = "Magnitude Difference -- Different Res",/xlog,yrange=[-5,0.5]
oplot,bm2001,imag1 - 4.08+2.5*ALOG10(s11/1.5),psym = 1, color = 50
oplot,bm2004,imag4 - 4.08+2.5*ALOG10(s14/1.5),psym = 1, color = 100
oplot,bm2002,imag2 - 4.08+2.5*ALOG10(s12/1.5),psym = 1, color = 200
oplot,bm2003,imag3 - 4.08+2.5*ALOG10(s13/1.5),psym = 1, color = 250
legend,['X2X2','X2X2nf','X3X3','X5X5','Goal'],color=[50,100,200,250,0],psym = [1,1,1,1,-3],linestyle=[0,0,0,0,0],/bottom,/left
set_plot,'ps'
device,filename='comp_mags_leo_ML2.eps',/color,bits_per_pixel=8
plot,[1e8,1e12],[0,0], xrange = [1e8,1e12],ytitle = "M/L = 1.5 - Leo I Mags", xtitle = "M200 -- Baryonic",title = "Magnitude Difference -- Different Res",/xlog,yrange=[-5,0.5]
oplot,bm2001,imag1 - 4.08+2.5*ALOG10(s11/1.5),psym = 1, color = 50
oplot,bm2004,imag4 - 4.08+2.5*ALOG10(s14/1.5),psym = 1, color = 100
oplot,bm2002,imag2 - 4.08+2.5*ALOG10(s12/1.5),psym = 1, color = 200
oplot,bm2003,imag3 - 4.08+2.5*ALOG10(s13/1.5) ,psym = 1, color = 250
legend,['X2X2','X2X2nf','X3X3','X5X5','Goal'],color=[50,100,200,250,0],psym = [1,1,1,1,-3],linestyle=[0,0,0,0,0],/bottom,/left
device,/close
stop




;*************************X2X2*******************************************************************
set_plot,'x'
plot,velocities,mags,title = 'Tully Fisher Relation - X2X2 -- Girardi',ytitle = 'MI - 5 Log(h)', xtitle = 'Vcirc',xrange=[100,350],yrange = [-18,-23],thick = 2,/xlog;yrange = [-12,-26]
;oplot,imag1 - 5*ALOG10(hubble),ALOG10(v2001),psym = 4,color = 50
oplot,v11,imag1 - 5*ALOG10(hubble),psym = 4,color = 50
oplot,v21,imag1 - 5*ALOG10(hubble),psym = 4,color = 240
legend,['20/220*v200/h','10/220*v200/h','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3],linestyle=[0,0,0]
set_plot, 'ps'
device, filename=outputfile+'X2X2.eps',/color,bits_per_pixel=8
plot,velocities,mags,title = 'Tully Fisher Relation - X2X2 -- Girardi',ytitle = 'MI - 5 Log(h)', xtitle = 'Log(Vcirc)',yrange = [-18,-23],thick = 2,/xlog,xrange=[100,350];xrange=[1.4,2.6],yrange = [-12,-26]
;oplot,imag1 - 5*ALOG10(hubble),ALOG10(v2001),psym = 4,color = 50
oplot,v11,imag1 - 5*ALOG10(hubble),psym = 4,color = 50
oplot,v21,imag1 - 5*ALOG10(hubble),psym = 4,color = 240
legend,['20/220*v200/h','10/220*v200/h','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3],linestyle=[0,0,0]
device,/close
stop

set_plot,'x'
plot,velocities,mags,title = 'Tully Fisher Relation - X2X2 -- M/L = 1.5',ytitle = 'MI - 5 Log(h)', xtitle = 'Vcirc',xrange=[100,350],yrange = [-18,-23],thick = 2,/xlog;yrange = [-12,-26]
;oplot,imag1 - 5*ALOG10(hubble),ALOG10(v2001),psym = 4,color = 50
oplot,v11,4.08-2.5*ALOG10(s11/1.5) - 5*ALOG10(hubble),psym = 4,color = 50
oplot,v21,4.08-2.5*ALOG10(s21/1.5) - 5*ALOG10(hubble),psym = 4,color = 240
legend,['20/220*v200/h','10/220*v200/h','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3],linestyle=[0,0,0]
set_plot, 'ps'
device, filename=outputfile+'X2X2stellar.eps',/color,bits_per_pixel=8
plot,velocities,mags,title = 'Tully Fisher Relation - X2X2 -- M/L = 1.5',ytitle = 'MI - 5 Log(h)', xtitle = 'Log(Vcirc)',yrange = [-18,-23],thick = 2,/xlog,xrange=[100,350];xrange=[1.4,2.6],yrange = [-12,-26]
;oplot,imag1 - 5*ALOG10(hubble),ALOG10(v2001),psym = 4,color = 50
oplot,v11,4.08-2.5*ALOG10(s11/1.5) - 5*ALOG10(hubble),psym = 4,color = 50
oplot,v21,4.08-2.5*ALOG10(s21/1.5) - 5*ALOG10(hubble),psym = 4,color = 240
legend,['20/220*v200/h','10/220*v200/h','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3],linestyle=[0,0,0]
device,/close
stop

set_plot,'x'
plot,[10,400],[7e5,1e12],title = 'Baryonic TF - X2X2',xtitle = 'Log(Vcirc)', ytitle = 'Baryonic Mass',/ylog,/xlog,xrange=[100,350],yrange = [1e9,1e12]
oplot,v11,bm11,psym = 4,color = 50
oplot,v21,bm21,psym = 4,color = 240
legend,['20/220*v200/h','10/220*v200/h','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3]
set_plot, 'ps'
device, filename='baryonicTFX2X2.eps',/color,bits_per_pixel=8
plot,[10,400],[7e5,1e12],title = 'Baryonic TF - X2X2',xtitle = 'Log(Vcirc)', ytitle = 'Baryonic Mass',/ylog,/xlog,xrange=[100,350],yrange = [1e9,1e12]
oplot,v11,bm11,psym = 4,color = 50
oplot,v21,bm21,psym = 4,color = 240
legend,['20/220*v200/h','10/220*v200/h','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3]
device,/close
stop

;***************X3X3***************************************************************************
set_plot, 'x'
plot,velocities,mags,title = 'Tully Fisher Relation - X3X3 -- Girardi',ytitle = 'MI - 5 Log(h)', xtitle = 'Log(Vcirc)',yrange = [-18,-23],thick = 2,/xlog,xrange=[100,350];xrange=[1.8,2.6],,yrange = [-18,-24]
;oplot,imag2 - 5*ALOG10(hubble),ALOG10(v2002),psym = 4,color = 50
oplot,v12,imag2 - 5*ALOG10(hubble),psym = 4,color = 50
oplot,v22,imag2 - 5*ALOG10(hubble),psym = 4,color = 240
legend,['20/220*v200/h','10/220*v200/h','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3],linestyle=[0,0,0]
set_plot, 'ps'
device, filename=outputfile+'X3X3.eps',/color,bits_per_pixel=8
plot,velocities,mags,title = 'Tully Fisher Relation - X3X3 -- Girardi',ytitle = 'MI - 5 Log(h)', xtitle = 'Log(Vcirc)',yrange = [-18,-23],thick = 2,/xlog,xrange=[100,350];xrange=[1.8,2.6],yrange = [-18,-24]
;oplot,imag2 - 5*ALOG10(hubble),ALOG10(v2002),psym = 4,color = 50
oplot,v12,imag2 - 5*ALOG10(hubble),psym = 4,color = 50
oplot,v22,imag2 - 5*ALOG10(hubble),psym = 4,color = 240
legend,['20/220*v200/h','10/220*v200/h','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3],linestyle=[0,0,0]
device,/close

stop
set_plot,'x'
plot,velocities,mags,title = 'Tully Fisher Relation - X3X3 -- M/L = 1.5',ytitle = 'MI - 5 Log(h)', xtitle = 'Vcirc',xrange=[100,350],yrange = [-18,-23],thick = 2,/xlog;yrange = [-12,-26]
;oplot,imag1 - 5*ALOG10(hubble),ALOG10(v2001),psym = 4,color = 50
oplot,v12,4.08-2.5*ALOG10(s12/1.5) - 5*ALOG10(hubble),psym = 4,color = 50
oplot,v22,4.08-2.5*ALOG10(s22/1.5) - 5*ALOG10(hubble),psym = 4,color = 240
legend,['20/220*v200/h','10/220*v200/h','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3],linestyle=[0,0,0]
set_plot, 'ps'
device, filename=outputfile+'X3X3stellar.eps',/color,bits_per_pixel=8
plot,velocities,mags,title = 'Tully Fisher Relation - X2X2 -- M/L = 1.5',ytitle = 'MI - 5 Log(h)', xtitle = 'Log(Vcirc)',yrange = [-18,-23],thick = 2,/xlog,xrange=[100,350];xrange=[1.4,2.6],yrange = [-12,-26]
;oplot,imag1 - 5*ALOG10(hubble),ALOG10(v2001),psym = 4,color = 50
oplot,v12,4.08-2.5*ALOG10(s12/1.5) - 5*ALOG10(hubble),psym = 4,color = 50
oplot,v22,4.08-2.5*ALOG10(s22/1.5) - 5*ALOG10(hubble),psym = 4,color = 240
legend,['20/220*v200/h','10/220*v200/h','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3],linestyle=[0,0,0]
device,/close
stop

set_plot,'x'
plot,[10,400],[7e5,1e12],title = 'Baryonic TF - X3X3',xtitle = 'Log(Vcirc)', ytitle = 'Baryonic Mass',/ylog,/xlog,xrange=[100,350],yrange = [1e9,1e12]
oplot,v12,bm12,psym = 4,color = 50
oplot,v22,bm22,psym = 4,color = 240
legend,['20/220*v200/h','10/220*v200/h','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3]
set_plot, 'ps'
device, filename='baryonicTFX3X3.eps',/color,bits_per_pixel=8
plot,[10,400],[7e5,1e12],title = 'Baryonic TF - X3X3',xtitle = 'Log(Vcirc)', ytitle = 'Baryonic Mass',/ylog,/xlog,xrange=[100,350],yrange = [1e9,1e12]
oplot,v12,bm12,psym = 4,color = 50
oplot,v22,bm22,psym = 4,color = 240
legend,['20/220*v200/h','10/220*v200/h','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3]
device,/close
stop

;************************X5X5*********************************************************************
set_plot, 'x'
plot,velocities,mags,title = 'Tully Fisher Relation - X5X5 -- Girardi',ytitle = 'MI - 5 Log(h)', xtitle = 'Log(Vcirc)',yrange = [-8,-23],thick = 2,/xlog,xrange=[10,350];,xrange=[1.4,2.4],,yrange = [-8,-24]
;oplot,imag2 - 5*ALOG10(hubble),ALOG10(v2003),psym = 4,color = 50
oplot,v13,imag3 - 5*ALOG10(hubble),psym = 4,color = 50
oplot,v23,imag3 - 5*ALOG10(hubble),psym = 4,color = 240
legend,['20/220*v200/h','10/220*v200/h','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3],linestyle=[0,0,0]
set_plot, 'ps'
device, filename=outputfile+'X5X5.eps',/color,bits_per_pixel=8
plot,velocities,mags,title = 'Tully Fisher Relation - X5X5 -- Girardi',ytitle = 'MI - 5 Log(h)', xtitle = 'Log(Vcirc)',yrange = [-8,-23],thick = 2,/xlog,xrange=[10,350];,xrange=[1.4,2.4],yrange = [-8,-24]
;oplot,imag3 - 5*ALOG10(hubble),ALOG10(v2003),psym = 4,color = 50
oplot,v13,imag3 - 5*ALOG10(hubble),psym = 4,color = 50
oplot,v23,imag3 - 5*ALOG10(hubble),psym = 4,color = 240
legend,['20/220*v200/h','10/220*v200/h','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3],linestyle=[0,0,0]
device,/close
set_plot, 'x'
stop

set_plot,'x'
plot,velocities,mags,title = 'Tully Fisher Relation - X5X5 -- M/L = 1.5',ytitle = 'MI - 5 Log(h)', xtitle = 'Vcirc',xrange=[10,350],yrange = [-8,-23],thick = 2,/xlog;yrange = [-12,-26]
;oplot,imag1 - 5*ALOG10(hubble),ALOG10(v2001),psym = 4,color = 50
oplot,v13,4.08-2.5*ALOG10(s13/1.5) - 5*ALOG10(hubble),psym = 4,color = 50
oplot,v23,4.08-2.5*ALOG10(s23/1.5) - 5*ALOG10(hubble),psym = 4,color = 240
legend,['20/220*v200/h','10/220*v200/h','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3],linestyle=[0,0,0]
set_plot, 'ps'
device, filename=outputfile+'X5X5stellar.eps',/color,bits_per_pixel=8
plot,velocities,mags,title = 'Tully Fisher Relation - X5X5 -- M/L = 1.5',ytitle = 'MI - 5 Log(h)', xtitle = 'Log(Vcirc)',yrange = [-8,-23],thick = 2,/xlog,xrange=[10,350];xrange=[1.4,2.6],yrange = [-12,-26]
;oplot,imag1 - 5*ALOG10(hubble),ALOG10(v2001),psym = 4,color = 50
oplot,v13,4.08-2.5*ALOG10(s13/1.5) - 5*ALOG10(hubble),psym = 4,color = 50
oplot,v23,4.08-2.5*ALOG10(s23/1.5) - 5*ALOG10(hubble),psym = 4,color = 240
legend,['20/220*v200/h','10/220*v200/h','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3],linestyle=[0,0,0]
device,/close
stop

set_plot,'x'
plot,[10,400],[7e5,1e12],title = 'Baryonic TF - X5X5',xtitle = 'Log(Vcirc)', ytitle = 'Baryonic Mass',/ylog,/xlog,xrange=[10,350],yrange = [7e5,1e12]
oplot,v13,bm13,psym = 4,color = 50
oplot,v23,bm23,psym = 4,color = 240
legend,['20/220*v200/h','10/220*v200/h','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3]
set_plot, 'ps'
device, filename='baryonicTFX5X5.eps',/color,bits_per_pixel=8
plot,[10,400],[7e5,1e12],title = 'Baryonic TF - X5X5',xtitle = 'Log(Vcirc)', ytitle = 'Baryonic Mass',/ylog,/xlog,xrange=[10,350],yrange = [7e5,1e12]
oplot,v13,bm13,psym = 4,color = 50
oplot,v23,bm23,psym = 4,color = 240
legend,['20/220*v200/h','10/220*v200/h','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3]
device,/close
stop

;**************************X2X2_NF**********************************************************
set_plot,'x'
plot,velocities,mags,title = 'No Feedback TF -- X2X2 -- Girardi',ytitle = 'MI - 5 Log(h)', xtitle = 'Log(Vcirc)',yrange = [-18 ,-23],thick = 2,/xlog,xrange=[100,350];,xrange=[1.4,2.6],yrange = [-8 ,-26]
;oplot,imag1 - 5*ALOG10(hubble),ALOG10(v2001),psym = 4,color = 50
oplot,v14,imag4 - 5*ALOG10(hubble),psym = 4,color = 50
oplot,v24,imag4 - 5*ALOG10(hubble),psym = 4,color = 240
legend,['20/220*v200/h','10/220*v200/h','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3],linestyle=[0,0,0]
set_plot, 'ps'
device, filename=outputfile+'X2X2nf.eps',/color,bits_per_pixel=8
plot,velocities,mags,title = 'No Feedback TF -- X2X2 -- Girardi',ytitle = 'MI - 5 Log(h)', xtitle = 'Log(Vcirc)',yrange = [-18 ,-23],thick = 2,/xlog,xrange=[100,350];,xrange=[1.4,2.6],yrange = [-8 ,-26]
;oplot,imag1 - 5*ALOG10(hubble),ALOG10(v2001),psym = 4,color = 50
oplot,v14,imag4 - 5*ALOG10(hubble),psym = 4,color = 50
oplot,v24,imag4 - 5*ALOG10(hubble),psym = 4,color = 240
legend,['20/220*v200/h','10/220*v200/h','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3],linestyle=[0,0,0]
device,/close
stop

set_plot,'x'
plot,velocities,mags,title = 'No Feedback TF - X2X2 -- M/L = 1.5',ytitle = 'MI - 5 Log(h)', xtitle = 'Vcirc',xrange=[100,350],yrange = [-18,-23],thick = 2,/xlog;yrange = [-12,-26]
;oplot,imag1 - 5*ALOG10(hubble),ALOG10(v2001),psym = 4,color = 50
oplot,v14,4.08-2.5*ALOG10(s14/1.5) - 5*ALOG10(hubble),psym = 4,color = 50
oplot,v24,4.08-2.5*ALOG10(s24/1.5) - 5*ALOG10(hubble),psym = 4,color = 240
legend,['20/220*v200/h','10/220*v200/h','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3],linestyle=[0,0,0]
set_plot, 'ps'
device, filename=outputfile+'X2X2nfstellar.eps',/color,bits_per_pixel=8
plot,velocities,mags,title = 'Tully Fisher Relation - X2X2nf -- M/L = 1.5',ytitle = 'MI - 5 Log(h)', xtitle = 'Log(Vcirc)',yrange = [-18,-23],thick = 2,/xlog,xrange=[100,350];xrange=[1.4,2.6],yrange = [-12,-26]
;oplot,imag1 - 5*ALOG10(hubble),ALOG10(v2001),psym = 4,color = 50
oplot,v14,4.08-2.5*ALOG10(s14/1.5) - 5*ALOG10(hubble),psym = 4,color = 50
oplot,v24,4.08-2.5*ALOG10(s24/1.5) - 5*ALOG10(hubble),psym = 4,color = 240
legend,['20/220*v200/h','10/220*v200/h','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3],linestyle=[0,0,0]
device,/close
stop

set_plot,'x'
plot,[10,400],[7e5,1e12],title = 'Baryonic TF - No Feedback X2X2',xtitle = 'Log(Vcirc)', ytitle = 'Baryonic Mass',/ylog,/xlog,xrange=[100,350],yrange = [1e9,1e12]
oplot,v14,bm14,psym = 4,color = 50
oplot,v24,bm24,psym = 4,color = 240
legend,['20/220*v200/h','10/220*v200/h','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3]
set_plot, 'ps'
device, filename='baryonicTFX2X2nf.eps',/color,bits_per_pixel=8
plot,[10,400],[7e5,1e12],title = 'Baryonic TF - No Feedback X2X2',xtitle = 'Log(Vcirc)', ytitle = 'Baryonic Mass',/ylog,/xlog,xrange=[100,350],yrange = [1e9,1e12]
oplot,v14,bm14,psym = 4,color = 50
oplot,v24,bm24,psym = 4,color = 240
legend,['20/220*v200/h','10/220*v200/h','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3]
device,/close
stop

;**********************Composit Feedback/No Feedback********************************
set_plot,'x'
plot,velocities,mags,title = 'X2X2 With/Without Feedback -- M/L = 1.5',ytitle = 'MI - 5 Log(h)', xtitle = 'Log(Vcirc)',yrange = [-18 ,-23],thick = 2,/xlog,xrange=[100,350]
oplot,v14,4.08-2.5*ALOG10(s14/1.5) - 5*ALOG10(hubble),psym = 4,color = 240
oplot,v11,4.08-2.5*ALOG10(s11/1.5) - 5*ALOG10(hubble),psym = 4,color = 50
legend,['Feedback','No Feedback','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3],linestyle=[0,0,0],/bottom,/righ
set_plot, 'ps'
device, filename='stellarTFX2X2f_nf.eps',/color,bits_per_pixel=8
plot,velocities,mags,title = 'X2X2 With/Without Feedback -- M/L = 1.5',ytitle = 'MI - 5 Log(h)', xtitle = 'Log(Vcirc)',yrange = [-18 ,-23],thick = 2,/xlog,xrange=[100,350]
oplot,v14,4.08-2.5*ALOG10(s14/1.5) - 5*ALOG10(hubble),psym = 4,color = 240
oplot,v11,4.08-2.5*ALOG10(s11/1.5) - 5*ALOG10(hubble),psym = 4,color = 50
legend,['Feedback -- 20/220*v200/h','No Feedback -- 20/220*v200/h','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3],linestyle=[0,0,0],/bottom,/righ
device,/close
stop

set_plot,'x'
plot,velocities,mags,title = 'X2X2 With/Without Feedback -- Girardi',ytitle = 'MI - 5 Log(h)', xtitle = 'Log(Vcirc)',yrange = [-18 ,-23],thick = 2,/xlog,xrange=[100,350]
oplot,v14,imag4 - 5*ALOG10(hubble),psym = 4,color = 240
oplot,v11,imag1 - 5*ALOG10(hubble),psym = 4,color = 50
legend,['Feedback','No Feedback','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3],linestyle=[0,0,0],/bottom,/righ
set_plot, 'ps'
device, filename=outputfile+'X2X2f_nf.eps',/color,bits_per_pixel=8
plot,velocities,mags,title = 'X2X2 With/Without Feedback -- Girardi',ytitle = 'MI - 5 Log(h)', xtitle = 'Log(Vcirc)',yrange = [-18 ,-23],thick = 2,/xlog,xrange=[100,350]
oplot,v14,imag4 - 5*ALOG10(hubble),psym = 4,color = 240
oplot,v11,imag1 - 5*ALOG10(hubble),psym = 4,color = 50
legend,['Feedback -- 20/220*v200/h','No Feedback -- 20/220*v200/h','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3],linestyle=[0,0,0],/bottom,/righ
device,/close
stop

set_plot,'x'
plot,[10,400],[7e5,1e12],title = 'Baryonic TF - With/Without Feedback X2X2',xtitle = 'Log(Vcirc)', ytitle = 'Baryonic Mass',/ylog,/xlog,xrange=[100,350],yrange = [1e9,1e12]
oplot,v14,bm14,psym = 4,color = 50
oplot,v11,bm11,psym = 4,color = 240
legend,['Feedback','No Feedback','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3],linestyle=[0,0,0],/bottom,/right
set_plot, 'ps'
device, filename='baryonicTFX2X2f_nf.eps',/color,bits_per_pixel=8
plot,[10,400],[7e5,1e12],title = 'Baryonic TF - With/Without Feedback X2X2',xtitle = 'Log(Vcirc)', ytitle = 'Baryonic Mass',/ylog,/xlog,xrange=[100,350],yrange = [1e9,1e12]
oplot,v14,bm14,psym = 4,color = 50
oplot,v11,bm11,psym = 4,color = 240
legend,['Feedback','No Feedback','Giovanelli et al'],color=[50,240,0],psym = [4,4,-3],linestyle=[0,0,0],/bottom,/right
device,/close
stop

;************************ Effect of Resolution *************************************

set_plot,'x'
plot,velocities,mags,title = 'Effect of Resolution -- Girardi',ytitle = 'MI - 5 Log(h)', xtitle = 'Log(Vcirc)',yrange = [-18 ,-23],thick = 2,/xlog,xrange=[100,350]
oplot,v11,imag1 - 5*ALOG10(hubble),psym = 4,color = 50
oplot,v12,imag2 - 5*ALOG10(hubble),psym = 4,color = 150
oplot,v13,imag3 - 5*ALOG10(hubble),psym = 4,color = 200
legend,['X2X2 -- 20/220*v200/h','X3X3 -- 20/220*v200/h','X5X5 -- 20/220*v200/h','Giovanelli et al'],color=[50,150,200,0],psym = [4,4,4,-3],linestyle=[0,0,0,0],/bottom,/righ
set_plot,'ps'
device, filename=outputfile+'res.eps',/color,bits_per_pixel=8
plot,velocities,mags,title = 'Effect of Resolution -- Girardi',ytitle = 'MI - 5 Log(h)', xtitle = 'Log(Vcirc)',yrange = [-18 ,-23],thick = 2,/xlog,xrange=[100,350]
oplot,v11,imag1 - 5*ALOG10(hubble),psym = 4,color = 50
oplot,v12,imag2 - 5*ALOG10(hubble),psym = 4,color = 150
oplot,v13,imag3 - 5*ALOG10(hubble),psym = 4,color = 200
legend,['X2X2 -- 20/220*v200/h','X3X3 -- 20/220*v200/h','X5X5 -- 20/220*v200/h','Giovanelli et al'],color=[50,150,200,0],psym = [4,4,4,-3],linestyle=[0,0,0,0],/bottom,/right
device,/close
stop

set_plot,'x'
plot,velocities,mags,title = 'Effect of Resolution -- M/L = 1.5',ytitle = 'MI - 5 Log(h)', xtitle = 'Log(Vcirc)',yrange = [-18 ,-23],thick = 2,/xlog,xrange=[100,350]
oplot,v11,4.08-2.5*ALOG10(s11/1.5) - 5*ALOG10(hubble),psym = 4,color = 50
oplot,v12,4.08-2.5*ALOG10(s12/1.5) - 5*ALOG10(hubble),psym = 4,color = 150
oplot,v13,4.08-2.5*ALOG10(s13/1.5) - 5*ALOG10(hubble),psym = 4,color = 200
legend,['X2X2 -- 20/220*v200/h','X3X3 -- 20/220*v200/h','X5X5 -- 20/220*v200/h','Giovanelli et al'],color=[50,150,200,0],psym = [4,4,4,-3],linestyle=[0,0,0,0],/bottom,/right
set_plot,'ps'
device, filename='stellarTFres.eps',/color,bits_per_pixel=8
plot,velocities,mags,title = 'Effect of Resolution -- M/L = 1.5',ytitle = 'MI - 5 Log(h)', xtitle = 'Log(Vcirc)',yrange = [-18 ,-23],thick = 2,/xlog,xrange=[100,350]
oplot,v11,4.08-2.5*ALOG10(s11/1.5) - 5*ALOG10(hubble),psym = 4,color = 50
oplot,v12,4.08-2.5*ALOG10(s12/1.5) - 5*ALOG10(hubble),psym = 4,color = 150
oplot,v13,4.08-2.5*ALOG10(s13/1.5) - 5*ALOG10(hubble),psym = 4,color = 200
legend,['X2X2 -- 20/220*v200/h','X3X3 -- 20/220*v200/h','X5X5 -- 20/220*v200/h','Giovanelli et al'],color=[50,150,200,0],psym = [4,4,4,-3],linestyle=[0,0,0,0],/bottom,/right
device,/close
stop

set_plot,'x'
plot,[10,400],[7e5,1e12],title = 'Effect of Resolution -- Baryonic TF',xtitle = 'Log(Vcirc)', ytitle = 'Baryonic Mass',/ylog,/xlog,xrange=[100,350],yrange = [1e9,1e12]
oplot,v11,bm11,psym = 4,color = 50
oplot,v12,bm12,psym = 4,color = 150
oplot,v13,bm13,psym = 4,color = 200
legend,['X2X2 -- 20/220*v200/h','X3X3 -- 20/220*v200/h','X5X5 -- 20/220*v200/h','Giovanelli et al'],color=[50,150,200,0],psym = [4,4,4,-3],linestyle=[0,0,0,0],/bottom,/right
set_plot,'ps'
device, filename='baryonicTFres.eps',/color,bits_per_pixel=8
plot,[10,400],[7e5,1e12],title = 'Effect of Resolution -- Baryonic TF',xtitle = 'Log(Vcirc)', ytitle = 'Baryonic Mass',/ylog,/xlog,xrange=[100,350],yrange = [1e9,1e12]
oplot,v11,bm11,psym = 4,color = 50
oplot,v12,bm12,psym = 4,color = 150
oplot,v13,bm13,psym = 4,color = 200
legend,['X2X2 -- 20/220*v200/h','X3X3 -- 20/220*v200/h','X5X5 -- 20/220*v200/h','Giovanelli et al'],color=[50,150,200,0],psym = [4,4,4,-3],linestyle=[0,0,0,0],/bottom,/right
device,/close
stop

end
