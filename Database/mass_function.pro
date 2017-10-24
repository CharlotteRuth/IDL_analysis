;This compares the stellar mass function for the X2X2 and X5X5 runs
;mass_func,'1'

PRO mass_func, outputfile
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
filename2 = 'Data_files/X5X5randv.dat'
;filename3 = 'Data_files/X3X3randv.dat'
;filename4 = 'Data_files/X2X2nfrandv.dat'

readcol,filename1,indicies1,r2001, v2001, s2001, bm2001, r11, v11, s11, bm11, r21, v21, s21, bm21,format = 'A,F,F,F,F,F,F,F,F,F,F,F,F',/silent
readcol,filename2,indicies2,r2002, v2002, s2002, bm2002, r12, v12, s12, bm12, r22, v22, s22, bm22,format = 'A,F,F,F,F,F,F,F,F,F,F,F,F',/silent
;readcol,filename3,indicies3,r2003, v2003, s2003, bm2003, r13, v13, s13, bm13, r23, v23, s23, bm23,format = 'A,F,F,F,F,F,F,F,F,F,F,F,F',/silent
;readcol,filename4,indicies4,r2004, v2004, s2004, bm2004, r14, v14, s14, bm14, r24, v24, s24, bm24,format = 'A,F,F,F,F,F,F,F,F,F,F,F,F',/silent



END
