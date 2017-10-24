pro format_grids
;This procedure takes in isochrone files and returns grids for each
;color that have metallicity as columns and age as rows


; Leo Girardi's Isochrones -----------------------------------------------------
;filenames = ['color0004.dat','color0005.dat','color0007.dat','color0010.dat','color0013.dat','color0017.dat','color0023.dat','color0030.dat','color0040.dat','color0054.dat','color0073.dat','color0098.dat','color0130.dat','color0234.dat','color0300.dat']
;filebase = 'data_leo'
;n_age = 42

;Starburst99 Isochrones --------------------------------------------------------
filenames = ['color0004.dat','color004.dat','color008.dat','color02.dat','color05.dat']
filebase = '/astro/users/christensen/code/Database/data_star99_MS'
;n_age = 43
n_age = 828

cd,filebase
nfiles = N_ELEMENTS(filenames)
print,'Number of Metallicities: ',nfiles
print,'Number of Ages:              ',n_age

V_grid = fltarr(nfiles,n_age)
bv_grid = fltarr(nfiles,n_age)
ri_grid = fltarr(nfiles,n_age)
rprime_from_bv_grid = fltarr(nfiles,n_age)
rprime_from_vr_grid = fltarr(nfiles,n_age)
ub_grid = fltarr(nfiles,n_age)
vk_grid = fltarr(nfiles,n_age)
vr_grid = fltarr(nfiles,n_age)
Mb_grid = fltarr(nfiles,n_age)

FOR ct = 0, nfiles-1 DO BEGIN
;    readcol, filenames[ct], metallicity, age, mbol, ux, bx, b, v, r,i, j,  h, k, l, lprime, m, format='d,d,d,d,d,d,d,d,d,d,d,d,d,d,d'    ;leo
;    v_grid[ct,*] = v
;    bv_grid[ct,*] = b - v ;leo
;    ri_grid[ct,*] = r - i;leo
;    rprime_from_bv_grid[ct,*] = r;leo
;    rprime_from_vr_grid[ct,*] = r ; leo
;    ub_grid[ct,*] = ux - b ; leo
;    vk_grid[ct,*] = v - k ; leo
;    vr_grid[ct,*] = v - r ; leo

    readcol, filenames[ct],age,v130,v210,ub,bv,vr,vi,vj,vh,vk,vl,v,b,mb,format='d,d,d,d,d,d,d,d,d,d,d,d,d,d' ;starburst99
    v_grid[ct,*] = v
    bv_grid[ct,*] = bv ;starburst99
    ri_grid[ct,*] = vi - vr ;starburst99
    rprime_from_bv_grid[ct,*] = v - vr ;starburst99
    rprime_from_vr_grid[ct,*] = v - vr ;starburst99
    ub_grid[ct,*] = ub ; starburst99
    vk_grid[ct,*] = vk ; starburst99
    vr_grid[ct,*] = vr ; starburst99
    Mb_grid[ct,*] = mb ; starburst99
ENDFOR
openw,1,'ages'
printf,1,TRANSPOSE(age)
close,1

openw,1,'V_grid'
printf,1,v_grid
;printf,1,v_grid,format = '(d,d,d,d,d,d,d,d,d,d,d,d,d,d,d)' ;leo
close,1

openw,1,'bv_grid'
printf,1,bv_grid
;printf,1,bv_grid,format = '(d,d,d,d,d,d,d,d,d,d,d,d,d,d,d)' ;leo
close,1

openw,1,'ri_grid'
printf,1,ri_grid
;printf,1,ri_grid,format = '(d,d,d,d,d,d,d,d,d,d,d,d,d,d,d)' ;leo
close,1

openw,1,'rprime_from_bv_grid'
printf,1,rprime_from_bv_grid
;printf,1,rprime_from_bv_grid,format = '(d,d,d,d,d,d,d,d,d,d,d,d,d,d,d)' ;leo
close,1

openw,1,'rprime_from_vr_grid'
printf,1,rprime_from_vr_grid
;printf,1,rprime_from_vr_grid,format = '(d,d,d,d,d,d,d,d,d,d,d,d,d,d,d)' ;leo
close,1

openw,1,'ub_grid'
printf,1,ub_grid
;printf,1,ub_grid,format = '(d,d,d,d,d,d,d,d,d,d,d,d,d,d,d)' ;leo
close,1

openw,1,'vk_grid'
printf,1,vk_grid
;printf,1,vk_grid,format = '(d,d,d,d,d,d,d,d,d,d,d,d,d,d,d)' ;leo
close,1

openw,1,'vr_grid'
printf,1,vr_grid
;printf,1,vr_grid,format = '(d,d,d,d,d,d,d,d,d,d,d,d,d,d,d)' ;leo
close,1

openw,1,'mb_grid'
printf,1,mb_grid
;printf,1,mb_grid,format = '(d,d,d,d,d,d,d,d,d,d,d,d,d,d,d)' ;leo
close,1

cd,'..'

END
