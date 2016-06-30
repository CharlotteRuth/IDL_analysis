;.r /astro/users/christensen/code/IDL/procedures/wtipsy.pro
;.r /astro/users/christensen/code/IDL/idl4tipsy/tipsySatsHI.pro
;.r /astro/users/christensen/code/IDL/idl4tipsy/alignHI.pro

;outarray =['FeMassFracdot','Metalsdot','OxMassFracdot','eCool','eCoolH2','eDot','eHeat','HI','H2','OxMassFrac']
;outarray =['amiga.grp','FeMassFrac','OxMassFrac','HI','H2','correL','coolontime','iord','lw','HeI','HeII']
;outarray = ['amiga.grp']
;tipsySatsHI,'h516.cosmo25cmb.3072g14HBWK'      ,1,25000,2.310e15,  /cutout_rad,outarray= outarray
;tipsySatsHI,'h799.cosmo25cmb.3072g14HBWK.00084',2,25000,2.310e15,  /cutout_rad,outarray= outarray
;tipsySatsHI,'h986.cosmo50cmb.3072g14HBWK.00092',2,50000,1.84793e16,/cutout_rad,outarray= outarray
pro tipsysatshi, infile, groups, dist_units, mass_unit, cutout_rad=cutout_rad, outarray = outarray, $
                 outfile=outfile,num = num,gas = gas,dark = dark,stars = stars, gethi=gethi, ascii=ascii,  writemacro = writemacro, $
                 notipsy = notipsy, dwarfv = dwarfv, spinaxes = spinaxes,cx = cx, cy = cy, cz = cz,ind = ind, valign = valign,debug = debug

if not keyword_set(outfile) then outfile = infile

grpfile = infile+'.amiga.grp'
gtpfile = infile+'.amiga.gtp'
grp = read_ascii_array(grpfile)

;cmp_file = filename+'.cmp'

rtipsy, infile, h,g,d,s

if keyword_set(gethi) then hi = (read_ascii_array(infile+'.HI')) else hi = fltarr(N_ELEMENTS(h.ngas)) + 1.0

;*******************************************
;INDEX GAS,DM & STARS FROM PARTICULAR HALO
;*******************************************

for m = 0, n_elements(groups)-1 do begin

    i = groups[m]

    nostar = 0
    ind = where(grp eq i,comp=indcomp)
    IF ((where(ind ge h.ngas+h.ndark))[0] ne -1) THEN inds = ind(where(ind ge h.ngas+h.ndark)) ELSE nostar = 1
    indg = ind(where(ind lt h.ngas))
    indd = ind(where(ind ge h.ngas and ind lt h.ngas+h.ndark))
;tempdata = read_ascii_array(cmp_file)
;components = tempdata[h.ngas+h.ndark:h.ndark+h.nstar+h.ngas-1]

    IF NOT nostar THEN stars = s[inds-h.ngas-h.ndark]
    dark = d[indd-h.ngas] 
    gas = g[indg]


;***************************************************************************
; get CENTRE OF MASS from amiga (first entry in gtp file): REPOSITION
;**************************************************************************
;rtipsy, gtpfile, h1,g1,d1,s1

;cx= s1[i-1].x
;cy= s1[i-1].y
;cz= s1[i-1].z

    IF NOT nostar THEN BEGIN
        foo = min(stars.phi, cmindex)
        cx = stars[cmindex].x
        cy = stars[cmindex].y
        cz = stars[cmindex].z
        
        stars.x =  stars.x - cx
        stars.y =  stars.y - cy
        stars.z =  stars.z - cz
    ENDIF ELSE BEGIN
        foo = min(dark.phi, cmindex)
        cx = dark[cmindex].x
        cy = dark[cmindex].y
        cz = dark[cmindex].z    
    ENDELSE

    dark.x = dark.x - cx
    dark.y = dark.y - cy
    dark.z = dark.z - cz
    
    gas.x = gas.x - cx
    gas.y = gas.y - cy
    gas.z = gas.z - cz
    IF keyword_set(cutout_rad) THEN BEGIN
        IF NOT nostar THEN rads = sqrt(stars.x*stars.x+stars.y*stars.y+stars.z*stars.z)
        radg = sqrt(gas.x*gas.x+gas.y*gas.y+gas.z*gas.z)
        radd = sqrt(dark.x*dark.x+dark.y*dark.y+dark.z*dark.z)
        IF cutout_rad EQ 1 THEN $
          IF NOT nostar THEN maxrad = MAX([rads,radg,radd]) ELSE maxrad = MAX([radg,radd]) $
          ELSE maxrad = cutout_rad/dist_units
        IF h.nstar NE 0 THEN BEGIN
            s.x =  s.x - cx
            s.y =  s.y - cy
            s.z =  s.z - cz
        ENDIF
        d.x = d.x - cx
        d.y = d.y - cy
        d.z = d.z - cz
        g.x = g.x - cx
        g.y = g.y - cy
        g.z = g.z - cz
        IF h.nstar NE 0 THEN rads = sqrt(s.x*s.x+s.y*s.y+s.z*s.z)
        radg = sqrt(g.x*g.x+g.y*g.y+g.z*g.z)
        radd = sqrt(d.x*d.x+d.y*d.y+d.z*d.z)
        IF h.nstar NE 0 THEN nostar = 0
        IF h.nstar NE 0 THEN $
          IF (where(rads LT maxrad))[0] ne -1 THEN inds = where(rads LT maxrad) ELSE nostar = 1 ;note that these are in comoving kpc!
        indg = where(radg LT maxrad) ;note that these are in comoving kpc!
        indd = where(radd LT maxrad) ;note that these are in comoving kpc!
        gas = g[indg]
        IF NOT nostar THEN stars = s[inds]
        dark = d[indd]

        indd = indd + h.ngas
        IF NOT nostar THEN inds = inds + h.ngas + h.ndark 
        IF NOT nostar THEN ind = [indg,indd,inds] ELSE ind = [indg,indd]
    ENDIF

;***********************************************
;eliminate proper motion of the halo/galaxy
;***********************************************
    IF NOT nostar THEN BEGIN
        propx=mean(stars.vx)
        propy=mean(stars.vy)
        propz=mean(stars.vz)
    
        stars.vx=stars.vx-propx
        stars.vy=stars.vy-propy
        stars.vz=stars.vz-propz
    ENDIF ELSE BEGIN
        propx=mean(dark.vx)
        propy=mean(dark.vy)
        propz=mean(dark.vz)
    ENDELSE
    gas.vx=gas.vx-propx
    gas.vy=gas.vy-propy
    gas.vz=gas.vz-propz
    dark.vx=dark.vx-propx
    dark.vy=dark.vy-propy
    dark.vz=dark.vz-propz
    valign = [propx,propy,propz]
;    hi = hi[ind]

;**************************************************************
; find rotation curve and transtarslate in XY plane using align.pro
;*************************************************************
    limit=5/h.time/dist_units ;use stars within this limit (kpc/dist_units)
    print,'h.time: ',h.time
    IF keyword_set (gethi) THEN alignHI,stars,dark,gas,limit,hi = hi[ind],spinaxes = spinaxes ELSE align,gas,dark,stars,limit,spinaxes = spinaxes 
;sout = where(components ge 1 and components le 4)
;radius = sqrt(stars.x^2+stars.y^2+stars.z^2)
;limit = max(radius)
;gout = where (sqrt(gas.x*gas.x + gas.y*gas.y +gas.z*gas.z) lt limit)
;dout = where (sqrt(dark.x*dark.x + dark.y*dark.y +dark.z*dark.z) lt limit)
    
    num = h 
    IF NOT nostar THEN num.nstar = n_elements(stars) ELSE num.nstar = 0
    num.ngas = n_elements(gas)
    num.ndark = n_elements(dark)
    IF NOT nostar THEN num.n = num.nstar + num.ngas + num.ndark ELSE num.n = num.ngas + num.ndark 
    IF keyword_set(cutout_rad) THEN $
      IF cutout_rad NE 1.0 THEN fileout = outfile + '.halo.' + strtrim(string(i),2) + '.' + strtrim(cutout_rad,2) ELSE fileout = outfile + '.halo.' + strtrim(string(i),2)

    IF NOT keyword_set(notipsy) THEN BEGIN
    IF (1) THEN BEGIN
        if keyword_set(ascii) then $
          wtipsy, fileout + '.std', num, gas, dark, stars $
        else wtipsy,fileout + '.std',num,gas,dark,stars,/standard
    ENDIF ELSE wtipsy,fileout + '.std',num,gas,dark,stars
    ENDIF


    if keyword_set(outarray) then begin
        FOR k = 0, N_ELEMENTS(outarray) - 1 DO BEGIN
;            array = (read_ascii_array(infile+'.' +outarray[k])).c
            IF (outarray[k] eq 'decomp' OR outarray[k] eq 'iord' ) THEN read_tipsy_arr,infile+'.' +outarray[k],h,array,type = 'long' ELSE read_tipsy_arr,infile+'.' +outarray[k],h,array,type = 'float'
            openw,lunhi,fileout + '.' + outarray[k],/get_lun
            printf,lunhi,num.n
            for j=0L,n_elements(ind)-1 do begin
;                printf,lunhi,strtrim(array[ind[j]],2)
                IF (outarray[k] eq 'decomp' OR outarray[k] eq 'iord' ) THEN printf,lunhi,array[ind[j]],format = '(I)'  ELSE printf,lunhi,array[ind[j]],format = '(g)'
            endfor
            close,lunhi
        ENDFOR
    endif
    IF file_test(infile+'.H2') THEN molecularH = 1 ELSE molecularH = 0
    IF keyword_set(writemacro) THEN writeHImacro,outfile + '.halo.' + strtrim(string(i),2),dist_unit = dist_units,mass_unit = mass_unit,molecularH = molecularH,dwarfv = dwarfv

    IF keyword_set(debug) THEN stop

endfor
;close,/all
end

