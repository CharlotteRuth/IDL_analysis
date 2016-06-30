;.r /astro/users/christensen/code/procedures/wtipsy.pro
;.r /astro/users/christensen/code/idl4tipsy/tipsySatsHI.pro
;.r /astro/users/christensen/code/idl4tipsy/alignHI.pro

pro tipsySatsArrays, infile, groups, dist_units, mass_unit, outfile=outfile, h=h,g=g,d=d,s=s,  $
               ascii=ascii,  cutout_rad=cutout_rad, outarray = outarray

if not keyword_set(outfile) then outfile = infile

grpfile = infile+'.amiga.grp'
gtpfile = infile+'.amiga.gtp'
grp = read_ascii_array(grpfile)
rtipsy, infile, h,g,d,s


;*******************************************
;INDEX GAS,DM & STARS FROM PARTICULAR HALO
;*******************************************

for m = 0, n_elements(groups)-1 do begin

    i = groups[m]

    ind = where(grp eq i,comp=indcomp)
    inds = ind(where(ind ge h.ngas+h.ndark))
    indg = ind(where(ind lt h.ngas))
    indd = ind(where(ind ge h.ngas and ind lt h.ngas+h.ndark))
;tempdata = read_ascii_array(cmp_file)
;components = tempdata[h.ngas+h.ndark:h.ndark+h.nstar+h.ngas-1]

    stars = s[inds-h.ngas-h.ndark]
    dark = d[indd-h.ngas] 
    gas = g[indg]

;***************************************************************************
; get CENTRE OF MASS from amiga (first entry in gtp file): REPOSITION
;**************************************************************************
;rtipsy, gtpfile, h1,g1,d1,s1

;cx= s1[i-1].x
;cy= s1[i-1].y
;cz= s1[i-1].z

    foo = min(stars.phi, cmindex)
    cx = stars[cmindex].x
    cy = stars[cmindex].y
    cz = stars[cmindex].z
    
    stars.x =  stars.x - cx
    stars.y =  stars.y - cy
    stars.z =  stars.z - cz
    
    dark.x = dark.x - cx
    dark.y = dark.y - cy
    dark.z = dark.z - cz
    
    gas.x = gas.x - cx
    gas.y = gas.y - cy
    gas.z = gas.z - cz

    IF keyword_set(cutout_rad) THEN BEGIN
        rads = sqrt(stars.x*stars.x+stars.y*stars.y+stars.z*stars.z)
        radg = sqrt(gas.x*gas.x+gas.y*gas.y+gas.z*gas.z)
        radd = sqrt(dark.x*dark.x+dark.y*dark.y+dark.z*dark.z)
        cutout_rad = MAX([rads,radg,radd])
        s.x =  s.x - cx
        s.y =  s.y - cy
        s.z =  s.z - cz
        d.x = d.x - cx
        d.y = d.y - cy
        d.z = d.z - cz
        g.x = g.x - cx
        g.y = g.y - cy
        g.z = g.z - cz
        rads = sqrt(s.x*s.x+s.y*s.y+s.z*s.z)
        radg = sqrt(g.x*g.x+g.y*g.y+g.z*g.z)
        radd = sqrt(d.x*d.x+d.y*d.y+d.z*d.z)
        inds = where(rads lt cutout_rad) ;note that these are in comoving kpc!
        indg = where(radg lt cutout_rad) ;note that these are in comoving kpc!
        indd = where(radd lt cutout_rad) ;note that these are in comoving kpc!
        gas = g[indg]
        stars = s[inds]
        dark = d[indd]

        indd = indd + h.ngas
        inds = inds + h.ngas + h.ndark 
        ind = [indg,indd,inds]
    ENDIF


    
    num = h 
    num.nstar = n_elements(stars)
    num.ngas = n_elements(gas)
    num.ndark = n_elements(dark)
    num.n = num.nstar + num.ngas + num.ndark
    fileout = outfile + '.halo.' + strtrim(string(i),2) + '.std'

    IF (1) THEN BEGIN
        if keyword_set(ascii) then $
          wtipsy, fileout, num, gas, dark, stars $
        else wtipsy,fileout,num,gas,dark,stars,/standard
    ENDIF ELSE wtipsy,fileout,num,gas,dark,stars


    if keyword_set(outarray) then begin
        FOR k = 0, N_ELEMENTS(outarray) - 1 DO BEGIN
;            array = (read_ascii_array(infile+'.' +outarray[k])).c
            IF (outarray[k] eq 'decomp' OR outarray[k] eq 'iord' ) THEN read_tipsy_arr,infile+'.' +outarray[k],h,array,type = 'long' ELSE read_tipsy_arr,infile+'.' +outarray[k],h,array,type = 'float'
            fileout = outfile + '.halo.' + strtrim(string(i),2) + '.' +outarray[k];+'.HI'
            openw,lunhi,fileout,/get_lun
            printf,lunhi,num.n
            for j=0L,n_elements(ind)-1 do begin
;                printf,lunhi,strtrim(array[ind[j]],2)
                IF (outarray[k] eq 'decomp' OR outarray[k] eq 'iord' ) THEN printf,lunhi,array[ind[j]],format = '(I)'  ELSE printf,lunhi,array[ind[j]],format = '(g)'
            endfor
            close,lunhi
        ENDFOR
    endif

endfor
end

