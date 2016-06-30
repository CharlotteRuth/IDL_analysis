FUNCTION starmass,filename,stars,inds,ms_units,index,s

;pro starmass,tipsyfile,halonumb,dMsolUnit
;if (N_PARAMS() eq 0) then begin
;    print, 'Format iput as follows'
;    print, 'starmass, tipsyfile, Amiga halo number, dMsolUnit(from the param file)'
;endif

;components = filename+'.cmp'
;group = filename+'.amiga.grp'

;Read in the important files
;rtipsy,tipsyfile,h,g,d,s
;readcol,components,comp,format='F',/silent
;readcol,group,galaxy,format='F',/silent

;galaxy = grp
comp = index
dMsolUnit = ms_units


;break up the and grp
;ind = where(galaxy eq halo,complement=indcomp)
;inds = ind(where(ind ge h.ngas+h.ndark)) 
;indg = ind(where(ind lt h.ngas))
;indd = ind(where(ind ge h.ngas and ind lt h.ngas+h.ndark))
;indsats= indcomp(where(indcomp ge h.ngas+h.ndark))


;inde=lonarr(h.nstar)
;for i=0L,h.nstar-1 do inde[i] = i
;indstars=inde[inds-h.ngas-h.ndark]

;stars = s[inds-h.ngas-h.ndark] 

starcomp = comp[inds-1]

disk = where(starcomp eq 1)
halo = where(starcomp eq 2)
bulge = where(starcomp eq 3)
thick = where(starcomp eq 4)
psbulge = where(starcomp eq 5)

diskstars = s[starcomp[disk]]
halostars = s[starcomp[halo]]
bulgestars = s[starcomp[bulge]]
thickstars = s[starcomp[thick]]
psbulgestars = s[starcomp[psbulge]]


diskstarmass=total(diskstars.mass)*dMsolUnit
halostarmass=total(halostars.mass)*dMsolUnit
bulgestarmass=total(bulgestars.mass)*dMsolUnit
thickstarmass=total(thickstars.mass)*dMsolUnit
psbulgestarmass=total(psbulgestars.mass)*dMsolUnit

BD = (bulgestarmass) / (diskstarmass + thickstarmass + psbulgestarmass)


openw,1,filename+'.masses'
printf,1,'Masses in Solar Units '
printf,1,'For Tipsy Units divide by dMsolUnit'
printf,1,' '
printf,1,'B/D = '+string(BD)
printf,1,'For this ratio the disk is considered to include the'
printf,1,'disk, thickdisk, and pseudobulge'
printf,1,' '
printf,1,'The Total Simulation'
printf,1,'Stellar mass = '+string(total(s.mass)*dMsolUnit)
printf,1,' '
printf,1,'The Whole Galaxy'
printf,1,'Stellar mass = '+string(total(stars.mass)*dMsolUnit)
printf,1,' '
printf,1,'Disk mass = '+string(diskstarmass)
printf,1,' '
printf,1,'Thickdisk mass = '+string(thickstarmass)
printf,1,' '
printf,1,'Bulge mass = '+string(bulgestarmass)
printf,1,' '
printf,1,'Pseudobulge mass = '+string(psbulgestarmass)
printf,1,' '
printf,1,'Halo mass = '+string(halostarmass)

close,1

RETURN,BD

END
