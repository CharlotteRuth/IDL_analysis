pro stdsats, filename, h=h,g=g,d=d,s=s

grpfile = filename+'.amiga.grp'
gtpfile = filename+'.amiga.gtp'
statfile = filename+'.amiga.stat'
grp = read_lon_array(grpfile)

;***SARAH'S CHANGES***
FMT = 'X,X,F,F,F,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,A'
readcol, statfile, F=FMT, Number_Gas, Number_Star, Number_Dark, Contam
extension_name = ' '
Read, 'Please input the file prefix ', extension_name

galaxy_index = intarr(1) ;dummy value of 0

;making dummy values for i=0 (want i=1 to correspond largest halo in stat file)
Number_Gas=[0.,Number_Gas]
Number_Star=[0.,Number_Star]
Number_Dark=[0.,Number_Dark]
Contam=['contam',Contam]

;*****END CHANGES****

;cmp_file = filename+'.cmp'
rtipsy, filename, h,g,d,s
;*******************************************
;INDEX GAS,DM & STARS FROM PARTICULAR HALO
;*******************************************

;***MORE SARAH'S CHANGES***
;before the users specified the number of groups, but i just wanted it
;to run through all of them

ngroups = n_elements(contam)

for i = 1,ngroups-1 do begin     ;before went to ngroups instead of ngroups-1

if (Contam[i] eq 'clean') and (Number_Dark[i] ge 200) and (Number_Gas[i] gt 32) and $
(Number_Star[i] gt 0) then begin
    
;want to select clean halos w/ more than 200 dark matter particles.  
;***END SARAH'S CHANGES***

ind = (where(grp eq i,comp=indcomp))
inds = ind[where(ind ge h.ngas+h.ndark)]-h.ngas-h.ndark
indg = ind[where(ind lt h.ngas)]
indd = ind[where(ind ge h.ngas and ind lt h.ngas+h.ndark)]-h.ngas
;tempdata = read_ascii_array(cmp_file)
;components = tempdata[h.ngas+h.ndark:h.ndark+h.nstar+h.ngas-1]


stars = s[inds]
gas = g[indg]
dark = d[indd] 

;***************************************************************************
; get CENTRE OF MASS from amiga (first entry in gtp file): REPOSITION
;**************************************************************************
rtipsy, gtpfile, h1,g1,d1,s1

cx= s1[i-1].x
cy= s1[i-1].y
cz= s1[i-1].z

gas.x = gas.x - cx
gas.y = gas.y - cy
gas.z = gas.z - cz
dark.x = dark.x - cx
dark.y = dark.y - cy
dark.z = dark.z - cz

stars.x =  stars.x - cx
stars.y =  stars.y - cy
stars.z =  stars.z - cz


;***********************************************
;eliminate proper motion of the halo/galaxy
;***********************************************
propx=mean(stars.vx)
propy=mean(stars.vy)
propz=mean(stars.vz)

stars.vx=stars.vx-propx
stars.vy=stars.vy-propy
stars.vz=stars.vz-propz
gas.vx=gas.vx-propx
gas.vy=gas.vy-propy
gas.vz=gas.vz-propz
dark.vx=dark.vx-propx
dark.vy=dark.vy-propy
dark.vz=dark.vz-propz

;**************************************************************
; find rotation curve and transtarslate in XY plane using align.pro
;*************************************************************
;dist_units=1e5*h.time
;dist_units=2.85714e4*h.time
;limit=5./dist_units            ;use stars within this limit (kpc/dist_units)
;align,stars,dark,gas,limit



num = h 
num.nstar = n_elements(stars)
num.ngas = n_elements(gas)
num.ndark = n_elements(dark)
num.n = num.nstar + num.ngas + num.ndark

;fileout='X5.084'+strmid(filename,5,6,/reverse_offset)+'_'+strtrim(string(i),2)+'.asc'

;fileout='Cosmo50.036'+'_'+strtrim(string(i),2)+'.asc'

fileout = extension_name + '_' + strtrim(string(i),2) + '.asc'
galaxy_index = [galaxy_index, i]

wtipsy,fileout,num,gas,dark,stars,/standard

endif 
endfor

galaxy_index = galaxy_index[1:*] ;getting ride of dummy value 

fileout2 = 'gal_list'
openw, 1, fileout2

for k=0,n_elements(galaxy_index)-1 do begin
printf,1, strtrim(string(galaxy_index[k]),2)
endfor

close,1

end

