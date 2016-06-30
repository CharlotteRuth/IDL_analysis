
;dir ='/astro/net/scratch2/christensen/MolecH/Cosmo/h603.cosmo50cmb.1536g/h603.cosmo50cmb.1536g3HBWK/'
;cd,dir
;makemovie,'steps/h603.cosmo50cmb.1536g3HBWK.00300.dir/h603.cosmo50cmb.1536g3HBWK.00300','h603.cosmo50cmb.1536g3HBWK.mark_1_gd',1e-3,pfile='h603.cosmo50cmb.1536g3HBWK',dfile = 'h603.cosmo50cmb.1536g3HBWK',/rotate


;This function returns a tipsy particle whoes position has been
;shifted over by the center vector
FUNCTION center_particle,particle,center

particle.x = particle.x - center[0]
particle.y = particle.y - center[1]
particle.z = particle.z - center[2]

RETURN,particle
END


;This procedure will make a mark file and output only the marks for
;the gas and dark particle
pro select_gd_mark,mfile,halo = halo
if not keyword_set(halo) then halo = '1'
readcol,mfile+'.' + halo + '.mark',n,ngas,nstar,format = 'L,L,L',/silent
n = n[0]
ngas = ngas[0]
nstar = nstar[0]
ndark = n - ngas - nstar

readcol,mfile+'.' + halo + '.mark',iord,format = 'L'
gas_iord = iord[where(iord LT ngas)]
dark_iord = iord[where(iord LT ndark + ngas AND iord GE ngas)]
IF nstar NE 0 THEN star_iord = iord[where(iord LT n AND iord GE ndark + ngas)]
openw,1,mfile+'.' + halo + '.mark_gd'
;printf,1,ngas+ndark,ngas,0,format='(3I)'
printf,1,n,ngas,nstar,format='(3I12)'
printf,1,TRANSPOSE([gas_iord,dark_iord]),format='(I12)'
close,1
end

pro select_d_mark,mfile,halo = halo
if not keyword_set(halo) then halo = '1'
readcol,mfile+'.' + halo + '.mark',n,ngas,nstar,format = 'L,L,L',/silent
n = n[0]
ngas = ngas[0]
nstar = nstar[0]
ndark = n - ngas - nstar

readcol,mfile+'.' + halo + '.mark',iord,format = 'L'
gas_iord = iord[where(iord LT ngas)]
dark_iord = iord[where(iord LT ndark + ngas AND iord GE ngas)]
IF nstar NE 0 THEN star_iord = iord[where(iord LT n AND iord GE ndark + ngas)]
openw,1,mfile+'.' + halo + '.mark_d'
;printf,1,ngas+ndark,ngas,0,format='(3I)'
printf,1,n,ngas,nstar,format='(3I12)'
printf,1,TRANSPOSE([dark_iord]),format='(I12)'
close,1
end



;Given a mark file, this will output a photogenic file
pro make_photogenic,mfile,filename = filename
IF NOT KEYWORD_SET(filename) THEN filename = mfile + '.photogenic'

openr,1,mfile
readf,1,n,ngas,nstar;,format='(3I)'
ndark = n - ngas - nstar
close,1
readcol,mfile,iord,format='(d)'
iord = LONG(iord)
ngashalo = 1
IF (where(iord LT ngas))[0] ne -1 THEN gas_iord = iord[where(iord LT ngas)] - 1 ELSE ngashalo = 0;Mark file iorders are +1
dark_iord = iord[where(iord LT ndark + ngas AND iord GE ngas)] - 1
IF ngashalo NE 0 THEN writecol,filename,[gas_iord,dark_iord],format='(I)' ELSE writecol,filename,[dark_iord],format='(I)' 
END

;Given a mark file, this will output a photogenic file
pro make_mark,pfile,tfile,filename = filename
IF NOT KEYWORD_SET(filename) THEN filename = mfile + '.mrk'

rtipsy,tfile,h,g,d,s,/justhead
readcol,pfile,iord,format='(d)'
iord = LONG(iord)
miord = iord + 1

openr,filename,1
printf,1,h.n,h.ngas,h.nstar
printf,1,miord
close,1


;readf,1,n,ngas,nstar;,format='(3I)'
;ndark = n - ngas - nstar
;close,1
;readcol,mfile,iord,format='(d)'
;iord = LONG(iord)

;ngashalo = 1
;IF (where(iord LT ngas))[0] ne -1 THEN gas_iord = iord[where(iord LT ngas)] - 1 ELSE ngashalo = 0;Mark file iorders are +1
;dark_iord = iord[where(iord LT ndark + ngas AND iord GE ngas)] - 1
;IF ngashalo NE 0 THEN writecol,filename,[gas_iord,dark_iord],format='(I)' ELSE writecol,filename,[dark_iord],format='(I)' 

END

;this function will write a director file
pro writedirector, filename, number, suffix, mass_scale, fov, eye, zeye, up,colorscheme,comment = comment
    openw,1,filename+'.director' + number
    printf,1,'# size: ',strtrim(2.0*zeye*TAN(FOV*!PI/360.0),2)
    printf,1,'file movie/',filename+'.' + suffix
    printf,1,'encode ppm'
    printf,1,'size 800 800' ;Could use to make larger movie but this is probably fine
    printf,1,''
    printf,1,'clip 0.2 2.'   ;Can change to keep more gas in (0.3 2.0)
    printf,1,'render tsc'
    printf,1,'target photogenic'
    printf,1,'physical'
    printf,1,''
    printf,1,'project perspective'
    printf,1,'#softdark 0.2' ;0.5
    printf,1,'softgas 0.25'
    printf,1,'softstar 0.2'     ;0.125
    printf,1,'softgassph'
    printf,1,    '#coldark            0.00  0.00  1.00  1e-9'
    CASE colorscheme OF
        1: BEGIN ; My golden color scheme
            printf,1,'colgas              1.00  0.78  0.54 ',STRTRIM(mass_scale/250.0,2)
            printf,1,'colstarbricol       0.78  0.86  1.00 ',STRTRIM(mass_scale/1500.0,2) ;The smaller the number, the brighter the light
            printf,1,'logscalecoloursafe  10 50000' ;
        END
        2: BEGIN;Fabio scheme highlighting stars
            printf,1,'colgas              0.70  0.80  0.85 ',STRTRIM(mass_scale,2)
            printf,1,'colstarbricol       1.00  1.00  1.00 ',STRTRIM(mass_scale/1500.0,2) ;The smaller the number, the brighter the light
            printf,1,'logscalecoloursafe  1 5000' ;
        END
    ENDCASE
    printf,1,'dark off'
    printf,1,'FOV ',fov
    printf,1,'up ',up
    printf,1,'eye ',eye
    printf,1,'zeye ',STRTRIM(zeye,2)
    close,1
end

;This procedure has the ability, given a mark file, to write a
;photogenic file and two director file

;tfile is the name of a tipsy file
;radius_ppm should be the radius your sphere you would like included in your movie, in simulation units
;mfile is the name of a file in which the particles used for calculating the angle are marked

;makedirector,'h516.cosmo25cmb.2304g14HBWK.00512',6e-4,fov = 60,filename = 'h516.cosmo25cmb.2304g14HBWK',pfile = 'h516.cosmo25cmb.2304g14HBWK.photogenic',/angleup,origtfile = 'h516.cosmo25cmb.2304g.tbin',theta = [0,45,90],suffix = ['1FO','1AO','1EO']


;makedirector,'h516.cosmo25cmb.1536g14HBWK',6e-4,fov = 60,theta=[0,45,90],suffix = ['1FO','1AO','1EO'],dsuffix = ['t','t1','t2']
;makedirector,'h516.cosmo25cmb.2304g14HBWK',6e-4,fov = 60,theta =[45],suffix = ['AO'],dsuffix = ['t']

;makedirector,'h603.cosmo50cmb.2304g14HBWK',0.00128000,theta=[45],suffix = ['AO'],dsuffix = ['t']

;makedirector,'h516.cosmo25cmb.1536g14HBWK',2.3e-4,theta = [0,45,90],suffix =['FO','AO','EO'],dsuffix=['t','t1','t2'],angleup = [0.33030044 ,-0.59726007, -0.73087758]

;makedirector,'h516.cosmo25cmb.1536g16HBWK',4.8e-4,theta = [0,45,90],suffix =['FO','AO','EO'],angleup = [0.33030044 ,-0.59726007, -0.73087758]
;makedirector,'h516.cosmo25cmb.1536g16HBWK',4.8e-4,theta =[45],suffix=['AOstar'],angleup = [0.33030044 ,-0.59726007,-0.73087758],dsuffix= ['4'],colorscheme = 2

;cd,'/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK'
;makedirector,'h516.cosmo25cmb.3072g14HBWK',4.8e-4,theta=[0,45,90],suffix =['FO','AO','EO'],angleup = [-0.225629  ,-0.678328,   -0.699259],dsuffix=['3','4','5']
;makedirector,'h516.cosmo25cmb.3072g14HBWK',4.8e-4,theta=[45],suffix =['AOstar'],angleup = [-0.225629  ,-0.678328,   -0.699259],dsuffix=['6'],colorscheme = 2

;cd,'/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK'
;makedirector,'h603.cosmo50cmb.3072g14HBWK',0.00128000,theta=[0,45,90],suffix=['FO','AO','EO'],angleup= [0.241514,-0.730361,-0.638940],dsuffix = ['4','5','6']
;makedirector,'h603.cosmo50cmb.3072g14HBWK',0.00128000,theta=[45],suffix=['AOstar'],angleup= [0.241514,-0.730361,-0.638940],dsuffix = ['7']

;cd,'/astro/store/student-scratch1/christensen/MolecH/Cosmo/h937.cosmo25cmb.4096g/h937.cosmo25cmb.4096g14HBWK'
;makedirector,'h937.cosmo25cmb.4096g14HBWK',0.00256000,theta=[0,45,90],suffix=['FO','AO','EO']
;makedirector,'h937.cosmo25cmb.4096g14HBWK',0.00256000,theta=[45],suffix=['AOstar'],dsuffix= ['4'],colorscheme = 2

pro makedirector, filename, radius_ppm, paramfile = paramfile, fov = fov, tfile = tfile, mfile = mfile, pfile = pfile, origtfile = origtfile, angleup = angleup, theta = theta, suffix = suffix, colorscheme = colorscheme,dsuffixd = dsuffix
IF NOT KEYWORD_SET(FOV) THEN FOV = 40
distance = radius_ppm/TAN(FOV*!PI/180.0/2.0)
IF NOT KEYWORD_SET(tfile) THEN tfile = filename
IF NOT KEYWORD_SET(paramfile) THEN paramfile = filename + '.param'
IF NOT (KEYWORD_SET(colorscheme)) THEN colorscheme = 1

units = tipsyunits(paramfile)
print,'Gas Particle Mass: ',units.gasparmass*units.massunit
print,'Star Particle Mass: ',units.ISTARMASS*units.massunit
mass_scale = 20000.0/units.massunit ;units.istarmass

IF KEYWORD_SET(angleup) THEN BEGIN
    IF (N_ELEMENTS(angleup) EQ 3) THEN BEGIN
        spinvectorz = angleup/MAG(angleup)
    ENDIF ELSE BEGIN
        rtipsy,tfile,h,g,d,s

;give a markfile to select a halo for use in angleup 
        IF KEYWORD_SET(mfile) THEN BEGIN
            openr,1,mfile
            readf,1,n,ngas,nstar ;,format='(3I)'
            ndark = n - ngas - nstar
            close,1
            readcol,mfile,iord,format='(d)'
            iord = LONG(iord) - 1 ;Mark file iorders are +1
            ngashalo = 1
            IF (where(iord LT ngas))[0] ne -1 THEN gas_iord = iord[where(iord LT ngas)] ELSE ngashalo = 0 
            dark_iord = iord[where(iord LT ndark + ngas AND iord GE ngas)] 
            IF ngashalo NE 0 THEN ghalo = g[gas_iord]
            dhalo = d[dark_iord - ngas]
        ENDIF
        
        IF KEYWORD_SET(pfile) THEN BEGIN
;give a photogenic to select a halo for use in angleup 
            rtipsy,origtfile,horig,gorig,horig,sorig,/justhead
            close,/all
            nstar = h.nstar
            n = h.n
            readcol,pfile,iord,format='(d)'
            iord = LONG(iord) 
            ngashalo = 1
            IF (where(iord LT horig.ngas))[0] ne -1 THEN gas_iord = iord[where(iord LT horig.ngas)] ELSE ngashalo = 0 
            dark_iord = iord[where(iord LT horig.ndark + horig.ngas AND iord GE horig.ngas)] 
            IF ngashalo NE 0 THEN ghalo = g[gas_iord]
            dhalo = d[dark_iord - horig.ngas]  
        ENDIF
        
        IF KEYWORD_SET(pfile) OR KEYWORD_SET(mfile) THEN BEGIN
            IF ngashalo NE 0 THEN BEGIN
                cx = TOTAL([ghalo.x*ghalo.mass,dhalo.x*dhalo.mass])/TOTAL([ghalo.mass, dhalo.mass])
                cy = TOTAL([ghalo.y*ghalo.mass,dhalo.y*dhalo.mass])/TOTAL([ghalo.mass, dhalo.mass])
                cz = TOTAL([ghalo.z*ghalo.mass,dhalo.z*dhalo.mass])/TOTAL([ghalo.mass, dhalo.mass])
                center = [cx,cy,cz]
                radius = MAX(SQRT( ([ghalo.x,dhalo.x]-cx)*([ghalo.x,dhalo.x]-cx) + ([ghalo.y,dhalo.y]-cy)*([ghalo.y,dhalo.y]-cy) + ([ghalo.z,dhalo.z]-cz)*([ghalo.z,dhalo.z]-cz)))
            ENDIF ELSE BEGIN
                cx = TOTAL(dhalo.x*dhalo.mass)/TOTAL([dhalo.mass])
                cy = TOTAL(dhalo.y*dhalo.mass)/TOTAL([dhalo.mass])
                cz = TOTAL(dhalo.z*dhalo.mass)/TOTAL([dhalo.mass])
                center = [cx,cy,cz]
                radius = MAX(SQRT( (dhalo.x-cx)*(dhalo.x-cx) + (dhalo.y-cy)*(dhalo.y-cy) + (dhalo.z-cz)*(dhalo.z-cz)))   
            ENDELSE
        ENDIF ELSE BEGIN      ;IF you want to angleup the whole galaxy
            dhalo = d
            ghalo = g
            cx = 0
            cy = 0
            cz = 0
            center = [cx,cy,cz]    
            radius = radius_ppm/10.0
        ENDELSE
        
        
        dhalo_all = center_particle(d[WHERE(SQRT((d.x-cx)*(d.x-cx) + (d.y-cy)*(d.y-cy) + (d.z-cz)*(d.z-cz)) le radius_ppm)],center)
        ghalo_all = center_particle(g[WHERE(SQRT((g.x-cx)*(g.x-cx) + (g.y-cy)*(g.y-cy) + (g.z-cz)*(g.z-cz)) le radius_ppm)],center)
        IF nstar ne 0 THEN shalo_all = center_particle(s[WHERE(SQRT((s.x-cx)*(s.x-cx) + (s.y-cy)*(s.y-cy) + (s.z-cz)*(s.z-cz)) le radius_ppm)],center) ELSE shalo_all = [0]
        alignHI,shalo_all,dhalo_all,ghalo_all,radius_ppm/10,spinaxes = spinaxes
        
        set_plot,'x'
        loadct,39
        !p.multi = [0,2,1]
        plot,dhalo_all.x,dhalo_all.y,psym = 3,xtitle = 'X',ytitle = 'Y',xrange = [-1.0*radius_ppm,radius_ppm],yrange = [-1.0*radius_ppm,radius_ppm]
        oplot,ghalo_all.x,ghalo_all.y,psym = 3,color = 240
        IF nstar ne 0 THEN oplot,shalo_all.x,shalo_all.y,psym = 3,color = 100
        oplot,dhalo.x - cx,dhalo.y - cy,psym = 2,color = 140
        IF ngashalo NE 0 THEN  oplot,ghalo.x - cx,ghalo.y - cy,psym = 2,color = 140
        
        plot,dhalo_all.x,dhalo_all.z,psym = 3,xtitle = 'X',ytitle = 'Z',xrange = [-1.0*radius_ppm,radius_ppm],yrange = [-1.0*radius_ppm,radius_ppm]
        oplot,ghalo_all.x,ghalo_all.z,psym = 3,color = 240
        IF nstar ne 0 THEN oplot,shalo_all.x,shalo_all.z,psym = 3,color = 100
        oplot,dhalo.x - cx,dhalo.z - cz,psym = 2,color = 140
        IF ngashalo NE 0 THEN  oplot,ghalo.x - cx,ghalo.z - cz,psym = 2,color = 140
        !p.multi = 0
        spinvectorz = spinaxes
    ENDELSE
ENDIF ELSE spinvectorz = [0,0,1.0] ;spinaxes = [[0,1e-9,1],[1,1e-9,1e-9],[0.707107,1e-9,0.707107]]
x = spinvectorz[0]
y = spinvectorz[1]
z = spinvectorz[2]
spinvectorx = [z/sqrt(x*x + z*z),0,-1.0*x/sqrt(x*x + z*z)]
spinvectory = crossp(spinvectorz,spinvectorx);[0,-1.0*z/sqrt(y*y + z*z),y/sqrt(y*y + z*z)]
basis = [[spinvectorx],[spinvectory],[spinvectorz]]

IF NOT KEYWORD_SET(theta) THEN theta = [0.0,45.0,90.0]
theta = theta/180.0*!PI
IF NOT KEYWORD_SET(suffix) OR (N_ELEMENTS(suffix) NE N_ELEMENTS(theta)) THEN suffix = sindgen(N_ELEMENTS(theta))

FOR i = 0, N_ELEMENTS(theta) - 1 DO BEGIN
    eye = [SIN(theta[i]),0,COS(theta[i])]
    eye[where(eye eq 0)] = 1e-9
    up =  [SIN(theta[i] - !PI/2),0,COS(theta[i] - !PI/2)]
    up[where(up eq 0)] = 1e-9
    IF KEYWORD_SET (dsuffix) THEN number = dsuffix[i] ELSE BEGIN
       IF i eq 0 THEN number = '' ELSE  number = STRTRIM(i+1,2)
    ENDELSE
    writedirector, filename, number, STRTRIM(suffix[i],2), mass_scale, fov, basis#eye, distance, basis#up, colorscheme
ENDFOR
END




PRO junk
IF KEYWORD_SET(dfile) THEN BEGIN
    openw,1,dfile+'.director2'
    printf,1,'file movie/',dfile+'.FO'
    printf,1,'encode ppm'
    printf,1,'size 800 800' ;Could use to make larger movie but this is probably fine
    printf,1,''
    printf,1,'clip 0.2 2.'   ;Can change to keep more gas in (0.3 2.0)
    printf,1,'render tsc'
    printf,1,'target photogenic'
    printf,1,'physical'
    printf,1,''
    printf,1,'project perspective'
;printf,1,'softdark 0.2' ;0.5
    printf,1,'softgas 0.25'
    printf,1,'softstar 0.2'     ;0.125
    printf,1,'softgassph'
    printf,1,'logscalecoloursafe  10 50000' ;
;printf,1,'coldark 0 0 1 1e-9'
    printf,1,'colgas 1.0 0.78  0.54 ',STRTRIM(smass/250.0,2)
    printf,1,'colstarbricol 0.78 0.86 1.0 ',STRTRIM(smass/1500.0,2) ;
    printf,1,'dark off'
    printf,1,'FOV ',FOV
    printf,1,'up 0 1 0'
    printf,1,'eye ',spinaxes[*,1]
    printf,1,'zeye ',STRTRIM(distance,2)
    close,1

    openw,1,dfile+'.director3'
    printf,1,'file movie/',dfile+'.AO'
    printf,1,'encode ppm'
    printf,1,'size 800 800' ;Could use to make larger movie but this is probably fine
    printf,1,''
    printf,1,'clip 0.2 2.'   ;Can change to keep more gas in (0.3 2.0)
    printf,1,'render tsc'
    printf,1,'target photogenic'
    printf,1,'physical'
    printf,1,''
    printf,1,'project perspective'
;printf,1,'softdark 0.2' ;0.5
    printf,1,'softgas 0.25'
    printf,1,'softstar 0.2'     ;0.125
    printf,1,'softgassph'
    printf,1,'logscalecoloursafe  10 50000' ;
;printf,1,'coldark 0 0 1 1e-9'
    printf,1,'colgas 1.0 0.78  0.54 ',STRTRIM(smass/250.0,2)
    printf,1,'colstarbricol 0.78 0.86 1.0 ',STRTRIM(smass/1500.0,2) ;
    printf,1,'dark off'
    printf,1,'FOV ',FOV
    printf,1,'up 0 1 0'
    printf,1,'eye ',spinaxes[*,2]
    printf,1,'zeye ',STRTRIM(distance,2)
    close,1
ENDIF
end

;logscalecoloursafe 30 45000
;colgas 0.71 0.81  0.86 5.e-14
;colstarbricol 1 1 1 0.1e-14
;coldark 0 0 1 3e-13
;dark off
;FOV 25
;up 0 0 1
;eye -.00025 0.00025 -0.000925
;time 0.2

;h603
;logscalecoloursafe 50 50000
;colgas 0.65 0.9  0.75 5.5e-15
;colstarbricol 1 1 1 1e-14
;coldark 0 0 1 1e-13
;dark off

;FOV 30
;up 0 0 1
;eye  -.01 0.007 -0.0092
