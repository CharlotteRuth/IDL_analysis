;This file reads in a director file and output a macro for looking at
;a tipsy file

;cd,'/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g16HBWK'
;tfile ='/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g14HBWK/steps/h516.cosmo25cmb.1536g14HBWK.00512.dir/h516.cosmo25cmb.1536g14HBWK'
;origtfile = 'h516.cosmo25cmb.1536g.tbin'
;tipsymovie,tfile,origtfile = origtfile


pro tipsymovie,tfile,director = director,origtfile = origtfile,pfile = pfile
IF NOT KEYWORD_SET(director) THEN spawn,'ls *director*',director
        rtipsy,tfile,h,g,d,s
FOR i = 0, N_ELEMENTS(director) - 1 DO BEGIN
    print,director[i]
;READ IN DATA
    temphead = strsplit(tfile,'.',/extract)
    head = temphead[N_ELEMENTS(temphead) - 1]
    name = director[i]
    readcol,name,data,format='(A)',DELIMITER = ':'
    fovind = where(strcmp(data,'FOV',3,/fold_case))
    fov = DOUBLE((strsplit(data[fovind],' ',/extract))[1])
    eyeind = where(strcmp(data,'eye',3,/fold_case))
    eye =  DOUBLE((strsplit(data[eyeind],' ',/extract))[1:3])
    zeyeind = where(strcmp(data,'zeye',4,/fold_case))
    zeye = DOUBLE((strsplit(data[zeyeind],' ',/extract))[1])
    targetind = where(strcmp(data,'target',6,/fold_case))
    photogenic = 0
    IF strmatch(data[targetind],'*photogenic*',/fold_case) THEN photogenic = 1 $
      ELSE target = DOUBLE((strsplit(data[targetind],' ',/extract))[1:3])

;CONVERT DATA
    IF photogenic THEN BEGIN
;        rtipsy,tfile,h,g,d,s
        IF NOT KEYWORD_SET(pfile) THEN spawn,'ls *photogenic*',pfile
        rtipsy,origtfile,horig,gorig,horig,sorig,/justhead
        close,/all
        readcol,pfile[0],iord,format='(d)'
        iord = LONG(iord) 
        dark_iord = iord[where(iord LT horig.ndark + horig.ngas AND iord GE horig.ngas)] 
        dhalo = d[dark_iord - horig.ngas] 
        cx = TOTAL(dhalo.x*dhalo.mass)/TOTAL([dhalo.mass])
        cy = TOTAL(dhalo.y*dhalo.mass)/TOTAL([dhalo.mass])
        cz = TOTAL(dhalo.z*dhalo.mass)/TOTAL([dhalo.mass])
        target = [cx,cy,cz]
    ENDIF ELSE  rtipsy,h,g,d,s,/justhead
    r = sqrt(eye[0]*eye[0] + eye[1]*eye[1] + eye[2]*eye[2])
    distance = r * zeye
    radius_ppm = distance*TAN(FOV*!PI/180.0/2.0)
    iangle = [0,0,1.0]
    print,eye
 
    IF (where(abs(eye) lt 1e-5))[0] ne -1 THEN  eye[where(abs(eye) lt 1e-5)] = 0
    psi = acos(eye[2]/r)*180/!PI
;    if eye[2] lt 0 then psi = 180.0 + psi 
    IF (abs(eye[0]) gt 1e-5 AND abs(eye[1]) gt 1e-5) then theta = atan(eye[1]/eye[0])*180.0/!PI else begin
        if abs(eye[1]) gt 1e-5 then theta = asin(eye[1]/r)*180.0/!PI else theta = 0
    endelse
    IF (eye[0] lt 0 ) then theta= theta + 180
    if theta lt 0 then theta = theta + 360
    print,'Psi: ',psi
    print,'Theta: ',theta

;WRITE MACRO
    openw,1,'phot' +head + '.' + STRTRIM(i,2) + '.macro'
    printf,1,'displayppm'
    printf,1,' '
    printf,1,"openb "+tfile
    printf,1,"loads 1"
    printf,1,"xall"
    printf,1,"redshift ",STRTRIM(1.0/h.time - 1.0)," c 73 0.24 0.76 p 1"    
    printf,1,'setsphere 1 ',STRTRIM(target[0],2),' ',STRTRIM(target[1],2),' ',STRTRIM(target[2],2),' ',STRTRIM(radius_ppm,2)
    printf,1,'abox 1'
    printf,1,'zgas logrho 0 10'
    printf,1,'rotate c ',theta
    printf,1,'viewgas logrho 0 10'
 printf,1,'rotate left ',psi 
;    IF (eye[0] ge 1) THEN printf,1,'rotate left ',psi ELSE $
;      printf,1,'rotate right ',psi
;    IF (eye[2] gt 0) THEN BEGIN
;        if (theta gt 0) then printf,1,'rotate up ',theta else  printf,1,'rotate down ',-1.0*theta 
;    ENDIF ELSE BEGIN
;        if (theta gt 0) then printf,1,'rotate down ',theta else  printf,1,'rotate up ',-1.0*theta
;    ENDELSE
    printf,1,'viewgas logrho 0 10'
    printf,1,"closeb"
    printf,1,"end"
    close,1
ENDFOR

end
