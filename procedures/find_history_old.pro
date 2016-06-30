;******************************************************************************
;This function analyzes the history of gas and star particles, by
;taking a list of input particles and outputing a structure containing
;their properties as a function of time.
;
;REQUIREMENTS:
;     A list of simulation outputs must exist in a file with name
;     simoutputlist
;
;     This procedure requires two arrays: orig_iord, of the indices
;     of the objects you'd like to track, as given in the *iord files, 
;     and orig_igasorder, which has the igasorder indices of each
;     particle (for gas particles this will be zero).  If you have
;     only gas particles then you can choose not to specify
;     orig_igasorder. 
;
;     Use std_array_convert to convert all the .iord files into ascii
;     files called, for example, 50.iord.asc.  
;
;     Use tipsy's writebox command to dump binaries into ascii files,
;     called, for example, 50.writebox.
;
;     The program can also deal with SKID outputs, .stat and .grp
;     files.  To tell the program to use skid outputs run the program
;     with the /SKID flag set.  If this flag is set, the output
;     structure will contain the ID of all halos (history.haloid), the
;     max circular velocity in the halos (.vcmax) and the distance
;     from the center of the gal (.halodistance)
;
;     With the /WRITEMARK flag set, the program will write out a tipsy
;     mark file for each timestep using the name given in the
;     simoutputlist file <name>.mark
;******************************************************************************
FUNCTION find_history_old, simoutputlist, orig_iord, orig_igasorder,MARKFILE=markfile,SKID=skid,WRITEMARK=writemark


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; the next block contains information to translate into physical
; units.  CHECK carefully that they are correct for your simulation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
dKpcUnit=1.d5 ; to convert in kpc
massunit=13.6d16 ; to convert into solar masses
velunit=2.41846d3 ; to convert to km/sec
hubble=70./3.086d19 ; hubble constant in 1/sec
dRhoUnit=3.*hubble^2/(8.*!DPI*6.67259d-8) ; to convert to g/cc
timeunit=3.88d10 ;may not use this, because we want lookback and z
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;Read in the list of output filenames to search.
READCOL, simoutputlist, root, format = '(a)'

;get number of timesteps and desired particles
ntimesteps=N_ELEMENTS(root)
IF KEYWORD_SET(markfile) THEN BEGIN
  rdfloat,markfile,marked,skip=1
  io = rbarray(root[ntimesteps-1]+".iord")
  orig_iord=io[marked-1] 
ENDIF
norig=N_ELEMENTS(orig_iord)

print,ntimesteps

;create structure to hold information of particles over time:
;history structure is subscripted by number of desired particles, many
;elements within it are subscripted by number of timesteps
history=REPLICATE({iord:0l, igasorder:0l, type:STRARR(ntimesteps), mass:DBLARR(ntimesteps), x:DBLARR(ntimesteps), y:DBLARR(ntimesteps), z:DBLARR(ntimesteps), vx:DBLARR(ntimesteps), vy:DBLARR(ntimesteps), vz:DBLARR(ntimesteps), rho:DBLARR(ntimesteps), temp:DBLARR(ntimesteps), metallicity:DBLARR(ntimesteps), redshift:DBLARR(ntimesteps), lookback:DBLARR(ntimesteps), expfactor:DBLARR(ntimesteps), timestep:LONARR(ntimesteps), halovcmax:DBLARR(ntimesteps), haloid:LONARR(ntimesteps), halodistance:DBLARR(ntimesteps)}, norig)

history.iord=orig_iord
IF ARG_PRESENT(orig_igasorder) THEN history.igasorder=orig_igasorder

;For each output,
FOR i = 0l, ntimesteps - 1 DO BEGIN

    ;Read in the .iord file.
    ;iord=read_ascii_array(root[i]+'.iord.asc',/LONG)

    print,"Working on "+root[i]
    iord=rbarray(strcompress(root[i]+'.iord'))

    ;Read in the .writebox file and get its dimensions
    ;catalog = read_writebox_output(strcompress(root[i]+'.writebox'))
    ;dim=get_dim(strcompress(root[i]+'.writebox'))

   
    rtipsy,root[i],header,gas,dark,star
    print,"Finished reading "+root[i]

    ;if skid outputs exist then read in .grp and .stat files
    IF KEYWORD_SET(skid) THEN BEGIN
        haloind=read_ascii_array(root[i]+'.grp',/LONG)
        halo=read_stat(root[i]+'.stat')
        ;put lengths into physical units for later use.
        halo.x=halo.x*lengthunit
        halo.y=halo.y*lengthunit
        halo.z=halo.z*lengthunit
    ENDIF

    history[*].expfactor[i]=header.time
    ;now work out redshift and lookback time
    currentz=(1./header.time)-1
    currentlookback=lookback(currentz)
    history[*].redshift[i]=currentz
    history[*].lookback[i]=currentlookback
    ;find the timestep from the file name
    currentstep=LONG(STRMID(root[i],4,/REVERSE_OFFSET))
    history[*].timestep[i]=currentstep
    densityunit=dRhoUnit/(header.time^3.0)
    lengthunit=dKpcUnit*header.time

    ;IF writemark is set, open write mark file
    IF KEYWORD_SET(writemark) THEN BEGIN
        markfile=root[i]+'.mark'
        OPENW, 1, markfile
        outstring=STRTRIM(header.n,2)+' '+STRTRIM(header.ngas,2)+' '+STRTRIM(header.nstar,2)
        printf, 1, outstring
    ENDIF

    ;For each particle,
    print,"Matching particles from origiord with "+root[i]
    FOR j = 0l, n_elements(orig_iord) - 1 DO BEGIN
        ;Look for the particle's original index in the iord array.
        currentind=WHERE((orig_iord[j] EQ iord), ncurrentind)

        ;If it is not found, then it must be a
        ;star at timestep before it was born
        IF (ncurrentind EQ 0) THEN currentind=WHERE(orig_igasorder[j] EQ iord)

        ;Now fill up the history structure
        if ( currentind gt (header.ngas+header.ndark -1) ) then begin
            history[j].type[i]='star'
            history[j].mass[i]=star[currentind - ( header.ngas+header.ndark )].mass
            history[j].x[i]=star[currentind - ( header.ngas+header.ndark )].x
            history[j].y[i]=star[currentind - ( header.ngas+header.ndark )].y
            history[j].z[i]=star[currentind - ( header.ngas+header.ndark )].z
            history[j].vx[i]=star[currentind - ( header.ngas+header.ndark)].vx
            history[j].vy[i]=star[currentind - ( header.ngas+header.ndark )].vy
            history[j].vz[i]=star[currentind - ( header.ngas+header.ndark )].vz
            history[j].rho[i]=-1
            history[j].temp[i]=-1
            history[j].metallicity[i]=star[currentind - ( header.ngas+header.ndark )].metals
        endif else if (currentind gt (header.ngas-1) ) then begin
            history[j].type[i]='dark'
            history[j].mass[i]=dark[currentind - header.ngas ].mass
            history[j].x[i]=dark[currentind - header.ngas  ].x
            history[j].y[i]=dark[currentind - header.ngas].y
            history[j].z[i]=dark[currentind - header.ngas].z
            history[j].vx[i]=dark[currentind - header.ngas].vx
            history[j].vy[i]=dark[currentind - header.ngas].vy
            history[j].vz[i]=dark[currentind - header.ngas].vz
            history[j].rho[i]=-1
            history[j].temp[i]=-1
            history[j].metallicity[i]=-1
        endif else begin
            history[j].type[i]='gas'
            history[j].mass[i]=gas[currentind ].mass
            history[j].x[i]=gas[currentind ].x
            history[j].y[i]=gas[currentind ].y
            history[j].z[i]=gas[currentind ].z
            history[j].vx[i]=gas[currentind ].vx
            history[j].vy[i]=gas[currentind ].vy
            history[j].vz[i]=gas[currentind ].vz
            history[j].rho[i]=gas[currentind ].dens
            history[j].temp[i]=gas[currentind ].tempg
            history[j].metallicity[i]=gas[currentind ].zmetal
        endelse


        IF KEYWORD_SET(skid) THEN BEGIN
            ;find out which halo it is in
            haloid=haloind[currentind]
            history[j].haloid[i]=haloid                
            ;only calculate other quantities if particle is in halo
            IF (haloid GT 0) THEN BEGIN
                ;to get proper halo, the index haloid-1
                history[j].halovcmax[i]=halo.vcmax[haloid-1]*velunit
                distance=SQRT((history[j].x[i]-halo.x[haloid-1])^2+(history[j].y[i]-halo.y[haloid-1])^2+(history[j].z[i]-halo.z[haloid-1])^2)
                history[j].halodistance[i]=distance
            ENDIF
        ENDIF

        IF KEYWORD_SET(writemark) THEN BEGIN
            ; For now only print stars.  Change if you want gas
            IF (history[j].type[i] EQ 'star') THEN printf, 1, STRTRIM(currentind+1,2)
        ENDIF

    ENDFOR
    ;close the mark file for this timestep
    IF KEYWORD_SET(writemark) THEN BEGIN
        close, 1
        free_lun, 1
    ENDIF
ENDFOR

save,/compress,filename='find_history.sav'
RETURN, history

END
