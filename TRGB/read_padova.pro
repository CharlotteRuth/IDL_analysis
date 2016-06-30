FUNCTION READ_PADOVA,INTERPOL=interpol
; DESCRIPTION: Reads in all padova isochrone data into structure
; including ACS, WFPC2, Sloan, 2MASS, and Johnson magnitudes 
; NOTES: 
; 1) Creates a save file, padova.idl in the basedir.  If it is present
; already, This file is used rather than reading in the isochrones!
; 2) INTERPOL keyword should be set to a file downloaded from Leo's
; website after running proc_interp.pl (which adds a Z column to the data).

basedir='/Users/seth/astro/padova/'
;basedir='/pool/nebraska1/isochrones/padova/'
baseacs=basedir+'isoc_acs_wfc/'
base2mass=basedir+'isoc_2mass/'
base2sloan=basedir+'isoc_sloan/'
basewfpc2=basedir+'isoc_wfpc2vega/'

metaltag=['0001','0004','001','004','008','019']
ztag=[.0001,.0004,.001,.004,.008,0.019]
nmetals=N_ELEMENTS(ztag)

;to read in (ONLY) a set of interpolated ACS tracks
IF KEYWORD_SET(INTERPOL) THEN BEGIN    
    infile=baseacs+interpol
    IF (FILE_TEST(infile) EQ 0) THEN infile=baseacs+'isoc_10Gyr_proc.dat'
    READCOL,infile,z,logage,m_ini,m_act,f475mag,f606mag,f814mag,skip=2,FORMAT='F,F,F,F,X,X,X,X,X,F,X,X,F,X,X,X,X,F',/SILENT
    padova=REPLICATE({z:0.0, logage:0.0, m_ini:0.0, m_act:0.0,umag:0.0,bmag:0.0,vmag:0.0,rmag:0.0,imag:0.0,jmag:0.0,hmag:0.0,kmag:0.0,f475mag:0.0,f606mag:0.0,f814mag:0.0,j2mass:0.0,h2mass:0.0,k2mass:0.0,w300mag:0.0,w439mag:0.0,w555mag:0.0,w606mag:0.0,w814mag:0.0},N_ELEMENTS(z))
;    padova=REPLICATE({z:0.0, logage:0.0, m_ini:0.0, m_act:0.0,f475mag:0.0,f606mag:0.0,f814mag:0.0},N_ELEMENTS(z))
    padova.z=z
    padova.logage=logage 
    padova.m_ini=m_ini   
    padova.m_act=m_act   
    padova.f475mag=f475mag  
    padova.f606mag=f606mag  
    padova.f814mag=f814mag  
    return, padova
ENDIF ;interpolated tracks

savefile=basedir+'padova.idl'
IF (FILE_TEST(savefile)) THEN BEGIN
    restore, savefile
    return, padova
ENDIF

totallines=0l
FOR i=0,nmetals-1 DO BEGIN
    infile=basedir+'isocz'+metaltag[i]+'.dat'
    totallines=totallines+NUMLINES(infile)-2l
ENDFOR


padova=REPLICATE({z:0.0, logage:0.0, m_ini:0.0, m_act:0.0,umag:0.0,bmag:0.0,vmag:0.0,rmag:0.0,imag:0.0,jmag:0.0,hmag:0.0,kmag:0.0,f475mag:0.0,f606mag:0.0,f814mag:0.0,j2mass:0.0,h2mass:0.0,k2mass:0.0,usloan:0.0,gsloan:0.0,rsloan:0.0,isloan:0.0,zsloan:0.0,w300mag:0.0,w439mag:0.0,w555mag:0.0,w606mag:0.0,w814mag:0.0},totallines)

counter=0
FOR i=0,nmetals-1 DO BEGIN
    infile=basedir+'isocz'+metaltag[i]+'.dat'
    acsfile=baseacs+'isoc_z'+metaltag[i]+'.dat'
    massfile=base2mass+'isoc_z'+metaltag[i]+'.dat'
    sloanfile=base2sloan+'isoc_z'+metaltag[i]+'.dat'
    wfpc2file=basewfpc2+'isoc_z'+metaltag[i]+'.dat'
    print,'reading in data isoc_z',metaltag[i],'.dat'
    READCOL, infile, logage,m_ini,m_act,umag,bmag,vmag,rmag,imag,jmag,hmag,kmag,FORMAT='F,F,F,X,X,X,X,F,F,F,F,F,F,F,F',skip=2,/SILENT
    READCOL, acsfile,f475mag,f606mag,f814mag,skip=2,FORMAT='X,X,X,X,X,X,X,X,F,X,X,F,X,X,X,X,F',/SILENT
    READCOL, massfile,j2mass,h2mass,k2mass,skip=2,FORMAT='X,X,X,X,X,X,X,F,F,F',/SILENT
    READCOL,sloanfile,usloan,gsloan,rsloan,isloan,zsloan,FORMAT='X,X,X,X,X,X,X,F,F,F,F,F',/SILENT
    READCOL,wfpc2file,w300mag,w439mag,w555mag,w606mag,w814mag,FORMAT='X,X,X,X,X,X,X,X,X,X,F,X,F,X,F,F,X,F',/SILENT

    npoints=N_ELEMENTS(logage)
    z=REPLICATE(ztag[i],npoints)
    mini=counter
    maxi=counter+npoints-1
    padova[mini:maxi].logage=logage     
    padova[mini:maxi].m_ini=m_ini       
    padova[mini:maxi].m_act=m_act       
    padova[mini:maxi].umag=umag         
    padova[mini:maxi].bmag=bmag         
    padova[mini:maxi].vmag=vmag         
    padova[mini:maxi].rmag=rmag         
    padova[mini:maxi].imag=imag         
    padova[mini:maxi].jmag=jmag         
    padova[mini:maxi].hmag=hmag         
    padova[mini:maxi].kmag=kmag         
    padova[mini:maxi].f475mag=f475mag   
    padova[mini:maxi].f606mag=f606mag   
    padova[mini:maxi].f814mag=f814mag   
    padova[mini:maxi].j2mass=j2mass
    padova[mini:maxi].h2mass=h2mass
    padova[mini:maxi].k2mass=k2mass
    padova[mini:maxi].usloan=usloan          
    padova[mini:maxi].gsloan=gsloan          
    padova[mini:maxi].rsloan=rsloan         
    padova[mini:maxi].isloan=isloan          
    padova[mini:maxi].zsloan=zsloan
    padova[mini:maxi].w300mag=w300mag
    padova[mini:maxi].w439mag=w439mag
    padova[mini:maxi].w555mag=w555mag
    padova[mini:maxi].w606mag=w606mag
    padova[mini:maxi].w814mag=w814mag    
    padova[mini:maxi].z=z               
    counter=counter+npoints
ENDFOR
save, padova, file=savefile
return,padova

END













