;This program will read in .stat files for different runs and attempt
;to match the halos using position and mass


function distance,pos
;quick function to return the distance of a point from the origin
return, SQRT(pos[0]^2 + pos[0]^2 + pos[0]^2)
END


pro match_halos
mindis = 0.5
mass_range = 5.00

stat_file_1 = '/astro/net/scratch2/fabio/Xruns/h1201-h1755-X5X5g1bwK/OLD_AMIGA/h1201-h1755-X5X5g1bwK.00512.amiga.stat'
;stat_file_4 = '/astro/net/scratch2/fabio/Xruns/h1201-h1755-X2X2g3bwK/h1201-h1755-X2X2g3bwK.00512.amiga.stat'
;stat_file_3 = '/astro/net/scratch2/fabio/Xruns/h1201-h1755-X3X3g1bwK/h1201-h1755-X3X3g1bwK.00512.amiga.stat'
stat_file_2 = '/astro/net/scratch2/fabio/Xruns/h1201-h1755-X2X2g1bwK/OLD_AMIGA/h1201-h1755-X2X2g1bwK.00512.amiga.stat'

input_halos_1 = READ_STAT_STRUC_AMIGA(stat_file_1)
input_halos_2 = READ_STAT_STRUC_AMIGA(stat_file_2)
;input_halos_3 = READ_STAT_STRUC_AMIGA(stat_file_3)
;input_halos_4 = READ_STAT_STRUC_AMIGA(stat_file_4)
nhalos_1 = n_elements(input_halos_1)

close,1
openw,1,'halo_matches_new_full1.txt'
openw,2,'halo_matches_X5X5.txt'
openw,3,'halo_matches_X2X2.txt'
;matches = fltarr(10,nhalos_1)
matches = fltarr(12,nhalos_1)
matches1 = fltarr(5,nhalos_1)
matches2 = fltarr(5,nhalos_1)

FOR ct = 0, MIN([N_ELEMENTS(INPUT_HALOS_1),N_ELEMENTS(INPUT_HALOS_2) ])- 1 DO BEGIN
;    print,' '
;    print,'* Halo Number ',input_halos_1[ct].group,"*"
    matches[0,ct] = input_halos_1[ct].group
    matches1[0,ct] = input_halos_1[ct].group
    mb1 = input_halos_1.m_tot - input_halos_1.m_gas - input_halos_1.m_star
    mb2 = input_halos_2.m_tot - input_halos_2.m_gas - input_halos_2.m_star
    distances = SQRT((input_halos_2.xc - input_halos_1[ct].xc)^2 + (input_halos_2.yc - input_halos_1[ct].yc)^2 + (input_halos_2.zc - input_halos_1[ct].zc)^2)
    match = where(mb2 le (mb1*(1.0 + mass_range)) AND mb2 ge (mb1*(1.0 - mass_range)) AND distances le mindis)
    IF (match[0] NE -1) THEN BEGIN
        fracmass = input_halos_2[match].m_tot/input_halos_1[ct].m_tot
        print,input_halos_1[ct].group," Matches at :",input_halos_2[match].group
        print,"Fractional Mass: ", input_halos_2[match].m_tot/input_halos_1[ct].m_tot
        print,"Distance       : ",distances[match]
  ;      print,input_halos_1[ct].group,input_halos_2[match].group
;        print,distances[match],fracmass
;        temp_min = MIN(distances[match],min_ind)
        temp_min = MIN(ABS(1.0 - fracmass),min_ind)
        match = match[min_ind]
        matches[1,ct] = input_halos_1[ct].xc
        matches[2,ct] = input_halos_1[ct].yc
        matches[3,ct] = input_halos_1[ct].zc
        matches[4,ct] = input_halos_1[ct].RVIR
        matches[5,ct] = input_halos_2[match].group
        matches[6,ct] = input_halos_2[match].xc
        matches[7,ct] = input_halos_2[match].yc
        matches[8,ct] = input_halos_2[match].zc
        matches[9,ct] = input_halos_2[match].RVIR
        matches[10,ct] = mb2[match]/mb1[match]
        matches[11,ct] = distances[match]

        matches1[1,ct] = (input_halos_1[ct].xc-25.0)/50
        matches1[2,ct] = (input_halos_1[ct].yc-25.0)/50
        matches1[3,ct] = (input_halos_1[ct].zc-25.0)/50
        matches1[4,ct] = input_halos_1[ct].RVIR/50000

        matches2[0,ct] = input_halos_2[match].group
        matches2[1,ct] = (input_halos_2[match].xc-25.0)/50
        matches2[2,ct] = (input_halos_2[match].yc-25.0)/50
        matches2[3,ct] = (input_halos_2[match].zc-25.0)/50
        matches2[4,ct] = input_halos_2[match].RVIR/50000
        a = [input_halos_2[0:match-1],input_halos_2[match+1:N_ELEMENTS(INPUT_HALOS_2)-1]]
        input_halos_2 = a

        printf,2,matches1[*,ct]
        printf,3,matches2[*,ct]
    ENDIF; ELSE print,"For File 2, No match"

;    distances = SQRT((input_halos_3.xc - input_halos_1[ct].xc)^2 + (input_halos_3.yc - input_halos_1[ct].yc)^2 + (input_halos_3.zc - input_halos_1[ct].zc)^2)
;    match = where(input_halos_3.m_tot le (input_halos_1[ct].m_tot*1.1) AND input_halos_3.m_tot ge (input_halos_1[ct].m_tot*0.9) AND distances le 0.5)
;    IF (match[0] NE -1) THEN BEGIN
;        print,"For File 3, Matches at :",input_halos_3[match].group
;        print,"Fractional Mass: ", input_halos_3[match].m_tot/input_halos_1[ct].m_tot
;        print,"Distance       : ",distances[match]
        
;        temp_min = MIN(distances[match],min_ind)
;        match = match[min_ind]
;        matches[4,ct] = input_halos_3[match].group
;        matches[5,ct] = input_halos_3[match].m_tot/input_halos_1[ct].m_tot
;        matches[6,ct] = distances[match]
;    ENDIF ELSE print,"For File 2, No match"

;    distances = SQRT((input_halos_4.xc - input_halos_1[ct].xc)^2 + (input_halos_4.yc - input_halos_1[ct].yc)^2 + (input_halos_4.zc - input_halos_1[ct].zc)^2)
;    match = where(input_halos_4.m_tot le (input_halos_1[ct].m_tot*1.1) AND input_halos_4.m_tot ge (input_halos_1[ct].m_tot*0.9) AND distances le 0.5)
;    IF (match[0] NE -1) THEN BEGIN
;        print,"For File 4, Matches at :",input_halos_4[match].group
;        print,"Fractional Mass: ", input_halos_4[match].m_tot/input_halos_1[ct].m_tot
;        print,"Distance       : ",distances[match]
;        
;        temp_min = MIN(distances[match],min_ind)
;        match = match[min_ind]
;        matches[7,ct] = input_halos_4[match].group
;        matches[8,ct] = input_halos_4[match].m_tot/input_halos_1[ct].m_tot
;        matches[9,ct] = distances[match]
;    ENDIF ELSE print,"For File 4, No match"
    printf,1,matches[*,ct]

ENDFOR
close,1
close,2
close,3

indicies = WHERE(matches2[0,*] ne 0)
plot_match1 = fltarr(6,N_ELEMENTS(indicies))
plot_match2 = fltarr(6,N_ELEMENTS(indicies))
plot_match1[0,*] = findgen(N_ELEMENTS(indicies))/N_ELEMENTS(indicies)*224.
plot_match2[0,*] = plot_match1[0,*]
plot_match1[1:4,*] = matches1[1:4,indicies]
plot_match2[1:4,*] = matches2[1:4,indicies]
window,0
plot,plot_match1[2,*],plot_match1[3,*],psym = 1,title = 'X5X5 -- X',xrange=[-1,1],yrange=[-1,1]
FOR ct = 0, N_ELEMENTS(indicies)-1 DO plots,plot_match1[2,ct],plot_match1[3,ct],color = plot_match1[0,ct], psym = 1

window,3
plot,plot_match2[2,*],plot_match2[3,*], psym = 1,title = 'X2X2 -- X',xrange=[-1,1],yrange=[-1,1]
FOR ct = 0, N_ELEMENTS(indicies)-1 DO plots,plot_match2[2,ct],plot_match2[3,ct],color = plot_match2[0,ct], psym = 1
stop

window,0
plot,plot_match1[1,*],plot_match1[3,*],psym = 1,title = 'X5X5 -- Y',xrange=[-1,1],yrange=[-1,1]
FOR ct = 0, N_ELEMENTS(indicies)-1 DO plots,plot_match1[1,ct],plot_match1[3,ct],color = plot_match1[0,ct], psym = 1

window,3
plot,plot_match2[1,*],plot_match2[3,*], psym = 1,title = 'X2X2 -- Y',xrange=[-1,1],yrange=[-1,1]
FOR ct = 0, N_ELEMENTS(indicies)-1 DO plots,plot_match2[1,ct],plot_match2[3,ct],color = plot_match2[0,ct], psym = 1
stop

window,0
plot,plot_match1[1,*],plot_match1[2,*],psym = 1,title = 'X5X5 -- Z',xrange=[-1,1],yrange=[-1,1]
FOR ct = 0, N_ELEMENTS(indicies)-1 DO plots,plot_match1[1,ct],plot_match1[2,ct],color = plot_match1[0,ct], psym = 1

window,3
plot,plot_match2[1,*],plot_match2[2,*], psym = 1,title = 'X2X2 -- Z',xrange=[-1,1],yrange=[-1,1]
FOR ct = 0, N_ELEMENTS(indicies)-1 DO plots,plot_match2[1,ct],plot_match2[2,ct],color = plot_match2[0,ct], psym = 1
stop
END
