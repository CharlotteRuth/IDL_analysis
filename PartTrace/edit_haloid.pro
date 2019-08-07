;Uses the output of mergertree.pro (*halo_step.out) to edit the
;*haloid.dat file listing the main halo at each of the steps

PRO edit_haloid,filebase,finalid

openr, lun, filebase+'.grp'+finalid+'.halo_step.out', /GET_LUN
openw, lun2, filebase+'.grp'+finalid+'.haloid.dat',/GET_LUN 
line = ' '
previous_id = finalid
readf, lun, line
split = strsplit(line,' ',/extract)
previous_step = fix(split[2])
current_step = fix(split[1])
CASE 1 OF
    (previous_step GE 1000) AND (previous_step LT 10000): previous_step_str = '00'+strtrim(previous_step,2)
    (previous_step GT 100) AND (previous_step LT 1000): previous_step_str = '000'+strtrim(previous_step,2)
    (previous_step GT 10) AND (previous_step LT 100): previous_step_str = '0000'+strtrim(previous_step,2)
    (previous_step LT 10): previous_step_str = '0000'+strtrim(previous_step,2)
    ELSE: print,'Previous step, '+strtrim(previous_step,2)+', is not within range'
ENDCASE
print,filebase+'.'+previous_step_str+'/'+filebase+'.'+previous_step_str,' ',previous_id
printf,lun2,filebase+'.'+previous_step_str+'/'+filebase+'.'+previous_step_str,' ',previous_id,format='(80A,4A)'
n_previous_part = 0

WHILE NOT EOF(lun) DO BEGIN
    readf, lun, line
    split = strsplit(line,' ',/extract)
    IF n_elements(split) EQ 5 THEN BEGIN
        previous_step = fix(split[2])
        IF current_id NE '0' THEN BEGIN
            previous_id = current_id
            current_id = '0'
            IF previous_step NE current_step THEN print,'Missing a step: ', current_step
            current_step = fix(split[1])
            CASE 1 OF
                (previous_step GE 1000) AND (previous_step LT 10000): previous_step_str = '00'+strtrim(previous_step,2)
                (previous_step GT 100) AND (previous_step LT 1000): previous_step_str = '000'+strtrim(previous_step,2)
                (previous_step GT 10) AND (previous_step LT 100): previous_step_str = '0000'+strtrim(previous_step,2)
                (previous_step LT 10): previous_step_str = '0000'+strtrim(previous_step,2)
                ELSE: print,'Previous step, '+strtrim(previous_step,2)+', is not within range'
            ENDCASE
            print,filebase+'.'+previous_step_str+'/'+filebase+'.'+previous_step_str,' ',previous_id
            printf,lun2,filebase+'.'+previous_step_str+'/'+filebase+'.'+previous_step_str,' ',previous_id,format='(80A,4A)'
            n_previous_part = 0
        ENDIF ELSE print,'Untracked halo'
    ENDIF ELSE BEGIN
        IF n_elements(split) EQ 4 THEN BEGIN
;            print,split[1],' ',previous_id,long(split[3]),n_previous_part
            IF (split[1] EQ previous_id) AND (long(split[3]) GT n_previous_part) THEN BEGIN
                current_id = split[0]
                n_previous_part = long(split[3])
;                print,'valid'
            ENDIF
        ENDIF
    ENDELSE
ENDWHILE
IF current_id NE '0' THEN BEGIN
    CASE 1 OF
        (current_step GE 1000) AND (current_step LT 10000): current_step_str = '00'+strtrim(current_step,2)
        (current_step GT 100) AND (current_step LT 1000): current_step_str = '000'+strtrim(current_step,2)
        (current_step GT 10) AND (current_step LT 100): current_step_str = '0000'+strtrim(current_step,2)
        (current_step LT 10): current_step_str = '0000'+strtrim(current_step,2)
        ELSE: print,'Current step, '+strtrim(current_step,2)+', is not within range'
    ENDCASE
    print,filebase+'.'+current_step_str+'/'+filebase+'.'+current_step_str,' ',current_id
    printf,lun2,filebase+'.'+current_step_str+'/'+filebase+'.'+current_step_str,' ',current_id,format='(80A,4A)'
ENDIF ELSE print,'Untracked Halo'
close,lun
close,lun2
free_lun,lun,lun2
END
