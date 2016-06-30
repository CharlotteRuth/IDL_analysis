PRO writecol,name,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15, $
            v16,v17,v18,v19,v20,v21,v22,v23,v24,v25,v26,v27,v28,v29,$
             v30,v31,v32,v33,v34,v35,v36,v37,v38,v39,v40,v41,v42,v43,v44,v45,$
             format=format,update=update
close,11
if (keyword_set(update)) then openw,11,name,/append else openw,11,name
ncol = N_params() - 1           ;Number of columns of data expected
if (size(v1,/n_dimension) ge 2) then nel=n_elements(v1[0,*]) $
else nel=n_elements(v1)
vv = 'v' + strtrim( indgen(ncol)+1, 2)
cmd='for i=0L,nel-1 do printf,11'
i=0L
for j=0L,ncol-1 do begin
	r=execute('dim=size('+vv[j]+',/dimension)')
 	if (n_elements(dim) ge 2) then $
		for k=0,dim[0]-1 do cmd=cmd+','+vv[j]+'['+strtrim(string(k),2)+',i]' $
	else cmd=cmd+','+vv[j]+'[i]'
endfor
if (keyword_set(format)) then cmd=cmd+',format=format'
r=execute(cmd)
close,11
end

