function my_read_ascii, filename, template=template, _extra=_extra

  ;; the types: skip field, byte, integer, long integer, float, double,
  ;;            complex, string
  typelist = ['skip','0b','0','0L','0.0','0d','complex(0)',"''"]

  struct = read_ascii( filename, template=template, _extra=_extra )

  tags = template.fieldnames
  typeind = template.fieldtypes
  length = n_elements(struct.(0))

  ;; don't want skipped fields
  good = where( typeind ne 0 )
  tags = tags[good]

  ;; the types
  types = typelist[typeind[good]]

  ;; this is inefficient, but I wanted to do it this way in case there is a
  ;; ridiculously long structure, too long to make a string command to create
  ;; the structure

  ;; start with a throwaway
  temp = replicate( {_temp:0}, length )

  ;; now the rest of the tags
  add_tags, temporary(temp), tags, types, data

  ;; remove the temp tag
  remove_tags, temporary(data),'_temp',data

  for i = 0L, n_elements(tags)-1 do begin
    data.(i) = struct.(i)
  endfor

  return, data

end
