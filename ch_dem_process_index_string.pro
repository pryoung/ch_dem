
FUNCTION ch_dem_process_index_string, input

;+
; NAME:
;     CH_DEM_PROCESS_INDEX_STRING
;
; PURPOSE:
;     Processes the index string in the structure created by
;     ch_dem_read_line_ids.pro.
;
; CATEGORY:
;     CHIANTI; differential emission measure (DEM); string processing.
;
; CALLING SEQUENCE:
;     Result = CH_DEM_PROCESS_INDEX_STRING( Input )
;
; INPUTS:
;     Input:  A string containing the index information. For example,
;             '1-2,1-3' is interpreted as two transitions between
;             CHIANTI level indices 1 and 2, and 1 and 3.
;
; OUTPUTS:
;     A structure of N elements, where N is the number of transitions
;     (separated by commas in the input string). The tags of the
;     structure are:
;        .lower   Index of the lower level (integer).
;        .upper   Index of the upper level (integer).
;
; EXAMPLE:
;     IDL> indstr=ch_dem_process_index_string('1-2,1-3')
;
; MODIFICATION HISTORY:
;     Ver.1, 22-Jul-2019, Peter Young
;-


a=strsplit(input,',',/extract,count=n)

str={lower: 0, upper: 0}
output=replicate(str,n)

FOR i=0,n-1 DO BEGIN
  b=strsplit(a[i],'-',/extract)
  output[i].lower=fix(b[0])
  output[i].upper=fix(b[1])
ENDFOR 

return,output

END
