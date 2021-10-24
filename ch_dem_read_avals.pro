

FUNCTION ch_dem_read_avals, filename, ion_name=ion_name, lvl1=lvl1, lvl2=lvl2, $
                            status=status, wvl=wvl, aval=aval, $
                            line_data=line_data

;+
; NAME:
;     CH_DEM_READ_AVALS
;
; PURPOSE:
;     Read a file (pop_lookup_line_list.txt) that maps index
;     information to A-values. This file is automatically created as
;     part of the ch_dem procedures. 
;
; CATEGORY:
;     CHIANTI; differential emission measure (DEM); input/output.
;
; CALLING SEQUENCE:
;     Result = CH_DEM_READ_AVALS( Filename )
;
; INPUTS:
;     Filename: The name of the file to be read. The file is created
;               by the routine ch_dem_write_avals.
;
; OPTIONAL INPUTS:
;     Ion_Name: The name of an ion (CHIANTI format). Used in
;               conjunction with LVL1 and LVL2.
;     Lvl1:     The lower level index of a transition. Used in
;               conjunction with ION_NAME and LVL2.
;     Lvl1:     The upper level index of a transition. Used in
;               conjunction with ION_NAME and LVL1.
;     Line_Data: A structure in the format returned by
;                ch_dem_read_line_ids.pro. 
;
; OUTPUTS:
;     A structure (or array of structures) with the following tags:
;       .ion   Ion name (CHIANTI format)
;       .lvl1  The lower level index of a transition
;       .lvl2  The upper level index of a transition
;       .wvl   The wavelength of a transition
;       .aval  The A-value of a transition
;
; OPTIONAL OUTPUTS:
;     Status:   Either 0 or 1. A 1 indicates the transition specified
;               by ION_NAME, LVL1 and LVL2 is present in the output
;               structure.
;     Wvl:      If ION_NAME, LVL1 and LVL2 were input, and a match was
;               found, then WVL is the wavelength of the matched
;               transition.
;     Aval:     If ION_NAME, LVL1 and LVL2 were input, and a match was
;               found, then AVAL is the A-value of the matched
;               transition.
;
; EXAMPLE:
;     IDL> output=ch_dem_read_avals('pop_lookup_line_list.txt')
;     IDL> output=ch_dem_read_avals('pop_lookup_line_list.txt',ion_name='fe_12', 
;                       lvl1=1,lvl2=27,status=status,wvl=wvl,aval=aval)
;     IDL> output=ch_dem_read_avals('pop_lookup_line_list.txt',line_data=line_data)
;
; PROGRAMMING NOTES:
;     This routine has three functions:
;
;     1. Read the 'pop_lookup_line_list.txt' file.
;     2. Check if a specified transition is in the file, returning the
;        wavelength and A-value (this is used by
;        ch_dem_write_lookup_table.pro).
;     3. Filter the output structure to contain only those transitions
;        specified by LINE_DATA (this is used by
;        ch_dem_gauss_fit.pro). 
;
; MODIFICATION HISTORY:
;     Ver.1, 17-Jul-2019, Peter Young
;-


chck=file_info(filename)
IF chck.exists EQ 0 THEN BEGIN
  status=0
  return,-1
ENDIF 

openr,lin,filename,/get_lun

str={ion: '', lvl1: 0, lvl2: 0, wvl: 0., aval: 0.}
output=0
WHILE eof(lin) NE 1 DO BEGIN 
  readf,lin,format='(a6,2i7,f15.3,e15.3)',str
  str.ion=trim(str.ion)
  IF n_tags(output) EQ 0 THEN output=str ELSE output=[output,str]
ENDWHILE 

free_lun,lin

;
; The following code is intended as a check to see if the transition
; specified by ion_name, lvl1 and lvl2 is already present in the
; A-value structure: STATUS=1 if yes, 0 if no.
;
; Note WVL and AVAL are outputs. 
;
status=0
IF n_elements(ion_name) NE 0 AND n_elements(lvl1) NE 0 AND n_elements(lvl2) NE 0 THEN BEGIN
  k=where(output.ion EQ trim(ion_name) AND output.lvl1 EQ lvl1 AND output.lvl2 EQ lvl2,nk)
  IF nk NE 0 THEN BEGIN
    status=1
    wvl=output[k[0]].wvl
    aval=output[k[0]].aval
  ENDIF 
ENDIF 

;
; If LINE_DATA has been input, then it is assumed we only want to
; return a structure for the transitions specified in LINE_DATA.
;
IF n_tags(line_data) NE 0 THEN BEGIN
  out_index=-1
  n=n_elements(line_data)
  FOR i=0,n-1 DO BEGIN
    indstr=ch_dem_process_index_string(line_data[i].index)
    ntrans=n_elements(indstr)
    FOR j=0,ntrans-1 DO BEGIN
      k=where(output.ion EQ line_data[i].ion AND output.lvl1 EQ indstr[j].lower AND output.lvl2 EQ indstr[j].upper,nk)
      IF nk NE 0 THEN out_index=[out_index,k[0]]
    ENDFOR 
  ENDFOR 
  IF n_elements(out_index) GT 1 THEN BEGIN
    out_index=out_index[1:*]
    output=output[out_index]
  ENDIF 
ENDIF 



return,output

END

