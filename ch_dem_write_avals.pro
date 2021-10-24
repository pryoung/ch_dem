

PRO ch_dem_write_avals, ion_name, lvl1, lvl2, dir_lookup=dir_lookup, outfile=outfile, $
                        line_data=line_data, status=status

;+
; NAME:
;     CH_DEM_WRITE_AVALS
;
; PURPOSE:
;     Writes/appends to a file containing ion, level indices,
;     wavelength for the transitions to be used in the CHIANTI DEM
;     analysis. 
;
; CATEGORY:
;     CHIANTI; differential emission analysis (DEM); file output.
;
; CALLING SEQUENCE:
;     CH_DEM_WRITE_AVALS, Ion_Name, Lvl1, Lvl2
;
; INPUTS:
;     Ion_Name: Name of ion in CHIANTI format (e.g., 'fe_13' for Fe
;               XIII).
;     Lvl1:     CHIANTI index of transition's lower level.
;     Lvl2:     CHIANTI index of transition's upper level.
;
; OPTIONAL INPUTS:
;     Line_Data:  A single element structure in the format returned by
;                 ch_dem_read_line_ids. It is used to perform a check
;                 on the wavelength corresponding to the requested
;                 transition.
;     Dir_Lookup: If OUTFILE is not specified, then this specifies the
;                 directory where the output file is sent.
;     Outfile:    The name of the output text file. The default is
;                 pop_lookup_text_file.txt. The name can include an
;                 output sub-directory. 
;	
; OUTPUTS:
;     The transition information will be written to the output file
;     (see optional input OUTFILE). If the output file already exists,
;     then the new entry will be appended. If status is not equal to
;     one, then nothing will be written.
;     
; OPTIONAL OUTPUTS:
;     Status:  An integer with one of the values:
;               0 - no problems
;               1 - the transition wavelength is not consistent with
;                   WVL_CHECK in LINE_DATA.
;               2 - the transition lvl1-lvl2 was not found in the wgfa
;                   file. 
;
; EXAMPLE:
;     IDL> ch_dem_write_avals, 'fe_13', 1, 20
;     IDL> ch_dem_write_avals, 'fe_13', 1, 20, dir_lookup='pop_lookup'
;     IDL> ch_dem_write_avals, 'fe_13', 1, 20, line_data=line_data[0]
;
;     See the routine ch_dem_write_lookup_tables to see how this
;     routine is used.
;
; MODIFICATION HISTORY:
;     Ver.1, 22-Jul-2019, Peter Young
;-


IF n_params() LT 3 THEN BEGIN
  print,'Use:  IDL> ch_dem_write_avals, ion_name, lvl1, lvl2 [, dir_lookup=, outfile= '
  print,'                               line_data= ] )'
  return
ENDIF

;
; Note that LINE_DATA must be a single element structure.
;
IF n_tags(line_data) NE 0 AND n_elements(line_data) EQ 1 THEN BEGIN
  wvl_check=line_data.wvl
  dwvl_check=line_data.dwvl_check
ENDIF 

ch_dem_get_wvl_aval,ion_name,lvl1,lvl2,wvl,aval, status=status, /quiet, $
                    wvl_check=wvl_check, dwvl_check=dwvl_check
IF status NE 0 THEN return

IF n_elements(outfile) EQ 0 THEN BEGIN
  outfile='pop_lookup_line_list.txt'
  IF n_elements(dir_lookup) NE 0 THEN outfile=concat_dir(dir_lookup,outfile)
ENDIF 

chck=file_info(outfile)
IF chck.exists THEN BEGIN
  openu,lout,outfile,/append,/get_lun
ENDIF ELSE BEGIN
  openw,lout,outfile,/get_lun
ENDELSE
printf,lout,format='(a6,2i7,f15.3,e15.3)',ion_name,lvl1,lvl2,wvl,aval
free_lun,lout

END
