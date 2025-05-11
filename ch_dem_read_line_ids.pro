
FUNCTION ch_dem_read_line_ids, filename, int_file=int_file,  $
                               spec_gauss_int_file=spec_gauss_int_file, $
                               dwvl_check=dwvl_check, $
                               vshift=vshift, $
                               no_pop_delete=no_pop_delete

;+
; NAME:
;      CH_DEM_READ_LINE_IDS
;
; PURPOSE:
;      Read an emission line ID file in the format used by the CH_DEM
;      routines. 
;
; CATEGORY:
;      CHIANTI; differential emission measure; DEM.
;
; CALLING SEQUENCE:
;	Result = CH_DEM_READ_LINE_IDS( Filename )
;
; INPUTS:
;      Filename:  The name of an ID file to read. The Fortran format
;                 for the file is (a5,f10.0,a40) for ion name (CHIANTI
;                 format), wavelength, and transition information (see
;                 OUTPUTS section below).
;
; OPTIONAL INPUTS:
;      Int_File:  The name of a file containing line intensities. The
;                 file should have three columns for observed
;                 wavelength, intensity and intensity error. A free
;                 format can be used. Note that if an entry in
;                 FILENAME does not have an intensity, then it will be
;                 removed from the output. Matches between the two
;                 files are made by finding the lines in each file
;                 that are closest in wavelength, as long as the
;                 separation is less than 0.5 angstroms.
;      Spec_Gauss_Int_File:  The routine SPEC_GAUSS_WIDGET returns a
;                 text file containing line fit parameters in a
;                 standard format. This file can be input to the
;                 present routine with this keyword.
;      Dwvl_Check: When making matches in the line intensity files, the
;                 routine checks for lines that are within +/-
;                 DWVL_CHECK angstroms of the wavelength in the line ID
;                 file. The default is dwvl_check=0.5. If multiple
;                 lines are found, then the line nearest in wavelength
;                 is selected.
;      Vshift:    This is used if there's systematic Doppler
;                 shift of lines in the observed spectrum, which
;                 potentially affects the mapping between the lines in
;                 the intensity file and those in the ID file. For
;                 example, if the observed spectrum lines show a
;                 blueshift of -50 km/s, then set VSHIFT=-50. 
;
; KEYWORD PARAMETERS:
;      NO_POP_DELETE: If set, then the file pop_lookup_line_list.txt
;                 will not be deleted (this the software will use the
;                 existing file). You should only use this if you are
;                 sure the line list has not changed. For
;                 example, if you're recursively running the
;                 software on the same data-set.
;
; OUTPUTS:
;      A structure with the tags:
;        .ion   Ion name (CHIANTI format)
;        .wvl   Wavelength of line
;        .index A string identifying the CHIANTI transition indices
;               (see below).
;        .obs_wvl  Observed wavelength (only if INT_FILE specified).
;        .int   Intensity (only if INT_FILE specified).
;        .err   Intensity error (only if INT_FILE specified).
;        .model_int  The intensity computed from the model (inserted
;               by ch_dem_gauss_fit).
;        .logt_eff  The effective (log) temperature of the line
;               (inserted by ch_dem_gauss_fit).
;        .logt_max  The temperature at which the contribution
;               function peaks (inserted by ch_dem_add_contrib).
;        .init_ab_ind Integer that gets populated during the DEM
;               minimization process.
;        .ab_type   Integer that gets populated during the DEM
;               minimization process
;        .ab_ind Integer that gives the element index (populated
;                during DEM minimization process).
;        .abund  Float that will contain the element abundance after
;                the DEM minimization process.
;
;      The transition indices are of the form '1-27'. If the line
;      consists of self-blends, then they are specified as, e.g.,
;      '1-27,2-28,1-26'. The routine CH_DEM_PROCESS_INDEX_STRING
;      processes this string to extract the levels into integers.
;
;      Note that the output structure may not contain all the lines
;      listed in FILENAME. If an intensity file was specified (e.g.,
;      with INT_FILE) then only those transitions in FILENAME that
;      have an observed intensity will be included. If the intensity
;      file was not specified, then all lines in FILENAME are
;      included. 
;
; MODIFICATION HISTORY:
;      Ver.1, 17-Jul-2018, Peter Young.
;      Ver.2, 23-Jul-2018, Peter Young.
;         Added model_int and logt_eff to the output. These are filled
;         by ch_dem_gauss_fit after the fitting has been performed.
;      Ver.3, 5-Feb-2019, Peter Young
;         Added SPEC_GAUSS_INT_FILE and DWVL_CHECK optional inputs.
;      Ver.4, 12-Feb-2019, Peter Young
;         Added init_ab_ind (see routine ch_dem_init_abund.pro) and
;         ab_type to output.
;      Ver.5, 13-Feb-2019, Peter Young
;         Added abund to output.
;      Ver.6, 20-Jul-2019, Peter Young
;         Fixed error in how int_file is read.
;      Ver.7, 22-Jul-2019, Peter Young
;         Changed wvl_check to dwvl_check and added it to output
;         structure; added check on n_params.
;      Ver.8, 24-Jul-2019, Peter Young
;         Added a check to make sure a line intensity doesn't
;         get assigned to two different transitions; added ab_ind tag
;         to output structure.
;      Ver.9, 26-Jul-2019, Peter Young
;         Added BLEND_FILE input that tells the routine to add
;         together separate lines within the FILENAME file.
;      Ver.10, 2-Aug-2019, Peter Young
;         Added VSHIFT optional input.
;      Ver.11, 14-Aug-2019, Peter Young
;         Now prints wavelength check list and warning messages;
;         removed the blend_file optional input.
;      Ver.12, 20-Aug-2019, Peter Young
;         Now deletes the pop_lookup_line_list.txt file (if it
;         exists), forcing the user to re-generate it using
;         ch_dem_write_lookup_tables.
;      Ver.13, 13-Dec-2019, Peter Young
;         Don't check velocities if the line intensities weren't
;         specified.
;      Ver.14, 19-May-2020, Peter Young
;         Added /no_pop_delete keyword; added logt_max to output structure.
;      Ver.15, 11-May-2025, Peter Young
;         Added logt_max to output structure.
;-


IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> output = ch_dem_read_line_ids( filename [, spec_gauss_int_file= '
  print,'                                          int_file=, dwvl_chck=, vshift= '
  print,'                                          /no_pop_delete ] )'
  return,-1
ENDIF 

IF n_elements(dwvl_check) EQ 0 THEN dwvl_check=0.5

chck=file_search(filename,count=count)
IF count EQ 0 THEN BEGIN
  print,'% CH_DEM_READ_LINE_IDS: file not found. Returning...'
  return,-1
ENDIF

;
; Create the tags for the LINE_DATA output. Note that contrib is a
; list as we don't know how big the array is at this point.
;
str={ion: '', $
     element_num: 0, $
     wvl: 0., $
     index: '', $
     obs_wvl: -1., $
     int: -1., $
     err: -1., $
     model_int: -1., $    ; to be filled after DEM computed
     logt_eff: -1., $     ; to be filled after DEM computed
     logt_max: -1., $     ; to be filled after contrib. fn. computed.
     ab_ind: -1, $
     dwvl_check: dwvl_check, $
     blend_tag: '', $
     label: ''}
output=0
a=''
b=''
c=0.
d=''

openr,lin,filename,/get_lun
WHILE eof(lin) NE 1 DO BEGIN
  readf,lin,format='(a6,a2,f10.0,a40)',a,b,c,d
  str.ion=trim(a)
  str.blend_tag=trim(b)
  convertname,str.ion,iz
  str.element_num=iz
  str.wvl=c
  str.label=trim(string(format='(f12.2)',str.wvl))
  str.index=trim(d)
  IF n_tags(output) EQ 0 THEN output=str ELSE output=[output,str]
ENDWHILE
free_lun,lin


;
; Do a check on 'blend_tag' in case two lines from different species
; have been flagged.
;
k=where(output.blend_tag NE '',nk)
IF nk NE 0 THEN BEGIN
  blend_tag=output[k].blend_tag
  tags=blend_tag[uniq(blend_tag,sort(blend_tag))]
  nt=n_elements(tags)
  FOR i=0,nt-1 DO BEGIN
    k=where(output.blend_tag EQ tags[i])
    ions=output[k].ion
    ions=ions[uniq(ions,sort(ions))]
    IF n_elements(ions) GT 1 THEN BEGIN
      print,'% CH_DEM_READ_LINE_IDS: a blend tag has been assigned to two different ion species. This is not'
      print,'                        allowed. Please revise the LINE_LIST file. Returning...'
      print,format='("                        blend_tag:",a2,"  ions: ",6a6)',tags[i],ions
      return,-1
    ENDIF 
  ENDFOR
ENDIF


;
;
; If int_file was specified, then we load it into 'int_struc', which has
; the same tags as the output from read_line_fits (see below).
;
str={wvl: 0., int: 0., sint: 0.}
IF keyword_set(int_file) THEN BEGIN
  int_struc=0
  openr,lin,int_file,/get_lun
  WHILE eof(lin) NE 1 DO BEGIN
    readf,lin,str
    IF n_tags(int_struc) EQ 0 THEN int_struc=str ELSE int_struc=[int_struc,str]
  ENDWHILE
  free_lun,lin
 ;
 ; This handles the special case of sint (error on intensity) being
 ; -1. In this case, it is assumed that int is an upper limit, and
 ; that the intensity is actually int/2 +/- int/2.
 ; For example, if int=6, then the intensity used for the minimization
 ; process will be 3 +/- 3.
 ; 
  k=where(int_struc.sint EQ -1.,nk)
  IF nk NE 0 THEN BEGIN
    FOR i=0,nk-1 DO BEGIN
      int_struc[k[i]].int=int_struc[k[i]].int/2.
      int_struc[k[i]].sint=int_struc[k[i]].int
    ENDFOR 
  ENDIF 
ENDIF


;
; This handles the spec_gauss line fit parameter files.
; If int_struc already exists, then I merge the line list from the two
; files. 
;
IF keyword_set(spec_gauss_int_file) THEN BEGIN
  read_line_fits,spec_gauss_int_file,sg_struc
  n=n_elements(sg_struc)
  m=n_elements(int_struc)
  int_struc2=replicate(str,n+m)
  int_struc2[0:n-1].wvl=sg_struc.wvl
  int_struc2[0:n-1].int=sg_struc.int
  int_struc2[0:n-1].sint=sg_struc.sint
 ;
  IF m NE 0 THEN BEGIN
    int_struc2[n:n+m-1].wvl=int_struc.wvl
    int_struc2[n:n+m-1].int=int_struc.int
    int_struc2[n:n+m-1].sint=int_struc.sint
  ENDIF
  junk=temporary(int_struc)
  int_struc=temporary(int_struc2)
ENDIF


;
; If int_struc has been created, then assign intensities to lines in
; 'output'.
;
; I use sg_flag to check if a line intensity has already been assigned
; to another transition. This can happen if two lines are very close
; together or if there's an error in one of the input files.
; However, if a line is a blend of two different species (e.g., Si VII
; 276.85+Si VIII 276.86) and they're flagged with the same
; 'blend_tag' then it's OK to assign the intensity to both
; lines as this is handled later in the software.
;
flag_warning=0b
IF n_tags(int_struc) NE 0 THEN BEGIN
  sg_flag=make_array(n_elements(int_struc),/integer,value=-1)
  n=n_elements(output)
  FOR i=0,n-1 DO BEGIN
    IF n_elements(vshift) NE 0 THEN shft=v2lamb(vshift,output[i].wvl) ELSE shft=0
    getmin=min(abs(output[i].wvl+shft-int_struc.wvl),imin)
    IF getmin LT dwvl_check THEN BEGIN
      output[i].obs_wvl=int_struc[imin].wvl
      output[i].int=int_struc[imin].int
      output[i].err=int_struc[imin].sint
      IF sg_flag[imin] GT -1 THEN BEGIN
        flag_warning=flag_warning+1b
        print,format='("WARNING: The same intensity has been assigned to ",a6,f8.2)',output[i].ion,output[i].wvl
        print,format='("                                             and ",a6,f8.2)',output[sg_flag[imin]].ion,output[sg_flag[imin]].wvl
      ENDIF 
      sg_flag[imin]=i
    ENDIF 
  ENDFOR
  k=where(output.obs_wvl GT 0.,nk)
  IF n NE nk THEN BEGIN
    print,'% CH_DEM_READ_LINE_IDS: '+trim(n-nk)+' of '+trim(n)+' lines do not have observed intensities.'
    print,'                        These lines were not placed in the output structure.'
  ENDIF 
  output=output[k]
ENDIF 



;
; Do a wavelength check (only if the line intensity file was
; specified, though).
;
IF n_tags(int_struc) NE 0 THEN BEGIN 
  nl=n_elements(output)
  i=sort(output.wvl)
  outputx=output[i]
  print,' Ion   Ref. wvl   Obs. wvl    Velocity'
  FOR i=0,nl-1 DO BEGIN
    v=lamb2v(outputx[i].obs_wvl-outputx[i].wvl,outputx[i].wvl)
    print,format='(a5,2f10.3,f8.1," km/s")',outputx[i].ion,outputx[i].wvl,outputx[i].obs_wvl,v
  ENDFOR
  print,'* check the velocities above as unusual values may indicate lines have been incorrectly identified.'
 ;
  IF keyword_set(flag_warning) THEN BEGIN
    print,'* WARNING: At least one intensity was assigned to two different lines in the line list.'
    print,'           See output above for more details.'
    print,'           Try using a smaller value of DWVL_CHECK (current value: '+trim(dwvl_check)+' angstroms).'
    print,'           Alternatively, check your input files.'
  ENDIF
ENDIF 

junk=temporary(outputx)

;
; Problems may occur if I read a modified line list but the ch_dem
; software still uses the previous pop_lookup_line_list.txt file
; (created by ch_dem_write_lookup_tables). I thus make sure to delete
; this file (if it exists) every time ch_dem_read_line_ids is
; called. To create the file again, you need to run
; ch_dem_write_lookup_tables.
;
; PRY, 19-May-2020: however, if you're recursively running the
; ch_dem software for the same data-set (e.g., each pixel in a
; raster), then you can safely keep this file. Therefore, I have the
; keyword /no_pop_delete for this situation. 
;
IF NOT keyword_set(no_pop_delete) THEN BEGIN 
  chck=file_info('pop_lookup_line_list.txt')
  IF chck.exists EQ 1 THEN file_delete,'pop_lookup_line_list.txt'
ENDIF 

return,output

END
