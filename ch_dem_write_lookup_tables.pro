

PRO ch_dem_write_lookup_tables, line_id_file, popdir=popdir,overwrite=overwrite, $
                                ldens_start=ldens_start, ldens_end=ldens_end, $
                                ldens_step=ldens_step, $
                                ltemp_start=ltemp_start, ltemp_end=ltemp_end, $
                                ltemp_step=ltemp_step, execute=execute, $
                                dir_lookup=dir_lookup, ioneq_file=ioneq_file

;+
; NAME:
;      CH_DEM_WRITE_LOOKUP_TABLES
;
; PURPOSE:
;      Write CHIANTI population lookup tables for a set of emission
;      lines identified in the input line identification file.
;
; CATEGORY:
;      CHIANTI; differential emission measure; DEM.
;
; CALLING SEQUENCE:
;      CH_DEM_WRITE_LOOKUP_TABLE, Line_ID_File
;
; INPUTS:
;      Line_ID_File:  The name of a text file containing emission line
;               identifications. It must be in the format that is read
;               by the routine ch_dem_read_line_ids.pro.
;
; OPTIONAL INPUTS:
;      Ldens_Start:  The starting value of the log density (units:
;                    cm^-3) range for which populations are
;                    required. Default is 8.0.
;      Ldens_End:    The ending value of the log density (units:
;                    cm^-3) range for which populations are
;                    required. Default is 12.0.
;      Ldens_Step:   The step size for log density values. Default is
;                    0.2. 
;      Ltemp_Start:  The starting value of log temperature (units:
;                    K) range for which populations are required. If
;                    not set, then the minimum value in the ioneq file
;                    is used.
;      Ltemp_End:    The ending value of log temperature (units:
;                    K) range for which populations are required. If
;                    not set, then the maximum value in the ioneq file
;                    is used.
;      Ltemp_Step:   The step size for log density values. Default is
;                    0.05.
;      Dir_Lookup:   The directory to which the lookup tables will be
;                    written. If not specified then the current
;                    working directory will be used.
;      PopDir:       **Obsolete**. Exactly the same as
;                    DIR_LOOKUP. Kept for backward compatibility purposes.
;      Ioneq_File:   The name of an ionization equilibrium file. The
;                    default is to use the CHIANTI file (!ioneq_file).
;
; KEYWORD PARAMETERS:
;      OVERWRITE:  By default the routine will not overwrite a lookup
;                  table that already exists. Setting this keyword
;                  will force an overwrite.
;      EXECUTE:    The lookup tables are only written if /EXECUTE is
;                  set. 
;
; OUTPUTS:
;      For each ion in the line ID file, creates a population lookup
;      table with populations only for the specific lines in the ID
;      file. The files will go in DIR_LOOKUP, $CH_DEM_LOOKUP or
;      $CHIANTI_LOOKUP, or the working directory depending how the
;      routine is called or what environment variables are set. See
;      CHIANTI Technical Report 22 for more details.
;
;      The file pop_lookup_line_list.txt is also created in the
;      working directory, which gives A-values for the atomic
;      transitions (this speeds up the DEM calculation).
;
;      The files are only created if the keyword /EXECUTE is set.
;
; CALLS:
;      CH_DEM_READ_LINE_IDS, CH_WRITE_POP_LOOKUP_TABLE,
;      CH_DEM_PROCESS_INDEX_STRING, CH_GET_VERSION, CH_DEM_READ_AVALS,
;      CH_DEM_WRITE_AVALS 
;
; MODIFICATION HISTORY:
;      Ver.1, 18-Jul-2018, Peter Young
;      Ver 2. 30-Jul-2018, Terry Kucera
;         Added POPDIR keyword.
;      Ver.3, 4-Feb-2019, Peter Young
;         Added /EXECUTE keyword; added a check on the CHIANTI version
;         number.
;      Ver.4, 15-Feb-2019, Peter Young
;         Added DIR_LOOKUP= keyword.
;      Ver.5, 22-Jul-2019, Peter Young
;         Added error messages for checks on wavelengths; now delete
;         the pop_lookup_line_list.txt file if it already exists;
;         added check on lookup table to see if it has all the levels
;         requested by the user.
;      Ver.6, 23-Jul-2019, Peter Young
;         Added check on CHIANTI version of existing lookup table;
;         added $CH_LOOKUP_DEM environment variable; modified output
;         messages; fixed bug related to the temperature arrays being
;         used.
;      Ver.7, 24-Jul-2019, Peter Young
;         Added 'extra_levels' to correctly handle the re-writing of
;         an existing file.
;      Ver.8, 13-Dec-2019, Peter Young
;         Routine now recognises the $CHIANTI_LOOKUP environment
;         variable. 
;         
;-


IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> ch_dem_write_lookup_tables, line_id_file [, dir_lookup=, /overwrite, '
  print,'                                       /execute, ioneq_file= ]'
  return
ENDIF 

IF n_elements(popdir) NE 0 AND n_elements(dir_lookup) EQ 0 THEN BEGIN
  dir_lookup=popdir
ENDIF

;
; The code below handles where the population lookup tables will be
; written. Priority is given to DIR_LOOKUP, then $CH_DEM_LOOKUP, then
; $CHIANTI_LOOKUP. 
;
chianti_lookup=getenv('CHIANTI_LOOKUP')
ch_dem_lookup=getenv('CH_DEM_LOOKUP')
IF n_elements(dir_lookup) EQ 0 THEN BEGIN
  IF ch_dem_lookup NE '' THEN BEGIN
    dir_lookup=ch_dem_lookup
  ENDIF ELSE BEGIN
    IF chianti_lookup NE '' THEN dir_lookup=chianti_lookup
  ENDELSE 
ENDIF
;
IF n_elements(dir_lookup) NE 0 THEN BEGIN
  chck=file_info(dir_lookup)
  IF NOT chck.exists THEN file_mkdir,dir_lookup
ENDIF


;
; If LINE_ID_FILE is not a string, then assume it is the structure
; created by ch_dem_read_line_ids. 
;
IF datatype(line_id_file) EQ 'STR' THEN BEGIN
  data=ch_dem_read_line_ids(line_id_file)
ENDIF ELSE BEGIN
  data=line_id_file
ENDELSE 

k=sort(data.ion)
data=data[k]

ions=data[uniq(data.ion)].ion
nions=n_elements(ions)

;
; If the ltemp parameters are not input, then I have to be careful
; when calling ch_write_pop_lookup_table as the ltemp parameters will
; be returned for the first ion, and then applied to all subsequent
; ions. I thus need to make use of ltmp parameters to stop this happening.
;
IF n_elements(ltemp_start) NE 0 THEN ltmp_start=ltemp_start
IF n_elements(ltemp_end) NE 0 THEN ltmp_end=ltemp_end
IF n_elements(ltemp_step) NE 0 THEN ltmp_step=ltemp_step


;
; The code below forces the user to think about the temperature and
; density ranges to be used for the calculation. Also checks the
; CHIANTI version number. 
;
IF not keyword_set(execute) THEN BEGIN
 ;
 ; Check the user's CHIANTI version to make sure it's up-to-date.
 ;
  chver=ch_get_version()
  print,'CHIANTI version number: ',chver
  IF n_elements(dir_lookup) EQ 0 THEN outdir=' (current working directory)' ELSE outdir=expand_path(dir_lookup)
  print,'Output directory for tables: '+outdir
  chver2=ch_get_version(/ssw)
  IF chver2 NE '' AND chver2 NE chver THEN BEGIN
    print,'% CH_DEM_WRITE_LOOKUP_TABLES: WARNING - your CHIANTI version is different from the one in Solarsoft. Please update your system.'
  ENDIF

  FOR i=0,nions-1 DO BEGIN
    ch_write_pop_lookup_table,ions[i],/no_execute, $
                              ldens_start=ldens_start, $
                              ldens_end=ldens_end, $
                              ldens_step=ldens_step, $
                              ltemp_start=ltmp_start, $
                              ltemp_end=ltmp_end, $
                              ltemp_step=ltmp_step, $
                              ioneq_file=ioneq_file
    IF i EQ 0 THEN BEGIN
      print,''
      print,'The log density parameters to be used for the calculation are:'
      print,'             Start      End    Step-size'
      print,format='(8x,3f10.2)',ldens_start,ldens_end,ldens_step
      print,''
      print,'The log temperature parameters to be used for each ion are:'
    ENDIF
    print,format='(3x,a5,3f10.2)',strpad(ions[i],5,fill=' ',/after), $
          ltmp_start,ltmp_end,ltmp_step
    IF n_elements(ltemp_start) EQ 0 THEN junk=temporary(ltmp_start)
    IF n_elements(ltemp_end) EQ 0 THEN junk=temporary(ltmp_start)
    IF n_elements(ltemp_step) EQ 0 THEN junk=temporary(ltmp_start)
  ENDFOR 
  print,''
  print,'Use the keywords ldens_start, ldens_end, ldens_step to change these values.'
  print,'                 ltemp_start, ltemp_end, ltemp_step'
  print,''
  print,'**LOOKUP TABLES NOT CALCULATED.**'
  print,'**Please run this routine again with the keyword /EXECUTE.**'
  return
ENDIF



;
; This gets written to the working directory
;
avalfile='pop_lookup_line_list.txt'
;
chck=file_search(avalfile,count=count)
IF count NE 0 THEN file_delete,avalfile

;
; Go through the list of ions (IONS) and write out the population
; lookup table. Also write the population line list file, adding
; entries for each transition in DATA.
;
FOR i=0,nions-1 DO BEGIN
  j=where(data.ion EQ ions[i],nj)
  junk=temporary(upper)
  FOR k=0,nj-1 DO BEGIN
    jk=j[k]
    levels=ch_dem_process_index_string(data[jk].index)
    FOR il=0,n_elements(levels)-1 DO BEGIN
      output=ch_dem_read_avals(avalfile,ion_name=data[jk].ion, $
                               lvl1=levels[il].lower,lvl2=levels[il].upper, $
                               status=status)
      IF status EQ 0 THEN BEGIN 
        ch_dem_write_avals,data[jk].ion,levels[il].lower,levels[il].upper, outfile=avalfile, $
                           line_data=data[jk],status=status
        IF status EQ 1 THEN BEGIN
          print,'% CH_DEM_WRITE_AVALS: inconsistency between transition indices and wavelength. Returning...'
          print,format='(20x,a5,2i5,f12.3)',data[jk].ion,levels[il].lower, $
                levels[il].upper, data[jk].wvl
          print,'  This is due to an error in the line list file. The output from the CHIANTI which_line'
          print,'  routine is shown below for the requested wavelength.'
          print,''
          which_line,data[jk].ion,data[jk].wvl,/narr
          return 
        ENDIF
        IF status EQ 2 THEN BEGIN
          print,'% CH_DEM_WRITE_AVALS: requested transition was not found in CHIANTI wgfa file. Returning...'
          print,format='(20x,a5,2i5,f12.3)',data[jk].ion,levels[il].lower, $
                levels[il].upper, data[jk].wvl
          return
        ENDIF 
      ENDIF 
    ENDFOR 
    IF n_elements(upper) EQ 0 THEN upper=levels.upper ELSE upper=[upper,levels.upper]
  ENDFOR
  upper=upper[sort(upper)]
  upper=upper[uniq(upper)]
  n=n_elements(upper)
  levstr=trim(upper[0])
  IF n GT 1 THEN BEGIN
    FOR j=1,n-1 DO levstr=levstr+', '+trim(upper[j])
  ENDIF 
  print,'Writing lookup tables...'
  print,'   Ion: '+trim(ions[i])+'  levels: '+levstr
 ;
  outfile='pop_lookup_'+trim(ions[i])+'.txt'
 ;
  IF n_elements(dir_lookup) NE 0 THEN outfile=concat_dir(dir_lookup,outfile)
 ;
 ; The following determines whether the lookup table needs to be
 ; written: swtch=1 means yes, swtch=0 means no.
 ;
 ; If the table already exists, then I check to see if it has all the
 ; levels contained in upper. If not, then the table needs to be
 ; re-written. I also check if the CHIANTI version of the lookup table
 ; matches the latest version on the user's computer.
 ;
 ; In the case the file is re-written, I want to keep any levels that
 ; were in the original file. For this reason I use 'extra_levels'.
 ;
  swtch=0
  extra_levels=-1
  chck=file_search(outfile,count=count)
  IF count EQ 1 THEN BEGIN
    pop_check=ch_read_pop_lookup_table(outfile)
    FOR j=0,n-1 DO BEGIN
      k=where(upper[j] EQ pop_check.levels,nk)
      IF nk EQ 0 THEN BEGIN
        swtch=1
        extra_levels=pop_check.levels
      ENDIF 
    ENDFOR
    IF pop_check.chianti_version NE ch_get_version() THEN swtch=1
  ENDIF ELSE BEGIN
    swtch=1
  ENDELSE
 ;
  IF keyword_set(overwrite) THEN swtch=1
 ;
  IF swtch EQ 1 THEN BEGIN
   ;
   ; This if statement adds extra_levels (if it exists) to upper,
   ; sorts it, and removes duplicate entries.
   ;
    IF extra_levels[0] NE -1 THEN BEGIN
      levels=[upper,extra_levels]
      levels=levels[uniq(levels,sort(levels))]
    ENDIF ELSE BEGIN
      levels=upper
    ENDELSE
   ;
   ; ch_write_pop_lookup_table has its own check for whether to
   ; overwrite a lookup table, but I over-ride this here by setting
   ; /overwrite.
   ;
    ch_write_pop_lookup_table,ions[i],levels=levels,outfile=outfile, $
                            ldens_start=ldens_start, $
                            ldens_end=ldens_end, $
                            ldens_step=ldens_step, $
                            ltemp_start=ltmp_start, $
                            ltemp_end=ltmp_end, $
                            ltemp_step=ltmp_step,/overwrite
    IF n_elements(ltemp_start) EQ 0 THEN junk=temporary(ltmp_start)
    IF n_elements(ltemp_end) EQ 0 THEN junk=temporary(ltmp_start)
    IF n_elements(ltemp_step) EQ 0 THEN junk=temporary(ltmp_start)
    print,'% CH_DEM_WRITE_LOOKUP_TABLES: the file '+outfile+' has been written.'
  ENDIF ELSE BEGIN
    print,'% CH_DEM_WRITE_LOOKUP_TABLES: current lookup file is OK. Not over-written.'
  ENDELSE 
ENDFOR 

END
