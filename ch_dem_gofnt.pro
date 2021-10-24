
FUNCTION ch_dem_gofnt, line_id_file, wvl, log_temp=log_temp, $
                       log_dens=log_dens, log_press=log_press, $
                       abund_file=abund_file, noabund=noabund, $
                       photons=photons, dir_lookup=dir_lookup, $
                       file_lookup=file_lookup

;+
; ****THIS ROUTINE HAS BEEN REPLACED BY CH_LOOKUP_GOFNT. DO NOT
; USE!!****
;
; PURPOSE:
;      This routine takes a population lookup table and computes the
;      contribution function for the emission line. The definition of
;      the contribution is the same as that of the standard CHIANTI
;      routine gofnt.pro.
;
; INPUTS:
;      Line_Id_File: The name of a file containing line IDs. See the
;               routine ch_dem_read_line_ids for the format. It can
;               also be the structure returned by
;               ch_dem_read_line_ids. 
;      Wvl:     The wavelength of the line for which the contribution
;               function is required.
;      Lookup_Table: The name of a file containing the population
;               lookup table for the ion.
;
; OPTIONAL INPUTS:
;      Log_Temp: A 1D array of Log10 temperatures. If not specified,
;                then the temperature array in the lookup table will be used.
;      Log_Dens: A scalar specifying the Log10 of the electron number
;                density. Either log_dens or log_press should be
;                defined. 
;      Log_Press: A scalar specifying the Log10 of the electron
;                 pressure (N_e*T). Either log_dens or log_press
;                 should be defined.
;      Abund_File: The name of  CHIANTI format element abundance
;                 file. If not specified then the default,
;                 !abund_file, is used.
;      File_Lookup: The name of the population lookup table file
;                 (produced by ch_write_pop_lookup_table.pro).
;      Dir_Lookup: The name of a directory containing population
;                 lookup table files. If set, then the routine will
;                 expect to find a filename of the form
;                 'pop_lookup_[ion_name].txt' in this directory.
;
; KEYWORD PARAMETERS:
;      NOABUND: If set, then the contribution function is not
;               multiplied by the element abundance.
;      PHOTONS: If set, then the contribution function is in photons

;
; OUTPUTS:
;      An IDL structure with the tags:
;      .ltemp  Log temperatures at which contrib. fn. defined.
;      .gofnt  The contribution function in erg cm^3 s^-1 sr^-1.
;      .ldens  Log densities at which contrib. fn. defined. Either a
;              scalar or a 1D array (if log_press defined).
;      .lpress Log pressure at which contrib. fn. defined. Set to -1
;              if log_press not specified.
;      .abund_file  Name of the element abundance file.
;      .ioneq_file  Name of the ionization balance file.
;      .chianti_version  The version of CHIANTI that was used.
;      .line_id_data  The structure returned by ch_dem_read_line_ids
;                     for the user-specified line.
;      .time_stamp  The time at which the structure was created.
;
; MODIFICATION HISTORY:
;      Ver.0.1, 5-Feb-2019, Peter Young
;      Ver.0.2, 8-Feb-2019, Peter Young
;        Allow LINE_ID_FILE to be the structure returned by
;        ch_dem_read_line_ids.pro.
;      Ver.0.3, 9-Apr-2019, Peter Young
;        Removed call to read_wgfa.
;       
;-


print,''
print,'**** THIS ROUTINE IS OBSOLETE. DO NOT USE!! ****'
print,''

;
; LINE_ID_FILE can be either the name of a line identification file,
; or it can be the structure created from the file by
; ch_dem_read_line_ids.pro. 
;
IF datatype(line_id_file) EQ 'STC' THEN BEGIN
  line_data=line_id_file
ENDIF ELSE BEGIN 
  line_data=ch_dem_read_line_ids(line_id_file)
ENDELSE 
;
getmin=min(abs(wvl-line_data.wvl),imin)
IF n_elements(wvl_chck) EQ 0 THEN wvl_chck=0.5
IF getmin LE wvl_chck THEN BEGIN
  line_data=line_data[imin]
  IF NOT keyword_set(quiet) THEN BEGIN
    print,'% CH_DEM_GOFNT: found the following wavelength match...'
    ion=strpad(line_data.ion,5,fill=' ',/after)
    trans=strpad(line_data.index,40,fill=' ',/after)
    print,format='(10x,a5,f10.2,3x,a40)',ion,wvl,trans
  ENDIF 
ENDIF ELSE BEGIN
  print,'% CH_DEM_GOFNT: no transitions found within '+trim(string(format='(f10.2)',wvl_chck))+' angstroms of WVL. Returning...'
  return,-1
ENDELSE
;
ion_name=line_data.ion

;
; Do some checks on the population lookup table file.
;
n1=n_elements(dir_lookup)
n2=n_elements(file_lookup)
CASE 1 OF
  n1 EQ 0 AND n2 EQ 0: BEGIN
    print,'% CH_DEM_GOFNT: please specify EITHER dir_lookup OR file_lookup. Returning...'
    return,-1
  END
  n1 NE 0 AND n2 EQ 0: lookup_file=concat_dir(dir_lookup,'pop_lookup_'+ion_name+'.txt')
  n1 EQ 0 AND n2 NE 0: lookup_file=file_lookup
  n1 NE 0 AND n2 NE 0: lookup_file=file_lookup
END
;
chck=file_info(lookup_file)
IF chck.exists EQ 0 THEN BEGIN
  print,'% CH_DEM_GOFNT: the population lookup table was not found. Returning...'
  return,-1
ENDIF 

;
; Process the index string into a structure.
;
indstr=ch_dem_process_index_string(line_data.index)
ntrans=n_elements(indstr)


;
; Look for the A-value file associated with the lookup tables.
;
avalfile='pop_lookup_line_list.txt'
IF n_elements(dir_lookup) NE 0 THEN avalfile=concat_dir(dir_lookup,avalfile)
chck=file_search(avalfile,count=count)




;
; If the A-val file exists, then we need to create aval_str, which
; gets sent to ch_lookup_gofnt.
;
IF count NE 0 THEN BEGIN 
  str={lvl1: 0, lvl2: 0, wvl: 0., aval: 0.}
  aval_str=replicate(str,ntrans)
  swtch=0
  FOR j=0,ntrans-1 DO BEGIN
    IF swtch EQ 0 THEN BEGIN 
      l1=indstr[j].lower
      l2=indstr[j].upper
 ;
 ; The following finds the wavelength and A-value for the
 ; transition. They should be in the "pop_lookup_line_list" file
 ; (which is the quickest way to get them), but otherwise they are
 ; extracted from the CHIANTI wgfa file. 
 ;
      output=ch_dem_read_avals(avalfile,ion_name=line_data.ion,lvl1=l1,lvl2=l2, $
                               status=status,wvl=wvl,aval=aval)
      IF status NE 0 THEN BEGIN
        aval_str[j].lvl1=l1
        aval_str[j].lvl2=l2
        aval_str[j].wvl=wvl
        aval_str[j].aval=aval
      ENDIF ELSE BEGIN
        swtch=1
        junk=temporary(aval_str)
      ENDELSE
    ENDIF
  ENDFOR
ENDIF 

;
; The call to ch_lookup_gofnt differs depending if aval_str or
; lower_levels/upper_levels are defined.
;
IF n_tags(aval_str) NE 0 THEN BEGIN
  print,'% CH_DEM_GOFNT: using pre-loaded wavelengths and A-values.'
  g=ch_lookup_gofnt(line_data.ion,aval_str=aval_str, $
                    log_temp=log_temp,log_dens=log_dens,log_press=log_press, $
                    abund_file=abund_file, noabund=noabund, $
                    photons=photons,file_lookup=lookup_file)
ENDIF ELSE BEGIN 
  lower_levels=indstr.lower
  upper_levels=indstr.upper
  print,'% CH_DEM_GOFNT: reading wavelenghts and A-values.'
  g=ch_lookup_gofnt(line_data.ion,lower_levels=lower_levels,upper_levels=upper_levels, $
                    log_temp=log_temp,log_dens=log_dens,log_press=log_press, $
                    abund_file=abund_file, noabund=noabund, $
                    photons=photons,file_lookup=lookup_file)
ENDELSE 


output=add_tag(g,line_data,'line_id_data')

;; output={ ltemp: log_temp, $
;;          gofnt: contrib_data, $
;;          ldens: log_dens, $
;;          lpress: lpress, $
;;          abund_file: abund_file, $
;;          ioneq_file: !ioneq_file, $
;;          chianti_version: popstr.chianti_version, $
;;          line_id_data: line_data, $
;;          time_stamp: systime() }


return, output

END
