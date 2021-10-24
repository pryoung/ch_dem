
PRO ch_dem_get_wvl_aval, ion_name, lvl1, lvl2, wvl, aval, wvl_check=wvl_check,  $
                         dwvl_check=dwvl_check, quiet=quiet, status=status

;+
; NAME:
;     CH_DEM_GET_WVL_AVAL
;
; PURPOSE:
;     Extracts the wavelength and A-value from CHIANTI for the
;     specified ion and atomic transition.
;
; CATEGORY:
;     CHIANTI; differential emission measure (DEM); data access.
;
; CALLING SEQUENCE:
;     CH_DEM_GET_WVL_AVAL, Ion_Name, Lvl1, Lvl2, Wvl, Aval
;
; INPUTS:
;     Ion_Name: Name of ion in CHIANTI format (e.g., 'fe_13' for Fe
;               XIII).
;     Lvl1:     CHIANTI index of transition's lower level.
;     Lvl2:     CHIANTI index of transition's upper level.
;
; OPTIONAL INPUTS:
;     Wvl_Check: A wavelength (angstroms). The routine checks if the
;                wavelength of transition lvl1-lvl2 is within DWVL_CHECK
;                angstroms of WVL_CHECK. If not, then status=1.
;     Dwvl_Check:When checking the wavelength, this specifies the
;                range to check, i.e., +/- DWVL_CHECK. The default is 0.5
;                angstroms.
;
; KEYWORD PARAMETERS:
;     QUIET:     If set, then warning messages are not printed. 
;
; OUTPUTS:
;     Wvl:   Wavelength (angstroms) for requested transition.
;     Aval:  A-value (s^-1) for requested transition.
;
; OPTIONAL OUTPUTS:
;     Status:  An integer with one of the values:
;               0 - no problems
;               1 - the output wavelength WVL is not consistent with
;                   WVL_CHECK.
;               2 - the transition lvl1-lvl2 was not found in the wgfa
;                   file. 
; 
; CALLS:
;     READ_WGFA_STR, CONVERTNAME, ZION2FILENAME
;
; EXAMPLE:
;     IDL> ch_dem_get_wvl_aval, 'fe_13', 1, 20, wvl, aval
;     IDL> ch_dem_get_wvl_aval, 'fe_13', 1, 20, wvl, aval, wvl_check=202.04
;
; MODIFICATION HISTORY:
;     Ver.1, 22-Jul-2019, Peter Young
;-

IF n_params() LT 3 THEN BEGIN
  print,'Use:  IDL> ch_dem_get_wvl_aval, ion_name, lvl1, lvl2, wvl, aval [, wvl_check=,'
  print,'                                dwvl_check=, /quiet, status= ]'
  return
ENDIF

IF n_elements(dwvl_check) EQ 0 THEN dwvl_check=0.5

;
; Read wgfa file
;
convertname,ion_name,iz,iion
zion2filename,iz,iion,fname
read_wgfa_str,fname+'.wgfa',wgfastr,ref

iw=where(wgfastr.lvl1 EQ lvl1 AND wgfastr.lvl2 EQ lvl2,niw)


;
; Extract outputs and perform wavelength check.
;
status=0
IF niw NE 0 THEN BEGIN
  aval=wgfastr[iw[0]].aval
  wvl=wgfastr[iw[0]].wvl
 ;
  IF n_elements(wvl_check) NE 0 THEN BEGIN
    IF abs(wvl-wvl_check) GT dwvl_check THEN status=1
  ENDIF 
ENDIF ELSE BEGIN
  status=2
ENDELSE 

;
; Print warning messages (if necessary).
;
IF NOT keyword_set(quiet) THEN BEGIN 
  CASE status OF
    2: print,'% CH_DEM_GET_WVL_AVAL: specified transition was not found in wgfa file.'
    1: print,'% CH_DEM_GET_WVL_AVAL: transition wavelength does not match WVL_CHECK.'
    ELSE: 
  ENDCASE 
ENDIF 

END
