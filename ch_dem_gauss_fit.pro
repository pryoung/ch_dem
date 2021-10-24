
FUNCTION ch_dem_gauss_fit, line_data, ldens=ldens, lpress=lpress, abund_file=abund_file, $
                           interr_scale=interr_scale, dir_lookup=dir_lookup, $
                           temp_log=temp_log, initial_params=initial_params,quiet=quiet, $
                           ab_elt_fix=ab_elt_fix, fixed_abund=fixed_abund, $
                           logt_min=logt_min, logt_max=logt_max, dlogt=dlogt, $
                           swtch_ab=swtch_ab

;+
; NAME:
;      CH_DEM_GAUSS_FIT
;
; PURPOSE:
;      Fits a Gaussian DEM to a set of observed line intensities using
;      the method of Warren & Brooks (2009, ApJ).
;
; CATEGORY:
;      CHIANTI; differential emission measure; DEM.
;
; CALLING SEQUENCE:
;      Result = CH_DEM_GAUSS_FIT(Line_Data)
;
; INPUTS:
;      Line_Data:   A structure containing details of the lines used
;                   to derive the DEM. Note that the model_int and
;                   logt_eff tags will be modified by
;                   ch_dem_gauss_int. 
;
; OPTIONAL INPUTS:
;      Ldens:  Specifies the logarithm of the electron number density
;              (units: cm^-3) to be used for the calculation. Default
;              is 9.0.
;      Abund_File:  The CHIANTI abundance file to be used for the
;                   calculation. The default is
;                   sun_coronal_1992_feldman_ext.abund.
;      Interr_Scale: A number between 0 and 1 that adds an additional
;                    amount to the lines' intensity errors that
;                    is interr_scale*intensity. For example,
;                    interr_scale=0.1 will add 10% of the intensity to
;                    the error (added in quadrature). This can be useful if
;                    the intensity errors are very small. 
;      Dir_Lookup:  Directory containing the pre-calculated lookup tables
;      Initial_Params:  A 3-element array specifying the initial parameters
;               for the fit. Only use this if the function fails to
;               return a solution. 
;      Quiet:   Suppresses some text produced by program.
;      Dlogt:   Specifies the logT intervals for the temperature
;               range. Default is 0.025
;      Logt_Min:  Specifies the minimum logT value for the temperature
;               range. Default is 5.0
;      Logt_Max:  Specifies the maximum logT value for the temperature
;               range. Default is 7.0
;      Ab_Elt_Fix:  Either an integer or a string specifying the
;                   reference element. For example, ab_elt_fix=26 or
;                   'fe' specifies iron. If not specified, then the
;                   routine uses the element with the most number of
;                   lines in LINE_DATA.
;      Swtch_Ab:  An array of same size as the number of unique
;                 elements in LINE_DATA. A zero indicates the
;                 element's abundance should be fixed, and a
;                 one indicates the abundance should be a
;                 variable. This is expected to be used only in
;                 special cases - see CHIANTI Technical Report
;                 No. 22. 
;
; KEYWORD PARAMETERS:
;      TEMP_LOG: If set, then the model Gaussian is defined in logT space
;                rather than T space.
;      QUIET: If set, text information to screen is not printed.
; 
; OUTPUTS:
;      A structure containing the results. The tags are:
;        .aa   3-element array containing fit params (EM0,T0,sigmaT).
;        .sigmaa Errors on fit params.
;        .chi2 Reduced chi^2 for fit.
;        .yfit The line intensities from the model.
;        .ltemp The log temperatures at which the DEM is defined.
;        .dem   The DEM function.
;        .lteff The effective temperatures for the lines.
;        .int_ratio  Ratios of observed to model intensities.
;        .int_ratio_err Errors on the ratios.
;        .temp_log    This is keyword_set(temp_log).
;        .coldepth  The column depth (cm) corresponding to the DEM. 
;        .abund  Structure containing information about the
;                abundances. 
;
; CALLS:
;      READ_ABUND, READ_IONEQ, 
;      CH_DEM_GOFNT, CH_DEM_PROCESS_INDEX_STRING, Z2ELEMENT,
;      CH_DEM_PROCESS_ABUND, CH_DEM_INIT_ABUND
;
; MODIFICATION HISTORY:
;      Ver.1, 17-Jul-2018, Peter Young.
;      Ver.2, 19-Jul-2018, Peter Young
;          Added INTERR_SCALE= keyword; output fit parameters are now
;          de-scaled back to their correct physical values; set a
;          lower limit for the Gaussian width to prevent
;          divide-by-zero problems.
;      Ver.3, 20-Jul-2018, Terry Kucera & Peter Young
;          Added POPDir Keyword; added /LOG keyword.
;      Ver.4, 23-Jul-2018, Peter Young.
;          Fixed bug if /log wasn't defined; modified some of
;          the text output; put T_eff and model intensity into the
;          line_data input.
;      Ver.5, 24-Jul-2018, Peter Young
;          Fixed bug in calculation of T_eff values.
;      Ver.6, 25-Jul-2018, Peter Young
;          Now use temporary on SCL to clear out the common block;
;          added INITIAL_PARAMS optional input; fixed another bug in
;          calculation of T_eff; modified initial parameters for /LOG
;          case.
;      Ver.7, 30-Jul-2018, Terry Kucera
;          Added keyword QUIET.
;      Ver.8, 13-Feb-2019, Peter Young
;          Major overhaul. If the lines come from multiple elements,
;          then the relative abundances are allowed to vary; now calls
;          a single function (ch_dem_gauss_fit_fn) for both the log=0
;          and log=1 cases; got rid of common block; replaced
;          "popdir=" with "dir_lookup=".
;      Ver.9, 17-Jul-2019, Peter Young
;          Now calls ch_lookup_gofnt to compute contribution
;          functions; added input parameters for specifying the
;          temperature range (dlogt, logt_min, logt_max).
;      Ver.10, 23-Jul-2019, Peter Young
;          Now checks for $CH_LOOKUP_DEM; changed default location of
;          the aval file (pop_lookup_line_list.txt).
;      Ver.11, 29-Jul-2019, Peter Young
;          Changed keyword /log to /temp_log due to conflict.
;      Ver.12, 19-May-2020, Peter Young
;          Removed rogue help statement.
;      Ver.13, 16-Jun-2021, Peter Young
;         /FIXED_ABUND is now implemented; the fit values are now
;         loaded into ld_fit; now checks to see if there are enough
;         free parameters.
;      Ver.14, 17-Jun-2021, Peter Young
;         Added SWTCH_AB optional input; now exits gracefully if
;         ch_dem_process_abund fails. 
;-


IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> output=ch_dem_gauss_fit(line_data [, ldens=, /temp_log, initial_params=, abund_file='
  print,'                                   interr_scale=, dir_lookup=, dlogt=, logt_min= '
  print,'                                   logt_max=, /fixed_abund, initial_params=, /quiet '
  print,'                                   swtch_ab= ])'
  return,-1
ENDIF 



IF n_elements(initial_params) GT 0 AND n_elements(initial_params) NE 3 THEN BEGIN
  print,'% CH_DEM_GAUSS_FIT: the input INITIAL_PARAMS must be a 3-element array. Returning...'
  return,-1
ENDIF 


; Check the ldens and lpress inputs. We need one or the other, but not
; both (or neither).
;
IF n_elements(ldens) NE 0 AND n_elements(lpress) NE 0 THEN BEGIN
  print,'% CH_DEM_GAUSS_FIT: Please specify either LDENS or LPRESS, but not both. Returning...'
  return,-1
ENDIF 
;
IF n_elements(ldens) EQ 0 AND n_elements(lpress) EQ 0 THEN BEGIN
  print,'% CH_DEM_GAUSS_FIT: Please specify either LDENS or LPRESS. Returning...'
  return,-1
ENDIF





;
; Set the temperature range.
; Note that Warren & Brooks (2009) used 0.025 dex intervals in logT.
;
IF n_elements(logt_min) EQ 0 THEN logt_min=5.0
IF n_elements(logt_max) EQ 0 THEN logt_max=7.0
IF n_elements(dlogt) EQ 0 THEN dlogt=0.025
nt=round((logt_max-logt_min)/dlogt)+1
ltemp=findgen(nt)*dlogt+logt_min


;
; Make a copy of LINE_DATA. From this point on, line_data is not
; used. 
;
ld_all=line_data

;
; Make sure ld_all has observed intensities.
;
nl=n_elements(ld_all)
k=where(ld_all.int NE -1.,nk)
IF nk EQ 0 THEN BEGIN
  print,'% CH_DEM_GAUSS_FIT: The input LINE_DATA structure does not contain observed intensities. Please run '
  print,'                     ch_dem_read_line_ids to load intensities from an input file. Returning...'
  return,-1
ENDIF 


;
; Note on the LINE_DATA structures.
;   LINE_DATA  Input to this routine; not modified.
;   LD_ALL     Identical to LINE_DATA, but is modified by this routine
;              and then placed in the output structure (tag:
;              line_data).
;   LD         This is obtained from LD_ALL by using
;              ch_dem_process_blends, which sums any blends (if
;              necessary).
;   LD_FIT     The sub-set of LD that will actually be fitted. For
;              example, if there's only one line from an
;              element then this line is not used to compute
;              the DEM (it's only used to derive the abundance) 
;
; The following call modifies LD_ALL (by adding the "contrib" tag
; containing the contribution functions), and then sums any blends to
; create the LD structure.
;
ld=ch_dem_add_contrib(ld_all, ltemp, avalfile=avalfile, $
                      log_press=lpress, log_dens=ldens, $
                      ioneq_file=ioneq_file, dir_lookup=dir_lookup)
IF n_tags(ld) EQ 0 THEN BEGIN
  print,'% CH_DEM_GAUSS_FIT: error found. Returning...'
  return,-1
ENDIF 



;
; INITIAL PARAMETER SETUP
; -----------------------
; Set the scaling factor (SCL) and the initial parameters for
; temp_log=0/temp_log=1 cases.
; The scaling factor is used to make the parameters that are actually
; fit by MPFIT to be close to unity.
;
IF keyword_set(temp_log) THEN BEGIN
 ;
 ; Set initial parameters for DEM.
 ;   init[0] -> EM_0/1e27
 ;   init[1] -> log T0
 ;   init[2] -> sigma_logT
 ;
  scl=[1d27,1.0,1.0]   ; contains scaling parameters
  IF n_elements(initial_params) EQ 0 THEN init=[1.0,6.0,0.1] ELSE init=initial_params/scl
ENDIF ELSE BEGIN 
 ;
 ; Set initial parameters for DEM.
 ;   init[0] -> EM_0/1e27
 ;   init[1] -> T_0/1e6
 ;   init[2] -> sigma_T/1e5
 ;
  scl=[1d27,1d6,1d5]   ; contains scaling parameters
  IF n_elements(initial_params) EQ 0 THEN init=[1.0,1.0,3.0] ELSE init=initial_params/scl
ENDELSE

;
; ELEMENT ABUNDANCE SETUP
; -----------------------
; - Element abundance information is stored in ABSTR and AB_REF.
; - A subset of LD is extracted to LD_FIT -> these are the lines that
;   will be used in the DEM fitting procedure
;
IF n_elements(abund_file) EQ 0 THEN abund_file=!abund_file
chck=file_info(abund_file)
IF NOT chck.exists THEN BEGIN
  print,'% CH_DEM_GAUSS_FIT: could not find the element abundance file. Returning...'
  return,-1
ENDIF 
abstr=ch_dem_process_abund(ld,abund_file=abund_file,ab_elt_fix=ab_elt_fix, $
                           init=init,out_init=out_init,fixed_abund=fixed_abund, $
                           swtch_ab=swtch_ab)
IF n_tags(abstr) EQ 0 THEN return,-1
init=out_init
k=where(abstr.type EQ 1)
ab_ref=abstr[k[0]].abund
;
ab_type=abstr[ld.ab_ind].type
k=where(ab_type GE 1)
ld_fit=ld[k]
IF n_tags(ld_fit) EQ 0 THEN return,-1
nfit=n_elements(ld_fit)


;
; Modify the line intensity errors with INTERR_SCALE (if ncessary). 
;
IF n_elements(interr_scale) NE 0 THEN BEGIN
  IF interr_scale LT 0 OR interr_scale GT 1 THEN BEGIN
    print,'% CH_DEM_GAUSS_FIT:  the input interr_scale should take a value between 0 and 1. Returning...'
    return,-1
  ENDIF ELSE BEGIN
    ld_fit.err=sqrt(ld_fit.err^2 + (ld_fit.int*interr_scale)^2 )
  ENDELSE
ENDIF
int=ld_fit.int
err=ld_fit.err




;
; SET PARAMETER LIMITS (PARINFO structure)
; --------------------
;  - the lower limit for sigma is necessary to prevent divide-by-zero problems.
;
np=n_elements(init)
parinfo=replicate({fixed: 0, limited: [0,0], limits:[0.d,0.d], step: 0.},np)
parinfo[0].limited=[1,0]
IF keyword_set(temp_log) THEN BEGIN
  parinfo[1].limited=[1,1]
  parinfo[1].limits=minmax(ltemp)
  parinfo[2].limited=[1,0]
  parinfo[2].limits[0]=1e-10
ENDIF ELSE BEGIN
  parinfo[1].limited=[1,1]
  parinfo[1].limits=minmax(10.^ltemp)/scl[1]
  parinfo[2].limited=[1,0]
  parinfo[2].limits[0]=1e-5
ENDELSE
IF np GT 3 THEN BEGIN
  FOR i=0,np-4 DO parinfo[i+3].limited=[1,0]
ENDIF



;
; I pass additional data to the fit function by using the structure
; 'other'. Note that 'line_data', 'abstr', 'ab_ref', 'ltemp' and 'dem'
; are essential. 
;
; The other tags are specific to the type of fit function.  
;
junk=temporary(other)
other={ line_data: ld_fit, $
        abstr: abstr, $
        ab_ref: ab_ref, $
        ltemp: ltemp, $
        scl: scl, $
        temp_log: keyword_set(temp_log) }


;
; Check the number of free parameters to see if the fit is possible.
;
IF np GE nfit THEN BEGIN
   print,'% CH_DEM_GAUSS_FIT: there are too many free parameters! Try using /fixed_abund to'
   print,'                    fix the element abundances.'
   print,'  No. free parameters: '+trim(np)
   print,'  No. lines:           '+trim(nfit)
   print,''
   return,-1
ENDIF 


;
; This is the minimization procedure, producing fit parameters AA.
;
aa=mpfitfun('ch_dem_gauss_fit_fn',findgen(nfit), int , err, init, $
            perror=sigmaa, /quiet, bestnorm=bestnorm,yfit=yfit, $
            parinfo=parinfo,status=status, functargs= other)


k=where(sigmaa EQ 0.,nk)
IF nk GT 0 THEN BEGIN
  print,'% CH_DEM_GAUSS_FIT: at least one of the parameter errors (sigmaa) is zero, therefore fit has failed.'
  print,'                    A possible solution is to adjust the initial parameters using the INIT optional input.'
  print,'                    The current values are:'
  print,format='(a27,3e12.3)',init*scl
  return,-1
ENDIF 


;
; Need to re-scale the DEM fit parameters
;
aa[0:2]=aa[0:2]*scl
sigmaa[0:2]=sigmaa[0:2]*scl

;
; Load fitted intensities into LD_FIT
;
ld_fit.model_int=yfit

;
; Compute the DEM from the fit parameters.
;
dem=ch_dem_gauss_compute_dem(ltemp,aa[0:2],temp_log=temp_log)

;
; This updates abstr with the optimized abundances.
;
ch_dem_abund_results,abstr,ld,aa,sigmaa,ltemp,dlogt,dem,interr_scale=interr_scale


IF NOT keyword_set(quiet) THEN ch_dem_write_results,ld_fit, abstr


chi2=bestnorm/float(nfit-np)
IF NOT keyword_set(quiet) THEN BEGIN 
  print,''
  print,'Gaussian params: '
  print,format='("    EM_0",e12.2," +/-",e12.2)',aa[0],sigmaa[0]
  IF keyword_set(temp_log) THEN BEGIN
    print,format='(" log T_0",e12.2," +/-",e12.2)',aa[1],sigmaa[1]
    print,format='("sig_logT",e12.2," +/-",e12.2)',aa[2],sigmaa[2]
  ENDIF ELSE BEGIN 
    print,format='("     T_0",e12.2," +/-",e12.2)',aa[1],sigmaa[1]
    print,format='("   sig_T",e12.2," +/-",e12.2)',aa[2],sigmaa[2]
  ENDELSE 
  print,''
  print,format='("Reduced chi^2: ",f8.2)',chi2
ENDIF 

;
; Compute the T_eff of each line in LD_ALL, and populate the ab_ind
; and model_int tags of LD_ALL.
;
lteff=fltarr(nl)
nab=n_elements(abstr)
FOR i=0,nl-1 DO BEGIN
  contrib_fn=ld_all[i].contrib
  func=contrib_fn*dem*10.^ltemp
  getmax=max(func,imax)
  lteff[i]=ltemp[imax]
  ld_all[i].logt_eff=lteff[i]
 ;
 ; Populate ld_all.ab_ind
 ;
  FOR j=0,nab-1 DO BEGIN
    k=where(ld_all.element_num EQ abstr[j].elt_num)
    ld_all[k].ab_ind=j
  ENDFOR
 ;
 ; Populate ld_all.model_int
 ;
  ab_i=abstr[ld_all[i].ab_ind].abund
  ld_all[i].model_int=ab_i*total(dem*contrib_fn*10.^ltemp*dlogt*alog(10.))
ENDFOR 



;
; Compute column depth distribution. 
;
IF n_elements(lpress) NE 0 THEN ldens=lpress - ltemp
h_to_e_ratio=proton_dens(ltemp,abund_file=abund_file,ioneq_file=ioneq_file,/hydrogen)
coldepth=10.^ltemp*alog(10.)*dem*dlogt/h_to_e_ratio/(10.^ldens)^2

;
; Create the output structure.
;
IF n_elements(ldens) EQ 0 THEN ldens=-1.
IF n_elements(lpress) EQ 0 THEN lpress=-1.
IF n_elements(interr_scale) EQ 0 THEN interr_scale=-1.
output={method: 'gauss', $
        aa: aa, $
        sigmaa: sigmaa, $
        chi2: chi2, $
        yfit: yfit, $
        ltemp: ltemp, $
        dem: dem, $
        coldepth: coldepth, $
        line_data: ld_all, $
        abstr: abstr, $
        interr_scale: interr_scale, $
        log_dens: ldens, $
        log_press: lpress, $
        temp_log: keyword_set(temp_log) }



return,output

END
