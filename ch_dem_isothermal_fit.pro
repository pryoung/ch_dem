
FUNCTION ch_dem_isothermal_fit, line_data, ldens=ldens, lpress=lpress, $
                            abund_file=abund_file, $
                            interr_scale=interr_scale, dir_lookup=dir_lookup, $
                            quiet=quiet, $
                            ab_elt_fix=ab_elt_fix, fixed_abund=fixed_abund, $
                            ioneq_file=ioneq_file, swtch_ab=swtch_ab

;+
; NAME:
;      CH_DEM_ISOTHERMAL_FIT
;
; PURPOSE:
;      Fits an isothermal model to the line intensities. That is, a
;      single value of temperature and emission measure. The intensity of
;      a line is given by I=Ab*G*EM, where Ab is the element abundance, G
;      is the contribution function and EM is the emission measure.
;
; CATEGORY:
;      CHIANTI; differential emission measure; DEM.
;
; CALLING SEQUENCE:
;      Result = CH_DEM_ISOTHERMAL_FIT(Line_Data)
;
; INPUTS:
;      Line_Data:   A structure containing details of the lines used
;                   to derive the DEM. 
;
; OPTIONAL INPUTS:
;      Ldens:  Specifies the logarithm of the electron number density
;              (units: cm^-3) to be used for the calculation.Either
;              ldens or lpress should be specified (not both).
;      Lpress:  Specifies the logarithm of the electron pressure
;              (units: K cm^-3) to be used for the calculation. Either
;              ldens or lpress should be specified (not both).
;      Abund_File:  The CHIANTI abundance file to be used for the
;                   calculation. The default is !abund_file.
;      Interr_Scale: A number between 0 and 1 that adds an additional
;                    amount to the lines' intensity errors that
;                    is interr_scale*intensity. For example,
;                    interr_scale=0.1 will add 10% of the intensity to
;                    the error (added in quadrature). This can be useful if
;                    the intensity errors are very small. 
;      Dir_Lookup:  Directory containing the pre-calculated lookup tables
;      Quiet:   Suppresses some text produced by program.
;      Ab_Elt_Fix:  Either an integer or a string specifying the
;                   reference element. For example, ab_elt_fix=26 or
;                   'fe' specifies iron. If not specified, then the
;                   routine uses the element with the most number of
;                   lines in LINE_DATA.
;      Ioneq_File:  The name of an ionization equilibrium file. The
;                   default is to use the CHIANTI file (!ioneq_file).
;      Swtch_Ab:  An array of same size as the number of unique
;                 elements in LINE_DATA. A zero indicates the
;                 element's abundance should be fixed, and a
;                 one indicates the abundance should be a
;                 variable. This is expected to be used only in
;                 special cases - see CHIANTI Technical Report
;                 No. 22. 
;
; KEYWORD PARAMETERS:
;      QUIET: If set, text information to screen is not printed.
;      FIXED_ABUND:  If set, then the abundances are not varied.
; 
; OUTPUTS:
;      A structure containing the results. The tags are:
;        .method  Set to 'isothermal'.
;        .aa   3-element array containing fit params (EM0,T0,sigmaT).
;        .sigmaa Errors on fit params.
;        .chi2 Reduced chi^2 for fit.
;        .yfit The line intensities from the model.
;        .log_t  Log of the plasma temperature (T).
;        .log_t_err  Error on logT.
;        .log_em  Log of the emission measure.
;        .log_em_err  Error on logEM.
;        .line_data  Contains the line_data input.
;        .abstr  Structure containing information about the
;                abundances.
;        .interr_scale Value of interr_scale input.
;        .log_dens Value of ldens input (if specified).
;        .log_press Value of lpress input (if specified).
;
; CALLS:
;      CH_DEM_ADD_CONTRIB, CH_DEM_PROCESS_ABUND, CH_ISOTHERMAL_ABUND_RESULTS
;
; MODIFICATION HISTORY:
;      Ver. 1, 11-May-2025, Peter Young
;      Ver. 2, 30-Jun-2025, Peter Young
;        Fixed bug in specification of INIT when abundances are varied.
;-


IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> output=ch_dem_isothermal_fit(line_data [, ldens=, lpress=, abund_file='
  print,'                                   interr_scale=, dir_lookup=, ioneq_file=  '
  print,'                                   /fixed_abund, /quiet, ab_elt_fix=, '
  print,'                                   swtch_ab= ])'
  return,-1
ENDIF 

; Check the ldens and lpress inputs. We need one or the other, but not
; both (or neither).
;
IF n_elements(ldens) NE 0 AND n_elements(lpress) NE 0 THEN BEGIN
  message,/info,/cont,'Please specify either LDENS or LPRESS, but not both. Returning...'
  return,-1
ENDIF 
;
IF n_elements(ldens) EQ 0 AND n_elements(lpress) EQ 0 THEN BEGIN
  message,/info,/cont,'Please specify either LDENS or LPRESS. Returning...'
  return,-1
ENDIF

nl=n_elements(line_data)
k=where(line_data.int NE -1.,nk)
IF nk EQ 0 THEN BEGIN
  message,/info,/cont,'The input LINE_DATA structure does not contain observed intensities. Please run ch_dem_read_line_ids to load intensities from an input file. Returning...'
  return,-1
ENDIF 

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

ltemp=findgen(81)/20.+4.

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
  message,/info,/cont,'error found. Returning...'
  return,-1
ENDIF 

;
; This is an initial definition of INIT that is modified later.
; SCL is not used in this routine, so just set to 1.
;
init=[1d0,1d0]
scl=init

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
  message,/info,/cont,'Could not find the element abundance file. Returning...'
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
; INITIAL PARAMETER SETUP
; -----------------------
; I make estimates of the temperature and EM using the line list and the
; contribution functions of the reference element.
;
k=where(abstr[ld_fit.ab_ind].type EQ 1,nk)
contrib_max=dblarr(nk)
FOR i=0,nk-1 DO BEGIN
  getmax=max(ld_fit[i].contrib)
  contrib_max[i]=getmax
ENDFOR
em_init=ld_fit[k].int/abstr[ld_fit[k].ab_ind].abund/contrib_max
;
init[0]=mean(ld_fit[k].logt_max)
init[1]=alog10(mean(em_init))

;
; Modify the line intensity errors with INTERR_SCALE (if ncessary). 
;
int=ld_fit.int
err=ld_fit.err
;
; Add an additional component to the error if interr_scale is given.
;
IF n_elements(interr_scale) NE 0 THEN BEGIN
  IF interr_scale LT 0 OR interr_scale GT 1 THEN BEGIN
    message,/info,/cont,'The input interr_scale should take a value between 0 and 1. Returning...'
    return,-1
  ENDIF ELSE BEGIN
    err=sqrt(err^2 + (int*interr_scale)^2 )
  ENDELSE
ENDIF


;
; SET PARAMETER LIMITS (PARINFO structure)
; --------------------
; Note that I'm fitting log(EM) and log(T). I'm not setting any parameter limits.
;
np=n_elements(init)
parinfo=replicate({fixed: 0, limited: [0,0], limits:[0.d,0.d], step: 0.},np)
;parinfo[*].limited=[1,0]
;parinfo[0:n_nodes-1].limits[0]=1e-20



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
        ltemp: ltemp }


;
; Check the number of free parameters to see if the fit is possible.
;
IF np GE nfit THEN BEGIN
   print,'% CH_DEM_LINEAR_FIT: there are too many free parameters! Try using /fixed_abund to'
   print,'                     fix the element abundances.'
   print,'  No. free parameters: '+trim(np)
   print,'  No. lines:           '+trim(nfit)
   print,''
   return,-1
ENDIF 

;
; This is the minimization procedure, producing fit parameters AA.
;
aa=mpfitfun('ch_dem_isothermal_fit_fn',findgen(nfit), int , err, init, $
            perror=sigmaa, /quiet, bestnorm=bestnorm,yfit=yfit, $
            parinfo=parinfo,status=status, functargs= other,errmsg=errmsg)


;
; This updates abstr with the optimized abundances.
;
ch_isothermal_abund_results,abstr,ld,aa,sigmaa,ltemp,interr_scale=interr_scale



chi2=bestnorm/float(nfit-np)

;
; Write out the following to the IDL window:
;  - the comparison between the observed intensity and model intensity
;    for each line
;  - the abundance results
;
ld_fit.model_int=yfit
IF NOT keyword_set(quiet) THEN ch_dem_write_results,ld_fit, abstr


chi2=bestnorm/float(nfit-np)
IF NOT keyword_set(quiet) THEN BEGIN 
  print,''
  print,'Fit parameters: '
  print,format='("  log T",f9.3," +/-",f7.3,"   ( T =",e9.2," )")',aa[0],sigmaa[0],10.^aa[0]
  print,format='(" log EM",f9.3," +/-",f7.3,"   ( EM =",e9.2," )")',aa[1],sigmaa[1],10.^aa[1]
  print,''
  print,format='("Reduced chi^2: ",f8.2)',chi2
ENDIF 


;
; Compute the T_eff of each line in LD_ALL and populate the ab_ind
; and model_int tags of LD_ALL
;
lteff=fltarr(nl)
nab=n_elements(abstr)
FOR i=0,nl-1 DO BEGIN
  ld_all[i].logt_eff=aa[0]
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
  contrib=ld_all[i].contrib
  k=where(contrib NE 0.)
  y2=spl_init(ltemp[k],alog10(contrib[k]))
  yi=spl_interp(ltemp[k],alog10(contrib[k]),y2,aa[0])
  ld_all[i].model_int=ab_i*yi*10.^aa[1]
ENDFOR 


;
; Compute column depth distribution. 
;
;; IF n_elements(lpress) NE 0 THEN ldens=lpress - ltemp
;; h_to_e_ratio=proton_dens(ltemp,abund_file=abund_file,ioneq_file=ioneq_file,/hydrogen)
;; coldepth=10.^ltemp*alog(10.)*dem*dlogt/h_to_e_ratio/(10.^ldens)^2

;
; Create the output structure.
;
IF n_elements(ldens) EQ 0 THEN ldens=-1.
IF n_elements(lpress) EQ 0 THEN lpress=-1.
IF n_elements(interr_scale) EQ 0 THEN interr_scale=-1.
output={method: 'isothermal', $
        aa: aa, $
        sigmaa: sigmaa, $
        chi2: chi2, $
        yfit: yfit, $
        ltemp: ltemp, $
        log_t: aa[0], $
        log_t_err: sigmaa[0], $ 
        log_em: aa[1], $
        log_em_err: sigmaa[1], $
;        coldepth: coldepth, $
        line_data: ld_all, $
        abstr: abstr, $
        interr_scale: interr_scale, $
        log_dens: ldens, $
        log_press: lpress}

return,output

END
