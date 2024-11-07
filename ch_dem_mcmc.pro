
FUNCTION ch_dem_mcmc, line_data, ltemp=ltemp, lpress=lpress, ldens=ldens, $
                      interr_scale=interr_scale, $
                      dem=dem, dlogt=dlogt, abund_file=abund_file, $
                      fixed_abund=fixed_abund, $
                      nsim=nsim, mcmc_savefile=mcmc_savefile

;+
; NAME:
;     CH_DEM_MCMC
;
; PURPOSE:
;     Compute a DEM using MCMC DEM method (PINTofALE).
;
; CATEGORY:
;     CHIANTI; differential emission measure (DEM).
;
; CALLING SEQUENCE:
;     Result = CH_DEM_MCMC ( Line_Data )
;
; INPUTS:
;     Line_Data:   A structure containing details of the lines used
;                  to derive the DEM. Should be in the format created
;                  by ch_dem_read_line_ids.pro.
;
; OPTIONAL INPUTS:
;     Ltemp:  An array of log10 temperatures at which the DEM will be
;             defined. If not specified, then the routine computes the
;             Tmax of each ion in LINE_DATA and then sets the range to
;             be -0.15 and +0.15 of the min and max temperatures.
;     Dlogt:  The bin size of the LTEMP array. Set to 0.10 dex by
;             default. 
;     Ldens:  Specifies the logarithm of the electron number density
;             (units: cm^-3) to be used for the calculation. Default
;             is 9.0.
;     Abund_File:  The CHIANTI abundance file to be used for the
;                  calculation. The default is
;                  sun_coronal_1992_feldman_ext.abund.
;     Interr_Scale: A number between 0 and 1 that adds an additional
;                   amount to the lines' intensity errors that
;                   is interr_scale*intensity. For example,
;                   interr_scale=0.1 will add 10% of the intensity to
;                   the error (added in quadrature). This can be useful if
;                   the intensity errors are very small.
;     Nsim:   The number of iterations to be performed. The default (set
;             in the mcmc_dem routine) is 100.
;     Mcmc_Savefile:  The name of the file where the MCMC output parameters
;             are saved. The default is 'mcmc.sav' in the working directory.
;	
; KEYWORD PARAMETERS:
;     FIXED_ABUND:  If set, then element abundances are fixed.
;
; OUTPUTS:
;     An IDL structure containing the tags:
;      .method  'mcmc'
;      .ltemp   Temperatures at which DEM tabulated.
;      .dem     The DEM values.
;      .line_data  The line_data input structure modified to include model
;                  intensities and T_eff values.
;      .abstr   Element abundance structure containing output abundances.
;      .interr_scale  The value of interr_scale (if set).
;      .log_dens   Log of density (if set).
;      .log_press  Log of pressure (if set).
;      .nsim    Value of nsim.
;      .simdem  Array of size (nT,nsim+1) containing all of the simulated
;               DEMs.
;      .demerr  Array of size (nT,2) giving the lower and upper bounds on the DEM.
;      .simprb  Array of size nsim giving the likelihood array (for checking
;               convergence).
;      .chi2_proxy  Proxy for the reduced chi^2 as suggested in MCMC software.
;
;     simdem[*,nsim] is the "best" DEM, which is also stored in the .dem tag.
;
; RESTRICTIONS:
;     Requires the PINTofALE IDL software to be installed (not
;     available in Solarsoft).
;
; CALLS:
;     CH_DEM_PROCESS_ABUND, MCMC_DEM, CH_DEM_ADD_CONTRIB, PROTON_DENS,
;     READ_ABUND, CH_DEM_WRITE_RESULTS
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;     Ver.0.1, 9-Aug-2019, Peter Young
;     Ver.0.2, 30-Oct-2024, Peter Young
;       Modified how the output results are printed to the IDL window.
;     Ver.0.3, 07-Nov-2024, Peter Young
;       Added nsim= and mcmc_savefile= inputs; added additional tags to output.
;-


chck=have_proc('mcmc_dem')
IF chck EQ 0 THEN BEGIN
  print,'% CH_DEM_MCMC:  This routine makes use of the PINTofALE routine "mcmc_dem.pro", but you do not have this installed.'
  print,'                Please check the website below:'
  print,'                   https://hea-www.harvard.edu/PINTofALE/'
  print,'                A tar file containing PINTofALE is available at :'
  print,'                   http://hea-www.harvard.edu/PINTofALE/PoA_current.tar.gz'
  return,-1
ENDIF


IF n_elements(lpress) EQ 0 THEN BEGIN
  IF n_elements(ldens) EQ 0 THEN ldens=9.0
ENDIF 

;
; This is where all of the MCMC output parameters are stored.
;
IF n_elements(mcmc_savefile) EQ 0 THEN mcmc_savefile='mcmc.sav'

;
; The following automatically works out the temperature range by
; considering +/- 0.15 either side of log Tmax of each ion.
;
IF n_elements(dlogt) EQ 0 THEN dlogt=0.10
IF n_elements(ltemp) EQ 0 THEN BEGIN 
  ltemp_min=7.0
  ltemp_max=5.0
  FOR i=0,n_elements(line_data)-1 DO BEGIN
    ltmax=ch_tmax(line_data[i].ion,/log)
    IF ltmax-0.15 LT ltemp_min THEN ltemp_min=ltmax-0.15
    IF ltmax+0.15 GT ltemp_max THEN ltemp_max=ltmax+0.15
  ENDFOR 
  nt=round((ltemp_max-ltemp_min)/dlogt)+1
  ltemp=findgen(nt)*dlogt + ltemp_min
ENDIF ELSE BEGIN
  nt=n_elements(ltemp)
ENDELSE 


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
  print,'% CH_DEM_MCMC: The input LINE_DATA structure does not contain observed intensities. Please run '
  print,'               ch_dem_read_line_ids to load intensities from an input file. Returning...'
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
  print,'% CH_DEM_MCMC: error found. Returning...'
  return,-1
ENDIF 


;
; Element abundance information is stored in ABSTR and AB_REF.
;
IF n_elements(abund_file) EQ 0 THEN abund_file=!abund_file
chck=file_info(abund_file)
IF NOT chck.exists THEN BEGIN
  print,'% CH_DEM_MCMC: could not find the element abundance file. Returning...'
  return,-1
ENDIF 
abstr=ch_dem_process_abund(ld,abund_file=abund_file,ab_elt_fix=ab_elt_fix,fixed_abund=fixed_abund)
k=where(abstr.type EQ 1)
ab_ref=abstr[k[0]].abund
;
ab_type=abstr[ld.ab_ind].type
ind_ld_fit=where(ab_type GE 1,n_ld)
IF n_ld EQ 0 THEN return,-1
ld_fit=ld[ind_ld_fit]
nfit=n_elements(ld_fit)


;
; Need to multiply by 10^23 for the MCMC routines.
;
emis=dblarr(nt,nfit)
FOR i=0,nfit-1 DO BEGIN
  emis[*,i]=ld_fit[i].contrib*1e23
ENDFOR 


z=ld_fit.element_num
wvl=ld_fit.wvl
flx=ld_fit.int
fsigma=ld_fit.err
;
; Add an additional component to the error if interr_scale is given.
;
IF n_elements(interr_scale) NE 0 THEN BEGIN
  IF interr_scale LT 0 OR interr_scale GT 1 THEN BEGIN
    print,'% CH_DEM_GAUSS_FIT:  the input interr_scale should_fit take a value between 0 and 1. Returning...'
    return,-1
  ENDIF ELSE BEGIN
    fsigma=sqrt(fsigma^2 + (flx*interr_scale)^2 )
  ENDELSE
ENDIF


nhne=proton_dens(ltemp,/hydrogen,abund_file=abund_file,ioneq_file=ioneq_file)


;
; Set an initial guess for the DEM. I assume a Gaussian centered on
; the mid-point of ltemp, and with a sigma set to the ltemp range
; divided by four. Note that MCMC uses DEM*T rather than DEM.
;
IF n_elements(dem) EQ 0 THEN BEGIN
  dem=dblarr(nt)+1e22
  ;; lt0=mean(ltemp)
  ;; sig=(max(ltemp)-min(ltemp))/4.
  ;; diffem=1e22*exp(-(ltemp-lt0)^2/sig^2)*10.^ltemp
ENDIF ELSE BEGIN
  IF n_elements(dem) NE nt THEN BEGIN
    print,'% CH_DEM_MCMC: The size of the DEM input does not match the temperature array size.'
    print,'               If you did not specify LTEMP=, then run the routine again with the '
    print,'               LTEMP= optional output to see the temperatures for which the DEM needs'
    print,'               to be defined.'
    return,-1
  ENDIF 
ENDELSE
diffem=dem*10.^ltemp/nhne


;
; Need the abundance array (ab) for mcmc_dem, and set up the "abrng"
; input that allows abundances to vary. Only the "type 2" elements are
; allowed to vary.
;
; Note that I'm truncating to 30 elements (CHIANTI gives 50 by default).
;
read_abund,abund_file,ab,ref
ab=ab[0:29]
nab=n_elements(ab)
abrng=fltarr(nab,2)
abrng[*,0]=ab
abrng[*,1]=ab
k=where(abstr.type EQ 2,nk)
IF nk NE 0 THEN BEGIN
  j=abstr[k].elt_num-1
  abrng[j,0]=ab[j]/10.
  abrng[j,1]=ab[j]*10.
ENDIF 


;
; This is the call to the PINTofALE routine mcmc_dem. 
;
result=mcmc_dem(wvl, flx, emis, fsigma=fsigma, z=z, $
                chidir=!xuvtop, logt=ltemp, /noph, nhne=nhne, $
                diffem=diffem, abund=ab, savfil=mcmc_savefile, $
                abrng=abrng,aberr=aberr, nsim=nsim, simdem=simdem, $
                simprb=simprb,lscal=lscal,demerr=demerr)

dem=result/10.^ltemp*nhne
FOR i=0,1 DO demerr[*,i]=demerr[*,i]/10.^ltemp*nhne


;
; 'ab' contains updated abundances for 'type 2' elements, so load
; these into abstr
;
; I'm not updating the error as mcmc gives lower and upper
; confidence bounds (in the aberr output) rather than 1-sigma
; errors. Instead I just print the bounds to the screen.
;
k=where(abstr.type EQ 2,nk)
IF nk NE 0 THEN BEGIN
  abstr[k].abund=ab[abstr[k].elt_num-1]
  err_minus=abs(abstr[k].abund-aberr[abstr[k].elt_num-1,0])
  err_plus=abs(abstr[k].abund-aberr[abstr[k].elt_num-1,1])
  abstr[k].error=mean([err_minus,err_plus])
  abstr[k].ratio=abstr[k].ratio/ab_ref
ENDIF

;
; 'type 0' elements are not part of the minimization procedure but
; their updated abundances are obtained algebraically from the DEM.
;
; Note that I have to use ld here (not ld_all or ld_fit)
;
k=where(abstr.type EQ 0,nk)
IF nk NE 0 THEN BEGIN
  FOR i=0,nk-1 DO BEGIN
    j=k[i]
    ii=where(ld.element_num EQ abstr[j].elt_num)
    contrib=ld[ii].contrib
    line_int=ld[ii].int
    line_err=ld[ii].err
    IF n_elements(interr_scale) NE 0 THEN line_err=sqrt(line_err^2 + (line_int*interr_scale)^2)
    abstr[j].abund=line_int / $
                   total( contrib*dem*10.^ltemp*dlogt*alog(10.) )
    abstr[j].error=abstr[j].abund * line_err / line_int
    abstr[j].ratio=abstr[j].abund/ab_ref
  ENDFOR 
ENDIF





;
; Compute the T_eff of each line in LINE_DATA, and populate the ab_ind
; and model_int tags of LINE_DATA.
;
nab=n_elements(abstr)
FOR i=0,nl-1 DO BEGIN
  contrib_fn=ld_all[i].contrib
  func=contrib_fn*dem*10.^ltemp
  getmax=max(func,imax)
  ld_all[i].logt_eff=ltemp[imax]
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
; Write the list of lines to the IDL window with observed and model intensities.
;
ld_fit=ld_all[ind_ld_fit]
IF NOT keyword_set(quiet) THEN ch_dem_write_results,ld_fit, abstr

;
; I follow the prescription in Sect. 5 of the MCMC online manual for obtaining
; the reduced chi-square proxy.
;
jnk=varsmooth(fltarr(nt),lscal,nueff=nueff)
chi2_proxy=2*simprb[nsim]/(float(nfit)-nueff)


;
; Create the output structure.
;
IF n_elements(ldens) EQ 0 THEN ldens=-1.
IF n_elements(lpress) EQ 0 THEN lpress=-1.
IF n_elements(interr_scale) EQ 0 THEN interr_scale=-1.
output={method: 'mcmc', $
        ltemp: ltemp, $
        dem: dem, $
        line_data: ld_all, $
        abstr: abstr, $
        interr_scale: interr_scale, $
        log_dens: ldens, $
        log_press: lpress, $
        nsim: nsim, $
        simdem: simdem, $
        demerr: demerr, $
        simprb: simprb, $
        chi2_proxy: chi2_proxy, $
        time_stamp: systime()}




return,output

END
