
FUNCTION ch_dem_linear_compute_dem, ltemp, ltemp_nodes, node_values, scl=scl, em=em

dlogt=ltemp[1]-ltemp[0]
t=10.^ltemp

log_f=make_array(n_elements(ltemp),/double,value=-1.)
phi=dblarr(n_elements(ltemp))

nn=n_elements(ltemp_nodes)
p=node_values

IF n_elements(scl) EQ 0 THEN scl=1.0

;
; Here I linearly interpolate between the nodes to give the function log_f
;
FOR i=0,nn-2 DO BEGIN
  x1=ltemp_nodes[i]
  x2=ltemp_nodes[i+1]
  y1=alog10(p[i])
  y2=alog10(p[i+1])
  k=where(ltemp GE x1 AND ltemp LE x2)
  log_f[k]=ltemp[k]*(y2-y1)/(x2-x1) + (x2*y1-x1*y2)/(x2-x1)
ENDFOR 

;
; If /em was set, then the log_f values are log(EM) values, so I need
; to convert them to DEM (phi) values.
;
k=where(log_f NE -1.)
IF keyword_set(em) THEN BEGIN
  phi[k]=10.^log_f[k] * scl / (alog(10.) * t[k] * dlogt)
ENDIF ELSE BEGIN
  phi[k]=10.^log_f[k] * scl
ENDELSE

return,phi

END
