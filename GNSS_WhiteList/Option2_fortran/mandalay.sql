CREATE VIEW mandalay AS
SELECT
  statid,date,time,an_sens_obs,timeslot,lat,lon,stalt,orography,obsvalue,fg_depar,datum_status.active@body,datum_status.blacklisted@body,an_depar
FROM hdr,body,errstat,index,modsurf
WHERE (obstype == 1 && varno=128 )
