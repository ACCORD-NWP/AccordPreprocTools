SELECT
statid,rad2deg(lat),rad2deg(lon),stalt,obsvalue,fg_depar,datum_status.active,datum_status.blacklisted,datum_status.passive,datum_status.rejected
FROM  hdr,body
WHERE (codetype = 110 ) AND (fg_depar is not NULL)

