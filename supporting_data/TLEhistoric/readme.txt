README file for TLE data repository 

Space-track.org TLE as CSV file from query:

https://www.space-track.org/basicspacedata/query/class/gp_history/EPOCH/2023-01-01--2023-01-03/DECAY_DATE/null-val/MEAN_MOTION/>3/orderby/NORAD_CAT_ID/format/csv

NOTE: 
* Login via browser before trying query
* Any more than 3 days result in error to download (too big a file?)
* Above query explained: history of TLE's (GP) for epoch between 1/1/23 - 1/3/23, non-decayed objects, mean motion > 3 revs/day
* The CSV files in this folder is the result of the above query with the years changed


API info: https://www.space-track.org/documentation#/api
API info on gp_history object: https://www.space-track.org/basicspacedata/modeldef/class/gp_history/format/html
Data definition follows CCSDS Recommended Standard 502.0-B-2: https://public.ccsds.org/Pubs/502x0b2c1e2.pdf



To merge with DISCOSweb data (object type, size, mass, etc), see https://github.mit.edu/arclab/orbitalrisk_MC/tree/master/DISCOSweb/
DISCOSweb database: https://discosweb.esoc.esa.int/
