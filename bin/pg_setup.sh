export GDAL_DATA=/usr/local/opt/gdal2/share/gdal/

psql -c "DROP DATABASE \"nhd\"" postgresql://localhost:5432/

psql -c "CREATE DATABASE \"nhd\"" postgresql://localhost:5432/

psql -c "CREATE EXTENSION postgis; CREATE EXTENSION postgis_topology;" postgresql://localhost:5432/nhd

ogr2ogr -overwrite -progress -f "PostGreSQL" PG:"host=localhost dbname=nhd" -nlt PROMOTE_TO_MULTI -lco "GEOMETRY_NAME=the_geom" -a_srs EPSG:5070 out_nhdp_latest/wbd_viz.gpkg wbd_viz

ogr2ogr -overwrite -progress -f "PostGreSQL" PG:"host=localhost dbname=nhd" -nlt PROMOTE_TO_MULTI -lco "GEOMETRY_NAME=the_geom" -a_srs EPSG:5070 out_nhdp_latest/wbd_viz.gpkg net

psql -c "CREATE INDEX ON net (levelpathi);" postgresql://localhost:5432/nhd

ogr2ogr -overwrite -progress -f "PostGreSQL" PG:"host=localhost dbname=nhd" -nlt PROMOTE_TO_MULTI -lco "GEOMETRY_NAME=the_geom" -a_srs EPSG:5070 wbd_viz.gpkg net_grouped

ogr2ogr -overwrite -progress -f "PostGreSQL" PG:"host=localhost dbname=nhd" -nlt PROMOTE_TO_MULTI -lco "GEOMETRY_NAME=the_geom" -a_srs EPSG:5070 out_nhdp_latest/wbd_viz.gpkg wbd_matched

ogr2ogr -overwrite -progress -f "PostGreSQL" PG:"host=localhost dbname=nhd" -nlt PROMOTE_TO_MULTI -lco "GEOMETRY_NAME=the_geom" -a_srs EPSG:5070 out_nhdp_latest/wbd_viz.gpkg wbd_grouped

ogr2ogr -overwrite -progress -f "PostGreSQL" PG:"host=localhost dbname=nhd" -nlt PROMOTE_TO_MULTI -lco "GEOMETRY_NAME=the_geom" -nln wbd_hu12 -t_srs EPSG:5070 WBD_National_GDB.gdb WBDHU12

ogr2ogr -overwrite -progress -f "PostGreSQL" PG:"host=localhost dbname=nhd" -nlt PROMOTE_TO_MULTI -lco "GEOMETRY_NAME=the_geom" -nln nhdplus_hu12 -t_srs EPSG:5070 data/nhdplus/NHDPlusNationalData/NHDPlusV21_National_Seamless.gdb/ HUC12