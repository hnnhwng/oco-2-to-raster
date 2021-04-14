# oco-2-to-raster
Convert OCO-2 data to single band variable raster

Data available here: https://co2.jpl.nasa.gov/#mission=OCO-2

### Setup
1. Place data in the same directory as script.
2. To create environment:
```
conda create --name <env> --file requirements.txt
conda activate <env>
```
3. Run script with 
```
python .\toGeoTiff.py
```

### Usage
```
    -r <number> or --resolution <number> to specify the size of each pixel, in meters. Default: 1000
    -v <string> or --variable <dir> to specify a variable from the dataset to extract. Default "xco2"
    -e <EPSG code> or --epsg <EPSG code> to specify an output projection. Default: 4326
    -b <min_lon, min_lat, max_lon, max_lat> or --bounds min_lon, min_lat, max_lon, max_lat> to specify an raster bounds. Default: "-45,45,-180,180"
    -o <dir> or --output <dir> to specify a path to an output directory. Default: ""
    -f <filename> or --file <filename> to specify an input file. Omit to process all files. Default: ""
    -a or --average to generate yearly and monthly averages. Default: False
    -h to display this message
```
