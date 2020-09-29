
# Search for Ingress homogeneous fields

This program searches a portal export file (CSV) for homogeneous fields
of the given depth.

## Example

For example, how many H2s in this demo.csv file of five portals?

    Name,Latitude,Longitude
    "Aquatic West Off leash Area ",-27.481769,153.21646
    "A Place By The Bay ",-27.481458,153.218233
    "Annex Aquatica",-27.481773,153.217117
    "Aquatic Paradise Playground ",-27.481921,153.217877
    "North Annex Aquatica",-27.481626,153.217501

The h program finds these two solutions:

    $ ./h -d 2 demo.csv
    Found solution, area: 0.355092
    [{"type":"polygon","latLngs":[{"lat":-27.481458,"lng":153.218233},{"lat":-27.481769,"lng":153.21646},{"lat":-27.481921,"lng":153.217877}],"color":"#a24ac3"},{"type":"polygon","latLngs":[{"lat":-27.481773,"lng":153.217117},{"lat":-27.481458,"lng":153.218233},{"lat":-27.481769,"lng":153.21646}],"color":"#a24ac3"},{"type":"polygon","latLngs":[{"lat":-27.481773,"lng":153.217117},{"lat":-27.481769,"lng":153.21646},{"lat":-27.481921,"lng":153.217877}],"color":"#a24ac3"},{"type":"polygon","latLngs":[{"lat":-27.481773,"lng":153.217117},{"lat":-27.481921,"lng":153.217877},{"lat":-27.481458,"lng":153.218233}],"color":"#a24ac3"}]
    Found solution, area: 0.105709
    [{"type":"polygon","latLngs":[{"lat":-27.481458,"lng":153.218233},{"lat":-27.481769,"lng":153.21646},{"lat":-27.481773,"lng":153.217117}],"color":"#a24ac3"},{"type":"polygon","latLngs":[{"lat":-27.481626,"lng":153.217501},{"lat":-27.481458,"lng":153.218233},{"lat":-27.481769,"lng":153.21646}],"color":"#a24ac3"},{"type":"polygon","latLngs":[{"lat":-27.481626,"lng":153.217501},{"lat":-27.481769,"lng":153.21646},{"lat":-27.481773,"lng":153.217117}],"color":"#a24ac3"},{"type":"polygon","latLngs":[{"lat":-27.481626,"lng":153.217501},{"lat":-27.481773,"lng":153.217117},{"lat":-27.481458,"lng":153.218233}],"color":"#a24ac3"}]

## Input and output

The input CSV databases can be constructed using the CSV Export plugin installed into IITC:
 * https://iitc.app/
 * https://raw.githubusercontent.com/Zetaphor/IITC-Ingress-Portal-CSV-Export/master/ingress_export.js

The JSON output can be pasted into the DrawTools feature of IITC (DrawTools Opt →  Paste Drawn Items)

## Build options

When compiling the *h* program, you may try to make it faster with these build options:

 * NDEBUG - disable all debugging which makes structures smaller, and the program runs faster
 * OPENMP - if your compiler supports OpenMP, the program will run using as many CPUs as it can

## Other tools

 * csvtodrawjson - converts a CSV file into DrawTools markers (JSON)
 * csvview.html - a tool to visualise CSV and JSON polygons; used when developing h
 * maketri.py - a python script to generate a CSV file of portals of a perfect Hn field

