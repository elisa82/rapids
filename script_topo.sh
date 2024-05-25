gmt blockmean bedrock.xyz -R11.90/14.2/45.3/47.0 -I15s+e/15s+e > Friuli.xyz
./convert_lonlat2utm.pl Friuli.xyz 33 > Friuli.utm
/Applications/Coreform-Cubit-2023.1.app/Contents/MacOS/Coreform-Cubit-2023.1 -nographics python3 read_topo.py
