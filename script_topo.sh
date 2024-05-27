gmt blockmean bedrock.xyz -R11.70/14.8/45.2/47.3 -I15s+e/15s+e > Friuli.xyz
./convert_lonlat2utm.pl Friuli.xyz 33 > Friuli.utm
/Applications/Coreform-Cubit-2023.1.app/Contents/MacOS/Coreform-Cubit-2023.1 -nographics python3 read_topo.py
