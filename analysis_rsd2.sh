echo 'starting the Analysis - hope all the config files had the right paths!'
rm nohup.out
make -f Makefile_Digitizer_wfm
nohup nice ./Analysis_Digitizer_wfm
cd /home/daq/Desktop/Luca/DiscordBots/RSD2
python3 Analysis_alert.py
cd -
nohup nice root -lq Process_RSD2.cpp
cd /home/daq/Desktop/Luca/DiscordBots/RSD2
python3 Attach_results.py
cd -
rm *.pdf
mv RunXX.root Run XY.root

