echo 'starting the Analysis - I hope all the config files had the right paths!'
cd /home/daq/Desktop/Luca/Analisi_RSD/RSD2_Analysis
rm nohup.out
make -f Makefile_Digitizer_wfm
nohup nice ./Analysis_Digitizer_wfm
cd /home/daq/Desktop/Luca/DiscordBots/RSD2
python3 Analysis_alert.py
cd -

