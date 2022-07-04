echo 'starting the Analysis - hope all the config files had the right paths!'
cd /home/daq/Desktop/Luca/Analisi_RSD/RSD2_Analysis
rm nohup.out
make -f Makefile_Digitizer_wfm
nohup nice ./Analysis_Digitizer_wfm
cd /home/daq/Desktop/Luca/DiscordBots/RSD2
python3 Analysis_alert.py
cd /home/daq/Desktop/Luca/Analisi_RSD/RSD2_Analysis
nohup nice root -lq Process_RSD2.cpp
cd W3/
for i in {1..6}
do
   pdftoppm -png "c${i}.pdf" > "c${i}.png"
done
rm *pdf
cd /home/daq/Desktop/Luca/DiscordBots/RSD2
python3 Attach_results.py
cd /home/daq/Desktop/Luca/Analisi_RSD/RSD2_Analysis
mv RunXX.root RunXY.root

