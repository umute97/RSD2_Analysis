echo 'starting the Analysis - hope all the config files had the right paths!'
cd /home/daq/Desktop/Luca/Analisi_RSD/RSD2_Analysis
rm nohup.out
make -f Makefile_Digitizer_wfm
nohup nice ./Analysis_Digitizer_wfm
cd /home/daq/Desktop/Luca/DiscordBots/RSD2
python3 Analysis_alert.py
#nohup nice root -lq Process_RSD2.cpp #could be useful in the future
#cd W3/ 
cd /home/daq/Desktop/Luca/Analisi_RSD/RSD2_Analysis/test_newfiles
python3 test_rootfile.py
for i in {0..15}
do
   pdftoppm -png "c${i}.pdf" > "c${i}.png"
done
rm *pdf
cd /home/daq/Desktop/Luca/DiscordBots/RSD2
python3 Attach_results.py
cd /home/daq/Desktop/Luca/Analisi_RSD/RSD2_Analysis
mv RunXX.root RunXY.root

