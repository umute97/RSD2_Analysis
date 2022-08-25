echo 'starting the Analysis - hope all the config files had the right paths!'
cd /home/daq/Desktop/Luca/Analisi_RSD/RSD2_Analysis
rm nohup.out
make -f Makefile_Digitizer_wfm
nohup nice ./Analysis_Digitizer_wfm
cd /home/daq/Desktop/Luca/DiscordBots/RSD2
python3 Analysis_alert.py
#could be useful in the future
#nohup nice root -lq Process_RSD2.cpp 
#cd W3/ 
cd /home/daq/Desktop/Luca/Analisi_RSD/RSD2_Analysis/test_newfiles
#moves the file into the stats directory and renames it
python3 test_rootfile.py /home/daq/hdd8TB_bis/RSD2/stats_N_script/RunXX.root 1
#from here on no actions on the output root file are taken
for i in {0..15}
do
   pdftoppm -png "c${i}.pdf" > "c${i}.png"
done
rm *.pdf
cd /home/daq/Desktop/Luca/DiscordBots/RSD2
python3 Attach_results.py
cd /home/daq/Desktop/Luca/Analisi_RSD/RSD2_Analysis
echo 'done with the Analysis'
