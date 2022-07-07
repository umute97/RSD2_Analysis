cd /home/daq/Desktop/Luca/Analisi_RSD/RSD2_Analysis/test_newfiles
for i in {1..15}
do
   pdftoppm -png "c${i}.pdf" > "c${i}.png"
done
rm *pdf
cd /home/daq/Desktop/Luca/DiscordBots/RSD2
python3 Attach_results.py
cd /home/daq/Desktop/Luca/Analisi_RSD/RSD2_Analysis
