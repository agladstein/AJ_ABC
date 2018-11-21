cd /home/nu_agladstein/macsswig_simsaj/workflow/runs/macsswig_simsaj_genome_1513132007/workflow/macsswig_simsaj_genome_1513132007
for FILE in `pegasus-analyzer | grep "error file:" | sed 's/.*: //' | sed 's/\.err/.sub/'`; do cat $FILE | grep transfer_input_files | perl -p -i -e 's/.*(macsargs_[0-9]*.txt).*/$1/'; done | sort -u

echo "macsargs_1696.txt
macsargs_1722.txt
macsargs_1848.txt
macsargs_2375.txt
macsargs_2845.txt
macsargs_322.txt
macsargs_4394.txt
macsargs_4501.txt
macsargs_4922.txt" | cut -d "_" -f2 | cut -d "." -f1 | xargs -I % macss_env/bin/python gen_macsargs_AJmodel1.py % full prior 0 retry

scp retry/* nu_agladstein@submit-4.chtc.wisc.edu:/home/nu_agladstein/macsswig_simsaj/workflow/runs/macsswig_simsaj_genome_1513132007/scratch/macsswig_simsaj_genome_1513132007/

cd /home/nu_agladstein/macsswig_simsaj/workflow
pegasus-run /home/nu_agladstein/macsswig_simsaj/workflow/runs/macsswig_simsaj_genome_1513132007/workflow/macsswig_simsaj_genome_1513132007

###############################################

cd /home/nu_agladstein/macsswig_simsaj/workflow/runs/macsswig_simsaj_genome_1513297199/workflow/macsswig_simsaj_genome_1513297199
for FILE in `pegasus-analyzer | grep "error file:" | sed 's/.*: //' | sed 's/\.err/.sub/'`; do cat $FILE | grep transfer_input_files | perl -p -i -e 's/.*(macsargs_[0-9]*.txt).*/$1/'; done | sort -u

echo "macsargs_1221.txt
macsargs_145.txt
macsargs_154.txt
macsargs_2098.txt
macsargs_2315.txt
macsargs_315.txt
macsargs_4428.txt
macsargs_4581.txt
macsargs_925.txt
macsargs_949.txt" | cut -d "_" -f2 | cut -d "." -f1 | xargs -I % macss_env/bin/python gen_macsargs_AJmodel1.py % full prior 0 retry

scp retry/* nu_agladstein@submit-4.chtc.wisc.edu:/home/nu_agladstein/macsswig_simsaj/workflow/runs/macsswig_simsaj_genome_1513297199/scratch/macsswig_simsaj_genome_1513297199/

cd /home/nu_agladstein/macsswig_simsaj/workflow
pegasus-run /home/nu_agladstein/macsswig_simsaj/workflow/runs/macsswig_simsaj_genome_1513297199/workflow/macsswig_simsaj_genome_1513297199

###############################################

cd /home/nu_agladstein/macsswig_simsaj/workflow/runs/macsswig_simsaj_genome_1513454037/workflow/macsswig_simsaj_genome_1513454037
echo "macsargs_574.txt" | cut -d "_" -f2 | cut -d "." -f1 | xargs -I % macss_env/bin/python gen_macsargs_AJmodel1.py % full prior 0 retry
scp retry/* nu_agladstein@submit-4.chtc.wisc.edu:/home/nu_agladstein/macsswig_simsaj/workflow/runs/macsswig_simsaj_genome_1513454037/scratch/macsswig_simsaj_genome_1513454037/

###############################################

cd /home/nu_agladstein/macsswig_simsaj/workflow/runs/macsswig_simsaj_genome_1513535687/workflow/macsswig_simsaj_genome_1513535687
for FILE in `pegasus-analyzer | grep "error file:" | sed 's/.*: //' | sed 's/\.err/.sub/'`; do cat $FILE | grep transfer_input_files | perl -p -i -e 's/.*(macsargs_[0-9]*.txt).*/$1/'; done | sort -u

echo "macsargs_1263.txt
macsargs_4399.txt
macsargs_439.txt
macsargs_4492.txt
macsargs_4705.txt
macsargs_4998.txt
macsargs_79.txt
macsargs_2324.txt
macsargs_285.txt" | cut -d "_" -f2 | cut -d "." -f1 | xargs -I % macss_env/bin/python gen_macsargs_AJmodel1.py % full prior 0 retry
scp retry/* nu_agladstein@submit-4.chtc.wisc.edu:/home/nu_agladstein/macsswig_simsaj/workflow/runs/macsswig_simsaj_genome_1513535687/scratch/macsswig_simsaj_genome_1513535687

cd /home/nu_agladstein/macsswig_simsaj/workflow
pegasus-run /home/nu_agladstein/macsswig_simsaj/workflow/runs/macsswig_simsaj_genome_1513535687/workflow/macsswig_simsaj_genome_1513535687

###############################################

cd /home/nu_agladstein/macsswig_simsaj/workflow/runs/macsswig_simsaj_genome_1514587139/workflow/macsswig_simsaj_genome_1514587139
for FILE in `pegasus-analyzer | grep "error file:" | sed 's/.*: //' | sed 's/\.err/.sub/'`; do cat $FILE | grep transfer_input_files | perl -p -i -e 's/.*(macsargs_[0-9]*.txt).*/$1/'; done | sort -u

echo "macsargs_1036.txt
macsargs_1089.txt
macsargs_1643.txt
macsargs_288.txt" | cut -d "_" -f2 | cut -d "." -f1 | xargs -I % macss_env/bin/python gen_macsargs_AJmodel1.py % full prior 0 retry
scp retry/* nu_agladstein@submit-4.chtc.wisc.edu:/home/nu_agladstein/macsswig_simsaj/workflow/runs/macsswig_simsaj_genome_1514587139/scratch/macsswig_simsaj_genome_1514587139

cd /home/nu_agladstein/macsswig_simsaj/workflow
pegasus-run /home/nu_agladstein/macsswig_simsaj/workflow/runs/macsswig_simsaj_genome_1514587139/workflow/macsswig_simsaj_genome_1514587139

###############################################

cd /home/nu_agladstein/macsswig_simsaj/workflow/runs/macsswig_simsaj_genome_1514824713/workflow/macsswig_simsaj_genome_1514824713
for FILE in `pegasus-analyzer | grep "error file:" | sed 's/.*: //' | sed 's/\.err/.sub/'`; do cat $FILE | grep transfer_input_files | perl -p -i -e 's/.*(macsargs_[0-9]*.txt).*/$1/'; done | sort -u

echo "macsargs_1647.txt
macsargs_4973.txt
macsargs_938.txt
macsargs_957.txt" | cut -d "_" -f2 | cut -d "." -f1 | xargs -I % macss_env/bin/python gen_macsargs_AJmodel1.py % full prior 0 retry
scp retry/* nu_agladstein@submit-4.chtc.wisc.edu:/home/nu_agladstein/macsswig_simsaj/workflow/runs/macsswig_simsaj_genome_1514824713/scratch/macsswig_simsaj_genome_1514824713

cd /home/nu_agladstein/macsswig_simsaj/workflow
pegasus-run /home/nu_agladstein/macsswig_simsaj/workflow/runs/macsswig_simsaj_genome_1514824713/workflow/macsswig_simsaj_genome_1514824713

###############################################

cd /home/nu_agladstein/macsswig_simsaj/workflow/runs/macsswig_simsaj_genome_1515003738/workflow/macsswig_simsaj_genome_1515003738
for FILE in `pegasus-analyzer | grep "error file:" | sed 's/.*: //' | sed 's/\.err/.sub/'`; do cat $FILE | grep transfer_input_files | perl -p -i -e 's/.*(macsargs_[0-9]*.txt).*/$1/'; done | sort -u

echo "macsargs_1123.txt
macsargs_1413.txt
macsargs_1704.txt
macsargs_1753.txt
macsargs_1805.txt
macsargs_2351.txt
macsargs_2948.txt
macsargs_3209.txt
macsargs_3609.txt
macsargs_4168.txt
macsargs_4709.txt
macsargs_559.txt" | cut -d "_" -f2 | cut -d "." -f1 | xargs -I % macss_env/bin/python gen_macsargs_AJmodel1.py % full prior 0 retry
scp retry/* nu_agladstein@submit-4.chtc.wisc.edu:/home/nu_agladstein/macsswig_simsaj/workflow/runs/macsswig_simsaj_genome_1515003738/scratch/macsswig_simsaj_genome_1515003738

cd /home/nu_agladstein/macsswig_simsaj/workflow
pegasus-run /home/nu_agladstein/macsswig_simsaj/workflow/runs/macsswig_simsaj_genome_1515003738/workflow/macsswig_simsaj_genome_1515003738


###############################################

cd /home/nu_agladstein/macsswig_simsaj/workflow/runs/macsswig_simsaj_genome_1513297199/workflow/macsswig_simsaj_genome_1513297199
for FILE in `pegasus-analyzer | grep "error file:" | sed 's/.*: //' | sed 's/\.err/.sub/'`; do cat $FILE | grep transfer_input_files | perl -p -i -e 's/.*(macsargs_[0-9]*.txt).*/$1/'; done | sort -u


###############################################

cd /home/nu_agladstein/macsswig_simsaj/workflow/runs/macsswig_simsaj_genome_1515003738/workflow/macsswig_simsaj_genome_1515003738
for FILE in `pegasus-analyzer | grep "error file:" | sed 's/.*: //' | sed 's/\.err/.sub/'`; do cat $FILE | grep transfer_input_files | perl -p -i -e 's/.*(macsargs_[0-9]*.txt).*/$1/'; done | sort -u

echo "macsargs_1123.txt
macsargs_1413.txt
macsargs_1704.txt
macsargs_1753.txt
macsargs_1805.txt
macsargs_2351.txt
macsargs_2948.txt
macsargs_3209.txt
macsargs_3609.txt
macsargs_4168.txt
macsargs_4709.txt
macsargs_559.txt" | cut -d "_" -f2 | cut -d "." -f1 | xargs -I % macss_env/bin/python gen_macsargs_AJmodel1.py % full prior 0 retry
scp retry/* nu_agladstein@submit-4.chtc.wisc.edu:/home/nu_agladstein/macsswig_simsaj/workflow/runs/macsswig_simsaj_genome_1515003738/scratch/macsswig_simsaj_genome_1515003738

cd /home/nu_agladstein/macsswig_simsaj/workflow
pegasus-run /home/nu_agladstein/macsswig_simsaj/workflow/runs/macsswig_simsaj_genome_1515003738/workflow/macsswig_simsaj_genome_1515003738

###############################################

cd /home/nu_agladstein/macsswig_simsaj/workflow/runs/macsswig_simsaj_genome_1515190460/workflow/macsswig_simsaj_genome_1515190460
for FILE in `pegasus-analyzer | grep "error file:" | sed 's/.*: //' | sed 's/\.err/.sub/'`; do cat $FILE | grep transfer_input_files | perl -p -i -e 's/.*(macsargs_[0-9]*.txt).*/$1/'; done | sort -u


###############################################

cd /local-scratch/agladstein/workflows/macsswig_simsaj_genome_1512886866/workflow/macsswig_simsaj_genome_1512886866


###############################################

cd /local-scratch/agladstein/workflows/macsswig_simsaj_genome_1514911098/workflow/macsswig_simsaj_genome_1514911098
for SH_FILE in `pegasus-analyzer | grep "error file:" | sed 's/.*: //' | sed 's/\.err.*/.sh/'`; do cat $SH_FILE | grep PWD/macsargs_ | perl -p -i -e 's/.*(macsargs_[0-9]*.txt).*/$1/'; done | sort -u

BAD=$(echo "
macsargs_1093.txt
macsargs_1217.txt
macsargs_1297.txt
macsargs_1348.txt
macsargs_1797.txt
macsargs_1847.txt
macsargs_1957.txt
macsargs_204.txt
macsargs_2295.txt
macsargs_2701.txt
macsargs_2998.txt
macsargs_3261.txt
macsargs_3317.txt
macsargs_3444.txt
macsargs_3498.txt
macsargs_3745.txt
macsargs_4298.txt
macsargs_4406.txt
macsargs_453.txt
macsargs_4599.txt
macsargs_4793.txt
macsargs_4903.txt
macsargs_5102.txt
macsargs_5370.txt
macsargs_544.txt
macsargs_5599.txt
macsargs_6360.txt
macsargs_6648.txt
macsargs_6683.txt
macsargs_6898.txt
macsargs_6945.txt
macsargs_7000.txt
macsargs_7099.txt
macsargs_7137.txt
macsargs_7139.txt
macsargs_7140.txt
macsargs_752.txt
macsargs_773.txt
macsargs_7935.txt
macsargs_7992.txt
macsargs_8152.txt
macsargs_8571.txt
macsargs_8999.txt
macsargs_9298.txt
macsargs_9362.txt
macsargs_9603.txt
macsargs_9853.txt")
for FILE in $BAD; do find /stash2/user/agladstein/public/macsswig_simsaj_genome_1514911098 -name $FILE; done
/stash2/user/agladstein/public/macsswig_simsaj_genome_1514911098/00/1E/macsargs_1297.txt
/stash2/user/agladstein/public/macsswig_simsaj_genome_1514911098/00/3E/macsargs_8571.txt
/stash2/user/agladstein/public/macsswig_simsaj_genome_1514911098/00/14/macsargs_8999.txt
/stash2/user/agladstein/public/macsswig_simsaj_genome_1514911098/00/7E/macsargs_9853.txt

echo "macsargs_1297.txt
macsargs_8571.txt
macsargs_8999.txt
macsargs_9853.txt" | cut -d "_" -f2 | cut -d "." -f1 | xargs -I % macss_env/bin/python gen_macsargs_AJmodel1.py % full prior 0 retry

scp retry/macsargs_1297.txt agladstein@login02.osgconnect.net:/stash2/user/agladstein/public/macsswig_simsaj_genome_1514911098/00/1E/macsargs_1297.txt
scp retry/macsargs_8571.txt agladstein@login02.osgconnect.net:/stash2/user/agladstein/public/macsswig_simsaj_genome_1514911098/00/3E/macsargs_8571.txt
scp retry/macsargs_8999.txt agladstein@login02.osgconnect.net:/stash2/user/agladstein/public/macsswig_simsaj_genome_1514911098/00/14/macsargs_8999.txt
scp retry/macsargs_9853.txt agladstein@login02.osgconnect.net:/stash2/user/agladstein/public/macsswig_simsaj_genome_1514911098/00/7E/macsargs_9853.txt

pegasus-run /local-scratch/agladstein/workflows/macsswig_simsaj_genome_1514911098/workflow/macsswig_simsaj_genome_1514911098


cd /local-scratch/agladstein/workflows/macsswig_simsaj_genome_1514911098/workflow/macsswig_simsaj_genome_1514911098/
for SH_FILE in `pegasus-analyzer | grep "error file:" | sed 's/.*: //' | sed 's/\.err.*/.sh/'`; do cat $SH_FILE | grep PWD/macsargs_ | perl -p -i -e 's/.*(macsargs_[0-9]*.txt).*/$1/'; done | sort -u


BAD=$(echo "macsargs_1582.txt
macsargs_2037.txt
macsargs_2157.txt
macsargs_2323.txt
macsargs_392.txt
macsargs_5574.txt
macsargs_5726.txt
macsargs_6603.txt
macsargs_7212.txt
macsargs_8544.txt")
for FILE in $BAD; do find /stash2/user/agladstein/public/macsswig_simsaj_genome_1514911098 -name $FILE; done
/stash2/user/agladstein/public/macsswig_simsaj_genome_1514911098/00/69/macsargs_1582.txt
/stash2/user/agladstein/public/macsswig_simsaj_genome_1514911098/00/46/macsargs_2037.txt
/stash2/user/agladstein/public/macsswig_simsaj_genome_1514911098/00/87/macsargs_2157.txt
/stash2/user/agladstein/public/macsswig_simsaj_genome_1514911098/00/BC/macsargs_2323.txt
/stash2/user/agladstein/public/macsswig_simsaj_genome_1514911098/00/48/macsargs_392.txt
/stash2/user/agladstein/public/macsswig_simsaj_genome_1514911098/00/B6/macsargs_5574.txt
/stash2/user/agladstein/public/macsswig_simsaj_genome_1514911098/00/1B/macsargs_5726.txt
/stash2/user/agladstein/public/macsswig_simsaj_genome_1514911098/00/88/macsargs_6603.txt
/stash2/user/agladstein/public/macsswig_simsaj_genome_1514911098/00/BC/macsargs_7212.txt
/stash2/user/agladstein/public/macsswig_simsaj_genome_1514911098/00/3F/macsargs_8544.txt

echo ${BAD} | cut -d "_" -f2 | cut -d "." -f1 | xargs -I % macss_env/bin/python gen_macsargs_AJmodel1.py % full prior 0 retry

scp retry/macsargs_1582.txt agladstein@login02.osgconnect.net:/stash2/user/agladstein/public/macsswig_simsaj_genome_1514911098/00/69/macsargs_1582.txt
scp retry/macsargs_2037.txt agladstein@login02.osgconnect.net:/stash2/user/agladstein/public/macsswig_simsaj_genome_1514911098/00/46/macsargs_2037.txt
scp retry/macsargs_2157.txt agladstein@login02.osgconnect.net:/stash2/user/agladstein/public/macsswig_simsaj_genome_1514911098/00/87/macsargs_2157.txt
scp retry/macsargs_2323.txt agladstein@login02.osgconnect.net:/stash2/user/agladstein/public/macsswig_simsaj_genome_1514911098/00/BC/macsargs_2323.txt
scp retry/macsargs_392.txt agladstein@login02.osgconnect.net:/stash2/user/agladstein/public/macsswig_simsaj_genome_1514911098/00/48/macsargs_392.txt
scp retry/macsargs_5574.txt agladstein@login02.osgconnect.net:/stash2/user/agladstein/public/macsswig_simsaj_genome_1514911098/00/B6/macsargs_5574.txt
scp retry/macsargs_5726.txt agladstein@login02.osgconnect.net:/stash2/user/agladstein/public/macsswig_simsaj_genome_1514911098/00/1B/macsargs_5726.txt
scp retry/macsargs_6603.txt agladstein@login02.osgconnect.net:/stash2/user/agladstein/public/macsswig_simsaj_genome_1514911098/00/88/macsargs_6603.txt
scp retry/macsargs_7212.txt agladstein@login02.osgconnect.net:/stash2/user/agladstein/public/macsswig_simsaj_genome_1514911098/00/BC/macsargs_7212.txt
scp retry/macsargs_8544.txt agladstein@login02.osgconnect.net:/stash2/user/agladstein/public/macsswig_simsaj_genome_1514911098/00/3F/macsargs_8544.txt

pegasus-run /local-scratch/agladstein/workflows/macsswig_simsaj_genome_1514911098/workflow/macsswig_simsaj_genome_1514911098

###############################################

cd /local-scratch/agladstein/workflows/macsswig_simsaj_genome_1512886866/workflow/macsswig_simsaj_genome_1512886866/
for SH_FILE in `pegasus-analyzer | grep "error file:" | sed 's/.*: //' | sed 's/\.err.*/.sh/'`; do cat $SH_FILE | grep PWD/macsargs_ | perl -p -i -e 's/.*(macsargs_[0-9]*.txt).*/$1/'; done | sort -u

echo "macsargs_1030.txt
macsargs_1448.txt
macsargs_1597.txt
macsargs_1689.txt
macsargs_255.txt
macsargs_3206.txt
macsargs_3477.txt
macsargs_3488.txt
macsargs_3859.txt
macsargs_4371.txt
macsargs_4459.txt
macsargs_4972.txt
macsargs_5478.txt
macsargs_5887.txt
macsargs_5942.txt
macsargs_6692.txt
macsargs_6898.txt
macsargs_6940.txt
macsargs_6981.txt
macsargs_7187.txt
macsargs_7440.txt
macsargs_7630.txt
macsargs_844.txt
macsargs_9137.txt
macsargs_918.txt
macsargs_9313.txt"


###############################################


