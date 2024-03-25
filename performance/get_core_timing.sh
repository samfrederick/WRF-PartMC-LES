
ARCHIVE_PATH=/data/keeling/a/sf20/d/les_output/wrf-partmc
JOB=slurm-1972087

grep -o -P '.{0,19} on domain' $ARCHIVE_PATH/$JOB/rsl.error.00000 > $ARCHIVE_PATH/$JOB/clock.txt
grep -o -P '.{0,12}elapsed seconds' $ARCHIVE_PATH/$JOB/rsl.error.00000 > $ARCHIVE_PATH/$JOB/elapsed.txt

python core_timing.py $ARCHIVE_PATH/$JOB

rm $ARCHIVE_PATH/$JOB/clock.txt
rm $ARCHIVE_PATH/$JOB/elapsed.txt