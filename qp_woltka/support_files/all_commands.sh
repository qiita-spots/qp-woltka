cat ../S22205_S104_L001_R* > S22205_S104_L001.fastq.gz
cat ../S22282_S102_L001_R* > S22282_S102_L001.fastq.gz

bowtie2 -x ../../databases/wol/WoLmin -q S22282_S102_L001.fastq.gz -S S22282_S102_L001.sam --seed 42 --very-sensitive -k 16 --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-0.05" --no-head --no-unal
bowtie2 -x ../../databases/wol/WoLmin -q S22205_S104_L001.fastq.gz -S S22205_S104_L001.sam --seed 42 --very-sensitive -k 16 --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-0.05" --no-head --no-unal

xz -9 -c S22282_S102_L001.sam > S22282_S102_L001.sam.xz
xz -9 -c S22205_S104_L001.sam > S22205_S104_L001.sam.xz

woltka classify -i S22205_S104_L001.sam -o S22205_S104_L001.woltka-taxa --no-demux --lineage ../../databases/wol/WoLmin.tax --rank phylum,genus,species,free,none
woltka classify -i S22282_S102_L001.sam -o S22282_S102_L001.woltka-taxa --no-demux --lineage ../../databases/wol/WoLmin.tax --rank phylum,genus,species,free,none

woltka classify -i S22282_S102_L001.sam -o S22282_S102_L001.woltka-per-gene --no-demux -c ../../databases/wol/WoLmin.coords
woltka classify -i S22205_S104_L001.sam -o S22205_S104_L001.woltka-per-gene --no-demux -c ../../databases/wol/WoLmin.coords

woltka_merge --prep prep_information.txt --base $PWD --name species --glob "*.woltka-taxa/species.biom"
woltka_merge --prep prep_information.txt --base $PWD --name genus --glob "*.woltka-taxa/genus.biom"
woltka_merge --prep prep_information.txt --base $PWD --name phylum --glob "*.woltka-taxa/phylum.biom"
woltka_merge --prep prep_information.txt --base $PWD --name free --glob "*.woltka-taxa/free.biom"
woltka_merge --prep prep_information.txt --base $PWD --name none --glob "*.woltka-taxa/none.biom"
woltka_merge --prep prep_information.txt --base $PWD --name per-gene --glob "*.woltka-per-gene"

tar -cvf alignment.tar *.sam.xz
