echo 'Downloading embl files'
wget ftp://ftp.sanger.ac.uk/pub/project/pathogens/Plasmodium/falciparum/3D7/3D7.latest_version/version3.1/2018/May_2018/Pf3D7_{01..14}_v3.embl.gz
gunzip *

echo 'Converting embl to fasta'
for chr in {01..14}; do
	python3 embl2fasta.py Pf3D7_${chr}_v3.embl Pf3D7_${chr}_v3.fasta
done

echo 'Concatenating per-chromosome files'
cat Pf3D7_*_v3.fasta > Pf3D7_v3_May18.fasta

echo 'Cleaning contig names'
sed -i "s/_v3//g" Pf3D7_v3_May18.fasta

rm Pf3D7_{01..14}_v3.embl
rm Pf3D7_{01..14}_v3.fasta
