mkdir -p C_immitis/fastq
for n in $(cut -f1 Cimm_sra_samples.csv); do ln -s ../../fastq/${n}_1.fastq.gz C_immitis/fastq/${n}_1.fastq.gz; ln -s ../../fastq/${n}_2.fastq.gz C_immitis/fastq/${n}_2.fastq.gz; done

mkdir -p C_posadasii/fastq
for n in $(cut -f1 Cpos_sra_samples.csv); do ln -s ../../fastq/${n}_1.fastq.gz C_posadasii/fastq/${n}_1.fastq.gz; ln -s ../../fastq/${n}_2.fastq.gz C_posadasii/fastq/${n}_2.fastq.gz; done
