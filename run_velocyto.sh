export PATH=/home/hieunguyen/samtools/bin:$PATH;
maindir=/home/hieunguyen/CRC1382/outdir/GSE192742/cellranger_output;
path_to_gtf=/media/hieunguyen/HD0/storage/build-mm10/reference_sources/gencode.vM23.primary_assembly.annotation.gtf;

for file in Sample_GSM5764245 Sample_GSM5764259 Sample_GSM5764288;do \
path_to_cellranger_output=$maindir/$file
velocyto run10x $path_to_cellranger_output $path_to_gtf -@ 20;done
