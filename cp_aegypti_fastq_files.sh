# create parent directory where species data will be
export SPECIES_DIR=/media/jaxi/aedes_aegypti
mkdir "$SPECIES_DIR"

# within create raw_data and reference directory's
# for fastq files (reads) and gff3/fa (annotated genome and- 
# reference genome) files respectively
export SPECIES_DIR_RAW=$SPECIES_DIR/raw_data
export SPECIES_DIR_REF=$SPECIES_DIR/reference

mkdir "$SPECIES_DIR_RAW"
mkdir "$SPECIES_DIR_REF"

# create species named folders within the raw_data
# and reference folders
export SPECIES_DIR_RAW=$SPECIES_DIR_RAW/aedes_aegypti/
export SPECIES_DIR_REF=$SPECIES_DIR_REF/aedes_aegypti/

mkdir "$SPECIES_DIR_RAW"
mkdir "$SPECIES_DIR_REF"

# copy reference genome and genome annotation to species folder
cp /data/seqfiles_SW/index_files/aegypti.fa "$SPECIES_DIR_REF"
cp /media/jaxi/temp_bam/noriboaegypti.gff3 "$SPECIES_DIR_REF"

export PROJECT_BK=/data/Project_BROOKLYN/

cp $PROJECT_BK/Sample_JN01-050913/* "$SPECIES_DIR_RAW"
cp $PROJECT_BK/Sample_JN02-050913/* "$SPECIES_DIR_RAW"
