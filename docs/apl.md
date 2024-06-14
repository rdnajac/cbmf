# RNAseq for apl

This document describes the steps to analyze RNAseq data for the APL project.

## Download data from AWS

Locations:
* s3://lab-aaf-ngs-data-archive/RNAseq/

20220807_DNMT3A_RHOA_PreLymphoma_APL/
20231127_Tet2Rhoa-sgID2-3-cloneB_APL/
20240424_Tet2RhoaG17V-AGX51/

mv 20240424_Tet2RhoaG17V-AGX51/ to 20240424_Tet2Rhoa-AGX51_APL/
aws s3 mv s3://lab-aaf-ngs-data-archive/RNAseq/20240424_Tet2RhoaG17V-AGX51/ s3://lab-aaf-ngs-data-archive/RNAseq/20240424_Tet2Rhoa-AGX51_APL/

check move status:
aws s3 ls s3://lab-aaf-ngs-data-archive/RNAseq/20240424_Tet2Rhoa-AGX51_APL/
aws s3 sync s3://lab-aaf-ngs-data-archive/RNAseq/20240424_Tet2RhoaG17V-AGX51/ .

### rename files

rename files to remove hyphens 

copy all .sh files from here:
s3://lab-aaf-ngs-data-archive/RNAseq/20220807_DNMT3A_RHOA_PreLymphoma_APL/

aws s3 ls s3://lab-aaf-ngs-data-archive/RNAseq/Tet2Rhoa-AGX51_APL/

aws s3 copy all files in pwd to that folder
aws s3 cp . s3://lab-aaf-ngs-data-archive/RNAseq/Tet2Rhoa-AGX51_APL/ --recursive

# make a checksums .md file for every file in the folder
for f in ./*; do md5sum $f >> checksums.md; done

cat *.md5
for f in ./*; do md5sum $f >> checksums.md5; done 

ubuntu@ip-172-31-32-180:~/apl2/bam$ ls ballgown_input/*
scp -r my-ec2:/home/ubuntu/apl2/bam/ballgown_input/* .
