# ChIPseq for Kalay

This is a markdown file for the ChIPseq analysis of Kalay.

## Data

### BaseSpace
| run date | run name (basespace) | size (GB) | %Q30 | %PF | instrument | flow cell ID
|----------|----------------------|-----------|------|-----|------------|-------------
| 20230215 | HEL_H3K27ac_ChIPSeq1 | 27 | 94.49 | 92.21 | NB551203 | HC32JBGXN |
| 20230306 | HEL_STAT5_H3K27ac_2 | 37 | 92.74 | 89.86 | NS500289 | HW5Y3BGXN |
| 20230519 | CBP_p300 | 44 | 90.06 | 87.69 | NS500289 | HHLJNBGXT |
| 20240522 | KB_ChIPseq_4 | 10 | 94.47 | 93.50 | NB551203 | H7TMNAFX7 |

### AWS
Organize the data into a folder named `kalay/` with each seq date as a subfolder.

# make dirs for each seq date in aws s3 bucket s3://lab-aaf-ngs-data-archive/ChIPseq/kalay/
aws s3 ls s3://lab-aaf-ngs-data-archive/ChIPseq/kalay/
aws s3api put-object --bucket lab-aaf-ngs-data-archive --key ChIPseq/kalay/20230215/

# copy everything in pwd to s3://lab-aaf-ngs-data-archive/ChIPseq/kalay/20240522/
aws s3 cp . s3://lab-aaf-ngs-data-archive/ChIPseq/kalay/20240522/ --recursive

aws s3 cp . 3://lab-aaf-ngs-data-archive/RNAseq/20240424_Tet2Rhoa-AGX51_APL/ --recursive
# instead copy all .gz files 

save md5 of every ./ to file
for f in ./*.gz; do md5sum $f >> md5sum.txt; done
# copyeverything in this folder excluding .log files to s3
aws s3 cp . s3://lab-aaf-ngs-data-archive/RNAseq/20240424_Tet2Rhoa-AGX51_APL/ --recursive --exclude "*.log"
mv md5sum.txt checksums.txt
