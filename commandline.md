# Markdown file with shell commands (aws cli)
### `arn:aws:s3:::lab-aaf-ngs-data-archive`

``` bash
aws s3 ls lab-aaf-ngs-data-archive/ChIPseq/
#ssh ec2-35-175-217-48.compute-1.amazonaws.com
#vim scp://ubuntu@ec2-35-175-217-48.compute-1.amazonaws.com//home/ubuntu/README.md
```

Folders of interest:
* 20230216_ChIP-seq_HEL_H3K27ac_ruxolitinib_SGC-CBP30_CB/
* 20230308_ChIP-seq_HEL_H3K27ac_ruxolitinib_SGC-CBP30_rep2redo_CB/
* 20230308_ChIP-seq_HEL_STAT5_ruxolitinib_SGC-CBP30_CB/
* 20230315_mPTCL_CellLines_HistoneMarks/

Full path:
lab-aaf-ngs-data-archive/ChIPseq/20230216_ChIP-seq_HEL_H3K27ac_ruxolitinib_SGC-CBP30_CB/00_fastq/
lab-aaf-ngs-data-archive/ChIPseq/20230308_ChIP-seq_HEL_H3K27ac_ruxolitinib_SGC-CBP30_rep2redo_CB/00_fastq/
lab-aaf-ngs-data-archive/ChIPseq/20230308_ChIP-seq_HEL_STAT5_ruxolitinib_SGC-CBP30_CB/00_fastq

``bash
folder1=lab-aaf-ngs-data-archive/ChIPseq/20230216_ChIP-seq_HEL_H3K27ac_ruxolitinib_SGC-CBP30_CB/00_fastq/
folder2=lab-aaf-ngs-data-archive/ChIPseq/20230308_ChIP-seq_HEL_H3K27ac_ruxolitinib_SGC-CBP30_rep2redo_CB/00_fastq/
folder3=lab-aaf-ngs-data-archive/ChIPseq/20230308_ChIP-seq_HEL_STAT5_ruxolitinib_SGC-CBP30_CB/00_fastq/
folders=($folder1 $folder2 $folder3)

# make sure we got it
for folder in ${folders[@]}; do echo $folder; done

# output looks like this
# 2023-02-16 20:17:53  201831270 RUX-2-input_S6_L004_R2_001.fastq.gz
# iwant get a list of every unique sample upti the S##_
aws s3 ls $folder1 | awk '{print $4}' | cut -d'_' -f1-2 | sort | uniq
aws s3 ls $folder2 | awk '{print $4}' | cut -d'_' -f1-2 | sort | uniq
aws s3 ls $folder3 | awk '{print $4}' | cut -d'_' -f1-2 | sort | uniq | pbcopy

for folder in ${folders[@]}; do aws s3 ls $folder | awk '{print $4}' | cut -d'_' -f1-2 | sort | uniq; done

for folder in ${folders[@]}; do aws s3 ls $folder > bigls.txt; done
aws s3 ls $folder1 | pbcopy
 w

``` bash
```

``` bash

# get the list of files in the folder
aws s3 ls $folder3
# just the filenames
aws s3 ls $folder1 | awk '{print $4}' >> samples.txt
aws s3 ls $folder2 | awk '{print $4}' > samples.txt
aws s3 ls $folder3 | awk '{print $4}' > samples.txt
```

``` bash
```

``` bash
# get a list of every sample and save just the name of the sample to a file
for folder in ${folders[@]}; do echo "saving samples from $folder" > samples.txt; aws s3 ls $folder | awk '{print $4}' > samples.txt; done
```


## Basespace
``` bash
# https://github.com/basespace/homebrew-basespace
brew tap basespace/basespace
brew install bs-cli
bs authenticate
```
### `bs` commands

|Command|Description|
|---|---|
`accession`|Register a new item
`add`|Add lane or workflow settings
`archive`|Archive entities
`authenticate`|Make an authentication request to BaseSpace (aliases: auth)
`await`|Await completion or status change
`children`|Retrieve child entities (aliases: child)
`clear`|Clear lane or workflow settings
`content`|Show all files in an entity (aliases: contents, dir)
`copy`|Copy an entity to a new location (aliases: duplicate)
`create`|Create an entity
`delete`|Move entities to your trash (aliases: rm, move-to-trash)
`download`|Download files from an entity
`empty`|Empty a store
`export`|Export lane or workflow settings
`get`|Get detailed information about an entity (aliases: info, inspect)
`header`|List headers for an entity (aliases: headers)
`history`|Retrieve account activity records for a user or workgroup
`import`|Import lane or workflow settings
`kill`|Abort entities
`launch`|Execute a task
`link`|Get a direct link for an entity
`list`|List and filter entities (aliases: filter, list-all)
`load`|Load into your environment
`logs`|Retrieve log files (aliases: log)
`rename`|Rename entities
`restore`|Restore items
`revoke`|Invalidate a resource (aliases: expire)
`seqstats`|Sequencing stats
`set`|Set properties in an entity
`translate`|Translate v1 <-> v2 entity IDs
`unarchive`|Restore entities from archival storage
`unlock`|Unlock a locked entity
`update`|Update entities
`upload`|Upload files to an entity
`whoami`|Get information about selected user/configuration
`yield`|Return yield information for an entity (aliases: yields)

```
bs list project
```

|Name|Id|
|---|---|
| HEL_H3K27ac_ChIPSeq1 |381715669|
| Untitled from 230306_NS500289_1231_AHW5Y3BGXN |383461079|
| CBP_P300 |388989668| 

bs list biosample > biosamples.txt


see biosamples from a run

``` bash
bs list biosample --project-id 381715669
```

To list samples:
``` bash
bs list biosample --project-id 381715669
bs list biosample --project-id 383461079
bs list biosample --project-id 388989668
```

get a comma separated list of just the biosameple name and ids (two columns)
``` bash
bs list biosample --project-id 381715669 
bs list biosample --project-id 383461079 | awk '{print  $2 "," $4}' | pbcopy
bs list biosample --project-id 388989668 
```

,
BioSampleName,Id
,
CBP_DMSO_1,625657188
CBP_RUX_1,625657189
CBP_CBP_1,625657190
CBP_Combo_1,625657191
CBP_DMSO_2,625657192
CBP_RUX_2,625657193
CBP_CBP_2,625657194
CBP_Combo_2,625657195
P300_DMSO_1,625657196
P300_RUX_1,625657197
P300_CBP_1,625657198
P300_Combo_1,625657199
P300_DMSO_2,625657200
P300_RUX_2,625657201
P300_CBP_2,625657202
,

|Name|Id|ExperimentName|
|---|---|---|
| 230215_NB551203_0628_AHC32JBGXN |253069839|HEL_H3K27ac_ChIPSeq1|
| 230306_NS500289_1231_AHW5Y3BGXN |254559306|HEL_STAT5_H3K27ac_2|
| 230519_NS500289_1262_AHHLJNBGXT |258955746|CBP_p300|

