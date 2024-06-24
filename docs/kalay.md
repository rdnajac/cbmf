# ChIPseq for Kalay

This is a markdown file for the ChIPseq analysis of Kalay.

## Data

### BaseSpace

| run date | run name (basespace) | size (GB) | %Q30  | %PF   | instrument | flow cell ID |
| -------- | -------------------- | --------- | ----- | ----- | ---------- | ------------ |
| 20230215 | HEL_H3K27ac_ChIPSeq1 | 27        | 94.49 | 92.21 | NB551203   | HC32JBGXN    |
| 20230306 | HEL_STAT5_H3K27ac_2  | 37        | 92.74 | 89.86 | NS500289   | HW5Y3BGXN    |
| 20230519 | CBP_p300             | 44        | 90.06 | 87.69 | NS500289   | HHLJNBGXT    |
| 20240522 | KB_ChIPseq_4         | 10        | 94.47 | 93.50 | NB551203   | H7TMNAFX7    |
| 20240617 | KB_ChIPseq_5         | 281976710 |

Download the first three runs from BaseSpace.

```plaintext

ubuntu@ip-172-31-32-180:~/kalay$ bs list runs
+---------------------------------+-----------+----------------------+----------+
|              Name               |    Id     |    ExperimentName    |  Status  |
+---------------------------------+-----------+----------------------+----------+
| 230215_NB551203_0628_AHC32JBGXN | 253069839 | HEL_H3K27ac_ChIPSeq1 | Complete |
| 230306_NS500289_1231_AHW5Y3BGXN | 254559306 | HEL_STAT5_H3K27ac_2  | Complete |
| 230519_NS500289_1262_AHHLJNBGXT | 258955746 | CBP_p300             | Complete |
| 240522_NB551203_0777_AH7TMNAFX7 | 281093846 | KB_ChIPseq_4         | Complete |
| 240617_NB551203_0782_AH7NV5AFX7 | 281976710 | Kalay5               | Failed   |
+---------------------------------+-----------+----------------------+----------+
```

Download runs with `bs download run --id "$run_id"`

The first three runs should contain the correct sample sheets for demux.

```bash

```
