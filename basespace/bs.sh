#!/bin/bash
# basespace CLI commands

run_name=("230215_NB551203_0628_AHC32JBGXN" "230306_NS500289_1231_AHW5Y3BGXN" "230519_NS500289_1262_AHHLJNBGXT")
run_id=("253069839" "254559306" "258955746")
exp_name=("HEL_H3K27ac_ChIPSeq1" "HEL_STAT5_H3K27ac_2" "CBP_p300")

sample=( # {{{
1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
DMSO_1_STAT5
RUX_1_STAT5
CBP_1_STAT5
Combo_1_STAT5
DMSO_2_STAT5
RUX_2_STAT5
CBP_2_STAT5
Combo_2_STAT5
DMSO_2_input2
RUX_2_input2
CBP_2_input2
Combo_2_input2
DMSO_2_H3K27ac2
RUX_2_H3K27ac2
CBP_2_H3K27ac2
Combo_2_H3K27ac2
CBP_DMSO_1
CBP_RUX_1
CBP_CBP_1
CBP_Combo_1
CBP_DMSO_2
CBP_RUX_2
CBP_CBP_2
CBP_Combo_2
P300_DMSO_1
P300_RUX_1
P300_CBP_1
P300_Combo_1
P300_DMSO_2
P300_RUX_2
P300_CBP_2
) # }}}
ids=( # {{{
610706528
610706529
610706530
610706531
610706532
610706533
610706534
610706535
610706536
610706537
610706538
610706539
610706540
610706541
610706542
610706543
614276798
614276799
614276800
614276801
614276802
614276803
614276804
614276805
614276806
614276807
614276808
614276809
614276810
614276811
614276812
614276813
625657188
625657189
625657190
625657191
625657192
625657193
625657194
625657195
625657196
625657197
625657198
625657199
625657200
625657201
625657202
) # }}}
#bs get run --id "${run_id[0]}"

# get the biosample in the run
bs yield biosample $id[0]
# entites {{{
# entities: application - appresult - appsession - biosample - config
# dataset - file - lane - manifest - project - run  - trash - workflow
# }}}

exit 0
# install aws cli on ubuntu
sudo apt install awscli