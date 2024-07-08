#!/bin/bash
#
## This script generates phenodata.csv files for ballgown analyses

cat << EOF > "$HOME"/cbmf/data/agx/ballgown_input/phenodata.csv
ids,group,replicate
Tet2RhoaG17V_24hAGX51_1,treated,1
Tet2RhoaG17V_24hAGX51_2,treated,2
Tet2RhoaG17V_24hAGX51_3,treated,3
Tet2RhoaG17V_24hDMSO_1,untreated,1
Tet2RhoaG17V_24hDMSO_2,untreated,2
Tet2RhoaG17V_24hDMSO_3,untreated,3
EOF

cat << EOF > "$HOME"/cbmf/data/clone/ballgown_input/phenodata.csv
ids,group,replicate
Tet2Rhoa-sgID2-3-cloneB-induced-1,induced,1
Tet2Rhoa-sgID2-3-cloneB-induced-2,induced,2
Tet2Rhoa-sgID2-3-cloneB-induced-3,induced,3
Tet2Rhoa-sgID2-3-cloneB-noninduced-1,noninduced,1
Tet2Rhoa-sgID2-3-cloneB-noninduced-2,noninduced,2
Tet2Rhoa-sgID2-3-cloneB-noninduced-3,noninduced,3
EOF

cat << EOF > "$HOME"/cbmf/data/ra/ballgown_input/phenodata.csv
ids,treatment,replicate
DMSO1,DMSO,1
DMSO2,DMSO,2
DMSO3,DMSO,3
Fingolimod1,Fingolimod,1
Fingolimod2,Fingolimod,2
Fingolimod3,Fingolimod,3
Ozanimod1,Ozanimod,1
Ozanimod2,Ozanimod,2
Ozanimod3,Ozanimod,3
Ponesimod1,Ponesimod,1
Ponesimod2,Ponesimod,2
Ponesimod3,Ponesimod,3
EOF

"$()" and $(()): {{{1
# " $(..) is not supported by sh (Bourne shell).  However, apparently
# " some systems (HP?) have as their /bin/sh a (link to) Korn shell
# " (ie. Posix compliant shell).  /bin/ksh should work for those
# " systems too, however, so the following syntax will flag $(..) as
# " an Error under /bin/sh.  By consensus of vimdev'ers!

# example using $(..) instead of `..` (backticks)
echo "Today is $( "$(date)" )"
