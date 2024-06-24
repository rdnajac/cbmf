for f in ./*_L001_R1_001.fastq.gz; do
	# mv f "${f/_L001_R1_001.fastq.gz/_R1.fastq.gz}"
	# use sed to get hte new dilename
	mv "$f" "$(echo "$f" | sed 's/_L001_R1_001.fastq.gz/_R1.fastq.gz/')"
done

oneliner
for f in ./*_L001_R1_001.fastq.gz; do mv "$f" "$(echo "$f" | sed 's/_L001_R1_001.fastq.gz/_R1.fastq.gz/')"; done
for f in ./*_L001_R2_001.fastq.gz; do mv "$f" "$(echo "$f" | sed 's/_L001_R2_001.fastq.gz/_R2.fastq.gz/')"; done

for f1 in ./*_L001_R1_001.fastq.gz; do cat "$f1" >"$(echo "$f1" | sed 's/_L001_R1_001.fastq.gz/_R1.fastq.gz/')" && rm -f "$f1"; done
for f2 in ./*_L002_R1_001.fastq.gz; do cat "$f2" >"$(echo "$f2" | sed 's/_L002_R1_001.fastq.gz/_R1.fastq.gz/')" && rm -f "$f2"; done
for f2 in ./*_L002_R2_001.fastq.gz; do cat "$f2" >>"$(echo "$f2" | sed 's/_L002_R2_001.fastq.gz/_R2.fastq.gz/')" && rm -f "$f2"; done
for f3 in ./*_L003_R1_001.fastq.gz; do cat "$f3" >>"$(echo "$f3" | sed 's/_L003_R1_001.fastq.gz/_R1.fastq.gz/')" && rm -f "$f3"; done
for f3 in ./*_L003_R2_001.fastq.gz; do cat "$f3" >>"$(echo "$f3" | sed 's/_L003_R2_001.fastq.gz/_R2.fastq.gz/')" && rm -f "$f3"; done
for f4 in ./*_L004_R1_001.fastq.gz; do cat "$f4" >>"$(echo "$f4" | sed 's/_L004_R1_001.fastq.gz/_R1.fastq.gz/')" && rm -f "$f4"; done
for f4 in ./*_L004_R2_001.fastq.gz; do cat "$f4" >>"$(echo "$f4" | sed 's/_L004_R2_001.fastq.gz/_R2.fastq.gz/')" && rm -f "$f4"; done

# write protext all .fastq.gz files
chmod -w ./*.fastq.gz
