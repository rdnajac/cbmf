# StringTie

## Installation

```bash
git clone https://github.com/gpertea/stringtie
cd stringtie
make release -j4
```

### Add to PATH

``` bash
export PATH=$PATH:/path/to/stringtie

# or copy to /usr/local/bin
sudo cp stringtie /usr/local/bin
```


```

usage: stringtie [-o <output.gtf>] [other_options] <read_alignments.bam>
# using nproc threadds with -p ootions
stringtie -p $(nproc) -o output.gtf input.bam
# do it on these bbam files:

ubuntu@ip-172-31-32-180:~/rnaseq/ra/aligned$ ls
DMSO1.bam  Fingolimod1.bam  out            Ozanimod3.bam   Ponesimod3.bam
DMSO2.bam  Fingolimod2.bam  Ozanimod1.bam  Ponesimod1.bam  test
DMSO3.bam  Fingolimod3.bam  Ozanimod2.bam  Ponesimod2.bam
ubuntu@ip-172-31-32-180:~/rnaseq/ra/aligned$


for i in *.bam; do stringtie -p $(nproc) -o ${i}.gtf $i & ; done
# fix it so it runs in parallel
# zip all .gtf files into ra_rnaseq.zip
zip ra_rnaseq.zip *.gtf
```
