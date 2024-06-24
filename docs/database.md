# Database

## sql

- database stored on lab website server palomerolab.org
- make sure cPanel Remote MySQL has this Access Host: athens-public-ph1-ctrl%.cpmc.columbia.edu
  - # note the wildcard '%' (the IP address of the server is dynamic)

```sh
#mysql -u root -p123456 -e "create database test;"
#mysql --version
mysql -u labaf_ryan -p -h palomerolab.org --ssl-mode=REQUIRED
```

## [TOML (Tom's Obvious, Minimal Language)](https://toml.io/en/)

Use TOML to store database connection information.

```toml
[author]
name = "Teresa Palomero"

[[publications]]
title = "Genomic characterization of cutaneous T-cell lymphoma and Sezary syndrome"
journal = "Nature Genetics"
year = 2023
url = "https://doi.org/10.1038/s41588-023-01001-9"

[[publications]]
title = "Therapeutic targeting of NOTCH1 signaling in Peripheral T-cell Lymphomas"
journal = "Nature Reviews Cancer"
year = 2022
url = "https://doi.org/10.1038/s41586-022-05326-3"

[[publications]]
title = "The genomic landscape of angioimmunoblastic T-cell lymphoma"
journal = "Nature Communications"
year = 2021
url = "https://doi.org/10.1038/s41467-021-21410-3"

[[publications]]
title = "Genomic and epigenomic insights into the pathogenesis of cutaneous T-cell lymphoma"
journal = "Blood"
year = 2020
url = "https://doi.org/10.1182/blood.2020004834"

[[publications]]
title = "Role of Vav1 genomic alterations in T-Cell differentiation and transformation in Peripheral T-Cell Lymphoma"
journal = "Proceedings of the National Academy of Sciences"
year = 2019
url = "https://doi.org/10.1073/pnas.1901733116"

[[publications]]
title = "Mechanisms of therapeutic response to Tipifarnib in a mouse model of Angioimmunoblastic T-Cell Lymphoma"
journal = "Cancer Cell"
year = 2018
url = "https://doi.org/10.1016/j.ccell.2018.01.001"

[[publications]]
title = "Activating mutations and translocations in the guanine exchange factor VAV1 in peripheral T-cell lymphomas"
journal = "Proceedings of the National Academy of Sciences"
year = 2017
url = "https://www.pnas.org/doi/pdf/10.1073/pnas.1608839114"

[[publications]]
title = "Targeted cellular immunotherapy for T cell malignancies"
journal = "Nature Medicine"
year = 2017
url = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5848213/"

[[publications]]
title = "Targeted cellular immunotherapy for T cell malignancies"
journal = "Nature medicine"
year = 2017
url = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5848213/"

[[publications]]
title = "Activating mutations and translocations in the guanine exchange factor VAV1 in peripheral T-cell lymphomas"
journal = "Proceedings of the National Academy of Sciences"
year = 2017
url = "https://www.pnas.org/doi/pdf/10.1073/pnas.1608839114"

[[publications]]
title = "Mutational landscape, clonal evolution patterns, and role of RAS mutations in relapsed acute lymphoblastic leukemia"
journal = "Proceedings of the National Academy of Sciences"
year = 2016
url = "https://www.pnas.org/doi/pdf/10.1073/pnas.1608420113"

[[publications]]
title = "The curious origins of angioimmunoblastic T-cell lymphoma"
journal = "Current opinion in hematology"
year = 2016
url = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5499521/"

[[publications]]
title = "The mutational landscape of cutaneous T cell lymphoma and Sézary syndrome"
journal = "Nature genetics"
year = 2015
url = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4878831/"

[[publications]]
title = "A NOTCH1-driven MYC enhancer promotes T cell development, transformation and acute lymphoblastic leukemia"
journal = "Nature medicine"
year = 2014
url = "https://core.ac.uk/download/pdf/327342554.pdf"

[[publications]]
title = "N-Me, a long range oncogenic enhancer in T-cell acute lymphoblastic leukemia"
journal = "Nature medicine"
year = 2014
url = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4192073/"

[[publications]]
title = "Recurrent mutations in epigenetic regulators, RHOA and FYN kinase in peripheral T cell lymphomas"
journal = "Nature genetics"
year = 2014
url = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3963408/"

[[publications]]
title = "Activating mutations in the NT5C2 nucleotidase gene drive chemotherapy resistance in relapsed ALL"
journal = "Nature medicine"
year = 2013
url = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3594483/"

[[publications]]
title = "Therapeutic effect of γ-secretase inhibition in Kras G12V-driven non-small cell lung carcinoma by derepression of DUSP1 and inhibition of ERK"
journal = "Cancer cell"
year = 2012
url = "https://www.cell.com/cancer-cell/fulltext/S1535-6108(12)00263-2?idioma=galego"

[[publications]]
title = "Reverse engineering of TLX oncogenic transcriptional networks identifies RUNX1 as tumor suppressor in T-ALL"
journal = "Nature medicine"
year = 2012
url = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3298036/"

[[publications]]
title = "Whole exome sequencing to identify a novel gene (caveolin-1) associated with human pulmonary arterial hypertension"
journal = "Circulation: Cardiovascular Genetics"
year = 2012
url = "https://www.ahajournals.org/doi/full/10.1161/circgenetics.111.961888"

[[publications]]
title = "The TLX1 oncogene drives aneuploidy in T cell transformation"
journal = "Nature medicine"
year = 2010
url = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2974790/"

[[publications]]
title = "ParMap, an algorithm for the identification of small genomic insertions and deletions in nextgen sequencing data"
journal = "BMC research notes"
year = 2010
url = "https://link.springer.com/article/10.1186/1756-0500-3-147"

[[publications]]
title = "PHF6 mutations in T-cell acute lymphoblastic leukemia"
journal = "Nature genetics"
year = 2010
url = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2847364/"

[[publications]]
title = "Genome-wide RNA-mediated interference screen identifies miR-19 targets in Notch-induced T-cell acute lymphoblastic leukaemia"
journal = "Nature cell biology"
year = 2010
url = "http://fulltext.calis.edu.cn/nature/ncb/12/4/ncb2037.pdf"

[[publications]]
title = "Regulation of hematopoietic stem cell differentiation by a single ubiquitin ligase-substrate complex"
journal = "Nature immunology"
year = 2010
url = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2825759/"

[[publications]]
title = "ParMap, an algorithm for the identification of complex genomic variations in nextgen sequencing data"
journal = "Nature Publishing Group"
year = 2010
url = "https://www.nature.com/articles/npre.2010.4145.1.pdf"

[[publications]]
title = "Therapeutic targeting of NOTCH1 signaling in T-cell acute lymphoblastic leukemia"
journal = "Clinical Lymphoma and Myeloma"
year = 2009
url = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2814179/"

[[publications]]
title = "WT1 mutations in T-all"
journal = "Blood"
year = 2009
url = "https://www.sciencedirect.com/science/article/pii/S0006497120421767"

[[publications]]
title = "Genomic tools for dissecting oncogenic transcriptional networks in human leukemia"
journal = "Leukemia"
year = 2009
url = "https://www.researchgate.net/profile/Teresa-Palomero-2/publication/23807411_Genomic_tools_for_dissecting_oncogenic_transcriptional_networks_in_human_leukemia/links/565cb52d08ae1ef92981f254/Genomic-tools-for-dissecting-oncogenic-transcriptional-networks-in-human-leukemia.pdf"

[[publications]]
title = "ChIP-on-chip significance analysis reveals large-scale binding and regulation by human transcription factor oncogenes"
journal = "Proceedings of the National Academy of Sciences"
year = 2009
url = "https://www.pnas.org/doi/pdf/10.1073/pnas.0806445106"

[[publications]]
title = "γ-secretase inhibitors reverse glucocorticoid resistance in T cell acute lymphoblastic leukemia"
journal = "Nature medicine"
year = 2009
url = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2692090/"

[[publications]]
title = "Oncogenic NOTCH1 control of MYC and PI3K: challenges and opportunities for anti-NOTCH1 therapy in T-cell acute lymphoblastic leukemias and lymphomas"
journal = "Clinical cancer research"
year = 2008
url = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2577004/"

[[publications]]
title = "NOTCH1 extracellular juxtamembrane expansion mutations in T-ALL"
journal = "Blood"
year = 2008
url = "https://www.sciencedirect.com/science/article/pii/S0006497120604190"

[[publications]]
title = "The role of the PTEN/AKT Pathway in NOTCH1-induced leukemia"
journal = "Cell cycle"
year = 2008
url = "https://www.tandfonline.com/doi/pdf/10.4161/cc.7.8.5753"
```

### `import toml`

Use python to import the TOML file.

```python
import toml

# Read the TOML file
with open('palomero.toml', 'r') as file:
    data = toml.load(file)

# Extract author name
author_name = data.get('author', {}).get('name', 'Unknown Author')

# Extract publications
publications = data.get('publications', [])

# Generate Markdown content
markdown_content = f"# Publications by {author_name}\n\n"
for publication in publications:
    title = publication.get('title', 'No Title')
    journal = publication.get('journal', 'No Journal')
    year = publication.get('year', 'No Year')
    url = publication.get('url', '#')

    markdown_content += f"## {title}\n"
    markdown_content += f"*Journal:* {journal}\n"
    markdown_content += f"*Year:* {year}\n"
    markdown_content += f"*Link:* [{url}]({url})\n\n"

# Write the Markdown content to a file
with open('publications.md', 'w') as file:
    file.write(markdown_content)

print("Markdown file generated: publications.md")
```
