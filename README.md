# 16S Nanopore Workflow for the NEMADAPT (Nematode ADAPTation in response to Thermal stress) project

### Author: Rodgee Mae Guden  
**Marine Biology Research Group, Ghent University**  
**Contact:** rodgeemae.guden@ugent.be  

---

## Overview

This repository contains a modular, robust Bash workflow for processing **16S Nanopore sequencing data** within the framework of the [NEMADAPT project](https://www.marinebiology.ugent.be/). It automates the key steps for high-throughput microbial community analysis, from raw FASTQ files to taxonomic profiles and OTU tables.

The workflow includes:

✔️ Demultiplexing  
✔️ Quality control (FastQC, MultiQC)  
✔️ Adapter trimming (Porechop)  
✔️ Chimera detection & removal (yacrd)  
✔️ Length filtering (fastp)  
✔️ Post-filter QC (FastQC, MultiQC)  
✔️ Taxonomic classification (Kraken2)  
✔️ Abundance estimation (Bracken)  
✔️ OTU table generation  
✔️ Krona visualization  

---

## Repository Structure

```
📁 16S-Nanopore-Workflow/
├── 16S_workflow_nemadapt.sh    # Main Bash workflow script
├── README.md                  # Project documentation
├── LICENSE                    # License file (optional)
└── Example Output Folders Created by the Script:
    ├── concatenated/          # Concatenated FASTQ files per barcode
    ├── renamed/               # Renamed FASTQ files per sample
    ├── nanostat/              # Read statistics reports
    ├── qc/                    # Initial QC reports
    ├── trimmed/               # Adapter-trimmed reads
    ├── filtered/              # Length-filtered, non-chimeric reads
    ├── qc_filtered/           # Post-filter QC reports
    ├── taxonomy/              # Kraken2 classification outputs
    ├── bracken/               # Abundance estimation results
    ├── visualization/         # Krona charts
    └── otu_tables/            # OTU tables at multiple taxonomic levels
```

---

## Dependencies

The following tools must be installed and accessible in your `$PATH`:

| Tool         | Purpose                             | Installation link                           |
|--------------|-------------------------------------|----------------------------------------------|
| `fastqc`     | Read quality control               | [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) |
| `multiqc`    | Aggregate QC reports               | [MultiQC](https://multiqc.info/)            |
| `porechop`   | Adapter trimming for Nanopore      | [Porechop](https://github.com/rrwick/Porechop) |
| `minimap2`   | Overlap detection (chimera removal)| [Minimap2](https://github.com/lh3/minimap2) |
| `yacrd`      | Chimera detection & removal        | [yacrd](https://github.com/natir/yacrd)     |
| `fastp`      | Length filtering, quality trimming | [fastp](https://github.com/OpenGene/fastp)  |
| `kraken2`    | Taxonomic classification           | [Kraken2](https://github.com/DerrickWood/kraken2) |
| `bracken`    | Abundance estimation refinement    | [Bracken](https://github.com/jenniferlu717/Bracken) |
| `krona`      | Taxonomy visualization             | [Krona](https://github.com/marbl/Krona)     |
| `NanoStat`   | Read statistics summary            | [NanoStat](https://github.com/wdecoster/nanostat) |

You must also have a **Kraken2 database**, such as SILVA, properly set up.

---

## Usage

### 1. Make the Script Executable

```bash
chmod +x 16S_workflow_nemadapt.sh
```

### 2. Run the Full Workflow (All Steps in Order)

```bash
bash 16S_workflow_nemadapt.sh
```

### 3. Run Specific Steps Only

You can run selected steps by specifying them with the `-s` flag, separated by commas:

```bash
bash 16S_workflow_nemadapt.sh -s qc,trim_adapter,filter,classify
```

### 4. Available Workflow Steps

| Step                       | Description                                                        |
|----------------------------|--------------------------------------------------------------------|
| `concat_fastq`             | Concatenate FASTQ files per barcode                                |
| `rename_fastq_FC1/FC2/FC3` | Rename concatenated FASTQ files per sample (different sample sets) |
| `nanostat`                 | Generate read statistics with NanoStat                            |
| `qc`                       | Initial quality control (FastQC, MultiQC)                         |
| `trim_adapter`             | Remove adapters with Porechop                                     |
| `trim_chimeric`            | Remove chimeric reads (yacrd)                                     |
| `filter`                   | Filter reads by length and quality                                |
| `postqc`                   | Post-filter quality control                                       |
| `classify`                 | Classify reads with Kraken2                                       |
| `krona`                    | Create Krona visualizations                                       |
| `bracken_setup`            | Prepare Bracken database (only once)                              |
| `bracken`                  | Estimate species/genus abundances                                 |
| `extract_otu`              | Generate OTU tables for all taxonomic levels                      |

---

## Output Example

After running the full workflow, key outputs will be located in:

- `output/qc/multiqc_report/` — Initial QC report  
- `output/filtered/` — Clean, filtered reads  
- `output/taxonomy/` — Kraken2 classification reports  
- `output/bracken/` — Refined abundance estimates  
- `output/otu_tables/` — Final OTU tables (species, genus, etc.)  
- `output/visualization/multi-krona.html` — Interactive taxonomic chart  

---

## Notes

- The script is configured for **12 threads** by default. Adjust the `THREADS` parameter at the top of the script as needed.
- Input/output directories and database paths are defined in the script header. Modify them to match your environment.
- Ensure the Kraken2 database is properly built before running the classification step.

---

## Citation & License

© 2025 Rodgee Mae Guden  
**Marine Biology Research Group, Ghent University**  
All rights reserved.

This workflow and script are the intellectual property of the author. Reuse or redistribution requires prior written permission, except for brief academic quotations with proper attribution.

**For permissions, contact:** rodgeemae.guden@ugent.be  

---

## Acknowledgments

This workflow was developed within the **NEMADAPT** project, which focuses on developing molecular tools to study nematode adaptation in response to climate change using Nanopore sequencing technology.

---

## Questions & Issues

If you encounter problems or have suggestions, feel free to open an [Issue](https://github.com/yourusername/16S-Nanopore-Workflow/issues) or contact the author directly.
