# host-rm-bowtie2 (human read removal for paired-end FASTQ)

A lightweight, reproducible workflow to **remove human (host) reads** from paired-end Illumina FASTQs using **Bowtie2**.

✅ Works on macOS (Homebrew) and Conda/Mamba  
✅ Batch mode from a `fastq_list.txt` (2 columns: R1 R2)  
✅ Outputs host-removed paired FASTQs (`*_R1.clean.fastq.gz`, `*_R2.clean.fastq.gz`)  
✅ Logs Bowtie2 alignment summary per sample

---

## Repository name suggestion

Recommended repo name (short + descriptive):

**`host-rm-bowtie2-human`**

Alternative options:
- `human-host-read-removal-bowtie2`
- `host-decontam-bowtie2`
- `bowtie2-host-removal-batch`

---

## Repo layout (recommended)

```
host-rm-bowtie2-human/
├── README.md
├── LICENSE
├── .gitignore
├── .gitattributes
├── scripts/
│   └── host_remove_bt2_batch.sh
├── examples/
│   ├── fastq_list.example.txt
│   └── run_example.sh
└── reference/
    ├── README.md
    ├── human.fa.gz            # (optional) large file via Git LFS
    └── bt2_index/             # (optional) build outputs (NOT recommended to commit)
```

### Notes on storing the human reference
- **Best practice:** do *not* commit the full human genome FASTA or Bowtie2 index to GitHub.
- Instead, include **instructions** to download from an official source (e.g., GRCh38) and build the index locally.
- If you still choose to store it, use **Git LFS** and expect GitHub storage/bandwidth limits.

---

## Quick start

### 1) Install dependencies

**Option A (recommended on macOS): Homebrew**
```bash
brew install bowtie2 samtools
which bowtie2 && bowtie2 --version
which bowtie2-build && bowtie2-build --version
which samtools && samtools --version
```

**Option B: Conda/Mamba**
```bash
conda create -n hostrm -y -c bioconda -c conda-forge bowtie2 samtools
conda activate hostrm
which bowtie2 && bowtie2 --version
which samtools && samtools --version
```

---

### 2) Build the Bowtie2 index from `human.fa`

Pick a location for the index (example matches your paths):

```bash
mkdir -p /Volumes/MBOOWA/Uganda_Jan_2026/human_bt2_index
bowtie2-build \
  /Volumes/MBOOWA/Uganda_Jan_2026/human.fa \
  /Volumes/MBOOWA/Uganda_Jan_2026/human_bt2_index/human
```

This creates files like:
- `human.1.bt2` … `human.4.bt2`
- `human.rev.1.bt2`, `human.rev.2.bt2`
(or `.bt2l` for large genomes)

---

### 3) Prepare `fastq_list.txt`

Two columns per line (tab or space), **R1 then R2**:

```text
~/SampleA_1.fastq.gz ~/SampleA_2.fastq.gz
~/SampleB_1.fastq.gz ~/SampleB_2.fastq.gz
~/SampleC_1.fastq.gz ~/SampleC_2.fastq.gz
```

---

### 4) Run host removal (batch)

```bash
bash scripts/host_remove_bt2_batch.sh \
  -x /Volumes/MBOOWA/Uganda_Jan_2026/human_bt2_index/human \
  -l /Users/gmboowa/Desktop/Staphylococcus-TZ/fastq_list.txt \
  -o /Users/gmboowa/Desktop/Staphylococcus-TZ/host_removed_reads \
  -t 8
```

Example output:
```
▶ Processing: A55934
processed 5663646 reads
✓ Wrote: .../A55934_R1.clean.fastq.gz
✓ Wrote: .../A55934_R2.clean.fastq.gz
Done.
```

---

## The script

Place this at: `scripts/host_remove_bt2_batch.sh` and make executable:

```bash
chmod +x scripts/host_remove_bt2_batch.sh
```

Script (current working version):

```bash
#!/usr/bin/env bash
set -euo pipefail

IDX=""
LIST=""
OUT=""
THREADS=4

usage() {
  echo "Usage: $0 -x <bt2_index_prefix> -l <fastq_list.txt> -o <outdir> [-t <threads>]"
  echo "  fastq_list.txt: two columns per line: R1 <tab or space> R2"
  exit 1
}

while getopts "x:l:o:t:" opt; do
  case $opt in
    x) IDX="$OPTARG" ;;
    l) LIST="$OPTARG" ;;
    o) OUT="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    *) usage ;;
  esac
done

[[ -n "$IDX" && -n "$LIST" && -n "$OUT" ]] || usage
[[ -r "$LIST" ]] || { echo "ERROR: Cannot read list file: $LIST"; exit 1; }

if [[ ! -f "${IDX}.1.bt2" && ! -f "${IDX}.1.bt2l" ]]; then
  echo "ERROR: Bowtie2 index not found for prefix: $IDX"
  exit 1
fi

mkdir -p "$OUT"

echo "Index prefix : $IDX"
echo "List file    : $LIST"
echo "Output dir   : $OUT"
echo "Threads      : $THREADS"
echo

while IFS=$' \t' read -r R1 R2 REST || [[ -n "${R1:-}" ]]; do
  R1="${R1//$'\r'/}"
  R2="${R2//$'\r'/}"
  [[ -z "${R1:-}" ]] && continue
  [[ "${R1:0:1}" == "#" ]] && continue
  [[ -z "${R2:-}" ]] && { echo "WARNING: missing R2 for $R1"; continue; }

  [[ -f "$R1" ]] || { echo "WARNING: R1 not found: $R1"; continue; }
  [[ -f "$R2" ]] || { echo "WARNING: R2 not found: $R2"; continue; }

  base=$(basename "$R1")
  SAMPLE="${base%.fastq.gz}"
  SAMPLE="${SAMPLE%.fq.gz}"
  SAMPLE="${SAMPLE%_R1_001}"
  SAMPLE="${SAMPLE%_R1}"
  SAMPLE="${SAMPLE%_1}"

  echo "▶ Processing: $SAMPLE"

  bowtie2 -x "$IDX" \
    -1 "$R1" -2 "$R2" \
    -p "$THREADS" \
    --very-sensitive \
    --un-conc-gz "$OUT/${SAMPLE}.hostRemoved.%.fastq.gz" \
    -S /dev/null \
    2> "$OUT/${SAMPLE}.bowtie2.log"

  mv "$OUT/${SAMPLE}.hostRemoved.1.fastq.gz" "$OUT/${SAMPLE}_R1.clean.fastq.gz"
  mv "$OUT/${SAMPLE}.hostRemoved.2.fastq.gz" "$OUT/${SAMPLE}_R2.clean.fastq.gz"

  echo "  ✓ Wrote: $OUT/${SAMPLE}_R1.clean.fastq.gz"
  echo "  ✓ Wrote: $OUT/${SAMPLE}_R2.clean.fastq.gz"
  echo
done < "$LIST"

echo "Done."
```

---

## Git LFS setup for `human.fa` (if you insist on storing it)

### 1) Install Git LFS
```bash
git lfs install
```

### 2) Track the FASTA with LFS
From the repo root:
```bash
git lfs track "reference/*.fa"
git lfs track "reference/*.fa.gz"
```

### 3) Update `.gitattributes`
This will be created/updated automatically by `git lfs track`, but it should include lines like:

```gitattributes
reference/*.fa filter=lfs diff=lfs merge=lfs -text
reference/*.fa.gz filter=lfs diff=lfs merge=lfs -text
```

### 4) Add and commit
```bash
git add .gitattributes reference/human.fa.gz
git commit -m "Track human reference with Git LFS"
```

---

## Output files

For each sample, you get:
- `{SAMPLE}_R1.clean.fastq.gz`
- `{SAMPLE}_R2.clean.fastq.gz`
- `{SAMPLE}.bowtie2.log`

Validate output file integrity with:
```bash
gzip -t SampleA_R1.clean.fastq.gz
gzip -t SampleA_R2.clean.fastq.gz
```

---

## Citation

If you use this workflow in publications:
- Bowtie2: Langmead & Salzberg, Nat Methods (2012)
- SAMtools: Li et al., Bioinformatics (2009)

---

## License

MIT License

---

## Contact
Gerald Mboowa
