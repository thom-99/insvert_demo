# Task Plan: `benchmark_insert.sh` - Performance Benchmark for `inSVert insert`

## Goal
Create a reusable performance benchmark script `scripts/benchmark_insert.sh` dedicated specifically to measuring wall-clock time, peak RSS memory, and CPU usage for `inSVert insert`.

---

## 1. Features & Requirements
1. **Inputs:**
   - `<reference.fasta>` (Path to reference FASTA)
   - `<variants.vcf>` (Path to VCF file containing structural variants)
   - `[ploidy]` (Optional ploidy number, default: 2; if a `config.yaml` path is provided as ploidy argument or config, auto-extract ploidy)
   - `[output_dir]` (Optional output directory for FASTA and timing logs, default: `benchmark_insert_results`)
2. **Metrics:**
   - Wall-clock time (seconds / formatted)
   - Peak RSS Memory (MB and KB)
   - User CPU time (seconds)
   - System CPU time (seconds)
   - Average CPU usage %
3. **Execution Command:**
   - Uses `/usr/bin/time -v conda run -n demo inSVert insert ...`
4. **Output Report:**
   - Formatted console report outputting all stats for `inSVert insert`.

---

## 2. File Specifications (`scripts/benchmark_insert.sh`)
- `scripts/benchmark_insert.sh`: Executable bash script.
- Handles argument validation, directory creation, ploidy extraction (if a YAML file is supplied), process execution via `/usr/bin/time -v`, and formatted log reporting.

---

## 3. Implementation Plan
- [ ] **Step 1:** Create `scripts/benchmark_insert.sh`.
- [ ] **Step 2:** Make `scripts/benchmark_insert.sh` executable (`chmod +x scripts/benchmark_insert.sh`).
- [ ] **Step 3:** Test `scripts/benchmark_insert.sh` with a sample reference FASTA and VCF.
