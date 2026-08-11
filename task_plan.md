# Implementation Strategy: SNP Simulation and Insertion for inSVert

This document outlines the detailed implementation strategy for adding Single-Nucleotide Polymorphism (SNP) simulation and reference-insertion support to **inSVert**.

---

## Architecture Overview

SNPs differ from Structural Variants (SVs) in key ways:
1. **Coordinates**: A SNP represents a single-base substitution directly at position `POS` (no padding base or anchor base required).
2. **VCF Metadata**: Per VCF 4.2 conventions, SNPs do not use `SVTYPE`, `SVLEN`, or `END` tags. To maintain clean annotation while supporting 3rd-party VCF compatibility, inSVert will write `VT=SNP` in the VCF `INFO` column, while the insertion module will use structural detection (`svtype is None and len(REF) == 1 and len(ALT) == 1`).
3. **Simulation Strategy**: SNPs are simulated **after** all SVs are placed. This guarantees that SVs (which require large contiguous intervals and limited retries) claim their space first, while abundant single-base SNPs can easily route around existing variants.

---

## 1. `VariantObjects.py` — New `SNP` Class

Add a standalone `SNP` class (not inheriting from `StructuralVariant`, as SNPs are not structural variants):

```python
class SNP:
    """
    Represents a Single-Nucleotide Polymorphism (SNP) in VCF 4.2 format.
    
    Unlike StructuralVariant subclasses, SNPs have no padding base, no SVTYPE,
    no SVLEN, and no END field. The REF and ALT are both single nucleotides,
    and POS points directly to the base being substituted.
    """
    
    def __init__(self, chrom: str, pos: int, id: str, genotype: str, ref_base: str, alt_base: str):
        self.chrom = chrom
        self.pos = pos          # 1-indexed VCF position of the substituted base
        self.id = id
        self.genotype = genotype
        self.ref = ref_base     # Single reference nucleotide (A, C, G, or T)
        self.alt = alt_base     # Single alternate nucleotide
        self.qual = "."
        self.filter = "PASS"
    
    def get_end(self) -> int:
        """SNP occupies a single position."""
        return self.pos
    
    def format(self) -> str:
        """Formats the SNP as a VCF 4.2 line with VT=SNP in the INFO column."""
        return (f"{self.chrom}\t{self.pos}\t{self.id}\t{self.ref}\t{self.alt}\t"
                f"{self.qual}\t{self.filter}\tVT=SNP\tGT\t{self.genotype}")
```

---

## 2. `utils_sim.py` — Config Parsing, Ts/Tv Helper & Header Line

### 2.1 Config Parsing (`parse_config`)

Update `parse_config(config_path)` in `inSVert/utils_sim.py` to handle `SNP`:

```python
        # Iterate through the variants defined in YAML
        for sv_type, settings in config['variants'].items():
            # SNP: ratio-based, no length distribution
            if sv_type == 'SNP':
                snp_ratio = settings.get('ratio')
                if snp_ratio is None:
                    raise ValueError("Config Error: SNP requires a 'ratio' parameter (e.g., 0.0001)")
                if not (0.0 < snp_ratio < 1.0):
                    raise ValueError(f"Config Error: SNP 'ratio' must be between 0 and 1, got {snp_ratio}")
                tstv_ratio = settings.get('tstv_ratio', 2.0)
                if tstv_ratio < 0:
                    raise ValueError(f"Config Error: SNP 'tstv_ratio' must be >= 0, got {tstv_ratio}")
                sv_data['SNP'] = {
                    'ratio': snp_ratio,
                    'tstv_ratio': tstv_ratio
                }
                continue
```

### 2.2 Transition / Transversion Helper (`pick_snp_alt`)

Add `pick_snp_alt(ref_base, tstv_ratio)` to `inSVert/utils_sim.py`:

```python
def pick_snp_alt(ref_base: str, tstv_ratio: float):
    """
    Given a reference base and a transition/transversion ratio, returns a
    randomly selected ALT base respecting the Ts/Tv ratio.
    
    Transitions (Ts): A <-> G, C <-> T
    Transversions (Tv): A <-> {C, T}, G <-> {C, T}, C <-> {A, G}, T <-> {A, G}
    
    P(transition) = tstv_ratio / (tstv_ratio + 1)
    P(transversion) = 1 / (tstv_ratio + 1), split equally between 2 options.
    
    Returns None if ref_base is not A/C/G/T (e.g., N).
    """
    TRANSITIONS = {'A': 'G', 'G': 'A', 'C': 'T', 'T': 'C'}
    TRANSVERSIONS = {'A': ['C', 'T'], 'G': ['C', 'T'], 'C': ['A', 'G'], 'T': ['A', 'G']}
    
    ref = ref_base.upper()
    if ref not in TRANSITIONS:
        return None
    
    p_transition = tstv_ratio / (tstv_ratio + 1)
    if random.random() < p_transition:
        return TRANSITIONS[ref]
    else:
        return random.choice(TRANSVERSIONS[ref])
```

### 2.3 VCF Header Line (`buildheader`)

In `buildheader()`, add `VT` header definition:

```python
##INFO=<ID=VT,Number=1,Type=String,Description="Variant Type">
```

---

## 3. `simulate.py` — SNP Placement Loop

1. Exclude `SNP` from the standard SV length loop.
2. Run SNP generation after all SVs complete.
3. Use a 3-attempt retry loop per SNP (matching existing SV placement behavior).

```python
        # Exclude SNP from standard SV loop
        sv_types = [k for k in fakedict if k != "SNP"]
        for svtype in sv_types:
            # ... existing SV simulation loop ...

        # ---------------------------------------------------------------
        # SNP SIMULATION (processed after all SVs to avoid collisions)
        # ---------------------------------------------------------------
        if "SNP" in fakedict:
            snp_ratio = fakedict['SNP']['ratio']
            tstv_ratio = fakedict['SNP'].get('tstv_ratio', 2.0)
            total_genome_length = sum(lengths)
            n_snps = int(total_genome_length * snp_ratio)
            
            print(f"Generating {n_snps} SNPs (ratio={snp_ratio}, Ts/Tv={tstv_ratio})...")
            
            snp_count = 0
            for i in range(n_snps):
                placed = False
                for attempt in range(3):
                    chrom, chrom_length = utils_sim.select_chr(chroms, lengths)
                    pos = utils_sim.select_pos(chrom_length)
                    gt = utils_sim.generate_genotype(ploidy, heterozygosity)
                    
                    if (utils_sim.overlaps(chrom, pos, pos + 1, gt, sv_positions) or
                        utils_sim.overlaps_excluded_region(chrom, pos, pos + 1, excluded_regions)):
                        continue
                    
                    ref_base = utils_sim.fetch_ref_base(chrom, pos, ref_fasta).upper()
                    alt_base = utils_sim.pick_snp_alt(ref_base, tstv_ratio)
                    if alt_base is None:
                        continue  # Skip non-ACGT bases (N, etc.)
                    
                    placed = True
                    break
                
                if not placed:
                    print(f"SNP n: {i+1} could not be placed after 3 attempts, skipping")
                    continue
                
                snp_count += 1
                snp_id = f'inSVert.SNP.{snp_count}'
                snp = VariantObjects.SNP(chrom, pos, snp_id, gt, ref_base, alt_base)
                
                # Track single-base interval [pos, pos + 1)
                utils_sim.track_sv(sv_positions, chrom, pos, pos + 1, gt)
                vcf.write(snp.format() + '\n')
            
            print(f"Placed {snp_count}/{n_snps} SNPs")
```

---

## 4. `insert.py` — Single-Base Substitution Engine

In `inSVert/insert.py`, inspect `var` inside the variant stream loop. Detect SNPs using single-base structural signature (`svtype is None and len(REF) == 1 and len(ALT) == 1`):

```python
                        svtype = var.info.get("SVTYPE")
                        
                        # ── SNP HANDLING (Single-base substitution, no padding base) ──
                        if svtype is None:
                            ref_allele = var.ref
                            alt_alleles = var.alts
                            if (ref_allele and alt_alleles and 
                                len(ref_allele) == 1 and len(alt_alleles) == 1 and 
                                len(str(alt_alleles[0])) == 1):
                                
                                # Write ALT nucleotide in place of reference base
                                writer.write(str(alt_alleles[0]))
                                ref_pos = start + 1
                                continue
```

---

## 5. `config.yaml` — Configuration Example

Add example `SNP` settings to `config.yaml`:

```yaml
variants:
  # ... existing SVs ...

  # Single-Nucleotide Polymorphisms
  SNP:
    ratio: 0.0001          # SNP density: 1 per 10,000 bp
    tstv_ratio: 2.0        # Transition/Transversion ratio (~2.0 for human)
```

---

## 6. Verification and Unit Testing Plan

Create `tests/test_snp.py` containing:
1. `test_snp_formatting()`: Verifies `SNP.format()` outputs valid VCF 4.2 string with `VT=SNP`.
2. `test_pick_snp_alt()`: Verifies transition vs transversion outputs for A/C/G/T and returns `None` for N.
3. `test_parse_config_snp()`: Tests config parsing validation for `ratio` and `tstv_ratio`.
4. `test_end_to_end_snp_simulation_and_insertion()`: Runs `simulate` + `insert` on synthetic reference FASTA, confirming SNP substitution in output FASTA and zero coordinate clashes.
