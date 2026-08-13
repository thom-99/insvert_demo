import os 
import math
import pysam
import yaml


### INPUT CONFIG VALIDATION ###

def validate_config(config_path):
    """Validate a simulation YAML config and report defaults and ignored fields."""
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file not found at: {config_path}")
    if not os.path.isfile(config_path):
        raise ValueError(f"Config path is not a file: {config_path}")

    try:
        with open(config_path, 'r') as config_file:
            config = yaml.safe_load(config_file)
    except (OSError, yaml.YAMLError) as e:
        raise ValueError(f"Config Validation Failed: could not parse YAML. Details: {e}")

    errors = []
    warnings = []

    def is_number(value):
        return (
            not isinstance(value, bool)
            and isinstance(value, (int, float))
            and math.isfinite(value)
        )

    def warn_default(path, default):
        warnings.append(f"'{path}' is not specified; using default {default}.")

    if not isinstance(config, dict):
        raise ValueError("Config Validation Failed: the YAML root must be a mapping.")

    known_top_level = {'genome', 'variants'}
    for key in config.keys() - known_top_level:
        warnings.append(f"Unknown top-level setting '{key}' will be ignored.")

    genome = config.get('genome')
    if not isinstance(genome, dict):
        errors.append("'genome' must be a mapping containing genome settings.")
    else:
        for key in genome.keys() - {'ploidy', 'heterozygosity'}:
            warnings.append(f"Unknown setting 'genome.{key}' will be ignored.")

        if 'ploidy' not in genome:
            warn_default('genome.ploidy', 2)
        else:
            ploidy = genome['ploidy']
            if isinstance(ploidy, bool) or not isinstance(ploidy, int) or ploidy < 1:
                errors.append("'genome.ploidy' must be an integer greater than or equal to 1.")

        if 'heterozygosity' not in genome:
            warn_default('genome.heterozygosity', 0.5)
        else:
            heterozygosity = genome['heterozygosity']
            if not is_number(heterozygosity) or not 0.0 <= heterozygosity <= 1.0:
                errors.append("'genome.heterozygosity' must be a number between 0 and 1.")

    variants = config.get('variants')
    if not isinstance(variants, dict) or not variants:
        errors.append("'variants' must be a non-empty mapping.")
        variants = {}

    supported_variants = {
        'INS', 'DEL', 'INV', 'DUP', 'TRA_COPY', 'TRA_CUT', 'SNP', 'MNP'
    }

    for variant_type, settings in variants.items():
        path = f"variants.{variant_type}"

        if variant_type not in supported_variants:
            errors.append(
                f"Unsupported variant type '{variant_type}'. Supported types are: "
                f"{', '.join(sorted(supported_variants))}."
            )
            continue
        if not isinstance(settings, dict):
            errors.append(f"'{path}' must be a mapping.")
            continue

        if variant_type in {'SNP', 'MNP'}:
            for key in settings.keys() - {'ratio', 'tstv_ratio'}:
                warnings.append(f"Unknown setting '{path}.{key}' will be ignored.")

            ratio = settings.get('ratio')
            if not is_number(ratio) or not 0.0 < ratio < 1.0:
                errors.append(f"'{path}.ratio' is required and must be a number between 0 and 1.")

            if 'tstv_ratio' not in settings:
                warn_default(f'{path}.tstv_ratio', 2.0)
            else:
                tstv_ratio = settings['tstv_ratio']
                if not is_number(tstv_ratio) or tstv_ratio < 0:
                    errors.append(f"'{path}.tstv_ratio' must be a number greater than or equal to 0.")
            continue

        allowed_settings = {'count', 'distribution', 'parameters'}
        if variant_type == 'DUP':
            allowed_settings.add('copy_number')
        if variant_type in {'TRA_COPY', 'TRA_CUT'}:
            allowed_settings.add('reverse_ratio')
        for key in settings.keys() - allowed_settings:
            warnings.append(f"Unknown setting '{path}.{key}' will be ignored.")

        count = settings.get('count')
        if isinstance(count, bool) or not isinstance(count, int) or count < 0:
            errors.append(f"'{path}.count' is required and must be an integer greater than or equal to 0.")

        distribution = settings.get('distribution')
        if distribution not in {'normal', 'pareto'}:
            errors.append(f"'{path}.distribution' must be either 'normal' or 'pareto'.")

        parameters = settings.get('parameters')
        if not isinstance(parameters, dict):
            errors.append(f"'{path}.parameters' is required and must be a mapping.")
            parameters = {}

        allowed_parameters = {'median_length', 'min_length', 'max_length', 'sigma'}
        if variant_type in {'TRA_COPY', 'TRA_CUT'}:
            allowed_parameters.add('reverse_ratio')
        for key in parameters.keys() - allowed_parameters:
            warnings.append(f"Unknown setting '{path}.parameters.{key}' will be ignored.")

        lengths = {}
        for name in ('median_length', 'min_length', 'max_length'):
            value = parameters.get(name)
            if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
                errors.append(f"'{path}.parameters.{name}' is required and must be a positive integer.")
            else:
                lengths[name] = value

        if len(lengths) == 3:
            minimum = lengths['min_length']
            median = lengths['median_length']
            maximum = lengths['max_length']
            if minimum > maximum:
                errors.append(f"'{path}' requires min_length to be less than or equal to max_length.")
            if distribution == 'pareto' and not minimum < median <= maximum:
                errors.append(
                    f"'{path}' with a Pareto distribution requires min_length < median_length <= max_length."
                )
            elif distribution == 'normal' and not minimum <= median <= maximum:
                errors.append(
                    f"'{path}' with a normal distribution requires min_length <= median_length <= max_length."
                )

        sigma = parameters.get('sigma')
        if distribution == 'normal':
            if sigma is None:
                warnings.append(
                    f"'{path}.parameters.sigma' is not specified; using default 10% of median_length."
                )
            elif not is_number(sigma) or sigma < 0:
                errors.append(f"'{path}.parameters.sigma' must be a number greater than or equal to 0.")
        elif sigma is not None:
            warnings.append(f"'{path}.parameters.sigma' is ignored for a Pareto distribution.")

        if variant_type == 'DUP':
            copy_number = settings.get('copy_number')
            if not isinstance(copy_number, dict):
                errors.append(f"'{path}.copy_number' is required and must be a mapping.")
            else:
                for key in copy_number.keys() - {'min', 'max', 'weights'}:
                    warnings.append(f"Unknown setting '{path}.copy_number.{key}' will be ignored.")

                cn_min = copy_number.get('min')
                cn_max = copy_number.get('max')
                valid_min = isinstance(cn_min, int) and not isinstance(cn_min, bool) and cn_min >= 2
                valid_max = isinstance(cn_max, int) and not isinstance(cn_max, bool) and cn_max >= 2
                if not valid_min:
                    errors.append(f"'{path}.copy_number.min' is required and must be an integer >= 2.")
                if not valid_max:
                    errors.append(f"'{path}.copy_number.max' is required and must be an integer >= 2.")
                if valid_min and valid_max and cn_min > cn_max:
                    errors.append(f"'{path}.copy_number.min' cannot be greater than copy_number.max.")

                weights = copy_number.get('weights')
                if weights is None:
                    warnings.append(f"'{path}.copy_number.weights' is not specified; using uniform weights.")
                elif not isinstance(weights, list) or not weights:
                    errors.append(f"'{path}.copy_number.weights' must be a non-empty list.")
                elif not all(is_number(weight) and weight >= 0 for weight in weights):
                    errors.append(f"'{path}.copy_number.weights' must contain only non-negative numbers.")
                else:
                    if valid_min and valid_max and len(weights) != cn_max - cn_min + 1:
                        errors.append(
                            f"'{path}.copy_number.weights' must contain one value for each copy number "
                            f"from {cn_min} to {cn_max}."
                        )
                    if sum(weights) <= 0:
                        errors.append(f"'{path}.copy_number.weights' must contain at least one positive value.")

        if variant_type in {'TRA_COPY', 'TRA_CUT'}:
            top_level_ratio = settings.get('reverse_ratio')
            nested_ratio = parameters.get('reverse_ratio')
            if top_level_ratio is not None and nested_ratio is not None:
                warnings.append(
                    f"Both '{path}.reverse_ratio' and '{path}.parameters.reverse_ratio' are set; "
                    "the top-level value will be used."
                )
            reverse_ratio = top_level_ratio if top_level_ratio is not None else nested_ratio
            if reverse_ratio is None:
                warn_default(f'{path}.reverse_ratio', 0.0)
            elif not is_number(reverse_ratio) or not 0.0 <= reverse_ratio <= 1.0:
                errors.append(f"'{path}.reverse_ratio' must be a number between 0 and 1.")

    for warning in warnings:
        print(f"[Warning] {warning}")

    if errors:
        details = '\n'.join(f"- {error}" for error in errors)
        raise ValueError(f"Config Validation Failed:\n{details}")

    return True


### INPUT VCF VALIDATION ###

def validate_vcf(vcf_path, fasta_path, target_ploidy=None):
    """
    Validates the input VCF for structural integrity, compatibility with the reference FASTA,
    and genotype ploidy.
    """
    if not os.path.exists(vcf_path):
        raise FileNotFoundError(f"Input VCF file not found at: {vcf_path}")
    
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"Reference FASTA file not found at: {fasta_path}")

    try:
        # Load both files
        with pysam.VariantFile(vcf_path) as vcf, pysam.FastaFile(fasta_path) as fasta:
            vcf_contigs = vcf.header.contigs
            fasta_chroms = set(fasta.references)

            # Check 1: Do the VCF chromosomes exist in the FASTA?
            for vcf_chrom in vcf_contigs:
                if vcf_chrom not in fasta_chroms:
                    raise ValueError(f"Chromosome '{vcf_chrom}' found in VCF but missing from reference FASTA. "
                                     f"Are you sure the VCF was called against this exact FASTA?")

            # Check 2: Origin FASTA matching (Compare contig lengths)
            for vcf_chrom, contig_record in vcf_contigs.items():
                if contig_record.length is not None:
                    fasta_length = fasta.get_reference_length(vcf_chrom)
                    if contig_record.length != fasta_length:
                        raise ValueError(
                            f"Length mismatch for '{vcf_chrom}'. "
                            f"VCF expects {contig_record.length}bp, but FASTA has {fasta_length}bp. "
                            f"Are you sure the VCF was called against this exact FASTA?"
                        )

            if target_ploidy is not None:
                for rec in vcf:

                    # Check 3: presence of a sample 
                    if not rec.samples:
                        raise ValueError(f"VCF record at {rec.chrom}:{rec.pos} has no sample genotype.")

                    # Check 4: presence of a genotype
                    gt = rec.samples[0].get('GT')
                    if gt is None or len(gt) == 0:
                        raise ValueError(f"VCF record at {rec.chrom}:{rec.pos} has no GT value")

                    # Check 5: ploidy matching with target_ploidy
                    if len(gt) != target_ploidy:
                        raise ValueError(f"Ploidy mismatch at {rec.chrom}:{rec.pos}, GT={gt} has {len(gt)} alleles but --ploidy {target_ploidy} was specified")
                        

    except Exception as e:
        raise ValueError(f"VCF Validation Failed: {e}")




### INPUT FASTA VALIDATION ###

def validate_fasta(fasta_path):
    """Ensures the reference FASTA exists and can be parsed."""
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"Reference FASTA file not found at: {fasta_path}")
    try:
        with pysam.FastaFile(fasta_path) as fasta:
            if not fasta.references:
                raise ValueError(f"FASTA file '{fasta_path}' contains no sequences.")
    except Exception as e:
        raise ValueError(f"Failed to parse FASTA file. Details: {e}")
    return True


### INPUT BED VALIDATION ###

def validate_bed(bed_path, fasta_path=None):
    """Validates BED format (3 tab-separated columns) and alerts if chroms mismatch."""
    if not bed_path:
        return True
    if not os.path.exists(bed_path):
        raise FileNotFoundError(f"BED exclusion file not found at: {bed_path}")

    fasta_chroms = set()
    if fasta_path and os.path.exists(fasta_path):
        with pysam.FastaFile(fasta_path) as fasta:
            fasta_chroms = set(fasta.references)

    with open(bed_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#') or not line.strip():
                continue
            columns = line.strip().split('\t')
            if len(columns) < 3:
                raise ValueError(f"Malformed BED line {line_num}: Must have 3 tab-separated columns.")
            
            chrom, start, end = columns[0], columns[1], columns[2]
            if not start.isdigit() or not end.isdigit():
                raise ValueError(f"Malformed BED line {line_num}: Coordinates must be integers.")
            if int(start) >= int(end):
                raise ValueError(f"Malformed BED line {line_num}: Start coordinate cannot be >= End.")
            if fasta_chroms and chrom not in fasta_chroms:
                print(f"[Warning] BED line {line_num}: Chromosome '{chrom}' not in reference FASTA.")
    return True




def validate_output_path(reference_path, output_path):
    """
    Ensures that the output FASTA does not overwrite the input reference FASTA.
    """
    reference_real = os.path.realpath(os.path.abspath(reference_path))
    output_real = os.path.realpath(os.path.abspath(output_path))

    if reference_real == output_real:
        raise ValueError(
            "Output FASTA cannot be the same file as the reference FASTA. "
            "Please specify a different output path with -o/--output.")

    return True
