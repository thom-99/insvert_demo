import pysam

'''
parses a VCF file and gets relevant informations that are fed to simulator.py
it also allows to print some summary statistics of such informations
'''

# ADD TESTS FOR VALIDATING VCF 

def parse_vcf(vcf_path:str):
    #parses a vcf and inserts information into a dictonary
    vcf = pysam.VariantFile(vcf_path)
    
    sv_data = {
        'INS': {
            'lengths': [],
        },
        'DEL': {
            'lengths': [],
        },
        'DUP': {
            'lengths': [],
            'copy_numbers': [],  
        },
        'INV': {
            'lengths': [],
        }
    }

    for SV in vcf:
        svtype = SV.info.get("SVTYPE")
        if svtype in sv_data:
            # get SV length
            sv_length = SV.info.get("SVLEN")
            if sv_length is not None:
                # SVLEN can be negative for deletions, take absolute value
                sv_length = abs(sv_length)
                sv_data[svtype]['lengths'].append(sv_length)
            
            # for duplications, get also copy number
            if svtype == 'DUP':
                copy_number = SV.info.get("CN")
                if copy_number is not None:
                    sv_data[svtype]['copy_numbers'].append(copy_number)

    vcf.close()
    return sv_data


def print_summary(sv_data:dict):
    # optional summary statistics for parsed SV
    print("Structural Variant Summary:")
    print("-" * 30)
    for sv_type in sv_data:
        count = len(sv_data[sv_type]['lengths'])
        print(f"{sv_type}: {count} variants")
        if count > 0:
            lengths = sv_data[sv_type]['lengths']
            print(f"  Length range: {min(lengths)} - {max(lengths)} bp")
            print(f"  Mean length: {int(sum(lengths) / len(lengths))} bp")
            
            if sv_type == 'DUP' and sv_data[sv_type]['copy_numbers']:
                cn_values = sv_data[sv_type]['copy_numbers']
                print(f"  Copy numbers: {min(cn_values)} - {max(cn_values)}")
        print()
    print("-" * 30)




# ensuring vcfparser still works as a command line script when running it directly
if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python vcf_parser.py <vcf_file>")
        sys.exit(1)
    
    vcf_path = sys.argv[1]
    sv_data = parse_vcf(vcf_path)
    print_summary(sv_data)