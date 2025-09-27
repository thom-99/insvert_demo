import pysam
import sys

vcf_path = sys.argv[1]

# add checks on the extension on the path to verify it is a vcf and it is not zipped

vcf = pysam.VariantFile(vcf_path)

# for fast access, data will be stored in a dictonary
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


vcf = pysam.VariantFile(vcf_path)

for SV in vcf:
    svtype = SV.info.get("SVTYPE")
    if svtype in sv_data:
        sv_length = SV.info.get("SVLEN")
        if sv_length is not None: 
            #svlength is negative for deletions
            sv_length = abs(sv_length)
            sv_data[svtype]['lengths'].append(sv_length)

        if svtype == 'DUP':
            copy_number = SV.info.get("CN")
            if copy_number is not None:
                sv_data[svtype]['copy_numbers'].append(copy_number)




# SUMMARY STATISTICS

print("Structural Variant Summary:")
print("-" * 30)
for sv_type in sv_data:
    count = len(sv_data[sv_type]['lengths'])
    print(f"{sv_type}: {count} variants")
    if count > 0:
        lengths = sv_data[sv_type]['lengths']
        print(f"  Length range: {min(lengths)} - {max(lengths)} bp")
        print(f"  Mean length: {sum(lengths) / len(lengths):.1f} bp")
        
        if sv_type == 'DUP' and sv_data[sv_type]['copy_numbers']:
            cn_values = sv_data[sv_type]['copy_numbers']
            print(f"  Copy numbers: {min(cn_values)} - {max(cn_values)}")
    print()

vcf.close()