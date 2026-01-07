import pysam
import tempfile
import os

# Create a dummy unsorted VCF
with open("test_input.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n")
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    f.write("chr1\t200\t.\tA\tT\t.\t.\t.\n")
    f.write("chr1\t100\t.\tA\tT\t.\t.\t.\n")

tf = tempfile.NamedTemporaryFile(delete=False, suffix=".vcf.gz")
tf.close()
temp_out = tf.name

print(f"Pysam version: {pysam.__version__}")
print("Attempting sort...")

try:
    pysam.sort("-o", temp_out, "test_input.vcf")
    print("Sort successful.")
    
    pysam.tabix_index(temp_out, preset="vcf", force=True)
    print("Index successful.")
except Exception as e:
    print(f"Python Error: {e}")

# Cleanup
if os.path.exists("test_input.vcf"): os.remove("test_input.vcf")
if os.path.exists(temp_out): os.remove(temp_out)
if os.path.exists(temp_out + ".tbi"): os.remove(temp_out + ".tbi")
