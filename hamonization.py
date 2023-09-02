# %%
#import required modules
import gzip
import pandas as pd
from pyensembl import EnsemblRelease
import requests, sys
from liftover import get_lifter


# %%
#setup data - using pandas to leverage vectorization
#https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle
#https://ftp.ensembl.org/pub/release-110/variation/vcf/homo_sapiens/
def get_vcf_names(vcf_path):
    with gzip.open(vcf_path, "rt") as ifile:
          for line in ifile:
            if line.startswith("#CHROM"):
                  vcf_names = [x for x in line.split('\t')]
                  break
    ifile.close()
    return vcf_names
vcf_file = '1000GENOMES-phase_3.vcf.gz'
names = get_vcf_names(vcf_file)

df = pd.read_csv('ukb_phase1to3_heart_may_2022_pheno10.fastGWA', delimiter='\t') 
vcf = pd.read_csv(vcf_file, compression='gzip', comment='#',delim_whitespace=True, header=None, names=names) 
server = "https://rest.ensembl.org"
converter = get_lifter('hg19', 'hg38')

# %%
#for testing only
DF = df.copy()

# %% [markdown]
# Step 1: Mapping variant IDs to locations
# - Update base pair location value by mapping rsID using latest Ensembl release Q(‘hm_coordinate_conversion’ = ‘rs’); or
# - if above not possible, liftover base pair location to latest genome build (‘hm_coordinate_conversion’ = ‘lo’); or
# - if above not possible, remove variant from file.

# %%
def rsIDtoEnsembl(row):
    #attempt with local file
    local_vcf_query = vcf.loc[vcf['ID'] == row['SNP']]
    if not local_vcf_query.empty:
        return [local_vcf_query['#CHROM'].iloc[0], local_vcf_query['POS'].iloc[0]]
    #attempt with ensembl API
    ext = '/variation/human/' + row['SNP'] + '?'
    r = requests.get(server+ext, headers={ 'Content-Type' : 'application/json'})
    if not r.ok:
        return [None, None]
    decoded = r.json()
    if 'failed' in decoded:
        return [None, None]
    try:
        return list(map(int, decoded['mappings'][0]['location'].split('-', 1)[0].split(':', 1)))
    except:
        print(decoded)

def liftoverM(row):
    ret = converter.convert_coordinate(row['CHR'], row['POS'])
    if (ret == None or len(ret) == 0):
        return [None, None]
    return [ret[0][0][3], ret[0][1]]

def mapIDs(row):
    [chrom, pos] = rsIDtoEnsembl(row)
    if (chrom == None or pos == None):
        [chrom, pos] = liftoverM(row)
    if (chrom == None or pos == None):
        return None
    row['#CHROM'] = chrom
    row['POS'] = pos
    return row

DF = DF[DF.apply(mapIDs, axis=1) != None]

# %% [markdown]
# Step 2: Orientation, part(a): Infer the orientation of palindromic variants
# - Firstly, we randomly select 10% of sites. The effect and other alleles are compared with counterpart alternative and reference alleles in the Ensembl VCF references to identify the strand of the non-palindromic variants (forward or reverse).
# - The forward strand consensus can be calculated by forward/(forward+reverse) or reverse/(forward+reverse). To avoid any possibility of sampling bias:
#     - If the rate ≥ 0.995, the following harmonisation steps on the palindromic variants are inferred as on the forward (or reverse) strand;
#     - If the rate is the range of (0.995,0.9), this rate is recalculated by all non-palindromic variants in the data. The palindromic variants can be inferred as forward (or reverse) if the recalculated rate > 0.99, otherwise palindromic variants are dropped for harmonisation;
#     - If the rate ≤ 0.9, palindromic variants are dropped in the following harmonisation step.
# 
# 

# %%
class Counter: #create function class for vectorized implementation
    def __init__(self, seed):
        self.counter = 0
        self.valid = seed
    def fun(self, l):
        f = l[0]
        v = l[1]
        chrom = self['CHR']
        loc = self['POS']
        ref = self['A1']
        alt = self['A2']
        reference = vcf.loc[vcf['CHR'] == chrom and vcf['POS'] == loc]
        if reference.empty:
            return [self.counter + f, self.valid + v - 1]
        else:
            if ref == reference['REF'] and alt == reference['ALT']:
                return [self.counter + f + 1, self.valid + v]
            else:
                return [self.counter + f, self.valid + v]

def get_consensus(df):
    counter = Counter(len(df.index))
    count = df.apply(counter.fun)
    forward_consensus = count[0]/(count[0] + count[1])
    return forward_consensus

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement

def isPalindrome(variant):
    return variant == reverse_complement(variant)

def infer_orientation(df):
    sample_consenus = get_consensus(df.sample(frac=0.1))
    if forward_consensus >= 0.995:
        return 1 #forward
    elif forward_consensus <= 0.005:
        return -1 #reverse
    elif forward_consensus > 0.9 or forward_consensus < 0.1:
        sample_consenus = get_consensus(df[isPalindrome(df['A2'])])
        if forward_consensus >= 0.99:
            return 1 #forward
        elif forward_consensus <= 0.01:
            return -1
    else:
        return 0 #ambiguous

def harmonize_palindrome(row, orientation):
    chrom = row['CHR']
    loc = row['POS']
    ref = row['A1']
    alt = row['A2']
    reference = vcf.loc[vcf['CHR'] == chrom and vcf['POS'] == loc]
    if orientation == 1: 
        #assume on forward strand, no need to reverse complement
        if ref == alt:
            #alleles in correct orientation, no need to do anything
            return row
        else:
            #alleles in flipped orientation, need to flip strand
            return flip(row)
    elif orientation == -1:
        #assume on reverse strand, need to reverse complement
        row['A1'] = reverse_complement(ref)
        row['A2'] = reverse_complement(alt)
        if ref == reverse_complement(alt):
            #alleles in correct orientation, no need to do anything
            return row
        else:
            #alleles in flipped orientation, need to flip strand
            return flip(row)
    else:
        #ambiguous, do not harmonize
        return None
orientation = infer_orientation(DF)

# %% [markdown]
# Step 2: Orientation, part(b): Using chromosome, base pair location and the effect and other alleles, query each variant against the Ensembl VCF reference to harmonise as appropriate by either:
# - Variant harmonisation: Using chromosome, base pair location and the effect and other alleles, query each variant against the Ensembl VCF reference to harmonise as appropriate by either:
#     - keeping record as is because:
#         - it is already correctly orientated
#     - orientating to reference strand:
#         - reverse complement the effect and other alleles
#     - flipping the effect and other alleles
#         - because the effect and other alleles are flipped in the reference
#         - this also means the beta, odds ratio, 95% CI and effect allele frequency are inverted
#     - a combination of the orientating and flipping the alleles.
#     - replace with NA because:
#         - There is no counterpart record in the reference VCF file.

# %%
def flip(row):
    row['AF1'] = 1 - row['AF1']
    row['BETA'] = -row['BETA']
    row['SE'] = row['SE']
    row['P'] = row['P']
    return row

def harmonize(row):
    chrom = row['CHR']
    loc = row['POS']
    ref = row['A1']
    alt = row['A2']
    if isPalindrome(alt):
        return process_palindrome(row, orientation)
    reference = vcf.loc[vcf['CHR'] == chrom and vcf['POS'] == loc]
    if reference.empty:
        row['A1'] = 'NA'
        row['A2'] = 'NA'
        return row
    reference = reference.iloc[0]
    if ref == reference['A1'] and alt == reference['A2']:
        return row
    elif ref == reference['A2'] and alt == reference['A1']: 
        #alleles on forward strand, no need to reverse complement
        if ref == alt:
            #alleles in correct orientation, no need to do anything
            return row
        else:
            #alleles in flipped orientation, need to flip strand
            return flip(row)
    else:
        #alleles on reverse strand, need to reverse complement
        row['A1'] = reverse_complement(ref)
        row['A2'] = reverse_complement(alt)
        if ref == reverse_complement(alt):
            #alleles in correct orientation, no need to do anything
            return row
        else:
            #alleles in flipped orientation, need to flip strand
            return flip(row)

DF = DF.apply(harmonize, axis=1)  

