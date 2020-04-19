import os
import argparse
import re
import random
import math
import sys
import time
import datetime

#import helper functions
from helper import revcomp
from helper import motif_search_topstrand
from helper import motif_search_bothstrands
from helper import findmut
from helper import analyze_atg
from helper import analyze_ctg
from helper import analyze_kozak
from helper import gquad
from helper import motif_search_prte
from helper import fivetop
from helper import motif_search_jaspar
from helper import convert
from helper import getrefs
from helper import convert_back
from helper import num_to_triplet
from helper import triplet_to_num



cisbp_added = 0
cisbp_removed = 0
cisbp_affected = 0

rbpdb_added = 0
rbpdb_removed = 0
rbpdb_affected = 0

homer_added = 0
homer_removed = 0
homer_affected = 0

jaspar_added = 0
jaspar_removed = 0
jaspar_affected = 0

hughes_added = 0
hughes_removed = 0
hughes_affected = 0

atgs_added = 0
atgs_removed = 0

ctgs_added = 0
ctgs_removed = 0

kozaks_enhanced = 0
kozaks_weakened = 0

gquads = 0

prtes_added = 0
prtes_removed = 0
prtes_affected = 0

fivetops = 0

with open("homer_motifs_wmax.txt", "r") as homer_motif_file:
    homer_motifs = homer_motif_file.read().splitlines()
    homer_motif_file.close()

with open("jaspar_motifs_wmax.txt", "r") as jaspar_motif_file:
    jaspar_motifs = jaspar_motif_file.read().splitlines()
    jaspar_motif_file.close()

with open("compendium_motifs_wmax.txt", "r") as compendium_motif_file:
    compendium_motifs = compendium_motif_file.read().splitlines()
    compendium_motif_file.close()

with open("rbpdb_motifs_wmax.txt", "r") as rbpdb_motif_file:
    rbpdb_motifs = rbpdb_motif_file.read().splitlines()
    rbpdb_motif_file.close()

with open("cisbp_motifs_wmax.txt", "r") as cisbp_motif_file:
    cisbp_motifs = cisbp_motif_file.read().splitlines()
    cisbp_motif_file.close()

with open("prte_motif_wmax.txt", "r") as prte_pwm_file:
    prte_motif = prte_pwm_file.read().splitlines()
    prte_pwm_file.close()

#filename = input("UTR filename: ")
filename = "FINAL_MUFAM_MUTATIONS.txt"

with open(filename, "r") as UTR_file:
    UTRs = UTR_file.read().splitlines()
    UTR_file.close()


output_file = open("joint_element_analysis_{0}".format(sys.argv[1]), "w+")

##figure out probability of triplet mutation

ref_UTRs = getrefs(UTRs)


dic = {'____': 0, '___A': 0, '___C': 0, '___G': 0, '___T': 0, '__A_': 0, '__AA': 0.0, '__AC': 0.0, '__AG': 0.0, '__AT': 0.0, '__C_': 0, '__CA': 0.0, '__CC': 0.0, '__CG': 0.0, '__CT': 0.0, '__G_': 0, '__GA': 0.0, '__GC': 0.0, '__GG': 0.0, '__GT': 0.0, '__T_': 0, '__TA': 0.0, '__TC': 0.0, '__TG': 0.0, '__TT': 0.0, '_A__': 0, '_A_A': 0, '_A_C': 0, '_A_G': 0, '_A_T': 0, '_AA_': 0.0, '_AAA': 0.0, '_AAC': 0.0, '_AAG': 0.0, '_AAT': 0.0, '_AC_': 0.0, '_ACA': 0.0, '_ACC': 0.0, '_ACG': 0.0, '_ACT': 0.0, '_AG_': 0.0, '_AGA': 0.0, '_AGC': 0.0, '_AGG': 0.0, '_AGT': 0.0, '_AT_': 0.0, '_ATA': 0.0, '_ATC': 0.0, '_ATG': 0.0, '_ATT': 0.0, '_C__': 0, '_C_A': 0, '_C_C': 0, '_C_G': 0, '_C_T': 0, '_CA_': 0.0, '_CAA': 0.0, '_CAC': 0.0, '_CAG': 0.0, '_CAT': 0.0, '_CC_': 0.0, '_CCA': 0.0, '_CCC': 0.0, '_CCG': 0.0, '_CCT': 0.0, '_CG_': 0.0, '_CGA': 0.0, '_CGC': 0.0, '_CGG': 0.0, '_CGT': 0.0, '_CT_': 0.0, '_CTA': 0.0, '_CTC': 0.0, '_CTG': 0.0, '_CTT': 0.0, '_G__': 0, '_G_A': 0, '_G_C': 0, '_G_G': 0, '_G_T': 0, '_GA_': 0.0, '_GAA': 0.0, '_GAC': 0.0, '_GAG': 0.0, '_GAT': 0.0, '_GC_': 0.0, '_GCA': 0.0, '_GCC': 0.0, '_GCG': 0.0, '_GCT': 0.0, '_GG_': 0.0, '_GGA': 0.0, '_GGC': 0.0, '_GGG': 0.0, '_GGT': 0.0, '_GT_': 0.0, '_GTA': 0.0, '_GTC': 0.0, '_GTG': 0.0, '_GTT': 0.0, '_T__': 0, '_T_A': 0, '_T_C': 0, '_T_G': 0, '_T_T': 0, '_TA_': 0.0, '_TAA': 0.0, '_TAC': 0.0, '_TAG': 0.0, '_TAT': 0.0, '_TC_': 0.0, '_TCA': 0.0, '_TCC': 0.0, '_TCG': 0.0, '_TCT': 0.0, '_TG_': 0.0, '_TGA': 0.0, '_TGC': 0.0, '_TGG': 0.0, '_TGT': 0.0, '_TT_': 0.0, '_TTA': 0.0, '_TTC': 0.0, '_TTG': 0.0, '_TTT': 0.0, 'A___': 0, 'A__A': 0, 'A__C': 0, 'A__G': 0, 'A__T': 0, 'A_A_': 0, 'A_AA': 0.0, 'A_AC': 0.0, 'A_AG': 0.0, 'A_AT': 0.0, 'A_C_': 0, 'A_CA': 0.0, 'A_CC': 0.0, 'A_CG': 0.0, 'A_CT': 0.0, 'A_G_': 0, 'A_GA': 0.0, 'A_GC': 0.0, 'A_GG': 0.0, 'A_GT': 0.0, 'A_T_': 0, 'A_TA': 0.0, 'A_TC': 0.0, 'A_TG': 0.0, 'A_TT': 0.0, 'AA__': 0, 'AA_A': 0, 'AA_C': 0, 'AA_G': 0, 'AA_T': 0, 'AAA_': 0.0, 'AAAA': 0.0, 'AAAC': 0.0, 'AAAG': 0.0, 'AAAT': 0.0, 'AAC_': 0.0, 'AACA': 0.00037821482602118004, 'AACC': 0.0012350761630300535, 'AACG': 0.0008417508417508417, 'AACT': 0.0, 'AAG_': 0.0, 'AAGA': 0.0018367882445552348, 'AAGC': 0.001844532279314888, 'AAGG': 0.0007481296758104738, 'AAGT': 0.0008385744234800838, 'AAT_': 0.0, 'AATA': 0.0008038585209003215, 'AATC': 0.0, 'AATG': 0.0, 'AATT': 0.0009442870632672333, 'AC__': 0, 'AC_A': 0, 'AC_C': 0, 'AC_G': 0, 'AC_T': 0, 'ACA_': 0.0, 'ACAA': 0.0, 'ACAC': 0.0, 'ACAG': 0.0, 'ACAT': 0.0, 'ACC_': 0.0, 'ACCA': 0.003400396712949844, 'ACCC': 0.00041399296211964395, 'ACCG': 0.0031120331950207467, 'ACCT': 0.0009770395701025891, 'ACG_': 0.0, 'ACGA': 0.004395604395604396, 'ACGC': 0.006759040216289287, 'ACGG': 0.0035783994795055302, 'ACGT': 0.00897226753670473, 'ACT_': 0.0, 'ACTA': 0.0, 'ACTC': 0.004107272287992268, 'ACTG': 0.0006330449461911796, 'ACTT': 0.0003132832080200501, 'AG__': 0, 'AG_A': 0, 'AG_C': 0, 'AG_G': 0, 'AG_T': 0, 'AGA_': 0.0, 'AGAA': 0.0, 'AGAC': 0.0, 'AGAG': 0.0, 'AGAT': 0.0, 'AGC_': 0.0, 'AGCA': 0.0003259452411994785, 'AGCC': 0.0002304147465437788, 'AGCG': 0.0009545020680878142, 'AGCT': 0.00027739251040221914, 'AGG_': 0.0, 'AGGA': 0.0021342186388427793, 'AGGC': 0.0021151586368977674, 'AGGG': 0.0023835319609967496, 'AGGT': 0.000899685110211426, 'AGT_': 0.0, 'AGTA': 0.0008944543828264759, 'AGTC': 0.0, 'AGTG': 0.0006022282445046673, 'AGTT': 0.0, 'AT__': 0, 'AT_A': 0, 'AT_C': 0, 'AT_G': 0, 'AT_T': 0, 'ATA_': 0.0, 'ATAA': 0.0, 'ATAC': 0.0, 'ATAG': 0.0, 'ATAT': 0.0, 'ATC_': 0.0, 'ATCA': 0.0016168148746968471, 'ATCC': 0.0007926023778071334, 'ATCG': 0.005702066999287242, 'ATCT': 0.0005765350245027386, 'ATG_': 0.0, 'ATGA': 0.0003497726477789437, 'ATGC': 0.001275103602167676, 'ATGG': 0.000555247084952804, 'ATGT': 0.0025678650036683784, 'ATT_': 0.0, 'ATTA': 0.0005730659025787965, 'ATTC': 0.0, 'ATTG': 0.0, 'ATTT': 0.0004973887092762995, 'C___': 0, 'C__A': 0, 'C__C': 0, 'C__G': 0, 'C__T': 0, 'C_A_': 0, 'C_AA': 0.0, 'C_AC': 0.0, 'C_AG': 0.0, 'C_AT': 0.0, 'C_C_': 0, 'C_CA': 0.0, 'C_CC': 0.0, 'C_CG': 0.0, 'C_CT': 0.0, 'C_G_': 0, 'C_GA': 0.0, 'C_GC': 0.0, 'C_GG': 0.0, 'C_GT': 0.0, 'C_T_': 0, 'C_TA': 0.0, 'C_TC': 0.0, 'C_TG': 0.0, 'C_TT': 0.0, 'CA__': 0, 'CA_A': 0, 'CA_C': 0, 'CA_G': 0, 'CA_T': 0, 'CAA_': 0.0, 'CAAA': 0.0010655301012253596, 'CAAC': 0.000945179584120983, 'CAAG': 0.0006004202942059442, 'CAAT': 0.0, 'CAC_': 0.0, 'CACA': 0.0, 'CACC': 0.0, 'CACG': 0.0, 'CACT': 0.0, 'CAG_': 0.0, 'CAGA': 0.0018367882445552348, 'CAGC': 0.0, 'CAGG': 0.0004987531172069825, 'CAGT': 0.0004192872117400419, 'CAT_': 0.0, 'CATA': 0.001607717041800643, 'CATC': 0.000555247084952804, 'CATG': 0.0, 'CATT': 0.0, 'CC__': 0, 'CC_A': 0, 'CC_C': 0, 'CC_G': 0, 'CC_T': 0, 'CCA_': 0.0, 'CCAA': 0.0004154549231408392, 'CCAC': 0.003150157507875394, 'CCAG': 0.0008950548221078541, 'CCAT': 0.0, 'CCC_': 0.0, 'CCCA': 0.0, 'CCCC': 0.0, 'CCCG': 0.0, 'CCCT': 0.0, 'CCG_': 0.0, 'CCGA': 0.0014652014652014652, 'CCGC': 0.0013518080432578573, 'CCGG': 0.002277163305139883, 'CCGT': 0.0008156606851549756, 'CCT_': 0.0, 'CCTA': 0.0028469750889679717, 'CCTC': 0.0012080212611741967, 'CCTG': 0.0014771048744460858, 'CCTT': 0.0015664160401002505, 'CG__': 0, 'CG_A': 0, 'CG_C': 0, 'CG_G': 0, 'CG_T': 0, 'CGA_': 0.0, 'CGAA': 0.0006066120715802245, 'CGAC': 0.0, 'CGAG': 0.001254180602006689, 'CGAT': 0.0005558643690939411, 'CGC_': 0.0, 'CGCA': 0.0, 'CGCC': 0.0, 'CGCG': 0.0, 'CGCT': 0.0, 'CGG_': 0.0, 'CGGA': 0.00047427080863172874, 'CGGC': 0.0007050528789659225, 'CGGG': 0.00021668472372697725, 'CGGT': 0.001349527665317139, 'CGT_': 0.0, 'CGTA': 0.0008944543828264759, 'CGTC': 0.0, 'CGTG': 0.0006022282445046673, 'CGTT': 0.001443001443001443, 'CT__': 0, 'CT_A': 0, 'CT_C': 0, 'CT_G': 0, 'CT_T': 0, 'CTA_': 0.0, 'CTAA': 0.0006150061500615006, 'CTAC': 0.0, 'CTAG': 0.0, 'CTAT': 0.0007874015748031496, 'CTC_': 0.0, 'CTCA': 0.0, 'CTCC': 0.0, 'CTCG': 0.0, 'CTCT': 0.0, 'CTG_': 0.0, 'CTGA': 0.0003497726477789437, 'CTGC': 0.0, 'CTGG': 0.000277623542476402, 'CTGT': 0.0011005135730007337, 'CTT_': 0.0, 'CTTA': 0.001146131805157593, 'CTTC': 0.003493172435693871, 'CTTG': 0.0, 'CTTT': 0.002238249191743347, 'G___': 0, 'G__A': 0, 'G__C': 0, 'G__G': 0, 'G__T': 0, 'G_A_': 0, 'G_AA': 0.034482758620689655, 'G_AC': 0.0, 'G_AG': 0.0, 'G_AT': 0.037037037037037035, 'G_C_': 0, 'G_CA': 0.0, 'G_CC': 0.0, 'G_CG': 0.0, 'G_CT': 0.0, 'G_G_': 0, 'G_GA': 0.0, 'G_GC': 0.0, 'G_GG': 0.0, 'G_GT': 0.0, 'G_T_': 0, 'G_TA': 0.0, 'G_TC': 0.0, 'G_TG': 0.0, 'G_TT': 0.0, 'GA__': 0, 'GA_A': 0, 'GA_C': 0, 'GA_G': 0, 'GA_T': 0, 'GAA_': 0.0, 'GAAA': 0.0005327650506126798, 'GAAC': 0.0004725897920604915, 'GAAG': 0.0003002101471029721, 'GAAT': 0.0026680896478121665, 'GAC_': 0.0, 'GACA': 0.0018910741301059002, 'GACC': 0.0004116920543433512, 'GACG': 0.0008417508417508417, 'GACT': 0.0004347826086956522, 'GAG_': 0.0, 'GAGA': 0.0, 'GAGC': 0.0, 'GAGG': 0.0, 'GAGT': 0.0, 'GAT_': 0.0, 'GATA': 0.0040192926045016075, 'GATC': 0.000555247084952804, 'GATG': 0.004517221908526256, 'GATT': 0.00141643059490085, 'GC__': 0, 'GC_A': 0, 'GC_C': 0, 'GC_G': 0, 'GC_T': 0, 'GCA_': 0.0, 'GCAA': 0.0004154549231408392, 'GCAC': 0.0021001050052502626, 'GCAG': 0.0013425822331617813, 'GCAT': 0.0015228426395939086, 'GCC_': 0.018518518518518517, 'GCCA': 0.0, 'GCCC': 0.00020699648105982198, 'GCCG': 0.001037344398340249, 'GCCT': 0.0002442598925256473, 'GCG_': 0.0, 'GCGA': 0.0, 'GCGC': 0.0, 'GCGG': 0.0, 'GCGT': 0.0, 'GCT_': 0.0, 'GCTA': 0.0, 'GCTC': 0.0012080212611741967, 'GCTG': 0.0010550749103186326, 'GCTT': 0.0006265664160401002, 'GG__': 0, 'GG_A': 0, 'GG_C': 0, 'GG_G': 0, 'GG_T': 0, 'GGA_': 0.0, 'GGAA': 0.0, 'GGAC': 0.002514668901927913, 'GGAG': 0.0008361204013377926, 'GGAT': 0.0, 'GGC_': 0.0, 'GGCA': 0.000651890482398957, 'GGCC': 0.0004608294930875576, 'GGCG': 0.0012726694241170856, 'GGCT': 0.0005547850208044383, 'GGG_': 0.0, 'GGGA': 0.0, 'GGGC': 0.0, 'GGGG': 0.0, 'GGGT': 0.0, 'GGT_': 0.0, 'GGTA': 0.0, 'GGTC': 0.0014619883040935672, 'GGTG': 0.0027100271002710027, 'GGTT': 0.0, 'GT__': 0, 'GT_A': 0, 'GT_C': 0, 'GT_G': 0, 'GT_T': 0, 'GTA_': 0.0, 'GTAA': 0.0024600246002460025, 'GTAC': 0.0, 'GTAG': 0.0007112375533428165, 'GTAT': 0.0031496062992125984, 'GTC_': 0.0, 'GTCA': 0.0012126111560226355, 'GTCC': 0.0010568031704095112, 'GTCG': 0.0, 'GTCT': 0.0014413375612568463, 'GTG_': 0.0, 'GTGA': 0.0, 'GTGC': 0.0, 'GTGG': 0.0, 'GTGT': 0.0, 'GTT_': 0.0, 'GTTA': 0.0005730659025787965, 'GTTC': 0.00031756113051762465, 'GTTG': 0.0004, 'GTTT': 0.0004973887092762995, 'T___': 0, 'T__A': 0, 'T__C': 0, 'T__G': 0, 'T__T': 0, 'T_A_': 0, 'T_AA': 0.0, 'T_AC': 0.0, 'T_AG': 0.0, 'T_AT': 0.0, 'T_C_': 0, 'T_CA': 0.0, 'T_CC': 0.0, 'T_CG': 0.0, 'T_CT': 0.0, 'T_G_': 0, 'T_GA': 0.0, 'T_GC': 0.0, 'T_GG': 0.0, 'T_GT': 0.0, 'T_T_': 0, 'T_TA': 0.0, 'T_TC': 0.0, 'T_TG': 0.0, 'T_TT': 0.0, 'TA__': 0, 'TA_A': 0, 'TA_C': 0, 'TA_G': 0, 'TA_T': 0, 'TAA_': 0.0, 'TAAA': 0.0, 'TAAC': 0.0004725897920604915, 'TAAG': 0.0006004202942059442, 'TAAT': 0.0, 'TAC_': 0.0, 'TACA': 0.0015128593040847202, 'TACC': 0.0008233841086867024, 'TACG': 0.013468013468013467, 'TACT': 0.0008695652173913044, 'TAG_': 0.0, 'TAGA': 0.0013119916032537393, 'TAGC': 0.0007905138339920949, 'TAGG': 0.001745635910224439, 'TAGT': 0.0, 'TAT_': 0.0, 'TATA': 0.0, 'TATC': 0.0, 'TATG': 0.0, 'TATT': 0.0, 'TC__': 0, 'TC_A': 0, 'TC_C': 0, 'TC_G': 0, 'TC_T': 0, 'TCA_': 0.0, 'TCAA': 0.0, 'TCAC': 0.002450122506125306, 'TCAG': 0.0006712911165808906, 'TCAT': 0.0010152284263959391, 'TCC_': 0.0, 'TCCA': 0.0011334655709832814, 'TCCC': 0.0016559718484785758, 'TCCG': 0.008990318118948824, 'TCCT': 0.001709819247679531, 'TCG_': 0.0, 'TCGA': 0.003663003663003663, 'TCGC': 0.001013856032443393, 'TCGG': 0.0032530904359141183, 'TCGT': 0.0008156606851549756, 'TCT_': 0.0, 'TCTA': 0.0, 'TCTC': 0.0, 'TCTG': 0.0, 'TCTT': 0.0, 'TG__': 0, 'TG_A': 0, 'TG_C': 0, 'TG_G': 0, 'TG_T': 0, 'TGA_': 0.0, 'TGAA': 0.00030330603579011223, 'TGAC': 0.0, 'TGAG': 0.0006270903010033445, 'TGAT': 0.0005558643690939411, 'TGC_': 0.0, 'TGCA': 0.001303780964797914, 'TGCC': 0.0013824884792626728, 'TGCG': 0.0066815144766146995, 'TGCT': 0.0011095700416088765, 'TGG_': 0.0, 'TGGA': 0.00047427080863172874, 'TGGC': 0.0007050528789659225, 'TGGG': 0.0006500541711809318, 'TGGT': 0.001349527665317139, 'TGT_': 0.0, 'TGTA': 0.0, 'TGTC': 0.0, 'TGTG': 0.0, 'TGTT': 0.0, 'TT__': 0, 'TT_A': 0, 'TT_C': 0, 'TT_G': 0, 'TT_T': 0, 'TTA_': 0.0, 'TTAA': 0.0, 'TTAC': 0.0, 'TTAG': 0.0, 'TTAT': 0.0015748031496062992, 'TTC_': 0.0, 'TTCA': 0.002021018593371059, 'TTCC': 0.003170409511228534, 'TTCG': 0.005702066999287242, 'TTCT': 0.0037474776592678004, 'TTG_': 0.0, 'TTGA': 0.0027981811822315496, 'TTGC': 0.000956327701625757, 'TTGG': 0.001665741254858412, 'TTGT': 0.0022010271460014674, 'TTT_': 0.0, 'TTTA': 0.0, 'TTTC': 0.0, 'TTTG': 0.0, 'TTTT': 0.0}









for var in range(1):

    cisbp_added = 0
    cisbp_removed = 0
    cisbp_affected = 0

    rbpdb_added = 0
    rbpdb_removed = 0
    rbpdb_affected = 0

    homer_added = 0
    homer_removed = 0
    homer_affected = 0

    jaspar_added = 0
    jaspar_removed = 0
    jaspar_affected = 0

    hughes_added = 0
    hughes_removed = 0
    hughes_affected = 0

    atgs_added = 0
    atgs_removed = 0

    ctgs_added = 0
    ctgs_removed = 0

    kozaks_enhanced = 0
    kozaks_weakened = 0

    gquads = 0

    prtes_added = 0
    prtes_removed = 0
    prtes_affected = 0

    fivetops = 0
    
    total_mutations = 0    

    for i in ref_UTRs: 
        wildtype = i.split()[1]
        print("mutants counter: " + str(total_mutations))
        mutants_counter = 0
        mutants_list = []
        index_list = []
        length  = len(wildtype)
        
        #print("length is: " + str(length))

        for var in range(length):

        
            if (var == 0):
                if (random.uniform(0.00000000000000000000, 1.00000000000000000000) < dic[("A_" + wildtype[0:2])]):
                    mutants_list.append(str("A" + wildtype[1:length]))
                    index_list.append(var)
                    mutants_counter += 1
                if (random.uniform(0.00000000000000000000, 1.00000000000000000000) < dic[("G_" + wildtype[0:2])]):
                    mutants_list.append(str("G" + wildtype[1:length]))
                    index_list.append(var)
                    mutants_counter += 1
                if (random.uniform(0.00000000000000000000, 1.00000000000000000000) < dic[("C_" + wildtype[0:2])]):
                    mutants_list.append(str("C" + wildtype[1:length]))
                    index_list.append(var)
                    mutants_counter += 1
                if (random.uniform(0.00000000000000000000, 1.00000000000000000000) < dic[("T_" + wildtype[0:2])]):
                    mutants_list.append(str("T" + wildtype[1:length]))
                    index_list.append(var)
                    mutants_counter += 1

            elif (var == (len(wildtype) - 1)):
                if (random.uniform(0.00000000000000000000, 1.00000000000000000000) < dic[("A" + wildtype[var-1:var+1] + "_")]):
                    mutants_list.append(str(wildtype[0:var] + "A"))
                    index_list.append( var)
                    mutants_counter += 1
                if (random.uniform(0.00000000000000000000, 1.00000000000000000000) < dic[("G" + wildtype[var-1:var+1] + "_")]):
                    mutants_list.append( str(wildtype[0:var] + "G"))
                    index_list.append(var)
                    mutants_counter += 1
                if (random.uniform(0.00000000000000000000, 1.00000000000000000000) < dic[("C" + wildtype[var-1:var+1] + "_")]):
                    mutants_list.append(str(wildtype[0:var] + "C"))
                    index_list.append(var)
                    mutants_counter += 1
                if (random.uniform(0.00000000000000000000, 1.00000000000000000000) < dic[("G" + wildtype[var-1:var+1] + "_")]):
                    mutants_list.append(str(wildtype[0:var] + "G"))
                    index_list.append(var)
                    mutants_counter += 1

            else:
                if (random.uniform(0.00000000000000000000, 1.00000000000000000000) < dic[("A" + wildtype[var-1:var+2])]):
                    mutants_list.append(str(wildtype[0:var] + "A" + wildtype[var + 1 : length]))
                    index_list.append(var)
                    mutants_counter += 1
                if (random.uniform(0.00000000000000000000, 1.00000000000000000000) < dic[("T" + wildtype[var-1:var+2])]):
                    mutants_list.append(str(wildtype[0:var] + "T" + wildtype[var + 1 : length]))
                    index_list.append(var)
                    mutants_counter += 1
                if (random.uniform(0.00000000000000000000, 1.00000000000000000000) < dic[("G" + wildtype[var-1:var+2])]):
                    mutants_list.append(str(wildtype[0:var] + "G" + wildtype[var + 1 : length]))
                    index_list.append(var)
                    mutants_counter += 1
                if (random.uniform(0.00000000000000000000, 1.00000000000000000000) < dic[("C" + wildtype[var-1:var+2])]):
                    mutants_list.append(str(wildtype[0:var] + "C" + wildtype[var + 1 : length]))
                    index_list.append(var)
                    mutants_counter += 1

        #print(str(mutants_counter) + "mutations in this utr")
        total_mutations = total_mutations + mutants_counter

        for mutx in range(mutants_counter):
            x = index_list[mutx]
            mutant = mutants_list[mutx]

            #print("Wildtype: " + str(wildtype[x-5:x+6]))
            #print("Mutant: " + str(mutant[x-5:x+6]))

            #RBPDB
            indicator = motif_search_topstrand(rbpdb_motifs, wildtype, mutant, x)
            if indicator == 0:
                #do nothing
                p = 1
            elif indicator == 1:
                #added
                rbpdb_added = rbpdb_added + 1
                rbpdb_affected = rbpdb_affected + 1
            elif indicator == 2:
                #removed
                rbpdb_removed = rbpdb_removed + 1
                rbpdb_affected = rbpdb_affected + 1
            elif indicator == 3:
                rbpdb_affected = rbpdb_affected + 1

            #CISBP
            indicator = motif_search_topstrand(cisbp_motifs, wildtype, mutant, x)
            if indicator == 0:
                #do nothing
                p = 1
            elif indicator == 1:
                #added
                cisbp_added = cisbp_added + 1
                #print(cisbp_added)
                cisbp_affected = cisbp_affected + 1
            elif indicator == 2:
                #removed
                cisbp_removed = cisbp_removed + 1
                cisbp_affected = cisbp_affected + 1
            elif indicator == 3:
                cisbp_affected = cisbp_affected + 1

            #HUGHES
            indicator = motif_search_topstrand(compendium_motifs, wildtype, mutant, x)
            if indicator == 0:
                #do nothing
                p = 1
            elif indicator == 1:
                #added
                hughes_added = hughes_added + 1
                hughes_affected = hughes_affected + 1
            elif indicator == 2:
                #removed
                hughes_removed = hughes_removed + 1
                hughes_affected = hughes_affected + 1
            elif indicator == 3:
                hughes_affected = hughes_affected + 1

            #HOMER
            indicator = motif_search_bothstrands(homer_motifs, wildtype, mutant, x)
            if indicator == 0:
                #do nothing
                p = 1
            elif indicator == 1:
                #added
                homer_added = homer_added + 1
                homer_affected = homer_affected + 1
            elif indicator == 2:
                #removed
                homer_removed = homer_removed + 1
                homer_affected = homer_affected + 1
            elif indicator == 3:
                homer_affected = homer_affected + 1

            #JASPAR
            indicator = motif_search_jaspar(jaspar_motifs, wildtype, mutant, x)
            if indicator == 0:
                #do nothing
                p = 1
            elif indicator == 1:
                #added
                jaspar_added = jaspar_added + 1
                jaspar_affected = jaspar_affected + 1
            elif indicator == 2:
                #removed
                jaspar_removed = jaspar_removed + 1
                jaspar_affected = jaspar_affected + 1
            elif indicator == 3:
                jaspar_affected = jaspar_affected + 1

            #ATG
            indicator = analyze_atg(wildtype, mutant, x)
            if indicator == 0:
                #do nothing
                p = 1
            elif indicator == 1:
                atgs_added = atgs_added + 1
            elif indicator == 2:
                atgs_removed = atgs_removed + 1

            #CTG
            indicator = analyze_ctg(wildtype, mutant, x)
            if indicator == 0:
                #do nothing
                p = 1
            elif indicator == 1:
                ctgs_added = ctgs_added + 1
            elif indicator == 2:
                ctgs_removed = ctgs_removed + 1

            #KOZAK
            indicator = analyze_kozak(wildtype, mutant, x)
            if indicator == 0:
                #do nothing
                p = 1
            elif indicator == 1:
                kozaks_enhanced = kozaks_enhanced + 1
            elif indicator == 2:
                kozaks_weakened = kozaks_weakened + 1

            #G QUADRUPLEX
            indicator = gquad(wildtype, mutant, x)
            if indicator == 0:
                #do nothing
                p = 1
            elif indicator == 1:
                gquads = gquads + 1

            #PRTE
            indicator = motif_search_prte(prte_motif, wildtype, mutant, x)
            if indicator == 0:
                #do nothing
                p = 1
            elif indicator == 1:
                prtes_added = prtes_added + 1
                prtes_affected = prtes_affected + 1
            elif indicator == 2:
                prtes_removed = prtes_removed + 1
                prtes_affected = prtes_affected + 1
            elif indicator == 3:
                prtes_affected = prtes_affected + 1

            #FIVE PRIME TOPS
            indicator = fivetop(wildtype, mutant, x)
            if indicator == 0:
                #do nothing
                p = 1
            elif indicator == 1:
                fivetops = fivetops + 1

output_file.write(str(cisbp_added) + "\n" +
str(cisbp_removed) + "\n" +
str(cisbp_affected) + "\n" +

str(rbpdb_added) + "\n" +
str(rbpdb_removed) + "\n" +
str(rbpdb_affected) + "\n" +

str(homer_added) + "\n" +
str(homer_removed) + "\n" +
str(homer_affected) + "\n" +

str(jaspar_added) + "\n" +
str(jaspar_removed) + "\n" +
str(jaspar_affected) + "\n" +

str(hughes_added) + "\n" +
str(hughes_removed) + "\n" +
str(hughes_affected) + "\n" +

str(atgs_added) + "\n" +
str(atgs_removed) + "\n" +

str(ctgs_added) + "\n" +
str(ctgs_removed) + "\n" +


str(kozaks_enhanced) + "\n" +
str(kozaks_weakened) + "\n" +

str(gquads) + "\n" +

str(prtes_added) + "\n" +
str(prtes_removed) + "\n" +
str(prtes_affected) + "\n" +

str(fivetops) + "\n" +

str(total_mutations) + "\n" ) 


output_file.close()