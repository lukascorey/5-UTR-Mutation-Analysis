import os
import argparse
import re
import random
import math
import sys
import time
import os.path
import datetime

total_mutation_list = []

cisbp_added_list = []
cisbp_removed_list = []
cisbp_affected_list = []

rbpdb_added_list = []
rbpdb_removed_list = []
rbpdb_affected_list = []

homer_added_list = []
homer_removed_list = []
homer_affected_list = []

jaspar_added_list = []
jaspar_removed_list = []
jaspar_affected_list = []

hughes_added_list = []
hughes_removed_list = []
hughes_affected_list = []

atgs_added_list = []
atgs_removed_list = []

ctgs_added_list = []
ctgs_removed_list = []

kozaks_enhanced_list = []
kozaks_weakened_list = []

gquads_list = []

prtes_added_list = []
prtes_removed_list = []
prtes_affected_list = []

fivetops_list = []



for i in range(10000):
	
	if (os.path.exists("joint_element_analysis_{0}\r".format(i),)):
		with open("joint_element_analysis_{0}\r".format(i), "r") as datafile:
			print(i)
			data = datafile.read().splitlines()
			datafile.close()
			total_mutation = data[26]
			cisbp_added = data[0]
			cisbp_removed = data[1]
			cisbp_affected = data[2]

			rbpdb_added = data[3]
			rbpdb_removed = data[4]
			rbpdb_affected = data[5]

			homer_added = data[6]
			homer_removed = data[7]
			homer_affected = data[8]

			jaspar_added = data[9]
			jaspar_removed = data[10]
			jaspar_affected = data[11]

			hughes_added = data[12]
			hughes_removed = data[13]
			hughes_affected = data[14]

			atgs_added = data[15]
			atgs_removed = data[16]

			ctgs_added = data[17]
			ctgs_removed = data[18]

			kozaks_enhanced = data[19]
			kozaks_weakened = data[20]

			gquads = data[21]

			prtes_added = data[22]
			prtes_removed = data[23]
			prtes_affected = data[24]

			fivetops = data[25]

			total_mutation_list.append(str(total_mutation) + "\t")

			cisbp_added_list.append(str(cisbp_added) + "\t")
			cisbp_removed_list.append(str(cisbp_removed) + "\t")
			cisbp_affected_list.append(str(cisbp_affected) + "\t")

			rbpdb_added_list.append(str(rbpdb_added) + "\t")
			rbpdb_removed_list.append(str(rbpdb_removed) + "\t")
			rbpdb_affected_list.append(str(rbpdb_affected) + "\t")

			homer_added_list.append(str(homer_added) + "\t")
			homer_removed_list.append(str(homer_removed) + "\t")
			homer_affected_list.append(str(homer_affected) + "\t")

			jaspar_added_list.append(str(jaspar_added) + "\t")
			jaspar_removed_list.append(str(jaspar_removed) + "\t")
			jaspar_affected_list.append(str(jaspar_affected) + "\t")

			hughes_added_list.append(str(hughes_added) + "\t")
			hughes_removed_list.append(str(hughes_removed) + "\t")
			hughes_affected_list.append(str(hughes_affected) + "\t")

			atgs_added_list.append(str(atgs_added) + "\t")
			atgs_removed_list.append(str(atgs_removed) + "\t")

			ctgs_added_list.append(str(ctgs_added) + "\t")
			ctgs_removed_list.append(str(ctgs_removed) + "\t")

			kozaks_enhanced_list.append(str(kozaks_enhanced) + "\t")
			kozaks_weakened_list.append(str(kozaks_weakened) + "\t")

			gquads_list.append(str(gquads) + "\t")

			prtes_added_list.append(str(prtes_added) + "\t")
			prtes_removed_list.append(str(prtes_removed) + "\t")
			prtes_affected_list.append(str(prtes_affected) + "\t")

			fivetops_list.append(str(fivetops) + "\t")
	

output_file  = open("combined_data2.txt", "w+")

output_file.write(
	"total_mutation_list: " + "\n" + str("".join(total_mutation_list)) + "\n" +

    "cisbp_added_list: "  + "\n" + str("".join(cisbp_added_list)) + "\n" +
    "cisbp_removed_list: " + "\n" + str("".join(cisbp_removed_list)) + "\n" +
    "cisbp_affected_list: " + "\n" + str("".join(cisbp_affected_list)) + "\n" +

    "rbpdb_added_list: "  + "\n" + str("".join(rbpdb_added_list)) + "\n" +
    "rbpdb_removed_list: " + "\n" + str("".join(rbpdb_removed_list)) + "\n" +
    "rbpdb_affected_list: " + "\n" + str("".join(rbpdb_affected_list)) + "\n" +

    "homer_added_list: "  + "\n" + str("".join(homer_added_list)) + "\n" +
    "homer_removed_list: " + "\n" + str("".join(homer_removed_list)) + "\n" +
    "homer_affected_list: " + "\n" + str("".join(homer_affected_list)) + "\n" +

    "jaspar_added_list: "  + "\n" + str("".join(jaspar_added_list)) + "\n" +
    "jaspar_removed_list: " + "\n" + str("".join(jaspar_removed_list)) + "\n" +
    "jaspar_affected_list: " + "\n" + str("".join(jaspar_affected_list)) + "\n" +

    "hughes_added_list: "  + "\n" + str("".join(hughes_added_list)) + "\n" +
    "hughes_removed_list: " + "\n" + str("".join(hughes_removed_list)) + "\n" +
    "hughes_affected_list: " + "\n" + str("".join(hughes_affected_list)) + "\n" +


    "atgs_added_list: " + "\n" + str("".join(atgs_added_list)) + "\n" +
    "atgs_removed_list: " + "\n" + str("".join(atgs_removed_list)) + "\n" +

    "ctgs_added_list: " + "\n" + str("".join(ctgs_added_list)) + "\n" +
    "ctgs_removed_list: " + "\n" + str("".join(ctgs_removed_list)) + "\n" +

    "kozaks_enhanced_list: " + "\n" + str("".join(kozaks_enhanced_list)) + "\n" +
    "kozaks_weakened_list: " + "\n" + str("".join(kozaks_weakened_list)) + "\n" +

    "gquads_list: " + "\n" + str("".join(gquads_list)) + "\n" +

    "prtes_added_list: "  + "\n" + str("".join(prtes_added_list)) + "\n" +
    "prtes_removed_list: " + "\n" + str("".join(prtes_removed_list)) + "\n" +
    "prtes_affected_list: " + "\n" + str("".join(prtes_affected_list)) + "\n" +

    "fivetops_list: " + "\n" + str("".join(fivetops_list)))

output_file.close()