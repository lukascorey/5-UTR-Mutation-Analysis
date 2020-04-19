def revcomp(sequence):
    new = []
    for var in range(len(sequence)):
        if sequence[((len(sequence)) - 1 - var)] == "A":
            new.append("T")
        elif sequence[((len(sequence)) - 1 - var)] == "C":
            new.append("G")
        elif sequence[((len(sequence)) - 1 - var)] == "G":
            new.append("C")
        elif sequence[((len(sequence)) - 1 - var)] == "T" or sequence[((len(sequence)) - 1 -  var)] == "U":
            new.append("A")
        else:
            print("error")
            quit()
    return "".join(new)

def convert(letter):
    if (letter == "A"):
        return 1
    elif (letter == "C"):
        return 2
    elif (letter == "T"):
        return 3
    elif (letter == "G"):
        return 4
    elif (letter == "_"):
        return 0
    else: 
        print("error, __" + letter + "__ passed as letter")

def convert_back(number):
    if (number == 1):
        return "A"
    elif (number == 2):
        return "C"
    elif (number == 3):
        return "G"
    elif (number == 4):
        return "T"
    elif (number == 0):
        return "_"
    else: 
        print("error, __" + number + "__ passed as number")


def triplet_to_num(triplet, new_nuc):
	return ((convert(new_nuc) * 125) + (convert(triplet[0]) * 25) + (convert(triplet[1]) * 5) + (convert(triplet[2])))

def num_to_triplet(number): 
    returnv = "____"
    returnv_list = list(returnv)
    counter = 3
    while (number > 0):
        returnv_list[counter] = convert_back(number % 5)
        counter = counter - 1
        number = number // 5

    return "".join(returnv_list)
	

def motif_search_topstrand(motif_lines, wildtype, mutant, x):
    motif_in_wt = 0
    motif_in_mut = 0
    for line in range(len(motif_lines) // 6):
        marker = line * 6
        #print(line)
        #print(marker)
        name = motif_lines[(marker)].split()[0][1 : (len(motif_lines[(marker)].split()[0]))]
        length_of_motif = len(motif_lines[2 + (marker)].split())

        As = (motif_lines[marker + 2].split())
        #print(As)
        Cs = (motif_lines[marker + 3].split())
        #print(Cs)
        Gs = (motif_lines[marker + 4].split())
        #print(Gs)
        Ts = (motif_lines[marker + 5].split())
        #print(Ts)

        maximum = float(motif_lines[marker + 1])
        #print(maximum)


        for counter in range(length_of_motif):
            if motif_in_wt == 0:
                if (x - length_of_motif + counter + 1) >= 0 and (x + counter + 1) < len(wildtype) + 1:
                    eval_seq = wildtype[(x - length_of_motif + counter + 1) : (x + counter + 1)]
                    score = 0
                    for i in range(length_of_motif):
                        if eval_seq[i] == "A":
                            score = score + float(As[i])
                        elif eval_seq[i] == "C":
                            score = score + float(Cs[i])
                        elif eval_seq[i] == "G":
                            score = score + float(Gs[i])
                        elif eval_seq[i] == "T" or eval_seq[i] == "U":
                            score = score + float(Ts[i])
                        else:
                            print("error in scoring, terminating program")
                            quit()
                    if score >= .90 * maximum:
                        motif_in_wt = 1
        for counter in range(length_of_motif):
            if motif_in_mut == 0:
                if (x - length_of_motif + counter + 1) >= 0 and (x + counter + 1) < len(mutant) + 1:
                    eval_seq = mutant[(x - length_of_motif + counter + 1) : (x + counter + 1)]
                    score = 0
                    for i in range(length_of_motif):
                        if eval_seq[i] == "A":
                            score = score + float(As[i])
                        elif eval_seq[i] == "C":
                            score = score + float(Cs[i])
                        elif eval_seq[i] == "G":
                            score = score + float(Gs[i])
                        elif eval_seq[i] == "T" or eval_seq[i] == "U":
                            score = score + float(Ts[i])
                        else:
                            print("error in scoring, terminating program")
                            quit()
                    if score >= .90 * maximum:
                        motif_in_mut = 1
    if motif_in_mut + motif_in_wt == 2:
        #print(As + Cs + Gs + Ts)
        #print(wildtype[x-5:x+5])
        #print(mutant[x-5:x+5])
        return 3
    elif motif_in_mut + motif_in_wt == 0:
        return 0
    elif motif_in_mut == 1:
        return 1
    else:
        return 2


def motif_search_bothstrands(motif_lines, wildtype, mutant, x):
    #print("doing motif search on bothstrands")
    #print(len(motif_lines) // 5)
    motif_in_wt = 0
    motif_in_mut = 0
    for line in range(len(motif_lines) // 6):
        marker = line * 6
        #print(line)
        #print(marker)
        name = motif_lines[(marker)].split()[0][1 : (len(motif_lines[(marker)].split()[0]))]
        length_of_motif = len(motif_lines[2 + (marker)].split())

        As = (motif_lines[marker + 2].split())
        #print(As)
        Cs = (motif_lines[marker + 3].split())
        #print(Cs)
        Gs = (motif_lines[marker + 4].split())
        #print(Gs)
        Ts = (motif_lines[marker + 5].split())
        #print(Ts)

        maximum = float(motif_lines[marker + 1])
        #print("length of motif" + str(length_of_motif))
        #print(maximum)

        for counter in range(length_of_motif):
            if motif_in_wt == 0:
                if (x - length_of_motif + counter + 1) >= 0 and (x + counter + 1) < len(wildtype) + 1:
                    eval_seq = wildtype[(x - length_of_motif + counter + 1) : (x + counter + 1)]
                    score = 0
                    for i in range(length_of_motif):
                        if eval_seq[i] == "A":
                            score = score + float(As[i])
                        elif eval_seq[i] == "C":
                            score = score + float(Cs[i])
                        elif eval_seq[i] == "G":
                            score = score + float(Gs[i])
                        elif eval_seq[i] == "T" or eval_seq[i] == "U":
                            score = score + float(Ts[i])
                        else:
                            print("error in scoring, terminating program")
                            quit()
                    #print(score)
                    if score >= .90 * maximum:
                        motif_in_wt = 1
        for counter in range(length_of_motif):
            if motif_in_mut == 0:
                if (x - length_of_motif + counter + 1) >= 0 and (x + counter + 1) < len(mutant) + 1:
                    eval_seq = mutant[(x - length_of_motif + counter + 1) : (x + counter + 1)]
                    score = 0
                    for i in range(length_of_motif):
                        if eval_seq[i] == "A":
                            score = score + float(As[i])
                        elif eval_seq[i] == "C":
                            score = score + float(Cs[i])
                        elif eval_seq[i] == "G":
                            score = score + float(Gs[i])
                        elif eval_seq[i] == "T" or eval_seq[i] == "U":
                            score = score + float(Ts[i])
                        else:
                            print("error in scoring, terminating program")
                            quit()
                    if score >= .90 * maximum:
                        motif_in_mut = 1
        for counter in range(length_of_motif):
            if motif_in_wt == 0:
                if (x - length_of_motif + counter + 1) >= 0 and (x + counter + 1) < len(wildtype) + 1:
                    eval_seq = revcomp(wildtype[(x - length_of_motif + counter + 1) : (x + counter + 1)])
                    score = 0
                    for i in range(length_of_motif):
                        if eval_seq[i] == "A":
                            score = score + float(As[i])
                        elif eval_seq[i] == "C":
                            score = score + float(Cs[i])
                        elif eval_seq[i] == "G":
                            score = score + float(Gs[i])
                        elif eval_seq[i] == "T" or eval_seq[i] == "U":
                            score = score + float(Ts[i])
                        else:
                            print("error in scoring, terminating program")
                            quit()
                    if score >= .90 * maximum:
                        motif_in_wt = 1
        for counter in range(length_of_motif):
            if motif_in_mut == 0:
                if (x - length_of_motif + counter + 1) >= 0 and (x + counter + 1) < len(mutant) + 1:
                    eval_seq = revcomp(mutant[(x - length_of_motif + counter + 1) : (x + counter + 1)])
                    score = 0
                    for i in range(length_of_motif):
                        if eval_seq[i] == "A":
                            score = score + float(As[i])
                        elif eval_seq[i] == "C":
                            score = score + float(Cs[i])
                        elif eval_seq[i] == "G":
                            score = score + float(Gs[i])
                        elif eval_seq[i] == "T" or eval_seq[i] == "U":
                            score = score + float(Ts[i])
                        else:
                            print("error in scoring, terminating program")
                            quit()
                    if score >= .90 * maximum:
                        motif_in_mut = 1

    if motif_in_mut + motif_in_wt == 2:
        return 3
    elif motif_in_mut + motif_in_wt == 0:
        return 0
    elif motif_in_mut == 1:
        return 1
    else:
        return 2

def findmut(wt, mut):
    x = 0
    while wt[x] == mut[x]:
        x = x + 1
    return x

def analyze_atg(wildtype, mutant, x):
    reference = wildtype
    seq = mutant

    if ("ATG" == seq[x : (x + 3)]):
        return 1
    elif ("ATG" == seq[(x - 2) : (x + 1)]):
        return 1
    elif ("ATG" == seq[(x - 1) : (x + 2)]):
        return 1
    elif ("ATG" == reference[x : (x + 3)]):
        return 2
    elif ("ATG" == reference[(x - 2) : (x + 1)]):
        return 2
    elif ("ATG" == reference[(x - 1) : (x + 2)]):
        return 2
    else:
        return 0


def analyze_ctg(wildtype, mutant, x):
    reference = wildtype
    seq = mutant

    if ("CTG" == seq[x : (x + 3)]):
        return 1
    elif ("CTG" == seq[(x - 2) : (x + 1)]):
        return 1
    elif ("CTG" == seq[(x - 1) : (x + 2)]):
        return 1
    elif ("CTG" == reference[x : (x + 3)]):
        return 2
    elif ("CTG" == reference[(x - 2) : (x + 1)]):
        return 2
    elif ("CTG" == reference[(x - 1) : (x + 2)]):
        return 2
    else:
        return 0

def analyze_kozak(wildtype, mutant, x):
    reference = wildtype
    seq = mutant
    if ("ATG" == reference[(x + 3) : (x + 6)]):
        if reference[x] == "C" and seq[x] == "G":
            return 1
        elif reference[x] == "T" and seq[x] == "G":
            return 1
        elif reference[x] == "C" and seq[x] == "A":
            return 1
        elif reference[x] == "T" and seq[x] == "A":
            return 1

        elif reference[x] == "A" and seq[x] == "T":
            return 2
        elif reference[x] == "A" and seq[x] == "C":
            return 2
        elif reference[x] == "G" and seq[x] == "T":
            return 2
        elif reference[x] == "G" and seq[x] == "C":
            return 2
        else:
            return 0

    elif ("ATG" == reference[(x - 3): x]):
        if reference[x] == "C" and seq[x] == "G":
            return 1
        elif reference[x] == "T" and seq[x] == "G":
            return 1
        elif reference[x] == "A" and seq[x] == "G":
            return 1
        elif reference[x] == "G" and seq[x] == "T":
            return 2
        elif reference[x] == "G" and seq[x] == "C":
            return 2
        elif reference[x] == "G" and seq[x] == "A":
            return 2
        else:
            return 0

    elif ("ATG" == reference[(x + 2) : (x + 5)]) or ("ATG" == reference[(x + 1) : (x + 4)]):
        if reference[x] == "T" and seq[x] == "C":
            return 1
        elif reference[x] == "A" and seq[x] == "C":
            return 1
        elif reference[x] == "G" and seq[x] == "C":
            return 1

        elif reference[x] == "C" and seq[x] == "A":
            return 2
        elif reference[x] == "C" and seq[x] == "G":
            return 2
        elif reference[x] == "C" and seq[x] == "T":
            return 2
        else:
            return 0


    elif ("CTG" == reference[(x + 3) : (x + 6)]):
        if reference[x] == "C" and seq[x] == "G":
            return 1
        elif reference[x] == "T" and seq[x] == "G":
            return 1
        elif reference[x] == "C" and seq[x] == "A":
            return 1
        elif reference[x] == "T" and seq[x] == "A":
            return 1
        elif reference[x] == "A" and seq[x] == "T":
            return 2
        elif reference[x] == "A" and seq[x] == "C":
            return 2
        elif reference[x] == "G" and seq[x] == "T":
            return 2
        elif reference[x] == "G" and seq[x] == "C":
            return 2
        else:
            return 0

    elif ("CTG" == reference[(x - 3): x]):
        if reference[x] == "C" and seq[x] == "G":
            return 1
        elif reference[x] == "T" and seq[x] == "G":
            return 1
        elif reference[x] == "A" and seq[x] == "G":
            return 1
        elif reference[x] == "G" and seq[x] == "T":
            return 2
        elif reference[x] == "G" and seq[x] == "C":
            return 2
        elif reference[x] == "G" and seq[x] == "A":
            return 2
        else:
            return 0

    elif ("CTG" == reference[(x + 2) : (x + 5)]) or ("CTG" == reference[(x + 1) : (x + 4)]):
        if reference[x] == "T" and seq[x] == "C":
            return 1
        elif reference[x] == "A" and seq[x] == "C":
            return 1
        elif reference[x] == "G" and seq[x] == "C":
            return 1

        elif reference[x] == "C" and seq[x] == "A":
            return 2
        elif reference[x] == "C" and seq[x] == "G":
            return 2
        elif reference[x] == "C" and seq[x] == "T":
            return 2
        else:
            return 0
    else:
        return 0

def gquad(wildtype, mutant, x):
    present = 0
    reference = wildtype
    seq = mutant
    stop = len(reference) - 1
    for i in range(29):
        if not ((x - 29 + i + 1) > stop) and ((x - 29 + i + 1) >= 0):
            if reference[x - 29 + i] == "G" and reference[x - 29 + i + 1] == "G":
                for j in range(8):
                    if not ((x - 29 + i + 4 + j) > stop) and ((x - 29 + i + 4 + j) >= 0):
                        if reference[x - 29 + i + 3 + j] == "G" and reference[x - 29 + i + 4 + j] == "G":
                            for k in range(8):
                                if not ((x - 29 + i + 4 + j + 3 + k) > stop) and ((x - 29 + i + 4 + j + 3 + k) >= 0):
                                    if reference[x - 29 + i + 4 + j + 2 + k] == "G" and reference[x - 29 + i + 4 + j + 3 + k] == "G":
                                        for l in range(8):
                                            if not ((x - 29 + i + 4 + j + 3 + k + 3 + l) > stop) and ((x - 29 + i + 4 + j + 3 + k + 4 + l) > 0):
                                                if reference[x - 29 + i + 4 + j + 3 + k + 2 + l] == "G" and reference[x - 29 + i + 4 + j + 3 + k + 3 + l] == "G":
                                                    present = 1
    for i in range(29):
        if not ((x - 29 + i + 1) > stop) and ((x - 29 + i) >= 0):
            if seq[x - 29 + i] == "G" and seq[x - 29 + i + 1] == "G":
                for j in range(8):
                    if not ((x - 29 + i + 4 + j) > stop) and ((x - 29 + i + 4 + j) >= 0):
                        if seq[x - 29 + i + 3 + j] == "G" and seq[x - 29 + i + 4 + j] == "G":
                            for k in range(8):
                                if not ((x - 29 + i + 4 + j + 3 + k) > stop) and ((x - 29 + i + 4 + j + 3 + k) >= 0):
                                    if seq[x - 29 + i + 4 + j + 2 + k] == "G" and seq[x - 29 + i + 4 + j + 3 + k] == "G":
                                        for l in range(8):
                                            if not ((x - 29 + i + 4 + j + 3 + k + 3 + l) > stop) and ((x - 29 + i + 4 + j + 3 + k + 4 + l) > 0):
                                                if seq[x - 29 + i + 4 + j + 3 + k + 2 + l] == "G" and seq[x - 29 + i + 4 + j + 3 + k + 3 + l] == "G":
                                                    present = 1
    return present

def fivetop(wildtype, mutant, x):
    present1 = 0
    present2 = 0
    reference = wildtype
    seq = mutant


    if reference[0] == "C" or seq[0] == "C":

        for i in range(4):
            if ((1 + i) >= len(reference) - 1) or ((1 + i) > len(seq) - 1):
                present1 = 1
                present2 = 1

            else:
                if reference[1 + i] == "A" or reference[i + 1] == "G":
                    present1 = 1
                if seq[i + 1] == "A" or seq[1 + i] == "G":
                    present2 = 1

    else:
        present1 = 1
        present2 = 1
    if (present1 == 0 or present2 == 0) and x < 15:
        return 1
    else:
        return 0


def motif_search_prte(motif_lines, wildtype, mutant, x):
    motif_in_wt = 0
    motif_in_mut = 0
    for line in range(len(motif_lines) // 6):
        marker = line * 6
        #print(line)
        #print(marker)
        name = motif_lines[(marker)].split()[0][1 : (len(motif_lines[(marker)].split()[0]))]
        length_of_motif = len(motif_lines[2 + (marker)].split())

        As = (motif_lines[marker + 2].split())
        #print(As)
        Cs = (motif_lines[marker + 3].split())
        #print(Cs)
        Gs = (motif_lines[marker + 4].split())
        #print(Gs)
        Ts = (motif_lines[marker + 5].split())
        #print(Ts)

        maximum = float(motif_lines[marker + 1])


        for counter in range(length_of_motif):
            if motif_in_wt == 0:
                if (x - length_of_motif + counter + 1) >= 0 and (x + counter + 1) < len(wildtype) + 1:
                    eval_seq = wildtype[(x - length_of_motif + counter + 1) : (x + counter + 1)]
                    score = 0
                    for i in range(length_of_motif):
                        if eval_seq[i] == "A":
                            score = score + float(As[i])
                        elif eval_seq[i] == "C":
                            score = score + float(Cs[i])
                        elif eval_seq[i] == "G":
                            score = score + float(Gs[i])
                        elif eval_seq[i] == "T" or eval_seq[i] == "U":
                            score = score + float(Ts[i])
                        else:
                            print("error in scoring, terminating program")
                            quit()
                    if score >= .70 * maximum:
                        motif_in_wt = 1
        for counter in range(length_of_motif):
            if motif_in_mut == 0:
                if (x - length_of_motif + counter + 1) >= 0 and (x + counter + 1) < len(mutant) + 1:
                    eval_seq = mutant[(x - length_of_motif + counter + 1) : (x + counter + 1)]
                    score = 0
                    for i in range(length_of_motif):
                        if eval_seq[i] == "A":
                            score = score + float(As[i])
                        elif eval_seq[i] == "C":
                            score = score + float(Cs[i])
                        elif eval_seq[i] == "G":
                            score = score + float(Gs[i])
                        elif eval_seq[i] == "T" or eval_seq[i] == "U":
                            score = score + float(Ts[i])
                        else:
                            print("error in scoring, terminating program")
                            quit()
                    if score >= .70 * maximum:
                        motif_in_mut = 1
    if motif_in_mut + motif_in_wt == 2:
        #print(As + Cs + Gs + Ts)
        #print(wildtype[x-5:x+5])
        #print(mutant[x-5:x+5])
        return 3
    elif motif_in_mut + motif_in_wt == 0:
        return 0
    elif motif_in_mut == 1:
        return 1
    else:
        return 2

def getrefs(UTR_lines):
    refs = []
    for i in UTR_lines:
        if i.split()[2] == "ref":
            refs.append(i)
            refs.append("\n")
    return ("".join(refs)).splitlines()

def motif_search_jaspar(motif_lines, wildtype, mutant, x):
    #print("doing motif search on bothstrands")
    #print(len(motif_lines) // 5)
    motif_in_wt = 0
    motif_in_mut = 0
    for line in range(len(motif_lines) // 6):
        marker = line * 6
        #print(line)
        #print(marker)
        name = motif_lines[(marker)].split()[0][1 : (len(motif_lines[(marker)].split()[0]))]
        length_of_motif = len(motif_lines[2 + (marker)].split())

        As = (motif_lines[marker + 2].split())
        #print(As)
        Cs = (motif_lines[marker + 3].split())
        #print(Cs)
        Gs = (motif_lines[marker + 4].split())
        #print(Gs)
        Ts = (motif_lines[marker + 5].split())
        #print(Ts)

        maximum = float(motif_lines[marker + 1])
        #print("length of motif" + str(length_of_motif))
        #print(maximum)

        for counter in range(length_of_motif):
            if motif_in_wt == 0:
                if (x - length_of_motif + counter + 1) >= 0 and (x + counter + 1) < len(wildtype) + 1:
                    eval_seq = wildtype[(x - length_of_motif + counter + 1) : (x + counter + 1)]
                    score = 0
                    for i in range(length_of_motif):
                        if eval_seq[i] == "A":
                            score = score + float(As[i])
                        elif eval_seq[i] == "C":
                            score = score + float(Cs[i])
                        elif eval_seq[i] == "G":
                            score = score + float(Gs[i])
                        elif eval_seq[i] == "T" or eval_seq[i] == "U":
                            score = score + float(Ts[i])
                        else:
                            print("error in scoring, terminating program")
                            quit()
                    #print(score)
                    if score >= .65 * maximum:
                        motif_in_wt = 1
        for counter in range(length_of_motif):
            if motif_in_mut == 0:
                if (x - length_of_motif + counter + 1) >= 0 and (x + counter + 1) < len(mutant) + 1:
                    eval_seq = mutant[(x - length_of_motif + counter + 1) : (x + counter + 1)]
                    score = 0
                    for i in range(length_of_motif):
                        if eval_seq[i] == "A":
                            score = score + float(As[i])
                        elif eval_seq[i] == "C":
                            score = score + float(Cs[i])
                        elif eval_seq[i] == "G":
                            score = score + float(Gs[i])
                        elif eval_seq[i] == "T" or eval_seq[i] == "U":
                            score = score + float(Ts[i])
                        else:
                            print("error in scoring, terminating program")
                            quit()
                    if score >= .65 * maximum:
                        motif_in_mut = 1
        for counter in range(length_of_motif):
            if motif_in_wt == 0:
                if (x - length_of_motif + counter + 1) >= 0 and (x + counter + 1) < len(wildtype) + 1:
                    eval_seq = revcomp(wildtype[(x - length_of_motif + counter + 1) : (x + counter + 1)])
                    score = 0
                    for i in range(length_of_motif):
                        if eval_seq[i] == "A":
                            score = score + float(As[i])
                        elif eval_seq[i] == "C":
                            score = score + float(Cs[i])
                        elif eval_seq[i] == "G":
                            score = score + float(Gs[i])
                        elif eval_seq[i] == "T" or eval_seq[i] == "U":
                            score = score + float(Ts[i])
                        else:
                            print("error in scoring, terminating program")
                            quit()
                    if score >= .65 * maximum:
                        motif_in_wt = 1
        for counter in range(length_of_motif):
            if motif_in_mut == 0:
                if (x - length_of_motif + counter + 1) >= 0 and (x + counter + 1) < len(mutant) + 1:
                    eval_seq = revcomp(mutant[(x - length_of_motif + counter + 1) : (x + counter + 1)])
                    score = 0
                    for i in range(length_of_motif):
                        if eval_seq[i] == "A":
                            score = score + float(As[i])
                        elif eval_seq[i] == "C":
                            score = score + float(Cs[i])
                        elif eval_seq[i] == "G":
                            score = score + float(Gs[i])
                        elif eval_seq[i] == "T" or eval_seq[i] == "U":
                            score = score + float(Ts[i])
                        else:
                            print("error in scoring, terminating program")
                            quit()
                    if score >= .65 * maximum:
                        motif_in_mut = 1
    if motif_in_mut + motif_in_wt == 2:
        return 3
    elif motif_in_mut + motif_in_wt == 0:
        return 0
    elif motif_in_mut == 1:
        return 1
    else:
        return 2
