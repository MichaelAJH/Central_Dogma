from LectureCentralDogma import *

# pKa1= α-carboxyl group, pKa2 = α-ammonium ion, and pKa3 = side chain group
# "amino acid": (pKa1, pKa2, pKa3)
pKs = {
    "Gly":(2.34, 9.60, None),
    "Ala":(2.34, 9.69, None),
    "Val":(2.32, 9.62, None),
    "Leu":(2.36, 9.60, None),
    "Ile":(2.36, 9.60, None),
    "Met":(2.28, 9.21, None),
    "Pro":(1.99, 10.60, None),
    "Phe":(1.83, 9.13, None),
    "Trp":(2.83, 9.39, None),
    "Asp":(2.02, 8.80, None),
    "Gln":(2.17, 9.13, None),
    "Ser":(2.21, 9.15, None),
    "Thr":(2.09, 9.10, None),
    "Tyr":(2.20, 9.11, None),
    "Cys":(1.96, 8.18, None), 
    "Asp":(1.88, 9.60, 3.65),
    "Glu":(2.19, 9.67, 4.25),
    "Lys":(2.18, 8.95, 10.53),
    "Arg":(2.17, 9.04, 12.48),
    "His":(1.82, 9.17, 6.00)
}

def getpH(polypeptide):

    if "Lys" in polypeptide or "Arg" in polypeptide or "His" in polypeptide or "Asp" in polypeptide or "Glu" in polypeptide:
        charge = 0
        pK = []
        pKa_N = pKs["Met"][1]
        pKa_C = pKs[polypeptide[-1]][0]
        pK.extend([pKa_C, pKa_N])
        for amino in polypeptide:
            if amino=="Lys" or amino=="Arg" or amino=="His":
                charge += 1
                pK.append(pKs[amino][2])
            elif amino=="Asp" or amino=="Glu":
                charge -= 1
                pK.append(pKs[amino][2])
            
    else:
        pKa_N = pKs["Met"][1]
        pKa_C = pKs[polypeptide[-1]][0]
        pI = (pKa_C + pKa_N)/2
    return pI