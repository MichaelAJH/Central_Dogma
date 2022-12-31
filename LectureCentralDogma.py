COMP_BASE  = {'A':'T', 'T':'A', 'C':'G', 'G':'C'} 
CODON_DICT = { 
  "UUU":"Phe", "UUC":"Phe", "UUA":"Leu", "UUG":"Leu", 
  "CUU":"Leu", "CUC":"Leu", "CUA":"Leu", "CUG":"Leu", 
  "AUU":"Ile", "AUC":"Ile", "AUA":"Ile", "AUG":"Met", 
  "GUU":"Val", "GUC":"Val", "GUA":"Val", "GUG":"Val", 
  "UCU":"Ser", "UCC":"Ser", "UCA":"Ser", "UCG":"Ser", 
  "CCU":"Pro", "CCC":"Pro", "CCA":"Pro", "CCG":"Pro", 
  "ACU":"Thr", "ACC":"Thr", "ACA":"Thr", "ACG":"Thr", 
  "GCU":"Ala", "GCC":"Ala", "GCA":"Ala", "GCG":"Ala", 
  "UAU":"Tyr", "UAC":"Tyr", "UAA":"Stop", "UAG":"Stop", 
  "CAU":"His", "CAC":"His", "CAA":"Gln", "CAG":"Gln", 
  "AAU":"Asn", "AAC":"Asn", "AAA":"Lys", "AAG":"Lys", 
  "GAU":"Asp", "GAC":"Asp", "GAA":"Glu", "GAG":"Glu", 
  "UGU":"Cys", "UGC":"Cys", "UGA":"Stop", "UGG":"Trp", 
  "CGU":"Arg", "CGC":"Arg", "CGA":"Arg", "CGG":"Arg", 
  "AGU":"Ser", "AGC":"Ser", "AGA":"Arg", "AGG":"Arg", 
  "GGU":"Gly", "GGC":"Gly", "GGA":"Gly", "GGG":"Gly"}


class CentralDogma(str):
    def __init__(self, seq) -> None:
        super().__init__()
        self.seq = seq
    
    def complement(self):
        # 상보염기서열을 반환하시오.
        comp = ''
        for base in self.seq:
            comp += COMP_BASE[base]
        return comp
    
    def translate(self):
        # 상보염기서열의 순서를 뒤집고, T를 U로 치환하여 mRNA를 반환하시오.
        rev = self.complement()[:-1]
        mRNA = ''
        for base in rev:
            if base == 'T': mRNA += 'U'
            else: mRNA += base
        return mRNA
        '''
        return self.complement()[:-1].replace['T', 'U']
        '''

    def initiate(self, mRNA): 
        # mRNA에서 AUG를 찾아 개시 코돈 앞 부분의 염기 서열을 잘라내시오.
        for i in range(len(mRNA)):
            if i+2<len(mRNA) and mRNA[i] == 'A' and mRNA[i+1] == 'U' and mRNA[i+2] == 'G':
                mRNA = mRNA[i:]
        return mRNA[:(len(mRNA)//3)*3]  # 0~10까지 -> 012, 345, 678 -> [0:9]
        '''
        mRNA[self.translate().find('AUG'):]
        return mRNA[:(len(mRNA)//3)+3]
        '''
        
    def elongate(self, mRNA):
        # mRNA를 번역하여 아미노산 리스트로 반환하시오.
        amino_acid = []
        for i in range(len(mRNA)//3):
            codon = mRNA[i*3] + mRNA[i*3+1] + mRNA[i*3+2]
            amino_acid.append(CODON_DICT[codon])
        return amino_acid
        '''
        return [CODON_DICT[mRNA[i:i+3]] for i in range(0,len(mRNA),3)]
        '''

    def terminate(self):
        codonList = self.initiate(self.translate())
        aminoAcid = self.elongate(codonList)
        # 첫 번째 'Stop' 뒤의 아미노산을 삭제하여 반환하시오.
        for amino in aminoAcid:
            if amino == 'Stop': aminoAcid = aminoAcid[:aminoAcid.index(amino)+1]
        return aminoAcid
        '''
        if Stop in aminoAcid:
            aminoAcid = aminoAcid[:aminoAcid.index('Stop')+1]       # index는 오류 존재, find는 오류 X -> 조건문 사용했으므로 이 경우에는 오류 X
        '''
    def moving_after_translate(self):
        aminoAcid_list = self.terminate()
        #N말단에 있는 신호 서열 인식
        if ['Pro', 'Pro', 'Lys', 'Lys', 'Lys', 'Arg', 'Lys', 'Val'] in aminoAcid_list:
            return 'nucleus'
        elif len(aminoAcid_list) >= 39:
            if aminoAcid_list[0:39] ==  ['Met', 'Val', 'Ala', 'Met', 'Ala', 'Met', 'Ala', 'Ser', 'Leu', 'Gln', 'Ser', 'Ser', 'Met', 'Ser', 'Ser', 'Leu', 'Ser', 'Leu', 'Ser', 'Ser', 'Asn', 'Ser', 'Phe', 'Leu', 'Gly', 'Gln', 'Pro', 'Leu', 'Ser', 'Pro', 'Ile', 'Thr', 'Leu', 'Ser', 'Pro', 'Phe', 'Leu', 'Gln', 'Gly']:
                return 'plastid'
        elif len(aminoAcid_list) >= 26:
            if aminoAcid_list[0:26] == ['Met', 'Val', 'Leu', 'Ser', 'Leu', 'Arg', 'Gln', 'Ser', 'Ile', 'Arg', 'Phe', 'Phe', 'Lys', 'Pro', 'Ala', 'Thr', 'Arg', 'Thr', 'Leu', 'Cys', 'Ser', 'Ser', 'Arg', 'Tyr', 'Leu', 'Leu']:
                return 'mitochondria'
        elif len(aminoAcid_list) >= 29:
            if aminoAcid_list[0:29] == ['Met', 'Met', 'Ser', 'Phe', 'Val', 'Ser', 'Leu', 'Leu', 'Leu', 'Val', 'Gly', 'Ile', 'Leu', 'Phe', 'Trp', 'Ala', 'thr', 'Glu', 'Ala', 'Glu', 'Gln', 'Leu', 'Thr', 'Lys', 'Cys', 'Glu', 'Val', 'Phe', 'Gln']:
                return 'ER'
        else:
            return 'cytoplasm'
    def check_protein_location(self):
      if dna_temp.moving_after_translate() == 'ER':
        Protein_Glycosylation = input("단백질의 당화를 하시겠습니까? (Yes or No) : ")
        if Protein_Glycosylation == "Yes" :
          print("단백질이 골지체로 이동 중입니다.")
          Protein_Phosphorylation = input("단백질의 인산화를 하시겠습니까? (Yes or No) : ")
          if Protein_Phosphorylation == "Yes" :
            print("단백질이 리소좀으로 이동 중입니다.")
          else:
            print("단백질이 외부로 분비 되었습니다.")
    
    


temp = 'ATGTTAAAGAGCAGTCACAGACTTTAGCATTG'

dna_temp = CentralDogma(temp)
print("DNA 주형가닥 : ", "5'-", dna_temp , "-3'")
dna_comp = dna_temp.complement()
print("DNA 상보가닥 : ", "3'-", dna_comp , "-5'")
mRNA = dna_temp.translate()
print("mRNA         : ", "5'-", mRNA, "-3'")
print('initated     :', dna_temp.initiate(dna_temp.translate()))
print(dna_temp.terminate())
print('protein is moving to %s' %dna_temp.moving_after_translate())
dna_temp.check_protein_location()
    
