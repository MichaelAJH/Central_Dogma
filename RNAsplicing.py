# input: original RNA starting from the 5'-
# output: processed RNA starting from the 5'-
def RNAsplicing(rna_pre):
    exon_count = int(input("엑손의 수 입력: \t"))
    exon_info = list()
    
    for i in range(exon_count):
        exon_info.append((int(input(f"엑손{i+1}의 시작 위치: \t")),int(input(f"엑손{i+1}의 끝 위치: \t"))))
    
    rna_result = ""
    for i in range(exon_count):
        rna_result += rna_pre[exon_info[i][0]-1:exon_info[i][1]]
    
    return rna_result


if __name__ == "__main__":
    print(RNAsplicing(input("RNA 입력: ")))