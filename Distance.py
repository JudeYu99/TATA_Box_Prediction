#! /usr/bin/python

# 获取注释信息的列表结果
with open('Sc_gene.gff', 'r+') as genes:
    gene_loc_info = []
    for line in genes.readlines(): 
        temp =line.split(" ")
        temp[-1] = temp[-1].replace('\n', "")
        gene_loc_info.append(temp)


# 获取TATA box信息的列表结果
with open('TATA_result_p_0.01.gff', 'r+') as boxes:
    TATA_loc_info = []
    for line in boxes.readlines():
        temp =line.split("\t")
        temp[-1] = temp[-1].replace('\n', "")
        TATA_loc_info.append(temp)


# 预测结果加上终止位置
for item in TATA_loc_info:
    item.insert(2, str(int(item[1]) + 12))

    
# 按顺序获取所有染色体
all_chr = []
[all_chr.append(item[0]) for item in gene_loc_info if item[0] not in all_chr]


# 返回第n条染色体的信息
def chr_info(info, chrID):
    chr_info = []
    for item in info:
        if item[0] == all_chr[chrID]:
            chr_info.append(item)
    return chr_info


# 提取并区分同一个染色体的正负链，返回正链信息和负链信息
def ExtractPositiveNegativeChain(chr_info):
    positiveChain = []
    negativeChain = []
    for item in chr_info:
        if item[-1] == "+":
            positiveChain.append(item)
        else:
            negativeChain.append(item)
    return (positiveChain, negativeChain)


# 提取同一个染色体的正负链的起始位置，分别返回正负链的起始信息的列表
def Extract_Ini_pos(chr_info):
    positiveChainIni = []
    negativeChainIni = []
    temp = ExtractPositiveNegativeChain(chr_info)
    for item in temp[0]:
        positiveChainIni.append(int(item[1]))
    for item in temp[1]:
        negativeChainIni.append(int(item[1]))

    return (positiveChainIni, negativeChainIni)


'''
# 计算预测的TATA box距离下游最近的TSS的距离
def distance(TATA_loc, gene_loc):
    distance = []
    for i in TATA_loc:
        TATA_end = i + 11
        if TATA_end < gene_loc[-1]:
            gene_loc.append(TATA_end)
            gene_loc.sort()

            TATA_end_Ind = gene_loc.index(TATA_end)

            if TATA_end == gene_loc[TATA_end_Ind + 1]:
                k = 2
                while gene_loc[TATA_end_Ind + k] == TATA_end:
                    k = k + 1
                else:
                    distance.append(gene_loc[TATA_end_Ind + k] - i)
            else:
                distance.append(gene_loc[TATA_end_Ind + 1] - i)

    return distance
'''

# 正链计算距离; TATA_loc, gene_loc分别为起始位置
def positive_distance(TATA_loc, gene_loc):
    distance = []
    for i in TATA_loc:
        TATA_end = i + 11
        if TATA_end < gene_loc[-1]:
            gene_loc.append(TATA_end)
            gene_loc.sort()

            TATA_end_Ind = gene_loc.index(TATA_end)

            if TATA_end == gene_loc[TATA_end_Ind + 1]:
                k = 2
                while gene_loc[TATA_end_Ind + k] == TATA_end:
                    k = k + 1
                else:
                    distance.append(gene_loc[TATA_end_Ind + k] - i)
            else:
                distance.append(gene_loc[TATA_end_Ind + 1] - i)

    return distance

# 负链计算距离; TATA_loc, gene_loc分别为终止位置
def negative_distance(TATA_loc, gene_loc):
    distance = []
    for i in TATA_loc:
        TATA_begin = i - 11
        if TATA_begin > gene_loc[0]:
            gene_loc.append(TATA_begin)
            gene_loc.sort()

            TATA_begin_Ind = gene_loc.index(TATA_begin)
            
            TATA_begin_Ind = gene_loc.index(TATA_begin)
            distance.append(i - gene_loc[TATA_begin_Ind - 1])
    return distance


def TATA_positive(TATA):
    res = []
    for item in TATA:
        if item[-1] == "+":
            res.append(item)
    return res

def TATA_negative(TATA):
    res = []
    for item in TATA:
        if item[-1] == "-":
            res.append(item)
    return res


def ProcessSingleChr(ChrID):
    GENE = chr_info(gene_loc_info, ChrID)
    TATA = chr_info(TATA_loc_info, ChrID)

    DIS_positive = positive_distance(Extract_Ini_pos(TATA)[0], Extract_Ini_pos(GENE)[0])
    DIS_negative = negative_distance(Extract_Ini_pos(TATA)[1], Extract_Ini_pos(GENE)[1])

    positiveChain =TATA_positive(TATA)
    negativeChain = TATA_negative(TATA)

    valid_positive = positiveChain[:len(DIS_positive)]
    valid_negative = negativeChain[:len(DIS_negative)]

    for i in range(len(DIS_positive)):
        valid_positive[i].append(DIS_positive[i])   

    for j in range(len(DIS_negative)):
        valid_negative[j].append(DIS_negative[j])
    
    return (valid_positive, valid_negative) 


# 结果写入文件
def write_file(filename):
    with open(filename, "a+") as res:
        temp = []
        for item in valid_positive:
            temp.append(item[0] + "\t" + str(item[1]) + "\t" + str(item[2]) + "\t" + str(item[3]) + "\t" + item[4] + "\t" + str(item[5])+ "\t" + str(item[6]) + "\n")
        res.writelines(temp)

    with open(filename, "a+") as res:
        temp = []
        for item in valid_negative:
            temp.append(item[0] + "\t" + str(item[1]) + "\t" + str(item[2]) + "\t" + str(item[3]) + "\t" + item[4] + "\t" + str(item[5]) + "\t" + str(item[6]) + "\n")
        res.writelines(temp)


for i in range(len(all_chr)):
    singleChr = ProcessSingleChr(i)
    valid_positive = singleChr[0]
    valid_negative = singleChr[1]
    write_file("TATA_distance.gff")
