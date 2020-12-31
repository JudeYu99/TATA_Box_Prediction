#! /usr/bin/python

import math

# 调用处理HMM矩阵
hmm_900={}
with open('TATA_HMM.txt', 'r') as hmm:
    for line in hmm.readlines(): 
        #print(line)
        if line[0] == '%':
            percent = line.split('\t') 
            percent[-1] = percent[-1].replace('\n', "")
            #print(percent)
            for i in range(1, 13):
                hmm_900[str(i) + '_' + percent[0][-1]] = float(percent[i])/100

# 读取fna文件
def read_file(filename, filetype):
    from Bio import SeqIO
    records = list(SeqIO.parse(filename, filetype))
    return records

# 12个碱基的移动框打分(不跳过0的情况)
def hmm_score(seq):
    score = 1
    for i in range(len(seq)): 
        score *= float(hmm_900[str(i + 1) + "_" + seq[i]])
    return score

# 记录起始位点和HMM得分的字典
def ORIGIN_SCORE(seq):
    seq = seq.upper()

    ini_pos = []
    scores = []
    for i in range(len(seq) - 11):
        temp = seq[i:i+12]
        if len([nt for nt in temp if nt not in "AGCT"]) != 0 or hmm_score(temp) == 0:
            continue
        else:
            scores.append(hmm_score(temp))
            ini_pos.append(i+1)
    result = dict(zip(ini_pos, scores))
    
    return result

# 随机打乱
def RandomShuffle(subseq):
    import random
    temp = list(subseq)    
    random.shuffle(temp)
    subseq = "".join(temp)
    return subseq


# bootstrap第一次抽样, N取10，cutoff取0.05
def bootstrap_trial(origin_seq): 
    origin_states = ORIGIN_SCORE(origin_seq) #原始状态的起始位置和得分
    Ini_positions = list(origin_states.keys()) #原始状态的起始位置
    origin_scores = list(origin_states.values()) #原始状态的得分

    positions = [] #bootstrap 抽样评估的方法得到的新的起始位置
    p_values = []
    for j in range(len(Ini_positions)) :
        count = 0
        Ini_pos = Ini_positions[j]
        subseq = origin_seq[Ini_pos-1:Ini_pos+11]
        
        for i in range(10): # 抽样评估总次数为 10
            subseq = RandomShuffle(subseq) #随机打乱
            randomed_score = hmm_score(subseq) #计算随机打乱后的打分

            if randomed_score > origin_scores[j]: #每次评估的片段得分大于原始得分的次数为 n
                count += 1
                positions.append(Ini_pos)
                
        p_value = count/10 #计算p值
        p_values.append(p_value)
    
    result_positions = []
    positions_origin_score = []
    for i in positions: 
        if i not in result_positions:
            result_positions.append(i)
            positions_origin_score.append(origin_states[i])
            
    position_info = [] #result字典的键
    for k in range(len(result_positions)):
        temp = []
        temp.append(p_values[k])
        temp.append(positions_origin_score[k])
        position_info.append((temp))
    
    result = dict(zip(result_positions, position_info)) #result结果键为起始位置，值为p值，原始得分
    out = score_filter(result, 0.05) #只保留p值小于0.05的部分
    return out


# bootstrap抽样
def bootstrap(origin_seq, Ntimes): 
    origin_states = ORIGIN_SCORE(origin_seq) #原始状态的起始位置和得分
    temp = bootstrap_trial(origin_seq)

    Ini_positions = list(temp.keys()) #原始状态的起始位置
    origin_scores = []                 
    for item in list(temp.values()):
        origin_scores.append(item[1])    #原始状态的得分

    positions = [] #bootstrap 第二轮抽样评估的方法得到的新的起始位置
    p_values = []
    for j in range(len(Ini_positions)) :
        count = 0
        Ini_pos = Ini_positions[j]
        subseq = origin_seq[Ini_pos-1:Ini_pos+11]
        
        for i in range(Ntimes): # 抽样评估总次数为 N
            subseq = RandomShuffle(subseq) #随机打乱
            randomed_score = hmm_score(subseq) #计算随机打乱后的打分

            if randomed_score > origin_scores[j]: #每次评估的片段得分大于原始得分的次数为 n
                count += 1
                positions.append(Ini_pos)
                
        p_value = count/Ntimes #计算p值
        p_values.append(p_value)
    
    result_positions = []
    positions_origin_score = []
    for i in positions: 
        if i not in result_positions:
            result_positions.append(i)
            positions_origin_score.append(origin_states[i])
            
    position_info = [] #result字典的键
    for k in range(len(result_positions)):
        temp = []
        temp.append(p_values[k])
        temp.append(-math.log2(positions_origin_score[k]))
        position_info.append((temp))
    
    result = dict(zip(result_positions, position_info)) #result结果键为起始位置，值为p值，原始得分
    return result

def score_filter(bootResult, cutoff):
    dc = {key:value for key, value in bootResult.items() if value[0] < cutoff}
    return dc

# 反向互补序列
def reverse_comp(seq):
    comp = ""
    for nt in seq:
        if (nt == "A"): comp = "T" + comp
        elif (nt == "T"): comp = "A" + comp
        elif (nt == "G"): comp = "C" + comp
        elif (nt == "C"): comp = "G" + comp
    return comp

# 正义链
def write_file_template(output, ChrID, filename, mode):
    with open(filename, mode) as res:
        temp = []
        for i in range(len(output.items())):
            temp.append(ChrID + "\t" + str(list(output.keys())[i]) + "\t" + str(list(output.values())[i][1]) +  "\t"  + str(list(output.values())[i][0]) +  "\t"  + "+" + "\n")

        res.writelines(temp)

# 反义链
def write_file_complete(output, ChrID, filename, mode):
    with open(filename, mode) as res:
        temp = []
        for i in range(len(output.items())):
            temp.append(ChrID + "\t" + str(list(output.keys())[i]) + "\t" + str(list(output.values())[i][1]) +  "\t"  + str(list(output.values())[i][0]) +  "\t"  + "-" + "\n")

        res.writelines(temp)

# 主函数
def all_chromosome(inputFile, inputType, N_times, cut_off):
    records = read_file(inputFile, inputType)
    for i in range(len(records)):
        origin_seq = records[i].seq.upper()
        bootResult = bootstrap(origin_seq, N_times)
        output = score_filter(bootResult, cut_off)
        ChrID =  records[i].id
        write_file_template(output, ChrID, "TATA_result_p_0.01.gff", "a+")
        
        origin_seq_comp = reverse_comp(records[i].seq.upper())
        bootResult_comp = bootstrap(origin_seq_comp, N_times)
        output_comp = score_filter(bootResult_comp, cut_off)
        write_file_complete(output_comp, ChrID, "TATA_result_p_0.01.gff", "a+")

all_chromosome("Sc_YJM993.fna", "fasta", 1000, 0.01)
