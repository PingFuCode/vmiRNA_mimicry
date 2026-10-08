import csv
from collections import defaultdict
import os

def process_mirnas(mirna_list, source_name):
    seed_dict = defaultdict(list)
    invalid_count = 0
    
    for mirna in mirna_list:
        if len(mirna) < 7:
            invalid_count += 1
            continue
        seed = mirna[1:7].upper()  # 提取第2-7nt并转为大写
        seed_dict[seed].append(mirna)
    
    # 统计种子区域使用情况
    count_single = 0
    count_multi = 0
    
    for seed, mirnas in seed_dict.items():
        count = len(mirnas)
        if count == 1:
            count_single += 1
        else:
            count_multi += 1
    
    print(f"/n{source_name} 结果:")
    print(f"无效序列数量 (长度<7): {invalid_count}")
    print(f"唯一使用的种子区域数量: {count_single}")
    print(f"共享使用的种子区域数量: {count_multi}")
    print(f"总有效种子区域: {len(seed_dict)}")
    
    return seed_dict

# 处理所有病毒miRNA (FASTA文件)
all_mirnas = []
fasta_path = r"D:/0work/0wholetransctiptome/2sRNAminic/1identify/data/vmiRNA.fa"

if not os.path.exists(fasta_path):
    print(f"错误: FASTA文件不存在 {fasta_path}")
else:
    with open(fasta_path, 'r') as f:
        sequence = ""
        for line in f:
            if line.startswith('>'):
                if sequence:
                    all_mirnas.append(sequence)
                    sequence = ""
            else:
                sequence += line.strip()
        if sequence:
            all_mirnas.append(sequence)

    all_seed_dict = process_mirnas(all_mirnas, "所有病毒miRNA")

# 处理模拟人病毒miRNA (CSV文件)
human_mimic_mirnas = []
csv_path = r"D:/0work/0wholetransctiptome/2sRNAminic/3mRNATarget/R/High_LowMimic.csv"

if not os.path.exists(csv_path):
    print(f"错误: CSV文件不存在 {csv_path}")
else:
    with open(csv_path, 'r', newline='') as f:
        reader = csv.DictReader(f)
        if 'Virus_fa' not in reader.fieldnames:
            print("错误: CSV文件中缺少 'Virus_fa' 列")
        else:
            seen = set()
            for row in reader:
                if row['Mimic'].strip() == 'LowMimic':
                    mirna = row['Virus_fa'].strip()
                    if mirna and mirna not in seen:
                        human_mimic_mirnas.append(mirna)
                        seen.add(mirna)
    
    human_seed_dict = process_mirnas(human_mimic_mirnas, "模拟人病毒miRNA")

# 可选：保存详细结果到文件
def save_seed_details(seed_dict, filename):
    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['种子序列', 'miRNA列表', 'miRNA数量'])
        for seed, mirnas in seed_dict.items():
            writer.writerow([seed, '; '.join(mirnas), len(mirnas)])

# 保存所有病毒miRNA的种子区域详情
save_seed_details(all_seed_dict, 'D:/0work/0wholetransctiptome/2sRNAminic/1seedregionInmiRNA/all_virus_seed_details.csv')
# 保存模拟人病毒miRNA的种子区域详情
save_seed_details(human_seed_dict, 'D:/0work/0wholetransctiptome/2sRNAminic/1seedregionInmiRNA/Low_details.csv')

print("/n详细结果已保存到: all_virus_seed_details.csv 和 human_mimic_seed_details.csv")