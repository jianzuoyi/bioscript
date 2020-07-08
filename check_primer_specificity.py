#!/usr/bin/env python
import sys
# 程序功能：检查扩增子引物特异性
# 输入：引物blast到ncbi代表性基因组上的结果，列表格式（-outfmt 6）
lines = []
with open('primer-58-to-rep-all.m6', 'r') as fh:
    lines = [line.rstrip().split() for line in fh ]

amplicons = {}    # 用于保存Forward引物的配对引物，键=forword primer，值=[amplicon]
for F in lines:
    if F[0][-2] != 'F':
        continue
    # 第一个引物所在的链
    primer1_strand = "+" if int(F[8]) < int(F[9]) else "-"

    # 第一个引物 基因组_开始_结束_引物名
    primer1_pos = "_".join([F[1], F[0], F[8], F[9]])

    # 循环列表，找出与前一个引物配对的引物
    for R in lines:
        if any([R[0]==F[0], R[1]!=F[1]]):   # 引物名是自己，或者所在的基因组序列名不一样，都忽略
            continue

        primer2_strand = "+" if int(R[8]) < int(R[9]) else "-"
        # 条件1： 第1个引物在正链， 第2个引物在负链，第1个引物的起始位置小于第2个引物的起始位置（确保能扩出东西），扩增子长度 < 1000bp
        # 条件2： 第1个引物在负链， 第2个引物在正链，第1个引物的起始位置大于第2个引物的起始位置（确保能扩出东西），扩增子长度 < 1000bp
        cond1 = [primer1_strand=="+", primer2_strand=="-", int(R[8])>int(F[8]), int(R[8])-int(F[8])<1000]
        cond2 = [primer1_strand=="-", primer2_strand=="+", int(R[8])<int(F[8]), int(F[8])-int(R[8])<1000]
        if all(cond1) or all(cond2):
            amplicon_len = abs(int(R[8])-int(F[8]))
            primer2_pos = "_".join([R[0], R[8], R[9], primer1_strand+primer2_strand, str(amplicon_len)])
            print("\t".join([primer1_pos, primer2_pos]))
            print("\t".join(F))
            print("\t".join(R))

            if primer1_pos not in amplicons:
                amplicons[primer1_pos] = [primer2_pos]
            else:
                amplicons[primer1_pos].append(primer2_pos)

    #break

# 字典：用来记录 Forward 引物所在的基因组
genomes = {}
for pos1 in amplicons:
    primer_fp = pos1.split('_')[-3]
    if primer_fp not in genomes:
        genomes[primer_fp] = []
    for p in amplicons[pos1]:
        genome = '_'.join(pos1.split('_')[:-3])
        genomes[primer_fp].append(genome)

# 接上面，Forward 引物为键，值存储的是基因组列表，通过直接计算list的长度，以及经过set()过的list长度，就能知道list中是否有重复基因组        
amplicon_in_one_genome = {}
for p in genomes:
    flag = 'uniq_amplicon_in_one_genome' if len(genomes[p]) == len(set(genomes[p])) else 'multi_amplicon_in_one_genome'
    amplicon_in_one_genome[p] = flag

# 循环每个引物1
for pos1 in amplicons:
    primer_fp = pos1.split('_')[-3]
    # 循环与每个引物1配对的引物2
    for p in amplicons[pos1]:
        primer_tp = p.split('_')[0]
        primer_tp = primer_tp[:-2]+"F"+primer_tp[-1]
        flag = 'primer_to_self-primer' if primer_tp == primer_fp else 'primer_to_other_primer'
        # 输出结果，每一行为一个扩增子
        print('\t'.join([primer_fp, pos1, p, flag, amplicon_in_one_genome[primer_fp]]))
