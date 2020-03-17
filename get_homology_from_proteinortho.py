# -*- coding: utf-8 -*-
import os           #导入os模块
import numpy as np
import re
global now_dir
from Bio import SeqIO
now_dir = os.getcwd()           #得到当前目录  

def rename_proteinortho_output(proteinortho_output_name):  #该函数将proteinortho软件输出文件中的蛋白质文件名重新命名，命名结构遵循："species_name"++"protein_name" 

  with open(now_dir + "/renamed_" + proteinortho_output_name,"a") as write_file:
    write_file.write(",\t")
    with open(now_dir + "/" + proteinortho_output_name,"r") as read_file:
      first_line = read_file.readline()
      species_name_list = first_line.split("\t")[3:]
      species_name_list[-1] = species_name_list[-1][:-1]
      '''
      以上代码得到了一个包含所有物种名字的列表
      '''
      for each_line in read_file:
        for each_gene_clustr in each_line.split("\t")[:3]:
          write_file.write(each_gene_clustr + "\t")
        n = 0
        for each_gene_clustr in each_line.split("\t")[3:]:
          new_clustr_name = []
          for each_gene in each_gene_clustr.split(","):
            new_clustr_name.append(species_name_list[n] + "++" + each_gene + ",")
          write_file.write("".join(new_clustr_name))
          n = n + 1
          write_file.write("\t")
      '''
      以上代码将所有物种名字添加到了文件里面的蛋白质名字之前
      '''
  with open(now_dir + "/renamed_" + proteinortho_output_name + "_temp","a") as temp_write_file:
    with open(now_dir + "/renamed_" + proteinortho_output_name,"r") as temp_read_file:
      for each_line in temp_read_file:
        try:
          if int(each_line.split("\t")[1]) > 29:
            temp_write_file.write(each_line[2:-1] + "," + "\n")
        except:
          pass
  '''
  以上代码把刚刚生成的文件做一下整理，方便之后使用
  '''
  return(species_name_list)
      
      

    
file_dict = {}  
for each in os.listdir(os.getcwd()):        
    if ".pep" in each:
        file_dict[each] = SeqIO.to_dict(SeqIO.parse(each, "fasta"))
    elif ".cds" in each:
        file_dict[each] = SeqIO.to_dict(SeqIO.parse(each, "fasta"))
    else:
        pass

def get_protein(IDs): #该函数根据一个proteinortho输入文件中的基因ID号从基因库中提取对应的蛋白质序列
    return str(file_dict[IDs.split("++")[0][:-3] + "pep"][IDs.split("++")[1]].seq)      
        
       
def get_nucleotide(IDs): #该函数根据一个proteinortho输入文件中的基因ID号从基因库中提取对应的蛋白质序列
    return str(file_dict[IDs.split("++")[0][:-3] + "cds"][IDs.split("++")[1]].seq)               
       

def write_protein_result(protein_name):
    with open(now_dir + "/output/" + protein_name + "_pep.fas","a") as write_file:
        for each_species in species_list:
            if each_line.count(each_species) > 1:
                write_file.close()
                os.remove(now_dir + "/output/" + protein_name + "_pep.fas")
                break
            elif each_line.count(each_species) > 1:
                pattern_protein = re.compile(each_species + "\+\+.*?\,")
                pattern_protein_select1 = re.findall(pattern_protein,each_line)
                #print(pattern_protein_select1)
                a = "a"
                for each in pattern_protein_select1:
                    protein = get_protein(each[:-1])
                    if len(protein) > len(a):
                        a = protein
                        b = each
                write_file.write(">" + protein_name + "|" + b + "\n" + a + "\n")
            else:
                pattern_protein = re.compile(each_species + "\+\+.*?\,")
                pattern_protein_select3 = re.findall(pattern_protein,each_line)[0]
                if pattern_protein_select3[-2] == "*":
                    pass
                else:
                    pattern = get_protein(pattern_protein_select3[:-1])
                    write_file.write(">" + protein_name + "|" + pattern_protein_select3 + "\n" + pattern + "\n")

def write_nucleotide_result(protein_name):
    with open(now_dir + "/output/" + protein_name + "_cds.fasta","a") as write_file:
        for each_species in species_list:
            if each_line.count(each_species) > 1:
                write_file.close()
                os.remove(now_dir + "/output/" + protein_name + "_cds.fasta")
                break
            elif each_line.count(each_species) > 1:
                pattern_protein = re.compile(each_species + "\+\+.*?\,")
                pattern_protein_select1 = re.findall(pattern_protein,each_line)
                #print(pattern_protein_select1)
                a = "a"
                for each in pattern_protein_select1:
                    protein = get_nucleotide(each[:-1])
                    if len(protein) > len(a):
                        a = protein
                        b = each
                write_file.write(">" + protein_name + "|" + b + "\n" + a + "\n")
            else:
                pattern_protein = re.compile(each_species + "\+\+.*?\,")
                pattern_protein_select3 = re.findall(pattern_protein,each_line)[0]
                if pattern_protein_select3[-2] == "*":
                    pass
                else:
                    pattern = get_nucleotide(pattern_protein_select3[:-1])
                    write_file.write(">" + protein_name + "|" + pattern_protein_select3 + "\n" + pattern + "\n")


proteinortho_output_name = "myproject.proteinortho.tsv"
species_list = rename_proteinortho_output(proteinortho_output_name)



with open(now_dir + "/renamed_" + proteinortho_output_name + "_temp", "r") as ortho_file:
  n = 1
  for each_line in ortho_file:
      if each_line:
          protein_name = "ortho" + str(n)
          write_protein_result(protein_name)
          write_nucleotide_result(protein_name)
          n = n + 1
          

                    
