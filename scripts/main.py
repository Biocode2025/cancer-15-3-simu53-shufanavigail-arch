import random

def DNA_RNA_Cod(seq): #פונקציה שמקבלת רצף מקודד והופכת אותו לרנא
  RNA_seq = ''
  seq = seq.upper()
  for nuc in seq: # לולאה שמחליפה את נוקלאוטיד T לנוקלאוטיד U
    if nuc == "T":
      RNA_seq += 'U'
    else:
      RNA_seq += nuc  
  return RNA_seq

def Read_dict(file): #פונקציה שקוראת את הקובץ של הקודונים וחומצות האמינו והופפכת אותו למילון
  global RNA_codon_table
  for line in file:
    line = line.rstrip('\n')
    (Codon, sep, Amino_Acid) = line.partition('\t') #הרדה של הקודון מהחומצת אמינו בקובץ
    #print(Amino_Acid)
    RNA_codon_table[Codon] = Amino_Acid

def RNA_prot(seq): #פוקציה שמקבלת רצף רנא ומתרגמת אותו לרצף חלבונים
  prot_seq = ''
  for i in range(0, len(seq), 3): #לולאה שרצה על הרצף בקפיצות של 3
    codon = seq[i : i+3] #כל קודון שווה לשלישייה מהרצף
    Amino = RNA_codon_table[codon] #מציאת החומצת אמינו באמצעות המילון
    prot_seq += Amino
  return prot_seq

def Mutate_DNA(seq): #פונקציה שמקבלת רצף דנא ומחליפה נוקלאוטיד לנוקלאוטיד אחר במקום רנדומלי
  new_seq = ''
  base_list = ['A', 'C', 'G', 'T']
  seq = seq.upper()
  ran_place = random.randrange(0, len(seq)) #הגרלת מקום רנדומלי
  ran_nuc = seq[ran_place] #הנוקלאוטיד במקום הרנדומלי שנבחר
  base_list.remove(ran_nuc) #מחיקת הנוקלאוטיד שנבחר מרשימת הנוקלאוטיד
  ran_base = random.randrange(0,len(base_list)) #הגרלת בסיס רנדומלי
  new_nuc = base_list[ran_base] 

  seq_1 = seq[0 : ran_place]
  seq_2 = seq[ran_place + 1 :] 
  new_seq = seq_1 + new_nuc + seq_2 #יצירת הרצף עם הנוקלאוטיד החדש
  return new_seq

def Comp_seq(old,new): #פונקציה שמקבלת שני רצפים ומחזירה את מספר ההבדלים ביניהם
  diff_count = 0
  for i in range(0, len(old)): 
    if old[i] != new[i]: #אם הרצף במקום הi שונה מהרצף השני באותו מקום מוסיף לספירת ההבדלים 1
      diff_count += 1
  return diff_count

def Insert_DNA(seq):
    base_list = ['A', 'C', 'G', 'T']
    seq = seq.upper()
    nuc_num = random.randrange(1, 4)  # הגרלת מספר הנוקלאוטידים להוספה
    ran_place = random.randrange(0, len(seq))  # הגרלת מקום רנדומלי
    if nuc_num == 1:
        ran_base = random.randrange(0, len(base_list))
        new_nuc = base_list[ran_base]
    elif nuc_num == 2:
        new_nuc = ''
        for i in range(2):
            ran_base = random.randrange(0, len(base_list))
            new_nuc += base_list[ran_base]
    else:
        new_nuc = ''
        for i in range(3):
            ran_base = random.randrange(0, len(base_list))
            new_nuc += base_list[ran_base]

    seq_1 = seq[0: ran_place]
    seq_2 = seq[ran_place :]
    new_seq = seq_1 + new_nuc + seq_2
    return new_seq

def Delete_DNA(seq):
    seq = seq.upper()
    nuc_num = random.randrange(1, 4)  # הגרלת מספר הנוקלאוטידים למחיקה
    ran_place = random.randrange(0, len(seq) - nuc_num)  # הגרלת מקום רנדומלי
    seq_1 = seq[0: ran_place]
    seq_2 = seq[ran_place + nuc_num:] 
    new_seq = seq_1 + seq_2
    return new_seq

gen = 1000
RNA_codon_table = {}
gen_count = 1
diff_count = 0
overall_gen = 0
dna_file = open('data/human_p53_coding.txt', 'r')
codon_file = open('data/codon_AA.txt', 'r')
dna_seq = ''
#main#
for line in dna_file:
  if line[0] != ">":
    line = line.rstrip('\n')
    dna_seq += line
Read_dict(codon_file)

rna_seq = DNA_RNA_Cod(dna_seq)
prot_seq = RNA_prot(rna_seq)


for i in range(gen):
    diff_count = 0
    gen_count = 0
    mutated_dna = dna_seq  
    while diff_count == 0:
         gen_count += 1
         mutation = random.randrange(1, 10001) #הגרלת מספר בין 1 ל10000 כדי לקבוע אם תתבצע מוטציה
         if mutation == 1:       
             chance = random.randrange(1, 101) #הגרלת מספר בין 1 ל100 כדי לקבוע איזה סוג מוטציה תתבצע
             if chance <= 98:
                 mutated_dna = Mutate_DNA(mutated_dna)
             elif chance == 99:
                 mutated_dna = Insert_DNA(mutated_dna)
             else:
                 mutated_dna = Delete_DNA(mutated_dna)
                
             remain = len(mutated_dna) % 3 # מחזיר את השארית אם הרצף לא מתחלק בשלוש כדי להוות את הנוקלאוטידים שנשארים 
             mutated_dna = mutated_dna[0:len(mutated_dna) - remain]
             mutated_rna = DNA_RNA_Cod(mutated_dna)
             mutated_prot = RNA_prot(mutated_rna)
             if len(mutated_prot) != len(prot_seq):
                 diff_count += 1
             else:
                 diff_count = Comp_seq(prot_seq, mutated_prot)
    overall_gen += gen_count
        

print(overall_gen)
print("Average number of dna replications until mutation that resulted in protein change: %d" %(round(overall_gen / gen)))

dna_file.close()
codon_file.close()
