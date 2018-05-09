import allel
import matplotlib.pyplot as plt
import numpy as np
import plotly.plotly as py
from pattern_search import Trie, SuffixTrie

allmutations = {}

def getVariants(fileName):
	data = allel.read_vcf(fileName)
	refs = data['variants/REF']
	alts = data['variants/ALT']

	return refs,alts

def perg_a(fileName):

	refs,alts = getVariants(fileName)

	for r,a in zip(refs,alts):
		mutation = r+'>'+a[0]
		if mutation not in allmutations.keys():
			allmutations[mutation] = 1
		else:
			allmutations[mutation]+=1
		if a[1]!="": 
			mutation = r+'>'+a[1]
			if mutation not in allmutations.keys():
				allmutations[mutation] = 1
			else:
				allmutations[mutation]+=1

	return(str(len(allmutations)))
	#for m in (sorted(allmutations, key=allmutations.get,reverse = True)):
	#	print (m+" "+str(allmutations[m]))

def perg_b(fileName):

	refs,alts = getVariants(fileName)

	snp = {}

	for mut in allmutations.keys():
		if len(mut)==3:
			snp[mut] = allmutations[mut]

	labels = snp.keys()
	exp = snp.values()
	numb = list(snp)
	f = plt.figure()
	plt.bar(numb,exp, align = 'center')
	plt.title("Mutações SNP")
	plt.xticks(numb, labels);
	plt.xlabel("Mutação")
	plt.ylabel("Ocorrências")
	f.savefig("snp.pdf", bbox_inches='tight')

def perg_c(fileName):

	refs,alts = getVariants(fileName)

	snp=de=ad=0
	for r,a in zip(refs,alts):
		tamref = len(r)
		tamalt0 = len(a[0])
		if tamref==tamalt0:
			snp+=1
		elif tamref>tamalt0:
			de+=1
		else: ad+=1
		if a[1]!="":
			tamalt1 = len(a[1])
			if tamref==tamalt1:
				snp+=1
			elif tamref>tamalt1:
				de+=1
			else: ad+=1

	return(snp,de,ad)


def ex1():
	
	fileName = 'HG00154.chr21.raw.vcf'

	allmutations = {}

	print("--------------------------")
	mutacoes = perg_a(fileName)
	print("PERGUNTA A:\nMutações únicas - "+mutacoes)
	print("--------------------------")

	perg_b(fileName)
	print("PERGUNTA B:\nGráfico concluído e guardado na diretoria num ficheiro snp.pdf")
	print("--------------------------")


	(snp,de,ad) = perg_c(fileName)
	print("PERGUNTA C:")
	print("SNPs - "+str(snp))
	print("Deleções - "+str(de))
	print("Inserções - "+str(ad))
	print("--------------------------")


def compact(self):
	tam = len(self.nodes.keys())-1
	pos = 2
	string = ""
	while(pos!=tam):
		if(len(self.nodes[pos][1])==1):
			ini = list(self.nodes.keys())[pos]
			#print("ini "+str(ini))
			apagar = []
			while(self.nodes[pos][0]<0):
				if(self.nodes[pos][0]<0 and pos!=ini):
					apagar.append(pos)
					#print(apagar)
				string+=(list(self.nodes[pos][1].keys())[0])
				pos = list(self.nodes[pos][1].values())[0]
				#print("pos"+str(pos))
				#print(string)
				#print(self.nodes[pos][0])
			self.nodes[ini][1][string] = pos
			self.nodes[ini][1].pop(string[0],None)
			for i in apagar:
				#print(i)
				self.nodes.pop(i,None)
			#print(self.nodes[ini][1])
			print(self.nodes)
		else:
			#for k in self.nodes[pos][1].values():
			pos = 13
			#pos = 9

	

def ex2():
	
	seq = "ATAA"
	st = SuffixTrie()
	st.suffixTrieFromSeq(seq)
	st.print_tree()
	ct = compact(st)

if __name__=="__main__":

	#ex1()

	ex2()

