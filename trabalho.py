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

'''
def compact(self):
	tam = len(self.nodes.keys())-1 #13
	pos = 2
	apagados = 0
	string = ""
	while(pos!=tam):
		if(len(self.nodes[pos][1])==1):
			ini = list(self.nodes.keys())[pos] - apagados
			print("ini "+str(ini))
			apagar = []
			while(self.nodes[pos][0]<0):
				if(self.nodes[pos][0]<0 and pos!=ini):
					apagar.append(pos)
					print(apagar)
				string+=(list(self.nodes[pos][1].keys())[0])
				pos = list(self.nodes[pos][1].values())[0]
				print("pos"+str(pos))
				print(string)
				print(self.nodes[pos][0])
			self.nodes[ini][1][string] = pos
			self.nodes[ini][1].pop(string[0],None)
			for i in apagar:
				self.nodes.pop(i,None)
				apagados += 1
			print(self.nodes[ini][1])
			print(self.nodes)
			string = []
			pos+=1
		else:
			#for k in self.nodes[pos][1].values():
			pos = 13
			#pos = 9


def test(self):
	apagados = 0
	string = ""

	for i in self.nodes.keys():
		for k,v in self.nodes[i][1].items():
			symbol = k
			pos = v
			seg = (self.nodes[pos][1])
			print(pos)
			if len(seg)==1:
				ini = list(self.nodes.keys())[pos] - apagados
				print("ini "+str(ini))
				apagar = []
				while(self.nodes[pos][0]<0):
					if(self.nodes[pos][0]<0 and pos!=ini):
						apagar.append(pos)
						print(apagar)
					string+=(list(self.nodes[pos][1].keys())[0])
					pos = list(self.nodes[pos][1].values())[0]
					print("pos"+str(pos))
					print(string)
					print(self.nodes[pos][0])
				self.nodes[ini][1][string] = pos
				self.nodes[ini][1].pop(string[0],None)
				for i in apagar:
					self.nodes.pop(i,None)
					apagados += 1
				print(self.nodes[ini][1])
				print(self.nodes)
				string = []
				pos+=1
			else:
				print("NNN")
'''
def compactaux(self,x):
	apagados = 0
	it = 0
	for k,v in self.nodes[x][1].items():
		it+=1
		string = k
		pos = v
		seg = (self.nodes[pos][1])
		if len(seg)==1:
			apagar = []
			while(self.nodes[pos][0]<0):
				if(self.nodes[pos][0]<0 and pos!=x):
					apagar.append(pos)
				string+=(list(self.nodes[pos][1].keys())[0])
				pos = list(self.nodes[pos][1].values())[0]
			self.nodes[x][1][string] = pos
			self.nodes[x][1].pop(string[0],None)
			for i in apagar:
				self.nodes.pop(i,None)
				apagados += 1
		else:
			compactaux(self,pos)
		
def compact(self):
	print("----------------------------------------------------------")
	print("		ÁRVORE ANTES DA COMPACTAÇÃO")
	print("----------------------------------------------------------")
	self.print_tree()
	compactaux(self,0)
	print("----------------------------------------------------------")
	print("		ÁRVORE DEPOIS DA COMPACTAÇÃO")
	print("----------------------------------------------------------")
	self.print_tree()

def ex2():
	
	seq = "ATGAAC"
	st = SuffixTrie()
	st.suffixTrieFromSeq(seq)
	sim = compact(st)

if __name__=="__main__":

	#ex1()

	ex2()

