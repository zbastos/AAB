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


def compactaux(self,x):
	for string,pos in self.nodes[x][1].items():
		seg = (self.nodes[pos][1])
		if len(seg)==1:
			apagar = []
			while(self.nodes[pos][0]<0 and len(self.nodes[pos][1].values())==1):
				if(self.nodes[pos][0]<0 and pos!=x):
					apagar.append(pos)
				string+=(list(self.nodes[pos][1].keys())[0])
				pos = list(self.nodes[pos][1].values())[0]
			self.nodes[x][1][string] = pos
			self.nodes[x][1].pop(string[0],None)
			for i in apagar:
				self.nodes.pop(i,None)
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
'''
def findPattern(self, pattern):
    pos = 0
    node = 0
    while (pos != (len(pattern))):
        string = ""
        primelems = list(k[0] for k in list(self.nodes[node][1].keys()))
        elems = list(k for k in list(self.nodes[node][1].keys()))
        print("coiso"+pattern[pos])
        if pattern[pos] in primelems:
            indice = primelems.index(pattern[pos])
            string+=pattern[pos]
            print("aqui")
            print(elems[indice])
            if len(elems[indice])==1:
                print("1")
                node = self.nodes[node][1][pattern[pos]]
                pos += 1
            else:
                print("2")
                #print(elems[indice][pos])
                #print(pattern[pos])
                while(elems[indice][pos]==pattern[pos] and pos != len(elems[indice])-1):
                    print("3")
                    pos += 1
                    string+=pattern[pos]
                    #print(elems[indice][pos])
                    #print(pattern[pos])
                    #print(pos)
                    #print(len(pattern))
                if(pos==len(pattern)-1):
                    print("welelelelel")
                    print(node)
                    print(primelems)
                    print(string)
                    pos += 1
                    node = self.nodes[node][1][elems[indice]]
                    
                else:
                    print("4")
                    print(string)
                    print(node)
                    node = self.nodes[node][1][string]
                    pos += 1
                    print(node)
                    print(pattern[pos])                    
        else:
            return None
    return getLeafesBelow(self,node)
'''

def findPattern(self, pattern):
    pos = 0
    node = 0
    while (pos != (len(pattern))):
        string = ""
        primelems = list(k[0] for k in list(self.nodes[node][1].keys()))
        elems = list(k for k in list(self.nodes[node][1].keys()))

        if pattern[pos] in primelems:

            indice = primelems.index(pattern[pos])
            string+=pattern[pos]
            
            if len(elems[indice])==1:

                node = self.nodes[node][1][pattern[pos]]
                pos += 1

            else:
                while(elems[indice][pos]==pattern[pos] and pos != len(pattern)-1):

                    pos += 1
                    string+=pattern[pos]
                    

                if(pos==len(pattern)-1):
                    pos += 1
                    node = self.nodes[node][1][elems[indice]]
                    
                else:
                    node = self.nodes[node][1][string]
                    pos += 1                  
        else:
            return None
    return getLeafesBelow(self,node)

def getLeafesBelow(self, node):
        res = []
        if self.nodes[node][0] >= 0:
        	res.append(self.nodes[node][0])
        else:
            for k in self.nodes[node][1].keys():
            	newnode = self.nodes[node][1][k]
            	leafes = getLeafesBelow(self,newnode)
            	res.extend(leafes)
        return res

def repeats(st,seq,k,ocs):
	pos = 0
	res = []
	print(len(seq)-k)
	for pos in range(len(seq)-k):
		pat = seq[pos:pos+k]
		print(pat)
		find = findPattern(st,pat)
		print(find)
		if len(find) >= ocs:
			res.append(pat)
	print(set(res))





def ex2():
	
	seq = "GTAACTAAGAA"
	st = SuffixTrie()
	st.suffixTrieFromSeq(seq)
	compact(st)
	padrao = "AA"
	print("----------------------------------------------------------")
	print("		PROCURA DO PADRÃO *"+padrao+"* NA ÁRVORE")
	print("----------------------------------------------------------")
	posicoes = findPattern(st,"AA")
	print("POSIÇÕES: ")
	posicoes.sort()
	print(posicoes)
	print("||||||||||||||||||||||||||||")
	#repeats(st,seq,2,3)

if __name__=="__main__":

	#ex1()

	ex2()

