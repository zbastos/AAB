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
    outrpos = 0
    while (pos != (len(pattern))):
        string = ""
        primelems = list(k[0] for k in list(self.nodes[node][1].keys()))
        print("11")
        print(pos)
        print(len(pattern))
        elems = list(k for k in list(self.nodes[node][1].keys()))

        if pattern[pos] in primelems:
            print("22")

            indice = primelems.index(pattern[pos])
            string+=pattern[pos]
            outrpos = 0

            if len(elems[indice])==1:
                print("TTT")

                node = self.nodes[node][1][pattern[pos]]
                pos += 1
                

            else:
                while(elems[indice][outrpos]==pattern[pos] and pos != len(elems[indice])-1 and pos != len(pattern)-1):

                    print("93")
                    print(elems[indice][outrpos])
                    print(pattern[pos])
                    #print(len(elems[indice])-1)
                    pos += 1
                    outrpos += 1
                    string+=elems[indice][pos]
                    print(string)
                    
                if(elems[indice][outrpos]!=pattern[pos]):
                    print("aqui")
                    print(elems)
                    print(elems[indice][outrpos])
                    print(pattern[pos])
                    return None

                elif(pos==len(pattern)-1):
                    print("44")
                    print(elems[indice][outrpos])
                    print(pattern[pos])
                    pos += 1
                    outrpos += 1
                    print(pos)
                    node = self.nodes[node][1][elems[indice]]
                    print(node)
                    
                else:
                    print("55")
                    newpos = 0
                    if string in primelems:
                        ind = primelems.index(string)
                        string = elems[ind]
                    while(elems[indice][newpos]==pattern[pos] and newpos != len(elems[indice])-1 and pos != len(pattern)-1):
                        print("60")
                        pos += 1
                        newpos += 1
                    node = self.nodes[node][1][string]
                    print(node)
                    pos += 1
                    print(elems[indice][newpos])                  
        else:
            return None
    return getLeafesBelow(self,node)

'''
def findPattern(self, pattern):
    pos = 0
    node = 0
    outrpos = 0
    while (pos != (len(pattern))):
        string = ""
        primelems = list(k[0] for k in list(self.nodes[node][1].keys()))
        elems = list(k for k in list(self.nodes[node][1].keys()))

        if pattern[pos] in primelems:
            indice = primelems.index(pattern[pos])
            string+=pattern[pos]
            outrpos = 0

            if len(elems[indice])==1:

                node = self.nodes[node][1][pattern[pos]]
                pos += 1
                

            else:
                while(elems[indice][outrpos]==pattern[pos] and pos != len(elems[indice])-1 and pos != len(pattern)-1):

                    pos += 1
                    outrpos += 1
                    string+=elems[indice][pos]
                    
                if(elems[indice][outrpos]!=pattern[pos]):
                    return None

                elif(pos==len(pattern)-1):
                    pos += 1
                    outrpos += 1
                    node = self.nodes[node][1][elems[indice]]
                    
                else:
                    newpos = 0
                    if string in primelems:
                        ind = primelems.index(string)
                        string = elems[ind]
                    while(elems[indice][newpos]==pattern[pos] and newpos != len(elems[indice])-1 and pos != len(pattern)-1):
                        pos += 1
                        newpos += 1
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

def repeats(st,k,ocs):
	pos = 0
	res = []
	seq = getSeq(st,0,"")
	for pos in range(len(seq)-k+1):
		pat = seq[pos:pos+k]
		find = findPattern(st,pat)
		if find!=None and len(find) >= ocs:
			res.append(pat)
	if res:
		return set(res)
	else:
		return None


def getSeq(self,node,string):
    if self.nodes[node][0] < 0:
        valores = list(self.nodes[node][1].values())
        v = min(valores)
        for ke,va in list(self.nodes[node][1].items()):
        	if va==v:
        		k=ke
        string += k
        newnode = self.nodes[node][1][k]
        return getSeq(self,newnode,string)
    else:
        return string[:-1]


def ex2():
	
	seq = "TACTA"
	st = SuffixTrie()
	st.suffixTrieFromSeq(seq)
	compact(st)
	padrao = "TA"
	print("----------------------------------------------------------------------")
	print("		PROCURA DO PADRÃO *"+padrao+"* NA ÁRVORE")
	print("----------------------------------------------------------------------")
	posicoes = findPattern(st,padrao)
	print("POSIÇÕES: ")
	if posicoes!=None:
		posicoes.sort()
	print(posicoes)
	k=2
	ocs=2
	print("----------------------------------------------------------------------")
	print("		PADRÕES DE TAMANHO "+str(k)+" ENCONTRADOS "+str(ocs)+" OU MAIS VEZES")
	print("----------------------------------------------------------------------")
	res = repeats(st,k,ocs)
	print(res)

if __name__=="__main__":

	#ex1()

	ex2()

