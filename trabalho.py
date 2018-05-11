import allel
import matplotlib.pyplot as plt
import numpy as np
import plotly.plotly as py
from pattern_search import Trie, SuffixTrie

allmutations = {} #dicionário global que contém todas as mutações

def getVariants(fileName):
	data = allel.read_vcf(fileName)
	refs = data['variants/REF'] #obter lista das sequências de referência
	alts = data['variants/ALT'] #obter lista das sequências alternativas

	return refs,alts

def perg_a(fileName):

	refs,alts = getVariants(fileName)

	for r,a in zip(refs,alts): #para cada elemento em cada uma das listas
		mutation = r+'>'+a[0]
		if mutation not in allmutations.keys(): #se mutação não existe no dicionário, mete o valor a 1
			allmutations[mutation] = 1 
		else: 									#se já existe incrementa o valor
			allmutations[mutation]+=1
		if a[1]!="": 							#se existir outro sequência alternativa faz o mesmo
			mutation = r+'>'+a[1]
			if mutation not in allmutations.keys():
				allmutations[mutation] = 1
			else:
				allmutations[mutation]+=1

	return(str(len(allmutations))) #retorna o tamanho do dicionário

def perg_b(fileName):

	refs,alts = getVariants(fileName)

	snp = {} #dicionário que contem as SNPs e a quantas vezes aparecem

	for mut in allmutations.keys(): #para cada chave do dicionário
		if len(mut)==3:             #se o seu tamanho é 3 significa que é SNP
			snp[mut] = allmutations[mut] #logo, adiciona no dicionário de snps, o número de mutações dessa SNP

	labels = snp.keys() #eixo do x representado pelas chaves do dicionário
	exp = snp.values() #eixo do y representado pelos valores do dicionário
	numb = list(snp) 
	f = plt.figure()
	plt.bar(numb,exp, align = 'center') #criação do plot
	plt.title("Mutações SNP")
	plt.xticks(numb, labels);
	plt.xlabel("Mutação")
	plt.ylabel("Ocorrências")
	f.savefig("snp.pdf", bbox_inches='tight') #guarda o gráfico num ficheiro pdf na diretoria atual

def perg_c(fileName):

	refs,alts = getVariants(fileName)

	snp=de=ad=0
	for r,a in zip(refs,alts):
		tamref = len(r) #tamanho da seq de referência
		tamalt0 = len(a[0]) #tamanho da primeira seq alternativa
		if tamref==tamalt0:  #se for igual é uma snp
			snp+=1
		elif tamref>tamalt0: #se o tamanho da seq de referência, foi uma deleção
			de+=1
		else: ad+=1 #se não foi uma inserção
		if a[1]!="": #fazer o mesmo para a segunda sequência se esta não for vazia
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
	for string,pos in self.nodes[x][1].items(): #para cada key,value do segundo elemento do tuplo do nodo x
		seg = (self.nodes[pos][1])
		if len(seg)==1: #se o nodo contém apenas um elemento
			apagar = []
			while(self.nodes[pos][0]<0 and len(self.nodes[pos][1].values())==1): #enquanto nodo não for terminal e contiver apenas um elemento
				if(self.nodes[pos][0]<0 and pos!=x):
					apagar.append(pos)  #adicionar nodo à lista para apagar
				string+=(list(self.nodes[pos][1].keys())[0]) #retirar a key a adicionar à string
				pos = list(self.nodes[pos][1].values())[0] #retirar a posição
			self.nodes[x][1][string] = pos #atualizar o nodo inicial com a string e posição retirada
			self.nodes[x][1].pop(string[0],None) #apagar o elemento a mais
			for i in apagar: #apagar os nodos a mais da árvore
				self.nodes.pop(i,None)
		else:
			compactaux(self,pos) #se não tiver só um elemento fazer recursividade na árvore
		
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

def findPattern(self, pattern):
    pos = 0
    node = 0
    outrpos = 0
    while (pos != (len(pattern))):
        string = ""
        primelems = list(k[0] for k in list(self.nodes[node][1].keys())) #lista com os primeiros carateres dos elementos do nodo
        elems = list(k for k in list(self.nodes[node][1].keys())) #lista com os elementos do nodo

        if pattern[pos] in primelems: #se o carater do padrão está na lista dos primeiros elementos
            indice = primelems.index(pattern[pos]) #retirar o índice em que está contido na lista
            string+=pattern[pos] #adicionar o carater à string
            outrpos = 0

            if len(elems[indice])==1: #se o tamanho for apenas de 1

                node = self.nodes[node][1][pattern[pos]]
                pos += 1
                

            else:
                while(elems[indice][outrpos]==pattern[pos] and pos != len(elems[indice])-1 and pos != len(pattern)-1): 
                	#enquanto o padrão estiver igual à string do índice pretendido do nodo e enquanto houver padrão/string
                    pos += 1
                    outrpos += 1
                    string+=elems[indice][pos]
                    
                if(elems[indice][outrpos]!=pattern[pos]): #se for diferente retornar None
                    return None

                elif(pos==len(pattern)-1): #se acabar aí retirar o próximo nodo
                    pos += 1
                    outrpos += 1
                    node = self.nodes[node][1][elems[indice]]
                    
                else: #fazer o mesmo para o próximo nodo
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

def getLeafesBelow(self, node): #mesma função do pattern_search.py
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
	seq = getSeq(st,0,"") #função para obter a sequência original
	for pos in range(len(seq)-k+1): #técnica de sliding-window
		pat = seq[pos:pos+k]
		find = findPattern(st,pat)
		if find!=None and len(find) >= ocs: #Se no uso da função findPattern o output não for vazio e for maior que ocs, adicionar a lista de return
			res.append(pat)
	if res:
		return set(res)
	else:
		return None


def getSeq(self,node,string):
    if self.nodes[node][0] < 0: #verificar se não é nodo terminal
        valores = list(self.nodes[node][1].values()) #lista dos valores do segundo elemento do tuplo
        v = min(valores) #obter o menor valor da lista
        for ke,va in list(self.nodes[node][1].items()):
        	if va==v:
        		k=ke
        string += k #fazer append da chave à string
        newnode = self.nodes[node][1][k] #declaração doo próximo nodo
        return getSeq(self,newnode,string) #recursividade para executar a função para o próximo nodo, com a string calculada até ao momento
    else:
        return string[:-1] #eliminar o '$'


def ex2():
	
	seq = input("Sequência da árvore: ")
	st = SuffixTrie()
	st.suffixTrieFromSeq(seq)
	compact(st)


	padrao = input("Padrão: ")
	print("----------------------------------------------------------------------")
	print("		PROCURA DO PADRÃO *"+padrao+"* NA ÁRVORE")
	print("----------------------------------------------------------------------")
	posicoes = findPattern(st,padrao)
	print("POSIÇÕES: ")
	if posicoes!=None:
		posicoes.sort()
	print(posicoes)


	k=int(input("k: "))
	ocs=int(input("ocs: "))
	print("----------------------------------------------------------------------")
	print("		PADRÕES DE TAMANHO "+str(k)+" ENCONTRADOS "+str(ocs)+" OU MAIS VEZES")
	print("----------------------------------------------------------------------")
	res = repeats(st,k,ocs)
	print(res)

if __name__=="__main__":

	ex1()

	ex2()

