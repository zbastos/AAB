import allel
import matplotlib.pyplot as plt
import numpy as np
import plotly.plotly as py

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


if __name__=="__main__":

	ex1()

	#ex2()

