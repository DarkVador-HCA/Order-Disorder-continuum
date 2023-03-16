import pickle as pk
import sys

""" Creation of Datasets S2 and S3 from each database sequences and from the hcatk output """

def calcule_cov(bornes,taille):
        cov = 0
        for q in bornes:
            cov += (q[1]-q[0]+1)
        if taille==0:
            return 0
        cov=(cov/taille)
        return cov

class Protein:
	def __init__(self):
		self.size=0 #protein size
		self.hca_all=0 #HCA score over all protein
		self.sizes_dom=[] #foldable segment size
		self.hcas_dom=[] #foldable segment HCA score
		self.sizes_hc=[] #hydrophobic cluster size (1 list by foldable segment)
		self.bornes_hc=[] #hydrophobic cluster bounds (1 list by foldable segment)
		self.hc=[] # hydrophobic cluster binary code (1 list by foldable segment)
		self.nb_hc=[] # number of hydrophobic cluster  by foldable segment
		self.bornes_foldable=[] #foldable segment bounds
		self.cov_foldable=0 #coverage of amino acids in a foldable segments (computed over the whole protein sequence)
		self.sf_nbH=[] #number of hydrophobic amino acids by foldable segment
		self.covH=0 #% of hydrophibic residues in the whole protein sequence
		self.seq=[] #amino acid sequence

	def get_covs(self):
		self.cov_foldable=calcule_cov(self.bornes_foldable,self.size)

def parse_fasta(fasta,Prots):
    prot=None
    with open(fasta,"r") as f:
        for l in f:
            if l[0]==">":
                prot=l[1:].strip().split()[0].split("|")[-1]

                Prots[prot]=Protein()
            else:
                Prots[prot].seq+=list(l.strip())
    return Prots

def parse_segHCA(hca_file,Prots):
	prot=""
	with open(hca_file,"r")as file:
		for l in file:
			if l[0]!='#':
				if l[0]=='>':
					if prot!="":
						Prots[prot].covH/=Prots[prot].size
					l=l[1:].strip().split()
					prot=l[0].split('|')[-1]
					Prots[prot].size=int(l[1])
					Prots[prot].hca_all=float(l[-1])
					foldable=[]
					want=False

				elif l[:6]=='domain':
					l=l.strip().split()
					b1,b2=int(l[1]),int(l[2])
					Prots[prot].bornes_foldable.append([b1,b2])
					Prots[prot].hcas_dom.append(l[-1])
					Prots[prot].sizes_dom.append(b2-b1+1)
					Prots[prot].hc.append([])
					Prots[prot].bornes_hc.append([])
					Prots[prot].sizes_hc.append([])
					Prots[prot].nb_hc.append(0)
					Prots[prot].sf_nbH.append(0)
					want=True

				elif l[:7]=='cluster':
					l=l.strip().split()
					Prots[prot].covH+=l[-1].count("1")
					if want:
						i=0
						b1=int(l[1])
						b2=int(l[2])
						done=False
						while i<len(Prots[prot].bornes_foldable) and not done:
							if b1>=Prots[prot].bornes_foldable[i][0] and b2<=Prots[prot].bornes_foldable[i][1]:
								Prots[prot].hc[i].append(l[-1])
								Prots[prot].bornes_hc[i].append((int(l[1]),int(l[2])))
								Prots[prot].sizes_hc[i].append(len(l[-1]))
								Prots[prot].sf_nbH[i]+=l[-1].count("1")
								if l[-1] not in ["1","11"]:
									Prots[prot].nb_hc[i]+=1
								if prot=="d1dw0a_":
									done=True
							i+=1


	return Prots

def save_scores(Prots):
    l1=[Prots[p].hca_all for p in Prots.keys()]
    with open("all_seq.pkl","wb") as f:
        pk.dump(l1,f)


def get_ids_rm(db):
    with open("../data/"+db+"/trainID.lst","r") as f:
        train=[l.strip() for l in f]
    return train

def get_id_cat(infile):
    id_cat={}
    cats=[]
    with open(infile,'r') as f:
        for l in f:
            id_cat[l.strip().split(",")[1]]=l.split(",")[0]
            if l.split(",")[0] not in cats:
                cats.append(l.split(",")[0])
    return id_cat,cats

def get_rmaa_OPM():
	with open("../data/OPM/data/cut_data/nb_rm_aa30.csv","r") as f:
		rm_aa={l.split(',')[0]:l.strip().split(",")[1] for l in f}
	return rm_aa


def write_tbls(Prots,db,cats=False,rm_aa=False):
	tbl_all=[["database","class","sequence id","length (aa)","HCA score","% strong hydrophobic","sequence coverage by foldable segment",'number of removed aa']]
	tbl_rf=[["database","class","sequence id","start","length (aa)","HCA score","% strong hydrophobic","number  of hydrophobic clusters","start of hydrophobic clusters",'hydrophobic clusters (binary)']]
	cat=db
	for p in Prots.keys():
		if cats:
			if db=="SCOPe":
				cat=cats[p]
			elif db=="DisProt":
				cat=cats[p]#.split('r')[0]]
			elif  db=="OPM":
				cat=cats[p.split('_')[0]]

		if rm_aa:
			lessaa=rm_aa[p]
			if db!="OPM":
				print("oups "+db)
		else:
			lessaa="0"

		Prots[p].get_covs()
		tbl_all.append([db,cat,p,str(Prots[p].size),"%0.2f"%Prots[p].hca_all,"%0.2f"%Prots[p].covH,"%0.2f"%Prots[p].cov_foldable,lessaa])

		for i in range(len(Prots[p].bornes_foldable)):
			b=Prots[p].bornes_foldable[i]
			sizesHC=[]
			for h in Prots[p].sizes_hc[i]:
				sizesHC.append(str(h))
			startsHC=[]
			for v in Prots[p].bornes_hc[i]:
				startsHC.append(v[0])
			tbl_rf.append([db,cat,p,str(b[0]), str(Prots[p].sizes_dom[i]),"%0.2f"%Prots[p].hcas_dom[i],"%0.2f"%(Prots[p].sf_nbH[i]/Prots[p].sizes_dom[i]),str(Prots[p].nb_hc[i]),",".join([str(v) for v in startsHC]),",".join(Prots[p].hc[i])])


	with open("../../tables/"+db+"/datasetS2.csv","w") as f:
		for l in tbl_all:
			f.write("%s\n"%";".join(l))

	with open("../../tables/"+db+"/datasetS3.csv","w") as f:
		for l in tbl_rf:
			f.write("%s\n"%";".join(l))

if __name__=="__main__":
	fa=sys.argv[1]
	hca=sys.argv[2]
	db=sys.argv[3]

	Prots=parse_fasta(fa,{})
	Prots=parse_segHCA(hca,Prots)
	if db in ["SCOPe","DisProt"]:
		id_cat,cats=get_id_cat("../data/"+db+"/id_cats.csv")
		if db in ["SCOPe","DisProt"]:
			id_train=get_ids_rm(db)
		write_tbls(Prots,db,cats=id_cat,train=id_train)

	elif db=="OPM":
		rm_aa=get_rmaa_OPM() #residues in soluble segments removed to keep only membrane segments
		id_cat,cats=get_id_cat("../data/"+db+"/id_cats.csv")
		write_tbls(Prots,db,cats=id_cat,rm_aa=rm_aa)

	else:
		write_tbls(Prots,db)
