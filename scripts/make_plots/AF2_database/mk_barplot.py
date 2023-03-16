import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

"""
création des barplots (score HCA x pLDDT) à partir de la table : table_AF2DB.csv (script get_table.py)
"""

tableAF2DB="../../../tables/datasetS4.csv"

cats=["vl","l","c","vh"]

labels={"vl":"very low (pLDDT≤50)","l":"low (70≥pLDDT>50)","c":"confident (90≥pLDDT>70)","vh":"very high (pLDDT>90)"}
colors={"vl":"darkorange","l":"gold","c":"lightskyblue","vh":"darkblue"}

def get_barplot_data(df):
    """
    formatage des données pour les barplots
    """
    cats_nbaa={c:[0 for _ in range(-6,10)] for c in cats}
    allaa=[0 for _ in range(-6,10)]
    nb_allSF=[0 for _ in range(-6,10)]
    i=0
    for h in range(-6,10):
        dh=df[(df["HCA score"]<=h) & (df["HCA score"]>h-1)]
        allaa[i]=dh["contiguous regions length (aa)"].sum()
        nb_allSF[i]=len(dh.drop_duplicates(subset=['sequence id','whole segment start']))
        for c in cats:
            if allaa[i]==0:
                cats_nbaa[c][i]=0
            else:
                cats_nbaa[c][i]=100*dh[dh["pLDDT category"]==c]["contiguous regions length (aa)"].sum()/allaa[i]
        i+=1
    return cats_nbaa,allaa, nb_allSF



def mk_hca_barplot(cats_nbaa,nb_allaa,nb_allSF,figname,cats,labels,fig_title="",save=False):
    """
    Création du bar plot:
    cats_nbaa,nb_allaa,nb_allSF --> objets générés par la fonction get_barplot_data()
    figname --> pour sauvegarder la figure (adresse+nom du fichier + extension)
    cats --> nom des catégories (? pourquoi un argument ?)
    labels --> description à afficher dans la légende pour chaque catégorie
    fig_title="" --> titre à la figure (ex : code du protéome représenté)
    """
    f, (a0, a1) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [1, 3]},figsize=[6,7])
    # 1 : nb residues et SF par catégories HCA
    plt.sca(a0)
    plt.title(fig_title)
    a0.bar(np.arange(-6.5,9.5,1),nb_allaa,1,color="white",edgecolor="black")
    ymin,ymax=plt.ylim()
    a0.set_yscale("log")
    a0.set_ylabel("Residues")
    a0.vlines([-4.7,-1,3.5,6.7,7.7],ymin,ymax,linestyle="--",lw=1,color="grey")
    plt.xlim(-7,9)
    plt.xticks(range(-7,10,1),range(-7,10,1))
    a02=a0.twinx()
    # make a plot with different y-axis using second axis object
    a02.plot(np.arange(-6.5,9.5,1),nb_allSF,color="darkgray",marker=".")
    a02.set_ylabel("Foldable segments",color="darkgray")
    a02.tick_params(axis="y",colors="darkgray")
    a02.set_yscale("log")
    
    plt.sca(a1)
    y_offset=np.zeros(len(nb_allaa))
    for c in cats:
        a1.bar(np.arange(-6.5,9.5,1),cats_nbaa[c],1,bottom=y_offset,color=colors[c],label=labels[c])
        y_offset=y_offset+cats_nbaa[c]
    # a1.legend()
    plt.ylim((0,100))
    ymin,ymax=plt.ylim()
    plt.vlines([-4.7,-1,3.5,6.7,7.7],ymin,ymax,linestyle="--",lw=1,color="grey")
    plt.xlabel("HCA score")
    plt.ylabel("% of residues")
    plt.xlim(-7,9)
    plt.xticks(range(-7,10,1),range(-7,10,1))
    
    f.tight_layout()
    if save:
        f.savefig(figname)
        print("saved")
    plt.show()

def from_table_to_barplots(table_path,toplot):
    df=pd.read_csv(table_path,sep=";") 

    cats_nbaa,allaa,nb_allSF=get_barplot_data(df)

    if toplot=="all proteomes":
        mk_hca_barplot(cats_nbaa,allaa,nb_allSF,"barplot_all.pdf",cats,labels,fig_title="all proteomes") # tous les protéomes combinés

    else:
        o=toplot
        cats_nbaa,allaa,nb_allSF=get_barplot_data(df[df["proteome"]==o]) # Proteins of organism o 
        mk_hca_barplot(cats_nbaa,allaa,nb_allSF,o+"_barplot.pdf",cats,labels,fig_title=o)
    plt.close('all')

if __name__=="__main__":
    df=pd.read_csv(tableAF2DB,sep=";") 

    cats_nbaa,allaa,nb_allSF=get_barplot_data(df)
    mk_hca_barplot(cats_nbaa,allaa,nb_allSF,"barplot_all.pdf",cats,labels,fig_title="all proteomes",save=True) # tous les protéomes combinés


    ogns=df["proteome"].unique().tolist() # list of the 21 proteomes
    
    for o in ogns:
        cats_nbaa,allaa,nb_allSF=get_barplot_data(df[df["proteome"]==o]) # On ne considère que les protéines de l'organisme o
        mk_hca_barplot(cats_nbaa,allaa,nb_allSF,"figures/"+o+"_barplot.pdf",cats,labels,fig_title=o,save=True)
        plt.close('all')
