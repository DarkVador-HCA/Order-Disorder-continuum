import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

#### Figures parameters
#COLORS
ci={"SCOPe":5.5,"DisProt":9,"OPM":2.5,
    "a":3,"b":5,'c':6.8,'d':8,
    "1":3,"2":2,"11":0,
    'Disorder function':5,'Interaction partner':6.2,'Structural state':7.8, 'Structural transition':8.5,
    'Interaction partner':5,'Other':7.8,
    "OPM_1_11":1,"OPM_2":3
   }
cm = plt.cm.get_cmap('viridis')#,vmin=0,vmax=8)
norm=mpl.colors.Normalize(vmin=0,vmax=9)
cms={"SCOPe":plt.cm.get_cmap('viridis'),"OPM":plt.cm.get_cmap('viridis'),"DisProt":plt.cm.get_cmap('viridis')}


# LABELS
labels={"SCOPe":"SCOPe (a,b,c,d)","DisProt":"DisProt v8.0.2","OPM":"OPM",
    "a":"SCOPe a","b":"SCOPe b",'c':"SCOPe c",'d':"SCOPe d",
    "1":"Alpha-helical polytopic","2":"Beta-barrel transmembrane","11":"Bitopic proteins",
    'Structural state':'Structural state', 'Structural transition':'Structural transition', 'Disorder function':'Disorder function', 'Interaction partner':'Interaction partner', 'Other':'Other',
    "OPM_1_11":"OPM (alpha-helical polytopic and bitopic)","OPM_2":"OPM (beta-barrel)",
    }

### Function
def plot_density(distribs,labels,lab_dict,colors,xlabel,cm=plt.cm.get_cmap('viridis'),fix_xlim=False,bins=0.4,add_hist=False,legend=True,figsize=None,constrained_layout=True):
    plt.figure(figsize=figsize,constrained_layout=constrained_layout)
    if fix_xlim:
        plt.xlim(fix_xlim)
    # cm = plt.cm.get_cmap('viridis')#,vmin=0,vmax=8)
    norm=mpl.colors.Normalize(vmin=0,vmax=9)
    densities=[]
    for i in range(len(distribs)):
        # density = gaussian_kde(distribs[i],weights=[1/len(distribs[i]) for _ in range(len(distribs[i]))])
        if add_hist:
            h=[]
            for d in add_hist[0]:
                h+=d
            w=[]
            for d in add_hist[1]:
                w+=d
            plt.hist(h,color="red",weights=w,bins=5) #Train-Seqs=0.19 #,bins=3    ["red" for _ in range(len(add_hist[0]))]
            plt.hist(add_hist[0],color=cm(norm(colors[labels[i]])),weights=add_hist[1],bins=3) #Train-Seqs=0.19 #,bins=3

        density = gaussian_kde(distribs[i])
        xs = np.linspace(min(distribs[i]),max(distribs[i]),200) # v1 & v2 & v5 : 200 # v3 & v4:1000
        density.covariance_factor = lambda : .25 # v1 & v3 : .25 # v2 & v4 : Sans  # v5 : .05
        density._compute_covariance() #v1 & v3 : Avec # v2 & v4 : Sans
        densities.append([xs,density(xs)])
        plt.plot(xs,density(xs),color=cm(norm(colors[labels[i]])))#,color=colsD[i])
        plt.xticks(range(-7,10),range(-7,10))
        plt.ylabel("density")
        #
        plt.fill_between(xs,0,density(xs),color=cm(norm(colors[labels[i]])),alpha=0.3,label=lab_dict[labels[i]],zorder=5)
    if legend:
        leg=plt.legend()
        for i,lh in enumerate(leg.legendHandles):
            lh.set_edgecolor(cm(norm(colors[labels[i]])))
    plt.xlabel(xlabel)
    return densities

#### Prep data
def find_limit(distrib1,distrib2,bins=60,lim=[-10,10]):
    intervals = np.linspace(lim[0], lim[1], bins)
    valZ1, b1 = np.histogram(distrib1, bins=intervals)
    valZ2, b2 = np.histogram(distrib2, bins=intervals)
    valZ1=np.true_divide(valZ1,np.sum(valZ1))
    valZ2=np.true_divide(valZ2,np.sum(valZ2))
    i=0
    while valZ1[i]>=valZ2[i]:
        i+=1
    l=b1[i-1]
    return l

def all_fold(df):
    # DBs=df['DB'].unique()
    toplot=[]
    toplot.append(df[df["database"]=="SCOPe"]["HCA score"].values)
    toplot.append(df[(df["class"]=="1") | (df["class"]=="11")]["HCA score"].values) #OPM alpha
    toplot.append(df[df["class"]=="2"]["HCA score"].values) #OPM beta
    toplot.append(df[df["database"]=="DisProt"]["HCA score"].values)

    # #HCA scores to plot
    # print('Only segments from DisProt under : ',min(df[(df["database"]!="DisProt")]["HCA score"]))
    # print('5th percentile for SCOPe segments:',np.percentile(toplot[0],5))
    # print('OPM alpha and SCOPe distributions cross at :',find_limit(toplot[0],toplot[1],bins=150,lim=[-2.5,10]))
    # print('95th percentile for OPM alpha segments :',np.percentile(toplot[1],95))
    # print('HCA score minimum for segments with just 1 hydrophobic cluster :',min(df[df["number of hydrophobic clusters"]<2]["HCA score"]))
    # print('Only segments with just 1 hydrophobic cluster from :',max(df[df["number of hydrophobic clusters"]>2]["HCA score"]))

    groups=["SCOPe","OPM_1_11","OPM_2","DisProt"]
    densities=plot_density(toplot,groups,labels,ci,"HCA score",fix_xlim=(-7,9),legend=False,figsize=None,constrained_layout=True)
    leg=plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.2), ncol=1,framealpha=1)
    for i,lh in enumerate(leg.legendHandles):
        lh.set_edgecolor(cm(norm(ci[groups[i]])))
    ymin,ymax=plt.ylim()
    # print(np.argmax(densities[-1]))
    # print(len(densities[-1]))
    l=[-4.7,-1,3.5,6.7,7.7]

    plt.vlines(l,ymin,ymax,ls="--",color="grey",zorder=5,lw=1)
    plt.show()
    plt.close("all")
