import re, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, numpy as np
d={}; so={}
for line in open('../lis/output/job17_fig.txt'):
    t=line.split()
    if t[0]=='FIG': so[(t[1],t[2])]=float(t[6])
    if t[0]=='FIGSHARES':
        xs=[];ys=[]
        for kv in t[3:]:
            k,v=kv.split(':'); xs.append(float(k)); ys.append(np.nan if v=='NA' else float(v))
        d[(t[1],t[2])]=(np.array(xs),np.array(ys))
fig,axs=plt.subplots(1,2,figsize=(9,3.4),sharey=True)
FP={'dk10':133.272,'dk19':159.0}
for ax,(nm,yr) in zip(axs,[('dk10','2010'),('dk19','2019')]):
    x=d[(nm,'old')][0]
    parts=[('r2_nohb','transfer-reliant, no housing benefit','#4d4d4d'),('r2_hb','transfer-reliant, housing benefit','#9e9e9e'),('old_other','other','#d9d9d9')]
    bottom=np.zeros_like(x)
    for g,lab,c in parts:
        y=np.nan_to_num(d[(nm,g)][1])*so[(nm,g)]*100
        ax.bar(x,y,width=5,bottom=bottom,color=c,edgecolor='white',linewidth=.3,label=lab); bottom+=y
    ax.axvline(FP[nm],color='k',ls=':',lw=1)
    ax.set_xlim(85,260); ax.set_title(f'Denmark {yr}',fontsize=10); ax.set_xlabel('disposable income, thousand DKK',fontsize=9)
    ax.tick_params(labelsize=8)
axs[0].set_ylabel('per cent of people 65+ living alone',fontsize=9)
axs[1].legend(fontsize=7,frameon=False,loc='upper right')
plt.tight_layout(); plt.savefig('fig/dk_stack.pdf'); plt.savefig('fig/dk_stack.png',dpi=130)
