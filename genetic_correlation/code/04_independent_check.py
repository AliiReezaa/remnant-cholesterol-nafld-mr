"""Independent cross-check of LDSC: re-implements the LD score regression estimating equations with numpy
(weighted least squares, 22-chromosome leave-one-chromosome-out jackknife) and compares h2 and rg with ldsc.py.
   E[chi2_j]       = N * h2 * l_j / M + a
   E[z1_j * z2_j]  = sqrt(N1 N2) * rho_g * l_j / M + b          (alleles aligned to trait 1)
Regression weights: 1 / (l_j * (1 + N h2_hat l_j / M)^2) (Bulik-Sullivan 2015, simplified, no IRWLS iteration beyond one pass).
"""
import pandas as pd, numpy as np, glob, sys
R='/home/claude/rg_work'
ld=pd.concat([pd.read_csv(f,sep='\t') for f in sorted(glob.glob(f'{R}/ref/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.*.l2.ldscore.gz'))])[['CHR','SNP','L2']]
M=len(ld)
def load(n): return pd.read_csv(f'{R}/munged/{n}.sumstats.gz',sep='\t').dropna(subset=['Z'])
def wls(y,x,w):
    X=np.column_stack([np.ones_like(x),x]); W=w[:,None]
    b=np.linalg.solve(X.T@(W*X),X.T@(W[:,0]*y)); return b
def jack(y,x,w,chrom,N_scale):
    b=wls(y,x,w); slopes=[]
    for c in range(1,23):
        m=chrom!=c; slopes.append(wls(y[m],x[m],w[m])[1])
    slopes=np.array(slopes); k=len(slopes); pseudo=k*b[1]-(k-1)*slopes
    se=pseudo.std(ddof=1)/np.sqrt(k)
    return b,se,pseudo
def h2(df,N):
    y=df.Z.values**2; x=df.L2.values
    b=wls(y,x,1/x); hhat=max(b[1]*M/N,0.001)
    w=1/(x*(1+N*hhat*x/M)**2)
    return wls(y,x,w),w
def run(a,b):
    A=load(a); B=load(b); df=ld.merge(A[['SNP','A1','A2','Z','N']],on='SNP').merge(B[['SNP','A1','A2','Z','N']],on='SNP',suffixes=('1','2'))
    flip=(df.A11==df.A12)&(df.A21==df.A22)
    swap=(df.A11==df.A22)&(df.A21==df.A12)
    ok=flip|swap
    df=df[ok].copy(); df['Z2']=np.where(swap[ok],-df.Z2,df.Z2)   # same-strand data: swapped alleles => flip sign
    df=df[(df.Z1**2<max(80,0.001*df.N1.max()))&(df.Z2**2<max(80,0.001*df.N2.max()))]
    N1,N2=df.N1.mean(),df.N2.mean(); x=df.L2.values; chrom=df.CHR.values
    b1,w1=h2(df.assign(Z=df.Z1),N1); b2,w2=h2(df.assign(Z=df.Z2),N2)
    h21=b1[1]*M/N1; h22=b2[1]*M/N2
    y=(df.Z1*df.Z2).values
    w=1/(x*(1+N1*h21*x/M)**0.5*(1+N2*h22*x/M)**0.5)  # approximate gencov weights
    bg=wls(y,x,w); gcov=bg[1]*M/np.sqrt(N1*N2)
    # jackknife rg
    rgs=[]
    for c in range(1,23):
        m=chrom!=c
        bb1=wls(df.Z1.values[m]**2,x[m],w1[m]); bb2=wls(df.Z2.values[m]**2,x[m],w2[m]); bgc=wls(y[m],x[m],w[m])
        rgs.append((bgc[1]*M/np.sqrt(N1*N2))/np.sqrt(max(bb1[1]*M/N1,1e-6)*max(bb2[1]*M/N2,1e-6)))
    rgs=np.array(rgs); rg=gcov/np.sqrt(h21*h22); k=22; pseudo=k*rg-(k-1)*rgs
    print(f'{a} vs {b}: SNPs {len(df):,}  h2_1 {h21:.4f}  h2_2 {h22:.4f}  gcov {gcov:.4f}  rg {rg:.3f} (jackknife SE {pseudo.std(ddof=1)/np.sqrt(k):.3f})  intercepts: {b1[0]:.3f} {b2[0]:.3f} gencov {bg[0]:.4f}')
for a,b in [('RC','RC'),('RC','CAD'),('RC','Ghodsian'),('RC','Fairfield'),('Ghodsian','Fairfield')]: run(a,b)
