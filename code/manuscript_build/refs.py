"""Reference registry: key -> DOI. Metadata fetched from Crossref and cached (refs_cache.json); Vancouver strings built from it."""
import json, os, re, subprocess
HERE=os.path.dirname(__file__); CACHE=os.path.join(HERE,'refs_cache.json')
DOI={
'varbo13':'10.1016/j.jacc.2012.08.1026','nord16':'10.1161/CIRCRESAHA.115.306249','nordvarbo14':'10.1016/S0140-6736(14)61177-6',
'quispe21':'10.1093/eurheartj/ehab432','duran20':'10.1016/j.jacc.2020.02.059','richardson20':'10.1371/journal.pmed.1003062',
'younossi16':'10.1002/hep.28431','younossi23':'10.1097/HEP.0000000000000004','riazi22':'10.1016/S2468-1253(22)00165-0',
'rinella23':'10.1097/HEP.0000000000000520','eslam20':'10.1016/j.jhep.2020.03.039','chalasani18':'10.1002/hep.29367',
'romeo08':'10.1038/ng.257','kozlitina14':'10.1038/ng.2901','speliotes11':'10.1371/journal.pgen.1001324','fairfield22':'10.1002/hep4.1805',
'ghodsian21':'10.1016/j.xcrm.2021.100437','donnelly05':'10.1172/JCI23621','fabbrini10':'10.1002/hep.23280','fabbrini09':'10.1073/pnas.0904944106',
'zou21':'10.1186/s12944-021-01573-y','miao23':'10.1136/bmjopen-2022-069440','liu24':'10.1265/ehpm.23-00289','wang24':'10.1038/s41598-024-55282-4',
'guan23':'10.1186/s12933-023-01927-z','ko25':'10.3390/genes16020157',
'ds03':'10.1093/ije/dyg070','ds14':'10.1093/hmg/ddu328','sanderson22':'10.1038/s43586-021-00092-5','lawlor08':'10.1002/sim.3034',
'burgess13':'10.1002/gepi.21758','burgess24':'10.1186/s12916-024-03607-5','burgess20g':'10.12688/wellcomeopenres.15555.2',
'skriv21':'10.1001/jama.2021.18236','skrivbmj21':'10.1136/bmj.n2233',
'richardson22':'10.1371/journal.pbio.3001547','sudlow15':'10.1371/journal.pmed.1001779','bycroft18':'10.1038/s41586-018-0579-z',
'soininen15':'10.1161/CIRCGENETICS.114.000216','wurtz17':'10.1093/aje/kwx016','loh15':'10.1038/ng.3190','auton15':'10.1038/nature15393',
'kurki23':'10.1038/s41586-022-05473-8','sollis23':'10.1093/nar/gkac1010','nikpay15':'10.1038/ng.3396',
'hemani18':'10.7554/eLife.34408','hemani17':'10.1371/journal.pgen.1007081','bowden15':'10.1093/ije/dyv080','bowden16':'10.1002/gepi.21965',
'hartwig17':'10.1093/ije/dyx102','verbanck18':'10.1038/s41588-018-0099-7','bowden17':'10.1002/sim.7221','bowdenI2':'10.1093/ije/dyw220',
'burgess11':'10.1093/ije/dyr036','burgess16':'10.1002/gepi.21998','mounier23':'10.1002/gepi.22522','brion13':'10.1093/ije/dyt179',
'burgess14':'10.1093/ije/dyu005','slob20':'10.1002/gepi.22295','conmix20':'10.1038/s41467-019-14156-4','rees19':'10.1371/journal.pone.0222362',
'zhao20':'10.1214/19-AOS1866','ye21':'10.1214/20-AOS2027','yavorska17':'10.1093/ije/dyx034','bowden18':'10.1093/ije/dyy101',
'sanderson19':'10.1093/ije/dyy262','burgess15':'10.1093/aje/kwu283','burgess17s':'10.1097/EDE.0000000000000559','burgess17e':'10.1007/s10654-017-0255-x',
'stender24':'10.1186/s12944-024-02284-w','elsworth20':'10.1101/2020.08.10.244293','sanderson21':'10.1002/sim.9133',
}
JAB={'International Journal of Epidemiology':'Int J Epidemiol','Genetic Epidemiology':'Genet Epidemiol','Nature Genetics':'Nat Genet','PLoS Genetics':'PLoS Genet','PLOS Genetics':'PLoS Genet','PLOS ONE':'PLoS One','Journal of the American College of Cardiology':'J Am Coll Cardiol','European Heart Journal':'Eur Heart J','Circulation Research':'Circ Res','The Lancet':'Lancet','The Lancet Gastroenterology & Hepatology':'Lancet Gastroenterol Hepatol','Journal of Hepatology':'J Hepatol','Journal of Clinical Investigation':'J Clin Invest','Proceedings of the National Academy of Sciences':'Proc Natl Acad Sci U S A','Nature Communications':'Nat Commun','Nucleic Acids Research':'Nucleic Acids Res','Statistics in Medicine':'Stat Med','The Annals of Statistics':'Ann Stat','PLoS Biology':'PLoS Biol','PLOS Biology':'PLoS Biol','PLoS Medicine':'PLoS Med','PLOS Medicine':'PLoS Med','Cell Reports Medicine':'Cell Rep Med','Hepatology Communications':'Hepatol Commun','Lipids in Health and Disease':'Lipids Health Dis','Scientific Reports':'Sci Rep','Environmental Health and Preventive Medicine':'Environ Health Prev Med','BMC Medicine':'BMC Med','Cardiovascular Diabetology':'Cardiovasc Diabetol','Circulation: Cardiovascular Genetics':'Circ Cardiovasc Genet','American Journal of Epidemiology':'Am J Epidemiol','Wellcome Open Research':'Wellcome Open Res','Nature Reviews Methods Primers':'Nat Rev Methods Primers','European Journal of Epidemiology':'Eur J Epidemiol','Human Molecular Genetics':'Hum Mol Genet','BMJ Open':'BMJ Open','Journal of the American Medical Association':'JAMA','JAMA':'JAMA'}
def _fetch(doi):
    import time
    for i in range(6):
        out=subprocess.run(['curl','-sS','-m','25',f'https://api.crossref.org/works/{doi}'],capture_output=True,text=True).stdout
        try: return json.loads(out)['message']
        except Exception: time.sleep(2+i*2)
    raise RuntimeError('crossref failed '+doi)
def meta(key):
    c=json.load(open(CACHE)) if os.path.exists(CACHE) else {}
    if key not in c:
        m=_fetch(DOI[key]); c[key]={k:m.get(k) for k in ['author','title','container-title','issued','volume','issue','page','article-number','DOI','type','institution','subtitle']}
        json.dump(c,open(CACHE,'w'),indent=1)
    return c[key]
KEEP={'Mendelian','European','Biobank','UK','LDL','ApoB','NAFLD','NASH','MR','PNPLA3','TM6SF2','FinnGen','GWAS','I2','CAD','NMR','R','PLINK','Egger','MendelianRandomization','American','Japanese','Chinese','Taiwan'}
def sent(t):
    w=t.split(' ')
    caps=sum(1 for x in w if len(x)>3 and x[0].isupper())
    if caps<0.5*max(1,sum(1 for x in w if len(x)>3)): return t
    out=[]
    for i,x in enumerate(w):
        core=re.sub(r'[^A-Za-z0-9]','',x)
        if i==0 or core.isupper() or core in KEEP or re.search(r'\d',core) or x[:1]==x[:1].lower(): out.append(x)
        else:
            out.append(x[0].lower()+x[1:] if '-' not in x else '-'.join(p.lower() if j==0 or p not in KEEP else p for j,p in enumerate(x.split('-'))))
    return ' '.join(out).replace('‐','-')
import html
OVR={'nikpay15':dict(authors='Nikpay M, Goel A, Won HH, et al.; CARDIoGRAMplusC4D Consortium'),
'auton15':dict(authors='Auton A, Brooks LD, Durbin RM, et al.; 1000 Genomes Project Consortium'),
'bowdenI2':dict(page='1961–1974',vol='45',year=2016),'skrivbmj21':dict(vol='375',page='n2233'),
'brion13':dict(year=2013),'liu24':dict(page='5'),'fairfield22':dict(year=2022),'sollis23':dict(year=2023),'sanderson19':dict(year=2019),
'skriv21':dict(page='1614–1621'),'burgess11':dict(authors='Burgess S, Thompson SG'),'riazi22':dict(journal='Lancet Gastroenterol Hepatol'),
'chalasani18':dict(year=2018),'elsworth20':dict(year=2020)}
def fmt(key,nmax=3):
    m=meta(key); au=m.get('author') or []
    def nm(a):
        fam=a.get('family',''); g=a.get('given','')
        ini=''.join(w[0] for w in re.split(r'[\s\-\.]+',g) if w)
        return f'{fam} {ini}'.strip()
    names=[nm(a) for a in au[:nmax]]
    s=', '.join(names)+(', et al' if len(au)>nmax else '')
    s=OVR.get(key,{}).get('authors',s).rstrip('.')
    title=re.sub(r'<[^>]+>','',m['title'][0]).rstrip('.')
    if m.get('subtitle'): title+=': '+m['subtitle'][0].rstrip('.')
    title=html.unescape(title).replace('‐','-').replace('*','').strip()
    title=sent(title)
    for a,b in [('mendelian','Mendelian'),('Genome-Wide','Genome-wide'),('Triglyceride-Rich','Triglyceride-rich'),('Meta-analytic','meta-analytic'),('-Omic','-omic'),('mendelian','Mendelian'),('1000 Genomes–based','1000 Genomes-based'),('Mendelian randomisation','Mendelian randomisation')]: title=title.replace(a,b)
    o=OVR.get(key,{})
    j=html.unescape((m.get('container-title') or [''])[0]); j=JAB.get(j,j); j=OVR.get(key,{}).get('journal',j)
    o=OVR.get(key,{}); yr=o.get('year',m['issued']['date-parts'][0][0]); vol=o.get('vol',m.get('volume')); pg=o.get('page',m.get('page') or m.get('article-number'))
    if m.get('type')=='posted-content': j='bioRxiv (preprint)'
    title=title.replace('?.','?')
    ref=f'{s}. {title}. {j}. {yr}'.replace('. . ','. ').replace('M. F','M F')
    if vol: ref+=f';{vol}'
    if pg: ref+=f':{pg.replace("-","–")}'
    ref+=f'. doi:{DOI[key]}' if not (vol and pg) else '.'
    return ref
if __name__=='__main__':
    for k in DOI: print(k,'|',fmt(k))
