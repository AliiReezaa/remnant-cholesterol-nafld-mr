import re, os
from docx import Document
from docx.shared import Pt, Inches, RGBColor, Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_BREAK
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
import refs
FONT='Times New Roman'
class Doc:
    def __init__(s,size=10.5):
        s.d=Document(); s.order=[]; s.size=size
        sec=s.d.sections[0]; sec.left_margin=sec.right_margin=Inches(0.9); sec.top_margin=Inches(0.9); sec.bottom_margin=Inches(0.9)
        st=s.d.styles['Normal']; st.font.name=FONT; st.font.size=Pt(size); st.element.rPr.rFonts.set(qn('w:eastAsia'),FONT)
        st.paragraph_format.space_after=Pt(6); st.paragraph_format.line_spacing=1.15
        s.fig=0; s.tab=0
    # citation numbering
    def cite(s,keys):
        nums=[]
        for k in keys:
            if k not in s.order: s.order.append(k)
            nums.append(s.order.index(k)+1)
        nums=sorted(set(nums)); out=[]; i=0
        while i<len(nums):
            j=i
            while j+1<len(nums) and nums[j+1]==nums[j]+1: j+=1
            out.append(f'{nums[i]}–{nums[j]}' if j-i>=2 else ','.join(str(n) for n in nums[i:j+1])); i=j+1
        return '['+','.join(out)+']'
    def runs(s,p,text,size=None,bold=False,italic=False,color=None):
        """markup: {k1,k2} citation (superscript), **bold**, _italic_ via ^^x^^ for superscript, ~~x~~ subscript"""
        text=re.sub(r'(?:^|(?<=[\s(\[;,]))-(?=\d)','\u2212',text)
        tokens=re.split(r'(\{[A-Za-z0-9,]+\}|\*\*.+?\*\*|\^\^.+?\^\^|~~.+?~~|__.+?__)',text)
        for t in tokens:
            if not t: continue
            if t.startswith('{') and t.endswith('}'):
                r=p.add_run(s.cite(t[1:-1].split(','))); r.font.superscript=True
            elif t.startswith('**'): r=p.add_run(t[2:-2]); r.bold=True
            elif t.startswith('^^'): r=p.add_run(t[2:-2]); r.font.superscript=True
            elif t.startswith('~~'): r=p.add_run(t[2:-2]); r.font.subscript=True
            elif t.startswith('__'): r=p.add_run(t[2:-2]); r.italic=True
            else: r=p.add_run(t)
            if bold: r.bold=True
            if italic: r.italic=True
            if size: r.font.size=Pt(size)
            if color: r.font.color.rgb=RGBColor(*color)
            r.font.name=FONT
        return p
    def para(s,text,align='justify',size=None,bold=False,italic=False,after=6,before=0,indent=None,keep=False,color=None):
        p=s.d.add_paragraph(); s.runs(p,text,size,bold,italic,color)
        p.alignment={'justify':WD_ALIGN_PARAGRAPH.JUSTIFY,'left':WD_ALIGN_PARAGRAPH.LEFT,'center':WD_ALIGN_PARAGRAPH.CENTER}[align]
        p.paragraph_format.space_after=Pt(after); p.paragraph_format.space_before=Pt(before)
        if indent: p.paragraph_format.first_line_indent=Inches(indent)
        if keep: p.paragraph_format.keep_with_next=True
        return p
    def h1(s,t): return s.para(f'**{t}**',align='left',size=13,before=12,after=6,keep=True)
    def h2(s,t): return s.para(f'**{t}**',align='left',size=11,before=8,after=4,keep=True)
    def rule(s,color='7f7f7f'):
        p=s.d.add_paragraph(); pPr=p._p.get_or_add_pPr(); b=OxmlElement('w:pBdr'); e=OxmlElement('w:bottom')
        for k,v in [('val','single'),('sz','6'),('space','1'),('color',color)]: e.set(qn('w:'+k),v)
        b.append(e); pPr.append(b); p.paragraph_format.space_after=Pt(4); return p
    def image(s,path,width=6.4):
        p=s.d.add_paragraph(); p.alignment=WD_ALIGN_PARAGRAPH.CENTER; p.paragraph_format.keep_with_next=True; p.paragraph_format.space_after=Pt(2)
        p.add_run().add_picture(path,width=Inches(width))
    def figure(s,path,title,desc,width=6.4,label=None):
        if label is None: s.fig+=1; label=f'Figure {s.fig}'
        else: s.fig+=1
        s.image(path,width)
        s.para(f'**{label}. {title}** {desc}',size=9,align='justify',after=10)
        return s.fig
    def table(s,title,header,rows,widths=None,note=None,size=8.5,boldrows=(),align_first='left',label=None):
        s.tab+=1
        if isinstance(title,tuple): label,title=title
        if label is None: label=f'Table {s.tab}'
        p=s.para(f'**{label}.** {title}',align='left',size=9.5,after=3,before=6,keep=True)
        t=s.d.add_table(rows=1+len(rows),cols=len(header)); t.alignment=WD_TABLE_ALIGNMENT.CENTER
        def setcell(c,text,bold=False,first=False):
            c.text=''; pp=c.paragraphs[0]; s.runs(pp,str(text),size=size,bold=bold); pp.paragraph_format.space_after=Pt(1); pp.paragraph_format.line_spacing=1.0
            pp.alignment=WD_ALIGN_PARAGRAPH.LEFT if first else WD_ALIGN_PARAGRAPH.CENTER
        for j,h in enumerate(header): setcell(t.rows[0].cells[j],h,True,j==0 and align_first=='left')
        for i,r in enumerate(rows):
            for j,v in enumerate(r): setcell(t.rows[i+1].cells[j],v,i in boldrows,j==0 and align_first=='left')
        def border(cell,**kw):
            tcPr=cell._tc.get_or_add_tcPr(); b=tcPr.find(qn('w:tcBorders'))
            if b is None: b=OxmlElement('w:tcBorders'); tcPr.append(b)
            for edge,val in kw.items():
                e=OxmlElement('w:'+edge)
                for k,v in [('val','single'),('sz',str(val)),('color','000000')]: e.set(qn('w:'+k),v)
                b.append(e)
        for c in t.rows[0].cells: border(c,top=12,bottom=6)
        for c in t.rows[-1].cells: border(c,bottom=12)
        if widths:
            t.autofit=False
            tblPr=t._tbl.tblPr; lay=OxmlElement('w:tblLayout'); lay.set(qn('w:type'),'fixed'); tblPr.append(lay)
            for col,w in zip(t.columns,widths): col.width=Inches(w)
            for row in t.rows:
                for c,w in zip(row.cells,widths): c.width=Inches(w)
        if note: s.para(note,size=8,align='justify',after=10,before=2)
        else: s.para('',after=4)
        return s.tab
    def pagebreak(s): s.d.add_paragraph().add_run().add_break(WD_BREAK.PAGE)
    def reflist(s,title='References'):
        s.h1(title)
        for i,k in enumerate(s.order,1):
            p=s.para(f'[{i}] {refs.fmt(k)}',align='left',size=9,after=2); p.paragraph_format.left_indent=Inches(0.35); p.paragraph_format.first_line_indent=Inches(-0.35)
    def footer(s,text):
        f=s.d.sections[0].footer; p=f.paragraphs[0]; p.text=''; r=p.add_run(text); r.font.size=Pt(8); r.font.name=FONT; p.alignment=WD_ALIGN_PARAGRAPH.CENTER
        # page number
        r2=p.add_run('   |   Page '); r2.font.size=Pt(8)
        for tag,txt in [('begin',None),(None,'PAGE'),('end',None)]:
            rr=p.add_run(); rr.font.size=Pt(8)
            if tag: e=OxmlElement('w:fldChar'); e.set(qn('w:fldCharType'),tag); rr._r.append(e)
            else: e=OxmlElement('w:instrText'); e.set(qn('xml:space'),'preserve'); e.text=txt; rr._r.append(e)
    def save(s,path): s.d.save(path)
# number formatting
def f3(x): return f'{x:.3f}'
def f2(x): return f'{x:.2f}'
def pf(p):
    if p<0.001:
        e=int(f'{p:.1e}'.split('e')[1]); m=float(f'{p:.1e}'.split('e')[0]); return f'{m:.1f} × 10^^−{abs(e)}^^'
    return ('.'+f'{p:.3f}'.split('.')[1])
def orci(o,l,h,d=3): return f'{o:.{d}f} (95% CI: {l:.{d}f}–{h:.{d}f})'
