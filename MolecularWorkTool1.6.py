import tkinter.messagebox as msgbox
import tkinter.ttk as ttk
from tkinter import*

root=Tk()

root.title('Molecular Work Tool')
root.geometry("1250x600+50+50")

Label(root, text='Hello!').grid(row=0,column=2)
wframe=LabelFrame(root, text='Which work do you want to do?')
wframe.grid(row=1, column=0)
fframe=LabelFrame(root, text='Sequence (or subject sequence)')
fframe.grid(row=2, column=2)
sframe=LabelFrame(root, text='Sequence you want to find (or query sequence)')
sframe.grid(row=4, column=2)
tframe=LabelFrame(root, text='Result')
tframe.grid(row=2, column=5,rowspan=3, padx=10)
Label(root, text='Thank you :D').grid(row=8,column=2)

eframe=Frame(root)
aframe=Frame(root)
lframe=Frame(root)

codon3={'TTT':'Phe','TTC':'Phe','TTA':'Leu','TTG':'Leu','TCT':'Ser','TCC':'Ser','TCA':'Ser','TCG':'Ser','TAT':'Tyr','TAC':'Tyr','TAA':'___','TAG':'___','TGT':'Cys','TGC':'Cys','TGA':'___','TGG':'Trp','CTT':'Leu','CTC':'Leu','CTA':'Leu','CTG':'Leu','CCT':'Pro','CCC':'Pro','CCA':'Pro','CCG':'Pro','CAT':'His','CAC':'His','CAA':'Gln','CAG':'Gln','CGT':'Arg','CGC':'Arg','CGA':'Arg','CGG':'Arg','ATT':'Ile','ATC':'Ile','ATA':'Ile','ATG':'Met','ACT':'Thr','ACC':'Thr','ACA':'Thr','ACG':'Thr','AAT':'Asn','AAC':'Asn','AAA':'Lys','AAG':'Lys','AGT':'Ser','AGC':'Ser','AGA':'Arg','AGG':'Arg','GTT':'Val','GTC':'Val','GTA':'Val','GTG':'Val','GCT':'Ala','GCC':'Ala','GCA':'Ala','GCG':'Ala','GAT':'Asp','GAC':'Asp','GAA':'Glu','GAG':'Glu','GGT':'Gly','GGC':'Gly','GGA':'Gly','GGG':'Gly', 'CUN':'Leu', 'GUN':'Val', 'UCN':'Ser', 'CCN':'Pro', 'ACN':'Thr', 'GCN':'Ala', 'CGN':'Arg', 'GGN':'Gly'}
threeto1={'Ala':'A', 'Ser':'S','Cys':'C','Thr':'T','Val':'V','Met':'M','Trp':'W','Tyr':'Y','Asn':'N','Asp':'D','Gln':'Q','Glu':'E','His':'H','Gly':'G','Pro':'P','Lys':'K','Leu':'L','Ile':'I','Phe':'F','Arg':'R', 'Pyl':'O', 'Sec':'U', 'Asx':'B', 'Glx':'Z', 'Xle':'J', 'Xaa':'X', '___':'_'}
oneto3={'A':'Ala', 'S':'Ser','C':'Cys','T':'Thr','V':'Val','M':'Met','W':'Trp','Y':'Tyr','N':'Asn','D':'Asp','Q':'Gln','E':'Glu','H':'His','G':'Gly','P':'Pro','K':'Lys','L':'Leu','I':'Ile','F':'Phe','R':'Arg', 'O':'Pyl', 'U':'Sec', 'B':'Asx', 'Z':'Glx', 'J':'Xle', 'X':'Xaa', '_':'___'}
base=['A','T','U','C','G','N']
aaMW={'A':89.09, 'R':174.2, 'N': 132.12, 'D': 133.1, 'C':121.16, 'Q':146.14, 'E':147.13, 'Z':146.65, 'G':75.07, 'H':155.15, 'I':131.17, 'L':131.17, 'K':146.19, 'M':149.21, 'F':165.19, 'P':115.13, 'S':105.09, 'T':119.12, 'W':204.22, 'Y':181.19, 'V':117.15, 'B': 132.61, 'J':131.17, 'Z':146.635, 'O':255.31, 'U':167.06, 'X':110}

work=['Counting DNA length(+GC ratio)', 'Finding DNA seq location','RNA,DNA moles <-> mass','Generating complementray seq','Generating reverse seq','Translation', 'Simple DNA seq alignment', 'Counting a.a. seq length(+Mw)','Finding a.a. seq location', 'Converting a.a. abbreviation']
cb=ttk.Combobox(wframe, width=27, height=5, values=work, state='readonly')
cb.set('Please choose')
cb.grid(row=0,column=0)

letter=['1 letter','3 letters']
cb1=ttk.Combobox(wframe, height=2, values=letter, state='readonly')
cb1.set('a.a. abbreviation type')

DNAform=['separately', 'dsDNA']
cb2=ttk.Combobox(wframe, height=2, values=DNAform, state='readonly')
cb2.set('How to be shown')

findingmethod=['location of seq', 'seq of the location']
cb3=ttk.Combobox(wframe,height=2, values=findingmethod, state='readonly')
cb3.set('What to find')
esframe=LabelFrame(eframe,text='start')
esframe.grid(row=0,column=0)
etframe=LabelFrame(eframe,text='end')
etframe.grid(row=1,column=0)
es=ttk.Entry(esframe, width=8)
es.pack()
et=ttk.Entry(etframe, width=8)
et.pack()

caltype=['dsDNA','ssDNA','ssRNA']
calmethod=['mass -> moles','moles -> mass']
caltcb=ttk.Combobox(wframe,height=3, values=caltype, state='readonly')
caltcb.set('nucleic acid type')
calmcb=ttk.Combobox(wframe,height=2, values=calmethod, state='readonly')
calmcb.set('What to convert')
lnframe=LabelFrame(lframe,text='length(bp)')
lnframe.grid(row=0,column=0)
ln=ttk.Entry(lnframe,width=10)
ln.pack()
moleframe=LabelFrame(lframe,text='moles(pmol)')
moleframe.grid(row=2,column=0)
mole=ttk.Entry(moleframe,width=8)
mole.pack()
massframe=LabelFrame(lframe,text='mass(ng)')
massframe.grid(row=3,column=0)
mass=ttk.Entry(massframe,width=8)
mass.pack()


convert=['1 to 3', '3 to 1']
cb4=ttk.Combobox(wframe,height=2, values=convert, state='readonly')
cb4.set('What to convert')

matframe=LabelFrame(aframe,text='match score')
matframe.grid(row=0,column=0)
mismatframe=LabelFrame(aframe,text='mismatch score')
mismatframe.grid(row=1,column=0)
gapframe=LabelFrame(aframe,text='gap penalty')
gapframe.grid(row=2,column=0)
match=ttk.Entry(matframe,width=8)
match.pack()
mismatch=ttk.Entry(mismatframe,width=8)
mismatch.pack()
gap=ttk.Entry(gapframe,width=8)
gap.pack()

sb1=Scrollbar(fframe)
sb1.pack(side='right',fill='y')
t=Text(fframe, width=54, height=20, yscrollcommand=sb1.set)
t.insert(END,'Please enter DNA or RNA or a.a. seq here :D\n')
t.insert(END,'(case-insensitive)\n')
t.insert(END,'\n')
t.insert(END,"Numbers will be deleted automatically\n")
t.insert(END,'RNA seq will be converted to DNA seq automatically\n')
# t.insert(END,'a.a. seq 1 letter 약자로 입력해주세요\n')
t.insert(END,'\n')
t.insert(END,'\n')
# t.insert(END,'원하시는 작업을 고르신 후 선택 버튼을 누르고 Clear를 눌러주세요')
t.pack(side='left')
sb1.config(command=t.yview)

def clear():
    t.delete('1.0',END)
    tt.delete('1.0',END)

def a():
    global pseq
    global pseqrr
    p.config(state=NORMAL)
    p.delete('1.0',END)
    i=t.get('1.0',END)

    global seq1
    tem3=[z for z in i if z.isalpha()==True]
    i=''.join(tem3)
    
    if w!=9:
        i=i.upper()
    tem1=[]
    for z in range(len(i)):
        if i[z]=='U':
            tem1.append(z)
    tem2=[z for z in i]
    for z in tem1:
        tem2[z]='T'
    i=''.join(tem2)
    seq1=i.split()
    seq1=''.join(seq1)

    rtem=[]
    for z in range(len(i)):
        if i[z]=='A':
            rtem.append('T')
        if i[z]=='T':
            rtem.append('A')
        if i[z]=='G':
            rtem.append('C')
        if i[z]=='C':
            rtem.append('G')
        if i[z]=='N':
            rtem.append('N')
    rr=''.join(rtem)
    rtem=[]
    for z in range(len(rr)-1,-1,-1):
        rtem.append(rr[z])
    r=''.join(rtem)
    global seq1rr
    seq1rr=rr
    global seq1r
    seq1r=r

    global gcr1
    gn1=seq1.count('G')
    cn1=seq1.count('C')
    if len(seq1)!=0:
        gcr1=((gn1+cn1)/len(seq1))*100

    global seq2
    i=tt.get('1.0',END)
    
    tem3=[z for z in i if z.isalpha()==True]
    i=''.join(tem3)

    if w!=8:
        i=i.upper()
    tem1=[]
    for z in range(len(i)):
        if i[z]=='U':
            tem1.append(z)
    tem2=[z for z in i]
    for z in tem1:
        tem2[z]='T'
    i=''.join(tem2)
    seq2=i.split()
    seq2=''.join(seq2)

    rtem=[]
    for z in range(len(i)):
        if i[z]=='A':
            rtem.append('T')
        if i[z]=='T':
            rtem.append('A')
        if i[z]=='G':
            rtem.append('C')
        if i[z]=='C':
            rtem.append('G')
        if i[z]=='N':
            rtem.append('N')
    rr=''.join(rtem)
    rtem=[]
    for z in range(len(rr)-1,-1,-1):
        rtem.append(rr[z])
    r=''.join(rtem)
    seq2rr=rr
    seq2r=r

    global gcr2
    gn2=seq2.count('G')
    cn2=seq2.count('C')
    if len(seq2)!=0:
        gcr2=((gn2+cn2)/len(seq2))*100

    if w==0 or w==1 or w==2 or w==3 or w==4 or w==5 or w==6:
        error1=[]
        error2=[]
        for i in seq1:
            if i in base: continue
            else:
                error1.append(i)
        if len(error1)>0:
            msgbox.showerror('Error', 'There is character that is not base in your sequence!')
        for i in seq2:
            if i in base: continue
            else:
                error2.append(i)
        if len(error2)>0:
            msgbox.showerror('Error','There is character that is not base in the sequence you want to find!')
            
        
    if w==0:
        p.insert(END, 'The length of DNA seq\n')
        p.insert(END, '\n')
        p.insert(END, len(seq1))
        p.insert(END, '\n\n\n\n')
        p.insert(END, 'GC content\n')
        p.insert(END, '\n')
        p.insert(END, gcr1)
        p.insert(END, ' %')
        p.config(state=DISABLED)

    if w==1:
        l=cb3.get()
        l=findingmethod.index(l)
        if l==0:
            if seq2 in seq1:
                nos=seq1.count(seq2)
                p.insert(END, int(seq1.index(seq2))+1)
                p.insert(END, '\n')
                p.insert(END, '~')
                p.insert(END, '\n')
                p.insert(END, int(seq1.index(seq2))+len(seq2))
                p.insert(END, '\n')
                p.insert(END, '\n')
                underb=''
                for i in range(len(seq2)):
                    underb=underb+'_'
                temseq1=seq1.replace(seq2,underb,1)
                if nos>=2:
                    for _ in range(nos-1):
                        p.insert(END, int(temseq1.index(seq2))+1)
                        p.insert(END, '\n')
                        p.insert(END, '~')
                        p.insert(END, '\n')
                        p.insert(END, int(temseq1.index(seq2))+len(seq2))
                        p.insert(END, '\n')
                        p.insert(END, '\n')
                        temseq1=temseq1.replace(seq2,underb,1)
                p.config(state=DISABLED)
                
            else:
                p.insert(END, 'Cannot find')
                p.config(state=DISABLED)
        
        if l==1:
            esg=es.get()
            etg=et.get()
            fseq=seq1[int(esg)-1:int(etg)]
            pseq=[]
            i=0
            while 10*i+10<len(fseq):
                pseq.append(fseq[10*i:10*i+10])
                i+=1
            else:
                pseq.append(fseq[10*i:len(seq1)])
            j=0
            while j<len(pseq)//5:
                p.insert(END, pseq[5*j:5*j+5])
                p.insert(END, '\n')
                j+=1
            p.insert(END, pseq[5*j:len(pseq)])
            p.insert(END, '\n')
            p.config(state=DISABLED)

    if w==2:
        na=seq1.count('A')
        nt=seq1.count('T')
        nc=seq1.count('C')
        ng=seq1.count('G')
        l=ln.get()
        if l!='':
            l=float(ln.get())
        else: l=0
        calt=caltype.index(caltcb.get())
        calm=calmethod.index(calmcb.get())
        moleg=mole.get()
        massg=mass.get()
        if moleg!='':
            moleg=float(mole.get())
        if massg!='':
            massg=float(mass.get())
        
        if calt==0:
            if l>0:
                l=float(l)
                seqmw=617.96*l+36.04
            else:
                seqmw=(na+nt)*(313.23+304.21)+(nc+ng)*(329.23+289.2)+36.04
            
            if calm==0:
                mo=(massg/seqmw)*1000 #pmol
                p.insert(END, '\n')
                p.insert(END, 'Molecular weight\n')
                p.insert(END, '\n')
                p.insert(END, seqmw)
                p.insert(END, ' Da\n')
                p.insert(END, '\n')
                p.insert(END, 'Moles of dsDNA\n')
                p.insert(END, '\n')
                p.insert(END, mo)
                p.insert(END, ' pmol')
                p.config(state=DISABLED)

            if calm==1:
                ma=(moleg*seqmw)/1000 #ng
                p.insert(END, '\n')
                p.insert(END, 'Molecular weight\n')
                p.insert(END, '\n')
                p.insert(END, seqmw)
                p.insert(END, ' Da\n')
                p.insert(END, '\n')
                p.insert(END, 'dsDNA Mass\n')
                p.insert(END, '\n')
                p.insert(END, ma)
                p.insert(END, ' ng')
                p.config(state=DISABLED)

        if calt==1:
            if l>0:
                l=float(l)
                seqmw=308.97*l+18.02
            else:
                seqmw=na*313.23+nt*304.21+ng*329.23+nc*289.2+18.02
            if calm==0:
                mo=(massg/seqmw)*1000 #pmol
                p.insert(END, '\n')
                p.insert(END, 'Molecular weight\n')
                p.insert(END, '\n')
                p.insert(END, seqmw)
                p.insert(END, ' Da\n')
                p.insert(END, '\n')
                p.insert(END, 'Moles of ssDNA\n')
                p.insert(END, '\n')
                p.insert(END, mo)
                p.insert(END, ' pmol')
                p.config(state=DISABLED)

            if calm==1:
                ma=(moleg*seqmw)/1000 #ng
                p.insert(END, '\n')
                p.insert(END, 'Molecular weight\n')
                p.insert(END, '\n')
                p.insert(END, seqmw)
                p.insert(END, ' Da\n')
                p.insert(END, '\n')
                p.insert(END, 'ssDNA Mass\n')
                p.insert(END, '\n')
                p.insert(END, ma)
                p.insert(END, ' ng')
                p.config(state=DISABLED)

        if calt==2:
            if l>0:
                l=float(l)
                seqmw=321.47*l+18.02
            else:
                seqmw=na*329.2+nt*306.2+ng*345.2+nc*305.2+18.02
            if calm==0:
                mo=(massg/seqmw)*1000 #pmol
                p.insert(END, '\n')
                p.insert(END, 'Molecular weight\n')
                p.insert(END, '\n')
                p.insert(END, seqmw)
                p.insert(END, ' Da\n')
                p.insert(END, '\n')
                p.insert(END, 'Moles of ssRNA\n')
                p.insert(END, '\n')
                p.insert(END, mo)
                p.insert(END, ' pmol')
                p.config(state=DISABLED)

            if calm==1:
                ma=(moleg*seqmw)/1000 #ng
                p.insert(END, '\n')
                p.insert(END, 'Molecular weight\n')
                p.insert(END, '\n')
                p.insert(END, seqmw)
                p.insert(END, ' Da\n')
                p.insert(END, '\n')
                p.insert(END, 'ssRNA Mass\n')
                p.insert(END, '\n')
                p.insert(END, ma)
                p.insert(END, ' ng')
                p.config(state=DISABLED)

    if w==3:
        pseq=[]
        pseqrr=[]
        i=0
        while 10*i+10<len(seq1):
            pseq.append(seq1[10*i:10*i+10])
            pseqrr.append(seq1rr[10*i:10*i+10])
            i+=1
        else:
            pseq.append(seq1[10*i:len(seq1)])
            pseqrr.append(seq1rr[10*i:len(seq1)])
        l=cb2.get()
        l=DNAform.index(l)
        if l==0:
            j=0
            while j<len(pseq)//5:
                p.insert(END, pseq[5*j:5*j+5])
                p.insert(END, '\n')
                j+=1
            p.insert(END, pseq[5*j:len(pseq)])
            p.insert(END, '\n')
            p.insert(END, '\n')
            j=0
            while j<len(pseq)//5:
                p.insert(END, pseqrr[5*j:5*j+5])
                p.insert(END, '\n')
                j+=1
            p.insert(END, pseqrr[5*j:len(pseqrr)])
            p.config(state=DISABLED)

        if l==1:
            j=0
            while j<len(pseq)//5:
                p.insert(END, pseq[5*j:5*j+5])
                p.insert(END, '\n')
                p.insert(END, pseqrr[5*j:5*j+5])
                p.insert(END, '\n')
                p.insert(END, '\n')
                j+=1
            p.insert(END, pseq[5*j:len(pseq)])
            p.insert(END, '\n')
            p.insert(END, pseqrr[5*j:len(pseqrr)])      
        p.config(state=DISABLED)

    if w==4:
        rseq=[]
        for z in range(len(seq1)-1,-1,-1):
            rseq.append(seq1[z])
        rseq1=''.join(rseq)
        
        pseq=[]
        i=0
        while 10*i+10<len(rseq1):
            pseq.append(rseq1[10*i:10*i+10])
            i+=1
        else: pseq.append(rseq1[10*i:len(rseq1)])
        i=0
        while 5*i+5<len(pseq):
            p.insert(END, pseq[5*i:5*i+5])
            p.insert(END,'\n')
            i+=1
        else: p.insert(END, pseq[5*i:len(pseq)])
        p.config(state=DISABLED)
            
        
    if w==5:
        f1=[]
        f2=[]
        f3=[]

        for i in range(len(seq1)//3):
            f1.append(seq1[3*i:3*i+3])
            f2.append(seq1[3*i+1:3*i+4])
            f3.append(seq1[3*i+2:3*i+5])
            
        a3a1=[]
        a3a2=[]
        a3a3=[]
        a3a1r=[]
        a3a2r=[]
        a3a3r=[]

        for i in f1:
            if i in codon3.keys():
                a3a1.append(codon3[i])
            elif 'N' in i:
                a3a1.append('Xaa')
            else: continue
        for i in f2:
            if i in codon3.keys():
                a3a2.append(codon3[i])
            elif 'N' in i:
                a3a2.append('Xaa')
            else: continue
        for i in f3:
            if i in codon3.keys():
                a3a3.append(codon3[i])
            elif 'N' in i:
                a3a3.append('Xaa')
            else: continue

        f1=[]
        f2=[]
        f3=[]

        for i in range(len(seq1)//3):
            f1.append(seq1r[3*i:3*i+3])
            f2.append(seq1r[3*i+1:3*i+4])
            f3.append(seq1r[3*i+2:3*i+5])

        for i in f1:
            if i in codon3.keys():
                a3a1r.append(codon3[i])
            elif 'N' in i:
                a3a1r.append('Xaa')
            else: continue
        for i in f2:
            if i in codon3.keys():
                a3a2r.append(codon3[i])
            elif 'N' in i:
                a3a2r.append('Xaa')
            else: continue
        for i in f3:
            if i in codon3.keys():
                a3a3r.append(codon3[i])
            elif 'N' in i:
                a3a3r.append('Xaa')
            else: continue

        frame1=''.join(a3a1)
        frame2=''.join(a3a2)
        frame3=''.join(a3a3)
        frame1r=''.join(a3a1r)
        frame2r=''.join(a3a2r)
        frame3r=''.join(a3a3r)
        trfframe=[frame1,frame2,frame3]
        trrframe=[frame1r,frame2r,frame3r]
        
        l=cb1.get()
        l=letter.index(l)
        if l==1:
            p.insert(END, 'Forward\n')
            p.insert(END, '\n')
            p.insert(END, 'frame1\n')
            p.insert(END, frame1+ '\n')
            p.insert(END, '\n')
            p.insert(END, 'frame2\n')
            p.insert(END, frame2+ '\n')
            p.insert(END, '\n')
            p.insert(END, 'frame3\n')
            p.insert(END, frame3+ '\n')
            p.insert(END, '\n\n')
            p.insert(END, 'Reverse\n')
            p.insert(END, '\n')
            p.insert(END, 'frame1r\n')
            p.insert(END, frame1r+ '\n')
            p.insert(END, '\n')
            p.insert(END, 'frame2r\n')
            p.insert(END, frame2r+ '\n')
            p.insert(END, '\n')
            p.insert(END, 'frame3r\n')
            p.insert(END, frame3r)
            p.config(state=DISABLED)
                     
        if l==0:
            p.insert(END, 'Forward\n')
            p.insert(END, '\n')
            for b in trfframe:
                p.insert(END, 'frame'+ str(trfframe.index(b)+1)+'\n')
                tseq=[]
                i=0
                while 3*i+3<len(b):
                    tseq.append(b[3*i:3*i+3])
                    i+=1
                else: tseq.append(b[3*i:len(seq1)])

                tseq1=[]
                for z in tseq:
                    if z in threeto1:
                        tseq1.append(threeto1[z])
                tseq1=''.join(tseq1)

                i=0
                tseq=[]
                while 10*i+10<len(tseq1):
                    tseq.append(tseq1[10*i:10*i+10])                
                    i+=1
                else:
                    tseq.append(tseq1[10*i:len(seq1)])                

                i=0
                while 5*i+5<len(tseq):
                    p.insert(END, tseq[5*i:5*i+5])
                    p.insert(END,'\n')
                    i+=1
                else: p.insert(END, tseq[5*i:len(tseq)])
                p.insert(END, '\n\n')

            p.insert(END, 'Reverse\n')
            p.insert(END, '\n')
            for b in trrframe:
                p.insert(END, 'frame'+ str(trrframe.index(b)+1)+'r'+'\n')
                tseq=[]
                i=0
                while 3*i+3<len(b):
                    tseq.append(b[3*i:3*i+3])
                    i+=1
                else: tseq.append(b[3*i:len(seq1)])

                tseq1=[]
                for z in tseq:
                    if z in threeto1:
                        tseq1.append(threeto1[z])
                tseq1=''.join(tseq1)

                i=0
                tseq=[]
                while 10*i+10<len(tseq1):
                    tseq.append(tseq1[10*i:10*i+10])                
                    i+=1
                else:
                    tseq.append(tseq1[10*i:len(seq1)])                

                i=0
                while 5*i+5<len(tseq):
                    p.insert(END, tseq[5*i:5*i+5])
                    p.insert(END,'\n')
                    i+=1
                else: p.insert(END, tseq[5*i:len(tseq)])
                p.insert(END, '\n\n')

            p.config(state=DISABLED)

    if w==6:
        x=seq1
        y=seq2
        import random
        m=float(match.get())
        mi=float(mismatch.get())
        g=float(gap.get())

        a=[]
        for i in range(len(x)+1):
            a.append([])

        l=0
        for i in range(len(a)):
            a[i].append(l)
            l=l+g
        l=g
        for i in range(len(y)):
            a[0].append(l)
            l=l+g

        for l in range(len(y)):
            for i in range(1,len(x)+1):
                a[i].append(0)

        k=[]
        for i in range(len(x)+1):
            k.append([])

        for i in range(len(a)):
            k[i].append(3)

        for i in range(len(y)):
            k[0].append(2)

        for l in range(len(y)):
            for i in range(1,len(x)+1):
                k[i].append(0)
        k[0][0]=0

        i=0
        j=0
        for _ in range(len(x)):
            i+=1
            j=0
            for _ in range(len(y)):
                j+=1
                if x[i-1]==y[j-1]:
                    s=m
                else:s=mi
                b=a[i-1][j-1]+s
                c=a[i][j-1]+g
                d=a[i-1][j]+g
                a[i][j]=max(b,c,d)

                if max(b,c,d)==b:
                    if b!=c and b!=d:
                        k[i][j]=1
                    elif b==c and b!=d:
                        k[i][j]=12
                    elif b==d and b!=c:
                        k[i][j]=13
                    elif b==c==d:
                        k[i][j]=4
                elif max(b,c,d)==c:
                    if c!=b and c!=d:
                        k[i][j]=2
                    elif c==b and c!=d:
                        k[i][j]=21
                    elif c==d and c!=b:
                        k[i][j]=23
                    elif c==b==d:
                        k[i][j]=4
                elif max(b,c,d)==d:
                    if d!=b and d!=c:
                        k[i][j]=3
                    elif d==b and d!=c:
                        k[i][j]=31
                    elif d==c and d!=b:
                        k[i][j]=32
                    elif d==b==c:
                        k[i][j]=4

        seqxxx=[]
        seqyyy=[]
        for _ in range(10000):
            i=len(x)
            j=len(y)
            
            seqx=[]
            seqy=[]
            
            while i!=0 and j!=0:
                if k[i][j]==1:
                    seqx.append(x[i-1])
                    seqy.append(y[j-1])
                    i+=-1
                    j+=-1
                elif k[i][j]==2:
                    seqx.append('-')
                    seqy.append(y[j-1])
                    j+=-1
                elif k[i][j]==3:
                    seqy.append('-')
                    seqx.append(x[i-1])
                    i+=-1
                elif k[i][j]==4:
                    num=random.randint(1,3)
                    if num==1:
                        seqx.append(x[i-1])
                        seqy.append(y[j-1])
                        i+=-1
                        j+=-1
                    elif num==2:
                        seqx.append('-')
                        seqy.append(y[j-1])
                        j+=-1
                    elif num==3:
                        seqy.append('-')
                        seqx.append(x[i-1])
                        i+=-1
                elif k[i][j]==12:
                    
                    num=random.randint(1,2)
                    if num==1:
                        seqx.append(x[i-1])
                        seqy.append(y[j-1])
                        i+=-1
                        j+=-1
                    elif num==2:
                        seqx.append('-')
                        seqy.append(y[j-1])
                        j+=-1
                elif k[i][j]==13:
                    num=random.randrange(1,4,2)
                    if num==1:
                        seqx.append(x[i-1])
                        seqy.append(y[j-1])
                        i+=-1
                        j+=-1
                    elif num==3:
                        seqy.append('-')
                        seqx.append(x[i-1])
                        i+=-1
                elif k[i][j]==23:
                    num=random.randint(2,3)
                    if num==2:
                        seqx.append('-')
                        seqy.append(y[j-1])
                        j+=-1
                    elif num==3:
                        seqy.append('-')
                        seqx.append(x[i-1])
                        i+=-1
            if i!=0:
                while i!=0:
                    seqx.append(x[i-1])
                    seqy.append('-')
                    i+=-1
            if j!=0:
                while j!=0:
                    seqy.append(y[0])
                    seqx.append('-')
                    j+=-1
            seqxx=[]
            seqyy=[]
            for i in range(len(seqx)-1,-1,-1):
                seqxx.append(seqx[i])
            for i in range(len(seqy)-1,-1,-1):
                seqyy.append(seqy[i])

            seqxx=''.join(seqxx)
            seqyy=''.join(seqyy)
            seqxxx.append(seqxx)
            seqyyy.append(seqyy)



        SEQ=[]
        for i,j in zip(seqxxx,seqyyy):
            SEQ.append((i,j))

        SEQ=set(SEQ)
        SEQ=list(SEQ)

        if len(SEQ)>5:
            for i in range(5):
                p.insert(END,SEQ[i][0])
                p.insert(END, '\n')
                p.insert(END,SEQ[i][1])
                p.insert(END,'\n')
                p.insert(END, '\n')
            p.insert(END,'There are another '+str(len(SEQ)-5)+' sequences')
        else:
            for i in range(len(SEQ)):
                p.insert(END,SEQ[i][0])
                p.insert(END, '\n')
                p.insert(END,SEQ[i][1])
                p.insert(END,'\n')
                p.insert(END, '\n')
        p.config(state=DISABLED)

        
    if w==7:
        aaseqMW=0
        for i in seq1:
            aaseqMW=aaseqMW+aaMW[i]
        aaseqMW=aaseqMW-(18*(len(seq1)-1))
        p.insert(END, 'The length of a.a. seq\n')
        p.insert(END, '\n')
        p.insert(END, len(seq1))
        p.insert(END, '\n\n\n\n')
        p.insert(END, 'Molecular Weight\n')
        p.insert(END, '\n')
        p.insert(END, str(aaseqMW) +'   Da')
        p.config(state=DISABLED)

    if w==8:
        l=cb3.get()
        l=findingmethod.index(l)
        if l==0:
            if seq2 in seq1:
                nos=seq1.count(seq2)
                p.insert(END, int(seq1.index(seq2))+1)
                p.insert(END, '\n')
                p.insert(END, '~')
                p.insert(END, '\n')
                p.insert(END, int(seq1.index(seq2))+len(seq2))
                p.insert(END, '\n')
                p.insert(END, '\n')
                underb=''
                for i in range(len(seq2)):
                    underb=underb+'_'
                temseq1=seq1.replace(seq2,underb,1)
                if nos>=2:
                    for _ in range(nos-1):
                        p.insert(END, int(temseq1.index(seq2))+1)
                        p.insert(END, '\n')
                        p.insert(END, '~')
                        p.insert(END, '\n')
                        p.insert(END, int(temseq1.index(seq2))+len(seq2))
                        p.insert(END, '\n')
                        p.insert(END, '\n')
                        temseq1=temseq1.replace(seq2,underb,1)
                p.config(state=DISABLED)
            else:
                p.insert(END, 'Cannot find')
                p.config(state=DISABLED)
                       
        if l==1:
            esg=es.get()
            etg=et.get()
            fseq=seq1[int(esg)-1:int(etg)]
            pseq=[]
            i=0
            while 10*i+10<len(fseq):
                pseq.append(fseq[10*i:10*i+10])
                i+=1
            else:
                pseq.append(fseq[10*i:len(seq1)])
            j=0
            while j<len(pseq)//5:
                p.insert(END, pseq[5*j:5*j+5])
                p.insert(END, '\n')
                j+=1
            p.insert(END, pseq[5*j:len(pseq)])
            p.insert(END, '\n')
            p.config(state=DISABLED)

    if w==9:
        l=cb4.get()
        l=convert.index(l)
        if l==0:
            errora=[]
            seq1=seq1.upper()

            for i in seq1:
                if i in oneto3:continue
                else:
                    errora.append(i)
            if len(errora)>0:
                msgbox.showerror('Error','There is wrong a.a. name in your sequence!')
            
            tseq1=[]
            for z in seq1:
                tseq1.append(oneto3[z])
            
            i=0
            tseq=[]
            while 10*i+10<len(tseq1):
                tseq.append(tseq1[10*i:10*i+10])                
                i+=1
            else:
                tseq.append(tseq1[10*i:len(seq1)])                
            #if len(tseq1)%10==1:
                #tseq.pop()

            for z in tseq:
                p.insert(END, z)
                p.insert(END,'\n')
            
        if l==1:
            errora=[]
            tseq=[]
            i=0
            while 3*i+3<len(seq1):
                tseq.append(seq1[3*i:3*i+3])
                i+=1               
            else: tseq.append(seq1[3*i:len(seq1)])

            tseq1=[]
            for z in tseq:
                if z in threeto1:
                    tseq1.append(threeto1[z])
                else:
                    errora.append(z)
            if len(errora)>0:
                msgbox.showerror('Error','There is wrong a.a. name in your sequence!')
            tseq1=''.join(tseq1)

            i=0
            tseq=[]
            while 10*i+10<len(tseq1):
                tseq.append(tseq1[10*i:10*i+10])                
                i+=1
            else:
                tseq.append(tseq1[10*i:len(seq1)])                
            #if len(tseq1)%10==1:
                #tseq.pop()

            i=0
            while 5*i+5<len(tseq):
                p.insert(END, tseq[5*i:5*i+5])
                p.insert(END,'\n')
                i+=1
            else: p.insert(END, tseq[5*i:len(tseq)])
            p.config(state=DISABLED)
            

        
            
b1=Button(root, padx=8, pady=9, text='Clear', command=clear)
b1.grid(row=3, column=3, padx=5)

b2=Button(root, padx=8, pady=9, text='run', command=a)
b2.grid(row=3, column=4, padx=5)

sb2=Scrollbar(sframe)
sb2.pack(side='right',fill='y')
tt=Text(sframe, width=54, height=10,yscrollcommand=sb2.set)
tt.insert(END, 'Please enter the seq you want to find\n')
tt.insert(END, 'Or \n')
tt.insert(END,'Query seq for seq alignment')
tt.pack(side='left')
sb2.config(command=tt.yview)



def b():
    global w
    w=cb.get()
    w=work.index(w)
    if w==0 or w==4 or w==7:
        cb1.grid_forget()
        cb2.grid_forget()
        cb3.grid_forget()
        cb4.grid_forget()
        aframe.grid_forget()
        eframe.grid_forget()
        lframe.grid_forget()
        caltcb.grid_forget()
        calmcb.grid_forget()
    if w==1 or w==8:
        cb1.grid_forget()
        cb2.grid_forget()
        cb4.grid_forget()
        aframe.grid_forget()
        eframe.grid(row=2,column=0)
        cb3.grid(row=1,column=0)
        lframe.grid_forget()
        caltcb.grid_forget()
        calmcb.grid_forget()
    if w==3:
        cb1.grid_forget()
        cb3.grid_forget()
        cb4.grid_forget()
        eframe.grid_forget()
        aframe.grid_forget()
        cb2.grid(row=1,column=0)
        lframe.grid_forget()
        caltcb.grid_forget()
        calmcb.grid_forget()
    if w==5:
        cb2.grid_forget()
        cb3.grid_forget()
        cb4.grid_forget()
        eframe.grid_forget()
        aframe.grid_forget()
        cb1.grid(row=1, column=0)
        lframe.grid_forget()
        caltcb.grid_forget()
        calmcb.grid_forget()
    if w==6:
        cb1.grid_forget()
        cb2.grid_forget()
        cb3.grid_forget()
        cb4.grid_forget()
        eframe.grid_forget()
        aframe.grid(row=2,column=0)
        lframe.grid_forget()
        caltcb.grid_forget()
        calmcb.grid_forget()
    if w==9:
        cb1.grid_forget()
        cb2.grid_forget()
        cb3.grid_forget()
        eframe.grid_forget()
        aframe.grid_forget()
        cb4.grid(row=1, column=0)
        lframe.grid_forget()
        caltcb.grid_forget()
        calmcb.grid_forget()
    if w==2:
        cb1.grid_forget()
        cb2.grid_forget()
        cb3.grid_forget()
        cb4.grid_forget()
        aframe.grid_forget()
        eframe.grid_forget()
        caltcb.grid(row=1,column=0)
        calmcb.grid(row=2,column=0)
        lframe.grid(row=2,column=0)
        
    
b3=Button(wframe, padx=8, pady=5, text='choose', command=b)
b3.grid(row=0, column=1)

psb=Scrollbar(tframe)
psb.pack(side='right',fill='y')
p=Text(tframe, width=54, height=33,yscrollcommand=psb.set)
p.insert(END,'How to use this program\n')
p.insert(END,'\n')
p.insert(END,'Choose what you want to do\n')
p.insert(END,'->\n')
p.insert(END,'Press choose button\n')
p.insert(END,'->\n')
p.insert(END,'Choose option\n')
p.insert(END,'->\n')
p.insert(END,'Press clear button\n')
p.insert(END,'->\n')
p.insert(END,'Enter all required things\n')
p.insert(END,'->\n')
p.insert(END,'Press run button\n')
p.config(state=DISABLED)
p.pack(side='left')
psb.config(command=p.yview)

root.mainloop()
