#!/usr/bin/env python3

from Bio.KEGG.REST import kegg_get
from Bio.KEGG.REST import kegg_list
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas
from collections import defaultdict
import nsdm
import matplotlib.cm as cm
import re
import os
import scipy.special as ss
from scipy import interpolate
import seaborn as sn
import matplotlib.transforms as transforms
import numpy as np 
import math
import sys
import plotly.offline
import plotly.graph_objects as go
from time import sleep
from patsy import dmatrix
import statsmodels.api as sm


def main(option,variants,ppifile,gffile):
    gff = nsdm.fileparse.gff_read(gffile)
    gene_max = len(gff)

    if option == 13:
       stripro = {}
    else:
       labels = [ "neighborhood", "neighborhood_transferred" ,"fusion",
                  "cooccurence", "homology", "coexpression" ,
                  "coexpression_transferred" ,"experiments" ,"experiments_transferred" ,
                  "database" ,"database_transferred","textmining",
                  "textmining_transferred"]
       ebi = ""
       for i in option:
           ebi += labels[i] + "  "
       print('evidence type :',ebi,'\n')
       stripro = stringdict(option,ppifile)
    koumoku = "\tpathway number\tpathway name\tnumber of variant proteins\tThe number of proteins in each pathway\tp-value"
    result = []
    vtestsg = []
    up_res ={}
    up_res2 ={}
    accumulation_keys = set()
    accumulation = dict()

    for n  in range(len(variants)):
        f = open( variants[n],"r")
        varif = f.read()
        variresult = varif.split("\n")
        f.close()
        if n == 0:
           species =  variresult[0][0]+variresult[0][1]+variresult[0][2]
           pathways = pathdict(species)
           print("The number of the analyzed pathways is ", len(pathways), ".")
        n += 1
        variresult = set(variresult)
        print("Generation:%s" % (n))
        print(koumoku)
        res, keys,vtests = kgmlr(variresult, pathways, n, gene_max, stripro, option)
        vtestsg.append(vtests)
        result.append(res)
        [accumulation_keys.add(z) for z in keys]
    resultkeys = set()
    accumulation_keys = list(accumulation_keys)
    for i,v in zip(result,vtestsg):
        for k ,j in zip(i,v):
            if j == True:
               resultkeys.add(k)
    len_resultkeys = len(resultkeys)
    true_pathway = []
    for g_res in result:
        true_path = {}
        for name,p in g_res.items():
            if name in resultkeys:
               true_path[name] = p
        true_pathway.append(true_path)    
    ramudas = np.arange(0.01,0.90,0.01) 
    n = 0
    pf = open("q-value.tsv", "w")
    for g_res in true_pathway:#result:#generation
        n += 1
        up_res = sorted(g_res.items(), key=lambda x:-x[1]) 
        qi_plus1 = 1
        pais = []
        for ramuda in ramudas:#lambda
            big_pi = len([i for i in g_res.values() if i >= ramuda])
            pai = big_pi/(len_resultkeys * ( 1 - ramuda))
            pais.append(pai)
        cubic = interpolate.interp1d(ramudas,np.array(pais), kind="cubic",fill_value='extrapolate') 
        cs = interpolate.CubicSpline(ramudas,np.array(pais),bc_type = "natural")
        x3 = dmatrix("cr(ramu,df = 3)", {"ramu": ramudas}, return_type='dataframe')
        fit3 = sm.GLM(np.array(pais), x3).fit()
        x = np.arange(0,1.01,0.05)
        pred3 = fit3.predict(dmatrix("cr(x, df=3)", {"xp": x}, return_type='dataframe'))
        pai_0 = pred3[20]
        pf.write(f"Generation {n}:\n")
        j = 0
        for name in up_res:##Q-value(storey)
            if name[0] in resultkeys:
               i = len_resultkeys - j
               j += 1
               q = name[1]*len_resultkeys*pai_0/i #story
               if q > qi_plus1:
                  q = qi_plus1
               qi_plus1 = q
               pf.write(
                      f"\t{q}\t{name[0]}\n")
               accumulation.setdefault(name[0], []).append(q)
    print("INFO:finish! The number of pathways containing variant proteins is", len_resultkeys, ".")
    print("\t")
    ydata = dict()
    data = [] 
    trace_line = go.Scatter(
           x=[0, n+1],
           y=[0.05,0.05],
           mode='lines',
           line = dict(color="red"),
           name = '0.05'
           )
    data.append(trace_line)
    layout = dict(#title = dict(text = ''),
                  xaxis = dict(title = 'generation'),
                  yaxis = dict(title = 'q-value'),
              )
    x_generate = range(1,len(accumulation)+1)
    for k,v in accumulation.items():
        trace = go.Scatter(
              y = v,
	      x =list(x_generate),
              name = k,
              line_shape='linear',
              )
        data.append(trace)
    fig = dict(data=data, layout=layout)     
    plotly.offline.plot(fig, filename='q-value.html')
    pf.close()


def getkgml(pid):
    try:
        kgml = KGML_parser.read((kegg_get(pid, "kgml").read()))
        return kgml
    except KeyboardInterrupt:
        exit()
    except:
        sleep(1)
        return getkgml(pid)


def drowkgml(canvas, outfilename):
    try:
        canvas.draw(outfilename)
    except KeyboardInterrupt:
        exit()
    except:
        print("reconnect")
        drowkgml(canvas, outfilename)


def kgmlr(result, pathways, generation, gmax, stripro ,option):
    output = "./pathway_images/"
    os.makedirs(output, exist_ok=True)
    bbhdict = result
    generation = str(generation)
    vnum = len(bbhdict)
    print(vnum)
    ret = dict()
    accumulation = set()
    epathway = 0
    kari = 0 
    yuukou = 0
    vtests = []
    all_SCO = set()
    all_variSCO =set()
    all_kimuraSCO =set()
    all_kimuravariSCO =set()
    outfile = 'pathway_geneID'
    os.makedirs(outfile, exist_ok=True)
    sf = open("./pathway_geneID/geneID_list_"+generation+".txt","w")
    for pid in pathways.keys():
        kari +=1
        test3 = 0
        varigene = []
        rateadd = []
        variantgene_in_pathway = set()
        pmax = set()
        kgml = getkgml(pid)
        element = kgml.genes          
        judge = 0
        VnodeTgene = defaultdict(list) ##mutant
        NnodeTgene = defaultdict(list) ##no-mutant
        kumura_SCO = set()
        kumura_variSCO= set()
        for n, e in enumerate(element):#pathway node
            test2 = 0  
            genenum = 0
            kgml.genes[n].graphics[0].name = generation.zfill(
                4) + ":" + kgml.genes[n].graphics[0].name
            gene_name=set()
            for v in e._names: #gene
                pathway_gene = v.split(":")[1]
                pmax.add(pathway_gene)
                all_SCO.add(pathway_gene)
                gene_name.add(pathway_gene)		
                all_kimuraSCO.add(pathway_gene) 
                kumura_SCO.add(pathway_gene)
                interact = stripro.get(pathway_gene) 
                if interact !=None:
                   for pp in interact :##string data
                       pmax.add(pp)
                       all_SCO.add(pp)   
                       genenum +=1
                       if pp in bbhdict:
                          test2 += 1
                          all_variSCO.add(pp) 
                          variantgene_in_pathway.add(pp)
                genenum +=1
                if pathway_gene in bbhdict:
                    test2 +=1
                    all_variSCO.add(pathway_gene)    
                    kumura_variSCO.add(pathway_gene)   
                    all_kimuravariSCO.add(pathway_gene)
                    variantgene_in_pathway.add(pathway_gene) 
            if test2 == 0:
               kgml.genes[n].graphics[0].bgcolor = "#BCDDE5"
               continue
            test3 +=1
            rate = round(test2/genenum,3)##mutation rate
            varigene.append(test2)
            rateadd.append(rate)
            color = cm.Reds(rate/0.12)#node color
           # color = cm.Reds(math.log(test2+1,38+1))
           # color = cm.Reds(test2/5)
            kgml.genes[n].graphics[0].bgcolor = color
        if test3 == 0:
           vtest = False
        else:
           vtest = True
        canvas = KGMLCanvas(
            kgml,
            import_imagemap=True,
            fontsize=4,
            margins=(0.005, 0.005),
        )  	
        outfilename = output + pid + "_" + generation.zfill(4) + ".pdf"
        drowkgml(canvas, outfilename)
        spmax = str(pmax).replace("'", "").replace("}","").replace("{","").replace("set()","")
        svariantgene = str(variantgene_in_pathway).replace("'", "").replace("}","").replace("{","").replace("set()","")
        sf.write(f"{pid}\t{kgml.title}\n"\
                 f"proteins in each pathway:\n{spmax}\n"\
                 f"variant proteins:\n{svariantgene}\n\n")
        variantgene_in_pathway = len(variantgene_in_pathway)
        pmax = len(pmax)
        accumulation.add(kgml.title)	
        kumura_SCO= len(kumura_SCO)
        kumura_variSCO=len(kumura_variSCO)
        if option==13:
           sv = svalue(vnum, gmax, kumura_SCO, kumura_variSCO)
           if sv < 0.05:
              yuukou += 1
           if vtest == True: 
              print('\t', pid,'\t', kgml.title,'\t',kumura_variSCO , '\t', kumura_SCO,'\t',sv)   
        else:		
           sv = svalue(vnum, gmax, pmax, variantgene_in_pathway)   
           if sv < 0.05:
              yuukou += 1
           if vtest == True:
              print('\t', pid,'\t', kgml.title,'\t',variantgene_in_pathway , '\t', pmax,'\t',sv)   
        ret[kgml.title] = sv
        vtests.append(vtest)

    all_SCO =len(all_SCO)
    all_variSCO =len(all_variSCO) 
    all_kimuraSCO =len(all_kimuraSCO) 
    all_kimuravariSCO = len(all_kimuravariSCO)
    print("The number of proteins is",all_SCO,'\n'
          "The number of variant proteins is",all_variSCO ,'\n'
	  "The number of the analyzed pathways is ", yuukou, "\n")
    sf.close()

    return ret, list(accumulation) ,vtests



def svalue(vnum, gmax, n, variantgene_in_pathway): ##p-value	
    ret = []
    P = vnum / gmax
    co = 0
    comv = 0
    for i in range(variantgene_in_pathway,n+1): ##binomial distribution
        r = i  ##nCr
        if n - r < r: 
           r = n - r
        numerator = [n - r + k + 1 for k in range(r)]
        denominator = [k + 1 for k in range(r)]

        for p in range(2,r+1):
            pivot = denominator[p - 1]
            if pivot > 1:
                offset = (n - r) % p
                for k in range(p-1,r,p):
                    numerator[k - offset] /= pivot
                    denominator[k] /= pivot
        comv = (1-P)**(n-i)
        count =0
        for k in range(r):
            if numerator[k] > 1:
                comv *= numerator[k] 
                if comv > 1.0e+150 :
                   comv = comv * (P**100)
                   i = i - 100
            count += 1
        tmp = comv *(P**i)
        co +=1 
        ret.append(tmp)
    ret = sum(ret)
    return ret



def pathdict(hid):
    pathway = kegg_list("pathway", hid)
    ret = dict()
    repatter = re.compile("path:" + hid + "011.." +
                          "|" + "path:" + hid + "012..")
    for i in pathway:
        i = i.strip().split()
        key = i[0]
        if repatter.match(key):
            continue
        value = " ".join(i[1:])
        ret[key] = value
    return ret



def stringdict(typ,ppifile):#STRING
    w = open(ppifile,"r")
    string = w.read()
    w.close()

    protein1 = []
    protein2 = []
    protein11 = []
    protein22 = []

    evidence = []
    score = string.split('\n')
    score.pop(0)
    score.pop(-1)

    for i in score :
        i = re.split(' ',i)
        if int(i[-1]) < 400:##evidence score
           continue
        protein1 += re.findall('[A-Z0-9]{7}',i[0])
        protein2 += re.findall('[A-Z0-9]{7}',i[1])
        ty = []
        for k in typ :
            t = i[k+2]
            ty.append(t)
        evidence.append(ty) 

    stringdict = {}
    pro2ID = []
    #print(len(protein1))
    for i in range(len(protein1)):
        ev = evidence[i]
        if ev.count('0') != len(ev):
           protein11.append(protein1[i])
           protein22.append(protein2[i])

    for i in range(len(protein11)):
        pro2ID.append(protein22[i])
        if  i+1 == len(protein11) or protein11[i] != protein11[i+1] :
            stringdict[protein11[i]] = pro2ID
            pro2ID = []
            if i+1 == len(protein11) :
               break

    return stringdict


def parser():
    usage = 'Usage: python {} [--help]  [--ppi] [evidence_type]'\
            .format(__file__)
    command ="\n\n\t-h  or --help\tshow this help message"+\
               "\n\t-p  or --ppi\tversion ppi" +\
               "\n\t  evidence_type"+\
               "\n\t\t0.  neighborhood"+\
               "\n\t\t1.  neighborhood_transferred"+\
               "\n\t\t2.  fusion"+\
               "\n\t\t3.  cooccurence"+\
               "\n\t\t4.  homology"+\
               "\n\t\t5.  coexpression"+\
               "\n\t\t6.  coexpression_transferred"+\
               "\n\t\t7.  experiments"+\
               "\n\t\t8.  experiments_transferred"+\
               "\n\t\t9.  database"+\
               "\n\t\t10. database_transferred"+\
               "\n\t\t11. textmining"+\
               "\n\t\t12. textmining_transferred"+\
             "\n\n\t-pを加えるとppiのevidence_typeを指定出来ます。evidence_typeか上記の番号を引数に書いて下さい。"+\
	     "\n\t-pまたは--ppiのみで行うと全てのevidence_typeの情報を付与します。"+\
             "\n\tevidence_typeはSTRINGに準じています。 http://version10.string-db.org/help/getting_started/"+\
                "\n\t　ex) 'python {} -p 1  database'".format(__file__)+\
               "\n\t   　  はneighborhood_transferredとdatabaseのevidence_typeを持っているppiの情報をpathwayに付与します。"+\
              "\n\n\t引数なしで行うと、ppiを付与しないdefault versionになります。"
    
    ppitype = ["0",'neighborhood',"1",'neighborhood_transferred',"2",'fusion',"3",'cooccurence',
		"4",'homology',"5",'coexpression',"6",'coexpression_transferred',"7",'experiments',
		"8",'experiments_transferred',"9",'database',"10",'database_transferred',
		"11",'textmining',"12",'textmining_transferred']

    option = sys.argv
    vfile = [s for s in option if s.endswith('.txt')]
    gffile = [s for s in option if s.endswith('.gff3')]
    # option
    if '-h' in option or '--help' in option:
        print( usage, command)
        exit()
	
    if '-p' in option or '--ppi' in option:
        print("PPI VERSION")
        typ = []
        for i in range(0,len(ppitype),2):
            if ppitype[i] in option or ppitype[i+1]in option:
               typ.append(int(i/2))
        if len(typ) == 0:
            typ = [0,1,2,3,4,5,6,7,8,9,10,11,12] 
        ppifile = vfile[-1]
        vfile.pop(-1)
        main(typ,vfile,ppifile,gffile[0])
    else:
        print("DEFAULT VERSION")
        main(13,vfile,None,gffile[0])  



if __name__ == '__main__':
    parser()
