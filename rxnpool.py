# !/usr/bin/python
#  _*_ coding:utf-8 _*_
# 2016/08/04
# Author:LingWu
# Email:wu_l@tib.cas.cn  

import os 
import json
import itertools
import copy
import time
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import DataStructs
from rdkit.Chem import MolFromSmarts as mfsma 
from rdkit.Chem import MolFromSmiles as mfsmi 
from rdkit.Chem import MolToSmarts as mtsma
from rdkit.Chem import MolToSmiles as mtsmi
from rdkit.Chem.AllChem import ReactionFromSmarts as rxnfsma
from rdkit.Chem.Fingerprints import FingerprintMols


def BioReactor(queryRsSmi,queryPsSmi=False,Ec=False,CALSIMI=True,Draw=True):
    
    Rid_EcAndPattern = json.load(open("./Rid_EcAndPattern.json")) 
    rhea_ECAssigner = json.load(open("./rhea_ECAssigner.json"))
    ridPatRsPs = RidPatRsPs(Ec)

    prePsSmi_preRsSmiRxnid= defaultdict() 
    ridSubSimilarity = dict()
    
    queryRsMol= mfsmi(queryRsSmi) 
    (n,m,k) = (0,0,0) 

    for rid,Pat in ridPatRsPs.items():
        n += 1
        
        SamepatId = Rid_EcAndPattern[rid]["SamepatId"].keys()
        SamepatId.append(rid)

        hasSmirks = False  
        for _rid in SamepatId:
            try:
                smirks = rhea_ECAssigner[rid]["smirks"]
                hasSmirks = True
            except:
                continue
        if not hasSmirks:
            continue

        patRsSmi = [str(pat)for pat in Pat["patrs"]] 
        Num = int(len(patRsSmi)-1) 

        queryPosition = list()
        for a,pa_rs_smi in enumerate(patRsSmi):   
            has = pa_rs_smi.startswith("&")
            if has:
                pa_rs_smi = pa_rs_smi.split("(",1)[1].rsplit(")",1)[0].split(".")
                pa_rs_mol = [mfsma(smi) for smi in pa_rs_smi]
                
            else:
                pa_rs_mol= mfsma(pa_rs_smi)
            
            if (not has and queryRsMol.HasSubstructMatch(pa_rs_mol)) or \
                ( has and queryRsMol.HasSubstructMatch(pa_rs_mol[0]) and queryRsMol.HasSubstructMatch(pa_rs_mol[1])):
                m += 1
                queryPosition.append(a) 

        if queryPosition: 
            
            preRsSmi_patID_Dic = dict() 

            for patID in SamepatId: 
                patID = str(patID)
                try:
                    origRxnSmirks = rhea_ECAssigner[patID]["smirks"]  
                except:
                    continue

                _origRxnRsSmi =[str(smi) for smi in  origRxnSmirks.split(">>")[0].split(".")]
                origRxnRsSmi = copy.deepcopy(_origRxnRsSmi)
                
                for b,b_origrs in enumerate(_origRxnRsSmi):
                    if b_origrs == "[H+]": 
                        origRxnRsSmi.remove(b_origrs)
                    elif b_origrs == "[H]O[H]":
                        origRxnRsSmi[b] = "[OH2]"
                
                if len(origRxnRsSmi) != len(patRsSmi): 
                    continue
                
                origRsPosition = dict() 
                
                repeat = 0
                for c,c_origrs in enumerate(origRxnRsSmi):

                    if c_origrs not in origRsPosition.keys():
                        origRsPosition[c_origrs] = list()
                    else:
                        repeat += 1
                        c_origrs = c_origrs+'$%s' % repeat
                        origRsPosition[c_origrs] = list()

                    if c_origrs.count('$') == 0:
                        for d,d_sub in enumerate(patRsSmi):
                            has = d_sub.startswith('&')
                            if not has :
                                hasSub = mfsmi(c_origrs).HasSubstructMatch(mfsma(d_sub))
                                if hasSub:
                                    origRsPosition[c_origrs].append(d)
                            else:
                                d_sub_listsmi = d_sub.split('(',1)[1].rsplit(')',1)[0].split('.')
                                hasSub = True
                                for dsub in d_sub_listsmi:
                                    if not mfsmi(c_origrs).HasSubstructMatch(mfsma(dsub)):
                                        hasSub = False
                                if hasSub:
                                    origRsPosition[c_origrs].append(d)
                    else:
                        origRsPosition[c_origrs] = origRsPosition[c_origrs.split('$')[0]]
                
                for po in queryPosition: 
                    preRsSmi = range(len(patRsSmi)) 
                    preRsSmi[po] = queryRsSmi 

                    for origrs,polist in origRsPosition.items():
                        _polist = copy.deepcopy(polist)
                        for _po in polist:
                            if _po == po :
                                _polist.remove(po)
                        
                        if not _polist:
                            del(origRsPosition[origrs])

                            if CALSIMI : 
                                if patID not in ridSubSimilarity.keys():
                                    ridSubSimilarity[patID] = dict() 

                                if origrs.count('$') != 0:
                                    origrs = origrs.split('$')[0]
                                ridSubSimilarity[str(patID)][origrs] = similar(origrs,queryRsSmi)
                        else:
                            origRsPosition[origrs] = _polist

                    if len(origRsPosition.keys()) == Num:
                        
                        remainOrigRsPoCombin = [x for x in list(itertools.product(*origRsPosition.values())) if len(set(x)) == Num] #1

                        for remainOrignRsPo in remainOrigRsPoCombin:

                            preRsSmiNew = copy.deepcopy(preRsSmi)

                            for e,e_po in enumerate(remainOrignRsPo):

                                e_key_origRsPosition = origRsPosition.keys()[e] 

                                if e_key_origRsPosition.count('$') != 0:
                                    _e_key_origRsPosition = e_key_origRsPosition.split('$')[0]
                                    preRsSmiNew[e_po] = _e_key_origRsPosition
                                else:
                                    preRsSmiNew[e_po] = e_key_origRsPosition

                            preRsSmiNew = ".".join(preRsSmiNew)
                            if preRsSmiNew not in preRsSmi_patID_Dic.keys():
                                preRsSmi_patID_Dic[preRsSmiNew] = list()
                                preRsSmi_patID_Dic[preRsSmiNew].append(patID)
                        
                            else:
                                preRsSmi_patID_Dic[preRsSmiNew].append(patID)

                    else:
                
                        remainOrigRsCombin =list(itertools.combinations(origRsPosition.keys(),Num))
                        for remainOrigRs in remainOrigRsCombin:

                            remainOrigRsPoSet = [origRsPosition[x] for x in remainOrigRs]
                            if len(reduce(lambda x,y:set(x)|set(y),remainOrigRsPoSet)) != Num:
                                continue


                            _origRsPosition = copy.deepcopy(origRsPosition)
                            [_origRsPosition.pop(z) for z in origRsPosition.keys() if z not in remainOrigRs]

                            left = list(set(origRsPosition.keys())-set(remainOrigRs))[0]
                            
                            if CALSIMI :
                                if patID not in ridSubSimilarity.keys():
                                    ridSubSimilarity[patID] = dict() 

                                if left.count("$") != 0:
                                    left = left.split('$')[0]
                                
                                ridSubSimilarity[str(patID)][left] = similar(left,queryRsSmi)

                            remainOrigRsPoCombin = [q for q in list(itertools.product(*_origRsPosition.values())) if len(set(q)) == Num]
                            
                            for remainOrignRsPo in remainOrigRsPoCombin:

                                preRsSmiNew = copy.deepcopy(preRsSmi)

                                for f,f_po in enumerate(remainOrignRsPo):
                                    f_key_origRsPosition = _origRsPosition.keys()[f]
                                    
                                    if f_key_origRsPosition.count('$') != 0:
                                        _f_key_origRsPosition = f_key_origRsPosition.split('$')[0]
                                        preRsSmiNew[f_po] = _f_key_origRsPosition

                                    else:
                                        preRsSmiNew[f_po] = f_key_origRsPosition
                                
                                preRsSmiNew = ".".join(preRsSmiNew)
                                if preRsSmiNew not in preRsSmi_patID_Dic.keys():
                                    preRsSmi_patID_Dic[preRsSmiNew] = list()
                                    preRsSmi_patID_Dic[preRsSmiNew].append(patID)
                                else:
                                    preRsSmi_patID_Dic[preRsSmiNew].append(patID)

            reaction = rxnfsma(str(Rid_EcAndPattern[rid]["pat"]))
            
            preRsMol_list = [[mfsmi(y) for y in z.split(".")] for z in preRsSmi_patID_Dic.keys()]

            for g,preRsMol in enumerate(preRsMol_list):  
                prePsMol_tuple = reaction.RunReactants(preRsMol) 
                
                for prePsMol in prePsMol_tuple: 

                    if bool(prePsMol): 
                        if queryPsSmi:
                            pshassub = False
                            pssub = mfsma(queryPsSmi)
                            for psmol in prePsMol:
                                if psmol.HasSubstructMatch(pssub):
                                    pshassub = True
                        else:  
                            pshassub = True
                        if pshassub:
                            prePsSmi = sorted([NeutraliseCharges(mtsmi(h))[0] for h in list(prePsMol)])#得到smile格式的其中一种产物组合，list形式,并对产物进行去电荷处理，方便后期去重

                            prePsSmi = ".".join(prePsSmi)
                            if prePsSmi in prePsSmi_preRsSmiRxnid:
                                preRsSmiRxnid = prePsSmi_preRsSmiRxnid[prePsSmi]
                                _preRsSmi = ".".join(sorted(preRsSmi_patID_Dic.keys()[g].split(".")))
                                if _preRsSmi in preRsSmiRxnid.keys():
                                    preRsSmiRxnid[_preRsSmi].extend(preRsSmi_patID_Dic[preRsSmi_patID_Dic.keys()[g]])
                                    preRsSmiRxnid[_preRsSmi] = list(set(preRsSmiRxnid[_preRsSmi]))
                                else:
                                    preRsSmiRxnid[_preRsSmi] = preRsSmi_patID_Dic[preRsSmi_patID_Dic.keys()[g]]
                            else:
                                preRsSmiRxnid = dict()
                                _preRsSmi = ".".join(sorted(preRsSmi_patID_Dic.keys()[g].split(".")))
                                preRsSmiRxnid[_preRsSmi] = preRsSmi_patID_Dic[preRsSmi_patID_Dic.keys()[g]]
                                prePsSmi_preRsSmiRxnid[prePsSmi] = preRsSmiRxnid

                                #{'ps':{'rs':[id]}}
                        else:
                            continue
    if bool(prePsSmi_preRsSmiRxnid):

        currenttime = time.strftime(r'%Y-%m-%d-%H-%M-%S',time.localtime(time.time()))
    
        resultDir = "./preResult/%s/" % currenttime
        
        if not os.path.exists(resultDir):
            os.makedirs(resultDir)
        
        it = 0
        preResultsDic = defaultdict()

        f=open(os.path.join(resultDir+"result.txt"),"w")
        f.write('query : '+queryRsSmi+'\n')
        
        for ps,rsinfo in prePsSmi_preRsSmiRxnid.items():
            for rs,idlist in rsinfo.items():
                smirks = rs+">>"+ps
                it += 1

                line1 = "predicted reaction smirks : %s" % it
                f.write(line1+"\n"+"Smirks:"+smirks+"\n")

                with open(os.path.join(resultDir+str(it))+".smi","w") as smif:
                    smif.write(smirks)
                if Draw:
                    rxnSmiFile = os.path.join(resultDir+str(it)+".smi")
                    rxnImagefile = os.path.join(resultDir+str(it)+".png")
                    DrawImage(rxnSmiFile,rxnImagefile)

                preResultsDic[it] = dict()

                preResultsDic[it]["Smikrs"] = smirks
                
                ref = 0
                for _id in idlist:
                    ref+=1

                    preResultsDic[it]["Ref Rxn|Ec|Similarity %s" % ref] = str([_id,rhea_ECAssigner[_id]["ecnumber"],max(ridSubSimilarity[_id].values())])
                    line2 = "Ref%s Rxnid:%s    Ec:%s    Similarity:%s " % (ref,_id,str(rhea_ECAssigner[_id]["ecnumber"]),max(ridSubSimilarity[_id].values()))
                    f.write(line2+"\n")
                    
                    if CALSIMI:
                        preResultsDic[it].update(ridSubSimilarity[_id])
                        for sub,val in ridSubSimilarity[_id].items():
                            f.write("--Substrate:"+sub+"\n"+"--Similarity:"+str(val)+"\n")
                f.write("-"*100+"\n")          
        f.close()
        
        rs_smilarity = dict()
        
        for _id,rssimi in ridSubSimilarity.items():
            rs_smilarity[_id] = str(max(rssimi.values()))
        
        with open(os.path.join(resultDir,'result.json'),"w") as f1:
            json.dump(preResultsDic,f1,indent = 2)
        with open(os.path.join(resultDir,'similarity.txt'),"w") as f2:
            f2.writelines(":".join(list(x))+"\n" for x in arrange(rs_smilarity))

        print "总迭代次数   = ",n
        
    else:
        print "No ReSult"

def splitPat(rid,item,ridPatRsPs):

    ridPatRsPs[rid] = dict()

    _ridPatRs = item["pat"].split(">>")[0].split(".")
    _ridPatPs = item["pat"].split(">>")[1].split(".")
    
    ridPatRs = list()
    ridPatPs = list()
    
    _PatDic = {"rs":_ridPatRs,"ps":_ridPatPs}
    PatDic = {"rs":ridPatRs,"ps":ridPatPs}
    
    for t,pat in _PatDic.items():
        for i,s in enumerate(pat):
            num1 = s.count("(")
            num2 = s.count(")")
            if num1 == num2:
                PatDic[t].append(s)
            elif num1 > num2:
                PatDic[t].append("&"+pat[i]+"."+pat[i+1]) 
            else:                                       
                pass
    ridPatRsPs[rid]["patrs"] = ridPatRs
    ridPatRsPs[rid]["patps"] = ridPatPs    

def RidPatRsPs(Ec = False):
    
    ridPatRsPs = defaultdict() 

    Rid_EcAndPattern = json.load(open("./Rid_EcAndPattern.json")) 
    
    for rid,item in Rid_EcAndPattern.items(): 
        
        if Ec:
            contain = False
            ecList = item["Ec"]
            for ec in ecList:
                if str(ec).startswith(str(Ec)):
                    contain = True
            
            if contain:
                splitPat(rid,item,ridPatRsPs)
            else:
                continue
        else:
            splitPat(rid,item,ridPatRsPs)

    with open("./ridPatRsPs.json","w") as f:
        json.dump(ridPatRsPs,f,indent = 2)
    
    return ridPatRsPs


def InitialiseNeutralisationReactions(): 

    patts= (
    # Imidazoles
    ('[n+;H]','n'),
    # Amines
    ('[N+;!H0]','N'),
    # Carboxylic acids and alcohols ('[$([O-]);!$([O-][#7])]','O'), # Thiols
    ('[S-;X1]','S'),
    # Sulfonamides 
    ('[$([N-;X2]S(=O)=O)]','N'), 
    # Enamines 
    ('[$([N-;X2][C,N]=C)]','N'), 
    # Tetrazoles 
    ('[n-]','[nH]'),
    # Sulfoxides 
    ('[$([S-]=O)]','S'),
    # Amides
    ('[$([N-]C=O)]','N'),
    #
    ('[O-;X1]',"O"),
    #
    ('[$([O-]=C)]','O')
    )
    return [(Chem.MolFromSmarts(x),Chem.MolFromSmiles(y,False)) for x,y in patts]

_reactions=None

def NeutraliseCharges(smiles, reactions=None):
    global _reactions
    if reactions is None: 
        if _reactions is None:
            _reactions=InitialiseNeutralisationReactions()
        reactions=_reactions 
    mol = Chem.MolFromSmiles(smiles)

    if not mol:
        mol = Chem.MolFromSmarts(smiles)

    replaced = False
    for i,(reactant, product) in enumerate(reactions):
        while mol.HasSubstructMatch(reactant): 
            replaced = True
            rms = AllChem.ReplaceSubstructs(mol, reactant, product) 
            mol = rms[0] 
    if replaced:
        return (Chem.MolToSmiles(mol,True), True) 
    else:
        return (smiles, False)

def similar(smiles1,smiles2):
    smiles1_mol = mfsmi(smiles1)
    if not bool(smiles1_mol):
        smiles1_mol = mfsma(smiles1)
    
    smiles2_mol = mfsmi(smiles2)
    if not bool(smiles2_mol):
        smiles2_mol = mfsma(smiles2)
    
    refmolfp = FingerprintMols.FingerprintMol(smiles1_mol)
    molfp = FingerprintMols.FingerprintMol(smiles2_mol)
    similarity = DataStructs.FingerprintSimilarity(refmolfp,molfp)
    return similarity
    
def arrange(dic):
    array = sorted(dic.iteritems(),key=lambda t:t[1],reverse=True)
    return array

def DrawImage(rxnSmiFile,rxnImagefile):
    
    command='/Applications/ChemAxon/MarvinBeans/bin/molconvert png:w1600,h800 '+ rxnSmiFile+ ' -o '+ rxnImagefile
    os.popen(command)
    return rxnImagefile

def testSmi(sm='',rxn =''):
    if sm:
        mol = mfsma(sm)
        if mol:
            print "smi correct"
        else:
            print "smi error"
    if rxn:
        rule = AllChem.ReactionFromSmarts(rxn)
        if rule:
            print "pattern correct"
        else:
            print "pattern error"

def main():
    start = time.strftime(r'%Y-%m-%d-%H-%M-%S',time.localtime(time.time()))
    
    queryRsSmi = 'CSC'
    BioReactor(queryRsSmi,queryPsSmi=False,Ec=False,CALSIMI=True,Draw=False)
    
    end = time.strftime(r'%Y-%m-%d-%H-%M-%S',time.localtime(time.time()))
    print start
    print end
        
if __name__ == "__main__":
    main()





