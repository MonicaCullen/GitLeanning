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
                                c_origrs_mol = mfsmi(c_origrs) ####modified
                                if not bool(c_origrs_mol): ####modified
                                    c_origrs_mol = mfsma(c_origrs) ####modified
                                hasSub = c_origrs_mol.HasSubstructMatch(mfsma(d_sub))

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
                            prePsSmi = sorted([NeutraliseCharges(mtsmi(h))[0] for h in list(prePsMol)])

                            prePsSmi = ".".join(prePsSmi)
                            if prePsSmi in prePsSmi_preRsSmiRxnid:
                                preRsSmiRxnid = prePsSmi_preRsSmiRxnid[prePsSmi]
                                _preRsSmi = ".".join(sorted(preRsSmi_patID_Dic.keys()[g].split(".")))
                                
                                if _preRsSmi in preRsSmiRxnid.keys():
                                    _preRsSmi_idList = preRsSmi_patID_Dic[preRsSmi_patID_Dic.keys()[g]]

                                    for _preRsSmi_id in _preRsSmi_idList:
                                        if _preRsSmi_id not in preRsSmiRxnid[_preRsSmi]:
                                            preRsSmiRxnid[_preRsSmi].append(_preRsSmi_id)
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
        
                    
    #存结果
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

                #文本文件保存
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

def customize(Rs,Pattern,queryPsSmi=False):

    patOK = True
    for pat in Pattern:
        try:
            pat = rxnfsma(pat)
        except Exception,e :
            patOK = False
            return ('Input Rules error',pat,Exception,':',e)
            break

    if patOK:
        preRsSmiList = list(itertools.product(*Rs))
        preRsMolList = [[mfsmi(i) for i in j] for j in preRsSmiList]
        
        prePsSmi_preRsSmi_pat = dict()

        for pat in Pattern:
            
            patRs = pat.split('>>')[0].split('.')
            if len(patRs) == len(Rs):
                reaction = rxnfsma(pat)
                
                for a,preRsMol in enumerate(preRsMolList):
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
                                preRsSmi = ".".join(sorted(preRsSmiList[a]))

                                if prePsSmi in prePsSmi_preRsSmi_pat.keys():
                                    preRsSmi_pat = prePsSmi_preRsSmi_pat[prePsSmi]
                                    
                                    if preRsSmi in preRsSmi_pat.keys():
                                        preRsSmi_pat[preRsSmi].append(pat)

                                    else:
                                        patList = list()
                                        patList.append(pat)
                                        preRsSmi_pat[preRsSmi] = patList

                                else:
                                    prePsSmi_preRsSmi_pat[prePsSmi] = dict()
                                    patList = list()
                                    patList.append(pat)
                                    prePsSmi_preRsSmi_pat[prePsSmi][preRsSmi] = patList

            else:
                return('Input reactants not equal Rule\'s','Rule:',pat)

        with open('./customize.json','w') as fn:
            json.dump(prePsSmi_preRsSmi_pat,fn,indent = 2)

        f = open('./customize.txt','w')
        num = 0
        for ps,rsitems in prePsSmi_preRsSmi_pat.items():
            for rs,patList in rsitems.items():
                smirks = rs + '>>' + ps
                num += 1
                f.write('PredictedSmirks %s:'  % num + '\t' + smirks + '\n')
                for pat in patList:
                    f.write('Reference Rule:'+ '\t' + pat + '\n')
                    f.write('-'*100 + '\n')
            f.flush()
        f.close()

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
                PatDic[t].append("&"+pat[i]+"."+pat[i+1]) #此处对于(CCC.CCC)一个分子内含有两个片段的rs 做一个标记，方便后期使用识别
            else:                                       #当num1<num2时,前一个s 一定是num1>num2,已经被添加进去
                pass
    ridPatRsPs[rid]["patrs"] = ridPatRs
    ridPatRsPs[rid]["patps"] = ridPatPs    

def RidPatRsPs(Ec = False):
    
    '''#构建一个{rid1:{patrs:list,patps:list},rid2:{patrs:list,patps:list},}的字典'''

    ridPatRsPs = defaultdict() 

    Rid_EcAndPattern = json.load(open("./Rid_EcAndPattern.json")) #加载反应ID对应的反应Patterns，Smirks,Ec等信息
    
    for rid,item in Rid_EcAndPattern.items(): #将patterns拆分为反应物和产物两部分的列表
        
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

    #构建需要replace的带电原子类型与其对应的中性原子的pair对
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
    ('[$([O-]=C)]','O'),
    )
    return [(Chem.MolFromSmarts(x),Chem.MolFromSmiles(y,False)) for x,y in patts]

_reactions=None

def NeutraliseCharges(smiles, reactions=None):
    global _reactions
    if reactions is None: #默认不输入
        if _reactions is None:
            _reactions=InitialiseNeutralisationReactions()
        reactions=_reactions #reactions引入带电原子类型与其对应的中性原子的pair对表单函数
    mol = Chem.MolFromSmiles(smiles)

    #判断Chem.MolFromSmiles是否成功读取smile，如没有换成smarts读取
    #例如“NC(=O)c1cccn(c1)C1OC(COP(=O)([O-])OP(=O)([O-])OCC2OC(C(O)C2O)n2cnc3c(N)ncnc32)C(O)C1O”
    
    if not mol:
        mol = Chem.MolFromSmarts(smiles)
    if not mol:
        return None

    else:
        replaced = False
        for i,(reactant, product) in enumerate(reactions):
            while mol.HasSubstructMatch(reactant): #一直循环至不含该带电原子
                replaced = True
                rms = AllChem.ReplaceSubstructs(mol, reactant, product) #ReplaceSubstructures可选择Replacement = True（默认为False）一步替换所有
                mol = rms[0] #rms是一个tuple,内含多个重复的mol，原因不明
        if replaced:
            return (Chem.MolToSmiles(mol,True), True) 
        else:
            return (Chem.MolToSmiles(mol,True), False) # modeified@20161121

def similar(smiles1,smiles2):
    smiles1_mol = mfsmi(smiles1) ####modified
    if not bool(smiles1_mol):####modified
        smiles1_mol = mfsma(smiles1)####modified
    
    smiles2_mol = mfsmi(smiles2)####modified
    if not bool(smiles2_mol):####modified
        smiles2_mol = mfsma(smiles2)####modified
    
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

def Testsmi(sm='',rxn =''):
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

def smirks2rxnName(smirks): 

    newCPSmiId = json.load(open('./newCPSmiId.json'))

    try:
        StandSmiName = json.load(open('./StandSmi_NameId.json')) # modeified@20161121
    except:
        ChEBISmiName = json.load(open('./ChEBISmi_NameId.json')) # modeified@20161121
        StandSmiName = dict()
                                            # modeified@20161121
        for smi,name in ChEBISmiName.items(): 
            # try:
            #     smiMol = mfsmi(smi)
            # except:
            #     smiMol = mfsma(smi)

            # if smiMol:
            standSmi = NeutraliseCharges(smi)
            if not standSmi:
                standSmi = smi #如果标准化不成功，则还是用原来的smiles
            else:
                standSmi = standSmi[0]
            StandSmiName[standSmi] = name
            # else:
            #     continue
        with open('./StandSmi_NameId.json','w') as f1:
            json.dump(StandSmiName,f1,indent = 2)

    rs = [i.strip() for i in smirks.split('>>')[0].split('.')]
    ps = [j.strip() for j in smirks.split('>>')[1].split('.')]

   
    _rs = list(set(rs))
    _ps = list(set(ps))
    coefRs = list();coefRsName = list();coefRsId = list()
    coefPs = list();coefPsName = list();coefPsId = list()

    index = {'0':rs,'1':ps}
    index1 = {'0':coefRs,'1':coefPs}
    index2 = {'0':coefRsName,'1':coefPsName}
    index3 = {'0':coefRsId,'1':coefPsId}

    for m,it in enumerate([_rs,_ps]):

        for s in it:

            num = index[str(m)].count(s)
            
            # if s == '[OH2]': # modeified@20161121
            #     s = 'O'

            stand_s = NeutraliseCharges(s)

            if not stand_s: #如果标准化不成功，则还是用原来的smiles
                stand_s = s
            else:
                stand_s = stand_s[0] # modeified@20161121
            
            try:
                sName = StandSmiName[stand_s]["ChEBI Name"]
                sId = StandSmiName[stand_s]["ChEBI ID"][0] # modeified@20161121
            except:
                if stand_s in newCPSmiId.keys():
                    sName = newCPSmiId[stand_s]
                    sId = newCPSmiId[stand_s]
                else:
                    maxid = max([int(i.split('W')[1]) for i in newCPSmiId.values()])
                    sName = 'W' + str(maxid+1)
                    newCPSmiId[stand_s] = sName
                    sId = sName    # modeified@20161121
            
            if num != 1:
                index1[str(m)].append(str(num) + ' ' + stand_s)
                index2[str(m)].append(str(num) + ' ' + str(sName))
                index3[str(m)].append(str(num) + ' ' + str(sId)) # modeified@20161121

            else:
                index1[str(m)].append(stand_s)
                index2[str(m)].append(str(sName))
                index3[str(m)].append(str(sId)) # modeified@20161121


    coefSmirks =' + '.join(coefRs)+ ' >> ' + ' + '.join(coefPs)
    coefrxnName =' + '.join(coefRsName) + ' >> ' +' + '.join(coefPsName)
    coefrxnId =' + '.join(coefRsId)+ ' >> ' + ' + '.join(coefPsId) # modeified@20161121

    with open('./newCPSmiId.json','w') as fn:
        json.dump(newCPSmiId,fn,indent = 2)
    return (coefSmirks,coefrxnName,coefrxnId)


def main():
    start = time.strftime(r'%Y-%m-%d-%H-%M-%S',time.localtime(time.time()))
    queryRsSmi ="CSC"
    queryPsSmi =''#"c1cc(C)ccc1C(C)C" #raw_input('''Product Input:''')
    BioReactor(queryRsSmi,queryPsSmi,Ec = False,CALSIMI=True,Draw=False)
    end = time.strftime(r'%Y-%m-%d-%H-%M-%S',time.localtime(time.time()))
    
    print start
    print end

        
if __name__ == "__main__":
    main()
   
    




