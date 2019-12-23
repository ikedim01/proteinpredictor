
import collections
import csv
import itertools
import pickle
import random
import re

# CODE TO PARSE .DAT FILES (Swissprot flat file format):
def iterDat(fPath) :
    with open(fPath,'r') as f :
        for (isEnd,it) in itertools.groupby(f,lambda ln : ln.strip()=='//') :
            if not isEnd :
                yield(list(it))
def scanDat(datFPath,fn,**kwargs) :
    res = []
    for datGroup in iterDat(datFPath) :
        item = fn(datGroup,**kwargs)
        if isinstance(item,list) :
            res.extend(item)
        elif item is not None :
            res.append(item)
    return res
def countDat(datFPath,fn,**kwargs) :
    return collections.Counter(scanDat(datFPath,fn,**kwargs))
def filterDatForLMSet(datGroup, skipUndet=True, skipPyl=True,
            minLen=50, maxLen=400, maxPE=None,
            requireInOS=None,
            requireInName='', excludeStrs=[]) :
    seq = getDatSeq(datGroup)[1]
    # print(seq)
    if ((skipPyl and 'O' in seq)
        or (skipUndet and any(c in seq for c in 'BJZX*-'))
        or (minLen is not None and len(seq)<minLen)
        or (maxLen is not None and len(seq)>maxLen)
        or (maxPE is not None and getDatPE(datGroup)>maxPE)
        or (requireInOS is not None
            and requireInOS not in getDatLns(datGroup,'OS')[0].lower())) :
        return None
    name = getDatName(datGroup)
    lName = name.lower()
    excludeStrs = listIfStr(excludeStrs)
    if ((requireInName not in lName)
        or any (excludeStr in lName for excludeStr in excludeStrs)) :
        return None
    return (getDatAC(datGroup),seq)
def getLMDatasetFromDatFile(datFPath, **kwargs) :
    return dict(scanDat(datFPath,filterDatForLMSet,**kwargs))
def makeLMDataset(datFPath, defKwArgs, outFPath=None, clusterFPath=None, clusterLen=30,
                    **kwOverrideArgs) :
    kwargs = dict(defKwArgs)
    kwargs.update(kwOverrideArgs)
    seqM = getLMDatasetFromDatFile(datFPath,**kwargs)
    allLetters = ''.join(allCs(seqM.values()))
    print('total length {:,} in {:,} sequences; {} letters: {}'
          .format(sum(len(ss) for ss in seqM.values()),len(seqM.values()),
                  len(allLetters),allLetters))
    if outFPath is not None :
        writeSeqMCsv(outFPath, seqM)
    if clusterFPath is None :
        return seqM
    clusters = clusterRepeats(seqM,clusterLen)
    if clusterFPath != True :
        with open(clusterFPath,'wb') as f :
            pickle.dump(clusters,f)
    return seqM,clusters
def getDatSeq(datGroup) :
    seqL = []
    for ln in reversed(datGroup) :
        ln = ln.strip()
        if ln.startswith('SQ ') :
            break
        seqL.append(ln.replace(' ',''))
    seq = ''.join(reversed(seqL))
    return ln,seq,len(seq)
def getDatLns(datGroup,code) :
    res = []
    for ln in datGroup :
        ln = ln.strip()
        if ln.startswith(code+' ') :
            res.append(ln[len(code):].lstrip())
    return res
def listIfStr(lOrStr) :
    if isinstance(lOrStr,str) :
        return [lOrStr]
    return lOrStr
numPatStr = r'\d+(?:\.\d*)?'
optimPHPat = re.compile(
    r'optimum\s+ph\s+is\s+(?:between\s+)?('+numPatStr+r')(?:\s*(?:-|and)\s*('+numPatStr+'))?',
    re.IGNORECASE)
def phStrToFloat(phStr) :
    res = float(phStr)
    if res == 70.0 :
        res = 7.0
    return res
def printGroupSummary(msg,datGroup,lns) :
    print(msg)
    for ln in datGroup[:2] :
        print(ln,end='')
    for ln in listIfStr(lns) :
        print(ln)
def getDatGoTerm(datGroup,goTerm) :
    drLns = getDatLns(datGroup,'DR')
    for drLn in drLns :
        if goTerm in drLn :
            return 'pos'
    return 'neg'
def getDatAtpBinding(datGroup) : return getDatGoTerm(datGroup,'GO:0005524;')
def getDatGtpBinding(datGroup) : return getDatGoTerm(datGroup,'GO:0005525;')
def getDatMetalBinding(datGroup) : return getDatGoTerm(datGroup,'GO:0046872;')
def getDatOptimPHRange(datGroup) :
    phRange = getDatOptimPH(datGroup)
    if phRange is None :
        return None
    if phRange < 5.51 :
        return 'acidic'
    if phRange < 8.49 :
        return 'neutral'
    return 'basic'
def getDatOptimPH(datGroup) :
    ccLns = getDatLns(datGroup,'CC')
    for ccLn in ccLns :
        m = optimPHPat.match(ccLn)
        if m is not None :
            res = phStrToFloat(m.group(1))
            if m.group(2) :
                res2 = phStrToFloat(m.group(2))
                if res<res2<=14.0 :
                    res = (res + res2) / 2.0
                else :
                    printGroupSummary("couldn't parse pH value:",datGroup,ccLn)
                    return None
            if (not 2.0<=res<=12.0) :
                printGroupSummary('extreme pH value {:.2f}:'.format(res),datGroup,ccLn)
            return res
def getDatPE(datGroup) :
    peLns = getDatLns(datGroup,'PE')
    return int(peLns[0][0])
def getDatAC(datGroup) :
    acLns = getDatLns(datGroup,'AC')
    return acLns[0].split(';')[0]
def getDatGO(datGroup) :
    drLns = getDatLns(datGroup,'DR')
    return [[drLn for drLn in drLns if drLn.startswith('GO;')]]
def getDatName(datGroup) :
    deLns = getDatLns(datGroup,'DE')
    name = deLns[0].split('=')[1].rstrip(';')
    flagSet = set()
    for deLn in deLns[1:] :
        if deLn.startswith('Flags') :
            flagSet.update(flag.rstrip(';') for flag in deLn.split()[1:])
    if len(flagSet) > 0 :
        name += ' (' + '; '.join(sorted(flagSet)) +')'
    return name

def readCsv(fPath) :
    with open(fPath,'r') as f :
        return list(csv.reader(f))
def writeCsv(fPath,rows,headers=None) :
    with open(fPath,'w',newline='') as f :
        csvw = csv.writer(f)
        if headers is not None :
            csvw.writerow(headers)
        csvw.writerows(rows)
def writeSeqMCsv(fPath,m) :
    writeCsv(fPath,
            [[k,m[k]] for k in sorted(m.keys())],
            ['identifier','sequence'])
def readSeqMCsv(fPath) :
    return dict(readCsv(fPath)[1:])

def nSubstrs(s,m) :
    return len(s)+1-m
def itSubstrs(s,m) :
    return (s[i:i+m] for i in range(nSubstrs(s,m)))
def totNSubstrs(ss,m) :
    return sum(nSubstrs(s,m) for s in ss)
def getUniqueSubstrs(ss,m) :
    return set(itertools.chain.from_iterable(itSubstrs(s,m) for s in ss))
def totNUniqueSubstrs(ss,m) :
    return len(getUniqueSubstrs(ss,m))
def findLongestRepeat(s,sourceStrs) :
    for m in range(2,30) :
        sSet = set(itSubstrs(s,m))
        sourceSet = set(itertools.chain.from_iterable(itSubstrs(sourceS,m) for sourceS in sourceStrs))
        if len(sSet&sourceSet) == 0 :
            return m-1
    return m
def clusterRepeats(identToSeq,m) :
    subseqRep = {}
    edges = collections.defaultdict(set)
    for ident,seq in identToSeq.items() :
        for subseq in itSubstrs(seq,m) :
            if subseq not in subseqRep :
                subseqRep[subseq] = ident
            else :
                edges[ident].add(subseqRep[subseq])
                edges[subseqRep[subseq]].add(ident)
    print ('length-{} subsequences: total {:,}, unique {:,}'.format(
            m, totNSubstrs(identToSeq.values(),m), len(subseqRep)))
    subseqRep = None
    identToComponentNo = {}
    components = []
    for ident in identToSeq :
        if ident in identToComponentNo :
            continue
        curComponent = set()
        curComponentNo = len(components)
        components.append(curComponent)
        stack = [ident]
        while len(stack) > 0 :
            ident = stack.pop()
            if ident in identToComponentNo :
                continue
            curComponent.add(ident)
            identToComponentNo[ident] = curComponentNo
            stack.extend(identNeighbor for identNeighbor in edges[ident]
                         if identNeighbor not in identToComponentNo)
    print('{:,} clusters; identifier counts: {:,}, {:,}'.format(
            len(components),len(identToSeq),len(identToComponentNo)))
    return components,identToComponentNo

def dupPropVal(nCopies,propVal) :
    return lambda prop,datGroup : nCopies if prop==propVal else 1
def dupPos(nCopies) : return dupPropVal(nCopies,'pos')
def makeClasDataset(datFPath, clusterFPath, propFn, propKwArgs={},
                    outFPath=None, testFPath=None,
                    split_pct=0.25, test_split_pct=0.4,
                    nCopiesFn=lambda prop,datGroup : 1, seed=3141) :
    with open(clusterFPath,'rb') as f :
        _,identToComponentNo = pickle.load(f)
    if seed is not None :
        random.seed(seed)
    trainAndValSet,testSet = [],[]
    componentNoToSplit = {}
    nTrainExamples = nValidExamples = nTestExamples = nSkipped = nPropClusters = 0
    for datGroup in iterDat(datFPath) :
        ident = getDatAC(datGroup)
        if ident not in identToComponentNo :
            continue
        componentNo = identToComponentNo[ident]
        prop = propFn(datGroup,**propKwArgs)
        if prop is None :
            continue
        example = (ident,getDatName(datGroup),getDatSeq(datGroup)[1],prop)
        if componentNo not in componentNoToSplit :
            nPropClusters += 1
            if random.random() >= split_pct :
                componentNoToSplit[componentNo] = 'train'
            elif random.random() < test_split_pct :
                componentNoToSplit[componentNo] = 'test'
                nTestExamples += 1
                testSet.append(example)
                continue
            else :
                componentNoToSplit[componentNo] = 'valid'
                nValidExamples += 1
                trainAndValSet.append(example + (True,))
                continue
        if componentNoToSplit[componentNo] == 'train' :
            nTrainExamples += 1
            trainAndValSet.extend(nCopiesFn(prop,datGroup)*[example + (False,)])
        else :
            nSkipped += 1
    print('train:',nTrainExamples,'valid:',nValidExamples,'test:',nTestExamples,
            'skipped:',nSkipped,'clusters:',nPropClusters)
    if outFPath is not None :
        writeCsv(outFPath,trainAndValSet,['identifier','name','sequence','label','is_valid'])
    if testFPath is not None :
        writeCsv(testFPath,testSet,['identifier','name','sequence','label'])
    return trainAndValSet,testSet

def allCs(strs) :
    return sorted(set(itertools.chain.from_iterable(strs)))
def cCounts(strs) :
    res = collections.defaultdict(lambda : 0)
    for c in itertools.chain.from_iterable(strs) :
        res[c] += 1
    return res

def markov(strs,n,allC=None) :
    mmod = collections.defaultdict(lambda : collections.defaultdict(lambda : 0))
    for s in strs :
        for i in range(n,len(s)) :
            mmod[s[i-n:i]][s[i]] += 1
    res = {}
    for k,mm in mmod.items() :
        res[k] = max((mmk for mmk in mm), key=lambda x : mm[x])
    if allC is None :
        allC = allCs(strs)
    for tup in itertools.product(allC,repeat=n) :
        k = ''.join(tup)
        if k not in res :
            res[k] = tup[-1]
    return res
def markovAcc(mm,strs) :
    tot = nCorrect = 0
    n = len(mm.keys()[0])
    for s in strs :
        tot += len(s)-n
        for i in range(n,len(s)) :
            if mm[s[i-n:i]] == s[i] :
                nCorrect += 1
    return nCorrect,tot,nCorrect/float(tot)

# OLD CODE:
def randSplit(it,p,seed=None) :
    if seed is not None :
        random.seed(seed)
    l = sorted(it)
    random.shuffle(l)
    res = ([],[])
    for i in it :
        res[0 if random.random()<p else 1].append(i)
    return res
def addSeq(seqM, curName, curSeq,
            skipUndet=True, skipPyl=True,
            minLen=50, maxLen=400, requireInName='', excludeStrs=[]) :
    curSeq = "".join(curSeq)
    if curName is None :
        if len(curSeq) > 0 :
            print('missing name for',len(curSeq),' - length sequence')
        return
    excludeStrs = listIfStr(excludeStrs)
    lCurName = curName.lower()
    if ((skipUndet and any(c in curSeq for c in 'BJZX*-'))
            or (skipPyl and 'O' in curSeq)
            or (minLen is not None and len(curSeq)<minLen)
            or (maxLen is not None and len(curSeq)>maxLen)
            or (not requireInName in lCurName)
            or any(excludeStr in lCurName for excludeStr in excludeStrs)) :
        return
    if curName in seqM :
        print('repeated name',curName)
        print('seq1',seqM[curName])
        print('seq2',curSeq)
    if len(curSeq) == 0 :
        print('empty sequence for',curName)
    seqM[curName] = curSeq
def readpseq(fPath, **kwargs) :
    res = {}
    with open(fPath,'r') as f :
        curName,curSeq = None,[]
        for ln in f :
            ln = ln.strip()
            if ln.startswith('>') :
                addSeq(res,curName,curSeq,**kwargs)
                curName,curSeq = ln,[]
            else :
                curSeq.append(ln)
        addSeq(res,curName,curSeq,**kwargs)
    return res
