
from fastai import *
from fastai.text import *
import random

def randomSeedForTraining(seed) :
    random.seed(seed)
    torch.manual_seed(seed)
    torch.backends.cudnn.deterministic = True
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    # torch.cuda.manual_seed(seed_value)
    # torch.backends.cudnn.benchmark = False

def getModelFNameFromHyperPars(pref,*hyperParsList) :
    res = pref
    for hyperPars in hyperParsList :
        if isinstance(hyperPars,tuple) or isinstance(hyperPars,list) :
            hyperPars, defHyperPars = hyperPars
        else :
            defHyperPars = {}
        for k in sorted(hyperPars.keys()) :
            if hyperPars[k] != defHyperPars.get(k) :
                res += '_{}_{}'.format(k,hyperPars[k]).replace('.','_')
    return res
def getModelVersions(modelFName) :
    i = 0
    res = [modelFName]
    while (modelDirPath/(res[-1]+'.pth')).exists() :
        i += 1
        res.append(modelFName + '_' + str(i))
    return res
def loadModelVersion(learn,modelFName,i) :
    modelVersions = getModelVersions(modelFName)[:-1]
    if i >= 0 :
        versionAvailable = i < len(modelVersions)
    else :
        versionAvailable = -i <= len(modelVersions)
    if versionAvailable :
        print('Loading',modelVersions[i])
        learn.load(modelVersions[i])
    else :
        print('Version not available to load!')
def loadLatestModelVersion(learn,modelFName) :
    loadModelVersion(learn,modelFName,-1)
def removeLatestModelVersion(learn,modelFName) :
    modelVersions = getModelVersions(modelFName)
    if len(modelVersions) >= 2 :
        print('Removing',modelVersions[-2])
        (modelDirPath/(modelVersions[-2]+'.pth')).unlink()
    else :
        print('No version available to load!')
def saveNextModelVersion(learn,modelFName) :
    modelVersions = getModelVersions(modelFName)
    print('Saving',modelVersions[-1])
    learn.save(modelVersions[-1])

class SeqTokenizer(BaseTokenizer):
    def __init__(self, lang='en', ngram=1, stride=1):
        self.lang, self.ngram, self.stride = lang, ngram, stride
    def tokenizer(self, t):
        res = []
        pos = 0
        while pos+self.ngram <= len(t) :
            if t.startswith(BOS,pos) :
                res.append(BOS)
                pos += len(BOS)+1
            else :
                res.append(t[pos:pos+self.ngram])
                pos += self.stride
        return res
    def add_special_cases(self, toks):
        pass

aaLetters = 'ACDEFGHIKLMNPQRSTVWY'
aaLettersSet = set(aaLetters)
tokenizer = Tokenizer(SeqTokenizer,pre_rules=[],post_rules=[])
vocab = Vocab.create(aaLetters,max_vocab=60000,min_freq=0)
processor=[TokenizeProcessor(tokenizer=tokenizer,mark_fields=False),
           NumericalizeProcessor(vocab=vocab)]
