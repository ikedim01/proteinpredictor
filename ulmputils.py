
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

aaLetters = 'ACDEFGHIKLMNPQRSTUVWY'
tokenizer = Tokenizer(SeqTokenizer,pre_rules=[],post_rules=[])
vocab = Vocab.create(aaLetters,max_vocab=60000,min_freq=0)
processor=[TokenizeProcessor(tokenizer=tokenizer,mark_fields=False),
           NumericalizeProcessor(vocab=vocab)]
