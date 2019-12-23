# Predicting protein properties (Gene Ontology/GO) using ULMFiT

This project is an attempt to use ULMFiT to predict functional properties of a protein from its one-dimensional amino acid sequence. First results seem promising - using the Swissprot database as a source, we predicted the ATP-binding and GTP-binding properties with high accuracy.

Since we're using ULMFit on proteins, the first two things to cover are:
**What is ULMFiT?** and, **What are proteins?**

## What is ULMFiT?
An algorithm with breakthrough results on NLP; a clever application of transfer learning - [the original paper](https://arxiv.org/abs/1801.06146).

I'll give a quick walkthrough of applying ULMFiT to a particular NLP problem, as presented in fastai 2019 part1. For more detail, see the IMDB notebook from lesson 3.

The problem: **sentiment analysis on IMDB movie reviews**. A classification problem - given a movie review written in English, is the review positive or negative? Training set of 25K labeled reviews.

**The difficulty:** This classification is quite a complex function of the input text, so we probably need a powerful model like a deep RNN to get good results. But, we don't have enough data to train a powerful model! And this is typical of NLP tasks (labeled training data may be hard to get).

**The clever idea:** use a form of transfer learning.
- Train a language model on English text (given a text segment, predict the word that comes next). We have plenty of plain text data to train a powerful model like a multilayer RNN for this! Note we'd usually use a pretrained model for this step, since these models can take a long time to train. For example, the fastai IMDB notebook from lesson 3 uses a model with the AWD_LSTM architecture that was pretrained on a subset of wikipedia called wikitext-103 (size approximately 100M words/tokens).
- Fine-tune the language model for the kind of text in our problem (movie reviews). Note in this step we can use all of our problem's data, including the training, validation and test sets!	
- Use the internal state of this language model as an encoder; train the classifier based on the encoded values of the positive/negative examples in our training set.

**More intuition about how ULMFiT works:** To do a really good job of predicting the next word, our language model has to have a pretty deep knowledge of English (or whatever language it was trained on). This knowledge is encoded into the internal state of the model, and so can be used to do other tasks like classifying text.

## What are proteins?

**The basics:** There are assorted complications on top of this, but the basic picture is:
- Chromosomes are 1D DNA sequences of 4 possible "tokens" (nucleotides).
- Chromosomes contain genes, which are **transcribed** into 1D mRNA sequences, again of 4 possible tokens (nucleotides).
- These mRNA transcripts are then **translated** into proteins, which start out as 1D sequences of 20 possible tokens (amino acids - each mRNA triplet specifies an amino acid, with some redundancy).
- These proteins then fold into a 3D functional shape that enable them to do their thing! This last part is a very complex process.

**Why is this interesting?** Of course, proteins ultimately determine most of what goes on in biological organisms, and so are crucial to figuring out what causes diseases, to developing drugs, and so on.

But also: proteins are incredible molecular machines! They are basically working nanotechnology. For example, consider enzymes (proteins that catalyze chemical reactions). Enzymes are capable of chemistry that is far beyond that designed by humans.
Chemistry designed by humans:
- often works at extreme temperature/pressure; we have to bash the molecules together at very high speed to get a few lucky collisions that give us our wanted reaction.
- is relatively nonspecific - will often produce a high proportion of unwanted "side reactions" that must then be purified out.
- can generally only do one reaction at a time.
Enzyme chemistry:
- generally works in very mild conditions (body temperature, neutral pH)
- is very specific - side reactions are often only a few percent or less; the enzyme is basically grabbing the reactant molecules and holding them in exactly the right position to react!
- many reactions are proceeding at the same time in a very small space!

**So, being able to predict a protein's functions and properties** from its 1D sequence is very interesting! And, as mentioned, the function mapping this 1D sequence to its functional shape (and therefore to its functions) is quite complex. And, experiments on protein function are very difficult and expensive compared to sequencing - this starts to sound like the NLP situation! So, I thought it would be worth trying a version of ULMFiT on this. I did some searching and found one existing project applying ULMFiT to genomic data at https://github.com/kheyer/Genomic-ULMFiT - however, the experiments in this project used ULMFiT on DNA and RNA sequence data, so I decided to try it on protein sequence data, with the objective of usefully predicting some protein functional properties.

## Experimental methods/Results

**Building the data sets:** I used the Swiss-Prot database. This apparently aims to be a comprehensive and up-to-date database of protein candidates that have been reviewed by human experts. For each protein, we have (of course) the amino acid sequence, but we also have lots of other helpful annotations, including various functional information, and the level of evidence for the protein's existence (this is on a scale of 1 to 5; 1=protein level, 2=transcript, 3=homology, 4=predicted, 5=uncertain). The latest version of Swiss-Prot can be downloaded from https://www.uniprot.org/downloads - I used the flat-file (.dat) format as of Dec. 3, 2019; the format of this file is given at https://web.expasy.org/docs/userman.html 

Note: there's also a much larger database of protein sequences at the same site, but as I understand it almost all of these are predicted automatically from genomic (DNA) data without expert review. Since there seemed to be enough data in Swiss-Prot (see below) and these proteins seem more likely to actually exist, I decided to restrict training to sequences in Swiss-Prot.

I restricted to protein sequences with evidence level <= 3. I also restricted to sequences with lengths between 40 and 500, since I thought the "rules" for very short or very long sequences might be different. Finally, I restricted to sequences with no unknown residues, and to sequences that weren't flagged as fragments. This resulted in a total of about 110M amino acids in 424K sequences; according to Jeremy Howard, 100M tokens is generally enough to build a good language model.

For details, see the [dataproc notebook](nb/dataproc.ipynb).

**Building the language model:** I used single amino acids as tokens, and used fastai to train a language model with an AWD_LSTM architecture with default parameters using fit_one_cycle with a learning rate of 3e-3 for 40 epochs (other hyperparameters: moms=(0.8,0.7), drop_mult=0.2, random split of 0.1 for the validation set).

**Language model results:** The model took about 2 days to train on my desktop and reached an accuracy of 53.2% on the validation set. I haven't yet tried to fine-tune the model as in the full ULMFiT algorithm. I also did the obligatory cool PCA graph on the amino acid embedding values in the model; this seemed to do a pretty good job of clustering the amino acids according to their chemical structure.

For details, see the [training notebook](nb/ulmptrain.ipynb).

**Testing property prediction:** I used the GO (Gene Ontology) annotations in the Swiss-Prot database to generate a couple of binary classification problems. Specifically, I picked the ATP-binding and GTP-binding terms (GO:0005524 and GO:0005525, respectively), and attempted to predict whether these were present based on the protein sequence. To try to avoid a potential problem with repeated proteins, I clustered all sequences that shared any common 30-aa sequence, and randomly assigned each cluster to only one of the training, validation, and test sets; also, if a cluster was assigned to the validation or test sets I only used one sequence from that cluster. I randomly split clusters 75:15:10 into training, validation, and test.

**Property prediction results:**
For ATP binding, the classifier reached 98.7% accuracy after 5 epochs of training. Precision was .945, recall .918, and F1 score .932. (all with a classification threshold of 0.5).
For GTP binding,  the classifier reached 99.8% accuracy after 5 epochs of training. Precision was .989, recall .893, and F1 score .938 (all with a classification threshold of 0.5).
I haven't yet tried to fine-tune the classifiers by unfreezing shallower layers.

For details, see the [classification notebook](nb/ulmpclas.ipynb).

## Tentative conclusion:

While much remains to be done, initial results are promising.

**Stuff to do:**
- Improve the language model (try different hyperparameters, more epochs).
- Try fine-tuning as in the full ULMFiT (ex. fine tune the language model by organism; fine-tune classifiers by unfreezing shallower layers).
- Try more properties.
- Use evidence tags for the GO terms (requires accessing other databases?).
