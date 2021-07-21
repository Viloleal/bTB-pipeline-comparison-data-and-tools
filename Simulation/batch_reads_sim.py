#!/usr/bin/env python

import glob
import re
import subprocess


fasta = glob.glob("*fasta")
prefixes = [re.sub(r'.fasta','',file) for file in fasta]
for (files,prefix) in zip(fasta,prefixes):
	print ("fasta: ", files, "; prefix: ", prefix)
	listToStr1 = ''.join(map(str, files))
	listToStr2 = ''.join(map(str, prefix))
	file = open("fasta_prefixes.txt", "a")
	file.write("List of files \n fasta \n %s \n Prefix \n %s \n ------------------- \n" % (listToStr1,listToStr2))
	file.close()
	subprocess.call(['bash', '-c', 'java -jar ~/bin/ArtificialFastqGenerator/ArtificialFastqGenerator.jar -O ./%s -R %s -S \'>\' -RL 250 -TLM 650 -TLSD 60 -CMP 250 -CSD 0.1 -N 100000 SE true -URQS true' % (prefix,files)])
