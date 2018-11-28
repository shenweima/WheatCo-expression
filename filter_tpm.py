#!/usr/bin/python
# -*- coding: utf-8 -*-
__author__ = 'shengwei ma'
__author_email__ = 'shengweima@icloud.com'

import math

with open('PRJEB25639_tpm.tsv', 'r') as f:
    for line in f:
        max = 0
        lin = line.strip().split('\t')
        if line.startswith('GeneID'):
            print('\t'.join(lin))
        else:
            for i in range(1,len(lin)):
                if float(lin[i]) >= 5:
                    max += 1
                if float(lin[i]) < 0.5:
                    lin[i] = 0
                else:
                    lin[i] = math.log(float(lin[i]) + 1, 2)
                
        if max >=3:
            print('\t'.join('%s' % t for t in lin))