#!/bin/bash

# homer software needs to be installed for this (http://homer.ucsd.edu/homer/)

findMotifsGenome.pl cluster1_tc-FC-based-for-homer.bed hg38 cluster1_tc/ -size given -bg all-consensus-FC-based-for-homer.bed

findMotifsGenome.pl cluster2_tc-FC-based-for-homer.bed hg38 cluster2_tc/ -size given -bg all-consensus-FC-based-for-homer.bed

findMotifsGenome.pl cluster3_tc-FC-based-for-homer.bed hg38 cluster3_tc/ -size given -bg all-consensus-FC-based-for-homer.bed

findMotifsGenome.pl cluster3_KO-FC-based-for-homer.bed hg38 cluster3_KO/ -size given -bg all-consensus-FC-based-for-homer.bed

findMotifsGenome.pl cluster1_cluster2_KO-FC-based-for-homer.bed hg38 cluster1_cluster2_KO/ -size given -bg all-consensus-FC-based-for-homer.bed

