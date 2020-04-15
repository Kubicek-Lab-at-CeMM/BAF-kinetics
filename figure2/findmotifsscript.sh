#!/bin/bash

# homer software needs to be installed for this (http://homer.ucsd.edu/homer/)

findMotifsGenome.pl ../data/cluster1_tc-FC-based-for-homer.bed hg38 ../data/cluster1_tc/ -size given -bg ../data/all-consensus-FC-based-for-homer.bed

findMotifsGenome.pl ../data/cluster2_tc-FC-based-for-homer.bed hg38 ../data/cluster2_tc/ -size given -bg ../data/all-consensus-FC-based-for-homer.bed

findMotifsGenome.pl ../data/cluster3_tc-FC-based-for-homer.bed hg38 ../data/cluster3_tc/ -size given -bg ../data/all-consensus-FC-based-for-homer.bed

findMotifsGenome.pl ../data/cluster3_KO-FC-based-for-homer.bed hg38 ../data/cluster3_KO/ -size given -bg ../data/all-consensus-FC-based-for-homer.bed

findMotifsGenome.pl ../data/cluster1_cluster2_KO-FC-based-for-homer.bed hg38 ../data/cluster1_cluster2_KO/ -size given -bg ../data/all-consensus-FC-based-for-homer.bed

