#!/usr/bin/bash

head -1 dataDir/0106/G-6-4A22.1.cu.bo.hg19.coToBa.coToBe.unSo.coBeToSiFr.slBeb6.coToFiRa10.soBe.coBeToFa.diFr.txt >dataDir/mergedDinuc.txt
tail -q --lines=+2 dataDir/0106/*cu.bo.hg19.coToBa.coToBe.unSo.coBeToSiFr.slBeb6.coToFiRa10.soBe.coBeToFa.diFr.txt >>dataDir/mergedDinuc.txt