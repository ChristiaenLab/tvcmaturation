URL = http://ghost.zool.kyoto-u.ac.jp/datas/
GFF = HT.Gene.gff
ZIP = HT.KYGene.gff3.zip

install: peakToGene dirfns BSgenome.Crobusta.HT.KY CrobustaTFs tfenrichr moreComplexHeatmap
	make -C peakToGene
	make -C dirfns
	make -C BSgenome.Crobusta.HT.KY
	make -C CrobustaTFs
	make -C tfenrichr
	make -C moreComplexHeatmap

$(GFF):
	wget -U firefox $(URL)$(ZIP)
	unzip -o $(ZIP)

DAmotifs:
	gunzip DAmotifs.csv.gz

peakToGene:
	git clone https://github.com/kewiechecki/peakToGene

dirfns:
	git clone https://github.com/kewiechecki/dirfns

BSgenome.Crobusta.HT.KY:
	git clone https://github.com/kewiechecki/BSgenome.Crobusta.HT.KY

CrobustaTFs:
	git clone https://github.com/kewiechecki/CrobustaTFs

tfenrichr:
	git clone https://github.com/kewiechecki/tfenrichr

moreComplexHeatmap:
	git clone https://github.com/kewiechecki/moreComplexHeatmap

clean:
	rm -rf peakToGene dirfns BSgenome.Crobusta.HT.KY CrobustaTFs tfenrichr moreComplexHeatmap
	rm -f $(ZIP)
